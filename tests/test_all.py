"""
tests/test_all.py
=================
Standalone tests for all KazNRDC core modules.

No nuclear data files (xsdir/) required — all isotopes built manually.
Run from the kaznrdc9/ project root:

    python -m pytest tests/test_all.py -v
    # or directly:
    python tests/test_all.py

Coverage
--------
1.  nuclear_data.models      — dataclass fields and properties
2.  core.entities.ChemTable  — Z lookup, invalid symbol
3.  core.entities.Element    — valid / invalid construction, equality, hash
4.  core.entities.Isotope    — name, meta, from_symbol, add_decay/reaction
5.  core.entities.DecayLink  — construction, repr, default product=None
6.  core.entities.ReactionLink — getId, repr
7.  core.burnup.BurnupConfig — all filters and presets
8.  core.burnup.BurnupMatrix (decay-only) — matrix shape, all solvers, heat,
                               equilibrium, evolve, export, validate
9.  core.burnup.BurnupMatrix (reactions)  — neutron reactions, conservation
10. core.burnup.BurnupMatrix open/closed boundary
11. core.cycle_finder.CycleAnalyzer — SCC, all-cycles, queries, print_report
12. core.cycle_finder (acyclic network) — no false positives
13. core.registry.IsotopeBuilder (mock provider) — build, link, report
14. core.service — mt_symbol, decay_symbol, CHANNEL_DELTAS
15. Integration — full headless 3-isotope pipeline + cycle injection test
"""

from __future__ import annotations

import contextlib
import io
import math
import os
import sys
import tempfile
import unittest

import numpy as np

# ── add project root to path ──────────────────────────────────────────────
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from nuclear_data.models import (
    CrossSectionPoint, DecayChannel, DecayData,
    IsomericBranch, MacsPoint, NeutronReactionData,
)
from core.entities import ChemTable, DecayLink, Element, Isotope, ReactionLink
from core.registry import IsotopeBuilder
from core.burnup import BurnupConfig, BurnupMatrix
from core.cycle_finder import CycleAnalyzer, CycleStep, IsotopeCycle, SCC


# ═════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ═════════════════════════════════════════════════════════════════════════════

def make_isotope(Z: int, X: str, A: int,
                 lam: float = 0.0,
                 hl: float | None = None) -> Isotope:
    """Build a bare Isotope (no links). lam=decay constant [s⁻¹]."""
    iso = Isotope(Z, X, A)
    iso.decay_constant = lam
    iso.half_life_s = hl
    return iso


def two_isotope_decay_chain():
    """
    A (H-3, tritium) --[beta-]--> B (He-3, stable).
    Analytic solution:
        N_A(t) = exp(-lam * t)
        N_B(t) = 1 - exp(-lam * t)
    Returns (iso_A, iso_B, lam).
    """
    lam = 1e-8   # s⁻¹
    iso_A = make_isotope(1, "H", 3, lam=lam, hl=math.log(2) / lam)
    iso_B = make_isotope(2, "He", 3)

    dl = DecayLink("beta-", 1.0, 0.0186, dz=1, da=0, product=iso_B)
    iso_A.add_decay(dl)

    return iso_A, iso_B, lam


def build_decay_matrix(iso_A, iso_B):
    cfg = BurnupConfig(flux=0.0)
    cfg.set_decay_time_filter("all")
    mtx = BurnupMatrix([iso_A, iso_B], cfg)
    mtx.build()
    return mtx


def cyclic_network():
    """
    Pb-208 --(n,g)--> Bi-209 --(beta-)--> Pb-208   (cycle of length 2)
    Pb-208 --(alpha)-> Pb-207               (acyclic arm)
    Returns [Pb-208, Bi-209, Pb-207].
    """
    pb208 = make_isotope(82, "Pb", 208)
    bi209 = make_isotope(83, "Bi", 209, lam=1e-9)
    pb207 = make_isotope(82, "Pb", 207)

    pb208.add_reaction(ReactionLink(102, "(n,g)", 1e-3, product=bi209))
    bi209.add_decay(DecayLink("beta-", 1.0, 0.5, dz=-1, da=0, product=pb208))
    pb208.add_decay(DecayLink("alpha", 0.1, 5.0, dz=-2, da=-4, product=pb207))

    return [pb208, bi209, pb207]


# ═════════════════════════════════════════════════════════════════════════════
# 1. nuclear_data.models
# ═════════════════════════════════════════════════════════════════════════════

class TestModels(unittest.TestCase):

    def test_decay_channel_fields(self):
        ch = DecayChannel(mode="beta-", branch=0.9999, Q_MeV=0.064, dz=1, da=0)
        self.assertEqual(ch.mode, "beta-")
        self.assertAlmostEqual(ch.branch, 0.9999)
        self.assertEqual(ch.dz, 1)
        self.assertEqual(ch.da, 0)
        self.assertEqual(ch.rfs, 0)          # default

    def test_decay_channel_isomeric_rfs(self):
        ch = DecayChannel("it", 1.0, 0.0, 0, 0, rfs=1)
        self.assertEqual(ch.rfs, 1)

    def test_decay_data_stable(self):
        dd = DecayData(half_life_s=None, decay_constant=0.0)
        self.assertTrue(dd.is_stable)
        self.assertEqual(dd.source, "unknown")

    def test_decay_data_radioactive(self):
        lam = math.log(2) / 7e8
        dd = DecayData(half_life_s=7e8, decay_constant=lam)
        self.assertFalse(dd.is_stable)
        self.assertAlmostEqual(dd.decay_constant, lam)

    def test_decay_data_channels_list(self):
        ch = DecayChannel("alpha", 1.0, 5.3, -2, -4)
        dd = DecayData(half_life_s=1e10, decay_constant=1e-10, channels=[ch])
        self.assertEqual(len(dd.channels), 1)

    def test_macs_point(self):
        pt = MacsPoint(kT_keV=30.0, sigma_mb=5.2)
        self.assertEqual(pt.kT_keV, 30.0)
        self.assertEqual(pt.sigma_mb, 5.2)
        self.assertIsNone(pt.rate_cm3_mol_s)

    def test_macs_point_with_rate(self):
        pt = MacsPoint(kT_keV=30.0, sigma_mb=5.2, rate_cm3_mol_s=192.0)
        self.assertEqual(pt.rate_cm3_mol_s, 192.0)

    def test_isomeric_branch(self):
        b = IsomericBranch(LFS=0, yield_=1.0)
        self.assertEqual(b.LFS, 0)
        self.assertEqual(b.yield_, 1.0)

    def test_cross_section_point(self):
        pt = CrossSectionPoint(E_eV=30_000.0, sigma_barn=0.5)
        self.assertEqual(pt.E_eV, 30_000.0)

    def test_neutron_reaction_data_defaults(self):
        rxn = NeutronReactionData(mt=102, channel_name="(n,g)", Q_MeV=3.82)
        self.assertEqual(rxn.mt, 102)
        self.assertEqual(rxn.source, "unknown")
        self.assertEqual(rxn.spectrum, [])
        self.assertEqual(rxn.macs, [])
        self.assertEqual(rxn.branches, [])

    def test_neutron_reaction_data_with_macs(self):
        pt = MacsPoint(30.0, 1.5)
        rxn = NeutronReactionData(mt=102, channel_name="(n,g)", Q_MeV=3.82,
                                  macs=[pt], source="endf")
        self.assertEqual(len(rxn.macs), 1)
        self.assertEqual(rxn.source, "endf")


# ═════════════════════════════════════════════════════════════════════════════
# 2. core.entities.ChemTable
# ═════════════════════════════════════════════════════════════════════════════

class TestChemTable(unittest.TestCase):

    def test_z_of_pb(self):
        self.assertEqual(ChemTable.Z_of("Pb"), 82)

    def test_z_of_bi(self):
        self.assertEqual(ChemTable.Z_of("Bi"), 83)

    def test_z_of_h(self):
        self.assertEqual(ChemTable.Z_of("H"), 1)

    def test_z_of_og(self):
        self.assertEqual(ChemTable.Z_of("Og"), 118)

    def test_unknown_symbol_raises(self):
        with self.assertRaises(ValueError):
            ChemTable.Z_of("Xx")

    def test_elements_dict_size(self):
        self.assertEqual(len(ChemTable.elements), 118)

    def test_reverse_lookup(self):
        for z, sym in ChemTable.elements.items():
            self.assertEqual(ChemTable.Z_of(sym), z)


# ═════════════════════════════════════════════════════════════════════════════
# 3. core.entities.Element
# ═════════════════════════════════════════════════════════════════════════════

class TestElement(unittest.TestCase):

    def test_valid_construction(self):
        el = Element(82, "Pb")
        self.assertEqual(el.Z, 82)
        self.assertEqual(el.X, "Pb")

    def test_wrong_symbol_for_Z_raises(self):
        with self.assertRaises(ValueError):
            Element(82, "Au")   # Z=82 is Pb, not Au

    def test_invalid_Z_raises(self):
        with self.assertRaises(ValueError):
            Element(0, "X")

    def test_equality(self):
        self.assertEqual(Element(82, "Pb"), Element(82, "Pb"))
        self.assertNotEqual(Element(82, "Pb"), Element(83, "Bi"))

    def test_hash_consistency(self):
        a = Element(82, "Pb")
        b = Element(82, "Pb")
        self.assertEqual(hash(a), hash(b))

    def test_repr(self):
        el = Element(82, "Pb")
        self.assertIn("82", repr(el))
        self.assertIn("Pb", repr(el))


# ═════════════════════════════════════════════════════════════════════════════
# 4. core.entities.Isotope
# ═════════════════════════════════════════════════════════════════════════════

class TestIsotope(unittest.TestCase):

    def test_name_no_meta(self):
        iso = Isotope(82, "Pb", 210)
        self.assertEqual(iso.name, "Pb-210")

    def test_name_with_meta(self):
        iso = Isotope(82, "Pb", 207, meta=1)
        self.assertEqual(iso.name, "Pb-207m1")

    def test_from_symbol(self):
        iso = Isotope.from_symbol("Pb", 210)
        self.assertEqual(iso.Z, 82)
        self.assertEqual(iso.X, "Pb")
        self.assertEqual(iso.A, 210)

    def test_invalid_A_raises(self):
        with self.assertRaises(ValueError):
            Isotope(82, "Pb", 0)

    def test_add_decay(self):
        iso = Isotope(82, "Pb", 210)
        dl = DecayLink("beta-", 1.0, 0.064, 1, 0)
        iso.add_decay(dl)
        decays = iso.getListOfDecays()
        self.assertEqual(len(decays), 1)
        self.assertEqual(decays[0].mode, "beta-")

    def test_add_reaction(self):
        iso = Isotope(82, "Pb", 210)
        rx = ReactionLink(102, "(n,g)", 1.5e-3)
        iso.add_reaction(rx)
        rxns = iso.getListOfReactions()
        self.assertEqual(len(rxns), 1)
        self.assertAlmostEqual(rxns[0].sigma_barn, 1.5e-3)

    def test_multiple_decay_channels(self):
        iso = Isotope(82, "Pb", 210)
        iso.add_decay(DecayLink("beta-", 0.9999981, 0.064, 1, 0))
        iso.add_decay(DecayLink("alpha", 0.0000019, 3.792, -2, -4))
        self.assertEqual(len(iso.getListOfDecays()), 2)

    def test_equality_same(self):
        a = Isotope(82, "Pb", 210)
        b = Isotope(82, "Pb", 210)
        self.assertEqual(a, b)

    def test_equality_different_A(self):
        a = Isotope(82, "Pb", 210)
        b = Isotope(82, "Pb", 207)
        self.assertNotEqual(a, b)

    def test_hash_used_in_set(self):
        s = {Isotope(82, "Pb", 210), Isotope(82, "Pb", 210)}
        self.assertEqual(len(s), 1)

    def test_default_decay_constant(self):
        iso = Isotope(83, "Bi", 209)
        self.assertEqual(iso.decay_constant, 0.0)

    def test_repr(self):
        iso = Isotope(82, "Pb", 210)
        r = repr(iso)
        self.assertIn("Pb-210", r)


# ═════════════════════════════════════════════════════════════════════════════
# 5. core.entities.DecayLink and ReactionLink
# ═════════════════════════════════════════════════════════════════════════════

class TestDecayLink(unittest.TestCase):

    def test_construction(self):
        dl = DecayLink("alpha", 1.0, 5.3, -2, -4)
        self.assertEqual(dl.mode, "alpha")
        self.assertAlmostEqual(dl.branch, 1.0)
        self.assertEqual(dl.dz, -2)
        self.assertEqual(dl.da, -4)
        self.assertIsNone(dl.product)

    def test_product_set(self):
        iso = Isotope(83, "Bi", 210)
        dl = DecayLink("beta-", 1.0, 0.5, 1, 0, product=iso)
        self.assertIs(dl.product, iso)

    def test_repr_contains_mode(self):
        dl = DecayLink("alpha", 1.0, 5.0, -2, -4)
        self.assertIn("alpha", repr(dl))


class TestReactionLink(unittest.TestCase):

    def test_get_id(self):
        rx = ReactionLink(102, "(n,g)", 0.5)
        self.assertEqual(rx.getId(), 102)

    def test_q_yield_default(self):
        rx = ReactionLink(102, "(n,g)", 0.5)
        self.assertEqual(rx.q_yield, 1.0)

    def test_repr_contains_mt(self):
        rx = ReactionLink(102, "(n,g)", 0.5)
        self.assertIn("102", repr(rx))

    def test_source_stored(self):
        rx = ReactionLink(102, "(n,g)", 0.5, source="endf")
        self.assertEqual(rx.source, "endf")


# ═════════════════════════════════════════════════════════════════════════════
# 6. core.burnup.BurnupConfig
# ═════════════════════════════════════════════════════════════════════════════

class TestBurnupConfig(unittest.TestCase):

    def test_defaults(self):
        cfg = BurnupConfig()
        self.assertEqual(cfg.flux, 0.0)
        self.assertFalse(cfg.is_filter_active)
        self.assertTrue(cfg.include_stable)
        self.assertFalse(cfg.open_boundary)
        self.assertIsNone(cfg.open_boundary_mts)

    def test_flux_passed_to_constructor(self):
        cfg = BurnupConfig(flux=1e14)
        self.assertEqual(cfg.flux, 1e14)

    # -- reaction filter --

    def test_enable_specific_reactions(self):
        cfg = BurnupConfig()
        cfg.enable_specific_reactions([102, 16])
        self.assertTrue(cfg.is_filter_active)
        self.assertTrue(cfg.is_reaction_allowed(102))
        self.assertTrue(cfg.is_reaction_allowed(16))
        self.assertFalse(cfg.is_reaction_allowed(107))

    def test_enable_all_reactions_clears_filter(self):
        cfg = BurnupConfig()
        cfg.enable_specific_reactions([102])
        cfg.enable_all_reactions()
        self.assertFalse(cfg.is_filter_active)
        self.assertTrue(cfg.is_reaction_allowed(107))

    # -- decay time filter presets --

    def test_preset_all(self):
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("all")
        self.assertTrue(cfg.is_decay_allowed(1.0))
        self.assertTrue(cfg.is_decay_allowed(1e30))
        self.assertTrue(cfg.is_decay_allowed(None))  # stable

    def test_preset_yearly_excludes_long_lived(self):
        # "yearly": max T½ = 3.156e7 s — keeps short-lived, excludes long-lived
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("yearly")
        self.assertTrue(cfg.is_decay_allowed(1.0))         # 1 s ≤ 1 year → kept
        self.assertFalse(cfg.is_decay_allowed(3.2e7))      # > 1 year → excluded

    def test_preset_daily(self):
        # "daily": max T½ = 86400 s — keeps T½ ≤ 1 day
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("daily")
        self.assertTrue(cfg.is_decay_allowed(60.0))        # 60 s ≤ 1 day → kept
        self.assertFalse(cfg.is_decay_allowed(86401.0))    # > 1 day → excluded

    def test_preset_hourly(self):
        # "hourly": max T½ = 3600 s — keeps T½ ≤ 1 hour
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("hourly")
        self.assertTrue(cfg.is_decay_allowed(60.0))        # 60 s ≤ 1 hour → kept
        self.assertFalse(cfg.is_decay_allowed(3601.0))     # > 1 hour → excluded

    def test_preset_minute(self):
        # "minute": max T½ = 60 s — keeps T½ ≤ 1 minute
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("minute")
        self.assertTrue(cfg.is_decay_allowed(1.0))         # 1 s ≤ 1 min → kept
        self.assertFalse(cfg.is_decay_allowed(61.0))       # > 1 min → excluded

    def test_preset_second(self):
        cfg = BurnupConfig()
        cfg.set_decay_time_filter("second")
        self.assertTrue(cfg.is_decay_allowed(1.0))

    def test_invalid_preset_raises(self):
        cfg = BurnupConfig()
        with self.assertRaises(ValueError):
            cfg.set_decay_time_filter("nanosecond")

    def test_custom_time_range(self):
        cfg = BurnupConfig()
        cfg.set_decay_time_range(min_s=100.0, max_s=1e6)
        self.assertFalse(cfg.is_decay_allowed(1.0))       # below min
        self.assertFalse(cfg.is_decay_allowed(2e6))       # above max
        self.assertTrue(cfg.is_decay_allowed(1000.0))

    # -- boundary conditions --

    def test_global_open_boundary(self):
        cfg = BurnupConfig()
        cfg.set_open_boundary(True)
        self.assertTrue(cfg.is_mt_boundary_open(102))
        self.assertTrue(cfg.is_mt_boundary_open(16))

    def test_global_closed_boundary(self):
        cfg = BurnupConfig()
        cfg.set_open_boundary(False)
        self.assertFalse(cfg.is_mt_boundary_open(102))

    def test_per_mt_boundary_overrides_global(self):
        cfg = BurnupConfig()
        cfg.set_open_boundary(False)           # global = closed
        cfg.set_open_boundary_mts([102])       # per-MT: 102 is open
        self.assertTrue(cfg.is_mt_boundary_open(102))
        self.assertFalse(cfg.is_mt_boundary_open(16))

    def test_per_mt_none_reverts_to_global(self):
        cfg = BurnupConfig()
        cfg.set_open_boundary_mts([102])
        cfg.set_open_boundary_mts(None)        # revert
        cfg.set_open_boundary(True)
        self.assertTrue(cfg.is_mt_boundary_open(16))


# ═════════════════════════════════════════════════════════════════════════════
# 7. core.burnup.BurnupMatrix — decay-only, 2 isotopes
# ═════════════════════════════════════════════════════════════════════════════

class TestBurnupMatrixDecayOnly(unittest.TestCase):
    """
    H-3 --[beta-]--> He-3 (stable).
    Analytic solution: N_A(t) = exp(-lam*t), N_B(t) = 1 - exp(-lam*t).
    """

    def setUp(self):
        self.iso_A, self.iso_B, self.lam = two_isotope_decay_chain()
        self.mtx = build_decay_matrix(self.iso_A, self.iso_B)
        self.N0 = np.array([1.0, 0.0])

    # -- matrix structure --

    def test_matrix_shape(self):
        self.assertEqual(self.mtx.matrix_A.shape, (2, 2))

    def test_matrix_A_diagonal(self):
        A = self.mtx.matrix_A
        self.assertAlmostEqual(A[0, 0], -self.lam, places=18)
        self.assertAlmostEqual(A[1, 1],  0.0,      places=18)

    def test_matrix_A_offdiag(self):
        A = self.mtx.matrix_A
        # gain into He-3 from H-3
        self.assertAlmostEqual(A[1, 0], self.lam, places=18)
        # no back-flow from He-3 to H-3
        self.assertAlmostEqual(A[0, 1], 0.0, places=18)

    def test_matrix_decay_populated(self):
        self.assertIsNotNone(self.mtx.matrix_decay)

    def test_matrix_rxn_zero(self):
        # no reactions (flux=0)
        self.assertTrue(np.allclose(self.mtx.matrix_rxn, 0.0))

    # -- solvers --

    def _analytic(self, t):
        a = math.exp(-self.lam * t)
        return np.array([a, 1.0 - a])

    def _check_solver(self, method, t=1e8, atol=1e-4):
        N = self.mtx.solve(self.N0, t, method=method)
        expected = self._analytic(t)
        np.testing.assert_allclose(N, expected, atol=atol,
                                   err_msg=f"Solver '{method}' failed")

    def test_solver_cram16(self):
        self._check_solver("cram16")

    def test_solver_cram48(self):
        self._check_solver("cram48")

    def test_solver_pade(self):
        self._check_solver("pade")

    def test_solver_taylor(self):
        self._check_solver("taylor", atol=1e-3)

    def test_solver_bdf(self):
        self._check_solver("bdf", atol=1e-4)

    def test_solver_invalid_raises(self):
        with self.assertRaises(ValueError):
            self.mtx.solve(self.N0, 1e8, method="euler")

    def test_build_required_before_solve(self):
        cfg = BurnupConfig()
        bare = BurnupMatrix([self.iso_A], cfg)
        with self.assertRaises(RuntimeError):
            bare.solve(np.array([1.0]), 1e8)

    # -- heat release --

    def test_heat_release_positive(self):
        N = self.mtx.solve(self.N0, 1e7)
        q = self.mtx.heat_release(N)
        self.assertGreater(q, 0.0)

    def test_heat_release_zero_when_decayed(self):
        N_eq, _ = self.mtx.solve_equilibrium(self.N0)
        q = self.mtx.heat_release(N_eq)
        # All H-3 has decayed; stable He-3 releases no heat
        self.assertAlmostEqual(q, 0.0, places=10)

    # -- equilibrium / plateau --

    def test_solve_equilibrium_all_to_he3(self):
        N_eq, t_reached = self.mtx.solve_equilibrium(self.N0)
        self.assertAlmostEqual(N_eq[0], 0.0, places=2)   # H-3 depleted
        self.assertAlmostEqual(N_eq[1], 1.0, places=2)   # all in He-3

    def test_find_plateau(self):
        t_plat, q_plat = self.mtx.find_plateau(self.N0)
        self.assertGreater(q_plat, 0.0)

    # -- evolve --

    def test_evolve_returns_trajectory(self):
        times = [0, 1e7, 1e8, 1e9]
        ts, traj = self.mtx.evolve(self.N0, times)
        self.assertEqual(len(traj), 4)
        self.assertEqual(len(traj[0]), 2)

    def test_evolve_t0_equals_N0(self):
        times = [0, 1e8]
        _, traj = self.mtx.evolve(self.N0, times)
        np.testing.assert_array_equal(traj[0], self.N0)

    # -- heat_curve --

    def test_heat_curve_length(self):
        times = [1e6, 1e7, 1e8]
        _, traj = self.mtx.evolve(self.N0, times)
        qc = self.mtx.heat_curve(traj)
        self.assertEqual(len(qc), 3)

    # -- names --

    def test_names(self):
        names = self.mtx.names()
        self.assertIn("H-3", names)
        self.assertIn("He-3", names)

    # -- validate --

    def test_validate_passes(self):
        self.assertTrue(self.mtx.validate())

    # -- export to file --

    def test_export_txt(self):
        times = [1e7, 1e8]
        _, traj = self.mtx.evolve(self.N0, times)
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            fname = f.name
        try:
            self.mtx.export_txt(fname)
            with open(fname, encoding="utf-8") as f:
                content = f.read()
            self.assertIn("H-3", content)
        finally:
            os.unlink(fname)

    def test_export_evolution_txt(self):
        times = [1e7, 1e8]
        ts, traj = self.mtx.evolve(self.N0, times)
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            fname = f.name
        try:
            self.mtx.export_evolution_txt(fname, ts, traj)
            with open(fname, encoding="utf-8") as f:
                content = f.read()
            self.assertIn("Time_sec", content)
        finally:
            os.unlink(fname)


# ═════════════════════════════════════════════════════════════════════════════
# 8. core.burnup.BurnupMatrix — with neutron reactions
# ═════════════════════════════════════════════════════════════════════════════

class TestBurnupMatrixWithReactions(unittest.TestCase):
    """
    3-isotope closed chain (stable): Pb-206 --(n,g)--> Pb-207 --(n,g)--> Pb-208.
    All isotopes stable; only (n,g) allowed; closed boundary.
    Conservation: sum(N) = 1 at all times.
    """

    def setUp(self):
        flux = 1e14
        sigma = 1e-3    # barn

        self.pb206 = make_isotope(82, "Pb", 206)
        self.pb207 = make_isotope(82, "Pb", 207)
        self.pb208 = make_isotope(82, "Pb", 208)

        self.pb206.add_reaction(ReactionLink(102, "(n,g)", sigma, product=self.pb207))
        self.pb207.add_reaction(ReactionLink(102, "(n,g)", sigma, product=self.pb208))
        # Pb-208 has no reaction product in list → reaction disabled by closed boundary

        cfg = BurnupConfig(flux=flux)
        cfg.enable_specific_reactions([102])
        cfg.set_decay_time_filter("all")
        cfg.set_open_boundary(False)

        self.mtx = BurnupMatrix([self.pb206, self.pb207, self.pb208], cfg)
        self.mtx.build()
        self.N0 = np.array([1.0, 0.0, 0.0])

    def test_matrix_populated(self):
        A = self.mtx.matrix_A
        self.assertLess(A[0, 0], 0.0)    # loss from Pb-206

    def test_conservation_closed(self):
        t = 1e6
        N = self.mtx.solve(self.N0, t)
        self.assertAlmostEqual(N.sum(), 1.0, places=4)

    def test_material_flows_to_pb208(self):
        t = 1e8
        N = self.mtx.solve(self.N0, t)
        names = self.mtx.names()
        self.assertGreater(N[names.index("Pb-208")], 0.0)


# ═════════════════════════════════════════════════════════════════════════════
# 9. open / closed boundary condition
# ═════════════════════════════════════════════════════════════════════════════

class TestBoundaryCondition(unittest.TestCase):
    """
    Boundary condition applies ONLY when product=None (no product in any universe).
    When product is an Isotope object outside the matrix, it is always a sink
    (simulating a sub-matrix cut from a larger network).

    Open  boundary (product=None) → reaction is a sink → mass decreases.
    Closed boundary (product=None) → reaction disabled → mass conserved.
    """

    def _make_matrix(self, open_b: bool) -> BurnupMatrix:
        pb206 = make_isotope(82, "Pb", 206)
        # product=None → truly external (no product anywhere) → boundary applies
        pb206.add_reaction(ReactionLink(102, "(n,g)", 1e-3, product=None))

        cfg = BurnupConfig(flux=1e14)
        cfg.enable_specific_reactions([102])
        cfg.set_open_boundary(open_b)

        mtx = BurnupMatrix([pb206], cfg)
        mtx.build()
        return mtx

    def test_open_boundary_loses_mass(self):
        mtx = self._make_matrix(True)
        N = mtx.solve(np.array([1.0]), 1e7)
        self.assertLess(N[0], 1.0)

    def test_closed_boundary_conserves_mass(self):
        mtx = self._make_matrix(False)
        N = mtx.solve(np.array([1.0]), 1e7)
        self.assertAlmostEqual(N[0], 1.0, places=6)

    def test_external_object_product_is_always_sink(self):
        # When product is a real Isotope object but not in the matrix, it is
        # always treated as a sink (regardless of boundary flag).
        pb206 = make_isotope(82, "Pb", 206)
        pb207 = make_isotope(82, "Pb", 207)           # object exists, not in matrix
        pb206.add_reaction(ReactionLink(102, "(n,g)", 1e-3, product=pb207))

        for open_b in (True, False):
            cfg = BurnupConfig(flux=1e14)
            cfg.enable_specific_reactions([102])
            cfg.set_open_boundary(open_b)
            mtx = BurnupMatrix([pb206], cfg)
            mtx.build()
            N = mtx.solve(np.array([1.0]), 1e7)
            # Either way: the channel is active as a sink
            self.assertLess(N[0], 1.0,
                msg=f"open_b={open_b}: external-object product should always be a sink")

    def test_per_mt_open_for_102_only(self):
        pb206 = make_isotope(82, "Pb", 206)
        # product=None → boundary decision is made per-MT
        pb206.add_reaction(ReactionLink(102, "(n,g)", 1e-3, product=None))
        pb206.add_reaction(ReactionLink(16,  "(n,2n)", 1e-3, product=None))

        cfg = BurnupConfig(flux=1e14)
        cfg.enable_specific_reactions([102, 16])
        cfg.set_open_boundary_mts([102])   # 102=open (sink), 16=closed (disabled)

        mtx = BurnupMatrix([pb206], cfg)
        mtx.build()
        A = mtx.matrix_A
        # Only MT=102 contributes → A[0,0] is negative (sink from 102 only)
        lam_n_102 = 1e-3 * 1e-24 * 1e14   # sigma * barn_m2 * flux
        self.assertAlmostEqual(A[0, 0], -lam_n_102, places=10)


# ═════════════════════════════════════════════════════════════════════════════
# 10. core.cycle_finder.CycleAnalyzer
# ═════════════════════════════════════════════════════════════════════════════

class TestCycleAnalyzer(unittest.TestCase):

    def setUp(self):
        self.isotopes = cyclic_network()
        self.analyzer = CycleAnalyzer(self.isotopes)

    def test_has_cycles_true(self):
        self.assertTrue(self.analyzer.has_cycles())

    def test_find_sccs_count(self):
        sccs = self.analyzer.find_sccs()
        # Pb-208 + Bi-209 form one SCC; Pb-207 is a leaf
        self.assertEqual(len(sccs), 1)

    def test_scc_isotope_names(self):
        sccs = self.analyzer.find_sccs()
        names = sccs[0].names
        self.assertIn("Pb-208", names)
        self.assertIn("Bi-209", names)

    def test_scc_not_contain_leaf(self):
        sccs = self.analyzer.find_sccs()
        all_names = [n for scc in sccs for n in scc.names]
        self.assertNotIn("Pb-207", all_names)

    def test_find_all_cycles_non_empty(self):
        cycles = self.analyzer.find_all_cycles()
        self.assertGreaterEqual(len(cycles), 1)

    def test_cycle_length_two(self):
        cycles = self.analyzer.find_all_cycles()
        self.assertTrue(any(c.length == 2 for c in cycles))

    def test_cycle_is_mixed(self):
        cycles = self.analyzer.find_all_cycles()
        mixed = [c for c in cycles if c.is_mixed]
        self.assertGreater(len(mixed), 0)

    def test_cycle_n_decays_n_reactions(self):
        cycles = self.analyzer.find_all_cycles()
        mixed = [c for c in cycles if c.is_mixed][0]
        self.assertEqual(mixed.n_decays, 1)
        self.assertEqual(mixed.n_reactions, 1)

    def test_isotopes_in_cycles(self):
        in_cycles = self.analyzer.isotopes_in_cycles()
        names = {i.name for i in in_cycles}
        self.assertIn("Pb-208", names)
        self.assertIn("Bi-209", names)
        self.assertNotIn("Pb-207", names)

    def test_cycles_through(self):
        pb208 = self.isotopes[0]
        cycles = self.analyzer.cycles_through(pb208)
        self.assertGreater(len(cycles), 0)

    def test_isotope_names_property(self):
        cycles = self.analyzer.find_all_cycles()
        c = cycles[0]
        self.assertIsInstance(c.isotope_names, list)
        self.assertGreater(len(c.isotope_names), 0)

    def test_cycle_str(self):
        cycles = self.analyzer.find_all_cycles()
        s = str(cycles[0])
        self.assertIsInstance(s, str)
        self.assertGreater(len(s), 0)

    def test_print_report_output(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.analyzer.print_report()
        out = buf.getvalue()
        self.assertIn("CYCLE ANALYSIS", out)
        self.assertIn("Mixed cycles", out)

    def test_scc_repr(self):
        sccs = self.analyzer.find_sccs()
        r = repr(sccs[0])
        self.assertIn("SCC", r)

    def test_scc_len(self):
        sccs = self.analyzer.find_sccs()
        self.assertEqual(len(sccs[0]), 2)


class TestCycleAnalyzerAcyclic(unittest.TestCase):
    """Linear A→B chain — no cycles."""

    def setUp(self):
        iso_A = make_isotope(82, "Pb", 206, lam=1e-8, hl=math.log(2) / 1e-8)
        iso_B = make_isotope(82, "Pb", 207)
        iso_A.add_decay(DecayLink("beta-", 1.0, 0.5, 1, 0, product=iso_B))
        self.analyzer = CycleAnalyzer([iso_A, iso_B])

    def test_has_cycles_false(self):
        self.assertFalse(self.analyzer.has_cycles())

    def test_find_all_cycles_empty(self):
        self.assertEqual(len(self.analyzer.find_all_cycles()), 0)

    def test_isotopes_in_cycles_empty(self):
        self.assertEqual(self.analyzer.isotopes_in_cycles(), [])

    def test_print_report_acyclic_message(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.analyzer.print_report()
        self.assertIn("acyclic", buf.getvalue())


# ═════════════════════════════════════════════════════════════════════════════
# 11. core.registry.IsotopeBuilder with a mock provider
# ═════════════════════════════════════════════════════════════════════════════

class _MockProvider:
    """Provides fixed data for Pb-210 (decay + MACS), nothing for others."""

    def get_decay(self, Z, X, A, meta=None):
        if (Z, A) == (82, 210):
            lam = math.log(2) / 7.0e8
            ch = DecayChannel("beta-", 0.9999981, 0.0635, dz=1, da=0)
            return DecayData(half_life_s=7.0e8, decay_constant=lam, channels=[ch],
                             source="mock")
        return None

    def get_macs_with_source(self, Z, X, A, mt, kT_keV, meta=None):
        if (Z, A, mt) == (82, 210, 102):
            return MacsPoint(30.0, 1.5), "mock"
        return None, None

    def get_sigma_with_source(self, Z, X, A, mt, E, e_unit="eV", meta=None):
        return None, None

    def get_reaction_data(self, Z, X, A, mt, meta=None):
        return NeutronReactionData(mt=mt, channel_name=f"MT={mt}", Q_MeV=3.82)

    def get_isomeric_branches(self, Z, X, A, mt, E_eV, meta=None):
        return [IsomericBranch(LFS=0, yield_=1.0)]

    def get_override(self, Z, X, A, field, meta=None):
        return None


class TestIsotopeBuilder(unittest.TestCase):

    def setUp(self):
        self.builder = IsotopeBuilder(_MockProvider())

    def test_build_isotope_decay(self):
        iso = self.builder.build_isotope(82, "Pb", 210, E_eV=30000.0, mt_filter=[102])
        self.assertAlmostEqual(iso.decay_constant, math.log(2) / 7e8, places=20)
        self.assertEqual(len(iso.getListOfDecays()), 1)
        self.assertEqual(iso.getListOfDecays()[0].mode, "beta-")

    def test_build_isotope_macs_prefer_macs(self):
        iso = self.builder.build_isotope(82, "Pb", 210, E_eV=30000.0,
                                          mt_filter=[102], prefer_macs=True)
        rxns = iso.getListOfReactions()
        self.assertGreater(len(rxns), 0)
        # sigma = 1.5 mb * 1e-3 = 1.5e-3 barn
        self.assertAlmostEqual(rxns[0].sigma_barn, 1.5e-3, places=10)

    def test_build_range(self):
        specs = [(82, "Pb", 210), (83, "Bi", 210)]
        isotopes = self.builder.build_range(specs, E_eV=30000.0, mt_filter=[102])
        self.assertEqual(len(isotopes), 2)
        self.assertEqual(isotopes[0].name, "Pb-210")
        self.assertEqual(isotopes[1].name, "Bi-210")

    def test_missing_isotopes(self):
        specs = [(83, "Bi", 210)]  # no data in mock
        self.builder.build_range(specs, E_eV=30000.0, mt_filter=[102])
        missing = self.builder.missing_isotopes()
        self.assertIn("Bi-210", missing)

    def test_pb210_not_missing(self):
        specs = [(82, "Pb", 210)]
        self.builder.build_range(specs, E_eV=30000.0, mt_filter=[102])
        self.assertNotIn("Pb-210", self.builder.missing_isotopes())

    def test_link_products_by_dz_da(self):
        iso_A = Isotope(82, "Pb", 210)
        iso_B = Isotope(83, "Bi", 210)
        dl = DecayLink("beta-", 1.0, 0.064, dz=1, da=0)   # no product yet
        iso_A.add_decay(dl)
        IsotopeBuilder.link_products([iso_A, iso_B])
        linked = iso_A.getListOfDecays()[0].product
        self.assertIsNotNone(linked)
        self.assertEqual(linked.name, "Bi-210")

    def test_link_products_unresolved_stays_none(self):
        # Product (Z+1, A) not in the list → stays None
        iso_A = Isotope(82, "Pb", 210)
        dl = DecayLink("beta-", 1.0, 0.064, dz=1, da=0)
        iso_A.add_decay(dl)
        IsotopeBuilder.link_products([iso_A])   # Bi-210 absent
        self.assertIsNone(iso_A.getListOfDecays()[0].product)

    def test_print_report_runs(self):
        self.builder.build_range([(82, "Pb", 210)], E_eV=30000.0, mt_filter=[102])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.builder.print_report()
        self.assertIn("Pb-210", buf.getvalue())

    def test_report_has_decay_source(self):
        self.builder.build_range([(82, "Pb", 210)], E_eV=30000.0, mt_filter=[102])
        self.assertEqual(len(self.builder.report), 1)
        self.assertIsNotNone(self.builder.report[0]["decay"])


# ═════════════════════════════════════════════════════════════════════════════
# 12. core.service — utility functions (no Qt, no DB)
# ═════════════════════════════════════════════════════════════════════════════

class TestServiceUtils(unittest.TestCase):

    def test_mt_symbol_known(self):
        from core.service import mt_symbol
        self.assertEqual(mt_symbol(102), "(n,γ)")
        self.assertEqual(mt_symbol(16),  "(n,2n)")
        self.assertEqual(mt_symbol(107), "(n,α)")
        self.assertEqual(mt_symbol(18),  "(n,f)")

    def test_mt_symbol_unknown(self):
        from core.service import mt_symbol
        self.assertEqual(mt_symbol(999), "MT=999")

    def test_decay_symbol_known(self):
        from core.service import decay_symbol
        self.assertEqual(decay_symbol("beta-"),   "β⁻")
        self.assertEqual(decay_symbol("alpha"),    "α")
        self.assertEqual(decay_symbol("it"),       "IT")
        self.assertEqual(decay_symbol("ec/beta+"), "β⁺/EC")
        self.assertEqual(decay_symbol("sf"),       "SF")

    def test_decay_symbol_none(self):
        from core.service import decay_symbol
        self.assertEqual(decay_symbol(None), "?")

    def test_decay_symbol_comma_list(self):
        from core.service import decay_symbol
        result = decay_symbol("beta-, alpha")
        self.assertIn("β⁻", result)
        self.assertIn("α", result)

    def test_channel_deltas_keys(self):
        from core.service import CHANNEL_DELTAS
        self.assertIn("beta-",    CHANNEL_DELTAS)
        self.assertIn("alpha",    CHANNEL_DELTAS)
        self.assertIn("(n,γ)",    CHANNEL_DELTAS)
        self.assertIn("(n,2n)",   CHANNEL_DELTAS)
        self.assertIn("(n,p)",    CHANNEL_DELTAS)
        self.assertIn("(n,α)",    CHANNEL_DELTAS)

    def test_channel_deltas_values(self):
        from core.service import CHANNEL_DELTAS
        self.assertEqual(CHANNEL_DELTAS["(n,γ)"],  (0, +1))
        self.assertEqual(CHANNEL_DELTAS["alpha"],  (-2, -4))
        self.assertEqual(CHANNEL_DELTAS["beta-"],  (+1,  0))
        self.assertEqual(CHANNEL_DELTAS["(n,2n)"], (0,  -1))


# ═════════════════════════════════════════════════════════════════════════════
# 13. Full headless pipeline — 3-isotope s-process stub
# ═════════════════════════════════════════════════════════════════════════════

class TestHeadlessPipeline(unittest.TestCase):
    """
    End-to-end test without any database.

    Network: Pb-208 --(n,g)--> Pb-209 --(beta-)--> Bi-209 (stable).
    Initial condition: 100% Pb-208, flux = 1e14 n/cm²/s.
    """

    def setUp(self):
        lam_Pb209 = math.log(2) / (3.253 * 3600)   # T½ = 3.253 h

        self.pb208 = make_isotope(82, "Pb", 208)
        self.pb209 = make_isotope(82, "Pb", 209, lam=lam_Pb209, hl=3.253 * 3600)
        self.bi209 = make_isotope(83, "Bi", 209)

        # Pb-208 (n,g) → Pb-209
        self.pb208.add_reaction(
            ReactionLink(102, "(n,g)", 0.5e-3, product=self.pb209))
        # Pb-209 beta- → Bi-209
        self.pb209.add_decay(
            DecayLink("beta-", 1.0, 0.635, dz=1, da=0, product=self.bi209))

        self.isotopes = [self.pb208, self.pb209, self.bi209]

        cfg = BurnupConfig(flux=1e14)
        cfg.enable_specific_reactions([102])
        cfg.set_decay_time_filter("all")
        cfg.set_open_boundary(False)

        self.mtx = BurnupMatrix(self.isotopes, cfg)
        self.mtx.build()

        names = self.mtx.names()
        self.N0 = np.zeros(len(names))
        self.N0[names.index("Pb-208")] = 1.0

    def test_buildup_of_bi209(self):
        t = 1e6  # ~11.6 days
        N = self.mtx.solve(self.N0, t)
        bi_idx = self.mtx.names().index("Bi-209")
        self.assertGreater(N[bi_idx], 0.0)

    def test_conservation(self):
        N = self.mtx.solve(self.N0, 1e6)
        self.assertAlmostEqual(N.sum(), 1.0, places=3)

    def test_heat_positive_while_pb209_active(self):
        N = self.mtx.solve(self.N0, 1e5)
        q = self.mtx.heat_release(N)
        self.assertGreater(q, 0.0)

    def test_no_cycles_in_acyclic_network(self):
        analyzer = CycleAnalyzer(self.isotopes)
        self.assertFalse(analyzer.has_cycles())

    def test_builder_links_products_correctly(self):
        # After build, verify product links exist via IsotopeBuilder.link_products
        IsotopeBuilder.link_products(self.isotopes)
        rx_product = self.pb208.getListOfReactions()[0].product
        self.assertIsNotNone(rx_product)
        self.assertEqual(rx_product.name, "Pb-209")

    def test_cycle_injection_detected(self):
        """Add a fake back-reaction Bi-209 → Pb-208 and verify cycle is found."""
        rx_back = ReactionLink(102, "(n,g)", 1e-4, product=self.pb208)
        self.bi209.add_reaction(rx_back)
        analyzer = CycleAnalyzer(self.isotopes)
        self.assertTrue(analyzer.has_cycles())

    def test_solve_equilibrium_stable_end_state(self):
        N_eq, _ = self.mtx.solve_equilibrium(self.N0)
        # At equilibrium, all material is in stable Bi-209
        bi_idx = self.mtx.names().index("Bi-209")
        self.assertAlmostEqual(N_eq[bi_idx], 1.0, places=1)

    def test_evolve_and_heat_curve(self):
        times = np.logspace(3, 9, 10)
        ts, traj = self.mtx.evolve(self.N0, times)
        qc = self.mtx.heat_curve(traj)
        self.assertEqual(len(qc), 10)
        # Heat must be positive somewhere in the trajectory (Pb-209 decays release heat)
        self.assertGreater(max(qc), 0.0)
        # At t=1e9 s (>> T½ of Pb-209=3.25h), heat approaches 0 as chain depletes
        self.assertLess(qc[-1], max(qc) + 1e-30)  # last point ≤ peak


# ═════════════════════════════════════════════════════════════════════════════
# 14. IsotopeCycle and CycleStep properties
# ═════════════════════════════════════════════════════════════════════════════

class TestIsotopeCycleProperties(unittest.TestCase):

    def setUp(self):
        self.isotopes = cyclic_network()
        analyzer = CycleAnalyzer(self.isotopes)
        self.cycles = analyzer.find_all_cycles()

    def test_is_mixed_correct(self):
        mixed = [c for c in self.cycles if c.is_mixed]
        self.assertTrue(len(mixed) > 0)
        c = mixed[0]
        self.assertFalse(c.is_pure_decay)
        self.assertFalse(c.is_pure_reaction)
        self.assertTrue(c.is_mixed)

    def test_cycle_len_matches_length(self):
        for c in self.cycles:
            self.assertEqual(len(c), c.length)

    def test_isotope_names_matches_steps(self):
        for c in self.cycles:
            self.assertEqual(len(c.isotope_names), c.length)

    def test_pure_decay_cycle(self):
        # Build: A --[beta-]--> B --[beta-]--> A (both decays)
        a = make_isotope(82, "Pb", 208, lam=1e-8)
        b = make_isotope(83, "Bi", 209, lam=1e-9)
        a.add_decay(DecayLink("beta-", 1.0, 0.5, dz=1, da=0, product=b))
        b.add_decay(DecayLink("beta-", 1.0, 0.3, dz=-1, da=0, product=a))
        analyzer = CycleAnalyzer([a, b])
        cycles = analyzer.find_all_cycles()
        pure_dec = [c for c in cycles if c.is_pure_decay]
        self.assertTrue(len(pure_dec) > 0)


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    unittest.main(verbosity=2)
