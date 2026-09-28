"""
core/registry.py
================
IsotopeBuilder — builds isotopes from data (via NuclearDataProvider)
and links them into a transition graph.

This is the ONLY place where domain entities meet the databases.
The burnup matrix does not use this layer — it works with a ready-made list.
"""

from __future__ import annotations
from typing import List, Optional

from .entities import Isotope, DecayLink, ReactionLink, ChemTable


class IsotopeBuilder:
    """
    Builds a list of isotopes and wires up decay/reaction products.

    provider — an instance of nuclear_data.NuclearDataProvider.
    """

    def __init__(self, provider):
        self.provider = provider
        # data-provenance report, filled in by build_range
        self.report: List[dict] = []

    # ── build a single isotope with data ──
    def build_isotope(self, Z, X, A, meta=None,
                      E_eV: float = 0.0253,
                      mt_filter: Optional[List[int]] = None,
                      audit: Optional[dict] = None,
                      prefer_macs: bool = False) -> Isotope:
        """
        Creates an isotope and populates it with decays and reactions from provider.
        E_eV       — energy at which cross sections and isomeric branches are taken.
        mt_filter  — if given, only these MT channels are loaded.
        prefer_macs — if True, cross section is taken from MACS (kT=E_eV) FIRST,
                      with point-wise sigma(E) from ENDF as fallback. This is the
                      correct order for astrophysics (s-process at kT), where the
                      Maxwell-averaged quantity is required.
        audit      — if a dict is passed, data provenance is written into it.
        """
        iso = Isotope(Z, X, A, meta)

        # --- Decay ---
        dd = self.provider.get_decay(Z, X, A, meta)
        if dd is not None:
            iso.half_life_s = dd.half_life_s
            iso.decay_constant = dd.decay_constant
            for ch in dd.channels:
                iso.add_decay(DecayLink(ch.mode, ch.branch, ch.Q_MeV,
                                        ch.dz, ch.da, ch.rfs))
        if audit is not None:
            audit["decay"] = dd.source if dd is not None else None
            audit["reactions"] = {}
            audit["missing"] = []

        # --- Neutron reactions ---
        mts = mt_filter if mt_filter else [102, 16, 18, 103, 107]
        kT_keV = E_eV / 1e3
        for mt in mts:
            src = None
            sigma = None
            try:
                if prefer_macs:
                    # 1) MACS (kT), 2) point-wise sigma(E) from ENDF
                    pt, msrc = self.provider.get_macs_with_source(Z, X, A, mt, kT_keV, meta)
                    if pt:
                        sigma = pt.sigma_mb * 1e-3   # mb -> barn
                        src = msrc or "macs"
                    if sigma is None or sigma <= 0:
                        sigma, ssrc = self.provider.get_sigma_with_source(Z, X, A, mt, E_eV, "eV", meta)
                        if sigma and sigma > 0:
                            src = ssrc or "endf"
                else:
                    # 1) point-wise sigma(E), 2) MACS
                    sigma, ssrc = self.provider.get_sigma_with_source(Z, X, A, mt, E_eV, "eV", meta)
                    if sigma is not None and sigma > 0:
                        src = ssrc or "endf"
                    else:
                        pt, msrc = self.provider.get_macs_with_source(Z, X, A, mt, kT_keV, meta)
                        if pt:
                            sigma = pt.sigma_mb * 1e-3
                            src = msrc or "macs"
            except Exception as _exc:
                print(f"[registry] sigma lookup error {iso.name} MT={mt}: {_exc}")
                sigma = None
            if sigma is None or sigma <= 0:
                if audit is not None:
                    audit["missing"].append(mt)
                continue

            rxn_data = self.provider.get_reaction_data(Z, X, A, mt, meta)
            q_mev = rxn_data.Q_MeV if (rxn_data and rxn_data.Q_MeV is not None) else 0.0

            branches = self.provider.get_isomeric_branches(Z, X, A, mt, E_eV, meta)
            name = rxn_data.channel_name if rxn_data else f"MT={mt}"
            for b in branches:
                iso.add_reaction(ReactionLink(
                    mt=mt, name=name,
                    sigma_barn=sigma * b.yield_,
                    q_yield=b.yield_, Q_MeV=q_mev, lfs=b.LFS, source=src))
            if audit is not None:
                audit["reactions"][mt] = src
        return iso

    # ── build a set of isotopes ──
    def build_range(self, specs, E_eV: float = 0.0253,
                    mt_filter: Optional[List[int]] = None,
                    prefer_macs: bool = False) -> List[Isotope]:
        """
        specs       — list of tuples (Z, X, A) or (Z, X, A, meta).
        prefer_macs — see build_isotope (set True for astrophysics).
        Returns a linked list of isotopes.
        Also fills self.report with data provenance for each isotope.
        """
        self.report = []
        isotopes = []
        for s in specs:
            Z, X, A = s[0], s[1], s[2]
            meta = s[3] if len(s) > 3 else None
            audit = {"isotope": None}
            iso = self.build_isotope(Z, X, A, meta, E_eV, mt_filter, audit, prefer_macs)
            audit["isotope"] = iso.name
            self.report.append(audit)
            isotopes.append(iso)
        self.link_products(isotopes)
        # Print diagnostics for isotopes with no sigma data at all
        for a in self.report:
            missing = a.get("missing", [])
            if missing and not a.get("reactions"):
                print(f"[registry] {a['isotope']}: no sigma data found for MTs {missing}")
        return isotopes

    # ── source diagnostics ──
    def missing_isotopes(self) -> List[str]:
        """Isotopes for which neither decay nor reaction data was found anywhere."""
        out = []
        for a in self.report:
            no_decay = a.get("decay") is None
            no_rxn = not a.get("reactions")
            if no_decay and no_rxn:
                out.append(a["isotope"])
        return out

    def print_report(self):
        """Prints a table: isotope | decay source | reaction sources."""
        print(f"{'isotope':10s} {'decay':8s} reactions (MT:source)")
        print("-" * 60)
        for a in self.report:
            dec = a.get("decay") or "—"
            rxns = a.get("reactions", {})
            rxn_str = ", ".join(f"{mt}:{src}" for mt, src in rxns.items()) or "—"
            miss = a.get("missing", [])
            miss_str = f"  [missing: {miss}]" if miss else ""
            print(f"{a['isotope']:10s} {dec:8s} {rxn_str}{miss_str}")


    # ── product linking ──
    @staticmethod
    def link_products(isotopes: List[Isotope]):
        """Sets the product references on decay and reaction links by Z/A/meta."""
        index = {(i.Z, i.A, i.meta): i for i in isotopes}

        def find(Z, A, lfs):
            meta = None if lfs == 0 else lfs
            return index.get((Z, A, meta)) or index.get((Z, A, None))

        for iso in isotopes:
            for d in iso.getListOfDecays():
                d.product = find(iso.Z + d.dz, iso.A + d.da, d.rfs)
            for r in iso.getListOfReactions():
                dZ, dA = IsotopeBuilder._reaction_shift(r.mt)
                r.product = find(iso.Z + dZ, iso.A + dA, r.lfs)

    @staticmethod
    def _reaction_shift(mt):
        """(dZ, dA) change for the main neutron channels."""
        return {
            102: (0, 1),    # (n,g): +1 neutron
            16:  (0, -1),   # (n,2n): -1 neutron
            103: (-1, 0),   # (n,p)
            107: (-2, -3),  # (n,alpha)
            18:  (0, 0),    # fission — fragments handled separately
        }.get(mt, (0, 0))
