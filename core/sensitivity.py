"""
core/sensitivity.py
===================
One-At-a-Time (OAT) sensitivity analysis for the burnup network.

A base run is computed once; then each input parameter is perturbed in turn
(everything else fixed), the matrix is rebuilt and re-solved, and the response
of a chosen target is compared with the base. All physics/maths lives here —
the GUI only feeds a SensitivityConfig and renders the SensitivityResult.

Four sensitivity types (spec §3.2):
  sigma        — per-(isotope, MT) cross-section σ_j  → S = (ΔN/N)/(Δσ/σ)
  halflife     — per-isotope decay constant λ_j       → S = (ΔN/N)/(Δλ/λ)
  flux_energy  — neutron flux Φ and/or evaluation energy E
  library      — discrete: same run on each data library, spread of N_i and Q

Target (spec §3.3): a single isotope concentration N_i, or total heat Q,
evaluated at a fixed time t or at the equilibrium plateau.

Defaults (spec §3.8): two-sided central difference (±δ), δ = 1 %, library
deviation taken about the mean (or a chosen reference library).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from .service import mt_symbol

ProgressCb = Optional[Callable[[float, str], None]]
CancelCb = Optional[Callable[[], bool]]

# MACS / σ-source names that are selected through the provider priority list
# rather than the ENDF-6 decay-library folder switch.
_MACS_SOURCES = {"talys", "eaf2010", "endfb71", "rawmacs", "endf_computed", "user"}
_MACS_DEFAULT_ORDER = ["user", "endfb71", "eaf2010", "rawmacs", "talys", "endf_computed"]


# ──────────────────────────────────────────────────────────────────────────
#  Configuration & results
# ──────────────────────────────────────────────────────────────────────────
@dataclass
class SensitivityConfig:
    initial_by_label: Dict[str, float] = field(default_factory=dict)

    # target -------------------------------------------------------------
    target_kind: str = "concentration"      # "concentration" | "heat"
    target_label: Optional[str] = None       # isotope for the concentration target
    eval_mode: str = "time"                  # "time" | "plateau"
    t_eval_s: float = 1.0e9

    # which analyses to run (multi-select) -------------------------------
    types: List[str] = field(default_factory=lambda: ["sigma"])

    # perturbation -------------------------------------------------------
    delta: float = 0.01
    two_sided: bool = True

    # scope (None → every channel / isotope in the working set) -----------
    scope_channels: Optional[List[Tuple[str, int]]] = None   # for sigma: (iso, MT)
    scope_isotopes: Optional[List[str]] = None               # for halflife

    # run parameters (same meaning as service.run_evolution) -------------
    flux: float = 0.0
    energy_eV: float = 30_000.0
    mt_filter: Optional[List[int]] = None
    method: str = "cram16"                   # fast CRAM by default (spec §3.7)
    hl_preset: str = "all"
    open_boundary_mts: Optional[List[int]] = None
    prefer_macs: bool = False                # astro (MACS σ) vs reactor (point-wise σ)

    # library type -------------------------------------------------------
    libraries: List[str] = field(default_factory=list)
    reference_library: Optional[str] = None  # None → deviation about the mean


@dataclass
class SensRow:
    parameter: str
    type: str
    base: float
    perturbed: float
    S: float
    rank: int = 0


@dataclass
class LibRow:
    library: str
    N_i: Optional[float]
    Q: float
    deviation_pct: float


@dataclass
class SensitivityResult:
    target_desc: str
    base_value: float
    rows: List[SensRow] = field(default_factory=list)          # σ / λ / flux / energy
    library_rows: List[LibRow] = field(default_factory=list)   # library type
    lib_stats: Dict[str, float] = field(default_factory=dict)  # min/max/mean/std

    def ranked(self) -> List[SensRow]:
        """Rows sorted by |S| descending, with rank filled in."""
        ordered = sorted(self.rows, key=lambda r: abs(r.S), reverse=True)
        for i, r in enumerate(ordered, 1):
            r.rank = i
        return ordered


# ──────────────────────────────────────────────────────────────────────────
#  Engine
# ──────────────────────────────────────────────────────────────────────────
class SensitivityEngine:
    """Runs OAT sensitivity over the current NuclearService universe."""

    def __init__(self, service):
        self.svc = service

    # -- public ----------------------------------------------------------
    def run(self, cfg: SensitivityConfig,
            progress: ProgressCb = None,
            should_cancel: CancelCb = None) -> SensitivityResult:
        if not self.svc.initialized:
            raise RuntimeError("Service not initialised. Call init() first.")
        if cfg.target_kind == "concentration" and not cfg.target_label:
            raise ValueError("A target isotope is required for a concentration target.")

        # Remember how the universe was assembled so energy/library rebuilds keep
        # the same Z-range (or label list) instead of falling back to a default.
        self.cfg = cfg
        a = self.svc._build_args
        self._orig_args = dict(a) if a else None
        self._orig_labels = None if a else list(self.svc.universe_labels())
        self._prefer_macs = bool(cfg.prefer_macs)

        # Build the base universe at the requested energy / σ-mode, then the subset.
        self._rebuild_universe(cfg.energy_eV)
        subset = self._working_subset(cfg)

        base_mtx = self._build(subset, cfg)
        base_value = self._evaluate(base_mtx, cfg)

        result = SensitivityResult(target_desc=self._target_desc(cfg),
                                   base_value=base_value)

        # progress accounting over the total perturbed solves
        self._total = self._count_runs(cfg, subset)
        self._done = 0
        self._progress = progress
        self._cancel = should_cancel

        if "sigma" in cfg.types:
            self._sigma(cfg, subset, base_value, result)
        if "halflife" in cfg.types:
            self._halflife(cfg, subset, base_value, result)
        if "flux_energy" in cfg.types:
            self._flux_energy(cfg, subset, base_value, result)
        if "library" in cfg.types:
            self._library(cfg, base_value, result)

        result.rows = result.ranked()
        return result

    # -- universe (re)build, preserving the original range / labels ------
    def _rebuild_universe(self, energy_eV: float, force: bool = False) -> None:
        """Assemble the universe at `energy_eV` with the configured σ-mode, using
        the same Z-range / label list the caller originally built. `force` drops
        the build cache (needed when only the data source changed, not the args).
        """
        if force:
            self.svc._build_args = None
        if self._orig_args is not None:
            o = self._orig_args
            self.svc.build_universe(
                o["z_from"], o["z_to"], o["n_from"], o["n_to"],
                E_eV=energy_eV, mt_filter=self.cfg.mt_filter,
                prefer_macs=self._prefer_macs)
        elif self._orig_labels:
            self.svc.build_universe_from_labels(
                self._orig_labels, E_eV=energy_eV, mt_filter=self.cfg.mt_filter,
                prefer_macs=self._prefer_macs)
        else:
            self.svc._ensure_universe(energy_eV, self.cfg.mt_filter, None)

    # -- working set -----------------------------------------------------
    def _working_subset(self, cfg: SensitivityConfig):
        """Isotopes the target actually depends on (keeps OAT cheap, §3.7).

        Concentration target → minimal source→target sub-network. Heat target →
        the whole universe (every nuclide can contribute heat)."""
        if cfg.target_kind == "concentration" and cfg.target_label:
            return self.svc.minimal_subnetwork(
                cfg.initial_by_label.keys(), [cfg.target_label], cfg.mt_filter)
        return list(self.svc.isotopes)

    # -- matrix build + evaluate ----------------------------------------
    def _build(self, subset, cfg: SensitivityConfig):
        return self.svc.build_matrix(
            flux=cfg.flux, energy_eV=cfg.energy_eV, mt_filter=cfg.mt_filter,
            open_boundary_mts=cfg.open_boundary_mts, hl_preset=cfg.hl_preset,
            isotopes=subset)

    def _solve_vector(self, mtx, cfg: SensitivityConfig):
        """Solve the matrix and return (N_vector, name→index). Name-indexed so it
        survives universe rebuilds."""
        names = mtx.names()
        idx = {n: i for i, n in enumerate(names)}
        N0 = np.zeros(len(names))
        for lab, conc in cfg.initial_by_label.items():
            nm = self.svc.normalize_label(lab)
            if nm in idx:
                N0[idx[nm]] += float(conc)
        if cfg.eval_mode == "plateau":
            N, _ = mtx.solve_equilibrium(N0, method=cfg.method)
        else:
            N = mtx.solve(N0, float(cfg.t_eval_s), cfg.method)
        return N, idx

    def _evaluate(self, mtx, cfg: SensitivityConfig) -> float:
        """Scalar target value (N_i or Q) for a built matrix."""
        N, idx = self._solve_vector(mtx, cfg)
        if cfg.target_kind == "heat":
            return float(mtx.heat_release(N))
        nm = self.svc.normalize_label(cfg.target_label)
        return float(N[idx[nm]]) if nm in idx else 0.0

    def _evaluate_both(self, mtx, cfg: SensitivityConfig):
        """One solve → (N_i or None, Q). Used by the library type, which reports
        both columns per library."""
        N, idx = self._solve_vector(mtx, cfg)
        q = float(mtx.heat_release(N))
        ni = None
        if cfg.target_label:
            nm = self.svc.normalize_label(cfg.target_label)
            ni = float(N[idx[nm]]) if nm in idx else 0.0
        return ni, q

    # -- progress / cancellation ----------------------------------------
    def _tick(self, msg: str) -> None:
        self._done += 1
        if self._cancel and self._cancel():
            raise RuntimeError("Sensitivity analysis cancelled.")
        if self._progress and self._total:
            self._progress(min(1.0, self._done / self._total),
                           f"{msg}  ({self._done}/{self._total})")

    @staticmethod
    def _coeff(base: float, v_plus: float, v_minus: Optional[float], delta: float) -> float:
        """Dimensionless sensitivity S = (ΔN/N)/(Δp/p). Undefined when the base
        response is ~0 → reported as 0."""
        if base == 0 or delta == 0:
            return 0.0
        if v_minus is None:
            return ((v_plus - base) / base) / delta            # one-sided
        return ((v_plus - v_minus) / base) / (2.0 * delta)     # central difference

    # -- σ sensitivity ---------------------------------------------------
    def _sigma(self, cfg, subset, base_value, result):
        for iso, rx in self._sigma_channels(cfg, subset):
            base_sigma = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
            if base_sigma <= 0:
                continue
            mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
            label = f"{iso.name} {mt_symbol(mt)}"

            rx.sigma_barn = base_sigma * (1.0 + cfg.delta)
            v_plus = self._evaluate(self._build(subset, cfg), cfg)
            self._tick(f"σ {label}")
            v_minus = None
            if cfg.two_sided:
                rx.sigma_barn = base_sigma * (1.0 - cfg.delta)
                v_minus = self._evaluate(self._build(subset, cfg), cfg)
                self._tick(f"σ {label}")
            rx.sigma_barn = base_sigma   # restore

            S = self._coeff(base_value, v_plus, v_minus, cfg.delta)
            result.rows.append(SensRow(label, "σ", base_value, v_plus, S))

    def _sigma_channels(self, cfg, subset):
        """Yield (isotope, ReactionLink) for the σ scope."""
        wanted = None
        if cfg.scope_channels:
            wanted = {(self.svc.normalize_label(n), int(m)) for n, m in cfg.scope_channels}
        for iso in subset:
            for rx in iso.getListOfReactions():
                mt = int(rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1))
                if wanted is not None and (iso.name, mt) not in wanted:
                    continue
                yield iso, rx

    # -- λ (half-life) sensitivity --------------------------------------
    def _halflife(self, cfg, subset, base_value, result):
        wanted = None
        if cfg.scope_isotopes:
            wanted = {self.svc.normalize_label(n) for n in cfg.scope_isotopes}
        for iso in subset:
            if wanted is not None and iso.name not in wanted:
                continue
            base_lam = float(getattr(iso, "decay_constant", 0.0) or 0.0)
            if base_lam <= 0:
                continue   # stable → no decay constant to perturb
            label = f"{iso.name} λ"

            iso.decay_constant = base_lam * (1.0 + cfg.delta)
            v_plus = self._evaluate(self._build(subset, cfg), cfg)
            self._tick(f"λ {label}")
            v_minus = None
            if cfg.two_sided:
                iso.decay_constant = base_lam * (1.0 - cfg.delta)
                v_minus = self._evaluate(self._build(subset, cfg), cfg)
                self._tick(f"λ {label}")
            iso.decay_constant = base_lam   # restore

            S = self._coeff(base_value, v_plus, v_minus, cfg.delta)
            result.rows.append(SensRow(label, "λ", base_value, v_plus, S))

    # -- flux / energy sensitivity --------------------------------------
    def _flux_energy(self, cfg, subset, base_value, result):
        # Flux Φ — only the matrix rebuilds (σ unchanged, rate = σ·Φ).
        if cfg.flux > 0:
            v_plus = self._eval_flux(cfg, subset, cfg.flux * (1.0 + cfg.delta))
            self._tick("Φ flux")
            v_minus = None
            if cfg.two_sided:
                v_minus = self._eval_flux(cfg, subset, cfg.flux * (1.0 - cfg.delta))
                self._tick("Φ flux")
            S = self._coeff(base_value, v_plus, v_minus, cfg.delta)
            result.rows.append(SensRow("Flux Φ", "Φ", base_value, v_plus, S))

        # Energy E — the universe is re-evaluated (σ re-interpolated at new E).
        v_plus = self._eval_energy(cfg, cfg.energy_eV * (1.0 + cfg.delta))
        self._tick("E energy")
        v_minus = None
        if cfg.two_sided:
            v_minus = self._eval_energy(cfg, cfg.energy_eV * (1.0 - cfg.delta))
            self._tick("E energy")
        # restore the universe at the base energy
        self._rebuild_universe(cfg.energy_eV)
        S = self._coeff(base_value, v_plus, v_minus, cfg.delta)
        result.rows.append(SensRow("Energy E", "E", base_value, v_plus, S))

    def _eval_flux(self, cfg, subset, flux):
        mtx = self.svc.build_matrix(
            flux=flux, energy_eV=cfg.energy_eV, mt_filter=cfg.mt_filter,
            open_boundary_mts=cfg.open_boundary_mts, hl_preset=cfg.hl_preset,
            isotopes=subset)
        return self._evaluate(mtx, cfg)

    def _eval_energy(self, cfg, energy_eV):
        # Rebuild the universe at the perturbed energy, then a fresh subset.
        self._rebuild_universe(energy_eV)
        subset = self._working_subset(cfg)
        mtx = self.svc.build_matrix(
            flux=cfg.flux, energy_eV=energy_eV, mt_filter=cfg.mt_filter,
            open_boundary_mts=cfg.open_boundary_mts, hl_preset=cfg.hl_preset,
            isotopes=subset)
        return self._evaluate(mtx, cfg)

    # -- library sensitivity (discrete) ---------------------------------
    def _library(self, cfg, base_value, result):
        base_decay = getattr(self.svc, "decay_library", None)
        base_macs = (self.svc.provider.get_macs_priority()
                     if self.svc.provider else list(_MACS_DEFAULT_ORDER))
        decay_libs = set(self.svc.decay_libraries())

        ni_vals: List[Optional[float]] = []
        q_vals: List[float] = []
        rows: List[Tuple[str, Optional[float], float]] = []

        for lib in cfg.libraries:
            self._apply_library(lib, decay_libs, base_decay, base_macs)
            self._rebuild_universe(cfg.energy_eV, force=True)   # σ re-evaluated
            subset = self._working_subset(cfg)
            mtx = self.svc.build_matrix(
                flux=cfg.flux, energy_eV=cfg.energy_eV, mt_filter=cfg.mt_filter,
                open_boundary_mts=cfg.open_boundary_mts, hl_preset=cfg.hl_preset,
                isotopes=subset)
            ni, q = self._evaluate_both(mtx, cfg)
            rows.append((lib, ni, q))
            if ni is not None:
                ni_vals.append(ni)
            q_vals.append(q)
            self._tick(f"library {lib}")

        # restore base library state + universe
        self._restore_library(base_decay, base_macs)
        self._rebuild_universe(cfg.energy_eV, force=True)

        # reference for the deviation column
        series = ni_vals if (cfg.target_label and ni_vals) else q_vals
        ref = self._reference_value(cfg, rows, series)
        for lib, ni, q in rows:
            val = ni if (cfg.target_label and ni is not None) else q
            dev = 100.0 * (val - ref) / ref if ref else 0.0
            result.library_rows.append(LibRow(lib, ni, q, dev))
        if series:
            arr = np.array(series, dtype=float)
            result.lib_stats = {"min": float(arr.min()), "max": float(arr.max()),
                                "mean": float(arr.mean()), "std": float(arr.std())}

    def _reference_value(self, cfg, rows, series) -> float:
        if cfg.reference_library:
            for lib, ni, q in rows:
                if lib == cfg.reference_library:
                    return (ni if (cfg.target_label and ni is not None) else q) or 0.0
        return float(np.mean(series)) if series else 0.0

    def _apply_library(self, lib, decay_libs, base_decay, base_macs):
        """Reset to base, then apply the one library under test."""
        self._restore_library(base_decay, base_macs)
        if lib in decay_libs:
            self.svc.switch_decay_library(lib)
        elif lib.lower() in _MACS_SOURCES:
            key = lib.lower()
            self.svc.set_macs_priority([key] + [s for s in _MACS_DEFAULT_ORDER if s != key])

    def _restore_library(self, base_decay, base_macs):
        if base_decay and base_decay != getattr(self.svc, "decay_library", None):
            self.svc.switch_decay_library(base_decay)
        if base_macs:
            self.svc.set_macs_priority(list(base_macs))

    # -- misc ------------------------------------------------------------
    def _target_desc(self, cfg: SensitivityConfig) -> str:
        when = "plateau" if cfg.eval_mode == "plateau" else f"t={cfg.t_eval_s:g}s"
        if cfg.target_kind == "heat":
            return f"Heat Q @ {when}"
        return f"N[{cfg.target_label}] @ {when}"

    def _count_runs(self, cfg: SensitivityConfig, subset) -> int:
        mult = 2 if cfg.two_sided else 1
        total = 0
        if "sigma" in cfg.types:
            total += mult * sum(1 for _ in self._sigma_channels(cfg, subset))
        if "halflife" in cfg.types:
            wanted = ({self.svc.normalize_label(n) for n in cfg.scope_isotopes}
                      if cfg.scope_isotopes else None)
            total += mult * sum(1 for iso in subset
                                if (wanted is None or iso.name in wanted)
                                and float(getattr(iso, "decay_constant", 0.0) or 0.0) > 0)
        if "flux_energy" in cfg.types:
            total += mult * (2 if cfg.flux > 0 else 1)
        if "library" in cfg.types:
            total += len(cfg.libraries)
        return max(1, total)
