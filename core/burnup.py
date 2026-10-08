"""
core/burnup.py
==============
BurnupConfig + BurnupMatrix.

KEY POINT: BurnupMatrix accepts ANY list of isotope objects.
Each isotope must provide the interface:
    .Z, .A, .name
    .decay_constant          -> float [s⁻¹]
    .getListOfDecays()       -> objects with .branch, .product, (.Q_MeV)
    .getListOfReactions()    -> objects with .getId()/.mt, .sigma_barn,
                                .q_yield, .product, (.Q_MeV)

This is exactly the core.entities.Isotope interface, but any custom
class with the same attributes works — the matrix is database-agnostic.

Filters in BurnupConfig:
  - neutron reactions — by MT list;
  - decays — by half-life (all / yearly / daily / ...).

Instead of input() for missing cross sections, the matrix collects a list
of "holes" and optionally invokes a callback. The GUI/script decides what to substitute.
"""

from __future__ import annotations
import math
from typing import Callable, Dict, List, Optional, Set

import numpy as np
import scipy.linalg as sla
from scipy.integrate import solve_ivp
from scipy.linalg import expm

import warnings
try:
    from scipy.linalg import LinAlgWarning
    warnings.filterwarnings("ignore", category=LinAlgWarning)
except Exception:
    pass


class BurnupConfig:
    """Calculation parameters: flux, reaction filter, decay time filter."""

    # half-life threshold presets [seconds]
    _PRESETS = {
        "all":     (0.0, None),
        "yearly":  (0.0, 3.156e7),
        "daily":   (0.0, 86400.0),
        "hourly":  (0.0, 3600.0),
        "minute":  (0.0, 60.0),
        "second":  (0.0, 1.0),

    }

    def __init__(self, flux: float = 0.0):
        self.flux = float(flux)
        # reaction filter
        self.is_filter_active = False
        self.allowed_mts: Set[int] = set()
        # decay time filter
        self.min_half_life_s: float = 0.0
        self.max_half_life_s: Optional[float] = None
        self.include_stable: bool = True
        # decay mode filter: None = all allowed; set of strings = only those modes
        self.decay_mode_filter: Optional[Set[str]] = None
        # manual decay links: list of dicts with keys:
        #   parent (str), daughter (str or ""), mode (str), half_life_s (float), branch (float)
        self.manual_decay_links: list = []
        # global boundary condition for neutron reactions (fallback when
        # open_boundary_mts is None):
        # True  — open: reactions of boundary isotopes (product outside list) = sinks
        # False — closed: such reactions are disabled; decays are always active
        self.open_boundary: bool = False
        # per-MT boundary override: if set, each MT is checked individually;
        # MTs in this set are open (sink), others are closed.
        # None means: use the global open_boundary flag for all MTs.
        self.open_boundary_mts: Optional[Set[int]] = None

    def set_open_boundary(self, open_b: bool = True) -> None:
        """
        True  — open: neutron reactions with product outside the list = sinks
                (material flows out of the chain).
        False — closed: neutron reactions of boundary isotopes are disabled.
                Decays are always active regardless of this flag.
        """
        self.open_boundary = bool(open_b)

    def set_open_boundary_mts(self, mts) -> None:
        """Per-MT open boundary control.

        Only the listed MT reactions act as sinks when their product is absent
        from the network; every other MT is treated as closed.  Passing None
        falls back to the global open_boundary flag.
        """
        self.open_boundary_mts = {int(m) for m in mts} if mts is not None else None

    def is_mt_boundary_open(self, mt: int) -> bool:
        """True if reactions of this MT type use open boundary.

        Returns per-MT setting when available, otherwise falls back to the
        global open_boundary flag.
        """
        if self.open_boundary_mts is not None:
            return int(mt) in self.open_boundary_mts
        return bool(self.open_boundary)

    # ── reactions ──
    # ── decay mode filter ──
    def set_decay_mode_filter(self, modes) -> None:
        """Allow only the listed decay modes (e.g. ['beta-']). None = allow all."""
        self.decay_mode_filter = (
            {str(m).lower() for m in modes} if modes is not None else None)

    def is_decay_mode_allowed(self, mode) -> bool:
        if self.decay_mode_filter is None:
            return True
        return str(mode).lower() in self.decay_mode_filter

    # ── reactions ──
    def enable_all_reactions(self):
        self.is_filter_active = False
        self.allowed_mts.clear()

    def enable_specific_reactions(self, mt_list):
        self.is_filter_active = True
        self.allowed_mts = {int(m) for m in mt_list}

    def is_reaction_allowed(self, mt) -> bool:
        if not self.is_filter_active:
            return True
        try:
            return int(mt) in self.allowed_mts
        except (ValueError, TypeError):
            return False

    # ── decay time filter ──
    def set_decay_time_filter(self, preset: str = "all"):
        """preset: all | yearly | daily | hourly | minute | second"""
        if preset not in self._PRESETS:
            raise ValueError(f"Unknown preset '{preset}'. Available: {list(self._PRESETS)}")
        self.min_half_life_s, self.max_half_life_s = self._PRESETS[preset]

    def set_decay_time_range(self, min_s: float = 0.0, max_s: Optional[float] = None):
        self.min_half_life_s = float(min_s)
        self.max_half_life_s = max_s

    def is_decay_allowed(self, half_life_s: Optional[float]) -> bool:
        """Returns True if a decay with this half-life passes the filter."""
        if half_life_s is None:        # stable
            return self.include_stable
        if half_life_s < self.min_half_life_s:
            return False
        if self.max_half_life_s is not None and half_life_s > self.max_half_life_s:
            return False
        return True


class BurnupMatrix:
    """
    Assembly and solution of the burnup equation system.

    isotopes — ANY list of isotopes (from a database or assembled manually).
    """

    def __init__(self, isotopes: List, config: BurnupConfig):
        self.isotopes = list(isotopes)
        self.config = config
        self.flux = config.flux
        self.index = {iso: i for i, iso in enumerate(self.isotopes)}

        self.matrix_A: Optional[np.ndarray] = None
        self.matrix_decay: Optional[np.ndarray] = None
        self.matrix_rxn: Optional[np.ndarray] = None
        self.matrix_heat: Optional[np.ndarray] = None   # heat release [MeV/s]

        self.missing: List[dict] = []   # entries without a cross section

    # ── product index ──
    def _find_index(self, product) -> Optional[int]:
        if product is None:
            return None
        if product in self.index:
            return self.index[product]
        # fallback: search by name
        name = getattr(product, "name", str(product))
        for iso, idx in self.index.items():
            if getattr(iso, "name", None) == name:
                return idx
        return None

    # ── assembly ──
    def build(self, missing_callback: Optional[Callable[[dict], float]] = None):
        """
        Builds matrix_A.

        Boundary condition (config.open_boundary):
          True  — open: neutron reactions of "boundary" isotopes (product outside
                  the list) are sinks: material is lost, heat is counted.
          False — closed: neutron reactions with product outside the list are
                  disabled (reaction does not occur). Decays are always active.

        missing_callback(info) -> float|None: called when sigma is absent.
        """
        manual_keys = set()
        by_name = {iso.name: idx for idx, iso in enumerate(self.isotopes)}
        for link in self.config.manual_decay_links:
            parent = link.get("parent", "")
            daughter = link.get("daughter", "")
            hl = float(link.get("half_life_s", 0))
            branch = float(link.get("branch", 1))
            if parent not in by_name or (daughter and daughter not in by_name):
                raise ValueError(f"Manual decay requires isotopes in the matrix: {parent} -> {daughter}")
            if not np.isfinite(hl) or hl <= 0 or not np.isfinite(branch) or not 0 < branch <= 1:
                raise ValueError("Manual decay requires a positive half-life and 0 < branch <= 1.")
            key = (parent, str(link.get("mode", "")).lower(), daughter)
            if key in manual_keys:
                raise ValueError(f"Duplicate manual decay: {parent} -> {daughter}")
            manual_keys.add(key)
        N = len(self.isotopes)
        self.matrix_decay = np.zeros((N, N))
        self.matrix_rxn   = np.zeros((N, N))
        # matrix_heat is diagonal: matrix_heat[j,j] = total heat release
        # coefficient for isotope j [MeV·s⁻¹ per unit concentration].
        # heat_release(N) = sum_j matrix_heat[j,j] * N[j]
        self.matrix_heat  = np.zeros((N, N))
        self.missing.clear()

        for j, parent in enumerate(self.isotopes):

            # 1) Decays (with time filter)
            # Decays are always active — regardless of open_boundary.
            # If product is outside the list: sink (material is lost, heat counted).
            hl = getattr(parent, "half_life_s", None)
            if self.config.is_decay_allowed(hl):
                lam = float(getattr(parent, "decay_constant", 0.0))
                if lam > 0:
                    for d in parent.getListOfDecays():
                        key = (parent.name, str(d.mode).lower(),
                               getattr(getattr(d, "product", None), "name", ""))
                        if key in manual_keys:
                            continue
                        if not self.config.is_decay_mode_allowed(getattr(d, "mode", "")):
                            continue
                        br    = float(getattr(d, "branch", 1.0))
                        lam_b = lam * br
                        q     = float(getattr(d, "Q_MeV", 0.0))
                        idx   = self._find_index(getattr(d, "product", None))
                        if idx is not None:
                            # product is in the chain
                            self.matrix_decay[j,   j] -= lam_b
                            self.matrix_decay[idx, j] += lam_b
                            self.matrix_heat[j, j]    += lam_b * q
                        else:
                            # product outside the chain: sink (always, in both modes)
                            self.matrix_decay[j, j] -= lam_b
                            self.matrix_heat[j, j]  += lam_b * q

            # 2) Neutron reactions (with MT filter)
            # open_boundary=True:  reactions of "boundary" isotopes are sinks.
            # open_boundary=False: reactions with product outside the list are disabled.
            if self.flux > 0:
                for rx in parent.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if not self.config.is_reaction_allowed(mt):
                        continue

                    sigma = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                    if sigma <= 0:
                        info = {"isotope": getattr(parent, "name", str(parent)), "mt": mt}
                        if missing_callback:
                            val = missing_callback(info)
                            if val:
                                sigma = float(val)
                                rx.sigma_barn = sigma
                        if sigma <= 0:
                            self.missing.append(info)
                            continue

                    lam_n = sigma * 1e-24 * self.flux   # [s⁻¹]
                    q     = float(getattr(rx, "Q_MeV", 0.0))
                    prod  = getattr(rx, "product", None)
                    idx   = self._find_index(prod)
                    if idx is not None:
                        # product in chain — always active regardless of boundary mode
                        self.matrix_rxn[j,   j] -= lam_n
                        self.matrix_rxn[idx, j] += lam_n
                        self.matrix_heat[j, j]  += lam_n * q
                    elif prod is not None:
                        # Product exists in the full universe but lies outside this
                        # (sub-)matrix. In the full network this channel was active
                        # and removed material from the parent, so keep that loss as
                        # a sink. This is what makes a minimal sub-network solve
                        # numerically identical to the full run for the kept nuclei.
                        self.matrix_rxn[j, j]  -= lam_n
                        self.matrix_heat[j, j] += lam_n * q
                    elif self.config.is_mt_boundary_open(mt):
                        # product outside the full universe: honour the boundary mode —
                        # open boundary, so the isotope acts as a sink.
                        self.matrix_rxn[j, j]  -= lam_n
                        self.matrix_heat[j, j] += lam_n * q
                    # closed boundary for a truly-external product: reaction disabled — skip

        # Manual entries override matching library channels and bypass filters.
        for link in self.config.manual_decay_links:
            p_idx = by_name[link["parent"]]
            rate = np.log(2) / float(link["half_life_s"]) * float(link.get("branch", 1))
            daughter = link.get("daughter", "")
            self.matrix_decay[p_idx, p_idx] -= rate
            if daughter:
                self.matrix_decay[by_name[daughter], p_idx] += rate
            q = link.get("Q_MeV")
            if q is None:
                q = next((d.Q_MeV for d in self.isotopes[p_idx].getListOfDecays()
                          if str(d.mode).lower() == str(link.get("mode", "")).lower()
                          and getattr(getattr(d, "product", None), "name", "") == daughter), 0.0)
            self.matrix_heat[p_idx, p_idx] += rate * float(q or 0.0)

        self.matrix_A = self.matrix_decay + self.matrix_rxn
        return self.matrix_A

    # ── solvers ──
    def solve(self, N0, t_sec: float, method: str = "cram16") -> np.ndarray:
        if self.matrix_A is None:
            raise RuntimeError("Call build() first.")
        N0 = np.asarray(N0, dtype=float)
        m = method.lower()
        if m == "cram16":
            from .cram import cram16_step
            if not np.isfinite(t_sec) or t_sec < 0:
                raise ValueError("Evolution duration must be finite and non-negative.")
            if t_sec == 0:
                return N0.copy()
            return cram16_step(self.matrix_A, N0, t_sec)
        if m == "cram48":
            from .cram import cram48_step
            return cram48_step(self.matrix_A, N0, t_sec)
        if m == "cram16_adaptive":
            from .cram import adaptive_cram16
            return adaptive_cram16(self.matrix_A, N0, t_sec)
        At = self.matrix_A * t_sec
        if m == "pade":
            return expm(At).dot(N0)
        if m == "taylor":
            return self._taylor(N0, At)
        if m == "bdf":
            sol = solve_ivp(lambda t, y: self.matrix_A.dot(y),
                            [0, t_sec], N0, method="BDF")
            return sol.y[:, -1]
        raise ValueError(f"Method '{method}' is not supported.")

    def _taylor(self, N0, At, order=20):
        res = N0.astype(float).copy()
        term = N0.astype(float).copy()
        for i in range(1, order + 1):
            term = At.dot(term) / i
            res = res + term
        return res

    def _cram16(self, N0, At):
        from .cram import cram16_step
        return cram16_step(At, N0, 1.0)

    def _cram48(self, N0, At):
        """Order-48 incomplete partial factorization."""
        from .cram import cram48_step
        return cram48_step(At, N0, 1.0)

    # ── heat release ──
    def heat_release(self, N) -> float:
        """Total heat release [MeV/s] for concentration vector N."""
        if self.matrix_heat is None:
            raise RuntimeError("Call build() first.")
        N = np.asarray(N, dtype=float)
        # matrix_heat is diagonal: heat_j = matrix_heat[j,j] * N[j]
        return float(np.dot(np.diag(self.matrix_heat), N))

    # ── equilibrium (plateau) ──
    def solve_equilibrium(self, N0, method: str = "cram16",
                          rtol: float = 1e-4, t_start: float = None,
                          t_max: float = 1e25, max_steps: int = 60):
        """
        Finds the equilibrium (plateau) concentration slice — solves at
        increasing t until the N vector stops changing (relative dN < rtol).

        Returns (N_eq, t_reached). If plateau is not reached within t_max,
        returns the last slice and t_max.
        """
        if self.matrix_A is None:
            raise RuntimeError("Call build() first.")
        N0 = np.asarray(N0, dtype=float)

        # Automatically choose starting time: 10 half-lives of the slowest
        # active process in the matrix. Without this, at small flux (lambda_n << 1)
        # the comparison N(t_start) ~ N(3*t_start) ~ N0 gives false convergence
        # long before reaching the true plateau (100-300 years away).
        if t_start is None:
            diag_rates = np.abs(np.diag(self.matrix_A))
            nonzero = diag_rates[diag_rates > 0]
            if len(nonzero) > 0:
                t_start = min(10.0 / nonzero.min(), t_max / 10.0)
            else:
                t_start = 1e3

        t = t_start
        prev = self.solve(N0, t, method)
        for _ in range(max_steps):
            t *= 3.0
            if t > t_max:
                t = t_max
            cur = self.solve(N0, t, method)
            denom = np.maximum(np.abs(cur), 1e-30)
            rel = np.max(np.abs(cur - prev) / denom)
            if rel < rtol:
                return cur, t
            prev = cur
            if t >= t_max:
                break
        return prev, t

    def find_plateau(self, N0, method: str = "cram16",
                     n_scan: int = 80,
                     t_min: float = 1e-3, t_max: float = 1e25
                     ) -> tuple:
        """
        Finds the time and value of the heat plateau.

        For closed chains Q grows monotonically to steady-state.
        For open chains — Q passes through a maximum (secular equilibrium),
        then falls to zero as material flows out.
        In both cases max Q(t) = Q at the plateau.

        Returns (t_plateau, Q_plateau).
        """
        if self.matrix_heat is None:
            raise RuntimeError("Call build() first.")
        N0 = np.asarray(N0, dtype=float)
        times = np.logspace(np.log10(t_min), np.log10(t_max), n_scan)
        q_best, t_best = 0.0, times[0]
        for t in times:
            q = self.heat_release(self.solve(N0, float(t), method))
            if q > q_best:
                q_best, t_best = q, t
        return t_best, q_best

    def equilibrium_heat(self, N0, method: str = "cram16",
                         t_plateau: float = None,
                         n_scan: int = 80,
                         t_min: float = 1e-3, t_max: float = 1e25) -> float:
        """
        Heat density [MeV/s] at the equilibrium plateau.

        Strategy depends on config.open_boundary:
          False (closed chain) — finds true steady-state via
              solve_equilibrium(): N(t) ~ N(3t), adaptive t_start.
              For closed chains Q grows monotonically to the plateau.
          True (open chain) — finds max Q(t) via find_plateau():
              Q grows -> peak (secular equilibrium) -> falls to zero.

        t_plateau — specify time manually [s], overrides auto-search.
        """
        if self.matrix_heat is None:
            raise RuntimeError("Call build() first.")
        N0 = np.asarray(N0, dtype=float)

        if t_plateau is not None:
            return self.heat_release(self.solve(N0, float(t_plateau), method))

        open_b = getattr(self.config, "open_boundary", True)
        if not open_b:
            # Closed chain: find true steady-state
            N_eq, _ = self.solve_equilibrium(N0, method)
            return self.heat_release(N_eq)
        else:
            # Open chain: peak Q(t) = secular equilibrium
            _, q = self.find_plateau(N0, method=method, n_scan=n_scan,
                                      t_min=t_min, t_max=t_max)
            return q

    # ── time-grid evolution ──
    def evolve(self, N0, times, method: str = "cram16"):
        """
        Runs the solution over a list of times.
        Returns (times, trajectories), where trajectories[k] is the N vector at times[k].
        Convenient for plotting and export.
        """
        N0 = np.asarray(N0, dtype=float)
        traj = []
        for t in times:
            traj.append(self.solve(N0, float(t), method) if t > 0 else N0.copy())
        return list(times), traj

    def heat_curve(self, trajectories):
        """Heat release curve [MeV/s] for each moment in the trajectory."""
        return [self.heat_release(N) for N in trajectories]

    # ── output / export ──
    def names(self) -> List[str]:
        return [getattr(i, "name", str(i)) for i in self.isotopes]

    def export_txt(self, filename: str, mode: str = "total"):
        names = self.names()
        N = len(names)
        with open(filename, "w", encoding="utf-8") as f:
            f.write("\t" + "\t".join(names) + "\n")
            for i in range(N):
                row = [names[i]]
                for j in range(N):
                    if mode == "rxn":
                        v = self.matrix_rxn[i, j]
                    elif mode == "decay":
                        v = self.matrix_decay[i, j]
                    else:
                        v = self.matrix_A[i, j]
                    row.append(f"{v:.4e}" if abs(v) > 1e-45 else "0")
                f.write("\t".join(row) + "\n")

    def export_evolution_txt(self, filename, times, trajectories):
        """
        Saves N(t): first column Time_sec, then isotope columns.
        times — list of times; trajectories — list of N vectors of equal length.
        """
        names = self.names()
        with open(filename, "w", encoding="utf-8") as f:
            f.write("Time_sec\t" + "\t".join(names) + "\n")
            for t, vec in zip(times, trajectories):
                f.write(f"{t:.6e}\t" + "\t".join(f"{v:.6e}" for v in vec) + "\n")

    def validate(self) -> bool:
        if self.matrix_A is None:
            return False
        diag_ok = np.all(np.diag(self.matrix_A) <= 1e-12)
        col_ok = np.all(np.sum(self.matrix_A, axis=0) <= 1e-9)
        return bool(diag_ok and col_ok)


# =====================================================================
# Equilibrium heat release curve Q_eq(Phi)
# =====================================================================

def equilibrium_heat_vs_flux(isotopes, config_template, N0, fluxes,
                             method: str = "cram16", n_scan: int = 80,
                             t_plateau=None):
    """
    Builds the equilibrium heat density as a function of neutron flux.

    For EACH flux the matrix is rebuilt (reaction rates proportional to Phi),
    then Q at the plateau is computed.

    Parameters
    ----------
    isotopes        — isotope list;
    config_template — BurnupConfig (flux is ignored, taken from fluxes);
    N0              — initial concentrations;
    fluxes          — list of fluxes Phi [n/cm²/s];
    method          — solver;
    n_scan          — log-grid points for auto plateau search;
    t_plateau       — plateau time [s], one value for all fluxes.
                      None -> auto-search via find_plateau() per flux.
                      May be float or a list/array of length len(fluxes).

    Returns
    -------
    (fluxes, heats) — fluxes and heat densities [MeV/s].
    """

    import copy
    fluxes = list(fluxes)
    # t_plateau may be: None / a single number / a list of numbers
    if t_plateau is None:
        t_list = [None] * len(fluxes)
    elif np.ndim(t_plateau) == 0:
        t_list = [float(t_plateau)] * len(fluxes)
    else:
        t_list = [float(t) for t in t_plateau]
        if len(t_list) != len(fluxes):
            raise ValueError("len(t_plateau) must equal len(fluxes).")

    heats = []
    for phi, tp in zip(fluxes, t_list):
        cfg = copy.copy(config_template)
        cfg.flux = float(phi)
        mtx = BurnupMatrix(isotopes, cfg)
        mtx.build()
        q = mtx.equilibrium_heat(N0, method=method, t_plateau=tp, n_scan=n_scan)
        heats.append(q)
    return fluxes, heats
