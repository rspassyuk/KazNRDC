"""
core/service.py
===============
Application facade that wires the GUI to the domain + data layers.

This is the SINGLE place the PyQt5 models talk to. It owns:
  * path / provider initialisation (xsdir, ENDF libraries, MACS, TALYS, user db);
  * per-isotope inspection for the CoreNucleo panel (decay + neutron + all MACS);
  * the burnup workspace (isotope universe → matrix → evolution) shared by the
    IsotopeChart, Results and Graph panels.

It is deliberately Qt-free: it returns plain dicts / core dataclasses so the
models stay pure UI. A module-level singleton (init_service / get_service) lets
every dynamically-loaded model reach the same instance without threading a
reference through the launcher's model-loading contract.

Mirrors the initialisation used by examples/sprocess_termination.py and test.py.
"""

from __future__ import annotations


import pickle
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from .entities import ChemTable, Isotope
from .registry import IsotopeBuilder
from .burnup import BurnupConfig, BurnupMatrix

from nuclear_data import EndfPathHolder, NuclearDataProvider
from nuclear_data.endf_reader import EndfReader, _HAS_ENDF
from nuclear_data.macs_reader import MacsReader, MultiMacsReader
from nuclear_data.talys_reader import NreacTalys


Meta = Optional[Union[int, str]]
ProgressCb = Optional[Callable[[float, str], None]]

# dec-082_Pb_206.endf  /  dec-082_Pb_206m1.endf
_DECAY_RE = re.compile(
    r"^dec-(?P<Z>\d{3})_(?P<X>[A-Za-z]{1,2})_(?P<A>\d{3})(?P<meta>m\d+)?\.endf$"
)
# "Pb-206", "Pb206", "U-235m1", "U 235 m1"
_LABEL_RE = re.compile(
    r"^\s*(?P<X>[A-Za-z]{1,2})\s*-?\s*(?P<A>\d{1,3})\s*(?:m(?P<meta>\d+))?\s*$"
)
# "Tl206..Tl216" / "Tl-206..Tl-216" — inclusive mass-number range, same element.
_RANGE_RE = re.compile(
    r"^\s*(?P<X1>[A-Za-z]{1,2})\s*-?\s*(?P<A1>\d{1,3})\s*\.\.\s*"
    r"(?P<X2>[A-Za-z]{1,2})\s*-?\s*(?P<A2>\d{1,3})\s*$"
)

# Libraries that carry a rawmacs.txt column and a MacsReader alias.
_RAWMACS_ALIASES = ["Endfb7", "Endfb6", "Jeff", "Jendl"]
# MT channels that actually appear in rawmacs.txt.
_RAWMACS_MTS = [102, 16, 18, 103, 107]
# kT at which "the" MACS is reported throughout the app (s-process peak).
MACS_KT_KEV = 30.0




# Reader decay-mode string  ->  conventional symbol.
_DECAY_SYMBOLS = {
    "beta-": "β⁻",
    "beta+": "β⁺",
    "ec/beta+": "β⁺/EC",
    "ec": "EC",
    "it": "IT",
    "alpha": "α",
    "n": "n",
    "p": "p",
    "sf": "SF",
    "unknown": "?",
    "stable": "stable",
}

# MT  ->  conventional neutron-reaction symbol (superset of EndfReader.MT_NAMES).
_MT_SYMBOLS = {
    1: "(n,total)", 2: "(n,el)", 4: "(n,inl)", 16: "(n,2n)", 17: "(n,3n)",
    18: "(n,f)", 91: "(n,n')", 102: "(n,γ)", 103: "(n,p)", 104: "(n,d)",
    105: "(n,t)", 106: "(n,³He)", 107: "(n,α)",
}

# Net (dZ, dA) one transition of each channel applies to a parent (Z, A) to
# reach its product — used for chart reachability (Task 5). Pure nuclide
# bookkeeping over the chart of nuclides, with no cross-section evaluation, so
# the selection shape can be computed before any burnup universe is built.
# Channels with no single daughter (sf, n-induced fission, unknown) are absent.
CHANNEL_DELTAS = {
    # decay modes (chart UI names)
    "beta-":    (+1,  0),
    "ec/beta+": (-1,  0),
    "alpha":    (-2, -4),
    "n":        (0,  -1),   # neutron emission
    "p":        (-1, -1),   # proton emission
    "it":       (0,   0),   # isomeric transition (same nuclide)
    # neutron-induced reactions
    "(n,γ)":    (0,  +1),   # radiative capture
    "(n,α)":    (-2, -3),   # absorb n (A+1), emit α (Z-2, A-4)
    "(n,p)":    (-1,  0),   # absorb n (A+1), emit p (Z-1, A-1)
    "(n,2n)":   (0,  -1),   # absorb n (A+1), emit 2n (A-2)
}


def decay_symbol(mode: Optional[str]) -> str:
    """Map a reader decay-mode string ('beta-', 'beta-, alpha') to symbols."""
    if not mode:
        return "?"
    parts = [p.strip().lower() for p in str(mode).split(",")]
    return ", ".join(_DECAY_SYMBOLS.get(p, p) for p in parts)


def mt_symbol(mt: int) -> str:
    return _MT_SYMBOLS.get(int(mt), f"MT={mt}")


@dataclass
class RunResult:
    """Output of a burnup evolution — shared by Results and Graph panels."""
    names: List[str]
    times: List[float]                 # seconds
    trajectories: List[np.ndarray]     # one N-vector per time
    heat: List[float]                  # MeV/s per initial nucleus, per time
    isotopes: List[Isotope] = field(default_factory=list)
    initial: Dict[str, float] = field(default_factory=dict)
    flux: float = 0.0
    energy_eV: float = 0.0253
    t_end_s: float = 0.0
    method: str = "cram16"
    # Requested target isotopes (normalized labels). Empty → no restriction;
    # otherwise the run used the minimal source→target sub-network and the UI
    # should display only these isotopes.
    targets: List[str] = field(default_factory=list)

    exposure_metadata: dict = field(default_factory=dict)
    capture_sigma_barn: Optional[np.ndarray] = None
    exposure_result: object = None

    def final(self) -> Dict[str, float]:
        if not self.trajectories:
            return dict(self.initial)
        last = self.trajectories[-1]
        return {n: float(last[i]) for i, n in enumerate(self.names)}


@dataclass
class EqResult:
    """Output of an equilibrium scan (flux sweep) — shared by Results and Graph panels."""
    names: List[str]
    flux_list: List[float]             # neutron flux points [n/cm²/s]
    time_list: List[float]             # irradiation time per flux point [s]
    N_eq: List[np.ndarray]            # equilibrium concentration vector per flux point
    Q_eq: List[float]                  # total heat at equilibrium [MeV/s per initial nucleus]
    isotopes: List[Isotope] = field(default_factory=list)
    initial: Dict[str, float] = field(default_factory=dict)
    energy_eV: float = 0.0253
    method: str = "cram16"


class NuclearService:
    """High-level facade over nuclear_data + core. One instance per app run."""

    def __init__(self):
        self.project_root: Optional[Path] = None
        self.xsdir: Optional[Path] = None
        self.endf6_root: Optional[Path] = None
        self.provider: Optional[NuclearDataProvider] = None

        self.decay_library = "ENDFB-VIII.0"
        self.neutron_library = "ENDFB-VIII.0"
        self.macs_library = "Endfb7"

        self._dec_dir: Optional[str] = None
        self._neu_dir: Optional[str] = None
        self._macs_file: Optional[str] = None
        self._endfb71_dir: Optional[str] = None
        self._eaf2010_dir: Optional[str] = None
        self._talys_dir: Optional[str] = None

        # workspace (burnup)
        self.isotopes: List[Isotope] = []
        self.matrix: Optional[BurnupMatrix] = None
        self.config: Optional[BurnupConfig] = None
        self._build_args: Optional[dict] = None     # last build_universe arguments
        self._build_energy_eV: Optional[float] = None
        self._manual_prefer_macs = False
        self.last_run: Optional[RunResult] = None
        self.last_eq: Optional[EqResult] = None
        self.last_rho: float = 0.0
        self.last_cycles: list = []
        self.last_sensitivity = None   # core.sensitivity.SensitivityResult

        # Isotopes the user excluded from calculations (Task 6); dropped from
        # every range/label build. Keyed by (Z, A), folding metastable states.
        self.excluded_zas: set = set()

        # Manual decay links added by the user (persist across matrix rebuilds)
        self._manual_decay_links: list = []

        # cross-model event bus
        self._matrix_built_callbacks: List[Callable] = []
        self._cycles_found_callbacks: List[Callable] = []
        self._theme_changed_callbacks: List[Callable] = []
        self.theme: str = "dark"   # current global UI theme ("dark" | "light")

        # nuclide-chart caches (full chart of nuclides)
        self._chart_cells: Optional[List[dict]] = None
        self._chart_modes: Optional[Dict[Tuple[int, int], Optional[str]]] = None

    # ──────────────────────────────────────────────────────────────────
    # Initialisation
    # ──────────────────────────────────────────────────────────────────
    @property
    def initialized(self) -> bool:
        return self.provider is not None

    def init(self, project_root, decay_library="ENDFB-VIII.0",
             neutron_library="ENDFB-VIII.0", macs_library="Endfb7",
             progress: ProgressCb = None) -> None:
        """
        Discover xsdir/user_db under project_root and build the provider with
        every available source. `progress(frac, message)` drives the splash.
        """
        def step(frac, msg):
            if progress:
                progress(frac, msg)

        # The 'endf' library warns about covariance sections it cannot parse
        # (MF=32/33, MT=151). These are irrelevant here and would spam the
        # in-app terminal on every isotope read.
        import warnings
        warnings.filterwarnings("ignore", category=UserWarning, module=r"endf\..*")

        self.project_root = Path(project_root).expanduser().resolve()
        self.decay_library = decay_library
        self.neutron_library = neutron_library
        self.macs_library = macs_library

        step(0.05, "Initializing databases…")
        info = EndfPathHolder.init_xsdir(str(self.project_root / "xsdir"))
        self.xsdir = Path(info["xsdir"])
        self.endf6_root = Path(info["endf-6"])

        step(0.15, "Resolving ENDF libraries…")
        self._dec_dir = EndfPathHolder.set_DEFAULT_DEC_DIR(decay_library)
        self._neu_dir = EndfPathHolder.set_DEFAULT_NEU_DIR(neutron_library)
        self._macs_file = EndfPathHolder.set_DEFAULT_MACS_DIR(macs_library)

        # Optional multi-kT MACS + TALYS (auto-discovered by init_xsdir).
        self._endfb71_dir = info.get("ENDFB71_MACS")
        self._eaf2010_dir = info.get("EAF2010_MACS")
        self._talys_dir = info.get("TALYS")

        user_db = self.project_root / "user_db"
        user_dir = EndfPathHolder.set_USER_DB_DIR(str(user_db))

        step(0.30, "Building data provider…")
        self.provider = NuclearDataProvider(
            user_dir=user_dir,
            dec_dir=self._dec_dir,
            neu_dir=self._neu_dir,
            macs_file=self._macs_file,
            macs_lib=macs_library,
            endfb71_macs_dir=self._endfb71_dir,
            eaf2010_macs_dir=self._eaf2010_dir,
            talys_dir=self._talys_dir,
            dec_priority=["user", "endf6"],
            macs_priority=["user", "endfb71", "eaf2010", "rawmacs", "talys", "endf_computed"],
        )

        step(0.45, "Indexing MACS tables…")
        # Touch the multi-kT readers so their (cached) parse happens here, behind
        # the splash, instead of on first isotope click.
        try:
            if self._endfb71_dir:
                MultiMacsReader(self._endfb71_dir)
            if self._eaf2010_dir:
                MultiMacsReader(self._eaf2010_dir)
            if self._macs_file:
                MacsReader(self._macs_file, macs_library)
        except Exception:
            pass

        step(0.55, "Databases ready.")

    def info(self) -> Dict[str, object]:
        return {
            "initialized": self.initialized,
            "decay_library": self.decay_library,
            "neutron_library": self.neutron_library,
            "macs_library": self.macs_library,
            "xsdir": str(self.xsdir) if self.xsdir else None,
            "has_endf": _HAS_ENDF,
        }

    def decay_libraries(self) -> List[str]:
        """ENDF-6 libraries available for decay (sub-folders of xsdir/endf-6)."""
        if not self.endf6_root or not self.endf6_root.is_dir():
            return []
        return sorted(p.name for p in self.endf6_root.iterdir()
                      if p.is_dir() and (p / "decay").is_dir())

    def switch_decay_library(self, library: str) -> None:
        """Switch the active ENDF-6 decay library and invalidate cached universe.

        Rebuilds the data provider so the new library is used on the next
        build_universe() / build_matrix() call.
        """
        if not self.endf6_root:
            raise RuntimeError("Service not initialised. Call init() first.")
        dec_path = self.endf6_root / library / "decay"
        if not dec_path.is_dir():
            raise ValueError(
                f"Library '{library}' has no decay sub-folder under {self.endf6_root}")
        self._dec_dir = str(dec_path)
        EndfPathHolder.set_DEFAULT_DEC_DIR(library)
        self.decay_library = library

        user_db = (self.project_root / "user_db") if self.project_root else None
        user_dir = EndfPathHolder.set_USER_DB_DIR(str(user_db)) if user_db else None
        self.provider = NuclearDataProvider(
            user_dir=user_dir,
            dec_dir=self._dec_dir,
            neu_dir=self._neu_dir,
            macs_file=self._macs_file,
            macs_lib=self.macs_library,
            endfb71_macs_dir=self._endfb71_dir,
            eaf2010_macs_dir=self._eaf2010_dir,
            talys_dir=self._talys_dir,
            dec_priority=["user", "endf6"],
            macs_priority=["user", "endfb71", "eaf2010", "rawmacs",
                           "talys", "endf_computed"],
        )
        # Invalidate everything derived from the decay library
        self.isotopes = []
        self.matrix = None
        self.config = None
        self._build_args = None
        self._chart_cells = None
        self._chart_modes = None
        print(f"[service] Decay library switched to: {library}")

    # ──────────────────────────────────────────────────────────────────
    # CoreNucleo inspection  (Task 1) — multi-library, all MACS sources
    # ──────────────────────────────────────────────────────────────────
    def _lib_dir(self, library: str, kind: str) -> Optional[Path]:
        if not self.endf6_root:
            return None
        d = self.endf6_root / library / kind
        return d if d.is_dir() else None

    def scan_isotopes_for_element(self, library: str, Z: int, X: str
                                  ) -> List[Tuple[int, Meta]]:
        """[(A, meta), …] for an element from a library's dec-*.endf files."""
        dec = self._lib_dir(library, "decay")
        if dec is None:
            return []
        z3 = f"{int(Z):03d}"
        Xn = str(X).strip()
        out: List[Tuple[int, Meta]] = []
        for fp in dec.iterdir():
            m = _DECAY_RE.match(fp.name)
            if not m or m.group("Z") != z3 or m.group("X") != Xn:
                continue
            meta_s = m.group("meta")
            meta: Meta = int(meta_s[1:]) if meta_s else None
            out.append((int(m.group("A")), meta))
        out.sort(key=lambda t: (t[0], 0 if t[1] in (None, 0, "") else int(t[1])))
        return out

    def read_decay(self, library: str, Z: int, X: str, A: int, meta: Meta = None
                   ) -> Optional[dict]:
        """Decay parameters + channels for one isotope from ONE ENDF-6 library."""
        dec = self._lib_dir(library, "decay")
        if dec is None or not _HAS_ENDF:
            return None
        try:
            r = EndfReader(Z, X, A, meta, dec_dir=str(dec))
        except Exception:
            return None
        if not r.has_decay:
            return None
        dd = r.get_decay_data()
        mass = r.get_mass(to_kg=False)
        channels = [{
            "mode": ch.mode,
            "symbol": decay_symbol(ch.mode),
            "branch": ch.branch,
            "Q_MeV": ch.Q_MeV,
            "dz": ch.dz, "da": ch.da, "rfs": ch.rfs,
        } for ch in dd.channels]
        return {
            "library": library,
            "half_life_s": dd.half_life_s,
            "decay_constant": dd.decay_constant,
            "is_stable": dd.is_stable,
            "mass_amu": mass,
            "channels": channels,
            "source": dd.source,
        }

    def read_neutron(self, library: str, Z: int, X: str, A: int, meta: Meta = None
                     ) -> List[dict]:
        """
        Per-MT neutron-reaction data for one isotope from ONE library:
        channel, σ(E) summary, reaction Q, isomeric branch, MACS@30 keV (from
        that library's own spectrum). Empty when the library has no n-file
        (e.g. JEFF here) — callers fall back to the MACS sources below.
        """
        neu = self._lib_dir(library, "neutrons")
        if neu is None or not _HAS_ENDF:
            return []
        try:
            r = EndfReader(Z, X, A, meta, neu_dir=str(neu))
        except Exception:
            return []
        if not r.has_neutron:
            return []
        out: List[dict] = []
        for mt in r.available_mts():
            rd = r.get_reaction_data(mt)
            if rd is None:
                continue
            macs30 = r.compute_macs(mt, MACS_KT_KEV, "keV")
            branches = r.get_isomeric_branches(mt, MACS_KT_KEV * 1e3)
            out.append({
                "mt": mt,
                "name": mt_symbol(mt),
                "Q_MeV": rd.Q_MeV,
                "n_points": len(rd.spectrum),
                "sigma_thermal_barn": r.get_sigma_at(mt, 0.0253, "eV"),
                "sigma_30keV_barn": r.get_sigma_at(mt, MACS_KT_KEV, "keV"),
                "sigma_1MeV_barn": r.get_sigma_at(mt, 1.0, "MeV"),
                "macs30_mb": macs30.sigma_mb if macs30 else None,
                "branches": [{"LFS": b.LFS, "yield": b.yield_} for b in branches],
                "spectrum": rd.spectrum,
            })
        return out

    def read_macs_sources(self, Z: int, X: str, A: int, meta: Meta = None) -> dict:
        """
        ALL MACS available in xsdir/MACS for an isotope, independent of the
        ENDF-6 libraries selected for decay:
          rawmacs   — every library/MT present in rawmacs.txt
          endfb71   — ENDF/B-VII.1 multi-kT (MT102)
          eaf2010   — EAF-2010 multi-kT (MT102)
          talys     — TALYS σ(E) spectra (na/ng/np) + single MACS@30 keV (MT102)
        """
        sym = str(X).capitalize()
        iso_key = f"{int(Z)}-{sym}-{int(A)}"
        result: dict = {"rawmacs": {}, "endfb71": None, "eaf2010": None, "talys": None}

        # rawmacs — per library, per MT
        if self._macs_file:
            for alias in _RAWMACS_ALIASES:
                try:
                    rdr = MacsReader(self._macs_file, alias)
                except Exception:
                    continue
                lib_out: Dict[int, dict] = {}
                for mt in _RAWMACS_MTS:
                    pts = rdr.get_all_points(iso_key, mt)
                    if not pts:
                        continue
                    at30 = rdr.get_macs_at(iso_key, mt, MACS_KT_KEV)
                    lib_out[mt] = {
                        "macs30_mb": at30.sigma_mb if at30 else None,
                        "points": [(p.kT_keV, p.sigma_mb) for p in pts],
                    }
                if lib_out:
                    result["rawmacs"][rdr._library] = lib_out

        # multi-kT tables (MT102)
        for key, directory in (("endfb71", self._endfb71_dir),
                               ("eaf2010", self._eaf2010_dir)):
            if not directory:
                continue
            try:
                mr = MultiMacsReader(directory)
            except Exception:
                continue
            pts = mr.get_all_points(iso_key, 102)
            if pts:
                at30 = mr.get_macs_at(iso_key, 102, MACS_KT_KEV)
                result[key] = {
                    "macs30_mb": at30.sigma_mb if at30 else None,
                    "points": [(p.kT_keV, p.sigma_mb) for p in pts],
                }

        # TALYS — recommended spectra + single 30 keV MACS
        if self._talys_dir:
            try:
                t = NreacTalys(self._talys_dir, Z, sym, A, meta)
            except Exception:
                t = None
            if t is not None:
                spectra: Dict[int, dict] = {}
                for mt in t.available_mts():
                    spec = t.get_cross_section(mt)
                    if spec:
                        spectra[mt] = {"name": mt_symbol(mt), "n_points": len(spec),
                                       "spectrum": spec}
                m30 = t.get_macs_30keV(102)
                if spectra or m30:
                    result["talys"] = {
                        "spectra": spectra,
                        "macs30_mb": m30.sigma_mb if m30 else None,
                    }
        return result

    # ──────────────────────────────────────────────────────────────────
    # Burnup workspace  (Tasks 2–4)
    # ──────────────────────────────────────────────────────────────────
    def _scan_specs(self, z_from: int, z_to: int,
                    n_from: Optional[int] = None, n_to: Optional[int] = None
                    ) -> List[tuple]:
        """
        (Z, X, A, meta) specs from the default library's decay files.

        Optional n_from/n_to restrict the neutron number N = A - Z, so the two
        ranges together describe a bounding box on the chart of nuclides.
        """
        if not self._dec_dir:
            return []
        specs: List[tuple] = []
        for fp in Path(self._dec_dir).iterdir():
            m = _DECAY_RE.match(fp.name)
            if not m:
                continue
            Z = int(m.group("Z"))
            if not (z_from <= Z <= z_to):
                continue
            A = int(m.group("A"))
            N = A - Z
            if n_from is not None and N < n_from:
                continue
            if n_to is not None and N > n_to:
                continue
            if (Z, A) in self.excluded_zas:          # Task 6: drop excluded
                continue
            meta_s = m.group("meta")
            meta: Meta = int(meta_s[1:]) if meta_s else None
            specs.append((Z, m.group("X"), A, meta))
        specs.sort(key=lambda s: (s[0], s[2], 0 if s[3] in (None, 0) else int(s[3])))
        return specs

    def count_specs(self, z_from: int, z_to: int,
                    n_from: Optional[int] = None, n_to: Optional[int] = None) -> int:
        """How many isotopes a build over this range would assemble (fast)."""
        return len(self._scan_specs(z_from, z_to, n_from, n_to))

    def build_universe(self, z_from: int = 80, z_to: int = 84,
                       n_from: Optional[int] = None, n_to: Optional[int] = None,
                       E_eV: float = MACS_KT_KEV * 1e3,
                       mt_filter: Optional[List[int]] = None,
                       prefer_macs: bool = False,
                       progress: ProgressCb = None) -> List[Isotope]:
        """
        Build the isotope set (and wire decay/reaction products) over a Z window,
        optionally restricted to an N band (bounding box on the chart).
        prefer_macs=True uses MACS tables for σ (astrophysics mode).
        Cached: identical arguments return the previously built universe.
        """
        if not self.initialized:
            raise RuntimeError("Service not initialised. Call init() first.")
        mt_filter = mt_filter or [102]
        args = {
            "z_from": int(z_from), "z_to": int(z_to),
            "n_from": None if n_from is None else int(n_from),
            "n_to": None if n_to is None else int(n_to),
            "E_eV": float(E_eV), "mt_filter": tuple(mt_filter),
            "prefer_macs": bool(prefer_macs),
            "excluded": frozenset(self.excluded_zas),   # Task 6: invalidate cache
        }
        if args == self._build_args and self.isotopes:
            return self.isotopes

        specs = self._scan_specs(z_from, z_to, n_from, n_to)
        if progress:
            progress(0.05, f"Assembling {len(specs)} isotopes…")
        builder = IsotopeBuilder(self.provider)
        self.isotopes = builder.build_range(
            specs, E_eV=E_eV, mt_filter=mt_filter, prefer_macs=prefer_macs)
        self._build_args = args
        self._build_energy_eV = float(E_eV)
        # workspace changed → invalidate downstream artefacts
        self.matrix = None
        return self.isotopes

    def build_universe_from_labels(self, labels: List[str],
                                   E_eV: float = MACS_KT_KEV * 1e3,
                                   mt_filter: Optional[List[int]] = None,
                                   prefer_macs: bool = False,
                                   progress: ProgressCb = None) -> List[Isotope]:
        """
        Build the isotope set from an explicit list of isotope labels
        (e.g. ['Tl-206', 'Pb-207', 'Pb-208']). Used for manual range input.
        """
        if not self.initialized:
            raise RuntimeError("Service not initialised. Call init() first.")
        specs: List[tuple] = []
        for label in labels:
            parsed = self.parse_label(label)
            if parsed:
                specs.append(parsed)
        specs = [s for s in specs if (s[0], s[2]) not in self.excluded_zas]  # Task 6
        if not specs:
            raise ValueError("No valid isotope labels in the provided list.")
        if progress:
            progress(0.05, f"Assembling {len(specs)} isotopes…")
        builder = IsotopeBuilder(self.provider)
        self.isotopes = builder.build_range(
            specs, E_eV=E_eV, mt_filter=mt_filter, prefer_macs=prefer_macs)
        self._build_args = None          # manual build — do not cache by range
        self._build_energy_eV = float(E_eV)
        self.matrix = None
        return self.isotopes

    def edit_workspace_isotopes(self, add_labels=(), remove_labels=()):
        """Edit membership independently of decay and reaction channel filters."""
        if not self.initialized:
            raise RuntimeError("Load nuclear data before editing the isotope set.")
        additions = []
        for label in add_labels:
            spec = self.parse_label(label)
            if spec is None:
                raise ValueError(f"Invalid isotope label: {label}")
            additions.append((self.normalize_label(label), spec))
        removed = {self.normalize_label(label) for label in remove_labels}
        labels = [iso.name for iso in self.isotopes if iso.name not in removed]
        for label, spec in additions:
            if label not in labels:
                labels.append(label)
        if not labels:
            raise ValueError("Keep at least one isotope in the calculation set.")
        prefer = (self._build_args or {}).get("prefer_macs", self._manual_prefer_macs)
        current = {iso.name: iso for iso in self.isotopes}
        specs = [self.parse_label(label) for label in labels if label not in current]
        if specs:
            builder = IsotopeBuilder(self.provider)
            loaded = builder.build_range(
                specs, E_eV=self._build_energy_eV or MACS_KT_KEV * 1e3,
                mt_filter=None, prefer_macs=prefer)
            current.update((iso.name, iso) for iso in loaded)
        updated = [current[label] for label in labels]
        IsotopeBuilder.link_products(updated)
        self.isotopes = updated
        for label, spec in additions:
            self.excluded_zas.discard((spec[0], spec[2]))
        self._manual_prefer_macs = prefer
        self._build_args = None
        self.matrix = None
        return self.isotopes

    def get_universe(self, default_range=(80, 84)) -> List[Isotope]:
        """Current isotope set; builds a sensible default on first access."""
        if not self.isotopes and self.initialized:
            self.build_universe(default_range[0], default_range[1])
        return self.isotopes

    def universe_labels(self) -> List[str]:
        return [iso.name for iso in self.isotopes]

    # ──────────────────────────────────────────────────────────────────
    # Excluded isotopes (Task 6) — dropped from every range/label build
    # ──────────────────────────────────────────────────────────────────
    def get_excluded(self) -> set:
        """Copy of the excluded (Z, A) set."""
        return set(self.excluded_zas)

    def is_excluded(self, za) -> bool:
        return tuple(za) in self.excluded_zas

    def add_excluded(self, za) -> None:
        self.excluded_zas.add(tuple(za))

    def remove_excluded(self, za) -> None:
        self.excluded_zas.discard(tuple(za))

    def set_excluded(self, zas) -> None:
        self.excluded_zas = {tuple(za) for za in zas}

    def clear_excluded(self) -> None:
        self.excluded_zas.clear()

    # ──────────────────────────────────────────────────────────────────
    # Full chart of nuclides (display layer — independent of the burnup set)
    # ──────────────────────────────────────────────────────────────────
    def chart_nuclides(self) -> List[dict]:
        """
        Every available nuclide as a chart cell, one per (Z, A) — fast (built
        from decay filenames, no ENDF parsing). Metastable states are folded
        into their ground-state cell (has_meta flag).
        Each cell: {Z, X, A, N, metas:set, has_meta:bool}.
        """
        if self._chart_cells is not None:
            return self._chart_cells
        cells: Dict[Tuple[int, int], dict] = {}
        for (Z, X, A, meta) in self._scan_specs(1, 118):
            key = (Z, A)
            c = cells.get(key)
            if c is None:
                c = {"Z": Z, "X": X, "A": A, "N": A - Z, "metas": set(), "has_meta": False}
                cells[key] = c
            if meta not in (None, 0, ""):
                c["metas"].add(meta)
                c["has_meta"] = True
        self._chart_cells = sorted(cells.values(), key=lambda c: (c["Z"], c["A"]))
        return self._chart_cells

    def reachable_cells(self, seed, channels, transitive: bool = True,
                        max_cells: int = 20000) -> set:
        """
        (Z, A) cells reachable from `seed` by repeatedly applying the active
        `channels`' (dZ, dA) transforms, restricted to nuclides that actually
        exist on the chart of nuclides (Task 5).

        `seed`     iterable of (Z, A) start cells (the selected range).
        `channels` iterable of channel names — decay UI names ("beta-",
                   "alpha", …) and/or neutron-reaction names ("(n,γ)", …).

        Returns a set containing the seed cells present on the chart plus every
        cell reachable through the active channels. This is pure graph
        reachability (no σ evaluation), so the selection outline can be drawn
        before any burnup universe is built. With transitive=False only the
        immediate one-step products are added.
        """
        chart = {(c["Z"], c["A"]) for c in self.chart_nuclides()}
        deltas = [CHANNEL_DELTAS[c] for c in channels if c in CHANNEL_DELTAS]
        result = {za for za in seed if za in chart}
        if not deltas or not result:
            return result
        frontier = list(result)
        while frontier:
            nxt = []
            for (Z, A) in frontier:
                for (dZ, dA) in deltas:
                    nb = (Z + dZ, A + dA)
                    if nb in chart and nb not in result:
                        result.add(nb)
                        nxt.append(nb)
                        if len(result) >= max_cells:
                            return result
            if not transitive:
                break
            frontier = nxt
        return result

    def _chart_cache_path(self) -> Optional[Path]:
        if not self.project_root:
            return None
        return self.project_root / ".cache" / f"nuclide_modes_{self.decay_library}.pkl"

    def _chart_signature(self):
        n = len(self.chart_nuclides())
        return (self.decay_library, n)

    def chart_decay_modes(self, progress: ProgressCb = None,
                          should_cancel: Optional[Callable[[], bool]] = None
                          ) -> Dict[Tuple[int, int], Optional[str]]:
        """
        {(Z, A): primary_decay_mode}, where the value is None for a stable
        nuclide and a reader mode string otherwise ("beta-", "alpha", …).

        Heavy on first run (~5 s over the whole chart); cached to a pickle and
        returned instantly afterwards. Intended to run on a worker thread.
        """
        if self._chart_modes is not None:
            return self._chart_modes

        cache = self._chart_cache_path()
        sig = self._chart_signature()
        if cache and cache.is_file():
            try:
                pl = pickle.loads(cache.read_bytes())
                if pl.get("sig") == list(sig) or pl.get("sig") == sig:
                    self._chart_modes = {tuple(k): v for k, v in pl["modes"].items()}
                    return self._chart_modes
            except Exception:
                pass

        cells = self.chart_nuclides()
        modes: Dict[Tuple[int, int], Optional[str]] = {}
        total = max(1, len(cells))
        cancelled = False
        for i, c in enumerate(cells):
            if should_cancel and should_cancel():
                cancelled = True
                break
            if progress and (i % 100 == 0 or i == total - 1):
                progress(i / total, f"Indexing nuclide chart… {i}/{total}")
            d = self.read_decay(self.decay_library, c["Z"], c["X"], c["A"], None)
            if d is None:
                modes[(c["Z"], c["A"])] = "unknown"
            elif d["is_stable"]:
                modes[(c["Z"], c["A"])] = None
            else:
                chans = d.get("channels") or []
                if chans:
                    top = max(chans, key=lambda ch: ch.get("branch", 0.0))
                    modes[(c["Z"], c["A"])] = top["mode"]
                else:
                    modes[(c["Z"], c["A"])] = "unknown"

        self._chart_modes = modes
        if cache and not cancelled:
            try:
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_bytes(pickle.dumps({"sig": sig, "modes": modes},
                                               protocol=pickle.HIGHEST_PROTOCOL))
            except Exception:
                pass
        return modes

    def _ensure_universe(self, energy_eV: float,
                         mt_filter: Optional[List[int]] = None,
                         progress: ProgressCb = None) -> None:
        """Make sure self.isotopes exists and is evaluated at `energy_eV`.
        Cross sections are evaluated at build time, so a different energy forces
        a rebuild over the same range."""
        if not self.isotopes:
            self.build_universe(E_eV=energy_eV, mt_filter=mt_filter, progress=progress)
        elif self._build_energy_eV != float(energy_eV) and self._build_args:
            a = self._build_args
            self.build_universe(a["z_from"], a["z_to"], a["n_from"], a["n_to"],
                                E_eV=energy_eV, mt_filter=mt_filter,
                                prefer_macs=a.get("prefer_macs", False),
                                progress=progress)
        elif self._build_energy_eV != float(energy_eV) and self.isotopes:
            self.build_universe_from_labels(
                self.universe_labels(), E_eV=energy_eV, mt_filter=mt_filter,
                prefer_macs=self._manual_prefer_macs, progress=progress)

    def build_matrix(self, flux: float = 0.0, energy_eV: float = MACS_KT_KEV * 1e3,
                     mt_filter: Optional[List[int]] = None,
                     decay_mode_filter: Optional[List[str]] = None,
                     open_boundary: bool = False,
                     open_boundary_mts: Optional[List[int]] = None,
                     hl_preset: str = "all",
                     isotopes: Optional[List[Isotope]] = None,
                     manual_decay_links: Optional[list] = None,
                     progress: ProgressCb = None) -> BurnupMatrix:
        """Assemble a BurnupMatrix over the current universe.

        open_boundary_mts — list of MT codes whose reactions with out-of-network
            products are treated as sinks (open boundary). All other MTs are
            closed (reaction is disabled for edge isotopes). When None, the
            global open_boundary bool is used instead.
        isotopes — explicit isotope subset to assemble the matrix over (e.g. a
            minimal sub-network). When None, the full current universe is used.
            A subset built from objects of the current universe stays numerically
            exact: reactions whose product is outside the subset become sinks.
        """
        if progress:
            progress(0.05, "Preparing isotopes…")
        self._ensure_universe(energy_eV, mt_filter, progress)

        isos = isotopes if isotopes is not None else self.isotopes
        if progress:
            progress(0.45, f"Assembling {len(isos)}×{len(isos)} matrix…")
        cfg = BurnupConfig(flux=float(flux))
        if mt_filter is not None:
            cfg.enable_specific_reactions(mt_filter)
        if decay_mode_filter is not None:
            cfg.set_decay_mode_filter(decay_mode_filter)
        cfg.set_decay_time_filter(hl_preset if hl_preset in BurnupConfig._PRESETS else "all")
        if open_boundary_mts is not None:
            cfg.set_open_boundary_mts(open_boundary_mts)
        else:
            cfg.set_open_boundary(open_boundary)
        cfg.manual_decay_links = list(
            self._manual_decay_links if manual_decay_links is None else manual_decay_links)
        cfg.manual_decay_links = [
            dict(link, parent=self.normalize_label(link["parent"]),
                 daughter=self.normalize_label(link["daughter"]) if link.get("daughter") else "")
            for link in cfg.manual_decay_links
        ]
        mtx = BurnupMatrix(isos, cfg)
        mtx.build()
        self.config = cfg
        # Only cache as the canonical matrix when it spans the full universe.
        if isotopes is None:
            self.matrix = mtx
        return mtx

    def minimal_subnetwork(self, source_labels, target_labels,
                           mt_filter: Optional[List[int]] = None) -> List[Isotope]:
        """Minimal correct sub-network: the isotopes lying on at least one path
        from a source (start inventory) to a target.

            keep = (reachable forward from sources) ∩ (can reach a target)

        Edges are decays (always) plus the selected reaction channels, restricted
        to products that exist in the current universe. Every precursor that
        influences a target is retained, so a burnup solved over this subset is
        numerically identical to the full network for the target isotopes, while
        branches that never reach a target are dropped (the speed-up).

        Returns isotopes in universe order. Falls back to the full universe when
        sources or targets cannot be resolved.
        """
        isos = self.isotopes
        if not isos:
            return []
        by_name = {i.name: i for i in isos}
        src = {by_name[n] for n in (self.normalize_label(l) for l in source_labels)
               if n in by_name}
        tgt = {by_name[n] for n in (self.normalize_label(l) for l in target_labels)
               if n in by_name}
        if not src or not tgt:
            return list(isos)

        mtf = set(mt_filter) if mt_filter else None
        iso_set = set(isos)
        fwd: Dict[Isotope, list] = {i: [] for i in isos}
        rev: Dict[Isotope, list] = {i: [] for i in isos}
        for i in isos:
            nbrs = set()
            for d in i.getListOfDecays():
                p = getattr(d, "product", None)
                if p in iso_set:
                    nbrs.add(p)
            for rx in i.getListOfReactions():
                mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                if mtf is not None and mt not in mtf:
                    continue
                p = getattr(rx, "product", None)
                if p in iso_set:
                    nbrs.add(p)
            for p in nbrs:
                fwd[i].append(p)
                rev[p].append(i)

        def _reach(starts, adj):
            seen = set(starts)
            stack = list(starts)
            while stack:
                n = stack.pop()
                for m in adj[n]:
                    if m not in seen:
                        seen.add(m)
                        stack.append(m)
            return seen

        forward = _reach(src, fwd)
        backward = _reach(tgt, rev)
        keep = (forward & backward) | tgt   # always report the requested targets
        return [i for i in isos if i in keep]

    def run_evolution(self, initial_by_label: Dict[str, float], t_end_s: float,
                      flux: float = 0.0, energy_eV: float = MACS_KT_KEV * 1e3,
                      method: str = "cram16",
                      mt_filter: Optional[List[int]] = None,
                      decay_mode_filter: Optional[List[str]] = None,
                      hl_preset: str = "all",
                      n_points: int = 400, t_start_s: float = 1e-6,
                      open_boundary: bool = False,
                      open_boundary_mts: Optional[List[int]] = None,
                      target_labels: Optional[List[str]] = None,
                      progress: ProgressCb = None,
                      should_cancel: Optional[Callable[[], bool]] = None) -> RunResult:
        """
        Full pipeline: (re)build matrix → evolve over a geometric time grid →
        store + return a RunResult. Drives the Results table and Graph panel.

        open_boundary_mts — per-MT open boundary (see build_matrix).  When set,
            it takes precedence over the global open_boundary flag.
        target_labels — when given, the run is restricted to the minimal
            source→target sub-network (faster, same target results). The result's
            `targets` field is set so the UI can display only those isotopes.
        `progress(frac, msg)` and `should_cancel()` let a GUI worker thread show
        a progress bar and cancel a long run.
        """
        # Build the sub-network from the full universe (built/refreshed here).
        # Targets are resolved against the universe; unknown labels are dropped,
        # and if none resolve we fall back to the full network (no restriction).
        subset = None
        targets: List[str] = []
        if target_labels:
            self._ensure_universe(energy_eV, mt_filter, progress)
            present = {i.name for i in self.isotopes}
            targets = [t for t in (self.normalize_label(l) for l in target_labels)
                       if t in present]
            if targets and not self._manual_decay_links:
                subset = self.minimal_subnetwork(initial_by_label.keys(),
                                                 targets, mt_filter)

        mtx = self.build_matrix(flux=flux, energy_eV=energy_eV,
                                mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
                                open_boundary=open_boundary,
                                open_boundary_mts=open_boundary_mts,
                                hl_preset=hl_preset, isotopes=subset, progress=progress)
        names = mtx.names()
        index = {n: i for i, n in enumerate(names)}

        N0 = np.zeros(len(names))
        applied: Dict[str, float] = {}
        for label, conc in initial_by_label.items():
            nm = self.normalize_label(label)
            if nm in index:
                N0[index[nm]] += float(conc)
                applied[nm] = applied.get(nm, 0.0) + float(conc)

        from .exposure import capture_snapshot
        capture_sigma, capture_sources = capture_snapshot(names, mtx.isotopes)
        preference = (self._build_args or {}).get("prefer_macs", self._manual_prefer_macs)
        metadata = dict(
            constant_conditions=True, energy_eV=float(energy_eV), method=method,
            sigma_mode="MACS preferred with fallback" if preference else "Pointwise preferred with fallback",
            capture_sources=capture_sources, decay_library=self.decay_library,
            neutron_library=self.neutron_library,
            macs_priority=list(self.provider.get_macs_priority()) if self.provider else [],
            open_boundary=bool(open_boundary), open_boundary_mts=open_boundary_mts,
            mt_filter=mt_filter, decay_mode_filter=decay_mode_filter, hl_preset=hl_preset,
            boundary_columns=[
                names[i] for i, value in enumerate(mtx.matrix_A.sum(axis=0))
                if value < -1e-12 * max(abs(mtx.matrix_A[i, i]), 1e-300)],
            manual_decay_links=[dict(link) for link in mtx.config.manual_decay_links],
            sigma_provenance="Branch-weighted capture sum at Evolution run time")

        metadata["solver_implementation"] = (
            "Adaptive sparse CRAM-16 with step doubling"
            if method.lower() in ("cram16", "cram48") else method)
        if method.lower() == "cram48":
            metadata["solver_implementation"] += " (legacy cram48 alias)"

        # Manual evolve loop (mirrors BurnupMatrix.evolve) so we can report
        # per-step progress and honour cancellation.
        times = self._geo_times(t_end_s, n_points, t_start=t_start_s)
        traj = []
        total = max(1, len(times))
        for i, t in enumerate(times):
            if should_cancel and should_cancel():
                raise RuntimeError("Calculation cancelled.")
            if progress and (i % 4 == 0 or i == total - 1):
                progress(0.5 + 0.48 * (i / total), f"Solving time step {i + 1}/{total}…")
            if method.lower() in ("cram16", "cram48") and t > 0:
                from .cram import adaptive_cram16
                traj.append(adaptive_cram16(mtx.matrix_A, N0, float(t),
                                            should_cancel=should_cancel))
            else:
                traj.append(mtx.solve(N0, float(t), method) if t > 0 else N0.copy())
        heat = mtx.heat_curve(traj)

        self.last_run = RunResult(
            names=names, times=list(times), trajectories=traj, heat=heat,
            isotopes=list(mtx.isotopes), initial=applied, flux=float(flux),
            energy_eV=float(energy_eV), t_end_s=float(t_end_s), method=method,
            targets=targets, exposure_metadata=metadata, capture_sigma_barn=capture_sigma,
        )
        return self.last_run

    def run_sensitivity(self, config, progress: ProgressCb = None,
                        should_cancel: Optional[Callable[[], bool]] = None):
        """Run an OAT sensitivity analysis (see core.sensitivity).

        `config` is a SensitivityConfig; returns a SensitivityResult. Imported
        lazily to avoid a core import cycle.
        """
        from .sensitivity import SensitivityEngine
        self.last_sensitivity = SensitivityEngine(self).run(
            config, progress=progress, should_cancel=should_cancel)
        return self.last_sensitivity

    def run_equilibrium_scan(self, initial_by_label: Dict[str, float],
                             flux_list: List[float], time_list: List[float],
                             energy_eV: float = MACS_KT_KEV * 1e3,
                             method: str = "cram16",
                             mt_filter: Optional[List[int]] = None,
                             decay_mode_filter: Optional[List[str]] = None,
                             hl_preset: str = "all",
                             open_boundary: bool = False,
                             open_boundary_mts: Optional[List[int]] = None,
                             progress: ProgressCb = None,
                             should_cancel: Optional[Callable[[], bool]] = None) -> EqResult:
        """
        Equilibrium scan over flux_list.  For each flux[i] the system is irradiated
        for time[i] seconds and the final (equilibrium) state + heat are recorded.
        Stores result in self.last_eq and returns it.

        open_boundary_mts — per-MT open boundary (see build_matrix).
        """
        if not self.initialized:
            raise RuntimeError("Service not initialised. Call init() first.")
        mt_filter = mt_filter or [102]
        if not self.isotopes:
            raise RuntimeError("No burnup universe. Build the matrix first.")
        if len(flux_list) != len(time_list):
            raise ValueError("flux_list and time_list must have the same length.")

        # Build the initial N0 vector (consistent with current universe names).
        names = [iso.name for iso in self.isotopes]
        index = {n: i for i, n in enumerate(names)}
        N0 = np.zeros(len(names))
        applied: Dict[str, float] = {}
        for label, conc in initial_by_label.items():
            nm = self.normalize_label(label)
            if nm in index:
                N0[index[nm]] += float(conc)
                applied[nm] = applied.get(nm, 0.0) + float(conc)

        total = max(1, len(flux_list))
        N_eq_list: List[np.ndarray] = []
        Q_eq_list: List[float] = []

        for i, (flux, t_irr) in enumerate(zip(flux_list, time_list)):
            if should_cancel and should_cancel():
                raise RuntimeError("Equilibrium scan cancelled.")
            if progress and (i % max(1, total // 20) == 0 or i == total - 1):
                progress(0.1 + 0.88 * i / total,
                         f"Equilibrium scan {i + 1}/{total}  φ={flux:.2e} n/cm²/s…")
            mtx = self.build_matrix(flux=flux, energy_eV=energy_eV,
                                    mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
                                    open_boundary=open_boundary,
                                    open_boundary_mts=open_boundary_mts,
                                    hl_preset=hl_preset)
            N_i = mtx.solve(N0, float(t_irr), method) if t_irr > 0 else N0.copy()
            Q_i = float(mtx.heat_release(N_i))
            N_eq_list.append(N_i)
            Q_eq_list.append(Q_i)

        self.last_eq = EqResult(
            names=names, flux_list=list(flux_list), time_list=list(time_list),
            N_eq=N_eq_list, Q_eq=Q_eq_list,
            isotopes=list(self.isotopes), initial=applied,
            energy_eV=float(energy_eV), method=method,
        )
        return self.last_eq

    # ──────────────────────────────────────────────────────────────────
    # Cross-model event bus
    # ──────────────────────────────────────────────────────────────────

    def set_dec_priority(self, libs: List[str]) -> None:
        if self.provider:
            self.provider.set_dec_priority(libs)

    def set_macs_priority(self, libs: List[str]) -> None:
        if self.provider:
            self.provider.set_macs_priority(libs)

    def set_manual_decay_links(self, links: list) -> None:
        """Store manual decay links that are injected into every subsequent build_matrix()."""
        self._manual_decay_links = list(links)

    def get_manual_decay_links(self) -> list:
        return list(self._manual_decay_links)

    def register_matrix_built_callback(self, cb: Callable) -> None:
        if cb not in self._matrix_built_callbacks:
            self._matrix_built_callbacks.append(cb)

    def _fire_matrix_built(self) -> None:
        for cb in list(self._matrix_built_callbacks):
            try:
                cb()
            except Exception as exc:
                print(f"[service] matrix_built callback error: {exc}")

    def register_cycles_found_callback(self, cb: Callable) -> None:
        if cb not in self._cycles_found_callbacks:
            self._cycles_found_callbacks.append(cb)

    def _fire_cycles_found(self, cycles: list) -> None:
        for cb in list(self._cycles_found_callbacks):
            try:
                cb(cycles)
            except Exception as exc:
                print(f"[service] cycles_found callback error: {exc}")

    # -- global UI theme (one switch drives every loaded model page) -----

    def register_theme_changed_callback(self, cb: Callable) -> None:
        if cb not in self._theme_changed_callbacks:
            self._theme_changed_callbacks.append(cb)

    def set_theme(self, theme: str) -> None:
        self.theme = theme
        for cb in list(self._theme_changed_callbacks):
            try:
                cb(theme)
            except Exception as exc:
                print(f"[service] theme_changed callback error: {exc}")

    # ──────────────────────────────────────────────────────────────────
    # Auto-plateau equilibrium scan
    # ──────────────────────────────────────────────────────────────────

    def run_equilibrium_auto(self, initial_by_label: Dict[str, float],
                             flux_list: List[float],
                             energy_eV: float = MACS_KT_KEV * 1e3,
                             method: str = "cram16",
                             mt_filter: Optional[List[int]] = None,
                             decay_mode_filter: Optional[List[str]] = None,
                             hl_preset: str = "all",
                             open_boundary: bool = False,
                             open_boundary_mts: Optional[List[int]] = None,
                             progress: ProgressCb = None,
                             should_cancel: Optional[Callable[[], bool]] = None) -> EqResult:
        """Equilibrium scan using solve_equilibrium for auto-plateau detection per flux point."""
        if not self.initialized:
            raise RuntimeError("Service not initialised. Call init() first.")
        mt_filter = mt_filter or [102]
        if not self.isotopes:
            raise RuntimeError("No burnup universe. Build the matrix first.")

        names = [iso.name for iso in self.isotopes]
        index = {n: i for i, n in enumerate(names)}
        N0 = np.zeros(len(names))
        applied: Dict[str, float] = {}
        for label, conc in initial_by_label.items():
            nm = self.normalize_label(label)
            if nm in index:
                N0[index[nm]] += float(conc)
                applied[nm] = applied.get(nm, 0.0) + float(conc)

        total = max(1, len(flux_list))
        N_eq_list: List[np.ndarray] = []
        Q_eq_list: List[float] = []
        t_reached_list: List[float] = []

        for i, flux in enumerate(flux_list):
            if should_cancel and should_cancel():
                raise RuntimeError("Equilibrium scan cancelled.")
            if progress:
                progress(0.1 + 0.88 * i / total,
                         f"Auto-plateau {i + 1}/{total}  φ={flux:.2e} n/cm²/s…")
            mtx = self.build_matrix(flux=flux, energy_eV=energy_eV,
                                    mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
                                    open_boundary=open_boundary,
                                    open_boundary_mts=open_boundary_mts,
                                    hl_preset=hl_preset)
            try:
                N_i, t_reached = mtx.solve_equilibrium(N0, method=method)
            except Exception:
                N_i = N0.copy()
                t_reached = 0.0
            Q_i = float(mtx.heat_release(N_i))
            N_eq_list.append(N_i)
            Q_eq_list.append(Q_i)
            t_reached_list.append(float(t_reached))

        self.last_eq = EqResult(
            names=names, flux_list=list(flux_list), time_list=t_reached_list,
            N_eq=N_eq_list, Q_eq=Q_eq_list,
            isotopes=list(self.isotopes), initial=applied,
            energy_eV=float(energy_eV), method=method,
        )
        return self.last_eq

    @staticmethod
    def _geo_times(t_end_s: float, n: int = 400, t_start: float = 1e-6) -> List[float]:
        if t_end_s <= 0.0 or n <= 1:
            return [0.0, max(0.0, t_end_s)]
        t_min = max(t_start, 1e-9)
        return [0.0] + np.geomspace(t_min, t_end_s, n - 1).tolist()

    # ──────────────────────────────────────────────────────────────────
    # Label helpers
    # ──────────────────────────────────────────────────────────────────
    @staticmethod
    def parse_label(label: str) -> Optional[tuple]:
        """'Pb-206' / 'U-235m1' → (Z, X, A, meta) or None."""
        m = _LABEL_RE.match(label or "")
        if not m:
            return None
        X = m.group("X").capitalize()
        try:
            Z = ChemTable.Z_of(X)
        except ValueError:
            return None
        meta = int(m.group("meta")) if m.group("meta") else None
        return (Z, X, int(m.group("A")), meta)

    @staticmethod
    def expand_label_token(token: str) -> List[str]:
        """
        Expand one manual-list token. Plain isotope labels ('Pb-206') pass
        through unchanged. An inclusive mass-number range for the same
        element ('Tl206..Tl216') expands into every ground-state label in
        that range. Mismatched elements or unparsable tokens pass through
        unchanged (parse_label will drop them later, same as before).
        """
        m = _RANGE_RE.match(token or "")
        if not m:
            return [token]
        x1, x2 = m.group("X1").capitalize(), m.group("X2").capitalize()
        if x1 != x2:
            return [token]
        a1, a2 = int(m.group("A1")), int(m.group("A2"))
        lo, hi = min(a1, a2), max(a1, a2)
        return [f"{x1}-{a}" for a in range(lo, hi + 1)]

    @classmethod
    def normalize_label(cls, label: str) -> str:
        """Canonical isotope name matching Isotope.name ('Pb-206', 'Pb-209m1')."""
        spec = cls.parse_label(label)
        if spec is None:
            return str(label).strip()
        Z, X, A, meta = spec
        return f"{X}-{A}" + (f"m{meta}" if meta else "")

    def find_isotope(self, label: str) -> Optional[Isotope]:
        nm = self.normalize_label(label)
        for iso in self.isotopes:
            if iso.name == nm:
                return iso
        return None

    @staticmethod
    def isotope_label(iso) -> str:
        return getattr(iso, "name", str(iso))


# ──────────────────────────────────────────────────────────────────────
# Module-level singleton
# ──────────────────────────────────────────────────────────────────────
_SERVICE: Optional[NuclearService] = None


def get_service() -> NuclearService:
    """The shared service instance (created lazily, initialised by the launcher)."""
    global _SERVICE
    if _SERVICE is None:
        _SERVICE = NuclearService()
    return _SERVICE


def init_service(project_root, progress: ProgressCb = None, **kwargs) -> NuclearService:
    """Create (if needed) and initialise the shared service."""
    svc = get_service()
    svc.init(project_root, progress=progress, **kwargs)
    return svc


def is_initialized() -> bool:
    return _SERVICE is not None and _SERVICE.initialized
