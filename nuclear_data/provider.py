"""
nuclear_data/provider.py
========================
Facade for nuclear data access.

Default decay data source priority:
    user > endf6

Users may specify a custom order via the dec_priority parameter:

    nd = NuclearDataProvider(
        ...
        dec_priority=["user", "endf6"]
    )

Available decay source names:
    "user"   — user database (user_db/isotopes.json)
    "endf6"  — ENDF-6 format files (dec_dir, MF=8/MT=457)

Default MACS source priority:
    user > rawmacs > endfb71 > eaf2010 > talys > endf_computed

Users may specify a custom order via the macs_priority parameter:

    nd = NuclearDataProvider(
        ...
        macs_priority=["user", "endfb71", "rawmacs", "talys", "eaf2010", "endf_computed"]
    )

Available MACS source names:
    "user"          — user database (user_db)
    "rawmacs"       — rawmacs.txt (MacsReader)
    "endfb71"       — multi-kT ENDF/B-VII.1 tables (MultiMacsReader)
    "eaf2010"       — multi-kT EAF-2010 tables    (MultiMacsReader)
    "talys"         — TALYS 1.96, kT=30 keV only, MT=102 only
    "endf_computed" — sigma(E) from ENDF file integrated over Maxwell-Boltzmann

sigma(E) priority: user > ENDF > TALYS.
"""

from __future__ import annotations
from typing import Dict, List, Optional

from .models import DecayData, MacsPoint, NeutronReactionData
from .macs_reader import MacsReader, MultiMacsReader
from .user_reader import UserDataReader
from .talys_reader import NreacTalys

try:
    from .endf_reader import EndfReader, _HAS_ENDF
except Exception:
    EndfReader = None   # type: ignore
    _HAS_ENDF = False

# Default decay source priority
DEC_DEFAULT_PRIORITY: List[str] = ["user", "endf6"]

# All known MACS sources in recommended order
MACS_DEFAULT_PRIORITY: List[str] = [
    "user", "rawmacs", "endfb71", "eaf2010", "talys", "endf_computed"
]


class NuclearDataProvider:

    def __init__(
        self,
        dec_dir=None,
        neu_dir=None,
        macs_file=None,
        macs_lib=None,
        user_dir=None,
        talys_dir=None,
        endfb71_macs_dir=None,
        eaf2010_macs_dir=None,
        macs_priority: Optional[List[str]] = None,
        dec_priority: Optional[List[str]] = None,
    ):
        """
        Parameters
        ----------
        dec_dir          : folder with ENDF decay files (dec-*.endf)
        neu_dir          : folder with ENDF neutron cross-section files (n-*.endf)
        macs_file        : path to rawmacs.txt
        macs_lib         : library name in rawmacs.txt (e.g. "Endfb7")
        user_dir         : user database folder
        talys_dir        : TALYS root directory
        endfb71_macs_dir : ENDF/B-VII.1 multi-kT MACS directory
        eaf2010_macs_dir : EAF-2010 multi-kT MACS directory
        macs_priority    : MACS source list in descending priority order.
                           Default: MACS_DEFAULT_PRIORITY
        dec_priority     : decay source list in descending priority order.
                           Default: DEC_DEFAULT_PRIORITY
        """
        self._dec_dir   = dec_dir
        self._neu_dir   = neu_dir
        self._user      = UserDataReader(user_dir) if user_dir else None
        self._talys_dir = talys_dir

        # rawmacs
        self._macs: Optional[MacsReader] = (
            MacsReader(macs_file, macs_lib) if (macs_file and macs_lib) else None
        )
        # multi-temperature MACS
        self._macs_endfb71: Optional[MultiMacsReader] = (
            MultiMacsReader(endfb71_macs_dir) if endfb71_macs_dir else None
        )
        self._macs_eaf2010: Optional[MultiMacsReader] = (
            MultiMacsReader(eaf2010_macs_dir) if eaf2010_macs_dir else None
        )

        self._macs_priority: List[str] = (
            list(macs_priority) if macs_priority is not None else list(MACS_DEFAULT_PRIORITY)
        )
        self._dec_priority: List[str] = (
            list(dec_priority) if dec_priority is not None else list(DEC_DEFAULT_PRIORITY)
        )

        # unified ENDF reader cache (one object per isotope covers both decay and neutron)
        self._endf_cache: Dict[tuple, object] = {}
        self._talys_cache: Dict[tuple, object] = {}

    # ── keys ──

    @staticmethod
    def _key(Z, X, A, meta):
        return (int(Z), str(X).capitalize(), int(A), meta)

    @staticmethod
    def _iso_str(Z, X, A, meta) -> str:
        sym = str(X).capitalize()
        if meta is None or meta == 0 or meta == "":
            m = ""
        elif isinstance(meta, int):
            m = f"m{meta}"
        else:
            s = str(meta).strip().lower()
            m = s if s.startswith("m") else f"m{s}"
        return f"{Z}-{sym}-{A}{m}"

    # ── readers (lazy, cached) ──

    def _endf_reader(self, Z, X, A, meta) -> Optional["EndfReader"]:
        """
        Returns an EndfReader for the given isotope.
        Loads both files (decay + neutron) if present.
        Result is cached — one object per isotope.
        """
        if not (_HAS_ENDF and (self._dec_dir or self._neu_dir)):
            return None
        k = self._key(Z, X, A, meta)
        if k not in self._endf_cache:
            try:
                self._endf_cache[k] = EndfReader(
                    Z, X, A, meta,
                    dec_dir=self._dec_dir,
                    neu_dir=self._neu_dir,
                )
            except (FileNotFoundError, ImportError):
                self._endf_cache[k] = None
        return self._endf_cache[k]

    def _talys_reader(self, Z, X, A, meta):
        if not self._talys_dir:
            return None
        k = self._key(Z, X, A, meta)
        if k not in self._talys_cache:
            self._talys_cache[k] = NreacTalys(self._talys_dir, Z, X, A, meta)
        return self._talys_cache[k]

    # ── decay ──

    def get_decay(self, Z, X, A, meta=None) -> Optional[DecayData]:
        """Source query order is determined by self._dec_priority."""
        for source in self._dec_priority:
            data = self._decay_from_source(source, Z, X, A, meta)
            if data is not None:
                return data
        return None

    def _decay_from_source(self, source: str, Z, X, A, meta) -> Optional[DecayData]:
        if source == "user":
            if self._user and self._user.has_decay(Z, X, A, meta):
                data = self._user.get_decay_data(Z, X, A, meta)
                if data:
                    import math
                    ov = self._user.get_override(Z, X, A, "half_life_s", meta)
                    if ov is not None:
                        t = float(ov)
                        data.half_life_s = t
                        data.decay_constant = math.log(2) / t if t > 0 else 0.0
                    return data
        elif source == "endf6":
            r = self._endf_reader(Z, X, A, meta)
            if r and r.has_decay:
                return r.get_decay_data()
        return None

    def set_dec_priority(self, priority: List[str]) -> None:
        self._dec_priority = list(priority)

    def get_dec_priority(self) -> List[str]:
        return list(self._dec_priority)

    # ── MACS (priority controlled by user) ──

    def get_macs(self, Z, X, A, mt, kT_keV, meta=None) -> Optional[MacsPoint]:
        """
        Returns MACS for (Z, X, A, MT) at the given kT [keV].
        Source query order is determined by self._macs_priority.
        """
        iso_key = self._iso_str(Z, X, A, meta)
        for source in self._macs_priority:
            pt = self._macs_from_source(source, Z, X, A, mt, kT_keV, meta, iso_key)
            if pt is not None:
                return pt
        return None

    def get_macs_with_source(self, Z, X, A, mt, kT_keV, meta=None):
        """Like get_macs, but also returns the source name that supplied the
        value (e.g. "eaf2010", "talys"). Returns (None, None) when nothing
        is found. Lets callers record data provenance for the UI."""
        iso_key = self._iso_str(Z, X, A, meta)
        for source in self._macs_priority:
            pt = self._macs_from_source(source, Z, X, A, mt, kT_keV, meta, iso_key)
            if pt is not None:
                return pt, source
        return None, None

    def _macs_from_source(
        self, source: str,
        Z, X, A, mt, kT_keV, meta, iso_key: str
    ) -> Optional[MacsPoint]:
        if source == "user":
            if self._user and self._user.has_reaction(Z, X, A, mt, meta):
                rxn = self._user.get_reaction_data(Z, X, A, mt, meta)
                if rxn and rxn.macs:
                    return self._interp_macs(rxn.macs, kT_keV)

        elif source == "rawmacs":
            if self._macs:
                return self._macs.get_macs_at(iso_key, mt, kT_keV)

        elif source == "endfb71":
            if self._macs_endfb71:
                return self._macs_endfb71.get_macs_at(iso_key, mt, kT_keV)

        elif source == "eaf2010":
            if self._macs_eaf2010:
                return self._macs_eaf2010.get_macs_at(iso_key, mt, kT_keV)

        elif source == "talys":
            t = self._talys_reader(Z, X, A, meta)
            if t:
                return t.get_macs_30keV(mt)

        elif source == "endf_computed":
            r = self._endf_reader(Z, X, A, meta)
            if r and r.has_neutron:
                return r.compute_macs(mt, kT_keV, "keV")

        return None

    @staticmethod
    def _interp_macs(pts, kT: float) -> Optional[MacsPoint]:
        pts = sorted(pts, key=lambda p: p.kT_keV)
        if kT <= pts[0].kT_keV:
            return pts[0]
        if kT >= pts[-1].kT_keV:
            return pts[-1]
        for i in range(len(pts) - 1):
            a, b = pts[i], pts[i + 1]
            if a.kT_keV <= kT <= b.kT_keV:
                w = (kT - a.kT_keV) / (b.kT_keV - a.kT_keV)
                return MacsPoint(kT, a.sigma_mb * (1 - w) + b.sigma_mb * w)
        return None

    # ── sigma(E) ──

    def get_sigma(self, Z, X, A, mt, E, e_unit="eV", meta=None) -> Optional[float]:
        """Point-wise σ(E). Priority: user > ENDF > TALYS."""
        v, _ = self.get_sigma_with_source(Z, X, A, mt, E, e_unit, meta)
        return v

    def get_sigma_with_source(self, Z, X, A, mt, E, e_unit="eV", meta=None):
        """get_sigma + the source name ("user"/"endf"/"talys").

        A zero/negative value from a higher-priority source means that
        evaluation carries no real cross section at this energy (e.g. the
        ENDF Po-210 (n,γ) evaluation is empty below ~285 keV), so fall through
        to the next source instead of reporting 0. Returns (None, None) when no
        source has data.
        """
        if self._user and self._user.has_reaction(Z, X, A, mt, meta):
            rxn = self._user.get_reaction_data(Z, X, A, mt, meta)
            if rxn and rxn.spectrum:
                conv  = {"ev": 1.0, "kev": 1e3, "mev": 1e6}
                E_eV  = float(E) * conv[e_unit.lower()]
                pairs = [(p.E_eV, p.sigma_barn) for p in rxn.spectrum]
                v     = self._interp_pairs(pairs, E_eV)
                if v is not None and v > 0:
                    return v, "user"
        r = self._endf_reader(Z, X, A, meta)
        if r and r.has_neutron:
            v = r.get_sigma_at(mt, E, e_unit)
            if v is not None and v > 0:
                return v, "endf"
        t = self._talys_reader(Z, X, A, meta)
        if t:
            v = t.get_sigma_at(mt, E, e_unit)
            if v is not None and v > 0:
                return v, "talys"
        return None, None

    @staticmethod
    def _interp_pairs(pairs, x: float) -> Optional[float]:
        pairs = sorted(pairs)
        if not pairs:
            return None
        if x <= pairs[0][0]:
            return pairs[0][1]
        if x >= pairs[-1][0]:
            return pairs[-1][1]
        for i in range(len(pairs) - 1):
            x0, y0 = pairs[i]; x1, y1 = pairs[i + 1]
            if x0 <= x <= x1:
                return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        return None

    # ── isomeric branching ──

    def get_isomeric_branches(self, Z, X, A, mt, E_eV, meta=None):
        """Q(ij) coefficients. ENDF only; falls back to [LFS=0, yield=1.0]."""
        from .models import IsomericBranch
        r = self._endf_reader(Z, X, A, meta)
        if r and r.has_neutron:
            return r.get_isomeric_branches(mt, E_eV)
        return [IsomericBranch(0, 1.0)]

    # ── full reaction data ──

    def get_reaction_data(self, Z, X, A, mt, meta=None) -> Optional[NeutronReactionData]:
        """Priority: user > ENDF > TALYS."""
        if self._user and self._user.has_reaction(Z, X, A, mt, meta):
            d = self._user.get_reaction_data(Z, X, A, mt, meta)
            if d:
                return d
        r = self._endf_reader(Z, X, A, meta)
        if r and r.has_neutron:
            d = r.get_reaction_data(mt)
            if d:
                return d
        t = self._talys_reader(Z, X, A, meta)
        return t.get_reaction_data(mt) if t else None

    # ── utilities ──

    def get_override(self, Z, X, A, field, meta=None):
        if self._user:
            return self._user.get_override(Z, X, A, field, meta)
        return None

    def reload_user_db(self):
        if self._user:
            self._user.reload()

    def set_macs_priority(self, priority: List[str]) -> None:
        self._macs_priority = list(priority)

    def get_macs_priority(self) -> List[str]:
        return list(self._macs_priority)

    def __repr__(self):
        active = []
        if self._user:                          active.append("user")
        if _HAS_ENDF and self._dec_dir:         active.append("endf-dec")
        if _HAS_ENDF and self._neu_dir:         active.append("endf-neu")
        if self._macs:                          active.append("rawmacs")
        if self._macs_endfb71:                  active.append("endfb71")
        if self._macs_eaf2010:                  active.append("eaf2010")
        if self._talys_dir:                     active.append("talys")
        return (f"NuclearDataProvider(sources=[{', '.join(active)}], "
                f"dec_priority={self._dec_priority}, "
                f"macs_priority={self._macs_priority})")
