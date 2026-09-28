"""
nuclear_data/macs_reader.py
===========================
All MACS data sources in one module (data resides in xsdir/MACS).

MacsReader      — reads rawmacs.txt (semicolon-delimited CSV)
MultiMacsReader — reads ENDFB71 / EAF2010 directories
                  (files {LIB}-Max-MT{mt}-{kT_eV}.txt, MACS in barn)

Both classes use a pkl cache for fast repeated runs.
"""

from __future__ import annotations
import os
import pickle
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .models import MacsPoint

_Cache = Dict[str, Dict[str, Dict[int, List[Tuple[float, float, float]]]]]


class MacsReader:
    _cache: Dict[Path, _Cache] = {}

    MACS_LIB_MAP = {
        "EndfbMacs": "ENDF/B-VII.0", "Endfb7": "ENDF/B-VII.0", "Endfb70": "ENDF/B-VII.0",
        "Endfb6Macs": "ENDF/B-VI.8", "Endfb6": "ENDF/B-VI.8", "Endfb68": "ENDF/B-VI.8",
        "JeffMacs": "JEFF-3.1", "Jeff": "JEFF-3.1", "Jeff31": "JEFF-3.1",
        "JendlMacs": "JENDL-3.3", "Jendl": "JENDL-3.3", "Jendl33": "JENDL-3.3",
        "ENDF/B-VII.0": "ENDF/B-VII.0", "ENDF/B-VI.8": "ENDF/B-VI.8",
        "JEFF-3.1": "JEFF-3.1", "JENDL-3.3": "JENDL-3.3",
    }
    _ALLOWED = {"ENDF/B-VII.0", "ENDF/B-VI.8", "JEFF-3.1", "JENDL-3.3"}

    def __init__(self, macs_file, library_alias):
        self._file = Path(macs_file).resolve()
        if not self._file.is_file():
            raise FileNotFoundError(f"MACS file not found: {self._file}")
        lib = self.MACS_LIB_MAP.get(library_alias.strip(), library_alias.strip())
        if lib not in self._ALLOWED:
            raise ValueError(f"Unknown MACS library '{library_alias}'.")
        self._library = lib
        self._ensure_loaded()

    def _signature(self):
        st = self._file.stat()
        return (st.st_size, int(st.st_mtime))

    def _ensure_loaded(self):
        if self._file in MacsReader._cache:
            return
        pkl = self._file.with_suffix(self._file.suffix + ".pkl")
        try:
            if pkl.is_file():
                pl = pickle.loads(pkl.read_bytes())
                if isinstance(pl, dict) and tuple(pl.get("sig", ())) == self._signature():
                    MacsReader._cache[self._file] = pl["data"]
                    return
        except Exception:
            pass
        data = self._parse()
        MacsReader._cache[self._file] = data
        try:
            tmp = pkl.with_suffix(".tmp")
            tmp.write_bytes(pickle.dumps({"sig": self._signature(), "data": data},
                                         protocol=pickle.HIGHEST_PROTOCOL))
            os.replace(tmp, pkl)
        except Exception:
            pass

    def _parse(self) -> _Cache:
        res: _Cache = {}
        with open(self._file, encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#") or s.lower().startswith("library"):
                    continue
                parts = [p.strip() for p in s.split(";")]
                if len(parts) < 6:
                    continue
                try:
                    lib, iso = parts[0], self._norm_key(parts[1])
                    mt = int(parts[2]); kT = self._num(parts[3])
                    mb = self._num(parts[4]); rate = self._num(parts[5])
                except Exception:
                    continue
                res.setdefault(lib, {}).setdefault(iso, {}).setdefault(mt, []).append((kT, mb, rate))
        for lib in res:
            for iso in res[lib]:
                for mt in res[lib][iso]:
                    res[lib][iso][mt].sort(key=lambda t: t[0])
        return res

    @staticmethod
    def _num(s):
        t = str(s).strip().replace(",", ".")
        if not t:
            raise ValueError("empty")
        if "e" not in t.lower():
            for i in range(len(t)-1, 0, -1):
                if t[i] in "+-":
                    try:
                        return float(t[:i] + "e" + t[i:])
                    except Exception:
                        pass
        return float(t)

    @staticmethod
    def _norm_key(s):
        m = re.match(r"^\s*(\d+)\s*-\s*([A-Za-z]{1,2})\s*-\s*(\d+)\s*([Mm]\d*)?\s*$", s.strip())
        if not m:
            return s.strip()
        meta = (m.group(4) or "").strip()
        if meta.upper().startswith("M"):
            meta = "M"
        return f"{m.group(1)}-{m.group(2).capitalize()}-{m.group(3)}{meta}"

    @staticmethod
    def _adapt(iso_key):
        return re.sub(r"m\d+$", "M", iso_key)

    # ── public ──
    def get_all_points(self, iso_key, mt) -> List[MacsPoint]:
        key = self._adapt(iso_key)
        series = MacsReader._cache.get(self._file, {}).get(self._library, {}).get(key, {}).get(mt, [])
        return [MacsPoint(kT, mb, r) for kT, mb, r in series]

    def get_macs_at(self, iso_key, mt, kT_keV) -> Optional[MacsPoint]:
        pts = self.get_all_points(iso_key, mt)
        if not pts:
            return None
        kTs = [p.kT_keV for p in pts]
        if kT_keV <= kTs[0]:
            return pts[0]
        if kT_keV >= kTs[-1]:
            return pts[-1]
        for i in range(len(kTs)-1):
            if kTs[i] <= kT_keV <= kTs[i+1]:
                w = (kT_keV - kTs[i]) / (kTs[i+1] - kTs[i])
                mb = pts[i].sigma_mb*(1-w) + pts[i+1].sigma_mb*w
                r0, r1 = pts[i].rate_cm3_mol_s, pts[i+1].rate_cm3_mol_s
                rate = (r0*(1-w)+r1*w) if (r0 is not None and r1 is not None) else None
                return MacsPoint(kT_keV, mb, rate)
        return None

    def __repr__(self):
        return f"MacsReader({self._file.name}, {self._library})"


# ═══════════════════════════════════════════════════════════════════════════
#  MultiMacsReader  — multi-temperature tables (ENDFB71 / EAF2010)
# ═══════════════════════════════════════════════════════════════════════════

_CacheData = Dict[str, Dict[int, List[Tuple[float, float, Optional[float]]]]]
# {iso_key: {mt: [(kT_keV, macs_mb, rate_cm3_mol_s), ...]}}


class MultiMacsReader:
    """
    MACS reader for multi-temperature tabulated files.
    One instance = one directory (ENDFB71 or EAF2010).

    File format: {LIB}-Max-MT{mt}-{kT_eV}.txt
    Columns (comma CSV): #, library, MT, MAT, isotope, kT_eV, MACS_barn, [unc], rate, [rate_unc]
    MACS stored in barn in the files; the API returns millibarn.

    Reference: Pritychenko & Mughabghab, NDS 113, 3120 (2012).
    """

    _file_cache: Dict[Path, _CacheData] = {}

    def __init__(self, directory):
        self._dir = Path(directory).resolve()
        if not self._dir.is_dir():
            raise FileNotFoundError(
                f"MultiMacsReader: directory not found: {self._dir}")
        self._data = self._ensure_loaded(self._dir)

    # ── cache ──

    @classmethod
    def _ensure_loaded(cls, directory: Path) -> _CacheData:
        if directory in cls._file_cache:
            return cls._file_cache[directory]
        pkl = directory / "_multi_macs_cache.pkl"
        try:
            sig = cls._dir_signature(directory)
            if pkl.is_file():
                pl = pickle.loads(pkl.read_bytes())
                if isinstance(pl, dict) and pl.get("sig") == sig:
                    cls._file_cache[directory] = pl["data"]
                    return cls._file_cache[directory]
        except Exception:
            pass
        data = cls._load_directory(directory)
        cls._file_cache[directory] = data
        try:
            tmp = pkl.with_suffix(".tmp")
            tmp.write_bytes(pickle.dumps(
                {"sig": cls._dir_signature(directory), "data": data},
                protocol=pickle.HIGHEST_PROTOCOL))
            os.replace(tmp, pkl)
        except Exception:
            pass
        return data

    @staticmethod
    def _dir_signature(directory: Path):
        files = list(directory.glob("*-Max-MT*-*.txt"))
        return (len(files), sum(f.stat().st_size for f in files if f.is_file()))

    # ── loading ──

    @classmethod
    def _load_directory(cls, directory: Path) -> _CacheData:
        data: _CacheData = {}
        for fpath in sorted(directory.glob("*-Max-MT*-*.txt")):
            cls._parse_file(fpath, data)
        for iso in data:
            for mt in data[iso]:
                data[iso][mt].sort(key=lambda t: t[0])
        return data

    @classmethod
    def _parse_file(cls, path: Path, data: _CacheData) -> None:
        m = re.search(r"-MT(\d+)-(.+)\.txt$", path.name, re.IGNORECASE)
        if not m:
            return
        try:
            mt = int(m.group(1))
            kT_keV = float(m.group(2)) / 1000.0
        except (ValueError, TypeError):
            return
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            return
        for line in text.splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(",")]
            if len(parts) < 7:
                continue
            try:
                mat = int(parts[3])
                macs_mb = float(parts[6]) * 1000.0  # barn → mb
            except (ValueError, IndexError):
                continue
            key = cls._parse_iso_name(parts[4], mat)
            if key is None:
                continue
            rate: Optional[float] = None
            if len(parts) > 8 and parts[8]:
                try:
                    rate = float(parts[8])
                except ValueError:
                    pass
            data.setdefault(key, {}).setdefault(mt, []).append((kT_keV, macs_mb, rate))

    @staticmethod
    def _parse_iso_name(raw: str, mat: int) -> Optional[str]:
        raw = raw.strip()
        m = re.match(r"^(\d+)\s*-\s*([A-Za-z]{1,2})\s*-\s*(\d+)\s*$", raw)
        if m:
            Z, sym, A = int(m.group(1)), m.group(2).capitalize(), int(m.group(3))
            return None if A == 0 else f"{Z}-{sym}-{A}"
        m = re.match(r"^([A-Za-z]{1,2})-(\d+)\s*$", raw)
        if m:
            sym = m.group(1).capitalize()
            A   = int(m.group(2))
            Z   = mat // 100
            return None if (A == 0 or Z <= 0) else f"{Z}-{sym}-{A}"
        return None

    @staticmethod
    def _norm_key(iso_key: str) -> str:
        m = re.match(r"^(\d+)-([A-Za-z]{1,2})-(\d+)([Mm]\d*)?$",
                     str(iso_key).strip())
        if m:
            meta = (m.group(4) or "").strip()
            return f"{int(m.group(1))}-{m.group(2).capitalize()}-{int(m.group(3))}{meta}"
        return iso_key.strip()

    # ── interpolation ──

    @staticmethod
    def _interp(pts, kT: float) -> Optional[MacsPoint]:
        if not pts:
            return None
        kTs = [p[0] for p in pts]
        if kT <= kTs[0]:
            return MacsPoint(kT, pts[0][1], pts[0][2])
        if kT >= kTs[-1]:
            return MacsPoint(kT, pts[-1][1], pts[-1][2])
        for i in range(len(kTs) - 1):
            if kTs[i] <= kT <= kTs[i + 1]:
                w  = (kT - kTs[i]) / (kTs[i + 1] - kTs[i])
                mb = pts[i][1] * (1 - w) + pts[i + 1][1] * w
                r0, r1 = pts[i][2], pts[i + 1][2]
                rate = (r0 * (1 - w) + r1 * w
                        if r0 is not None and r1 is not None else None)
                return MacsPoint(kT, mb, rate)
        return None

    # ── public API ──

    def get_all_points(self, iso_key: str, mt: int = 102) -> List[MacsPoint]:
        """All MACS(kT) points [mb] for an isotope and MT."""
        key    = self._norm_key(iso_key)
        series = self._data.get(key, {}).get(mt, [])
        return [MacsPoint(kT, mb, r) for kT, mb, r in series]

    def get_macs_at(self, iso_key: str, mt: int,
                    kT_keV: float) -> Optional[MacsPoint]:
        """Interpolated MACS [mb] at kT [keV]."""
        key = self._norm_key(iso_key)
        pts = self._data.get(key, {}).get(mt, [])
        return self._interp(pts, kT_keV)

    def available_isotopes(self) -> List[str]:
        return list(self._data.keys())

    def __repr__(self):
        n_iso = len(self._data)
        return f"MultiMacsReader({self._dir.name}, {n_iso} isotopes)"
