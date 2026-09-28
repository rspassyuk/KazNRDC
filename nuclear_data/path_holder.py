"""
nuclear_data/path_holder.py
===========================
Central path configuration point.

Expected xsdir structure:
    xsdir/
      endf-6/<LIB>/decay/
      endf-6/<LIB>/neutrons/
      MACS/
        rawmacs.txt
        ENDFB71-EAF2010/
          ENDFB71/   (ENDF/B-VII.1, multi-temperature MACS)
          EAF2010/   (EAF-2010, multi-temperature MACS)
      TALYS/
        recommended-ng-cross-sections/
        recommended-na-cross-sections/
        recommended-np-cross-sections/
        talys_macs_30_keV/
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional




class EndfPathHolder:
    _xsdir: Optional[Path] = None
    _default_dec_dir: Optional[str] = None
    _default_neu_dir: Optional[str] = None
    _macs_file_path: Optional[Path] = None
    _current_macs_library: Optional[str] = None
    _user_db_dir: Optional[Path] = None
    _talys_dir: Optional[Path] = None
    _endfb71_macs_dir: Optional[Path] = None
    _eaf2010_macs_dir: Optional[Path] = None

    MACS_LIB_MAP: Dict[str, str] = {
        "EndfbMacs": "ENDF/B-VII.0", "Endfb7": "ENDF/B-VII.0", "Endfb70": "ENDF/B-VII.0",
        "Endfb6Macs": "ENDF/B-VI.8", "Endfb6": "ENDF/B-VI.8", "Endfb68": "ENDF/B-VI.8",
        "JeffMacs": "JEFF-3.1", "Jeff": "JEFF-3.1", "Jeff31": "JEFF-3.1",
        "JendlMacs": "JENDL-3.3", "Jendl": "JENDL-3.3", "Jendl33": "JENDL-3.3",
        "ENDF/B-VII.0": "ENDF/B-VII.0", "ENDF/B-VI.8": "ENDF/B-VI.8",
        "JEFF-3.1": "JEFF-3.1", "JENDL-3.3": "JENDL-3.3",
    }
    _ALLOWED = {"ENDF/B-VII.0", "ENDF/B-VI.8", "JEFF-3.1", "JENDL-3.3"}

    # ── helpers ──
    @staticmethod
    def _require_dir(p: Path, err: str) -> Path:
        if not p.is_dir():
            raise FileNotFoundError(f"{err}: {p}")
        return p

    @staticmethod
    def _require_file(p: Path, err: str) -> Path:
        if not p.is_file():
            raise FileNotFoundError(f"{err}: {p}")
        return p

    @staticmethod
    def _require_xsdir() -> Path:
        if EndfPathHolder._xsdir is None:
            raise RuntimeError("xsdir not initialized. Call init_xsdir(path) first.")
        return EndfPathHolder._xsdir

    # ── init ──
    @staticmethod
    def init_xsdir(xsdir_path: str) -> Dict[str, str]:
        xsdir = Path(xsdir_path).expanduser().resolve()
        EndfPathHolder._require_dir(xsdir, "xsdir not found")
        endf6 = EndfPathHolder._require_dir(xsdir / "endf-6", "endf-6 not found")
        macs = EndfPathHolder._require_dir(xsdir / "MACS", "MACS not found")
        EndfPathHolder._xsdir = xsdir
        EndfPathHolder._default_dec_dir = None
        EndfPathHolder._default_neu_dir = None
        EndfPathHolder._macs_file_path = None
        EndfPathHolder._current_macs_library = None
        result: Dict[str, str] = {"xsdir": str(xsdir), "endf-6": str(endf6), "MACS": str(macs)}
        # Auto-discover TALYS (optional directory)
        talys = xsdir / "TALYS"
        if talys.is_dir():
            EndfPathHolder._talys_dir = talys
            result["TALYS"] = str(talys)
        else:
            EndfPathHolder._talys_dir = None
        # Auto-discover multi-temperature MACS directories
        endfb71 = xsdir / "MACS" / "ENDFB71-EAF2010" / "ENDFB71"
        if endfb71.is_dir():
            EndfPathHolder._endfb71_macs_dir = endfb71
            result["ENDFB71_MACS"] = str(endfb71)
        else:
            EndfPathHolder._endfb71_macs_dir = None
        eaf2010 = xsdir / "MACS" / "ENDFB71-EAF2010" / "EAF2010"
        if eaf2010.is_dir():
            EndfPathHolder._eaf2010_macs_dir = eaf2010
            result["EAF2010_MACS"] = str(eaf2010)
        else:
            EndfPathHolder._eaf2010_macs_dir = None
        return result

    # ── ENDF dirs ──
    @staticmethod
    def set_DEFAULT_DEC_DIR(library: str) -> str:
        xs = EndfPathHolder._require_xsdir()
        d = EndfPathHolder._require_dir(xs / "endf-6" / library / "decay", "decay folder not found")
        EndfPathHolder._default_dec_dir = str(d)
        return str(d)

    @staticmethod
    def get_DEFAULT_DEC_DIR() -> str:
        if EndfPathHolder._default_dec_dir is None:
            raise RuntimeError("Decay path not set. Call set_DEFAULT_DEC_DIR(library) first.")
        return EndfPathHolder._default_dec_dir

    @staticmethod
    def set_DEFAULT_NEU_DIR(library: str) -> str:
        xs = EndfPathHolder._require_xsdir()
        d = EndfPathHolder._require_dir(xs / "endf-6" / library / "neutrons", "neutrons folder not found")
        EndfPathHolder._default_neu_dir = str(d)
        return str(d)

    @staticmethod
    def get_DEFAULT_NEU_DIR() -> str:
        if EndfPathHolder._default_neu_dir is None:
            raise RuntimeError("Neutrons path not set. Call set_DEFAULT_NEU_DIR(library) first.")
        return EndfPathHolder._default_neu_dir

    # ── MACS ──
    @staticmethod
    def set_DEFAULT_MACS_DIR(library_alias: str, filename: str = "rawmacs.txt") -> str:
        xs = EndfPathHolder._require_xsdir()
        macs_dir = EndfPathHolder._require_dir(xs / "MACS", "MACS folder not found")
        macs_file = (macs_dir / filename).resolve()
        EndfPathHolder._require_file(macs_file, "MACS file not found")
        lib = EndfPathHolder.MACS_LIB_MAP.get(library_alias.strip(), library_alias.strip())
        if lib not in EndfPathHolder._ALLOWED:
            raise ValueError(f"Unknown MACS library '{library_alias}'.")
        EndfPathHolder._macs_file_path = macs_file
        EndfPathHolder._current_macs_library = lib
        return str(macs_file)

    @staticmethod
    def get_macs_file() -> Path:
        if EndfPathHolder._macs_file_path is None:
            raise RuntimeError("MACS file not set.")
        return EndfPathHolder._macs_file_path

    @staticmethod
    def get_macs_library() -> str:
        if EndfPathHolder._current_macs_library is None:
            raise RuntimeError("MACS library not selected.")
        return EndfPathHolder._current_macs_library

    # ── Multi-temperature MACS (ENDFB71 / EAF2010) ──

    @staticmethod
    def set_ENDFB71_MACS_DIR(path: Optional[str] = None) -> str:
        """
        Sets the ENDF/B-VII.1 multi-temperature MACS directory.
        If path=None — searches for xsdir/MACS/ENDFB71-EAF2010/ENDFB71.
        """
        if path is None:
            xs = EndfPathHolder._require_xsdir()
            p = xs / "MACS" / "ENDFB71-EAF2010" / "ENDFB71"
        else:
            p = Path(path).expanduser().resolve()
        EndfPathHolder._require_dir(p, "ENDFB71 MACS directory not found")
        EndfPathHolder._endfb71_macs_dir = p
        return str(p)

    @staticmethod
    def get_ENDFB71_MACS_DIR() -> Path:
        if EndfPathHolder._endfb71_macs_dir is None:
            raise RuntimeError("ENDFB71 MACS not set. Call set_ENDFB71_MACS_DIR() first.")
        return EndfPathHolder._endfb71_macs_dir

    @staticmethod
    def set_EAF2010_MACS_DIR(path: Optional[str] = None) -> str:
        """
        Sets the EAF-2010 multi-temperature MACS directory.
        If path=None — searches for xsdir/MACS/ENDFB71-EAF2010/EAF2010.
        """
        if path is None:
            xs = EndfPathHolder._require_xsdir()
            p = xs / "MACS" / "ENDFB71-EAF2010" / "EAF2010"
        else:
            p = Path(path).expanduser().resolve()
        EndfPathHolder._require_dir(p, "EAF2010 MACS directory not found")
        EndfPathHolder._eaf2010_macs_dir = p
        return str(p)

    @staticmethod
    def get_EAF2010_MACS_DIR() -> Path:
        if EndfPathHolder._eaf2010_macs_dir is None:
            raise RuntimeError("EAF2010 MACS not set. Call set_EAF2010_MACS_DIR() first.")
        return EndfPathHolder._eaf2010_macs_dir

    # ── TALYS ──
    @staticmethod
    def set_TALYS_DIR(path: Optional[str] = None) -> str:
        """
        Sets the TALYS root directory.
        If path=None — searches for xsdir/TALYS (requires prior init_xsdir).
        """
        if path is None:
            xs = EndfPathHolder._require_xsdir()
            p = xs / "TALYS"
        else:
            p = Path(path).expanduser().resolve()
        EndfPathHolder._require_dir(p, "TALYS directory not found")
        EndfPathHolder._talys_dir = p
        return str(p)

    @staticmethod
    def get_TALYS_DIR() -> Path:
        if EndfPathHolder._talys_dir is None:
            raise RuntimeError("TALYS path not set. Call set_TALYS_DIR() first.")
        return EndfPathHolder._talys_dir

    # ── user db ──
    @staticmethod
    def set_USER_DB_DIR(path: str) -> str:
        p = Path(path).expanduser().resolve()
        p.mkdir(parents=True, exist_ok=True)
        (p / "spectra").mkdir(exist_ok=True)
        EndfPathHolder._user_db_dir = p
        return str(p)


    @staticmethod
    def get_USER_DB_DIR() -> Path:
        if EndfPathHolder._user_db_dir is None:
            raise RuntimeError("User database path not set. Call set_USER_DB_DIR(path) first.")
        return EndfPathHolder._user_db_dir
