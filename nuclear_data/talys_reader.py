"""
nuclear_data/talys_reader.py
============================
Reader for recommended TALYS 1.96 cross sections (PSI, Switzerland).

TALYS directory structure:
    TALYS/
      recommended-ng-cross-sections/   MT=102  sigma(E), E[MeV], sigma[mb]
      recommended-na-cross-sections/   MT=107  sigma(E), E[MeV], sigma[mb]
      recommended-np-cross-sections/   MT=103  sigma(E), E[MeV], sigma[mb]
      talys_macs_30_keV/               MT=102  MACS at kT=30 keV, sigma[mb]

Filename: {Sym}{A:03d}.txt  (e.g. Xe135.txt, Ag088.txt).
Data line format: E_MeV  sigma_mb  uncertainty_mb

Instance class: one object = one isotope.
Constructor always succeeds; methods return empty data / None
when the required file is absent.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .models import CrossSectionPoint, MacsPoint, NeutronReactionData


class NreacTalys:
    """TALYS cross-section reader for one isotope."""

    # MT -> (subdirectory, channel label)
    _MT_DIR: Dict[int, Tuple[str, str]] = {
        102: ("recommended-ng-cross-sections", "(n,γ)"),
        103: ("recommended-np-cross-sections", "(n,p)"),
        107: ("recommended-na-cross-sections", "(n,α)"),
    }
    _MACS_SUBDIR = "talys_macs_30_keV"
    _MACS_KT_KEV = 30.0

    def __init__(self, talys_root, Z, X, A, meta=None):
        self._root = Path(talys_root).resolve()
        self._Z = int(Z)
        self._X = str(X).capitalize()
        self._A = int(A)
        self._fname = f"{self._X}{self._A:03d}.txt"

    # ── internal ──

    @staticmethod
    def _interp(pairs: List[Tuple[float, float]], x: float) -> Optional[float]:
        if not pairs:
            return None
        pairs = sorted(pairs, key=lambda t: t[0])
        if x <= pairs[0][0]:
            return pairs[0][1]
        if x >= pairs[-1][0]:
            return pairs[-1][1]
        for i in range(len(pairs) - 1):
            x0, y0 = pairs[i]
            x1, y1 = pairs[i + 1]
            if x0 <= x <= x1:
                return y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        return None

    def _read_pairs(self, mt: int) -> Optional[List[Tuple[float, float]]]:
        """
        Reads the sigma(E) file for the given MT.
        Returns [(E_eV, sigma_barn), ...] or None if the file is absent.
        Conversion: E: MeV -> eV (x1e6),  sigma: mb -> barn (x1e-3).
        """
        info = self._MT_DIR.get(mt)
        if info is None:
            return None
        subdir, _ = info
        path = self._root / subdir / self._fname
        if not path.is_file():
            return None
        pairs: List[Tuple[float, float]] = []
        try:
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) < 2:
                    continue
                E_eV = float(parts[0]) * 1e6      # MeV → eV
                sigma_barn = float(parts[1]) * 1e-3  # mb → barn
                pairs.append((E_eV, sigma_barn))
        except Exception:
            return None
        return pairs if pairs else None

    # ── public ──

    def available_mts(self) -> List[int]:
        """List of MT channels for which cross-section files exist on disk."""
        return [
            mt for mt, (subdir, _) in self._MT_DIR.items()
            if (self._root / subdir / self._fname).is_file()
        ]

    def get_cross_section(self, mt: int) -> List[CrossSectionPoint]:
        """Full sigma(E) spectrum as a list of CrossSectionPoint (E in eV, sigma in barn)."""
        pairs = self._read_pairs(mt)
        if not pairs:
            return []
        return [CrossSectionPoint(E_eV, sigma_barn) for E_eV, sigma_barn in pairs]

    def get_sigma_at(self, mt: int, E, e_unit: str = "eV") -> Optional[float]:
        """Interpolated cross section [barn] at energy E."""
        pairs = self._read_pairs(mt)
        if not pairs:
            return None
        conv = {"ev": 1.0, "kev": 1e3, "mev": 1e6}
        E_eV = float(E) * conv[e_unit.lower()]
        return self._interp(pairs, E_eV)

    def get_macs_30keV(self, mt: int = 102) -> Optional[MacsPoint]:
        """
        MACS at kT = 30 keV from the talys_macs_30_keV directory.
        Available only for MT=102 (n,gamma). Returns None for other MT.
        """
        if mt != 102:
            return None
        path = self._root / self._MACS_SUBDIR / self._fname
        if not path.is_file():
            return None
        try:
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) >= 2:
                    return MacsPoint(kT_keV=self._MACS_KT_KEV, sigma_mb=float(parts[1]))
        except Exception:
            pass
        return None

    def get_reaction_data(self, mt: int) -> Optional[NeutronReactionData]:
        """
        Full NeutronReactionData for MT: sigma(E) spectrum + MACS at 30 keV (MT=102 only).
        Returns None if the cross-section file is absent.
        """
        spec = self.get_cross_section(mt)
        if not spec:
            return None
        info = self._MT_DIR.get(mt)
        channel = info[1] if info else f"MT={mt}"
        macs_pt = self.get_macs_30keV(mt) if mt == 102 else None
        return NeutronReactionData(
            mt=mt,
            channel_name=channel,
            Q_MeV=None,
            spectrum=spec,
            macs=[macs_pt] if macs_pt else [],
            source="talys",
        )

    def __repr__(self):
        return f"NreacTalys({self._Z}-{self._X}-{self._A}, root={self._root.name})"
