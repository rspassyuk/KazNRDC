"""
nuclear_data/user_reader.py
===========================
User nuclear data database.

    user_db/
      isotopes.json     <- decay data and MACS points
      override.json     <- point-wise overrides
      spectra/<key>_mt<MT>.csv   <- sigma(E) spectra
"""

from __future__ import annotations
import csv
import json
import math
from pathlib import Path
from typing import Any, List, Optional

from .models import CrossSectionPoint, DecayChannel, DecayData, MacsPoint, NeutronReactionData


class UserDataReader:

    def __init__(self, db_dir):
        self._dir = Path(db_dir).expanduser().resolve()
        if not self._dir.is_dir():
            raise FileNotFoundError(f"User DB not found: {self._dir}")
        self._iso: dict = {}
        self._ovr: dict = {}
        self._reload()

    def _reload(self):
        ip = self._dir / "isotopes.json"
        op = self._dir / "override.json"
        self._iso = json.loads(ip.read_text(encoding="utf-8")) if ip.is_file() else {}
        self._ovr = json.loads(op.read_text(encoding="utf-8")) if op.is_file() else {}
        names = sorted(self._iso.keys())
        print(f"[user_db] Loaded from: {ip}")
        print(f"[user_db] {len(names)} isotope(s) found: {', '.join(names) if names else '(none)'}")

    def reload(self):
        self._reload()

    @staticmethod
    def _key(Z, X, A, meta=None):
        sym = str(X).strip().capitalize()
        if meta is None or meta == 0 or meta == "":
            m = ""
        elif isinstance(meta, int):
            m = f"m{meta}"
        else:
            s = str(meta).strip().lower()
            m = s if s.startswith("m") else f"m{s}"
        return f"{Z}-{sym}-{A}{m}"

    # ── decay ──
    def get_decay_data(self, Z, X, A, meta=None) -> Optional[DecayData]:
        dec = self._iso.get(self._key(Z, X, A, meta), {}).get("decay")
        if not dec:
            return None
        t = dec.get("half_life_s")
        lam = (math.log(2)/t) if t else 0.0
        ch = [DecayChannel(c.get("mode", "unknown"), float(c.get("branch", 1.0)),
                           float(c.get("Q_MeV", 0.0)), int(c.get("dz", 0)),
                           int(c.get("da", 0)), int(c.get("rfs", 0)))
              for c in dec.get("channels", [])]
        return DecayData(t, lam, ch, "user")

    def has_decay(self, Z, X, A, meta=None):
        return bool(self._iso.get(self._key(Z, X, A, meta), {}).get("decay"))

    # ── reactions ──
    def get_reaction_data(self, Z, X, A, mt, meta=None) -> Optional[NeutronReactionData]:
        key = self._key(Z, X, A, meta)
        raw = self._iso.get(key, {}).get("neutrons", {}).get(str(mt), {})
        macs = [MacsPoint(float(p["kT_keV"]), float(p["sigma_mb"]), p.get("rate_cm3_mol_s"))
                for p in raw.get("macs", [])]
        spec = self._load_csv(key, mt)
        if not macs and not spec:
            return None
        q = raw.get("Q_MeV")
        return NeutronReactionData(mt, f"MT={mt} (user)",
                                   float(q) if q is not None else None,
                                   spectrum=spec, macs=macs, source="user")

    def has_reaction(self, Z, X, A, mt, meta=None):
        key = self._key(Z, X, A, meta)
        if self._iso.get(key, {}).get("neutrons", {}).get(str(mt)):
            return True
        return (self._dir / "spectra" / f"{key}_mt{mt}.csv").is_file()

    def _load_csv(self, key, mt) -> List[CrossSectionPoint]:
        p = self._dir / "spectra" / f"{key}_mt{mt}.csv"
        if not p.is_file():
            return []
        out = []
        with open(p, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    out.append(CrossSectionPoint(float(row["E_eV"]), float(row["sigma_barn"])))
                except (KeyError, ValueError):
                    continue
        return out

    # ── override ──
    def get_override(self, Z, X, A, field, meta=None) -> Optional[Any]:
        return self._ovr.get(self._key(Z, X, A, meta), {}).get(field)

    def list_isotopes(self):
        return list(self._iso.keys())

    def __repr__(self):
        return f"UserDataReader({self._dir.name}, isotopes={len(self._iso)})"
