"""
nuclear_data/endf_reader.py
===========================
Unified reader for a single isotope's ENDF data.

Combines decay reading (dec-*.endf, MF=8/MT=457) and neutron-reaction
reading (n-*.endf, MF=3 cross sections, MF=9/10 isomeric branchings)
into ONE instance class.

Design:
- One isotope = one EndfReader object. No global mutable state.
- The decay file and the neutron file are looked up independently. Either
  may be absent: short-lived nuclides often have a decay file but no
  evaluated neutron file. Missing files are tolerated, not fatal.
- Depends on the `endf` library. Import never fails; the error is raised
  only when an instance actually needs a file that requires `endf`.

Public API:
    Decay:
        get_decay_data()  -> DecayData
        get_mass(to_kg)   -> float | None
    Neutron reactions:
        available_mts()                 -> List[int]
        get_cross_section(mt, ...)      -> List[CrossSectionPoint]
        get_sigma_at(mt, E, e_unit)     -> float | None
        get_reaction_data(mt)           -> NeutronReactionData | None
        get_isomeric_branches(mt, E_eV) -> List[IsomericBranch]
        compute_macs(mt, kT, kT_unit)   -> MacsPoint | None
"""

from __future__ import annotations
import math
from pathlib import Path
from typing import Optional, Union, List, Tuple, Dict

from .models import (
    DecayData, DecayChannel,
    CrossSectionPoint, IsomericBranch, MacsPoint, NeutronReactionData,
)

try:
    import endf  # type: ignore
    _HAS_ENDF = True
except ImportError:
    endf = None
    _HAS_ENDF = False


class EndfReader:
    """Unified ENDF reader for one isotope (decay + neutron reactions)."""

    # Decay-mode decoding table for the RTYP field (MF=8/MT=457).
    _PARTICLE = {
        "1":  {"name": "beta-",    "dz": 1,  "da": 0},
        "2":  {"name": "ec/beta+", "dz": -1, "da": 0},
        "3":  {"name": "it",       "dz": 0,  "da": 0},
        "4":  {"name": "alpha",    "dz": -2, "da": -4},
        "5":  {"name": "n",        "dz": 0,  "da": -1},
        "6":  {"name": "sf",       "dz": 0,  "da": 0},
        "7":  {"name": "p",        "dz": -1, "da": -1},
        "10": {"name": "unknown",  "dz": 0,  "da": 0},
    }

    # MT code -> human-readable channel name.
    MT_NAMES: Dict[int, str] = {
        1: "(n,total)", 2: "(n,el)", 4: "(n,inl)", 16: "(n,2n)", 17: "(n,3n)",
        18: "fission", 91: "(n,n')", 102: "(n,γ)", 103: "(n,p)", 104: "(n,d)",
        105: "(n,t)", 106: "(n,3He)", 107: "(n,α)",
    }

    def __init__(self, Z, X, A, meta=None, dec_dir=None, neu_dir=None):
        """
        Z, X, A, meta : isotope identity.
        dec_dir       : folder with dec-*.endf files (decay). Optional.
        neu_dir       : folder with n-*.endf files (neutrons). Optional.

        At least one of dec_dir / neu_dir should be given. Each file is
        loaded lazily-tolerant: if the file is missing, that half of the
        data is simply unavailable (methods return empty/None), but the
        object is still constructed.
        """
        if not _HAS_ENDF:
            raise ImportError("The 'endf' library is not installed: pip install endf")
        self._Z, self._X, self._A, self._meta = int(Z), str(X).capitalize(), int(A), meta

        # Decay material (MF=8/457).
        self._dec_mat = None
        self._dec_path = None
        if dec_dir is not None:
            fname = self._decay_filename(Z, X, A, meta)
            full = (Path(dec_dir) / fname).resolve()
            if full.exists():
                self._dec_path = full
                self._dec_mat = endf.Material(str(full))

        # Neutron material (MF=3, 9, 10). Raw text is cached for manual parsing.
        self._neu_mat = None
        self._neu_path = None
        self._lines: List[str] = []
        if neu_dir is not None:
            fname = self._neutron_filename(Z, X, A, meta)
            full = (Path(neu_dir) / fname).resolve()
            if full.exists():
                self._neu_path = full
                self._neu_mat = endf.Material(str(full))
                self._lines = full.read_text(
                    encoding="utf-8", errors="replace").splitlines()

    # ------------------------------------------------------------------
    # Filename / identity helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _fmt_meta(meta):
        """Format the metastable suffix: 1 -> 'm1', 'g' -> '', None -> ''."""
        if meta is None or meta == "" or meta == 0:
            return ""
        if isinstance(meta, int):
            return f"m{meta}"
        s = str(meta).strip().lower()
        if s == "g":
            return ""
        return ("m" + s[1:]) if s.startswith("m") else s

    @staticmethod
    def _decay_filename(Z, X, A, meta):
        m = EndfReader._fmt_meta(meta)
        return f"dec-{int(Z):03d}_{str(X).capitalize()}_{int(A):03d}{m}.endf"

    @staticmethod
    def _neutron_filename(Z, X, A, meta):
        m = EndfReader._fmt_meta(meta)
        return f"n-{int(Z):03d}_{str(X).capitalize()}_{int(A):03d}{m}.endf"

    @property
    def has_decay(self) -> bool:
        return self._dec_mat is not None

    @property
    def has_neutron(self) -> bool:
        return self._neu_mat is not None

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _first_of(obj, names, default=None):
        """Return the first attribute/key found on obj from a name list."""
        for n in names:
            if hasattr(obj, n):
                return getattr(obj, n)
            if isinstance(obj, dict) and n in obj:
                return obj[n]
        return default

    @staticmethod
    def _scalar(raw):
        """Extract a number. endf 0.1.12 often stores (value, uncertainty)."""
        if isinstance(raw, (tuple, list)):
            return float(raw[0])
        return float(raw)

    @staticmethod
    def _endf_float(s):
        """Parse an 11-char ENDF float such as ' 1.000000-5' into float."""
        s = s.strip()
        if not s:
            return None
        if "e" not in s.lower():
            s = s[0] + s[1:].replace("+", "e+").replace("-", "e-")
        try:
            return float(s)
        except ValueError:
            return None

    @staticmethod
    def _fields6(line):
        """Split an 80-char ENDF line into 6 fields of 11 chars."""
        line = line.ljust(80)
        return [EndfReader._endf_float(line[i:i + 11]) for i in range(0, 66, 11)]

    @staticmethod
    def _mf_mt(line):
        """Extract (MF, MT) from the tail of an ENDF line (cols 70-75)."""
        line = line.ljust(80)
        try:
            return int(line[70:72]), int(line[72:75])
        except Exception:
            return None, None

    @staticmethod
    def _interp(pairs, x):
        """Linear interpolation over a list of (x, y) pairs."""
        if not pairs:
            return None
        pairs = sorted(pairs, key=lambda t: t[0])
        if x <= pairs[0][0]:
            return pairs[0][1]
        if x >= pairs[-1][0]:
            return pairs[-1][1]
        for i in range(len(pairs) - 1):
            x0, y0 = pairs[i]; x1, y1 = pairs[i + 1]
            if x0 <= x <= x1:
                return y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        return None

    # ==================================================================
    # DECAY DATA (MF=8 / MT=457)
    # ==================================================================
    def _decode_rtyp(self, rtyp):
        """Decode an RTYP float (e.g. 1.4 = beta- then alpha) into name, dZ, dA."""
        s = f"{float(rtyp):.6f}".rstrip("0").rstrip(".")
        parts = s.split(".")
        codes = [parts[0]] + list(parts[1] if len(parts) > 1 else "")
        names, dz, da = [], 0, 0
        for c in codes:
            p = self._PARTICLE.get(c)
            if p:
                names.append(p["name"]); dz += p["dz"]; da += p["da"]
        return (", ".join(names) or "unknown"), dz, da

    def get_decay_data(self) -> DecayData:
        """Read decay data. Returns a stable DecayData if no decay file/section."""
        if self._dec_mat is None:
            return DecayData(None, 0.0, source="endf")
        sec = self._dec_mat.section_data.get((8, 457))
        if sec is None:
            return DecayData(None, 0.0, source="endf")

        # Half-life: the real key in endf 0.1.12 is 'T1/2', value is a
        # (T, dT) tuple. NST=1 marks a stable nuclide.
        nst = self._first_of(sec, ["NST"], 0)
        try:
            nst = int(self._scalar(nst))
        except Exception:
            nst = 0
        raw_hl = self._first_of(sec, ["T1/2", "half_life", "T", "T12", "T_half"])
        try:
            t = self._scalar(raw_hl) if raw_hl is not None else None
        except Exception:
            t = None
        if nst == 1 or not t or t <= 0:
            return DecayData(None, 0.0, source="endf")

        lam = math.log(2) / t
        channels: List[DecayChannel] = []
        for m in self._first_of(sec, ["modes", "decay_modes", "decay"], []):
            rtyp = self._scalar(self._first_of(m, ["RTYP", "rtyp"], 0.0))
            rfs = int(self._scalar(self._first_of(m, ["RFS", "rfs"], 0.0)))
            q = self._scalar(self._first_of(m, ["Q", "q", "q_value"], 0.0))
            br = self._scalar(self._first_of(m, ["BR", "br", "branching_ratio"], 0.0))
            name, dz, da = self._decode_rtyp(rtyp)
            # Q is stored in eV in the file; convert to MeV.
            channels.append(DecayChannel(name, br, q / 1e6, dz, da, rfs))
        return DecayData(t, lam, channels, "endf")

    def get_mass(self, to_kg=False):
        """Atomic mass (AWR-based). In amu by default, kg if to_kg=True."""
        sec = None
        if self._dec_mat is not None:
            sec = self._dec_mat.section_data.get((1, 451)) \
                or self._dec_mat.section_data.get((8, 457))
        if sec is None and self._neu_mat is not None:
            sec = self._neu_mat.section_data.get((1, 451))
        raw = self._first_of(sec, ["AWR", "mass", "aw", "atomic_mass"]) if sec else None
        try:
            mass = self._scalar(raw)
        except Exception:
            return None
        return mass * 1.66053906660e-27 if to_kg else mass

    # ==================================================================
    # NEUTRON CROSS SECTIONS (MF=3)
    # ==================================================================
    def _try_section(self, mt):
        """Extract (E, sigma) arrays for an MF=3 MT section, several layouts."""
        if self._neu_mat is None:
            return None, None
        sec = self._neu_mat.section_data.get((3, mt))
        if sec is None:
            return None, None
        E = XS = None
        # Layout 1: raw ENDF text (list of strings).
        if isinstance(sec, (list, tuple)) and sec and isinstance(sec[0], str):
            try:
                t = self._fields6(sec[1]); NR, NP = int(t[4]), int(t[5])
                start = 2 + (NR * 2 + 5) // 6
                el, xl = [], []
                for ln in sec[start:]:
                    v = self._fields6(ln)
                    for i in range(0, 6, 2):
                        if v[i] is not None and v[i + 1] is not None:
                            el.append(v[i]); xl.append(v[i + 1])
                        if len(el) == NP:
                            break
                    if len(el) == NP:
                        break
                E, XS = el, xl
            except Exception:
                pass
        # Layout 2: dict with a Tabulated1D under 'sigma' (endf 0.1.12).
        elif isinstance(sec, dict):
            t1 = sec.get("sigma")
            if t1 is not None and hasattr(t1, "x"):
                E, XS = t1.x, t1.y
            if E is None:
                E = self._first_of(sec, ["E", "x", "energy"])
                XS = self._first_of(sec, ["xs", "sigma", "y"])
        # Layout 3: a bare Tabulated1D-like object.
        elif hasattr(sec, "x") and hasattr(sec, "y"):
            E, XS = sec.x, sec.y

        if E is None or XS is None or len(E) != len(XS) or not len(E):
            return None, None
        # Keep strictly increasing energies (drops duplicate-edge points
        # seen in some libraries such as JENDL).
        Ee, Xs, last = [], [], -1.0
        for e, s in zip(E, XS):
            ef = float(e)
            if ef > last:
                Ee.append(ef); Xs.append(float(s)); last = ef
        return Ee, Xs

    def _cs_pairs(self, mt, e_unit="eV", cs_unit="barn"):
        """Return [(E, sigma)] in requested units, or None."""
        ec = {"ev": 1.0, "kev": 1e-3, "mev": 1e-6}
        cc = {"barn": 1.0, "mbarn": 1e3, "mb": 1e3, "cm2": 1e-24}
        ef, sf = ec[e_unit.lower()], cc[cs_unit.lower()]
        E, XS = self._try_section(mt)
        if E is not None:
            return [(e * ef, s * sf) for e, s in zip(E, XS)]
        # Manual TAB1 fallback from cached raw lines.
        i, n = 0, len(self._lines)
        while i < n:
            mf, mtf = self._mf_mt(self._lines[i])
            if mf == 3 and mtf == mt:
                head = self._fields6(self._lines[i])
                NR = int(round(head[4] or 0)); NP = int(round(head[5] or 0))
                j = i + 1 + math.ceil((2 * NR) / 6)
                xy = []
                while len(xy) < 2 * NP and j < n:
                    xy.extend(self._fields6(self._lines[j])); j += 1
                return [(xy[k] * ef, xy[k + 1] * sf) for k in range(0, 2 * NP, 2)
                        if k + 1 < len(xy) and xy[k] is not None and xy[k + 1] is not None] or None
            i += 1
        return None

    def _q_value_ev(self, mt):
        """Reaction Q value in eV from the MF=3 section header."""
        if self._neu_mat is not None:
            sec = self._neu_mat.section_data.get((3, mt))
            if sec is not None:
                raw = self._first_of(sec, ["QI", "qi", "Q", "QM", "C2"])
                if raw is not None:
                    try:
                        return float(raw)
                    except Exception:
                        pass
        for ln in self._lines:
            mf, mtf = self._mf_mt(ln)
            if mf == 3 and mtf == mt:
                head = self._fields6(ln)
                for v in head[1::-1]:
                    if v is not None:
                        try:
                            return float(v)
                        except Exception:
                            pass
                break
        return None

    # ==================================================================
    # ISOMERIC BRANCHING (MF=9 / MF=10)  ->  the Q coefficient in the matrix
    # ==================================================================
    def _parse_mf9_10(self, mf, mt):
        """Parse MF=9 (multiplicities) or MF=10 (partial cross sections)."""
        if self._neu_mat is None:
            return None
        sec = self._neu_mat.section_data.get((mf, mt))
        if not sec:
            return None
        # Library already parsed it into a dict of levels.
        if isinstance(sec, dict) and "levels" in sec:
            out = []
            for lvl in sec["levels"]:
                LFS = int(round(self._first_of(lvl, ["LFS"]) or 0))
                so = self._first_of(lvl, ["sigma", "xs", "y"])
                E = self._first_of(so, ["x", "E"]); XS = self._first_of(so, ["y", "xs"])
                pr = []
                try:
                    if hasattr(E, "__len__") and len(E) == len(XS):
                        pr = [(float(e), float(s)) for e, s in zip(E, XS)]
                except Exception:
                    pass
                out.append({"LFS": LFS, "pairs": pr})
            return out
        # Raw ENDF text layout.
        if not (isinstance(sec, (list, tuple)) and sec and isinstance(sec[0], str)):
            return None
        out = []
        head = self._fields6(sec[0]); NS = int(round(head[4] if head[4] is not None else 1))
        i = 1
        for _ in range(NS):
            if i >= len(sec):
                break
            th = self._fields6(sec[i])
            LFS = int(round(th[3] if th[3] is not None else 0))
            NR = int(round(th[4] if th[4] is not None else 0))
            NP = int(round(th[5] if th[5] is not None else 0))
            i += 1 + (2 * NR + 5) // 6
            xy = []
            for _ in range((2 * NP + 5) // 6):
                if i < len(sec):
                    xy.extend(self._fields6(sec[i])); i += 1
            pr = [(xy[k], xy[k + 1]) for k in range(0, 2 * NP, 2)
                  if k + 1 < len(xy) and xy[k] is not None and xy[k + 1] is not None]
            out.append({"LFS": LFS, "pairs": pr})
        return out

    def get_isomeric_branches(self, mt, E_eV) -> List[IsomericBranch]:
        """
        Fraction of reactions going into each final state (yield = Q_ij).
        Falls back to [LFS=0, yield=1.0] when no MF=9/10 data is present.
        """
        for mf in (9, 10):
            blocks = self._parse_mf9_10(mf, mt)
            if not blocks:
                continue
            raw, tot = [], 0.0
            for b in blocks:
                y = self._interp(b["pairs"], E_eV)
                if y and y > 0:
                    raw.append((b["LFS"], y)); tot += y
            if tot > 0:
                return [IsomericBranch(L, y / tot) for L, y in raw]
        return [IsomericBranch(0, 1.0)]

    # ==================================================================
    # PUBLIC NEUTRON API
    # ==================================================================
    def available_mts(self):
        """List of known MT channels present as MF=3 sections."""
        if self._neu_mat is None:
            return []
        return sorted(mt for (mf, mt) in self._neu_mat.section_data
                      if mf == 3 and mt in self.MT_NAMES)

    def get_cross_section(self, mt, e_unit="eV", cs_unit="barn"):
        """Full sigma(E) spectrum as a list of CrossSectionPoint."""
        pairs = self._cs_pairs(mt, e_unit, cs_unit)
        return [CrossSectionPoint(e, s) for e, s in pairs] if pairs else []

    def get_sigma_at(self, mt, E, e_unit="eV"):
        """Interpolated sigma [barn] at a single energy E."""
        conv = {"ev": 1.0, "kev": 1e3, "mev": 1e6}
        pairs = self._cs_pairs(mt, e_unit="eV")
        return self._interp(pairs, float(E) * conv[e_unit.lower()]) if pairs else None

    def get_reaction_data(self, mt):
        """Full NeutronReactionData (spectrum + Q value) for one MT."""
        spec = self.get_cross_section(mt)
        if not spec:
            return None
        q = self._q_value_ev(mt)
        return NeutronReactionData(
            mt=mt, channel_name=self.MT_NAMES.get(mt, f"MT={mt}"),
            Q_MeV=(q / 1e6 if q is not None else None), spectrum=spec, source="endf")

    def compute_macs(self, mt, kT=30.0, kT_unit="keV"):
        """
        Maxwellian-averaged cross section by integrating the ENDF spectrum
        over a Maxwell-Boltzmann distribution at temperature kT.
        Returns a MacsPoint with sigma in mbarn.
        """
        conv = {"ev": 1.0, "kev": 1e3, "mev": 1e6}
        kT_eV = float(kT) * conv[kT_unit.lower()]
        if kT_eV <= 0:
            raise ValueError("kT must be positive")
        pairs = self._cs_pairs(mt, e_unit="eV", cs_unit="barn")
        if not pairs:
            return None
        inv = 1.0 / kT_eV
        pref = 2.0 * inv ** 1.5 / math.sqrt(math.pi)
        tw, tws = 0.0, 0.0
        for e, s in pairs:
            if e <= 0:
                continue
            w = pref * math.sqrt(e) * math.exp(-e * inv)
            tw += w; tws += w * s
        if tw <= 0:
            return None
        return MacsPoint(kT_eV / 1e3, (tws / tw) * 1e3)

    def __repr__(self):
        flags = []
        if self.has_decay:
            flags.append("decay")
        if self.has_neutron:
            flags.append("neutron")
        return f"EndfReader({self._Z}-{self._X}-{self._A}, {'+'.join(flags) or 'empty'})"


# ═══════════════════════════════════════════════════════════════════════════
#  Backward-compat wrappers (old API: separate folder-first constructors)
#
#  decay_reader.EndfDecayReader(folder, Z, X, A)  →  EndfDecayReader below
#  nreac_reader.EndfNReacReader(folder, Z, X, A)  →  EndfNReacReader below
# ═══════════════════════════════════════════════════════════════════════════

class EndfDecayReader:
    """Wrapper preserving the old decay_reader API for backward compatibility."""

    def __init__(self, folder, Z, X, A, meta=None):
        self._r = EndfReader(Z, X, A, meta, dec_dir=folder)

    def get_decay_data(self) -> DecayData:
        return self._r.get_decay_data()

    def get_mass(self, to_kg=False):
        return self._r.get_mass(to_kg)

    def __repr__(self):
        return repr(self._r).replace("EndfReader", "EndfDecayReader")


class EndfNReacReader:
    """Wrapper preserving the old nreac_reader API for backward compatibility."""

    MT_NAMES = EndfReader.MT_NAMES

    def __init__(self, folder, Z, X, A, meta=None):
        self._r = EndfReader(Z, X, A, meta, neu_dir=folder)

    def available_mts(self) -> List[int]:
        return self._r.available_mts()

    def get_cross_section(self, mt, e_unit="eV", cs_unit="barn") -> List[CrossSectionPoint]:
        return self._r.get_cross_section(mt, e_unit, cs_unit)

    def get_sigma_at(self, mt, E, e_unit="eV") -> Optional[float]:
        return self._r.get_sigma_at(mt, E, e_unit)

    def get_reaction_data(self, mt) -> Optional[NeutronReactionData]:
        return self._r.get_reaction_data(mt)

    def get_isomeric_branches(self, mt, E_eV) -> List[IsomericBranch]:
        return self._r.get_isomeric_branches(mt, E_eV)

    def compute_macs(self, mt, kT=30.0, kT_unit="keV") -> Optional[MacsPoint]:
        return self._r.compute_macs(mt, kT, kT_unit)

    def __repr__(self):
        return repr(self._r).replace("EndfReader", "EndfNReacReader")
