"""
core/entities.py
================
Domain entities: element, isotope, decay and reaction links.

Key idea: Isotope is a self-contained object. It stores the decay constant,
decay channels, and neutron reactions. Where the data came from
(ENDF / MACS / manual) is unknown and irrelevant to the isotope.

This allows building an isotope list MANUALLY and passing it to the burnup
matrix, bypassing any database entirely.
"""

from __future__ import annotations
from typing import List, Optional, Union


class ChemTable:
    elements = {
        1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
        11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
        21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
        31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
        41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
        51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
        61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
        71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
        81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
        91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
        101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt",
        110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og",
    }

    @staticmethod
    def Z_of(symbol: str) -> int:
        for z, s in ChemTable.elements.items():
            if s == symbol:
                return z
        raise ValueError(f"Unknown element symbol: {symbol}")


class Element:
    __slots__ = ("Z", "X")

    def __init__(self, Z: int, X: str):
        if Z not in ChemTable.elements:
            raise ValueError(f"Z={Z} not in the periodic table")
        if ChemTable.elements[Z] != X:
            raise ValueError(f"Z={Z} corresponds to {ChemTable.elements[Z]}, not {X}")
        self.Z, self.X = Z, X

    def __eq__(self, o):
        return isinstance(o, Element) and (self.Z, self.X) == (o.Z, o.X)

    def __hash__(self):
        return hash((self.Z, self.X))

    def __repr__(self):
        return f"Element({self.Z}-{self.X})"


class DecayLink:
    """Decay channel: target and branching fraction."""
    __slots__ = ("mode", "branch", "Q_MeV", "dz", "da", "rfs", "product")

    def __init__(self, mode, branch, Q_MeV, dz, da, rfs=0, product=None):
        self.mode = mode
        self.branch = float(branch)      # P(ij)
        self.Q_MeV = float(Q_MeV)
        self.dz = int(dz)
        self.da = int(da)
        self.rfs = int(rfs)
        self.product: Optional["Isotope"] = product

    def __repr__(self):
        return f"DecayLink({self.mode}, BR={self.branch:.3g}, dZ={self.dz}, dA={self.da})"


class ReactionLink:
    """
    Neutron reaction: MT, cross section [barn] at the working energy,
    isomeric yield Q (yield), reaction Q value [MeV], product.
    """
    __slots__ = ("mt", "name", "sigma_barn", "q_yield", "Q_MeV", "lfs", "product",
                 "source")

    def __init__(self, mt, name, sigma_barn, q_yield=1.0, Q_MeV=0.0, lfs=0, product=None,
                 source=None):
        self.mt = int(mt)
        self.name = name
        self.sigma_barn = float(sigma_barn)   # sigma at the selected energy
        self.q_yield = float(q_yield)         # Q(ij) — isomeric branching
        self.Q_MeV = float(Q_MeV)             # reaction energy (for heat)
        self.lfs = int(lfs)
        self.product: Optional["Isotope"] = product
        self.source = source                  # data source that supplied sigma

    def getId(self):
        return self.mt

    def __repr__(self):
        return f"ReactionLink(MT={self.mt} {self.name}, sigma={self.sigma_barn:.3g} b, Q_yield={self.q_yield:.2f})"


class Isotope(Element):
    """
    Self-contained isotope. Can be created manually:

        u235 = Isotope(92, "U", 235)
        u235.decay_constant = 3.1e-17
        u235.add_decay(DecayLink("alpha", 1.0, 4.68, -2, -4))
        u235.add_reaction(ReactionLink(102, "(n,g)", 98.5e-3))
    """
    __slots__ = ("A", "meta", "half_life_s", "decay_constant", "mass",
                 "_decays", "_reactions")

    def __init__(self, Z, X, A, meta=None):
        super().__init__(Z, X)
        if A <= 0:
            raise ValueError("A must be > 0")
        self.A = int(A)
        self.meta = meta
        self.half_life_s: Optional[float] = None
        self.decay_constant: float = 0.0
        self.mass: Optional[float] = None
        self._decays: List[DecayLink] = []
        self._reactions: List[ReactionLink] = []

    @classmethod
    def from_symbol(cls, X, A, meta=None):
        return cls(ChemTable.Z_of(X), X, A, meta)

    # links
    def add_decay(self, link: DecayLink):
        self._decays.append(link)

    def add_reaction(self, link: ReactionLink):
        self._reactions.append(link)

    def getListOfDecays(self) -> List[DecayLink]:
        return self._decays

    def getListOfReactions(self) -> List[ReactionLink]:
        return self._reactions

    @property
    def name(self) -> str:
        m = f"m{self.meta}" if self.meta else ""
        return f"{self.X}-{self.A}{m}"

    def __eq__(self, o):
        return isinstance(o, Isotope) and (self.Z, self.A, self.meta) == (o.Z, o.A, o.meta)

    def __hash__(self):
        return hash((self.Z, self.A, self.meta))

    def __repr__(self):
        return f"Isotope({self.name}, lambda={self.decay_constant:.3g})"
