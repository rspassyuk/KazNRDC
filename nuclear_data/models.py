"""
nuclear_data/models.py
======================
Standardized data structures returned by the readers.
Any consumer works with these objects without knowing the data source.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List


# ───────────────────────────── Decay ─────────────────────────────

@dataclass
class DecayChannel:
    """One radioactive decay channel."""
    mode: str            # "alpha", "beta-", "ec/beta+", "it", "sf", "n", "p", ...
    branch: float        # branching fraction [0..1] — P(ij) in the matrix
    Q_MeV: float         # decay Q value [MeV] (for heat release)
    dz: int              # change in Z (protons)
    da: int              # change in A (mass number)
    rfs: int = 0         # isomeric state of daughter nucleus (0 = ground)


@dataclass
class DecayData:
    """Complete decay data for one isotope."""
    half_life_s: Optional[float]              # None = stable
    decay_constant: float                     # lambda = ln2 / T½ [s⁻¹]; 0.0 = stable
    channels: List[DecayChannel] = field(default_factory=list)
    source: str = "unknown"                   # "endf" | "user"

    @property
    def is_stable(self) -> bool:
        return self.half_life_s is None or self.decay_constant == 0.0


# ──────────────────────── Neutron reactions ───────────────────────

@dataclass
class CrossSectionPoint:
    """One point of the sigma(E) table."""
    E_eV: float
    sigma_barn: float


@dataclass
class MacsPoint:
    """One point of the MACS(kT) table."""
    kT_keV: float
    sigma_mb: float                           # MACS [millibarn]
    rate_cm3_mol_s: Optional[float] = None    # reaction rate [cm³/(mol·s)]


@dataclass
class IsomericBranch:
    """
    Fraction of the reaction going into a specific isomeric state of the product.
    yield_ is the Q(ij) coefficient in the burnup matrix (NOT an energy).
    """
    LFS: int             # final state level (0 = ground state)
    yield_: float        # fraction [0..1]


@dataclass
class NeutronReactionData:
    """Complete neutron reaction data for one isotope and one MT."""
    mt: int
    channel_name: str
    Q_MeV: Optional[float]                                  # reaction Q value [MeV]
    spectrum: List[CrossSectionPoint] = field(default_factory=list)
    macs: List[MacsPoint]            = field(default_factory=list)
    branches: List[IsomericBranch]   = field(default_factory=list)
    source: str = "unknown"          # "endf" | "macs" | "user"
