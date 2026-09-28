"""
core — domain layer: entities, graph assembly, burnup calculation.

Key principle: BurnupMatrix accepts ANY list of isotopes —
from databases or assembled manually. The matrix does not touch files.
"""

from .entities import ChemTable, Element, Isotope, DecayLink, ReactionLink
from .registry import IsotopeBuilder
from .burnup import BurnupConfig, BurnupMatrix, equilibrium_heat_vs_flux
from .plotting import plot_concentrations, plot_heat, plot_heat_vs_flux
from .cycle_finder import CycleAnalyzer, IsotopeCycle, CycleStep, SCC

__all__ = [
    "ChemTable", "Element", "Isotope", "DecayLink", "ReactionLink",
    "IsotopeBuilder", "BurnupConfig", "BurnupMatrix", "equilibrium_heat_vs_flux",
    "plot_concentrations", "plot_heat", "plot_heat_vs_flux",
    "CycleAnalyzer", "IsotopeCycle", "CycleStep", "SCC",
]
