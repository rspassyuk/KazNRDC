"""
nuclear_data — nuclear data access layer.

Knows about file formats (ENDF, MACS, user database),
but has no knowledge of burnup physics. Returns standardized
dataclasses from models.py.
"""

from .models import (
    DecayChannel, DecayData,
    CrossSectionPoint, MacsPoint, IsomericBranch, NeutronReactionData,
)
from .path_holder import EndfPathHolder
from .endf_reader import EndfReader, _HAS_ENDF
from .macs_reader import MacsReader, MultiMacsReader
from .user_reader import UserDataReader
from .talys_reader import NreacTalys
from .provider import NuclearDataProvider, MACS_DEFAULT_PRIORITY

__all__ = [
    "DecayChannel", "DecayData",
    "CrossSectionPoint", "MacsPoint", "IsomericBranch", "NeutronReactionData",
    "EndfPathHolder",
    "EndfReader", "_HAS_ENDF",
    "MacsReader", "MultiMacsReader",
    "UserDataReader", "NreacTalys",
    "NuclearDataProvider", "MACS_DEFAULT_PRIORITY",
]
