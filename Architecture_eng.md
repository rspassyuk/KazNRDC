> For the maintained publication overview and current solver/averaging behavior, see [Architecture](docs/architecture.md) and [Methods](docs/methods.md). This earlier detailed API reference may describe older interfaces.

# NuMatRx Architecture Reference

Complete technical reference: every module, class, method, and call example.
Goal: understand the full system and run any method independently without the GUI.

---

## Table of Contents

1. [Layer map](#1-layer-map)
2. [nuclear_data/ — data layer](#2-nuclear_data--data-layer)
   - 2.1 [models.py — dataclasses](#21-modelspy--dataclasses)
   - 2.2 [path_holder.py — EndfPathHolder](#22-path_holderpy--endfpathholder)
   - 2.3 [endf_reader.py — EndfReader](#23-endf_readerpy--endfreader)
   - 2.4 [macs_reader.py — MacsReader / MultiMacsReader](#24-macs_readerpy--macsreader--multimacs-reader)
   - 2.5 [talys_reader.py — NreacTalys](#25-talys_readerpy--nreactalys)
   - 2.6 [user_reader.py — UserDataReader](#26-user_readerpy--userdatareader)
   - 2.7 [provider.py — NuclearDataProvider](#27-providerpy--nucleardataprovider)
3. [core/ — physics layer](#3-core--physics-layer)
   - 3.1 [entities.py — ChemTable, Element, Isotope, DecayLink, ReactionLink](#31-entitiespy)
   - 3.2 [registry.py — IsotopeBuilder](#32-registrypy--isotopebuilder)
   - 3.3 [burnup.py — BurnupConfig, BurnupMatrix](#33-burnuppy--burnupconfig-and-burnupmatrix)
   - 3.4 [cycle_finder.py — CycleAnalyzer](#34-cycle_finderpy--cycleanalyzer)
   - 3.5 [sensitivity.py — SensitivityEngine](#35-sensitivitypy--sensitivityengine)
   - 3.6 [plotting.py — plot helpers](#36-plottingpy--plot-helpers)
   - 3.7 [service.py — NuclearService](#37-servicepy--nuclearservice)
4. [GUI layer (model/)](#4-gui-layer-model)
5. [Interaction map](#5-interaction-map)
6. [Full call examples](#6-full-call-examples)
7. [Algorithm: cycle detection](#7-algorithm-cycle-detection)
8. [Algorithm: sensitivity analysis](#8-algorithm-sensitivity-analysis)

---

## 1. Layer map

```
┌─────────────────────────────────────────────────────────────────────┐
│  GUI  (model/)                                                      │
│  model_CoreNucleo · model_IsotopeChart · model_results · model_graph│
└────────────────────────────┬────────────────────────────────────────┘
                             │  calls
┌────────────────────────────▼────────────────────────────────────────┐
│  Service  (core/service.py)  — Qt-free singleton                    │
│  Bridges GUI ↔ domain; manages xsdir, universe, matrix, results     │
└──────┬────────────────────────────────────────┬─────────────────────┘
       │                                        │
┌──────▼──────────────────────┐  ┌─────────────▼──────────────────────┐
│  Core  (core/)              │  │  Nuclear data  (nuclear_data/)      │
│  entities · registry        │  │  EndfReader · MacsReader            │
│  BurnupMatrix · BurnupConfig│  │  MultiMacsReader · NreacTalys       │
│  CycleAnalyzer              │  │  UserDataReader                     │
│  SensitivityEngine          │  │  NuclearDataProvider                │
│  plotting                   │  │  EndfPathHolder                     │
└─────────────────────────────┘  └────────────────────────────────────┘
```

**Design rules**
- `core/` is DB-independent. `BurnupMatrix` accepts any list of `Isotope` objects.
- `nuclear_data/` is Qt-independent. It returns only standard dataclasses.
- `core/service.py` contains **no Qt imports** — returns plain Python objects.

---

## 2. nuclear_data/ — data layer

### 2.1 models.py — dataclasses

**File:** `nuclear_data/models.py`
Immutable dataclass objects passed between readers and consumers.

#### `DecayChannel`
```python
@dataclass
class DecayChannel:
    mode:   str    # "alpha" | "beta-" | "ec/beta+" | "it" | "sf" | "n" | "p"
    branch: float  # branching fraction [0..1]
    Q_MeV:  float  # decay energy [MeV]
    dz:     int    # change in Z
    da:     int    # change in A
    rfs:    int = 0  # daughter isomeric state (0 = ground)
```

#### `DecayData`
```python
@dataclass
class DecayData:
    half_life_s:    Optional[float]        # None = stable
    decay_constant: float                  # λ = ln2/T½ [s⁻¹]; 0 = stable
    channels:       List[DecayChannel]
    source:         str = "unknown"        # "endf6" | "user"

    @property
    def is_stable(self) -> bool: ...
```

#### `MacsPoint`
```python
@dataclass
class MacsPoint:
    kT_keV:        float
    sigma_mb:      float                   # MACS [millibarn]
    rate_cm3_mol_s: Optional[float] = None
```

#### `IsomericBranch`
```python
@dataclass
class IsomericBranch:
    LFS:    int    # final state level (0 = ground)
    yield_: float  # fraction of total σ going to this state
```

#### `NeutronReactionData`
```python
@dataclass
class NeutronReactionData:
    mt:           int
    channel_name: str                      # "(n,g)", "(n,2n)", ...
    Q_MeV:        Optional[float]
    spectrum:     List[CrossSectionPoint]  # σ(E) table
    macs:         List[MacsPoint]
    branches:     List[IsomericBranch]
    source:       str = "unknown"
```

#### `CrossSectionPoint`
```python
@dataclass
class CrossSectionPoint:
    E_eV:      float
    sigma_barn: float
```

**Usage examples:**
```python
from nuclear_data.models import DecayChannel, DecayData, MacsPoint

ch = DecayChannel("beta-", 0.9999981, 0.0635, dz=1, da=0)
dd = DecayData(half_life_s=7.0e8, decay_constant=9.9e-10, channels=[ch])
print(dd.is_stable)          # False
print(dd.decay_constant)     # 9.9e-10

pt = MacsPoint(kT_keV=30.0, sigma_mb=1.5)
print(pt.sigma_mb)           # 1.5
```

---

### 2.2 path_holder.py — EndfPathHolder

**File:** `nuclear_data/path_holder.py`
Class-level (static) path registry. All readers call it to get their directory.

#### Class attributes (all Optional, initially None)
| Attribute | Type | Purpose |
|-----------|------|---------|
| `_xsdir` | `Path` | Root of the xsdir/ tree |
| `_default_dec_dir` | `str` | ENDF decay files directory |
| `_default_neu_dir` | `str` | ENDF neutron files directory |
| `_macs_file_path` | `Path` | Path to rawmacs.txt |
| `_talys_dir` | `Path` | TALYS root |
| `_endfb71_macs_dir` | `Path` | ENDF/B-VII.1 multi-kT MACS |
| `_eaf2010_macs_dir` | `Path` | EAF-2010 multi-kT MACS |
| `_user_db_dir` | `Path` | User database folder |

#### Methods (all static)

```python
EndfPathHolder.init_xsdir(xsdir_path: str) -> Dict[str, str]
```
Auto-discovers subdirectories from `xsdir/`. Raises `FileNotFoundError` if
`endf-6/` or `MACS/` are missing. Returns a dict of discovered paths.

```python
EndfPathHolder.set_DEFAULT_DEC_DIR(library: str) -> str
# library: "ENDFB-VIII.0" | "JEFF" | "JENDL" | ...
# Returns path to xsdir/endf-6/<library>/decay/

EndfPathHolder.get_DEFAULT_DEC_DIR() -> str

EndfPathHolder.set_DEFAULT_NEU_DIR(library: str) -> str
EndfPathHolder.get_DEFAULT_NEU_DIR() -> str

EndfPathHolder.set_DEFAULT_MACS_DIR(library_alias: str, filename="rawmacs.txt") -> str
# library_alias: "Endfb7" | "Jeff" | "Jendl" | ...
# Returns path to rawmacs.txt

EndfPathHolder.get_macs_file() -> Path
EndfPathHolder.get_macs_library() -> str

EndfPathHolder.set_ENDFB71_MACS_DIR(path=None) -> str
EndfPathHolder.get_ENDFB71_MACS_DIR() -> Path

EndfPathHolder.set_EAF2010_MACS_DIR(path=None) -> str
EndfPathHolder.get_EAF2010_MACS_DIR() -> Path

EndfPathHolder.set_TALYS_DIR(path=None) -> str
EndfPathHolder.get_TALYS_DIR() -> Path

EndfPathHolder.set_USER_DB_DIR(path: str) -> str
# Creates the directory and spectra/ subdirectory if absent.
EndfPathHolder.get_USER_DB_DIR() -> Path
```

**Usage example:**
```python
from nuclear_data import EndfPathHolder

paths = EndfPathHolder.init_xsdir("xsdir")
# {'xsdir': '.../xsdir', 'endf-6': '...', 'MACS': '...', 'TALYS': '...'}

dec_dir = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0")
neu_dir = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0")
macs    = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7")
```

---

### 2.3 endf_reader.py — EndfReader

**File:** `nuclear_data/endf_reader.py`
Requires `pip install endf`. If absent, `_HAS_ENDF = False` and all methods return `None`.

#### Constructor
```python
EndfReader(Z, X, A, meta=None, dec_dir=None, neu_dir=None)
```
Locates files:
- `dec_dir/dec-{Z:03d}_{X}_{A:03d}[m{meta}].endf`
- `neu_dir/n{Z:03d}-{X}-{A:03d}[m{meta}].endf`

Sets `has_decay: bool` and `has_neutron: bool`.

#### Methods

```python
EndfReader.get_decay_data() -> DecayData | None
```
Reads MF=8/MT=457. Returns λ, half-life, decay channels (mode/branch/Q/dz/da/rfs).

```python
EndfReader.get_sigma_at(mt: int, E, e_unit="eV") -> float | None
```
Interpolates σ(E) from MF=3/MT tabulation using ENDF interpolation law.

```python
EndfReader.compute_macs(mt: int, kT, kT_unit="keV") -> MacsPoint | None
```
Maxwell-Boltzmann integration over the MF=3 σ(E) grid.

```python
EndfReader.get_reaction_data(mt: int) -> NeutronReactionData | None
```
Returns full reaction record with Q-value, spectrum, and isomeric branches.

```python
EndfReader.get_isomeric_branches(mt: int, E_eV: float) -> List[IsomericBranch]
```
Reads MF=9/10 branching fractions. Returns `[IsomericBranch(0, 1.0)]` if absent.

---

### 2.4 macs_reader.py — MacsReader / MultiMacsReader

#### `MacsReader` (rawmacs.txt)
```python
MacsReader(filepath: str, lib_tag: str = "Endfb7")
MacsReader.get_macs_at(iso_key: str, mt: int, kT_keV: float) -> MacsPoint | None
```
- `iso_key` format: `"82-Pb-210"` (from `NuclearDataProvider._iso_str()`).
- Caches parsed file as `rawmacs.pkl` next to the source for fast reload.
- Supported MTs: 102, 16, 18, 103, 107.

#### `MultiMacsReader` (ENDF/B-VII.1, EAF-2010)
```python
MultiMacsReader(directory: str, source_tag: str = "")
MultiMacsReader.get_macs_at(iso_key: str, mt: int, kT_keV: float) -> MacsPoint | None
```
- Reads per-isotope files with multi-kT columns (5, 10, 15…100 keV).
- Linear interpolation between kT points.

**Direct usage example:**
```python
from nuclear_data.macs_reader import MacsReader, MultiMacsReader

mr = MacsReader("xsdir/MACS/rawmacs.txt", "Endfb7")
pt = mr.get_macs_at("82-Pb-210", mt=102, kT_keV=30.0)
if pt:
    print(f"MACS = {pt.sigma_mb} mb")

mmr = MultiMacsReader("xsdir/MACS/ENDFB71-EAF2010/ENDFB71")
pt2 = mmr.get_macs_at("82-Pb-210", mt=102, kT_keV=30.0)
```

---

### 2.5 talys_reader.py — NreacTalys

```python
NreacTalys(talys_dir: str, Z: int, X: str, A: int, meta=None)
NreacTalys.get_sigma_at(mt: int, E, e_unit="eV") -> float | None
NreacTalys.get_macs_30keV(mt: int) -> MacsPoint | None
NreacTalys.get_reaction_data(mt: int) -> NeutronReactionData | None
```
Reads TALYS output files:
- `recommended-ng-cross-sections/` for MT=102
- `recommended-na-cross-sections/` for MT=107
- `recommended-np-cross-sections/` for MT=103
- `talys_macs_30_keV/` for MACS at kT=30 keV

---

### 2.6 user_reader.py — UserDataReader

```python
UserDataReader(user_dir: str)
```
Reads `user_dir/isotopes.json`. Key format: `"82-Pb-210"` or `"82-Pb-210m1"`.

#### Methods
```python
UserDataReader.has_decay(Z, X, A, meta=None)     -> bool
UserDataReader.get_decay_data(Z, X, A, meta=None) -> DecayData | None
UserDataReader.has_reaction(Z, X, A, mt, meta=None) -> bool
UserDataReader.get_reaction_data(Z, X, A, mt, meta=None) -> NeutronReactionData | None
UserDataReader.get_override(Z, X, A, field, meta=None)   -> Any | None
UserDataReader.reload()                                   -> None
```

**isotopes.json schema:**
```json
{
  "82-Pb-210": {
    "decay": {
      "half_life_s": 700535040.0,
      "channels": [
        {"mode": "beta-", "branch": 0.9999981, "Q_MeV": 0.0635, "dz": 1, "da": 0, "rfs": 0}
      ]
    },
    "neutrons": {
      "102": {
        "Q_MeV": 3.82,
        "macs": [
          {"kT_keV": 30.0, "sigma_mb": 1.5},
          {"kT_keV": 90.0, "sigma_mb": 0.8}
        ]
      }
    }
  }
}
```

---

### 2.7 provider.py — NuclearDataProvider

**File:** `nuclear_data/provider.py`
Single façade over all data sources. Priority lists are configurable.

#### Constructor
```python
NuclearDataProvider(
    dec_dir=None,           # ENDF decay directory
    neu_dir=None,           # ENDF neutron directory
    macs_file=None,         # path to rawmacs.txt
    macs_lib=None,          # rawmacs column alias ("Endfb7")
    user_dir=None,          # user_db/ directory
    talys_dir=None,         # TALYS root
    endfb71_macs_dir=None,  # ENDF/B-VII.1 multi-kT MACS
    eaf2010_macs_dir=None,  # EAF-2010 multi-kT MACS
    macs_priority=None,     # List[str], default: MACS_DEFAULT_PRIORITY
    dec_priority=None,      # List[str], default: DEC_DEFAULT_PRIORITY
)
```

Default priorities:
```python
DEC_DEFAULT_PRIORITY  = ["user", "endf6"]
MACS_DEFAULT_PRIORITY = ["user", "rawmacs", "endfb71", "eaf2010", "talys", "endf_computed"]
```

#### Core methods

```python
# Decay data (follows dec_priority)
NuclearDataProvider.get_decay(Z, X, A, meta=None) -> DecayData | None

# MACS (follows macs_priority)
NuclearDataProvider.get_macs(Z, X, A, mt, kT_keV, meta=None) -> MacsPoint | None
NuclearDataProvider.get_macs_with_source(Z, X, A, mt, kT_keV, meta=None)
    -> Tuple[MacsPoint | None, str | None]  # (point, source_name)

# σ(E) point-wise (priority: user > ENDF > TALYS)
NuclearDataProvider.get_sigma(Z, X, A, mt, E, e_unit="eV", meta=None) -> float | None
NuclearDataProvider.get_sigma_with_source(Z, X, A, mt, E, e_unit="eV", meta=None)
    -> Tuple[float | None, str | None]

# Isomeric branching (ENDF only, falls back to [(LFS=0, yield=1.0)])
NuclearDataProvider.get_isomeric_branches(Z, X, A, mt, E_eV, meta=None) -> List[IsomericBranch]

# Full reaction record
NuclearDataProvider.get_reaction_data(Z, X, A, mt, meta=None) -> NeutronReactionData | None

# Priority management
NuclearDataProvider.set_macs_priority(priority: List[str]) -> None
NuclearDataProvider.get_macs_priority() -> List[str]
NuclearDataProvider.set_dec_priority(priority: List[str]) -> None
NuclearDataProvider.get_dec_priority() -> List[str]

# Utilities
NuclearDataProvider.get_override(Z, X, A, field, meta=None) -> Any | None
NuclearDataProvider.reload_user_db() -> None
```

**Complete provider setup (matches examples/):**
```python
from nuclear_data import EndfPathHolder, NuclearDataProvider

EndfPathHolder.init_xsdir("xsdir")

nd = NuclearDataProvider(
    user_dir         = EndfPathHolder.set_USER_DB_DIR("user_db"),
    dec_dir          = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0"),
    neu_dir          = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0"),
    macs_file        = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7"),
    macs_lib         = "Endfb7",
    endfb71_macs_dir = str(EndfPathHolder.get_ENDFB71_MACS_DIR()),
    eaf2010_macs_dir = str(EndfPathHolder.get_EAF2010_MACS_DIR()),
    talys_dir        = str(EndfPathHolder.get_TALYS_DIR()),
    macs_priority    = ["user", "endfb71", "eaf2010", "rawmacs", "talys", "endf_computed"],
)

dd = nd.get_decay(82, "Pb", 210)
print(f"T½ = {dd.half_life_s:.2e} s, λ = {dd.decay_constant:.3e} s⁻¹")

pt, src = nd.get_macs_with_source(82, "Pb", 210, mt=102, kT_keV=30.0)
print(f"MACS = {pt.sigma_mb} mb  [{src}]")
```

---

## 3. core/ — physics layer

### 3.1 entities.py

**File:** `core/entities.py`

#### `ChemTable`
Static class. Stores `{Z: symbol}` for Z=1..118.

```python
ChemTable.elements: Dict[int, str]   # {1: "H", 2: "He", ..., 118: "Og"}

ChemTable.Z_of(symbol: str) -> int   # "Pb" → 82; raises ValueError if unknown
```

**Examples:**
```python
from core.entities import ChemTable

z = ChemTable.Z_of("Pb")   # 82
sym = ChemTable.elements[83]  # "Bi"
```

---

#### `Element`
```python
class Element:
    __slots__ = ("Z", "X")
    def __init__(self, Z: int, X: str)
    # Raises ValueError if Z is out of table or Z/X mismatch

    def __eq__(self, other) -> bool
    def __hash__(self) -> int
    def __repr__(self) -> str   # "Element(82-Pb)"
```

**Examples:**
```python
from core.entities import Element

pb = Element(82, "Pb")
bi = Element(83, "Bi")
print(pb == bi)   # False
print(hash(pb))
```

---

#### `DecayLink`
```python
class DecayLink:
    __slots__ = ("mode", "branch", "Q_MeV", "dz", "da", "rfs", "product")
    def __init__(self, mode, branch, Q_MeV, dz, da, rfs=0, product=None)
    # product: Optional[Isotope] — set by IsotopeBuilder.link_products()
```

**Examples:**
```python
from core.entities import DecayLink

# Alpha decay of Po-210: dZ=-2, dA=-4
dl = DecayLink("alpha", 1.0, 5.407, dz=-2, da=-4)
print(dl)  # DecayLink(alpha, BR=1.0, dZ=-2, dA=-4)
```

---

#### `ReactionLink`
```python
class ReactionLink:
    __slots__ = ("mt", "name", "sigma_barn", "q_yield", "Q_MeV", "lfs", "product", "source")
    def __init__(self, mt, name, sigma_barn, q_yield=1.0, Q_MeV=0.0, lfs=0,
                 product=None, source=None)

    def getId(self) -> int   # returns self.mt (SBML-compatible interface)
```

**Examples:**
```python
from core.entities import ReactionLink

rx = ReactionLink(102, "(n,g)", sigma_barn=1.5e-3, Q_MeV=3.82, source="endf")
print(rx.getId())         # 102
print(rx.sigma_barn)      # 0.0015
```

---

#### `Isotope`
```python
class Isotope(Element):
    __slots__ = ("A", "meta", "half_life_s", "decay_constant", "mass",
                 "_decays", "_reactions")

    def __init__(self, Z, X, A, meta=None)
    # Raises ValueError if A <= 0 or Z/X mismatch

    @classmethod
    def from_symbol(cls, X: str, A: int, meta=None) -> "Isotope"
    # Isotope.from_symbol("Pb", 210)  →  Isotope(82, "Pb", 210)

    def add_decay(self, link: DecayLink) -> None
    def add_reaction(self, link: ReactionLink) -> None
    def getListOfDecays() -> List[DecayLink]
    def getListOfReactions() -> List[ReactionLink]

    @property
    def name(self) -> str    # "Pb-210" or "Pb-207m1"
```

**Manual isotope (no database):**
```python
from core.entities import Isotope, DecayLink, ReactionLink
import math

pb210 = Isotope(82, "Pb", 210)
pb210.half_life_s    = 7.0e8                         # 22.3 years
pb210.decay_constant = math.log(2) / 7.0e8

pb210.add_decay(DecayLink("beta-", 0.9999981, 0.0635, dz=1, da=0))
pb210.add_decay(DecayLink("alpha", 0.0000019, 3.792,  dz=-2, da=-4))

pb210.add_reaction(ReactionLink(102, "(n,g)", 1.5e-3, Q_MeV=3.82))

print(pb210.name)                        # Pb-210
print(len(pb210.getListOfDecays()))      # 2
print(pb210.getListOfReactions()[0].sigma_barn)  # 0.0015
```

---

### 3.2 registry.py — IsotopeBuilder

**File:** `core/registry.py`
The **only** place where `core/` entities touch `nuclear_data/`.

```python
IsotopeBuilder(provider: NuclearDataProvider)
```

#### Methods

```python
IsotopeBuilder.build_isotope(
    Z, X, A, meta=None,
    E_eV: float = 0.0253,        # energy for σ(E) lookup
    mt_filter: Optional[List[int]] = None,  # None → [102, 16, 18, 103, 107]
    audit: Optional[dict] = None,           # written with provenance info
    prefer_macs: bool = False,              # True → MACS first (astrophysics)
) -> Isotope
```
Builds one isotope with decay and reaction data from provider.

```python
IsotopeBuilder.build_range(
    specs: List[Tuple],          # [(Z, X, A)] or [(Z, X, A, meta)]
    E_eV: float = 0.0253,
    mt_filter: Optional[List[int]] = None,
    prefer_macs: bool = False,
) -> List[Isotope]
```
Builds all isotopes and calls `link_products()` automatically.

```python
IsotopeBuilder.link_products(isotopes: List[Isotope]) -> None  # static
```
Resolves `DecayLink.product` and `ReactionLink.product` pointers by matching
`(Z+dZ, A+dA, lfs)` in the isotope list.

```python
IsotopeBuilder.missing_isotopes() -> List[str]
IsotopeBuilder.print_report() -> None
```

**Example with real database:**
```python
from core.registry import IsotopeBuilder
from nuclear_data import NuclearDataProvider, EndfPathHolder

EndfPathHolder.init_xsdir("xsdir")
nd = NuclearDataProvider(
    dec_dir  = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0"),
    neu_dir  = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0"),
    macs_file = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7"),
    macs_lib  = "Endfb7",
)
builder = IsotopeBuilder(nd)

specs = [(82,"Pb",209), (83,"Bi",209), (83,"Bi",210), (84,"Po",210)]
isotopes = builder.build_range(specs, E_eV=30_000.0, mt_filter=[102],
                               prefer_macs=True)
builder.print_report()
print(builder.missing_isotopes())
```

**Example without database (fully manual):**
```python
from core.entities import Isotope, DecayLink, ReactionLink
from core.registry import IsotopeBuilder
import math

pb208 = Isotope(82, "Pb", 208); pb208.decay_constant = 0.0
pb209 = Isotope(82, "Pb", 209)
lam = math.log(2) / (3.253 * 3600)
pb209.decay_constant = lam
pb209.half_life_s    = 3.253 * 3600
bi209 = Isotope(83, "Bi", 209); bi209.decay_constant = 0.0

pb208.add_reaction(ReactionLink(102, "(n,g)", 0.5e-3))
pb209.add_decay(DecayLink("beta-", 1.0, 0.635, dz=1, da=0))

IsotopeBuilder.link_products([pb208, pb209, bi209])

# Check links
print(pb208.getListOfReactions()[0].product)  # Isotope(Pb-209, ...)
print(pb209.getListOfDecays()[0].product)     # Isotope(Bi-209, ...)
```

---

### 3.3 burnup.py — BurnupConfig and BurnupMatrix

**File:** `core/burnup.py`

#### `BurnupConfig`

```python
BurnupConfig(flux: float = 0.0)
```

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `flux` | `float` | `0.0` | Neutron flux [n/cm²/s] |
| `is_filter_active` | `bool` | `False` | True when MT filter is set |
| `allowed_mts` | `Set[int]` | `{}` | Allowed MT reactions |
| `min_half_life_s` | `float` | `0.0` | Decay half-life lower limit |
| `max_half_life_s` | `Optional[float]` | `None` | Upper limit (None = no limit) |
| `include_stable` | `bool` | `True` | Include stable isotopes |
| `open_boundary` | `bool` | `False` | Global boundary condition |
| `open_boundary_mts` | `Optional[Set[int]]` | `None` | Per-MT boundary override |

**Half-life presets** (`set_decay_time_filter(preset)`):

| Preset | Min T½ | Max T½ |
|--------|--------|--------|
| `"all"` | 0 | None |
| `"yearly"` | 0 | 3.156×10⁷ s |
| `"daily"` | 0 | 86400 s |
| `"hourly"` | 0 | 3600 s |
| `"minute"` | 0 | 60 s |
| `"second"` | 0 | 1 s |

**Methods:**
```python
BurnupConfig.enable_all_reactions() -> None
BurnupConfig.enable_specific_reactions(mt_list: List[int]) -> None
BurnupConfig.is_reaction_allowed(mt: int) -> bool

BurnupConfig.set_decay_time_filter(preset: str) -> None   # raises ValueError for unknown preset
BurnupConfig.set_decay_time_range(min_s=0.0, max_s=None) -> None
BurnupConfig.is_decay_allowed(half_life_s: Optional[float]) -> bool

BurnupConfig.set_open_boundary(open_b: bool) -> None
BurnupConfig.set_open_boundary_mts(mts: Optional[Iterable[int]]) -> None
BurnupConfig.is_mt_boundary_open(mt: int) -> bool
```

**Usage example:**
```python
from core.burnup import BurnupConfig

cfg = BurnupConfig(flux=1e14)
cfg.enable_specific_reactions([102])          # only (n,γ)
cfg.set_decay_time_filter("all")              # all half-lives
cfg.set_open_boundary(True)                   # reactions to external isotopes = sinks

# Per-MT boundary: (n,g) is open, (n,2n) is closed
cfg2 = BurnupConfig(flux=1e14)
cfg2.enable_specific_reactions([102, 16])
cfg2.set_open_boundary_mts([102])

print(cfg2.is_mt_boundary_open(102))   # True
print(cfg2.is_mt_boundary_open(16))    # False
```

---

#### `BurnupMatrix`

```python
BurnupMatrix(isotopes: List[Isotope], config: BurnupConfig)
```

**Attributes after `build()`:**
| Attribute | Shape | Description |
|-----------|-------|-------------|
| `matrix_A` | (N,N) | Full burnup matrix |
| `matrix_decay` | (N,N) | Decay contribution |
| `matrix_rxn` | (N,N) | Neutron reaction contribution |
| `matrix_heat` | (N,N) | Diagonal: heat coefficients [MeV·s⁻¹] |
| `missing` | `List[dict]` | `{"isotope": name, "mt": mt}` |

**Methods:**

```python
BurnupMatrix.build(missing_callback=None) -> np.ndarray
# Assembles matrix_A = matrix_decay + matrix_rxn.
# missing_callback({"isotope": str, "mt": int}) -> float | None:
#   called when sigma is absent — can supply a fallback value.

BurnupMatrix.solve(N0: array, t_sec: float, method="cram16") -> np.ndarray
# Solvers: "cram16" | "cram48" | "pade" | "taylor" | "bdf"
# Returns N(t).

BurnupMatrix.heat_release(N: array) -> float
# Total heat [MeV/s] = sum_j matrix_heat[j,j] * N[j]

BurnupMatrix.solve_equilibrium(N0, method="cram16", rtol=1e-4,
                               t_start=None, t_max=1e25, max_steps=60)
    -> Tuple[np.ndarray, float]
# Returns (N_eq, t_reached). Adaptive doubling until dN/N < rtol.

BurnupMatrix.find_plateau(N0, method="cram16", n_scan=80,
                          t_min=1e-3, t_max=1e25)
    -> Tuple[float, float]
# Returns (t_plateau, Q_plateau). Scans Q(t) on log-time grid.

BurnupMatrix.equilibrium_heat(N0, method="cram16", t_plateau=None,
                              n_scan=80, t_min=1e-3, t_max=1e25) -> float
# Closed chain → solve_equilibrium; open chain → find_plateau.

BurnupMatrix.evolve(N0, times: Iterable[float], method="cram16")
    -> Tuple[List[float], List[np.ndarray]]
# Returns (times, trajectories). trajectories[k] = N at times[k].

BurnupMatrix.heat_curve(trajectories: List[np.ndarray]) -> List[float]
# [heat_release(N) for N in trajectories]

BurnupMatrix.names() -> List[str]
# [iso.name for iso in self.isotopes]

BurnupMatrix.validate() -> bool
# Checks: diagonal ≤ 0 and column sums ≤ 0 (conservation).

BurnupMatrix.export_txt(filename: str, mode="total") -> None
# mode: "total" | "rxn" | "decay"

BurnupMatrix.export_evolution_txt(filename, times, trajectories) -> None
# Writes N(t): Time_sec | isotope columns
```

**Standalone example:**
```python
import numpy as np
import math
from core.entities import Isotope, DecayLink
from core.burnup import BurnupConfig, BurnupMatrix

# A --[beta-]--> B (stable)
lam = 1e-8
iso_A = Isotope(1, "H", 3)
iso_A.decay_constant = lam
iso_A.half_life_s    = math.log(2) / lam
iso_B = Isotope(2, "He", 3)

dl = DecayLink("beta-", 1.0, 0.0186, dz=1, da=0, product=iso_B)
iso_A.add_decay(dl)

cfg = BurnupConfig(flux=0.0)
cfg.set_decay_time_filter("all")

mtx = BurnupMatrix([iso_A, iso_B], cfg)
mtx.build()

N0 = np.array([1.0, 0.0])
t  = 1e9   # seconds

N = mtx.solve(N0, t, method="cram16")
print(f"H-3:  {N[0]:.4f}")
print(f"He-3: {N[1]:.4f}")

q = mtx.heat_release(N)
print(f"Heat: {q:.3e} MeV/s")

N_eq, t_eq = mtx.solve_equilibrium(N0)
print(f"Equilibrium: H-3={N_eq[0]:.2e}, He-3={N_eq[1]:.4f} at t={t_eq:.1e}s")

times = np.logspace(6, 12, 30)
ts, traj = mtx.evolve(N0, times)
qc = mtx.heat_curve(traj)

mtx.export_txt("matrix.txt")
mtx.export_evolution_txt("evolution.txt", ts, traj)
```

#### `equilibrium_heat_vs_flux` (module-level function)

```python
from core.burnup import equilibrium_heat_vs_flux

fluxes = np.logspace(14, 28, 8)
fl, heats = equilibrium_heat_vs_flux(isotopes, cfg, N0, fluxes)
# Returns (fluxes_array, heats_array)
```

---

### 3.4 cycle_finder.py — CycleAnalyzer

**File:** `core/cycle_finder.py`

#### Data structures

```python
@dataclass
class CycleStep:
    isotope: Isotope
    label:   str           # "beta-", "(n,g) s=0.5b", ...
    kind:    str           # "decay" | "reaction"
    mt:      Optional[int] # None for decays

@dataclass
class IsotopeCycle:
    steps: List[CycleStep]

    @property length     -> int
    @property isotope_names -> List[str]
    @property is_pure_decay  -> bool
    @property is_pure_reaction -> bool
    @property is_mixed   -> bool
    @property n_decays   -> int
    @property n_reactions -> int
    def __len__()  -> int
    def __str__()  -> str  # "Pb-208 --[alpha]--> Po-212 --[...]-->"

@dataclass
class SCC:
    isotopes: List[Isotope]
    cycles:   List[IsotopeCycle]    # filled after find_all_cycles()

    @property names -> List[str]   # sorted isotope names
    def __len__() -> int
    def __repr__() -> str          # "SCC(['Bi-210', 'Po-210'])"
```

#### `CycleAnalyzer`

```python
CycleAnalyzer(isotopes: List[Isotope])
# isotopes must have .product references set (after link_products)
```

**Methods:**
```python
CycleAnalyzer.find_sccs() -> List[SCC]
# Tarjan's algorithm O(V+E). Returns only non-trivial SCCs (≥2 vertices).

CycleAnalyzer.find_all_cycles() -> List[IsotopeCycle]
# DFS within each SCC + canonical deduplication.
# Sorted by cycle length.

CycleAnalyzer.has_cycles() -> bool
# Quick check: bool(find_sccs())

CycleAnalyzer.isotopes_in_cycles() -> List[Isotope]
# All isotopes belonging to at least one SCC.

CycleAnalyzer.cycles_through(isotope: Isotope) -> List[IsotopeCycle]
# All cycles containing the given isotope.

CycleAnalyzer.print_report(cycles=None, max_show=30) -> None
# Full report: SCC groups, cycle categories, participating isotopes.
```

**Standalone example:**
```python
from core.entities import Isotope, DecayLink, ReactionLink
from core.cycle_finder import CycleAnalyzer

pb208 = Isotope(82, "Pb", 208); pb208.decay_constant = 0.0
bi209 = Isotope(83, "Bi", 209); bi209.decay_constant = 1e-9
po210 = Isotope(84, "Po", 210); po210.decay_constant = 5.8e-8

# Bi-209 (n,g) -> Po-210 -> [alpha] -> Pb-206? ... construct cycle manually
pb208.add_reaction(ReactionLink(102, "(n,g)", 1e-3, product=bi209))
bi209.add_decay(DecayLink("beta-", 1.0, 0.5, dz=-1, da=0, product=pb208))

analyzer = CycleAnalyzer([pb208, bi209, po210])
cycles = analyzer.find_all_cycles()
analyzer.print_report()

for c in cycles:
    print(f"Cycle L={c.length}: {c.isotope_names}, mixed={c.is_mixed}")

for scc in analyzer.find_sccs():
    print(f"SCC: {scc.names}")
```

---

### 3.5 sensitivity.py — SensitivityEngine

**File:** `core/sensitivity.py`

#### Data structures

```python
@dataclass
class SensitivityConfig:
    initial_by_label: Dict[str, float]   # {"Pb-208": 1.0}
    target_kind:      str     # "concentration" | "heat"
    target_label:     Optional[str]      # isotope name
    eval_mode:        str     # "time" | "plateau"
    t_eval_s:         float   # 1e9 [s]
    types:            List[str]          # ["sigma", "halflife", "flux_energy", "library"]
    delta:            float   # 0.01 (1%)
    two_sided:        bool    # True → central difference
    scope_channels:   Optional[List[Tuple[str,int]]]  # [(name, MT)]
    scope_isotopes:   Optional[List[str]]
    flux:             float
    energy_eV:        float   # 30000.0
    mt_filter:        Optional[List[int]]
    method:           str     # "cram16"
    hl_preset:        str     # "all"
    open_boundary_mts: Optional[List[int]]
    prefer_macs:      bool
    libraries:        List[str]
    reference_library: Optional[str]

@dataclass
class SensRow:
    parameter: str     # "Pb-210 (n,γ)" | "Pb-210 λ" | "Flux Φ" | "Energy E"
    type:      str     # "σ" | "λ" | "Φ" | "E"
    base:      float   # base target value
    perturbed: float   # value at +δ
    S:         float   # sensitivity coefficient
    rank:      int     # filled by ranked()

@dataclass
class LibRow:
    library:       str
    N_i:           Optional[float]
    Q:             float
    deviation_pct: float

@dataclass
class SensitivityResult:
    target_desc: str
    base_value:  float
    rows:        List[SensRow]     # OAT rows
    library_rows: List[LibRow]
    lib_stats:   Dict[str, float]  # min/max/mean/std

    def ranked() -> List[SensRow]  # sorted by |S| descending
```

#### `SensitivityEngine`

```python
SensitivityEngine(service: NuclearService)

SensitivityEngine.run(
    cfg: SensitivityConfig,
    progress: Optional[Callable[[float, str], None]] = None,
    should_cancel: Optional[Callable[[], bool]] = None,
) -> SensitivityResult
```

**OAT formula (two-sided):**
```
S_j = [(N(p_j + δ) - N(p_j - δ)) / N_base] / (2·δ)
```

---

### 3.6 plotting.py — plot helpers

**File:** `core/plotting.py`
Requires `pip install matplotlib`. All functions save to disk (Agg backend).

```python
plot_concentrations(
    names: List[str],
    times: List[float],
    trajectories: List[np.ndarray],
    filename="concentrations.png",
    mode="relative",           # "relative" | "absolute" | "sigma_n"
    sigma_weights=None,        # List[float] for σ·N mode
    stable_isotopes=None,
    title="Isotope Concentration Evolution",
    logy=False, logx=False,
) -> str                       # saved filename

plot_heat(
    times, heat_curve,
    filename="heat.png",
    title="Decay Heat",
    unit="MeV/s",              # "MeV/s" | "W/cm3"
    density=1.0, M=1.0,        # for W/cm3 conversion
) -> str

plot_heat_vs_flux(
    fluxes, heats,
    filename="heat_flux.png",
    title="Equilibrium Heat vs Flux",
) -> str
```

**Example:**
```python
from core.plotting import plot_concentrations, plot_heat

names = mtx.names()
times = list(np.logspace(4, 12, 50))
ts, traj = mtx.evolve(N0, times)

plot_concentrations(names, ts, traj, "conc.png", mode="relative")
qc = mtx.heat_curve(traj)
plot_heat(ts, qc, "heat.png", unit="MeV/s")
```

---

### 3.7 service.py — NuclearService

**File:** `core/service.py`
Qt-free application singleton. The GUI never calls `core/` directly — only through `NuclearService`.

#### Module-level utilities (no GUI needed)

```python
from core.service import mt_symbol, decay_symbol, CHANNEL_DELTAS

mt_symbol(mt: int) -> str
# 102 → "(n,γ)", 16 → "(n,2n)", 107 → "(n,α)", 999 → "MT=999"

decay_symbol(mode: Optional[str]) -> str
# "beta-" → "β⁻", "alpha" → "α", "it" → "IT", "ec/beta+" → "β⁺/EC"
# Comma-separated: "beta-, alpha" → "β⁻, α"

CHANNEL_DELTAS: Dict[str, Tuple[int,int]]
# Net (dZ, dA) for chart reachability:
# "beta-":(+1,0), "alpha":(-2,-4), "(n,γ)":(0,+1), "(n,2n)":(0,-1), etc.
```

#### Key regex patterns (module-level)
```python
_LABEL_RE   # parses "Pb-210", "Pb210", "U-235m1" → (X, A, meta)
_RANGE_RE   # parses "Tl206..Tl216" → element + mass range
_DECAY_RE   # parses ENDF filename "dec-082_Pb_210.endf"
```

---

## 4. GUI layer (model/)

| File | Panel | Purpose |
|------|-------|---------|
| `model_CoreNucleo.py` | Periodic Table | Element/isotope browser, data source diagnostics |
| `model_IsotopeChart.py` | Isotope Chart | Z vs N map; isotope selection; arrow filtering; Build Matrix |
| `model_results.py` | Results | 4 tabs: matrix view, evolution, equilibrium, sensitivity |
| `model_graph.py` | Graph | 6 plot types: concentrations, heat, abundance, eq-flux, sens. bar |
| `GUIcomp.py` | Shared | Reusable Qt widgets used by multiple panels |
| `core.py` | Launcher | `QMainWindow` entry point; loads panels dynamically |

GUI panels call only `NuclearService` methods — they never import `BurnupMatrix` or `NuclearDataProvider` directly.

---

## 5. Interaction map

```
User → GUI panel
         │
         └──► NuclearService.build_universe(z_from, z_to, n_from, n_to, E_eV, ...)
                  │
                  ├──► NuclearDataProvider.get_decay(Z, X, A)
                  │         → DecayData
                  ├──► NuclearDataProvider.get_macs_with_source(Z,X,A, mt, kT)
                  │         → (MacsPoint, source_name)
                  └──► IsotopeBuilder.build_range(specs, E_eV, mt_filter)
                            → List[Isotope]  (with product links)

         └──► NuclearService.build_matrix(flux, energy_eV, mt_filter, ...)
                  │
                  └──► BurnupMatrix(isotopes, cfg).build()
                            → matrix_A  (N×N numpy array)

         └──► NuclearService.run_evolution(N0, t_end, method, ...)
                  │
                  └──► BurnupMatrix.evolve(N0, times, method)
                            → (times, trajectories)

         └──► NuclearService.find_cycles()
                  │
                  └──► CycleAnalyzer(isotopes).find_all_cycles()
                            → List[IsotopeCycle]

         └──► NuclearService.run_sensitivity(cfg)
                  │
                  └──► SensitivityEngine(service).run(cfg)
                            → SensitivityResult
```

---

## 6. Full call examples

### 6.1 Provider → Builder → Matrix → Solve (with real xsdir)

```python
import sys, os, numpy as np
sys.path.insert(0, ".")

from nuclear_data import EndfPathHolder, NuclearDataProvider
from core.registry import IsotopeBuilder
from core.burnup import BurnupConfig, BurnupMatrix
from core.cycle_finder import CycleAnalyzer
from core.plotting import plot_concentrations, plot_heat

# 1. Paths
EndfPathHolder.init_xsdir("xsdir")
nd = NuclearDataProvider(
    dec_dir          = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0"),
    neu_dir          = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0"),
    macs_file        = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7"),
    macs_lib         = "Endfb7",
    user_dir         = EndfPathHolder.set_USER_DB_DIR("user_db"),
    endfb71_macs_dir = str(EndfPathHolder.get_ENDFB71_MACS_DIR()),
    eaf2010_macs_dir = str(EndfPathHolder.get_EAF2010_MACS_DIR()),
    talys_dir        = str(EndfPathHolder.get_TALYS_DIR()),
    macs_priority    = ["user","endfb71","eaf2010","rawmacs","talys","endf_computed"],
)

# 2. Build isotopes (s-process termination network, kT=30 keV)
specs = [
    (81,"Tl",206),(81,"Tl",207),
    (82,"Pb",206),(82,"Pb",207),(82,"Pb",208),
    (82,"Pb",209),(82,"Pb",210),(82,"Pb",211),
    (83,"Bi",209),(83,"Bi",210),(83,"Bi",211),
    (84,"Po",210),(84,"Po",211),
]
builder  = IsotopeBuilder(nd)
isotopes = builder.build_range(specs, E_eV=30_000.0, mt_filter=[102], prefer_macs=True)
builder.print_report()

# 3. Configure burnup
cfg = BurnupConfig(flux=1e14)
cfg.enable_specific_reactions([102])
cfg.set_decay_time_filter("all")
cfg.set_open_boundary(True)

# 4. Initial condition: all in Pb-206
names = [i.name for i in isotopes]
N0    = np.zeros(len(isotopes))
N0[names.index("Pb-206")] = 1.0

# 5. Build matrix and solve
mtx = BurnupMatrix(isotopes, cfg)
mtx.build()
print("Matrix valid:", mtx.validate())
print("Missing σ:", mtx.missing)

times = np.logspace(4, 18, 60)
ts, traj = mtx.evolve(N0, times)

# 6. Heat
qc = mtx.heat_curve(traj)
plot_heat(ts, qc, "heat.png")
plot_concentrations(names, ts, traj, "conc.png", mode="relative")

# 7. Equilibrium
N_eq, t_eq = mtx.solve_equilibrium(N0)
print(f"Equilibrium at t={t_eq:.1e} s")
for n, v in zip(names, N_eq):
    if v > 1e-4:
        print(f"  {n}: {v:.4f}")

# 8. Cycle detection
analyzer = CycleAnalyzer(isotopes)
cycles   = analyzer.find_all_cycles()
analyzer.print_report()
```

---

### 6.2 Fully manual — no database, no xsdir

```python
import numpy as np, math
from core.entities import Isotope, DecayLink, ReactionLink
from core.registry import IsotopeBuilder
from core.burnup import BurnupConfig, BurnupMatrix
from core.cycle_finder import CycleAnalyzer

# Build isotopes manually
def iso(Z, X, A, lam=0.0, hl=None):
    i = Isotope(Z, X, A)
    i.decay_constant = lam; i.half_life_s = hl
    return i

pb208 = iso(82, "Pb", 208)
pb209 = iso(82, "Pb", 209, lam=math.log(2)/(3.253*3600), hl=3.253*3600)
bi209 = iso(83, "Bi", 209)

pb208.add_reaction(ReactionLink(102, "(n,g)", 0.5e-3, product=pb209))
pb209.add_decay(DecayLink("beta-", 1.0, 0.635, dz=1, da=0, product=bi209))

isotopes = [pb208, pb209, bi209]

# Matrix
cfg = BurnupConfig(flux=1e14)
cfg.enable_specific_reactions([102])
cfg.set_decay_time_filter("all")
cfg.set_open_boundary(False)

mtx = BurnupMatrix(isotopes, cfg)
mtx.build()
print("A =\n", mtx.matrix_A)

N0 = np.array([1.0, 0.0, 0.0])
N  = mtx.solve(N0, 1e6, method="cram48")
print("N(1e6 s) =", N)

# Cycle check (no cycles in linear chain)
analyzer = CycleAnalyzer(isotopes)
print("Cycles:", analyzer.find_all_cycles())
```

---

### 6.3 Direct nuclear_data queries (no builder)

```python
# Query raw data directly from provider (no BurnupMatrix needed)
dd  = nd.get_decay(82, "Pb", 210)
lam = dd.decay_constant
print(f"Pb-210 λ = {lam:.3e} s⁻¹, T½ = {dd.half_life_s:.2e} s")
for ch in dd.channels:
    print(f"  {ch.mode}: BR={ch.branch:.4f}, Q={ch.Q_MeV} MeV, dZ={ch.dz}, dA={ch.da}")

pt, src = nd.get_macs_with_source(82, "Pb", 210, mt=102, kT_keV=30.0)
print(f"MACS = {pt.sigma_mb} mb  [{src}]")

sig, src2 = nd.get_sigma_with_source(82, "Pb", 210, mt=102, E=30000.0, e_unit="eV")
print(f"σ(30keV) = {sig:.3e} barn  [{src2}]")

branches = nd.get_isomeric_branches(82, "Pb", 210, mt=102, E_eV=30000.0)
for b in branches:
    print(f"  LFS={b.LFS}, yield={b.yield_}")

# Change MACS priority at runtime
nd.set_macs_priority(["user", "eaf2010", "endfb71", "rawmacs", "talys", "endf_computed"])
print(nd.get_macs_priority())
```

---

### 6.4 CycleAnalyzer — all public methods

```python
from core.cycle_finder import CycleAnalyzer

analyzer = CycleAnalyzer(isotopes)   # isotopes with product links

# Basic checks
print(analyzer.has_cycles())                    # bool

# SCC analysis
sccs = analyzer.find_sccs()                     # List[SCC]
for scc in sccs:
    print(f"SCC: {scc.names} ({len(scc)} isotopes)")
    print(f"  cycles: {len(scc.cycles)}")       # after find_all_cycles()

# Full enumeration
cycles = analyzer.find_all_cycles()             # List[IsotopeCycle]
for c in cycles:
    print(f"L={c.length}, mixed={c.is_mixed}, decays={c.n_decays}, rxn={c.n_reactions}")
    print(f"  isotopes: {c.isotope_names}")
    print(f"  {c}")                             # text path

# Targeted queries
in_cycles = analyzer.isotopes_in_cycles()
pb210 = next(i for i in isotopes if i.name == "Pb-210")
pb210_cycles = analyzer.cycles_through(pb210)

# Report
analyzer.print_report(cycles=cycles, max_show=20)
```

---

### 6.5 BurnupConfig — all filters

```python
from core.burnup import BurnupConfig

cfg = BurnupConfig(flux=1e14)

# Reaction filter
cfg.enable_specific_reactions([102, 16, 107])
print(cfg.is_reaction_allowed(102))   # True
print(cfg.is_reaction_allowed(18))    # False
cfg.enable_all_reactions()
print(cfg.is_reaction_allowed(18))    # True

# Decay filter
for preset in ["all", "yearly", "daily", "hourly", "minute", "second"]:
    cfg.set_decay_time_filter(preset)
    print(f"{preset}: 1s → {cfg.is_decay_allowed(1.0)}, 1yr → {cfg.is_decay_allowed(3.2e7)}")

cfg.set_decay_time_range(min_s=60.0, max_s=3.156e7)  # 1 min to 1 year
print(cfg.is_decay_allowed(30.0))    # False
print(cfg.is_decay_allowed(3600.0))  # True

# Boundary
cfg.set_open_boundary(True)
print(cfg.is_mt_boundary_open(102))  # True (global open)
cfg.set_open_boundary_mts([102])     # only 102 is open
print(cfg.is_mt_boundary_open(16))   # False
cfg.set_open_boundary_mts(None)      # revert to global
```

---

### 6.6 BurnupMatrix — all solvers compared

```python
import numpy as np, math, time
from core.entities import Isotope, DecayLink
from core.burnup import BurnupConfig, BurnupMatrix

lam = 1e-9
iso_A = Isotope(1, "H", 3);   iso_A.decay_constant = lam
iso_B = Isotope(2, "He", 3);  iso_B.decay_constant = 0.0
iso_A.add_decay(DecayLink("beta-", 1.0, 0.0186, 1, 0, product=iso_B))

cfg = BurnupConfig(); cfg.set_decay_time_filter("all")
mtx = BurnupMatrix([iso_A, iso_B], cfg)
mtx.build()

N0 = np.array([1.0, 0.0])
t  = 5e9   # 5 × 10⁹ s

analytic_A = math.exp(-lam * t)
print(f"Analytic: H-3={analytic_A:.6f}, He-3={1-analytic_A:.6f}")

for method in ["cram16", "cram48", "pade", "taylor", "bdf"]:
    t0 = time.perf_counter()
    N = mtx.solve(N0, t, method=method)
    dt = time.perf_counter() - t0
    err = abs(N[0] - analytic_A)
    print(f"  {method:8s}: H-3={N[0]:.6f}  err={err:.2e}  t={dt*1000:.2f}ms")
```

---

## 7. Algorithm: cycle detection

### Step 1 — Build adjacency graph `_build_adj()`

For each isotope, collect outgoing edges (decay and reaction links) where the
product is in the network:

```
for iso in isotopes:
    for d in iso.getListOfDecays():
        if d.product in iso_set:
            adj[iso].append((d.product, label, "decay", None))
    for r in iso.getListOfReactions():
        if r.product in iso_set:
            adj[iso].append((r.product, label, "reaction", r.mt))
```

### Step 2 — Tarjan's SCC `find_sccs()` — O(V+E)

Recursive DFS, tracking discovery index and low-link values:

```
visit(v):
    idx[v] = low[v] = counter++
    push(v)
    for w in adj[v]:
        if w not visited: visit(w); low[v] = min(low[v], low[w])
        elif w on stack:  low[v] = min(low[v], idx[w])
    if low[v] == idx[v]:
        pop SCC from stack until v
        if |SCC| > 1: add to result
```

### Step 3 — DFS cycle enumeration `find_all_cycles()`

For each non-trivial SCC, run DFS from each vertex looking for back-edges to the start:

```
dfs(v, steps, start):
    for w in sub[v]:
        if w == start and steps != []:
            record cycle(steps + [step])
        elif w not in path:
            dfs(w, steps + [step], start)
```

### Step 4 — Canonical deduplication `_canon_key()`

Rotate cycle so the lexicographically smallest isotope name comes first.
Two cycles with different starting vertices but identical sequence → same key:

```
key = pairs[min_index:] + pairs[:min_index]   # where pairs = [(name, label)]
```

---

## 8. Algorithm: sensitivity analysis

### OAT formula

Two-sided central difference (default, δ = 1%):
```
S_j = [(N_target(p_j + δ·p_j) - N_target(p_j - δ·p_j)) / N_base] / (2·δ)
```

One-sided:
```
S_j = [(N_target(p_j + δ·p_j) - N_base) / N_base] / δ
```

### Four analysis types

| Type | Perturbed parameter | Rebuild needed |
|------|---------------------|----------------|
| `sigma` | `rx.sigma_barn *= (1 ± δ)` on live `ReactionLink` | Matrix only |
| `halflife` | `iso.decay_constant *= (1 ± δ)` on live `Isotope` | Matrix only |
| `flux_energy` | Flux: new `cfg.flux`; Energy: `_rebuild_universe(E_new)` | Matrix / Universe |
| `library` | `switch_decay_library()` or `set_macs_priority()` | Full universe |

### Minimal subnetwork optimization

For concentration targets, only isotopes on paths source→target are perturbed:
```python
subset = service.minimal_subnetwork(initial_labels, [target_label], mt_filter)
```
Reduces OAT iterations proportionally to network size.

### Progress and cancellation

Each `_tick(msg)` call increments the done counter and checks `should_cancel()`:
```python
def _tick(self, msg):
    self._done += 1
    if self._cancel and self._cancel():
        raise RuntimeError("Sensitivity analysis cancelled.")
    if self._progress:
        self._progress(self._done / self._total, f"{msg} ({self._done}/{self._total})")
```

---

*References: Tarjan (1972), Johnson (1975), Pusa (2010, 2016 — CRAM coefficients), KADoNiS/BRUSLIB (MACS tables).*


## Exposure averaging (fixed-condition Evolution)

The Exposure averaging tab processes a stored Evolution trajectory. It does not
rebuild the burnup matrix or rerun CRAM when the mean exposure changes.

### GUI workflow

1. Run Evolution with fixed flux and energy/kT, or open an existing .kaz session.
2. In Results, select Exposure averaging and enter Mean exposure tau0 [mbarn^-1].
3. Click Average results. Review coverage, the uncovered tail, the numerical
   weight error, the inventory range, and negative-concentration diagnostics.
4. Inspect the Isotopes and Mass numbers tables. Plot averaged curves in the Graph module.
5. In the main Graph model choose Abundance, then Exposure-averaged. The time
   slider is disabled for this source; Absolute, Relative and Sigma N remain available.
6. Export averaged results writes both tables as TXT or CSV. Save averaged session
   includes the result and configuration in a .kaz file. Saved-session viewers
   expose the same averaging tab; old sessions remain readable.

### Definitions and finite range

For flux Phi in n/cm^2/s and time in seconds, tau = 1e-27 * Phi * t
in mbarn^-1. The exponential weight is exp(-tau/tau0)/tau0. This rho(tau)
is a probability density of exposures, not material mass density.

The module integrates the full computed isotope vector on the original exposure
grid using trapezoidal quadrature. It ignores display targets during integration.
A restricted source-to-target run still represents a restricted physical network;
its warning cannot be resolved by display settings.

The default output is the finite-range integral, without renormalization.
Analytical coverage is 1-exp(-tau_max/tau0); the tail is exp(-tau_max/tau0).
Tail > 0.001 produces a warning. The numerical integral of the weight is also
compared with analytical coverage; relative error > 0.005 warns about grid
resolution. A small tail alone does not certify convergence. Refine the original
saved time grid and compare averaged inventories. Extending the duration requires
a new Evolution run; changing tau0 alone does not.

Relative mode divides each averaged inventory by the sum of all averaged isotope
inventories, before any presentation filtering. It describes the composition of
the retained inventory, not recovered material, and must be read with the balance
diagnostics.

### Cross sections and aggregation

Registry ReactionLink.sigma_barn already includes the isomeric yield:
sigma_barn = total_capture_sigma * branch_yield. capture_snapshot sums all MT=102
links once; it must not multiply q_yield again. This convention also applies to
manually constructed links supplied to the exposure API.

New Evolution results store a numerical snapshot of these capture cross sections
and source provenance. Subsequent universe edits cannot alter the averaged Sigma N.
Legacy runs use their saved isotope objects and report limited provenance.
No capture channel means unavailable (NaN internally, NA on export), not a measured
zero cross section. A mass-number group with an unavailable contribution remains NA.

Absolute mass abundance is sum_i mean_N_i for equal A.
The mass Sigma N is sum_i sigma_i * mean_N_i, not a single sigma(A) multiplied by
sum_i mean_N_i. This aggregate is not automatically the classical stable s-only
isotope curve. Compare matching isotopes, units, normalization and physical histories.

### Diagnostics and limits

All times, concentrations, flux and tau0 must be finite. Times start at zero and
increase strictly. Flux and tau0 are positive. Variable-condition histories are
rejected. For legacy RunResult files, fixed conditions are inferred from the
fixed-matrix Evolution API and explicitly marked as legacy provenance.

The negative tolerance is negative_atol + negative_rtol * initial_inventory_scale,
with defaults 1e-14 and 1e-10 in input inventory units. More negative values block
averaging. Negative trajectory totals also block it. Smaller signed values are
retained and reported, not clipped. There is no automatic correction of boundary
losses, matrix defects, or poorly conditioned solver outputs.

Closed capture/one-heavy-daughter decay networks conserve the heavy-nucleus count;
open networks or other reaction stoichiometries need not. New run metadata records
columns with net loss, filters, manual channels, library names and capture sources.
These are diagnostics, not proof that every channel is physically complete.
An existing numerical issue remains outside this feature: the legacy method name cram48 delegates to adaptive CRAM-16, not genuine CRAM-48. It is not an independent
validation solver.

### API and data flow

core.exposure:
- ExposureConfig: tau0 and diagnostic tolerances.
- ExposureResult: full mean vector, capture snapshot, metadata, warnings and diagnostics.
- capture_snapshot(names, isotopes): sum branch-weighted MT=102 data and source names.
- average_exposure(...): validated finite-range trapezoidal integration.
- average_run(run, config): adapter reading only the original run.
- isotope_values(result, mode): absolute, relative or sigma_n vector.
- mass_values(result, mode): shared aggregation for tables, plots and export.
- export_exposure(result, path): TXT/CSV with provenance and both tables.

RunResult -> average_run -> ExposureResult -> Results tables / Graph / export.
RunResult.exposure_result is serialized in session format v3. The load_session
three-value return contract is unchanged. Legacy v1/v2 files may omit new fields.

Example:
```python
from core.exposure import ExposureConfig, average_run, mass_values, export_exposure

averaged = average_run(run, ExposureConfig(tau0=0.3))
run.exposure_result = averaged
mass_numbers, abundance = mass_values(averaged, "absolute")
_, capture_flow = mass_values(averaged, "sigma_n")
export_exposure(averaged, "exposure.csv")
```

The GUI adapter is model/exposure_panel.py. ExposurePanel is reused in live Results
and saved-session viewers. Calculations use the existing worker/progress mechanism;
model/modelgraph/exposure_plot.py provides plot_exposure for the main Graph model and the saved-session Graph tab. The integration works on private
array copies and never edits the original trajectory.

Tests: python -B tests/test_exposure.py


## Adaptive CRAM-16 correction

Evolution now uses sparse CRAM-16 with step doubling. Each requested output time
is still solved from the original inventory, but the solver refines internal
substeps until componentwise agreement (rtol=1e-9, absolute scale 1e-14 times
initial absolute inventory), signed-value checks and applicable inventory checks
pass. The stored time grid therefore remains a presentation grid; more plotted
points alone were not a remedy for inaccurate single-step CRAM.

No negative clipping or inventory renormalization is applied. Refinement failure
raises an explicit error; cancellation is checked between sparse solves.
The legacy cram48 name now delegates to the same adaptive CRAM-16 implementation.
It is not genuine CRAM-48 and is not an independent verification method.

Old .kaz trajectories are not recalculated when opened. Rerun Evolution after
updating the application, then use Average results. Invalid old trajectories now
report the isotope, time and concentration instead of displaying only a traceback.
Boundary losses remain part of the configured network and are not changed by this
numerical correction. Regression tests: tests/test_cram_adaptive.py.
