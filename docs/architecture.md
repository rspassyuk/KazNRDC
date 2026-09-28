# Architecture and source map

## Data flow

```text
ENDF / MACS / TALYS / user JSON and CSV
                  |
       NuclearDataProvider
                  |
          IsotopeBuilder
                  |
 Isotope / DecayLink / ReactionLink
                  |
    NuclearService -> BurnupMatrix -> adaptive CRAM / SciPy
          |                  |
          |            N(t), Q(t)
          |                  |
          +---- RunResult / EqResult / SensitivityResult
                             |
               ExposureResult (optional average)
                             |
                 GUI tables / graphs / exports
```

## Nuclear data

| Module | Responsibility |
|---|---|
| `nuclear_data/models.py` | Dataclasses for decay channels, cross sections, MACS and isomeric branching |
| `path_holder.py` | `EndfPathHolder`: database locations and discovery |
| `endf_reader.py` | ENDF-6 decay and neutron data, including relevant MF=8/MT=457 and MF=3/9/10 sections |
| `macs_reader.py`, `multi_macs_reader.py` | Single-temperature and multi-temperature MACS tables and caches |
| `talys_reader.py` | TALYS cross-section and recommended MACS input |
| `user_reader.py` | User JSON and CSV overrides |
| `provider.py` | Source priorities, fallback and unified query interfaces |

The default decay priority is user data, then ENDF-6. Default MACS priority is user, rawmacs, ENDF/B-VII.1, EAF-2010, TALYS, then ENDF-computed MACS. GUI configuration can change the order. Pointwise and MACS lookup are distinct operations.

## Calculation core

| Module | Main objects and role |
|---|---|
| `core/entities.py` | `ChemTable`, `Element`, `Isotope`, `DecayLink`, `ReactionLink`; linked domain objects |
| `registry.py` | `IsotopeBuilder`; obtains evaluated rates and links daughters |
| `burnup.py` | `BurnupConfig`, `BurnupMatrix`; filters, matrix and heat assembly, evolution and equilibrium |
| `cram.py` | Adaptive sparse CRAM-16 action with numerical checks |
| `service.py` | `NuclearService`, result dataclasses and GUI-facing orchestration |
| `sensitivity.py` | `SensitivityConfig`, `SensitivityEngine`, `SensRow`, `LibRow`, `SensitivityResult` |
| `cycle_finder.py` | `CycleAnalyzer`, `CycleStep`, `IsotopeCycle`, `SCC`; graph analysis |
| `exposure.py` | `ExposureConfig`, `ExposureResult`; physical-time/exposure quadrature, aggregation, export |
| `session.py` | Result-session serialization and loading |
| `plotting.py` | Matplotlib helpers for scripts |

`core/exposure.py` retains its module/class names to keep saved sessions and API callers compatible. Its public averaging workflow now supports physical-time and neutron-exposure distributions.

## GUI

`model/core.py` creates the application, loads themes and assembles the panels. `model/GUIcomp.py` supplies common widgets and worker/progress helpers.

- `model_CoreNucleo`: periodic table and nuclear-data inspection.
- `model_IsotopeChart`: isotope selection, displayed arrows, data configuration and network construction.
- `model_results`: network tables, manual channels, initial conditions, evolution, equilibrium, sensitivity and saved-result windows.
- `model_graph`: graph selection, units, saved-time selection and abundance source.
- `model/exposure_panel.py`: averaging configuration and tables; the graph remains in the graph module.
- `model/modelgraph/exposure_plot.py`: shared averaged-abundance renderer for live and saved results.

Graph display controls never rerun the solver. The averaging engine reads a saved run and its capture-cross-section snapshot. The model's selected energy or library may have changed since the run; those current controls must not replace the original provenance.

## Algorithms

Cycle detection builds an adjacency graph from linked transitions, finds strongly connected components with Tarjan's algorithm, enumerates simple cycles by depth-first search and removes equivalent rotations. This implementation is not described as a full Johnson implementation. The GUI graph analysis can include stored links that are disabled in the current matrix.

OAT sensitivity solves a baseline target, perturbs one parameter at a time, rebuilds/solves and ranks absolute dimensionless coefficients. A concentration target may use a source-to-target subnetwork; a heat target uses the wider network. Library comparisons change source selection and report deviations, rather than derivatives.

## Tests

- `test_all.py`: domain objects, assembly, solver interfaces, synthetic chains and core utilities.
- `test_workspace_membership.py`: selected-network membership and manual-channel behavior.
- `test_cram_adaptive.py`: solver accuracy and stiff-network regressions.
- `test_exposure.py`: exposure quadrature, aggregation, exports and GUI integration.
- `test_time_averaging.py`: temporal averages, endpoint interpolation and session compatibility.
- `test_graph_theme.py`: mode-button theme and selection behavior.

These are focused regression suites, not an assertion that every public/private method or every nuclear dataset is validated.
