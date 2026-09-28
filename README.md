# KazNRDC

**Kazakhstan Nuclear Reaction Data Code** is a Python/PyQt5 desktop application for nuclear-data inspection, nuclide-network construction and fixed-condition transmutation calculations. It supports reactor-physics and nucleosynthesis studies using user-supplied evaluated data.

## Start here

- [User guide (English)](docs/user-guide.md)
- [Calculation methods and limitations](docs/methods.md)
- [Architecture and source map](docs/architecture.md)
- [GitHub and Read the Docs publication](docs/publication.md)

## Installation

Python 3.11 or newer is required. Use a virtual environment. The application uses Qt 5; Linux systems also need the system libraries required by Qt's platform plugin.

```bash
git clone https://github.com/rspassyuk/KazNRDC.git
cd KazNRDC
python -m venv .venv
```

On Windows (PowerShell):

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-endf.txt
.\.venv\Scripts\python.exe model/core.py
```

On Linux/macOS:

```bash
.venv/bin/python -m pip install -r requirements-endf.txt
.venv/bin/python model/core.py
```

Prepare `xsdir/` before starting calculations; see the data layout in the user guide. The nuclear libraries are external datasets and are not included in the source repository. `requirements.txt` installs the numerical and GUI dependencies; `requirements-endf.txt` also installs the ENDF reader. Without `endf`, only data supported by other configured readers are available. The GUI still expects the configured directory structure.

## Workflow

1. Inspect decay and reaction data in **Periodic Table**.
2. Select isotopes, nuclear-data sources, energy and transition types in **Isotope Chart**.
3. Use **Build Burnup Matrix**, then inspect **Results → Isotopes & Matrix**. Add/remove isotopes or define manual decay channels as needed.
4. Set initial inventories, flux or neutron density, and irradiation duration. Run **Evolution**.
5. Inspect individual isotopes and mass-number sums with the **Time step** control.
6. In **Graph**, select concentration, heat, abundance, equilibrium or sensitivity plots.
7. Optionally use **Trajectory averaging** and select **Averaged results** in the abundance graph.
8. Export numerical tables and figures; save `.kaz` sessions to reopen results in separate windows.

## Capabilities

| Area | Implementation |
|---|---|
| Nuclear data | ENDF-6 decay/reaction readers, rawmacs, ENDF/B-VII.1 and EAF-2010 MACS tables, TALYS, JSON/CSV user data |
| Networks | Nuclide selection, decay-type and reaction-MT filters, half-life filters, manual decay channels, boundary diagnostics |
| Evolution | Adaptive sparse CRAM-16; Padé through SciPy; physical-time trajectories |
| Abundance | Isotope inventories and sums by mass number A; absolute, relative and capture-cross-section-weighted displays |
| Averaging | Uniform physical-time average; exponential irradiation-duration mixture; exponential neutron-exposure mixture |
| Analysis | Equilibrium search, flux sweeps, OAT sensitivity, strongly connected components and elementary cycles |
| Output | CSV/TXT tables, Matplotlib figures, saved result sessions |

ENDF/B, JEFF, JENDL and TENDL data can be configured through the ENDF-6 layout. Available channels and energies depend on the actual files. A library name does not guarantee complete network coverage.

## Examples and tests

`examples/sprocess_termination.py` constructs a small Pb–Bi–Po termination network and exports its evolution. `examples/heat_vs_flux.py` computes equilibrium heat across a flux grid. These scripts require local nuclear data; inspect their settings before using them as a scientific model.

```bash
python examples/sprocess_termination.py
python examples/heat_vs_flux.py
python -B -m unittest discover -s tests -p "test_*.py"
```

The regression tests use synthetic networks or temporary fixtures and do not require distributing `xsdir`. GUI tests use offscreen Qt.

## Documentation

The Sphinx site is built from the English documentation in `docs/`. It does not import the application or require nuclear data.

```bash
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` locally. `.readthedocs.yaml` configures hosted builds; the GitHub workflow checks the documentation. Hosting requires connecting the repository to a Read the Docs account.

## Scientific scope

Evolution solves a fixed-matrix model. It does not automatically model stellar structure, time-dependent neutron production, mixing between stellar zones, or stellar corrections to decay rates. Averaging combines stored trajectories; it does not repair missing reactions or boundary losses.

The accepted API name `cram48` currently calls the same adaptive CRAM-16 solver. It is **not** an independent CRAM-48 implementation. MMPA is not implemented. See the methods page before interpreting a method comparison.

## Project publications

Bibliographic details supplied by the project authors:

1. Kenzhebayev N., Khassanov M., Spassyuk R., Anarbek D., Aimuratov Y., Abishev M. *Just Beyond the S-Process Termination Point: Nucleosynthesis of Lead–Bismuth Cyclic Reactions.* Galaxies, 14, 46 (2026).
2. Spassyuk R., Anarbek D., Khassanov M., Aimuratov Y., Kenzhebayev N., Abishev M. *Sensitivity of the s-process termination point to neutron capture cross sections and irradiation parameters.* Advances in Nuclear Science and Applications, 2, 30–40 (2026).

The earlier [nr-tool project](https://github.com/Nurzat-89/nr-tool) provides historical context.

## Authors, funding and license

Developed at the Fesenkov Astrophysical Institute, Almaty, Kazakhstan, within grant AP23488136, *The Energy Release of Cyclic Reactions of the Kilonova GW170817–GRB170817A–AT2017gfo* (Principal Investigator: Yerlan Aimuratov).

Core team: N. Kenzhebayev, M. Khassanov, R. Spassyuk, D. Anarbek and M. Abishev.

This checkout has no project license file. A public repository alone does not grant an open-source license; the authors must select and add the appropriate license. Nuclear datasets retain their own terms and attribution requirements.
