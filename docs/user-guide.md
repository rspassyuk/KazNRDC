# User guide — NuMatRx

## Installation and data

Install Python 3.11+ and create an isolated environment in the repository root. Install `requirements-endf.txt` for GUI, numerical packages and ENDF parsing; `requirements.txt` alone omits the optional ENDF reader. Start with `python model/core.py` using that environment.

Nuclear libraries are supplied separately. The configured layout is:

```text
xsdir/
  endf-6/
    ENDFB-VIII.0/
      decay/
      neutrons/
    JEFF-library-name/
      decay/
      neutrons/
  MACS/
    rawmacs.txt
    ENDFB71-EAF2010/
      ENDFB71/
      EAF2010/
  TALYS/
    recommended-ng-cross-sections/
    recommended-na-cross-sections/
    recommended-np-cross-sections/
    talys_macs_30_keV/
user_db/
  isotopes.json
```

Use the actual library folder names offered by your data distribution. Startup expects `xsdir/endf-6` and `xsdir/MACS`; not every optional source needs to be installed. Empty folders do not provide evaluated data. For initialization errors, inspect **Terminal → output**, directory spelling and reader dependencies. No database download is performed automatically.

`user_db/isotopes.json` supplies overrides. The entry key includes Z, symbol and mass, such as `82-Pb-210`. Decay entries use seconds, branching fractions and MeV. MACS entries use `kT_keV` and `sigma_mb`; this differs from the internal barn unit. Cross-section CSV files use `E_eV` and `sigma_barn`. Follow the shipped user-data schema, keep an original copy and verify loaded values in Periodic Table after restarting/reloading.

## Application layout

The left navigation switches between **Periodic Table**, **Isotope Chart**, **Results** and **Graph**. Use **View** to select the light/dark theme or restore a dock through **Panels**. Plot color selection in Graph is separate from the application widget theme. Numeric fields support scientific notation such as `1e14`.

The application keeps the latest calculation results. Editing controls does not recalculate an existing result. Record the settings before launching another run.

## First calculation

1. Open Periodic Table and inspect a few isotopes and their decay/reaction data.
2. Open Isotope Chart. Select a small network using the map or **Manual** labels.
3. Choose the library, energy and reaction/decay channels. For astrophysical capture studies choose the MACS preference and document the evaluation energy.
4. Click **Build Burnup Matrix**.
5. Open Results → **Isotopes & Matrix** and inspect membership and rates.
6. Load the matrix isotopes into **Initial Conditions**. Set at least one positive inventory; Fe-56 = 1 defines a normalized seed inventory.
7. Set flux or neutron density, energy, final time and solver. Run **Evolution**.
8. Move **Time step** to inspect saved compositions. In Graph select **Abundance vs mass number A** with source **Time slice**.
9. Export numerical results and save a `.kaz` session before replacing the run.

Start small to confirm your data and units, then expand the network. The selected rectangular chart range is not necessarily a single reaction chain.

## Periodic Table

Select an element and isotope to inspect decay and neutron reaction data. Use the two information panes to compare sources, lifetimes and channels. Displayed source data do not by themselves confirm that a channel is enabled in the working burnup matrix.

## Isotope Chart

Selection coordinates use Z and neutron number N = A − Z. The rectangle limit is 5,000 candidate isotopes; practical speed and memory also depend on network connectivity. Manual labels define a selected list; same-element ranges such as `Tl206..Tl216` are supported.

Decay and neutron-channel controls determine which transitions are requested for calculation. Arrow visibility is a display aid. Some stored links can still appear in graph-analysis tools even when filtered out of the matrix.

Right-click exclusion and the excluded-isotope controls modify workspace membership. Reactor mode favors pointwise cross sections; astrophysical mode favors MACS. Missing input can trigger source fallback, so inspect the selected source rather than assuming one library provides every rate.

## Isotopes & Matrix

The main table contains isotope properties and transition information. It is not the complete numeric square matrix A. Use **Ctrl+F** to find/filter isotopes. Copy exports selected cells (or visible data when no cells are selected); the table CSV export contains isotope properties. Numerical matrix export is available through the core API.

Use the generic isotope input and **Add** to extend membership. Remove selected isotopes when appropriate. Isotopes referenced by manual channels can be protected or reintroduced so those transitions remain well-defined.

Use **Refresh Matrix View** after changing matrix settings. The preview's flux and the Evolution flux are separate controls. Data edits in a display table should not be treated as confirmed rate overrides; use the user database or manual decay editor and verify the resulting matrix.

### Decay types and half-lives

Select only the decay types and neutron MT channels required by your model. Disabling a decay does not imply deleting its parent isotope. The half-life presets keep eligible decays with half-life no greater than the selected threshold:

| Preset | Upper half-life, seconds |
|---|---:|
| all | No bound |
| yearly | 31,560,000 |
| daily | 86,400 |
| hourly | 3,600 |
| minute | 60 |
| second | 1 |

These are filters, not numerical integration steps. Excluding a physically relevant long-lived branch changes the model. Check boundary settings before using conservation as a diagnostic.

### Manual decay channels

Enter parent, daughter, mode, half-life in seconds and branch fraction in the manual-channel editor. A blank daughter represents a sink. Add the row and refresh the matrix. Manual transitions bypass global decay filters and replace matching evaluated branches. The generic mechanism can enable a selected alpha branch while other alpha channels remain disabled; there is no isotope-specific shortcut.

Validate the daughter, mass/charge change, branch fraction and lifetime yourself. If no matching evaluated energy is available, the manual channel may not contribute a known heat value.

## Evolution

Initial inventories may be entered manually or loaded from the supported two-column CSV. A direct reload of matrix isotopes into Initial Conditions can reset entries; check them before running.

Choose flux in neutrons/cm²/s or neutron density in neutrons/cm³. The density option converts through the entered energy; see [Methods](methods.md). Set the irradiation endpoint in seconds and use adaptive CRAM-16 for the primary calculation. Padé is an alternative for independent checks on manageable networks.

Targets can restrict the source-to-target calculation network. Leave the field empty for the full selected workspace. Inspect the resulting isotope list rather than treating Targets as only a visual filter.

Time step selects a saved time, normally among 400 outputs. The Evolution table shows individual isotope values and mass-number totals. Absolute preserves input units, Relative normalizes by the inventory total, and sigma N weights with capture cross sections. The extra A/N(A) rows group isotopes by mass number, not by chemical element.

CSV/TXT Evolution exports contain times, isotope inventories, fractions and mass-number summaries. They use their export-defined columns; changing a display mode does not necessarily change these files into sigma N exports. Read the header.

## Graph

Select a graph type and click **Plot / Refresh**. Use the Matplotlib toolbar for zoom, pan and figure saving. The mode buttons above the graph follow the application theme.

| Graph | Controls |
|---|---|
| Concentrations vs time | Isotope filter; Absolute / Relative / sigma N |
| Heat release vs time | MeV/s or W/cm³; density is required for volumetric conversion |
| Abundance vs mass number A | Time slice or Averaged results; Absolute / Relative / sigma N |
| Equilibrium heat vs flux | MeV/s or W/cm³ |
| Equilibrium N vs flux | Relative / Absolute / sigma N |
| Sensitivity Analysis | Latest sensitivity calculation |

For time-slice abundance, move the slider to change the selected saved time. For averaged abundance, the slider does not scan different averages; recompute averaging with different interval or distribution parameters instead. Minimum-value and isotope display filters do not rerun the physics calculation. Missing cross sections can create gaps in averaged sigma N curves.

## Trajectory averaging

First calculate Evolution or open a saved result. Choose one distribution:

| Selection | Inputs | Meaning |
|---|---|---|
| Time average over interval | Start time and End time, seconds | Uniform temporal mean over the selected interval |
| Exponential irradiation durations | Mean duration t0, seconds | Mixture of irradiation endpoints weighted by an exponential duration distribution |
| Exponential neutron exposure | Mean exposure tau0, mbarn⁻¹ | Equivalent exposure-distribution model for constant positive flux |

For the uniform mode, an empty End time uses the final saved time. The interval must satisfy 0 ≤ start < end ≤ final time. For the exponential time mode, integration runs from zero to the trajectory endpoint; t0 does not set that endpoint.

Click **Average results**. Inspect **Isotopes**, **Mass numbers** and the diagnostics. **Mean N(A)** sums the averaged inventories of equal-A isotopes. **Mean fraction(A)** divides that sum by the total averaged inventory. **Sum sigma_i N_i** sums isotope-wise weighted values, in barn times input inventory.

The exponential modes do not renormalize the uncovered tail. A high coverage does not guarantee adequate temporal resolution; inspect the weight error and compare finer/longer trajectories. Significant negative input concentrations block averaging. See the formulas and normalization in [Methods](methods.md).

Use **Export averaged results** for CSV/TXT with parameters and both tables. **Save averaged session** stores the run and its latest average. To plot it, select Graph → Abundance vs mass number A → **Averaged results**. Repeating averaging replaces the attached average; save separate sessions for parameter comparisons.

## Equilibrium

Configure the flux grid, initial composition, energy and either automatic plateau search or explicit evaluation times. Each flux is a separate calculation with the initial inventory. A flux sweep is not sequential irradiation with a changing flux.

Inspect returned times and inventories before calling a result equilibrium. Export its table and use the equilibrium graph types. Boundary losses and search limits can prevent a meaningful nonzero plateau.

## Sensitivity

Choose a concentration or heat target, fixed evaluation time or plateau, perturbation size and parameter types. The available scans cover cross sections, decay constants, flux/energy and library selection. Start with a limited scope; each perturbation requires another solve.

Inspect the baseline and rank by absolute sensitivity. A positive/negative sign describes the local target response, not whether a reaction is useful or harmful. Library deviations depend on the chosen reference and fallbacks. Verify configuration independently of the Evolution tab; not all manual overrides and filters are shared.

## Cycle analysis

**Find Cycle Reactions** searches the linked graph for elementary cycles and strongly connected components. A graph cycle shows a possible return path; it does not quantify the amount of material circulating. Compare its links with active matrix filters and rate magnitudes before drawing a physical conclusion.

## Sessions and exports

A `.kaz` file stores Evolution, the attached average, optional equilibrium results and density. Open it in a separate result window without replacing the live workspace. Multiple sessions can be inspected side by side.

Sessions are result archives, not a full editable GUI snapshot: sensitivity scans, unsaved widget settings and every workspace edit are not restored. Existing trajectories are not recalculated on opening. The format uses Python pickle, so open only trusted session files.

Save plot images from the Matplotlib toolbar. Numerical exports and figures are separate artifacts. Include their metadata and initial conditions when sharing a scientific result.

## Troubleshooting and reproducibility

| Symptom | Check |
|---|---|
| Nuclear-data initialization failed | xsdir directories, library name, installed endf reader and Terminal output |
| Empty average graph | Average results completed; Graph source is Averaged results; required cross sections exist |
| Averaging blocked by negative concentration | Recalculate Evolution; inspect solver errors and the reported isotope/time |
| High uncovered tail | Extend the trajectory or choose a justified smaller mean duration/exposure |
| Weight integration error | Refine the saved time grid through the API and repeat the convergence check |
| Unexpected inventory loss | Daughters outside the matrix, open channels, omitted reactions and stoichiometry |
| s-process curve disagrees with a reference | Seed composition, metallicity assumptions, MACS energy, density/flux, duration distribution, branching and normalization |

Record code version, data versions, source priorities, selected isotopes, filters, manual channels, boundaries, N0, flux/density, energy, solver, saved times, averaging parameters and display normalization. This model does not supply a universal flux/duration pair that reproduces every stellar s-process component.

## Shared nuclear-data directory

The application first looks for `xsdir/` inside the project. If it is absent, it also checks for `xsdir/` beside the project folder. Nuclear libraries and private `user_db/` overrides are not included in the public source distribution. Existing `.kaz` sessions retain their format.
