"""
examples/heat_density_flux.py
========================
Demonstration of CORRECT residual heat configuration:

  X axis — neutron flux Φ,
  Y axis — equilibrium decay heat density (plateau slice).

For each flux value, the matrix is rebuilt (reaction rates ∝ Φ),
the system reaches equilibrium, and total heat (reactions + decays)
is calculated based on the equilibrium concentration slice.

Run from NuMatRx folder:
    python examples/heat_vs_flux.py
"""


import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core import (
    BurnupConfig, BurnupMatrix, IsotopeBuilder,
    equilibrium_heat_vs_flux, plot_heat_vs_flux,
)

from nuclear_data import EndfPathHolder, NuclearDataProvider

DB_DIR = os.path.join(os.path.dirname(__file__), "user_db")


def main():
    # Locate xsdir relative to this file (examples/../xsdir)
    xsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "xsdir"))
    EndfPathHolder.init_xsdir(xsdir)  # also auto-discovers xsdir/TALYS

    # Locate user_db relative to this file (examples/user_db)
    user_db = os.path.abspath(os.path.join(os.path.dirname(__file__), "user_db"))

    # Set up the data provider with required libraries
    # MACS priority: user → ENDF/B-VII.1 tabulated → EAF-2010 → rawmacs → TALYS → ENDF computed
    nd = NuclearDataProvider(
        user_dir         = EndfPathHolder.set_USER_DB_DIR(user_db),                      # User modifications database
        dec_dir          = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0"),           # Decay data folder
        neu_dir          = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0"),           # Neutron cross-sections (σ(E) + endf_computed)
        macs_file        = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7"),                # rawmacs.txt
        macs_lib         = "Endfb7",
        endfb71_macs_dir = str(EndfPathHolder.get_ENDFB71_MACS_DIR()),                  # Multi-kT ENDF/B-VII.1 MACS
        eaf2010_macs_dir = str(EndfPathHolder.get_EAF2010_MACS_DIR()),                  # Multi-kT EAF-2010 MACS
        talys_dir        = str(EndfPathHolder.get_TALYS_DIR()),                          # TALYS σ(E) + MACS 30 keV
        macs_priority    = ["user", "endfb71", "eaf2010", "rawmacs", "talys", "endf_computed"],
    )






    builder = IsotopeBuilder(nd)

    # Define isotope chain specification
    specs = [
        (81, "Tl", 206), (81, "Tl", 207),
        (82, "Pb", 206), (82, "Pb", 207), (82, "Pb", 208),
        (82, "Pb", 209), (82, "Pb", 210), (82, "Pb", 211),
        (83, "Bi", 209), (83, "Bi", 210), (83, "Bi", 211),
        (84, "Po", 210), (84, "Po", 211),
    ]


    # Build isotope objects evaluating cross-sections at kT = 30 keV
    isotopes = builder.build_range(specs, E_eV=30_000.0, mt_filter=[102])



    # Configure burnup options
    cfg = BurnupConfig()
    cfg.enable_specific_reactions([102])
    cfg.set_decay_time_filter("all")
    cfg.set_open_boundary(True)   # open boundary: reactions outside the list are sinks

    # Set initial concentration: pure Pb-206
    names = [i.name for i in isotopes]
    N0 = np.zeros(len(isotopes))
    N0[names.index("Pb-206")] = 1.0

    # Define log-spaced neutron flux grid from 10^14 to 10^28 n/cm²/s



    fluxes = np.logspace(14, 28, 8)

    print("Calculating equilibrium heat for each flux level...")
    fl, heats = equilibrium_heat_vs_flux(isotopes, cfg, N0, fluxes)



    # t_plateaus = [1e6, 1e8, 1e10, 1e12, 1e14, 1e16, 1e18, 1e20]
    # fl, heats = equilibrium_heat_vs_flux(isotopes, cfg, N0, fluxes, t_plateau=t_plateaus)


    # t_plateaus = [1e6, 1e8, 1e10, 1e12, 1e14, 1e16, 1e18, 1e20]
    # fl, heats = equilibrium_heat_vs_flux(isotopes, cfg, N0, fluxes, t_plateau=t_plateaus)

    print(f"\n{'Flux, n/cm²/s':>14} {'Q_eq, MeV/s':>16}")
    for phi, q in zip(fl, heats):
        print(f"{phi:14.2e} {q:16.4e}")

    here = os.path.dirname(__file__)

    # ── Q vs Phi plot ─────────────────────────────────────────────────
    png = plot_heat_vs_flux(
        fl, heats,
        filename=os.path.join(here, "heat_density_flux.pdf"),
        title="Equilibrium Heat Release (closed chain)")
    print("\nQ vs Φ plot saved to →", png)

    # ── Q(t) for several fluxes — plateau verification ────────────────
    # Pick 3 representative fluxes: small, medium, large
    check_fluxes = [fluxes[0], fluxes[len(fluxes) // 2], fluxes[-1]]
    times_qt = np.logspace(-3, 25, 120)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#2980b9", "#27ae60", "#c0392b"]

    for phi, color in zip(check_fluxes, colors):
        cfg_phi = BurnupConfig()
        cfg_phi.enable_specific_reactions([102])
        cfg_phi.set_decay_time_filter("all")
        cfg_phi.set_open_boundary(False)
        cfg_phi.flux = float(phi)

        mtx = BurnupMatrix(isotopes, cfg_phi)
        mtx.build()

        qt = [mtx.heat_release(mtx.solve(N0, float(t))) for t in times_qt]
        ax.plot(times_qt, qt, color=color, linewidth=1.6, label=f"Φ = {phi:.1e}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Heat release Q(t)  [MeV/s per nucleus]")
    ax.set_title("Heat release vs time — plateau verification (closed chain)")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    qt_path = os.path.join(here, "heat_qt.pdf")
    fig.savefig(qt_path, dpi=140)
    plt.close(fig)
    print("Q(t) plot saved to →", qt_path)


if __name__ == "__main__":
    main()
