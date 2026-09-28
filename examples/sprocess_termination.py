"""
examples/sprocess_termination.py
================================
S-process termination cycle in the lead-bismuth region.

Physics: Neutron captures move the material flow up through the Pb isotopes to Pb-209,
which beta-decays into Bi-209. Capture on Bi-209 yields Bi-210 -> Po-210,
and Po-210 alpha-decays back to Pb-206 — closing the nucleosynthesis cycle.
Tl-206/207 appear from weak alpha branches of Bi-210/Bi-211.

Chain structure:
    Pb-206 →(n,γ) Pb-207 →(n,γ) Pb-208 →(n,γ) Pb-209 ─β⁻→ Bi-209
    Bi-209 →(n,γ) Bi-210 ─β⁻→ Po-210 ─α→ Pb-206   (LOOP CLOSURE)
    Bi-210 ─α(weak)→ Tl-206 ─β⁻→ Pb-206
    Pb-210 ─β⁻→ Bi-210                 (branching via Pb-210)
    Pb-211 ─β⁻→ Bi-211 ─α→ Tl-207 ─β⁻→ Pb-207
    Po-211 ─α→ Pb-207

Run from KazNRDC folder:
    python examples/sprocess_termination.py
"""

import os
import sys
import numpy as np

# Adjust path to import KazNRDC packages when executed from examples/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core import BurnupConfig, BurnupMatrix, IsotopeBuilder, plot_concentrations, plot_heat
from core import CycleAnalyzer
from nuclear_data import EndfPathHolder, NuclearDataProvider

DB_DIR = os.path.join(os.path.dirname(__file__), "sprocess_db")


def line(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def build_chain():
    """Assembles all 13 isotopes of the termination cluster using the specified data layers."""

    # Locate xsdir relative to this file (examples/../xsdir)
    xsdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "xsdir"))
    EndfPathHolder.init_xsdir(xsdir)  # also auto-discovers xsdir/TALYS

    # Locate user_db relative to the project root
    user_db = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "user_db"))

    # Configure the unified provider
    # Decay priority: user → ENDF-6 (MF=8/MT=457)
    # MACS priority:  user → ENDF/B-VII.1 tabulated → EAF-2010 → rawmacs → TALYS → ENDF computed

    nd = NuclearDataProvider(
        user_dir         = EndfPathHolder.set_USER_DB_DIR(user_db),                      # User modifications database
        dec_dir          = EndfPathHolder.set_DEFAULT_DEC_DIR("ENDFB-VIII.0"),           # Decay directory
        neu_dir          = EndfPathHolder.set_DEFAULT_NEU_DIR("ENDFB-VIII.0"),           # Neutron cross-sections (σ(E) + endf_computed)
        macs_file        = EndfPathHolder.set_DEFAULT_MACS_DIR("Endfb7"),                # rawmacs.txt
        macs_lib         = "Endfb7",
        endfb71_macs_dir = str(EndfPathHolder.get_ENDFB71_MACS_DIR()),                  # Multi-kT ENDF/B-VII.1 MACS
        eaf2010_macs_dir = str(EndfPathHolder.get_EAF2010_MACS_DIR()),                  # Multi-kT EAF-2010 MACS
        talys_dir        = str(EndfPathHolder.get_TALYS_DIR()),                          # TALYS σ(E) + MACS 30 keV
        dec_priority     = ["user", "endf6"],
        macs_priority    = ["user", "endfb71", "eaf2010", "rawmacs", "talys", "endf_computed"],
    )

    specs = [
        (81, "Tl", 206), (81, "Tl", 207),
        (82, "Pb", 206), (82, "Pb", 207), (82, "Pb", 208),
        (82, "Pb", 209), (82, "Pb", 210), (82, "Pb", 211),
        (83, "Bi", 209), (83, "Bi", 210), (83, "Bi", 211),
        (84, "Po", 210), (84, "Po", 211),
    ]

    builder = IsotopeBuilder(nd)
    # The s-process takes place at kT ~ 30 keV (stellar thermal peak),
    # so cross-sections are evaluated at this specific energy boundary
    isotopes = builder.build_range(specs, E_eV=30_000.0, mt_filter=[102])
    return nd, isotopes


def show_connections(isotopes):
    line("CHAIN STRUCTURE: Transmutation and decay pathways")
    for iso in isotopes:
        parts = []
        for d in iso.getListOfDecays():
            tgt = d.product.name if d.product else "(out of network)"
            parts.append(f"{d.mode}→{tgt} [BR={d.branch:.4g}]")
        for r in iso.getListOfReactions():
            tgt = r.product.name if r.product else "(out of network)"
            parts.append(f"{r.name}→{tgt} [σ={r.sigma_barn:.3g}b]")
        hl = "stable" if iso.half_life_s is None else f"T½={iso.half_life_s:.3g}s"
        print(f"  {iso.name:9s} {hl:16s} | " + "; ".join(parts) if parts
              else f"  {iso.name:9s} {hl:16s} | —")


def run(isotopes, flux, label):
    line(f"SIMULATION: {label}  (Φ = {flux:.1e} n/cm²/s)")

    cfg = BurnupConfig(flux=flux)
    cfg.enable_specific_reactions([102])   # Only (n,γ) reactions are included for the s-process
    cfg.set_decay_time_filter("all")       # All decay transitions are active

    mtx = BurnupMatrix(isotopes, cfg)
    mtx.build()
    print("Matrix assembled. Structural validation:", mtx.validate())
    if mtx.missing:
        print("Missing cross-section links for:", mtx.missing)

    # Set initial concentration: 100% Pb-206 (entry node into the cycle)
    names = mtx.names()
    N0 = np.zeros(len(names))
    N0[names.index("Pb-206")] = 1.0



    # Dense log-spaced time grid spanning from 1 ms to 10^15 seconds (approx. 31.7 million years)
    grid = [0.0] + list(np.logspace(-3, 15, 100))
    times, traj = mtx.evolve(N0, grid, method="cram16")

    # Define milestone points for reporting concentration snapshots to console
    print_times = [0, 1e-3, 1.0, 86400, 365*86400, 1e9, 1e15]
    print_labels = ["0", "1 ms", "1 s", "1 day", "1 year", "31.7 years", "31.7 Myr"]

    watch = ["Pb-206", "Pb-207", "Pb-208", "Pb-209", "Bi-209",
             "Bi-210", "Po-210", "Tl-206"]

    print(f"\n{'isotope':9s}", *[f"{l:>10s}" for l in print_labels])
    for w in watch:
        idx = names.index(w)
        vals = []
        for pt in print_times:
            # Reconstruct concentrations explicitly at the target timestamp milestone
            N = mtx.solve(N0, float(pt), "cram16") if pt > 0 else N0.copy()
            vals.append(N[idx])
        print(f"{w:9s} " + " ".join(f"{v:10.3e}" for v in vals))

    print(f"\nEnergy release at the end of the simulation: "
          f"{mtx.heat_release(traj[-1]):.4e} MeV/s per initial nucleus")

    return mtx, times, print_labels, traj


def export(mtx, times, labels, traj):
    line("DATA EXPORT")
    here = os.path.dirname(__file__)
    out_mtx = os.path.join(here, "sprocess_matrix.txt")
    out_evo = os.path.join(here, "sprocess_evolution.txt")
    mtx.export_txt(out_mtx, mode="total")
    mtx.export_evolution_txt(out_evo, times, traj)
    print("Burnup matrix file    →", out_mtx)
    print("Evolution pathway file →", out_evo)


    # Plot generation section
    try:
        png_c = plot_concentrations(
            mtx.names(), times, traj,
            filename=os.path.join(here, "sprocess_concentrations.png"),
            mode="relative", logy=True, logx=True,
            title="s-process: Isotope Concentration Evolution (Relative Fraction)")

        # Generates the mixture decay heat curve matching the exact logarithmic time step range
        png_h = plot_heat(
            times, mtx.heat_curve(traj),
            filename=os.path.join(here, "sprocess_heat.png"),
            logy=True, logx=True,  # Set logx=True to align time domain directly with concentrations
            title="s-process: Mixture Decay Heat Release Profile")

        print("Concentration plot N(t) →", png_c)
        print("Decay heat profile Q(t) →", png_h)
    except ImportError as e:
        print("Plotting routines skipped:", e)




if __name__ == "__main__":
    nd, isotopes = build_chain()
    show_connections(isotopes)

    # Automatic detection of cyclic reaction/decay loops in the network
    CycleAnalyzer(isotopes).print_report()

    # Scenario A: Pure radioactive decay simulation (neutron flux is turned off)
    run(isotopes, flux=0.0, label="PURE DECAY SCENARIO")

    # Scenario B: High-intensity stellar capture loop environment
    mtx, times, labels, traj = run(isotopes, flux=1e14, label="INTENSE CAPTURE SCENARIO (Active Loop)")

    export(mtx, times, labels, traj)
    line("Done")
