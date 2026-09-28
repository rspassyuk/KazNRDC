"""
core/plotting.py
================
Plotting utilities for burnup calculations. Saves results as PNG files.

matplotlib is an optional dependency. If it is not installed, these functions
will raise an informative ImportError, but the rest of the application remains functional.
"""

from __future__ import annotations
import re
from typing import List, Optional

try:
    import matplotlib
    matplotlib.use("Agg")  # Non-interactive backend — outputs directly to disk file
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    plt = None
    _HAS_MPL = False


def _check():
    if not _HAS_MPL:
        raise ImportError("Matplotlib is required for plotting utilities: pip install matplotlib")


def _seconds_to_unit(times):
    """Selects an optimal, human-readable time scale unit based on the maximum value."""
    tmax = max(times) if times else 1.0
    if tmax >= 3.156e7:
        return [t / 3.156e7 for t in times], "years"
    if tmax >= 86400:
        return [t / 86400 for t in times], "days"
    if tmax >= 3600:
        return [t / 3600 for t in times], "hours"
    return list(times), "seconds"


def _parse_isotope_base(name: str) -> str:
    """
    Extracts the base isotope name, stripping isomeric state tags.
    Example: 'Pb-209m1' or 'Pb-209m' -> 'Pb-209'
    """
    match = re.match(r"^([A-Za-z]+-\d+)", name.strip())
    return match.group(1) if match else name.strip()


def plot_concentrations(
    names: List[str],
    times: List[float],
    trajectories: List,
    filename: str = "concentrations.png",
    mode: str = "relative",
    sigma_weights: Optional[List[float]] = None,
    stable_isotopes: Optional[List[str]] = None,
    title: str = "Isotope Concentration Evolution",
    logy: bool = False,
    logx: bool = False,
    min_fraction: float = 1e-30,
    ymin: Optional[float] = 1e-31,  # Added minimum Y-axis limit parameter
):
    """
    Plots the structural changes of your network across the selected time grid.
    """
    _check()
    import numpy as np

    M = np.array(trajectories, dtype=float)  # shape (T, N)
    t_plot, t_unit = _seconds_to_unit(times)

    if logx:
        t_plot = list(times)
        t_unit = "seconds"

    if mode == "relative":
        sums = M.sum(axis=1, keepdims=True)
        sums[sums == 0] = 1.0
        Y = M / sums
        ylabel = "Relative Atomic Fraction"
    elif mode == "sigmaN":
        if sigma_weights is None:
            raise ValueError("sigma_weights parameter is mandatory for 'sigmaN' visualization mode.")
        w = np.asarray(sigma_weights, dtype=float)
        Y = M * w[None, :]
        ylabel = "σ·N Reaction Contribution Metric [barns]"
    else:  # absolute
        Y = M
        ylabel = "Absolute Concentration"

    fig, ax = plt.subplots(figsize=(10, 6))
    peak = Y.max(axis=0)
    drawn = 0

    if stable_isotopes is None:
        stable_isotopes = ["Pb-206", "Pb-207", "Pb-208", "Bi-209"]

    # --- COLOR MAP: One consistent color per base ISOTOPE ---
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Extract only pure isotopes without meta-states (e.g., 'Pb-209')
    unique_base_isotopes = sorted(list(set(_parse_isotope_base(name) for name in names)))

    # Assign a strict color to each base isotope
    isotope_color_map = {iso: color_cycle[idx % len(color_cycle)] for idx, iso in enumerate(unique_base_isotopes)}

    # Marker set for distinguishing different states/isomers of the same isotope
    markers = ["o", "s", "^", "v", "<", ">", "p", "*", "D"]
    # Tracker to iterate marker indices per base isotope
    isotope_marker_counter = {iso: 0 for iso in unique_base_isotopes}

    for j, name in enumerate(names):
        if peak[j] < min_fraction:
            continue

        base_iso = _parse_isotope_base(name)
        color = isotope_color_map.get(base_iso, "gray")

        # Pick the marker for the current state and advance the counter for the next isomer
        marker_idx = isotope_marker_counter[base_iso]
        marker = markers[marker_idx % len(markers)]
        isotope_marker_counter[base_iso] += 1

        # Line style depends on the stability of the isotope
        if name in stable_isotopes:
            linestyle = "-"     # Solid line for stable nuclei
            linewidth = 2.0
        else:
            linestyle = "--"    # Dashed line for radioactive nuclides
            linewidth = 1.3

        ax.plot(
            t_plot, Y[:, j],
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            marker=marker,
            markersize=4,
            markevery=max(1, len(t_plot) // 15),
            label=name
        )
        drawn += 1

    ax.set_xlabel(f"Time [{t_unit}]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if logy:
        ax.set_yscale("log")
        if ymin is not None:
            ax.set_ylim(bottom=ymin)  # Enforce bottom boundary to prevent zero-drop visually

    if logx:
        ax.set_xscale("log")
        valid_times = [t for t in t_plot if t > 0]
        ax.set_xlim(left=min(valid_times) if valid_times else 1e-3, right=max(t_plot))

    ax.grid(True, alpha=0.3, which="both" if logx or logy else "major")

    if drawn <= 25:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, borderaxespad=0.)

    fig.tight_layout()
    fig.savefig(filename, dpi=140)
    plt.close(fig)
    return filename


def plot_heat(
    times: List[float],
    heat_values: List[float],
    filename: str = "heat.png",
    title: str = "Mixture Decay Heat Profile",
    logy: bool = False,
    logx: bool = False,
    ymin: Optional[float] = 1e-31,  # Added minimum Y-axis limit parameter
):
    """
    Plots total energy release [MeV/s per initial nucleus] as a function of time.
    """
    _check()
    t_plot, t_unit = _seconds_to_unit(times)

    if logx:
        t_plot = list(times)
        t_unit = "seconds"

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(
        t_plot, heat_values,
        marker="s", markersize=4,
        markevery=max(1, len(t_plot) // 20),
        linewidth=1.6, color="#c0392b"
    )
    ax.set_xlabel(f"Time [{t_unit}]")
    ax.set_ylabel("Decay Heat Release [MeV/s per Nucleus]")
    ax.set_title(title)

    if logy:
        ax.set_yscale("log")
        if ymin is not None:
            ax.set_ylim(bottom=ymin)  # Enforce bottom boundary

    if logx:
        ax.set_xscale("log")
        valid_times = [t for t in t_plot if t > 0]
        ax.set_xlim(left=min(valid_times) if valid_times else 1e-3, right=max(t_plot))

    ax.grid(True, alpha=0.3, which="both" if logx or logy else "major")
    fig.tight_layout()
    fig.savefig(filename, dpi=140)
    plt.close(fig)
    return filename


def plot_heat_vs_flux(
    fluxes: List[float],
    heats: List[float],
    filename: str = "heat_vs_flux.png",
    title: str = "Equilibrium Plateau Heat vs. Neutron Flux",
    logx: bool = True,
    logy: bool = True,
):
    """
    Plots steady-state plateau decay heat density vs. incident neutron flux.
    """
    _check()
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(fluxes, heats, marker="o", markersize=5, linewidth=1.8, color="#2c3e9e")
    ax.set_xlabel("Neutron Flux Φ [n/(cm²·s)]")
    ax.set_ylabel("Equilibrium Decay Heat [MeV/s per Nucleus]")
    ax.set_title(title)

    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(filename, dpi=140)
    plt.close(fig)
    return filename
