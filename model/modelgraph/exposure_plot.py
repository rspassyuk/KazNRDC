"""Exposure-averaged abundance curves for graph views."""
import numpy as np
from core.exposure import mass_values


def plot_exposure(ax, result, mode="absolute", colors=None):
    colors = colors or dict(bg="#ffffff", fg="#16222e", accent="#2b7fc4", grid="#d2dde6")
    ax.set_facecolor(colors["bg"])
    ax.figure.set_facecolor(colors["bg"])
    aa, values = mass_values(result, mode)
    # Keep negative numerical noise visible; use a linear scale if present.
    ax.plot(aa, values, "-o", color=colors["accent"], markersize=2, linewidth=1.2)
    finite = values[np.isfinite(values)]
    if len(finite) and not np.any(finite < 0):
        positive = finite[finite > 0]
        if len(positive) and positive.max() / positive.min() > 100:
            ax.set_yscale("log")
    if not len(finite):
        ax.text(0.5, 0.5, "No finite values to plot.\nCapture cross sections are unavailable for these mass groups."
                if mode == "sigma_n" else "No finite abundance values to plot.",
                transform=ax.transAxes, ha="center", va="center", color=colors["fg"], usetex=False)
    elif np.any(~np.isfinite(values)):
        ax.text(0.02, 0.98, "Gaps indicate unavailable cross sections.",
                transform=ax.transAxes, va="top", color=colors["fg"], usetex=False)
    ax.set_xlabel("Mass number A", color=colors["fg"])
    ax.set_ylabel(dict(absolute="Mean N(A) [input inventory units]",
                       relative="Mean N(A) / total mean N",
                       sigma_n="Sum sigma_i * mean N_i [barn * inventory]")[mode],
                  color=colors["fg"])
    mode_name = getattr(result.config, "mode", "exposure")
    if mode_name == "time_uniform":
        title = (f"Time average: {result.metadata['time_start_s']:g} to "
                 f"{result.metadata['time_end_s']:g} s")
    elif mode_name == "time_exponential":
        title = (f"Exponential time mixture: t0={result.config.mean_time_s:g} s; "
                 f"end={result.metadata['time_end_s']:g} s")
    else:
        title = f"Exposure-averaged: tau0={result.config.tau0:g} mbarn^-1"
    ax.set_title(title, color=colors["fg"], fontsize=11)
    ax.tick_params(colors=colors["fg"])
    for spine in ax.spines.values():
        spine.set_color(colors["grid"])
    ax.grid(alpha=0.3, color=colors["grid"])

    # Use native Matplotlib text even when another graph enabled global LaTeX.
    from matplotlib.text import Text
    for text in ax.figure.findobj(Text):
        text.set_usetex(False)
