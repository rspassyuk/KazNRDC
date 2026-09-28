"""Physical-time averaging and legacy exposure averaging of saved trajectories."""
from dataclasses import dataclass, field
import csv
import json
import re
import numpy as np


@dataclass(frozen=True)
class ExposureConfig:
    tau0: float = 0.3
    tail_limit: float = 1e-3
    weight_rtol: float = 5e-3
    negative_atol: float = 1e-14
    negative_rtol: float = 1e-10
    mode: str = "exposure"
    start_s: float = 0.0
    end_s: float | None = None
    mean_time_s: float = 1e9


@dataclass
class ExposureResult:
    names: list
    mean: np.ndarray
    sigma_barn: np.ndarray
    config: ExposureConfig
    coverage: float
    tail: float
    weight_integral: float
    metadata: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


def mass_number(name):
    match = re.fullmatch(r"(?:\d+-)?[A-Za-z]+-?(\d+)(?:m\d*)?", name)
    if not match:
        raise ValueError(f"Cannot parse mass number: {name}")
    return int(match.group(1))


def capture_snapshot(names, isotopes):
    """Registry stores sigma_barn = total sigma * branch yield; sum once."""
    by_name = {iso.name: iso for iso in isotopes}
    values, sources = [], {}
    for name in names:
        iso = by_name.get(name)
        channels = [] if iso is None else [
            rx for rx in iso.getListOfReactions() if int(rx.getId()) == 102]
        values.append(sum(float(rx.sigma_barn) for rx in channels) if channels else np.nan)
        sources[name] = sorted({str(rx.source) for rx in channels if rx.source})
    return np.array(values, dtype=float), sources


def average_exposure(names, times, trajectories, flux, config,
                     sigma_barn=None, metadata=None, progress=None):
    """Finite-range trapezoid integral. Small negative noise is retained and reported."""
    names = list(names)
    t = np.array(times, dtype=float, copy=True)
    y = np.array(trajectories, dtype=float, copy=True)
    mode = getattr(config, "mode", "exposure")
    if mode not in ("exposure", "time_uniform", "time_exponential"):
        raise ValueError("Unknown averaging mode.")
    flux, tau0 = float(flux), float(config.tau0)
    if mode == "exposure" and (not np.isfinite(flux) or flux <= 0 or not np.isfinite(tau0) or tau0 <= 0):
        raise ValueError("Flux and mean exposure must be finite and positive.")
    for value in (config.tail_limit, config.weight_rtol,
                  config.negative_atol, config.negative_rtol):
        if not np.isfinite(value) or value < 0:
            raise ValueError("Diagnostic tolerances must be finite and non-negative.")
    if not names or len(set(names)) != len(names):
        raise ValueError("Isotope names must be nonempty and unique.")
    for name in names:
        mass_number(name)
    if t.ndim != 1 or len(t) < 2 or y.shape != (len(t), len(names)):
        raise ValueError("Expected at least two times and a time-by-isotope concentration array.")
    if not np.isfinite(t).all() or not np.isfinite(y).all():
        raise ValueError("Times and concentrations must be finite.")
    if t[0] != 0 or np.any(np.diff(t) <= 0):
        raise ValueError("Times must start at zero and increase strictly.")
    meta = dict(metadata or {})
    if mode == "exposure" and meta.get("constant_conditions") is False:
        raise ValueError("Exposure averaging supports constant conditions only.")
    if mode == "exposure" and meta.get("constant_conditions") is not True:
        raise ValueError("Confirm fixed-condition provenance before exposure averaging.")
    totals = y.sum(axis=1)
    if totals[0] <= 0 or np.any(totals < 0):
        raise ValueError("Non-positive initial inventory or negative trajectory totals.")
    tolerance = config.negative_atol + config.negative_rtol * max(abs(y[0]).max(), totals[0])
    if y.min() < -tolerance:
        step, isotope = np.unravel_index(np.argmin(y), y.shape)
        raise ValueError(
            f"Invalid Evolution result: {names[isotope]} has N={y[step, isotope]:.8g} "
            f"at t={t[step]:.8g} s (step {step}; tolerance {tolerance:.6g}). "
            "Rerun Evolution with the updated adaptive CRAM-16 solver, then average again. "
            "Opening an older session does not recalculate its concentrations.")
    sigma = (np.full(len(names), np.nan) if sigma_barn is None
             else np.array(sigma_barn, dtype=float, copy=True))
    if sigma.shape != (len(names),) or np.isinf(sigma).any() or np.any(sigma < 0):
        raise ValueError("Capture cross sections must be non-negative; use NaN for missing data.")
    if mode == "exposure":
        grid = t * (flux * 1e-27)
        scale = tau0
        if not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0):
            raise ValueError("Exposure grid overflow or underflow.")
    elif mode == "time_uniform":
        start = float(config.start_s)
        end = float(t[-1] if config.end_s is None else config.end_s)
        if not np.isfinite([start, end]).all() or not 0 <= start < end <= t[-1]:
            raise ValueError("Require 0 <= start < end <= the final saved time.")
        grid = np.concatenate(([start], t[(t > start) & (t < end)], [end]))
        y = np.column_stack([np.interp(grid, t, column) for column in y.T])
        scale = end - start
    else:
        grid = t
        scale = float(config.mean_time_s)
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("Mean duration t0 must be finite and positive.")
    density = (np.full_like(grid, 1.0 / scale) if mode == "time_uniform"
               else np.exp(-grid / scale) / scale)
    widths = np.diff(grid)
    weights = np.zeros_like(grid)
    weights[:-1] += widths / 2
    weights[1:] += widths / 2
    weights *= density
    mean = weights @ y
    coverage = 1.0 if mode == "time_uniform" else float(-np.expm1(-grid[-1] / scale))
    tail = 0.0 if mode == "time_uniform" else float(np.exp(-grid[-1] / scale))
    integral = float(weights.sum())
    if not np.isfinite(mean).all() or not np.isfinite(integral):
        raise ValueError("Non-finite quadrature result.")
    error = abs(integral - coverage) / max(coverage, np.finfo(float).tiny)
    warnings = []
    if tail > config.tail_limit:
        warnings.append("Uncovered distribution tail exceeds the configured tolerance; extend the run.")
    if error > config.weight_rtol:
        warnings.append("Saved grid is too coarse for the weight; refine the saved time grid.")
    if np.any(y < 0):
        warnings.append("Small negative numerical values were retained; no clipping was applied.")
    if abs(totals[-1] - totals[0]) > 1e-6 * totals[0]:
        warnings.append("Total inventory changed; inspect boundary losses and reaction stoichiometry.")
    if np.isnan(sigma).any():
        warnings.append("Missing capture cross sections are marked unavailable, not zero.")
    if meta.get("targets"):
        warnings.append("Input is a restricted calculation network; averaging uses all its computed isotopes.")
    if meta.get("legacy"):
        warnings.append("This saved run does not include complete library and boundary metadata.")
    if mean.sum() <= 0:
        raise ValueError("Averaged total is non-positive; relative abundance is undefined.")
    meta.update(flux_n_cm2_s=flux, time_start_s=float(grid[0] if mode == "time_uniform" else t[0]),
                time_end_s=float(grid[-1] if mode == "time_uniform" else t[-1]),
                averaging_mode=mode,
                distribution="Uniform" if mode == "time_uniform" else "Exponential",
                integration="Trapezoid on exposure grid" if mode == "exposure" else "Trapezoid on physical time grid",
                normalization="Interval duration" if mode == "time_uniform" else "None (finite-range integral)",
                concentration_units="Same as input inventory", sigma_units="barn")
    if mode == "exposure":
        meta["tau_max_mbarn_inv"] = float(grid[-1])
    elif mode == "time_exponential":
        meta["mean_duration_s"] = scale
    diagnostics = dict(initial_total=float(totals[0]), final_total=float(totals[-1]),
                       minimum_total=float(totals.min()), maximum_total=float(totals.max()),
                       minimum_concentration=float(y.min()), negative_count=int((y < 0).sum()),
                       negative_tolerance=float(tolerance), weight_relative_error=float(error),
                       missing_sigma=[name for name, v in zip(names, sigma) if np.isnan(v)])
    if progress:
        progress(1.0, "Averaging complete")
    return ExposureResult(names, mean, sigma, config, coverage, tail, integral, meta, diagnostics, warnings)


def average_run(run, config, progress=None):
    """Read only the original run and its stored provenance, never current GUI state."""
    meta = dict(getattr(run, "exposure_metadata", {}) or {})
    if not meta:
        # Legacy RunResult was produced by the fixed-matrix Evolution API.
        meta = dict(constant_conditions=True, legacy=True,
                    energy_eV=run.energy_eV, method=run.method)
    meta["targets"] = list(getattr(run, "targets", []) or [])
    sigma = getattr(run, "capture_sigma_barn", None)
    if sigma is None:
        sigma, sources = capture_snapshot(run.names, run.isotopes)
        meta["capture_sources"] = sources
        meta["sigma_provenance"] = "Saved isotope data"
    return average_exposure(run.names, run.times, run.trajectories, run.flux,
                            config, sigma, meta, progress)


def isotope_values(result, mode="absolute"):
    if mode == "relative":
        return result.mean / result.mean.sum()
    if mode == "sigma_n":
        return result.mean * result.sigma_barn
    if mode != "absolute":
        raise ValueError(f"Unknown abundance mode: {mode}")
    return result.mean.copy()


def mass_values(result, mode="absolute"):
    values = isotope_values(result, mode)
    grouped = {}
    for name, value in zip(result.names, values):
        a = mass_number(name)
        grouped[a] = grouped.get(a, 0.0) + float(value)
    aa = np.array(sorted(grouped), dtype=int)
    return aa, np.array([grouped[a] for a in aa])


def export_exposure(result, path):
    """Write both tables from the same aggregation used by the GUI."""
    delimiter = "," if str(path).lower().endswith(".csv") else "\t"
    header = dict(result.metadata, tau0_mbarn_inv=result.config.tau0,
                  exposure_coverage=result.coverage, uncovered_tail=result.tail,
                  numerical_weight_integral=result.weight_integral,
                  diagnostics=result.diagnostics, warnings=result.warnings)
    if getattr(result.config, "mode", "exposure") != "exposure":
        header.pop("tau0_mbarn_inv", None)
        header["distribution_coverage"] = header.pop("exposure_coverage")
    with open(path, "w", encoding="utf-8", newline="") as stream:
        stream.write("# " + json.dumps(header, ensure_ascii=True, allow_nan=False) + "\n")
        writer = csv.writer(stream, delimiter=delimiter)
        writer.writerow(["Isotope", "A", "Mean N", "Mean fraction", "Sigma N [barn * inventory]"])
        rel, sn = isotope_values(result, "relative"), isotope_values(result, "sigma_n")
        def fmt(v):
            return "NA" if not np.isfinite(v) else format(float(v), ".16g")
        for i, name in enumerate(result.names):
            writer.writerow([name, mass_number(name), fmt(result.mean[i]), fmt(rel[i]), fmt(sn[i])])
        writer.writerow([])
        writer.writerow(["A", "Mean N(A)", "Mean fraction(A)", "Sum sigma_i N_i [barn * inventory]"])
        aa, absolute = mass_values(result)
        _, relative = mass_values(result, "relative")
        _, sigma = mass_values(result, "sigma_n")
        for a, n, f, v in zip(aa, absolute, relative, sigma):
            writer.writerow([a, fmt(n), fmt(f), fmt(v)])
