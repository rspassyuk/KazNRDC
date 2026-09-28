"""Session serialisation — save/load RunResult + EqResult to a .kaz file (pickle)."""

import pickle

from core.service import RunResult, EqResult


def save_session(run: RunResult, eq, rho: float, path: str) -> None:
    """Pickle run + optional eq result + density to *path*."""
    payload = {
        "run":      run,
        "exposure_result": getattr(run, "exposure_result", None),
        "eq":       eq,       # EqResult or None
        "rho":      float(rho),
        "_version": 3,
    }
    with open(path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_session(path: str):
    """Return (RunResult, eq_or_None, rho) from a .kaz session file."""
    with open(path, "rb") as f:
        payload = pickle.load(f)
    if not isinstance(payload, dict) or "run" not in payload:
        raise ValueError("Not a valid .kaz session file.")
    run = payload["run"]
    eq  = payload.get("eq")
    rho = float(payload.get("rho", 0.0))
    if not isinstance(run, RunResult):
        raise ValueError("Session file does not contain a RunResult.")
    if "exposure_result" in payload:
        run.exposure_result = payload["exposure_result"]
    return run, eq, rho
