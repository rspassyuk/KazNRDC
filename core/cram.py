"""Sparse CRAM-16 with step doubling, signed-error checks and no clipping."""
from collections import OrderedDict
import numpy as np
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import splu

THETA = np.array([-10.843917078696988 + 19.27744616718165j, -5.264971343442647 + 16.22022147316793j, 5.948152268951177 + 3.587457362018322j, 3.509103608414918 + 8.436198985884374j, 6.416177699099435 + 1.1941223933701386j, 1.419375897185666 + 10.925363484496723j, 4.993174737717997 + 5.996881713603942j, -1.4139284624888862 + 13.497725698892745j], dtype=np.complex128)
ALPHA = np.array([-5.090152186522492e-07 - 2.422001765285229e-05j, 0.0002115174218246603 + 0.004389296964738067j, 113.3977517848393 + 101.94721704215857j, 15.059585270023467 - 5.751405277642182j, -64.50087802553965 - 224.59440762652096j, -1.4793007113557999 + 1.7686588323782937j, -62.51839246320792 - 11.190391094283228j, 0.04102313683541002 - 0.15743466173455467j], dtype=np.complex128)
ALPHA0 = 2.1248537104952236e-16

def _factor_step(matrix, step):
    identity = eye(matrix.shape[0], format="csc")
    return [splu(matrix * step - pole * identity) for pole in THETA]


def _apply_step(factors, vector):
    result = ALPHA0 * vector.copy()
    rhs = vector.astype(complex)
    for coefficient, factor in zip(ALPHA, factors):
        result += 2 * np.real(coefficient * factor.solve(rhs))
    return result


def cram16_step(matrix, vector, step):
    """One unchecked rational step, retained for diagnostics and regression tests."""
    return _apply_step(_factor_step(csc_matrix(matrix), step), np.asarray(vector, dtype=float))


def adaptive_cram16(matrix, initial, duration, *, rtol=1e-9, atol=1e-14,
                    max_attempts=4096, should_cancel=None):
    """Integrate a constant burnup matrix by accepted pairs of half steps.

    atol is relative to the initial absolute inventory; rtol is componentwise.
    Small signed roundoff is retained. Never clip or renormalize the output.
    """
    matrix = csc_matrix(matrix, dtype=float, copy=True)
    vector = np.array(initial, dtype=float, copy=True)
    duration = float(duration)
    if vector.ndim != 1 or matrix.shape != (len(vector), len(vector)):
        raise ValueError("Initial vector does not match the burnup matrix.")
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Evolution duration must be finite and non-negative.")
    if not np.isfinite(vector).all() or not np.isfinite(matrix.data).all():
        raise ValueError("Burnup matrix and initial vector must be finite.")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol <= 0 or atol <= 0:
        raise ValueError("CRAM tolerances must be finite and positive.")
    inventory = float(np.abs(vector).sum())
    if inventory == 0:
        return vector
    absolute = max(atol * inventory, np.finfo(float).tiny)
    if vector.min() < -absolute or vector.sum() < 0:
        raise ValueError("Initial inventory contains significant negative concentrations.")
    if duration == 0:
        return vector
    diagonal = matrix.diagonal()
    off = matrix.copy()
    off.setdiag(0)
    off.eliminate_zeros()
    if np.any(diagonal > 0) or np.any(off.data < 0):
        raise ValueError("Burnup matrix requires non-positive diagonal and non-negative transfers.")
    column_sum = np.asarray(matrix.sum(axis=0)).ravel()
    column_scale = np.asarray(abs(matrix).sum(axis=0)).ravel()
    column_tol = 1e-12 * column_scale
    conservative = np.all(abs(column_sum) <= column_tol)
    nonincreasing = np.all(column_sum <= column_tol)
    cache = OrderedDict()
    def step(v, h):
        if should_cancel and should_cancel():
            raise RuntimeError("Calculation cancelled.")
        if h not in cache:
            cache[h] = _factor_step(matrix, h)
            # At most 24 sparse factorizations; do not retain all rejected steps.
            while len(cache) > 3:
                cache.popitem(last=False)
        cache.move_to_end(h)
        return _apply_step(cache[h], v)

    remaining = duration
    h = duration
    attempts = 0
    while remaining > 0:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError("CRAM step refinement did not converge. "
                               "Check the network or use a verified independent solver.")
        h = min(h, remaining)
        if h <= 0 or h * .5 == 0 or remaining - h == remaining:
            raise RuntimeError("CRAM step size underflow during refinement.")
        coarse = step(vector, h)
        middle = step(vector, h * .5)
        fine = step(middle, h * .5)
        scale = absolute + rtol * np.maximum(abs(fine), abs(coarse))
        error = float(np.max(abs(fine - coarse) / scale))
        total = float(fine.sum())
        budget = 10 * (absolute * len(vector) + rtol * max(abs(vector.sum()), inventory))
        acceptable = (np.isfinite(fine).all() and np.isfinite(middle).all()
                      and np.isfinite(error) and error <= 1
                      and fine.min() >= -absolute and middle.min() >= -absolute
                      and total >= 0 and middle.sum() >= 0)
        if conservative:
            acceptable = acceptable and abs(total - vector.sum()) <= budget
        elif nonincreasing:
            acceptable = acceptable and total <= vector.sum() + budget
        if acceptable:
            vector = fine
            remaining = max(0., remaining - h)
            if error < 0.05:
                h *= 2
        else:
            h *= .5
    return vector
