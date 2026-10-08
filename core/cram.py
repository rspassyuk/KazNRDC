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


# IPF CRAM-48 coefficients: Pusa (2016), doi:10.13182/NSE15-26.
# Numerical table checked against openmc/deplete/cram.py (OpenMC).
# IPF residues require sequential updates, not the CRAM-16 PFD sum.
THETA48 = np.array([-44.65731934165702, -5.284616241568964, -8.867715667624458, 3.493013124279215, 15.64102508858634, 17.42097597385893, -28.34466755180654, 16.61569367939544, 8.011836167974721, -2.056267541998229, 14.49208170441839, 18.53807176907916, 9.932562704505182, -22.44223871767187, 0.8590014121680897, -12.86192925744479, 11.64596909542055, 18.06076684783089, 5.870672154659249, -35.42938819659747, 19.0132348906025, 18.85508331552577, -17.34689708174982, 13.1628423712519]) + 1j * np.array([62.33225190695437, 40.57499381311059, 43.25515754166724, 32.81615453173585, 15.58061616372237, 10.7662930571442, 54.92841024648724, 13.16994930024688, 27.8023211130941, 37.94824788914354, 17.99988210051809, 5.974332563100539, 25.32823409972962, 51.79633600312162, 35.3645619429435, 46.00304902833652, 22.87153304140217, 8.36820058009982, 30.29700159040121, 58.34381701800013, 1.194282058271408, 3.583428564427879, 48.83941101108207, 20.42951874827759])
ALPHA48 = np.array([638.7380733878774, 190.989617906573, 423.6195226571914, 464.5770595258726, 776.5163276752432, 1907.115136768522, 2909.892685603256, 194.477220662045, 138279.9786972332, 5628.442079602433, 215.168128379422, 1324.72024051442, 16175.48476343347, 111.2729040439685, 107.4624783191125, 88.35727765158191, 93.54078136054179, 94.18142823531574, 104.0012390717851, 68.61882624343235, 87.66654491283722, 105.600761938965, 77.38987569039419, 104.1366366475571]) + 1j * np.array([-674.3912502859256, -397.3203432721332, -2041.233768918671, -1652.917287299683, -17836.17639907328, -58870.68595142284, -9953.25534551456, -1427.131226068449, -3256885.197214938, -29242.84515884309, -1121.774011188224, -63700.88443140973, -1008798.413156542, -88.37109731680418, -145.724611640818, -63.8828618841936, -219.5424319460237, -671.9055740098034, -169.3747595553868, -11.77598523430493, -4596.464999363902, -1738.294585524067, -43.11715386228984, -277.7743732451969])

ALPHA048 = 2.258038182743983e-47


def cram48_step(matrix, initial, duration):
    """Single-step order-48 IPF CRAM; no adaptation, clipping or normalization."""
    matrix = csc_matrix(matrix, dtype=float, copy=True)
    state = np.array(initial, dtype=float, copy=True)
    duration = float(duration)
    if state.ndim != 1 or matrix.shape != (state.size, state.size):
        raise ValueError("Initial vector does not match the burnup matrix.")
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Evolution duration must be finite and non-negative.")
    if not np.isfinite(state).all() or not np.isfinite(matrix.data).all():
        raise ValueError("Burnup matrix and initial vector must be finite.")
    if duration == 0 or not np.any(state):
        return state
    scaled = matrix * duration
    identity = eye(state.size, format="csc")
    for pole, residue in zip(THETA48, ALPHA48):
        factor = splu(scaled - pole * identity)
        state += 2.0 * np.real(residue * factor.solve(state.astype(complex)))
    state *= ALPHA048
    if not np.isfinite(state).all():
        raise RuntimeError("CRAM-48 produced nonfinite concentrations.")
    return state
