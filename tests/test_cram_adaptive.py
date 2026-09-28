"""Regression tests for adaptive CRAM on long, positive burnup chains."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from scipy.sparse import diags, csc_matrix
from scipy.stats import poisson
from scipy.linalg import expm
from core.cram import adaptive_cram16, cram16_step
from core.burnup import BurnupConfig, BurnupMatrix
from core.entities import Isotope
from core.exposure import average_exposure, ExposureConfig


def chain(count=180):
    matrix = diags([-np.ones(count), np.ones(count-1)], [0, -1], format="csc")
    matrix[-1, -1] = 0
    matrix.eliminate_zeros()
    initial = np.zeros(count)
    initial[0] = 1
    return matrix, initial


class AdaptiveCramTests(unittest.TestCase):
    def test_long_chain_against_poisson(self):
        matrix, initial = chain()
        t = 70.
        reference = np.r_[poisson.pmf(np.arange(179), t), poisson.sf(178, t)]
        raw = cram16_step(matrix, initial, t)
        self.assertLess(raw.min(), -1e-10, "Fixture must expose the one-step regression")
        result = adaptive_cram16(matrix, initial, t)
        self.assertGreaterEqual(result.min(), -1e-14)
        np.testing.assert_allclose(result, reference, atol=2e-12, rtol=1e-8)
        self.assertAlmostEqual(result.sum(), 1., delta=1e-11)

    def test_closed_stiff_chain(self):
        matrix = csc_matrix([[-1e6, 0, 0], [1e6, -1, 0], [0, 1, 0]])
        result = adaptive_cram16(matrix, [1,0,0], 2.)
        b = 1e6/(1e6-1)*np.exp(-2.)
        np.testing.assert_allclose(result,[0,b,1-b],atol=1e-11,rtol=1e-9)

    def test_open_chain_loses_inventory_without_normalization(self):
        matrix = csc_matrix([[-1.,0],[1.,-2.]])
        expected = expm(matrix.toarray()*3) @ np.array([1.,0.])
        actual = adaptive_cram16(matrix,[1.,0.],3)
        np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-9)
        self.assertLess(actual.sum(),.1)

    def test_zero_time_and_zero_inventory(self):
        matrix, initial = chain(3)
        np.testing.assert_array_equal(adaptive_cram16(matrix,initial,0),initial)
        np.testing.assert_array_equal(adaptive_cram16(matrix,[0,0,0],10),[0,0,0])

    def test_display_grid_independent_at_same_final_time(self):
        matrix, initial = chain()
        direct = adaptive_cram16(matrix,initial,70)
        state = initial.copy()
        for _ in range(10):
            state = adaptive_cram16(matrix,state,7)
        np.testing.assert_allclose(direct,state,atol=2e-12,rtol=1e-8)

    def test_invalid_inputs_and_matrix(self):
        for matrix,initial,time in [
            ([[-1]],[1],-1), ([[-1]],[np.nan],1),
            ([[1]],[1],1), ([[-1,-1],[0,-1]],[1,0],1),
            ([[-1]],[-1],1), ([[-1]],[1],np.inf)]:
            with self.subTest(matrix=matrix,time=time):
                with self.assertRaises(ValueError):
                    adaptive_cram16(matrix,initial,time)

    def test_failure_is_explicit_not_clipped(self):
        matrix, initial=chain()
        with self.assertRaisesRegex(RuntimeError,"did not converge"):
            adaptive_cram16(matrix,initial,70,max_attempts=1)

    def test_cancellation_inside_refinement(self):
        matrix, initial=chain()
        with self.assertRaisesRegex(RuntimeError,"cancelled"):
            adaptive_cram16(matrix,initial,70,should_cancel=lambda:True)

    def test_public_solver_uses_adaptation(self):
        matrix, initial=chain()
        mt=BurnupMatrix([Isotope(26,"Fe",a) for a in range(1,181)],BurnupConfig())
        mt.matrix_A=matrix.toarray()
        actual=mt.solve(initial,70,"cram16")
        self.assertGreaterEqual(actual.min(),-1e-14)
        self.assertGreater(actual[100],0)
        np.testing.assert_array_equal(mt.matrix_A,matrix.toarray())

    def test_adaptive_trajectory_can_be_averaged(self):
        matrix, initial=chain(120)
        times=np.r_[0.,np.geomspace(1e-3,150,120)]
        trajectory=[adaptive_cram16(matrix,initial,t) for t in times]
        result=average_exposure([f"Fe-{a}" for a in range(1,121)],times,trajectory,1e27,
                               ExposureConfig(20),metadata={"constant_conditions":True})
        self.assertGreater(result.mean.sum(),.99)
        self.assertLess(result.diagnostics["minimum_concentration"],1e-12)


if __name__ == "__main__":
    unittest.main()
