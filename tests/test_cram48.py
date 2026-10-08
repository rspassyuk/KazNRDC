"""Independent numerical checks for the order-48 solver."""
import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.linalg import expm
from scipy.stats import poisson
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core.cram import cram48_step
from core.burnup import BurnupMatrix


class Cram48Tests(unittest.TestCase):
    def test_long_chain(self):
        a = np.diag(-np.ones(180)) + np.diag(np.ones(179), -1)
        a[-1,-1] = 0
        n = np.eye(180)[0]
        expected = np.r_[poisson.pmf(np.arange(179),70),poisson.sf(178,70)]
        actual = cram48_step(a,n,70)
        # A long nonnormal chain amplifies single-step roundoff; require 1e-9 absolute accuracy.
        np.testing.assert_allclose(actual,expected,atol=1e-9,rtol=1e-9)
        self.assertAlmostEqual(actual.sum(),1,places=12)

    def test_independent_exponential_and_dispatch(self):
        for a in (np.array([[-1.,0],[1.,-2.]]), np.array([[-1e6,0,0],[1e6,-1,0],[0,1,0]])):
            n=np.eye(len(a))[0]
            m=BurnupMatrix.__new__(BurnupMatrix); m.matrix_A=a
            for t in (0.,1e-6,1.,100.):
                expected=expm(a*t)@n
                np.testing.assert_allclose(m.solve(n,t,'cram48'),expected,atol=2e-12,rtol=1e-9)
            np.testing.assert_array_equal(n,np.eye(len(a))[0])

    def test_zero_and_invalid(self):
        np.testing.assert_array_equal(cram48_step([[-1]],[1],0),[1])
        np.testing.assert_array_equal(cram48_step([[-1]],[0],10),[0])
        for t in (-1,np.nan,np.inf):
            with self.assertRaises(ValueError): cram48_step([[-1]],[1],t)


if __name__=='__main__': unittest.main()
