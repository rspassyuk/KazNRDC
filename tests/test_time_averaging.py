"""Physical-time quadrature, compatibility and presentation checks."""
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[1] / "model"))
import json
import tempfile
import unittest
import numpy as np
from core.exposure import ExposureConfig, average_exposure, export_exposure

class TimeAveragingTests(unittest.TestCase):
    def calc(self, t, y, **kwargs):
        return average_exposure(["Fe-56"], t, np.asarray(y).reshape(-1, 1), 0,
                                ExposureConfig(**kwargs))

    def test_uniform_interpolated_endpoints(self):
        t = np.array([0., 1., 4., 10.])
        r = self.calc(t, 2 + 3*t, mode="time_uniform", start_s=2, end_s=7)
        self.assertAlmostEqual(r.mean[0], 15.5)
        self.assertAlmostEqual(r.weight_integral, 1)
        self.assertEqual(r.tail, 0)

    def test_constant_log_grid(self):
        t = np.r_[0., np.geomspace(1e-6, 1e9, 400)]
        r = self.calc(t, np.full_like(t, 7), mode="time_uniform")
        self.assertAlmostEqual(r.mean[0], 7)

    def test_decay_analytic(self):
        t = np.linspace(0, 10, 20001)
        r = self.calc(t, np.exp(-.2*t), mode="time_uniform", start_s=1, end_s=9)
        self.assertAlmostEqual(r.mean[0], (np.exp(-.2)-np.exp(-1.8))/1.6, delta=1e-8)
        r = self.calc(t, np.exp(-.2*t), mode="time_exponential", mean_time_s=2)
        self.assertAlmostEqual(r.mean[0], (1-np.exp(-7))/1.4, delta=1e-8)

    def test_exposure_equivalence(self):
        t = np.linspace(0, 10, 10001)
        y = np.exp(-t/3)
        r = self.calc(t, y, mode="time_exponential", mean_time_s=2)
        old = average_exposure(["Fe-56"], t, y[:, None], 2e26,
                              ExposureConfig(tau0=.4), metadata={"constant_conditions": True})
        np.testing.assert_allclose(r.mean, old.mean, rtol=1e-13)

    def test_tail_not_renormalized(self):
        t = np.linspace(0, 2, 10001)
        r = self.calc(t, np.ones_like(t), mode="time_exponential", mean_time_s=2)
        self.assertAlmostEqual(r.mean[0], 1-np.exp(-1), delta=1e-8)
        self.assertTrue(r.warnings)

    def test_invalid_parameters(self):
        for args in (dict(start_s=2, end_s=1), dict(start_s=-1), dict(end_s=11),
                     dict(end_s=float("nan")), dict(start_s=0, end_s=0)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                self.calc([0, 10], [1, 1], mode="time_uniform", **args)
        for value in (0, -1, float("inf"), float("nan")):
            with self.assertRaises(ValueError):
                self.calc([0, 10], [1, 1], mode="time_exponential", mean_time_s=value)

    def test_export_parameters(self):
        r = self.calc([0, 10], [1, 1], mode="time_uniform", start_s=2, end_s=8)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"data.txt"
            export_exposure(r, path)
            meta = json.loads(path.read_text().splitlines()[0][2:])
            self.assertEqual(meta["averaging_mode"], "time_uniform")
            self.assertEqual(meta["time_start_s"], 2)
            self.assertNotIn("tau0_mbarn_inv", meta)

    def test_panel_modes_and_session(self):
        from PyQt5.QtWidgets import QApplication
        from exposure_panel import ExposurePanel
        from core.session import save_session, load_session
        from test_exposure import make_run
        from model.modelgraph.exposure_plot import plot_exposure
        from matplotlib.figure import Figure
        app = QApplication.instance() or QApplication([])
        panel = ExposurePanel()
        self.assertEqual(panel.distribution.currentData(), "time_uniform")
        run = make_run()
        from core.exposure import average_run
        for mode in ("time_uniform", "time_exponential", "exposure"):
            run.exposure_result = average_run(run, ExposureConfig(mode=mode, mean_time_s=2))
            with tempfile.TemporaryDirectory() as folder:
                path = str(Path(folder)/"run.kaz")
                save_session(run, None, 0, path)
                loaded, _, _ = load_session(path)
                panel.set_run(loaded)
                self.assertEqual(panel.distribution.currentData(), mode)
                fig = Figure(); ax = fig.add_subplot()
                plot_exposure(ax, loaded.exposure_result)
                self.assertEqual(len(ax.lines), 1)
        panel.close()

if __name__ == "__main__":
    unittest.main()
