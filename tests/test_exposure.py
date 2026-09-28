"""Analytical and GUI regression tests; no external nuclear databases needed."""
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[1] / "model"))
import copy
import pickle
import tempfile
import time
import unittest
from unittest.mock import patch
import numpy as np
from core.exposure import *
from core.entities import Isotope, ReactionLink
from core.service import RunResult, NuclearService
from core.session import save_session, load_session

META = {"constant_conditions": True, "energy_eV": 30000}


def make_run():
    times = np.r_[0., np.geomspace(1e-6, 12, 1200)]
    iso = Isotope(26, "Fe", 56)
    iso.add_reaction(ReactionLink(102, "(n,g)", .2))
    return RunResult(["Fe-56"], times.tolist(), [np.array([1.]) for _ in times],
                     [0.] * len(times), isotopes=[iso], flux=1e27, energy_eV=30000,
                     exposure_metadata=META.copy(), capture_sigma_barn=np.array([.2]))


class ExposureTests(unittest.TestCase):
    def calc(self, times, y, tau0=1., **kwargs):
        return average_exposure(["Fe-56"], times, np.asarray(y).reshape(-1, 1),
                                1e27, ExposureConfig(tau0=tau0),
                                metadata=META, **kwargs)

    def test_constant_finite_interval(self):
        t = np.linspace(0, 2, 20001)
        r = self.calc(t, np.full_like(t, 3.))
        self.assertAlmostEqual(r.mean[0], 3 * (1-np.exp(-2)), delta=3e-9)
        self.assertAlmostEqual(r.coverage, 1-np.exp(-2))
        self.assertAlmostEqual(r.tail, np.exp(-2))

    def test_single_decay_analytic(self):
        t = np.linspace(0, 4, 30001)
        a, tau0 = .7, .6
        r = self.calc(t, 2*np.exp(-a*t), tau0)
        expected = 2*(1-np.exp(-(a+1/tau0)*t[-1]))/(1+a*tau0)
        self.assertAlmostEqual(r.mean[0], expected, delta=2e-8)

    def test_closed_parent_daughter(self):
        t = np.linspace(0, 8, 30001)
        y = np.c_[np.exp(-t), 1-np.exp(-t)]
        r = average_exposure(["H-3","He-3"], t, y, 1e27, ExposureConfig(1), metadata=META)
        self.assertAlmostEqual(r.mean.sum(), r.coverage, delta=1e-8)

    def test_log_grid_converges(self):
        errors = []
        for n in (100, 400, 1600):
            t = np.r_[0, np.geomspace(1e-6, 10, n)]
            r = self.calc(t, np.exp(-t))
            errors.append(abs(r.mean[0] - (1-np.exp(-20))/2))
        self.assertLess(errors[1], errors[0]/10)
        self.assertLess(errors[2], errors[1]/10)

    def test_coarse_weight_warning(self):
        r = self.calc([0, 10], [1, 1])
        self.assertGreater(r.diagnostics["weight_relative_error"], 1)
        self.assertTrue(any("coarse" in w for w in r.warnings))

    def test_tail_not_renormalized(self):
        t = np.linspace(0, .1, 1001)
        r = self.calc(t, np.ones_like(t))
        self.assertLess(r.mean[0], .1)
        self.assertTrue(any("tail" in w for w in r.warnings))

    def test_mass_and_sigma_branching(self):
        fe, co = Isotope(26,"Fe",60), Isotope(27,"Co",60,1)
        fe.add_reaction(ReactionLink(102,"capture",.3,q_yield=.3))
        fe.add_reaction(ReactionLink(102,"capture",.7,q_yield=.7,lfs=1))
        sig, _ = capture_snapshot([fe.name,co.name], [fe,co])
        self.assertEqual(sig[0], 1.)
        self.assertTrue(np.isnan(sig[1]))
        r = ExposureResult([fe.name,co.name],np.array([2.,3.]),np.array([1.,2.]),
                           ExposureConfig(),1,0,1)
        np.testing.assert_array_equal(mass_values(r)[0],[60])
        np.testing.assert_array_equal(mass_values(r)[1],[5])
        np.testing.assert_array_equal(mass_values(r,"sigma_n")[1],[8])
        np.testing.assert_array_equal(mass_values(r,"relative")[1],[1])

    def test_targets_do_not_filter_input(self):
        run=make_run()
        first=average_run(run,ExposureConfig())
        run.targets=["not-in-network"]
        second=average_run(run,ExposureConfig())
        np.testing.assert_array_equal(first.mean,second.mean)
        np.testing.assert_array_equal(run.trajectories,np.ones((len(run.times),1)))

    def test_bad_inputs(self):
        for flux,tau0,t,y in [
            (0,1,[0,1],[[1],[1]]), (np.inf,1,[0,1],[[1],[1]]),
            (1e27,0,[0,1],[[1],[1]]), (1e27,np.nan,[0,1],[[1],[1]]),
            (1e27,1,[0,0],[[1],[1]]), (1e27,1,[1,2],[[1],[1]]),
            (1e27,1,[0,1],[[1],[np.nan]]), (1e27,1,[0,1],[[1],[np.inf]]),
            (1e27,1,[0,np.inf],[[1],[1]]), (1e27,1,[0,1],[[1,2],[1,2]]),
            (1e27,1,[0,1],[[1],[-1e-4]])]:
            with self.subTest(flux=flux,tau0=tau0,t=t,y=y):
                with self.assertRaises(ValueError):
                    average_exposure(["Fe-56"],t,y,flux,ExposureConfig(tau0),metadata=META)

    def test_variable_conditions_rejected(self):
        with self.assertRaises(ValueError):
            average_exposure(["Fe-56"],[0,1],[[1],[1]],1e27,ExposureConfig(),
                             metadata={"constant_conditions":False})

    def test_small_negative_retained_and_reported(self):
        r=average_exposure(["H-3","He-3"],[0,1],[[1,0],[1,-1e-13]],1e27,
                          ExposureConfig(),metadata=META)
        self.assertLess(r.mean[1],0)
        self.assertEqual(r.diagnostics["negative_count"],1)

    def test_missing_sigma_is_not_zero(self):
        r=self.calc([0,1],[1,1])
        self.assertTrue(np.isnan(isotope_values(r,"sigma_n")[0]))
        self.assertEqual(r.diagnostics["missing_sigma"],["Fe-56"])

    def test_negative_total_blocked(self):
        with self.assertRaises(ValueError):
            self.calc([0,1],[1,-1e-16])

    def test_session_roundtrip_and_legacy(self):
        run=make_run()
        run.exposure_result=average_run(run,ExposureConfig(.4))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"new.kaz"
            save_session(run,None,2,p)
            loaded,_,rho=load_session(p)
            self.assertEqual(rho,2)
            np.testing.assert_array_equal(loaded.exposure_result.mean,run.exposure_result.mean)
            legacy=copy.deepcopy(run)
            for attr in ("exposure_result","capture_sigma_barn","exposure_metadata"):
                delattr(legacy,attr)
            with open(p,"wb") as f:pickle.dump({"run":legacy,"_version":1},f)
            loaded,_,rho=load_session(p)
            self.assertIsNone(getattr(loaded,"exposure_result",None))
            result=average_run(loaded,ExposureConfig())
            self.assertTrue(result.metadata["legacy"])

    def test_export_matches_aggregation(self):
        r=average_run(make_run(),ExposureConfig())
        with tempfile.TemporaryDirectory() as d:
            for ext in ("txt","csv"):
                p=Path(d)/("average."+ext);export_exposure(r,p)
                text=p.read_text()
                self.assertIn("tau0_mbarn_inv",text)
                self.assertIn("Mean N(A)",text)
                self.assertIn(format(r.mean[0],".16g"),text)

    def test_service_run_captures_provenance(self):
        from types import SimpleNamespace
        svc = NuclearService()
        svc.provider = SimpleNamespace(get_macs_priority=lambda: ["user"])
        svc._build_energy_eV = 30000.
        svc.isotopes = [Isotope(26, "Fe", 56)]
        run = svc.run_evolution({"Fe-56": 1.}, 2., flux=1e27,
                                n_points=200, energy_eV=30000., mt_filter=[])
        self.assertTrue(run.exposure_metadata["constant_conditions"])
        self.assertEqual(run.exposure_metadata["mt_filter"], [])
        self.assertTrue(np.isnan(run.capture_sigma_barn[0]))
        svc.isotopes.clear()
        result = average_run(run, ExposureConfig(1.))
        self.assertEqual(result.metadata["energy_eV"], 30000.)
        self.assertEqual(result.names, ["Fe-56"])

    def test_provenance_uses_snapshot(self):
        run=make_run()
        run.isotopes[0]._reactions[0].sigma_barn=99.
        r=average_run(run,ExposureConfig())
        self.assertEqual(r.sigma_barn[0],.2)


class ExposureGuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from PyQt5.QtWidgets import QApplication
        cls.app=QApplication.instance() or QApplication([])

    def test_worker_table_graph_export_and_saved_viewer(self):
        from exposure_panel import ExposurePanel
        from model.modelresults.model_results import ResultsViewerDialog
        panel=ExposurePanel();run=make_run();panel.set_run(run)
        panel.start()
        deadline=time.monotonic()+15
        while panel.busy and time.monotonic()<deadline:
            self.app.processEvents();time.sleep(.005)
        self.assertFalse(panel.busy)
        self.assertIsNotNone(panel.result,panel.status.text())
        self.assertEqual(panel.tables[0].model().rowCount(),1)
        self.assertEqual(panel.tabs.count(),2)
        self.assertFalse(hasattr(panel,"canvas"))
        with tempfile.TemporaryDirectory() as d:
            target=str(Path(d)/"average.csv")
            with patch("exposure_panel.QFileDialog.getSaveFileName",return_value=(target,"")):
                panel.export()
            self.assertTrue(Path(target).exists())
            path=str(Path(d)/"average.kaz");save_session(run,None,0,path)
            loaded,eq,rho=load_session(path)
            viewer=ResultsViewerDialog(loaded,eq,rho,path)
            self.assertIsNotNone(viewer.exposure_panel.result)
            viewer._graph_combo.setCurrentText("Abundance")
            viewer._abund_source.setCurrentIndex(1)
            viewer._graph_canvas.draw()
            self.assertIn("Time average",viewer._graph_fig.axes[0].get_title())
            viewer.close()
        panel.close()

    def test_main_graph_uses_separate_result(self):
        import model.modelgraph.model_graph as graph
        run=make_run();run.exposure_result=average_run(run,ExposureConfig())
        svc=NuclearService();svc.last_run=run
        with patch.object(graph,"get_service",return_value=svc):
            left=graph.Left()
            left.plot(graph.PLOT_ABUND,"Dark",exposure_averaged=True,abund_mode="sigma_n")
            self.assertIn("Exposure-averaged",left.ax.get_title())
            np.testing.assert_allclose(left.ax.lines[0].get_ydata(),
                                       mass_values(run.exposure_result,"sigma_n")[1])
            for theme in ("Dark", "Light"):
                for mode in ("absolute", "relative", "sigma_n"):
                    left.plot(graph.PLOT_ABUND,theme,exposure_averaged=True,abund_mode=mode)
                    left.canvas.draw()
                    self.assertEqual(len(left.ax.lines),1)
            run.exposure_result.sigma_barn[:] = np.nan
            left.plot(graph.PLOT_ABUND,"Dark",exposure_averaged=True,abund_mode="sigma_n")
            left.canvas.draw()
            self.assertTrue(any("unavailable" in t.get_text() for t in left.ax.texts))
            left.close()


if __name__ == "__main__":
    unittest.main()
