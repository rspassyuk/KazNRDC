"""Regression checks for isotope membership and independent channel filters."""
import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.append(str(ROOT / "model"))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from core.service import NuclearService
from core.entities import Isotope, DecayLink
from core.registry import IsotopeBuilder
from PyQt5.QtWidgets import QApplication
from model.modelresults.model_results import Left


class WorkspaceMembershipTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.svc = NuclearService()
        self.svc.provider = object()
        self.svc._build_energy_eV = 30000.0
        po = Isotope(84, "Po", 210)
        pb = Isotope(82, "Pb", 206)
        po.half_life_s = 10.0
        po.decay_constant = np.log(2) / 10.0
        po.add_decay(DecayLink("alpha", 1.0, 5.4, -2, -4, product=pb))
        self.svc.isotopes = [po, pb]

    def test_disabled_alpha_keeps_parent_and_zeroes_decay(self):
        matrix = self.svc.build_matrix(decay_mode_filter=["beta-"], mt_filter=[])
        self.assertIn("Po-210", matrix.names())
        np.testing.assert_allclose(matrix.matrix_decay, 0.0)
        with patch("model.modelresults.model_results.get_service", return_value=self.svc):
            widget = Left()
            widget.populate_unified_table(decay_mode_filter=["beta-"])
            names = [widget.unified_model.item(r, 0).text()
                     for r in range(widget.unified_model.rowCount())]
            self.assertIn("Po-210", names)
            row = names.index("Po-210")
            self.assertEqual(widget.unified_model.item(row, 3).text(), "disabled")
            self.assertEqual(widget.unified_model.item(row, 2).text(), "0")
            widget.close()

    def test_empty_filter_disables_every_decay(self):
        matrix = self.svc.build_matrix(decay_mode_filter=[])
        np.testing.assert_allclose(matrix.matrix_decay, 0.0)

    def test_alpha_enabled_has_parent_loss_and_daughter_gain(self):
        matrix = self.svc.build_matrix(decay_mode_filter=["alpha"])
        self.assertLess(matrix.matrix_decay[0, 0], 0)
        self.assertAlmostEqual(matrix.matrix_decay[1, 0], -matrix.matrix_decay[0, 0])

    def test_manual_membership_survives_energy_rebuild(self):
        def build(builder, specs, **kwargs):
            return [Isotope(*spec) for spec in specs]
        with patch.object(IsotopeBuilder, "build_range", autospec=True, side_effect=build):
            self.svc.edit_workspace_isotopes(["Bi-209"], ["Po-210"])
            self.assertEqual(set(self.svc.universe_labels()), {"Pb-206", "Bi-209"})
            self.svc.build_matrix(energy_eV=8000)
            self.assertEqual(set(self.svc.matrix.names()), {"Pb-206", "Bi-209"})
            self.svc.edit_workspace_isotopes(["Po210", "Po-210"])
            self.assertEqual(self.svc.universe_labels().count("Po-210"), 1)

    def test_invalid_label_does_not_change_membership(self):
        before = self.svc.universe_labels()
        with self.assertRaises(ValueError):
            self.svc.edit_workspace_isotopes(["not-an-isotope"])
        self.assertEqual(self.svc.universe_labels(), before)

    def test_cannot_remove_last_isotopes(self):
        with self.assertRaises(ValueError):
            self.svc.edit_workspace_isotopes(remove_labels=self.svc.universe_labels())

    def manual_link(self):
        return dict(parent="Po210", daughter="Pb206", mode="alpha",
                    half_life_s=10.0, branch=1.0)

    def test_manual_alpha_is_applied_once_and_solves_analytically(self):
        for filters in (["beta-"], ["alpha"], []):
            matrix = self.svc.build_matrix(
                decay_mode_filter=filters, mt_filter=[],
                manual_decay_links=[self.manual_link()])
            np.testing.assert_allclose(matrix.matrix_decay,
                                       [[-np.log(2)/10, 0], [np.log(2)/10, 0]])
            np.testing.assert_allclose(matrix.solve([1, 0], 10, "cram16"),
                                       [0.5, 0.5], atol=1e-12)
            self.assertGreater(matrix.heat_release(np.array([1., 0.])), 0)

    def test_manual_missing_endpoint_is_error(self):
        link = self.manual_link()
        link["daughter"] = "Pb-205"
        with self.assertRaisesRegex(ValueError, "requires isotopes"):
            self.svc.build_matrix(manual_decay_links=[link])

    def test_clearing_manual_channels_disables_alpha(self):
        self.svc.set_manual_decay_links([self.manual_link()])
        matrix = self.svc.build_matrix(decay_mode_filter=[], manual_decay_links=[])
        np.testing.assert_allclose(matrix.matrix_decay, 0)

    def test_adding_existing_isotope_does_not_reload_library(self):
        po = self.svc.isotopes[0]
        with patch.object(IsotopeBuilder, "build_range") as build:
            self.svc.edit_workspace_isotopes(["Po210"])
            build.assert_not_called()
        self.assertIs(self.svc.isotopes[0], po)

    def test_manual_channel_is_visible_in_results(self):
        link = dict(self.manual_link(), parent="Po-210", daughter="Pb-206")
        self.svc.set_manual_decay_links([link])
        self.svc.build_matrix(decay_mode_filter=["beta-"], mt_filter=[])
        with patch("model.modelresults.model_results.get_service", return_value=self.svc):
            widget = Left()
            widget.populate_unified_table(decay_mode_filter=["beta-"])
            row = next(r for r in range(widget.unified_model.rowCount())
                       if widget.unified_model.item(r, 0).text() == "Po-210")
            self.assertIn("alpha (manual)", widget.unified_model.item(row, 3).text())
            self.assertGreater(float(widget.unified_model.item(row, 2).text()), 0)
            widget.close()


if __name__ == "__main__":
    unittest.main()
