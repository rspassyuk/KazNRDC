"""Graph mode controls follow both application themes without changing values."""
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.append(str(ROOT / "model"))
import unittest
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QPalette
from model.modelgraph.model_graph import Left, PLOT_HEAT
from exposure_panel import ExposurePanel

class GraphThemeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_modes_in_both_themes(self):
        panel = Left()
        panel.set_plot_kind(PLOT_HEAT)
        panel.show()
        for theme, active, idle in (("dark", "#e8943a", "#2e2823"),
                                    ("light", "#2b7fc4", "#eaf0f6")):
            self.app.setStyleSheet((ROOT / "model" / (theme + ".css")).read_text())
            for buttons in (panel._heat_btns, panel._conc_btns, panel._eq_n_btns):
                for selected in buttons:
                    selected.click()
                    self.app.processEvents()
                    self.assertEqual(sum(b.isChecked() for b in buttons), 1)
                    for button in buttons:
                        self.assertTrue(button.property("graphMode"))
                        self.assertEqual(button.grab().toImage().pixelColor(8, button.height() // 2).name(),
                                         active if button.isChecked() else idle)
            panel._btn_mevs.click()
            self.assertEqual(panel.current_heat_unit(), "mev")
            panel._btn_wcm3.click()
            self.assertEqual(panel.current_heat_unit(), "wcm3")
        panel.close()
        self.app.setStyleSheet("")

    def test_averaging_labels(self):
        panel = ExposurePanel()
        labels = [panel.distribution.itemText(i) for i in range(panel.distribution.count())]
        self.assertIn("Exponential neutron exposure", labels)
        self.assertFalse(any("legacy" in label.lower() for label in labels))
        self.assertEqual(panel.distribution.itemData(2), "exposure")
        self.assertEqual(panel.tabs.count(), 2)
        panel.close()

if __name__ == "__main__":
    unittest.main()
