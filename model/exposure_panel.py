"""Reusable exposure-averaging panel for live and saved Evolution results."""
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QLabel, QLineEdit, QComboBox, QPushButton, QTabWidget, QTableView,
    QFileDialog, QMessageBox, QHeaderView)
import numpy as np
from core.exposure import (ExposureConfig, average_run, mass_values, isotope_values,
                           mass_number, export_exposure)
from GUIcomp import run_with_progress




class ExposurePanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.run = None
        self.result = None
        self.busy = False
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.distribution = QComboBox()
        self.distribution.addItem("Time average over interval", "time_uniform")
        self.distribution.addItem("Exponential irradiation durations", "time_exponential")
        self.distribution.addItem("Exponential neutron exposure", "exposure")
        self.tau_edit = QLineEdit("0.3")
        form.addRow("Distribution", self.distribution)
        form.addRow("Mean exposure tau0 [mbarn^-1]", self.tau_edit)
        self.form = form
        self.start_edit = QLineEdit("0")
        self.end_edit = QLineEdit()
        self.end_edit.setPlaceholderText("Final saved time")
        self.duration_edit = QLineEdit("1e9")
        form.addRow("Start time [s]", self.start_edit)
        form.addRow("End time [s]", self.end_edit)
        form.addRow("Mean duration t0 [s]", self.duration_edit)
        self.distribution.currentIndexChanged.connect(self._mode_changed)
        self._mode_changed()
        layout.addLayout(form)
        buttons = QHBoxLayout()
        self.average_button = QPushButton("Average results")
        self.export_button = QPushButton("Export averaged results")
        self.save_button = QPushButton("Save averaged session")
        self.average_button.setObjectName("exposure_average_pb")
        self.export_button.setObjectName("exposure_export_pb")
        self.save_button.setObjectName("exposure_save_pb")
        for button in (self.average_button, self.export_button, self.save_button):
            buttons.addWidget(button)
        layout.addLayout(buttons)
        self.status = QLabel("Run Evolution or open a saved session.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.PlainText)
        layout.addWidget(self.status)
        self.tabs = QTabWidget()
        self.tables = []
        for title in ("Isotopes", "Mass numbers"):
            table = QTableView()
            table.setModel(QStandardItemModel(table))
            table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tables.append(table)
            self.tabs.addTab(table, title)
        hint = QLabel("Plot in Graph: Abundance vs mass number A > Averaged results.\n"
                      "Saved sessions: Graph > Abundance > Averaged results.")
        hint.setWordWrap(True)
        layout.addWidget(hint)
        layout.addWidget(self.tabs)
        self.average_button.clicked.connect(self.start)
        self.export_button.clicked.connect(self.export)
        self.save_button.clicked.connect(self.save)
        self.export_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.average_button.setEnabled(False)

    def _mode_changed(self):
        mode = self.distribution.currentData()
        for widget, visible in ((self.tau_edit, mode == "exposure"),
                                (self.start_edit, mode == "time_uniform"),
                                (self.end_edit, mode == "time_uniform"),
                                (self.duration_edit, mode == "time_exponential")):
            widget.setVisible(visible)
            self.form.labelForField(widget).setVisible(visible)

    def set_run(self, run):
        self.run = run
        self.result = getattr(run, "exposure_result", None) if run else None
        self.average_button.setEnabled(run is not None and not self.busy)
        self.export_button.setEnabled(self.result is not None and not self.busy)
        self.save_button.setEnabled(self.result is not None and not self.busy)
        for table in self.tables:
            table.model().clear()
        if self.result:
            cfg = self.result.config
            self.tau_edit.setText(str(cfg.tau0))
            self.distribution.setCurrentIndex(self.distribution.findData(getattr(cfg, "mode", "exposure")))
            self.start_edit.setText(str(getattr(cfg, "start_s", 0)))
            end = getattr(cfg, "end_s", None)
            self.end_edit.setText("" if end is None else str(end))
            self.duration_edit.setText(str(getattr(cfg, "mean_time_s", 1e9)))
            self.render()
        else:
            self.status.setText("Ready to average the complete computed vector." if run else
                                "Run Evolution or open a saved session.")

    def start(self):
        if self.busy or self.run is None:
            return
        try:
            mode = self.distribution.currentData()
            params = dict(mode=mode)
            if mode == "exposure":
                params["tau0"] = float(self.tau_edit.text())
            elif mode == "time_uniform":
                params["start_s"] = float(self.start_edit.text())
                params["end_s"] = float(self.end_edit.text()) if self.end_edit.text().strip() else None
            else:
                params["mean_time_s"] = float(self.duration_edit.text())
            config = ExposureConfig(**params)
        except ValueError:
            QMessageBox.warning(self, "Trajectory averaging", "Enter valid numeric averaging parameters in seconds (neutron exposure in mbarn^-1).")
            return
        run = self.run
        self.busy = True
        for control in (self.average_button, self.export_button, self.save_button, self.tau_edit, self.distribution, self.start_edit, self.end_edit, self.duration_edit):
            control.setEnabled(False)
        self.status.setText("Averaging the stored trajectory...")
        def finish():
            self.busy = False
            self.average_button.setEnabled(self.run is not None)
            for control in (self.tau_edit, self.distribution, self.start_edit, self.end_edit, self.duration_edit):
                control.setEnabled(True)
            self.export_button.setEnabled(self.result is not None)
            self.save_button.setEnabled(self.result is not None)
        def done(result):
            run.exposure_result = result
            if self.run is run:
                self.result = result
                self.render()
            finish()
        def failed(message):
            if self.run is run:
                self.result = None
                run.exposure_result = None
                for table in self.tables:
                    table.model().clear()
                details = str(message).strip().splitlines()
                summary = details[-1] if details else "Unknown averaging error."
                self.status.setText("Averaging blocked:\n" + summary)
                self.status.setToolTip(str(message))
            finish()
        self.job = run_with_progress(self.window(), "Trajectory averaging",
                                     average_run, run, config,
                                     on_done=done, on_error=failed, cancelable=False)

    def render(self):
        result = self.result
        d = result.diagnostics
        text = (f"Distribution coverage: {result.coverage:.6%}; Uncovered tail: {result.tail:.6%}\n"
                f"Weight integral: {result.weight_integral:.8g}; "
                f"weight relative error: {d['weight_relative_error']:.3g}\n"
                f"Inventory: {d['initial_total']:.6g} -> {d['final_total']:.6g}; "
                f"range [{d['minimum_total']:.6g}, {d['maximum_total']:.6g}]; "
                f"minimum N: {d['minimum_concentration']:.3g}; negatives: {d['negative_count']}\n"
                f"Normalization: {result.metadata.get('normalization', 'None (finite-range integral)')}. Relative denominator: total mean N.")
        boundary = result.metadata.get("boundary_columns")
        if boundary is not None:
            text += f"\nNet-loss matrix columns: {len(boundary)}"
            if boundary:
                text += " (" + ", ".join(boundary[:8]) + (", ..." if len(boundary) > 8 else "") + ")"
        else:
            text += "\nBoundary metadata: unavailable for this saved run."
        text += "\nCross-section mode: " + str(result.metadata.get("sigma_mode", "unknown"))
        self.status.setText(text + ("\n" + "\n".join(result.warnings) if result.warnings else ""))
        rel = isotope_values(result, "relative")
        sn = isotope_values(result, "sigma_n")
        isotope_rows = [[name, mass_number(name), result.mean[i], rel[i], sn[i]]
                        for i, name in enumerate(result.names)]
        aa, nn = mass_values(result)
        _, rr = mass_values(result, "relative")
        _, ss = mass_values(result, "sigma_n")
        self._fill(self.tables[0], ["Isotope", "A", "Mean N", "Mean fraction", "Sigma N"],
                   isotope_rows)
        self._fill(self.tables[1], ["A", "Mean N(A)", "Mean fraction(A)", "Sum sigma_i N_i"],
                   zip(aa, nn, rr, ss))
        self.export_button.setEnabled(not self.busy)
        self.save_button.setEnabled(not self.busy)

    @staticmethod
    def _fill(table, headers, rows):
        model = table.model()
        model.clear()
        model.setHorizontalHeaderLabels(headers)
        for row in rows:
            items = []
            for v in row:
                text = str(v) if isinstance(v, (str, int, np.integer)) else (
                    f"{v:.8g}" if np.isfinite(v) else "NA")
                item = QStandardItem(text)
                item.setEditable(False)
                items.append(item)
            model.appendRow(items)

    def export(self):
        if self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export averaged results", "exposure.csv",
                                             "CSV (*.csv);;Text (*.txt)")
        if path:
            try:
                export_exposure(self.result, path)
            except Exception as exc:
                QMessageBox.warning(self, "Export failed", str(exc))

    def save(self):
        if self.run is None or self.result is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save averaged session", "exposure.kaz",
                                             "NuMatRx session (*.kaz)")
        if path:
            try:
                from core.session import save_session
                save_session(self.run, getattr(self, "session_eq", None),
                             getattr(self, "session_rho", 0.0), path)
            except Exception as exc:
                QMessageBox.warning(self, "Save failed", str(exc))
