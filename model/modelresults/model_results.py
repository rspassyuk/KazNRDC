"""
Results model: Build/filter burnup matrix → evolve or equilibrium scan → display.

Right panel (controls):
  Shared: IC table + Load Matrix + half-life filter + MT filter checkboxes
  Tab "Isotopes & Matrix": flux, open/closed channels, refresh
  Tab "Evolution"        : method, time start/end+units, flux, energy, density, display
  Tab "Equilibrium"      : auto-plateau OR manual flux/time, energy, density

Left panel (outputs):
  Tab "Isotopes & Matrix": Block-1 isotope cross-section table + Block-2 rate matrix
  Tab "Evolution"        : time-step slider, results table (N, fraction, Q), total Q
  Tab "Equilibrium"      : table (flux, time, Q_MeV, Q_Wcm3)
"""

import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt5.QtCore import Qt, QObject, QSortFilterProxyModel
from PyQt5.QtGui import QKeySequence, QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDialog, QFileDialog,
    QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QMessageBox, QPushButton, QRadioButton,
    QScrollArea, QShortcut, QSizePolicy, QSlider, QSplitter, QTabWidget,
    QTableView, QVBoxLayout, QWidget,
)

from core.service import get_service, EqResult, RunResult
from exposure_panel import ExposurePanel
from core.sensitivity import SensitivityConfig
from GUIcomp import run_with_progress


# ---------------------------------------------------------------------------
# Time-unit constants
# ---------------------------------------------------------------------------

_TIME_UNITS: List[Tuple[str, float]] = [
    ("s",   1.0),
    ("h",   3600.0),
    ("d",   86400.0),
    ("y",   365.25 * 86400.0),
]

_HL_PRESETS = ["all", "yearly", "daily", "hourly", "minute", "second"]

_MT_LABELS: List[Tuple[str, int]] = [
    ("(n,γ)", 102), ("(n,α)", 107), ("(n,p)", 103), ("(n,2n)", 16), ("(n,f)", 18),
]

_MT_NAMES = {102: "(n,γ)", 107: "(n,α)", 103: "(n,p)", 16: "(n,2n)", 18: "(n,f)"}

_DECAY_LABELS: List[Tuple[str, str]] = [
    ("β⁻", "beta-"), ("β⁺", "beta+"), ("α", "alpha"), ("EC", "ec"), ("IT", "it"), ("SF", "sf"),
]


def _unit_seconds(name: str) -> float:
    for u, f in _TIME_UNITS:
        if u == name:
            return f
    return 1.0


def _fmt_e(v) -> str:
    return "—" if v is None else f"{v:.5e}"


def _human_time(t: float) -> str:
    for div, unit in [(3.15576e7, "yr"), (86400.0, "d"), (3600.0, "h"),
                      (60.0, "min"), (1.0, "s"), (1e-3, "ms")]:
        if t >= div:
            return f"{t / div:.3g} {unit}"
    return f"{t:.3g} s"


def _parse_A(name: str):
    """Return integer mass number from 'El-NNN' or 'El-NNNm1'; None if unparseable."""
    try:
        core = name.split("-", 1)[1]
        digits = ""
        for ch in core:
            if ch.isdigit():
                digits += ch
            else:
                break
        return int(digits) if digits else None
    except (IndexError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Model container
# ---------------------------------------------------------------------------

class model_results(QObject):
    """Pages are loaded by core.py.  Exposes Left widget + Right layout."""

    def __init__(self):
        super().__init__()
        self.__right = Right()
        self.__left = Left()

        self.__right.calculate_pb.clicked.connect(self.on_calculate)
        self.__right.calc_eq_pb.clicked.connect(self.on_calculate_equilibrium)
        self.__right.sens_calc_pb.clicked.connect(self.on_calculate_sensitivity)
        self.__left.csv_pb.clicked.connect(self.on_export_csv)
        self.__left.txt_pb.clicked.connect(self.on_export_txt)
        self.__left.eq_csv_pb.clicked.connect(self.on_export_eq_csv)
        self.__left.save_session_pb.clicked.connect(self.on_save_session)
        self.__left.open_session_pb.clicked.connect(self.on_open_session)
        self.__right.refresh_matrix_pb.clicked.connect(self._on_refresh_matrix)
        self.__left.iso_add_pb.clicked.connect(self._add_workspace_isotopes)
        self.__left.iso_add_edit.returnPressed.connect(self._add_workspace_isotopes)
        self.__left.iso_remove_pb.clicked.connect(self._remove_workspace_isotopes)


        # sync right tab <-> left tab
        self.__right.mode_tabs.currentChanged.connect(
            self.__left.calc_tabs.setCurrentIndex)
        self.__left.calc_tabs.currentChanged.connect(
            lambda index: self.__right.mode_tabs.setCurrentIndex(index)
            if index < self.__right.mode_tabs.count() else None)

        # auto-refresh matrix view when matrix is built in IsotopeChart
        svc = get_service()
        svc.register_matrix_built_callback(self._on_matrix_built_external)
        svc.register_cycles_found_callback(self._on_cycles_found)

    def getLeftWidget(self):
        return self.__left

    def getRightLayout(self):
        return self.__right

    # -- matrix built callback -----------------------------------------------

    def _on_matrix_built_external(self):
        """Called when IsotopeChart fires _fire_matrix_built."""
        right = self.__right
        right._load_from_matrix()
        self._on_refresh_matrix()
        # Switch to Isotopes & Matrix tab (index 0)
        self.__left.calc_tabs.setCurrentIndex(0)
        self.__right.mode_tabs.setCurrentIndex(0)

    def _on_cycles_found(self, cycles: list) -> None:
        self.__left.populate_cycles(cycles)
        self.__left.calc_tabs.setCurrentIndex(0)

    def _edit_workspace_isotopes(self, additions=(), removals=()):
        self._on_refresh_matrix(additions=additions, removals=removals)

    def _add_workspace_isotopes(self):
        labels = [v for v in re.split(r"[\s,;]+", self.__left.iso_add_edit.text().strip()) if v]
        if labels:
            self._edit_workspace_isotopes(additions=labels)

    def _remove_workspace_isotopes(self):
        left = self.__left
        labels = [left.unified_proxy.index(index.row(), 0).data()
                  for index in left.unified_table.selectionModel().selectedRows()]
        if labels:
            self._edit_workspace_isotopes(removals=labels)

    def _on_refresh_matrix(self, checked=False, additions=(), removals=()):
        if getattr(self, "_matrix_busy", False):
            return
        svc, right, left = get_service(), self.__right, self.__left
        if not svc.initialized:
            return
        try:
            flux = float(right.matrix_flux_edit.text().strip() or "0")
            reactor_eV = float(right.reactor_energy_edit.text().strip() or "30000")
            astro_keV = float(right.astro_energy_edit.text().strip() or "30")
            manual_links = left.read_manual_decay_links()
        except ValueError as exc:
            QMessageBox.warning(left, "Matrix settings", str(exc))
            return
        open_mts = right.get_matrix_open_mts()
        decay_filter = right.get_decay_mode_filter()
        hl_preset = right.hl_filter_box.currentText()
        mt_filter = right.get_mt_filter()
        initial = {right.ic_model.item(row, 0).text(): right.ic_model.item(row, 1).text()
                   for row in range(right.ic_model.rowCount())
                   if right.ic_model.item(row, 0) and right.ic_model.item(row, 1)}
        self._matrix_busy = True
        left.setEnabled(False)

        def work(progress=None):
            required = [name for link in manual_links
                        for name in (link["parent"], link["daughter"]) if name]
            removed = {svc.normalize_label(name) for name in removals}
            if removed.intersection(required):
                raise ValueError("Remove the manual channel before removing its parent or daughter.")
            missing = [name for name in required if name not in svc.universe_labels()]
            if additions or removals or missing:
                svc.edit_workspace_isotopes(list(additions) + missing, removals)
            matrix = svc.build_matrix(
                flux=flux, energy_eV=reactor_eV, mt_filter=mt_filter,
                decay_mode_filter=decay_filter, open_boundary_mts=open_mts,
                hl_preset=hl_preset, manual_decay_links=manual_links, progress=progress)
            svc.set_manual_decay_links(manual_links)
            return matrix

        def finish():
            self._matrix_busy = False
            left.setEnabled(True)

        def done(matrix):
            finish()
            left.populate_unified_table(flux, open_mts, reactor_eV, astro_keV,
                                        decay_mode_filter=decay_filter)
            right._load_from_matrix()
            for row in range(right.ic_model.rowCount()):
                name = right.ic_model.item(row, 0).text()
                if name in initial:
                    right.ic_model.item(row, 1).setText(initial[name])
            left.iso_add_edit.clear()
            left.iso_search_edit.clear()

        def failed(message):
            finish()
            QMessageBox.warning(left, "Matrix rebuild failed", message)

        self._matrix_job = run_with_progress(
            left.window(), "Updating burnup matrix", work,
            on_done=done, on_error=failed, cancelable=False)

    # -- Evolution calculation -----------------------------------------------

    def on_calculate(self) -> None:
        svc = get_service()
        if not svc.initialized:
            QMessageBox.critical(None, "Burnup",
                                 "Nuclear-data environment is not initialised yet.")
            return
        right = self.__right

        try:
            initial = right.read_initial_conditions()
        except ValueError as e:
            QMessageBox.warning(None, "Initial conditions", str(e))
            return
        if not initial:
            QMessageBox.information(None, "Initial conditions",
                                    "Add at least one isotope with non-zero concentration.")
            return

        try:
            t_start_s = right.read_time_start_s()
            t_end_s   = right.read_time_end_s()
            flux      = right.read_flux()
            energy_eV = right.read_energy_eV() or 30_000.0
            method    = right.read_method()
            hl_preset         = right.hl_filter_box.currentText()
            mt_filter         = right.get_mt_filter()
            decay_mode_filter = right.get_decay_mode_filter()
            rho               = right.read_density()
        except ValueError as e:
            QMessageBox.warning(None, "Inputs", str(e))
            return

        left = self.__left
        parent = left
        requested_targets = right.read_targets()
        right.calculate_pb.setEnabled(False)

        def done(res: RunResult):
            right.calculate_pb.setEnabled(True)
            applied = set(res.initial.keys())
            ignored = [lab for lab in initial
                       if svc.normalize_label(lab) not in applied]
            if ignored:
                QMessageBox.warning(
                    parent, "Unknown isotopes",
                    "These rows are not in the current burnup matrix and were ignored:\n  - "
                    + "\n  - ".join(ignored))
            if requested_targets:
                resolved = set(res.targets)
                unknown = [t for t in requested_targets
                           if svc.normalize_label(t) not in resolved]
                if unknown:
                    QMessageBox.warning(
                        parent, "Unknown targets",
                        "These target isotopes are not in the current universe and "
                        "were ignored"
                        + (" (showing the full network):" if not resolved
                           else ":") + "\n  - " + "\n  - ".join(unknown))
            svc.last_rho = rho
            left.populate_evolution(res, rho)
            self.__left.calc_tabs.setCurrentIndex(1)

        def error(tb):
            right.calculate_pb.setEnabled(True)
            if "cancelled" in tb.lower():
                return
            QMessageBox.critical(parent, "Solver failed", tb)

        run_with_progress(
            parent, "Running burnup evolution…",
            svc.run_evolution, initial, t_end_s,
            flux=flux, energy_eV=energy_eV, method=method,
            mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
            hl_preset=hl_preset,
            t_start_s=max(t_start_s, 1e-9),
            open_boundary_mts=right.get_matrix_open_mts(),
            target_labels=right.read_targets(),
            on_done=done, on_error=error)

    # -- Equilibrium calculation ---------------------------------------------

    def on_calculate_equilibrium(self) -> None:
        svc = get_service()
        if not svc.initialized:
            QMessageBox.critical(None, "Burnup",
                                 "Nuclear-data environment is not initialised yet.")
            return
        right = self.__right
        try:
            initial   = right.read_initial_conditions()
            energy_eV = right.read_energy_eV() or 30_000.0
            method            = right.read_method()
            mt_filter         = right.get_mt_filter()
            decay_mode_filter = right.get_decay_mode_filter()
            hl_preset         = right.hl_filter_box.currentText()
            rho               = right.read_density()
        except ValueError as e:
            QMessageBox.warning(None, "Inputs", str(e))
            return
        if not initial:
            QMessageBox.information(None, "Initial conditions",
                                    "Add at least one isotope with non-zero concentration.")
            return

        auto_mode = right.eq_auto_rb.isChecked()

        try:
            flux_list = right.read_eq_flux_list()
        except ValueError as e:
            QMessageBox.warning(None, "Equilibrium", str(e))
            return
        if not flux_list:
            QMessageBox.warning(None, "Equilibrium", "Enter at least one flux value.")
            return

        if not auto_mode:
            try:
                time_list = right.read_eq_time_list()
            except ValueError as e:
                QMessageBox.warning(None, "Equilibrium", str(e))
                return
            if len(time_list) != len(flux_list):
                QMessageBox.warning(None, "Equilibrium",
                                    "Flux list and time list must have the same number of values.")
                return

        left = self.__left
        parent = left
        right.calc_eq_pb.setEnabled(False)

        def done(eq: EqResult):
            right.calc_eq_pb.setEnabled(True)
            svc.last_rho = rho
            left.populate_equilibrium(eq, rho)
            self.__left.calc_tabs.setCurrentIndex(2)

        def error(tb):
            right.calc_eq_pb.setEnabled(True)
            if "cancelled" in tb.lower():
                return
            QMessageBox.critical(parent, "Equilibrium scan failed", tb)

        open_mts = right.get_matrix_open_mts()
        if auto_mode:
            run_with_progress(
                parent, "Running auto-plateau equilibrium scan…",
                svc.run_equilibrium_auto, initial, flux_list,
                energy_eV=energy_eV, method=method,
                mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
                hl_preset=hl_preset,
                open_boundary_mts=open_mts,
                on_done=done, on_error=error)
        else:
            run_with_progress(
                parent, "Running equilibrium scan…",
                svc.run_equilibrium_scan, initial, flux_list, time_list,
                energy_eV=energy_eV, method=method,
                mt_filter=mt_filter, decay_mode_filter=decay_mode_filter,
                hl_preset=hl_preset,
                open_boundary_mts=open_mts,
                on_done=done, on_error=error)

    # -- Sensitivity analysis (OAT) -----------------------------------------

    def on_calculate_sensitivity(self) -> None:
        svc = get_service()
        if not svc.initialized:
            QMessageBox.critical(None, "Sensitivity",
                                 "Nuclear-data environment is not initialised yet.")
            return
        right = self.__right
        try:
            initial = right.read_initial_conditions()
        except ValueError as e:
            QMessageBox.warning(None, "Initial conditions", str(e))
            return
        if not initial:
            QMessageBox.information(None, "Initial conditions",
                                    "Add at least one isotope with non-zero concentration.")
            return

        cfg = right.read_sensitivity_config(initial)
        if not cfg.types:
            QMessageBox.information(None, "Sensitivity",
                                    "Select at least one sensitivity type.")
            return
        if cfg.target_kind == "concentration" and not cfg.target_label:
            QMessageBox.information(None, "Sensitivity",
                                    "Choose a target isotope, or switch the target to Heat Q.")
            return
        if "library" in cfg.types and not cfg.libraries:
            QMessageBox.information(None, "Sensitivity",
                                    "Select at least one library for the Library type.")
            return

        left = self.__left
        parent = left
        right.sens_calc_pb.setEnabled(False)

        def done(res):
            right.sens_calc_pb.setEnabled(True)
            left.populate_sensitivity(res)
            self.__left.calc_tabs.setCurrentIndex(3)

        def error(tb):
            right.sens_calc_pb.setEnabled(True)
            if "cancelled" in tb.lower():
                return
            QMessageBox.critical(parent, "Sensitivity failed", tb)

        run_with_progress(
            parent, "Running sensitivity analysis…",
            svc.run_sensitivity, cfg,
            on_done=done, on_error=error)

    # -- Export Evolution ---------------------------------------------------

    def on_export_csv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            None, "Export evolution (CSV)", "results.csv", "CSV (*.csv)")
        if not path:
            return
        run = get_service().last_run
        if not run:
            QMessageBox.warning(None, "Export", "No evolution results to export.")
            return
        try:
            target_set = set(getattr(run, "targets", []) or [])
            total_per_step = [float(np.sum(N)) or 1.0 for N in run.trajectories]
            with open(path, "w", newline="") as f:
                w = csv.writer(f)
                # --- per-isotope section ---
                w.writerow(["Time_s", "Isotope", "A", "N", "Fraction"])
                for step_i, (t, N_vec) in enumerate(zip(run.times, run.trajectories)):
                    total = total_per_step[step_i]
                    for i, name in enumerate(run.names):
                        if target_set and name not in target_set:
                            continue
                        n_val = float(N_vec[i])
                        if n_val > 0:
                            A = _parse_A(name)
                            w.writerow([f"{t:.6e}", name,
                                        str(A) if A is not None else "",
                                        f"{n_val:.8e}", f"{n_val/total:.6e}"])
                # --- N(A) section ---
                w.writerow([])
                w.writerow(["# N(A) by mass number"])
                w.writerow(["Time_s", "A", "N_A", "Fraction_A"])
                for step_i, (t, N_vec) in enumerate(zip(run.times, run.trajectories)):
                    total = total_per_step[step_i]
                    by_A: dict = {}
                    for i, name in enumerate(run.names):
                        if target_set and name not in target_set:
                            continue
                        n_val = float(N_vec[i])
                        if n_val <= 0:
                            continue
                        A = _parse_A(name)
                        if A is not None:
                            by_A[A] = by_A.get(A, 0.0) + n_val
                    for A in sorted(by_A):
                        na = by_A[A]
                        w.writerow([f"{t:.6e}", A, f"{na:.8e}", f"{na/total:.6e}"])
        except OSError as e:
            QMessageBox.critical(None, "CSV export failed", str(e))

    def on_export_txt(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            None, "Export evolution (TXT)", "results.txt", "Text (*.txt)")
        if not path:
            return
        run = get_service().last_run
        if not run:
            QMessageBox.warning(None, "Export", "No evolution results to export.")
            return
        try:
            target_set = set(getattr(run, "targets", []) or [])
            total_per_step = [float(np.sum(N)) or 1.0 for N in run.trajectories]
            with open(path, "w", encoding="utf-8") as f:
                f.write(f"# Burnup evolution — t_end = {run.t_end_s} s, "
                        f"flux = {run.flux} n/cm2/s, method = {run.method}\n")
                # --- per-isotope section ---
                f.write(f"{'Time_s':>14}  {'Isotope':<10}  {'A':>4}  "
                        f"{'N':>16}  {'Fraction':>12}\n")
                for step_i, (t, N_vec) in enumerate(zip(run.times, run.trajectories)):
                    total = total_per_step[step_i]
                    for i, name in enumerate(run.names):
                        if target_set and name not in target_set:
                            continue
                        n_val = float(N_vec[i])
                        if n_val > 0:
                            A = _parse_A(name)
                            a_str = f"{A:>4}" if A is not None else "   —"
                            f.write(f"{t:>14.6e}  {name:<10}  {a_str}  "
                                    f"{n_val:>16.8e}  {n_val/total:>12.6e}\n")
                # --- N(A) section ---
                f.write("\n# N(A) by mass number\n")
                f.write(f"{'Time_s':>14}  {'A':>4}  {'N_A':>16}  {'Fraction_A':>12}\n")
                for step_i, (t, N_vec) in enumerate(zip(run.times, run.trajectories)):
                    total = total_per_step[step_i]
                    by_A: dict = {}
                    for i, name in enumerate(run.names):
                        if target_set and name not in target_set:
                            continue
                        n_val = float(N_vec[i])
                        if n_val <= 0:
                            continue
                        A = _parse_A(name)
                        if A is not None:
                            by_A[A] = by_A.get(A, 0.0) + n_val
                    for A in sorted(by_A):
                        na = by_A[A]
                        f.write(f"{t:>14.6e}  {A:>4}  {na:>16.8e}  "
                                f"{na/total:>12.6e}\n")
        except OSError as e:
            QMessageBox.critical(None, "TXT export failed", str(e))

    def on_export_eq_csv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            None, "Export equilibrium (CSV)", "equilibrium.csv", "CSV (*.csv)")
        if not path:
            return
        try:
            with open(path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["Flux_n_cm2_s", "Time_s", "Q_MeV_s", "Q_Wcm3"])
                for row in self.__left.iter_eq_rows():
                    w.writerow(row)
        except OSError as e:
            QMessageBox.critical(None, "CSV export failed", str(e))

    # -- Session save / open -----------------------------------------------

    def on_save_session(self) -> None:
        svc = get_service()
        run = svc.last_run
        if run is None:
            QMessageBox.warning(None, "Save session",
                                "No evolution results to save.\nRun a calculation first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Save session", "session.kaz", "NuMatRx session (*.kaz);;All files (*)")
        if not path:
            return
        try:
            from core.session import save_session
            rho = getattr(svc, "last_rho", 0.0) or 0.0
            eq  = getattr(svc, "last_eq",  None)
            save_session(run, eq, rho, path)
        except Exception as e:
            QMessageBox.critical(None, "Save session failed", str(e))

    def on_open_session(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            None, "Open saved session", "", "NuMatRx session (*.kaz);;All files (*)")
        if not path:
            return
        try:
            from core.session import load_session
            run, eq, rho = load_session(path)
            dlg = ResultsViewerDialog(run, eq, rho, path)
            dlg.show()
            ResultsViewerDialog._open_viewers.append(dlg)
            dlg.destroyed.connect(
                lambda: ResultsViewerDialog._open_viewers.remove(dlg)
                if dlg in ResultsViewerDialog._open_viewers else None)
        except Exception as e:
            QMessageBox.critical(None, "Open session failed", str(e))


# ---------------------------------------------------------------------------
# ResultsViewerDialog — standalone non-modal window for saved sessions
# ---------------------------------------------------------------------------

class ResultsViewerDialog(QDialog):
    """Non-modal results viewer for a previously saved .kaz session file."""

    _open_viewers: list = []   # module-level refs so GC won't close the windows

    _THEMES = {
        "Dark":  {"bg": "#1a1714", "fg": "#f5ede4", "accent": "#e8943a", "grid": "#3d352e"},
        "Light": {"bg": "#f4f7fa", "fg": "#16222e", "accent": "#2b7fc4", "grid": "#d2dde6"},
    }

    def __init__(self, run: RunResult, eq, rho: float, source_path: str, parent=None):
        super().__init__(parent, Qt.Window)
        self._run  = run
        self._eq   = eq      # EqResult or None
        self._rho  = rho
        self._step = max(0, len(run.times) - 1)

        self.setWindowTitle(f"Saved session — {Path(source_path).name}")
        self.resize(920, 680)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self._build_ui()
        self._on_evo_slider(self._step)

    # ------------------------------------------------------------------
    def _build_ui(self):
        run  = self._run
        root = QVBoxLayout(self)

        info = (f"t_end = {run.t_end_s:.4g} s   |   "
                f"Φ = {run.flux:.4g} n/cm²/s   |   "
                f"E = {run.energy_eV:.4g} eV   |   "
                f"method = {run.method}   |   "
                f"isotopes = {len(run.names)}")
        hdr = QLabel(info)
        hdr.setWordWrap(True)
        root.addWidget(hdr)

        self._tabs = QTabWidget()
        root.addWidget(self._tabs)

        self._tabs.addTab(self._build_evo_tab(),   "Evolution")
        self._tabs.addTab(self._build_graph_tab(), "Graph")
        self.exposure_panel = ExposurePanel()
        self.exposure_panel.session_eq = self._eq
        self.exposure_panel.session_rho = self._rho
        self.exposure_panel.set_run(run)
        self._tabs.addTab(self.exposure_panel, "Trajectory averaging")

        if self._eq is not None and getattr(self._eq, "flux_list", None):
            self._tabs.addTab(self._build_eq_tab(), "Equilibrium")

    # ---- Evolution tab ---------------------------------------------------
    def _build_evo_tab(self):
        run = self._run
        w   = QWidget()
        v   = QVBoxLayout(w)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Display:"))
        self._abs_rb = QRadioButton("Absolute")
        self._rel_rb = QRadioButton("Relative")
        self._sn_rb  = QRadioButton("σN")
        self._abs_rb.setChecked(True)
        for rb in (self._abs_rb, self._rel_rb, self._sn_rb):
            rb.toggled.connect(self._refresh_evo)
            mode_row.addWidget(rb)
        mode_row.addStretch(1)
        v.addLayout(mode_row)

        slider_row = QHBoxLayout()
        self._evo_slider = QSlider(Qt.Horizontal)
        self._evo_slider.setMinimum(0)
        self._evo_slider.setMaximum(max(0, len(run.times) - 1))
        self._evo_slider.setValue(self._step)
        self._evo_slider.valueChanged.connect(self._on_evo_slider)
        self._evo_time_lbl = QLabel("t = —")
        slider_row.addWidget(QLabel("Time step:"))
        slider_row.addWidget(self._evo_slider, 1)
        slider_row.addWidget(self._evo_time_lbl)
        v.addLayout(slider_row)

        self._evo_mdl = QStandardItemModel()
        self._evo_mdl.setColumnCount(4)
        self._evo_mdl.setHorizontalHeaderLabels(["Isotope", "N", "Fraction", "N(A)"])
        self._evo_tbl = QTableView()
        self._evo_tbl.setModel(self._evo_mdl)
        self._evo_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._evo_tbl.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        v.addWidget(self._evo_tbl)

        self._heat_lbl = QLabel("Q = —")
        v.addWidget(self._heat_lbl)
        return w

    # ---- Graph tab -------------------------------------------------------
    def _build_graph_tab(self):
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        w = QWidget()
        v = QVBoxLayout(w)

        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel("Plot:"))
        self._graph_combo = QComboBox()
        self._graph_combo.addItems(["Concentrations", "Heat", "Abundance"])
        self._graph_combo.currentIndexChanged.connect(self._on_graph_kind_changed)
        ctrl_row.addWidget(self._graph_combo)
        ctrl_row.addWidget(QLabel("  Theme:"))
        self._graph_theme_combo = QComboBox()
        self._graph_theme_combo.addItems(["Dark", "Light"])
        self._graph_theme_combo.currentIndexChanged.connect(self._replot_graph)
        ctrl_row.addWidget(self._graph_theme_combo)
        ctrl_row.addStretch(1)
        v.addLayout(ctrl_row)

        # Abundance-specific controls
        self._abund_ctrl = QWidget()
        abh = QHBoxLayout(self._abund_ctrl)
        abh.setContentsMargins(0, 0, 0, 0)
        abh.addWidget(QLabel("Step:"))
        self._abund_slider = QSlider(Qt.Horizontal)
        self._abund_slider.setMinimum(0)
        self._abund_slider.setMaximum(max(0, len(self._run.times) - 1))
        self._abund_slider.setValue(max(0, len(self._run.times) - 1))
        self._abund_slider.valueChanged.connect(self._replot_graph)
        self._abund_time_lbl = QLabel("t = —")
        abh.addWidget(self._abund_slider, 1)
        abh.addWidget(self._abund_time_lbl)
        abh.addWidget(QLabel("  Mode:"))
        self._abund_source = QComboBox()
        self._abund_source.addItems(["Time slice", "Averaged results"])
        self._abund_source.currentIndexChanged.connect(self._replot_graph)
        abh.addWidget(self._abund_source)
        self._abund_mode_combo = QComboBox()
        self._abund_mode_combo.addItems(["Absolute", "Relative", "σN"])
        self._abund_mode_combo.currentIndexChanged.connect(self._replot_graph)
        abh.addWidget(self._abund_mode_combo)
        self._abund_ctrl.setVisible(False)
        v.addWidget(self._abund_ctrl)

        self._graph_fig    = Figure(figsize=(7, 5), dpi=100)
        self._graph_canvas = FigureCanvasQTAgg(self._graph_fig)
        v.addWidget(self._graph_canvas)

        self._replot_graph()
        return w

    # ---- Equilibrium tab -------------------------------------------------
    def _build_eq_tab(self):
        eq  = self._eq
        w   = QWidget()
        v   = QVBoxLayout(w)
        v.addWidget(QLabel(f"Equilibrium scan — {len(eq.flux_list)} flux points   |   "
                           f"method = {eq.method}"))
        mdl = QStandardItemModel()
        mdl.setColumnCount(4)
        mdl.setHorizontalHeaderLabels(
            ["Flux [n/cm²/s]", "Time [s]", "Q [MeV/s]", "Q [W/cm³]"])
        conv = _heat_conv_factor(self._rho, eq.names, eq.isotopes, eq.initial)
        for flux, t_irr, q_mev in zip(eq.flux_list, eq.time_list, eq.Q_eq):
            q_wcm3 = (q_mev * conv) if conv is not None else None
            mdl.appendRow([
                QStandardItem(f"{flux:.4e}"),
                QStandardItem(f"{t_irr:.4e}"),
                QStandardItem(f"{q_mev:.5e}"),
                QStandardItem(_fmt_e(q_wcm3)),
            ])
        tbl = QTableView()
        tbl.setModel(mdl)
        tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        tbl.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        v.addWidget(tbl)
        return w

    # ------------------------------------------------------------------
    def _current_evo_mode(self) -> str:
        if self._rel_rb.isChecked():
            return "relative"
        if self._sn_rb.isChecked():
            return "sigma_n"
        return "absolute"

    def _on_evo_slider(self, idx: int):
        self._step = idx
        run = self._run
        if run.times and idx < len(run.times):
            self._evo_time_lbl.setText(f"t = {_human_time(run.times[idx])}")
        self._refresh_evo()

    def _refresh_evo(self, *_):
        run  = self._run
        step = self._step
        if not run.trajectories or step >= len(run.trajectories):
            return

        N_vec = run.trajectories[step]
        total = float(np.sum(N_vec)) or 1.0
        mode  = self._current_evo_mode()

        sigma_by_name: Dict[str, float] = {}
        if mode == "sigma_n":
            for iso in run.isotopes:
                for rx in iso.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if mt == 102:
                        sigma_by_name[iso.name] = float(
                            getattr(rx, "sigma_barn", 0.0) or 0.0)
                        break

        target_set = set(getattr(run, "targets", []) or [])
        rows_data = []
        for i, name in enumerate(run.names):
            if target_set and name not in target_set:
                continue
            N_i = float(N_vec[i])
            if N_i <= 0 and name not in target_set:
                continue
            frac = N_i / total
            if mode == "relative":
                display_n = frac
            elif mode == "sigma_n":
                display_n = sigma_by_name.get(name, 0.0) * N_i
            else:
                display_n = N_i
            rows_data.append((name, N_i, frac, display_n))

        by_A: Dict[int, float] = {}
        for name, N_i, frac, display_n in rows_data:
            A = _parse_A(name)
            if A is not None:
                by_A[A] = by_A.get(A, 0.0) + display_n

        rows_data.sort(key=lambda r: r[1], reverse=True)
        self._evo_mdl.removeRows(0, self._evo_mdl.rowCount())
        for name, N_i, frac, display_n in rows_data:
            A  = _parse_A(name)
            na = by_A.get(A) if A is not None else None
            self._evo_mdl.appendRow([
                QStandardItem(name),
                QStandardItem(_fmt_e(display_n)),
                QStandardItem(f"{frac:.5e}"),
                QStandardItem(_fmt_e(na)),
            ])

        if by_A:
            sep = QStandardItem("── N(A) by mass number ──")
            sep.setEnabled(False)
            self._evo_mdl.appendRow([sep, QStandardItem(""),
                                     QStandardItem(""), QStandardItem("")])
            for A in sorted(by_A):
                self._evo_mdl.appendRow([
                    QStandardItem(f"A = {A}"),
                    QStandardItem(""),
                    QStandardItem(""),
                    QStandardItem(_fmt_e(by_A[A])),
                ])

        q_str = "Q = —"
        if run.heat and step < len(run.heat):
            q_mev = run.heat[step]
            q_str = f"Q = {_fmt_e(q_mev)} MeV/s"
            if self._rho > 0:
                conv = _heat_conv_factor(self._rho, run.names, run.isotopes, run.initial)
                if conv is not None:
                    q_str += f"  =  {_fmt_e(q_mev * conv)} W/cm³"
        self._heat_lbl.setText(q_str)

    # ---- Graph helpers ---------------------------------------------------
    def _on_graph_kind_changed(self, *_):
        is_abund = (self._graph_combo.currentText() == "Abundance")
        self._abund_ctrl.setVisible(is_abund)
        self._replot_graph()

    def _replot_graph(self, *_):
        run    = self._run
        theme  = self._graph_theme_combo.currentText()
        colors = self._THEMES.get(theme, self._THEMES["Dark"])
        kind   = self._graph_combo.currentText()

        fig = self._graph_fig
        fig.clear()
        fig.patch.set_facecolor(colors["bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(colors["bg"])
        for spine in ax.spines.values():
            spine.set_color(colors["fg"])
        ax.tick_params(colors=colors["fg"], which="both", labelsize=11)
        ax.xaxis.label.set_color(colors["fg"])
        ax.yaxis.label.set_color(colors["fg"])
        ax.title.set_color(colors["fg"])
        ax.grid(True, alpha=0.25, color=colors["grid"], which="both")

        try:
            if kind == "Heat":
                self._gplot_heat(ax, run, colors)
            elif kind == "Abundance":
                idx      = self._abund_slider.value()
                mode_txt = self._abund_mode_combo.currentText()
                mode     = {"Relative": "relative", "σN": "sigma_n"}.get(mode_txt, "absolute")
                if run.times and idx < len(run.times):
                    self._abund_time_lbl.setText(_human_time(run.times[idx]))
                averaged = self._abund_source.currentIndex() == 1
                self._abund_slider.setEnabled(not averaged)
                self._abund_time_lbl.setVisible(not averaged)
                if averaged:
                    from model.modelgraph.exposure_plot import plot_exposure
                    result = getattr(run, "exposure_result", None)
                    if result is None:
                        ax.text(0.5, 0.5, "Run Trajectory averaging first.",
                                transform=ax.transAxes, ha="center", color=colors["fg"], usetex=False)
                    else:
                        plot_exposure(ax, result, mode, colors)
                else:
                    self._gplot_abund(ax, run, colors, idx, mode)
            else:
                self._gplot_conc(ax, run, colors)
        except Exception as exc:
            ax.text(0.5, 0.5, f"Plot error:\n{exc}", ha="center", va="center",
                    transform=ax.transAxes, color=colors["fg"], fontsize=12)

        try:
            fig.tight_layout()
        except Exception:
            pass
        self._graph_canvas.draw()

    def _gplot_heat(self, ax, run, colors):
        xs = list(run.times)
        ys = list(run.heat)
        ax.plot(xs, ys, color=colors["accent"], linewidth=1.8,
                marker="s", markersize=4, markevery=max(1, len(xs) // 20))
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Decay heat [MeV/s per nucleus]")
        ax.set_title("Decay-heat profile")
        if xs and max(xs) / max(min(xs), 1e-300) > 100:
            ax.set_xscale("log")
        pos = [v for v in ys if v > 0]
        if pos and max(pos) / min(pos) > 100:
            ax.set_yscale("log")

    def _gplot_conc(self, ax, run, colors):
        M     = np.array(run.trajectories, dtype=float)
        peak  = np.max(np.abs(M), axis=0)
        thr   = 1e-6 * (peak.max() or 1.0)
        shown = 0
        for j, name in enumerate(run.names):
            if peak[j] < thr:
                continue
            ax.plot(run.times, M[:, j], linewidth=1.0, label=name)
            shown += 1
            if shown >= 40:
                break
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Atom count N")
        ax.set_title(f"Isotope concentrations (top {shown} by peak, absolute)")
        ax.set_xscale("log")
        ax.set_yscale("log")

    def _gplot_abund(self, ax, run, colors, step, mode):
        if not run.trajectories:
            return
        idx   = step if 0 <= step < len(run.trajectories) else len(run.trajectories) - 1
        N_vec = run.trajectories[idx]
        t_val = run.times[idx] if idx < len(run.times) else 0.0
        total = float(np.sum(N_vec)) or 1.0

        sigma_by_name: dict = {}
        if mode == "sigma_n":
            for iso in run.isotopes:
                for rx in iso.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if mt == 102:
                        sigma_by_name[iso.name] = float(
                            getattr(rx, "sigma_barn", 0.0) or 0.0)
                        break

        by_A: dict = {}
        for name, N_i_raw in zip(run.names, N_vec):
            N_i = float(N_i_raw)
            if N_i <= 0:
                continue
            val = (N_i / total         if mode == "relative"
                   else sigma_by_name.get(name, 0.0) * N_i if mode == "sigma_n"
                   else N_i)
            A = _parse_A(name)
            if A is not None:
                by_A[A] = by_A.get(A, 0.0) + val

        if not by_A:
            return
        As   = sorted(by_A)
        vals = [by_A[a] for a in As]
        ax.plot(As, vals, "-", color=colors["accent"], linewidth=1.2)
        ax.plot(As, vals, "o", color=colors["accent"], markersize=3)
        ax.set_xlabel("Mass number A")
        ylabels = {"relative": "N(A) / N_tot", "sigma_n": "σ·N(A) [barn·atoms]"}
        ax.set_ylabel(ylabels.get(mode, "N(A) [atoms]"))
        mode_label = {"absolute": "Absolute", "relative": "Relative", "sigma_n": "Capture-weighted abundance"}.get(mode, mode)
        ax.set_title(f"Abundance vs A   (t = {_human_time(t_val)},  {mode_label})")
        pos = [v for v in vals if v > 0]
        if pos and max(pos) / max(min(pos), 1e-300) > 100:
            ax.set_yscale("log")


# ---------------------------------------------------------------------------
# Left — results display (three tabs)
# ---------------------------------------------------------------------------

class Left(QWidget):

    def __init__(self):
        super().__init__()
        self.setMinimumSize(100, 100)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._evo_rho: float = 0.0
        self._evo_res: Optional[RunResult] = None
        self._eq_rho: float = 0.0
        self._eq_res: Optional[EqResult] = None

        root = QVBoxLayout(self)

        self.calc_tabs = QTabWidget()
        self.calc_tabs.setObjectName("result_calc_tabs")
        self.calc_tabs.setElideMode(Qt.ElideNone)
        self.calc_tabs.tabBar().setExpanding(False)
        self.calc_tabs.tabBar().setUsesScrollButtons(True)

        # ---- Tab 0: Isotopes & Matrix ----
        iso_w = QWidget()
        iso_v = QVBoxLayout(iso_w)

        self._user_overrides: Dict[str, dict] = {}
        self._unified_updating = False

        _unified_cols = [
            "Isotope", "T½", "λ [1/s]", "Decay mode",
            "σ(n,γ) [b]", "rate(n,γ) [1/s]",
            "σ(n,2n) [b]", "rate(n,2n) [1/s]",
            "σ(n,p) [b]", "rate(n,p) [1/s]",
            "σ(n,α) [b]", "rate(n,α) [1/s]",
            "σ(n,f) [b]", "rate(n,f) [1/s]",
            "Total rate [1/s]", "Source",
        ]
        self.unified_model = QStandardItemModel()
        self.unified_model.setColumnCount(len(_unified_cols))
        self.unified_model.setHorizontalHeaderLabels(_unified_cols)
        self.unified_model.itemChanged.connect(self._on_unified_cell_edited)

        self.unified_proxy = QSortFilterProxyModel()
        self.unified_proxy.setSourceModel(self.unified_model)
        self.unified_proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.unified_proxy.setFilterKeyColumn(0)

        # Search + action row
        iso_search_row = QHBoxLayout()
        self.iso_search_edit = QLineEdit()
        self.iso_search_edit.setObjectName("iso_search_edit")
        self.iso_search_edit.setPlaceholderText("Search isotope…  (Ctrl+F)")
        self.iso_search_edit.setClearButtonEnabled(True)
        self.iso_search_edit.textChanged.connect(self.unified_proxy.setFilterFixedString)
        iso_search_row.addWidget(self.iso_search_edit, 1)
        self.iso_copy_pb = QPushButton("Copy")
        self.iso_copy_pb.setObjectName("iso_copy_pb")
        self.iso_copy_pb.setToolTip("Copy selected rows to clipboard (Tab-separated)")
        self.iso_copy_pb.clicked.connect(self._on_copy_matrix)
        self.iso_export_csv_pb = QPushButton("Export CSV")
        self.iso_export_csv_pb.setObjectName("iso_export_csv_pb")
        self.iso_export_csv_pb.setToolTip("Save full matrix table to CSV file")
        self.iso_export_csv_pb.clicked.connect(self._on_export_matrix_csv)
        iso_search_row.addWidget(self.iso_copy_pb)
        iso_search_row.addWidget(self.iso_export_csv_pb)
        iso_v.addLayout(iso_search_row)

        membership_row = QHBoxLayout()
        self.iso_add_edit = QLineEdit()
        self.iso_add_edit.setPlaceholderText("Add isotopes: Fe-56 Ni-60")
        self.iso_add_pb = QPushButton("Add isotopes")
        self.iso_add_pb.setObjectName("iso_add_pb")
        self.iso_remove_pb = QPushButton("Remove selected")
        self.iso_remove_pb.setObjectName("iso_remove_pb")
        self.iso_remove_pb.setToolTip("Remove selected isotopes from the calculation network")
        membership_row.addWidget(self.iso_add_edit, 1)
        membership_row.addWidget(self.iso_add_pb)
        membership_row.addWidget(self.iso_remove_pb)
        iso_v.addLayout(membership_row)

        self.unified_table = QTableView()
        self.unified_table.setObjectName("unified_isotope_table")
        self.unified_table.setModel(self.unified_proxy)
        self.unified_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.unified_table.horizontalHeader().setStretchLastSection(True)
        self.unified_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        iso_v.addWidget(self.unified_table)

        # Cycle results sub-table
        iso_v.addWidget(QLabel("Cycle reactions:"))
        self.cycles_model = QStandardItemModel()
        self.cycles_model.setColumnCount(4)
        self.cycles_model.setHorizontalHeaderLabels(
            ["Cycle #", "Length", "Isotopes", "Type"])
        self.cycles_table = QTableView()
        self.cycles_table.setObjectName("cycles_table")
        self.cycles_table.setModel(self.cycles_model)
        self.cycles_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.cycles_table.setMaximumHeight(120)
        iso_v.addWidget(self.cycles_table)

        # Manual decay channels editor
        mdl_hdr = QHBoxLayout()
        mdl_hdr.addWidget(QLabel("Manual decay channels:"))
        self.mdl_add_pb = QPushButton("+ Add")
        self.mdl_add_pb.setObjectName("mdl_add_pb")
        self.mdl_del_pb = QPushButton("− Delete")
        self.mdl_del_pb.setObjectName("mdl_del_pb")
        mdl_hdr.addStretch(1)
        mdl_hdr.addWidget(self.mdl_add_pb)
        mdl_hdr.addWidget(self.mdl_del_pb)
        iso_v.addLayout(mdl_hdr)

        self.mdl_model = QStandardItemModel()
        self.mdl_model.setColumnCount(5)
        self.mdl_model.setHorizontalHeaderLabels(
            ["Parent", "Daughter (blank=sink)", "Mode", "T½ [s]", "Branch"])
        self.mdl_table = QTableView()
        self.mdl_table.setObjectName("manual_decay_table")
        self.mdl_table.setModel(self.mdl_model)
        self.mdl_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.mdl_table.horizontalHeader().setStretchLastSection(True)
        self.mdl_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.mdl_table.setMaximumHeight(130)
        iso_v.addWidget(self.mdl_table)

        self.mdl_add_pb.clicked.connect(self._mdl_add_row)
        self.mdl_del_pb.clicked.connect(self._mdl_del_row)

        # ---- Tab 1: Evolution ----
        evo_w = QWidget()
        evo_v = QVBoxLayout(evo_w)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Display:"))
        self.mode_abs_rb  = QRadioButton("Absolute")
        self.mode_rel_rb  = QRadioButton("Relative")
        self.mode_sn_rb   = QRadioButton("σN")
        self.mode_abs_rb.setObjectName("mode_abs_rb")
        self.mode_rel_rb.setObjectName("mode_rel_rb")
        self.mode_sn_rb.setObjectName("mode_sn_rb")
        self.mode_abs_rb.setChecked(True)
        for rb in (self.mode_abs_rb, self.mode_rel_rb, self.mode_sn_rb):
            rb.toggled.connect(self._on_mode_changed)
            mode_row.addWidget(rb)
        mode_row.addStretch(1)   # keep the radios packed left, not spread out
        evo_v.addLayout(mode_row)

        slider_row = QHBoxLayout()
        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setObjectName("time_slider")
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(0)
        self.time_slider.valueChanged.connect(self._on_slider_changed)
        self.time_label = QLabel("t = —")
        self.time_label.setObjectName("time_label_evo")
        slider_row.addWidget(QLabel("Time step:"))
        slider_row.addWidget(self.time_slider, 1)
        slider_row.addWidget(self.time_label)
        evo_v.addLayout(slider_row)

        self.evo_model = QStandardItemModel()
        self.evo_model.setColumnCount(6)
        self.evo_model.setHorizontalHeaderLabels(
            ["Isotope", "N", "Fraction", "N(A)", "Q [MeV/s]", "Q [W/cm³]"])
        self.evo_table = QTableView()
        self.evo_table.setObjectName("evo_result_table")
        self.evo_table.setModel(self.evo_model)
        self.evo_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        evo_v.addWidget(self.evo_table)

        self.total_q_label = QLabel("Total Q: —")
        self.total_q_label.setObjectName("total_q_label")
        evo_v.addWidget(self.total_q_label)

        export_row = QHBoxLayout()
        export_row.addWidget(QLabel("Export:"))
        self.csv_pb = QPushButton("csv")
        self.csv_pb.setObjectName("result_csv_pb")
        self.txt_pb = QPushButton("txt")
        self.txt_pb.setObjectName("result_txt_pb")
        export_row.addWidget(self.csv_pb)
        export_row.addWidget(self.txt_pb)
        evo_v.addLayout(export_row)

        session_row = QHBoxLayout()
        session_row.addWidget(QLabel("Session:"))
        self.save_session_pb = QPushButton("Save…")
        self.save_session_pb.setObjectName("result_save_session_pb")
        self.save_session_pb.setToolTip("Save current results to a .kaz file")
        self.open_session_pb = QPushButton("Open…")
        self.open_session_pb.setObjectName("result_open_session_pb")
        self.open_session_pb.setToolTip("Open a saved .kaz session in a separate window")
        session_row.addWidget(self.save_session_pb)
        session_row.addWidget(self.open_session_pb)
        evo_v.addLayout(session_row)

        # ---- Tab 2: Equilibrium ----
        eq_w = QWidget()
        eq_v = QVBoxLayout(eq_w)
        eq_v.addWidget(QLabel("Equilibrium results (flux scan):"))

        self.eq_model = QStandardItemModel()
        self.eq_model.setColumnCount(4)
        self.eq_model.setHorizontalHeaderLabels(
            ["Flux [n/cm²/s]", "Time [s]", "Q [MeV/s]", "Q [W/cm³]"])
        self.eq_table = QTableView()
        self.eq_table.setObjectName("eq_result_table")
        self.eq_table.setModel(self.eq_model)
        self.eq_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        eq_v.addWidget(self.eq_table)

        eq_export_row = QHBoxLayout()
        self.eq_csv_pb = QPushButton("Export CSV")
        self.eq_csv_pb.setObjectName("eq_csv_pb")
        eq_export_row.addWidget(self.eq_csv_pb)
        eq_v.addLayout(eq_export_row)

        # ---- Tab 3: Sensitivity ----
        sens_w = QWidget()
        sens_v = QVBoxLayout(sens_w)
        self.sens_target_label = QLabel("Run a sensitivity analysis from the controls panel.")
        self.sens_target_label.setObjectName("sens_target_label")
        self.sens_target_label.setWordWrap(True)
        sens_v.addWidget(self.sens_target_label)

        self.sens_model = QStandardItemModel()
        self.sens_model.setColumnCount(6)
        self.sens_model.setHorizontalHeaderLabels(
            ["Rank", "Parameter", "Type", "Base", "Perturbed", "S"])
        self.sens_table = QTableView()
        self.sens_table.setObjectName("sens_result_table")
        self.sens_table.setModel(self.sens_model)
        self.sens_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        sens_v.addWidget(self.sens_table)

        self.sens_lib_caption = QLabel("Library spread (Nᵢ / Q per library):")
        self.sens_lib_caption.setObjectName("sens_lib_caption")
        self.sens_lib_caption.setVisible(False)
        sens_v.addWidget(self.sens_lib_caption)
        self.sens_lib_model = QStandardItemModel()
        self.sens_lib_model.setColumnCount(4)
        self.sens_lib_model.setHorizontalHeaderLabels(
            ["Library", "Nᵢ", "Q [MeV/s]", "Deviation %"])
        self.sens_lib_table = QTableView()
        self.sens_lib_table.setObjectName("sens_lib_result_table")
        self.sens_lib_table.setModel(self.sens_lib_model)
        self.sens_lib_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.sens_lib_table.setMaximumHeight(170)
        self.sens_lib_table.setVisible(False)
        sens_v.addWidget(self.sens_lib_table)

        self.calc_tabs.addTab(iso_w, "Isotopes & Matrix")
        self.calc_tabs.addTab(evo_w, "Evolution")
        self.calc_tabs.addTab(eq_w, "Equilibrium")
        self.calc_tabs.addTab(sens_w, "Sensitivity")
        self.exposure_panel = ExposurePanel()
        self.calc_tabs.addTab(self.exposure_panel, "Trajectory averaging")

        root.addWidget(self.calc_tabs)

        # Ctrl+F → focus search bar
        sc = QShortcut(QKeySequence("Ctrl+F"), self)
        sc.activated.connect(self._focus_search)

    def _focus_search(self):
        self.calc_tabs.setCurrentIndex(0)
        self.iso_search_edit.setFocus()
        self.iso_search_edit.selectAll()

    # -- display mode toggle -----------------------------------------------

    def _on_mode_changed(self, *_):
        if self._evo_res is not None:
            self.populate_evolution(self._evo_res, self._evo_rho)

    def _current_mode(self) -> str:
        if self.mode_rel_rb.isChecked():
            return "relative"
        if self.mode_sn_rb.isChecked():
            return "sigma_n"
        return "absolute"

    # -- Isotopes & Matrix population — unified table ----------------------

    _EDITABLE_COLS = {2, 4, 6, 8, 10, 12}   # λ, σ columns (0-indexed)

    def populate_isotope_matrix(self, flux: float = 0.0,
                                open_mts: Optional[List[int]] = None,
                                reactor_energy_eV: float = 30000.0,
                                astro_energy_keV: float = 30.0) -> None:
        """Populate the unified table (delegates to populate_unified_table)."""
        self.populate_unified_table(flux, open_mts, reactor_energy_eV, astro_energy_keV)

    def _on_copy_matrix(self):
        """Copy selected (or all) rows to clipboard as tab-separated text."""
        sel = self.unified_table.selectionModel().selectedRows()
        if sel:
            rows = sorted(idx.row() for idx in sel)
        else:
            rows = range(self.unified_proxy.rowCount())
        cols = self.unified_proxy.columnCount()
        headers = [self.unified_model.horizontalHeaderItem(c).text() for c in range(cols)]
        lines = ["\t".join(headers)]
        for r in rows:
            line = []
            for c in range(cols):
                src_idx = self.unified_proxy.mapToSource(self.unified_proxy.index(r, c))
                item = self.unified_model.itemFromIndex(src_idx)
                line.append(item.text() if item else "")
            lines.append("\t".join(line))
        QApplication.clipboard().setText("\n".join(lines))

    def _on_export_matrix_csv(self):
        """Save full matrix table to a CSV file."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Isotope Matrix", "", "CSV files (*.csv)")
        if not path:
            return
        cols = self.unified_model.columnCount()
        rows = self.unified_model.rowCount()
        headers = [self.unified_model.horizontalHeaderItem(c).text() for c in range(cols)]
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(headers)
                for r in range(rows):
                    row_data = []
                    for c in range(cols):
                        item = self.unified_model.item(r, c)
                        row_data.append(item.text() if item else "")
                    w.writerow(row_data)
            QMessageBox.information(self, "Export", f"Saved {rows} rows to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def populate_unified_table(self, flux: float = 0.0,
                               open_mts: Optional[List[int]] = None,
                               reactor_energy_eV: float = 30000.0,
                               astro_energy_keV: float = 30.0,
                               decay_mode_filter: Optional[List[str]] = None) -> None:
        """Populate the unified 16-column isotope/rate table."""
        import math as _math
        svc = get_service()
        if not svc.isotopes:
            return

        self._unified_updating = True
        self.unified_model.removeRows(0, self.unified_model.rowCount())

        def _sort_hl(iso):
            hl = getattr(iso, "half_life_s", None)
            dc = getattr(iso, "decay_constant", 0.0)
            if hl is None and dc and dc > 0:
                hl = _math.log(2) / dc
            if hl is None or hl >= 1e30:
                return float("inf")
            return float(hl)

        self._matrix_name_indices = ({name: idx for idx, name in enumerate(svc.matrix.names())}
                                     if svc.matrix is not None else {})
        sorted_isotopes = sorted(svc.isotopes, key=_sort_hl)

        for iso in sorted_isotopes:
            name = iso.name
            ov = self._user_overrides.get(name, {})

            # Half-life
            hl = getattr(iso, "half_life_s", None)
            _dc = getattr(iso, "decay_constant", 0.0)
            if (hl is None) and _dc and _dc > 0:
                hl = _math.log(2) / _dc
            if hl is None or hl >= 1e30:
                hl_str = "stable"
                lam = 0.0
            else:
                hl_str = _human_time(float(hl))
                lam = _math.log(2) / float(hl) if float(hl) > 0 else 0.0

            lam = float(ov.get("lambda", lam))

            # Decay mode
            decays = iso.getListOfDecays() if hasattr(iso, "getListOfDecays") else []
            decay_str = decays[0].mode if decays else "stable"

            # Channel filters disable transitions, never isotope membership.
            active = [d for d in decays
                      if decay_mode_filter is None
                      or str(d.mode).lower() in decay_mode_filter]
            if svc.config is not None and not svc.config.is_decay_allowed(hl):
                active = []
            decay_str = ", ".join(dict.fromkeys(str(d.mode) for d in active))
            if not decay_str:
                decay_str = "disabled" if decays else "no decay channels"
            lam *= sum(float(getattr(d, "branch", 1.0)) for d in active)
            manual = [link for link in svc.get_manual_decay_links() if link["parent"] == name]
            if manual:
                decay_str = ", ".join(dict.fromkeys(
                    ([str(d.mode) for d in active] +
                     [str(link.get("mode", "")) + " (manual)" for link in manual])))
            if svc.matrix is not None:
                matrix_indices = getattr(self, "_matrix_name_indices", {})
                idx = matrix_indices.get(name)
                if idx is not None:
                    lam = max(0.0, -float(svc.matrix.matrix_decay[idx, idx]))


            # Cross-sections per MT
            sigma_by_mt: Dict[int, float] = {}
            boundary_mts: set = set()  # MTs whose product is outside the burnup network
            open_mts_set = set(open_mts or [])
            source = "—"
            for rx in (iso.getListOfReactions() if hasattr(iso, "getListOfReactions") else []):
                mt = int(rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1))
                sigma = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                sigma_by_mt[mt] = sigma_by_mt.get(mt, 0.0) + sigma  # accumulate isomeric branches
                if getattr(rx, "product", None) is None:
                    boundary_mts.add(mt)
                if sigma > 0 and source == "—":
                    src = getattr(rx, "source", None)
                    source = src or "—"

            # Zero sigma for closed-boundary reactions (product outside network, MT not open)
            for _bmt in boundary_mts:
                if _bmt not in open_mts_set:
                    sigma_by_mt[_bmt] = 0.0

            mts = [102, 16, 103, 107, 18]
            total_rxn_rate = 0.0
            row_cells = []

            def _sig(mt):
                v = float(ov.get(f"sigma_{mt}", sigma_by_mt.get(mt, 0.0)))
                return v

            for mt in mts:
                sig = _sig(mt)
                rate = sig * 1e-24 * flux if sig > 0 else 0.0
                total_rxn_rate += rate
                row_cells.append((_fmt_e(sig) if sig > 0 else "—", _fmt_e(rate) if rate > 0 else "—"))

            total_rate = lam + total_rxn_rate

            cells = [
                (name, False),
                (hl_str, False),
                (_fmt_e(lam) if lam > 0 else "0", True),
                (str(decay_str), False),
            ]
            for sig_str, rate_str in row_cells:
                cells.append((sig_str, True))
                cells.append((rate_str, False))
            cells.append((_fmt_e(total_rate) if total_rate > 0 else "0", False))
            cells.append((source, False))

            row = []
            for val, editable in cells:
                item = QStandardItem(val)
                item.setEditable(editable)
                row.append(item)
            self.unified_model.appendRow(row)

        self._unified_updating = False

    def _on_unified_cell_edited(self, item: "QStandardItem") -> None:
        if self._unified_updating:
            return
        col = item.column()
        row = item.row()
        name_item = self.unified_model.item(row, 0)
        if not name_item:
            return
        name = name_item.text()
        col_keys = {2: "lambda", 4: "sigma_102", 6: "sigma_16", 8: "sigma_103",
                    10: "sigma_107", 12: "sigma_18"}
        if col in col_keys:
            try:
                val = float(item.text())
            except ValueError:
                return
            if name not in self._user_overrides:
                self._user_overrides[name] = {}
            self._user_overrides[name][col_keys[col]] = val

    def populate_cycles(self, cycles: list) -> None:
        """Populate the cycle results sub-table."""
        self.cycles_model.removeRows(0, self.cycles_model.rowCount())
        for i, cyc in enumerate(cycles, 1):
            iso_names = getattr(cyc, "isotope_names", None)
            if iso_names is None:
                iso_names = [s.isotope.name for s in cyc.steps if hasattr(s, "isotope")]
            length = getattr(cyc, "length", len(cyc.steps))
            ctype = "mixed"
            if getattr(cyc, "is_pure_decay", False):
                ctype = "pure decay"
            elif getattr(cyc, "is_pure_reaction", False):
                ctype = "pure reaction"
            row = [
                QStandardItem(str(i)),
                QStandardItem(str(length)),
                QStandardItem(", ".join(iso_names) if iso_names else "—"),
                QStandardItem(ctype),
            ]
            for item in row:
                item.setEditable(False)
            self.cycles_model.appendRow(row)

    # -- Manual decay channels editor ---------------------------------------

    def _mdl_add_row(self):
        row = [
            QStandardItem(""),   # Parent
            QStandardItem(""),   # Daughter
            QStandardItem("beta-"),  # Mode
            QStandardItem(""),   # T½ [s]
            QStandardItem("1.0"),  # Branch
        ]
        self.mdl_model.appendRow(row)
        new_row = self.mdl_model.rowCount() - 1
        self.mdl_table.scrollTo(self.mdl_model.index(new_row, 0))
        self.mdl_table.setCurrentIndex(self.mdl_model.index(new_row, 0))

    def _mdl_del_row(self):
        sel = self.mdl_table.selectionModel().selectedRows()
        for idx in sorted(sel, key=lambda i: i.row(), reverse=True):
            self.mdl_model.removeRow(idx.row())

    def read_manual_decay_links(self) -> list:
        """Read all rows from the manual decay table and return as list of dicts."""
        links = []
        for r in range(self.mdl_model.rowCount()):
            parent   = (self.mdl_model.item(r, 0) or QStandardItem()).text().strip()
            daughter = (self.mdl_model.item(r, 1) or QStandardItem()).text().strip()
            mode     = (self.mdl_model.item(r, 2) or QStandardItem()).text().strip() or "beta-"
            hl_str   = (self.mdl_model.item(r, 3) or QStandardItem()).text().strip()
            br_str   = (self.mdl_model.item(r, 4) or QStandardItem()).text().strip() or "1.0"
            if not parent and not daughter and not hl_str:
                continue
            svc = get_service()
            if not svc.parse_label(parent) or (daughter and not svc.parse_label(daughter)):
                raise ValueError(f"Manual channel row {r + 1}: invalid isotope label.")
            parent = svc.normalize_label(parent)
            daughter = svc.normalize_label(daughter) if daughter else ""
            mode = mode.lower()
            try:
                hl, br = float(hl_str), float(br_str)
            except ValueError:
                raise ValueError(f"Manual channel row {r + 1}: enter numeric half-life and branch.")
            if not np.isfinite(hl) or hl <= 0 or not np.isfinite(br) or not 0 < br <= 1:
                raise ValueError(f"Manual channel row {r + 1}: require half-life > 0 and 0 < branch <= 1.")
            links.append({"parent": parent, "daughter": daughter,
                          "mode": mode, "half_life_s": hl, "branch": br})
        return links

    # -- Evolution population -----------------------------------------------

    def populate_evolution(self, res: RunResult, rho: float) -> None:
        self.exposure_panel.set_run(res)
        self.exposure_panel.session_rho = rho
        self.exposure_panel.session_eq = get_service().last_eq
        self._evo_res = res
        self._evo_rho = rho
        n_steps = len(res.times)
        self.time_slider.setMaximum(max(0, n_steps - 1))
        self.time_slider.setValue(n_steps - 1)
        self._fill_evo_table(res, n_steps - 1, rho)

    def _on_slider_changed(self, idx: int):
        res = self._evo_res
        if res is None:
            return
        self._fill_evo_table(res, idx, self._evo_rho)

    def _fill_evo_table(self, res: RunResult, step: int, rho: float):
        if not res.trajectories or step >= len(res.trajectories):
            return
        N_vec = res.trajectories[step]
        t = res.times[step]
        self.time_label.setText(f"t = {_human_time(t)}")

        sigma_by_name: Dict[str, float] = {}
        for iso in res.isotopes:
            for rx in iso.getListOfReactions():
                mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                if mt == 102:
                    sigma_by_name[iso.name] = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                    break

        total = float(np.sum(N_vec)) or 1.0
        mode = self._current_mode()

        svc = get_service()
        heat_diag = None
        if svc.matrix is not None and svc.matrix.matrix_heat is not None:
            try:
                heat_diag = np.diag(svc.matrix.matrix_heat)
            except Exception:
                pass

        conv = _heat_conv_factor(rho, res.names, res.isotopes, res.initial)

        target_set = set(getattr(res, "targets", []) or [])

        self.evo_model.removeRows(0, self.evo_model.rowCount())
        rows_data = []
        for i, name in enumerate(res.names):
            if target_set and name not in target_set:
                continue   # restricted run: show only requested targets
            N_i = float(N_vec[i])
            if N_i <= 0 and name not in target_set:
                continue   # keep zero-valued targets visible, hide other zeros
            frac = N_i / total
            if mode == "relative":
                display_n = frac
            elif mode == "sigma_n":
                display_n = sigma_by_name.get(name, 0.0) * N_i
            else:
                display_n = N_i

            q_mev = (float(heat_diag[i]) * N_i) if heat_diag is not None else None
            q_wcm3 = (q_mev * conv) if (q_mev is not None and conv is not None) else None
            rows_data.append((name, N_i, frac, display_n, q_mev, q_wcm3))

        # Compute N(A): sum of display_n for all isotopes sharing the same mass number A
        by_A: Dict[int, float] = {}
        for name, N_i, frac, display_n, q_mev, q_wcm3 in rows_data:
            A = _parse_A(name)
            if A is not None:
                by_A[A] = by_A.get(A, 0.0) + display_n

        rows_data.sort(key=lambda r: r[1], reverse=True)
        for name, N_i, frac, display_n, q_mev, q_wcm3 in rows_data:
            A = _parse_A(name)
            na = by_A.get(A) if A is not None else None
            self.evo_model.appendRow([
                QStandardItem(name),
                QStandardItem(_fmt_e(display_n)),
                QStandardItem(f"{frac:.5e}"),
                QStandardItem(_fmt_e(na)),
                QStandardItem(_fmt_e(q_mev)),
                QStandardItem(_fmt_e(q_wcm3)),
            ])

        # Separator + per-A summary block
        if by_A:
            sep = QStandardItem("── N(A) by mass number ──")
            sep.setEnabled(False)
            self.evo_model.appendRow([sep] + [QStandardItem("") for _ in range(5)])
            for A in sorted(by_A):
                na_val = by_A[A]
                self.evo_model.appendRow([
                    QStandardItem(f"A = {A}"),
                    QStandardItem(""),
                    QStandardItem(""),
                    QStandardItem(_fmt_e(na_val)),
                    QStandardItem(""),
                    QStandardItem(""),
                ])

        step_q_mev = res.heat[step] if step < len(res.heat) else None
        step_q_wcm3 = (step_q_mev * conv) if (step_q_mev is not None and conv is not None) else None
        q_str = f"Q = {_fmt_e(step_q_mev)} MeV/s"
        if step_q_wcm3 is not None:
            q_str += f"  =  {_fmt_e(step_q_wcm3)} W/cm³"
        self.total_q_label.setText(q_str)

    def populate_equilibrium(self, eq: EqResult, rho: float) -> None:
        self._eq_res = eq
        self._eq_rho = rho
        conv = _heat_conv_factor(rho, eq.names, eq.isotopes, eq.initial)
        self.eq_model.removeRows(0, self.eq_model.rowCount())
        for flux, t_irr, q_mev in zip(eq.flux_list, eq.time_list, eq.Q_eq):
            q_wcm3 = (q_mev * conv) if (conv is not None) else None
            self.eq_model.appendRow([
                QStandardItem(f"{flux:.4e}"),
                QStandardItem(f"{t_irr:.4e}"),
                QStandardItem(f"{q_mev:.5e}"),
                QStandardItem(_fmt_e(q_wcm3)),
            ])

    def populate_sensitivity(self, res) -> None:
        """Fill the Sensitivity tab from a core.sensitivity.SensitivityResult."""
        self.sens_target_label.setText(
            f"Target: {res.target_desc}    |    base value = {_fmt_e(res.base_value)}")

        self.sens_model.removeRows(0, self.sens_model.rowCount())
        for r in res.ranked():
            self.sens_model.appendRow([
                QStandardItem(str(r.rank)),
                QStandardItem(r.parameter),
                QStandardItem(r.type),
                QStandardItem(_fmt_e(r.base)),
                QStandardItem(_fmt_e(r.perturbed)),
                QStandardItem(f"{r.S:+.4e}"),
            ])

        has_lib = bool(res.library_rows)
        self.sens_lib_caption.setVisible(has_lib)
        self.sens_lib_table.setVisible(has_lib)
        self.sens_lib_model.removeRows(0, self.sens_lib_model.rowCount())
        if has_lib:
            stats = res.lib_stats or {}
            cap = "Library spread (Nᵢ / Q per library)"
            if stats:
                cap += (f"   —   min {_fmt_e(stats.get('min'))}, "
                        f"max {_fmt_e(stats.get('max'))}, "
                        f"mean {_fmt_e(stats.get('mean'))}, "
                        f"σ {_fmt_e(stats.get('std'))}")
            self.sens_lib_caption.setText(cap)
            for lr in res.library_rows:
                self.sens_lib_model.appendRow([
                    QStandardItem(lr.library),
                    QStandardItem(_fmt_e(lr.N_i) if lr.N_i is not None else "—"),
                    QStandardItem(_fmt_e(lr.Q)),
                    QStandardItem(f"{lr.deviation_pct:+.2f}"),
                ])

    def iter_evo_rows(self):
        for r in range(self.evo_model.rowCount()):
            yield [self.evo_model.item(r, c).text() if self.evo_model.item(r, c) else ""
                   for c in range(5)]

    def iter_eq_rows(self):
        for r in range(self.eq_model.rowCount()):
            yield [self.eq_model.item(r, c).text() if self.eq_model.item(r, c) else ""
                   for c in range(4)]


# ---------------------------------------------------------------------------
# Heat conversion helper
# ---------------------------------------------------------------------------

def _heat_conv_factor(rho: float, names: List[str], isotopes: list,
                      initial: Dict[str, float]) -> Optional[float]:
    """Returns f such that Q[W/cm³] = Q[MeV/s/nucleus] * f."""
    if rho <= 0:
        return None
    dominant_name = max(initial.items(), key=lambda kv: kv[1])[0] if initial else None
    M = 1.0
    if dominant_name:
        m = re.search(r"-(\d+)", dominant_name)
        if m:
            M = float(m.group(1))
    N_A = 6.02214076e23
    return 1.602e-13 * rho * N_A / M


# ---------------------------------------------------------------------------
# Right — inputs
# ---------------------------------------------------------------------------

class Right(QVBoxLayout):

    def __init__(self):
        super().__init__()

        # ---- Shared: IC table ----
        self.init_cond_label = QLabel('Initial Conditions')
        self.init_cond_label.setObjectName('init_cond_label')

        self.ic_model = QStandardItemModel()
        self.ic_model.setColumnCount(2)
        self.ic_model.setHorizontalHeaderLabels(['Isotope', 'Concentration'])
        self.ic_model.appendRow([QStandardItem(""), QStandardItem("")])

        self.ic_table = QTableView()
        self.ic_table.setObjectName('init_cond_table')
        self.ic_table.setModel(self.ic_model)
        self.ic_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ic_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.ic_table.setMaximumHeight(150)

        row_btns = QHBoxLayout()
        self.load_matrix_pb = QPushButton('⟳ Load matrix isotopes')
        self.load_matrix_pb.setObjectName('load_matrix_pb')
        self.load_matrix_pb.setToolTip("Fill the table from the current burnup matrix.")
        self.load_matrix_pb.clicked.connect(self._load_from_matrix)
        self.add_row_pb   = QPushButton('+ row')
        self.add_row_pb.setObjectName('add_row_pb')
        self.add_row_pb.clicked.connect(self._add_row)
        self.del_row_pb   = QPushButton('- row')
        self.del_row_pb.setObjectName('del_row_pb')
        self.del_row_pb.clicked.connect(self._del_row)
        self.clear_pb     = QPushButton('clear')
        self.clear_pb.setObjectName('clear_rows_pb')
        self.clear_pb.clicked.connect(self._clear_rows)
        self.load_csv_pb  = QPushButton('load csv')
        self.load_csv_pb.setObjectName('load_csv_pb')
        self.load_csv_pb.clicked.connect(self._load_csv)
        for w in (self.load_matrix_pb, self.add_row_pb, self.del_row_pb,
                  self.clear_pb, self.load_csv_pb):
            row_btns.addWidget(w)

        # ---- Shared: Filters ----
        filter_group = QGroupBox('Filters')
        filter_group.setObjectName('filter_group')
        fg = QVBoxLayout(filter_group)
        fg.setContentsMargins(4, 4, 4, 4)

        hl_row = QHBoxLayout()
        hl_row.addWidget(QLabel('Half-life:'))
        self.hl_filter_box = QComboBox()
        self.hl_filter_box.setObjectName('hl_filter_box')
        self.hl_filter_box.addItems(_HL_PRESETS)
        hl_row.addWidget(self.hl_filter_box)
        fg.addLayout(hl_row)

        mt_row = QHBoxLayout()
        mt_row.addWidget(QLabel('Reactions:'))
        self._mt_checks: Dict[str, QCheckBox] = {}
        for name, mt in _MT_LABELS:
            cb = QCheckBox(name)
            cb.setObjectName(f"mt_{mt}_cb")
            cb.setChecked(True)
            mt_row.addWidget(cb)
            self._mt_checks[name] = cb
        fg.addLayout(mt_row)

        decay_row = QHBoxLayout()
        decay_row.addWidget(QLabel('Decays:'))
        self._decay_checks: Dict[str, QCheckBox] = {}
        for symbol, mode in _DECAY_LABELS:
            cb = QCheckBox(symbol)
            cb.setObjectName(f"decay_{mode.replace('-', '_')}_cb")
            cb.setChecked(True)
            decay_row.addWidget(cb)
            self._decay_checks[mode] = cb
        fg.addLayout(decay_row)

        # ---- Shared parameters (method, energy, density) ----
        param_group = QGroupBox('Parameters')
        param_group.setObjectName('param_group')
        pg = QVBoxLayout(param_group)
        pg.setContentsMargins(4, 4, 4, 4)

        meth_row = QHBoxLayout()
        meth_row.addWidget(QLabel('Method:'))
        self.method_box = QComboBox()
        self.method_box.setObjectName('method_box')
        self.method_box.addItem('CRAM-16', 'cram16')
        self.method_box.addItem('CRAM-16 (adaptive)', 'cram16_adaptive')
        self.method_box.setToolTip('CRAM-16: one rational step per output time.\nCRAM-16 (adaptive): internal step refinement and error checks.')
        self.method_box.addItem('CRAM-48', 'cram48')
        self.method_box.addItem('Padé (expm)', 'pade')
        meth_row.addWidget(self.method_box)
        pg.addLayout(meth_row)

        energy_row = QHBoxLayout()
        energy_row.addWidget(QLabel('Energy [eV]:'))
        self.energy_lineedit = QLineEdit('30000')
        self.energy_lineedit.setObjectName('energy_lineedit')
        self.energy_lineedit.setPlaceholderText('e.g. 30000')
        energy_row.addWidget(self.energy_lineedit)
        pg.addLayout(energy_row)

        rho_row = QHBoxLayout()
        rho_row.addWidget(QLabel('Density ρ [g/cm³]:'))
        self.rho_lineedit = QLineEdit()
        self.rho_lineedit.setObjectName('rho_lineedit')
        self.rho_lineedit.setPlaceholderText('0 → Q in MeV/s only')
        rho_row.addWidget(self.rho_lineedit)
        pg.addLayout(rho_row)

        # ---- Mode tabs ----
        self.mode_tabs = QTabWidget()
        self.mode_tabs.setObjectName('mode_tabs')
        self.mode_tabs.setElideMode(Qt.ElideNone)
        self.mode_tabs.tabBar().setExpanding(False)
        self.mode_tabs.tabBar().setUsesScrollButtons(True)

        # -- Tab 0: Isotopes & Matrix (right controls) --
        iso_ctrl_w = QWidget()
        icv = QVBoxLayout(iso_ctrl_w)
        icv.setContentsMargins(4, 4, 4, 4)

        flux_iso_row = QHBoxLayout()
        flux_iso_row.addWidget(QLabel('Flux Φ [n/cm²/s]:'))
        self.matrix_flux_edit = QLineEdit('0')
        self.matrix_flux_edit.setObjectName('matrix_flux_edit')
        self.matrix_flux_edit.setPlaceholderText('0 = decay only')
        flux_iso_row.addWidget(self.matrix_flux_edit)
        icv.addLayout(flux_iso_row)

        reactor_e_row = QHBoxLayout()
        reactor_e_row.addWidget(QLabel('Reactor energy [eV]:'))
        self.reactor_energy_edit = QLineEdit('30000')
        self.reactor_energy_edit.setObjectName('reactor_energy_edit')
        reactor_e_row.addWidget(self.reactor_energy_edit)
        icv.addLayout(reactor_e_row)

        astro_e_row = QHBoxLayout()
        astro_e_row.addWidget(QLabel('Astro energy [keV]:'))
        self.astro_energy_edit = QLineEdit('30')
        self.astro_energy_edit.setObjectName('astro_energy_edit')
        astro_e_row.addWidget(self.astro_energy_edit)
        icv.addLayout(astro_e_row)

        chan_group = QGroupBox('Open channels')
        chan_group.setObjectName('chan_open_group')
        cg = QVBoxLayout(chan_group)
        cg.setContentsMargins(4, 4, 4, 4)
        self._chan_checks: Dict[int, QCheckBox] = {}
        for label, mt in _MT_LABELS:
            cb = QCheckBox(f"{label}  open")
            cb.setObjectName(f"chan_{mt}_cb")
            cb.setChecked(True)
            cg.addWidget(cb)
            self._chan_checks[mt] = cb
        icv.addWidget(chan_group)

        dec_lib_group = QGroupBox('Decay library')
        dec_lib_group.setObjectName('dec_lib_group')
        dlg = QVBoxLayout(dec_lib_group)
        dlg.setContentsMargins(4, 4, 4, 4)
        self.dec_lib_combo = QComboBox()
        self.dec_lib_combo.setObjectName('dec_lib_combo')
        self.dec_lib_apply_pb = QPushButton('Apply')
        self.dec_lib_apply_pb.setObjectName('dec_lib_apply_pb')
        self.dec_lib_apply_pb.clicked.connect(self._on_switch_library)
        dlg.addWidget(self.dec_lib_combo)
        dlg.addWidget(self.dec_lib_apply_pb)
        icv.addWidget(dec_lib_group)

        self.refresh_matrix_pb = QPushButton('Refresh Matrix View')
        self.refresh_matrix_pb.setObjectName('refresh_matrix_pb')
        icv.addWidget(self.refresh_matrix_pb, alignment=Qt.AlignCenter)
        icv.addStretch(1)

        # -- Tab 1: Evolution --
        evo_w = QWidget()
        ev = QVBoxLayout(evo_w)
        ev.setContentsMargins(4, 4, 4, 4)

        tstart_row = QHBoxLayout()
        tstart_row.addWidget(QLabel('Time start:'))
        self.t_start_edit = QLineEdit('0')
        self.t_start_edit.setObjectName('t_start_edit')
        self.t_start_unit = QComboBox()
        self.t_start_unit.setObjectName('t_start_unit')
        self.t_start_unit.addItems([u for u, _ in _TIME_UNITS])
        tstart_row.addWidget(self.t_start_edit)
        tstart_row.addWidget(self.t_start_unit)
        ev.addLayout(tstart_row)

        tend_row = QHBoxLayout()
        tend_row.addWidget(QLabel('Time end:'))
        self.t_end_edit = QLineEdit('1e15')
        self.t_end_edit.setObjectName('t_end_edit')
        self.t_end_edit.setPlaceholderText('e.g. 3600')
        self.t_end_unit = QComboBox()
        self.t_end_unit.setObjectName('t_end_unit')
        self.t_end_unit.addItems([u for u, _ in _TIME_UNITS])
        tend_row.addWidget(self.t_end_edit)
        tend_row.addWidget(self.t_end_unit)
        ev.addLayout(tend_row)

        # Neutron input mode: flux Φ [n/cm²/s] or neutron density n [n/cm³]
        flux_mode_row = QHBoxLayout()
        self.flux_phi_rb = QRadioButton('Flux  Φ  [n/cm²/s]')
        self.flux_phi_rb.setObjectName('flux_phi_rb')
        self.flux_phi_rb.setChecked(True)
        self.flux_n_rb = QRadioButton('Neutron density  n  [n/cm³]')
        self.flux_n_rb.setObjectName('flux_n_rb')
        flux_mode_row.addWidget(self.flux_phi_rb)
        flux_mode_row.addWidget(self.flux_n_rb)
        ev.addLayout(flux_mode_row)

        flux_row = QHBoxLayout()
        self.flux_mode_label = QLabel('Φ [n/cm²/s]:')
        self.flux_mode_label.setObjectName('flux_mode_label')
        flux_row.addWidget(self.flux_mode_label)
        self.flux_lineedit = QLineEdit()
        self.flux_lineedit.setObjectName('flux_lineedit')
        self.flux_lineedit.setPlaceholderText('0 = decay only')
        flux_row.addWidget(self.flux_lineedit)
        ev.addLayout(flux_row)

        self.flux_phi_rb.toggled.connect(self._update_flux_label)
        self.flux_n_rb.toggled.connect(self._update_flux_label)

        # Target isotopes: when set, only the minimal source→target sub-network
        # is solved (much faster, identical target results) and only these
        # isotopes are shown. Empty → full network, every isotope shown.
        ev.addWidget(QLabel('Target isotopes (optional):'))
        self.targets_lineedit = QLineEdit()
        self.targets_lineedit.setObjectName('targets_lineedit')
        self.targets_lineedit.setPlaceholderText('e.g.  Po-210 Bi-209   (empty = all)')
        self.targets_lineedit.setToolTip(
            "Restrict the run to the minimal precursor sub-network that feeds "
            "these isotopes from the initial inventory. Faster, with identical "
            "target concentrations. Leave empty to compute the whole network.")
        ev.addWidget(self.targets_lineedit)

        self.calculate_pb = QPushButton('Calculate Evolution')
        self.calculate_pb.setObjectName('calculate_pb')
        ev.addWidget(self.calculate_pb, alignment=Qt.AlignCenter)
        ev.addStretch(1)   # pack controls to the top (no spread-out gaps)

        # -- Tab 2: Equilibrium --
        eq_w = QWidget()
        eqv = QVBoxLayout(eq_w)
        eqv.setContentsMargins(4, 4, 4, 4)

        # Mode toggle: auto-plateau vs manual
        eq_mode_row = QHBoxLayout()
        eq_mode_row.addWidget(QLabel('Mode:'))
        self.eq_auto_rb   = QRadioButton('Auto plateau')
        self.eq_manual_rb = QRadioButton('Manual times')
        self.eq_auto_rb.setObjectName('eq_auto_rb')
        self.eq_manual_rb.setObjectName('eq_manual_rb')
        self.eq_auto_rb.setChecked(True)
        eq_mode_row.addWidget(self.eq_auto_rb)
        eq_mode_row.addWidget(self.eq_manual_rb)
        eqv.addLayout(eq_mode_row)
        self.eq_auto_rb.toggled.connect(self._on_eq_mode_changed)

        eqv.addWidget(QLabel('Flux list [n/cm²/s]:'))
        self.eq_flux_edit = QLineEdit()
        self.eq_flux_edit.setObjectName('eq_flux_edit')
        self.eq_flux_edit.setPlaceholderText('e.g.  1e12 1e13 1e14')
        eqv.addWidget(self.eq_flux_edit)

        # Manual times (hidden in auto mode)
        self.eq_manual_widget = QWidget()
        emv = QVBoxLayout(self.eq_manual_widget)
        emv.setContentsMargins(0, 0, 0, 0)
        emv.addWidget(QLabel('Irradiation times:'))
        eq_time_row = QHBoxLayout()
        self.eq_time_edit = QLineEdit()
        self.eq_time_edit.setObjectName('eq_time_edit')
        self.eq_time_edit.setPlaceholderText('e.g.  1e7 1e8 1e9')
        self.eq_time_unit = QComboBox()
        self.eq_time_unit.setObjectName('eq_time_unit')
        self.eq_time_unit.addItems([u for u, _ in _TIME_UNITS])
        eq_time_row.addWidget(self.eq_time_edit, 1)
        eq_time_row.addWidget(self.eq_time_unit)
        emv.addLayout(eq_time_row)
        eqv.addWidget(self.eq_manual_widget)
        self.eq_manual_widget.setVisible(False)

        eq_note = QLabel("Auto plateau: uses solve_equilibrium per flux point to find steady state.")
        eq_note.setObjectName("eq_auto_note")
        eq_note.setWordWrap(True)
        eqv.addWidget(eq_note)
        self.eq_auto_note = eq_note

        self.calc_eq_pb = QPushButton('Calculate Equilibrium')
        self.calc_eq_pb.setObjectName('calc_eq_pb')
        eqv.addWidget(self.calc_eq_pb, alignment=Qt.AlignCenter)
        eqv.addStretch(1)   # pack controls to the top (no spread-out gaps)

        # -- Tab 3: Sensitivity (OAT) --
        sens_w = self._build_sensitivity_tab()
        sens_scroll = QScrollArea()
        sens_scroll.setObjectName('sens_scroll')
        sens_scroll.setWidgetResizable(True)
        sens_scroll.setWidget(sens_w)

        self.mode_tabs.addTab(iso_ctrl_w, 'Isotopes & Matrix')
        self.mode_tabs.addTab(evo_w, 'Evolution')
        self.mode_tabs.addTab(eq_w, 'Equilibrium')
        self.mode_tabs.addTab(sens_scroll, 'Sensitivity')

        # -- Assemble inside a scroll area so the panel works on small screens --
        inner_w = QWidget()
        inner_v = QVBoxLayout(inner_w)
        inner_v.setContentsMargins(4, 4, 4, 4)
        inner_v.setSpacing(4)
        inner_v.addWidget(self.init_cond_label)
        inner_v.addWidget(self.ic_table)
        inner_v.addLayout(row_btns)
        inner_v.addWidget(filter_group)
        inner_v.addWidget(param_group)
        inner_v.addWidget(self.mode_tabs)
        inner_v.addStretch(1)

        right_scroll = QScrollArea()
        right_scroll.setObjectName('right_panel_scroll')
        right_scroll.setWidgetResizable(True)
        right_scroll.setWidget(inner_w)
        right_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.addWidget(right_scroll)

        # Auto-populate if matrix already built
        svc = get_service()
        self._populate_dec_lib_combo(svc)
        self._populate_sens_lists(svc)
        if svc.initialized and svc.isotopes:
            self._load_from_matrix()

    # -- equilibrium mode toggle -------------------------------------------

    def _on_eq_mode_changed(self, auto_checked: bool):
        self.eq_manual_widget.setVisible(not auto_checked)
        self.eq_auto_note.setVisible(auto_checked)

    # -- decay library selector -------------------------------------------

    def _populate_dec_lib_combo(self, svc=None):
        if svc is None:
            svc = get_service()
        libs = svc.decay_libraries() if svc.initialized else []
        self.dec_lib_combo.clear()
        for lib in libs:
            self.dec_lib_combo.addItem(lib)
        current = getattr(svc, "decay_library", None)
        if current and self.dec_lib_combo.findText(current) >= 0:
            self.dec_lib_combo.setCurrentText(current)

    def _on_switch_library(self):
        lib = self.dec_lib_combo.currentText().strip()
        if not lib:
            return
        svc = get_service()
        try:
            svc.switch_decay_library(lib)
            QMessageBox.information(None, "Decay library",
                                    f"Switched to: {lib}\nRebuild the matrix to apply.")
        except Exception as exc:
            QMessageBox.critical(None, "Decay library error", str(exc))

    # -- open channels for matrix view ------------------------------------

    def get_matrix_open_mts(self) -> List[int]:
        return [mt for mt, cb in self._chan_checks.items() if cb.isChecked()]

    def read_targets(self) -> List[str]:
        """Target isotope labels from the Evolution tab (space/comma separated).
        Empty list → no restriction (solve the whole network)."""
        raw = self.targets_lineedit.text().strip()
        if not raw:
            return []
        return [tok for tok in raw.replace(",", " ").split() if tok]

    # -- Sensitivity (OAT) tab --------------------------------------------

    _SENS_MACS_SOURCES = ['eaf2010', 'talys', 'endfb71', 'rawmacs']

    def _build_sensitivity_tab(self) -> QWidget:
        """Construct the Sensitivity controls (returned widget goes into a tab)."""
        w = QWidget()
        sv = QVBoxLayout(w)
        sv.setContentsMargins(4, 4, 4, 4)

        # Target -----------------------------------------------------------
        tgt_row = QHBoxLayout()
        tgt_row.addWidget(QLabel('Target:'))
        self.sens_target_kind = QComboBox()
        self.sens_target_kind.setObjectName('sens_target_kind')
        self.sens_target_kind.addItem('Concentration Nᵢ', 'concentration')
        self.sens_target_kind.addItem('Heat Q', 'heat')
        tgt_row.addWidget(self.sens_target_kind, 1)
        sv.addLayout(tgt_row)

        iso_row = QHBoxLayout()
        iso_row.addWidget(QLabel('Isotope:'))
        self.sens_target_iso = QComboBox()
        self.sens_target_iso.setObjectName('sens_target_iso')
        self.sens_target_iso.setEditable(True)
        iso_row.addWidget(self.sens_target_iso, 1)
        sv.addLayout(iso_row)

        # Evaluate at ------------------------------------------------------
        eval_row = QHBoxLayout()
        eval_row.addWidget(QLabel('Evaluate at:'))
        self.sens_eval_mode = QComboBox()
        self.sens_eval_mode.setObjectName('sens_eval_mode')
        self.sens_eval_mode.addItem('Time t', 'time')
        self.sens_eval_mode.addItem('Plateau', 'plateau')
        eval_row.addWidget(self.sens_eval_mode, 1)
        self.sens_time_edit = QLineEdit('1e9')
        self.sens_time_edit.setObjectName('sens_time_edit')
        self.sens_time_edit.setPlaceholderText('t [s]')
        eval_row.addWidget(self.sens_time_edit)
        sv.addLayout(eval_row)

        # Flux + σ-mode ----------------------------------------------------
        sflux_row = QHBoxLayout()
        sflux_row.addWidget(QLabel('Flux Φ [n/cm²/s]:'))
        self.sens_flux_edit = QLineEdit('1e15')
        self.sens_flux_edit.setObjectName('sens_flux_edit')
        sflux_row.addWidget(self.sens_flux_edit)
        sv.addLayout(sflux_row)

        smode_row = QHBoxLayout()
        smode_row.addWidget(QLabel('σ mode:'))
        self.sens_sigma_mode = QComboBox()
        self.sens_sigma_mode.setObjectName('sens_sigma_mode')
        self.sens_sigma_mode.addItem('Reactor (point-wise)', False)
        self.sens_sigma_mode.addItem('Astro (MACS)', True)
        smode_row.addWidget(self.sens_sigma_mode, 1)
        sv.addLayout(smode_row)

        # Types ------------------------------------------------------------
        types_group = QGroupBox('Sensitivity types')
        types_group.setObjectName('sens_types_group')
        tgl = QVBoxLayout(types_group)
        tgl.setContentsMargins(4, 4, 4, 4)
        self.sens_type_sigma = QCheckBox('σ  cross-section')
        self.sens_type_sigma.setChecked(True)
        self.sens_type_lambda = QCheckBox('λ  half-life')
        self.sens_type_flux = QCheckBox('Flux / Energy')
        self.sens_type_lib = QCheckBox('Library')
        for cb, name in ((self.sens_type_sigma, 'sens_type_sigma'),
                         (self.sens_type_lambda, 'sens_type_lambda'),
                         (self.sens_type_flux, 'sens_type_flux'),
                         (self.sens_type_lib, 'sens_type_lib')):
            cb.setObjectName(name)
            tgl.addWidget(cb)
        sv.addWidget(types_group)

        # Perturbation -----------------------------------------------------
        delta_row = QHBoxLayout()
        delta_row.addWidget(QLabel('δ (perturbation):'))
        self.sens_delta_edit = QLineEdit('0.01')
        self.sens_delta_edit.setObjectName('sens_delta_edit')
        delta_row.addWidget(self.sens_delta_edit)
        sv.addLayout(delta_row)
        self.sens_two_sided = QCheckBox('Two-sided (±δ central difference)')
        self.sens_two_sided.setObjectName('sens_two_sided')
        self.sens_two_sided.setChecked(True)
        sv.addWidget(self.sens_two_sided)

        # Libraries (for the Library type) ---------------------------------
        lib_group = QGroupBox('Libraries (Library type)')
        lib_group.setObjectName('sens_lib_group')
        lgl = QVBoxLayout(lib_group)
        lgl.setContentsMargins(4, 4, 4, 4)
        self.sens_lib_list = QListWidget()
        self.sens_lib_list.setObjectName('sens_lib_list')
        self.sens_lib_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.sens_lib_list.setMaximumHeight(110)
        lgl.addWidget(self.sens_lib_list)
        ref_row = QHBoxLayout()
        ref_row.addWidget(QLabel('Reference:'))
        self.sens_ref_combo = QComboBox()
        self.sens_ref_combo.setObjectName('sens_ref_combo')
        ref_row.addWidget(self.sens_ref_combo, 1)
        lgl.addLayout(ref_row)
        sv.addWidget(lib_group)

        self.sens_calc_pb = QPushButton('Run Sensitivity')
        self.sens_calc_pb.setObjectName('sens_calc_pb')
        sv.addWidget(self.sens_calc_pb, alignment=Qt.AlignCenter)
        sv.addStretch(1)
        return w

    def _populate_sens_lists(self, svc=None):
        """Fill the target-isotope, library and reference selectors from the
        current universe (called on build and whenever the matrix changes)."""
        if svc is None:
            svc = get_service()
        labels = svc.universe_labels() if svc.initialized else []
        cur_iso = self.sens_target_iso.currentText()
        self.sens_target_iso.clear()
        self.sens_target_iso.addItems(labels)
        if cur_iso:
            self.sens_target_iso.setEditText(cur_iso)

        decay = svc.decay_libraries() if svc.initialized else []
        libs = decay + self._SENS_MACS_SOURCES
        self.sens_lib_list.clear()
        for lib in libs:
            self.sens_lib_list.addItem(QListWidgetItem(lib))
        self.sens_ref_combo.clear()
        self.sens_ref_combo.addItem('Mean (default)', None)
        for lib in libs:
            self.sens_ref_combo.addItem(lib, lib)

    def read_sensitivity_config(self, initial: Dict[str, float]) -> SensitivityConfig:
        """Assemble a SensitivityConfig from the Sensitivity-tab widgets."""
        def _f(text, default):
            try:
                return float(str(text).strip())
            except (ValueError, TypeError):
                return default

        types: List[str] = []
        if self.sens_type_sigma.isChecked():
            types.append('sigma')
        if self.sens_type_lambda.isChecked():
            types.append('halflife')
        if self.sens_type_flux.isChecked():
            types.append('flux_energy')
        if self.sens_type_lib.isChecked():
            types.append('library')

        iso = self.sens_target_iso.currentText().strip()
        return SensitivityConfig(
            initial_by_label=initial,
            target_kind=self.sens_target_kind.currentData(),
            target_label=(iso or None),
            eval_mode=self.sens_eval_mode.currentData(),
            t_eval_s=_f(self.sens_time_edit.text(), 1e9),
            types=types,
            delta=_f(self.sens_delta_edit.text(), 0.01),
            two_sided=self.sens_two_sided.isChecked(),
            flux=_f(self.sens_flux_edit.text(), 0.0),
            energy_eV=self.read_energy_eV() or 30_000.0,
            mt_filter=self.get_mt_filter(),
            method=self.read_method(),
            hl_preset=self.hl_filter_box.currentText(),
            open_boundary_mts=self.get_matrix_open_mts(),
            prefer_macs=bool(self.sens_sigma_mode.currentData()),
            libraries=[i.text() for i in self.sens_lib_list.selectedItems()],
            reference_library=self.sens_ref_combo.currentData(),
        )

    # -- IC table helpers --------------------------------------------------

    def _load_from_matrix(self):
        svc = get_service()
        if not svc.initialized:
            QMessageBox.warning(None, "Initial conditions",
                                "Nuclear-data environment is not initialised.")
            return
        labels = svc.universe_labels() or [iso.name for iso in svc.get_universe()]
        if not labels:
            QMessageBox.information(
                None, "Initial conditions",
                "No burnup matrix yet — build one in the Isotope-chart panel.")
            return
        self.ic_model.removeRows(0, self.ic_model.rowCount())
        for lab in labels:
            iso_item = QStandardItem(lab)
            iso_item.setEditable(False)
            self.ic_model.appendRow([iso_item, QStandardItem("0")])
        self._populate_sens_lists(svc)   # refresh sensitivity selectors too

    def _add_row(self):
        self.ic_model.appendRow([QStandardItem(""), QStandardItem("")])

    def _del_row(self):
        idx = self.ic_table.currentIndex()
        row = idx.row() if idx.isValid() else self.ic_model.rowCount() - 1
        if 0 <= row < self.ic_model.rowCount():
            self.ic_model.removeRow(row)

    def _clear_rows(self):
        self.ic_model.removeRows(0, self.ic_model.rowCount())
        self.ic_model.appendRow([QStandardItem(""), QStandardItem("")])

    def _load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            None, "Load initial conditions (CSV)", "", "CSV (*.csv)")
        if not path:
            return
        try:
            with open(path, "r") as f:
                reader = csv.reader(f)
                rows: List[Tuple[str, str]] = []
                for row in reader:
                    if not row:
                        continue
                    label = row[0].strip()
                    conc  = row[1].strip() if len(row) > 1 else ""
                    if not label or label.lower() in ("isotope", "iso", "name"):
                        continue
                    rows.append((label, conc))
        except OSError as e:
            QMessageBox.critical(None, "CSV load failed", str(e))
            return
        self.ic_model.removeRows(0, self.ic_model.rowCount())
        for label, conc in rows:
            self.ic_model.appendRow([QStandardItem(label), QStandardItem(conc)])
        if not rows:
            self.ic_model.appendRow([QStandardItem(""), QStandardItem("")])

    # -- parsing -----------------------------------------------------------

    def read_initial_conditions(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for r in range(self.ic_model.rowCount()):
            label_item = self.ic_model.item(r, 0)
            conc_item  = self.ic_model.item(r, 1)
            label = label_item.text().strip() if label_item else ""
            conc_text = conc_item.text().strip() if conc_item else ""
            if not label or not conc_text:
                continue
            try:
                conc = float(conc_text)
            except ValueError:
                raise ValueError(f"Row {r + 1}: concentration '{conc_text}' is not a number.")
            if conc == 0.0:
                continue
            out[label] = out.get(label, 0.0) + conc
        return out

    def read_time_start_s(self) -> float:
        text = self.t_start_edit.text().strip() or "0"
        try:
            t = float(text)
        except ValueError:
            raise ValueError(f"Time start '{text}' is not a number.")
        return max(0.0, t) * _unit_seconds(self.t_start_unit.currentText())

    def read_time_end_s(self) -> float:
        text = self.t_end_edit.text().strip()
        if not text:
            return 0.0
        try:
            t = float(text)
        except ValueError:
            raise ValueError(f"Time end '{text}' is not a number.")
        if t < 0:
            raise ValueError("Time must be ≥ 0.")
        return t * _unit_seconds(self.t_end_unit.currentText())

    def _update_flux_label(self):
        if self.flux_phi_rb.isChecked():
            self.flux_mode_label.setText('Φ [n/cm²/s]:')
            self.flux_lineedit.setPlaceholderText('0 = decay only')
        else:
            self.flux_mode_label.setText('n [n/cm³]:')
            self.flux_lineedit.setPlaceholderText('0 = decay only')

    def read_flux(self) -> float:
        import math
        text = self.flux_lineedit.text().strip()
        if not text:
            return 0.0
        try:
            val = float(text)
        except ValueError:
            raise ValueError(f"Flux/density '{text}' is not a number.")
        if val < 0:
            raise ValueError("Flux/density must be ≥ 0.")
        if getattr(self, 'flux_n_rb', None) and self.flux_n_rb.isChecked():
            # Convert neutron density n [n/cm³] → flux Φ [n/cm²/s]
            # Φ = n · v,  v = sqrt(2·E / m_n)
            try:
                e_eV = float(self.energy_lineedit.text().strip() or '30000')
            except (ValueError, AttributeError):
                e_eV = 30_000.0
            E_J   = e_eV * 1.602176634e-19       # J
            m_n   = 1.67492749804e-27             # kg
            v_cms = math.sqrt(2.0 * E_J / m_n) * 100.0  # cm/s
            return val * v_cms
        return val

    def read_energy_eV(self) -> Optional[float]:
        text = self.energy_lineedit.text().strip()
        if not text:
            return None
        try:
            e = float(text)
        except ValueError:
            raise ValueError(f"Energy '{text}' is not a number.")
        if e <= 0:
            raise ValueError("Energy must be > 0.")
        return e

    def read_method(self) -> str:
        return self.method_box.currentData() or "cram16"

    def read_density(self) -> float:
        text = self.rho_lineedit.text().strip()
        if not text:
            return 0.0
        try:
            return float(text)
        except ValueError:
            return 0.0

    def get_mt_filter(self) -> List[int]:
        mts = []
        for name, cb in self._mt_checks.items():
            if cb.isChecked():
                for n2, mt in _MT_LABELS:
                    if n2 == name:
                        mts.append(mt)
        return mts if mts else [102]

    def get_decay_mode_filter(self) -> Optional[List[str]]:
        checked = [mode for mode, cb in self._decay_checks.items() if cb.isChecked()]
        all_modes = [mode for _, mode in _DECAY_LABELS]
        if set(checked) == set(all_modes):
            return None  # all selected → no filter (faster)
        return checked if checked else []

    def read_eq_flux_list(self) -> List[float]:
        text = self.eq_flux_edit.text().strip()
        if not text:
            return []
        try:
            return [float(v) for v in text.split()]
        except ValueError:
            raise ValueError("Flux list must be space-separated numbers (e.g. 1e12 1e13 1e14).")

    def read_eq_time_list(self) -> List[float]:
        text = self.eq_time_edit.text().strip()
        if not text:
            return []
        unit_s = _unit_seconds(self.eq_time_unit.currentText())
        try:
            return [float(v) * unit_s for v in text.split()]
        except ValueError:
            raise ValueError("Time list must be space-separated numbers.")
