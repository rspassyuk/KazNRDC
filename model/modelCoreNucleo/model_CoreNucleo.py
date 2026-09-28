# Qt________________________________________
from PyQt5.QtCore import Qt, QObject, QThread, QSize, pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QComboBox, QTextEdit, QCheckBox, QGroupBox,
    QSizePolicy, QButtonGroup, QProgressBar, QMessageBox,
    QFileDialog, QSplitter
)
from PyQt5.QtGui import QFont

# OTHERS____________________________________
import traceback
from typing import Optional, List, Dict, Union

# internal GUI comps (for styling objectNames etc.)
from GUIcomp import *  # noqa: F401,F403

# nuclear core — ALL data access goes through the shared service facade.
from core.service import get_service, decay_symbol, mt_symbol

Meta = Optional[Union[int, str]]


# =========================
# Small formatting helpers (display only)
# =========================
def _format_isotope_label(X: str, A: int, Meta: Meta) -> str:
    if Meta is None or Meta == "" or Meta == 0:
        return f"{X}-{A}"
    if isinstance(Meta, int):
        return f"{X}-{A}m{Meta}"
    return f"{X}-{A}{str(Meta).strip()}"


def _human_time(t: Optional[float]) -> str:
    if t is None:
        return "stable"
    for div, unit in [(3.15576e7, "yr"), (86400.0, "d"), (3600.0, "h"),
                      (60.0, "min"), (1.0, "s"), (1e-3, "ms"), (1e-6, "µs")]:
        if t >= div:
            return f"{t / div:.3g} {unit}"
    return f"{t:.3g} s"


def _fmt_hl(t: Optional[float]) -> str:
    if t is None:
        return "stable"
    return f"{t:.4e} s   ({_human_time(t)})"


def _fmt_sigma(v) -> str:
    return "—" if v is None else f"{v:.3e} b"


def _fmt_mb(v) -> str:
    return "—" if v is None else f"{v:.4g} mb"


def _fmt_mev(v) -> str:
    return "—" if v is None else f"{v:.4g} MeV"


# =========================
# UI Overlay (blocking loader while ENDF is read on a worker thread)
# =========================
class LoadingOverlay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background-color: rgba(0, 0, 0, 120);")

        layout = QVBoxLayout(self)
        self.loader = QProgressBar()
        self.loader.setRange(0, 0)  # indeterminate
        self.loader.setTextVisible(False)
        self.loader.setFixedWidth(200)
        self.loader.setStyleSheet(
            "QProgressBar { background-color: white; border: 1px solid grey; border-radius: 3px; }")
        self.label = QLabel("Reading nuclear data…")
        self.label.setStyleSheet(
            "color: white; font-size: 14px; font-weight: bold; background: transparent;")
        layout.addWidget(self.label, alignment=Qt.AlignCenter)
        layout.addWidget(self.loader, alignment=Qt.AlignCenter)

    def mousePressEvent(self, event):
        event.accept()


class RightContainer(QWidget):
    def __init__(self):
        super().__init__()
        self.overlay = LoadingOverlay(self)
        self.overlay.hide()

    def resizeEvent(self, event):
        self.overlay.resize(event.size())
        super().resizeEvent(event)


# =========================
# Worker thread: gather decay + neutron + all-MACS for one isotope
# via the core service (no ENDF logic lives in the model).
# =========================
class IsotopeLoadWorker(QObject):
    finished = pyqtSignal(dict)
    failed = pyqtSignal(str)

    def __init__(self, libs: List[str], Z: int, X: str, A: int, meta: Meta):
        super().__init__()
        self.libs, self.Z, self.X, self.A, self.meta = libs, Z, X, A, meta

    def run(self):
        try:
            svc = get_service()
            decay_by_lib: Dict[str, dict] = {}
            neutron_by_lib: Dict[str, list] = {}
            for lib in self.libs:
                d = svc.read_decay(lib, self.Z, self.X, self.A, self.meta)
                if d:
                    decay_by_lib[lib] = d
                n = svc.read_neutron(lib, self.Z, self.X, self.A, self.meta)
                neutron_by_lib[lib] = n  # may be [] (e.g. JEFF has no n-files)
            macs = svc.read_macs_sources(self.Z, self.X, self.A, self.meta)
            self.finished.emit({
                "Z": self.Z, "X": self.X, "A": self.A, "meta": self.meta,
                "libs": list(self.libs),
                "decay_by_lib": decay_by_lib,
                "neutron_by_lib": neutron_by_lib,
                "macs": macs,
            })
        except Exception:
            self.failed.emit(traceback.format_exc())


# =========================
# Model class (loaded by core.py)
# =========================
class model_CoreNucleo(object):
    """
    Left: big periodic table.
    Right: database selector + isotope selector + Decaydata / NReactiondata.
    """

    def __init__(self, parent=None):
        super().__init__()
        self.main_parent = parent
        self.win_w = parent.window_width
        self.win_h = parent.window_height

        self.__left = Left(self)
        self.__right = Right(self)

        # Left -> Right
        self.__left.elementSelected.connect(self.__right.on_element_selected)

    def getLeftWidget(self):
        return self.__left

    def getRightLayout(self):
        return self.__right

    def getMainParent(self):
        return self.main_parent


# =========================
# Left side: periodic table (big)
# =========================
class Left(QWidget):
    elementSelected = pyqtSignal(int, str)  # Z, X

    def __init__(self, core_model):
        super().__init__()
        self.core_model = core_model

        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(5)
        # Don't let the (fixed-size) grid pin a large minimum on this widget —
        # the table must be free to shrink as well as grow.
        root.setSizeConstraint(QVBoxLayout.SetNoConstraint)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.pt_widget = QWidget()
        self.pt_grid = QGridLayout(self.pt_widget)
        self.pt_grid.setSpacing(2)
        self.pt_grid.setContentsMargins(0, 0, 0, 0)

        self.element_group = QButtonGroup(self)
        self.element_group.setExclusive(True)

        self.element_buttons: Dict[str, QPushButton] = {}
        for (sym, row, col, Z, category) in self._periodic_layout():
            b = QPushButton(sym)
            b.setObjectName(category)
            b.setFixedSize(40, 40)              # square; recomputed in resizeEvent

            b.setCheckable(True)
            self.element_group.addButton(b)

            b.clicked.connect(lambda _, z=Z, x=sym: self.elementSelected.emit(int(z), str(x)))

            self.pt_grid.addWidget(b, row, col)
            self.element_buttons[sym] = b

        # Center the square-celled grid; the cell side is recomputed on resize
        # so the table stays square and scales proportionally at any size.
        root.addStretch(1)
        root.addWidget(self.pt_widget, 0, Qt.AlignHCenter)
        root.addStretch(1)
        self._cell = None

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._apply_cell_size()

    def _apply_cell_size(self):
        """Uniform square cell = min(availW/cols, availH/rows) (Qt has no
        native aspect-ratio for layouts)."""
        cols = max(1, self.pt_grid.columnCount())
        rows = max(1, self.pt_grid.rowCount())
        spacing = self.pt_grid.spacing()
        margins = self.layout().contentsMargins()
        avail_w = self.width() - margins.left() - margins.right()
        avail_h = self.height() - margins.top() - margins.bottom()
        if avail_w <= 0 or avail_h <= 0:
            return
        cell = int(max(8, min((avail_w - (cols - 1) * spacing) / cols,
                              (avail_h - (rows - 1) * spacing) / rows)))
        if cell == self._cell:
            return
        self._cell = cell
        font = QFont("Arial")
        font.setPixelSize(max(6, int(cell * 0.34)))
        for b in self.element_buttons.values():
            b.setFixedSize(cell, cell)
            b.setFont(font)

    # Decoupled from the grid's content size so the page can shrink freely;
    # the cell side is always recomputed from the actual size in resizeEvent.
    def minimumSizeHint(self):
        return QSize(160, 100)

    def sizeHint(self):
        return QSize(760, 400)

    def _periodic_layout(self):
        return [
            ("H", 0, 0, 1, "nonmetal"), ("He", 0, 17, 2, "noble_gas"),
            ("Li", 1, 0, 3, "alkali_metal"), ("Be", 1, 1, 4, "alkaline_earth_metal"),
            ("B", 1, 12, 5, "metalloid"), ("C", 1, 13, 6, "nonmetal"),
            ("N", 1, 14, 7, "nonmetal"), ("O", 1, 15, 8, "nonmetal"),
            ("F", 1, 16, 9, "halogen"), ("Ne", 1, 17, 10, "noble_gas"),
            ("Na", 2, 0, 11, "alkali_metal"), ("Mg", 2, 1, 12, "alkaline_earth_metal"),
            ("Al", 2, 12, 13, "post_transition_metal"), ("Si", 2, 13, 14, "metalloid"),
            ("P", 2, 14, 15, "nonmetal"), ("S", 2, 15, 16, "nonmetal"),
            ("Cl", 2, 16, 17, "halogen"), ("Ar", 2, 17, 18, "noble_gas"),
            ("K", 3, 0, 19, "alkali_metal"), ("Ca", 3, 1, 20, "alkaline_earth_metal"),
            ("Sc", 3, 2, 21, "transition_metal"), ("Ti", 3, 3, 22, "transition_metal"),
            ("V", 3, 4, 23, "transition_metal"), ("Cr", 3, 5, 24, "transition_metal"),
            ("Mn", 3, 6, 25, "transition_metal"), ("Fe", 3, 7, 26, "transition_metal"),
            ("Co", 3, 8, 27, "transition_metal"), ("Ni", 3, 9, 28, "transition_metal"),
            ("Cu", 3, 10, 29, "transition_metal"), ("Zn", 3, 11, 30, "transition_metal"),
            ("Ga", 3, 12, 31, "post_transition_metal"), ("Ge", 3, 13, 32, "metalloid"),
            ("As", 3, 14, 33, "metalloid"), ("Se", 3, 15, 34, "nonmetal"),
            ("Br", 3, 16, 35, "halogen"), ("Kr", 3, 17, 36, "noble_gas"),
            ("Rb", 4, 0, 37, "alkali_metal"), ("Sr", 4, 1, 38, "alkaline_earth_metal"),
            ("Y", 4, 2, 39, "transition_metal"), ("Zr", 4, 3, 40, "transition_metal"),
            ("Nb", 4, 4, 41, "transition_metal"), ("Mo", 4, 5, 42, "transition_metal"),
            ("Tc", 4, 6, 43, "transition_metal"), ("Ru", 4, 7, 44, "transition_metal"),
            ("Rh", 4, 8, 45, "transition_metal"), ("Pd", 4, 9, 46, "transition_metal"),
            ("Ag", 4, 10, 47, "transition_metal"), ("Cd", 4, 11, 48, "transition_metal"),
            ("In", 4, 12, 49, "post_transition_metal"), ("Sn", 4, 13, 50, "post_transition_metal"),
            ("Sb", 4, 14, 51, "metalloid"), ("Te", 4, 15, 52, "metalloid"),
            ("I", 4, 16, 53, "halogen"), ("Xe", 4, 17, 54, "noble_gas"),
            ("Cs", 5, 0, 55, "alkali_metal"), ("Ba", 5, 1, 56, "alkaline_earth_metal"),
            ("La", 7, 3, 57, "lanthanide"),
            ("Hf", 5, 3, 72, "transition_metal"), ("Ta", 5, 4, 73, "transition_metal"),
            ("W", 5, 5, 74, "transition_metal"), ("Re", 5, 6, 75, "transition_metal"),
            ("Os", 5, 7, 76, "transition_metal"), ("Ir", 5, 8, 77, "transition_metal"),
            ("Pt", 5, 9, 78, "transition_metal"), ("Au", 5, 10, 79, "transition_metal"),
            ("Hg", 5, 11, 80, "transition_metal"), ("Tl", 5, 12, 81, "post_transition_metal"),
            ("Pb", 5, 13, 82, "post_transition_metal"), ("Bi", 5, 14, 83, "post_transition_metal"),
            ("Po", 5, 15, 84, "metalloid"), ("At", 5, 16, 85, "halogen"), ("Rn", 5, 17, 86, "noble_gas"),
            ("Fr", 6, 0, 87, "alkali_metal"), ("Ra", 6, 1, 88, "alkaline_earth_metal"),
            ("Ac", 8, 3, 89, "actinide"),
            ("Rf", 6, 3, 104, "transition_metal"), ("Db", 6, 4, 105, "transition_metal"),
            ("Sg", 6, 5, 106, "transition_metal"), ("Bh", 6, 6, 107, "transition_metal"),
            ("Hs", 6, 7, 108, "transition_metal"), ("Mt", 6, 8, 109, "transition_metal"),
            ("Ds", 6, 9, 110, "transition_metal"), ("Rg", 6, 10, 111, "transition_metal"),
            ("Cn", 6, 11, 112, "transition_metal"), ("Nh", 6, 12, 113, "post_transition_metal"),
            ("Fl", 6, 13, 114, "post_transition_metal"), ("Mc", 6, 14, 115, "post_transition_metal"),
            ("Lv", 6, 15, 116, "post_transition_metal"), ("Ts", 6, 16, 117, "halogen"), ("Og", 6, 17, 118, "noble_gas"),
            ("Ce", 7, 4, 58, "lanthanide"), ("Pr", 7, 5, 59, "lanthanide"),
            ("Nd", 7, 6, 60, "lanthanide"), ("Pm", 7, 7, 61, "lanthanide"),
            ("Sm", 7, 8, 62, "lanthanide"), ("Eu", 7, 9, 63, "lanthanide"),
            ("Gd", 7, 10, 64, "lanthanide"), ("Tb", 7, 11, 65, "lanthanide"),
            ("Dy", 7, 12, 66, "lanthanide"), ("Ho", 7, 13, 67, "lanthanide"),
            ("Er", 7, 14, 68, "lanthanide"), ("Tm", 7, 15, 69, "lanthanide"),
            ("Yb", 7, 16, 70, "lanthanide"), ("Lu", 7, 17, 71, "lanthanide"),
            ("Th", 8, 4, 90, "actinide"), ("Pa", 8, 5, 91, "actinide"),
            ("U", 8, 6, 92, "actinide"), ("Np", 8, 7, 93, "actinide"),
            ("Pu", 8, 8, 94, "actinide"), ("Am", 8, 9, 95, "actinide"),
            ("Cm", 8, 10, 96, "actinide"), ("Bk", 8, 11, 97, "actinide"),
            ("Cf", 8, 12, 98, "actinide"), ("Es", 8, 13, 99, "actinide"),
            ("Fm", 8, 14, 100, "actinide"), ("Md", 8, 15, 101, "actinide"),
            ("No", 8, 16, 102, "actinide"), ("Lr", 8, 17, 103, "actinide"),
        ]


# =========================
# Right side: database selector + isotope selector + info tabs
# =========================
class Right(QVBoxLayout):

    def __init__(self, parent=None):
        super().__init__(None)
        self.model = parent
        self.w = self.model.win_w
        self.h = self.model.win_h

        self.svc = get_service()

        # ── container (so the loading overlay can cover everything) ──
        self.container = RightContainer()
        self.container_layout = QVBoxLayout(self.container)
        self.container_layout.setContentsMargins(0, 0, 0, 0)

        # ── database multi-select (the 3 ENDF-6 decay libraries; no TENDL) ──
        self.db_group = QGroupBox("Databases")
        self.db_group.setObjectName("db_group")
        db_lay = QHBoxLayout(self.db_group)
        db_lay.setContentsMargins(8, 4, 8, 4)
        self.db_checks: Dict[str, QCheckBox] = {}
        libs = self.svc.decay_libraries() or ["ENDFB-VIII.0", "JEFF", "JENDL"]
        for lib in libs:
            cb = QCheckBox(lib)
            cb.setObjectName("db_check")
            cb.setChecked(True)  # work with ALL selected bases by default
            cb.stateChanged.connect(self._on_databases_changed)
            db_lay.addWidget(cb)
            self.db_checks[lib] = cb

        # ── isotope selector ──
        self.isotope_cb = QComboBox()
        self.isotope_cb.setObjectName('select_cb')
        self.isotope_cb.setFixedHeight(self.h // 40)
        self.isotope_cb.addItem("Choose element on the left…")
        self.isotope_cb.currentIndexChanged.connect(self._on_isotope_changed)

        # ── info panels (split view: Decay Data above, NReaction Data below) ──
        self.info_main = QTextEdit()
        self.info_main.setObjectName('info_text_main')
        self.info_main.setReadOnly(True)
        self.info_main.setLineWrapMode(QTextEdit.NoWrap)

        self.info_chanels = QTextEdit()
        self.info_chanels.setObjectName('info_text_chanels')
        self.info_chanels.setReadOnly(True)
        self.info_chanels.setLineWrapMode(QTextEdit.NoWrap)

        decay_box = QGroupBox("Decay Data")
        decay_box.setObjectName("decay_data_box")
        decay_lay = QVBoxLayout(decay_box)
        decay_lay.setContentsMargins(4, 4, 4, 4)
        decay_lay.addWidget(self.info_main)

        nreac_box = QGroupBox("NReaction Data")
        nreac_box.setObjectName("nreac_data_box")
        nreac_lay = QVBoxLayout(nreac_box)
        nreac_lay.setContentsMargins(4, 4, 4, 4)
        nreac_lay.addWidget(self.info_chanels)

        self.info_splitter = QSplitter(Qt.Vertical)
        self.info_splitter.setObjectName("info_splitter")
        self.info_splitter.addWidget(decay_box)
        self.info_splitter.addWidget(nreac_box)
        self.info_splitter.setSizes([300, 300])
        self.info_splitter.setHandleWidth(6)

        # ── export buttons ──
        export_row = QHBoxLayout()
        self.btn_export_decay = QPushButton("Save DecayData")
        self.btn_export_decay.setObjectName("export_decay_btn")
        self.btn_export_decay.setEnabled(False)
        self.btn_export_decay.clicked.connect(self._export_decay)
        self.btn_export_neutron = QPushButton("Save NReactionData")
        self.btn_export_neutron.setObjectName("export_neutron_btn")
        self.btn_export_neutron.setEnabled(False)
        self.btn_export_neutron.clicked.connect(self._export_neutron)
        export_row.addWidget(self.btn_export_decay)
        export_row.addWidget(self.btn_export_neutron)

        self.container_layout.addWidget(self.db_group)
        self.container_layout.addSpacing(self.h // 90)
        self.container_layout.addWidget(self.isotope_cb)
        self.container_layout.addSpacing(8)
        self.container_layout.addWidget(self.info_splitter)
        self.container_layout.addSpacing(4)
        self.container_layout.addLayout(export_row)

        self.addWidget(self.container)
        self.setAlignment(Qt.AlignTop)

        # thread handles
        self._thread: Optional[QThread] = None
        self._worker: Optional[IsotopeLoadWorker] = None
        self._last_payload: Optional[dict] = None

        # current element context
        self._cur_Z: Optional[int] = None
        self._cur_X: Optional[str] = None

        self.info_main.setPlainText("Pick an element from the periodic table.")

    # ── selected database names ──
    def selected_libs(self) -> List[str]:
        return [lib for lib, cb in self.db_checks.items() if cb.isChecked()]

    # ── element chosen on the periodic table ──
    def on_element_selected(self, Z: int, X: str):
        self._cur_Z = int(Z)
        self._cur_X = str(X)
        self._populate_isotopes()

    def _on_databases_changed(self, *_):
        # Re-scan the isotope list for the current element under the new selection.
        if self._cur_Z is not None:
            self._populate_isotopes()

    def _populate_isotopes(self):
        libs = self.selected_libs()
        if not libs:
            self.isotope_cb.blockSignals(True)
            self.isotope_cb.clear()
            self.isotope_cb.addItem("Select at least one database")
            self.isotope_cb.blockSignals(False)
            self.info_main.setPlainText("Select at least one database above.")
            return

        # Union of isotopes available in ANY selected library.
        seen = {}
        for lib in libs:
            for (A, meta) in self.svc.scan_isotopes_for_element(lib, self._cur_Z, self._cur_X):
                seen[(A, meta)] = (A, meta)
        isotopes = sorted(seen.values(),
                          key=lambda t: (t[0], 0 if t[1] in (None, 0, "") else int(t[1])))

        self.isotope_cb.blockSignals(True)
        self.isotope_cb.clear()
        if not isotopes:
            self.isotope_cb.addItem("No decay files for this element")
            self.isotope_cb.blockSignals(False)
            self.info_main.setPlainText("No decay data found for this element in the selected databases.")
            self.info_chanels.setPlainText("")
            return
        self.isotope_cb.addItem("Choose isotope…")
        for (A, meta) in isotopes:
            self.isotope_cb.addItem(_format_isotope_label(self._cur_X, A, meta), (A, meta))
        self.isotope_cb.blockSignals(False)
        self.info_main.setPlainText(
            f"{self._cur_X} (Z={self._cur_Z}) — choose an isotope from the list above.")

    def _on_isotope_changed(self, idx: int):
        if idx <= 0 or self._cur_Z is None:
            return
        data = self.isotope_cb.itemData(idx)
        if not data:
            return
        A, meta = data
        self._start_load(self._cur_Z, self._cur_X, int(A), meta)

    # ── threaded load ──
    def _stop_thread(self):
        if self._thread is not None and self._thread.isRunning():
            self._thread.requestInterruption()
            self._thread.quit()
            self._thread.wait(3000)
        self._thread = None
        self._worker = None

    def _start_load(self, Z, X, A, meta):
        self._stop_thread()
        self.info_main.setPlainText("Reading nuclear data…")
        self.info_chanels.setPlainText("Reading nuclear data…")
        self.container.overlay.show()
        self.model.getLeftWidget().setEnabled(False)

        self._thread = QThread()
        self._worker = IsotopeLoadWorker(self.selected_libs(), Z, X, A, meta)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_loaded)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._on_thread_finished)
        self._worker.finished.connect(self._worker.deleteLater)
        self._worker.failed.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_thread_finished(self):
        self._thread = None
        self._worker = None

    def _restore_ui(self):
        self.container.overlay.hide()
        self.model.getLeftWidget().setEnabled(True)

    # ── rendering ──
    def _on_loaded(self, payload: dict):
        try:
            self._last_payload = payload
            header = _format_isotope_label(payload["X"], payload["A"], payload["meta"])
            self.info_main.setPlainText(self._render_decay(header, payload))
            self.info_chanels.setPlainText(self._render_neutron(header, payload))
            self.btn_export_decay.setEnabled(True)
            self.btn_export_neutron.setEnabled(True)
        except Exception as e:
            tb = traceback.format_exc()
            self.info_main.setPlainText(f"Render error: {type(e).__name__}: {e}\n\n{tb}")
            self.info_chanels.setPlainText("—")
        finally:
            self._restore_ui()

    def _on_failed(self, msg: str):
        self._restore_ui()
        self.info_main.setPlainText("Error reading isotope:\n\n" + msg)
        self.info_chanels.setPlainText("—")

    def _render_decay(self, header: str, payload: dict) -> str:
        lines = [header, f"Z = {payload['Z']}    A = {payload['A']}", ""]
        decay_by_lib = payload["decay_by_lib"]
        for lib in payload["libs"]:
            d = decay_by_lib.get(lib)
            lines.append(f"── {lib} ──")
            if not d:
                lines.append("   (no decay data in this database)")
                lines.append("")
                continue
            if d["is_stable"]:
                lines.append("   Stable nuclide")
            else:
                lines.append(f"   Half-life      : {_fmt_hl(d['half_life_s'])}")
                lines.append(f"   Decay constant : {d['decay_constant']:.4e} 1/s")
            mass = d.get("mass_amu")
            lines.append(f"   Atomic mass    : {mass:.4f} amu" if mass else "   Atomic mass    : —")
            chans = d.get("channels") or []
            if chans:
                lines.append("   Decay channels :")
                for ch in chans:
                    lines.append(
                        f"      • {ch['symbol']:<8}  BR = {ch['branch'] * 100:7.3f} %"
                        f"   Q = {_fmt_mev(ch['Q_MeV'])}")
            elif not d["is_stable"]:
                lines.append("   Decay channels : (none listed)")
            lines.append("")
        return "\n".join(lines)

    def _render_neutron(self, header: str, payload: dict) -> str:
        lines = [header, f"Z = {payload['Z']}    A = {payload['A']}", ""]

        # (A) per selected library — channels, σ(E), Q, isomeric branch, MACS@30 keV
        lines.append("═══ Neutron reactions (selected databases) ═══")
        neutron_by_lib = payload["neutron_by_lib"]
        for lib in payload["libs"]:
            chans = neutron_by_lib.get(lib) or []
            lines.append("")
            lines.append(f"── {lib} ──")
            if not chans:
                lines.append("   (no neutron evaluation in this database)")
                continue
            for ch in chans:
                lines.append(f"   {ch['name']:<9} MT{ch['mt']:<3}   Q = {_fmt_mev(ch['Q_MeV'])}")
                lines.append(
                    f"      σ:  thermal {_fmt_sigma(ch['sigma_thermal_barn'])}"
                    f" | 30 keV {_fmt_sigma(ch['sigma_30keV_barn'])}"
                    f" | 1 MeV {_fmt_sigma(ch['sigma_1MeV_barn'])}   [{ch['n_points']} pts]")
                lines.append(
                    f"      MACS(30 keV) = {_fmt_mb(ch['macs30_mb'])}"
                    "   (Maxwellian average of this database's spectrum)")
                br = ch.get("branches") or []
                if br and not (len(br) == 1 and br[0]["LFS"] == 0):
                    parts = ", ".join(f"LFS{b['LFS']} → {b['yield'] * 100:.1f} %" for b in br)
                    lines.append(f"      isomeric branch: {parts}")

        # (B) ALL MACS sources from xsdir/MACS — independent of selected libraries
        lines.append("")
        lines.append("═══ MACS sources — xsdir/MACS (all databases) ═══")
        macs = payload["macs"]

        raw = macs.get("rawmacs") or {}
        lines.append("")
        lines.append(" rawmacs.txt:")
        if raw:
            for lib, by_mt in raw.items():
                lines.append(f"   {lib}:")
                for mt, info in sorted(by_mt.items()):
                    pts = info.get("points") or []
                    lines.append(
                        f"      {mt_symbol(mt)}  MACS@30 keV = {_fmt_mb(info['macs30_mb'])}"
                        f"   [{len(pts)} kT pts]")
                    for kT, sig in pts:
                        lines.append(f"         kT={kT:8.3f} keV   σ={sig:10.5g} mb")
        else:
            lines.append("   (this isotope is not in rawmacs.txt)")

        e71 = macs.get("endfb71")
        e10 = macs.get("eaf2010")
        lines.append("")
        lines.append(" Multi-kT tables (n,γ):")
        for label, data in [("ENDF/B-VII.1", e71), ("EAF-2010", e10)]:
            if data:
                pts = data.get("points") or []
                lines.append(
                    f"   {label} : MACS@30 keV = {_fmt_mb(data['macs30_mb'])}"
                    f"   [{len(pts)} kT pts]")
                for kT, sig in pts:
                    lines.append(f"      kT={kT:8.3f} keV   σ={sig:10.5g} mb")
            else:
                lines.append(f"   {label} : —")

        talys = macs.get("talys")
        lines.append("")
        lines.append(" TALYS (recommended):")
        if talys:
            spectra = talys.get("spectra") or {}
            lines.append(f"   MACS@30 keV (n,γ) = {_fmt_mb(talys.get('macs30_mb'))}")
            if spectra:
                for mt, spec_info in sorted(spectra.items()):
                    name = spec_info.get("name", f"MT={mt}")
                    pts = spec_info.get("spectrum") or []
                    lines.append(f"   {name} — {len(pts)} points:")
                    for pt in pts:
                        E_MeV = pt.E_eV / 1e6
                        sigma_mb = pt.sigma_barn * 1e3
                        lines.append(f"      {E_MeV:10.5g} MeV   {sigma_mb:10.5g} mb")
            else:
                lines.append("   spectra: —")
        else:
            lines.append("   (no TALYS data for this isotope)")

        return "\n".join(lines)

    def _export_decay(self):
        if not self._last_payload:
            return
        p = self._last_payload
        default = f"decay_{p['Z']}-{p['X']}-{p['A']}.txt"
        path, _ = QFileDialog.getSaveFileName(
            self.container, "Save Decay Data", default, "Text files (*.txt);;All files (*)")
        if not path:
            return
        try:
            header = _format_isotope_label(p["X"], p["A"], p["meta"])
            with open(path, "w", encoding="utf-8") as f:
                f.write(self._render_decay(header, p))
        except Exception as e:
            QMessageBox.critical(self.container, "Export error", str(e))

    def _export_neutron(self):
        if not self._last_payload:
            return
        p = self._last_payload
        default = f"neutron_{p['Z']}-{p['X']}-{p['A']}.txt"
        path, _ = QFileDialog.getSaveFileName(
            self.container, "Save NReaction Data", default, "Text files (*.txt);;All files (*)")
        if not path:
            return
        try:
            header = _format_isotope_label(p["X"], p["A"], p["meta"])
            with open(path, "w", encoding="utf-8") as f:
                f.write(self._render_neutron(header, p))
        except Exception as e:
            QMessageBox.critical(self.container, "Export error", str(e))

    def stop_loading(self):
        self._stop_thread()
