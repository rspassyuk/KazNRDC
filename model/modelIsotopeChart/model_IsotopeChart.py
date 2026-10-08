# Qt________________________________________
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import math

# nuclear core — all data + burnup logic comes from the shared service.
from core.entities import ChemTable
from core.service import get_service, decay_symbol, mt_symbol

from GUIcomp import *  # noqa: F401,F403  (brings Worker / run_with_progress too)


CELL = 22  # scene units per chart cell


# =====================================================================
# Single source of truth for how each reaction/decay channel is drawn:
# colour + line style. Shared by BOTH the chart arrows (Left.draw_arrow) and
# the legend previews beside the reaction checkboxes (ArrowLegend), so the two
# can never drift apart. Channels are told apart by colour AND dash pattern
# (not colour alone) so they stay distinguishable for colour-blind readers.
# =====================================================================
REACTION_STYLES = {
    # decay modes
    "beta-":    ("#e5534b", Qt.SolidLine),
    "ec/beta+": ("#4dabf7", Qt.DashLine),
    "it":       ("#c0c0c0", Qt.DotLine),
    "alpha":    ("#7ed957", Qt.DashDotLine),
    "n":        ("#f4a261", Qt.DashDotDotLine),
    "p":        ("#ffd166", Qt.DotLine),
    "sf":       ("#9aa0a6", Qt.DashLine),
    "unknown":  ("#5a5f66", Qt.SolidLine),
    # neutron reactions
    "(n,γ)":    ("#1a6fa8", Qt.SolidLine),
    "(n,α)":    ("#2e8b57", Qt.DashDotLine),
    "(n,p)":    ("#e0883c", Qt.DashLine),
    "(n,2n)":   ("#9b59b6", Qt.DashDotDotLine),
    "(n,f)":    ("#c0392b", Qt.DotLine),
}

# Thinner arrows (Task 4). Chart width is in scene units; the legend width is
# in device pixels (a few px so the dash pattern reads in the tiny preview).
ARROW_WIDTH = 0.8
ARROW_HEAD = 6.0
ARROW_LEGEND_WIDTH = 2.0


def reaction_style(name):
    """(colour_hex, Qt.PenStyle) for a channel — the shared drawing recipe."""
    return REACTION_STYLES.get(name, ("#5a5f66", Qt.SolidLine))


class ArrowLegend(QWidget):
    """Tiny fixed-size preview of the chart arrow for one reaction channel,
    painted from the same colour + line style the chart itself uses."""

    def __init__(self, name, parent=None):
        super().__init__(parent)
        self._name = name
        self.setFixedSize(34, 14)
        self.setToolTip(f"Isotope-chart arrow for {name}")

    def paintEvent(self, event):
        color_str, style = reaction_style(self._name)
        color = QColor(color_str)
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing, True)
        y = self.height() / 2.0
        x1, x2 = 2.0, self.width() - 7.0
        pen = QPen(color, ARROW_LEGEND_WIDTH)
        pen.setStyle(style)
        p.setPen(pen)
        p.drawLine(QPointF(x1, y), QPointF(x2, y))
        # filled arrowhead (solid, like the chart)
        head = QPolygonF([
            QPointF(self.width() - 2.0, y),
            QPointF(x2, y - 3.5),
            QPointF(x2, y + 3.5),
        ])
        p.setPen(QPen(color, 1.0))
        p.setBrush(QBrush(color))
        p.drawPolygon(head)
        p.end()


class model_IsotopeChart(QObject):

    def __init__(self, parent=None):
        super().__init__()
        self.main_parent = parent
        self.win_w = parent.window_width
        self.win_h = parent.window_height

        self.__left = Left()
        self.__right = Right(self)

        self.__right.set_left_widget(self.__left)
        self.__left.set_reaction_panel(self.__right)

        # Two-way sync: chart click <-> dropdowns.
        self.__left.selectionChanged.connect(self.__right.on_chart_selection)
        self.__left.colorProgress.connect(self.__right.set_status)

        # Colour the full chart in the background (cached after first run).
        self.__left.start_color_loading()
        self.__left.update_arrows()

    def getLeftWidget(self):
        return self.__left

    def getRightLayout(self):
        return self.__right

    def getMainParent(self):
        return self.main_parent


# =====================================================================
# Single graphics item painting the WHOLE chart of nuclides (fast).
# =====================================================================
class ChartItem(QGraphicsItem):
    def __init__(self, chart):
        super().__init__()
        self.chart = chart

    def boundingRect(self):
        return QRectF(0, 0, self.chart.cols * CELL, self.chart.rows * CELL)

    def paint(self, painter, option, widget=None):
        self.chart._paint_chart(painter, option.exposedRect)


# =====================================================================
# Left — the full chart of nuclides (N along x, Z along y)
# =====================================================================
class Left(QGraphicsView):

    selectionChanged = pyqtSignal(object, object)  # (initial_label|None, last_label|None)
    colorProgress = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setObjectName("foo")
        self.svc = get_service()

        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setMinimumSize(10, 10)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMouseTracking(True)
        self.setRenderHints(self.renderHints() | QPainter.Antialiasing | QPainter.TextAntialiasing)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        self.setDragMode(QGraphicsView.ScrollHandDrag)

        self._zoom = 0
        self._zoom_step = 1.2
        self._cycle_izas: set = set()      # (Z,A) tuples in detected cycles
        self.setFocusPolicy(Qt.StrongFocus)

        # all nuclides (one cell per Z,A)
        self.cells = self.svc.chart_nuclides()
        self.cell_by_za = {(c["Z"], c["A"]): c for c in self.cells}
        self.modes = {}                      # (Z,A) -> mode, filled in background
        if self.cells:
            self.zmin = min(c["Z"] for c in self.cells)
            self.zmax = max(c["Z"] for c in self.cells)
            self.nmin = min(c["N"] for c in self.cells)
            self.nmax = max(c["N"] for c in self.cells)
        else:
            self.zmin, self.zmax, self.nmin, self.nmax = 0, 1, 0, 1
        self.rows = self.zmax - self.zmin + 1
        self.cols = self.nmax - self.nmin + 1

        # range selection
        self.sel_initial = None              # (Z, A)
        self.sel_last = None
        self._click_stage = 0                # 0: next click = initial, 1: = last, 2: restart
        self._sel_cells: set = set()         # cells enclosed by the selection shape (Task 5)

        self.reaction_panel = None
        self._arrow_items = []
        self._color_worker = None
        self._fitted = False

        # Arrow colours + line styles live in the module-level REACTION_STYLES
        # (single source of truth, shared with the legend previews).

        self.chart_item = ChartItem(self)
        self.scene.addItem(self.chart_item)
        self.scene.setSceneRect(self.chart_item.boundingRect())
        # NB: do NOT self.show() here. This widget is built during the launcher's
        # model-loading phase (behind the splash); showing it now pops it up as a
        # stray top-level window. It is parented into the main layout and shown
        # via the page mechanism once the main window is up.

    # -----------------------------------------------------------------
    # coordinate mapping
    # -----------------------------------------------------------------
    def cell_rect(self, Z, A) -> QRectF:
        N = A - Z
        x = (N - self.nmin) * CELL
        y = (self.zmax - Z) * CELL
        return QRectF(x, y, CELL, CELL)

    def _center(self, za) -> QPointF:
        return QPointF(self.cell_rect(*za).center())

    def za_at_scene(self, pt: QPointF):
        col = int(pt.x() // CELL)
        row = int(pt.y() // CELL)
        if col < 0 or row < 0 or col >= self.cols or row >= self.rows:
            return None
        N = col + self.nmin
        Z = self.zmax - row
        za = (Z, Z + N)
        return za if za in self.cell_by_za else None

    # -----------------------------------------------------------------
    # colours
    # -----------------------------------------------------------------
    def _mode_color(self, mode) -> QColor:
        if mode == "__loading__":
            return QColor("#454a51")
        if mode is None:                        # stable
            return QColor("#1b1b1b")
        m = str(mode).lower()
        if "," in m:
            return QColor("#8a2be2")
        return {
            "beta-": QColor("#e5534b"),
            "beta+": QColor("#4dabf7"), "ec/beta+": QColor("#4dabf7"),
            "ec": QColor("#7fd1a6"),
            "alpha": QColor("#7ed957"),
            "it": QColor("#c0c0c0"),
            "p": QColor("#ffd166"),
            "n": QColor("#f4a261"),
            "sf": QColor("#9aa0a6"),
        }.get(m, QColor("#5a5f66"))

    @staticmethod
    def _text_color(bg: QColor) -> QColor:
        return QColor("#101010") if bg.lightnessF() > 0.6 else QColor("#f0f0f0")

    # -----------------------------------------------------------------
    # painting
    # -----------------------------------------------------------------
    def _paint_chart(self, painter, exposed: QRectF):
        scale = painter.worldTransform().m11()
        show_labels = scale * CELL >= 26
        if show_labels:
            f = QFont("Arial")
            f.setPixelSize(int(CELL * 0.40))
            painter.setFont(f)

        for c in self.cells:
            r = self.cell_rect(c["Z"], c["A"])
            if not r.intersects(exposed):
                continue
            color = self._mode_color(self.modes.get((c["Z"], c["A"]), "__loading__"))
            painter.fillRect(r, color)
            if c["has_meta"]:
                painter.setPen(QPen(QColor("#ffd700"), 0))
                painter.setBrush(Qt.NoBrush)
                painter.drawRect(r.adjusted(0.5, 0.5, -1, -1))
            if show_labels:
                painter.setPen(self._text_color(color))
                painter.drawText(r, Qt.AlignCenter, f"{c['A']}{c['X']}")

        self._paint_selection(painter)

    def _paint_selection(self, painter):
        # Selection region: an arbitrary (non-rectangular) shape that follows
        # the active reaction channels (Task 5). _sel_cells is the set of cells
        # reachable from the selected range via those channels (computed in
        # core). Fill them and stroke only the *outer* cell borders, so the
        # boundary hugs the reachable set rather than a min/max rectangle.
        if self._sel_cells:
            fill = QColor(255, 205, 0, 45)
            for za in self._sel_cells:
                if za in self.cell_by_za:
                    painter.fillRect(self.cell_rect(*za), fill)
            painter.setBrush(Qt.NoBrush)
            painter.setPen(QPen(QColor("#ffcc00"), 2))
            for seg in self._selection_boundary_segments():
                painter.drawLine(seg)

        # Excluded isotopes (Task 6): translucent red fill + red border + a
        # diagonal cross, so excluded cells are unmistakable. They are also
        # dropped from the inputs the core build consumes.
        excluded = self.svc.get_excluded()
        if excluded:
            painter.setBrush(QColor(255, 59, 48, 70))
            painter.setPen(QPen(QColor("#ff3b30"), 2))
            for za in excluded:
                if za in self.cell_by_za:
                    r = self.cell_rect(*za)
                    painter.drawRect(r.adjusted(1, 1, -1, -1))
                    painter.drawLine(r.topLeft(), r.bottomRight())
                    painter.drawLine(r.topRight(), r.bottomLeft())

        # Cycle highlights — yellow border per cell
        if self._cycle_izas:
            painter.setBrush(Qt.NoBrush)
            painter.setPen(QPen(QColor("#ffff00"), 3))
            for za in self._cycle_izas:
                if za in self.cell_by_za:
                    painter.drawRect(self.cell_rect(*za).adjusted(2, 2, -2, -2))

        painter.setBrush(Qt.NoBrush)
        if self.sel_initial:
            painter.setPen(QPen(QColor("#19e0ff"), 3))
            painter.drawRect(self.cell_rect(*self.sel_initial).adjusted(1.5, 1.5, -1.5, -1.5))
        if self.sel_last:
            painter.setPen(QPen(QColor("#ff4081"), 3))
            painter.drawRect(self.cell_rect(*self.sel_last).adjusted(1.5, 1.5, -1.5, -1.5))

    # -----------------------------------------------------------------
    # selection shape (Task 5)
    # -----------------------------------------------------------------
    def recompute_selection_shape(self):
        """Recompute the cells enclosed by the selection outline.

        Seed = chart cells inside the selected rectangular range; the reachable
        closure via the *active* reaction channels is computed in core
        (svc.reachable_cells). The result (an arbitrary, non-rectangular set)
        drives _paint_selection. Called whenever the range or the active
        channels change — never per-paint.
        """
        if not (self.sel_initial and self.sel_last):
            if self._sel_cells:
                self._sel_cells = set()
                self.chart_item.update()
                self.update_arrows()
            return
        (z1, a1), (z2, a2) = self.sel_initial, self.sel_last
        n1, n2 = a1 - z1, a2 - z2
        zlo, zhi = sorted((z1, z2))
        nlo, nhi = sorted((n1, n2))
        seed = [(c["Z"], c["A"]) for c in self.cells
                if zlo <= c["Z"] <= zhi and nlo <= c["N"] <= nhi]
        channels = []
        if self.reaction_panel and hasattr(self.reaction_panel, "get_enabled_reactions"):
            channels = self.reaction_panel.get_enabled_reactions()
        self._sel_cells = self.svc.reachable_cells(seed, channels)
        self.chart_item.update()
        self.update_arrows()

    def _selection_boundary_segments(self):
        """Outer borders of the reachable selection set: cell edges that are
        NOT shared with another selected cell (the union-of-cell-borders
        outline). Works for concave shapes and holes."""
        segs = []
        cells = self._sel_cells
        for (Z, A) in cells:
            if (Z, A) not in self.cell_by_za:
                continue
            r = self.cell_rect(Z, A)
            x, y, w, h = r.x(), r.y(), r.width(), r.height()
            if (Z, A - 1) not in cells:        # left edge  (lower N)
                segs.append(QLineF(x, y, x, y + h))
            if (Z, A + 1) not in cells:        # right edge (higher N)
                segs.append(QLineF(x + w, y, x + w, y + h))
            if (Z + 1, A + 1) not in cells:    # top edge   (one row up, higher Z)
                segs.append(QLineF(x, y, x + w, y))
            if (Z - 1, A - 1) not in cells:    # bottom edge(one row down, lower Z)
                segs.append(QLineF(x, y + h, x + w, y + h))
        return segs

    # -----------------------------------------------------------------
    # selection
    # -----------------------------------------------------------------
    def _label(self, za):
        if za is None:
            return None
        c = self.cell_by_za.get(za)
        return f"{c['X']}-{c['A']}" if c else None

    def za_from_label(self, label):
        spec = self.svc.parse_label(label or "")
        if not spec:
            return None
        za = (spec[0], spec[2])
        return za if za in self.cell_by_za else None

    def _emit_selection(self):
        self.selectionChanged.emit(self._label(self.sel_initial), self._label(self.sel_last))

    def set_cycle_highlights(self, izas: set):
        self._cycle_izas = set(izas)
        self.chart_item.update()

    def set_selection(self, initial_za, last_za, emit=True):
        self.sel_initial = initial_za if initial_za in self.cell_by_za else None
        self.sel_last = last_za if last_za in self.cell_by_za else None
        self._click_stage = 2 if (self.sel_initial and self.sel_last) else (1 if self.sel_initial else 0)
        self.recompute_selection_shape()
        self.chart_item.update()
        if emit:
            self._emit_selection()

    def _handle_click(self, za):
        if self._click_stage == 0:
            self.sel_initial, self.sel_last, self._click_stage = za, None, 1
        elif self._click_stage == 1:
            self.sel_last, self._click_stage = za, 2
        else:                                  # third click → new selection
            self.sel_initial, self.sel_last, self._click_stage = za, None, 1
        self.recompute_selection_shape()
        self.chart_item.update()
        self._emit_selection()

    # -----------------------------------------------------------------
    # events
    # -----------------------------------------------------------------
    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            za = self.za_at_scene(self.mapToScene(event.pos()))
            if za is not None:
                self._handle_click(za)
                event.accept()
                return
        elif event.button() == Qt.RightButton:
            za = self.za_at_scene(self.mapToScene(event.pos()))
            if za is not None:
                self.toggle_excluded(za)   # right-click excludes / un-excludes (Task 6)
                event.accept()
                return
        super().mousePressEvent(event)

    def toggle_excluded(self, za):
        """Right-click toggles whether an isotope is excluded from calculations
        (Task 6): updates the service (so the build drops it), repaints the
        chart overlay, and refreshes the reaction panel's excluded list."""
        if self.svc.is_excluded(za):
            self.svc.remove_excluded(za)
        else:
            self.svc.add_excluded(za)
        self.chart_item.update()
        if self.reaction_panel and hasattr(self.reaction_panel, "refresh_excluded_list"):
            self.reaction_panel.refresh_excluded_list()

    def sync_excluded(self):
        """Repaint after the excluded set was changed elsewhere (e.g. via the
        reaction panel's list widget)."""
        self.chart_item.update()

    def mouseMoveEvent(self, event):
        za = self.za_at_scene(self.mapToScene(event.pos()))
        if za is not None:
            c = self.cell_by_za[za]
            mode = self.modes.get(za, "__loading__")
            if mode == "__loading__":
                mtxt = "…"
            elif mode is None:
                mtxt = "stable"
            elif mode == "unknown":
                mtxt = "—"
            else:
                mtxt = decay_symbol(mode)
            meta = f"  (+meta {sorted(c['metas'])})" if c["has_meta"] else ""
            QToolTip.showText(
                event.globalPos(),
                f"{c['X']}-{c['A']}   Z={c['Z']}  N={c['N']}{meta}\nprimary decay: {mtxt}")
        else:
            QToolTip.hideText()
        super().mouseMoveEvent(event)

    def wheelEvent(self, event):
        if event.angleDelta().y() > 0:
            factor, self._zoom = self._zoom_step, self._zoom + 1
        else:
            factor, self._zoom = 1 / self._zoom_step, self._zoom - 1
        if self._zoom < -6:
            self._zoom = -6
            return
        if self._zoom > 25:
            self._zoom = 25
            return
        self.scale(factor, factor)

    def showEvent(self, event):
        super().showEvent(event)
        if not self._fitted and self.cells:
            self.fitInView(self.scene.sceneRect(), Qt.KeepAspectRatio)
            self._fitted = True

    # -----------------------------------------------------------------
    # background colouring
    # -----------------------------------------------------------------
    def start_color_loading(self):
        worker = Worker(self.svc.chart_decay_modes)
        worker.setAutoDelete(False)          # we hold the reference (self._color_worker)
        worker.signals.progress.connect(lambda f, m: self.colorProgress.emit(m))
        worker.signals.finished.connect(self._on_modes_ready)
        worker.signals.error.connect(lambda tb: print("[chart] colour load failed:\n", tb))
        self._color_worker = worker
        QThreadPool.globalInstance().start(worker)

    def _on_modes_ready(self, modes):
        self.modes = modes or {}
        self.colorProgress.emit("")
        self.chart_item.update()

    # -----------------------------------------------------------------
    # arrows (built universe only — they need wired products)
    # -----------------------------------------------------------------
    def set_reaction_panel(self, reaction_panel):
        self.reaction_panel = reaction_panel

    def decay_ui_names(self, mode):
        """Map a raw decay-mode string from the reader to the display checkbox key(s)."""
        if not mode:
            return []
        m = str(mode).strip().lower()
        # composite modes (e.g. "beta-, alpha") → treat as unknown
        if "," in m:
            return ["unknown"]
        return {
            "beta-": ["beta-"],
            "beta+": ["ec/beta+"], "ec": ["ec/beta+"], "ec/beta+": ["ec/beta+"],
            "it": ["it"],
            "alpha": ["alpha"],
            "n": ["n"],
            "p": ["p"],
            "sf": ["sf"],
            "unknown": ["unknown"],
        }.get(m, ["unknown"])

    def reaction_id_to_name(self, mt):
        return {102: "(n,γ)", 107: "(n,α)", 103: "(n,p)",
                16: "(n,2n)", 18: "(n,f)"}.get(int(mt), None)

    def _clear_arrows(self):
        for it in self._arrow_items:
            self.scene.removeItem(it)
        self._arrow_items = []

    def update_arrows(self):
        self._clear_arrows()
        if not self.reaction_panel or not hasattr(self.reaction_panel, "arrows_cb"):
            return
        if not self.reaction_panel.arrows_cb.isChecked():
            return
        allowed = set(self.reaction_panel.get_enabled_reactions())
        sel = self._sel_cells  # non-empty → restrict to selected isotopes only
        for iso in self.svc.isotopes:
            start = (iso.Z, iso.A)
            if sel and start not in sel:
                continue   # source outside selection: skip entirely
            for d in iso.getListOfDecays():
                hit = next((n for n in self.decay_ui_names(d.mode) if n in allowed), None)
                if hit and d.product is not None:
                    end = (d.product.Z, d.product.A)
                    outside = bool(sel) and end not in sel
                    self.draw_arrow(start, end, hit, outside=outside)
            for r in iso.getListOfReactions():
                rn = self.reaction_id_to_name(r.getId())
                if rn and rn in allowed and r.product is not None:
                    end = (r.product.Z, r.product.A)
                    outside = bool(sel) and end not in sel
                    self.draw_arrow(start, end, rn, outside=outside)

    def draw_arrow(self, start_za, end_za, name, outside: bool = False):
        if start_za not in self.cell_by_za or end_za not in self.cell_by_za:
            return
        p1, p2 = self._center(start_za), self._center(end_za)
        dx, dy = p2.x() - p1.x(), p2.y() - p1.y()
        dist = math.hypot(dx, dy)
        if dist < 1e-6:
            return
        ux, uy = dx / dist, dy / dist

        if outside:
            # Destination is outside the selected set: draw as a faint dotted
            # gray arrow so the user can see "where it would go" without it
            # competing visually with the in-selection arrows.
            color = QColor(160, 160, 160, 100)
            pen = QPen(color, ARROW_WIDTH * 0.7)
            pen.setStyle(Qt.DotLine)
        else:
            color_str, style = reaction_style(name)
            color = QColor(color_str)
            pen = QPen(color, ARROW_WIDTH)
            pen.setStyle(style)

        margin = CELL * 0.5
        end = QPointF(p2.x() - ux * margin, p2.y() - uy * margin)
        line = QGraphicsLineItem(p1.x(), p1.y(), end.x(), end.y())
        line.setPen(pen)
        line.setZValue(5 if not outside else 4)
        self.scene.addItem(line)
        self._arrow_items.append(line)

        size = ARROW_HEAD * (0.6 if outside else 1.0)
        ang = math.atan2(uy, ux)
        left = QPointF(end.x() - size * math.cos(ang - math.radians(25)),
                       end.y() - size * math.sin(ang - math.radians(25)))
        right = QPointF(end.x() - size * math.cos(ang + math.radians(25)),
                        end.y() - size * math.sin(ang + math.radians(25)))
        head = QGraphicsPolygonItem(QPolygonF([end, left, right]))
        head.setBrush(QBrush(color))
        head.setPen(QPen(color, ARROW_WIDTH * (0.7 if outside else 1.0)))
        head.setZValue(6 if not outside else 4)
        self.scene.addItem(head)
        self._arrow_items.append(head)


# Reaction MT values for the reaction checkboxes
_REACTION_MTS = {
    "(n,γ)": 102, "(n,α)": 107, "(n,p)": 103, "(n,2n)": 16, "(n,f)": 18,
}


# =====================================================================
# Right — isotope range, reaction filter, calc type, build
# =====================================================================
class Right(QVBoxLayout):

    def __init__(self, parent=None):
        super().__init__(None)
        self.model = parent
        self.w = self.model.win_w
        self.h = self.model.win_h
        self.svc = get_service()
        self.left_widget = None

        # -- Range header --
        self.set_label = QLabel('Range (click two cells, use dropdowns, or type manually)')
        self.set_label.setObjectName('set_label')
        self.set_label.setWordWrap(True)

        # -- Autocomplete dropdowns (chart-click driven) --
        self.from_box = self._make_isotope_combo('from_box')
        self.to_box = self._make_isotope_combo('to_box')
        self._fill_isotope_combos()
        self.from_box.currentTextChanged.connect(self._on_combo_changed)
        self.to_box.currentTextChanged.connect(self._on_combo_changed)

        # -- Manual text input --
        self.manual_label = QLabel('Manual:')
        self.manual_label.setObjectName('manual_list_label')
        self.manual_edit = QLineEdit()
        self.manual_edit.setObjectName('manual_list_edit')
        self.manual_edit.setPlaceholderText('e.g.  Tl-206 Tl-207 Pb-206 Pb-207')

        self.status_label = QLabel('')
        self.status_label.setObjectName('chart_status_label')
        self.status_label.setWordWrap(True)

        # -- Arrows toggle --
        self.arrows_cb = QCheckBox('Show arrows')
        self.arrows_cb.setObjectName('arrows_cb')
        self.arrows_cb.setChecked(True)
        self.arrows_cb.stateChanged.connect(self._on_arrows_changed)

        # -- Reaction-config tabs (Task 3) --
        # These controls (decay types, neutron reactions, calculation type) used
        # to stack vertically as three group boxes and overflowed the dashboard.
        # They now live in a 3-tab QTabWidget so only one group shows at a time
        # and nothing clips or needs scrolling.
        self.config_tabs = QTabWidget()
        self.config_tabs.setObjectName('reaction_tabs')
        self.config_tabs.setElideMode(Qt.ElideNone)
        self.config_tabs.tabBar().setExpanding(False)
        self.config_tabs.tabBar().setUsesScrollButtons(True)

        # Tab 1 — Decay types
        decay_page = QWidget()
        decay_page.setObjectName('reaction_tab_page')
        decay_layout = QVBoxLayout(decay_page)
        decay_layout.setContentsMargins(6, 6, 6, 6)
        decay_layout.setSpacing(2)
        decay_layout.setAlignment(Qt.AlignTop)
        self._decay_checks: dict = {}
        for key in ["beta-", "ec/beta+", "it", "alpha", "n", "sf", "p", "unknown"]:
            row, cb = self._make_reaction_row(key, f"decay_{key.replace('/', '_')}_cb")
            decay_layout.addWidget(row)
            self._decay_checks[key] = cb
        self.config_tabs.addTab(decay_page, 'Decay')

        # Tab 2 — Neutron reactions
        reac_page = QWidget()
        reac_page.setObjectName('reaction_tab_page')
        reac_layout = QVBoxLayout(reac_page)
        reac_layout.setContentsMargins(6, 6, 6, 6)
        reac_layout.setSpacing(2)
        reac_layout.setAlignment(Qt.AlignTop)
        self._reac_checks: dict = {}
        for key in ["(n,γ)", "(n,α)", "(n,p)", "(n,2n)", "(n,f)"]:
            row, cb = self._make_reaction_row(
                key, f"reac_{key.replace(',', '').replace('(', '').replace(')', '')}_cb")
            reac_layout.addWidget(row)
            self._reac_checks[key] = cb
        self.config_tabs.addTab(reac_page, 'Neutron')

        # Tab 3 — Calculation type: Reactor physics / Astrophysics (§2.3)
        calc_page = QWidget()
        calc_page.setObjectName('reaction_tab_page')
        calc_layout = QVBoxLayout(calc_page)
        calc_layout.setContentsMargins(6, 6, 6, 6)
        calc_layout.setAlignment(Qt.AlignTop)
        self.calc_reactor_rb = QRadioButton('Reactor physics')
        self.calc_reactor_rb.setObjectName('calc_reactor_rb')
        self.calc_reactor_rb.setChecked(True)
        self.calc_astro_rb = QRadioButton('Astrophysics')
        self.calc_astro_rb.setObjectName('calc_astro_rb')
        calc_layout.addWidget(self.calc_reactor_rb)
        calc_layout.addWidget(self.calc_astro_rb)
        self.calc_reactor_rb.toggled.connect(self._on_calc_type_changed)

        # -- Reactor physics sub-panel --
        self.reactor_widget = QWidget()
        self.reactor_widget.setObjectName('reactor_widget')
        rv = QVBoxLayout(self.reactor_widget)
        rv.setContentsMargins(0, 0, 0, 0)

        rv.addWidget(QLabel('Database (multi-select):'))
        self.db_list = QListWidget()
        self.db_list.setObjectName('db_list')
        self.db_list.setSelectionMode(QAbstractItemView.NoSelection)
        _dec_libs = ["user"] + list(self.svc.decay_libraries() or ["ENDFB-VIII.0"])
        if "TALYS" not in _dec_libs:
            _dec_libs.append("TALYS")
        _dec_libs.append("All")
        for lib in _dec_libs:
            item = QListWidgetItem(lib)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if lib in ("user", "All") else Qt.Unchecked)
            self.db_list.addItem(item)
        self.db_list.setMaximumHeight(110)
        self.db_list.itemChanged.connect(self._on_db_list_changed)
        rv.addWidget(self.db_list)

        energy_row = QHBoxLayout()
        energy_row.addWidget(QLabel('Energy:'))
        self.energy_box = QComboBox()
        self.energy_box.setObjectName('energy_box')
        for label, ev in [("0.0253 eV (thermal)", 0.0253), ("1 eV", 1.0), ("1 keV", 1e3),
                          ("30 keV", 3e4), ("100 keV", 1e5), ("1 MeV", 1e6),
                          ("14.1 MeV (DT)", 1.41e7), ("Manual…", None)]:
            self.energy_box.addItem(label, ev)
        self.energy_box.setCurrentIndex(3)
        self.energy_box.currentIndexChanged.connect(self._on_energy_combo_changed)
        energy_row.addWidget(self.energy_box)
        rv.addLayout(energy_row)

        self.energy_manual_edit = QLineEdit()
        self.energy_manual_edit.setObjectName('energy_manual_edit')
        self.energy_manual_edit.setPlaceholderText('Energy in eV')
        self.energy_manual_edit.setVisible(False)
        rv.addWidget(self.energy_manual_edit)

        # -- Astrophysics sub-panel --
        self.astro_widget = QWidget()
        self.astro_widget.setObjectName('astro_widget')
        av = QVBoxLayout(self.astro_widget)
        av.setContentsMargins(0, 0, 0, 0)
        av.addWidget(QLabel('MACS source (multi-select):'))
        self.macs_source_list = QListWidget()
        self.macs_source_list.setObjectName('macs_source_list')
        self.macs_source_list.setSelectionMode(QAbstractItemView.NoSelection)
        for s in ["user", "rawmacs", "ENDFB71", "EAF2010", "TALYS_30keV", "All"]:
            item = QListWidgetItem(s)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if s in ("user", "All") else Qt.Unchecked)
            self.macs_source_list.addItem(item)
        self.macs_source_list.setMaximumHeight(110)
        av.addWidget(self.macs_source_list)

        astro_energy_row = QHBoxLayout()
        astro_energy_row.addWidget(QLabel('Energy:'))
        self.astro_energy_box = QComboBox()
        self.astro_energy_box.setObjectName('astro_energy_box')
        for label, ev in [("0.0253 eV (thermal)", 0.0253), ("1 eV", 1.0), ("1 keV", 1e3),
                          ("30 keV", 3e4), ("100 keV", 1e5), ("1 MeV", 1e6),
                          ("14.1 MeV (DT)", 1.41e7), ("Manual…", None)]:
            self.astro_energy_box.addItem(label, ev)
        self.astro_energy_box.setCurrentIndex(3)  # 30 keV default
        self.astro_energy_box.currentIndexChanged.connect(self._on_astro_energy_combo_changed)
        astro_energy_row.addWidget(self.astro_energy_box)
        av.addLayout(astro_energy_row)

        self.astro_energy_manual = QLineEdit()
        self.astro_energy_manual.setObjectName('astro_energy_manual')
        self.astro_energy_manual.setPlaceholderText('Energy in eV')
        self.astro_energy_manual.setVisible(False)
        av.addWidget(self.astro_energy_manual)

        self.astro_widget.setVisible(False)

        calc_layout.addWidget(self.reactor_widget)
        calc_layout.addWidget(self.astro_widget)
        self.config_tabs.addTab(calc_page, 'Calculation')

        # -- Build button --
        self.burnup_pb = QPushButton('Build Burnup Matrix')
        self.burnup_pb.setObjectName('burnup_pb')
        self.burnup_pb.clicked.connect(self._on_build_burnup)

        self.find_cycle_pb = QPushButton('Find Cycle Reactions')
        self.find_cycle_pb.setObjectName('find_cycle_pb')
        self.find_cycle_pb.clicked.connect(self._on_find_cycles)

        # -- Excluded isotopes (Task 6) --
        # Right-click a chart cell to exclude it from calculations; this list
        # mirrors the excluded set (sourced from core) and offers an un-exclude
        # path. Excluded isotopes are dropped from the build inputs in core.
        self.excluded_label = QLabel('Excluded from calculations (right-click a cell):')
        self.excluded_label.setObjectName('excluded_label')
        self.excluded_label.setWordWrap(True)
        self.excluded_list = QListWidget()
        self.excluded_list.setObjectName('excluded_list')
        self.excluded_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.excluded_list.setMaximumHeight(90)
        self._excl_btn_row = QHBoxLayout()
        self.unexclude_pb = QPushButton('Un-exclude selected')
        self.unexclude_pb.setObjectName('unexclude_pb')
        self.unexclude_pb.clicked.connect(self._on_unexclude_selected)
        self.clear_excluded_pb = QPushButton('Clear')
        self.clear_excluded_pb.setObjectName('clear_excluded_pb')
        self.clear_excluded_pb.clicked.connect(self._on_clear_excluded)
        self._excl_btn_row.addWidget(self.unexclude_pb)
        self._excl_btn_row.addWidget(self.clear_excluded_pb)

        # -- Assemble inside a scroll area so the panel works on small screens --
        inner_w = QWidget()
        inner_v = QVBoxLayout(inner_w)
        inner_v.setContentsMargins(4, 4, 4, 4)
        inner_v.setSpacing(4)
        inner_v.addWidget(self.set_label)
        inner_v.addSpacing(4)
        inner_v.addWidget(QLabel('From isotope:'))
        inner_v.addWidget(self.from_box)
        inner_v.addWidget(QLabel('To isotope:'))
        inner_v.addWidget(self.to_box)
        inner_v.addSpacing(4)
        inner_v.addWidget(self.manual_label)
        inner_v.addWidget(self.manual_edit)
        inner_v.addWidget(self.status_label)
        inner_v.addSpacing(6)
        inner_v.addWidget(self.arrows_cb)
        inner_v.addWidget(self.config_tabs)
        inner_v.addSpacing(6)
        inner_v.addWidget(self.excluded_label)
        inner_v.addWidget(self.excluded_list)
        inner_v.addLayout(self._excl_btn_row)
        inner_v.addSpacing(8)
        inner_v.addWidget(self.burnup_pb, alignment=Qt.AlignCenter)
        inner_v.addWidget(self.find_cycle_pb, alignment=Qt.AlignCenter)
        inner_v.addStretch(1)

        right_scroll = QScrollArea()
        right_scroll.setObjectName('right_panel_scroll')
        right_scroll.setWidgetResizable(True)
        right_scroll.setWidget(inner_w)
        right_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.addWidget(right_scroll)
        self.refresh_excluded_list()

    # -- reaction rows (legend + checkbox) ---------------------------------
    def _make_reaction_row(self, key, object_name):
        """A reaction checkbox preceded by an arrow-style legend (Task 4): the
        legend previews the exact colour + line style the chart draws."""
        row = QWidget()
        row.setObjectName('reaction_tab_page')   # inherit themed page bg/text
        rl = QHBoxLayout(row)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(8)
        legend = ArrowLegend(key)
        cb = QCheckBox(key)
        cb.setObjectName(object_name)
        cb.stateChanged.connect(self._on_arrows_changed)
        rl.addWidget(legend)
        rl.addWidget(cb)
        rl.addStretch(1)
        return row, cb

    # -- isotope combos ----------------------------------------------------
    def _make_isotope_combo(self, name):
        box = QComboBox()
        box.setObjectName(name)
        box.setEditable(True)
        box.setInsertPolicy(QComboBox.NoInsert)
        box.setMaxVisibleItems(20)
        comp = box.completer()
        comp.setCompletionMode(QCompleter.PopupCompletion)
        comp.setFilterMode(Qt.MatchContains)
        comp.setCaseSensitivity(Qt.CaseInsensitive)
        return box

    def _fill_isotope_combos(self):
        for box in (self.from_box, self.to_box):
            box.blockSignals(True)
            box.clear()
            box.addItem("")
            for c in self.svc.chart_nuclides():
                box.addItem(f"{c['X']}-{c['A']}", (c["Z"], c["A"]))
            box.setCurrentIndex(0)
            box.blockSignals(False)

    def _za_from_combo(self, box):
        spec = self.svc.parse_label(box.currentText())
        if not spec:
            return None
        za = (spec[0], spec[2])
        return za if (self.left_widget and za in self.left_widget.cell_by_za) else None

    def _on_combo_changed(self, *_):
        if not self.left_widget:
            return
        iz = self._za_from_combo(self.from_box)
        lz = self._za_from_combo(self.to_box)
        self.left_widget.set_selection(iz, lz)

    def on_chart_selection(self, initial_label, last_label):
        for box, lab in ((self.from_box, initial_label), (self.to_box, last_label)):
            box.blockSignals(True)
            box.setCurrentText(lab or "")
            box.blockSignals(False)

    def set_status(self, text):
        self.status_label.setText(text or "")

    # -- calc type switch --------------------------------------------------
    def _on_calc_type_changed(self, reactor_checked):
        self.reactor_widget.setVisible(reactor_checked)
        self.astro_widget.setVisible(not reactor_checked)

    def _on_energy_combo_changed(self, idx):
        is_manual = self.energy_box.currentData() is None
        self.energy_manual_edit.setVisible(is_manual)

    def _on_astro_energy_combo_changed(self, idx):
        is_manual = self.astro_energy_box.currentData() is None
        self.astro_energy_manual.setVisible(is_manual)

    def _read_energy_eV(self) -> float:
        if self.calc_astro_rb.isChecked():
            data = self.astro_energy_box.currentData()
            if data is None:
                txt = self.astro_energy_manual.text().strip()
                try:
                    return float(txt)
                except ValueError:
                    return 30_000.0
            return float(data)
        data = self.energy_box.currentData()
        if data is None:
            txt = self.energy_manual_edit.text().strip()
            try:
                return float(txt)
            except ValueError:
                return 30_000.0
        return float(data)

    def _get_mt_filter(self) -> list:
        mts = []
        for name, cb in self._reac_checks.items():
            if cb.isChecked():
                mts.append(_REACTION_MTS[name])
        return mts if mts else [102]

    _UI_TO_BURNUP_MODES = {
        "beta-":    ["beta-"],
        "ec/beta+": ["ec", "beta+"],
        "it":       ["it"],
        "alpha":    ["alpha"],
        "n":        ["n"],
        "sf":       ["sf"],
        "p":        ["p"],
        "unknown":  ["unknown"],
    }

    def _get_decay_mode_filter(self):
        """Return None if all decay types checked (no filter), else list of burnup mode strings."""
        checked = [k for k, cb in self._decay_checks.items() if cb.isChecked()]
        if set(checked) == set(self._decay_checks.keys()):
            return None
        modes = []
        for key in checked:
            modes.extend(self._UI_TO_BURNUP_MODES.get(key, []))
        return modes

    def _get_prefer_macs(self) -> bool:
        return self.calc_astro_rb.isChecked()

    # -- arrows / reactions ------------------------------------------------
    def _on_arrows_changed(self, *_):
        if self.left_widget:
            self.left_widget.update_arrows()
            # The selection shape follows the active channels, so refresh it too.
            self.left_widget.recompute_selection_shape()

    def get_enabled_reactions(self):
        enabled = []
        for key, cb in self._decay_checks.items():
            if cb.isChecked():
                enabled.append(key)
        for key, cb in self._reac_checks.items():
            if cb.isChecked():
                enabled.append(key)
        return enabled

    def set_left_widget(self, left_widget):
        self.left_widget = left_widget

    # -- excluded isotopes (Task 6) ----------------------------------------
    def _za_label(self, Z, A):
        c = self.left_widget.cell_by_za.get((Z, A)) if self.left_widget else None
        return f"{c['X']}-{c['A']}" if c else f"Z{Z}·A{A}"

    def refresh_excluded_list(self):
        """Mirror the service's excluded (Z, A) set into the list widget."""
        self.excluded_list.blockSignals(True)
        self.excluded_list.clear()
        for (Z, A) in sorted(self.svc.get_excluded()):
            item = QListWidgetItem(self._za_label(Z, A))
            item.setData(Qt.UserRole, (Z, A))
            self.excluded_list.addItem(item)
        self.excluded_list.blockSignals(False)

    def _on_unexclude_selected(self):
        for it in self.excluded_list.selectedItems():
            za = it.data(Qt.UserRole)
            if za is not None:
                self.svc.remove_excluded(tuple(za))
        self.refresh_excluded_list()
        if self.left_widget:
            self.left_widget.sync_excluded()

    def _on_clear_excluded(self):
        self.svc.clear_excluded()
        self.refresh_excluded_list()
        if self.left_widget:
            self.left_widget.sync_excluded()

    # -- build burnup matrix (threaded) ------------------------------------
    def _dialog_parent(self):
        return self.left_widget

    def _read_checked_libs(self, list_widget, all_libs):
        """Return checked items; if 'All' checked, return all_libs."""
        checked = []
        for i in range(list_widget.count()):
            item = list_widget.item(i)
            if item.checkState() == Qt.Checked:
                checked.append(item.text())
        if "All" in checked:
            return list(all_libs)
        return [c for c in checked if c != "All"] or list(all_libs)

    def _on_db_list_changed(self, _item=None):
        """Different libraries have different branching ratios — re-run the
        build automatically so the chart's decay-transition tree (arrows)
        reflects the newly checked decay base right away, without requiring
        a manual rebuild click. Only acts once a universe already exists;
        toggling checkboxes before the first build does nothing extra."""
        if self.svc.isotopes:
            self._on_build_burnup()

    def _on_build_burnup(self):
        if not self.left_widget:
            return

        _all_dec = ["user", "endf6", "TALYS", "endf_computed"]
        dec_libs_raw = self._read_checked_libs(self.db_list, _all_dec)
        # A checked item may be a specific ENDF library name (e.g.
        # "ENDFB-VIII.0") rather than one of the generic priority tokens
        # above — the provider only understands "endf6". Point the active
        # ENDF decay directory at the checked library and translate its
        # name to "endf6" so decay lookups still reach it.
        endf_names = [lib for lib in dec_libs_raw if lib in self.svc.decay_libraries()]
        if endf_names:
            self.svc.switch_decay_library(endf_names[0])
        seen = set()
        dec_libs = []
        for lib in dec_libs_raw:
            token = "endf6" if lib in endf_names else lib
            if token not in seen:
                seen.add(token)
                dec_libs.append(token)
        self.svc.set_dec_priority(dec_libs)

        _all_macs = ["user", "rawmacs", "ENDFB71", "EAF2010", "TALYS_30keV"]
        macs_libs_raw = self._read_checked_libs(self.macs_source_list, _all_macs)
        # The checklist shows display names ("ENDFB71", "TALYS_30keV") but the
        # provider's priority tokens are lowercase ("endfb71", "talys").
        _macs_token_map = {"ENDFB71": "endfb71", "EAF2010": "eaf2010",
                           "TALYS_30keV": "talys"}
        seen = set()
        macs_libs = []
        for lib in macs_libs_raw:
            token = _macs_token_map.get(lib, lib)
            if token not in seen:
                seen.add(token)
                macs_libs.append(token)
        self.svc.set_macs_priority(macs_libs)

        manual_text = self.manual_edit.text().strip()
        E_eV = self._read_energy_eV()
        mt_filter = self._get_mt_filter()
        decay_mode_filter = self._get_decay_mode_filter()
        prefer_macs = self._get_prefer_macs()

        if manual_text:
            # Manual list mode — tokens may be single labels ('Pb-206') or
            # inclusive mass-number ranges ('Tl206..Tl216').
            labels = []
            for tok in manual_text.split():
                labels.extend(self.svc.expand_label_token(tok))
            n_specs = len(labels)
            if n_specs == 0:
                QMessageBox.warning(self._dialog_parent(), "Range",
                                    "No isotopes in the manual list.")
                return

            def build_manual(progress=None, should_cancel=None):
                self.svc.build_universe_from_labels(
                    labels, E_eV=E_eV, mt_filter=mt_filter,
                    prefer_macs=prefer_macs, progress=progress)
                return self.svc.build_matrix(
                    flux=0.0, energy_eV=E_eV, mt_filter=mt_filter,
                    decay_mode_filter=decay_mode_filter, progress=progress)

            self.burnup_pb.setEnabled(False)
            self._build_meta = (None, None, 0.0, E_eV)
            run_with_progress(
                self._dialog_parent(),
                f"Building burnup matrix ({n_specs} isotopes)…",
                build_manual,
                on_done=self._on_built,
                on_error=self._on_build_error,
            )
            return

        # Range mode (chart selection or dropdowns)
        iz = self.left_widget.sel_initial
        lz = self.left_widget.sel_last
        if not iz or not lz:
            QMessageBox.information(
                self._dialog_parent(), "Range",
                "Select a range: click an initial then a last cell on the chart, "
                "pick both dropdowns, or type isotope names in the manual list field.")
            return

        (z1, a1), (z2, a2) = iz, lz
        n1, n2 = a1 - z1, a2 - z2
        z_from, z_to = sorted((z1, z2))
        n_from, n_to = sorted((n1, n2))

        n_specs = self.svc.count_specs(z_from, z_to, n_from, n_to)
        if n_specs == 0:
            QMessageBox.warning(self._dialog_parent(), "Range",
                                "No isotopes found in the selected range.")
            return
        if n_specs > 5000:
            QMessageBox.warning(
                self._dialog_parent(), "Range too large",
                f"The selected range contains {n_specs} isotopes. "
                "Please narrow the range to ≲ 5000 isotopes.")
            return

        def build(progress=None, should_cancel=None):
            self.svc.build_universe(
                z_from, z_to, n_from, n_to,
                E_eV=E_eV, mt_filter=mt_filter,
                prefer_macs=prefer_macs, progress=progress)
            return self.svc.build_matrix(
                flux=0.0, energy_eV=E_eV, mt_filter=mt_filter,
                decay_mode_filter=decay_mode_filter, progress=progress)

        self.burnup_pb.setEnabled(False)
        self._build_meta = (z_from, z_to, 0.0, E_eV)
        run_with_progress(
            self._dialog_parent(),
            f"Building burnup matrix ({n_specs} isotopes)…",
            build,
            on_done=self._on_built,
            on_error=self._on_build_error,
        )

    def _on_built(self, mtx):
        self.burnup_pb.setEnabled(True)
        if self.left_widget:
            self.left_widget.update_arrows()
        # Notify all registered listeners (e.g. model_results) that the matrix is fresh
        self.svc._fire_matrix_built()
        A = mtx.matrix_A
        n = A.shape[0]
        nnz = int((A != 0).sum())
        z_from, z_to, _, E_eV = self._build_meta
        mode = "Reactor physics" if self.calc_reactor_rb.isChecked() else "Astrophysics"
        QMessageBox.information(
            self._dialog_parent(), "Burnup",
            "Burnup matrix built.\n"
            f"  Mode     : {mode}\n"
            f"  Isotopes : {n}\n"
            f"  Size     : {n} × {n}\n"
            f"  Non-zeros: {nnz}\n"
            f"  Energy   : {E_eV:g} eV\n\n"
            "Results → Isotopes & Matrix tab is now populated.")

    def _on_build_error(self, tb):
        self.burnup_pb.setEnabled(True)
        QMessageBox.critical(self._dialog_parent(), "Burnup build failed", tb)

    def _on_find_cycles(self):
        if not self.svc.isotopes:
            QMessageBox.information(self._dialog_parent(), "Find Cycles",
                                    "Build burnup matrix first.")
            return
        try:
            from core.cycle_finder import CycleAnalyzer
        except ImportError:
            QMessageBox.critical(self._dialog_parent(), "Find Cycles",
                                 "core/cycle_finder.py not found.")
            return
        if getattr(self, "_cycle_search_running", False):
            return
        analyzer = CycleAnalyzer(list(self.svc.isotopes))
        matrix = getattr(self.svc, "matrix", None)
        config = getattr(matrix, "config", None)
        allowed_mt = {_REACTION_MTS[name] for name, cb in self._reac_checks.items() if cb.isChecked()}
        allowed_decay = self._get_decay_mode_filter()
        # Snapshot the adjacency before starting the worker; no GUI access there.
        for iso in analyzer._isotopes:
            edges = []
            for d in iso.getListOfDecays():
                if d.product not in analyzer._iso_set or d.branch <= 0:
                    continue
                if allowed_decay is not None and d.mode not in allowed_decay and not (d.mode == "ec/beta+" and ("ec" in allowed_decay or "beta+" in allowed_decay)):
                    continue
                if config and not config.is_decay_allowed(iso.half_life_s):
                    continue
                edges.append((d.product, d.mode, "decay", None))
            for r in iso.getListOfReactions():
                if r.product in analyzer._iso_set and r.mt in allowed_mt and r.sigma_barn > 0:
                    edges.append((r.product, mt_symbol(r.mt), "reaction", r.mt))
            analyzer._adj[iso] = edges
        by_name = {iso.name: iso for iso in analyzer._isotopes}
        for link in getattr(config, "manual_decay_links", []):
            parent = by_name.get(link.get("parent"))
            daughter = by_name.get(link.get("daughter"))
            if parent is not None and daughter is not None and link.get("branch", 1) > 0:
                analyzer._adj[parent].append((daughter, link.get("mode", "manual"), "decay", None))
        self._cycle_search_running = True
        self.find_cycle_pb.setEnabled(False)

        def search(progress=None, should_cancel=None):
            cycles = analyzer.find_all_cycles(max_cycles=200, max_length=30,
                timeout_s=5.0, should_cancel=should_cancel, progress=progress)
            return cycles, analyzer.search_status

        def finish(result):
            self._cycle_search_running = False
            self.find_cycle_pb.setEnabled(True)
            cycles, status = result
            cells = {(step.isotope.Z, step.isotope.A) for c in cycles for step in c.steps}
            if self.left_widget:
                self.left_widget.set_cycle_highlights(cells)
            self.svc.last_cycles = cycles
            self.svc._fire_cycles_found(cycles)
            QMessageBox.information(self._dialog_parent(), "Find Cycles",
                f"Found {len(cycles)} cycles involving {len(cells)} isotopes.\n"
                + (status or "Search completed.")
                + "\nLimits: 200 cycles, 30 transitions per cycle, 5 seconds.")

        def failed(message):
            self._cycle_search_running = False
            self.find_cycle_pb.setEnabled(True)
            QMessageBox.critical(self._dialog_parent(), "Find Cycles", message)

        run_with_progress(self._dialog_parent(), "Finding cycles", search,
                          on_done=finish, on_error=failed, cancelable=True)
