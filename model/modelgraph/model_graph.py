from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

import shutil

import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

import numpy as np

# Reuse the core time-axis helper so the Graph panel agrees with core.plotting.
from core.plotting import _seconds_to_unit


def _fmt_t(t_s: float) -> str:
    """Format seconds into a compact human-readable string."""
    if t_s == 0:
        return "0 s"
    for factor, unit in [(3.156e7, "yr"), (86400, "d"), (3600, "h"), (60, "min")]:
        if t_s >= factor:
            return f"{t_s / factor:.3g} {unit}"
    return f"{t_s:.3g} s"
from core.service import get_service, EqResult, RunResult


# Astrophysics-style global matplotlib settings (from the task spec).
# text.usetex needs a LaTeX install — guarded below so the app still runs
# without texlive.
_HAS_LATEX = bool(shutil.which("latex") and shutil.which("dvipng"))


def apply_astro_rcparams(use_tex: bool):
    matplotlib.rcParams.update({
        "text.usetex": bool(use_tex),
        "font.family": "serif",
        "font.serif": ["Times"],
        "xtick.direction": "in",
        "ytick.direction": "in",
        "axes.labelsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
        "text.latex.preamble": r"\usepackage{amsmath}",
    })
    matplotlib.rcParams['axes.titlesize'] = 18
    matplotlib.rcParams['xtick.major.size'] = 7
    matplotlib.rcParams['xtick.major.width'] = 1
    matplotlib.rcParams['ytick.major.size'] = 7
    matplotlib.rcParams['ytick.major.width'] = 1
    matplotlib.rcParams['xtick.minor.size'] = 4
    matplotlib.rcParams['xtick.minor.width'] = 0.5
    matplotlib.rcParams['ytick.minor.size'] = 4
    matplotlib.rcParams['ytick.minor.width'] = 0.5


# Plot kinds (also used to populate the selector).
PLOT_HEAT   = "Heat release vs time"
PLOT_CONC   = "Concentrations vs time"
PLOT_ABUND  = "Abundance vs mass number A"
PLOT_EQ_Q   = "Equilibrium heat vs flux"
PLOT_EQ_N   = "Equilibrium N_eq vs flux"
PLOT_SENSITIVITY = "Sensitivity Analysis"
PLOT_KINDS  = [PLOT_HEAT, PLOT_CONC, PLOT_ABUND, PLOT_EQ_Q, PLOT_EQ_N, PLOT_SENSITIVITY]

# Per-theme colours (figure/axes/text/ticks) — same AP.CORE tokens as
# dark.css/light.css (bg-panel/text-primary/border/accent) so the chart
# canvas matches the rest of the app. Typography + tick style stay in
# rcParams; only colours change with the theme.
_THEMES = {
    "Dark":  {"bg": "#241f1b", "fg": "#f5ede4", "grid": "#3d352e", "accent": "#e8943a"},
    "Light": {"bg": "#ffffff", "fg": "#16222e", "grid": "#d2dde6", "accent": "#2b7fc4"},
}


_GREEK_TEX = {"γ": r"\gamma", "α": r"\alpha", "β": r"\beta", "ν": r"\nu",
              "δ": r"\delta", "Φ": r"\Phi", "λ": r"\lambda", "σ": r"\sigma"}


def _latex_safe(s: str) -> str:
    """Make a label safe for matplotlib's LaTeX renderer: Greek letters (from
    reaction symbols like (n,γ)) are wrapped in math mode, and a couple of
    Unicode punctuation marks are replaced. Harmless when usetex is off."""
    out = str(s)
    for ch, tex in _GREEK_TEX.items():
        out = out.replace(ch, f"${tex}$")
    return out.replace("—", "-").replace("≈", r"$\approx$")


def _parse_iso_filter(text: str, all_names: list):
    """Return a set of allowed names from a space-separated filter string, or None=all."""
    if not text:
        return None
    tokens = text.split()
    if not tokens:
        return None
    allowed = set()
    for t in tokens:
        t_norm = t.strip().capitalize()
        for n in all_names:
            if n.lower() == t_norm.lower() or n.lower().replace("-", "") == t_norm.lower():
                allowed.add(n)
    return allowed if allowed else None


class model_graph(QObject):

    def __init__(self):
        super().__init__()
        self.__left = Left()
        self.__right = Right()

        # Right controls -> Left canvas
        self.__right.graphtype_combo.currentIndexChanged.connect(self._on_kind_changed)
        self.__right.themes_combo.currentIndexChanged.connect(self._replot)
        self.__right.plot_pb.clicked.connect(self._replot)
        self.__right.iso_filter_edit.returnPressed.connect(self._replot)
        self.__right.min_conc_edit.returnPressed.connect(self._replot)
        self.__right.abund_slider.valueChanged.connect(self._on_abund_step_changed)
        self.__right.abund_source.currentIndexChanged.connect(self._replot)
        for rb in (self.__right.abund_abs_rb, self.__right.abund_rel_rb, self.__right.abund_sn_rb):
            rb.toggled.connect(self._replot)

        # Left mode-switcher buttons -> replot
        for btn in self.__left._conc_btns + self.__left._heat_btns + self.__left._eq_n_btns:
            btn.clicked.connect(self._replot)

        # Follow the app-wide theme switch (View menu) so the plot colours
        # change together with the rest of the UI from one switch.
        svc = get_service()
        svc.register_theme_changed_callback(self._on_global_theme_changed)
        self._sync_theme_combo(svc.theme)

        self._replot()

    def getLeftWidget(self):
        return self.__left

    def getRightLayout(self):
        return self.__right

    def _on_kind_changed(self, *_):
        kind = self.__right.graphtype_combo.currentText()
        self.__left.set_plot_kind(kind)
        svc = get_service()
        self.__right.sync_abund_controls(kind, svc)
        self._replot()

    def _on_abund_step_changed(self, _):
        svc = get_service()
        self.__right._update_abund_time_label(svc)
        self._replot()

    def _on_global_theme_changed(self, theme: str) -> None:
        self._sync_theme_combo(theme)
        self._replot()

    def _sync_theme_combo(self, theme: str) -> None:
        target = "Light" if theme == "light" else "Dark"
        combo = self.__right.themes_combo
        if combo.currentText() != target:
            combo.blockSignals(True)
            combo.setCurrentText(target)
            combo.blockSignals(False)

    def _replot(self, *_):
        kind  = self.__right.graphtype_combo.currentText()
        theme = self.__right.themes_combo.currentText()
        iso_filter = self.__right.iso_filter_edit.text().strip()
        conc_mode = self.__left.current_conc_mode()
        heat_unit = self.__left.current_heat_unit()
        eq_n_mode = self.__left.current_eq_n_mode()
        min_conc_text = self.__right.min_conc_edit.text().strip()
        svc = get_service()
        self.__right.sync_abund_controls(kind, svc)
        abund_step = self.__right.abund_current_step()
        abund_mode = self.__right.abund_current_mode()
        self.__left.plot(kind, theme, iso_filter, conc_mode=conc_mode,
                         heat_unit=heat_unit, eq_n_mode=eq_n_mode,
                         min_conc_text=min_conc_text,
                         abund_step=abund_step, abund_mode=abund_mode,
                         exposure_averaged=self.__right.abund_source.currentIndex() == 1)


class CreateCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.figure = Figure(figsize=(7, 5), dpi=100)
        super().__init__(self.figure)


class Left(QWidget):
    def __init__(self):
        super().__init__()

        # LaTeX is opt-in: enable only if available, and remember if it fails.
        self._use_tex = _HAS_LATEX
        apply_astro_rcparams(self._use_tex)

        self.canvas = CreateCanvas(self)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.ax = self.canvas.figure.add_subplot(111)

        # -- Concentration mode switcher (shown for PLOT_CONC only) --
        self._conc_mode_row = QWidget()
        ch = QHBoxLayout(self._conc_mode_row)
        ch.setContentsMargins(2, 2, 2, 2)
        ch.addStretch(0)
        self._btn_abs = QPushButton("Absolute")
        self._btn_rel = QPushButton("Relative")
        self._btn_sn  = QPushButton("σ·N")
        self._conc_btns = [self._btn_abs, self._btn_rel, self._btn_sn]
        for b in self._conc_btns:
            b.setProperty("graphMode", True)
            b.setCheckable(True)
            b.setFlat(False)
            ch.addWidget(b)
        self._btn_abs.setChecked(True)
        self._conc_mode_row.setVisible(False)

        # -- Heat unit switcher (shown for PLOT_HEAT and PLOT_EQ_Q) --
        self._heat_unit_row = QWidget()
        hh = QHBoxLayout(self._heat_unit_row)
        hh.setContentsMargins(2, 2, 2, 2)
        hh.addStretch(0)
        self._btn_mevs = QPushButton("MeV/s")
        self._btn_wcm3 = QPushButton("W/cm³")
        self._heat_btns = [self._btn_mevs, self._btn_wcm3]
        for b in self._heat_btns:
            b.setProperty("graphMode", True)
            b.setCheckable(True)
            b.setFlat(False)
            hh.addWidget(b)
        self._btn_mevs.setChecked(True)
        self._heat_unit_row.setVisible(False)

        # -- Equilibrium-N mode switcher (shown for PLOT_EQ_N only) --
        self._eq_n_mode_row = QWidget()
        en = QHBoxLayout(self._eq_n_mode_row)
        en.setContentsMargins(2, 2, 2, 2)
        en.addStretch(0)
        self._btn_eq_rel  = QPushButton("N/N_tot")
        self._btn_eq_abs  = QPushButton("Absolute")
        self._btn_eq_sn   = QPushButton("σ·N")
        self._eq_n_btns = [self._btn_eq_rel, self._btn_eq_abs, self._btn_eq_sn]
        for b in self._eq_n_btns:
            b.setProperty("graphMode", True)
            b.setCheckable(True)
            b.setFlat(False)
            en.addWidget(b)
        self._btn_eq_rel.setChecked(True)
        self._eq_n_mode_row.setVisible(False)

        # Exclusive toggle for conc buttons
        self._btn_abs.clicked.connect(lambda: self._set_conc_mode("absolute"))
        self._btn_rel.clicked.connect(lambda: self._set_conc_mode("relative"))
        self._btn_sn.clicked.connect(lambda: self._set_conc_mode("sigma_n"))
        self._btn_mevs.clicked.connect(lambda: self._set_heat_unit("mev"))
        self._btn_wcm3.clicked.connect(lambda: self._set_heat_unit("wcm3"))
        self._btn_eq_rel.clicked.connect(lambda: self._set_eq_n_mode("relative"))
        self._btn_eq_abs.clicked.connect(lambda: self._set_eq_n_mode("absolute"))
        self._btn_eq_sn.clicked.connect(lambda: self._set_eq_n_mode("sigma_n"))

        self._conc_mode_val = "absolute"
        self._heat_unit_val = "mev"
        self._eq_n_mode_val = "relative"

        # Scrollable, multi-column legend panel (handles 60+ isotope curves).
        self.legend_widget = QWidget()
        self.legend_widget.setObjectName('graph_legend_widget')
        self.legend_layout = QGridLayout(self.legend_widget)
        self.legend_layout.setAlignment(Qt.AlignTop)
        self.legend_layout.setContentsMargins(4, 4, 4, 4)
        self.legend_layout.setHorizontalSpacing(10)
        self.legend_layout.setVerticalSpacing(2)

        self.legend_scroll = QScrollArea()
        self.legend_scroll.setObjectName('graph_legend_scroll')
        self.legend_scroll.setWidgetResizable(True)
        self.legend_scroll.setWidget(self.legend_widget)
        self.legend_scroll.setMinimumWidth(160)
        self.legend_scroll.setMaximumWidth(280)
        self.legend_scroll.setVisible(False)   # only shown when there's a legend to draw

        body_splitter = QSplitter(Qt.Horizontal)
        body_splitter.addWidget(self.canvas)
        body_splitter.addWidget(self.legend_scroll)
        body_splitter.setStretchFactor(0, 4)
        body_splitter.setStretchFactor(1, 1)

        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.addWidget(self._conc_mode_row)
        root.addWidget(self._heat_unit_row)
        root.addWidget(self._eq_n_mode_row)
        root.addWidget(self.toolbar)
        root.addWidget(body_splitter, stretch=1)

    def set_plot_kind(self, kind: str):
        self._conc_mode_row.setVisible(kind == PLOT_CONC)
        self._heat_unit_row.setVisible(kind in (PLOT_HEAT, PLOT_EQ_Q))
        self._eq_n_mode_row.setVisible(kind == PLOT_EQ_N)

    def _update_legend(self, entries, colors, n_cols: int = 2):
        # Qt-native legend, not matplotlib's — matplotlib's can't scroll.
        layout = self.legend_layout
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()

        self.legend_scroll.setVisible(bool(entries))
        if not entries:
            return

        fg = colors["fg"]
        for i, (name, color) in enumerate(entries):
            row, col = divmod(i, n_cols)
            label = QLabel(
                f'<span style="color:{color};">■</span> '
                f'<span style="color:{fg};">{name}</span>')
            label.setStyleSheet("background: transparent;")
            layout.addWidget(label, row, col)

    def _set_conc_mode(self, mode: str):
        self._conc_mode_val = mode
        for btn, m in zip(self._conc_btns, ("absolute", "relative", "sigma_n")):
            btn.setChecked(m == mode)

    def _set_heat_unit(self, unit: str):
        self._heat_unit_val = unit
        for btn, u in zip(self._heat_btns, ("mev", "wcm3")):
            btn.setChecked(u == unit)

    def _set_eq_n_mode(self, mode: str):
        self._eq_n_mode_val = mode
        for btn, m in zip(self._eq_n_btns, ("relative", "absolute", "sigma_n")):
            btn.setChecked(m == mode)

    def current_conc_mode(self) -> str:
        return self._conc_mode_val

    def current_heat_unit(self) -> str:
        return self._heat_unit_val

    def current_eq_n_mode(self) -> str:
        return self._eq_n_mode_val

    # -- public ------------------------------------------------------------
    def plot(self, kind: str, theme: str, iso_filter: str = "",
             conc_mode: str = "absolute", heat_unit: str = "mev",
             eq_n_mode: str = "relative", min_conc_text: str = "",
             abund_step: int = -1, abund_mode: str = "absolute",
             exposure_averaged: bool = False):
        colors = _THEMES.get(theme, _THEMES["Dark"])
        fig = self.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self.ax = ax
        svc = get_service()
        # Default to no side legend; concentration/eq-N plots populate it below.
        self._update_legend([], colors)

        # Equilibrium plots read from last_eq
        if kind in (PLOT_EQ_Q, PLOT_EQ_N):
            eq = svc.last_eq
            if eq is None or not eq.flux_list:
                self._placeholder(ax, colors,
                                  "No equilibrium results yet.\n"
                                  "Run an Equilibrium scan in the Results panel.")
                self._finish(fig, ax, colors)
                return
            try:
                if kind == PLOT_EQ_Q:
                    self._plot_eq_heat(ax, eq, colors, heat_unit)
                else:
                    self._plot_eq_n(ax, eq, colors, iso_filter, min_conc_text, eq_n_mode)
            except Exception as e:
                self._placeholder(ax, colors, f"Plot error:\n{e}")
            self._finish(fig, ax, colors)
            return

        # Sensitivity analysis / library spread reads from last_sensitivity
        if kind == PLOT_SENSITIVITY:
            sens = getattr(svc, "last_sensitivity", None)
            if sens is None or (not sens.rows and not sens.library_rows):
                self._placeholder(ax, colors,
                                  "No sensitivity results yet.\n"
                                  "Run a Sensitivity analysis in the Results panel.")
                self._finish(fig, ax, colors)
                return
            try:
                self._plot_sensitivity(ax, sens, colors)
            except Exception as e:
                self._placeholder(ax, colors, f"Plot error:\n{e}")
            self._finish(fig, ax, colors)
            return

        res = svc.last_run
        if res is None or not res.times:
            self._placeholder(ax, colors,
                              "No results yet.\nRun a calculation in the Results panel.")
            self._finish(fig, ax, colors)
            return

        try:
            if kind == PLOT_HEAT:
                self._plot_heat(ax, res, colors, heat_unit)
            elif kind == PLOT_ABUND:
                if exposure_averaged:
                    from model.modelgraph.exposure_plot import plot_exposure
                    result = getattr(res, "exposure_result", None)
                    if result is None:
                        self._placeholder(ax, colors, "Run Trajectory averaging in Results first.")
                    else:
                        plot_exposure(ax, result, abund_mode, colors)
                else:
                    self._plot_abundance(ax, res, colors, abund_step, abund_mode)
            else:
                self._plot_concentrations(ax, res, colors, iso_filter, conc_mode, min_conc_text)
        except Exception as e:
            self._placeholder(ax, colors, f"Plot error:\n{e}")

        self._finish(fig, ax, colors)

    # -- individual plots --------------------------------------------------
    def _plot_heat(self, ax, res, colors, heat_unit: str = "mev"):
        x = list(res.times)
        y = list(res.heat)

        if heat_unit == "wcm3":
            svc = get_service()
            rho = getattr(svc, "last_rho", 0.0)
            if rho > 0 and res.isotopes and res.initial:
                import re as _re
                dominant = max(res.initial.items(), key=lambda kv: kv[1])[0] if res.initial else None
                M = 1.0
                if dominant:
                    m = _re.search(r"-(\d+)", dominant)
                    if m:
                        M = float(m.group(1))
                factor = 1.602e-13 * rho * 6.02214076e23 / M
                y = [v * factor for v in y]
                ylabel = r"Decay heat [W/cm³]"
            else:
                ylabel = r"Decay heat [MeV/s per nucleus]  (set ρ in Results for W/cm³)"
        else:
            ylabel = r"Decay heat [MeV/s per nucleus]"

        ax.plot(x, y, color=colors["accent"], linewidth=1.8, marker="s",
                markersize=4, markevery=max(1, len(x) // 20))
        ax.set_xlabel(r"Time [s]")
        ax.set_ylabel(ylabel)
        ax.set_title("Mixture decay-heat profile")
        self._maybe_log(ax, x, y)

    def _plot_concentrations(self, ax, res, colors, iso_filter: str = "",
                             conc_mode: str = "absolute", min_conc_text: str = ""):
        M = np.array(res.trajectories, dtype=float)        # (T, N)
        x = list(res.times)
        allowed = _parse_iso_filter(iso_filter, res.names)
        # Restricted run: show only the requested target isotopes.
        targets = set(getattr(res, "targets", []) or [])
        if targets:
            allowed = targets if allowed is None else (allowed & targets)

        # Build sigma lookup for σ·N mode
        sigma_by_name = {}
        if conc_mode == "sigma_n":
            for iso in res.isotopes:
                for rx in iso.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if mt == 102:
                        sigma_by_name[iso.name] = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                        break

        if conc_mode == "relative":
            sums = M.sum(axis=1, keepdims=True)
            sums[sums == 0] = 1.0
            Y = M / sums
            ylabel = r"Relative atomic fraction"
        elif conc_mode == "sigma_n":
            sigma_arr = np.array([sigma_by_name.get(n, 0.0) for n in res.names])
            Y = M * sigma_arr[np.newaxis, :]
            ylabel = r"$\sigma \cdot N$ [barn · atoms]"
        else:
            Y = M
            ylabel = r"Atom count $N$"

        peak = np.max(np.abs(Y), axis=0)
        try:
            threshold = float(min_conc_text) if min_conc_text else 1e-6 * (peak.max() or 1.0)
        except ValueError:
            threshold = 1e-6 * (peak.max() or 1.0)
        legend_entries = []
        for j, name in enumerate(res.names):
            if allowed is not None and name not in allowed:
                continue
            if peak[j] < threshold and not (targets and name in targets):
                continue   # always draw an explicitly requested target
            line, = ax.plot(x, Y[:, j], linewidth=1.2, label=name)
            legend_entries.append((name, line.get_color()))

        ax.set_xlabel(r"Time [s]")
        ax.set_ylabel(ylabel)
        mode_label = {"absolute": "Absolute", "relative": "Relative", "sigma_n": "Capture-weighted abundance"}.get(conc_mode, conc_mode)
        ax.set_title(f"Isotope concentration evolution ({mode_label})")
        ax.set_xscale("log")
        if conc_mode != "sigma_n":
            ax.set_yscale("log")
            ax.set_ylim(bottom=1e-12 * (np.max(Y) or 1.0))
        self._update_legend(legend_entries, colors)

    def _plot_eq_heat(self, ax, eq: EqResult, colors, heat_unit: str = "mev"):
        import re as _re
        x = eq.flux_list
        y = list(eq.Q_eq)

        if heat_unit == "wcm3":
            svc = get_service()
            rho = getattr(svc, "last_rho", 0.0)
            if rho > 0 and eq.isotopes and eq.initial:
                dominant = max(eq.initial.items(), key=lambda kv: kv[1])[0] if eq.initial else None
                M = 1.0
                if dominant:
                    m = _re.search(r"-(\d+)", dominant)
                    if m:
                        M = float(m.group(1))
                factor = 1.602e-13 * rho * 6.02214076e23 / M
                y = [v * factor for v in y]
                ylabel = r"Equilibrium heat [W/cm³]"
            else:
                ylabel = r"Equilibrium heat [MeV/s per nucleus]  (set ρ in Results for W/cm³)"
        else:
            ylabel = r"Equilibrium heat $Q_{eq}$ [MeV/s per nucleus]"

        ax.plot(x, y, color=colors["accent"], linewidth=1.8, marker="o", markersize=5)
        ax.set_xlabel(r"Neutron flux $\Phi$ [n/cm$^2$/s]")
        ax.set_ylabel(ylabel)
        ax.set_title("Equilibrium heat release vs neutron flux")
        self._maybe_log(ax, x, y)

    def _plot_eq_n(self, ax, eq: EqResult, colors, iso_filter: str = "",
                   min_conc_text: str = "", eq_n_mode: str = "relative"):
        x = eq.flux_list
        allowed = _parse_iso_filter(iso_filter, eq.names)
        M_eq = np.array(eq.N_eq, dtype=float)       # (n_flux, n_iso)

        if eq_n_mode == "sigma_n":
            sigma_by_name = {}
            for iso in eq.isotopes:
                for rx in iso.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if mt == 102:
                        sigma_by_name[iso.name] = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                        break
            sigma_arr = np.array([sigma_by_name.get(n, 0.0) for n in eq.names])
            Y = M_eq * sigma_arr[np.newaxis, :]
            ylabel = r"$\sigma \cdot N_{eq}$ [barn · atoms]"
        elif eq_n_mode == "absolute":
            Y = M_eq
            ylabel = r"Equilibrium concentration $N_{eq}$"
        else:  # relative
            sums = M_eq.sum(axis=1, keepdims=True)
            sums[sums == 0] = 1.0
            Y = M_eq / sums
            ylabel = r"Equilibrium fraction $N_{eq}/N_{tot}$"

        peak = np.max(np.abs(Y), axis=0)
        try:
            threshold = float(min_conc_text) if min_conc_text else 1e-8 * (peak.max() or 1.0)
        except ValueError:
            threshold = 1e-8 * (peak.max() or 1.0)
        legend_entries = []
        for j, name in enumerate(eq.names):
            if allowed is not None and name not in allowed:
                continue
            if peak[j] < threshold:
                continue
            line, = ax.plot(x, Y[:, j], linewidth=1.2, label=name, marker=".")
            legend_entries.append((name, line.get_color()))
        ax.set_xlabel(r"Neutron flux $\Phi$ [n/cm$^2$/s]")
        ax.set_ylabel(ylabel)
        mode_label = {"absolute": "Absolute", "relative": "Relative", "sigma_n": "Capture-weighted abundance"}.get(eq_n_mode, eq_n_mode)
        ax.set_title(f"Equilibrium concentration vs neutron flux ({mode_label})")
        self._maybe_log(ax, x, [v for col in Y.T for v in col if v > 0] or [1])
        if ax.get_xscale() != "log" and max(x) / max(min(x), 1e-30) > 10:
            ax.set_xscale("log")
        if eq_n_mode != "sigma_n":
            ax.set_yscale("log")
        self._update_legend(legend_entries, colors)

    def _plot_abundance(self, ax, res, colors, step: int = -1, mode: str = "absolute"):
        if not res.trajectories:
            self._placeholder(ax, colors, "No trajectories.")
            return
        idx = step if 0 <= step < len(res.trajectories) else len(res.trajectories) - 1
        N_vec = res.trajectories[idx]
        t_val = res.times[idx] if idx < len(res.times) else 0.0

        sigma_by_name: dict = {}
        if mode == "sigma_n":
            for iso in res.isotopes:
                for rx in iso.getListOfReactions():
                    mt = rx.getId() if hasattr(rx, "getId") else getattr(rx, "mt", -1)
                    if mt == 102:
                        sigma_by_name[iso.name] = float(getattr(rx, "sigma_barn", 0.0) or 0.0)
                        break

        total = float(np.sum(N_vec)) or 1.0

        by_A: dict = {}
        for name, N_i_raw in zip(res.names, N_vec):
            N_i = float(N_i_raw)
            if N_i <= 0:
                continue
            if mode == "relative":
                val = N_i / total
            elif mode == "sigma_n":
                val = sigma_by_name.get(name, 0.0) * N_i
            else:
                val = N_i
            A = self._mass_of(name)
            if A is not None:
                by_A[A] = by_A.get(A, 0.0) + val

        if not by_A:
            self._placeholder(ax, colors, "No non-zero abundances to show.")
            return

        As   = sorted(by_A)
        vals = [by_A[a] for a in As]

        ax.plot(As, vals, '-', color=colors["accent"], linewidth=1.2)
        ax.plot(As, vals, 'o', color=colors["accent"], markersize=3)
        ax.set_xlabel(r"Mass number $A$")
        _ylabels = {
            "relative": r"Relative abundance  $N(A)/N_{tot}$",
            "sigma_n":  r"$\sigma \cdot N(A)$  [barn · atoms]",
        }
        ax.set_ylabel(_ylabels.get(mode, r"$N(A)$  [atoms]"))
        mode_label = {"absolute": "Absolute", "relative": "Relative", "sigma_n": "Capture-weighted abundance"}.get(mode, mode)
        ax.set_title(f"Abundance vs $A$   (t = {_fmt_t(t_val)},  {mode_label})")
        pos_vals = [v for v in vals if v > 0]
        if pos_vals and max(pos_vals) / max(min(pos_vals), 1e-300) > 100:
            ax.set_yscale("log")

    def _plot_sensitivity(self, ax, sens, colors):
        """Sensitivity bars for σ/λ/flux/energy (sorted by |S|); for a library-only
        run, a bar chart of N_i (or Q) per library with the mean line."""
        rows = sens.ranked()
        if rows:
            top = rows[:15][::-1]   # largest |S| ends up at the top
            labels = [_latex_safe(r.parameter) for r in top]
            vals = [r.S for r in top]
            ypos = list(range(len(top)))
            bar_colors = [colors["accent"] if v >= 0 else "#c0504d" for v in vals]
            ax.barh(ypos, vals, color=bar_colors)
            ax.set_yticks(ypos)
            ax.set_yticklabels(labels, fontsize=8)
            ax.axvline(0, color=colors["fg"], linewidth=0.8)
            ax.set_xlabel("Sensitivity coefficient  S")
            ax.set_title(_latex_safe(f"Sensitivity Analysis — {sens.target_desc}"))
        elif sens.library_rows:
            libs = [_latex_safe(lr.library) for lr in sens.library_rows]
            has_ni = any(lr.N_i is not None for lr in sens.library_rows)
            vals = [(lr.N_i if (has_ni and lr.N_i is not None) else lr.Q)
                    for lr in sens.library_rows]
            xs = list(range(len(libs)))
            ax.bar(xs, vals, color=colors["accent"])
            ax.set_xticks(xs)
            ax.set_xticklabels(libs, rotation=30, ha="right", fontsize=8)
            ax.set_ylabel(r"$N_i$" if has_ni else "Q [MeV/s]")
            ax.set_title(_latex_safe(f"Library spread — {sens.target_desc}"))
            if vals:
                ax.axhline(float(np.mean(vals)), color="#c0504d",
                           linewidth=1.0, linestyle="--", label="mean")
                leg = ax.legend(fontsize=9)
                self._style_legend(leg, colors)

    # -- helpers -----------------------------------------------------------
    @staticmethod
    def _mass_of(name: str):
        # "Pb-206" / "Pb-209m1" -> 206 / 209
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

    @staticmethod
    def _maybe_log(ax, x, y):
        xs = [v for v in x if v > 0]
        ys = [v for v in y if v > 0]
        if xs and max(xs) / min(xs) > 100:
            ax.set_xscale("log")
        if ys and max(ys) / max(min(ys), 1e-300) > 100:
            ax.set_yscale("log")

    def _placeholder(self, ax, colors, text):
        ax.text(0.5, 0.5, text, ha="center", va="center",
                transform=ax.transAxes, color=colors["fg"], fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])

    def _style_legend(self, leg, colors):
        if leg is None:
            return
        leg.get_frame().set_facecolor(colors["bg"])
        leg.get_frame().set_edgecolor(colors["grid"])
        for txt in leg.get_texts():
            txt.set_color(colors["fg"])

    def _style_axes(self, ax, colors):
        ax.set_facecolor(colors["bg"])
        for spine in ax.spines.values():
            spine.set_color(colors["fg"])
        ax.tick_params(colors=colors["fg"], which="both", labelsize=13)
        ax.xaxis.label.set_color(colors["fg"])
        ax.yaxis.label.set_color(colors["fg"])
        ax.xaxis.label.set_fontsize(14)
        ax.yaxis.label.set_fontsize(14)
        ax.title.set_color(colors["fg"])
        ax.title.set_fontsize(16)
        ax.title.set_fontweight("semibold")
        ax.grid(True, alpha=0.25, color=colors["grid"], which="both")

    def _finish(self, fig, ax, colors):
        fig.patch.set_facecolor(colors["bg"])
        self._style_axes(ax, colors)
        try:
            fig.tight_layout()
        except Exception:
            pass
        # Draw, with a graceful fallback if a LaTeX render fails at runtime.
        try:
            self.canvas.draw()
        except Exception:
            if self._use_tex:
                self._use_tex = False
                apply_astro_rcparams(False)
                from matplotlib.text import Text
                for text in fig.findobj(Text):
                    text.set_usetex(False)
                try:
                    self.canvas.draw()
                except Exception:
                    pass


class Right(QVBoxLayout):

    def __init__(self):
        super().__init__()

        # ---- which data to display
        self.type_hlayout = QHBoxLayout()
        self.graphtype_label = QLabel('Graph data')
        self.graphtype_label.setObjectName('graphtype_label')
        self.graphtype_combo = QComboBox()
        self.graphtype_combo.setObjectName('graphtype_combo')
        self.graphtype_combo.addItems(PLOT_KINDS)
        self.type_hlayout.addWidget(self.graphtype_label)
        self.type_hlayout.addWidget(self.graphtype_combo)

        # ---- theme (figure colours)
        self.themes_hlayout = QHBoxLayout()
        self.themes_label = QLabel('Theme')
        self.themes_label.setObjectName('themes_label')
        self.themes_combo = QComboBox()
        self.themes_combo.setObjectName('themes_combo')
        self.themes_combo.addItems(['Dark', 'Light'])
        self.themes_hlayout.addWidget(self.themes_label)
        self.themes_hlayout.addWidget(self.themes_combo)

        # ---- isotope filter (for Concentrations and Equilibrium N plots)
        self.iso_filter_label = QLabel('Isotope filter:')
        self.iso_filter_label.setObjectName('iso_filter_label')
        self.iso_filter_edit = QLineEdit()
        self.iso_filter_edit.setObjectName('iso_filter_edit')
        self.iso_filter_edit.setPlaceholderText('empty = show all  |  e.g. Pb-208 Tl-208')

        # ---- minimum concentration/value to display (drop near-zero curves)
        self.min_conc_label = QLabel('Min value to display (peak, empty = auto):')
        self.min_conc_label.setObjectName('min_conc_label')
        self.min_conc_edit = QLineEdit()
        self.min_conc_edit.setObjectName('min_conc_edit')
        self.min_conc_edit.setPlaceholderText('e.g. 1e-10')

        # ---- (re)plot
        self.plot_pb = QPushButton('Plot / Refresh')
        self.plot_pb.setObjectName('plot_pb')

        # ---- abundance-specific controls (time step slider + mode radios)
        self.abund_widget = QWidget()
        self.abund_widget.setObjectName('abund_controls_widget')
        abund_v = QVBoxLayout(self.abund_widget)
        abund_v.setContentsMargins(0, 0, 0, 0)
        abund_v.setSpacing(4)

        abund_step_row = QHBoxLayout()
        abund_step_row.addWidget(QLabel("Time step:"))
        self.abund_slider = QSlider(Qt.Horizontal)
        self.abund_slider.setObjectName("abund_slider")
        self.abund_slider.setMinimum(0)
        self.abund_slider.setMaximum(0)
        self.abund_time_label = QLabel("t = —")
        self.abund_time_label.setObjectName("abund_time_label")
        abund_step_row.addWidget(self.abund_slider, 1)
        abund_step_row.addWidget(self.abund_time_label)
        abund_v.addLayout(abund_step_row)

        self.abund_source = QComboBox()
        self.abund_source.addItems(["Time slice", "Averaged results"])
        abund_v.addWidget(QLabel("Abundance source"))
        abund_v.addWidget(self.abund_source)
        abund_mode_row = QHBoxLayout()
        abund_mode_row.addWidget(QLabel("Mode:"))
        self.abund_abs_rb = QRadioButton("Absolute")
        self.abund_rel_rb = QRadioButton("Relative")
        self.abund_sn_rb  = QRadioButton("σN")
        self.abund_abs_rb.setObjectName("abund_abs_rb")
        self.abund_rel_rb.setObjectName("abund_rel_rb")
        self.abund_sn_rb.setObjectName("abund_sn_rb")
        self.abund_abs_rb.setChecked(True)
        for rb in (self.abund_abs_rb, self.abund_rel_rb, self.abund_sn_rb):
            abund_mode_row.addWidget(rb)
        abund_v.addLayout(abund_mode_row)

        self.abund_widget.setVisible(False)

        self.note_label = QLabel(
            "Evolution plots: source = last run (Results → Evolution).\n"
            "Equilibrium plots: source = last scan (Results → Equilibrium)."
            + ("" if _HAS_LATEX else "\nLaTeX not found — using mathtext fallback."))
        self.note_label.setObjectName('graph_note_label')
        self.note_label.setWordWrap(True)

        self.addLayout(self.type_hlayout)
        self.addSpacing(8)
        self.addLayout(self.themes_hlayout)
        self.addSpacing(8)
        self.addWidget(self.iso_filter_label)
        self.addWidget(self.iso_filter_edit)
        self.addSpacing(8)
        self.addWidget(self.min_conc_label)
        self.addWidget(self.min_conc_edit)
        self.addSpacing(8)
        self.addWidget(self.plot_pb)
        self.addSpacing(4)
        self.addWidget(self.abund_widget)
        self.addSpacing(8)
        self.addWidget(self.note_label)
        self.setAlignment(Qt.AlignTop)

    def sync_abund_controls(self, kind: str, svc) -> None:
        is_abund = (kind == PLOT_ABUND)
        self.abund_widget.setVisible(is_abund)
        averaged = self.abund_source.currentIndex() == 1
        self.abund_slider.setEnabled(not averaged)

        if is_abund:
            res = svc.last_run
            if res and res.times:
                n = len(res.times)
                self.abund_slider.blockSignals(True)
                self.abund_slider.setMaximum(n - 1)
                cur = self.abund_slider.value()
                if cur >= n:
                    self.abund_slider.setValue(n - 1)
                self.abund_slider.blockSignals(False)
                self._update_abund_time_label(svc)

    def _update_abund_time_label(self, svc) -> None:
        res = svc.last_run
        if res and res.times:
            idx = self.abund_slider.value()
            if 0 <= idx < len(res.times):
                self.abund_time_label.setText(_fmt_t(res.times[idx]))

    def abund_current_step(self) -> int:
        return self.abund_slider.value()

    def abund_current_mode(self) -> str:
        if self.abund_rel_rb.isChecked():
            return "relative"
        if self.abund_sn_rb.isChecked():
            return "sigma_n"
        return "absolute"
