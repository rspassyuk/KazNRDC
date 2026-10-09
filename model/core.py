# core.py
# Application entry point.


# All nuclear-data initialisation (xsdir discovery, ENDF library selection,
# provider build) lives in core.service.init_service (the shared facade).
# This file only wires the GUI: splash, window, dashboard, toolbar, dynamic
# model loading, and the application icon.

import importlib.util
import inspect
import os
import sys
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import (
    QAction,
    QApplication,
    QFileDialog,
    QFrame,
    QGroupBox,
    QMainWindow,
    QMessageBox,
    QSizePolicy,
    QSplitter,
    QToolBar,
    QVBoxLayout,
)

# --- Make the project root + this folder importable, so the dynamically
#     loaded models can `from core.service import ...` and `from GUIcomp import *`
#     regardless of the current working directory (Ubuntu, forward slashes). ---
from pathlib import Path

MODELS_DIR = Path(__file__).resolve().parent          # …/NuMatRx/model
PROJECT_ROOT = MODELS_DIR.parent                       # …/NuMatRx  (xsdir, icons)
# Order matters: PROJECT_ROOT must precede MODELS_DIR so the `core/` package
# wins over this file (which is itself named core.py) for `import core.*`.
for _p in (str(MODELS_DIR), str(PROJECT_ROOT)):
    if _p in sys.path:
        sys.path.remove(_p)
sys.path.insert(0, str(MODELS_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

from GUIcomp import *  # noqa: F401,F403
from core.service import init_service, get_service


# -------------------------
# Asset resolution (Ubuntu-safe, pathlib)
# -------------------------
def app_icon_path() -> str:
    return logo_path("dark")


def logo_path(theme="dark") -> str:
    p = MODELS_DIR / "resources" / "branding" / f"numatrx_logo_{theme}.png"
    return str(p) if p.is_file() else ""


def splash_image_path() -> str:
    p = Path(logo_path("dark")) if logo_path("dark") else PROJECT_ROOT / "start_image.png"
    return str(p) if p.is_file() else ""


# -------------------------
# Startup splash with a real-progress bar (Task 5)
# -------------------------
class SplashScreen(QWidget):
    """Frameless splash: start_image.png + status line + progress bar.

    Progress is driven by the actual initialisation stages (database load,
    provider build, model loading, interface prep) via set_progress().
    """

    def __init__(self, image_path: str = ""):
        super().__init__(None, Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint | Qt.SplashScreen)
        self.setObjectName("SplashScreen")
        self.setStyleSheet(
            "QWidget#SplashScreen { background-color: #1c1f22; border: 1px solid orange; }"
            "QWidget#SplashBar { background-color: #1c1f22; }"
            "QLabel#SplashStatus { color: #e6e6e6; font-size: 13px; font-family: 'Italiana'; }"
            "QProgressBar#SplashProgress { background-color: #2c2f33; border: 1px solid #444;"
            " border-radius: 3px; }"
            "QProgressBar#SplashProgress::chunk { background-color: orange; border-radius: 3px; }"
        )

        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        if image_path:
            pm = QPixmap(image_path)
            if not pm.isNull():
                pm = pm.scaled(QSize(680, 440), Qt.KeepAspectRatio, Qt.SmoothTransformation)
                self.image_label.setPixmap(pm)
        v.addWidget(self.image_label)

        bar_box = QWidget()
        bar_box.setObjectName("SplashBar")
        b = QVBoxLayout(bar_box)
        b.setContentsMargins(18, 8, 18, 14)
        b.setSpacing(6)
        self.status = QLabel("Starting…")
        self.status.setObjectName("SplashStatus")
        self.progress = QProgressBar()
        self.progress.setObjectName("SplashProgress")
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setTextVisible(False)
        self.progress.setFixedHeight(6)   # thin bar (Task 2)
        b.addWidget(self.status)
        b.addWidget(self.progress)
        v.addWidget(bar_box)

        self.adjustSize()
        screen = QApplication.primaryScreen()
        if screen is not None:
            self.move(screen.availableGeometry().center() - self.rect().center())

    def set_progress(self, frac: float, message: str = None):
        self.progress.setValue(max(0, min(100, int(frac * 100))))
        if message:
            self.status.setText(message)
        QApplication.processEvents()

    def finish_with(self, window):
        self.set_progress(1.0, "Ready.")
        self.close()


# -------------------------
# Main window
# -------------------------
class MainWindow(QMainWindow):
    def __init__(self, splash=None):
        super().__init__()
        self._splash = splash
        self.init_window()

    def _progress(self, frac, message=None):
        """Forward init progress to the splash screen, if one is present."""
        if self._splash is not None:
            self._splash.set_progress(frac, message)

    def init_window(self):
        # Theme first (optional)
        self.apply_theme("dark.css")

        self.setWindowTitle("NuMatRx")
        icon_path = app_icon_path()
        if icon_path:
            self.setWindowIcon(QIcon(icon_path))

        # Screen size
        screen = QApplication.primaryScreen()
        geom = screen.availableGeometry()
        self.window_width = int(geom.width())
        self.window_height = int(geom.height())
        self.resize(self.window_width, self.window_height)

        # Central widget
        self.center_widget = QGroupBox("___________________")
        self.center_widget.setObjectName("The_center_widget")
        self.center_widget.setContentsMargins(0, 0, 0, 0)

        self.center_layout = QVBoxLayout()
        self.center_layout.setContentsMargins(1, 4, 4, 0)
        self.center_widget.setLayout(self.center_layout)
        self.setCentralWidget(self.center_widget)

        # Main widget area
        self.MainWidget = QGroupBox("T")
        self.MainWidget.setObjectName("MainWidget")
        self.MainLayout = QVBoxLayout()
        self.MainWidget.setLayout(self.MainLayout)

        # Map for pages
        self.MainWidgetsMap = dict()

        # Dashboard
        self.dashboard = None
        self.init_dashboard()

        # Terminal (singleton from GUIcomp). Wire stdout/stderr now so the
        # nuclear-data init prints land in the in-app terminal.
        self.terminal = GUITerminal(self)
        self.terminal.setMinimumHeight(100)
        self.terminal.hide()
        sys.stdout = self.terminal
        sys.stderr = self.terminal

        # Splitter
        self.mainhor_splitter = QSplitter(Qt.Vertical)
        self.mainhor_splitter.setHandleWidth(6)
        self.mainhor_splitter.setObjectName("mainhorsplitter")
        self.mainhor_splitter.addWidget(self.MainWidget)
        self.mainhor_splitter.addWidget(self.terminal)
        self.center_layout.addWidget(self.mainhor_splitter)

        # Menu / line
        self.init_menu()
        self.add_menu_line()

        # ---- CRITICAL: bring up the nuclear-data environment BEFORE models load.
        # The shared core.service is the single backend facade every model calls.
        try:
            init_service(PROJECT_ROOT, progress=self._progress)
        except Exception as e:
            QMessageBox.critical(
                self,
                "Nuclear-data init failed",
                "Nuclear-data environment failed to initialise.\n\n"
                f"{e}\n\n"
                "Verify that the 'xsdir' folder sits at the project root or beside it "
                "and contains 'endf-6/' and 'MACS/'."
            )
            # Continue: the GUI is still usable for non-data tasks.

        # Toolbar (scans `model*` folders and instantiates them)
        self._progress(0.6, "Loading models…")
        self.init_toolbar()

        # Hide all pages
        self._progress(0.95, "Preparing interface…")
        self.HidePages()

        # Apply theme again (optional)
        self.apply_theme("dark.css")

        print("The software is ready to use.")
        print(self.window_width, self.window_height)

    # ---------------- dashboard
    def init_dashboard(self):
        if self.dashboard is None:
            self.dashboard = QDashboard(self)
            self.dashboard.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)
            self.addDockWidget(Qt.RightDockWidgetArea, self.dashboard)

            self.dashboard.setFloating(False)
            self.dashboard.setMinimumWidth(300)
            self.dashboard.setMaximumWidth(self.window_width // 2)
            self.dashboard.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            # Give the panel a comfortable default width when the window first appears.
            self.resizeDocks([self.dashboard], [400], Qt.Horizontal)
            # When undocked, resize to a usable floating window.
            self.dashboard.topLevelChanged.connect(self._on_dashboard_float_changed)
        else:
            self.dashboard.Show()
        self.dashboard.hide()

    def _on_dashboard_float_changed(self, floating: bool):
        if floating:
            self.dashboard.resize(420, 680)

    # ---------------- menu
    def init_menu(self):
        menu_bar = self.menuBar()
        self._sync_branding(getattr(self, "_theme_key", "dark"))

        view_menu = menu_bar.addMenu("View")
        view_light_action = QAction("Light", self)
        view_dark_action = QAction("Dark", self)
        view_light_action.triggered.connect(lambda: self.apply_theme("light.css"))
        view_dark_action.triggered.connect(lambda: self.apply_theme("dark.css"))
        view_menu.addAction(view_light_action)
        view_menu.addAction(view_dark_action)
        view_menu.addSeparator()

        # Panels submenu — each panel can be shown/hidden here if lost.
        self._panels_menu = view_menu.addMenu("Panels")
        right_action = self.dashboard.toggleViewAction()
        right_action.setText("Right Panel")
        self._panels_menu.addAction(right_action)

        terminal_menu = menu_bar.addMenu("Terminal")
        terminal_output_action = QAction("output", self)
        terminal_output_action.triggered.connect(lambda: self.show_output())
        terminal_menu.addAction(terminal_output_action)

        extensions_menu = menu_bar.addMenu("Extensions")
        extension_rootfolder_action = QAction("In root folder", self)
        extension_rootfolder_action.triggered.connect(lambda: self.loadExtensionfromRoot())
        extensions_menu.addAction(extension_rootfolder_action)

    def add_menu_line(self):
        self.menu_line = QFrame(self)
        self.menu_line.setObjectName("menuline")
        self.menu_line.setFrameShape(QFrame.HLine)
        self.menu_line.setFrameShadow(QFrame.Sunken)
        self.menu_line.setGeometry(0, self.menuBar().height(), self.width(), 2)
        self.resizeEvent = self.update_menu_line

    def update_menu_line(self, event):
        self.menu_line.setGeometry(0, self.menuBar().height(), self.width(), 2)
        super().resizeEvent(event)

    # ---------------- actions
    def show_output(self):
        try:
            self.terminal.show()
        except AttributeError:
            pass

    # ---------------- theme
    def apply_theme(self, theme_file):
        def resource_path(filename):
            base = MODELS_DIR
            p1 = base / filename
            p2 = base / "resources" / filename
            if p1.is_file():
                return p1
            if p2.is_file():
                return p2
            return p1  # will fall through to FileNotFoundError

        try:
            real_theme_path = resource_path(theme_file)
            with open(real_theme_path, "r", encoding="windows-1252") as file:
                qapp = QApplication.instance()
                if qapp is not None:
                    qapp.setStyleSheet(file.read())
                    print(f"Stylesheet: {real_theme_path}")
                else:
                    print("Stylesheet skipped (no QApplication).")
        except FileNotFoundError:
            print(f"Fail: {theme_file} not found.")
            return

        # Broadcast to every loaded model page (e.g. the Graph panel's own
        # matplotlib colour theme) so one switch updates the whole app at once.
        theme_key = "light" if "light" in theme_file.lower() else "dark"
        self._theme_key = theme_key
        self._sync_branding(theme_key)
        get_service().set_theme(theme_key)

        # Dark theme → orange (OR) icons; light theme → green (GR) icons.
        # Icon color signals the theme palette; active state is shown via
        # CSS :checked (border/background), not by icon switching.
        self._sync_toolbar_icons(theme_key)

    def _sync_branding(self, theme_key: str):
        pixmap = QPixmap(logo_path(theme_key))
        if pixmap.isNull():
            return
        # Use the emblem area of the horizontal artwork for a square window icon.
        emblem = pixmap.copy(int(pixmap.width() * 0.07),
                             int(pixmap.height() * 0.23),
                             int(pixmap.width() * 0.25),
                             int(pixmap.height() * 0.51))
        square = QPixmap(max(emblem.width(), emblem.height()),
                         max(emblem.width(), emblem.height()))
        square.fill(Qt.transparent)
        from PyQt5.QtGui import QPainter
        painter = QPainter(square)
        painter.drawPixmap((square.width() - emblem.width()) // 2,
                           (square.height() - emblem.height()) // 2, emblem)
        painter.end()
        icon = QIcon(square)
        self.setWindowIcon(icon)
        app = QApplication.instance()
        if app is not None:
            app.setWindowIcon(icon)

    def _sync_toolbar_icons(self, theme_key: str):
        if not hasattr(self, "models_info") or not self.models_info:
            return
        for info in self.models_info:
            btn = info["button"]
            icon = btn.active_icon if theme_key == "dark" else btn.default_icon
            btn_size = btn.width()
            btn.setIcon(QIcon(icon.pixmap(btn_size, btn_size)))
            btn.setIconSize(QSize(btn_size, btn_size))

    # ---------------- models: extension loading
    def loadExtensionfromRoot(self):
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select model folder (model_Xxx)",
            "",
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks,
        )
        if not folder:
            print("Extension loading cancelled.")
            return

        print("Selected extension folder:", folder)

        model_data = self.scan_folder(folder)
        if model_data is None:
            QMessageBox.warning(
                self,
                "Invalid structure",
                "Could not load a valid model from the selected folder.\n"
                "Make sure it contains model_*.py plus two .png files (GR/OR)."
            )
            return

        self.register_model(model_data)
        QMessageBox.information(
            self,
            "Model added",
            f"Model '{model_data['name']}' has been installed."
        )
        print(f"Extension '{model_data['name']}' loaded from {folder}.")

    def scan_folder(self, folder_path):
        if not os.path.isdir(folder_path):
            print("Not a directory:", folder_path)
            return None

        folder_name = os.path.basename(folder_path)
        if not folder_name.startswith("model"):
            print("Folder does not start with 'model':", folder_name)
            return None

        model_name = folder_name[5:]
        if not model_name:
            print("Invalid model folder name (empty model_name).")
            return None

        py_path = os.path.join(folder_path, f"model_{model_name}.py")
        icon_gr_path = os.path.join(folder_path, f"FIATb-GR-{model_name}.png")
        icon_or_path = os.path.join(folder_path, f"FIATb-OR-{model_name}.png")

        for f in [py_path, icon_gr_path, icon_or_path]:
            if not os.path.isfile(f):
                print("Missing required file:", f)
                return None

        try:
            module_name = f"model_{model_name}"
            spec = importlib.util.spec_from_file_location(module_name, py_path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
        except Exception as e:
            print(f"Failed to import model {model_name}:", e)
            return None

        class_model = getattr(module, module_name, None)
        left_cls = getattr(module, "Left", None)
        right_cls = getattr(module, "Right", None)

        if not (class_model and left_cls and right_cls):
            print(f"Invalid module {model_name}: missing class(es).")
            return None
        if not (issubclass(left_cls, QWidget) and issubclass(right_cls, QVBoxLayout)):
            print(f"Invalid widget types in model {model_name}.")
            return None

        sig = inspect.signature(class_model.__init__)
        if "parent" in sig.parameters:
            model_instance = class_model(self)
        else:
            model_instance = class_model()
        if hasattr(model_instance, "set_dialog_parent"):
            model_instance.set_dialog_parent(self)

        icon_gr = QIcon(icon_gr_path)
        icon_or = QIcon(icon_or_path)

        return {
            "name": model_name,
            "model_instance": model_instance,
            "icon_gr": icon_gr,
            "icon_or": icon_or
        }

    def register_model(self, model_data):
        model_name = model_data["name"]
        model_instance = model_data["model_instance"]
        icon_gr = model_data["icon_gr"]
        icon_or = model_data["icon_or"]

        new_index = len(self.models_info)
        btn_size = int(self.window_width / 27)

        button = AnimatedToolButton(
            icon_gr,
            icon_or,
            "",
            self.window_width,
            self.window_height,
            self
        )
        button.setFixedSize(btn_size, btn_size)
        button.setToolButtonStyle(Qt.ToolButtonIconOnly)

        pixmap_gr = icon_gr.pixmap(btn_size, btn_size)
        pixmap_or = icon_or.pixmap(btn_size, btn_size)
        icon_gr_scaled = QIcon(pixmap_gr)
        icon_or_scaled = QIcon(pixmap_or)

        button.setIcon(icon_gr_scaled)
        button.setIconSize(QSize(btn_size, btn_size))
        if hasattr(button, "setHoverIcon"):
            button.setHoverIcon(icon_or_scaled)

        self.Tb.addWidget(button)
        self.buttons.append(button)

        self.TbBtnMap[button] = lambda _checked, i=new_index: self.on_toolbar_button_clicked(i)
        button.clicked.connect(lambda checked, b=button: self.TbBtnMap[b](checked))

        left_widget = model_instance.getLeftWidget()
        right_layout = model_instance.getRightLayout()
        self.dashboard.AddLayout(new_index, right_layout)
        self.AddMainWidget(new_index, left_widget)

        self.models_info.append({
            "index": new_index,
            "name": model_name,
            "model_instance": model_instance,
            "button": button
        })

        print(f"Model '{model_name}' registered.")

    # ---------------- toolbar & model scan
    def init_toolbar(self):
        if hasattr(self, "Tb") and self.Tb:
            self.removeToolBar(self.Tb)

        self.Tb = QToolBar("TOOLS")
        self.Tb.setObjectName("left_nav_toolbar")
        self.addToolBar(Qt.LeftToolBarArea, self.Tb)
        self.Tb.setAllowedAreas(Qt.LeftToolBarArea | Qt.RightToolBarArea)

        self.Tb.setMinimumWidth(int(self.window_width / 27))
        self.Tb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.Tb.setMovable(True)

        # Allow the toolbar to be restored from Panels menu if accidentally closed.
        if hasattr(self, "_panels_menu"):
            tb_action = self.Tb.toggleViewAction()
            tb_action.setText("Left Toolbar")
            self._panels_menu.addAction(tb_action)

        found_models = self.Scan_models()

        self.models_info = []
        self.MainWidgetsMap = {}
        self.buttons = []
        self.TbBtnMap = {}

        btn_size = int(self.window_width / 27)

        for idx, model_data in enumerate(found_models):
            model_name = model_data["name"]
            model_instance = model_data["model_instance"]

            button = AnimatedToolButton(
                model_data["icon_gr"],
                model_data["icon_or"],
                "",
                self.window_width,
                self.window_height,
                self
            )
            button.setFixedSize(btn_size, btn_size)
            button.setToolButtonStyle(Qt.ToolButtonIconOnly)

            pixmap_gr = model_data["icon_gr"].pixmap(btn_size, btn_size)
            pixmap_or = model_data["icon_or"].pixmap(btn_size, btn_size)
            icon_gr_scaled = QIcon(pixmap_gr)
            icon_or_scaled = QIcon(pixmap_or)

            button.setIcon(icon_gr_scaled)
            button.setIconSize(QSize(btn_size, btn_size))
            if hasattr(button, "setHoverIcon"):
                button.setHoverIcon(icon_or_scaled)

            self.Tb.addWidget(button)
            self.buttons.append(button)

            self.TbBtnMap[button] = lambda _checked, i=idx: self.on_toolbar_button_clicked(i)

            left_widget = model_instance.getLeftWidget()
            right_layout = model_instance.getRightLayout()
            self.dashboard.AddLayout(idx, right_layout)
            self.AddMainWidget(idx, left_widget)

            print(f"Model {model_name} has been installed.")

            self.models_info.append({
                "index": idx,
                "name": model_name,
                "model_instance": model_instance,
                "button": button
            })

        for button in self.buttons:
            button.clicked.connect(lambda checked, b=button: self.TbBtnMap[b](checked))

    def on_toolbar_button_clicked(self, index):
        for info in self.models_info:
            info["button"].setActive(info["index"] == index)
        self.SetPage(index)

    def Scan_models(self, base_path=None):
        # Models live next to this launcher (…/NuMatRx/model/model*).
        base_path = str(MODELS_DIR)
        print("Scanning models in:", base_path)

        # Candidate model folders (so splash progress reflects REAL load stages).
        # Left-nav order is fixed by spec — CoreNucleo, Isotope Chart, Result,
        # Graph — so order by model name explicitly rather than relying on the
        # (case-sensitive) alphabetical directory listing. Button order and the
        # stacked-widget page indices both follow this list (indices are assigned
        # sequentially below). Unknown folders (e.g. loaded extensions) keep a
        # stable alphabetical order after the known ones.
        _NAV_ORDER = ["CoreNucleo", "IsotopeChart", "results", "graph"]
        _order_rank = {name: i for i, name in enumerate(_NAV_ORDER)}
        entries = [e for e in os.listdir(base_path)
                   if e.startswith("model")
                   and os.path.isdir(os.path.join(base_path, e))]
        entries.sort(key=lambda e: (_order_rank.get(e[5:], len(_NAV_ORDER)), e))
        total = max(1, len(entries))

        found_models = []
        for i, entry in enumerate(entries):
            self._progress(0.6 + 0.32 * (i / total), f"Loading model: {entry[5:] or entry}…")
            folder_path = os.path.join(base_path, entry)

            model_name = entry[5:]
            py_path = os.path.join(folder_path, f"model_{model_name}.py")
            icon_gr_path = os.path.join(folder_path, f"FIATb-GR-{model_name}.png")
            icon_or_path = os.path.join(folder_path, f"FIATb-OR-{model_name}.png")

            if not (os.path.isfile(py_path) and os.path.isfile(icon_gr_path) and os.path.isfile(icon_or_path)):
                continue

            module_name = f"model_{model_name}"
            spec = importlib.util.spec_from_file_location(module_name, py_path)
            if spec is None:
                continue
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            try:
                spec.loader.exec_module(module)
            except Exception as e:
                print(f"Failed to import model {model_name}: {e}")
                continue

            class_model = getattr(module, module_name, None)
            left_cls = getattr(module, "Left", None)
            right_cls = getattr(module, "Right", None)
            if not (class_model and left_cls and right_cls):
                continue
            if not (issubclass(left_cls, QWidget) and issubclass(right_cls, QVBoxLayout)):
                continue

            sig = inspect.signature(class_model.__init__)
            try:
                if "parent" in sig.parameters:
                    model_instance = class_model(self)
                else:
                    model_instance = class_model()
            except Exception as e:
                print(f"Failed to instantiate model {model_name}: {e}")
                continue
            if hasattr(model_instance, "set_dialog_parent"):
                model_instance.set_dialog_parent(self)

            icon_gr = QIcon(icon_gr_path)
            icon_or = QIcon(icon_or_path)

            print(f"Model {model_name} has been found.")
            found_models.append({
                "name": model_name,
                "model_instance": model_instance,
                "icon_gr": icon_gr,
                "icon_or": icon_or
            })

        return found_models

    # ---------------- page control
    def SetPage(self, index):
        self.SetMainWidget(index)
        self.dashboard.SetLayout(index)

    def AddMainWidget(self, num, w):
        if num not in self.MainWidgetsMap:
            self.MainWidgetsMap[num] = w
            self.MainLayout.addWidget(w)

    def SetMainWidget(self, num):
        if num in self.MainWidgetsMap:
            for widget in self.MainWidgetsMap.values():
                widget.hide()
            self.MainWidgetsMap[num].show()

    def HidePages(self):
        self.dashboard.hide()
        for widget in self.MainWidgetsMap.values():
            widget.hide()

    def Qprint(self, message):
        print("FROM modules...")
        print(message)


# -------------------------
# __main__
# -------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Application icon (also used by the OS task bar).
    _icon = app_icon_path()
    if _icon:
        app.setWindowIcon(QIcon(_icon))

    # Splash first; its progress bar is driven by the real init stages.
    splash = SplashScreen(splash_image_path())
    splash.show()
    app.processEvents()

    win = MainWindow(splash=splash)
    splash.finish_with(win)
    win.show()
    sys.exit(app.exec_())
