#!/usr/bin/python3

import sys
import os
import json
import re
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QListWidget, QListWidgetItem, QTreeWidget, 
                             QTreeWidgetItem, QPushButton, QFileDialog, QTabWidget, 
                             QLineEdit, QLabel, QSplitter, QMessageBox, QInputDialog, 
                             QComboBox, QTableWidget, QTableWidgetItem, QHeaderView, 
                             QDialog, QDialogButtonBox, QMenu, QAbstractItemView, 
                             QFormLayout)
from PyQt6.QtCore import Qt, pyqtSignal, QObject
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

class ExportDialog(QDialog):
    def __init__(self, signals, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export Signals to CSV")
        self.resize(300, 400)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select signals to include in CSV:"))
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
        self.list_widget.addItems(sorted(signals.keys()))
        layout.addWidget(self.list_widget)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def get_selected(self):
        return [item.text() for item in self.list_widget.selectedItems()]

class FormulaDialog(QDialog):
    def __init__(self, signal_manager, current_expr="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Formula Editor")
        self.resize(400, 350)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Double-click a signal to insert it:"))
        self.sig_list = QListWidget()
        for name in sorted(signal_manager.signals.keys()):
            item = QListWidgetItem(signal_manager.display_name(name))
            item.setData(Qt.ItemDataRole.UserRole, name)
            self.sig_list.addItem(item)
        layout.addWidget(self.sig_list)
        self.expr_edit = QLineEdit(current_expr)
        layout.addWidget(self.expr_edit)
        self.sig_list.itemDoubleClicked.connect(
            lambda item: self.expr_edit.insert(item.data(Qt.ItemDataRole.UserRole)))
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

class SignalManager(QObject):
    data_changed = pyqtSignal()
    def __init__(self):
        super().__init__()
        self.signals = {}
        self.expressions = {}
        self.files = {}
        self.source_map = {}

    def add_imported_signal(self, name, data, filename, varname, comp_idx=None, comp_label=None):
        self.signals[name] = data
        self.source_map[name] = {"file": filename, "var": varname, "idx": comp_idx, "lbl": comp_label}
        self.data_changed.emit()

    def add_expression_signal(self, name, data, expression):
        self.signals[name] = data
        self.expressions[name] = expression
        self.data_changed.emit()

    def delete_signal(self, name):
        if name in self.signals:
            self.signals.pop(name)
            self.expressions.pop(name, None)
            self.source_map.pop(name, None)
            self.data_changed.emit()

    def rename_signal(self, old_name, new_name):
        if old_name == new_name or not new_name:
            return
        new_name = re.sub(r'[^a-zA-Z0-9_]', '_', new_name)
        if new_name in self.signals:
            raise ValueError(f"Name '{new_name}' already exists.")
        
        self.signals[new_name] = self.signals.pop(old_name)
        if old_name in self.expressions:
            self.expressions[new_name] = self.expressions.pop(old_name)
        if old_name in self.source_map:
            self.source_map[new_name] = self.source_map.pop(old_name)
        
        # Update other formulas that might use this signal
        for name, expr in self.expressions.items():
            self.expressions[name] = re.sub(rf'\b{old_name}\b', new_name, expr)
        self.data_changed.emit()

    def evaluate_expression(self, name, expr):
        try:
            context = {**np.__dict__, **self.signals}
            result = eval(expr, {"__builtins__": None}, context)
            self.add_expression_signal(name, result, expr)
            return True
        except Exception as e:
            return str(e)

    def display_name(self, name):
        """Return 'name  [source]' for UI display."""
        src = self.source_map.get(name)
        if src:
            return f"{name}  [{src['file']}]"
        elif name in self.expressions:
            return f"{name}  [formula]"
        return name

class PlotSettingsDialog(QDialog):
    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Plot Settings")
        layout = QVBoxLayout(self)
        form = QWidget()
        fl = QFormLayout(form)
        self.title_edit  = QLineEdit(config.get("title",  ""))
        self.xlabel_edit = QLineEdit(config.get("xlabel", ""))
        self.ylabel_edit = QLineEdit(config.get("ylabel", ""))
        self.xmin_edit   = QLineEdit(config.get("xmin",  ""))
        self.xmax_edit   = QLineEdit(config.get("xmax",  ""))
        self.ymin_edit   = QLineEdit(config.get("ymin",  ""))
        self.ymax_edit   = QLineEdit(config.get("ymax",  ""))
        fl.addRow("Title:",   self.title_edit)
        fl.addRow("X Label:", self.xlabel_edit)
        fl.addRow("Y Label:", self.ylabel_edit)
        fl.addRow("X Min (blank = auto):", self.xmin_edit)
        fl.addRow("X Max (blank = auto):", self.xmax_edit)
        fl.addRow("Y Min (blank = auto):", self.ymin_edit)
        fl.addRow("Y Max (blank = auto):", self.ymax_edit)
        layout.addWidget(form)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def get_config(self):
        return {
            "title":  self.title_edit.text(),
            "xlabel": self.xlabel_edit.text(),
            "ylabel": self.ylabel_edit.text(),
            "xmin":   self.xmin_edit.text(),
            "xmax":   self.xmax_edit.text(),
            "ymin":   self.ymin_edit.text(),
            "ymax":   self.ymax_edit.text(),
        }

class LoadSessionFileDialog(QDialog):
    """Dialog shown when loading a session to confirm or remap .nc file paths."""
    def __init__(self, saved_files, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Load Session – File Mapping")
        self.resize(700, 180 + 40 * max(len(saved_files), 1))
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "The session references the following NetCDF files.\n"
            "Keep the saved paths or browse to select replacements."))
        self.table = QTableWidget(len(saved_files), 3)
        self.table.setHorizontalHeaderLabels(["File", "Path", ""])
        hdr = self.table.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.file_map = {}
        for row, (basename, filepath) in enumerate(saved_files.items()):
            self.file_map[basename] = filepath
            ni = QTableWidgetItem(basename)
            ni.setFlags(ni.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 0, ni)
            pi = QTableWidgetItem(filepath)
            pi.setFlags(pi.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, 1, pi)
            btn = QPushButton("Browse…")
            btn.clicked.connect(lambda checked, r=row, bn=basename: self._browse(r, bn))
            self.table.setCellWidget(row, 2, btn)
        layout.addWidget(self.table)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def _browse(self, row, basename):
        path, _ = QFileDialog.getOpenFileName(
            self, f"Replace '{basename}'", "", "NetCDF (*.nc)")
        if path:
            self.file_map[basename] = path
            self.table.item(row, 1).setText(path)

    def get_file_map(self):
        return dict(self.file_map)

class MBDynPlotCanvas(FigureCanvasQTAgg):
    def __init__(self):
        self.fig = Figure(tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)

class MBDynToolbar(NavigationToolbar2QT):
    """Custom toolbar that captures line styles after the figure options dialog."""
    def __init__(self, canvas, parent, plot_tab):
        super().__init__(canvas, parent)
        self.plot_tab = plot_tab

    def edit_parameters(self):
        super().edit_parameters()
        self.plot_tab._capture_line_styles()
        self.plot_tab._capture_plot_config()

    def home(self, *args):
        # Clear forced axis limits so autoscaling works
        cfg = self.plot_tab.plot_config
        cfg["xmin"] = ""
        cfg["xmax"] = ""
        cfg["ymin"] = ""
        cfg["ymax"] = ""
        self.plot_tab.update_plot()

class PlotTab(QWidget):
    def __init__(self, signal_manager):
        super().__init__()
        self.mgr = signal_manager
        self.plot_config = {"title": "", "xlabel": "", "ylabel": "",
                            "xmin": "", "xmax": "", "ymin": "", "ymax": ""}
        self.curve_styles = {}   # row_index -> dict of Line2D properties
        self._active_rows = []   # order of row indices that produced ax.lines
        layout = QVBoxLayout(self)
        self.splitter = QSplitter(Qt.Orientation.Vertical)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Active", "Label", "X Signal", "Y Signal"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._table_context_menu)
        self.table.itemChanged.connect(self._on_item_changed)

        btn_row = QWidget()
        btn_row_lay = QHBoxLayout(btn_row)
        btn_row_lay.setContentsMargins(0, 0, 0, 0)
        btn_add = QPushButton("+ Add Curve")
        btn_add.clicked.connect(lambda: self.add_curve_row())
        btn_del = QPushButton("- Delete Selected")
        btn_del.clicked.connect(self.delete_selected_curves)
        btn_row_lay.addWidget(btn_add)
        btn_row_lay.addWidget(btn_del)

        self.canvas = MBDynPlotCanvas()
        self.toolbar = MBDynToolbar(self.canvas, self, plot_tab=self)
        top_ui = QWidget()
        top_lay = QVBoxLayout(top_ui)
        top_lay.addWidget(self.table)
        top_lay.addWidget(btn_row)
        self.splitter.addWidget(top_ui)
        self.splitter.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.splitter)
        btn_refresh = QPushButton("Refresh Plot")
        btn_refresh.clicked.connect(self.update_plot)
        btn_settings = QPushButton("⚙ Plot Settings")
        btn_settings.clicked.connect(self.open_plot_settings)
        bottom_row = QWidget()
        bottom_lay = QHBoxLayout(bottom_row)
        bottom_lay.setContentsMargins(0, 0, 0, 0)
        bottom_lay.addWidget(btn_refresh)
        bottom_lay.addWidget(btn_settings)
        layout.addWidget(bottom_row)

    def _populate_signal_combo(self, cb, include_auto=False):
        """Fill a combobox with signal display names, storing real names as UserRole."""
        cb.clear()
        if include_auto:
            cb.addItem("Index (Auto)", "Index (Auto)")
        for name in sorted(self.mgr.signals.keys()):
            cb.addItem(self.mgr.display_name(name), name)

    def _combo_signal_name(self, cb):
        """Get the real signal name from a combobox (stored in UserRole)."""
        data = cb.currentData(Qt.ItemDataRole.UserRole)
        return data if data else cb.currentText()

    def _set_combo_by_name(self, cb, sig_name):
        """Select a combobox item by its real signal name."""
        for i in range(cb.count()):
            if cb.itemData(i, Qt.ItemDataRole.UserRole) == sig_name:
                cb.setCurrentIndex(i)
                return
        cb.setCurrentText(sig_name)

    def add_curve_row(self, x_sig="Index (Auto)", y_sig="", label="", active=True):
        row = self.table.rowCount()
        self.table.blockSignals(True)
        self.table.insertRow(row)

        chk = QTableWidgetItem()
        chk.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
        chk.setCheckState(Qt.CheckState.Checked if active else Qt.CheckState.Unchecked)
        self.table.setItem(row, 0, chk)
        self.table.setItem(row, 1, QTableWidgetItem(label if label else f"Curve {row+1}"))

        x_cb = QComboBox()
        y_cb = QComboBox()
        self._populate_signal_combo(x_cb, include_auto=True)
        self._populate_signal_combo(y_cb)

        if x_sig:
            self._set_combo_by_name(x_cb, x_sig)
        if y_sig:
            self._set_combo_by_name(y_cb, y_sig)
        self.table.setCellWidget(row, 2, x_cb)
        self.table.setCellWidget(row, 3, y_cb)
        self.table.blockSignals(False)

        # Connect AFTER setting values to avoid spurious update_plot calls
        x_cb.currentIndexChanged.connect(lambda: self.update_plot())
        y_cb.currentIndexChanged.connect(lambda: self.update_plot())

    def _line_to_style(self, line):
        import matplotlib.colors as mcolors
        alpha = line.get_alpha()
        def _hex(c):
            try:
                return str(mcolors.to_hex(c, keep_alpha=False))
            except Exception:
                return None
        return {
            "color":            _hex(line.get_color()) or "#000000",
            "linestyle":        str(line.get_linestyle()),
            "linewidth":        float(line.get_linewidth()),
            "marker":           str(line.get_marker()),
            "markersize":       float(line.get_markersize()),
            "markeredgecolor":  _hex(line.get_markeredgecolor()),
            "markerfacecolor":  _hex(line.get_markerfacecolor()),
            "markeredgewidth":  float(line.get_markeredgewidth()),
            "alpha":            float(alpha) if alpha is not None else 1.0,
        }

    def _capture_line_styles(self):
        """Read current line styles from ax.lines using GID → row mapping.
        Each row may produce multiple lines (multi-component signals),
        so curve_styles[r] is a list of style dicts."""
        collected = {}  # r -> [style, style, ...]
        for line in self.canvas.ax.lines:
            gid = line.get_gid()
            if gid is None:
                continue
            parts = gid.split(".")
            if parts[0].isdigit():
                r = int(parts[0])
                if r < self.table.rowCount():
                    collected.setdefault(r, []).append(self._line_to_style(line))
        for r, styles in collected.items():
            self.curve_styles[r] = styles

    def _capture_plot_config(self):
        """Read current title, labels, and axis limits from the axes."""
        ax = self.canvas.ax
        self.plot_config["title"]  = ax.get_title()
        self.plot_config["xlabel"] = ax.get_xlabel()
        self.plot_config["ylabel"] = ax.get_ylabel()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        self.plot_config["xmin"] = str(xlim[0])
        self.plot_config["xmax"] = str(xlim[1])
        self.plot_config["ymin"] = str(ylim[0])
        self.plot_config["ymax"] = str(ylim[1])

    def _apply_line_style(self, line, style):
        if not style:
            return
        if "color"           in style: line.set_color(style["color"])
        if "linestyle"       in style: line.set_linestyle(style["linestyle"])
        if "linewidth"       in style: line.set_linewidth(style["linewidth"])
        if "marker"          in style: line.set_marker(style["marker"])
        if "markersize"      in style: line.set_markersize(style["markersize"])
        if "markeredgecolor" in style and style["markeredgecolor"]:
            line.set_markeredgecolor(style["markeredgecolor"])
        if "markerfacecolor" in style and style["markerfacecolor"]:
            line.set_markerfacecolor(style["markerfacecolor"])
        if "markeredgewidth" in style:
            line.set_markeredgewidth(style["markeredgewidth"])
        if "alpha"           in style: line.set_alpha(style["alpha"])

    def open_plot_settings(self):
        dlg = PlotSettingsDialog(self.plot_config, self)
        if dlg.exec():
            self.plot_config = dlg.get_config()
            self.update_plot()

    def _apply_plot_config(self):
        cfg = self.plot_config
        ax = self.canvas.ax
        if cfg.get("title"):
            ax.set_title(cfg["title"])
        if cfg.get("xlabel"):
            ax.set_xlabel(cfg["xlabel"])
        if cfg.get("ylabel"):
            ax.set_ylabel(cfg["ylabel"])
        def _f(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
        xmin, xmax = _f(cfg.get("xmin")), _f(cfg.get("xmax"))
        ymin, ymax = _f(cfg.get("ymin")), _f(cfg.get("ymax"))
        if xmin is not None or xmax is not None:
            cur_xl = ax.get_xlim()
            ax.set_xlim(xmin if xmin is not None else cur_xl[0],
                        xmax if xmax is not None else cur_xl[1])
        if ymin is not None or ymax is not None:
            cur_yl = ax.get_ylim()
            ax.set_ylim(ymin if ymin is not None else cur_yl[0],
                        ymax if ymax is not None else cur_yl[1])

    def _on_item_changed(self, item):
        if item.column() == 0:
            self.update_plot()

    def _table_context_menu(self, pos):
        row = self.table.rowAt(pos.y())
        if row < 0:
            return
        menu = QMenu()
        del_act = menu.addAction("Delete Curve")
        toggle_act = menu.addAction("Toggle Active")
        act = menu.exec(self.table.mapToGlobal(pos))
        if act == del_act:
            self._remove_style_for_row(row)
            self.table.removeRow(row)
            self.update_plot()
        elif act == toggle_act:
            chk = self.table.item(row, 0)
            if chk:
                new_state = (Qt.CheckState.Unchecked
                             if chk.checkState() == Qt.CheckState.Checked
                             else Qt.CheckState.Checked)
                chk.setCheckState(new_state)

    def _remove_style_for_row(self, row):
        self.curve_styles.pop(row, None)
        # shift keys above the deleted row down by 1
        self.curve_styles = {
            (k - 1 if k > row else k): v
            for k, v in self.curve_styles.items()
        }

    def delete_selected_curves(self):
        rows = sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            self._remove_style_for_row(r)
            self.table.removeRow(r)
        if rows:
            self.update_plot()

    def update_plot(self):
        self.canvas.ax.clear()
        self._active_rows = []
        for r in range(self.table.rowCount()):
            chk = self.table.item(r, 0)
            if chk and chk.checkState() != Qt.CheckState.Checked:
                continue
            lbl_item = self.table.item(r, 1)
            lbl = lbl_item.text() if lbl_item else f"Curve {r+1}"
            x_widget = self.table.cellWidget(r, 2)
            y_widget = self.table.cellWidget(r, 3)
            if x_widget is None or y_widget is None:
                continue
            x_n = self._combo_signal_name(x_widget)
            y_n = self._combo_signal_name(y_widget)
            y_d = self.mgr.signals.get(y_n)
            if y_d is not None:
                x_d = self.mgr.signals.get(x_n) if x_n != "Index (Auto)" else None
                x_vals = x_d if x_d is not None else range(len(y_d))
                lines = self.canvas.ax.plot(x_vals, y_d)
                comp_labels = ["X", "Y", "Z"] if len(lines) == 3 else \
                              [f"C{j}" for j in range(len(lines))]
                saved = self.curve_styles.get(r, [])
                if isinstance(saved, dict):
                    saved = [saved]
                for j, line in enumerate(lines):
                    line.set_gid(f"{r}.{j}")
                    if len(lines) == 1:
                        line.set_label(lbl)
                    else:
                        line.set_label(f"{lbl} [{comp_labels[j]}]")
                    if j < len(saved):
                        self._apply_line_style(line, saved[j])
                self._active_rows.append(r)
        
        # Only update legend if we have labeled artists
        if self.canvas.ax.lines:
            self.canvas.ax.legend()
        
        self.canvas.ax.grid(True, alpha=0.3)
        self._apply_plot_config()
        self.canvas.draw()

class PostProcessor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MBDyn Post-Processor")
        self.resize(1400, 900)
        self.mgr = SignalManager()
        self.mgr.data_changed.connect(self.sync_ui)
        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        side = QWidget()
        s_lay = QVBoxLayout(side)
        
        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("NetCDF Variables")
        self.tree.itemDoubleClicked.connect(self.import_signal)
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._tree_context_menu)
        
        self.var_search = QLineEdit()
        self.var_search.setPlaceholderText("Filter variables…")
        self.var_search.textChanged.connect(self._filter_tree)
        
        self.sig_list = QListWidget()
        self.sig_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.sig_list.customContextMenuRequested.connect(self.show_signal_menu)
        self.sig_list.itemDoubleClicked.connect(self.edit_signal)
        
        btn_open = QPushButton("Open NetCDF")
        btn_open.clicked.connect(self.open_file)
        btn_refresh_files = QPushButton("Refresh Files Data")
        btn_refresh_files.clicked.connect(self.refresh_all_data)
        btn_calc = QPushButton("New Formula")
        btn_calc.clicked.connect(self.create_expression)
        btn_export = QPushButton("Export CSV")
        btn_export.clicked.connect(self.export_to_csv)
        btn_save = QPushButton("Save Session")
        btn_save.clicked.connect(self.save_session)
        btn_load = QPushButton("Load Session")
        btn_load.clicked.connect(self.load_session)
        
        for w in [btn_open, btn_refresh_files, QLabel("<b>Variables</b>"), self.var_search, self.tree, QLabel("<b>Signals (Dbl-Click to Edit)</b>"), 
                  self.sig_list, btn_calc, btn_export, btn_save, btn_load]:
            s_lay.addWidget(w)
            
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(lambda i: self.tabs.removeTab(i))
        self.tabs.tabBarDoubleClicked.connect(self._rename_tab)
        
        main_box = QWidget()
        m_lay = QVBoxLayout(main_box)
        btn_tab = QPushButton("+ New Plot")
        btn_tab.clicked.connect(lambda: self.add_plot_tab())
        m_lay.addWidget(btn_tab)
        m_lay.addWidget(self.tabs)
        
        splitter.addWidget(side)
        splitter.addWidget(main_box)
        splitter.setSizes([467, 933])  # ~1/3 and ~2/3 of 1400
        layout.addWidget(splitter)

    def _filter_tree(self, text):
        text = text.lower()
        for i in range(self.tree.topLevelItemCount()):
            root = self.tree.topLevelItem(i)
            any_visible = False
            for j in range(root.childCount()):
                var_item = root.child(j)
                match = text in var_item.text(0).lower() if text else True
                var_item.setHidden(not match)
                if match:
                    any_visible = True
            root.setHidden(not any_visible and bool(text))

    def show_signal_menu(self, pos):
        item = self.sig_list.itemAt(pos)
        if item:
            name = item.data(Qt.ItemDataRole.UserRole) or item.text()
            menu = QMenu()
            del_act = menu.addAction("Delete Signal")
            ren_act = menu.addAction("Rename Signal")
            act = menu.exec(self.sig_list.mapToGlobal(pos))
            if act == del_act:
                self._delete_signal_with_deps(name)
            elif act == ren_act:
                self.rename_signal_prompt(name)

    def _delete_signal_with_deps(self, name):
        """Delete a signal, warning about and deactivating dependent expressions/curves."""
        missing = {name}

        # Find dependent expressions transitively
        dep_exprs = set()
        changed = True
        while changed:
            changed = False
            for ename, expr_str in self.mgr.expressions.items():
                if ename in missing:
                    continue
                for sig in list(missing):
                    if re.search(rf'\b{re.escape(sig)}\b', expr_str):
                        missing.add(ename)
                        dep_exprs.add(ename)
                        changed = True
                        break

        # Find affected curves
        affected_curves = []
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            for r in range(tab.table.rowCount()):
                x_w = tab.table.cellWidget(r, 2)
                y_w = tab.table.cellWidget(r, 3)
                if x_w and y_w:
                    x_n = tab._combo_signal_name(x_w)
                    y_n = tab._combo_signal_name(y_w)
                    if x_n in missing or y_n in missing:
                        lbl = tab.table.item(r, 1)
                        affected_curves.append(
                            (i, r, self.tabs.tabText(i),
                             lbl.text() if lbl else f"Curve {r+1}"))

        # Build warning if anything is affected
        if dep_exprs or affected_curves:
            msg = f"Deleting signal '{name}' will also affect:\n"
            if dep_exprs:
                msg += "\nExpressions (will be deleted):\n"
                msg += "\n".join(f"  • {e}" for e in sorted(dep_exprs))
            if affected_curves:
                msg += "\nCurves (will be deleted):\n"
                msg += "\n".join(f"  • '{c[2]}' → '{c[3]}'" for c in affected_curves)
            msg += "\n\nProceed?"
            reply = QMessageBox.warning(
                self, "Dependent Items", msg,
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply != QMessageBox.StandardButton.Yes:
                return

        # Delete affected curves (reverse order to keep row indices valid)
        for tab_idx, row, _, _ in sorted(affected_curves, key=lambda c: c[1], reverse=True):
            tab = self.tabs.widget(tab_idx)
            tab._remove_style_for_row(row)
            tab.table.removeRow(row)

        # Delete dependent expressions first, then the signal itself
        for ename in dep_exprs:
            self.mgr.delete_signal(ename)
        self.mgr.delete_signal(name)

    def rename_signal_prompt(self, old_name):
        new_name, ok = QInputDialog.getText(self, "Rename", "New Name:", text=old_name)
        if ok and new_name:
            try:
                self.mgr.rename_signal(old_name, new_name)
            except ValueError as e:
                QMessageBox.warning(self, "Error", str(e))

    def edit_signal(self, item):
        name = item.data(Qt.ItemDataRole.UserRole) or item.text()
        if name in self.mgr.expressions:
            # Edit Formula
            dlg = FormulaDialog(self.mgr, current_expr=self.mgr.expressions[name], parent=self)
            if dlg.exec():
                res = self.mgr.evaluate_expression(name, dlg.expr_edit.text())
                if res is not True:
                    QMessageBox.critical(self, "Formula Error", res)
        else:
            # For imported signals, just offer rename
            self.rename_signal_prompt(name)

    def open_file(self, path=None):
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Open NetCDF", "", "NetCDF (*.nc)")
        if path:
            fname = os.path.basename(path)
            self.mgr.files[fname] = path
            with Dataset(path) as ds:
                root = QTreeWidgetItem([fname])
                for v_name in ds.variables:
                    v = ds.variables[v_name]
                    var_item = QTreeWidgetItem([v_name])
                    root.addChild(var_item)
                    if len(v.dimensions) > 1:
                        sz = ds.dimensions[v.dimensions[1]].size
                        for i in range(sz):
                            lbl = ["X", "Y", "Z"][i] if sz == 3 else f"C{i}"
                            ci = QTreeWidgetItem([lbl])
                            ci.setData(0, Qt.ItemDataRole.UserRole, i)
                            var_item.addChild(ci)
                self.tree.addTopLevelItem(root)
                root.setExpanded(True)

    def refresh_all_data(self):
        if not self.mgr.files:
            return
        for name, meta in list(self.mgr.source_map.items()):
            filepath = self.mgr.files.get(meta["file"])
            if filepath and os.path.exists(filepath):
                with Dataset(filepath) as ds:
                    vn = meta["var"]
                    idx = meta["idx"]
                    if vn in ds.variables:
                        data = ds.variables[vn][:].copy()
                        if idx is not None:
                            data = data[:, idx]
                        self.mgr.signals[name] = data

        for name, expr in self.mgr.expressions.items():
            self.mgr.evaluate_expression(name, expr)
        
        self.mgr.data_changed.emit()
        for i in range(self.tabs.count()):
            self.tabs.widget(i).update_plot()
        QMessageBox.information(self, "Refresh", "All signals reloaded from disk.")

    def import_signal(self, item):
        parent = item.parent()
        if not parent:
            return
        idx = item.data(0, Qt.ItemDataRole.UserRole)
        if parent.parent():
            fn = parent.parent().text(0)
            vn = parent.text(0)
        else:
            fn = parent.text(0)
            vn = item.text(0)
            
        with Dataset(self.mgr.files[fn]) as ds:
            data = ds.variables[vn][:].copy()
            if idx is not None:
                data = data[:, idx]
        
        name = f"{vn.replace('.', '_')}_{item.text(0)}_{len(self.mgr.signals)}"
        self.mgr.add_imported_signal(name, data, fn, vn, idx, item.text(0))

    def create_expression(self):
        dlg = FormulaDialog(self.mgr, parent=self)
        if dlg.exec():
            expr = dlg.expr_edit.text()
            name = f"calc_{len(self.mgr.signals)}"
            res = self.mgr.evaluate_expression(name, expr)
            if res is not True:
                QMessageBox.critical(self, "Error", res)

    def export_to_csv(self):
        dlg = ExportDialog(self.mgr.signals, self)
        if dlg.exec():
            sel = dlg.get_selected()
            if not sel:
                return
            path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV (*.csv)")
            if path:
                pd.DataFrame({s: self.mgr.signals[s] for s in sel}).to_csv(path, index=False)

    def save_session(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Session", "", "JSON (*.json)")
        if not path:
            return
        session = {
            "files": self.mgr.files, 
            "imports": self.mgr.source_map, 
            "exprs": self.mgr.expressions, 
            "plots": []
        }
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            tab._capture_line_styles()  # snapshot current ax state before serialising
            tab._capture_plot_config()
            curves = []
            for r in range(tab.table.rowCount()):
                chk = tab.table.item(r, 0)
                lbl_item = tab.table.item(r, 1)
                curves.append({
                    "lbl": lbl_item.text() if lbl_item else f"Curve {r+1}",
                    "x": tab._combo_signal_name(tab.table.cellWidget(r, 2)),
                    "y": tab._combo_signal_name(tab.table.cellWidget(r, 3)),
                    "active": chk.checkState() == Qt.CheckState.Checked if chk else True,
                    "style": tab.curve_styles.get(r, {})
                })
            session["plots"].append({
                "title": self.tabs.tabText(i),
                "curves": curves,
                "plot_config": tab.plot_config
            })
        with open(path, 'w') as f:
            json.dump(session, f)

    # ---- helpers for file compatibility checks ----

    def _check_file_vars(self, filepath, imports, basename):
        """Return set of signal names whose variables are missing from *filepath*."""
        missing = set()
        with Dataset(filepath) as ds:
            for sig_name, meta in imports.items():
                if meta["file"] != basename:
                    continue
                vn = meta["var"]
                if vn not in ds.variables:
                    missing.add(sig_name)
                elif meta["idx"] is not None:
                    v = ds.variables[vn]
                    if (len(v.dimensions) <= 1
                            or meta["idx"] >= ds.dimensions[v.dimensions[1]].size):
                        missing.add(sig_name)
        return missing

    def _find_dependent_exprs(self, missing_signals, exprs):
        """Return set of expression names that transitively depend on *missing_signals*."""
        all_missing = set(missing_signals)
        changed = True
        while changed:
            changed = False
            for name, expr_str in exprs.items():
                if name in all_missing:
                    continue
                for sig in list(all_missing):
                    if re.search(rf'\b{re.escape(sig)}\b', expr_str):
                        all_missing.add(name)
                        changed = True
                        break
        return all_missing - missing_signals

    def _find_affected_curves(self, missing_signals, plots):
        """Return list of descriptive strings for curves that use *missing_signals*."""
        affected = []
        for plot in plots:
            for curve in plot["curves"]:
                if curve["x"] in missing_signals or curve["y"] in missing_signals:
                    affected.append(f"'{plot['title']}' → '{curve['lbl']}'")
        return affected

    # ---- right-click on tree file names ----

    def _tree_context_menu(self, pos):
        item = self.tree.itemAt(pos)
        if not item or item.parent():
            return
        menu = QMenu()
        change_act = menu.addAction("Change File…")
        delete_act = menu.addAction("Delete File")
        act = menu.exec(self.tree.mapToGlobal(pos))
        if act == change_act:
            self._change_tree_file(item)
        elif act == delete_act:
            self._delete_tree_file(item)

    def _delete_tree_file(self, tree_item):
        """Remove a loaded .nc file and all its signals, dependent expressions, and curves."""
        bn = tree_item.text(0)

        # Collect signals sourced from this file
        file_signals = {n for n, m in self.mgr.source_map.items() if m["file"] == bn}

        # Transitively find dependent expressions
        all_to_remove = set(file_signals)
        dep_exprs = set()
        changed = True
        while changed:
            changed = False
            for ename, expr_str in self.mgr.expressions.items():
                if ename in all_to_remove:
                    continue
                for sig in list(all_to_remove):
                    if re.search(rf'\b{re.escape(sig)}\b', expr_str):
                        all_to_remove.add(ename)
                        dep_exprs.add(ename)
                        changed = True
                        break

        # Find affected curves
        affected_curves = []
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            for r in range(tab.table.rowCount()):
                x_w = tab.table.cellWidget(r, 2)
                y_w = tab.table.cellWidget(r, 3)
                if x_w and y_w:
                    x_n = tab._combo_signal_name(x_w)
                    y_n = tab._combo_signal_name(y_w)
                    if x_n in all_to_remove or y_n in all_to_remove:
                        lbl = tab.table.item(r, 1)
                        affected_curves.append(
                            (i, r, self.tabs.tabText(i),
                             lbl.text() if lbl else f"Curve {r+1}"))

        # Confirm
        msg = f"Delete file '{bn}' and all its data?\n"
        if file_signals:
            msg += "\nSignals (will be deleted):\n"
            msg += "\n".join(f"  • {s}" for s in sorted(file_signals))
        if dep_exprs:
            msg += "\nExpressions (will be deleted):\n"
            msg += "\n".join(f"  • {e}" for e in sorted(dep_exprs))
        if affected_curves:
            msg += "\nCurves (will be deleted):\n"
            msg += "\n".join(f"  • '{c[2]}' → '{c[3]}'" for c in affected_curves)

        reply = QMessageBox.warning(
            self, "Delete File", msg,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if reply != QMessageBox.StandardButton.Yes:
            return

        # Delete affected curves (reverse order to keep row indices valid)
        for tab_idx, row, _, _ in sorted(affected_curves, key=lambda c: c[1], reverse=True):
            tab = self.tabs.widget(tab_idx)
            tab._remove_style_for_row(row)
            tab.table.removeRow(row)

        # Delete signals and expressions
        for sig in all_to_remove:
            self.mgr.delete_signal(sig)

        # Remove file registration
        self.mgr.files.pop(bn, None)

        # Remove tree item
        idx = self.tree.indexOfTopLevelItem(tree_item)
        if idx >= 0:
            self.tree.takeTopLevelItem(idx)

        # Refresh plots
        for i in range(self.tabs.count()):
            self.tabs.widget(i).update_plot()

    def _change_tree_file(self, tree_item):
        """Replace a loaded .nc file via right-click, checking compatibility."""
        old_bn = tree_item.text(0)
        new_path, _ = QFileDialog.getOpenFileName(
            self, f"Replace '{old_bn}'", "", "NetCDF (*.nc)")
        if not new_path or not os.path.exists(new_path):
            return
        new_bn = os.path.basename(new_path)

        # Signals sourced from the old file
        affected_imports = {n: m for n, m in self.mgr.source_map.items()
                           if m["file"] == old_bn}

        # Check which variables are missing in the new file
        missing_sigs = set()
        with Dataset(new_path) as ds:
            for sig_name, meta in affected_imports.items():
                vn = meta["var"]
                if vn not in ds.variables:
                    missing_sigs.add(sig_name)
                elif meta["idx"] is not None:
                    v = ds.variables[vn]
                    if (len(v.dimensions) <= 1
                            or meta["idx"] >= ds.dimensions[v.dimensions[1]].size):
                        missing_sigs.add(sig_name)

        if missing_sigs:
            # Transitively find dependent expressions
            missing_exprs = set()
            all_missing = set(missing_sigs)
            changed = True
            while changed:
                changed = False
                for name, expr_str in self.mgr.expressions.items():
                    if name in all_missing:
                        continue
                    for sig in list(all_missing):
                        if re.search(rf'\b{re.escape(sig)}\b', expr_str):
                            all_missing.add(name)
                            missing_exprs.add(name)
                            changed = True
                            break

            # Find affected curves
            affected_curves = []
            for i in range(self.tabs.count()):
                tab = self.tabs.widget(i)
                for r in range(tab.table.rowCount()):
                    x_w = tab.table.cellWidget(r, 2)
                    y_w = tab.table.cellWidget(r, 3)
                    if x_w and y_w:
                        x_n = tab._combo_signal_name(x_w)
                        y_n = tab._combo_signal_name(y_w)
                        if x_n in all_missing or y_n in all_missing:
                            lbl = tab.table.item(r, 1)
                            affected_curves.append(
                                (i, r, self.tabs.tabText(i),
                                 lbl.text() if lbl else f"Curve {r+1}"))

            msg = (f"The replacement file '{new_bn}' is missing variables "
                   f"needed by these signals:\n\n")
            msg += "\n".join(f"  • {s}" for s in sorted(missing_sigs))
            if missing_exprs:
                msg += "\n\nDependent expressions:\n"
                msg += "\n".join(f"  • {e}" for e in sorted(missing_exprs))
            if affected_curves:
                msg += "\n\nCurves that will be deleted:\n"
                msg += "\n".join(f"  • '{c[2]}' → '{c[3]}'" for c in affected_curves)
            msg += ("\n\nProceed (affected items will be removed),"
                    " or cancel the file change?")

            reply = QMessageBox.warning(
                self, "Missing Variables", msg,
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply != QMessageBox.StandardButton.Yes:
                return

            # Delete affected curves (reverse order to keep row indices valid)
            for tab_idx, row, _, _ in sorted(affected_curves, key=lambda c: c[1], reverse=True):
                tab = self.tabs.widget(tab_idx)
                tab._remove_style_for_row(row)
                tab.table.removeRow(row)

            # Delete missing signals and expressions
            for sig in all_missing:
                self.mgr.delete_signal(sig)

        # Update file registration
        self.mgr.files.pop(old_bn, None)
        self.mgr.files[new_bn] = new_path

        # Update source_map references
        for meta in self.mgr.source_map.values():
            if meta["file"] == old_bn:
                meta["file"] = new_bn

        # Reload data for remaining signals from new file
        with Dataset(new_path) as ds:
            for name, meta in self.mgr.source_map.items():
                if meta["file"] == new_bn:
                    data = ds.variables[meta["var"]][:].copy()
                    if meta["idx"] is not None:
                        data = data[:, meta["idx"]]
                    self.mgr.signals[name] = data

        # Re-evaluate surviving expressions
        for name, expr in list(self.mgr.expressions.items()):
            self.mgr.evaluate_expression(name, expr)

        # Rebuild tree item
        tree_item.setText(0, new_bn)
        tree_item.takeChildren()
        with Dataset(new_path) as ds:
            for v_name in ds.variables:
                v = ds.variables[v_name]
                var_item = QTreeWidgetItem([v_name])
                tree_item.addChild(var_item)
                if len(v.dimensions) > 1:
                    sz = ds.dimensions[v.dimensions[1]].size
                    for i in range(sz):
                        lbl = ["X", "Y", "Z"][i] if sz == 3 else f"C{i}"
                        ci = QTreeWidgetItem([lbl])
                        ci.setData(0, Qt.ItemDataRole.UserRole, i)
                        var_item.addChild(ci)

        self.mgr.data_changed.emit()
        for i in range(self.tabs.count()):
            self.tabs.widget(i).update_plot()

    # ---- session load / save ----

    def load_session(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load Session", "", "JSON (*.json)")
        if not path:
            return
        with open(path, 'r') as f:
            s = json.load(f)

        # Show file-mapping dialog
        dlg = LoadSessionFileDialog(s["files"], self)
        if not dlg.exec():
            return
        file_map = dlg.get_file_map()          # old_basename → new_full_path

        # Derive basename remapping
        basename_remap = {old: os.path.basename(fp) for old, fp in file_map.items()}

        # Check compatibility for every changed file
        skipped_signals = set()
        for old_bn, new_fp in file_map.items():
            if new_fp == s["files"][old_bn]:
                continue  # unchanged
            if not os.path.exists(new_fp):
                QMessageBox.warning(self, "File Not Found",
                                    f"File not found:\n{new_fp}\n\nReverting to saved path.")
                file_map[old_bn] = s["files"][old_bn]
                basename_remap[old_bn] = old_bn
                continue

            missing_imports = self._check_file_vars(new_fp, s["imports"], old_bn)
            if not missing_imports:
                continue

            missing_exprs = self._find_dependent_exprs(
                missing_imports, s.get("exprs", {}))
            all_missing = missing_imports | missing_exprs
            affected_curves = self._find_affected_curves(
                all_missing, s.get("plots", []))

            new_bn = os.path.basename(new_fp)
            msg = (f"The replacement file '{new_bn}' is missing variables "
                   f"needed by these signals:\n\n")
            msg += "\n".join(f"  • {sig}" for sig in sorted(missing_imports))
            if missing_exprs:
                msg += "\n\nDependent expressions:\n"
                msg += "\n".join(f"  • {e}" for e in sorted(missing_exprs))
            if affected_curves:
                msg += "\n\nCurves that will be deleted:\n"
                msg += "\n".join(f"  • {c}" for c in affected_curves)
            msg += ("\n\nClick Yes to proceed (affected items will be removed).\n"
                    "Click No to keep the original file.")

            reply = QMessageBox.warning(
                self, "Missing Variables", msg,
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                skipped_signals |= all_missing
            else:
                file_map[old_bn] = s["files"][old_bn]
                basename_remap[old_bn] = old_bn

        # Open files
        for new_fp in file_map.values():
            self.open_file(new_fp)

        # Import signals (skip missing ones)
        for n, m in s["imports"].items():
            if n in skipped_signals:
                continue
            new_bn = basename_remap.get(m["file"], m["file"])
            if new_bn in self.mgr.files:
                with Dataset(self.mgr.files[new_bn]) as ds:
                    data = ds.variables[m["var"]][:].copy()
                    if m["idx"] is not None:
                        data = data[:, m["idx"]]
                    self.mgr.add_imported_signal(
                        n, data, new_bn, m["var"], m["idx"], m["lbl"])

        # Evaluate expressions (skip affected ones)
        for n, e in s.get("exprs", {}).items():
            if n in skipped_signals:
                continue
            self.mgr.evaluate_expression(n, e)

        # Build plots – skip curves that reference skipped signals
        for p in s["plots"]:
            tab = self.add_plot_tab(p["title"])
            if "plot_config" in p:
                tab.plot_config = p["plot_config"]
            row = 0
            for c in p["curves"]:
                if c["x"] in skipped_signals or c["y"] in skipped_signals:
                    continue
                tab.add_curve_row(c["x"], c["y"], c["lbl"], c.get("active", True))
                if c.get("style"):
                    style = c["style"]
                    if isinstance(style, dict):
                        style = [style]
                    tab.curve_styles[row] = style
                row += 1
            tab.update_plot()


    def add_plot_tab(self, title=None):
        tab = PlotTab(self.mgr)
        self.tabs.addTab(tab, title if title else f"Plot {self.tabs.count()+1}")
        return tab

    def _rename_tab(self, index):
        if index < 0:
            return
        old = self.tabs.tabText(index)
        new_name, ok = QInputDialog.getText(self, "Rename Tab", "New name:", text=old)
        if ok and new_name:
            self.tabs.setTabText(index, new_name)

    def sync_ui(self):
        self.sig_list.clear()
        names = sorted(self.mgr.signals.keys())
        for name in names:
            item = QListWidgetItem(self.mgr.display_name(name))
            item.setData(Qt.ItemDataRole.UserRole, name)
            src = self.mgr.source_map.get(name)
            item.setToolTip(f"Source: {src['file']}" if src else "")
            self.sig_list.addItem(item)
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            for r in range(tab.table.rowCount()):
                for c_idx in [2, 3]:
                    cb = tab.table.cellWidget(r, c_idx)
                    curr = tab._combo_signal_name(cb)
                    cb.blockSignals(True)
                    tab._populate_signal_combo(cb, include_auto=(c_idx == 2))
                    tab._set_combo_by_name(cb, curr)
                    cb.blockSignals(False)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = PostProcessor()
    win.show()
    sys.exit(app.exec())

