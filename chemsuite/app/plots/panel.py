# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.
"""A plot panel: the plot, its options, selection linking and hover cards.

The options (kind, columns, bins, fits…) live in a popup, :class:`PlotOptions`, opened from the
settings button on the panel's dock title bar (see ``PanelDock``), so the plot keeps the room.

The panel follows the table: it redraws (debounced) when cells, columns or filters change, and
plots only the rows the table shows. Selecting in the plot sets the shared selection; selecting
elsewhere highlights the matching points here.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QEvent, QPoint, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QActionGroup, QColor, QCursor, QImage, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QLabel,
    QSpinBox,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.medchem_space import GUIDES
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table.columns import NUMBER, TEXT
from chemsuite.plotting.data import colorable_columns, number_columns, shown_oids
from chemsuite.plotting.stats import CURVES

from ..flow_layout import compact
from .glyphs import glyph
from .renderers import RENDERERS, Hover, PlotSpec
from .view_box import BOX, LASSO, ZOOM, SelectViewBox

NONE = "(none)"
KINDS = {"scatter": "Scatter", "histogram": "Histogram", "box": "Box Plot"}
REBUILD_DELAY_MS = 150
ImageProvider = Callable[[int], QImage | None]


class HoverCard(QFrame):
    """A tooltip-like card: structure drawing and values for the point under the mouse."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.ToolTip)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setStyleSheet(
            "HoverCard { background: palette(base); border: 1px solid palette(mid); }"
        )
        self.image = QLabel()
        self.text = QLabel()
        self.text.setTextFormat(Qt.TextFormat.PlainText)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.image, 0, Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.text)

    def show_at(self, pos: QPoint, image: QImage | None, lines: list[str]) -> None:
        self.image.setVisible(image is not None)
        if image is not None:
            self.image.setPixmap(QPixmap.fromImage(image))
        self.text.setText("\n".join(lines))
        self.adjustSize()
        self.move(pos + QPoint(16, 16))
        self.show()


class PlotOptions(QFrame):
    """The plot's options, in a popup that closes when you click elsewhere."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.Popup)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.form = QFormLayout(self)
        self.form.setContentsMargins(10, 10, 10, 10)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

    def popup(self, anchor: QWidget) -> None:
        """Show below *anchor* (the settings button), kept on its screen."""
        self.adjustSize()
        below = anchor.mapToGlobal(QPoint(anchor.width() - self.width(), anchor.height()))
        screen = anchor.screen().availableGeometry()
        x = min(max(below.x(), screen.left()), screen.right() - self.width())
        y = below.y()
        if y + self.height() > screen.bottom():  # no room below: open above
            y = anchor.mapToGlobal(QPoint(0, 0)).y() - self.height()
        self.move(x, max(y, screen.top()))
        self.show()


class PlotPanel(QWidget):
    titleChanged = Signal(str)

    def __init__(
        self,
        table: Table,
        selection: Selection,
        images: ImageProvider | None = None,
        kind: str = "scatter",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self.selection = selection
        self._images = images
        self._origin = f"plot:{id(self)}"
        self._last_axes: tuple | None = None

        self.kind = QComboBox()
        for key, title in KINDS.items():
            self.kind.addItem(title, key)
        self.kind.setCurrentIndex(max(0, self.kind.findData(kind)))  # unknown kind: scatter
        self.x = QComboBox()
        self.y = QComboBox()
        self.color = QComboBox()
        self.size = QComboBox()
        self.bins = QSpinBox()
        self.bins.setRange(0, 200)
        self.bins.setSpecialValueText("Auto")
        self.bins.setToolTip("Number of bins (Auto picks a width from the data)")
        self.fit = QCheckBox("Fit line")
        self.fit.setToolTip("Least-squares line; one per category when coloured by text")
        self.curve = QComboBox()
        self.curve.addItem("None", None)
        for key, title in CURVES.items():
            self.curve.addItem(title, key)
        self.guide = QComboBox()
        self.guide.addItem("None", None)
        for key, guide in GUIDES.items():
            self.guide.addItem(guide.title, key)
            self.guide.setItemData(
                self.guide.count() - 1,
                f"Expects X: {guide.x_hint}, Y: {guide.y_hint}",
                Qt.ItemDataRole.ToolTipRole,
            )
        self.guide.setToolTip("Reference regions or lines drawn behind the points")
        compact(self.x, self.y, self.color, self.size, self.curve, self.guide, chars=16)
        #: Settings button on the dock title bar opens this (``PanelDock`` looks for ``options``)
        self.options = PlotOptions(self)
        form = self.options.form
        form.addRow("Plot:", self.kind)
        self._labels: dict[QWidget, QLabel] = {}
        for title, widget in (
            ("X", self.x),
            ("Y", self.y),
            ("Color", self.color),
            ("Size", self.size),
            ("Bins", self.bins),
            ("Curve", self.curve),
            ("Guide", self.guide),
        ):
            label = QLabel(title + ":")
            self._labels[widget] = label
            form.addRow(label, widget)
        form.addRow("", self.fit)
        #: what was drawn: point counts, rows left out, fitted curves and their parameters
        self.summary = QLabel(" ")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.summary.setStyleSheet("color: palette(placeholder-text);")
        form.addRow(self.summary)

        # A slim toolbar joined to the top of the plot: icon buttons, the mode ones exclusive.
        self.toolbar = QToolBar()
        self.toolbar.setIconSize(QSize(16, 16))
        self.toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.toolbar.setMovable(False)
        self.mode_buttons: dict[str, QAction] = {}
        self._glyphs: dict[QAction, str] = {}
        modes = QActionGroup(self)
        for mode, icon, text, tip in (
            (ZOOM, "pan", "Pan and Zoom",
             "Drag to pan, scroll or right-drag to zoom; click a point to select it"),
            (BOX, "box", "Box Select", "Drag a rectangle to select (Shift adds, ⌘/Ctrl toggles)"),
            (LASSO, "lasso", "Lasso Select",
             "Draw around points to select (Shift adds, ⌘/Ctrl toggles)"),
        ):  # fmt: skip
            action = self.toolbar.addAction(text)
            action.setToolTip(f"{text}: {tip}")
            action.setCheckable(True)
            modes.addAction(action)
            action.toggled.connect(lambda on, m=mode: on and self.set_mode(m))
            self.mode_buttons[mode] = action
            self._glyphs[action] = icon
        self.toolbar.addSeparator()
        reset = self.toolbar.addAction("Reset View", self.reset_view)
        reset.setToolTip("Reset View: fit the plot to the data")
        export = self.toolbar.addAction("Export…", self.export_dialog)
        export.setToolTip("Export: save the plot as PNG or SVG")
        self._glyphs[reset] = "fit"
        self._glyphs[export] = "export"

        self.view_box = SelectViewBox()
        self.plot_widget = pg.PlotWidget(viewBox=self.view_box)
        self.plot = self.plot_widget.getPlotItem()
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.setMenuEnabled(False)
        self._apply_palette()
        self.card = HoverCard(self)

        self.plot_widget.setFrameShape(QFrame.Shape.NoFrame)
        frame = QFrame()  # one bordered box round the toolbar and the plot
        frame.setObjectName("PlotFrame")
        inside = QVBoxLayout(frame)
        inside.setContentsMargins(1, 1, 1, 1)
        inside.setSpacing(0)
        inside.addWidget(self.toolbar)
        inside.addWidget(self.plot_widget, 1)
        self._frame = frame
        self._style_toolbar()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(frame, 1)

        self.renderer = None
        self._rebuild_timer = QTimer(self)
        self._rebuild_timer.setSingleShot(True)
        self._rebuild_timer.setInterval(REBUILD_DELAY_MS)
        self._rebuild_timer.timeout.connect(self.rebuild)
        self._fill_columns(initial=True)
        for combo in (self.kind, self.x, self.y, self.color, self.size):
            combo.currentIndexChanged.connect(self._controls_changed)
        self.bins.valueChanged.connect(self._controls_changed)
        self.curve.currentIndexChanged.connect(self._controls_changed)
        self.guide.currentIndexChanged.connect(self._controls_changed)
        self.fit.toggled.connect(self._controls_changed)
        self.view_box.clicked.connect(self._on_click)
        self.view_box.boxSelected.connect(
            lambda x0, y0, x1, y1, mods: self._select(self.renderer.hit_box(x0, y0, x1, y1), mods)
        )
        self.view_box.lassoSelected.connect(
            lambda polygon, mods: self._select(self.renderer.hit_lasso(polygon), mods)
        )
        self._hover_proxy = pg.SignalProxy(
            self.plot_widget.scene().sigMouseMoved, rateLimit=30, slot=self._on_mouse_moved
        )
        self._unsubscribe_table = table.subscribe(self._on_table_event)
        self._unsubscribe_selection = selection.subscribe(self._on_selection)
        self.mode_buttons[ZOOM].setChecked(True)
        self.rebuild()

    # -- controls -------------------------------------------------------------------------------
    def spec(self) -> PlotSpec:
        def name(combo: QComboBox) -> str | None:
            text = combo.currentText()
            return text if text and text != NONE else None

        kind = self.kind.currentData()
        return PlotSpec(
            kind=kind,
            x=name(self.x),
            y=name(self.y) if kind in ("scatter", "box") else None,
            color=name(self.color) if kind == "scatter" else None,
            size=name(self.size) if kind == "scatter" else None,
            bins=self.bins.value() or None,
            fit=self.fit.isChecked() and kind == "scatter",
            curve=self.curve.currentData() if kind == "histogram" else None,
            guide=self.guide.currentData() if kind == "scatter" else None,
        )

    def set_spec(self, **values) -> None:
        """Set controls by name (x, y, color, size, kind, bins, fit, curve, guide); redraw now."""
        widgets = {"x": self.x, "y": self.y, "color": self.color, "size": self.size}
        controls = (
            self.kind,
            self.x,
            self.y,
            self.color,
            self.size,
            self.bins,
            self.fit,
            self.curve,
            self.guide,
        )
        for control in controls:  # one redraw at the end, not one per control
            control.blockSignals(True)
        try:
            if "kind" in values:
                self.kind.setCurrentIndex(self.kind.findData(values.pop("kind")))
                self._fill_columns()
            for key, value in values.items():
                if key in widgets:
                    widgets[key].setCurrentText(value or NONE)
                elif key == "bins":
                    self.bins.setValue(value or 0)
                elif key == "fit":
                    self.fit.setChecked(bool(value))
                elif key == "curve":
                    self.curve.setCurrentIndex(max(0, self.curve.findData(value)))
                elif key == "guide":
                    self.guide.setCurrentIndex(max(0, self.guide.findData(value)))
        finally:
            for control in controls:
                control.blockSignals(False)
        self.rebuild()

    def state(self) -> dict:
        """The plot's settings, for saving with the document (JSON-friendly)."""
        spec = self.spec()
        return {
            "kind": spec.kind,
            "x": spec.x,
            "y": spec.y,
            "color": spec.color,
            "size": spec.size,
            "bins": spec.bins,
            "fit": spec.fit,
            "curve": spec.curve,
            "guide": spec.guide,
        }

    def set_state(self, state: dict) -> None:
        """Apply saved settings; names of columns that no longer exist are ignored."""
        known = ("kind", "x", "y", "color", "size", "bins", "fit", "curve", "guide")
        values = {key: state[key] for key in known if key in state}
        if values.get("kind") not in KINDS:
            values.pop("kind", None)
        self.set_spec(**values)

    def _fill_columns(self, *, initial: bool = False) -> None:
        table = self.table
        numbers = number_columns(table)
        kind = self.kind.currentData()
        texts = [n for n in table.column_names if table.kind(n) == TEXT]
        x_choices = {
            "scatter": numbers,
            "histogram": [n for n in table.column_names if table.kind(n) in (NUMBER, TEXT)],
            "box": [NONE, *texts],  # the grouping column (optional)
        }[kind]
        lists = {
            self.x: x_choices,
            self.y: numbers,
            self.color: [NONE, *colorable_columns(table)],
            self.size: [NONE, *numbers],
        }
        for combo, names in lists.items():
            current = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(names)
            if current in names:
                combo.setCurrentText(current)
            elif initial and combo is self.y and len(names) > 1:
                combo.setCurrentIndex(1)
            combo.blockSignals(False)
        shown = {
            "scatter": (self.y, self.color, self.size, self.fit, self.guide),
            "histogram": (self.bins, self.curve),
            "box": (self.y,),
        }[kind]
        for widget in (self.y, self.color, self.size, self.fit, self.bins, self.curve, self.guide):
            self.options.form.setRowVisible(widget, widget in shown)
        self._labels[self.x].setText("Group by:" if kind == "box" else "X:")
        self._labels[self.y].setText("Values:" if kind == "box" else "Y:")

    def _controls_changed(self, *_args) -> None:
        if self.sender() is self.kind:
            self._fill_columns()
        self.rebuild()

    def set_mode(self, mode: str) -> None:
        self.view_box.set_mode(mode)
        if not self.mode_buttons[mode].isChecked():
            self.mode_buttons[mode].setChecked(True)

    # -- drawing --------------------------------------------------------------------------------
    def rebuild(self) -> None:
        self._rebuild_timer.stop()
        spec = self.spec()
        renderer_type = RENDERERS[spec.kind]
        if not isinstance(self.renderer, renderer_type):
            if self.renderer is not None:
                self.renderer.clear()
            self.renderer = renderer_type(self.plot)
        self.summary.setText(self.renderer.build(self.table, spec, shown_oids(self.table)))
        self.renderer.show_selection(self.selection.oids)
        axes = (spec.kind, spec.x, spec.y, spec.bins)
        if axes != self._last_axes:
            self._last_axes = axes
            self.reset_view()
        self.titleChanged.emit(self.title())

    def title(self) -> str:
        spec = self.spec()
        if spec.kind == "scatter" and spec.x and spec.y:
            return f"{spec.y} vs {spec.x}"
        if spec.kind == "box" and spec.y:
            return f"{spec.y} by {spec.x}" if spec.x else f"Box Plot: {spec.y}"
        if spec.kind == "histogram" and spec.x:
            return f"Histogram: {spec.x}"
        return KINDS[spec.kind]

    def reset_view(self) -> None:
        self.plot.enableAutoRange()
        self.plot.autoRange()

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "plot"):
            self._apply_palette()
            self._style_toolbar()
            self._rebuild_timer.start()  # guide labels and marker outlines use the palette

    def _style_toolbar(self) -> None:
        """Toolbar and frame colours from the palette (called again when the theme changes)."""
        palette = self.palette()
        ink = palette.text().color()
        for action, name in self._glyphs.items():
            action.setIcon(glyph(name, ink))
        base = palette.base().color().name()
        line = palette.mid().color().name()
        checked = QColor(palette.highlight().color())
        checked.setAlphaF(0.28)
        hover = QColor(ink)
        hover.setAlphaF(0.08)
        rgba = "rgba({}, {}, {}, {:.2f})".format
        self._frame.setStyleSheet(
            f"QFrame#PlotFrame {{ border: 1px solid {line}; border-radius: 4px;"
            f" background: {base}; }}"
        )
        self.toolbar.setStyleSheet(
            f"QToolBar {{ background: {base}; border: none;"
            f" border-bottom: 1px solid {line}; padding: 2px 4px; spacing: 2px; }}"
            " QToolButton { border: none; border-radius: 4px; padding: 3px; }"
            f" QToolButton:hover {{ background: {rgba(*hover.getRgb()[:3], hover.alphaF())}; }}"
            " QToolButton:checked {"
            f" background: {rgba(*checked.getRgb()[:3], checked.alphaF())}; }}"
        )

    def _apply_palette(self) -> None:
        palette = self.palette()
        self.plot_widget.setBackground(palette.base().color())
        pen = pg.mkPen(palette.text().color())
        for axis in ("left", "bottom"):
            self.plot.getAxis(axis).setPen(pen)
            self.plot.getAxis(axis).setTextPen(pen)

    # -- following the table and the selection --------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        if event.structural:
            self._fill_columns()
        self._rebuild_timer.start()

    def _on_selection(self, oids: frozenset[int], origin: str) -> None:
        if origin != self._origin and self.renderer is not None:
            self.renderer.show_selection(oids)

    # -- selecting ------------------------------------------------------------------------------
    def _select(self, oids: np.ndarray, modifiers) -> None:
        mods = modifiers
        if mods & Qt.KeyboardModifier.ShiftModifier:
            self.selection.add(oids.tolist(), self._origin)
        elif mods & Qt.KeyboardModifier.ControlModifier:
            current = self.selection.oids
            hit = {int(o) for o in oids}
            self.selection.set(current ^ hit, self._origin)
        else:
            self.selection.set(oids.tolist(), self._origin)
        self.renderer.show_selection(self.selection.oids)

    def _on_click(self, x: float, y: float, modifiers) -> None:
        sx, sy = self.view_box.viewPixelSize()
        self._select(self.renderer.hit_point(x, y, sx, sy), modifiers)

    # -- hovering -------------------------------------------------------------------------------
    def _on_mouse_moved(self, args) -> None:
        pos = args[0]
        if self.renderer is None or not self.view_box.sceneBoundingRect().contains(pos):
            self.card.hide()
            return
        point = self.view_box.mapSceneToView(pos)
        sx, sy = self.view_box.viewPixelSize()
        hover = self.renderer.hover(self.table, point.x(), point.y(), sx, sy)
        if hover is None:
            self.card.hide()
            return
        self.show_hover(hover)

    def show_hover(self, hover: Hover) -> None:
        lines = list(hover.lines)
        title = hover.title
        if hover.oid is not None:
            label = next((n for n in self.table.column_names if self.table.kind(n) == TEXT), None)
            if label:
                title = self.table.cell(label, hover.oid)
        if title:
            lines.insert(0, title)
        image = self._images(hover.oid) if self._images and hover.oid is not None else None
        self.card.show_at(QCursor.pos(), image, lines)

    def leaveEvent(self, event) -> None:
        self.card.hide()
        super().leaveEvent(event)

    # -- export and closing ---------------------------------------------------------------------
    def export_dialog(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Plot", f"{self.title()}.png", "PNG image (*.png);;SVG drawing (*.svg)"
        )
        if path:
            self.export(Path(path))

    def export(self, path: Path) -> None:
        from pyqtgraph.exporters import ImageExporter, SVGExporter

        exporter = SVGExporter if path.suffix.lower() == ".svg" else ImageExporter
        exporter(self.plot).export(str(path))

    def close_panel(self) -> None:
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self._rebuild_timer.stop()
        self.card.hide()
        self.options.hide()
