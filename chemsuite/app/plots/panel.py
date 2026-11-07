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

from dataclasses import replace
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QEvent, QPoint, QPointF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QActionGroup, QColor, QCursor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMenu,
    QPushButton,
    QSpinBox,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.medchem_space import GUIDES
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table.columns import NUMBER, STRUCTURE, TEXT
from chemsuite.plotting.data import colorable_columns, number_columns, shown_oids
from chemsuite.plotting.stats import CURVES

from ..flow_layout import compact
from ..glyphs import glyph
from .cards import (
    MAX_PINNED,
    SIZES,
    WHEN,
    CardContent,
    CardOptions,
    HoverCard,
    ImageProvider,
    PinnedCards,
)
from .renderers import RENDERERS, Hover, PlotSpec
from .view_box import BOX, LASSO, ZOOM, SelectViewBox

NONE = "(none)"
KINDS = {"scatter": "Scatter", "histogram": "Histogram", "box": "Box Plot"}
REBUILD_DELAY_MS = 150
#: left, top, right, bottom space round the plot, so the axis titles do not touch the edges
AXIS_MARGINS = (6, 6, 10, 6)


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
        #: what was drawn: point counts, rows left out, fitted curves and their parameters.
        #: Not shown in the options (no explanatory text there); the fit and curve controls
        #: carry it as their tooltip.
        self.summary = QLabel(" ")
        self._build_card_options(form)

        # A slim toolbar joined to the top of the plot: icon buttons, the mode ones exclusive.
        self.toolbar = QToolBar()
        self.toolbar.setIconSize(QSize(15, 15))
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
        # breathing room between the axis titles and the pane's edges
        self.plot.layout.setContentsMargins(*AXIS_MARGINS)
        self._apply_palette()
        self.card = HoverCard(self)
        self.pinned = PinnedCards(self.plot_widget, self._data_to_widget)
        self.view_box.sigRangeChanged.connect(lambda *_: self.pinned.place())
        self.plot_widget.installEventFilter(self)

        self.plot_widget.setFrameShape(QFrame.Shape.NoFrame)
        frame = QFrame()  # one bordered box round the toolbar and the plot
        frame.setObjectName("PlotFrame")
        inside = QVBoxLayout(frame)
        inside.setContentsMargins(1, 1, 1, 1)
        inside.setSpacing(0)
        inside.addWidget(self.toolbar)
        inside.addWidget(self.plot_widget, 1)
        self._frame = frame
        self._header_attached = False
        #: a title the user typed (None: the automatic "pIC50 vs MW")
        self.custom_title: str | None = None
        self._style_toolbar()

        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.addWidget(frame, 1)

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
            "title": self.custom_title,
            "cards": self.card_options.state(),
        }

    def set_state(self, state: dict) -> None:
        """Apply saved settings; names of columns that no longer exist are ignored."""
        title = state.get("title")
        self.custom_title = title.strip() or None if isinstance(title, str) else None
        if "cards" in state:
            self.set_card_options(CardOptions.from_state(state["cards"]))
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
        for control in (self.fit, self.curve):
            control.setToolTip(self.summary.text())
        self.renderer.show_selection(self.selection.oids)
        self.update_pinned()
        axes = (spec.kind, spec.x, spec.y, spec.bins)
        if axes != self._last_axes:
            self._last_axes = axes
            self.reset_view()
        self.titleChanged.emit(self.title())

    def set_custom_title(self, text: str | None) -> None:
        """Name the plot; empty or None goes back to the automatic title."""
        self.custom_title = (text or "").strip() or None
        self.titleChanged.emit(self.title())

    def attach_header(self, attached: bool) -> None:
        """The dock's header now holds the toolbar and meets the top of the plot (one card)."""
        self._header_attached = attached
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._style_toolbar()

    def title(self) -> str:
        return self.custom_title or self.automatic_title()

    def automatic_title(self) -> str:
        """The title from the plotted columns: "pIC50 vs MW", "Histogram: MW"."""
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
        if self._header_attached:  # the toolbar sits in the pane's title bar
            bar = "background: transparent; border: none; padding: 0px; spacing: 2px;"
        else:
            bar = f"background: {base}; border: none; border-bottom: 1px solid {line};"
            bar += " padding: 2px 4px; spacing: 2px;"
        # Flush with the window: no box of its own, the dock divider and title bar frame it.
        self._frame.setStyleSheet(f"QFrame#PlotFrame {{ border: none; background: {base}; }}")
        self.toolbar.setStyleSheet(
            f"QToolBar {{ {bar} }}"
            " QToolButton { border: none; border-radius: 4px; padding: 2px; }"
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
        self.update_pinned()

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
        self.update_pinned()

    def _on_click(self, x: float, y: float, modifiers) -> None:
        sx, sy = self.view_box.viewPixelSize()
        self._select(self.renderer.hit_point(x, y, sx, sy), modifiers)

    # -- hovering -------------------------------------------------------------------------------
    def _on_mouse_moved(self, args) -> None:
        pos = args[0]
        if (
            self.renderer is None
            or not self.card_options.hovering
            or not self.view_box.sceneBoundingRect().contains(pos)
        ):
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
        if hover.oid is None:  # a group of rows (a histogram bar, a box)
            content = CardContent(hover.title, list(hover.lines))
        else:
            content = self.card_content(hover.oid)
        if content.empty():
            self.card.hide()
            return
        self.card.show_at(QCursor.pos(), content)

    # -- point cards ----------------------------------------------------------------------------
    def _build_card_options(self, form: QFormLayout) -> None:
        self.card_options = CardOptions()
        self.card_when = QComboBox()
        for key, title in WHEN.items():
            self.card_when.addItem(title, key)
        self.card_when.setToolTip(
            "When point cards show. Selected: each selected point keeps a card beside it "
            f"(up to {MAX_PINNED} points)."
        )
        compact(self.card_when, chars=16)
        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setFrameShadow(QFrame.Shadow.Plain)
        rule.setEnabled(False)
        form.addRow(rule)
        form.addRow("Point cards:", self.card_when)
        shows = QWidget()
        line = QHBoxLayout(shows)
        line.setContentsMargins(0, 0, 0, 0)
        line.setSpacing(8)
        self.card_structure = QCheckBox("Structure")
        self.card_label = QCheckBox("Label")
        self.card_label.setToolTip("The row's first text column (its name or ID)")
        self.card_values = QCheckBox("Values")
        self.card_values.setToolTip("The plotted columns: X, Y, colour and size")
        for box in (self.card_structure, self.card_label, self.card_values):
            line.addWidget(box)
        line.addStretch(1)
        form.addRow("Show:", shows)
        self.card_columns_button = QPushButton("None")
        self.card_columns_button.setToolTip("More columns to show on the cards")
        self.card_columns_menu = QMenu(self.card_columns_button)
        self.card_columns_menu.aboutToShow.connect(self._fill_card_columns)
        self.card_columns_button.setMenu(self.card_columns_menu)
        self.card_columns_button.setStyleSheet("text-align: left; padding: 3px 8px;")
        form.addRow("More columns:", self.card_columns_button)
        self.card_size = QComboBox()
        for key, (title, _size) in SIZES.items():
            self.card_size.addItem(title, key)
        compact(self.card_size, chars=16)
        form.addRow("Drawing size:", self.card_size)
        self._show_card_options(self.card_options)
        for box in (self.card_structure, self.card_label, self.card_values):
            box.toggled.connect(self._card_controls_changed)
        for combo in (self.card_when, self.card_size):
            combo.currentIndexChanged.connect(self._card_controls_changed)

    def _show_card_options(self, options: CardOptions) -> None:
        controls = (
            self.card_when,
            self.card_size,
            self.card_structure,
            self.card_label,
            self.card_values,
        )
        for control in controls:
            control.blockSignals(True)
        self.card_when.setCurrentIndex(self.card_when.findData(options.when))
        self.card_size.setCurrentIndex(self.card_size.findData(options.size))
        self.card_structure.setChecked(options.structure)
        self.card_label.setChecked(options.label)
        self.card_values.setChecked(options.values)
        for control in controls:
            control.blockSignals(False)
        count = len(options.columns)
        self.card_columns_button.setText(
            "None" if not count else options.columns[0] if count == 1 else f"{count} columns"
        )

    def _card_controls_changed(self, *_args) -> None:
        self.set_card_options(
            CardOptions(
                when=self.card_when.currentData(),
                structure=self.card_structure.isChecked(),
                label=self.card_label.isChecked(),
                values=self.card_values.isChecked(),
                columns=self.card_options.columns,
                size=self.card_size.currentData(),
            )
        )

    def _fill_card_columns(self) -> None:
        menu = self.card_columns_menu
        menu.clear()
        chosen = set(self.card_options.columns)
        for name in self.table.column_names:
            action = menu.addAction(name)
            action.setCheckable(True)
            action.setChecked(name in chosen)
            action.toggled.connect(lambda on, n=name: self._toggle_card_column(n, on))

    def _toggle_card_column(self, name: str, on: bool) -> None:
        columns = [c for c in self.card_options.columns if c != name]
        if on:
            columns.append(name)
        self.set_card_options(replace(self.card_options, columns=tuple(columns)))

    def set_card_options(self, options: CardOptions) -> None:
        """When cards show and what they hold; redraws the pinned cards."""
        self.card_options = options
        self._show_card_options(options)
        if not options.hovering:
            self.card.hide()
        self.update_pinned()

    def card_content(self, oid: int) -> CardContent:
        """The card for one row, as the card options say."""
        options = self.card_options
        table = self.table
        content = CardContent()
        if options.label:
            label = next((n for n in table.column_names if table.kind(n) == TEXT), None)
            if label:
                content.title = table.cell(label, oid)
        shown: list[str] = []
        if options.values:
            spec = self.spec()
            shown = [n for n in (spec.x, spec.y, spec.color, spec.size) if n and n in table]
        for name in options.columns:
            if name in table and name not in shown and table.kind(name) != STRUCTURE:
                shown.append(name)
        content.lines = [f"{name}: {table.cell(name, oid)}" for name in dict.fromkeys(shown)]
        if options.structure and self._images is not None:
            content.image = self._images(oid, SIZES[options.size][1])
        return content

    def _data_to_widget(self, x: float, y: float) -> QPointF:
        scene = self.view_box.mapViewToScene(QPointF(x, y))
        return QPointF(self.plot_widget.mapFromScene(scene))

    def update_pinned(self) -> None:
        """Cards beside the selected points (when the options ask for them)."""
        renderer = self.renderer
        if (
            not hasattr(self, "pinned")
            or renderer is None
            or not self.card_options.pinning
            or not hasattr(renderer, "x")
            or not len(self.selection)
        ):
            if hasattr(self, "pinned"):
                self.pinned.clear()
            return
        mask = np.isin(renderer.oids, self.selection.as_array())
        points = [
            (int(o), float(x), float(y))
            for o, x, y in zip(renderer.oids[mask], renderer.x[mask], renderer.y[mask], strict=True)
        ]
        self.pinned.show(points, self.card_content)

    def eventFilter(self, watched, event) -> bool:
        if watched is self.plot_widget and event.type() == QEvent.Type.Resize:
            self.pinned.place()
        return super().eventFilter(watched, event)

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
        self.pinned.clear()
        self.options.hide()
