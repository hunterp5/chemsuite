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

2D kinds draw with pyqtgraph; 3D kinds (3D Scatter, 3D Landscape) on the GPU in the same panel
(``scene3d``), so switching between any two kinds keeps the pane, its toolbar, title, columns,
cards and selection, and only swaps the drawing.

The panel follows the table: it redraws (debounced) when cells, columns or filters change, and
plots only the rows the table shows. Selecting in the plot sets the shared selection; selecting
elsewhere highlights the matching points here.
"""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import Executor
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QEvent, QPoint, QPointF, QSize, Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QColor, QCursor
from PySide6.QtWidgets import (
    QFileDialog,
    QFrame,
    QLabel,
    QStackedWidget,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table.columns import STRUCTURE, TEXT
from chemsuite.plotting.data import shown_oids
from chemsuite.plotting.stats import AGGREGATES
from chemsuite.sar.constellation import SUMMARIES

from ..shell import graphics
from ..shell.glyphs import glyph
from ..widgets.panel import OptionsPopup, ToolPanel
from .cards import (
    SIZES,
    CardContent,
    CardControls,
    CardOptions,
    HoverCard,
    ImageProvider,
    PinnedCards,
    pinned_size,
)
from .constellation import ConstellationRenderer, CoreFinder, Runner, core_image
from .controls import CONSTELLATION, KINDS, KINDS_3D, STATE_KEYS, PlotControls
from .renderers import RENDERERS
from .spec import Hover, PlotSpec
from .view_box import BOX, LASSO, ZOOM, SelectViewBox

#: what the first mode button does and looks like, in 2D and in 3D
MOVE_MODE = {
    False: (
        "pan",
        "Pan and Zoom",
        "Drag to pan, scroll or right-drag to zoom; click a point to select it",
    ),
    True: ("rotate", "Turn", "Drag to turn, scroll to zoom; click a point to select it"),
}
REBUILD_DELAY_MS = 150
#: left, top, right, bottom space round the plot, so the axis titles do not touch the edges
AXIS_MARGINS = (6, 6, 10, 6)


#: what draws each 2D kind
DRAWN_BY = {**RENDERERS, CONSTELLATION: ConstellationRenderer}
PlotOptions = OptionsPopup  # the plot's options popup


class PlotPanel(ToolPanel):
    renamable = True

    def __init__(
        self,
        table: Table,
        selection: Selection,
        images: ImageProvider | None = None,
        kind: str = "scatter",
        parent: QWidget | None = None,
        *,
        run: Runner | None = None,
        workers: Callable[[int], Executor | None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self.selection = selection
        self._images = images
        self._origin = f"plot:{id(self)}"
        self._last_axes: tuple | None = None
        #: finds a constellation's cores (in the background with *run*)
        self.cores = CoreFinder(table, run, workers, self)
        self.cores.ready.connect(self.rebuild)

        #: Settings button on the dock title bar opens this (``PanelDock`` looks for ``options``)
        self.options = PlotOptions(self)
        self.controls = PlotControls(table, self.options.form, kind, self)
        controls = self.controls
        self.kind, self.x, self.y, self.z = controls.kind, controls.x, controls.y, controls.z
        self.color, self.size, self.bins = controls.color, controls.size, controls.bins
        self.fit, self.curve = controls.fit, controls.curve
        if reason := graphics.problem():
            self._disable_3d(reason)
        #: what was drawn: point counts, rows left out, fitted curves and their parameters.
        #: Not shown in the options (no explanatory text there); the fit and curve controls
        #: carry it as their tooltip.
        self.summary = QLabel(" ")
        self.card_options = CardOptions()
        self.card_controls = CardControls(self.options.form, lambda: self.table.column_names, self)
        self.card_controls.changed.connect(self.set_card_options)

        # A slim toolbar joined to the top of the plot: icon buttons, the mode ones exclusive.
        self.toolbar = QToolBar()
        self.toolbar.setIconSize(QSize(15, 15))
        self.toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.toolbar.setMovable(False)
        self.mode_buttons: dict[str, QAction] = {}
        self._glyphs: dict[QAction, str] = {}
        modes = QActionGroup(self)
        for mode, icon, text, tip in (
            (ZOOM, *MOVE_MODE[False]),
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
        self.plot_widget.setFrameShape(QFrame.Shape.NoFrame)
        #: the 2D plot, and the GPU view once a 3D kind is first shown (made lazily)
        self.stack = QStackedWidget()
        self.stack.addWidget(self.plot_widget)
        self.view3d = None
        self.scene3d = None
        self._camera: dict | None = None  # a saved camera, applied when the view is made
        self.pinned = PinnedCards(self.stack, self._data_to_widget)
        self.view_box.sigRangeChanged.connect(lambda *_: self.pinned.place())
        self.plot_widget.installEventFilter(self)

        frame = QFrame()  # one bordered box round the toolbar and the plot
        frame.setObjectName("PlotFrame")
        inside = QVBoxLayout(frame)
        inside.setContentsMargins(1, 1, 1, 1)
        inside.setSpacing(0)
        inside.addWidget(self.toolbar)
        inside.addWidget(self.stack, 1)
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
        controls.changed.connect(self.rebuild)
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
        return self.controls.spec()

    def is_3d(self) -> bool:
        return self.kind.currentData() in KINDS_3D

    def set_spec(self, **values) -> None:
        """Set controls by name (:data:`STATE_KEYS`, and ``camera``); redraw now."""
        camera = values.pop("camera", None)
        self.controls.set_values(values)
        self.rebuild()
        if isinstance(camera, dict):
            self.set_camera(camera)

    def state(self) -> dict:
        """The plot's settings, for saving with the document (JSON-friendly)."""
        return {
            **self.controls.state(),
            "camera": self.camera(),
            "title": self.custom_title,
            "cards": self.card_options.state(),
        }

    def set_state(self, state: dict) -> None:
        """Apply saved settings; names of columns that no longer exist are ignored."""
        title = state.get("title")
        self.custom_title = title.strip() or None if isinstance(title, str) else None
        if "cards" in state:
            self.set_card_options(CardOptions.from_state(state["cards"]))
        values = {key: state[key] for key in (*STATE_KEYS, "camera") if key in state}
        if values.get("kind") not in KINDS:
            values.pop("kind", None)
        self.set_spec(**values)

    def set_mode(self, mode: str) -> None:
        self.view_box.set_mode(mode)
        if self.view3d is not None:
            self.view3d.mode = mode
        if not self.mode_buttons[mode].isChecked():
            self.mode_buttons[mode].setChecked(True)

    # -- drawing --------------------------------------------------------------------------------
    def rebuild(self) -> None:
        self._rebuild_timer.stop()
        spec = self.spec()
        three = spec.kind in KINDS_3D and self._ensure_3d()
        self._show_page(three)
        if three:
            if self.renderer is not None:  # leave nothing of a 2D kind behind
                self.renderer.clear()
                self.renderer = None
            self.summary.setText(self.scene3d.build(self.table, spec, shown_oids(self.table)))
        else:
            renderer_type = DRAWN_BY.get(spec.kind, RENDERERS["scatter"])
            if not isinstance(self.renderer, renderer_type):
                if self.renderer is not None:
                    self.renderer.clear()
                self.renderer = renderer_type(self.plot)
            oids = shown_oids(self.table)
            if isinstance(self.renderer, ConstellationRenderer):
                self.renderer.data = self.cores.result(spec, oids)
                self.renderer.waiting = self.cores.message
            self.summary.setText(self.renderer.build(self.table, spec, oids))
        for control in (self.fit, self.curve):
            control.setToolTip(self.summary.text())
        self._active().show_selection(self.selection.oids)
        self.update_pinned()
        axes = (spec.kind, spec.x, spec.y, spec.bins, spec.aggregate, spec.columns)
        if isinstance(self.renderer, ConstellationRenderer) and spec.kind == CONSTELLATION:
            axes += (id(self.renderer.data),)  # a new map fits the view to it
        if axes != self._last_axes:
            # a new 2D layout fits the view to it; turning a 3D plot's columns keeps the camera
            moved_in_3d = three and self._last_axes is not None and self._last_axes[0] in KINDS_3D
            self._last_axes = axes
            if not moved_in_3d:
                self.reset_view()
        self.titleChanged.emit(self.title())

    def _active(self):
        """What draws the current kind: a 2D renderer or the 3D scene."""
        return self.scene3d if self.is_3d() and self.scene3d is not None else self.renderer

    def _ensure_3d(self) -> bool:
        """Make the GPU view the first time a 3D kind is shown; False when 3D cannot draw."""
        if self.view3d is not None:
            return True
        if graphics.problem():
            return False
        from .scene3d import Scene3D
        from .view3d import View3D

        self.view3d = View3D()
        self.view3d.failed.connect(self._3d_failed)
        self.view3d.mode = self.view_box.mode
        self.scene3d = Scene3D(self.view3d)
        self.stack.addWidget(self.view3d)
        self.view3d.clicked.connect(
            lambda pos, mods: self._select(self.scene3d.hit_point(pos), mods)
        )
        self.view3d.boxed.connect(lambda rect, mods: self._select(self.scene3d.hit_box(rect), mods))
        self.view3d.lassoed.connect(
            lambda polygon, mods: self._select(self.scene3d.hit_lasso(polygon), mods)
        )
        self.view3d.hovered.connect(self._on_hover_3d)
        self.view3d.turning.connect(self.pinned.clear)
        self.view3d.moved.connect(self.update_pinned)
        self.view3d.installEventFilter(self)
        if self._camera is not None:
            self.view3d.set_camera(self._camera)
            self._camera = None
        else:
            self.view3d.reset_camera()
        return True

    def _disable_3d(self, reason: str) -> None:
        for key in KINDS_3D:
            item = self.kind.model().item(self.kind.findData(key))
            item.setEnabled(False)
            item.setToolTip(f"⚠ {reason}")

    def _3d_failed(self) -> None:
        """The graphics cannot draw 3D after all: say why on the 3D kinds, and draw flat."""
        graphics.set_problem(graphics.NO_GPU)
        self._disable_3d(graphics.NO_GPU)
        if self.is_3d():
            self.set_spec(kind="scatter")

    def _show_page(self, three: bool) -> None:
        """The 2D plot or the 3D view, and the toolbar's first button to match."""
        page = self.view3d if three else self.plot_widget
        if self.stack.currentWidget() is not page:
            self.card.hide()
            self.pinned.clear()
            self.stack.setCurrentWidget(page)
        icon, text, tip = MOVE_MODE[three]
        move = self.mode_buttons[ZOOM]
        move.setText(text)
        move.setToolTip(f"{text}: {tip}")
        self._glyphs[move] = icon
        move.setIcon(glyph(icon, self.palette().text().color()))

    def camera(self) -> dict | None:
        """Where a 3D plot is looked at from (None before a 3D kind was shown)."""
        if self.view3d is not None:
            return self.view3d.camera_state()
        return dict(self._camera) if self._camera else None

    def set_camera(self, camera: object) -> None:
        if not isinstance(camera, dict):
            return
        if self.view3d is not None:
            self.view3d.set_camera(camera)
        else:
            self._camera = dict(camera)

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
        if spec.kind == "scatter3d" and spec.x and spec.y and spec.z:
            return f"3D: {spec.x}, {spec.y}, {spec.z}"
        if spec.kind == "landscape" and spec.x and spec.y and spec.z:
            return f"Landscape: {spec.z} over {spec.x}, {spec.y}"
        if spec.kind == "box" and spec.y:
            return f"{spec.y} by {spec.x}" if spec.x else f"Box Plot: {spec.y}"
        if spec.kind == "histogram" and spec.x:
            return f"Histogram: {spec.x}"
        if spec.kind == "line" and spec.x and spec.y:
            return f"{spec.y} by {spec.x}"
        if spec.kind == "density" and spec.x and spec.y:
            return f"Density: {spec.y} vs {spec.x}"
        if spec.kind == "bar" and spec.x:
            if spec.y:
                return f"{AGGREGATES[spec.aggregate]} {spec.y} by {spec.x}"
            return f"Rows by {spec.x}"
        if spec.kind == CONSTELLATION and spec.structure:
            if spec.color:
                return f"Constellation: {SUMMARIES[spec.summary]} {spec.color}"
            return f"Constellation: {spec.structure}"
        return KINDS[spec.kind]

    def reset_view(self) -> None:
        """Fit the plot to the data (3D: look from the start position). A 2D view keeps fitting
        as the pane is resized until you pan or zoom, so a plot drawn while its pane was still
        tiny never stays tiny."""
        if self.is_3d() and self.view3d is not None:
            self.view3d.reset_camera()
            return
        self.plot.enableAutoRange()

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
            self.controls.fill_columns()
        self._rebuild_timer.start()

    def _on_selection(self, oids: frozenset[int], origin: str) -> None:
        active = self._active()
        if origin != self._origin and active is not None:
            active.show_selection(oids)
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
        self._active().show_selection(self.selection.oids)
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

    def _on_hover_3d(self, pos: QPoint) -> None:
        i = self.scene3d.hit_index(pos) if self.card_options.hovering else None
        if i is None:
            self.card.hide()
            return
        content = self.card_content(int(self.scene3d.oids[i]))
        if content.empty():
            self.card.hide()
            return
        self.card.show_at(self.view3d.mapToGlobal(pos), content)

    def show_hover(self, hover: Hover) -> None:
        if hover.oid is None:  # a group of rows (a histogram bar, a box)
            content = CardContent(hover.title, list(hover.lines))
            if hover.smiles and self.card_options.structure:
                size = SIZES[self.card_options.size][1]
                dark = self.palette().base().color().lightness() < 128
                ratio = self.devicePixelRatioF()
                content.image = core_image(hover.smiles, size.width(), size.height(), ratio, dark)
        else:
            content = self.card_content(hover.oid)
        if content.empty():
            self.card.hide()
            return
        self.card.show_at(QCursor.pos(), content)

    # -- point cards ----------------------------------------------------------------------------
    def set_card_options(self, options: CardOptions) -> None:
        """When cards show and what they hold; redraws the pinned cards."""
        self.card_options = options
        self.card_controls.show(options)
        if not options.hovering:
            self.card.hide()
        self.update_pinned()

    def card_content(self, oid: int, size: QSize | None | bool = True) -> CardContent:
        """The card for one row, as the card options say. *size* is the largest drawing (True:
        the options' size; None: no drawing)."""
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
            plotted = (spec.x, spec.y, spec.z, spec.color, spec.size)
            shown = [n for n in plotted if n and n in table]
        for name in options.columns:
            if name in table and name not in shown and table.kind(name) != STRUCTURE:
                shown.append(name)
        content.lines = [f"{name}: {table.cell(name, oid)}" for name in dict.fromkeys(shown)]
        if size is True:
            size = SIZES[options.size][1]
        if options.structure and self._images is not None and size is not None:
            content.image = self._images(oid, size)
        return content

    def _data_to_widget(self, x: float, y: float) -> QPointF:
        if self.is_3d():  # 3D cards are placed at their points' screen pixels already
            return QPointF(x, y)
        scene = self.view_box.mapViewToScene(QPointF(x, y))
        return QPointF(self.plot_widget.mapFromScene(scene))

    def update_pinned(self) -> None:
        """Cards beside the selected points (when the options ask for them)."""
        renderer = self._active() if hasattr(self, "stack") else None
        three = renderer is not None and renderer is self.scene3d
        if (
            not hasattr(self, "pinned")
            or renderer is None
            or not self.card_options.pinning
            or not (three or hasattr(renderer, "x"))
            or not len(self.selection)
        ):
            if hasattr(self, "pinned"):
                self.pinned.clear()
            return
        mask = np.isin(renderer.oids, self.selection.as_array())
        if three:  # where the points are on screen; behind the camera: no card
            xy, front = renderer.screen_points()
            mask &= front
            xs, ys = xy[:, 0], xy[:, 1]
        else:
            xs, ys = renderer.x, renderer.y
        points = [
            (int(o), float(x), float(y))
            for o, x, y in zip(renderer.oids[mask], xs[mask], ys[mask], strict=True)
        ]
        chosen = SIZES[self.card_options.size][1]
        size = pinned_size(chosen, len(points), self.plot_widget.size())
        self.pinned.show(points, lambda oid: self.card_content(oid, size))

    def eventFilter(self, watched, event) -> bool:
        if watched is self.plot_widget and event.type() == QEvent.Type.Resize:
            self.update_pinned()  # a smaller plot may need smaller cards
        if watched is self.view3d and event.type() == QEvent.Type.Leave:
            self.card.hide()
        return super().eventFilter(watched, event)

    def leaveEvent(self, event) -> None:
        self.card.hide()
        super().leaveEvent(event)

    # -- export and closing ---------------------------------------------------------------------
    def export_dialog(self) -> None:
        kinds = "PNG image (*.png)" if self.is_3d() else "PNG image (*.png);;SVG drawing (*.svg)"
        path, _ = QFileDialog.getSaveFileName(self, "Export Plot", f"{self.title()}.png", kinds)
        if path:
            self.export(Path(path))

    def export(self, path: Path) -> None:
        if self.is_3d() and self.view3d is not None:
            self.view3d.export(path.with_suffix(".png"))
            return
        from pyqtgraph.exporters import ImageExporter, SVGExporter

        exporter = SVGExporter if path.suffix.lower() == ".svg" else ImageExporter
        exporter(self.plot).export(str(path))

    def close_panel(self) -> None:
        self.cores.close()
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self._rebuild_timer.stop()
        self.card.hide()
        self.pinned.clear()
        self.options.hide()
