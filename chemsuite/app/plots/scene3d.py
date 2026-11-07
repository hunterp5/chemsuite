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
"""The 3D plots of a plot panel: a cube you turn, drawn with OpenGL (pyqtgraph.opengl).

* **3D Scatter**: three number columns, coloured and sized by others, with optional drop
  lines to the floor so heights read at a glance.
* **3D Landscape**: an activity landscape, a smooth surface of one column (pIC50) over two
  others (usually a chemical-space map), with the compounds standing on it.

Each axis is scaled to the cube's side, so columns with very different ranges share the space.
The panel's toolbar drives it: Turn (drag; scroll zooms), Box and Lasso select what you see,
Reset View, Export. Selection and hover go through the panel exactly as for the 2D plots. Always
on a white background, like a printed figure.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QMatrix4x4, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget

from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import BLANK_LABEL, OTHER_LABEL, categories, numbers
from chemsuite.plotting.geometry import lasso_mask, project
from chemsuite.plotting.stats import landscape
from chemsuite.plotting.style import (
    DEFAULT_COLOR,
    GREY,
    category_color,
    marker_sizes,
    scale_bins,
    scale_colors,
)

from .view_box import BOX, LASSO, ZOOM

POINT_SIZE = 7.0
HIGHLIGHT = (230, 70, 20)
DIMMED = 0.25  # opacity of unselected points while something is selected
HOVER_RADIUS = 9.0
CLICK_SLOP = 4  # pixels the mouse may move and still count as a click
SCALE_STEPS = 24
BACKGROUND = (255, 255, 255)
INK = (70, 70, 70)  # the cube's edges and labels
LABEL_OFFSET = 0.14  # how far labels sit outside the cube
FLOOR_LINES = 4  # grid lines across the floor, each way
DROP_MAX = 20_000  # more points than this: drop lines only for the selected ones
SPECIAL = (OTHER_LABEL, BLANK_LABEL)  # categories drawn grey
CAMERA = {"distance": 4.6, "elevation": 22.0, "azimuth": -55.0}
KINDS_3D = ("scatter3d", "landscape")
#: the cube's corners are at ±1 on every axis
CORNERS = np.array(
    [
        [a, b, c]
        for a, b, c in (
            (-1, -1, -1), (1, -1, -1), (-1, 1, -1), (1, 1, -1),
            (-1, -1, 1), (1, -1, 1), (-1, 1, 1), (1, 1, 1),
        )
    ],
    dtype=float,
)  # fmt: skip
EDGE_PAIRS = (
    (0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3), (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7),
)  # fmt: skip


def scaled(values: np.ndarray, lo: float | None = None, hi: float | None = None):
    """*values* mapped onto -1..1 along [lo, hi] (default: their range; a constant column sits
    at 0): ``(mapped, lo, hi)``."""
    finite = values[np.isfinite(values)]
    if lo is None:
        lo = float(finite.min()) if len(finite) else 0.0
    if hi is None:
        hi = float(finite.max()) if len(finite) else 0.0
    span = hi - lo
    mapped = 2 * (values - lo) / span - 1 if span > 1e-12 else np.zeros(len(values))
    return mapped, lo, hi


def _matrix(m: QMatrix4x4) -> np.ndarray:
    return np.array(m.data(), dtype=float).reshape(4, 4).T  # Qt stores columns first


def _rgba(rgb, alpha: float = 1.0) -> tuple[float, float, float, float]:
    return (rgb[0] / 255, rgb[1] / 255, rgb[2] / 255, alpha)


try:  # PyOpenGL is a dependency, but a broken install should not stop the app
    import pyqtgraph.opengl as gl
except ImportError:  # pragma: no cover - only without PyOpenGL
    gl = None

AVAILABLE = gl is not None


class _Trace(QWidget):
    """The box or lasso being drawn, painted over the OpenGL view."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.points: list[QPointF] = []
        self.box: QRect | None = None
        self.hide()

    def paintEvent(self, _event) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        pen = QPen(QColor(*HIGHLIGHT), 1.2, Qt.PenStyle.DashLine)
        painter.setPen(pen)
        painter.setBrush(QColor(*HIGHLIGHT, 30))
        if self.box is not None:
            painter.drawRect(self.box)
        elif len(self.points) > 1:
            path = QPainterPath(self.points[0])
            for point in self.points[1:]:
                path.lineTo(point)
            path.closeSubpath()
            painter.drawPath(path)


if gl is not None:

    class View3D(gl.GLViewWidget):
        """The OpenGL view: turns on drag (pyqtgraph), reports clicks, boxes, lassos, hovers and
        camera moves (so the panel can move the cards beside selected points)."""

        clicked = Signal(QPoint, object)
        boxed = Signal(QRect, object)
        lassoed = Signal(object, object)  # (n × 2 widget pixels, modifiers)
        hovered = Signal(QPoint)
        turning = Signal()  # a drag started
        moved = Signal()  # the camera settled (drag released or zoomed)

        def __init__(self) -> None:
            super().__init__()
            self.mode = ZOOM
            self._press: QPoint | None = None
            self._trace = _Trace(self)
            self.setMouseTracking(True)
            self.setBackgroundColor(QColor(*BACKGROUND))

        def resizeEvent(self, event) -> None:
            super().resizeEvent(event)
            self._trace.resize(self.size())
            self.moved.emit()

        def mousePressEvent(self, event) -> None:
            self._press = event.position().toPoint()
            if event.button() == Qt.MouseButton.LeftButton and self.mode in (BOX, LASSO):
                self._trace.box = QRect(self._press, QSize()) if self.mode == BOX else None
                self._trace.points = [QPointF(self._press)]
                self._trace.show()
                self._trace.raise_()
                return
            self.turning.emit()
            super().mousePressEvent(event)

        def mouseMoveEvent(self, event) -> None:
            pos = event.position().toPoint()
            if event.buttons() == Qt.MouseButton.NoButton:
                self.hovered.emit(pos)
                return
            if self._trace.isVisible() and self._press is not None:
                if self._trace.box is not None:
                    self._trace.box = QRect(self._press, pos).normalized()
                else:
                    self._trace.points.append(QPointF(pos))
                self._trace.update()
                return
            super().mouseMoveEvent(event)

        def mouseReleaseEvent(self, event) -> None:
            pos = event.position().toPoint()
            press, self._press = self._press, None
            if self._trace.isVisible():
                self._trace.hide()
                box, points = self._trace.box, self._trace.points
                self._trace.box, self._trace.points = None, []
                if box is not None and (box.width() > CLICK_SLOP or box.height() > CLICK_SLOP):
                    self.boxed.emit(box, event.modifiers())
                    return
                if box is None and len(points) > 2:
                    polygon = np.array([[p.x(), p.y()] for p in points])
                    span = polygon.max(axis=0) - polygon.min(axis=0)
                    if span.max() > CLICK_SLOP:
                        self.lassoed.emit(polygon, event.modifiers())
                        return
            if press is not None and (pos - press).manhattanLength() <= CLICK_SLOP:
                self.clicked.emit(pos, event.modifiers())
            super().mouseReleaseEvent(event)
            self.moved.emit()

        def wheelEvent(self, event) -> None:
            super().wheelEvent(event)
            self.moved.emit()

        def camera(self) -> dict:
            params = self.cameraParams()
            return {key: float(params[key]) for key in CAMERA}

        def set_camera(self, camera: object) -> None:
            if not isinstance(camera, dict):
                return
            values = {
                k: float(v) for k, v in camera.items() if k in CAMERA and isinstance(v, int | float)
            }
            self.setCameraPosition(**values)
            self.moved.emit()

        def reset_camera(self) -> None:
            from pyqtgraph import Vector

            self.setCameraPosition(pos=Vector(0, 0, 0), **CAMERA)
            self.moved.emit()

        def export(self, path) -> bool:
            return self.grabFramebuffer().save(str(path))

else:  # pragma: no cover - only without PyOpenGL
    View3D = None


class Scene3D:
    """What a 3D plot draws in a ``View3D``; the 3D counterpart of the 2D renderers.

    ``oids``/``pos`` are the plotted points (``pos`` in cube coordinates, ±1); the hit tests
    work on where those points are on screen, so they select what you see.
    """

    def __init__(self, view) -> None:
        self.view = view
        self.kind = "scatter3d"
        self.oids = np.empty(0, dtype=np.int64)
        self.pos = np.empty((0, 3))
        self.colors = np.empty((0, 4))
        self.sizes: np.ndarray | float = POINT_SIZE
        self.ranges: list[tuple[float, float]] = []
        self.names: tuple[str | None, ...] = (None, None, None)
        self.drop = False
        self._labels: list = []
        self._frame = gl.GLLinePlotItem(mode="lines", antialias=True)
        self._floor = gl.GLLinePlotItem(mode="lines", antialias=True)
        self._drops = gl.GLLinePlotItem(mode="lines", antialias=True)
        self._surface = gl.GLMeshItem(smooth=True, shader=None, glOptions="opaque")
        self._points = gl.GLScatterPlotItem(pxMode=True)
        self._points.setGLOptions("translucent")
        self._chosen = gl.GLScatterPlotItem(pxMode=True)
        self._chosen.setGLOptions("translucent")
        for item in (
            self._floor,
            self._frame,
            self._surface,
            self._drops,
            self._points,
            self._chosen,
        ):
            view.addItem(item)

    # -- drawing --------------------------------------------------------------------------------
    def build(self, table: Table, spec, oids: np.ndarray) -> str:
        """Draw *spec* (a ``PlotSpec`` of a 3D kind); returns the status line."""
        self.kind = spec.kind
        self.names = (spec.x, spec.y, spec.z)
        self.drop = bool(spec.drop) and spec.kind == "scatter3d"
        names = self.names
        self._surface.setVisible(False)
        if not all(names) or any(n not in table for n in names):
            self._empty()
            return "Choose number columns for X, Y and Z."
        values = [numbers(table, n, oids) for n in names]
        ok = np.isfinite(values[0]) & np.isfinite(values[1]) & np.isfinite(values[2])
        self.oids = oids[ok]
        if not ok.any():
            self._empty()
            return "No rows have all three values."
        axes = [scaled(v[ok]) for v in values]
        self.pos = np.column_stack([a[0] for a in axes])
        self.ranges = [(a[1], a[2]) for a in axes]
        missing = int((~ok).sum())
        status = f"{len(self.oids):,} points" + (
            f" · {missing:,} rows without all three values" if missing else ""
        )
        if spec.kind == "landscape":
            self.colors = self._height_colours(self.pos[:, 2])
            self.sizes = POINT_SIZE - 2
            cells = self._draw_surface(values, ok, spec.bins)
            status += f" · surface on a {cells}×{cells} grid" if cells else ""
        else:
            self.colors = self._colours(table, spec.color)
            self.sizes = self._sizes(table, spec.size)
        self._draw()
        return status

    def _empty(self) -> None:
        self.oids = np.empty(0, np.int64)
        self.pos = np.empty((0, 3))
        self.colors = np.empty((0, 4))
        self.ranges = []
        self._draw()

    def _colours(self, table: Table, column: str | None) -> np.ndarray:
        n = len(self.oids)
        rgb = np.tile(np.array(DEFAULT_COLOR, dtype=float), (n, 1))
        if column and column in table and n:
            if table.kind(column) == NUMBER:
                bins, _lo, _hi = scale_bins(numbers(table, column, self.oids), SCALE_STEPS)
                steps = np.array(scale_colors(SCALE_STEPS), dtype=float)
                rgb = np.where((bins >= 0)[:, None], steps[np.clip(bins, 0, None)], GREY)
            else:
                cats = categories(table, column, self.oids)
                palette = np.array(
                    [
                        category_color(c, label, special=SPECIAL)
                        for c, label in enumerate(cats.labels)
                    ],
                    dtype=float,
                )
                rgb = palette[cats.codes]
        return np.column_stack([rgb / 255.0, np.full(n, 0.85)])

    def _height_colours(self, heights: np.ndarray) -> np.ndarray:
        bins, _lo, _hi = scale_bins(heights, SCALE_STEPS, -1.0, 1.0)
        steps = np.array(scale_colors(SCALE_STEPS), dtype=float) / 255.0
        rgb = np.where((bins >= 0)[:, None], steps[np.clip(bins, 0, None)], np.array(GREY) / 255)
        return np.column_stack([rgb * 0.75, np.full(len(heights), 0.95)])  # darker than the sheet

    def _sizes(self, table: Table, column: str | None) -> np.ndarray | float:
        if column and column in table and table.kind(column) == NUMBER and len(self.oids):
            return marker_sizes(numbers(table, column, self.oids), 4.0, 16.0)
        return POINT_SIZE

    def _draw_surface(self, values, ok, bins: int | None) -> int:
        """The activity landscape mesh, in the points' cube coordinates; returns its grid size."""
        x, y, z = (v[ok] for v in values)
        surface = landscape(x, y, z, bins)
        if surface is None:
            return 0
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = self.ranges
        xs = scaled(surface.xs, xlo, xhi)[0]
        ys = scaled(surface.ys, ylo, yhi)[0]
        heights = scaled(surface.heights.ravel(), zlo, zhi)[0].reshape(surface.heights.shape)
        nx, ny = heights.shape
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        known = np.isfinite(heights)
        vertices = np.column_stack([gx.ravel(), gy.ravel(), np.where(known, heights, -1).ravel()])
        index = np.arange(nx * ny).reshape(nx, ny)
        a, b = index[:-1, :-1].ravel(), index[1:, :-1].ravel()
        c, d = index[1:, 1:].ravel(), index[:-1, 1:].ravel()
        whole = known.ravel()
        quad = whole[a] & whole[b] & whole[c] & whole[d]
        faces = np.concatenate([np.column_stack([a, b, c])[quad], np.column_stack([a, c, d])[quad]])
        bins_, _lo, _hi = scale_bins(vertices[:, 2], SCALE_STEPS, -1.0, 1.0)
        steps = np.array(scale_colors(SCALE_STEPS), dtype=float) / 255.0
        colours = np.column_stack([steps[np.clip(bins_, 0, None)], np.ones(len(vertices))])
        mesh = gl.MeshData(vertexes=vertices, faces=faces, vertexColors=colours)
        self._surface.setMeshData(meshdata=mesh, drawEdges=False)
        self._surface.setVisible(len(faces) > 0)
        return nx

    def _draw(self) -> None:
        ink = _rgba(INK, 0.35)
        self._frame.setData(pos=CORNERS[np.array(EDGE_PAIRS).ravel()], color=ink, width=1.0)
        ticks = np.linspace(-1, 1, FLOOR_LINES + 2)[1:-1]
        floor = []
        for t in ticks:
            floor += [(t, -1, -1), (t, 1, -1), (-1, t, -1), (1, t, -1)]
        self._floor.setData(pos=np.array(floor, dtype=float), color=_rgba(INK, 0.12), width=1.0)
        for label in self._labels:
            self.view.removeItem(label)
        self._labels = []
        if len(self.ranges) == 3 and len(self.oids):
            font = QFont()
            font.setPointSizeF(max(8.0, font.pointSizeF() - 2))
            ink_colour = QColor(*INK)
            start = np.array([-1.0, -1.0, -1.0])
            for axis, (name, (lo, hi)) in enumerate(zip(self.names, self.ranges, strict=True)):
                along = np.zeros(3)
                along[axis] = 2.0
                outward = np.full(3, -LABEL_OFFSET)
                outward[axis] = 0.0  # away from the cube, square to its edge
                for text, t in ((name, 0.5), (f"{lo:.3g}", 0.12), (f"{hi:.3g}", 0.96)):
                    if text is None:
                        continue
                    label = gl.GLTextItem(
                        pos=start + t * along + outward, text=text, color=ink_colour, font=font
                    )
                    self.view.addItem(label)
                    self._labels.append(label)
        self._points.setData(pos=self.pos, color=self.colors, size=self.sizes)
        self.show_selection(frozenset())

    def show_selection(self, selected) -> None:
        mask = self.selected_mask(selected)
        colors = self.colors.copy()
        if mask.any():
            colors[:, 3] = DIMMED
        self._points.setData(color=colors)
        n = int(mask.sum())
        sizes = self.sizes[mask] + 3 if isinstance(self.sizes, np.ndarray) else self.sizes + 3
        self._chosen.setData(
            pos=self.pos[mask] if n else np.empty((0, 3)),
            color=np.tile(_rgba(HIGHLIGHT), (n, 1)) if n else np.empty((0, 4)),
            size=sizes if n else POINT_SIZE,
        )
        self._draw_drops(mask)

    def _draw_drops(self, mask: np.ndarray) -> None:
        if not self.drop or not len(self.pos):
            self._drops.setVisible(False)
            return
        rows = self.pos if len(self.pos) <= DROP_MAX else self.pos[mask]
        tops = rows
        feet = rows.copy()
        feet[:, 2] = -1.0
        lines = np.empty((2 * len(rows), 3))
        lines[0::2], lines[1::2] = tops, feet
        self._drops.setData(pos=lines, color=_rgba(INK, 0.22), width=1.0)
        self._drops.setVisible(len(rows) > 0)

    def selected_mask(self, selected) -> np.ndarray:
        chosen = np.fromiter((int(o) for o in selected), dtype=np.int64, count=len(selected))
        if not len(chosen) or not len(self.oids):
            return np.zeros(len(self.oids), dtype=np.bool_)
        return np.isin(self.oids, chosen)

    # -- hits -----------------------------------------------------------------------------------
    def screen_points(self) -> tuple[np.ndarray, np.ndarray]:
        """Where each point is on the view, in widget pixels, and whether it is in front."""
        view = self.view
        if not len(self.pos):
            return np.empty((0, 2)), np.zeros(0, dtype=np.bool_)
        width, height = max(view.width(), 1), max(view.height(), 1)
        viewport = (0, 0, width, height)
        mvp = view.projectionMatrix(viewport, viewport) * view.viewMatrix()
        return project(self.pos, _matrix(mvp), width, height)

    def hit_box(self, rect: QRect) -> np.ndarray:
        xy, front = self.screen_points()
        inside = (
            front
            & (xy[:, 0] >= rect.left())
            & (xy[:, 0] <= rect.right())
            & (xy[:, 1] >= rect.top())
            & (xy[:, 1] <= rect.bottom())
        )
        return self.oids[inside]

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        xy, front = self.screen_points()
        if not len(xy):
            return self.oids[:0]
        return self.oids[front & lasso_mask(xy[:, 0], xy[:, 1], polygon)]

    def hit_index(self, pos: QPoint) -> int | None:
        """The index of the point under *pos* (the nearest, within ``HOVER_RADIUS``)."""
        xy, front = self.screen_points()
        if not front.any():
            return None
        distance = np.hypot(xy[:, 0] - pos.x(), xy[:, 1] - pos.y())
        distance[~front] = np.inf
        i = int(np.argmin(distance))
        return i if distance[i] <= HOVER_RADIUS else None

    def hit_point(self, pos: QPoint) -> np.ndarray:
        i = self.hit_index(pos)
        return self.oids[[i]] if i is not None else self.oids[:0]

    def title(self) -> str | None:
        x, y, z = self.names
        if not (x and y and z):
            return None
        return f"Landscape: {z} over {x}, {y}" if self.kind == "landscape" else f"3D: {x}, {y}, {z}"
