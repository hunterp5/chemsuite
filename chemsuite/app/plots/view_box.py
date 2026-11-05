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
"""A pyqtgraph ViewBox with selection gestures: click, box and lasso.

In *zoom* mode dragging pans (right-drag zooms, the wheel zooms); in *box* and *lasso* mode a
left drag draws the selection shape instead. Clicks are reported in every mode.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QColor, QPainterPath, QPen
from PySide6.QtWidgets import QGraphicsPathItem

ZOOM, BOX, LASSO = "zoom", "box", "lasso"


class SelectViewBox(pg.ViewBox):
    #: (x, y, modifiers) in data coordinates
    clicked = Signal(float, float, object)
    #: (x0, y0, x1, y1, modifiers)
    boxSelected = Signal(float, float, float, float, object)
    #: (k × 2 polygon, modifiers)
    lassoSelected = Signal(object, object)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.mode = ZOOM
        self._points: list[QPointF] = []
        self._shape: QGraphicsPathItem | None = None

    def set_mode(self, mode: str) -> None:
        self.mode = mode
        self.setMouseMode(pg.ViewBox.PanMode)

    def mouseClickEvent(self, ev) -> None:
        if ev.button() == Qt.MouseButton.LeftButton:
            ev.accept()
            pos = self.mapSceneToView(ev.scenePos())
            self.clicked.emit(pos.x(), pos.y(), ev.modifiers())
            return
        super().mouseClickEvent(ev)

    def mouseDragEvent(self, ev, axis=None) -> None:
        if self.mode == ZOOM or ev.button() != Qt.MouseButton.LeftButton:
            super().mouseDragEvent(ev, axis)
            return
        ev.accept()
        point = self.mapSceneToView(ev.scenePos())
        if ev.isStart():
            start = self.mapSceneToView(ev.buttonDownScenePos())
            self._points = [start]
            self._shape = QGraphicsPathItem()
            pen = QPen(QColor(230, 90, 20))
            pen.setCosmetic(True)
            pen.setWidthF(1.5)
            pen.setStyle(Qt.PenStyle.DashLine)
            self._shape.setPen(pen)
            self._shape.setBrush(QColor(230, 90, 20, 30))
            self.addItem(self._shape, ignoreBounds=True)
        if self.mode == LASSO:
            self._points.append(point)
        else:
            self._points = [self._points[0], point]
        self._draw_shape()
        if ev.isFinish():
            modifiers = ev.modifiers()
            points = self._points
            self._clear_shape()
            if self.mode == LASSO and len(points) >= 3:
                polygon = np.array([(p.x(), p.y()) for p in points])
                self.lassoSelected.emit(polygon, modifiers)
            elif self.mode == BOX:
                a, b = points[0], points[-1]
                self.boxSelected.emit(a.x(), a.y(), b.x(), b.y(), modifiers)

    def _draw_shape(self) -> None:
        if self._shape is None:
            return
        path = QPainterPath()
        if self.mode == LASSO:
            path.moveTo(self._points[0])
            for p in self._points[1:]:
                path.lineTo(p)
            path.closeSubpath()
        else:
            a, b = self._points[0], self._points[-1]
            path.addRect(
                min(a.x(), b.x()), min(a.y(), b.y()), abs(b.x() - a.x()), abs(b.y() - a.y())
            )
        self._shape.setPath(path)

    def _clear_shape(self) -> None:
        if self._shape is not None:
            self.removeItem(self._shape)
            self._shape = None
        self._points = []
