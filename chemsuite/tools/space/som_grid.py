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
"""The self-organizing map's cells, drawn as squares or hexagons.

Hexagons are pointy-topped, odd rows shifted half a cell right and rows √3/2 of a cell apart,
the layout ``analysis.som`` trains on, so neighbors on screen are neighbors on the map.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from PySide6.QtCore import QEvent, QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen, QPolygonF
from PySide6.QtWidgets import QSizePolicy, QToolTip, QWidget

from chemsuite.analysis.som import HEX_PITCH

EMPTY = QColor(235, 235, 235)
MARGIN = 4.0


@dataclass
class Cell:
    text: str = ""
    color: QColor | None = None  # None: empty
    oids: list[int] = field(default_factory=list)
    tip: str = ""


class SomGrid(QWidget):
    """*rows* × *columns* cells; a click on one emits :attr:`cellClicked` (row, column)."""

    cellClicked = Signal(int, int)

    def __init__(self) -> None:
        super().__init__()
        self._rows = self._columns = 0
        self._cells: list[Cell] = []
        self.hexagonal = False
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumSize(120, 120)

    # -- contents ------------------------------------------------------------------------------
    def set_cells(self, rows: int, columns: int, cells: list[Cell], hexagonal: bool) -> None:
        if len(cells) != rows * columns:
            raise ValueError("one cell per grid position")
        self._rows, self._columns, self._cells = rows, columns, cells
        self.hexagonal = hexagonal
        self.update()

    def clear(self) -> None:
        self.set_cells(0, 0, [], self.hexagonal)

    def rowCount(self) -> int:
        return self._rows

    def columnCount(self) -> int:
        return self._columns

    def cell(self, row: int, column: int) -> Cell | None:
        if 0 <= row < self._rows and 0 <= column < self._columns:
            return self._cells[row * self._columns + column]
        return None

    def text(self, row: int, column: int) -> str:
        cell = self.cell(row, column)
        return cell.text if cell is not None else ""

    # -- geometry --------------------------------------------------------------------------
    def _layout(self) -> tuple[float, float, float]:
        """``(cell size, left, top)``: a cell's width in pixels and where the grid starts."""
        width, height = self.width() - 2 * MARGIN, self.height() - 2 * MARGIN
        if not self._rows or not self._columns:
            return 0.0, MARGIN, MARGIN
        if self.hexagonal:  # in cell widths: columns + ½ across, rows·pitch + (1/√3 − pitch)
            across = self._columns + (0.5 if self._rows > 1 else 0.0)
            down = (self._rows - 1) * HEX_PITCH + 2 / math.sqrt(3)
        else:
            across, down = float(self._columns), float(self._rows)
        size = max(1.0, min(width / across, height / down))
        left = MARGIN + (width - size * across) / 2
        top = MARGIN + (height - size * down) / 2
        return size, left, top

    def center(self, row: int, column: int) -> QPointF:
        size, left, top = self._layout()
        if self.hexagonal:
            x = left + size * (column + 0.5 + 0.5 * (row % 2))
            y = top + size * (row * HEX_PITCH + 1 / math.sqrt(3))
        else:
            x, y = left + size * (column + 0.5), top + size * (row + 0.5)
        return QPointF(x, y)

    def shape(self, row: int, column: int) -> QPolygonF:
        size, _left, _top = self._layout()
        c = self.center(row, column)
        if not self.hexagonal:
            h = size / 2
            corners = [(-h, -h), (h, -h), (h, h), (-h, h)]
        else:
            r = size / math.sqrt(3)  # corner radius of a hexagon one cell wide
            corners = [(r * math.cos(a), r * math.sin(a))
                       for a in (math.radians(30 + 60 * k) for k in range(6))]  # fmt: skip
        return QPolygonF([QPointF(c.x() + dx, c.y() + dy) for dx, dy in corners])

    def cell_at(self, point: QPointF) -> tuple[int, int] | None:
        """The cell under *point*, or None."""
        for row in range(self._rows):
            for column in range(self._columns):
                if self.shape(row, column).containsPoint(point, Qt.FillRule.OddEvenFill):
                    return row, column
        return None

    # -- drawing and events ----------------------------------------------------------------
    def sizeHint(self) -> QSize:
        return QSize(360, 300)

    def paintEvent(self, _event) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        size, _left, _top = self._layout()
        font = painter.font()
        font.setPixelSize(max(7, min(13, int(size * 0.3))))
        painter.setFont(font)
        outline = QPen(self.palette().mid().color(), 1)
        for row in range(self._rows):
            for column in range(self._columns):
                cell = self._cells[row * self._columns + column]
                path = QPainterPath()
                path.addPolygon(self.shape(row, column))
                path.closeSubpath()
                color = cell.color if cell.color is not None else EMPTY
                painter.fillPath(path, color)
                painter.setPen(outline)
                painter.drawPath(path)
                if cell.text and size >= 14:
                    painter.setPen(
                        QColor(0, 0, 0) if color.lightness() > 120 else QColor(255, 255, 255)
                    )
                    c = self.center(row, column)
                    box = QRectF(c.x() - size / 2, c.y() - size / 2, size, size)
                    painter.drawText(box, Qt.AlignmentFlag.AlignCenter, cell.text)
        painter.end()

    def mousePressEvent(self, event) -> None:
        found = self.cell_at(event.position())
        if found is not None and event.button() == Qt.MouseButton.LeftButton:
            self.cellClicked.emit(*found)

    def event(self, event) -> bool:
        if event.type() == QEvent.Type.ToolTip:
            found = self.cell_at(QPointF(event.pos()))
            cell = self.cell(*found) if found is not None else None
            if cell is not None and cell.tip:
                QToolTip.showText(event.globalPos(), cell.tip, self)
            else:
                QToolTip.hideText()
            return True
        return super().event(event)
