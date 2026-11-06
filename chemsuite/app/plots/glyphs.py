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
"""Small line icons for the plot toolbar, drawn as vectors in the theme's text colour.

Drawn rather than taken from a font, so they look the same on every OS and follow light and
dark themes (the panel redraws them when the palette changes).
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPen, QPixmap

#: drawn on a 16 × 16 grid; pixmaps are made at 2× for sharp high-DPI screens
GRID = 16.0
SCALE = 2

Draw = Callable[[QPainter], None]


def _pan(p: QPainter) -> None:
    """Four arrows from the centre: pan and zoom."""
    path = QPainterPath()
    c = GRID / 2
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        tip = QPointF(c + dx * 6.5, c + dy * 6.5)
        path.moveTo(c + dx * 1.5, c + dy * 1.5)
        path.lineTo(tip)
        side = QPointF(dy * 2.2, dx * 2.2)
        back = QPointF(tip.x() - dx * 2.2, tip.y() - dy * 2.2)
        path.moveTo(back + side)
        path.lineTo(tip)
        path.lineTo(back - side)
    p.drawPath(path)


def _box(p: QPainter) -> None:
    """A dashed selection rectangle."""
    pen = QPen(p.pen())
    pen.setDashPattern([2.0, 1.6])
    p.setPen(pen)
    p.drawRect(QRectF(2.5, 3.5, 11, 9))


def _lasso(p: QPainter) -> None:
    """A free-drawn loop with a tail."""
    path = QPainterPath(QPointF(5.0, 12.5))
    path.cubicTo(1.0, 11.0, 1.0, 3.0, 8.0, 2.8)
    path.cubicTo(14.5, 2.6, 15.0, 9.5, 9.5, 10.8)
    path.cubicTo(7.0, 11.4, 5.5, 11.0, 5.0, 12.5)
    path.cubicTo(4.6, 13.8, 5.6, 14.8, 7.2, 14.6)
    pen = QPen(p.pen())
    pen.setDashPattern([2.0, 1.4])
    p.setPen(pen)
    p.drawPath(path)


def _fit(p: QPainter) -> None:
    """Four corners: fit the view to the data."""
    path = QPainterPath()
    lo, hi, arm = 2.5, 13.5, 4.0
    for x, y, dx, dy in ((lo, lo, 1, 1), (hi, lo, -1, 1), (lo, hi, 1, -1), (hi, hi, -1, -1)):
        path.moveTo(x + dx * arm, y)
        path.lineTo(x, y)
        path.lineTo(x, y + dy * arm)
    p.drawPath(path)
    p.drawEllipse(QPointF(8, 8), 1.4, 1.4)


def _export(p: QPainter) -> None:
    """An arrow down into a tray: save the plot as an image."""
    path = QPainterPath(QPointF(8, 2))
    path.lineTo(8, 10)
    path.moveTo(4.8, 7)
    path.lineTo(8, 10.2)
    path.lineTo(11.2, 7)
    path.moveTo(2.5, 10.5)
    path.lineTo(2.5, 13.5)
    path.lineTo(13.5, 13.5)
    path.lineTo(13.5, 10.5)
    p.drawPath(path)


GLYPHS: dict[str, Draw] = {
    "pan": _pan,
    "box": _box,
    "lasso": _lasso,
    "fit": _fit,
    "export": _export,
}


def glyph(name: str, colour: QColor) -> QIcon:
    """The icon *name* drawn in *colour*, with a fainter disabled state."""
    icon = QIcon()
    faint = QColor(colour)
    faint.setAlphaF(0.35)
    for mode, ink in ((QIcon.Mode.Normal, colour), (QIcon.Mode.Disabled, faint)):
        icon.addPixmap(_pixmap(GLYPHS[name], ink), mode)
    return icon


def _pixmap(draw: Draw, colour: QColor) -> QPixmap:
    size = int(GRID * SCALE)
    pixmap = QPixmap(size, size)
    pixmap.setDevicePixelRatio(SCALE)
    pixmap.fill(Qt.GlobalColor.transparent)
    p = QPainter(pixmap)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(colour, 1.3)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen)
    p.setBrush(Qt.BrushStyle.NoBrush)
    draw(p)
    p.end()
    return pixmap
