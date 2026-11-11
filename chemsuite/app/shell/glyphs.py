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
"""Small line icons (plot toolbar, panel headers), drawn as vectors in the theme's text colour.

Drawn rather than taken from a font, so they look the same on every OS and follow light and
dark themes (the panel redraws them when the palette changes).
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QIcon, QPainter, QPainterPath, QPalette, QPen, QPixmap

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


def _rotate(p: QPainter) -> None:
    """An arrow round a circle: turn the 3D view."""
    path = QPainterPath()
    path.arcMoveTo(QRectF(2.5, 2.5, 11, 11), 60)
    path.arcTo(QRectF(2.5, 2.5, 11, 11), 60, 290)
    tip = path.currentPosition()
    path.moveTo(tip + QPointF(-2.6, -0.6))
    path.lineTo(tip)
    path.lineTo(tip + QPointF(0.4, -2.7))
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


def _gear(p: QPainter) -> None:
    """A cog: settings."""
    import math

    c = GRID / 2
    path = QPainterPath()
    for i in range(8):
        a = i * math.pi / 4
        inner, outer = 4.6, 6.6
        path.moveTo(c + inner * math.cos(a), c + inner * math.sin(a))
        path.lineTo(c + outer * math.cos(a), c + outer * math.sin(a))
    p.drawPath(path)
    p.drawEllipse(QPointF(c, c), 4.4, 4.4)
    p.drawEllipse(QPointF(c, c), 1.8, 1.8)


def _undock(p: QPainter) -> None:
    """A window with an arrow out of its corner: float the panel (or dock it again)."""
    path = QPainterPath(QPointF(7, 3))
    path.lineTo(3, 3)
    path.lineTo(3, 13)
    path.lineTo(13, 13)
    path.lineTo(13, 9)
    path.moveTo(8, 8)
    path.lineTo(13.5, 2.5)
    path.moveTo(9.5, 2.5)
    path.lineTo(13.5, 2.5)
    path.lineTo(13.5, 6.5)
    p.drawPath(path)


def _close(p: QPainter) -> None:
    path = QPainterPath(QPointF(4, 4))
    path.lineTo(12, 12)
    path.moveTo(12, 4)
    path.lineTo(4, 12)
    p.drawPath(path)


def _arrange(p: QPainter) -> None:
    """A pane split in two, the right half split again: arrange panes."""
    p.drawRoundedRect(QRectF(2.5, 3, 11, 10), 1.5, 1.5)
    path = QPainterPath(QPointF(8, 3))
    path.lineTo(8, 13)
    path.moveTo(8, 8)
    path.lineTo(13.5, 8)
    p.drawPath(path)


def _plus(p: QPainter) -> None:
    """A plus: add."""
    path = QPainterPath(QPointF(8, 3))
    path.lineTo(8, 13)
    path.moveTo(3, 8)
    path.lineTo(13, 8)
    p.drawPath(path)


def _sets(p: QPainter) -> None:
    """Three stacked lines with a bookmark: saved sets."""
    path = QPainterPath(QPointF(2.5, 4.5))
    path.lineTo(8.5, 4.5)
    path.moveTo(2.5, 8)
    path.lineTo(8.5, 8)
    path.moveTo(2.5, 11.5)
    path.lineTo(8.5, 11.5)
    path.moveTo(10.5, 3)  # the bookmark
    path.lineTo(13.5, 3)
    path.lineTo(13.5, 13)
    path.lineTo(12, 11.2)
    path.lineTo(10.5, 13)
    path.closeSubpath()
    p.drawPath(path)


def _trash(p: QPainter) -> None:
    """A bin: remove all."""
    path = QPainterPath(QPointF(3, 4.5))
    path.lineTo(13, 4.5)
    path.moveTo(6.5, 4.5)
    path.lineTo(6.5, 2.8)
    path.lineTo(9.5, 2.8)
    path.lineTo(9.5, 4.5)
    path.moveTo(4.3, 4.5)
    path.lineTo(5.1, 13.2)
    path.lineTo(10.9, 13.2)
    path.lineTo(11.7, 4.5)
    path.moveTo(7, 7)
    path.lineTo(7, 11)
    path.moveTo(9, 7)
    path.lineTo(9, 11)
    p.drawPath(path)


def _history(p: QPainter) -> None:
    """A clock whose rim turns back into an arrow: recent items."""
    path = QPainterPath()
    face = QRectF(2.5, 2.5, 11, 11)
    path.arcMoveTo(face, 160)
    path.arcTo(face, 160, -310)
    start = QPointF(2.5 + 5.5 - 5.5 * 0.94, 8 - 5.5 * 0.342)  # the rim point at 160°
    path.moveTo(start + QPointF(-1.6, -1.8))
    path.lineTo(start)
    path.lineTo(start + QPointF(2.3, -0.6))
    path.moveTo(8, 5)
    path.lineTo(8, 8)
    path.lineTo(10.2, 9.4)
    p.drawPath(path)


def _search(p: QPainter) -> None:
    """A magnifying glass."""
    p.drawEllipse(QRectF(2.5, 2.5, 8, 8))
    p.drawLine(QPointF(9.4, 9.4), QPointF(13.5, 13.5))


def _chevron(p: QPainter, flip: bool, bar: bool) -> None:
    """A chevron pointing left (or right when *flip*), with a bar at its tip for first/last."""
    x = (lambda v: GRID - v) if flip else (lambda v: v)
    path = QPainterPath()
    shift = 1.5 if bar else 0.0
    path.moveTo(x(10 + shift), 3.5)
    path.lineTo(x(5.5 + shift), 8)
    path.lineTo(x(10 + shift), 12.5)
    if bar:
        path.moveTo(x(4), 3.5)
        path.lineTo(x(4), 12.5)
    p.drawPath(path)


def _letters(p: QPainter, text: str, y: float = 12.0, size: float = 9.5) -> None:
    """*text* filled in the glyph's colour, centred across the grid."""
    font = QFont()
    font.setPixelSize(round(size))
    font.setWeight(QFont.Weight.DemiBold)
    path = QPainterPath()
    path.addText(0, 0, font, text)
    box = path.boundingRect()
    path.translate((GRID - box.width()) / 2 - box.left(), y)
    p.save()
    p.setBrush(p.pen().color())
    p.setPen(Qt.PenStyle.NoPen)
    p.drawPath(path)
    p.restore()


def _case(p: QPainter) -> None:
    """ "Aa": match case."""
    _letters(p, "Aa")


def _whole(p: QPainter) -> None:
    """ "ab" over a bracket: the whole cell."""
    _letters(p, "ab", y=10.0, size=9)
    path = QPainterPath()
    path.moveTo(2, 10.5)
    path.lineTo(2, 13.5)
    path.lineTo(14, 13.5)
    path.lineTo(14, 10.5)
    p.drawPath(path)


def _regex(p: QPainter) -> None:
    """ ".*": a regular expression."""
    p.save()
    p.setBrush(p.pen().color())
    p.drawEllipse(QPointF(4.0, 12.0), 1.1, 1.1)
    p.restore()
    centre = QPointF(10.5, 6.0)
    path = QPainterPath()
    for dx, dy in ((0, 3.4), (2.95, 1.7), (2.95, -1.7)):
        path.moveTo(centre - QPointF(dx, dy))
        path.lineTo(centre + QPointF(dx, dy))
    p.drawPath(path)


def _not(p: QPainter) -> None:
    """ "¬": logical not (keep what the filter would hide)."""
    path = QPainterPath()
    path.moveTo(3, 6.5)
    path.lineTo(13, 6.5)
    path.lineTo(13, 10.5)
    p.drawPath(path)


def _tick_box(p: QPainter, ticked: bool) -> None:
    p.drawRoundedRect(QRectF(2.5, 2.5, 11, 11), 2.5, 2.5)
    if ticked:
        path = QPainterPath()
        path.moveTo(5.2, 8.2)
        path.lineTo(7.2, 10.2)
        path.lineTo(10.8, 5.8)
        p.drawPath(path)


def _pencil(p: QPainter) -> None:
    """A pencil: draw a structure."""
    path = QPainterPath()
    path.moveTo(3, 13)
    path.lineTo(3.6, 10.2)
    path.lineTo(10.8, 3)
    path.lineTo(13, 5.2)
    path.lineTo(5.8, 12.4)
    path.closeSubpath()
    path.moveTo(9.2, 4.6)
    path.lineTo(11.4, 6.8)
    p.drawPath(path)


def _ruler(p: QPainter) -> None:
    """A ruler at a slant: measure distances and angles."""
    p.save()
    p.translate(8, 8)
    p.rotate(-45)
    p.drawRoundedRect(QRectF(-7, -2.6, 14, 5.2), 1, 1)
    path = QPainterPath()
    for i, x in enumerate((-4.5, -2.0, 0.5, 3.0)):
        path.moveTo(x, -2.6)
        path.lineTo(x, -0.2 if i % 2 else -1.0)
    p.drawPath(path)
    p.restore()


def _overlay(p: QPainter) -> None:
    """Two offset rings: show conformers on top of each other."""
    p.drawEllipse(QPointF(6.2, 8), 3.8, 3.8)
    p.drawEllipse(QPointF(9.8, 8), 3.8, 3.8)


def _play(p: QPainter) -> None:
    """A triangle pointing right: play."""
    path = QPainterPath()
    path.moveTo(5, 3.5)
    path.lineTo(12, 8)
    path.lineTo(5, 12.5)
    path.closeSubpath()
    p.drawPath(path)


def _pause(p: QPainter) -> None:
    """Two upright bars: pause."""
    p.drawLine(QPointF(6, 4), QPointF(6, 12))
    p.drawLine(QPointF(10, 4), QPointF(10, 12))


GLYPHS: dict[str, Draw] = {
    "pencil": _pencil,
    "not": _not,
    "tick_all": lambda p: _tick_box(p, True),
    "tick_none": lambda p: _tick_box(p, False),
    "search": _search,
    "first": lambda p: _chevron(p, False, True),
    "previous": lambda p: _chevron(p, False, False),
    "next": lambda p: _chevron(p, True, False),
    "last": lambda p: _chevron(p, True, True),
    "play": _play,
    "pause": _pause,
    "case": _case,
    "whole": _whole,
    "regex": _regex,
    "history": _history,
    "plus": _plus,
    "sets": _sets,
    "trash": _trash,
    "arrange": _arrange,
    "gear": _gear,
    "undock": _undock,
    "close": _close,
    "pan": _pan,
    "rotate": _rotate,
    "box": _box,
    "lasso": _lasso,
    "fit": _fit,
    "export": _export,
    "ruler": _ruler,
    "overlay": _overlay,
}


def glyph(name: str, colour: QColor) -> QIcon:
    """The icon *name* drawn in *colour*, with a fainter disabled state."""
    icon = QIcon()
    faint = QColor(colour)
    faint.setAlphaF(0.35)
    for mode, ink in ((QIcon.Mode.Normal, colour), (QIcon.Mode.Disabled, faint)):
        icon.addPixmap(_pixmap(GLYPHS[name], ink), mode)
    return icon


def flat_buttons(palette: QPalette, selector: str = "QToolButton") -> str:
    """Style sheet for glyph buttons as the plot toolbars draw them: no frame, a faint wash
    on hover and a tint of the highlight colour when checked."""
    ink = palette.text().color()
    hover = QColor(ink)
    hover.setAlphaF(0.08)
    checked = QColor(palette.highlight().color())
    checked.setAlphaF(0.28)

    def rgba(c: QColor) -> str:
        return f"rgba({c.red()}, {c.green()}, {c.blue()}, {c.alphaF():.2f})"

    return (
        f"{selector} {{ border: none; border-radius: 4px; padding: 2px; }}"
        f" {selector}:hover {{ background: {rgba(hover)}; }}"
        f" {selector}:checked {{ background: {rgba(checked)}; }}"
        f" {selector}::menu-indicator {{ image: none; }}"
    )


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
