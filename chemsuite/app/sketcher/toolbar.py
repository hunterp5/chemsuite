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

"""The sketcher's tool bar: one checkable button per tool, icons drawn to follow the theme.

Every icon shares one style: a 20-unit grid with the drawing kept inside 3..17, a 1.6-unit
round-capped stroke in the button text colour, and one font (the system UI font, semi-bold,
12 px) for every element label, so C and Cl sit at the same size and weight. Buttons are all
the same square size with the same rounded hover and checked backgrounds.
"""

from __future__ import annotations

import math
from collections.abc import Callable

from PySide6.QtCore import QEvent, QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontDatabase,
    QIcon,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QPolygonF,
)
from PySide6.QtWidgets import QButtonGroup, QHBoxLayout, QToolButton, QVBoxLayout, QWidget

from chemsuite.app.sketcher import tools
from chemsuite.sketch.model import Order, Stereo
from chemsuite.sketch.templates import RING_BAR

GRID = 20.0
STROKE = 1.6
BUTTON = 30  # px, every tool button
LABEL_PX = 12
BUTTON_STYLE = (
    "QToolButton { border: 1px solid transparent; border-radius: 6px; padding: 0; }"
    " QToolButton:hover { background: palette(midlight); }"
    " QToolButton:checked { background: palette(highlight); border-color: palette(highlight); }"
)


def _line(p: QPainter, n: int, gap: float = 4.0) -> None:
    for i in range(n):
        off = (i - (n - 1) / 2) * gap
        p.drawLine(QPointF(4, 15 + off), QPointF(16, 5 + off))


def _wedge(p: QPainter) -> None:
    p.setBrush(p.pen().color())
    p.drawPolygon(QPolygonF([QPointF(4, 16), QPointF(15, 3), QPointF(17.5, 6)]))


def _hash(p: QPainter) -> None:
    for i in range(1, 6):
        t = i / 5
        cx, cy = 4 + 12 * t, 16 - 11.5 * t
        h = 2.0 * t
        p.drawLine(QPointF(cx - h, cy - h), QPointF(cx + h, cy + h))


def _wavy(p: QPainter) -> None:
    path = QPainterPath(QPointF(3, 15))
    for i in range(4):
        x = 3 + i * 3.5
        path.quadTo(
            QPointF(x + 1.75, 15 - i * 3 - 3 + (2.5 if i % 2 else -2.5)),
            QPointF(x + 3.5, 15 - (i + 1) * 3),
        )
    p.drawPath(path)


def _chain(p: QPainter) -> None:
    path = QPainterPath(QPointF(2, 13))
    for x, y in ((7, 7), (12, 13), (17, 7)):
        path.lineTo(QPointF(x, y))
    p.drawPath(path)


def _select(p: QPainter) -> None:
    p.setBrush(p.pen().color())
    p.drawPolygon(
        QPolygonF([QPointF(5, 3), QPointF(5, 16), QPointF(8.5, 12.5), QPointF(11, 17.5),
                   QPointF(13, 16.5), QPointF(10.5, 11.5), QPointF(15, 11)])
    )  # fmt: skip


def _lasso(p: QPainter) -> None:
    path = QPainterPath(QPointF(6, 15.5))
    path.cubicTo(QPointF(1, 12), QPointF(2, 3.5), QPointF(10, 3.5))
    path.cubicTo(QPointF(18, 3.5), QPointF(19, 11), QPointF(11, 12))
    path.cubicTo(QPointF(8, 12.3), QPointF(6.5, 13.5), QPointF(6, 15.5))
    path.lineTo(QPointF(7.5, 18))
    p.drawPath(path)


def _erase(p: QPainter) -> None:
    p.drawRoundedRect(QRectF(4, 6, 12, 8), 2, 2)
    p.drawLine(QPointF(9, 6), QPointF(9, 14))


def _arrow(p: QPainter) -> None:
    p.drawLine(QPointF(2.5, 10), QPointF(17, 10))
    p.drawLine(QPointF(17, 10), QPointF(13, 7))
    p.drawLine(QPointF(17, 10), QPointF(13, 13))


def _rxn_plus(p: QPainter) -> None:
    p.drawLine(QPointF(10, 4.5), QPointF(10, 15.5))
    p.drawLine(QPointF(4.5, 10), QPointF(15.5, 10))


def _map(p: QPainter) -> None:
    p.drawEllipse(QPointF(5, 13), 2.2, 2.2)
    p.drawEllipse(QPointF(15, 7), 2.2, 2.2)
    pen = QPen(p.pen())
    pen.setStyle(Qt.PenStyle.DashLine)
    pen.setDashPattern([1.5, 1.5])
    p.setPen(pen)
    p.drawLine(QPointF(7, 12), QPointF(13, 8))


def label_font() -> QFont:
    font = QFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.GeneralFont))
    font.setPixelSize(LABEL_PX)
    font.setWeight(QFont.Weight.DemiBold)
    font.setHintingPreference(QFont.HintingPreference.PreferNoHinting)
    return font


def _text(text: str) -> Callable[[QPainter], None]:
    def draw(p: QPainter) -> None:
        p.setFont(label_font())
        p.drawText(QRectF(0, 0, GRID, GRID), Qt.AlignmentFlag.AlignCenter, text)

    return draw


def _charge(sign: int) -> Callable[[QPainter], None]:
    """A circled plus or minus, drawn as strokes like the other tools."""

    def draw(p: QPainter) -> None:
        p.drawEllipse(QPointF(10, 10), 6.5, 6.5)
        p.drawLine(QPointF(7, 10), QPointF(13, 10))
        if sign > 0:
            p.drawLine(QPointF(10, 7), QPointF(10, 13))

    return draw


def icon(draw: Callable[[QPainter], None], colour: QColor) -> QIcon:
    return QIcon(_pixmap(draw, colour))


def _pixmap(draw: Callable[[QPainter], None], colour: QColor) -> QPixmap:
    pixmap = QPixmap(int(GRID * 2), int(GRID * 2))
    pixmap.setDevicePixelRatio(2)
    pixmap.fill(Qt.GlobalColor.transparent)
    p = QPainter(pixmap)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(colour, STROKE)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen)
    draw(p)
    p.end()
    return pixmap


#: (key, tooltip, icon, tool factory). Shortcuts in the tooltips are the hover keys.
TOOLS: list[tuple[str, str, Callable[[QPainter], None], Callable[[], tools.Tool]]] = [
    ("select", "Select and move; drag the handle to rotate (Esc)", _select, tools.SelectTool),
    ("lasso", "Lasso select: draw round atoms to select them", _lasso,
     lambda: tools.SelectTool(lasso=True)),
    ("erase", "Erase (Delete)", _erase, tools.EraseTool),
    ("single", "Single bond (1); click a bond to cycle", lambda p: _line(p, 1),
     lambda: tools.BondTool(Order.SINGLE)),
    ("double", "Double bond (2)", lambda p: _line(p, 2), lambda: tools.BondTool(Order.DOUBLE)),
    ("triple", "Triple bond (3)", lambda p: _line(p, 3, 3.5), lambda: tools.BondTool(Order.TRIPLE)),
    ("wedge", "Wedge: towards the viewer", _wedge,
     lambda: tools.BondTool(Order.SINGLE, Stereo.WEDGE)),
    ("hash", "Hash: away from the viewer", _hash,
     lambda: tools.BondTool(Order.SINGLE, Stereo.HASH)),
    ("wavy", "Wavy: unknown configuration", _wavy,
     lambda: tools.BondTool(Order.SINGLE, Stereo.EITHER)),
    ("chain", "Chain: drag to draw a zig-zag chain", _chain, tools.ChainTool),
    ("plus", "Charge + (=)", _charge(+1), lambda: tools.ChargeTool(+1)),
    ("minus", "Charge − (-)", _charge(-1), lambda: tools.ChargeTool(-1)),
    ("arrow", "Reaction arrow: drag from tail to head; right-click it for other arrows",
     _arrow, tools.ArrowTool),
    ("rxn-plus", "Reaction plus sign", _rxn_plus, tools.PlusTool),
    ("map", "Map atoms: drag from an atom to the same atom on the other side", _map,
     tools.MapTool),
]  # fmt: skip


def _ring(size: int, aromatic: bool) -> Callable[[QPainter], None]:
    def draw(p: QPainter) -> None:
        r = 7.0 if size > 4 else 7.5
        pts = [
            QPointF(10 + r * math.cos(math.pi / 2 + 2 * math.pi * i / size),
                    10 - r * math.sin(math.pi / 2 + 2 * math.pi * i / size))
            for i in range(size)
        ]  # fmt: skip
        p.drawPolygon(QPolygonF(pts))
        if aromatic:
            p.drawEllipse(QPointF(10, 10), r * 0.55, r * 0.55)

    return draw


ATOMS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H")
ATOM_KEYS = {"C": "c", "N": "n", "O": "o", "S": "s", "P": "p", "F": "f", "Cl": "l", "Br": "b",
             "I": "i", "H": "h"}  # fmt: skip


#: R-group and wildcard buttons under the elements (label, tooltip)
QUERY_TOOLS = (
    ("R", "R-group: each new one takes the next number (R1, R2, …)"),
    ("A", "Any atom except hydrogen (A); more queries in the atom's right-click menu"),
)


def tool_key(tool: tools.Tool) -> str:
    """The toolbar key of a tool instance (used to check the matching button)."""
    if isinstance(tool, tools.AtomTool):
        return f"atom:{tool.label}"
    if isinstance(tool, tools.BondTool):
        if tool.stereo is Stereo.WEDGE:
            return "wedge"
        if tool.stereo is Stereo.HASH:
            return "hash"
        if tool.stereo is Stereo.EITHER:
            return "wavy"
        return {Order.SINGLE: "single", Order.DOUBLE: "double", Order.TRIPLE: "triple"}.get(
            tool.order, ""
        )
    if isinstance(tool, tools.TemplateTool):
        return f"template:{tool.template.name}"
    if isinstance(tool, tools.ChargeTool):
        return "plus" if tool.delta > 0 else "minus"
    if isinstance(tool, tools.SelectTool):
        return "lasso" if tool.lasso else "select"
    return {
        tools.EraseTool: "erase",
        tools.ChainTool: "chain",
        tools.ArrowTool: "arrow",
        tools.PlusTool: "rxn-plus",
        tools.MapTool: "map",
    }.get(type(tool), "")


class SketchToolBar(QWidget):
    """The tool buttons, in three strips as in Ketcher: drawing tools (``top``), rings
    (``bottom``) and elements (``side``, a column). One exclusive group spans them all;
    ``picked`` carries a new tool instance."""

    picked = Signal(object)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons: dict[str, QToolButton] = {}
        self._factories: dict[str, Callable[[], tools.Tool]] = {}
        self._draws: dict[str, Callable[[QPainter], None]] = {}
        self.top, self.bottom, self.side = QWidget(), QWidget(), QWidget()
        strips = (QHBoxLayout(self.top), QHBoxLayout(self.bottom), QVBoxLayout(self.side))
        for strip in strips:
            strip.setContentsMargins(0, 0, 0, 0)
            strip.setSpacing(3)
        top, bottom, side = strips
        for i, (key, tip, draw, factory) in enumerate(TOOLS):
            self._add(top, key, tip, draw, factory)
            if i in (2, 5, 8, 9, 11):
                top.addSpacing(10)
        top.addStretch(1)
        for template in RING_BAR:
            self._add(
                bottom,
                f"template:{template.name}",
                f"{template.name} (t cycles the rings)",
                _ring(len(template.atoms), template.aromatic),
                lambda t=template: tools.TemplateTool(t),
            )
        bottom.addStretch(1)
        for symbol in ATOMS:
            self._add(
                side,
                f"atom:{symbol}",
                f"{symbol} ({ATOM_KEYS[symbol]} over an atom)",
                _text(symbol),
                lambda s=symbol: tools.AtomTool(s),
            )
        side.addSpacing(10)
        for symbol, tip in QUERY_TOOLS:
            self._add(
                side, f"atom:{symbol}", tip, _text(symbol), lambda s=symbol: tools.AtomTool(s)
            )
        side.addStretch(1)
        self.buttons["select"].setChecked(True)

    def _add(self, layout, key, tip, draw, factory) -> None:
        button = QToolButton()
        button.setCheckable(True)
        button.setAutoRaise(True)
        button.setFixedSize(BUTTON, BUTTON)
        button.setStyleSheet(BUTTON_STYLE)
        button.setToolTip(tip)
        button.setIconSize(QSize(int(GRID), int(GRID)))
        button.setObjectName(f"sketch-tool-{key}")
        button.clicked.connect(lambda _c=False, k=key: self.picked.emit(self._factories[k]()))
        self.group.addButton(button)
        layout.addWidget(button)
        self.buttons[key] = button
        self._factories[key] = factory
        self._draws[key] = draw
        self._paint_icons_for(button, key)

    def _paint_icons_for(self, button: QToolButton, key: str) -> None:
        palette = self.palette()
        normal = _pixmap(self._draws[key], palette.buttonText().color())
        on = _pixmap(self._draws[key], palette.highlightedText().color())
        both = QIcon()
        for mode in (QIcon.Mode.Normal, QIcon.Mode.Active):
            both.addPixmap(normal, mode, QIcon.State.Off)
            both.addPixmap(on, mode, QIcon.State.On)
        button.setIcon(both)

    def show_tool(self, tool: tools.Tool) -> None:
        button = self.buttons.get(tool_key(tool))
        if button is not None:
            button.setChecked(True)

    def changeEvent(self, event) -> None:
        if event.type() == QEvent.Type.PaletteChange:
            for key, button in self.buttons.items():
                self._paint_icons_for(button, key)
        super().changeEvent(event)
