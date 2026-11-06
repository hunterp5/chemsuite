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

"""The sketcher's tool bar: one checkable button per tool, icons drawn to follow the theme."""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import QEvent, QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPen, QPixmap, QPolygonF
from PySide6.QtWidgets import QButtonGroup, QHBoxLayout, QToolButton, QWidget

from chemsuite.app.sketcher import tools
from chemsuite.sketch.model import Order, Stereo

GRID = 20.0


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


def _erase(p: QPainter) -> None:
    p.drawRoundedRect(QRectF(4, 6, 12, 8), 2, 2)
    p.drawLine(QPointF(9, 6), QPointF(9, 14))


def _text(text: str) -> Callable[[QPainter], None]:
    def draw(p: QPainter) -> None:
        font = p.font()
        font.setPixelSize(13 if len(text) < 2 else 11)
        font.setBold(True)
        p.setFont(font)
        p.drawText(QRectF(0, 0, GRID, GRID), Qt.AlignmentFlag.AlignCenter, text)

    return draw


def icon(draw: Callable[[QPainter], None], colour: QColor) -> QIcon:
    pixmap = QPixmap(int(GRID * 2), int(GRID * 2))
    pixmap.setDevicePixelRatio(2)
    pixmap.fill(Qt.GlobalColor.transparent)
    p = QPainter(pixmap)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(colour, 1.5)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    p.setPen(pen)
    draw(p)
    p.end()
    return QIcon(pixmap)


#: (key, tooltip, icon, tool factory). Shortcuts in the tooltips are the hover keys.
TOOLS: list[tuple[str, str, Callable[[QPainter], None], Callable[[], tools.Tool]]] = [
    ("select", "Select and move (Esc)", _select, tools.SelectTool),
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
    ("plus", "Charge + (=)", _text("+"), lambda: tools.ChargeTool(+1)),
    ("minus", "Charge − (-)", _text("−"), lambda: tools.ChargeTool(-1)),
]  # fmt: skip
ATOMS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H")
ATOM_KEYS = {"C": "c", "N": "n", "O": "o", "S": "s", "P": "p", "F": "f", "Cl": "l", "Br": "b",
             "I": "i", "H": "h"}  # fmt: skip


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
    if isinstance(tool, tools.ChargeTool):
        return "plus" if tool.delta > 0 else "minus"
    return {
        tools.SelectTool: "select",
        tools.EraseTool: "erase",
        tools.ChainTool: "chain",
    }.get(type(tool), "")


class SketchToolBar(QWidget):
    """Tool buttons; ``picked`` carries a new tool instance."""

    picked = Signal(object)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons: dict[str, QToolButton] = {}
        self._factories: dict[str, Callable[[], tools.Tool]] = {}
        self._draws: dict[str, Callable[[QPainter], None]] = {}
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        for i, (key, tip, draw, factory) in enumerate(TOOLS):
            self._add(layout, key, tip, draw, factory)
            if i in (1, 7, 8):
                layout.addSpacing(8)
        layout.addSpacing(8)
        for symbol in ATOMS:
            self._add(
                layout,
                f"atom:{symbol}",
                f"{symbol} ({ATOM_KEYS[symbol]} over an atom)",
                _text(symbol),
                lambda s=symbol: tools.AtomTool(s),
            )
        layout.addStretch(1)
        self.buttons["select"].setChecked(True)

    def _add(self, layout, key, tip, draw, factory) -> None:
        button = QToolButton()
        button.setCheckable(True)
        button.setAutoRaise(True)
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
        button.setIcon(icon(self._draws[key], self.palette().buttonText().color()))

    def show_tool(self, tool: tools.Tool) -> None:
        button = self.buttons.get(tool_key(tool))
        if button is not None:
            button.setChecked(True)

    def changeEvent(self, event) -> None:
        if event.type() == QEvent.Type.PaletteChange:
            for key, button in self.buttons.items():
                self._paint_icons_for(button, key)
        super().changeEvent(event)
