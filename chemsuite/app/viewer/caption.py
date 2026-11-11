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
"""The Structure Viewer's caption: chosen column values drawn over the structure, the same
in 2D and 3D (and on exported pictures)."""

from __future__ import annotations

import re

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPalette
from PySide6.QtWidgets import QWidget

from chemsuite.core.table import Table
from chemsuite.core.table.columns import REACTION, STRUCTURE

MARGIN = 8.0  # from the view's corner
PAD = 6.0  # inside the plate
_ID = re.compile(r"(^|[\s_\-])(id|name)($|[\s_\-])", re.IGNORECASE)


def caption_columns(table: Table) -> list[str]:
    """The columns a caption can show: everything but structures and reactions."""
    return [n for n in table.column_names if table.kind(n) not in (STRUCTURE, REACTION)]


def default_caption(table: Table) -> list[str]:
    """The first ID- or name-like column, or nothing."""
    return next(([n] for n in caption_columns(table) if _ID.search(n)), [])


def paint_caption(
    painter: QPainter, lines: list[tuple[str, str]], width: float, font: QFont, palette: QPalette
) -> None:
    """*lines* (name, value) in the top-left corner of a *width*-wide view, on a plate."""
    if not lines:
        return
    metrics = QFontMetrics(font)
    room = max(width - 2 * (MARGIN + PAD), 40.0)
    texts = [
        metrics.elidedText(f"{name}: {value}", Qt.TextElideMode.ElideRight, int(room))
        for name, value in lines
    ]
    line = metrics.height()
    plate = QRectF(
        MARGIN,
        MARGIN,
        max(metrics.horizontalAdvance(t) for t in texts) + 2 * PAD,
        line * len(texts) + 2 * PAD - 2,
    )
    fill = QColor(palette.base().color())
    fill.setAlphaF(0.82)
    border = QColor(palette.mid().color())
    border.setAlphaF(0.6)
    painter.save()
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.setPen(border)
    painter.setBrush(fill)
    painter.drawRoundedRect(plate, 5, 5)
    painter.setFont(font)
    painter.setPen(palette.text().color())
    for i, text in enumerate(texts):
        box = QRectF(plate.left() + PAD, plate.top() + PAD - 1 + i * line, room, line)
        painter.drawText(box, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, text)
    painter.restore()


class CaptionOverlay(QWidget):
    """Draws the caption over whatever is under it; clicks go through."""

    def __init__(self) -> None:
        super().__init__()
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self.lines: list[tuple[str, str]] = []

    def set_lines(self, lines: list[tuple[str, str]]) -> None:
        if lines != self.lines:
            self.lines = list(lines)
            self.update()

    def paintEvent(self, _event) -> None:
        if not self.lines:
            return
        painter = QPainter(self)
        paint_caption(painter, self.lines, self.width(), self.font(), self.palette())
        painter.end()
