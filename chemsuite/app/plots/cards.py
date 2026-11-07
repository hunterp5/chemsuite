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
"""Point cards: the structure and values of a plotted row, shown on hover or kept beside
selected points.

:class:`CardOptions` says when cards show (``WHEN``) and what they hold. A hover card follows
the mouse; with ``when`` "selected" or "both", each selected point (up to ``MAX_PINNED``) gets a
card pinned beside it that moves with the plot as you pan and zoom, until the point is
deselected.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QSize, Qt
from PySide6.QtGui import QImage, QPixmap
from PySide6.QtWidgets import QFrame, QLabel, QVBoxLayout, QWidget

#: when cards show: key -> menu title
WHEN = {
    "hover": "On Hover",
    "both": "On Hover and Selected",
    "selected": "Selected Points Only",
    "off": "Off",
}
#: structure drawing sizes: key -> (title, size in pixels)
SIZES = {
    "small": ("Small", QSize(150, 110)),
    "medium": ("Medium", QSize(220, 160)),
    "large": ("Large", QSize(320, 230)),
}
#: at most this many pinned cards; more selected points show none (they would cover the plot)
MAX_PINNED = 12
GAP = 8  # pixels between a pinned card and its point, and between cards

ImageProvider = Callable[[int, QSize], QImage | None]


@dataclass(frozen=True)
class CardOptions:
    when: str = "hover"
    structure: bool = True
    label: bool = True
    values: bool = True  # the plotted columns (x, y, colour, size)
    columns: tuple[str, ...] = ()  # more columns to show
    size: str = "medium"

    @property
    def hovering(self) -> bool:
        return self.when in ("hover", "both")

    @property
    def pinning(self) -> bool:
        return self.when in ("selected", "both")

    def state(self) -> dict:
        return {
            "when": self.when,
            "structure": self.structure,
            "label": self.label,
            "values": self.values,
            "columns": list(self.columns),
            "size": self.size,
        }

    @classmethod
    def from_state(cls, state: object) -> CardOptions:
        """Saved options; unknown or damaged values fall back to the defaults."""
        if not isinstance(state, dict):
            return cls()
        default = cls()
        columns = state.get("columns", ())
        return replace(
            default,
            when=state.get("when") if state.get("when") in WHEN else default.when,
            structure=bool(state.get("structure", default.structure)),
            label=bool(state.get("label", default.label)),
            values=bool(state.get("values", default.values)),
            columns=tuple(c for c in columns if isinstance(c, str))
            if isinstance(columns, list | tuple)
            else (),
            size=state.get("size") if state.get("size") in SIZES else default.size,
        )


@dataclass
class CardContent:
    title: str = ""
    lines: list[str] = field(default_factory=list)
    image: QImage | None = None

    def text(self) -> str:
        return "\n".join(([self.title] if self.title else []) + self.lines)

    def empty(self) -> bool:
        return not self.title and not self.lines and self.image is None


class CardBody(QFrame):
    """The look shared by hover and pinned cards: a drawing over lines of text."""

    def __init__(self, parent: QWidget | None, flags=Qt.WindowType.Widget) -> None:
        super().__init__(parent, flags)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setObjectName("PointCard")
        self.setStyleSheet(
            "QFrame#PointCard { background: palette(base); border: 1px solid palette(mid);"
            " border-radius: 4px; }"
        )
        self.image = QLabel()
        self.text = QLabel()
        self.text.setTextFormat(Qt.TextFormat.PlainText)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.image, 0, Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.text)

    def fill(self, content: CardContent) -> None:
        self.image.setVisible(content.image is not None)
        if content.image is not None:
            self.image.setPixmap(QPixmap.fromImage(content.image))
        text = content.text()
        self.text.setVisible(bool(text))
        self.text.setText(text)
        self.adjustSize()


class HoverCard(CardBody):
    """A tooltip-like card for the point under the mouse."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.ToolTip)

    def show_at(self, pos: QPoint, content: CardContent) -> None:
        self.fill(content)
        self.move(pos + QPoint(16, 16))
        self.show()


class PinnedCards:
    """Cards kept beside selected points, inside the plot (they pan and zoom with it)."""

    def __init__(self, host: QWidget, to_widget: Callable[[float, float], QPointF]) -> None:
        self.host = host  # the plot widget
        self._to_widget = to_widget  # data coordinates -> host pixels
        self.cards: dict[int, tuple[CardBody, float, float]] = {}

    def show(self, points: Sequence[tuple[int, float, float]], content: Callable) -> None:
        """Pin a card at each ``(oid, x, y)``; *content(oid)* fills it. Others are removed."""
        if len(points) > MAX_PINNED:
            points = []
        wanted = {oid for oid, _x, _y in points}
        for oid in [o for o in self.cards if o not in wanted]:
            self.cards.pop(oid)[0].deleteLater()
        for oid, x, y in points:
            card = self.cards[oid][0] if oid in self.cards else self._new()
            card.fill(content(oid))
            self.cards[oid] = (card, x, y)
        self.place()

    def _new(self) -> CardBody:
        card = CardBody(self.host)
        # clicks and drags go through to the plot underneath
        card.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        return card

    def place(self) -> None:
        """Move each card beside its point (after a pan, zoom or resize), never on another."""
        width, height = self.host.width(), self.host.height()
        taken: list[QRect] = []
        for card, x, y in sorted(self.cards.values(), key=lambda c: (c[1], c[2])):
            if not (np.isfinite(x) and np.isfinite(y)):
                card.hide()
                continue
            at = self._to_widget(x, y)
            point = QPoint(round(at.x()), round(at.y()))
            if not (0 <= point.x() <= width and 0 <= point.y() <= height):
                card.hide()  # the point is out of view
                continue
            rect = _free_spot(point, card.size(), QRect(0, 0, width, height), taken)
            taken.append(rect)
            card.move(rect.topLeft())
            card.show()
            card.raise_()

    def clear(self) -> None:
        for card, _x, _y in self.cards.values():
            card.deleteLater()
        self.cards = {}


def _free_spot(point: QPoint, size: QSize, bounds: QRect, taken: list[QRect]) -> QRect:
    """Where a card for *point* goes: beside it (above right, above left, below right, below
    left, in that order), inside *bounds* and clear of the *taken* cards. When every corner
    is blocked the card slides down from the first one until it is clear."""
    w, h = size.width(), size.height()
    corners = [
        QPoint(point.x() + GAP, point.y() - GAP - h),
        QPoint(point.x() - GAP - w, point.y() - GAP - h),
        QPoint(point.x() + GAP, point.y() + GAP),
        QPoint(point.x() - GAP - w, point.y() + GAP),
    ]

    def fit(corner: QPoint) -> QRect:
        x = min(max(corner.x(), bounds.left()), max(bounds.left(), bounds.right() - w))
        y = min(max(corner.y(), bounds.top()), max(bounds.top(), bounds.bottom() - h))
        return QRect(x, y, w, h)

    def clear(rect: QRect) -> bool:
        padded = rect.adjusted(-GAP // 2, -GAP // 2, GAP // 2, GAP // 2)
        return not any(padded.intersects(other) for other in taken)

    for corner in corners:
        rect = fit(corner)
        if clear(rect):
            return rect
    rect = fit(corners[0])
    for _ in range(200):  # slide down, then wrap to the top, until clear
        rect.moveTop(rect.top() + GAP)
        if rect.bottom() > bounds.bottom():
            rect.moveTop(bounds.top())
        if clear(rect):
            return rect
    return fit(corners[0])
