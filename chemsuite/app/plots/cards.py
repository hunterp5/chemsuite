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

Cards are compact: the drawing is sized to its molecule (``SIZES`` are the largest boxes), at
the normal bond length so it stays legible. Pinned cards are laid out together: each takes the
nearest free spot around its point, clear of the other cards and of the other selected points,
with a thin leader line back to its point when it had to move away.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace

import numpy as np
from PySide6.QtCore import QLineF, QPoint, QPointF, QRect, QSize, Qt
from PySide6.QtGui import QFontMetrics, QImage, QPainter, QPen, QPixmap
from PySide6.QtWidgets import QFrame, QLabel, QVBoxLayout, QWidget

#: when cards show: key -> menu title
WHEN = {
    "hover": "On Hover",
    "both": "On Hover and Selected",
    "selected": "Selected Points Only",
    "off": "Off",
}
#: largest structure drawing: key -> (title, size in pixels); smaller molecules get less
SIZES = {
    "small": ("Small", QSize(100, 72)),
    "medium": ("Medium", QSize(140, 100)),
    "large": ("Large", QSize(200, 144)),
}
#: at most this many pinned cards; more selected points show none (they would cover the plot)
MAX_PINNED = 12
GAP = 6  # pixels between a pinned card and its point, and between cards
TEXT_WIDTH = 180  # longer lines are cut short with …
POINT_CLEARANCE = 5  # cards keep this far from the other pinned points
#: where a card may go relative to its point, in order of preference
DIRECTIONS = ((1, -1), (-1, -1), (1, 1), (-1, 1), (1, 0), (-1, 0), (0, -1), (0, 1))
RINGS = 8  # how far out (in steps of the card's half size) to look for a free spot
LEADER = 3 * GAP  # a card further than this from its point gets a leader line
#: pinned cards may cover at most this share of the plot; beyond it their drawings shrink
PINNED_SHARE = 0.3
TINY = QSize(72, 52)  # the smallest drawing a pinned card falls back to
TEXT_HEIGHT = 44  # room a card's text lines take, for the estimate

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


def pinned_size(chosen: QSize, count: int, plot: QSize) -> QSize | None:
    """The drawing size for *count* pinned cards on a *plot* of this size: the chosen one when
    they fit, else smaller, so the cards can sit side by side; None when even the smallest
    would crowd the plot (the cards then show text only)."""
    if count == 0:
        return chosen
    budget = PINNED_SHARE * plot.width() * plot.height() / count
    smaller = sorted(
        (size for _title, size in SIZES.values() if size.width() <= chosen.width()),
        key=lambda size: -size.width(),
    )
    for size in [*smaller, TINY]:
        if (size.width() + 2 * GAP) * (size.height() + TEXT_HEIGHT) <= budget:
            return size
    return None


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
        font = self.text.font()
        font.setPointSizeF(max(8.0, font.pointSizeF() - 2))
        self.text.setFont(font)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 3, 4, 3)
        layout.setSpacing(1)
        layout.addWidget(self.image, 0, Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.text)

    def fill(self, content: CardContent) -> None:
        self.image.setVisible(content.image is not None)
        if content.image is not None:
            self.image.setPixmap(QPixmap.fromImage(content.image))
        metrics = QFontMetrics(self.text.font())
        lines = [
            metrics.elidedText(line, Qt.TextElideMode.ElideRight, TEXT_WIDTH)
            for line in content.text().splitlines()
        ]
        self.text.setVisible(bool(lines))
        self.text.setText("\n".join(lines))
        self.adjustSize()


class HoverCard(CardBody):
    """A tooltip-like card for the point under the mouse."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.ToolTip)

    def show_at(self, pos: QPoint, content: CardContent) -> None:
        self.fill(content)
        self.move(pos + QPoint(16, 16))
        self.show()


class Leaders(QWidget):
    """Thin lines from moved pinned cards back to their points (under the cards)."""

    def __init__(self, host: QWidget) -> None:
        super().__init__(host)
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.lines: list[QLineF] = []

    def paintEvent(self, _event) -> None:
        if not self.lines:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(self.palette().mid().color(), 1))
        for line in self.lines:
            painter.drawLine(line)


class PinnedCards:
    """Cards kept beside selected points, inside the plot (they pan and zoom with it)."""

    def __init__(self, host: QWidget, to_widget: Callable[[float, float], QPointF]) -> None:
        self.host = host  # the plot widget
        self._to_widget = to_widget  # data coordinates -> host pixels
        self.cards: dict[int, tuple[CardBody, float, float]] = {}
        self.leaders = Leaders(host)
        self.leaders.hide()

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
        """Lay the cards out beside their points (after a pan, zoom or resize), never on
        another card or another pinned point."""
        width, height = self.host.width(), self.host.height()
        bounds = QRect(0, 0, width, height)
        shown: list[tuple[CardBody, QPoint]] = []
        for card, x, y in self.cards.values():
            at = self._to_widget(x, y) if np.isfinite(x) and np.isfinite(y) else None
            point = QPoint(round(at.x()), round(at.y())) if at is not None else None
            if point is None or not bounds.contains(point):
                card.hide()  # the point is out of view
                continue
            shown.append((card, point))
        anchors = [point for _card, point in shown]
        taken: list[QRect] = []
        lines: list[QLineF] = []
        # top to bottom, left to right: neighbours end up in a tidy order
        for card, point in sorted(shown, key=lambda c: (c[1].y(), c[1].x())):
            others = [p for p in anchors if p != point]
            rect = _free_spot(point, card.size(), bounds, taken, others)
            taken.append(rect)
            card.move(rect.topLeft())
            card.show()
            card.raise_()
            near = _nearest_on(rect, point)
            if (QPointF(near) - QPointF(point)).manhattanLength() > LEADER:
                lines.append(QLineF(QPointF(point), QPointF(near)))
        self.leaders.lines = lines
        self.leaders.setGeometry(bounds)
        self.leaders.setVisible(bool(lines))
        if lines:
            self.leaders.raise_()
            for card, _point in shown:
                card.raise_()
            self.leaders.update()

    def clear(self) -> None:
        for card, _x, _y in self.cards.values():
            card.deleteLater()
        self.cards = {}
        self.leaders.lines = []
        self.leaders.hide()


def _nearest_on(rect: QRect, point: QPoint) -> QPoint:
    return QPoint(
        min(max(point.x(), rect.left()), rect.right()),
        min(max(point.y(), rect.top()), rect.bottom()),
    )


def _free_spot(
    point: QPoint,
    size: QSize,
    bounds: QRect,
    taken: list[QRect],
    points: Sequence[QPoint] = (),
) -> QRect:
    """Where a card for *point* goes: the nearest spot around it (above right first, then the
    other corners and sides, then further out) inside *bounds*, clear of the *taken* cards and
    of the other *points*. When nothing is clear, the spot that overlaps least."""
    w, h = size.width(), size.height()

    def fit(x: int, y: int) -> QRect:
        x = min(max(x, bounds.left()), max(bounds.left(), bounds.right() - w + 1))
        y = min(max(y, bounds.top()), max(bounds.top(), bounds.bottom() - h + 1))
        return QRect(x, y, w, h)

    def candidate(dx: int, dy: int, reach: int) -> QRect:
        x = point.x() + GAP + reach if dx > 0 else point.x() - GAP - reach - w
        y = point.y() + GAP + reach if dy > 0 else point.y() - GAP - reach - h
        if dx == 0:
            x = point.x() - w // 2
        if dy == 0:
            y = point.y() - h // 2
        return fit(x, y)

    def overlap(rect: QRect) -> int:
        padded = rect.adjusted(-GAP // 2, -GAP // 2, GAP // 2, GAP // 2)
        area = 0
        for other in taken:
            common = padded.intersected(other)
            if not common.isEmpty():
                area += common.width() * common.height()
        c = POINT_CLEARANCE
        for p in [point, *points]:
            if rect.intersects(QRect(p.x() - c, p.y() - c, 2 * c, 2 * c)):
                area += 4 * c * c * (50 if p == point else 4)  # never cover its own point
        return area

    step = max(GAP, min(w, h) // 2)
    best: tuple[int, QRect] | None = None
    for ring in range(RINGS):
        for dx, dy in DIRECTIONS:
            rect = candidate(dx, dy, ring * step)
            cost = overlap(rect)
            if cost == 0:
                return rect
            if best is None or cost < best[0]:
                best = (cost, rect)
    return best[1] if best is not None else fit(point.x() + GAP, point.y() - GAP - h)
