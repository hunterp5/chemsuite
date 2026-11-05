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
"""A layout that wraps its widgets onto more lines when there is not room for one.

Panels docked beside the table use it for their control rows. A plain row is as wide as all its
controls together, and that width becomes the dock's minimum, so the dock could not be made
narrow. Wrapped, the minimum is only the widest single control.
"""

from __future__ import annotations

from PySide6.QtCore import QMargins, QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QLayout, QLayoutItem, QWidget


class FlowLayout(QLayout):
    def __init__(self, parent: QWidget | None = None, spacing: int = 6) -> None:
        super().__init__(parent)
        self._items: list[QLayoutItem] = []
        self.setContentsMargins(0, 0, 0, 0)
        self.setSpacing(spacing)

    # -- QLayout ------------------------------------------------------------------------------
    def addItem(self, item: QLayoutItem) -> None:
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int) -> QLayoutItem | None:
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index: int) -> QLayoutItem | None:
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self) -> Qt.Orientation:
        return Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._arrange(QRect(0, 0, width, 0), apply=False)

    def setGeometry(self, rect: QRect) -> None:
        super().setGeometry(rect)
        self._arrange(rect, apply=True)

    def sizeHint(self) -> QSize:
        """Everything on one line."""
        width = height = 0
        for item in self._visible():
            hint = item.sizeHint()
            width += hint.width() + (self.spacing() if width else 0)
            height = max(height, hint.height())
        return QSize(width, height).grownBy(self.contentsMargins())

    def minimumSize(self) -> QSize:
        """The widest single item: everything else can wrap."""
        size = QSize()
        for item in self._visible():
            size = size.expandedTo(item.minimumSize())
        return size.grownBy(self.contentsMargins())

    # -- placing ------------------------------------------------------------------------------
    def _visible(self) -> list[QLayoutItem]:
        return [item for item in self._items if not item.isEmpty()]

    def _arrange(self, rect: QRect, *, apply: bool) -> int:
        """Lay items out left to right, wrapping; returns the height used."""
        margins: QMargins = self.contentsMargins()
        area = rect.marginsRemoved(margins)
        spacing = self.spacing()
        x, y, line_height = area.x(), area.y(), 0
        line: list[tuple[QLayoutItem, QSize]] = []

        def place_line() -> None:
            # Centre each item vertically on its line, so labels and fields line up.
            if apply:
                left = area.x()
                for item, size in line:
                    top = y + (line_height - size.height()) // 2
                    item.setGeometry(QRect(QPoint(left, top), size))
                    left += size.width() + spacing

        for item in self._visible():
            size = item.sizeHint().boundedTo(QSize(area.width(), 10**6))
            size = size.expandedTo(item.minimumSize())
            if line and x + size.width() > area.right() + 1:
                place_line()
                x = area.x()
                y += line_height + spacing
                line, line_height = [], 0
            line.append((item, size))
            x += size.width() + spacing
            line_height = max(line_height, size.height())
        place_line()
        return y + line_height - rect.y() + margins.bottom()


def field(label: str | QLabel, *widgets: QWidget) -> QWidget:
    """*widgets* with a label in front, kept together on one line of a :class:`FlowLayout`."""
    holder = QWidget()
    row = QHBoxLayout(holder)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(4)
    row.addWidget(label if isinstance(label, QLabel) else QLabel(label))
    for widget in widgets:
        row.addWidget(widget)
    return holder


def compact(*combos: QComboBox, chars: int = 8) -> None:
    """Size combo boxes by a few characters instead of their longest item (often a long column
    name), so a docked panel can be narrow; the popup list still shows whole names."""
    for combo in combos:
        combo.setMinimumContentsLength(chars)
        combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
