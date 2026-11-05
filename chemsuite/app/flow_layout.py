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
"""A row of widgets that wraps onto more lines when there is not room for one.

Panels docked beside the table use it for their control rows. A plain row is as wide as all its
controls together, and that width becomes the dock's minimum, so the dock could not be made
narrow. Wrapped, the minimum is only the widest single control.

It is a widget that places its children itself, not a ``QLayout`` subclass: Qt calls a layout's
pure virtual methods while it destroys it, and a Python-implemented layout can be called after
its Python half is gone, which aborts the process at exit.
"""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QSizePolicy, QWidget

#: child events that can change how the row should be laid out
_RELAYOUT_EVENTS = (QEvent.Type.Show, QEvent.Type.Hide, QEvent.Type.LayoutRequest)


class FlowRow(QWidget):
    def __init__(self, parent: QWidget | None = None, spacing: int = 6) -> None:
        super().__init__(parent)
        self._widgets: list[QWidget] = []
        self._spacing = spacing
        policy = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred)
        policy.setHeightForWidth(True)
        self.setSizePolicy(policy)

    def addWidget(self, widget: QWidget) -> None:
        hidden = widget.isHidden() and widget.testAttribute(
            Qt.WidgetAttribute.WA_WState_ExplicitShowHide
        )
        widget.setParent(self)
        widget.installEventFilter(self)
        self._widgets.append(widget)
        if not hidden:
            widget.show()
        self._relayout()

    def widgets(self) -> list[QWidget]:
        return list(self._widgets)

    # -- size ---------------------------------------------------------------------------------
    def _visible(self) -> list[QWidget]:
        return [w for w in self._widgets if not w.isHidden()]

    def sizeHint(self) -> QSize:
        """Everything on one line."""
        width = height = 0
        for widget in self._visible():
            hint = widget.sizeHint()
            width += hint.width() + (self._spacing if width else 0)
            height = max(height, hint.height())
        return QSize(width, height)

    def minimumSizeHint(self) -> QSize:
        """The widest single widget: everything else can wrap."""
        size = QSize(0, 0)
        for widget in self._visible():
            size = size.expandedTo(widget.minimumSizeHint())
        return QSize(size.width(), self.heightForWidth(max(size.width(), 1)))

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._arrange(width, apply=False)

    # -- placing ------------------------------------------------------------------------------
    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._arrange(self.width(), apply=True)

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        if event.type() in _RELAYOUT_EVENTS and watched in self._widgets:
            self._relayout()
        return False

    def event(self, event: QEvent) -> bool:
        if event.type() == QEvent.Type.LayoutRequest:  # a child's size hint changed
            self._relayout()
        return super().event(event)

    def _relayout(self) -> None:
        self.updateGeometry()
        self._arrange(self.width(), apply=True)

    def _arrange(self, width: int, *, apply: bool) -> int:
        """Lay widgets out left to right, wrapping at *width*; returns the height used."""
        spacing = self._spacing
        x = y = line_height = 0
        line: list[tuple[QWidget, QSize]] = []

        def place_line() -> None:
            # Centre each widget vertically on its line, so labels and fields line up.
            if apply:
                left = 0
                for widget, size in line:
                    top = y + (line_height - size.height()) // 2
                    widget.setGeometry(QRect(QPoint(left, top), size))
                    left += size.width() + spacing

        for widget in self._visible():
            size = widget.sizeHint().boundedTo(QSize(max(width, 1), 10**6))
            size = size.expandedTo(widget.minimumSizeHint())
            if line and x + size.width() > width:
                place_line()
                x = 0
                y += line_height + spacing
                line, line_height = [], 0
            line.append((widget, size))
            x += size.width() + spacing
            line_height = max(line_height, size.height())
        place_line()
        return y + line_height


def field(label: str | QLabel, *widgets: QWidget) -> QWidget:
    """*widgets* with a label in front, kept together on one line of a :class:`FlowRow`."""
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
