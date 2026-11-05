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

"""A two-handle range slider that can draw a histogram of the values behind it."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QKeyEvent, QMouseEvent, QPainter, QPaintEvent
from PySide6.QtWidgets import QSizePolicy, QWidget

HANDLE_RADIUS = 6.0
HISTOGRAM_HEIGHT = 28


class RangeSlider(QWidget):
    """Integer positions ``0..steps``; ``rangeChanged(low, high)`` while dragging or typing."""

    rangeChanged = Signal(int, int)
    released = Signal()

    def __init__(self, steps: int = 1000, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._steps = steps
        self._low = 0
        self._high = steps
        self._active: str | None = None  # "low" | "high" while dragging
        self._last = "high"  # the handle the arrow keys move
        self._counts: np.ndarray | None = None
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)

    # -- values -------------------------------------------------------------------------------
    @property
    def steps(self) -> int:
        return self._steps

    def values(self) -> tuple[int, int]:
        return self._low, self._high

    def set_values(self, low: int, high: int, *, emit: bool = True) -> None:
        low = max(0, min(int(low), self._steps))
        high = max(0, min(int(high), self._steps))
        low, high = min(low, high), max(low, high)
        if (low, high) == (self._low, self._high):
            return
        self._low, self._high = low, high
        self.update()
        if emit:
            self.rangeChanged.emit(low, high)

    def set_histogram(self, counts: np.ndarray | None) -> None:
        """Bar heights across the full range (any number of bins), or None for no histogram."""
        self._counts = None if counts is None or not len(counts) else np.asarray(counts, float)
        self.updateGeometry()
        self.update()

    def sizeHint(self) -> QSize:
        return QSize(160, self._track_y() + int(HANDLE_RADIUS) + 3)

    def minimumSizeHint(self) -> QSize:
        return QSize(60, self.sizeHint().height())

    # -- geometry -----------------------------------------------------------------------------
    def _track_y(self) -> int:
        return (HISTOGRAM_HEIGHT + 4 if self._counts is not None else 0) + int(HANDLE_RADIUS) + 2

    def _left(self) -> float:
        return HANDLE_RADIUS + 1

    def _width(self) -> float:
        return max(1.0, self.width() - 2 * (HANDLE_RADIUS + 1))

    def _x(self, value: int) -> float:
        return self._left() + value / self._steps * self._width()

    def _value_at(self, x: float) -> int:
        t = (x - self._left()) / self._width()
        return round(min(max(t, 0.0), 1.0) * self._steps)

    def _nearest(self, x: float) -> str:
        low, high = self._x(self._low), self._x(self._high)
        if abs(x - low) == abs(x - high):  # handles on top of each other: pick by direction
            return "low" if x < low else "high"
        return "low" if abs(x - low) < abs(x - high) else "high"

    # -- painting -----------------------------------------------------------------------------
    def paintEvent(self, _event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        palette = self.palette()
        highlight = palette.highlight().color()
        muted = palette.mid().color()
        x_low, x_high = self._x(self._low), self._x(self._high)

        if self._counts is not None:
            counts = self._counts
            peak = counts.max() or 1.0
            bin_width = self._width() / len(counts)
            painter.setPen(Qt.PenStyle.NoPen)
            for i, count in enumerate(counts.tolist()):
                if not count:
                    continue
                height = max(1.0, count / peak * HISTOGRAM_HEIGHT)
                left = self._left() + i * bin_width
                inside = left + bin_width > x_low and left < x_high
                color = QColor(highlight if inside else muted)
                color.setAlpha(150 if inside else 90)
                painter.setBrush(color)
                painter.drawRect(
                    QRectF(left + 0.5, HISTOGRAM_HEIGHT - height + 2, bin_width - 1, height)
                )

        y = self._track_y()
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(muted)
        painter.drawRoundedRect(QRectF(self._left(), y - 1.5, self._width(), 3), 1.5, 1.5)
        painter.setBrush(highlight)
        painter.drawRoundedRect(QRectF(x_low, y - 1.5, max(1.0, x_high - x_low), 3), 1.5, 1.5)
        painter.setPen(palette.dark().color())
        painter.setBrush(palette.button())
        for x in (x_low, x_high):
            painter.drawEllipse(QPointF(x, y), HANDLE_RADIUS, HANDLE_RADIUS)
        if self.hasFocus():
            painter.setPen(highlight)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            x = x_low if self._last == "low" else x_high
            painter.drawEllipse(QPointF(x, y), HANDLE_RADIUS + 2, HANDLE_RADIUS + 2)
        painter.end()

    # -- interaction --------------------------------------------------------------------------
    def _move(self, handle: str, value: int) -> None:
        if handle == "low":
            self.set_values(min(value, self._high), self._high)
        else:
            self.set_values(self._low, max(value, self._low))

    def mousePressEvent(self, event: QMouseEvent) -> None:
        if event.button() != Qt.MouseButton.LeftButton:
            return
        x = event.position().x()
        self._active = self._last = self._nearest(x)
        self._move(self._active, self._value_at(x))

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        if self._active is not None:
            self._move(self._active, self._value_at(event.position().x()))

    def mouseReleaseEvent(self, _event: QMouseEvent) -> None:
        if self._active is not None:
            self._active = None
            self.released.emit()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        step = self._steps // 100 if event.modifiers() & Qt.KeyboardModifier.ShiftModifier else 1
        key = event.key()
        if key in (Qt.Key.Key_Left, Qt.Key.Key_Down):
            delta = -step
        elif key in (Qt.Key.Key_Right, Qt.Key.Key_Up):
            delta = step
        else:
            super().keyPressEvent(event)
            return
        current = self._low if self._last == "low" else self._high
        self._move(self._last, current + delta)
        self.released.emit()
