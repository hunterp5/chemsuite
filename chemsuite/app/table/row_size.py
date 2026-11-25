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

"""How tall the table's rows are when it draws structures, and the ways to change it.

Every row has the same height (Qt keeps a table of uniform rows cheap at any size), so each
gesture sets them all:

- dragging any row divider on the row numbers (a double-click puts the Settings height back);
- Ctrl+scroll wheel (⌘ on macOS), or Larger / Smaller Structures, which zoom: the rows and the
  structure and reaction columns grow or shrink together, so the drawings grow both ways.

The height starts at the Settings value; once the user changes it, it belongs to the document
(saved with its layout) until Settings sets a new value.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import QEvent, QObject, Qt
from PySide6.QtGui import QMouseEvent, QWheelEvent
from PySide6.QtWidgets import QHeaderView

from chemsuite.core.settings import DEFAULT_ROW_HEIGHT, MAX_ROW_HEIGHT, MIN_ROW_HEIGHT
from chemsuite.core.table.columns import REACTION, STRUCTURE

if TYPE_CHECKING:
    from .table_view import TableView

#: how much one zoom step (a wheel notch, Ctrl+= or Ctrl+-) scales the drawings
ZOOM_STEP = 1.15
#: how close to a row divider, in pixels, the pointer grabs it
GRIP = 4
_NOTCH = 120  # a wheel notch in angleDelta units


def clamp(height: float) -> int:
    return max(MIN_ROW_HEIGHT, min(MAX_ROW_HEIGHT, round(height)))


class RowSizer(QObject):
    """The structure row height of one :class:`TableView`, and the gestures that change it."""

    def __init__(self, view: TableView) -> None:
        super().__init__(view)
        self.view = view
        self.height = DEFAULT_ROW_HEIGHT
        #: set by the user (a drag, a zoom or the saved document) rather than Settings
        self.custom = False
        self._default = DEFAULT_ROW_HEIGHT
        # a zoom scales from where it started, so stepping in and out comes back exactly
        self._zoom: tuple[int, dict[int, int]] | None = None
        self._zoom_level = 0
        self._zooming = False
        self._wheel = 0
        # a divider drag: the pointer's start, the height at the start and how many rows, from
        # the top of the header down to the dragged divider, grow with it
        self._drag: tuple[float, int, int] | None = None
        view.horizontalHeader().sectionResized.connect(self._column_resized)
        self.watch(view.verticalHeader())

    # -- the height ---------------------------------------------------------------------------
    def drawn(self) -> bool:
        """Whether the table draws structures, so its rows are structure rows."""
        table = self.view.table_model.table
        return any(table.kind(name) in (STRUCTURE, REACTION) for name in table.column_names)

    def set_default(self, height: int) -> bool:
        """Settings' height. A new value replaces the document's; True when the rows change."""
        height = clamp(height)
        if height == self._default:
            return False
        self._default = height
        changed = height != self.height
        self.height, self.custom, self._zoom = height, False, None
        return changed

    def set_height(self, height: float) -> None:
        """Every row at *height* (clamped), chosen by the user."""
        self._zoom = None
        self._apply(clamp(height))

    def reset(self) -> None:
        """Back to the Settings height."""
        self._zoom = None
        self._apply(self._default)
        self.custom = False

    def restore(self, saved: object) -> None:
        """The height a document was saved with (ignored unless a whole number in range)."""
        whole = isinstance(saved, int) and not isinstance(saved, bool)
        if whole and MIN_ROW_HEIGHT <= saved <= MAX_ROW_HEIGHT:
            self.set_height(saved)

    def saved(self) -> int | None:
        """What the document keeps: the user's height, or None for the Settings one."""
        return self.height if self.custom else None

    def fit_header(self, text_height: int) -> None:
        """Size the view's rows: the structure height, or *text_height* without drawings.
        Quietly, as Qt signals each row's resize, which is slow on a large table; the view
        updates its geometry and overlays itself (``TableView.fit_rows``)."""
        header = self.view.verticalHeader()
        header.blockSignals(True)
        header.setDefaultSectionSize(self.height if self.drawn() else text_height)
        header.blockSignals(False)
        header.viewport().update()
        self.view.viewport().update()

    def _apply(self, height: int) -> None:
        self.custom = True
        if height != self.height:
            self.height = height
            self.view.fit_rows()

    # -- zoom ---------------------------------------------------------------------------------
    def zoom(self, steps: int) -> None:
        """Scale the rows and the drawing columns by *steps* zoom steps (negative: smaller)."""
        if not steps or not self.drawn():
            return
        view = self.view
        if self._zoom is None:
            model = view.table_model
            table = model.table
            widths = {
                col: view.columnWidth(col)
                for col in range(model.columnCount())
                if table.kind(model.column_name(col)) in (STRUCTURE, REACTION)
            }
            self._zoom, self._zoom_level = (self.height, widths), 0
        base_height, widths = self._zoom
        level = self._zoom_level + steps
        height = clamp(base_height * ZOOM_STEP**level)
        if height == self.height:
            return  # at a limit: the columns stay with the rows
        self._zoom_level = level
        self._zooming = True
        try:
            self._apply(height)
            scale = height / base_height
            for col, width in widths.items():
                if col < view.table_model.columnCount():
                    view.setColumnWidth(col, max(20, round(width * scale)))
        finally:
            self._zooming = False

    def wheel(self, event: QWheelEvent) -> bool:
        """Ctrl+wheel zooms (True: the event is used). Without Ctrl the table scrolls."""
        if not event.modifiers() & Qt.KeyboardModifier.ControlModifier or not self.drawn():
            return False
        self._wheel += event.angleDelta().y()
        steps = int(self._wheel / _NOTCH)  # toward zero: a trackpad's small deltas add up
        if steps:
            self._wheel -= steps * _NOTCH
            self.zoom(steps)
        event.accept()
        return True

    def _column_resized(self, *_args) -> None:
        if not self._zooming:
            self._zoom = None  # a column the user sized: the next zoom starts from it

    # -- dragging a row divider ---------------------------------------------------------------
    def watch(self, header: QHeaderView) -> None:
        """Let the row dividers of *header* (the table's row numbers, or an overlay's) be
        dragged."""
        header.setMouseTracking(True)
        header.viewport().setMouseTracking(True)
        header.viewport().installEventFilter(self)

    def divider_at(self, header: QHeaderView, y: float) -> int:
        """The row whose bottom edge is within :data:`GRIP` of *y* (-1: none)."""
        for probe in (y - GRIP, y + GRIP):
            row = header.logicalIndexAt(round(probe))
            if row < 0:
                continue
            bottom = header.sectionViewportPosition(row) + header.sectionSize(row)
            if abs(bottom - y) <= GRIP:
                return row
        return -1

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        header = watched.parent()
        if not isinstance(header, QHeaderView) or not isinstance(event, QMouseEvent):
            return False
        kind = event.type()
        y = event.position().y()
        left = event.button() == Qt.MouseButton.LeftButton
        if self._drag is not None:
            if kind == QEvent.Type.MouseMove:
                start_y, start, rows = self._drag
                self.set_height(start + (y - start_y) / rows)  # the divider follows the pointer
                return True
            if kind == QEvent.Type.MouseButtonRelease and left:
                self._drag = None
                return True
            return False
        row = self.divider_at(header, y) if self.drawn() else -1
        if kind == QEvent.Type.MouseMove and not event.buttons():
            if row >= 0:
                header.viewport().setCursor(Qt.CursorShape.SplitVCursor)
            else:
                header.viewport().unsetCursor()
            return False
        if row < 0 or not left:
            return False
        if kind == QEvent.Type.MouseButtonDblClick:
            self.reset()
            return True
        if kind == QEvent.Type.MouseButtonPress:
            rows = max(1, row - header.logicalIndexAt(0) + 1)
            self._drag = (y, self.height, rows)
            return True
        return False
