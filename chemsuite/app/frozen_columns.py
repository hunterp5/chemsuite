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

"""Frozen columns: the first columns of the table stay put while the rest scroll sideways.

``QTableView`` cannot do this itself, so the frozen columns are drawn by a second view laid
over the table's left edge (the arrangement of Qt's "frozen column" example). The two views
share one model, one selection model and one cell delegate, so a cell is drawn, selected and
highlighted the same in either, and they keep in step:

- scrolling up or down in either moves both;
- row heights and the widths of the frozen columns follow the table (and a frozen column
  resized on the overlay's header resizes the table's too);
- a cell or header right-click, a structure double-click and typing into a cell go to the
  table's own handlers, so the frozen cells behave like any other.

The overlay never takes keyboard focus: clicking a frozen cell makes it current in the shared
selection model and the table keeps the keyboard, so arrow keys and copy work as before.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import QEvent, QModelIndex, Qt
from PySide6.QtWidgets import QAbstractItemView, QFrame, QHeaderView, QTableView

if TYPE_CHECKING:
    from .table_view import TableView


class FrozenView(QTableView):
    """The overlay that draws the first *count* columns of *main*."""

    def __init__(self, main: TableView) -> None:
        super().__init__(main)
        self.main = main
        self.count = 0
        self.setModel(main.model())
        self.setSelectionModel(main.selectionModel())
        self.setItemDelegate(main.itemDelegate())
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setSelectionMode(main.selectionMode())
        self.setSelectionBehavior(main.selectionBehavior())
        self.setEditTriggers(main.editTriggers())
        self.setWordWrap(False)
        self.setAlternatingRowColors(True)
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollMode(main.verticalScrollMode())
        self.verticalHeader().hide()
        self.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        header = self.horizontalHeader()
        header.setSectionsClickable(True)
        header.setResizeContentsPrecision(200)
        header.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        header.customContextMenuRequested.connect(self._header_menu)
        header.sectionResized.connect(self._overlay_resized)
        header.viewport().installEventFilter(self)
        # a hairline marks where the frozen part ends
        self.setStyleSheet("QTableView { border-right: 1px solid palette(mid); }")

        main.horizontalHeader().sectionResized.connect(self._main_resized)
        main.verticalHeader().sectionResized.connect(
            lambda row, _old, new: self.setRowHeight(row, new)
        )
        main.verticalScrollBar().valueChanged.connect(self.verticalScrollBar().setValue)
        self.verticalScrollBar().valueChanged.connect(main.verticalScrollBar().setValue)
        model = main.model()
        for signal in (model.modelReset, model.columnsInserted, model.columnsRemoved):
            signal.connect(self.refresh)
        model.layoutChanged.connect(self.refresh)

    # -- layout -------------------------------------------------------------------------------
    def set_count(self, count: int) -> None:
        self.count = max(0, count)
        self.refresh()

    def refresh(self, *_args) -> None:
        """Show the first :attr:`count` columns at the table's widths and row heights."""
        main = self.main
        columns = main.model().columnCount()
        for col in range(columns):
            self.setColumnHidden(col, col >= self.count)
            if col < self.count:
                self.setColumnWidth(col, main.columnWidth(col))
        vertical = self.verticalHeader()
        vertical.setDefaultSectionSize(main.verticalHeader().defaultSectionSize())
        self.horizontalHeader().setFixedHeight(main.horizontalHeader().height())
        self.verticalScrollBar().setValue(main.verticalScrollBar().value())
        self.place()

    def frozen_width(self) -> int:
        main = self.main
        columns = min(self.count, main.model().columnCount())
        return sum(main.columnWidth(col) for col in range(columns))

    def place(self) -> None:
        """Lie over the table's left edge, from the header down to the bottom of its cells."""
        main = self.main
        frame = main.frameWidth()
        width = self.frozen_width()
        self.setVisible(self.count > 0 and width > 0)
        self.setGeometry(
            main.verticalHeader().width() + frame,
            frame,
            width + 1,  # the hairline
            main.viewport().height() + main.horizontalHeader().height(),
        )
        self.raise_()

    def _main_resized(self, col: int, _old: int, new: int) -> None:
        if col < self.count and self.columnWidth(col) != new:
            self.setColumnWidth(col, new)
        self.place()

    def _overlay_resized(self, col: int, _old: int, new: int) -> None:
        if col < self.count and self.main.columnWidth(col) != new:
            self.main.setColumnWidth(col, new)

    # -- passing events to the table ----------------------------------------------------------
    def eventFilter(self, watched, event) -> bool:
        if (
            watched is self.horizontalHeader().viewport()
            and event.type() == QEvent.Type.MouseButtonPress
        ):
            self.main.start_column_selection()  # a header click selects whole columns
        return super().eventFilter(watched, event)

    def _header_menu(self, pos) -> None:
        column = self.horizontalHeader().logicalIndexAt(pos)
        if column >= 0:
            self.main.headerMenuRequested.emit(column, self.horizontalHeader().mapToGlobal(pos))

    def mousePressEvent(self, event) -> None:
        self.main.start_cell_selection()
        super().mousePressEvent(event)

    def contextMenuEvent(self, event) -> None:
        index = self.indexAt(event.pos())
        if index.isValid():
            self.main.cellMenuRequested.emit(index, event.globalPos())

    def wheelEvent(self, event) -> None:
        self.main.wheelEvent(event)  # scroll the table (and so this view) as one

    def edit(self, index: QModelIndex, trigger=None, event=None):  # type: ignore[override]
        """Structures go to the table's sketcher handling; other cells edit in place here."""
        if trigger is None:
            return super().edit(index)
        from .table_view import DRAWN, KIND_ROLE

        if index.isValid() and index.data(KIND_ROLE) in DRAWN:
            return self.main.edit(index, trigger, event)
        return super().edit(index, trigger, event)

    def scrollTo(self, index: QModelIndex, hint=QAbstractItemView.ScrollHint.EnsureVisible):
        if index.isValid():
            self.main.scrollTo(index, hint)  # the table scrolls; this view follows it
