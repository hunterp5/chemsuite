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

"""Frozen columns and rows: the first columns stay put while the rest scroll sideways, and the
first rows stay put while the rest scroll up and down.

``QTableView`` cannot do this itself, so the frozen parts are drawn by views laid over the
table (the arrangement of Qt's "frozen column" example):

- :class:`FrozenView` over the left edge shows the frozen columns and scrolls up and down with
  the table;
- :class:`FrozenRowsView` over the top shows the frozen rows (with their row numbers) and
  scrolls sideways with the table;
- :class:`FrozenCorner`, where the two meet, shows the frozen rows of the frozen columns and
  never scrolls.

Every overlay shares the table's model, selection model and cell delegate, so a cell is drawn,
selected and highlighted the same wherever it shows, and:

- row heights and column widths follow the table (and a frozen column resized on the column
  overlay's header resizes the table's too);
- a cell, column-header or row-header right-click, a structure double-click and typing into a
  cell go to the table's own handlers, so frozen cells behave like any other.

The overlays never take keyboard focus: clicking a frozen cell makes it current in the shared
selection model and the table keeps the keyboard, so arrow keys and copy work as before.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import QEvent, QModelIndex, Qt
from PySide6.QtWidgets import QAbstractItemView, QFrame, QHeaderView, QTableView

if TYPE_CHECKING:
    from .table_view import TableView


class _Overlay(QTableView):
    """What every frozen overlay shares: the table's model, selection and delegate, and its
    events passed on to the table."""

    #: which parts the overlay scrolls with the table
    follows_vertical = True
    follows_horizontal = False
    #: whether it shows the row numbers (the overlays at the top cover the table's own)
    shows_row_numbers = False

    def __init__(self, main: TableView) -> None:
        super().__init__(main)
        self.main = main
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
        self.setHorizontalScrollMode(main.horizontalScrollMode())
        vertical = self.verticalHeader()
        vertical.setSectionResizeMode(QHeaderView.ResizeMode.Fixed)
        if self.shows_row_numbers:
            source = main.verticalHeader()
            vertical.setStyleSheet(source.styleSheet())
            vertical.setDefaultAlignment(source.defaultAlignment())
            vertical.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            vertical.customContextMenuRequested.connect(self._row_header_menu)
            vertical.viewport().installEventFilter(self)
        else:
            vertical.hide()

        main.horizontalHeader().sectionResized.connect(self._main_resized)
        main.verticalHeader().sectionResized.connect(
            lambda row, _old, new: self.setRowHeight(row, new)
        )
        if self.follows_vertical:
            main.verticalScrollBar().valueChanged.connect(self.verticalScrollBar().setValue)
            self.verticalScrollBar().valueChanged.connect(main.verticalScrollBar().setValue)
        if self.follows_horizontal:
            main.horizontalScrollBar().valueChanged.connect(self.horizontalScrollBar().setValue)
            self.horizontalScrollBar().valueChanged.connect(main.horizontalScrollBar().setValue)
        model = main.model()
        for signal in (model.modelReset, model.columnsInserted, model.columnsRemoved):
            signal.connect(self.refresh)
        model.layoutChanged.connect(self.refresh)

    # -- layout (each overlay places itself) --------------------------------------------------
    def refresh(self, *_args) -> None:
        main = self.main
        for col in range(main.model().columnCount()):
            hidden = self.column_hidden(col)
            self.setColumnHidden(col, hidden)
            if not hidden:
                self.setColumnWidth(col, main.columnWidth(col))
        # rows past the frozen ones are simply clipped by the overlay's height
        self.verticalHeader().setDefaultSectionSize(main.verticalHeader().defaultSectionSize())
        if self.shows_row_numbers:
            self.verticalHeader().setFont(main.verticalHeader().font())
            self.verticalHeader().setFixedWidth(main.verticalHeader().width())
        else:
            self.horizontalHeader().setFixedHeight(main.horizontalHeader().height())
        if self.follows_vertical:
            self.verticalScrollBar().setValue(main.verticalScrollBar().value())
        if self.follows_horizontal:
            self.horizontalScrollBar().setValue(main.horizontalScrollBar().value())
        self.place()

    def column_hidden(self, col: int) -> bool:
        return False

    def place(self) -> None:
        raise NotImplementedError

    def _main_resized(self, col: int, _old: int, new: int) -> None:
        if not self.column_hidden(col) and self.columnWidth(col) != new:
            self.setColumnWidth(col, new)
        self.main.place_frozen()

    # -- passing events to the table ----------------------------------------------------------
    def eventFilter(self, watched, event) -> bool:
        if event.type() == QEvent.Type.MouseButtonPress:
            if watched is self.horizontalHeader().viewport():
                self.main.start_column_selection()  # a header click selects whole columns
            elif watched is self.verticalHeader().viewport():
                self.main.start_cell_selection()
        return super().eventFilter(watched, event)

    def _row_header_menu(self, pos) -> None:
        header = self.verticalHeader()
        row = header.logicalIndexAt(pos)
        if row >= 0:
            self.main.request_row_menu(row, header.mapToGlobal(pos))

    def mousePressEvent(self, event) -> None:
        self.main.start_cell_selection()
        super().mousePressEvent(event)

    def contextMenuEvent(self, event) -> None:
        index = self.indexAt(event.pos())
        if index.isValid():
            self.main.cellMenuRequested.emit(index, event.globalPos())

    def wheelEvent(self, event) -> None:
        self.main.wheelEvent(event)  # scroll the table (and so the overlays) as one

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
            self.main.scrollTo(index, hint)  # the table scrolls; the overlays follow it


class FrozenView(_Overlay):
    """The overlay that draws the first *count* columns of *main* (scrolls up and down)."""

    def __init__(self, main: TableView) -> None:
        self.count = 0
        super().__init__(main)
        header = self.horizontalHeader()
        header.setSectionsClickable(True)
        header.setResizeContentsPrecision(200)
        header.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        header.customContextMenuRequested.connect(self._header_menu)
        header.sectionResized.connect(self._overlay_resized)
        header.viewport().installEventFilter(self)
        # a hairline marks where the frozen part ends
        self.setStyleSheet("QTableView { border-right: 1px solid palette(mid); }")

    def set_count(self, count: int) -> None:
        self.count = max(0, count)
        self.refresh()

    def column_hidden(self, col: int) -> bool:
        return col >= self.count

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

    def _overlay_resized(self, col: int, _old: int, new: int) -> None:
        if col < self.count and self.main.columnWidth(col) != new:
            self.main.setColumnWidth(col, new)

    def _header_menu(self, pos) -> None:
        column = self.horizontalHeader().logicalIndexAt(pos)
        if column >= 0:
            self.main.headerMenuRequested.emit(column, self.horizontalHeader().mapToGlobal(pos))


class FrozenRowsView(_Overlay):
    """The overlay that draws the first *count* rows of *main*, with their row numbers
    (scrolls sideways)."""

    follows_vertical = False
    follows_horizontal = True
    shows_row_numbers = True

    def __init__(self, main: TableView) -> None:
        self.count = 0
        super().__init__(main)
        self.horizontalHeader().hide()
        # a hairline marks where the frozen part ends
        self.setStyleSheet("QTableView { border-bottom: 1px solid palette(mid); }")

    def set_count(self, count: int) -> None:
        self.count = max(0, count)
        self.refresh()

    def frozen_height(self) -> int:
        main = self.main
        rows = min(self.count, main.model().rowCount())
        return sum(main.rowHeight(row) for row in range(rows))

    def place(self) -> None:
        """Lie over the table's top rows, under the column header, row numbers included."""
        main = self.main
        frame = main.frameWidth()
        height = self.frozen_height()
        self.setVisible(self.count > 0 and height > 0)
        self.verticalScrollBar().setValue(0)
        self.setGeometry(
            frame,
            frame + main.horizontalHeader().height(),
            main.verticalHeader().width() + main.viewport().width(),
            height + 1,  # the hairline
        )
        self.raise_()


class FrozenCorner(_Overlay):
    """Where frozen rows meet frozen columns: those cells, never scrolling."""

    follows_vertical = False
    follows_horizontal = False
    shows_row_numbers = True

    def __init__(self, main: TableView) -> None:
        super().__init__(main)
        self.horizontalHeader().hide()
        self.setStyleSheet(
            "QTableView { border-right: 1px solid palette(mid);"
            " border-bottom: 1px solid palette(mid); }"
        )

    def column_hidden(self, col: int) -> bool:
        return col >= self.main.frozen_count

    def place(self) -> None:
        main = self.main
        frame = main.frameWidth()
        width = main._frozen.frozen_width() if main._frozen is not None else 0
        height = main._frozen_rows.frozen_height() if main._frozen_rows is not None else 0
        self.setVisible(main.frozen_count > 0 and main.frozen_row_count > 0)
        self.verticalScrollBar().setValue(0)
        self.horizontalScrollBar().setValue(0)
        self.setGeometry(
            frame,
            frame + main.horizontalHeader().height(),
            main.verticalHeader().width() + width + 1,
            height + 1,
        )
        self.raise_()
