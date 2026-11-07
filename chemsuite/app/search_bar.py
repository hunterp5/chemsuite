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
"""The quick search box above the table: shows only rows with a matching cell.

A column menu at the left chooses where to look: **All Columns** (every text and number
column) or one column. Choosing a structure column turns the box into a substructure search
(SMILES or SMARTS), run in the background like the structure filter card; matching atoms are
highlighted in the drawings.

The bar starts hidden. Edit ▸ Search Rows (or View ▸ Search Bar) opens it. Escape clears the
search, and a second Escape (or the close button) hides the bar. Hiding the bar always clears its
search, so it never leaves rows hidden by a search you cannot see.
"""

from __future__ import annotations

import re

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QAction, QKeyEvent
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLineEdit, QMenu, QToolButton, QWidget

from chemsuite.analysis.structure_search import QueryError, StructureQuery
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import REACTION, STRUCTURE
from chemsuite.core.table.filters import OidFilter, SearchFilter
from chemsuite.core.table.search import Query

FILTER_KEY = "search"
DELAY_MS = 200
#: structure searches wait a little longer: each keystroke would start a background job
STRUCTURE_DELAY_MS = 500
ALL_COLUMNS = "All Columns"
INVALID_STYLE = "border: 1px solid #d33"


class _SearchEdit(QLineEdit):
    escapeWhenEmpty = Signal()

    def keyPressEvent(self, event: QKeyEvent) -> None:
        if event.key() == Qt.Key.Key_Escape:
            if self.text():
                self.clear()
            else:
                self.escapeWhenEmpty.emit()
            return
        super().keyPressEvent(event)


class SearchBar(QWidget):
    #: the bar was opened (True) or dismissed (False)
    openChanged = Signal(bool)

    def __init__(self, table: Table, parent: QWidget | None = None, searcher=None) -> None:
        super().__init__(parent)
        self.table = table
        #: runs substructure searches (``StructureSearcher``); None: structure columns are not
        #: offered
        self.searcher = searcher
        self.job = None  # the structure search running, if any
        self._highlighted: str | None = None  # the column whose drawings show matches
        self.setVisible(False)
        self.column = QComboBox()
        self.column.setToolTip("Where to search: every text and number column, or one column")
        self.column.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.edit = _SearchEdit()
        self.edit.setClearButtonEnabled(True)
        self.case = QAction("Match Case", self, checkable=True)
        self.whole = QAction("Whole Cell", self, checkable=True)
        self.regex = QAction("Regular Expression", self, checkable=True)
        options = QMenu(self)
        for action in (self.case, self.whole, self.regex):
            options.addAction(action)
            action.toggled.connect(lambda _on: self.apply_now())
        self.whole.toggled.connect(lambda on: on and self.regex.setChecked(False))
        self.regex.toggled.connect(lambda on: on and self.whole.setChecked(False))
        self.options_button = QToolButton()
        self.options_button.setText("Options")
        self.options_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.options_button.setMenu(options)

        close = QToolButton()
        close.setText("✕")
        close.setAutoRaise(True)
        close.setToolTip("Close the search bar and clear the search (Esc)")
        close.clicked.connect(self.dismiss)
        self.edit.escapeWhenEmpty.connect(self.dismiss)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.column)
        layout.addWidget(self.edit, 1)
        layout.addWidget(self.options_button)
        layout.addWidget(close)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self.apply_now)
        self.edit.textChanged.connect(lambda _t: self._timer.start(self._delay()))
        self.fill_columns()
        self.column.currentIndexChanged.connect(lambda _i: self._column_changed())
        self._unsubscribe = table.subscribe(self._on_table_event)

    # -- the column -------------------------------------------------------------------------
    def fill_columns(self) -> None:
        """List the searchable columns again, keeping the choice while it still exists."""
        current = self.column_name()
        self.column.blockSignals(True)
        self.column.clear()
        self.column.addItem(ALL_COLUMNS, None)
        for name in self.table.column_names:
            kind = self.table.kind(name)
            if kind == REACTION or (kind == STRUCTURE and self.searcher is None):
                continue
            self.column.addItem(name, name)
        index = self.column.findData(current) if current is not None else 0
        self.column.setCurrentIndex(max(0, index))
        self.column.blockSignals(False)
        self._show_mode()
        if current is not None and index < 0:  # the column went: search everything instead
            self.apply_now()

    def column_name(self) -> str | None:
        """The chosen column (None: all columns)."""
        return self.column.currentData() if self.column.count() else None

    def set_column(self, name: str | None) -> None:
        index = self.column.findData(name) if name is not None else 0
        self.column.setCurrentIndex(max(0, index))

    def structure_mode(self) -> bool:
        """The chosen column holds structures: the text is a substructure query."""
        name = self.column_name()
        return name is not None and name in self.table and self.table.kind(name) == STRUCTURE

    def _show_mode(self) -> None:
        structure = self.structure_mode()
        self.options_button.setEnabled(not structure)
        name = self.column_name()
        self.edit.setPlaceholderText(
            f"Substructure in {name} (SMILES or SMARTS)"
            if structure
            else f"Search {name}"
            if name
            else "Search all columns"
        )

    def _column_changed(self) -> None:
        self._show_mode()
        self.apply_now()

    def _delay(self) -> int:
        return STRUCTURE_DELAY_MS if self.structure_mode() else DELAY_MS

    # -- searching --------------------------------------------------------------------------
    def query(self) -> Query:
        mode = (
            "regex" if self.regex.isChecked() else "whole" if self.whole.isChecked() else "contains"
        )
        return Query(self.edit.text(), mode, self.case.isChecked())

    def structure_query(self) -> StructureQuery:
        return StructureQuery(self.edit.text().strip(), "substructure")

    def apply_now(self) -> None:
        self._timer.stop()
        self._cancel()
        if self.structure_mode():
            self._search_structures()
            return
        self._highlight(None)
        query = self.query()
        problem = ""
        try:
            query.pattern()
        except re.error as err:
            problem = f"Invalid regular expression: {err}"
        self._show_problem(problem)
        flt = SearchFilter(query, column=self.column_name()) if query.text and not problem else None
        self._set_filter(flt)

    def _search_structures(self) -> None:
        query = self.structure_query()
        column = self.column_name()
        if not query.text:
            self._show_problem("")
            self._highlight(None)
            self._set_filter(None)
            return
        try:
            self.job = self.searcher.start(
                column, query, lambda oids, result: self._structures_found(query, oids, result)
            )
        except QueryError as err:
            self._show_problem(str(err))
            return
        self._show_problem("")

    def _structures_found(self, query: StructureQuery, oids: np.ndarray, result) -> None:
        self.job = None
        if not self.structure_mode() or self.structure_query() != query:
            return  # the search changed while this one ran
        column = self.column_name()
        hits = frozenset(oids[result.hits].tolist())
        self._set_filter(OidFilter(hits, label=f"Substructure in {column}"))
        self._highlight(column)

    def _highlight(self, column: str | None) -> None:
        if self.searcher is None:
            return
        from .structure_search import highlight_for

        if self._highlighted is not None and self._highlighted != column:
            self.searcher.highlight(self._highlighted, None)
        if column is not None:
            self.searcher.highlight(column, highlight_for(self.structure_query()))
        self._highlighted = column

    def _cancel(self) -> None:
        if self.job is not None and self.job.running:
            self.job.cancel()
        self.job = None

    def _set_filter(self, flt) -> None:
        if flt is not None or FILTER_KEY in self.table.filters:
            self.table.set_filter(FILTER_KEY, flt)

    def _show_problem(self, problem: str) -> None:
        self.edit.setStyleSheet(INVALID_STYLE if problem else "")
        self.edit.setToolTip(problem)

    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        if any(
            isinstance(
                change,
                ev.ColumnInserted | ev.ColumnRemoved | ev.ColumnRenamed | ev.ColumnKindChanged,
            )
            for change in event.changes
        ):
            renamed = {c.old: c.new for c in event.changes if isinstance(c, ev.ColumnRenamed)}
            current = self.column_name()
            if current in renamed:
                self.column.setItemData(self.column.currentIndex(), renamed[current])
            self.fill_columns()

    # -- opening and closing ----------------------------------------------------------------
    @property
    def is_open(self) -> bool:
        return not self.isHidden()

    def focus(self) -> None:
        """Open the bar and put the cursor in it."""
        self.open()
        self.edit.setFocus(Qt.FocusReason.ShortcutFocusReason)
        self.edit.selectAll()

    def open(self) -> None:
        if not self.is_open:
            self.setVisible(True)
            self.openChanged.emit(True)

    def dismiss(self) -> None:
        """Clear the search and hide the bar."""
        if self.edit.text():
            self.edit.clear()
        self.apply_now()
        if self.is_open:
            self.setVisible(False)
            self.openChanged.emit(False)

    def close_bar(self) -> None:
        """The document is closing: stop listening to the table."""
        self._timer.stop()
        self._cancel()
        self._unsubscribe()
