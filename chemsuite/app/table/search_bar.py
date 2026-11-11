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

In its normal mode the box understands field syntax
(:func:`chemsuite.core.table.search.parse_terms`): ``ID:CPD-12``, ``pIC50>7``,
``"exact phrase"``, ``-Series:B``; Whole Cell and Regular Expression take the text as typed. A
count of matching rows sits beside the box, and the clock button lists the last searches (saved
with the document).
"""

from __future__ import annotations

import re

import numpy as np
from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QIcon, QKeyEvent
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QToolButton,
    QWidget,
)

from chemsuite.analysis.structure_search import QueryError, StructureQuery
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import REACTION, STRUCTURE
from chemsuite.core.table.filters import OidFilter, SearchFilter
from chemsuite.core.table.search import Query, parse_terms

FILTER_KEY = "search"
DELAY_MS = 200
#: structure searches wait a little longer: each keystroke would start a background job
STRUCTURE_DELAY_MS = 500
ALL_COLUMNS = "All Columns"
INVALID_STYLE = "border: 1px solid #d33"
HISTORY_LIMIT = 10
#: a search kept this long (or confirmed with Return) goes into the history
HISTORY_AFTER_MS = 1500
SYNTAX_TIP = (
    'Words must all be found. "exact phrase" · Column:text · pIC50>7 (also >=, <, <=, =, !=) · '
    '-word or -Column:text excludes. Quote column names with spaces: "Mol Wt"<400.'
)


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
        self.edit.setToolTip(SYNTAX_TIP)
        self.history: list[str] = []
        self._search_icon = self.edit.addAction(QIcon(), QLineEdit.ActionPosition.LeadingPosition)
        self._search_icon.setToolTip(SYNTAX_TIP)
        self.history_button = QToolButton()
        self.history_button.setAutoRaise(True)
        self.history_button.setToolTip("Recent searches")
        self.history_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.history_menu = QMenu(self)
        self.history_menu.aboutToShow.connect(self._fill_history)
        self.history_button.setMenu(self.history_menu)
        self.count = QLabel()
        self.count.setEnabled(False)  # the muted colour of a hint
        self.count.hide()
        self.case = QAction("Match Case", self, checkable=True)
        self.whole = QAction("Whole Cell", self, checkable=True)
        self.regex = QAction("Regular Expression", self, checkable=True)
        self.case.setToolTip("Match Case: capitals and lower case differ")
        self.whole.setToolTip("Whole Cell: the cell is exactly the text")
        self.regex.setToolTip("Regular Expression: the text is a Python regular expression")
        #: the three toggles (one widget, so structure searches can turn them all off)
        self.options_button = QWidget()
        toggles = QHBoxLayout(self.options_button)
        toggles.setContentsMargins(0, 0, 0, 0)
        toggles.setSpacing(1)
        for action in (self.case, self.whole, self.regex):
            button = QToolButton()
            button.setDefaultAction(action)
            button.setAutoRaise(True)
            toggles.addWidget(button)
            action.toggled.connect(lambda _on: self.apply_now())
        self.whole.toggled.connect(lambda on: on and self.regex.setChecked(False))
        self.regex.toggled.connect(lambda on: on and self.whole.setChecked(False))

        self.close_button = QToolButton()
        self.close_button.setAutoRaise(True)
        self.close_button.setToolTip("Close the search bar and clear the search (Esc)")
        self.close_button.clicked.connect(self.dismiss)
        self.edit.escapeWhenEmpty.connect(self.dismiss)
        self.setObjectName("SearchBar")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._paint_glyphs()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 5, 8, 5)
        layout.setSpacing(6)
        layout.addWidget(self.column)
        layout.addWidget(self.history_button)
        layout.addWidget(self.edit, 1)
        layout.addWidget(self.count)
        layout.addWidget(self.options_button)
        layout.addSpacing(4)
        layout.addWidget(self.close_button)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self.apply_now)
        self.edit.textChanged.connect(lambda _t: self._timer.start(self._delay()))
        self._remember_timer = QTimer(self)
        self._remember_timer.setSingleShot(True)
        self._remember_timer.setInterval(HISTORY_AFTER_MS)
        self._remember_timer.timeout.connect(self.remember)
        self.edit.textChanged.connect(lambda _t: self._remember_timer.start())
        self.edit.returnPressed.connect(self.remember)
        self._count_timer = QTimer(self)
        self._count_timer.setSingleShot(True)
        self._count_timer.setInterval(DELAY_MS)
        self._count_timer.timeout.connect(self.update_count)
        self.fill_columns()
        self.column.currentIndexChanged.connect(lambda _i: self._column_changed())
        self._unsubscribe = table.subscribe(self._on_table_event)

    def _paint_glyphs(self) -> None:
        """Glyphs and the flat toolbar look, in the palette's colours (again on theme changes)."""
        from ..shell.glyphs import flat_buttons, glyph

        palette = self.palette()
        ink = palette.text().color()
        muted = palette.placeholderText().color()
        self._search_icon.setIcon(glyph("search", muted))
        self.history_button.setIcon(glyph("history", ink))
        self.close_button.setIcon(glyph("close", ink))
        for action, name in ((self.case, "case"), (self.whole, "whole"), (self.regex, "regex")):
            action.setIcon(glyph(name, ink))
        line = palette.mid().color().name()
        self.setStyleSheet(
            f"QWidget#SearchBar {{ border-bottom: 1px solid {line}; }} " + flat_buttons(palette)
        )

    def changeEvent(self, event) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange:
            self._paint_glyphs()  # follow light and dark themes

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
        self.edit.setToolTip("" if structure else SYNTAX_TIP)
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

    def search_filter(self) -> SearchFilter | None:
        """The filter for the box's text; raises :class:`re.error` for a bad expression."""
        query = self.query()
        if not query.text.strip():
            return None
        if query.mode != "contains":  # Whole Cell and Regular Expression: the text as typed
            query.pattern()
            return SearchFilter(query, column=self.column_name())
        terms = parse_terms(query.text, self.table.column_names)
        if not terms:
            return None
        return SearchFilter(query, column=self.column_name(), terms=terms)

    def structure_query(self) -> StructureQuery:
        return StructureQuery(self.edit.text().strip(), "substructure")

    def apply_now(self) -> None:
        self._timer.stop()
        self._cancel()
        if self.structure_mode():
            self._search_structures()
            return
        self._highlight(None)
        problem = ""
        flt = None
        try:
            flt = self.search_filter()
        except re.error as err:
            problem = f"Invalid regular expression: {err}"
        self._show_problem(problem)
        self._set_filter(flt)
        self.update_count()

    def _search_structures(self) -> None:
        query = self.structure_query()
        column = self.column_name()
        if not query.text:
            self._show_problem("")
            self._highlight(None)
            self._set_filter(None)
            self.update_count()
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
        self.update_count()

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
        self.edit.setToolTip(problem or ("" if self.structure_mode() else SYNTAX_TIP))

    # -- the count ----------------------------------------------------------------------------
    def matching(self) -> int | None:
        """Rows the search keeps (before the other filters), or None with no search."""
        mask = self.table.filter_mask(FILTER_KEY)
        return None if mask is None else int(np.count_nonzero(mask))

    def update_count(self) -> None:
        self._count_timer.stop()
        found = self.matching()
        self.count.setVisible(found is not None)
        if found is not None:
            self.count.setText(f"{found:,} {'row' if found == 1 else 'rows'}")
            shown = len(self.table.view())
            self.count.setToolTip(
                f"{found:,} of {len(self.table):,} rows match; {shown:,} are shown with "
                "the other filters"
            )

    # -- history ------------------------------------------------------------------------------
    def remember(self, text: str | None = None) -> None:
        """Put *text* (default: the box's text) at the top of the recent searches."""
        self._remember_timer.stop()
        text = (self.edit.text() if text is None else text).strip()
        if not text:
            return
        self.history = [text, *[item for item in self.history if item != text]][:HISTORY_LIMIT]

    def _fill_history(self) -> None:
        menu = self.history_menu
        menu.clear()
        for text in self.history:
            action = menu.addAction(text)
            action.triggered.connect(lambda _c=False, t=text: self.edit.setText(t))
        if not self.history:
            menu.addAction("No Recent Searches").setEnabled(False)
            return
        menu.addSeparator()
        menu.addAction("Clear Recent Searches").triggered.connect(self.history.clear)

    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        if FILTER_KEY in self.table.filters:
            self._count_timer.start()  # edits change what the search matches
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

    def focus(self, column: str | None = None) -> None:
        """Open the bar and put the cursor in it; ``column`` also chooses where to look."""
        if column is not None:
            self.set_column(column)
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
        self._remember_timer.stop()
        self._count_timer.stop()
        self._cancel()
        self._unsubscribe()
