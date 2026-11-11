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

"""Find and Replace: a non-modal window that walks the table's matching cells.

**In** sets where it looks: all columns, the selected columns, the selected cells, or columns
chosen from a list. The selection is read once, when the search starts, so moving the cursor
from match to match does not shrink it. In regular-expression mode the replacement can use the
groups (``$1`` or ``\\1``, ``${name}`` or ``\\g<name>``); a reference to a group the expression
does not have is reported before anything is replaced.
"""

from __future__ import annotations

import re

import numpy as np
from PySide6.QtCore import QModelIndex, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLayout,
    QLineEdit,
    QMenu,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.table import Changed
from chemsuite.core.table.search import (
    Cells,
    Matches,
    Query,
    find,
    replace_all,
    replace_cell,
    searchable_columns,
)

from ..widgets.dialogs import Section, inline
from ..widgets.dialogs.kit import muted
from .table_view import TableView

MODES = {"contains": "Contains", "whole": "Whole cell", "regex": "Regular expression"}


class FindReplaceDialog(QDialog):
    def __init__(self, view: TableView, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Find and Replace")
        self.setWindowFlag(Qt.WindowType.Tool)
        self.setModal(False)
        self.view = view
        self.table = view.table_model.table
        self._matches: Matches | None = None
        self._scope_columns: list[str] | None = None
        self._scope_cells: dict[str, np.ndarray] | None = None
        self.chosen: list[str] = []

        self.find_edit = QLineEdit()
        self.replace_edit = QLineEdit()
        self.mode = QComboBox()
        for key, title in MODES.items():
            self.mode.addItem(title, key)
        self.case = QCheckBox("Match case")
        self.scope = QComboBox()
        self.scope.addItem("All columns", "all")
        self.scope.addItem("Selected columns", "selected")
        self.scope.addItem("Selected cells", "cells")
        self.scope.addItem("Chosen columns", "chosen")
        self.columns_button = QToolButton()
        self.columns_button.setText("Columns…")
        self.columns_button.setToolTip("The columns to search in")
        self.columns_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.columns_menu = QMenu(self)
        self.columns_menu.aboutToShow.connect(self._fill_columns_menu)
        self.columns_button.setMenu(self.columns_menu)
        self.columns_button.setVisible(False)
        self.replace_edit.setToolTip(
            "In regular-expression mode, $1 or \\1 is the first group, ${name} or \\g<name> "
            "a named one, and $$ a dollar sign."
        )
        self.hidden = QCheckBox("Hidden rows too")
        self.case.setToolTip("Upper and lower case must match")
        self.status = QLabel(" ")
        muted(self.status)

        find = Section("Find")
        find.row("Find", self.find_edit)
        find.row("Replace with", self.replace_edit, help=self.replace_edit.toolTip())
        find.row("Match", self.mode)
        find.row(
            "In",
            inline(self.scope, self.columns_button, stretch=0),
            help="Where to look. The selection is read when the search starts, so moving from "
            "match to match does not change it.",
        )
        find.row("Options", inline(self.case, self.hidden))
        find.set_help(self.hidden, "Replace All also changes rows hidden by filters.")

        self.previous_button = QPushButton("Previous")
        self.next_button = QPushButton("Next")
        self.replace_button = QPushButton("Replace")
        self.replace_all_button = QPushButton("Replace All")
        self.next_button.setDefault(True)
        buttons = QHBoxLayout()
        buttons.setSpacing(8)
        buttons.addWidget(self.status, 1)
        for button in (
            self.previous_button,
            self.next_button,
            self.replace_button,
            self.replace_all_button,
        ):
            buttons.addWidget(button)
        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setFrameShadow(QFrame.Shadow.Plain)
        rule.setEnabled(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 14)
        layout.setSpacing(6)
        layout.addWidget(find)
        layout.addSpacing(4)
        layout.addWidget(rule)
        layout.addLayout(buttons)
        layout.setSizeConstraint(QLayout.SizeConstraint.SetFixedSize)

        for signal in (
            self.find_edit.textChanged,
            self.mode.currentIndexChanged,
            self.case.toggled,
        ):
            signal.connect(self._invalidate)
        self.scope.currentIndexChanged.connect(self.rescope)
        self.scope.currentIndexChanged.connect(
            lambda _i: self.columns_button.setVisible(self.scope.currentData() == "chosen")
        )
        self.previous_button.clicked.connect(lambda: self.find_next(backwards=True))
        self.next_button.clicked.connect(lambda: self.find_next())
        self.replace_button.clicked.connect(self.replace_one)
        self.replace_all_button.clicked.connect(self.replace_everything)
        self._unsubscribe = self.table.subscribe(self._on_table_event)

    def close_dialog(self) -> None:
        self._unsubscribe()
        self.close()
        self.deleteLater()

    # -- query --------------------------------------------------------------------------------
    def query(self) -> Query | None:
        query = Query(self.find_edit.text(), self.mode.currentData(), self.case.isChecked())
        if not query.text:
            self._say("")
            return None
        try:
            query.pattern()
        except re.error as err:
            self._say(f"Invalid regular expression: {err}")
            return None
        return query

    def _columns(self) -> list[str] | None:
        scope = self.scope.currentData()
        if scope == "chosen":
            return [name for name in self.table.column_names if name in self.chosen]
        if scope != "selected":
            return None
        if self._scope_columns is None:  # freeze the scope while walking the matches
            model = self.view.table_model
            selected = {model.column_name(i.column()) for i in self.view.selectedIndexes()}
            self._scope_columns = [n for n in self.table.column_names if n in selected]
        return self._scope_columns

    def _cells(self) -> Cells | None:
        """The selected cells, by column, when the scope is the selection (read once)."""
        if self.scope.currentData() != "cells":
            return None
        if self._scope_cells is None:
            self._scope_cells = self.view.selected_cells_by_column()
        return self._scope_cells

    def matches(self) -> Matches | None:
        if self._matches is None:
            query = self.query()
            if query is None:
                return None
            self._matches = find(self.table, query, self._columns(), within=self._cells())
        return self._matches

    def _invalidate(self, *_args) -> None:
        """The query changed. The scope stays: walking the matches moves the selection."""
        self._matches = None
        self._say(" ")

    def rescope(self, *_args) -> None:
        """Read the selection afresh (the In choice changed, or the window was reopened)."""
        self._scope_columns = None
        self._scope_cells = None
        self._invalidate()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self.rescope()

    def set_chosen(self, names: list[str]) -> None:
        """Search in *names* (the Chosen columns scope)."""
        self.chosen = [name for name in names if name in self.table]
        self.scope.setCurrentIndex(self.scope.findData("chosen"))
        self._invalidate()

    def _fill_columns_menu(self) -> None:
        menu = self.columns_menu
        menu.clear()
        for name in searchable_columns(self.table):
            action = menu.addAction(name)
            action.setCheckable(True)
            action.setChecked(name in self.chosen)
            action.toggled.connect(lambda on, n=name: self._toggle_column(n, on))

    def _toggle_column(self, name: str, on: bool) -> None:
        chosen = [n for n in self.chosen if n != name] + ([name] if on else [])
        self.set_chosen(chosen)

    def _on_table_event(self, event: object) -> None:
        if isinstance(event, Changed):
            self._matches = None  # positions or contents moved; recompute on next use

    def _say(self, text: str) -> None:
        self.status.setText(text or " ")

    # -- walking ------------------------------------------------------------------------------
    def _positions(self, matches: Matches) -> tuple[np.ndarray, np.ndarray]:
        """(display row, model column) per match."""
        model = self.view.table_model
        rows = model.rows_of_oids(matches.oids)
        cols = np.array([model.column_index(n) for n in matches.columns], dtype=np.int64)
        return rows, cols[matches.column_index]

    def _current_match(self, matches: Matches) -> int | None:
        current = self.view.currentIndex()
        if not current.isValid() or not len(matches):
            return None
        rows, cols = self._positions(matches)
        hit = np.flatnonzero((rows == current.row()) & (cols == current.column()))
        return int(hit[0]) if len(hit) else None

    def find_next(self, *, backwards: bool = False) -> bool:
        matches = self.matches()
        if matches is None:
            return False
        if not len(matches):
            self._say("No matches")
            return False
        rows, cols = self._positions(matches)
        current = self.view.currentIndex()
        here = (current.row(), current.column()) if current.isValid() else (-1, -1)
        order = rows * (cols.max() + 1) + cols  # display order as one number
        point = here[0] * (cols.max() + 1) + here[1]
        if backwards:
            before = np.flatnonzero(order < point)
            index = int(before[-1]) if len(before) else len(matches) - 1  # wrap around
        else:
            after = np.flatnonzero(order > point)
            index = int(after[0]) if len(after) else 0
        self._go(int(rows[index]), int(cols[index]))
        self._say(f"{index + 1:,} of {len(matches):,}")
        return True

    def _go(self, row: int, column: int) -> None:
        model = self.view.table_model
        target: QModelIndex = model.index(row, column)
        self.view.setCurrentIndex(target)
        self.view.scrollTo(target, QAbstractItemView.ScrollHint.EnsureVisible)

    # -- replacing ----------------------------------------------------------------------------
    def _replacement_ok(self, query: Query) -> bool:
        problem = query.replacement_problem(self.replace_edit.text())
        if problem:
            self._say(problem)
        return not problem

    def replace_one(self) -> None:
        query = self.query()
        matches = self.matches()
        if query is None or matches is None or not self._replacement_ok(query):
            return
        index = self._current_match(matches)
        if index is None:  # not on a match yet: go to the next one first
            self.find_next()
            return
        oid, column = matches.cell(index)
        replace_cell(self.table, oid, column, query, self.replace_edit.text())
        self.find_next()

    def replace_everything(self) -> int:
        query = self.query()
        if query is None or not self._replacement_ok(query):
            return 0
        changed = replace_all(
            self.table,
            query,
            self.replace_edit.text(),
            self._columns(),
            shown_only=not self.hidden.isChecked(),
            within=self._cells(),
        )
        self._say(f"Replaced in {changed:,} cell(s)" if changed else "Nothing to replace")
        return changed
