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

"""Find and Replace: a non-modal window that walks the table's matching cells."""

from __future__ import annotations

import re

import numpy as np
from PySide6.QtCore import QModelIndex, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.table import Changed
from chemsuite.core.table.search import Matches, Query, find, replace_all, replace_cell

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

        self.find_edit = QLineEdit()
        self.replace_edit = QLineEdit()
        self.mode = QComboBox()
        for key, title in MODES.items():
            self.mode.addItem(title, key)
        self.case = QCheckBox("Match case")
        self.scope = QComboBox()
        self.scope.addItem("All columns", "all")
        self.scope.addItem("Selected columns", "selected")
        self.hidden = QCheckBox("Replace All includes hidden rows")
        self.status = QLabel(" ")

        form = QFormLayout()
        form.addRow("Find:", self.find_edit)
        form.addRow("Replace with:", self.replace_edit)
        form.addRow("Match:", self.mode)
        form.addRow("In:", self.scope)
        options = QHBoxLayout()
        options.addWidget(self.case)
        options.addWidget(self.hidden)
        options.addStretch(1)

        self.previous_button = QPushButton("Previous")
        self.next_button = QPushButton("Next")
        self.replace_button = QPushButton("Replace")
        self.replace_all_button = QPushButton("Replace All")
        self.next_button.setDefault(True)
        buttons = QHBoxLayout()
        for button in (
            self.previous_button,
            self.next_button,
            self.replace_button,
            self.replace_all_button,
        ):
            buttons.addWidget(button)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addLayout(options)
        layout.addWidget(self.status)
        layout.addLayout(buttons)

        for signal in (
            self.find_edit.textChanged,
            self.mode.currentIndexChanged,
            self.case.toggled,
            self.scope.currentIndexChanged,
        ):
            signal.connect(self._invalidate)
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
        if self.scope.currentData() != "selected":
            return None
        if self._scope_columns is None:  # freeze the scope while walking the matches
            model = self.view.table_model
            selected = {model.column_name(i.column()) for i in self.view.selectedIndexes()}
            self._scope_columns = [n for n in self.table.column_names if n in selected]
        return self._scope_columns

    def matches(self) -> Matches | None:
        if self._matches is None:
            query = self.query()
            if query is None:
                return None
            self._matches = find(self.table, query, self._columns())
        return self._matches

    def _invalidate(self, *_args) -> None:
        self._matches = None
        self._scope_columns = None
        self._say(" ")

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
    def replace_one(self) -> None:
        query = self.query()
        matches = self.matches()
        if query is None or matches is None:
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
        if query is None:
            return 0
        changed = replace_all(
            self.table,
            query,
            self.replace_edit.text(),
            self._columns(),
            shown_only=not self.hidden.isChecked(),
        )
        self._say(f"Replaced in {changed:,} cell(s)" if changed else "Nothing to replace")
        return changed
