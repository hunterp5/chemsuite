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

"""The quick search box above the table: shows only rows with a matching cell."""

from __future__ import annotations

import re

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QKeyEvent
from PySide6.QtWidgets import QHBoxLayout, QLineEdit, QMenu, QToolButton, QWidget

from chemsuite.core.table import Table
from chemsuite.core.table.filters import SearchFilter
from chemsuite.core.table.search import Query

FILTER_KEY = "search"
DELAY_MS = 200


class _SearchEdit(QLineEdit):
    def keyPressEvent(self, event: QKeyEvent) -> None:
        if event.key() == Qt.Key.Key_Escape and self.text():
            self.clear()
            return
        super().keyPressEvent(event)


class SearchBar(QWidget):
    def __init__(self, table: Table, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.edit = _SearchEdit()
        self.edit.setPlaceholderText("Search all columns")
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
        button = QToolButton()
        button.setText("Options")
        button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        button.setMenu(options)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.edit, 1)
        layout.addWidget(button)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(DELAY_MS)
        self._timer.timeout.connect(self.apply_now)
        self.edit.textChanged.connect(lambda _t: self._timer.start())

    def query(self) -> Query:
        mode = (
            "regex" if self.regex.isChecked() else "whole" if self.whole.isChecked() else "contains"
        )
        return Query(self.edit.text(), mode, self.case.isChecked())

    def apply_now(self) -> None:
        self._timer.stop()
        query = self.query()
        valid = True
        try:
            query.pattern()
        except re.error as err:
            valid = False
            self.edit.setToolTip(f"Invalid regular expression: {err}")
        self.edit.setStyleSheet("" if valid else "border: 1px solid #d33")
        if valid:
            self.edit.setToolTip("")
        flt = SearchFilter(query) if query.text and valid else None
        if flt is not None or FILTER_KEY in self.table.filters:
            self.table.set_filter(FILTER_KEY, flt)

    def focus(self) -> None:
        self.edit.setFocus(Qt.FocusReason.ShortcutFocusReason)
        self.edit.selectAll()
