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

The bar starts hidden. Edit ▸ Search Rows (or View ▸ Search Bar) opens it. Escape clears the
search, and a second Escape (or the close button) hides the bar. Hiding the bar always clears its
search, so it never leaves rows hidden by a search you cannot see.
"""

from __future__ import annotations

import re

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QAction, QKeyEvent
from PySide6.QtWidgets import QHBoxLayout, QLineEdit, QMenu, QToolButton, QWidget

from chemsuite.core.table import Table
from chemsuite.core.table.filters import SearchFilter
from chemsuite.core.table.search import Query

FILTER_KEY = "search"
DELAY_MS = 200


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

    def __init__(self, table: Table, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.setVisible(False)
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

        close = QToolButton()
        close.setText("✕")
        close.setAutoRaise(True)
        close.setToolTip("Close the search bar and clear the search (Esc)")
        close.clicked.connect(self.dismiss)
        self.edit.escapeWhenEmpty.connect(self.dismiss)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.edit, 1)
        layout.addWidget(button)
        layout.addWidget(close)

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
