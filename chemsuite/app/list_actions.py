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

"""Data ▸ Lists: save selections under a name, bring them back, and manage the names."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QInputDialog,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.core.lists import ListError, RowLists
from chemsuite.core.selection import Selection
from chemsuite.core.table import Table
from chemsuite.core.table.filters import OidFilter

FILTER_KEY = "list"
ORIGIN = "lists"


class ListActions:
    """The commands behind the Lists menu, for one workspace."""

    def __init__(
        self,
        parent: QWidget,
        table: Table,
        selection: Selection,
        lists: RowLists,
        show_message: Callable[[str], None],
    ) -> None:
        self.parent = parent
        self.table = table
        self.selection = selection
        self.lists = lists
        self.show_message = show_message

    # -- commands -----------------------------------------------------------------------------
    def add_selected(self) -> str | None:
        """Ask for a list name (new or existing) and add the selected rows to it."""
        oids = self.selection.as_array()
        if not len(oids):
            self.show_message("Select some rows first")
            return None
        names = self.lists.names
        name, ok = QInputDialog.getItem(
            self.parent,
            "Add to List",
            f"Add {len(oids):,} selected row(s) to list (pick one or type a new name):",
            names,
            0,
            True,
        )
        if not ok:
            return None
        try:
            added = self.lists.add(name, oids)
        except ListError as err:
            QMessageBox.warning(self.parent, "Add to List", str(err))
            return None
        name = name.strip()
        self.show_message(f"Added {added:,} row(s) to “{name}” ({self.lists.count(name):,} rows)")
        return name

    def remove_selected(self) -> str | None:
        oids = self.selection.as_array()
        if not len(oids) or not len(self.lists):
            self.show_message("Select some rows first" if len(self.lists) else "No lists yet")
            return None
        name, ok = QInputDialog.getItem(
            self.parent,
            "Remove from List",
            f"Remove {len(oids):,} selected row(s) from list:",
            self.lists.names,
            0,
            False,
        )
        if not ok or name not in self.lists:
            return None
        removed = self.lists.remove(name, oids)
        self.show_message(f"Removed {removed:,} row(s) from “{name}”")
        return name

    def select(self, name: str, *, extend: bool = False) -> None:
        oids = self.lists.oids(name)
        if extend:
            self.selection.add(oids.tolist(), ORIGIN)
        else:
            self.selection.set(oids.tolist(), ORIGIN)
        shown = self.table.oids[self.table.view()]
        hidden = len(oids) - int(np.isin(oids, shown).sum())
        message = f"Selected {len(oids):,} row(s) from “{name}”"
        self.show_message(message + (f"; {hidden:,} hidden by filters" if hidden else ""))

    def show_only(self, name: str) -> None:
        oids = frozenset(self.lists.oids(name).tolist())
        self.table.set_filter(FILTER_KEY, OidFilter(oids, label=name))
        self.show_message(f"Showing only “{name}”")

    def show_all(self) -> None:
        if FILTER_KEY in self.table.filters:
            self.table.set_filter(FILTER_KEY, None)

    def manage(self) -> ManageListsDialog:
        dialog = ManageListsDialog(self, self.parent)
        dialog.open()
        return dialog

    # -- menu ---------------------------------------------------------------------------------
    def fill_menu(self, menu: QMenu, add_action, remove_action) -> None:
        """Rebuild *menu* (called just before it shows)."""
        menu.clear()
        menu.addAction(add_action)
        menu.addAction(remove_action)
        remove_action.setEnabled(bool(len(self.lists)))
        menu.addSeparator()
        if not len(self.lists):
            menu.addAction("No lists").setEnabled(False)
        for name in self.lists.names:
            sub = menu.addMenu(f"{name} ({self.lists.count(name):,})")
            sub.addAction("Select Rows").triggered.connect(lambda _c=False, n=name: self.select(n))
            sub.addAction("Add to Selection").triggered.connect(
                lambda _c=False, n=name: self.select(n, extend=True)
            )
            sub.addAction("Show Only These Rows").triggered.connect(
                lambda _c=False, n=name: self.show_only(n)
            )
        menu.addSeparator()
        if FILTER_KEY in self.table.filters:
            menu.addAction("Stop Showing Only a List").triggered.connect(
                lambda _c=False: self.show_all()
            )
        manage = menu.addAction("Manage Lists…")
        manage.setEnabled(bool(len(self.lists)))
        manage.triggered.connect(lambda _c=False: self.manage())


class ManageListsDialog(QDialog):
    def __init__(self, actions: ListActions, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage Lists")
        self.actions = actions
        self.lists = actions.lists
        self.view = QListWidget()
        self.rename_button = QPushButton("Rename…")
        self.delete_button = QPushButton("Delete")
        self.select_button = QPushButton("Select Rows")
        close = QPushButton("Close")
        buttons = QVBoxLayout()
        for button in (self.select_button, self.rename_button, self.delete_button):
            buttons.addWidget(button)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout = QHBoxLayout(self)
        layout.addWidget(self.view, 1)
        layout.addLayout(buttons)
        self.resize(420, 300)

        self.rename_button.clicked.connect(self.rename_current)
        self.delete_button.clicked.connect(self.delete_current)
        self.select_button.clicked.connect(self.select_current)
        self.view.itemDoubleClicked.connect(lambda _i: self.select_current())
        self.view.currentRowChanged.connect(lambda _r: self._update_buttons())
        close.clicked.connect(self.accept)
        self._unsubscribe = self.lists.subscribe(self.reload)
        self.finished.connect(lambda _r: self._unsubscribe())
        self.reload()

    def reload(self) -> None:
        current = self.current_name()
        self.view.clear()
        for name in self.lists.names:
            item = QListWidgetItem(f"{name}  —  {self.lists.count(name):,} rows")
            item.setData(Qt.ItemDataRole.UserRole, name)
            self.view.addItem(item)
            if name == current:
                self.view.setCurrentItem(item)
        if self.view.currentRow() < 0 and self.view.count():
            self.view.setCurrentRow(0)
        self._update_buttons()

    def current_name(self) -> str | None:
        item = self.view.currentItem()
        return None if item is None else item.data(Qt.ItemDataRole.UserRole)

    def _update_buttons(self) -> None:
        on = self.current_name() is not None
        for button in (self.rename_button, self.delete_button, self.select_button):
            button.setEnabled(on)

    def rename_current(self) -> None:
        old = self.current_name()
        if old is None:
            return
        new, ok = QInputDialog.getText(self, "Rename List", "New name:", text=old)
        if not ok:
            return
        try:
            self.lists.rename(old, new)
        except ListError as err:
            QMessageBox.warning(self, "Rename List", str(err))
            return
        self._select_name(new.strip())

    def delete_current(self) -> None:
        name = self.current_name()
        if name is None:
            return
        answer = QMessageBox.question(
            self, "Delete List", f"Delete the list “{name}”? Its rows stay in the table."
        )
        if answer == QMessageBox.StandardButton.Yes:
            self.lists.delete(name)

    def select_current(self) -> None:
        name = self.current_name()
        if name is not None:
            self.actions.select(name)

    def _select_name(self, name: str) -> None:
        for row in range(self.view.count()):
            if self.view.item(row).data(Qt.ItemDataRole.UserRole) == name:
                self.view.setCurrentRow(row)
