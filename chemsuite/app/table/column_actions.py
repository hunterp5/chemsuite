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

"""Column operations behind the table header's context menu.

Each operation asks for what it needs (a name, a confirmation) and then makes one undoable change
to the table. The dialogs are plain ``QInputDialog``/``QMessageBox`` calls so tests can answer
them by monkeypatching (Colour Column… has its own dialog, :mod:`.column_colour`).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from PySide6.QtGui import QAction, QActionGroup
from PySide6.QtWidgets import QInputDialog, QMenu, QMessageBox, QWidget

from chemsuite.core.selection import Selection
from chemsuite.core.table import NUMBER, TEXT, Table
from chemsuite.core.table import picks as row_picks
from chemsuite.core.table.columns import REACTION, STRUCTURE, to_number
from chemsuite.core.wording import plural

KIND_TITLES = {TEXT: "Text", NUMBER: "Number", STRUCTURE: "Structure", REACTION: "Reaction"}
#: the Select submenu: title, pick, tooltip (None: a separator)
SELECT_PICKS = (
    ("Empty", row_picks.EMPTY, "Select the shown rows where this cell is empty."),
    ("Not Empty", row_picks.NOT_EMPTY, "Select the shown rows where this cell has a value."),
    None,
    (
        "First Occurrence",
        row_picks.FIRST,
        "Select the first shown row of each value, top to bottom. Structures count as the "
        "same when they are the same molecule.",
    ),
    (
        "Duplicates",
        row_picks.DUPLICATES,
        "Select every shown row whose value appears in a row above it; deleting these keeps "
        "one row of each value.",
    ),
    (
        "All Duplicated Values",
        row_picks.ALL_DUPLICATED,
        "Select every shown row whose value appears more than once, the first included.",
    ),
)
SELECT_ORIGIN = "column-menu"


class ColumnActions:
    def __init__(
        self,
        parent: QWidget,
        table: Table,
        show_message: Callable[[str], None] = lambda _m: None,
        add_filter: Callable[[str], object] | None = None,
        search: Callable[[str], object] | None = None,
        selection: Selection | None = None,
        freeze: Callable[[int], object] | None = None,
        frozen: Callable[[], int] = lambda: 0,
        align: Callable[[str], object] | None = None,
        unalign: Callable[[str], object] | None = None,
        aligned: Callable[[str], bool] = lambda _name: False,
    ) -> None:
        self._parent = parent
        self._table = table
        self._show = show_message
        self._add_filter = add_filter
        self._search = search
        self._selection = selection
        self._freeze = freeze
        self._frozen = frozen
        self._align = align
        self._unalign = unalign
        self._aligned = aligned

    # -- sorting ------------------------------------------------------------------------------
    def sort(self, name: str, ascending: bool) -> None:
        self._table.set_sort([(name, ascending)])

    def add_to_sort(self, name: str, ascending: bool) -> None:
        """Sort by *name* within the current sort (as the least significant key)."""
        keys = [k for k in self._table.sort_keys if k[0] != name]
        self._table.set_sort([*keys, (name, ascending)])

    def clear_sort(self) -> None:
        self._table.set_sort([])

    # -- selecting rows by value --------------------------------------------------------------
    def select_rows(self, name: str, pick: str) -> int:
        """Select the shown rows that ``pick`` (:mod:`chemsuite.core.table.picks`) finds in
        column *name*; returns how many."""
        if self._selection is None:
            return 0
        oids = row_picks.PICKS[pick](self._table, name, self._table.visible_oids())
        self._selection.set(oids.tolist(), SELECT_ORIGIN)
        self._show(f"Selected {plural(len(oids), 'row')}")
        return len(oids)

    # -- structure ----------------------------------------------------------------------------
    def rename(self, name: str) -> bool:
        new, ok = QInputDialog.getText(self._parent, "Rename Column", "Column name:", text=name)
        new = new.strip()
        if not ok or not new or new == name:
            return False
        if new in self._table:
            QMessageBox.warning(self._parent, "Rename Column", f"“{new}” already exists.")
            return False
        self._table.rename_column(name, new)
        return True

    def delete(self, names: Sequence[str]) -> None:
        if names:
            self._table.remove_columns(list(names))
            self._show(
                f"Deleted column {names[0]}" if len(names) == 1 else f"Deleted {len(names)} columns"
            )

    def duplicate(self, name: str) -> str:
        new = self._table.duplicate_column(name)
        self._show(f"Duplicated {name} as {new}")
        return new

    # -- numbers ------------------------------------------------------------------------------
    def log_problem(self, name: str) -> str:
        """Why Logarithmic cannot be switched on or off for *name* ("" when it can)."""
        table = self._table
        if table.kind(name) != NUMBER:
            return "Only Number columns can be made logarithmic."
        column = table.column(name)
        if column.style.log10:
            return ""
        values = column.values()
        values = values[~np.isnan(values)]
        if not len(values):
            return "The column has no numbers."
        if (values <= 0).any():
            return "log₁₀ needs every value above zero."
        return ""

    def toggle_log(self, name: str) -> bool:
        """Convert *name*'s values to log₁₀, or back if they already are. One undo step."""
        problem = self.log_problem(name)
        if problem:
            QMessageBox.warning(self._parent, "Logarithmic", problem)
            return False
        on = not self._table.column(name).style.log10
        count = self._table.set_log10(name, on)
        self._show(
            f"Converted {plural(count, 'value')} in “{name}” to log₁₀"
            if on
            else f"Converted {plural(count, 'value')} in “{name}” back from log₁₀"
        )
        return True

    def number_format(self, name: str):
        """Open Number Format… for *name* (display only; one undoable step on OK)."""
        from .number_format import NumberFormatDialog

        dialog = NumberFormatDialog(self._parent, self._table, name)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    def colour(self, name: str):
        from .column_colour import ColumnColourDialog

        dialog = ColumnColourDialog(self._parent, self._table, name)
        dialog.finished.connect(lambda _r: dialog.deleteLater())
        dialog.open()
        return dialog

    # -- type ---------------------------------------------------------------------------------
    def convert(self, name: str, kind: str) -> bool:
        table = self._table
        if table.kind(name) == kind:
            return False
        if kind != NUMBER:
            table.convert_column(name, kind)
            return True
        column = table.column(name)
        if to_number(column) is not None:
            table.convert_column(name, NUMBER)
            return True
        forced = to_number(column, force=True)
        lost = forced[1] if forced is not None else len(column)
        answer = QMessageBox.question(
            self._parent,
            "Convert to Number",
            f"{lost:,} value(s) in “{name}” are not plain numbers and will be cleared.\n\n"
            "Convert anyway? (You can undo this.)",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return False
        blanked = table.convert_column(name, NUMBER, force=True)
        self._show(f"Converted “{name}” to numbers; cleared {plural(blanked, 'value')}")
        return True

    # -- menu ---------------------------------------------------------------------------------
    def build_menu(self, name: str, selected: Sequence[str], position: int) -> QMenu:
        """Context menu for the header of column *name* (at *position*).

        *selected* are the fully selected columns; delete applies to them when *name* is one.
        """
        table = self._table
        menu = QMenu(self._parent)
        menu.setToolTipsVisible(True)
        sortable = table.column(name).sortable

        def add(
            title: str, slot: Callable[[], object], enabled: bool = True, into: QMenu = menu
        ) -> QAction:
            action = into.addAction(title)
            action.triggered.connect(lambda _checked=False: slot())
            action.setEnabled(enabled)
            return action

        drawn = table.kind(name) in (STRUCTURE, REACTION)
        if self._selection is not None:
            chooser = menu.addMenu("Select")
            chooser.setToolTipsVisible(True)
            for item in SELECT_PICKS:
                if item is None:
                    chooser.addSeparator()
                    continue
                title, pick, tip = item
                add(title, lambda p=pick: self.select_rows(name, p), True, chooser).setToolTip(tip)
        sort = menu.addMenu("Sort")
        sort.setToolTipsVisible(True)
        add("Ascending", lambda: self.sort(name, True), sortable, sort)
        add("Descending", lambda: self.sort(name, False), sortable, sort)
        if any(k[0] != name for k in table.sort_keys):
            then = add("Then Ascending", lambda: self.add_to_sort(name, True), sortable, sort)
            then.setToolTip("Sort by this column within the current sort")
            then = add("Then Descending", lambda: self.add_to_sort(name, False), sortable, sort)
            then.setToolTip("Sort by this column within the current sort")
        sort.addSeparator()
        add("Clear", self.clear_sort, bool(table.sort_keys), sort)
        colour = add("Colour Column…", lambda: self.colour(name), not drawn)
        colour.setToolTip(
            "Colour the cells by value (gradient) or by category."
            if not drawn
            else f"{KIND_TITLES[table.kind(name)]} columns cannot be coloured."
        )
        menu.addSeparator()

        add("Rename…", lambda: self.rename(name))
        add("Duplicate Column", lambda: self.duplicate(name))
        targets = list(selected) if name in selected else [name]
        add(
            "Delete Column" if len(targets) == 1 else f"Delete {len(targets)} Columns",
            lambda: self.delete(targets),
        )
        menu.addSeparator()

        if self._add_filter is not None:
            add_filter = self._add_filter
            add("Filter", lambda: add_filter(name))
        if self._search is not None:
            search = self._search
            add("Search", lambda: search(name)).setToolTip(
                "Open the search bar to search this column"
            )
        if self._freeze is not None:
            freeze = self._freeze
            frozen = self._frozen()
            through = add("Freeze Columns up to Here", lambda: freeze(position + 1))
            through.setToolTip(
                "Keep this column and the ones left of it in place while the rest scroll "
                "sideways. Drag columns such as ID and Structure to the front first."
            )
            through.setEnabled(frozen != position + 1)
            if frozen:
                add("Unfreeze Columns", lambda: freeze(0))
        if table.kind(name) == STRUCTURE and self._align is not None:
            align = self._align
            aligning = add("Align Drawings…", lambda: align(name))
            aligning.setToolTip(
                "Draw every structure that contains a core with the core in the same place "
                "and orientation; the structures themselves are not changed."
            )
            if self._aligned(name) and self._unalign is not None:
                unalign = self._unalign
                add("Stop Aligning Drawings", lambda: unalign(name))
        kinds = menu.addMenu("Column Type")
        group = QActionGroup(kinds)
        for kind, title in KIND_TITLES.items():
            action = kinds.addAction(title)
            action.setCheckable(True)
            action.setChecked(table.kind(name) == kind)
            group.addAction(action)
            action.triggered.connect(lambda _checked=False, k=kind: self.convert(name, k))
        menu.addSeparator()

        number = table.kind(name) == NUMBER
        log = add("Logarithmic (log₁₀)", lambda: self.toggle_log(name), not self.log_problem(name))
        log.setCheckable(True)
        log.setChecked(table.column(name).style.log10)
        log.setToolTip(
            self.log_problem(name)
            or "Convert the values to log₁₀; choose again to convert them back. You can undo it."
        )
        formatting = add("Number Format…", lambda: self.number_format(name), number)
        formatting.setToolTip(
            "Decimals, thousands separator, scientific, significant figures or percent; the "
            "stored values are not rounded."
            if number
            else "Only Number columns have a number format."
        )
        return menu
