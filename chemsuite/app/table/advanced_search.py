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
"""Edit ▸ Advanced Search…: a panel of conditions combined with "all" or "any".

Each condition is a column (or any column), an operator and a value, with Not. Conditions sit
in groups: the search is a group, and **Add Group** nests another, with its own "all / any" and
Not, so searches such as ``(A or B) and not (C and D)`` can be built. Text and number conditions are
evaluated at once; structure conditions run as background searches first. The matches can be
shown alone (one filter, ANDed with the Filters panel), selected or added to a list. The search,
whether it is showing, and the recent searches are saved with the document.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.analysis.structure_search import QueryError, StructureQuery
from chemsuite.core.lists import ListError, RowLists
from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import STRUCTURE
from chemsuite.core.table.criteria import (
    ANY_COLUMN,
    CriteriaFilter,
    CriteriaSearch,
    Criterion,
    CriterionError,
    ops_for,
)

from ..widgets.flow_layout import FlowRow, compact
from ..widgets.pane import menu_button, pane_layout, section
from .structure_search import QueryInput, StructureSearcher

FILTER_KEY = "advanced"
EXTRAS_KEY = "advanced_search"
MAX_RECENT = 10
ANY_LABEL = "Any column"
ORIGIN = "advanced search"


class CriterionRow(QFrame):
    """One condition: column, operator, value(s), match case, Not, remove."""

    changed = Signal()
    removeRequested = Signal()

    def __init__(self, table: Table, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        section(self, "CriterionRow")
        self.column = QComboBox()
        self.op = QComboBox()
        compact(self.column, self.op, chars=8)
        self.value = QLineEdit()
        self.value.setPlaceholderText("value")
        self.value.setClearButtonEnabled(True)
        self.and_label = QLabel("and")
        self.value2 = QLineEdit()
        self.value2.setPlaceholderText("value")
        self.structure = QueryInput()
        self.case = QCheckBox("Aa")
        self.case.setToolTip("Match case")
        self.negate = QCheckBox("Not")
        self.negate.setToolTip("Keep the rows this condition does not match")
        remove = QToolButton()
        remove.setText("✕")
        remove.setAutoRaise(True)
        remove.setToolTip("Remove condition")
        self.error = QLabel()
        self.error.setWordWrap(True)
        self.error.setStyleSheet("color: #d33")
        self.error.hide()

        flow = FlowRow()
        for widget in (
            self.column,
            self.op,
            self.value,
            self.and_label,
            self.value2,
            self.structure,
            self.case,
            self.negate,
            remove,
        ):
            flow.addWidget(widget)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 8, 0, 10)
        layout.setSpacing(4)
        layout.addWidget(flow)
        layout.addWidget(self.error)

        self.refresh_columns()
        self.column.currentIndexChanged.connect(lambda _i: self._column_changed())
        self.op.currentIndexChanged.connect(lambda _i: self._op_changed())
        for edit in (self.value, self.value2):
            edit.textChanged.connect(lambda _t: self.changed.emit())
        self.structure.changed.connect(self.changed.emit)
        self.case.toggled.connect(lambda _on: self.changed.emit())
        self.negate.toggled.connect(lambda _on: self.changed.emit())
        remove.clicked.connect(self.removeRequested.emit)

    # -- choices ------------------------------------------------------------------------------
    def column_name(self) -> str:
        data = self.column.currentData()
        return ANY_COLUMN if data is None else data

    def _kind(self) -> str | None:
        name = self.column_name()
        return self.table.kind(name) if name and name in self.table else None

    def refresh_columns(self) -> None:
        """Columns changed: list them again, keeping the choice (a missing one stays listed)."""
        current = self.column_name() if self.column.count() else ANY_COLUMN
        self.column.blockSignals(True)
        self.column.clear()
        self.column.addItem(ANY_LABEL, None)
        for name in self.table.column_names:
            self.column.addItem(name, name)
        if current and self.column.findData(current) < 0:
            self.column.addItem(f"{current} (missing)", current)
        index = self.column.findData(current) if current else 0
        self.column.setCurrentIndex(max(0, index))
        self.column.blockSignals(False)
        self._fill_ops()

    def _fill_ops(self) -> None:
        current = self.op.currentData()
        ops = ops_for(self._kind())
        self.op.blockSignals(True)
        self.op.clear()
        for key, title in ops.items():
            self.op.addItem(title, key)
        index = self.op.findData(current)
        self.op.setCurrentIndex(index if index >= 0 else 0)
        self.op.blockSignals(False)
        self._show_value_fields()

    def _column_changed(self) -> None:
        self._fill_ops()
        self.changed.emit()

    def _op_changed(self) -> None:
        self._show_value_fields()
        self.changed.emit()

    def _show_value_fields(self) -> None:
        op = self.op.currentData()
        structure = op in ("substructure", "exact")
        self.value.setVisible(op != "blank" and not structure)
        self.and_label.setVisible(op == "between")
        self.value2.setVisible(op == "between")
        self.structure.setVisible(structure)
        self.case.setVisible(self._kind() not in ("number", STRUCTURE) and op != "blank")

    # -- the condition ------------------------------------------------------------------------
    def criterion(self) -> Criterion:
        op = self.op.currentData() or "contains"
        structure = op in ("substructure", "exact")
        return Criterion(
            self.column_name(),
            op,
            self.structure.text() if structure else self.value.text(),
            self.value2.text() if op == "between" else "",
            case_sensitive=self.case.isChecked() and self.case.isVisible(),
            negate=self.negate.isChecked(),
        )

    def set_criterion(self, criterion: Criterion) -> None:
        widgets = (self.column, self.op, self.value, self.value2, self.case, self.negate)
        for widget in widgets:
            widget.blockSignals(True)
        try:
            index = self.column.findData(criterion.column) if criterion.column else 0
            if index < 0:
                self.column.addItem(f"{criterion.column} (missing)", criterion.column)
                index = self.column.count() - 1
            self.column.setCurrentIndex(index)
            self._fill_ops()
            self.op.setCurrentIndex(max(0, self.op.findData(criterion.op)))
            if criterion.op in ("substructure", "exact"):
                self.structure.set_text(criterion.value)
            else:
                self.value.setText(criterion.value)
            self.value2.setText(criterion.value2)
            self.case.setChecked(criterion.case_sensitive)
            self.negate.setChecked(criterion.negate)
        finally:
            for widget in widgets:
                widget.blockSignals(False)
        self._show_value_fields()

    def show_error(self, text: str) -> None:
        self.error.setText(text)
        self.error.setVisible(bool(text))


class GroupBox(QFrame):
    """A group of conditions and nested groups: "Match all/any of", optionally "Not".

    The top-level group is the whole search; nested groups are indented with a rule down their
    left side, and have Not and a remove button.
    """

    changed = Signal()
    removeRequested = Signal()

    def __init__(self, table: Table, *, top: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.table = table
        self.top = top
        self.items: list[CriterionRow | GroupBox] = []
        self.setObjectName("SearchGroup")
        if not top:
            self.setStyleSheet(
                "QFrame#SearchGroup { border: none; border-left: 2px solid palette(mid); }"
            )
        self.negate = QCheckBox("Not")
        self.negate.setToolTip("Keep the rows this group does not match")
        self.negate.setVisible(not top)
        self.combine = QComboBox()
        self.combine.addItem("all", "all")
        self.combine.addItem("any", "any")
        self.combine.setToolTip("all: every condition must match (and). any: one is enough (or).")
        remove = QToolButton()
        remove.setText("✕")
        remove.setAutoRaise(True)
        remove.setToolTip("Remove this group")
        remove.setVisible(not top)
        header_row = QHBoxLayout()
        header_row.setContentsMargins(0, 0, 0, 0)
        header_row.setSpacing(6)
        header_row.addWidget(self.negate)
        if not top:
            header_row.addSpacing(6)
        header_row.addWidget(QLabel("Match"))
        header_row.addWidget(self.combine)
        header_row.addWidget(QLabel("of these conditions" if top else "of these"))
        header_row.addStretch(1)
        header_row.addWidget(remove)
        header = QWidget()
        header.setLayout(header_row)

        self.add_condition_button = QPushButton("Add Condition")
        self.add_group_button = QPushButton("Add Group")
        self.add_group_button.setToolTip(
            "A group of conditions with its own all / any and Not, e.g. not (A and B)"
        )
        buttons = FlowRow()
        buttons.addWidget(self.add_condition_button)
        buttons.addWidget(self.add_group_button)

        self._items_box = QVBoxLayout()
        self._items_box.setContentsMargins(0, 0, 0, 0)
        self._items_box.setSpacing(0)  # rows divide themselves with hairlines
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0 if top else 10, 4, 0, 4)
        layout.setSpacing(2)
        layout.addWidget(header)
        layout.addLayout(self._items_box)
        layout.addWidget(buttons)

        self.negate.toggled.connect(lambda _on: self.changed.emit())
        self.combine.currentIndexChanged.connect(lambda _i: self.changed.emit())
        remove.clicked.connect(self.removeRequested.emit)
        self.add_condition_button.clicked.connect(lambda: self.add_row())
        self.add_group_button.clicked.connect(lambda: self.add_group())

    # -- building -------------------------------------------------------------------------------
    def add_row(self, criterion: Criterion | None = None) -> CriterionRow:
        row = CriterionRow(self.table)
        if criterion is not None:
            row.set_criterion(criterion)
        self._add(row)
        return row

    def add_group(self, search: CriteriaSearch | None = None) -> GroupBox:
        """Add a nested group (by default holding one empty condition)."""
        group = GroupBox(self.table)
        if search is not None:
            group.set_search(search)
        else:
            group.add_row()
        self._add(group)
        return group

    def _add(self, item: CriterionRow | GroupBox) -> None:
        self.items.append(item)
        self._items_box.addWidget(item)
        item.changed.connect(self.changed.emit)
        item.removeRequested.connect(lambda: self.remove(item))
        self.changed.emit()

    def remove(self, item: CriterionRow | GroupBox) -> None:
        if item in self.items:
            self.items.remove(item)
            item.setParent(None)
            item.deleteLater()
        if self.top and not self.items:
            self.add_row()  # the search always has a condition to fill in
        self.changed.emit()

    def clear(self) -> None:
        for item in list(self.items):
            self.items.remove(item)
            item.setParent(None)
            item.deleteLater()

    # -- the search ------------------------------------------------------------------------------
    def search(self) -> CriteriaSearch:
        return CriteriaSearch(
            tuple(
                item.search() if isinstance(item, GroupBox) else item.criterion()
                for item in self.items
            ),
            self.combine.currentData(),
            self.negate.isChecked() and not self.top,
        )

    def set_search(self, search: CriteriaSearch) -> None:
        self.blockSignals(True)
        try:
            self.clear()
            self.combine.setCurrentIndex(max(0, self.combine.findData(search.combine)))
            self.negate.setChecked(search.negate and not self.top)
            for item in search.criteria:
                if isinstance(item, CriteriaSearch):
                    self.add_group(item)
                else:
                    self.add_row(item)
        finally:
            self.blockSignals(False)

    def rows(self) -> list[CriterionRow]:
        """Every condition row, depth first (the order of ``CriteriaSearch.leaves``)."""
        out: list[CriterionRow] = []
        for item in self.items:
            out.extend(item.rows() if isinstance(item, GroupBox) else [item])
        return out

    def groups(self) -> list[GroupBox]:
        """The nested groups, depth first."""
        out: list[GroupBox] = []
        for item in self.items:
            if isinstance(item, GroupBox):
                out += [item, *item.groups()]
        return out

    def refresh_columns(self) -> None:
        for row in self.rows():
            row.refresh_columns()


class AdvancedSearchPanel(QWidget):
    """The conditions, what to do with the matches, and recent searches."""

    def __init__(
        self,
        table: Table,
        selection: Selection,
        lists: RowLists,
        searcher: StructureSearcher,
        show_message: Callable[[str], None],
        saved: dict[str, Any] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self.selection = selection
        self.lists = lists
        self.searcher = searcher
        self.show_message = show_message
        self.recent: list[CriteriaSearch] = []
        self.showing: CriteriaSearch | None = None  # the search the table is filtered by
        self.job = None  # the structure search running, if any
        self._pending: tuple | None = None

        #: the whole search: the top-level group
        self.root = GroupBox(table, top=True)
        self.combine = self.root.combine
        self.recent_button, self.recent_menu = menu_button("Recent")
        self.recent_menu.aboutToShow.connect(self._fill_recent)
        host = QWidget()
        host_box = QVBoxLayout(host)
        host_box.setContentsMargins(0, 0, 0, 0)
        host_box.addWidget(self.root)
        host_box.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(host)

        self.show_button = QPushButton("Show Only")
        self.show_button.setToolTip("Hide the rows that do not match (with the other filters)")
        self.select_button = QPushButton("Select")
        self.list_button = QPushButton("Add to List:")
        self.list_name = QComboBox()
        self.list_name.setEditable(True)
        self.list_name.lineEdit().setPlaceholderText("List name")
        compact(self.list_name, chars=10)
        self.clear_button = QPushButton("Show All")
        self.clear_button.setToolTip("Stop filtering by this search")
        to_list = QWidget()  # the button and its list name stay on one line
        to_list_row = QHBoxLayout(to_list)
        to_list_row.setContentsMargins(0, 0, 0, 0)
        to_list_row.setSpacing(4)
        to_list_row.addWidget(self.list_button)
        to_list_row.addWidget(self.list_name)
        actions = FlowRow()
        for widget in (self.show_button, self.select_button, self.clear_button, to_list):
            actions.addWidget(widget)
        tools = FlowRow()
        tools.addWidget(self.recent_button)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: palette(placeholder-text)")

        layout = pane_layout(QVBoxLayout(self))
        layout.addWidget(scroll, 1)
        layout.addWidget(tools)
        layout.addWidget(actions)
        layout.addWidget(self.status)
        self.setMinimumWidth(220)

        self._rerun = QTimer(self)
        self._rerun.setSingleShot(True)
        self._rerun.setInterval(400)
        self._rerun.timeout.connect(self._refresh_showing)
        self.root.changed.connect(self._edited)
        self.show_button.clicked.connect(lambda: self.run("show"))
        self.select_button.clicked.connect(lambda: self.run("select"))
        self.list_button.clicked.connect(lambda: self.run("list"))
        self.clear_button.clicked.connect(self.show_all)
        self._unsubscribe = table.subscribe(self._on_table_event)
        self._unsubscribe_lists = lists.subscribe(lambda *_: self._fill_list_names())
        self._fill_list_names()
        self.list_name.editTextChanged.connect(lambda _t: self._update())

        self.load_state(saved or {})
        if not self.rows:
            self.add_row()
        self._update()

    # -- rows ---------------------------------------------------------------------------------
    @property
    def rows(self) -> list[CriterionRow]:
        """Every condition row, depth first through the groups."""
        return self.root.rows()

    def add_row(self, criterion: Criterion | None = None) -> CriterionRow:
        """Add a condition to the top-level group."""
        return self.root.add_row(criterion)

    def add_group(self, search: CriteriaSearch | None = None) -> GroupBox:
        """Add a nested group to the top-level group."""
        return self.root.add_group(search)

    def remove_row(self, row: CriterionRow) -> None:
        for group in [self.root, *self.root.groups()]:
            if row in group.items:
                group.remove(row)
                return

    def current_search(self) -> CriteriaSearch:
        return self.root.search()

    def set_search(self, search: CriteriaSearch) -> None:
        self.root.set_search(search)
        if not self.root.items:
            self.root.add_row()
        self._update()

    def _edited(self) -> None:
        self._update()

    def _problems(self) -> list[str]:
        """One message per condition row ("" where the row is fine), depth first."""
        problems = []
        for row in self.rows:
            try:
                criterion = row.criterion()
                criterion.check(self.table)
                if criterion.needs_structure_search:
                    from chemsuite.analysis.structure_search import parse_query

                    parse_query(StructureQuery(criterion.value, criterion.op))
            except (CriterionError, QueryError) as err:
                problems.append(str(err))
            else:
                problems.append("")
        return problems

    def _update(self) -> None:
        problems = self._problems()
        for row, problem in zip(self.rows, problems, strict=True):
            # Do not shout about a row the user has not filled in yet.
            row.show_error(
                problem if row.criterion().value or row.criterion().op == "blank" else ""
            )
        empty_group = any(not g.items for g in self.root.groups())
        ok = not any(problems) and not empty_group and self.job is None
        self.show_button.setEnabled(ok)
        self.select_button.setEnabled(ok)
        self.list_button.setEnabled(ok and bool(self.list_name.currentText().strip()))
        self.clear_button.setEnabled(self.showing is not None)
        self.recent_button.setEnabled(bool(self.recent))

    # -- running ------------------------------------------------------------------------------
    def run(self, action: str) -> bool:
        """Search, then show / select / list the matches; False if the search is not valid."""
        search = self.current_search()
        try:
            search.check(self.table)
        except CriterionError as err:
            self.status.setText(str(err))
            return False
        if any(self._problems()):
            self.status.setText("Fix the conditions marked in red.")
            return False
        list_name = self.list_name.currentText().strip()
        self._remember(search)
        structure = [i for i, c in enumerate(search.leaves()) if c.needs_structure_search]
        self._search_structures(
            search, structure, {}, lambda hits: self._finish(action, search, hits, list_name)
        )
        return True

    def _search_structures(self, search, todo: list[int], hits: dict, then) -> None:
        """Run the structure conditions one after another, then call *then(hits)*."""
        if not todo:
            self.job = None
            self._update()
            then(hits)
            return
        index, rest = todo[0], todo[1:]
        criterion = search.leaves()[index]
        query = StructureQuery(criterion.value, criterion.op)
        self.status.setText(f"Searching structures ({len(hits) + 1} of {len(hits) + len(todo)})…")

        def done(oids: np.ndarray, result) -> None:
            hits[index] = oids[result.hits]
            self._search_structures(search, rest, hits, then)

        self.job = self.searcher.start(criterion.column, query, done)
        self.job.cancelled.connect(self._cancelled)
        self.job.failed.connect(lambda _e: self._cancelled())
        self._update()

    def _cancelled(self) -> None:
        self.job = None
        self.status.setText("Search stopped.")
        self._update()

    def _finish(self, action: str, search: CriteriaSearch, hits: dict, list_name: str) -> None:
        mask = search.mask(self.table, hits)
        oids = self.table.oids[mask]
        if action == "show":
            frozen = tuple((i, frozenset(o.tolist())) for i, o in sorted(hits.items()))
            self.table.set_filter(FILTER_KEY, CriteriaFilter(search, frozen))
            self.showing = search
            message = f"{len(oids):,} row(s) match: showing only those"
        elif action == "select":
            self.selection.set(oids.tolist(), ORIGIN)
            message = f"Selected {len(oids):,} matching row(s)"
        else:
            try:
                self.lists.add(list_name, oids)
            except ListError as err:
                self.status.setText(str(err))
                return
            message = f"Added {len(oids):,} matching row(s) to “{list_name}”"
        self.status.setText(message)
        self.show_message(message)
        self._update()

    def show_all(self) -> None:
        if FILTER_KEY in self.table.filters:
            self.table.set_filter(FILTER_KEY, None)
        self.showing = None
        self.status.setText("")
        self._update()

    # -- following the table ------------------------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        columns_changed = False
        structures_changed = False
        used = set(self.showing.columns) if self.showing is not None else set()
        for change in event.changes:
            if isinstance(
                change,
                ev.ColumnInserted | ev.ColumnRemoved | ev.ColumnRenamed | ev.ColumnKindChanged,
            ):
                columns_changed = True
            if isinstance(change, ev.RowsInserted):
                structures_changed = True
            if isinstance(change, ev.CellsChanged) and change.column in used:
                structures_changed |= self.table.kind(change.column) == STRUCTURE
        if columns_changed:
            self.root.refresh_columns()
            self._update()
        if self.showing is not None and FILTER_KEY not in self.table.filters:
            # The table dropped the filter (one of its columns was removed or converted).
            self.showing = None
            self.status.setText("The search stopped: a column it used changed.")
            self._update()
        elif (
            self.showing is not None
            and structures_changed
            and any(c.needs_structure_search for c in self.showing.leaves())
        ):
            self._rerun.start()  # structure matches were found before the edit: search again

    def _refresh_showing(self) -> None:
        search = self.showing
        if search is None or self.job is not None:
            return
        structure = [i for i, c in enumerate(search.leaves()) if c.needs_structure_search]
        self._search_structures(
            search, structure, {}, lambda hits: self._finish("show", search, hits, "")
        )

    # -- recent and saved state -----------------------------------------------------------------
    def _remember(self, search: CriteriaSearch) -> None:
        self.recent = [search] + [s for s in self.recent if s != search]
        del self.recent[MAX_RECENT:]

    def _fill_recent(self) -> None:
        self.recent_menu.clear()
        for search in self.recent:
            text = search.describe()
            action = self.recent_menu.addAction(text if len(text) <= 80 else text[:79] + "…")
            action.setToolTip(text)
            action.triggered.connect(lambda _c=False, s=search: self.set_search(s))

    def _fill_list_names(self) -> None:
        current = self.list_name.currentText()
        self.list_name.blockSignals(True)
        self.list_name.clear()
        self.list_name.addItems(self.lists.names)
        self.list_name.setCurrentText(current)
        self.list_name.blockSignals(False)

    def state(self) -> dict[str, Any]:
        """What to save with the document ({} when there is nothing worth keeping)."""
        search = self.current_search()
        filled = any(c.value or c.op == "blank" for c in search.leaves())
        if not filled and not self.recent:
            return {}
        return {
            "search": search.to_json(),
            "showing": self.showing is not None,
            "recent": [s.to_json() for s in self.recent],
        }

    def load_state(self, saved: dict[str, Any]) -> None:
        recent = []
        for entry in saved.get("recent") or []:
            try:
                recent.append(CriteriaSearch.from_json(entry))
            except (ValueError, TypeError, AttributeError):
                continue
        self.recent = recent[:MAX_RECENT]
        if "search" in saved:
            try:
                self.set_search(CriteriaSearch.from_json(saved["search"]))
            except (ValueError, TypeError, AttributeError):
                return
            if saved.get("showing"):
                QTimer.singleShot(0, lambda: self.run("show"))

    def close_panel(self) -> None:
        self._rerun.stop()
        if self.job is not None and self.job.running:
            self.job.cancel()
        self._unsubscribe()
        self._unsubscribe_lists()
