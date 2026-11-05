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

"""Qt item model over a :class:`~chemsuite.core.table.Table`.

The model is a *mirror*: it caches the column names and the visible rows (as oids) that Qt last
saw, and answers ``rowCount``/``data`` from that mirror. When the table reports a change, the
model compares the mirror with the table and tells Qt precisely what happened:

* only cells changed           -> ``dataChanged`` for the affected rows of each column
* same rows, new order (sort)  -> ``layoutChanged``, with selections and the current cell
                                  carried to the rows' new positions
* rows or columns differ       -> a model reset (the view restores selection and scroll
                                  position from oids; see :mod:`chemsuite.app.table_view`)

All edits go through the table, so they are undoable; the model never stores cell values.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPersistentModelIndex, Qt

from chemsuite.core.table import NUMBER, Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import STRUCTURE, StructureColumn

OID_ROLE = Qt.ItemDataRole.UserRole + 1
KIND_ROLE = Qt.ItemDataRole.UserRole + 2
#: ``(StructureColumn, slot)`` for the structure delegate.
STRUCTURE_ROLE = Qt.ItemDataRole.UserRole + 3
#: ``"blank" | "pending" | "ok" | "invalid"`` for structure cells.
STATUS_ROLE = Qt.ItemDataRole.UserRole + 4

_RIGHT = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
_LEFT = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
_TOOLTIP_CHARS = 400

Index = QModelIndex | QPersistentModelIndex


class TableModel(QAbstractTableModel):
    def __init__(self, table: Table, parent=None) -> None:
        super().__init__(parent)
        self._table = table
        self._columns: list[str] = table.column_names
        self._slots: np.ndarray = table.view()
        self._oids: np.ndarray = table.oids[self._slots]
        self._row_of_slot: np.ndarray | None = None
        self._unsubscribe = table.subscribe(self._on_table_event)

    @property
    def table(self) -> Table:
        return self._table

    def close(self) -> None:
        self._unsubscribe()

    # -- mapping between rows and oids --------------------------------------------------------
    @property
    def visible_oids(self) -> np.ndarray:
        """Oids in display order (read-only)."""
        return self._oids

    def oid_at(self, row: int) -> int:
        return int(self._oids[row])

    def column_name(self, column: int) -> str:
        return self._columns[column]

    def column_index(self, name: str) -> int:
        return self._columns.index(name)

    def rows_of_oids(self, oids: np.ndarray | list[int]) -> np.ndarray:
        """Display rows of *oids* (-1 where the row is hidden or gone), vectorised."""
        oids = np.asarray(oids, dtype=np.int64)
        if not len(oids):
            return np.empty(0, dtype=np.int64)
        table_oids = self._table.oids
        slots = np.searchsorted(table_oids, oids)
        clipped = np.minimum(slots, max(len(table_oids) - 1, 0))
        present = (slots < len(table_oids)) & (
            table_oids[clipped] == oids if len(table_oids) else False
        )
        rows = np.full(len(oids), -1, dtype=np.int64)
        rows[present] = self._row_lookup()[slots[present]]
        return rows

    def _row_lookup(self) -> np.ndarray:
        if self._row_of_slot is None or len(self._row_of_slot) != len(self._table):
            lookup = np.full(len(self._table), -1, dtype=np.int64)
            lookup[self._slots] = np.arange(len(self._slots), dtype=np.int64)
            self._row_of_slot = lookup
        return self._row_of_slot

    # -- QAbstractTableModel ------------------------------------------------------------------
    def rowCount(self, parent: Index = QModelIndex()) -> int:  # noqa: B008 (Qt signature)
        return 0 if parent.isValid() else len(self._slots)

    def columnCount(self, parent: Index = QModelIndex()) -> int:  # noqa: B008
        return 0 if parent.isValid() else len(self._columns)

    def data(self, index: Index, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        name = self._columns[index.column()]
        column = self._table.column(name)
        slot = int(self._slots[index.row()])
        if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
            return column.display(slot)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return _RIGHT if column.kind == NUMBER else _LEFT
        if role == OID_ROLE:
            return int(self._oids[index.row()])
        if role == KIND_ROLE:
            return column.kind
        if isinstance(column, StructureColumn):
            if role == STRUCTURE_ROLE:
                return column, slot
            if role == STATUS_ROLE:
                return _status(column, slot)
            if role == Qt.ItemDataRole.ToolTipRole:
                text = column.display(slot)
                if not text:
                    return None
                prefix = "Invalid structure:\n" if _status(column, slot) == "invalid" else ""
                return prefix + text[:_TOOLTIP_CHARS]
        elif role == Qt.ItemDataRole.ToolTipRole:
            text = column.display(slot)
            return text[:_TOOLTIP_CHARS] if len(text) > 40 or "\n" in text else None
        return None

    def headerData(self, section: int, orientation, role: int = Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal:
            if not 0 <= section < len(self._columns):
                return None
            name = self._columns[section]
            if role == Qt.ItemDataRole.DisplayRole:
                return name
            if role == Qt.ItemDataRole.ToolTipRole:
                return f"{name}\n{self._table.kind(name)} column"
            if role == KIND_ROLE:
                return self._table.kind(name)
            return None
        if role == Qt.ItemDataRole.DisplayRole and 0 <= section < len(self._slots):
            return str(section + 1)
        return None

    def flags(self, index: Index) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemFlag.NoItemFlags
        return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable

    def setData(self, index: Index, value: Any, role: int = Qt.ItemDataRole.EditRole) -> bool:
        if role != Qt.ItemDataRole.EditRole or not index.isValid():
            return False
        name = self._columns[index.column()]
        text = "" if value is None else str(value)
        if text == self._table.cell(name, self.oid_at(index.row())):
            return False
        self._table.set_cells(name, [self.oid_at(index.row())], [text], label=f"Edit {name}")
        return True

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        if not 0 <= column < len(self._columns):
            self._table.set_sort([])
            return
        name = self._columns[column]
        if self._table.column(name).sortable:
            self._table.set_sort([(name, order == Qt.SortOrder.AscendingOrder)])

    # -- following the table --------------------------------------------------------------------
    def _on_table_event(self, event: ev.TableEvent) -> None:
        if not isinstance(event, Changed):
            return
        if event.structural:
            self._resync()
            return
        for change in event.changes:
            if isinstance(change, ev.CellsChanged) and change.column in self._columns:
                rows = self.rows_of_oids(change.oids)
                rows = rows[rows >= 0]
                if len(rows):
                    col = self._columns.index(change.column)
                    self.dataChanged.emit(
                        self.index(int(rows.min()), col), self.index(int(rows.max()), col)
                    )

    def _resync(self) -> None:
        table = self._table
        columns = table.column_names
        slots = table.view()
        oids = table.oids[slots]
        if columns != self._columns or len(oids) != len(self._oids):
            self.beginResetModel()
            self._set_mirror(columns, slots, oids)
            self.endResetModel()
        elif np.array_equal(oids, self._oids):
            # Same rows in the same order; slots may have shifted and kinds may have changed.
            self._set_mirror(columns, slots, oids)
            if len(oids) and columns:
                self.dataChanged.emit(self.index(0, 0), self.index(len(oids) - 1, len(columns) - 1))
            if columns:
                self.headerDataChanged.emit(Qt.Orientation.Horizontal, 0, len(columns) - 1)
        else:
            self.layoutAboutToBeChanged.emit()
            persistent = self.persistentIndexList()
            old_oids = [int(self._oids[i.row()]) for i in persistent]
            self._set_mirror(columns, slots, oids)
            new_rows = self.rows_of_oids(old_oids).tolist()
            self.changePersistentIndexList(
                persistent,
                [
                    self.index(row, i.column()) if row >= 0 else QModelIndex()
                    for i, row in zip(persistent, new_rows, strict=True)
                ],
            )
            self.layoutChanged.emit()

    def _set_mirror(self, columns: list[str], slots: np.ndarray, oids: np.ndarray) -> None:
        self._columns = columns
        self._slots = slots
        self._oids = oids
        self._row_of_slot = None


def _status(column: StructureColumn, slot: int) -> str:
    if not column.display(slot):
        return "blank"
    blob = column.blobs(np.array([slot]))[0]
    if blob is None:
        return "pending"
    return "ok" if blob else "invalid"


__all__ = ["KIND_ROLE", "OID_ROLE", "STATUS_ROLE", "STRUCTURE", "STRUCTURE_ROLE", "TableModel"]
