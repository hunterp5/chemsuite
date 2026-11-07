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
from PySide6.QtGui import QColor

from chemsuite.core.table import NUMBER, Changed, Table
from chemsuite.core.table import events as ev
from chemsuite.core.table.columns import REACTION, STRUCTURE, NumberColumn, ParsedColumn
from chemsuite.core.table.summary import column_summary, format_title

from .structure_preview import preview_html

OID_ROLE = Qt.ItemDataRole.UserRole + 1
KIND_ROLE = Qt.ItemDataRole.UserRole + 2
#: ``(StructureColumn or ReactionColumn, slot)`` for the delegate that draws them.
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
        self._text_align = _LEFT
        self._number_align = _RIGHT
        #: column name -> (min, max) of its numbers, for gradient colours; cleared on any change
        self._ranges: dict[str, tuple[float, float] | None] = {}
        self._unsubscribe = table.subscribe(self._on_table_event)

    def set_alignment(self, horizontal: str, vertical: str) -> None:
        """Text cells take *horizontal*; numbers stay right-aligned; both take *vertical*."""
        h = {
            "left": Qt.AlignmentFlag.AlignLeft,
            "center": Qt.AlignmentFlag.AlignHCenter,
            "right": Qt.AlignmentFlag.AlignRight,
        }.get(horizontal, Qt.AlignmentFlag.AlignLeft)
        v = {
            "top": Qt.AlignmentFlag.AlignTop,
            "center": Qt.AlignmentFlag.AlignVCenter,
            "bottom": Qt.AlignmentFlag.AlignBottom,
        }.get(vertical, Qt.AlignmentFlag.AlignVCenter)
        text, number = h | v, Qt.AlignmentFlag.AlignRight | v
        if (text, number) != (self._text_align, self._number_align):
            self._text_align, self._number_align = text, number
            if self.rowCount() and self.columnCount():
                self.dataChanged.emit(
                    self.index(0, 0),
                    self.index(self.rowCount() - 1, self.columnCount() - 1),
                    [Qt.ItemDataRole.TextAlignmentRole],
                )

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
        if role == Qt.ItemDataRole.EditRole and isinstance(column, NumberColumn):
            return column.plain(slot)  # "45.7%" or "1,234" would edit as text
        if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
            return column.display(slot)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return self._number_align if column.kind == NUMBER else self._text_align
        if role == OID_ROLE:
            return int(self._oids[index.row()])
        if role == KIND_ROLE:
            return column.kind
        if role == Qt.ItemDataRole.BackgroundRole:
            return self._background(name, column, slot) if column.style.colour else None
        if isinstance(column, ParsedColumn):
            if role == STRUCTURE_ROLE:
                return column, slot
            if role == STATUS_ROLE:
                return _status(column, slot)
            if role == Qt.ItemDataRole.ToolTipRole:
                text = column.display(slot)
                if not text:
                    return None
                what = "reaction" if column.kind == REACTION else "structure"
                status = _status(column, slot)
                if status == "invalid":
                    return f"Invalid {what}:\n" + text[:_TOOLTIP_CHARS]
                if column.kind == STRUCTURE and status == "ok":  # parsed: drawn on screen
                    mol = column.obj(slot)
                    if mol is not None and mol.GetNumAtoms():
                        return preview_html(mol, text)
                return text[:_TOOLTIP_CHARS]
        elif role == Qt.ItemDataRole.ToolTipRole:
            text = column.display(slot)
            return text[:_TOOLTIP_CHARS] if len(text) > 40 or "\n" in text else None
        return None

    def headerData(self, section: int, orientation, role: int = Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal:
            if not 0 <= section < len(self._columns):
                return None
            name = self._columns[section]
            if name not in self._table:  # asked mid-change, before the mirror catches up
                return name if role == Qt.ItemDataRole.DisplayRole else None
            style = self._table.column(name).style
            if role == Qt.ItemDataRole.DisplayRole:
                return f"{name} (log₁₀)" if style.log10 else name
            if role == Qt.ItemDataRole.ToolTipRole:
                lines = [name, f"{self._table.kind(name)} column"]
                if style.log10:
                    lines.append("Values converted to log₁₀ (Logarithmic)")
                if style.formats_numbers:
                    lines.append(f"Shown as {format_title(style)}")
                lines += ["", *column_summary(self._table.column(name))]
                return "\n".join(lines)
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

    def _background(self, name: str, column, slot: int) -> QColor | None:
        rule = column.style.colour
        if rule.numeric:
            if not isinstance(column, NumberColumn):
                return None
            value = column.value(slot)
            if value is None:
                return None
            if name not in self._ranges:
                self._ranges[name] = column.finite_range()
            bounds = rule.bounds(self._ranges[name])
            if bounds is None:
                return None
            return QColor.fromRgba(rule.number_colour(value, bounds))
        text = column.display(slot)
        return QColor.fromRgba(rule.text_colour(text)) if text else None

    # -- following the table --------------------------------------------------------------------
    def _on_table_event(self, event: ev.TableEvent) -> None:
        if not isinstance(event, Changed):
            return
        self._ranges.clear()
        if event.structural:
            self._resync()
            return
        for change in event.changes:
            if isinstance(change, ev.CellsChanged) and change.column in self._columns:
                rows = self.rows_of_oids(change.oids)
                rows = rows[rows >= 0]
                col = self._columns.index(change.column)
                rule = self._table.column(change.column).style.colour
                if rule is not None and rule.numeric and self._slots.size:
                    # A gradient that follows the data recolours the whole column.
                    rows = np.array([0, len(self._slots) - 1])
                if len(rows):
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


def _status(column: ParsedColumn, slot: int) -> str:
    if not column.display(slot):
        return "blank"
    blob = column.blobs(np.array([slot]))[0]
    if blob is None:
        return "pending"
    return "ok" if blob else "invalid"


__all__ = ["KIND_ROLE", "OID_ROLE", "STATUS_ROLE", "STRUCTURE", "STRUCTURE_ROLE", "TableModel"]
