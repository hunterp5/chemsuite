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

"""The table's undoable operations (the :class:`~chemsuite.core.history.Op` subclasses).

Each captures exactly the data it needs to run in both directions (the old and new cells it
touched, the rows it removed, the column it replaced) and nothing else. Because commands work on
the table core alone, undo never needs a reference to a window or dialog. ``Command`` and
``History`` live in :mod:`chemsuite.core.history`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from chemsuite.core.history import Op

from .columns import Column
from .style import ColumnStyle

if TYPE_CHECKING:
    from .table import Table


@dataclass(eq=False)
class InsertRows(Op):
    oids: np.ndarray  # ascending
    cells: dict[str, Column]  # one entry per table column, len == len(oids)

    def apply(self, table: Table) -> None:
        table._insert_rows(self.oids, self.cells)

    def inverse(self) -> Op:
        return DeleteRows(self.oids, self.cells)

    @property
    def nbytes(self) -> int:
        return self.oids.nbytes + sum(c.nbytes for c in self.cells.values())


@dataclass(eq=False)
class DeleteRows(Op):
    oids: np.ndarray  # ascending
    cells: dict[str, Column]  # snapshot of the removed rows, for the inverse

    def apply(self, table: Table) -> None:
        table._delete_rows(self.oids)

    def inverse(self) -> Op:
        return InsertRows(self.oids, self.cells)

    @property
    def nbytes(self) -> int:
        return self.oids.nbytes + sum(c.nbytes for c in self.cells.values())


@dataclass(eq=False)
class PutCells(Op):
    column: str
    oids: np.ndarray
    new: Column
    old: Column
    structural: bool = field(default=False)

    def apply(self, table: Table) -> None:
        table._put_cells(self.column, self.oids, self.new)

    def inverse(self) -> Op:
        return PutCells(self.column, self.oids, self.old, self.new, self.structural)

    @property
    def nbytes(self) -> int:
        return self.oids.nbytes + self.new.nbytes + self.old.nbytes


@dataclass(eq=False)
class InsertColumn(Op):
    name: str
    position: int
    column: Column

    def apply(self, table: Table) -> None:
        table._insert_column(self.name, self.position, self.column)

    def inverse(self) -> Op:
        return RemoveColumn(self.name, self.position, self.column)

    @property
    def nbytes(self) -> int:
        return self.column.nbytes


@dataclass(eq=False)
class RemoveColumn(Op):
    name: str
    position: int
    column: Column

    def apply(self, table: Table) -> None:
        table._remove_column(self.name)

    def inverse(self) -> Op:
        return InsertColumn(self.name, self.position, self.column)

    @property
    def nbytes(self) -> int:
        return self.column.nbytes


@dataclass(eq=False)
class MoveColumn(Op):
    name: str
    old: int
    new: int

    def apply(self, table: Table) -> None:
        table._move_column(self.name, self.new)

    def inverse(self) -> Op:
        return MoveColumn(self.name, self.new, self.old)


@dataclass(eq=False)
class RenameColumn(Op):
    old: str
    new: str

    def apply(self, table: Table) -> None:
        table._rename_column(self.old, self.new)

    def inverse(self) -> Op:
        return RenameColumn(self.new, self.old)


@dataclass(eq=False)
class ReplaceColumn(Op):
    """Swap a whole column object, e.g. a number column converted to text."""

    name: str
    old: Column
    new: Column

    def apply(self, table: Table) -> None:
        table._replace_column(self.name, self.new)

    def inverse(self) -> Op:
        return ReplaceColumn(self.name, self.new, self.old)

    @property
    def nbytes(self) -> int:
        return self.old.nbytes + self.new.nbytes


@dataclass(eq=False)
class SetStyle(Op):
    """Change how a column is shown (see :mod:`.style`)."""

    name: str
    old: ColumnStyle
    new: ColumnStyle

    def apply(self, table: Table) -> None:
        table._set_style(self.name, self.new)

    def inverse(self) -> Op:
        return SetStyle(self.name, self.new, self.old)
