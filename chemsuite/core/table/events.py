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

"""Change notifications emitted by :class:`~chemsuite.core.table.table.Table`.

Every operation is bracketed by exactly one :class:`AboutToChange` and one :class:`Changed`.
``structural`` is True when the visible rows, their order or the columns may differ afterwards
(rows added or removed, columns added/removed/renamed, a sort or filter changed, or cells changed
in a column the view sorts or filters on). Otherwise only cell contents changed, and
``Changed.changes`` lists exactly which ones. The Qt adapter maps the former to a layout change and
the latter to ``dataChanged`` on the affected cells.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, eq=False)
class RowsInserted:
    oids: np.ndarray


@dataclass(frozen=True, eq=False)
class RowsRemoved:
    oids: np.ndarray


@dataclass(frozen=True, eq=False)
class CellsChanged:
    column: str
    oids: np.ndarray


@dataclass(frozen=True)
class ColumnInserted:
    name: str
    position: int


@dataclass(frozen=True)
class ColumnRemoved:
    name: str
    position: int


@dataclass(frozen=True)
class ColumnMoved:
    name: str
    old: int
    new: int


@dataclass(frozen=True)
class ColumnRenamed:
    old: str
    new: str


@dataclass(frozen=True)
class ColumnKindChanged:
    name: str
    kind: str


@dataclass(frozen=True)
class ViewChanged:
    """Filters or sort order changed (or a column they used went away)."""


Change = (
    RowsInserted
    | RowsRemoved
    | CellsChanged
    | ColumnInserted
    | ColumnRemoved
    | ColumnMoved
    | ColumnRenamed
    | ColumnKindChanged
    | ViewChanged
)


@dataclass(frozen=True)
class AboutToChange:
    structural: bool


@dataclass(frozen=True)
class Changed:
    structural: bool
    changes: tuple[Change, ...]


TableEvent = AboutToChange | Changed
