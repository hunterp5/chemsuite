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
"""Rows to select by one column's values (the header menu's Select submenu).

Each takes the rows to look at, in order (usually the shown rows, top to bottom) and returns
the oids to select, in that order:

- :func:`empty` and :func:`not_empty`;
- :func:`first_occurrences`: the first row of each value;
- :func:`duplicates`: every later row of a value seen above, so deleting them keeps one of
  each (with :func:`first_occurrences`, they split the non-empty rows in two);
- :func:`all_duplicated`: every row whose value appears more than once.

Values compare as shown for text and numbers; structures compare as canonical SMILES, so the
same molecule written two ways counts once. Empty cells are never a value.
"""

from __future__ import annotations

import numpy as np

from .columns import STRUCTURE
from .table import Table

EMPTY = "empty"
NOT_EMPTY = "not_empty"
FIRST = "first"
DUPLICATES = "duplicates"
ALL_DUPLICATED = "all_duplicated"


def empty(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    oids = np.asarray(oids, np.int64)
    blank = table.column(name).blank_mask()[table.slots_of(oids)]
    return oids[blank]


def not_empty(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    oids = np.asarray(oids, np.int64)
    blank = table.column(name).blank_mask()[table.slots_of(oids)]
    return oids[~blank]


def first_occurrences(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    oids, codes = _codes(table, name, oids)
    present = codes >= 0
    _, first = np.unique(codes[present], return_index=True)
    return oids[present][np.sort(first)]


def duplicates(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    oids, codes = _codes(table, name, oids)
    present = codes >= 0
    later = np.ones(int(present.sum()), bool)
    _, first = np.unique(codes[present], return_index=True)
    later[first] = False
    return oids[present][later]


def all_duplicated(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    oids, codes = _codes(table, name, oids)
    present = codes >= 0
    _, inverse, counts = np.unique(codes[present], return_inverse=True, return_counts=True)
    return oids[present][counts[inverse] > 1]


PICKS = {
    EMPTY: empty,
    NOT_EMPTY: not_empty,
    FIRST: first_occurrences,
    DUPLICATES: duplicates,
    ALL_DUPLICATED: all_duplicated,
}


def _codes(table: Table, name: str, oids) -> tuple[np.ndarray, np.ndarray]:
    """The oids and a code per row: equal values share a code, empty cells are -1."""
    oids = np.asarray(oids, np.int64)
    column = table.column(name)
    slots = table.slots_of(oids)
    blank = column.blank_mask()[slots]
    keys = column.values(slots)
    if column.kind == STRUCTURE:
        canonical = _canonical(column, slots)
        keys = np.array([c or (None if t is None else f"text:{t}")
                         for c, t in zip(canonical, keys, strict=True)], dtype=object)  # fmt: skip
    codes = np.full(len(oids), -1, np.int64)
    keep = ~blank & np.array([k is not None and k == k for k in keys], bool)  # no NaN, None
    if keep.any():
        values = keys[keep]
        if values.dtype == object:
            values = values.astype(str)
        _, codes[keep] = np.unique(values, return_inverse=True)
    return oids, codes


def _canonical(column, slots: np.ndarray) -> list[str | None]:
    """Canonical SMILES of each cell's molecule (None where there is none)."""
    from chemsuite.core.mols import smiles_of

    out: list[str | None] = []
    for slot in slots.tolist():
        mol = column.obj(slot)
        out.append(smiles_of(mol) or None if mol is not None else None)
    return out
