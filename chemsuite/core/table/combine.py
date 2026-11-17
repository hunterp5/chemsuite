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
"""Another table's rows brought into this one (File ▸ Append and File ▸ Merge), each as one
undoable step.

Append adds every row; columns are matched by name and the new ones added. Merge matches rows
by a key column in each table (text as typed, structures by canonical SMILES): a matched row
takes the other table's values (blank cells never overwrite), and rows with no match are added
or left out. The other table's assets (proteins its rows link to) come too, unless this table
already holds one with the same id.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from chemsuite.core.assets import PutAsset

from .columns import STRUCTURE, Column, blank_column
from .history import InsertColumn, InsertRows, Op, PutCells
from .table import Table


def asset_ops(table: Table, other: Table) -> list[Op]:
    """Add the other table's assets that this one does not hold (an id in both keeps ours)."""
    return [PutAsset(i, None, a) for i, a in other.assets.items() if i not in table.assets]


def append(table: Table, other: Table, *, label: str = "Append rows") -> np.ndarray:
    """Add every row of *other* to *table*; returns the new rows' oids."""
    if not len(other):
        return np.empty(0, dtype=np.int64)
    data = {name: other.column(name) for name in other.column_names}
    return table.add_rows(data, label=label, extra_ops=asset_ops(table, other))


def keys_of(column: Column) -> list[str]:
    """Each row's key ("" for none): a structure's canonical SMILES, else the text shown."""
    slots = np.arange(len(column))
    if column.kind != STRUCTURE:
        return [text.strip() for text in column.displays(slots)]
    from chemsuite.core.mols import smiles_of

    keys = []
    for slot in slots:
        mol = column.obj(int(slot))
        keys.append(smiles_of(mol) if mol is not None else "")
    return keys


@dataclass
class Merged:
    """What a merge did: rows given new values, rows added, columns added."""

    updated: int = 0
    added: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    new_columns: list[str] = field(default_factory=list)
    unmatched: int = 0  # rows of the other table with no match (added, or left out)


def merge(
    table: Table,
    other: Table,
    key: str,
    other_key: str,
    *,
    replace: bool = True,
    add_unmatched: bool = True,
    keys: Sequence[str] | None = None,
    other_keys: Sequence[str] | None = None,
    label: str = "Merge rows",
) -> Merged:
    """Merge *other* into *table*, matching *key* in one to *other_key* in the other.

    A matched row takes the other table's non-blank values; with ``replace=False`` it keeps
    the values it has and only fills its blanks. When the other table has a key more than
    once, its first row is matched and the rest count as unmatched. *keys* and *other_keys*
    (:func:`keys_of` of the two key columns, in storage order) may be passed when worked out
    already, for example in a background job.
    """
    keys = keys_of(table.column(key)) if keys is None else keys
    other_keys = keys_of(other.column(other_key)) if other_keys is None else other_keys
    first: dict[str, int] = {}
    for slot, k in enumerate(other_keys):
        if k and k not in first:
            first[k] = slot
    t_slots = np.array([s for s, k in enumerate(keys) if k in first], dtype=np.int64)
    o_slots = np.array([first[keys[s]] for s in t_slots], dtype=np.int64)
    used = set(o_slots.tolist())
    left = np.array([s for s in range(len(other)) if s not in used], dtype=np.int64)
    new_rows = left if add_unmatched else np.empty(0, dtype=np.int64)
    result = Merged(unmatched=len(left))

    # The other key column lands in ours; any other column of that name is left out.
    sources = {
        (key if name == other_key else name): other.column(name)
        for name in other.column_names
        if name == other_key or name != key
    }
    ops: list[Op] = asset_ops(table, other)
    fitted: dict[str, Column] = {}
    current: dict[str, Column] = {}
    for name, source in sources.items():
        if name in table.column_names:
            fitted[name], current[name] = table._fit(name, source, ops)

    n, k = len(table), len(new_rows)
    if k:
        result.added = np.arange(table.next_oid, table.next_oid + k, dtype=np.int64)
        cells = {
            name: fitted[name].take(new_rows)
            if name in fitted
            else blank_column(table.column(name).kind, k)
            for name in table.column_names
        }
        ops.append(InsertRows(result.added, cells))

    touched: set[int] = set()
    position = len(table.column_names)
    for name, source in sources.items():
        if name in table.column_names:
            continue
        column = blank_column(source.kind, n)
        column.style = source.style
        given = ~source.blank_mask()[o_slots] if len(o_slots) else np.zeros(0, dtype=bool)
        if given.any():
            column.put(t_slots[given], source.take(o_slots[given]))
            touched |= set(t_slots[given].tolist())
        if k:
            column.insert(np.full(k, n, dtype=np.int64), source.take(new_rows))
        ops.append(InsertColumn(name, position, column))
        position += 1
        result.new_columns.append(name)

    oids = table.oids
    for name in fitted:
        if name == key or not len(t_slots):
            continue
        write = ~fitted[name].blank_mask()[o_slots]
        if not replace:
            write &= current[name].blank_mask()[t_slots]
        same = np.array(fitted[name].displays(o_slots)) == np.array(current[name].displays(t_slots))
        write &= ~same  # a value it has already is not a change
        if write.any():
            rows = t_slots[write]
            new, old = fitted[name].take(o_slots[write]), current[name].take(rows)
            ops.append(PutCells(name, oids[rows], new, old))
            touched |= set(rows.tolist())

    result.updated = len(touched)
    table.edit(label, ops)
    return result
