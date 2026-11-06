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

"""Text search over a table: quick search, find and find/replace.

Everything matches the cells' displayed text, the same text the user sees and copies.
Structure columns are skipped: their text is SMILES or a molfile, which a text search would
match in surprising ways (substructure search is the tool for those).
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from .columns import REACTION, STRUCTURE

if TYPE_CHECKING:
    from .table import Table

SearchMode = Literal["contains", "whole", "regex"]


@dataclass(frozen=True)
class Query:
    """What to look for. ``whole`` matches the entire cell; ``regex`` uses Python syntax."""

    text: str
    mode: SearchMode = "contains"
    case_sensitive: bool = False

    def pattern(self) -> re.Pattern[str]:
        """The compiled pattern; raises :class:`re.error` for an invalid regular expression."""
        flags = 0 if self.case_sensitive else re.IGNORECASE
        if self.mode == "regex":
            return re.compile(self.text, flags)
        escaped = re.escape(self.text)
        return re.compile(f"^{escaped}$" if self.mode == "whole" else escaped, flags)

    def predicate(self) -> Callable[[str], bool]:
        pattern = self.pattern()
        if self.mode == "whole":
            return lambda text: pattern.fullmatch(text) is not None
        return lambda text: pattern.search(text) is not None

    def replace(self, text: str, replacement: str) -> str:
        """*text* with every match replaced (``\\1`` and ``\\g<name>`` work in regex mode)."""
        pattern = self.pattern()
        if self.mode == "whole":
            return replacement if pattern.fullmatch(text) else text
        if self.mode == "regex":
            return pattern.sub(replacement, text)
        return pattern.sub(lambda _m: replacement, text)  # literal: no escapes in replacement


def searchable_columns(table: Table, columns: Sequence[str] | None = None) -> list[str]:
    names = table.column_names if columns is None else [c for c in columns if c in table]
    return [name for name in names if table.kind(name) not in (STRUCTURE, REACTION)]


def match_mask(table: Table, query: Query, columns: Sequence[str] | None = None) -> np.ndarray:
    """``bool[rows, columns]`` in storage order for the searchable *columns*."""
    names = searchable_columns(table, columns)
    if not query.text or not names:
        return np.zeros((len(table), len(names)), dtype=np.bool_)
    predicate = query.predicate()
    return np.column_stack([table.column(name).match_text(predicate) for name in names])


@dataclass(frozen=True)
class Matches:
    """Matching cells in display order (row by row, left to right)."""

    oids: np.ndarray  # int64, one per match
    columns: tuple[str, ...]  # the searched columns
    column_index: np.ndarray  # index into *columns*, one per match

    def __len__(self) -> int:
        return len(self.oids)

    def cell(self, i: int) -> tuple[int, str]:
        return int(self.oids[i]), self.columns[int(self.column_index[i])]


def find(
    table: Table,
    query: Query,
    columns: Sequence[str] | None = None,
    *,
    shown_only: bool = True,
) -> Matches:
    """All matching cells, in the order the table shows them (hidden rows last if included)."""
    names = tuple(searchable_columns(table, columns))
    mask = match_mask(table, query, names)
    if shown_only:
        slots = table.view()
    else:
        shown = table.view()
        hidden = np.setdiff1d(np.arange(len(table)), shown, assume_unique=True)
        slots = np.concatenate([shown, hidden])
    rows, cols = np.nonzero(mask[slots]) if len(names) else (np.empty(0, int), np.empty(0, int))
    return Matches(table.oids[slots[rows]], names, cols.astype(np.int64))


def replace_all(
    table: Table,
    query: Query,
    replacement: str,
    columns: Sequence[str] | None = None,
    *,
    shown_only: bool = True,
) -> int:
    """Replace in every matching cell as one undoable step; returns how many cells changed."""
    if not query.text:
        return 0
    matches = find(table, query, columns, shown_only=shown_only)
    updates: dict[str, tuple[np.ndarray, list[str]]] = {}
    changed = 0
    for index, name in enumerate(matches.columns):
        oids = matches.oids[matches.column_index == index]
        if not len(oids):
            continue
        old = table.column(name).displays(table.slots_of(oids))
        new = [query.replace(text, replacement) for text in old]
        keep = np.array([a != b for a, b in zip(old, new, strict=True)], dtype=np.bool_)
        if keep.any():
            updates[name] = (oids[keep], [t for t, k in zip(new, keep, strict=True) if k])
            changed += int(keep.sum())
    if updates:
        table.set_block(updates, label=f"Replace “{query.text}”")
    return changed


def replace_cell(table: Table, oid: int, column: str, query: Query, replacement: str) -> bool:
    """Replace within one cell (Replace button); returns whether the cell changed."""
    old = table.cell(column, oid)
    new = query.replace(old, replacement)
    if new == old:
        return False
    table.set_cells(column, [oid], [new], label=f"Replace “{query.text}”")
    return True
