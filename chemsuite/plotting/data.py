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
"""Turning table columns into plot data.

Everything is keyed by row oid so the GUI can map a hit point back to a row (and to the shared
selection) no matter how the table is sorted or filtered.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER, REACTION, STRUCTURE, TEXT, TextColumn
from chemsuite.core.table.filters import number_values

BLANK_LABEL = "(blank)"
OTHER_LABEL = "Other"
MAX_CATEGORIES = 12


def number_columns(table: Table) -> list[str]:
    return [n for n in table.column_names if table.kind(n) == NUMBER]


def colorable_columns(table: Table) -> list[str]:
    """Columns that can colour points: numbers (a scale) and text (categories)."""
    return [n for n in table.column_names if table.kind(n) in (NUMBER, TEXT)]


def numbers(table: Table, name: str, oids: np.ndarray) -> np.ndarray:
    """float64 per oid; NaN where the cell is blank or not a number."""
    return number_values(table.column(name))[table.slots_of(oids)]


@dataclass(frozen=True)
class Categories:
    codes: np.ndarray  # int per point, index into labels
    labels: list[str]
    counts: list[int]


def categories(
    table: Table, name: str, oids: np.ndarray, *, limit: int = MAX_CATEGORIES
) -> Categories:
    """The most frequent *limit* values keep their own code; the rest share "Other".

    Labels are ordered by frequency; blank cells get their own "(blank)" category, listed last.
    """
    column = table.column(name)
    slots = table.slots_of(oids)
    if isinstance(column, TextColumn):  # group by dictionary code: no string per row
        codes, strings = column.coded(slots)
        counts_all = np.bincount(codes, minlength=len(strings))
        present = np.flatnonzero(counts_all)
        values = np.array([strings[i] for i in present], dtype=object)
        counts = counts_all[present]
        lookup = np.full(len(strings), -1, dtype=np.int64)
        lookup[present] = np.arange(len(present))
        inverse = lookup[codes]
    else:
        texts = np.array(column.displays(slots), dtype=object)
        values, inverse, counts = np.unique(texts, return_inverse=True, return_counts=True)
    texts = inverse  # only its length is used below
    order = sorted(range(len(values)), key=lambda i: (values[i] == "", -counts[i], values[i]))
    keep = [i for i in order if values[i] != ""][:limit]
    remap = np.full(len(values), -1, dtype=np.int64)
    labels: list[str] = []
    for code, i in enumerate(keep):
        remap[i] = code
        labels.append(str(values[i]))
    others = [i for i in order if values[i] != "" and remap[i] < 0]
    if others:
        remap[others] = len(labels)
        labels.append(OTHER_LABEL)
    blank = np.flatnonzero(values == "")
    if len(blank):
        remap[blank] = len(labels)
        labels.append(BLANK_LABEL)
    codes = remap[inverse] if len(texts) else np.empty(0, dtype=np.int64)
    counts_out = np.bincount(codes, minlength=len(labels)).tolist() if len(codes) else []
    return Categories(codes, labels, counts_out)


@dataclass(frozen=True)
class XY:
    """Points with both coordinates present, plus how many rows were left out."""

    oids: np.ndarray
    x: np.ndarray
    y: np.ndarray
    missing: int

    def __len__(self) -> int:
        return len(self.oids)


def xy(table: Table, x: str, y: str, oids: np.ndarray) -> XY:
    xs, ys = numbers(table, x, oids), numbers(table, y, oids)
    ok = np.isfinite(xs) & np.isfinite(ys)
    return XY(oids[ok], xs[ok], ys[ok], int((~ok).sum()))


def shown_oids(table: Table) -> np.ndarray:
    """Rows the table shows now (filters applied), in display order."""
    return table.oids[table.view()]


def scope_oids(table: Table, selected: np.ndarray, scope: str) -> np.ndarray:
    """The rows in *scope*: ``"selected"`` (*selected*), ``"shown"`` (passing the filters, in
    display order) or ``"all"`` (storage order)."""
    if scope == "selected":
        return selected
    if scope == "shown":
        return shown_oids(table)
    return table.oids.copy()


def is_plottable_axis(table: Table, name: str | None) -> bool:
    return bool(name) and name in table and table.kind(name) not in (STRUCTURE, REACTION)
