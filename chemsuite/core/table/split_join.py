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

"""Split one column into several on a delimiter, or join several columns into one.

Both work on the cells' displayed text and write their output as one undoable step. New
columns infer their kind from the values, so splitting ``"IC50; 5.2"`` gives a number column
for the second field.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from .table import Table

SPLIT_DELIMITERS = {"comma": ",", "semicolon": ";", "tab": "\t", "pipe": "|", "space": " "}
JOIN_DELIMITERS = {
    "comma": ",",
    "comma_space": ", ",
    "semicolon": "; ",
    "tab": "\t",
    "pipe": "|",
    "space": " ",
    "dot": ".",  # joins SMILES into one multi-component structure
    "none": "",
}
_AUTO = (",", ";", "\t", "|")
Keep = Literal["all", "first", "last", "largest", "smallest"]
MAX_SPLIT_COLUMNS = 256


def unescape(text: str) -> str:
    r"""A delimiter typed by the user: ``\t`` means a tab."""
    return text.replace("\\t", "\t").replace("\\n", "\n")


def detect_delimiter(texts: Iterable[str]) -> str:
    """The punctuation that appears most (comma, semicolon, tab, pipe), else whitespace."""
    counts = dict.fromkeys(_AUTO, 0)
    spaced = False
    for text in texts:
        for delimiter in _AUTO:
            counts[delimiter] += text.count(delimiter)
        spaced = spaced or len(text.split()) > 1
    best = max(counts, key=counts.__getitem__)
    if counts[best]:
        return best
    return " " if spaced else ","


def split_text(text: str, delimiter: str) -> list[str]:
    """Fields of one cell. Quotes protect delimiters (``"a, b", c``); ``" "`` splits on runs."""
    if not text.strip():
        return []
    if delimiter == " ":
        return text.split()
    if len(delimiter) != 1 or '"' not in text:
        return [part.strip() for part in text.split(delimiter)]
    try:
        row = next(csv.reader(io.StringIO(text), delimiter=delimiter, skipinitialspace=True))
    except (StopIteration, csv.Error):
        return [part.strip() for part in text.split(delimiter)]
    return [part.strip() for part in row]


def _number(text: str) -> float | None:
    try:
        value = float(text)
    except ValueError:
        return None
    return None if value != value else value  # NaN is not a usable number


def keep_fields(fields: list[str], keep: Keep) -> list[str]:
    """All fields, or the one *keep* picks (largest/smallest prefer numbers over text)."""
    present = [f for f in fields if f]
    if keep == "all" or not present:
        return fields if keep == "all" else [""]
    if keep == "first":
        return [present[0]]
    if keep == "last":
        return [present[-1]]
    numbers = [(n, f) for f in present if (n := _number(f)) is not None]
    pool = numbers or [(f.casefold(), f) for f in present]
    pick = max if keep == "largest" else min
    return [pick(pool, key=lambda pair: pair[0])[1]]


@dataclass(frozen=True)
class SplitResult:
    delimiter: str
    oids: np.ndarray
    fields: list[list[str]]  # one list per oid

    @property
    def width(self) -> int:
        return max((len(f) for f in self.fields), default=0)

    def column(self, i: int) -> list[str]:
        return [f[i] if i < len(f) else "" for f in self.fields]


def plan_split(
    table: Table,
    column: str,
    delimiter: str | None = None,
    *,
    keep: Keep = "all",
    oids: np.ndarray | None = None,
) -> SplitResult:
    """Split without writing (for previews). ``delimiter=None`` detects it from the cells."""
    oids = table.oids if oids is None else np.asarray(oids, dtype=np.int64)
    texts = table.column(column).displays(table.slots_of(oids))
    if delimiter is None:
        delimiter = detect_delimiter(dict.fromkeys(texts))
    if not delimiter:
        raise ValueError("Enter a delimiter.")
    memo: dict[str, list[str]] = {}  # dictionary-encoded columns repeat values a lot
    fields = []
    for text in texts:
        parts = memo.get(text)
        if parts is None:
            parts = memo[text] = keep_fields(split_text(text, delimiter), keep)
        fields.append(parts)
    return SplitResult(delimiter, oids, fields)


def split_names(prefix: str, width: int) -> list[str]:
    return [f"{prefix}_{i + 1}" for i in range(width)]


def split_column(
    table: Table,
    column: str,
    delimiter: str | None = None,
    *,
    prefix: str | None = None,
    keep: Keep = "all",
    oids: np.ndarray | None = None,
) -> list[str]:
    """Write the fields into new columns right of *column*; returns their names.

    Rows outside *oids* (when given) are left blank in the new columns.
    """
    result = plan_split(table, column, delimiter, keep=keep, oids=oids)
    width = result.width
    if not width:
        return []
    if width > MAX_SPLIT_COLUMNS:
        raise ValueError(
            f"Splitting on {result.delimiter!r} would make {width} columns "
            f"(the limit is {MAX_SPLIT_COLUMNS}). Check the delimiter."
        )
    prefix = (prefix or "").strip() or column
    requested = split_names(prefix, width) if keep == "all" else [prefix + "_" + keep]
    names = table.apply_results(
        {name: (result.oids, result.column(i)) for i, name in enumerate(requested)},
        label=f"Split {column}",
        position=table.column_names.index(column) + 1,
    )
    return list(names.values())


def join_texts(rows: Iterable[Sequence[str]], delimiter: str, *, skip_blank: bool) -> list[str]:
    joined = []
    for parts in rows:
        if skip_blank:
            parts = [p for p in parts if p.strip()]
        elif not any(p.strip() for p in parts):
            parts = []  # an all-blank row stays blank, not ",,"
        joined.append(delimiter.join(parts))
    return joined


def plan_join(
    table: Table,
    columns: Sequence[str],
    delimiter: str,
    *,
    skip_blank: bool = True,
    oids: np.ndarray | None = None,
) -> tuple[np.ndarray, list[str]]:
    oids = table.oids if oids is None else np.asarray(oids, dtype=np.int64)
    slots = table.slots_of(oids)
    cells = [table.column(name).displays(slots) for name in columns]
    return oids, join_texts(zip(*cells, strict=True), delimiter, skip_blank=skip_blank)


def join_columns(
    table: Table,
    columns: Sequence[str],
    delimiter: str,
    *,
    name: str | None = None,
    skip_blank: bool = True,
    oids: np.ndarray | None = None,
) -> str:
    """Write the joined text into a new column right of the last source; returns its name."""
    if len(columns) < 2:
        raise ValueError("Choose at least two columns to join.")
    oids, joined = plan_join(table, columns, delimiter, skip_blank=skip_blank, oids=oids)
    requested = (name or "").strip() or "_".join(columns)
    position = max(table.column_names.index(c) for c in columns) + 1
    names = table.apply_results(
        {requested: (oids, joined)}, label=f"Join {len(columns)} columns", position=position
    )
    return names[requested]


__all__ = [
    "JOIN_DELIMITERS",
    "MAX_SPLIT_COLUMNS",
    "SPLIT_DELIMITERS",
    "SplitResult",
    "detect_delimiter",
    "join_columns",
    "keep_fields",
    "plan_join",
    "plan_split",
    "split_column",
    "split_text",
    "unescape",
]
