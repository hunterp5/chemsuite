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

"""Row filters. Each filter is an immutable description that produces a boolean mask.

Filters are evaluated with whole-column numpy operations; text predicates run once per unique
value of a dictionary-encoded column rather than once per row. The table combines all active
filters with AND (see :meth:`Table.view`).
"""

from __future__ import annotations

import dataclasses
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np

from .columns import Column, NumberColumn, TextColumn, to_number
from .search import Query, match_mask

if TYPE_CHECKING:
    from .table import Table


class Filter(ABC):
    @property
    @abstractmethod
    def columns(self) -> tuple[str, ...]:
        """Columns this filter reads (renames are applied, removals drop the filter)."""

    @abstractmethod
    def mask(self, table: Table) -> np.ndarray:
        """``bool`` per row (physical order): True keeps the row."""

    def matches(self, table: Table) -> np.ndarray:
        """The mask, inverted when the filter has ``negate`` set."""
        keep = self.mask(table)
        return ~keep if getattr(self, "negate", False) else keep

    def renamed(self, old: str, new: str) -> Filter:
        if getattr(self, "column", None) == old:
            return dataclasses.replace(self, column=new)  # type: ignore[type-var]
        return self


@dataclass(frozen=True)
class RangeFilter(Filter):
    """Keep rows whose number lies in ``[low, high]`` (either bound may be open).

    ``negate`` keeps numbers *outside* the range; blanks are decided by ``include_blank`` in
    both cases ("not 6-8, no blanks" means numbers below 6 or above 8).
    """

    column: str
    low: float | None = None
    high: float | None = None
    include_blank: bool = False
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return (self.column,)

    def mask(self, table: Table) -> np.ndarray:
        column = table.column(self.column)
        if not isinstance(column, NumberColumn):
            column = to_number(column, force=True)[0]  # non-numbers count as blank
        values = column.values()
        numbers = ~np.isnan(values)
        with np.errstate(invalid="ignore"):
            inside = numbers.copy()
            if self.low is not None:
                inside &= values >= self.low
            if self.high is not None:
                inside &= values <= self.high
        keep = numbers & ~inside if self.negate else inside
        if self.include_blank:
            keep |= column.blank_mask()
        return keep

    def matches(self, table: Table) -> np.ndarray:
        return self.mask(table)  # negate is handled inside mask


TextMode = Literal["contains", "equals", "starts", "ends", "regex"]


@dataclass(frozen=True)
class TextFilter(Filter):
    """Keep rows whose displayed text matches *query*."""

    column: str
    query: str
    mode: TextMode = "contains"
    case_sensitive: bool = False
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return (self.column,)

    def mask(self, table: Table) -> np.ndarray:
        return table.column(self.column).match_text(self._predicate())

    def _predicate(self):
        if self.mode == "regex":
            pattern = re.compile(self.query, 0 if self.case_sensitive else re.IGNORECASE)
            return lambda text: pattern.search(text) is not None
        fold = (lambda s: s) if self.case_sensitive else str.casefold
        query = fold(self.query)
        tests = {
            "contains": lambda text: query in fold(text),
            "equals": lambda text: fold(text) == query,
            "starts": lambda text: fold(text).startswith(query),
            "ends": lambda text: fold(text).endswith(query),
        }
        return tests[self.mode]


@dataclass(frozen=True)
class CategoryFilter(Filter):
    """Keep rows whose displayed value is one of *values* (``""`` selects blanks)."""

    column: str
    values: frozenset[str]
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return (self.column,)

    def mask(self, table: Table) -> np.ndarray:
        return table.column(self.column).match_text(self.values.__contains__)


@dataclass(frozen=True)
class OidFilter(Filter):
    """Keep a fixed set of rows (substructure hits, a named list, a plot selection)."""

    oids: frozenset[int]
    label: str = ""
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return ()

    def mask(self, table: Table) -> np.ndarray:
        wanted = np.fromiter(self.oids, dtype=np.int64, count=len(self.oids))
        return np.isin(table.oids, wanted)


@dataclass(frozen=True)
class SearchFilter(Filter):
    """Keep rows where any text or number column matches *query* (the quick search box)."""

    query: Query
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return ()  # every searchable column at evaluation time; never dropped with a column

    def mask(self, table: Table) -> np.ndarray:
        if not self.query.text:
            return np.ones(len(table), dtype=np.bool_)
        return match_mask(table, self.query).any(axis=1)


# -- summaries for filter widgets ------------------------------------------------------------


def number_values(column: Column) -> np.ndarray:
    """float64 per row; NaN for blanks and for text that is not a plain number."""
    if not isinstance(column, NumberColumn):
        column = to_number(column, force=True)[0]
    return column.values()


def value_counts(column: Column) -> list[tuple[str, int]]:
    """Distinct displayed values with their row counts, most frequent first ("" = blank)."""
    if isinstance(column, TextColumn):
        pairs = column.value_counts()
    else:
        displays = column.displays(np.arange(len(column)))
        values, counts = np.unique(np.array(displays, dtype=object), return_counts=True)
        pairs = [(str(v), int(c)) for v, c in zip(values, counts, strict=True)]
    pairs.sort(key=lambda pair: (-pair[1], pair[0]))
    return pairs
