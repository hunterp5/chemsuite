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
"""Advanced search: rows of conditions combined with "all" or "any".

Each :class:`Criterion` is one condition (column, operator, value). Text and number conditions
are evaluated here with whole-column operations. Structure conditions need RDKit and may take a
while, so the app runs them as background jobs and passes their matching rows in
(``structure_hits``); this module stays free of chemistry.
"""

from __future__ import annotations

import dataclasses
import math
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from .columns import NUMBER, STRUCTURE
from .filters import Filter, number_values

if TYPE_CHECKING:
    from .table import Table

ANY_COLUMN = ""  # Criterion.column for "any column"
Combine = Literal["all", "any"]

TEXT_OPS: dict[str, str] = {
    "contains": "contains",
    "equals": "equals",
    "starts": "starts with",
    "ends": "ends with",
    "regex": "matches regex",
    "blank": "is blank",
}
NUMBER_OPS: dict[str, str] = {
    "eq": "=",
    "ne": "≠",
    "lt": "<",
    "le": "≤",
    "gt": ">",
    "ge": "≥",
    "between": "between",
    "blank": "is blank",
}
STRUCTURE_OPS: dict[str, str] = {
    "substructure": "contains substructure",
    "exact": "is exactly",
    "blank": "is blank",
}
#: operators whose value is a structure (evaluated by the app, not here)
STRUCTURE_QUERY_OPS = frozenset({"substructure", "exact"})


def ops_for(kind: str | None) -> dict[str, str]:
    """Operators offered for a column kind (None: any column)."""
    if kind == NUMBER:
        return NUMBER_OPS
    if kind == STRUCTURE:
        return STRUCTURE_OPS
    return TEXT_OPS


class CriterionError(ValueError):
    """A condition that cannot be evaluated (a bad number or regular expression)."""


@dataclass(frozen=True)
class Criterion:
    column: str  # ANY_COLUMN for any text or number column
    op: str
    value: str = ""
    value2: str = ""  # the upper bound for "between"
    case_sensitive: bool = False
    negate: bool = False

    @property
    def needs_structure_search(self) -> bool:
        return self.op in STRUCTURE_QUERY_OPS

    def check(self, table: Table) -> None:
        """Raise :class:`CriterionError` if this condition cannot be evaluated on *table*."""
        if self.column != ANY_COLUMN and self.column not in table:
            raise CriterionError(f"There is no column “{self.column}”.")
        kind = None if self.column == ANY_COLUMN else table.kind(self.column)
        if self.op not in ops_for(kind):
            raise CriterionError(f"“{self.op}” does not apply to this column.")
        if self.op == "blank":
            return
        if self.needs_structure_search:
            if not self.value.strip():
                raise CriterionError("Enter a structure.")
            return
        if self.op in NUMBER_OPS and kind == NUMBER:
            self._number(self.value)
            if self.op == "between":
                self._number(self.value2)
            return
        if self.op == "regex":
            try:
                re.compile(self.value)
            except re.error as err:
                raise CriterionError(f"Invalid regular expression: {err}") from err
        elif not self.value:
            raise CriterionError("Enter some text.")

    @staticmethod
    def _number(text: str) -> float:
        try:
            value = float(text.strip())
        except ValueError:
            raise CriterionError(f"“{text}” is not a number.") from None
        if not math.isfinite(value):
            raise CriterionError(f"“{text}” is not a number.")
        return value

    # -- evaluating ---------------------------------------------------------------------------
    def mask(self, table: Table, structure_hits: np.ndarray | None = None) -> np.ndarray:
        """``bool`` per row (storage order), before ``negate``.

        Structure conditions use *structure_hits* (the oids that matched), found by the app.
        """
        self.check(table)
        n = len(table)
        if self.needs_structure_search:
            if structure_hits is None:
                raise CriterionError("The structure search has not run.")
            return np.isin(table.oids, structure_hits)
        columns = self._columns(table)
        out = np.zeros(n, dtype=np.bool_)
        for name in columns:
            out |= self._column_mask(table, name)
        return out

    def matches(self, table: Table, structure_hits: np.ndarray | None = None) -> np.ndarray:
        keep = self.mask(table, structure_hits)
        return ~keep if self.negate else keep

    def _columns(self, table: Table) -> list[str]:
        if self.column != ANY_COLUMN:
            return [self.column]
        return [n for n in table.column_names if table.kind(n) != STRUCTURE]

    def _column_mask(self, table: Table, name: str) -> np.ndarray:
        column = table.column(name)
        if self.op == "blank":
            return column.blank_mask()
        if table.kind(name) == NUMBER and self.op in NUMBER_OPS:
            return self._number_mask(number_values(column))
        return column.match_text(self._text_predicate())

    def _number_mask(self, values: np.ndarray) -> np.ndarray:
        a = self._number(self.value)
        with np.errstate(invalid="ignore"):
            if self.op == "eq":
                out = np.isclose(values, a, rtol=1e-12, atol=0.0)
            elif self.op == "ne":
                out = ~np.isclose(values, a, rtol=1e-12, atol=0.0)
            elif self.op == "lt":
                out = values < a
            elif self.op == "le":
                out = values <= a
            elif self.op == "gt":
                out = values > a
            elif self.op == "ge":
                out = values >= a
            else:  # between, in either order
                b = self._number(self.value2)
                low, high = min(a, b), max(a, b)
                out = (values >= low) & (values <= high)
        return out & ~np.isnan(values)  # a blank is never a number match

    def _text_predicate(self) -> Callable[[str], bool]:
        if self.op == "regex":
            pattern = re.compile(self.value, 0 if self.case_sensitive else re.IGNORECASE)
            return lambda text: pattern.search(text) is not None
        fold = (lambda s: s) if self.case_sensitive else str.casefold
        query = fold(self.value)
        tests = {
            "contains": lambda text: query in fold(text),
            "equals": lambda text: fold(text) == query,
            "starts": lambda text: fold(text).startswith(query),
            "ends": lambda text: fold(text).endswith(query),
        }
        return tests[self.op]

    # -- describing and saving ----------------------------------------------------------------
    def describe(self) -> str:
        column = self.column or "any column"
        kind_ops = {**TEXT_OPS, **NUMBER_OPS, **STRUCTURE_OPS}
        op = kind_ops.get(self.op, self.op)
        if self.op == "blank":
            text = f"{column} {op}"
        elif self.op == "between":
            text = f"{column} between {self.value} and {self.value2}"
        else:
            value = self.value if len(self.value) <= 30 else self.value[:29] + "…"
            text = f"{column} {op} {value}"
        return f"not ({text})" if self.negate else text

    def to_json(self) -> dict[str, Any]:
        return dataclasses.asdict(self)

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> Criterion:
        names = {f.name for f in dataclasses.fields(cls)}
        values = {k: v for k, v in data.items() if k in names}
        if not isinstance(values.get("column"), str) or not isinstance(values.get("op"), str):
            raise ValueError(f"not a saved condition: {data!r}")
        return cls(**values)


@dataclass(frozen=True)
class CriteriaSearch:
    criteria: tuple[Criterion, ...]
    combine: Combine = "all"

    def check(self, table: Table) -> None:
        if not self.criteria:
            raise CriterionError("Add at least one condition.")
        for criterion in self.criteria:
            criterion.check(table)

    def mask(
        self, table: Table, structure_hits: Mapping[int, np.ndarray] | None = None
    ) -> np.ndarray:
        """``bool`` per row; *structure_hits* maps a criterion's index to its matching oids."""
        hits = structure_hits or {}
        masks = [c.matches(table, hits.get(i)) for i, c in enumerate(self.criteria)]
        if not masks:
            return np.ones(len(table), dtype=np.bool_)
        combined = np.logical_and if self.combine == "all" else np.logical_or
        return combined.reduce(masks)

    def describe(self) -> str:
        joiner = " and " if self.combine == "all" else " or "
        return joiner.join(c.describe() for c in self.criteria)

    def renamed(self, old: str, new: str) -> CriteriaSearch:
        criteria = tuple(
            dataclasses.replace(c, column=new) if c.column == old else c for c in self.criteria
        )
        return dataclasses.replace(self, criteria=criteria)

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(c.column for c in self.criteria if c.column != ANY_COLUMN))

    def to_json(self) -> dict[str, Any]:
        return {"combine": self.combine, "criteria": [c.to_json() for c in self.criteria]}

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> CriteriaSearch:
        combine = data.get("combine", "all")
        if combine not in ("all", "any"):
            raise ValueError(f"unknown combination {combine!r}")
        criteria = data.get("criteria")
        if not isinstance(criteria, Sequence):
            raise ValueError("a saved search needs a list of conditions")
        return cls(tuple(Criterion.from_json(c) for c in criteria), combine)


@dataclass(frozen=True)
class CriteriaFilter(Filter):
    """Keep the rows an advanced search matches.

    Structure conditions were searched when the filter was made: their matching oids are kept
    here (``structure_hits``), so editing a structure afterwards needs a new search.
    """

    search: CriteriaSearch
    structure_hits: tuple[tuple[int, frozenset[int]], ...] = field(default=())
    negate: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return self.search.columns

    def mask(self, table: Table) -> np.ndarray:
        hits = {
            i: np.fromiter(oids, dtype=np.int64, count=len(oids)) for i, oids in self.structure_hits
        }
        return self.search.mask(table, hits)

    def renamed(self, old: str, new: str) -> Filter:
        return dataclasses.replace(self, search=self.search.renamed(old, new))
