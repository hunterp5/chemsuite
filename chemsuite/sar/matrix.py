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
"""The SAR matrix: activity pivoted on two substituent positions (a third makes pages).

Each cell holds the compounds with that pair of substituents. A cell nobody has made shows
the Free-Wilson prediction when the additive model can estimate both groups.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .free_wilson import FreeWilson, free_wilson
from .rgroups import HYDROGEN, Decomposition

AGGREGATES = {"mean": "Mean", "median": "Median", "max": "Maximum", "min": "Minimum"}
MAX_LABELS = 60
BLANKS = frozenset({"", "n/a", "na", "none", "nan"})


@dataclass
class Cell:
    rows: list[int] = field(default_factory=list)  # input row indices
    values: list[float] = field(default_factory=list)  # their finite activities

    def value(self, aggregate: str) -> float:
        if not self.values:
            return math.nan
        return float(getattr(np, aggregate)(self.values))


@dataclass
class SarMatrix:
    row_labels: list[str]
    column_labels: list[str]
    page_labels: list[str]  # [""] without a page axis
    cells: dict[tuple[str, str, str], Cell]  # (page, row, column)
    model: FreeWilson | None
    dropped: int  # rows left out: blank labels, or past the label limit

    def cell(self, page: str, row: str, column: str) -> Cell | None:
        return self.cells.get((page, row, column))

    def predicted(self, page: str, row: str, column: str) -> float:
        if self.model is None:
            return math.nan
        choice = [row, column] + ([page] if len(self.model.positions) == 3 else [])
        return self.model.predict(choice)

    def values(self, aggregate: str) -> np.ndarray:
        return np.array([c.value(aggregate) for c in self.cells.values()])


def _order(labels: Counter[str]) -> list[str]:
    """Hydrogen first, then the most frequent."""
    return sorted(labels, key=lambda lab: (lab != HYDROGEN, -labels[lab], lab))


def build_matrix(
    rows: Sequence[str],
    columns: Sequence[str],
    activity: np.ndarray,
    pages: Sequence[str] | None = None,
    *,
    predict: bool = True,
    max_labels: int = MAX_LABELS,
) -> SarMatrix:
    """Pivot one label per input row on each axis (blank labels are left out)."""
    n = len(rows)
    page_values = list(pages) if pages is not None else [""] * n
    keep = [
        i
        for i in range(n)
        if rows[i].strip().lower() not in BLANKS
        and columns[i].strip().lower() not in BLANKS
        and (pages is None or page_values[i].strip().lower() not in BLANKS)
    ]
    row_counts = Counter(rows[i] for i in keep)
    column_counts = Counter(columns[i] for i in keep)
    page_counts = Counter(page_values[i] for i in keep)
    for counts, name in ((row_counts, "rows"), (column_counts, "columns")):
        if len(counts) > max_labels:
            raise ValueError(
                f"{len(counts)} different substituents on the {name} axis; the matrix shows at "
                f"most {max_labels}. Filter to a series first."
            )
    cells: dict[tuple[str, str, str], Cell] = defaultdict(Cell)
    activity = np.asarray(activity, dtype=np.float64)
    for i in keep:
        cell = cells[(page_values[i], rows[i], columns[i])]
        cell.rows.append(i)
        if math.isfinite(activity[i]):
            cell.values.append(float(activity[i]))
    model = None
    if predict:
        positions = ("rows", "columns") + (("pages",) if pages is not None else ())
        labels = {"rows": list(rows), "columns": list(columns), "pages": page_values}
        matched = np.zeros(n, dtype=np.bool_)
        matched[keep] = True
        decomposition = Decomposition(positions, {p: labels[p] for p in positions}, matched, "")
        try:
            model = free_wilson(decomposition, activity, positions=positions, min_count=1)
        except ValueError:
            model = None
    return SarMatrix(
        _order(row_counts),
        _order(column_counts),
        _order(page_counts) if pages is not None else [""],
        dict(cells),
        model,
        n - len(keep),
    )
