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

"""Table-core benchmarks at the target size (docs/PLAN.md §2.6). Run with ``pytest -m bench -s``.

Timings are the best of several runs so a noisy CI machine does not flap the gate; the limits
leave headroom over what a laptop measures.
"""

from __future__ import annotations

import time
from collections.abc import Callable

import numpy as np
import pytest

from chemsuite.core.table import NumberColumn, RangeFilter, Table, TextColumn, TextFilter

ROWS = 200_000
NUMBER_COLUMNS = 40
TEXT_COLUMNS = 10

pytestmark = pytest.mark.bench


def best_of(fn: Callable[[], object], repeat: int = 5) -> float:
    best = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - start)
    return best


@pytest.fixture(scope="module")
def big_table() -> Table:
    rng = np.random.default_rng(0)
    table = Table()
    table.add_rows(count=ROWS, undoable=False)
    for i in range(NUMBER_COLUMNS):
        values = rng.normal(5, 2, ROWS)
        values[rng.random(ROWS) < 0.05] = np.nan  # 5% blanks
        table._insert_column(f"num{i}", len(table.column_names), NumberColumn.from_floats(values))
    series = [f"Series-{k}" for k in range(40)]
    for i in range(TEXT_COLUMNS):
        if i == 0:
            texts = [f"CPD-{k:06d}" for k in rng.permutation(ROWS)]  # unique IDs
        else:
            texts = [series[k] for k in rng.integers(0, len(series), ROWS)]
        table._insert_column(f"text{i}", len(table.column_names), TextColumn.from_texts(texts))
    return table


def report(name: str, seconds: float, limit: float) -> None:
    print(f"\n  {name:<34} {seconds * 1000:8.1f} ms   (limit {limit * 1000:.0f} ms)")
    assert seconds < limit, f"{name}: {seconds * 1000:.1f} ms exceeds {limit * 1000:.0f} ms"


def test_memory_footprint(big_table: Table) -> None:
    mb = big_table.nbytes / 1024**2
    print(f"\n  table memory: {mb:.1f} MB for {ROWS:,} rows x {NUMBER_COLUMNS + TEXT_COLUMNS} cols")
    assert mb < 150


def test_sort_numeric(big_table: Table) -> None:
    def run() -> None:
        big_table.set_sort([("num0", False)])
        big_table.view()

    report("sort numeric", best_of(run), 0.050)


def test_sort_text_repeated_values(big_table: Table) -> None:
    def run() -> None:
        big_table.set_sort([("text1", True)])
        big_table.view()

    report("sort text (40 distinct)", best_of(run), 0.150)


def test_sort_text_unique_ids(big_table: Table) -> None:
    big_table.set_sort([("text0", True)])
    big_table.view()  # first sort builds the natural-order ranks once

    def run() -> None:
        big_table.set_sort([("text0", False)])
        big_table.view()

    report("sort text (unique IDs, cached)", best_of(run), 0.150)


def test_multi_column_sort(big_table: Table) -> None:
    def run() -> None:
        big_table.set_sort([("text1", True), ("num1", False)])
        big_table.view()

    report("sort text + numeric", best_of(run), 0.150)


def test_range_filter(big_table: Table) -> None:
    big_table.set_sort([])
    lows = iter(np.linspace(0, 8, 100))

    def run() -> None:
        big_table.set_filter("r", RangeFilter("num2", low=float(next(lows))))
        big_table.view()

    report("range filter", best_of(run), 0.030)
    big_table.clear_filters()


def test_text_filter_on_unique_ids(big_table: Table) -> None:
    def run() -> None:
        big_table.set_filter("t", TextFilter("text0", "0042"))
        big_table.view()

    report("text filter (unique IDs)", best_of(run), 0.150)
    big_table.clear_filters()


def test_edit_cells_and_undo(big_table: Table) -> None:
    oids = big_table.oids[::10]
    values = list(range(len(oids)))

    def run() -> None:
        big_table.set_cells("num3", oids, values)
        big_table.undo()

    report("set 20k cells + undo", best_of(run, 3), 0.300)


def test_delete_rows_and_undo(big_table: Table) -> None:
    oids = big_table.oids[::4]

    def run() -> None:
        big_table.remove_rows(oids)
        big_table.undo()

    report("delete 50k rows + undo", best_of(run, 3), 0.500)
