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

from __future__ import annotations

import time

import numpy as np
import pytest

from chemsuite.core.table import Table
from chemsuite.plotting.data import BLANK_LABEL, OTHER_LABEL, categories, numbers, xy
from chemsuite.plotting.geometry import box_mask, lasso_mask, nearest
from chemsuite.plotting.stats import bin_of, histogram_edges, linear_fit
from chemsuite.plotting.style import marker_sizes, scale_bins, scale_colors


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "x": ["1", "2", "", "4", "5"],
            "y": ["2", "4", "6", "x", "10"],
            "series": ["a", "b", "a", "", "c"],
        },
        undoable=False,
    )
    return t


def test_xy_drops_rows_missing_a_coordinate(table: Table) -> None:
    data = xy(table, "x", "y", table.oids)
    assert data.x.tolist() == [1, 2, 5] and data.y.tolist() == [2, 4, 10]
    assert data.missing == 2
    assert np.isnan(numbers(table, "y", table.oids)[3])


def test_categories_by_frequency_with_other_and_blank(table: Table) -> None:
    cats = categories(table, "series", table.oids, limit=1)
    assert cats.labels == ["a", OTHER_LABEL, BLANK_LABEL]
    assert cats.codes.tolist() == [0, 1, 0, 2, 1]
    assert cats.counts == [2, 2, 1]


def test_box_and_lasso() -> None:
    x = np.array([0.5, 1.5, 0.9, -1.0])
    y = np.array([0.5, 0.5, 0.9, 0.0])
    assert box_mask(x, y, 1, 0, 0, 1).tolist() == [True, False, True, False]
    triangle = np.array([[0, 0], [1, 0], [0, 1]])
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    assert lasso_mask(x, y, square).tolist() == [True, False, True, False]
    assert lasso_mask(np.array([0.2]), np.array([0.2]), triangle).tolist() == [True]
    assert lasso_mask(np.array([0.8]), np.array([0.8]), triangle).tolist() == [False]
    assert not lasso_mask(x, y, square[:2]).any()  # not a polygon


def test_lasso_is_fast_enough() -> None:
    rng = np.random.default_rng(0)
    x, y = rng.random(200_000), rng.random(200_000)
    angles = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    circle = np.column_stack([0.5 + 0.4 * np.cos(angles), 0.5 + 0.4 * np.sin(angles)])
    start = time.perf_counter()
    inside = lasso_mask(x, y, circle)
    assert time.perf_counter() - start < 1.0
    assert inside.mean() == pytest.approx(np.pi * 0.16, abs=0.01)


def test_nearest_uses_screen_distance() -> None:
    x, y = np.array([0.0, 10.0]), np.array([0.0, 0.1])
    # 1 data unit per pixel in x, 0.01 in y: from (7, 0.1) the first point is 12 px away
    # (7 across, 10 up) and the second only 3 px, though in data units y barely differs.
    assert nearest(x, y, 7.0, 0.1, 1.0, 0.01, radius=12) == 1
    assert nearest(x, y, 5.0, 0.5, 1.0, 0.01, radius=5) is None


def test_histogram_edges_and_bins() -> None:
    values = np.array([0.0, 1.0, 2.0, 3.0, np.nan])
    edges = histogram_edges(values, 3)
    assert edges.tolist() == [0.0, 1.0, 2.0, 3.0]
    assert bin_of(values, edges).tolist() == [0, 1, 2, 2, -1]
    assert histogram_edges(np.array([5.0, 5.0])).tolist() == [4.5, 5.5]
    assert len(histogram_edges(np.random.default_rng(1).random(10**6))) - 1 <= 100


def test_linear_fit() -> None:
    x = np.arange(10, dtype=float)
    fit = linear_fit(x, 2 * x + 1)
    assert fit.slope == pytest.approx(2) and fit.intercept == pytest.approx(1)
    assert fit.r2 == pytest.approx(1) and fit.n == 10
    assert "R² = 1.000" in fit.describe()
    assert linear_fit(np.ones(5), x[:5]) is None


def test_scales_and_sizes() -> None:
    colors = scale_colors(5)
    assert colors[0] == (68, 1, 84) and colors[-1] == (253, 231, 37)
    bins, lo, hi = scale_bins(np.array([0.0, 5.0, 10.0, np.nan]), 4)
    assert bins.tolist() == [0, 2, 3, -1] and (lo, hi) == (0.0, 10.0)
    sizes = marker_sizes(np.array([1.0, 3.0, np.nan]), 4, 16)
    assert sizes.tolist() == [4, 16, 4]


def test_box_stats_and_jitter() -> None:
    from chemsuite.plotting.geometry import jitter
    from chemsuite.plotting.stats import box_stats

    values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 100, np.nan], dtype=float)
    box = box_stats(values)
    assert (box.n, box.median, box.q1, box.q3) == (9, 5.0, 3.0, 7.0)
    assert (box.low, box.high) == (1.0, 8.0)  # 100 lies beyond 1.5 IQR: an outlier
    assert box_stats(np.array([np.nan])) is None
    offsets = jitter(np.arange(1000), 0.5)
    assert offsets.min() >= -0.25 and offsets.max() <= 0.25
    assert np.array_equal(offsets, jitter(np.arange(1000), 0.5))
    assert abs(offsets.mean()) < 0.03
