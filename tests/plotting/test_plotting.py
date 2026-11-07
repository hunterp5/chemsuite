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


def test_distribution_fits() -> None:
    from chemsuite.plotting.stats import fit_distribution

    rng = np.random.default_rng(7)
    normal = fit_distribution(rng.normal(5, 2, 20_000), "normal")
    assert normal.mu == pytest.approx(5, abs=0.05) and normal.sigma == pytest.approx(2, abs=0.05)
    edges = np.linspace(-5, 15, 41)
    xs, counts = normal.curve(edges)
    # The curve's area matches the number of values (counts per bin of width 0.5).
    assert np.trapezoid(counts, xs) / 0.5 == pytest.approx(20_000, rel=0.01)
    lognormal = fit_distribution(rng.lognormal(1.0, 0.5, 20_000), "lognormal")
    assert lognormal.mu == pytest.approx(1.0, abs=0.02)
    assert "median = 2.7" in lognormal.describe()
    assert lognormal.pdf(np.array([-1.0, 0.0])).tolist() == [0.0, 0.0]
    assert fit_distribution(np.array([1.0, 1.0, 1.0]), "normal") is None
    assert fit_distribution(np.array([-1.0, -2.0, 3.0]), "lognormal") is None
    with pytest.raises(ValueError):
        fit_distribution(np.ones(5), "gamma")


def test_aggregate_per_group() -> None:
    from chemsuite.plotting.stats import aggregate

    values = np.array([1.0, 3.0, np.nan, 10.0, 2.0])
    codes = np.array([0, 0, 0, 1, -1])
    mean = aggregate(values, codes, 3, "mean")
    assert mean.value[:2].tolist() == [2.0, 10.0] and np.isnan(mean.value[2])
    assert mean.spread[0] == pytest.approx(np.sqrt(2)) and np.isnan(mean.spread[1])
    assert mean.count.tolist() == [2, 1, 0]
    assert aggregate(None, codes, 3, "count").value.tolist() == [3.0, 1.0, 0.0]
    assert aggregate(values, codes, 2, "sum").value.tolist() == [4.0, 10.0]
    with pytest.raises(ValueError):
        aggregate(values, codes, 2, "mode")


def test_density_grid() -> None:
    from chemsuite.plotting.stats import density_grid

    x = np.array([0.0, 0.1, 0.9, 1.0, np.nan])
    y = np.array([0.0, 0.0, 1.0, 1.0, 5.0])
    grid = density_grid(x, y, bins=2)
    assert grid.counts.tolist() == [[2.0, 0.0], [0.0, 2.0]]
    assert grid.cell[:4].tolist() == [[0, 0], [0, 0], [1, 1], [1, 1]]
    assert grid.cell[4].tolist() == [-1, -1]
    assert grid.cell_at(1.0, 1.0) == (1, 1) and grid.cell_at(2.0, 0.0) is None
    assert density_grid(np.array([np.nan]), np.array([1.0])) is None


def test_unit_scale() -> None:
    from chemsuite.plotting.stats import unit_scale

    scaled, lows, highs = unit_scale([np.array([1.0, 3.0, np.nan]), np.array([5.0, 5.0, 5.0])])
    assert scaled[:, 0][:2].tolist() == [0.0, 1.0] and np.isnan(scaled[2, 0])
    assert scaled[:, 1].tolist() == [0.5, 0.5, 0.5]
    assert lows.tolist() == [1.0, 5.0] and highs.tolist() == [3.0, 5.0]


def test_mean_by_x() -> None:
    from chemsuite.plotting.stats import mean_by_x

    xs, ys = mean_by_x(np.array([2.0, 1.0, 2.0, np.nan]), np.array([4.0, 1.0, 6.0, 9.0]))
    assert xs.tolist() == [1.0, 2.0] and ys.tolist() == [1.0, 5.0]


def test_project_to_the_screen() -> None:
    from chemsuite.plotting.geometry import project

    # a camera 2 units back looking down -z, 90° field of view
    near, far = 0.1, 100.0
    perspective = np.array(
        [
            [1.0, 0, 0, 0],
            [0, 1.0, 0, 0],
            [0, 0, -(far + near) / (far - near), -2 * far * near / (far - near)],
            [0, 0, -1.0, 0],
        ]
    )
    view = np.eye(4)
    view[2, 3] = -2.0
    points = np.array([[0.0, 0, 0], [1.0, 1.0, 0], [0, 0, 5.0]])
    xy, in_front = project(points, perspective @ view, 200, 100)
    assert xy[0].tolist() == [100.0, 50.0]  # the centre
    assert xy[1].tolist() == pytest.approx([150.0, 25.0])  # right and up
    assert in_front.tolist() == [True, True, False]  # the last is behind the camera


def test_activity_landscape_surface() -> None:
    from chemsuite.plotting.stats import landscape

    rng = np.random.default_rng(0)
    # two islands in the plane: one active (z ≈ 8), one inactive (z ≈ 5), far apart
    x = np.concatenate([rng.normal(0, 0.3, 200), rng.normal(10, 0.3, 200), [np.nan]])
    y = np.concatenate([rng.normal(0, 0.3, 200), rng.normal(10, 0.3, 200), [1.0]])
    z = np.concatenate([np.full(200, 8.0), np.full(200, 5.0), [9.0]])
    surface = landscape(x, y, z, bins=20)
    assert surface.heights.shape == (20, 20) and len(surface.xs) == len(surface.ys) == 20
    assert surface.heights[0, 0] == pytest.approx(8.0)  # a weighted mean, not a sum
    assert surface.heights[-1, -1] == pytest.approx(5.0)
    assert np.isnan(surface.heights[10, 10])  # nothing near the middle: no height made up
    finite = surface.heights[np.isfinite(surface.heights)]
    assert finite.min() >= 5.0 - 1e-9 and finite.max() <= 8.0 + 1e-9
    assert landscape(x[:2], y[:2], z[:2]) is None
