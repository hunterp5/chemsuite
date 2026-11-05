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
"""Histogram binning and simple statistics for plot overlays."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MAX_BINS = 100


def histogram_edges(values: np.ndarray, bins: int | None = None) -> np.ndarray:
    """Bin edges for the finite *values*: *bins* equal bins, or numpy's "auto" (capped)."""
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.array([0.0, 1.0])
    lo, hi = float(finite.min()), float(finite.max())
    if lo == hi:
        return np.array([lo - 0.5, hi + 0.5])
    if bins is None:
        edges = np.histogram_bin_edges(finite, bins="auto")
        if len(edges) - 1 <= MAX_BINS:
            return edges
        bins = MAX_BINS
    return np.linspace(lo, hi, max(1, int(bins)) + 1)


def bin_of(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Bin index per value (the last bin includes its right edge); -1 outside or NaN."""
    index = np.searchsorted(edges, values, side="right") - 1
    index = np.where(values == edges[-1], len(edges) - 2, index)
    valid = np.isfinite(values) & (index >= 0) & (index < len(edges) - 1)
    return np.where(valid, index, -1)


@dataclass(frozen=True)
class LinearFit:
    slope: float
    intercept: float
    r: float
    n: int

    @property
    def r2(self) -> float:
        return self.r * self.r

    def describe(self) -> str:
        sign = "+" if self.intercept >= 0 else "−"
        return (
            f"y = {self.slope:.3g}x {sign} {abs(self.intercept):.3g} · "
            f"R² = {self.r2:.3f} · r = {self.r:.3f} · n = {self.n:,}"
        )


def linear_fit(x: np.ndarray, y: np.ndarray) -> LinearFit | None:
    """Least-squares line through the finite pairs; None for fewer than 3 or no spread."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or np.ptp(x) == 0:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    r = float(np.corrcoef(x, y)[0, 1]) if np.ptp(y) > 0 else 0.0
    return LinearFit(float(slope), float(intercept), r, len(x))


@dataclass(frozen=True)
class BoxStats:
    """Tukey box: quartiles, whiskers at the furthest values within 1.5 IQR of the box."""

    n: int
    q1: float
    median: float
    q3: float
    low: float
    high: float
    mean: float

    def describe(self) -> str:
        return f"n = {self.n:,} · median {self.median:.3g} · IQR {self.q1:.3g} to {self.q3:.3g}"


def box_stats(values: np.ndarray) -> BoxStats | None:
    finite = values[np.isfinite(values)]
    if not len(finite):
        return None
    q1, median, q3 = (float(v) for v in np.percentile(finite, [25, 50, 75]))
    reach = 1.5 * (q3 - q1)
    inside = finite[(finite >= q1 - reach) & (finite <= q3 + reach)]
    return BoxStats(
        len(finite), q1, median, q3, float(inside.min()), float(inside.max()), float(finite.mean())
    )
