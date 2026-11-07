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


CURVES = {"normal": "Normal", "lognormal": "Log-normal"}


@dataclass(frozen=True)
class DistributionFit:
    """A fitted density, drawn over a histogram scaled to its counts."""

    kind: str
    mu: float  # mean (normal) or mean of log values (log-normal)
    sigma: float  # standard deviation, of the values or of their logs
    n: int

    def pdf(self, x: np.ndarray) -> np.ndarray:
        if self.kind == "normal":
            z = (x - self.mu) / self.sigma
            return np.exp(-0.5 * z * z) / (self.sigma * np.sqrt(2 * np.pi))
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (np.log(x) - self.mu) / self.sigma
            density = np.exp(-0.5 * z * z) / (x * self.sigma * np.sqrt(2 * np.pi))
        return np.where(x > 0, density, 0.0)

    def curve(self, edges: np.ndarray, points: int = 200) -> tuple[np.ndarray, np.ndarray]:
        """``(x, expected count per bin)`` across the histogram (equal-width bins)."""
        lo, hi = float(edges[0]), float(edges[-1])
        if self.kind == "lognormal":
            lo = max(lo, np.finfo(float).tiny)
        xs = np.linspace(lo, hi, points)
        width = float(np.mean(np.diff(edges)))
        return xs, self.pdf(xs) * self.n * width

    def describe(self) -> str:
        if self.kind == "normal":
            return f"Normal fit: μ = {self.mu:.4g}, σ = {self.sigma:.4g}"
        return f"Log-normal fit: median = {np.exp(self.mu):.4g}, σ(ln) = {self.sigma:.4g}"


def fit_distribution(values: np.ndarray, kind: str) -> DistributionFit | None:
    """Maximum-likelihood fit (log-normal uses the positive values only); None if degenerate."""
    finite = values[np.isfinite(values)]
    if kind == "lognormal":
        finite = finite[finite > 0]
        sample = np.log(finite)
        sigma = float(sample.std()) if len(sample) > 1 else 0.0
    elif kind == "normal":
        sample = finite
        sigma = float(sample.std(ddof=1)) if len(sample) > 1 else 0.0
    else:
        raise ValueError(f"unknown curve {kind!r}")
    if len(sample) < 3 or sigma <= 1e-12:
        return None
    return DistributionFit(kind, float(sample.mean()), sigma, len(finite))


# -- bar charts ---------------------------------------------------------------------------------
#: how a bar chart sums up each group: key -> title
AGGREGATES = {
    "mean": "Mean",
    "median": "Median",
    "sum": "Sum",
    "min": "Minimum",
    "max": "Maximum",
    "count": "Count",
}


@dataclass(frozen=True)
class GroupValues:
    value: np.ndarray  # one per group (NaN: no numbers in the group)
    spread: np.ndarray  # standard deviation for "mean", else NaN
    count: np.ndarray  # numbers per group (rows per group for "count")


def aggregate(values: np.ndarray | None, codes: np.ndarray, groups: int, how: str) -> GroupValues:
    """Sum up *values* per group (*codes* index the groups; -1 is left out). For "count",
    *values* may be None and every row counts."""
    if how not in AGGREGATES:
        raise ValueError(f"unknown aggregate {how!r}")
    value = np.full(groups, np.nan)
    spread = np.full(groups, np.nan)
    count = np.zeros(groups, dtype=np.int64)
    for group in range(groups):
        member = codes == group
        if how == "count":
            count[group] = int(member.sum())
            value[group] = float(count[group])
            continue
        sample = values[member]
        sample = sample[np.isfinite(sample)]
        count[group] = len(sample)
        if not len(sample):
            continue
        value[group] = {
            "mean": np.mean,
            "median": np.median,
            "sum": np.sum,
            "min": np.min,
            "max": np.max,
        }[how](sample)
        if how == "mean" and len(sample) > 1:
            spread[group] = float(np.std(sample, ddof=1))
    return GroupValues(value, spread, count)


def mean_by_x(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The mean Y at each distinct X, in X order (finite pairs only)."""
    ok = np.isfinite(x) & np.isfinite(y)
    xs, inverse = np.unique(x[ok], return_inverse=True)
    if not len(xs):
        return xs, xs
    sums = np.bincount(inverse, weights=y[ok], minlength=len(xs))
    return xs, sums / np.bincount(inverse, minlength=len(xs))


# -- density ------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Grid:
    counts: np.ndarray  # (x bins, y bins)
    x_edges: np.ndarray
    y_edges: np.ndarray
    cell: np.ndarray  # (x bin, y bin) per point; -1 where the point has no finite value

    def cell_at(self, x: float, y: float) -> tuple[int, int] | None:
        i = int(np.searchsorted(self.x_edges, x, side="right")) - 1
        j = int(np.searchsorted(self.y_edges, y, side="right")) - 1
        nx, ny = self.counts.shape
        if x == self.x_edges[-1]:
            i = nx - 1
        if y == self.y_edges[-1]:
            j = ny - 1
        return (i, j) if 0 <= i < nx and 0 <= j < ny else None


def density_grid(x: np.ndarray, y: np.ndarray, bins: int | None = None) -> Grid | None:
    """Points counted on a regular grid (a 2D histogram); None without finite points."""
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return None
    n = bins or int(np.clip(np.sqrt(ok.sum() / 4), 10, 80))
    x_edges = _edges(x[ok], n)
    y_edges = _edges(y[ok], n)
    counts, _, _ = np.histogram2d(x[ok], y[ok], bins=[x_edges, y_edges])
    cell = np.full((len(x), 2), -1, dtype=np.int64)
    cell[ok, 0] = np.clip(np.searchsorted(x_edges, x[ok], side="right") - 1, 0, n - 1)
    cell[ok, 1] = np.clip(np.searchsorted(y_edges, y[ok], side="right") - 1, 0, n - 1)
    return Grid(counts, x_edges, y_edges, cell)


def _edges(values: np.ndarray, n: int) -> np.ndarray:
    lo, hi = float(values.min()), float(values.max())
    if hi - lo < 1e-12:
        lo, hi = lo - 0.5, hi + 0.5
    return np.linspace(lo, hi, n + 1)


# -- radar --------------------------------------------------------------------------------------
def unit_scale(columns: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Each column scaled to 0..1 by its own range: ``(scaled (rows, columns), lows, highs)``.
    A constant column sits at 0.5; blanks stay NaN."""
    data = np.column_stack(columns) if columns else np.empty((0, 0))
    lows = np.full(data.shape[1], np.nan)
    highs = np.full(data.shape[1], np.nan)
    scaled = np.full(data.shape, np.nan)
    for k in range(data.shape[1]):
        values = data[:, k]
        finite = values[np.isfinite(values)]
        if not len(finite):
            continue
        lows[k], highs[k] = float(finite.min()), float(finite.max())
        span = highs[k] - lows[k]
        scaled[:, k] = (values - lows[k]) / span if span > 1e-12 else 0.5
    return scaled, lows, highs


# -- activity landscape -------------------------------------------------------------------------
@dataclass(frozen=True)
class Surface:
    """A smooth height over a regular grid: ``heights[i, j]`` at ``(xs[i], ys[j])``; NaN where
    no point is near enough to say."""

    xs: np.ndarray
    ys: np.ndarray
    heights: np.ndarray


LANDSCAPE_REACH = 0.02  # cells whose kernel weight is below this share of the peak are empty


def landscape(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, bins: int | None = None, smooth: float = 1.2
) -> Surface | None:
    """The activity landscape of *z* over the plane of *x* and *y* (e.g. a chemical-space map):
    each grid node is the Gaussian-weighted mean of the points round it (Nadaraya-Watson
    kernel regression, width *smooth* grid cells), computed on binned sums so 200k points
    cost no more than a few hundred. Rough regions are activity cliffs; smooth plateaus are
    regions where structure changes little alter activity (Peltason & Bajorath 2007).

    Returns None with fewer than three complete points.
    """
    from scipy.ndimage import gaussian_filter

    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 3:
        return None
    x, y, z = x[ok], y[ok], z[ok]
    n = bins or int(np.clip(np.sqrt(len(x)) * 1.5, 16, 64))
    x_edges, y_edges = _edges(x, n), _edges(y, n)
    counts, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    sums, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=z)
    weight = gaussian_filter(counts, smooth, mode="constant")
    total = gaussian_filter(sums, smooth, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        heights = total / weight
    reach = LANDSCAPE_REACH * gaussian_filter(np.pad([[1.0]], 3), smooth).max()
    heights[weight < reach] = np.nan
    return Surface((x_edges[:-1] + x_edges[1:]) / 2, (y_edges[:-1] + y_edges[1:]) / 2, heights)
