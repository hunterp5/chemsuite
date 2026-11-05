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
