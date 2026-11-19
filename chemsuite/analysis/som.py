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

"""Self-organizing maps (batch training): compounds binned onto a grid of cells, with
similar compounds in the same or neighboring cells.

Cells are squares or hexagons (odd rows shifted half a cell right, rows √3/2 apart, so every
cell is one apart from its six neighbors), optionally wrapped at the edges into a torus so
the border cells have as many neighbors as the middle ones.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from .embedding import Features, reduce

Progress = Callable[[int, int], None]
MAX_SIDE = 40
CHUNK_ROWS = 8192
NEIGHBORHOODS = {"gaussian": "Gaussian", "bubble": "Bubble"}
INITS = {"pca": "PCA plane", "random": "Random compounds"}
END_RADIUS = 0.7
HEX_PITCH = math.sqrt(3.0) / 2.0  # the distance between hexagon rows


@dataclass
class Som:
    rows: int
    columns: int
    weights: np.ndarray  # (rows × columns) × d
    unit: np.ndarray  # per input row: cell index (row-major), -1 without features
    error: float  # mean distance from each compound to its cell
    hexagonal: bool = False
    wrap: bool = False

    def cell(self) -> tuple[np.ndarray, np.ndarray]:
        """1-based ``(grid row, grid column)`` per input row (NaN without features)."""
        ok = self.unit >= 0
        r = np.where(ok, self.unit // self.columns + 1, np.nan)
        c = np.where(ok, self.unit % self.columns + 1, np.nan)
        return r, c

    def counts(self) -> np.ndarray:
        units = self.unit[self.unit >= 0]
        return np.bincount(units, minlength=self.rows * self.columns).reshape(
            self.rows, self.columns
        )


def default_size(n: int) -> tuple[int, int]:
    """About 5·√n cells (Vesanto's rule of thumb), a little wider than tall."""
    cells = max(4, 5 * math.sqrt(max(n, 1)))
    rows = max(2, min(MAX_SIDE, round(math.sqrt(cells / 1.25))))
    columns = max(2, min(MAX_SIDE, round(cells / rows)))
    return rows, columns


def _best_units(x: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest cell and its squared distance per row, chunked."""
    w2 = (weights**2).sum(axis=1)
    best = np.empty(len(x), dtype=np.int64)
    dist = np.empty(len(x))
    for start in range(0, len(x), CHUNK_ROWS):
        part = x[start : start + CHUNK_ROWS]
        d = w2[None, :] - 2.0 * part @ weights.T
        best[start : start + len(part)] = d.argmin(axis=1)
        dist[start : start + len(part)] = d[np.arange(len(part)), best[start : start + len(part)]]
    dist += (x**2).sum(axis=1)
    return best, np.maximum(dist, 0.0)


def grid_positions(rows: int, columns: int, hexagonal: bool = False) -> np.ndarray:
    """``(rows × columns) × 2`` cell centers, row-major: squares one apart, or hexagons with
    odd rows shifted half a cell right and rows :data:`HEX_PITCH` apart."""
    r, c = np.divmod(np.arange(rows * columns), columns)
    if not hexagonal:
        return np.column_stack([r, c]).astype(np.float64)
    return np.column_stack([r * HEX_PITCH, c + 0.5 * (r % 2)])


def grid_distances2(
    rows: int, columns: int, *, hexagonal: bool = False, wrap: bool = False
) -> np.ndarray:
    """Squared distances between every pair of cells; with *wrap*, the shorter way round
    each edge (a torus)."""
    at = grid_positions(rows, columns, hexagonal)
    delta = np.abs(at[:, None, :] - at[None, :, :])
    if wrap:
        period = np.array([rows * (HEX_PITCH if hexagonal else 1.0), float(columns)])
        delta = np.minimum(delta, period - delta)
    return (delta**2).sum(axis=2)


def train(
    x: np.ndarray,
    rows: int,
    columns: int,
    *,
    epochs: int = 20,
    seed: int = 0,
    hexagonal: bool = False,
    wrap: bool = False,
    neighborhood: str = "gaussian",
    start_radius: float | None = None,
    end_radius: float = END_RADIUS,
    init: str = "pca",
    fit: np.ndarray | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Batch SOM on *x*: ``(weights, cell per row, squared distance per row)``.

    The cells start spread over the plane of the first two principal components (*init*
    "pca"), so the result does not depend on a random start, or on randomly drawn compounds
    ("random", by *seed*). The neighborhood (Gaussian, or a "bubble" that takes every cell
    within the radius alike) shrinks from *start_radius* (None: half the grid) to
    *end_radius*. With *fit* (row positions) only those rows train the map; every row is then
    placed on its nearest cell.
    """
    if neighborhood not in NEIGHBORHOODS:
        raise ValueError(f"unknown neighborhood {neighborhood!r}")
    if init not in INITS:
        raise ValueError(f"unknown SOM start {init!r}")
    every = np.asarray(x, dtype=np.float64)
    x = every if fit is None else every[fit]
    units = rows * columns
    rng = np.random.default_rng(seed)
    mean = x.mean(axis=0)
    if init == "random":
        weights = x[rng.choice(len(x), units, replace=len(x) < units)].copy()
    else:
        centred = x - mean
        if len(x) > 2 and x.shape[1] > 1:
            sample = centred[rng.permutation(len(x))[:20_000]]
            _u, s, vt = np.linalg.svd(sample, full_matrices=False)
            scale = s[:2] / math.sqrt(max(1, len(sample) - 1))
            two = len(scale) == 2
            axes = vt[:2] * scale[:, None] if two else np.vstack([vt[:1], vt[:1] * 0])
        else:
            axes = np.zeros((2, x.shape[1]))
        at = grid_positions(rows, columns, hexagonal)
        a = (at[:, 0] / max(1e-12, at[:, 0].max()) - 0.5) * 4.0 if rows > 1 else at[:, 0] * 0
        b = (at[:, 1] / max(1e-12, at[:, 1].max()) - 0.5) * 4.0 if columns > 1 else at[:, 1] * 0
        weights = mean + a[:, None] * axes[0] + b[:, None] * axes[1]
    grid_d2 = grid_distances2(rows, columns, hexagonal=hexagonal, wrap=wrap)
    # A bubble that covers the whole grid makes every cell the same mean, every compound then
    # ties to the first cell and the map never unfolds: each bubble leaves out its farthest cells.
    farthest = grid_d2.max(axis=1, keepdims=True) * (1.0 - 1e-9)
    start_sigma = max(rows, columns) / 2.0 if start_radius is None else float(start_radius)
    end_sigma = min(float(end_radius), start_sigma)
    if end_sigma <= 0:
        raise ValueError("the neighborhood radius must be above 0")
    for epoch in range(epochs):
        check(cancel)
        sigma = start_sigma * (end_sigma / start_sigma) ** (epoch / max(1, epochs - 1))
        best, _ = _best_units(x, weights)
        counts = np.bincount(best, minlength=units).astype(np.float64)
        sums = np.zeros((units, x.shape[1]))
        np.add.at(sums, best, x)
        if neighborhood == "bubble":
            h = (grid_d2 <= np.minimum(sigma * sigma, farthest)).astype(np.float64)
        else:
            h = np.exp(-grid_d2 / (2.0 * sigma * sigma))
        numerator = h @ sums
        denominator = h @ counts
        moved = denominator > 1e-12
        weights[moved] = numerator[moved] / denominator[moved, None]
        if progress is not None:
            progress(epoch + 1, epochs)
    best, d2 = _best_units(every, weights)
    return weights, best, d2


def som(
    features: Features,
    rows: int | None = None,
    columns: int | None = None,
    *,
    epochs: int = 20,
    fit: np.ndarray | None = None,
    hexagonal: bool = False,
    wrap: bool = False,
    neighborhood: str = "gaussian",
    start_radius: float | None = None,
    end_radius: float = END_RADIUS,
    init: str = "pca",
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Som:
    """A map of *features* on a *rows* × *columns* grid (None: :func:`default_size`); *fit*
    (input rows) trains it on those rows alone. See :func:`train` for the rest."""
    if len(features) < 2:
        raise ValueError("A self-organizing map needs at least two compounds with features.")
    fit = features.positions(fit)
    if fit is not None and len(fit) < 2:
        raise ValueError("A self-organizing map needs at least two compounds to train on.")
    if rows is None or columns is None:
        rows, columns = default_size(len(features) if fit is None else len(fit))
    x = reduce(features, cancel=cancel)
    weights, best, d2 = train(
        x, rows, columns, epochs=epochs, seed=seed, hexagonal=hexagonal, wrap=wrap,
        neighborhood=neighborhood, start_radius=start_radius, end_radius=end_radius,
        init=init, fit=fit, progress=progress, cancel=cancel,
    )  # fmt: skip
    unit = np.full(features.n, -1, dtype=np.int64)
    unit[features.rows] = best
    return Som(rows, columns, weights, unit, float(np.sqrt(d2).mean()), hexagonal, wrap)
