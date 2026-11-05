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

"""Self-organising maps (batch training): compounds binned onto a grid of cells, with
similar compounds in the same or neighbouring cells."""

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


@dataclass
class Som:
    rows: int
    columns: int
    weights: np.ndarray  # (rows × columns) × d
    unit: np.ndarray  # per input row: cell index (row-major), -1 without features
    error: float  # mean distance from each compound to its cell

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


def train(
    x: np.ndarray,
    rows: int,
    columns: int,
    *,
    epochs: int = 20,
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Batch SOM on *x*: ``(weights, cell per row, squared distance per row)``.

    Cells start spread over the plane of the first two principal components, so the result
    does not depend on a random start; the neighbourhood shrinks from half the grid to one cell.
    """
    x = np.asarray(x, dtype=np.float64)
    units = rows * columns
    grid = np.array([(r, c) for r in range(rows) for c in range(columns)], dtype=np.float64)
    mean = x.mean(axis=0)
    centred = x - mean
    if len(x) > 2 and x.shape[1] > 1:
        sample = centred[np.random.default_rng(seed).permutation(len(x))[:20_000]]
        _u, s, vt = np.linalg.svd(sample, full_matrices=False)
        scale = s[:2] / math.sqrt(max(1, len(sample) - 1))
        axes = vt[:2] * scale[:, None] if len(scale) == 2 else np.vstack([vt[:1], vt[:1] * 0])
    else:
        axes = np.zeros((2, x.shape[1]))
    a = (grid[:, 0] / max(1, rows - 1) - 0.5) * 4.0
    b = (grid[:, 1] / max(1, columns - 1) - 0.5) * 4.0
    weights = mean + a[:, None] * axes[0] + b[:, None] * axes[1]
    grid_d2 = ((grid[:, None, :] - grid[None, :, :]) ** 2).sum(axis=2)
    start_sigma, end_sigma = max(rows, columns) / 2.0, 0.7
    for epoch in range(epochs):
        check(cancel)
        sigma = start_sigma * (end_sigma / start_sigma) ** (epoch / max(1, epochs - 1))
        best, _ = _best_units(x, weights)
        counts = np.bincount(best, minlength=units).astype(np.float64)
        sums = np.zeros((units, x.shape[1]))
        np.add.at(sums, best, x)
        h = np.exp(-grid_d2 / (2.0 * sigma * sigma))
        numerator = h @ sums
        denominator = h @ counts
        moved = denominator > 1e-12
        weights[moved] = numerator[moved] / denominator[moved, None]
        if progress is not None:
            progress(epoch + 1, epochs)
    best, d2 = _best_units(x, weights)
    return weights, best, d2


def som(
    features: Features,
    rows: int | None = None,
    columns: int | None = None,
    *,
    epochs: int = 20,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Som:
    if len(features) < 2:
        raise ValueError("A self-organising map needs at least two compounds with features.")
    if rows is None or columns is None:
        rows, columns = default_size(len(features))
    x = reduce(features, cancel=cancel)
    weights, best, d2 = train(x, rows, columns, epochs=epochs, progress=progress, cancel=cancel)
    unit = np.full(features.n, -1, dtype=np.int64)
    unit[features.rows] = best
    return Som(rows, columns, weights, unit, float(np.sqrt(d2).mean()))
