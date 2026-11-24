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
"""Distances between rows, accurately and without running out of memory.

Three things in this package measure how far apart two rows are: the nearest-neighbor model,
the local difficulty an interval is scaled by, and the applicability domain. They all need the
same two properties, and getting either wrong is a real bug rather than a rounding detail.

**Accuracy.** The quick way to a distance matrix is ``‖a‖² − 2a·b + ‖b‖²``, which BLAS does in
one pass. It also cancels catastrophically when two rows are nearly equal: for a row compared
against *itself* it returns something around 1e-8 instead of zero. That is not a rounding
nuisance — scikit-learn's inverse-distance weighting gives a coincident training row all of
the weight, and a distance of 1e-8 instead of 0 silently takes a different branch. So the
difference is taken directly, ``a − b`` and then its norm, which is exact for equal rows.

**Memory.** Taking the difference directly needs rows × columns × features numbers at once, so
the rows are done in tiles sized from the actual width — a few million numbers at a time,
whatever the shape of the problem.

Numpy only: :mod:`chemsuite.qsar.apply` imports this, and nothing on the apply path may reach
for scikit-learn.
"""

from __future__ import annotations

import numpy as np

#: Numbers held at once while a tile of distances is computed, about 64 MB as float64.
TILE_NUMBERS = 8_000_000


def tile_rows(against: int, features: int) -> int:
    """How many rows to measure at a time, given what each one costs."""
    per_row = max(against * max(features, 1), 1)
    return max(1, min(4096, TILE_NUMBERS // per_row))


def distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Every distance between the rows of *a* and the rows of *b*, as (len(a), len(b)).

    Exact for equal rows, which is the property the expanded form gives away.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if not len(a) or not len(b):
        return np.zeros((len(a), len(b)))
    out = np.empty((len(a), len(b)))
    step = tile_rows(len(b), a.shape[1])
    for start in range(0, len(a), step):
        stop = min(start + step, len(a))
        out[start:stop] = np.linalg.norm(a[start:stop, None, :] - b[None, :, :], axis=2)
    return out


def nearest(a: np.ndarray, b: np.ndarray, k: int, skip: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """For each row of *a*, the *k* nearest rows of *b*: their indices and distances, nearest
    first, ties by row order.

    *skip* drops that many of the very nearest, which is how a row measured against the set it
    belongs to stops being its own neighbor.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    take = min(k + skip, len(b))
    if not len(a) or take <= skip:
        return np.zeros((len(a), 0), dtype=np.int64), np.zeros((len(a), 0))
    indices = np.empty((len(a), take - skip), dtype=np.int64)
    gaps = np.empty((len(a), take - skip))
    step = tile_rows(len(b), a.shape[1])
    for start in range(0, len(a), step):
        stop = min(start + step, len(a))
        block = distances(a[start:stop], b)
        close = np.argpartition(block, take - 1, axis=1)[:, :take] if take < len(b) else None
        if close is None:
            close = np.tile(np.arange(len(b)), (stop - start, 1))
        rows = np.arange(stop - start)[:, None]
        # lexsort: by distance, and by index where two distances are equal
        order = np.lexsort((close, block[rows, close]), axis=1)
        close = close[rows, order][:, skip:]
        indices[start:stop] = close
        gaps[start:stop] = block[rows, close]
    return indices, gaps


def mean_distance(a: np.ndarray, b: np.ndarray, k: int, skip: int = 0) -> np.ndarray:
    """Each row of *a*: the mean distance to its *k* nearest rows of *b*."""
    if not len(b) or not len(a) or k + skip > len(b):
        return np.zeros(len(a))
    return nearest(a, b, k, skip)[1].mean(axis=1)


def mean_of_neighbors(
    a: np.ndarray, b: np.ndarray, values: np.ndarray, k: int, skip: int = 0
) -> np.ndarray:
    """Each row of *a*: the mean *value* of its *k* nearest rows of *b*."""
    if not len(b) or not len(a) or k + skip > len(b):
        return np.zeros(len(a))
    indices = nearest(a, b, k, skip)[0]
    return np.asarray(values, dtype=np.float64)[indices].mean(axis=1)
