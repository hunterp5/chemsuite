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
"""Diverse subsets by MaxMin picking: each pick is the compound farthest from all earlier ones.

Distances are 1 − similarity for fingerprints, or Euclidean on z-scored number columns.
Picking can start from compounds already chosen (*seeds*), so a new batch complements them.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import Cancelled, CancelToken

from .clustering import standardize
from .fingerprints import Fingerprints, max_similarity, similarity

Progress = Callable[[int, int], None]


@dataclass(frozen=True)
class Picks:
    rows: np.ndarray  # picked row indices, in pick order
    distances: np.ndarray  # each pick's distance to the nearest earlier pick or seed

    def __len__(self) -> int:
        return len(self.rows)


def _maxmin(
    nearest: np.ndarray,
    count: int,
    distance_to: Callable[[int], np.ndarray],
    first: int | None,
    progress: Progress | None,
    cancel: CancelToken | None,
) -> Picks:
    """*nearest*: each row's distance to the seeds (inf without seeds, -inf if not pickable).

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the picks made so far
    (each as far from the earlier ones as a full run would make it), or None if there are none.
    """
    rows: list[int] = []
    distances: list[float] = []
    while len(rows) < count:
        if cancel is not None and cancel.cancelled:
            raise Cancelled(_picks(rows, distances) if rows else None)
        if first is not None and not rows:
            pick = first
        else:
            pick = int(nearest.argmax())
            if not nearest[pick] > 0:  # only duplicates of earlier picks remain
                break
        rows.append(pick)
        distances.append(float(nearest[pick]))
        nearest = np.minimum(nearest, distance_to(pick))
        nearest[pick] = -np.inf
        if progress is not None:
            progress(len(rows), count)
    return _picks(rows, distances)


def _picks(rows: list[int], distances: list[float]) -> Picks:
    return Picks(np.array(rows, dtype=np.int64), np.array(distances))


def maxmin_fingerprints(
    fps: Fingerprints,
    count: int,
    *,
    seeds: Sequence[int] = (),
    metric: str = "tanimoto",
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Picks:
    """Pick *count* rows (fewer if the rest duplicate the picks). Seeds are never picked."""
    nearest = np.where(fps.valid, np.inf, -np.inf)
    seeds = np.asarray(list(seeds), dtype=np.int64)
    first = None
    if len(seeds):
        try:
            best, _ = max_similarity(fps, fps.take(seeds), metric, cancel=cancel)
        except Cancelled:
            raise Cancelled from None  # no picks yet; its partial is not a Picks
        nearest = np.where(fps.valid, 1.0 - np.nan_to_num(best, nan=0.0), -np.inf)
        nearest[seeds] = -np.inf
    else:
        candidates = np.flatnonzero(fps.valid)
        if not len(candidates):
            return Picks(np.empty(0, np.int64), np.empty(0))
        first = int(np.random.default_rng(seed).choice(candidates))

    def distance_to(row: int) -> np.ndarray:
        return 1.0 - np.nan_to_num(similarity(fps, fps.take(np.array([row])), metric), nan=1.0)

    return _maxmin(nearest, count, distance_to, first, progress, cancel)


def maxmin_numeric(
    columns: Sequence[np.ndarray],
    count: int,
    *,
    seeds: Sequence[int] = (),
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Picks:
    """MaxMin on z-scored columns; rows missing a value are never picked."""
    x, rows = standardize(columns)
    n = len(columns[0])
    full = np.full((n, x.shape[1]), np.nan)
    full[rows] = x
    ok = np.zeros(n, dtype=np.bool_)
    ok[rows] = True
    nearest = np.where(ok, np.inf, -np.inf)
    seeds = [s for s in seeds if ok[s]]
    first = None
    for s in seeds:
        nearest = np.minimum(
            nearest, np.where(ok, np.sqrt(((full - full[s]) ** 2).sum(axis=1)), -np.inf)
        )
    if len(seeds):
        nearest[np.asarray(seeds)] = -np.inf
    elif len(rows):
        first = int(np.random.default_rng(seed).choice(rows))
    else:
        return Picks(np.empty(0, np.int64), np.empty(0))

    def distance_to(row: int) -> np.ndarray:
        return np.where(ok, np.sqrt(((full - full[row]) ** 2).sum(axis=1)), -np.inf)

    return _maxmin(nearest, count, distance_to, first, progress, cancel)
