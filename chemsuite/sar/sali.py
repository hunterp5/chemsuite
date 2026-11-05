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
"""The Structure-Activity Landscape Index: SALI = |Δactivity| / (1 − similarity).

Large values flag activity cliffs: very similar compounds with very different activity.
Only pairs at least *min_similarity* alike are considered (the interesting ones, and it keeps
the pair count manageable).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from chemsuite.analysis.fingerprints import Fingerprints, from_counts, neighbours
from chemsuite.core.jobs import CancelToken

#: Similarity at or above 1 − EPSILON counts as identical (finite SALI, still ranked top).
EPSILON = 1e-3


@dataclass(frozen=True)
class Cliffs:
    a: np.ndarray  # row index per pair (a < b)
    b: np.ndarray
    similarity: np.ndarray
    delta: np.ndarray  # activity(b) − activity(a)
    sali: np.ndarray

    def __len__(self) -> int:
        return len(self.a)

    def per_row(self, n: int) -> tuple[np.ndarray, np.ndarray]:
        """``(highest SALI, partner row)`` per row (NaN and -1 without a pair)."""
        rows = np.concatenate([self.a, self.b])
        others = np.concatenate([self.b, self.a])
        values = np.concatenate([self.sali, self.sali])
        best = np.full(n, -np.inf)
        np.maximum.at(best, rows, values)
        partner = np.full(n, -1, dtype=np.int64)
        hit = values == best[rows]
        partner[rows[hit]] = others[hit]  # ties: any of the partners
        return np.where(np.isfinite(best), best, np.nan), partner


def sali_pairs(
    fps: Fingerprints,
    activity: np.ndarray,
    *,
    min_similarity: float = 0.6,
    metric: str = "tanimoto",
    skip_identical: bool = False,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> Cliffs:
    """Every pair at least *min_similarity* alike where both activities are numbers,
    highest SALI first."""
    activity = np.asarray(activity, dtype=np.float64)
    usable = Fingerprints(fps.kind, fps.packed, fps.valid & np.isfinite(activity))
    indptr, indices = neighbours(usable, min_similarity, metric, progress=progress, cancel=cancel)
    rows = np.repeat(np.arange(len(fps)), np.diff(indptr))
    keep = rows < indices
    a, b = rows[keep], indices[keep]
    words = fps.words
    common = np.bitwise_count(words[a] & words[b]).sum(axis=1)
    sim = from_counts(common, fps.counts[a], fps.counts[b], metric)
    if skip_identical:  # salt forms, stereoisomers, duplicates: SALI would only reflect EPSILON
        different = sim < 1.0
        a, b, sim = a[different], b[different], sim[different]
    delta = activity[b] - activity[a]
    sali = np.abs(delta) / np.maximum(1.0 - sim, EPSILON)
    order = np.lexsort((-sim, -np.abs(delta), -sali))
    return Cliffs(a[order], b[order], sim[order], delta[order], sali[order])
