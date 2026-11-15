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

"""Chemical Space Networks: compounds as nodes, joined when they are alike.

Maggiora, G. M. and Bajorath, J. (2014). "Chemical space networks: a powerful new paradigm
for the description of chemical space." *J. Comput.-Aided Mol. Des.* 28, 795–802; Zwierzyna,
M., Vogt, M., Maggiora, G. M. and Bajorath, J. (2015). "Design and characterization of chemical
space networks for different compound data sets." *J. Comput.-Aided Mol. Des.* 29, 113–125.

- **Edges** join every pair whose Tanimoto similarity reaches a threshold (threshold CSNs), or
  every matched molecular pair (MMP CSNs, :func:`chemsuite.sar.mmp.find_pairs`).
- **Layout.** Each connected part is laid out by force (:mod:`.layout`) and the parts packed,
  largest first; compounds with no edge (singletons) sit in rows at the end, or are left out.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import numpy as np

from chemsuite.analysis.fingerprints import Fingerprints, from_counts, neighbours
from chemsuite.core.jobs import CancelToken

from .knn import SEED
from .layout import layout
from .trees import components

Progress = Callable[[int, int], None]
THRESHOLD = 0.65
#: more edges than this make a hairball that draws slowly: the threshold should go up
MAX_EDGES = 50_000


@dataclass(frozen=True)
class Network:
    """A network of *n* rows: edges ``(a, b)`` (``a < b``) with their similarity (NaN for a
    matched pair), and a position per row (NaN for rows not drawn)."""

    a: np.ndarray
    b: np.ndarray
    similarity: np.ndarray
    coords: np.ndarray  # n × 2
    component: np.ndarray  # per row; -1 when not drawn

    @property
    def placed(self) -> np.ndarray:
        return np.isfinite(self.coords[:, 0])

    def degree(self) -> np.ndarray:
        n = len(self.coords)
        return np.bincount(self.a, minlength=n) + np.bincount(self.b, minlength=n)


def similarity_edges(
    fps: Fingerprints,
    threshold: float = THRESHOLD,
    *,
    max_edges: int = MAX_EDGES,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(a, b, similarity)``: every pair at least *threshold* alike, ``a < b``. Raises
    :class:`~chemsuite.analysis.fingerprints.TooManyPairs` past *max_edges*."""
    indptr, indices = neighbours(
        fps, threshold, max_pairs=max_edges, progress=progress, cancel=cancel
    )
    rows = np.repeat(np.arange(len(fps)), np.diff(indptr))
    keep = rows < indices
    a, b = rows[keep], indices[keep]
    common = np.bitwise_count(fps.words[a] & fps.words[b]).sum(axis=1)
    return a, b, from_counts(common, fps.counts[a], fps.counts[b], "tanimoto")


def mmp_edges(
    cells: Sequence[bytes | str | None],
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(a, b)``: every matched molecular pair, ``a < b``."""
    from chemsuite.sar.mmp import find_pairs

    pairs = find_pairs(
        cells,
        np.full(len(cells), np.nan),
        executor=executor,
        progress=progress,
        cancel=cancel,
    )
    a = np.array([p.a for p in pairs], dtype=np.int64)
    b = np.array([p.b for p in pairs], dtype=np.int64)
    return a, b


def network(
    n: int,
    a: np.ndarray,
    b: np.ndarray,
    similarity: np.ndarray | None = None,
    *,
    valid: np.ndarray | None = None,
    singletons: bool = True,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Network:
    """Lay out the network of *n* rows with edges ``(a, b)``. Rows not *valid* are not drawn;
    nor are rows without an edge unless *singletons*."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    sim = np.full(len(a), np.nan) if similarity is None else np.asarray(similarity, np.float64)
    shown = np.ones(n, bool) if valid is None else np.asarray(valid, bool).copy()
    if not singletons:
        linked = np.zeros(n, bool)
        linked[a] = linked[b] = True
        shown &= linked
    keep = shown[a] & shown[b]
    a, b, sim = a[keep], b[keep], sim[keep]
    rows = np.flatnonzero(shown)
    coords = np.full((n, 2), np.nan)
    label = np.full(n, -1, dtype=np.int64)
    if len(rows):
        local = np.full(n, -1, dtype=np.int64)
        local[rows] = np.arange(len(rows))
        la, lb = local[a], local[b]
        labels = components(len(rows), la, lb)
        coords[rows] = layout(
            len(rows), la, lb, labels, seed=seed, progress=progress, cancel=cancel
        )
        label[rows] = labels
    return Network(a, b, sim, coords, label)
