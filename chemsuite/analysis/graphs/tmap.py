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

"""TMAP: a large set drawn as the minimum spanning tree of its nearest-neighbour graph.

Probst, D. and Reymond, J.-L. (2020). "Visualization of very large high-dimensional data sets
as minimum spanning trees." *J. Cheminf.* 12, 12. The steps are the paper's:

1. the *k* nearest neighbours of every compound (:mod:`.knn`: exact up to 10,000 compounds,
   MinHash and a prefix forest above);
2. the graph joining each compound to them, weighted by 1 − Tanimoto;
3. its minimum spanning tree (a forest when the graph falls apart);
4. a layout of the tree (:mod:`.layout`, multilevel force), its components packed.

Each compound is a point; branches hold its nearest relatives. This is chemsuite's own code
(no part of the ``tmap`` package or OGDF): the neighbours rank by the same Tanimoto every other
chemsuite tool uses.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from chemsuite.analysis.fingerprints import Fingerprints
from chemsuite.core.jobs import CancelToken

from .knn import SEED, knn, knn_edges
from .layout import layout
from .trees import components, spanning_forest

Progress = Callable[[int, int], None]
NEIGHBOURS = 10


@dataclass(frozen=True)
class Tree:
    """A TMAP of *n* rows: tree edges ``(a, b)`` with their similarity, and a position per row
    (NaN for rows without a valid structure)."""

    a: np.ndarray
    b: np.ndarray
    similarity: np.ndarray
    coords: np.ndarray  # n × 2
    component: np.ndarray  # per row; -1 without a structure

    @property
    def placed(self) -> np.ndarray:
        return np.isfinite(self.coords[:, 0])

    def degree(self) -> np.ndarray:
        n = len(self.coords)
        return np.bincount(self.a, minlength=n) + np.bincount(self.b, minlength=n)


def tmap(
    fps: Fingerprints,
    k: int = NEIGHBOURS,
    *,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Tree:
    """The TMAP of *fps*. Progress runs over the neighbour search, then the layout."""
    n = len(fps)
    ok = np.flatnonzero(fps.valid)
    coords = np.full((n, 2), np.nan)
    label = np.full(n, -1, dtype=np.int64)
    empty = np.empty(0, dtype=np.int64)
    if not len(ok):
        return Tree(empty, empty.copy(), np.empty(0), coords, label)

    def stage(offset: int):
        if progress is None:
            return None
        return lambda done, total: progress(offset * 100 + round(100 * done / max(total, 1)), 200)

    sub = fps.take(ok)
    index, sim = knn(sub, k, seed=seed, progress=stage(0), cancel=cancel)
    a, b, s = knn_edges(index, sim)
    ta, tb, tw = spanning_forest(len(ok), a, b, 1.0 - s)
    labels = components(len(ok), ta, tb)
    coords[ok] = layout(len(ok), ta, tb, labels, seed=seed, progress=stage(1), cancel=cancel)
    label[ok] = labels
    return Tree(ok[ta], ok[tb], 1.0 - tw, coords, label)
