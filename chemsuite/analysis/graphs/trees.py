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

"""Minimum spanning forests and connected components of compound graphs."""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree

#: added to every weight: scipy reads a weight of 0 as no edge, and identical compounds are
#: 0 apart
OFFSET = 1e-9


def spanning_forest(
    n: int, a: np.ndarray, b: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The minimum spanning forest of the undirected graph ``(a, b)`` with weights *w*:
    ``(a, b, w)`` with ``a < b``, one tree per connected component (n − components edges)."""
    if not len(a):
        empty = np.empty(0, dtype=np.int64)
        return empty, empty.copy(), np.empty(0)
    graph = coo_matrix((np.asarray(w, dtype=np.float64) + OFFSET, (a, b)), shape=(n, n))
    tree = minimum_spanning_tree(graph.tocsr()).tocoo()
    ta, tb = tree.row.astype(np.int64), tree.col.astype(np.int64)
    lo, hi = np.minimum(ta, tb), np.maximum(ta, tb)
    order = np.lexsort((hi, lo))
    return lo[order], hi[order], tree.data[order] - OFFSET


def components(n: int, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """A component label per node (0 the largest component, then by size; ties by first node)."""
    if n == 0:
        return np.empty(0, dtype=np.int64)
    graph = coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    first = np.full(len(sizes), n)
    np.minimum.at(first, labels, np.arange(n))
    rank = np.empty(len(sizes), dtype=np.int64)
    rank[np.lexsort((first, -sizes))] = np.arange(len(sizes))
    return rank[labels]
