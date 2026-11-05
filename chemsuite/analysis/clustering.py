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
"""Clustering compounds by fingerprint (Tanimoto) or by number columns (scaled Euclidean).

Clusters are numbered by size: cluster 1 is the largest. Each has one representative (its
"centre"): Butina's centroid, the member nearest a k-means centre, or a hierarchical
cluster's medoid. Rows that cannot be clustered (no structure, a missing number) get none.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from .fingerprints import Fingerprints, from_counts, neighbours

METHODS = {
    "butina": "Butina (similarity cutoff)",
    "kmeans": "K-means (number of clusters)",
    "hierarchical": "Hierarchical (number of clusters)",
}
LINKAGES = {"average": "Average", "complete": "Complete", "single": "Single", "ward": "Ward"}
MAX_HIERARCHICAL = 8000  # the full distance matrix grows with the square of the rows
BUTINA_WARN_ROWS = 50_000
Progress = Callable[[int, int], None]


@dataclass(frozen=True)
class Clusters:
    labels: np.ndarray  # 0 = largest cluster, -1 = not clustered
    centres: np.ndarray  # bool: the representative of its cluster

    @property
    def count(self) -> int:
        return int(self.labels.max()) + 1 if (self.labels >= 0).any() else 0

    def sizes(self) -> np.ndarray:
        """Cluster size per row (0 where not clustered)."""
        per_cluster = np.bincount(self.labels[self.labels >= 0], minlength=self.count)
        return np.where(self.labels >= 0, per_cluster[np.maximum(self.labels, 0)], 0)

    @property
    def singletons(self) -> int:
        return int(np.sum(np.bincount(self.labels[self.labels >= 0]) == 1))


def finish(labels: np.ndarray, centres: Sequence[int]) -> Clusters:
    """Renumber clusters largest first (ties: first row first) and flag the centres."""
    labels = np.asarray(labels, dtype=np.int64)
    out = np.full(len(labels), -1, dtype=np.int64)
    member = labels >= 0
    if member.any():
        ids, first, counts = np.unique(labels[member], return_index=True, return_counts=True)
        order = np.lexsort((first, -counts))
        rank = np.empty(len(ids), dtype=np.int64)
        rank[order] = np.arange(len(ids))
        out[member] = rank[np.searchsorted(ids, labels[member])]
    flags = np.zeros(len(labels), dtype=np.bool_)
    flags[np.asarray(list(centres), dtype=np.int64)] = True
    return Clusters(out, flags & member)


# -- fingerprints ------------------------------------------------------------------------------
def butina(
    fps: Fingerprints,
    cutoff: float = 0.6,
    metric: str = "tanimoto",
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Clusters:
    """Taylor-Butina: the compound with the most neighbours (similarity ≥ *cutoff*) founds a
    cluster with its unclustered neighbours; repeat. Remaining compounds are singletons."""
    indptr, indices = neighbours(fps, cutoff, metric, progress=progress, cancel=cancel)
    n = len(fps)
    counts = np.diff(indptr)
    labels = np.full(n, -1, dtype=np.int64)
    centres: list[int] = []
    order = np.lexsort((np.arange(n), -counts))  # most neighbours first, then row order
    for i in order[fps.valid[order]].tolist():
        if labels[i] >= 0:
            continue
        cluster = len(centres)
        labels[i] = cluster
        members = indices[indptr[i] : indptr[i + 1]]
        free = members[labels[members] < 0]
        labels[free] = cluster
        centres.append(i)
    return finish(labels, centres)


def _tiles(n: int, size: int = 4096):
    return [(s, min(n, s + size)) for s in range(0, n, size)]


def kmeans_fingerprints(
    fps: Fingerprints,
    k: int,
    *,
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Clusters:
    """Mini-batch k-means on the fingerprint bits, streamed in tiles to keep memory flat."""
    from sklearn.cluster import MiniBatchKMeans

    rows = np.flatnonzero(fps.valid)
    if len(rows) < k:
        raise ValueError(f"{k} clusters need at least {k} structures (there are {len(rows)}).")
    valid = fps.take(rows)
    model = MiniBatchKMeans(n_clusters=k, random_state=seed, batch_size=4096, n_init=3)
    tiles = _tiles(len(valid))
    epochs = 3
    steps = epochs * len(tiles) + len(tiles)
    step = 0
    rng = np.random.default_rng(seed)
    for _ in range(epochs):
        for start, stop in [tiles[i] for i in rng.permutation(len(tiles))]:
            check(cancel)
            block = valid.unpacked(start, stop)
            if len(block) >= k:
                model.partial_fit(block)
            step += 1
            if progress is not None:
                progress(step, steps)
    if not hasattr(model, "cluster_centers_"):  # every tile smaller than k: fit at once
        model.partial_fit(valid.unpacked())
    labels = np.empty(len(valid), dtype=np.int64)
    distance = np.empty(len(valid))
    for start, stop in tiles:
        check(cancel)
        block = valid.unpacked(start, stop)
        labels[start:stop] = model.predict(block)
        centre = model.cluster_centers_[labels[start:stop]]
        distance[start:stop] = ((block - centre) ** 2).sum(axis=1)
        step += 1
        if progress is not None:
            progress(step, steps)
    return _with_nearest_centres(rows, labels, distance, len(fps))


def _with_nearest_centres(rows: np.ndarray, labels: np.ndarray, distance: np.ndarray, n: int):
    order = np.lexsort((distance, labels))
    first = order[np.r_[True, labels[order][1:] != labels[order][:-1]]]
    full = np.full(n, -1, dtype=np.int64)
    full[rows] = labels
    return finish(full, rows[first].tolist())


def _check_hierarchical(n: int) -> None:
    if n > MAX_HIERARCHICAL:
        raise ValueError(
            f"Hierarchical clustering is limited to {MAX_HIERARCHICAL:,} rows (this has {n:,}): "
            "use Butina or k-means, or select fewer rows."
        )


def _cut(z: np.ndarray, k: int) -> np.ndarray:
    from scipy.cluster.hierarchy import fcluster

    return fcluster(z, t=k, criterion="maxclust").astype(np.int64) - 1


def hierarchical_fingerprints(
    fps: Fingerprints,
    k: int,
    linkage: str = "average",
    metric: str = "tanimoto",
    *,
    cancel: CancelToken | None = None,
) -> Clusters:
    from scipy.cluster.hierarchy import linkage as build
    from scipy.spatial.distance import squareform

    if linkage == "ward":
        raise ValueError("Ward linkage needs number columns; use average or complete here.")
    rows = np.flatnonzero(fps.valid)
    _check_hierarchical(len(rows))
    if len(rows) < max(2, k):
        raise ValueError(f"{k} clusters need at least {max(2, k)} structures.")
    valid = fps.take(rows)
    bits = valid.unpacked()
    sim = from_counts(bits @ bits.T, valid.counts[:, None], valid.counts[None, :], metric)
    check(cancel)
    distance = 1.0 - sim
    np.fill_diagonal(distance, 0.0)
    labels = _cut(build(squareform(distance, checks=False), method=linkage), k)
    centres = []
    for cluster in np.unique(labels):
        members = np.flatnonzero(labels == cluster)
        within = sim[np.ix_(members, members)].sum(axis=1)
        centres.append(int(rows[members[within.argmax()]]))  # the medoid
    full = np.full(len(fps), -1, dtype=np.int64)
    full[rows] = labels
    return finish(full, centres)


# -- number columns ----------------------------------------------------------------------------
def standardize(columns: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """``(z-scores of complete rows, row indices)``: rows missing any value are left out.

    A constant column carries no information and is dropped.
    """
    data = np.column_stack([np.asarray(c, dtype=np.float64) for c in columns])
    rows = np.flatnonzero(np.isfinite(data).all(axis=1))
    x = data[rows]
    if not len(rows):
        return x, rows
    sd = x.std(axis=0)
    keep = sd > 0
    return (x[:, keep] - x[:, keep].mean(axis=0)) / sd[keep], rows


def kmeans_numeric(columns: Sequence[np.ndarray], k: int, *, seed: int = 0) -> Clusters:
    from sklearn.cluster import KMeans

    x, rows = standardize(columns)
    if len(rows) < k:
        raise ValueError(f"{k} clusters need at least {k} complete rows (there are {len(rows)}).")
    if not x.shape[1]:
        raise ValueError("The chosen columns do not vary.")
    model = KMeans(n_clusters=k, random_state=seed, n_init="auto").fit(x)
    distance = ((x - model.cluster_centers_[model.labels_]) ** 2).sum(axis=1)
    return _with_nearest_centres(rows, model.labels_.astype(np.int64), distance, len(columns[0]))


def hierarchical_numeric(columns: Sequence[np.ndarray], k: int, linkage: str = "ward") -> Clusters:
    from scipy.cluster.hierarchy import linkage as build

    x, rows = standardize(columns)
    _check_hierarchical(len(rows))
    if len(rows) < max(2, k):
        raise ValueError(f"{k} clusters need at least {max(2, k)} complete rows.")
    if not x.shape[1]:
        raise ValueError("The chosen columns do not vary.")
    labels = _cut(build(x, method=linkage), k)
    centres = []
    for cluster in np.unique(labels):
        members = np.flatnonzero(labels == cluster)
        centre = x[members].mean(axis=0)
        centres.append(int(rows[members[((x[members] - centre) ** 2).sum(axis=1).argmin()]]))
    full = np.full(len(columns[0]), -1, dtype=np.int64)
    full[rows] = labels
    return finish(full, centres)
