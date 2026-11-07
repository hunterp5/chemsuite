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

"""Chemical-space maps: PCA, UMAP, t-SNE and PaCMAP of fingerprints or number columns, in two
or three dimensions, and the reduction to a few dozen dimensions that t-SNE, PaCMAP and
self-organising maps start from.

Fingerprints are streamed in tiles, so 200k compounds never need a dense float matrix.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from .clustering import standardize
from .fingerprints import Fingerprints

Progress = Callable[[int, int], None]

METHODS = {"pca": "PCA", "umap": "UMAP", "tsne": "t-SNE", "pacmap": "PaCMAP"}
#: PaCMAP distances. Jaccard and Hamming need fingerprints without PCA (they are measured on
#: the bits themselves); the others work on z-scored columns, bits or PCA dimensions.
PACMAP_METRICS = {
    "jaccard": "Jaccard (Tanimoto)",
    "hamming": "Hamming",
    "euclidean": "Euclidean",
    "manhattan": "Manhattan",
    "angular": "Angular",
}
DIMENSIONS = (2, 3)
#: UMAP distances: metric name -> title. The bit metrics need fingerprints.
UMAP_METRICS = {
    "jaccard": "Jaccard (Tanimoto)",
    "dice": "Dice",
    "rogerstanimoto": "Rogers-Tanimoto",
    "hamming": "Hamming",
    "euclidean": "Euclidean",
    "manhattan": "Manhattan",
    "cosine": "Cosine",
    "correlation": "Correlation",
}
BIT_METRICS = frozenset({"jaccard", "dice", "rogerstanimoto", "hamming"})
REDUCED_DIMS = 50  # PCA dimensions handed to the SOM (and to t-SNE / PaCMAP when asked)
TILE_ROWS = 4096
UMAP_FIT_MAX = 50_000  # larger sets: fit on a sample, then place the rest onto that map


@dataclass
class Embedding:
    coords: np.ndarray  # n × 2 or n × 3 (NaN for rows without features)
    method: str
    explained: tuple[float, ...] = ()  # PCA: share of variance on each axis
    notes: list[str] = field(default_factory=list)

    @property
    def placed(self) -> int:
        return int(np.isfinite(self.coords[:, 0]).sum())


@dataclass
class Features:
    """What a map is built from: fingerprints, or z-scored number columns (complete rows)."""

    n: int  # input rows
    rows: np.ndarray  # input rows that have features
    fps: Fingerprints | None = None  # of the rows in *rows*
    matrix: np.ndarray | None = None  # len(rows) × d

    @classmethod
    def from_fingerprints(cls, fps: Fingerprints) -> Features:
        rows = np.flatnonzero(fps.valid)
        return cls(len(fps), rows, fps=fps.take(rows))

    @classmethod
    def from_columns(cls, columns: Sequence[np.ndarray]) -> Features:
        x, rows = standardize(columns)
        n = len(columns[0]) if columns else 0
        return cls(n, rows, matrix=np.ascontiguousarray(x, dtype=np.float64))

    @classmethod
    def combined(
        cls,
        fps: Fingerprints,
        columns: Sequence[np.ndarray],
        *,
        dims: int = REDUCED_DIMS,
        progress: Progress | None = None,
        cancel: CancelToken | None = None,
    ) -> Features:
        """Fingerprints and number columns together, for rows that have both.

        The fingerprints are reduced to at most *dims* principal components (streamed, so no
        dense bit matrix), and the columns are z-scored. Each block is then scaled to the same
        total variance, so a 2048-bit fingerprint does not drown out three property columns,
        nor the reverse. Distances on the result are Euclidean.
        """
        x, column_rows = standardize(columns)
        rows = np.intersect1d(np.flatnonzero(fps.valid), column_rows)
        n = len(fps)
        if not len(rows) or not x.shape[1]:
            return cls(n, rows, matrix=np.empty((len(rows), 0)))
        bits = cls(n, rows, fps=fps.take(rows))
        scores = reduce(bits, dims, progress=progress, cancel=cancel)
        numbers = x[np.searchsorted(column_rows, rows)]

        def unit(block: np.ndarray) -> np.ndarray:
            total = float(block.var(axis=0).sum())
            return block / np.sqrt(total) if total > 1e-12 else block

        matrix = np.column_stack([unit(np.asarray(scores, dtype=np.float64)), unit(numbers)])
        return cls(n, rows, matrix=np.ascontiguousarray(matrix))

    def __len__(self) -> int:
        return len(self.rows)

    @property
    def dims(self) -> int:
        return self.fps.kind.bits if self.fps is not None else self.matrix.shape[1]

    def tiles(self) -> list[tuple[int, int]]:
        n = len(self)
        bounds = list(range(0, n, TILE_ROWS)) + [n]
        out = list(zip(bounds[:-1], bounds[1:], strict=True))
        if len(out) > 1 and out[-1][1] - out[-1][0] < REDUCED_DIMS:  # a tiny last tile
            out[-2] = (out[-2][0], out[-1][1])
            out.pop()
        return out

    def dense(self, start: int = 0, stop: int | None = None, dtype=np.float32) -> np.ndarray:
        if self.fps is not None:
            return self.fps.unpacked(start, stop).astype(dtype, copy=False)
        return self.matrix[start:stop].astype(dtype, copy=False)

    def bits(self, index: np.ndarray) -> np.ndarray:
        """Boolean fingerprint rows (for Jaccard distances)."""
        part = self.fps.take(index)
        return np.unpackbits(part.packed, axis=1, count=self.fps.kind.bits).astype(bool)


def _scatter(features: Features, values: np.ndarray, method: str, dims: int = 2, **kw) -> Embedding:
    coords = np.full((features.n, dims), np.nan)
    coords[features.rows] = values[:, :dims]
    return Embedding(coords, method, **kw)


def _dims(dims: int) -> int:
    if dims not in DIMENSIONS:
        raise ValueError(f"a map has 2 or 3 dimensions, not {dims}")
    return dims


def _need(features: Features, minimum: int) -> None:
    if len(features) < minimum:
        raise ValueError(f"A map needs at least {minimum} compounds with features.")
    if features.dims < 2:
        raise ValueError("A map needs at least two features (columns or fingerprint bits).")


def pca(
    features: Features,
    dims: int = 2,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(scores, explained variance ratios)`` for the rows that have features."""
    dims = max(1, min(dims, len(features) - 1, features.dims))
    if features.fps is None:
        from sklearn.decomposition import PCA

        model = PCA(n_components=dims, svd_solver="full" if features.dims <= 200 else "auto")
        scores = model.fit_transform(features.matrix)
        return scores, model.explained_variance_ratio_
    # Exact PCA from the covariance, summed tile by tile: d × d (2048² for most fingerprints)
    # is small, and the sums are matrix products, so this streams 200k rows quickly.
    tiles = features.tiles()
    d = features.dims
    gram = np.zeros((d, d), dtype=np.float64)
    total = np.zeros(d, dtype=np.float64)
    steps = 2 * len(tiles) + 1
    for i, (start, stop) in enumerate(tiles):
        check(cancel)
        block = features.dense(start, stop)
        gram += block.T @ block
        total += block.sum(axis=0)
        if progress is not None:
            progress(i + 1, steps)
    n = len(features)
    mean = total / n
    covariance = gram / n - np.outer(mean, mean)
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1][:dims]
    values = np.clip(values[order], 0.0, None)
    vectors = vectors[:, order]
    # Fix each axis's sign (largest loading positive), so a rerun draws the same map.
    signs = np.sign(vectors[np.abs(vectors).argmax(axis=0), np.arange(dims)])
    vectors *= np.where(signs == 0, 1.0, signs)
    variance = max(float(np.clip(np.trace(covariance), 1e-12, None)), 1e-12)
    if progress is not None:
        progress(len(tiles) + 1, steps)
    scores = np.empty((n, dims), dtype=np.float64)
    shift = mean @ vectors
    for i, (start, stop) in enumerate(tiles):
        check(cancel)
        scores[start:stop] = features.dense(start, stop, np.float64) @ vectors - shift
        if progress is not None:
            progress(len(tiles) + 2 + i, steps)
    return scores, values / variance


def bit_neighbours(
    fps: Fingerprints,
    k: int,
    metric: str = "jaccard",
    *,
    scores: bool = False,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """The *k* nearest rows of every row (self excluded), nearest first, as an ``n × k`` index
    array: by Tanimoto similarity (*metric* "jaccard") or by the number of differing bits
    ("hamming"). Exact, tile by tile, so memory stays small. With *scores*, also returns the
    matching similarities (Jaccard) or minus the differing bits (Hamming)."""
    if metric not in ("jaccard", "hamming"):
        raise ValueError(f"unknown bit distance {metric!r}")
    n = len(fps)
    k = max(1, min(k, n - 1))
    counts = fps.counts.astype(np.float32)
    bounds = list(range(0, n, TILE_ROWS)) + [n]
    tiles = list(zip(bounds[:-1], bounds[1:], strict=True))
    best = np.empty((n, k), dtype=np.int64)
    best_sim = np.empty((n, k), dtype=np.float32)
    for i, (a0, a1) in enumerate(tiles):
        block_a = fps.unpacked(a0, a1)
        top_sim = np.empty((a1 - a0, 0), dtype=np.float32)
        top_idx = np.empty((a1 - a0, 0), dtype=np.int64)
        for b0, b1 in tiles:
            check(cancel)
            common = block_a @ fps.unpacked(b0, b1).T
            union = counts[a0:a1, None] + counts[None, b0:b1] - common
            if metric == "jaccard":
                sim = np.where(union > 0, common / np.maximum(union, 1), 1.0)
            else:  # fewer differing bits is nearer
                sim = common - union
            sim = sim.astype(np.float32)
            if a0 == b0:
                np.fill_diagonal(sim, -np.inf)  # never its own neighbour
            sims = np.concatenate([top_sim, sim], axis=1)
            idx = np.concatenate([top_idx, np.broadcast_to(np.arange(b0, b1), sim.shape)], axis=1)
            keep = np.argpartition(-sims, k - 1, axis=1)[:, :k] if sims.shape[1] > k else None
            if keep is not None:
                sims = np.take_along_axis(sims, keep, axis=1)
                idx = np.take_along_axis(idx, keep, axis=1)
            top_sim, top_idx = sims, idx
        order = np.argsort(-top_sim, axis=1, kind="stable")
        best[a0:a1] = np.take_along_axis(top_idx, order, axis=1)
        best_sim[a0:a1] = np.take_along_axis(top_sim, order, axis=1)
        if progress is not None:
            progress(i + 1, len(tiles))
    return (best, best_sim) if scores else best


def _features_for(
    features: Features,
    reduce_to: int | None,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """What t-SNE and PaCMAP run on: *reduce_to* PCA dimensions, or (None or 0) the features
    themselves (fingerprints as 0/1 bits)."""
    if reduce_to:
        x = reduce(features, max(2, reduce_to), progress=progress, cancel=cancel)
    else:
        x = features.dense(0, None) if features.fps is not None else features.matrix
    return np.ascontiguousarray(x, dtype=np.float32)


def reduce(
    features: Features,
    dims: int = REDUCED_DIMS,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """At most *dims* columns per row: PCA scores, or the z-scores themselves if fewer."""
    if features.fps is None and features.dims <= dims:
        return features.matrix
    return pca(features, dims, progress=progress, cancel=cancel)[0]


def pca_map(
    features: Features, *, components: tuple[int, ...] = (1, 2), dims: int = 2, **kw
) -> Embedding:
    """PCA scores on two components (1-based, e.g. ``(1, 3)`` for PC1 against PC3); a 3D map
    (*dims* 3) is always PC1, PC2 and PC3."""
    _need(features, 3)
    dims = _dims(dims)
    if dims == 3:
        components = (1, 2, 3)
    if (
        len(components) != dims
        or any(not 1 <= a < b for a, b in zip(components, components[1:], strict=False))
        or components[0] < 1
    ):
        raise ValueError(f"components must be {dims} increasing numbers from 1, not {components}")
    last = components[-1]
    scores, ratio = pca(features, last, **kw)
    if scores.shape[1] < last:  # fewer components than asked: the missing axes are flat
        pad = last - scores.shape[1]
        scores = np.column_stack([scores, np.zeros((len(scores), pad))])
        ratio = np.concatenate([np.asarray(ratio, dtype=float), np.zeros(pad)])
    picked = [c - 1 for c in components]
    explained = tuple(float(ratio[i]) for i in picked)
    parts = [f"PC{c} {e:.0%}" for c, e in zip(components, explained, strict=True)]
    note = ", ".join(parts[:-1]) + f" and {parts[-1]} of the variance"
    return _scatter(features, scores[:, picked], "pca", dims, explained=explained, notes=[note])


def umap_map(
    features: Features,
    *,
    neighbours: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    metric: str | None = None,
    epochs: int | None = None,
    seed: int | None = 42,
    fit_max: int = UMAP_FIT_MAX,
    dims: int = 2,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """UMAP in *dims* (2 or 3) dimensions. The distance (*metric*) defaults to Jaccard on
    fingerprint bits and Euclidean on z-scored columns; a fingerprint metric needs
    fingerprints (see ``UMAP_METRICS``).

    Above *fit_max* rows the map is fitted on a random sample and the rest are placed onto it.
    A *seed* makes the map repeatable, at the cost of running on one core. *epochs* None lets
    UMAP choose (200 for large sets, 500 for small).
    """
    _need(features, 5)
    import umap  # slow to import (numba): only when a map is made

    n = len(features)
    neighbours = max(2, min(neighbours, n - 1))
    metric = metric or ("jaccard" if features.fps is not None else "euclidean")
    if metric in BIT_METRICS and features.fps is None:
        raise ValueError(f"{UMAP_METRICS[metric]} distance needs fingerprints")
    if min_dist > spread:
        raise ValueError("Min. distance must not be larger than the spread")
    rng = np.random.default_rng(seed if seed is not None else 0)
    fit_rows = np.sort(rng.choice(n, fit_max, replace=False)) if n > fit_max else np.arange(n)

    def take(index: np.ndarray) -> np.ndarray:
        return features.bits(index) if features.fps is not None else features.matrix[index]

    dims = _dims(dims)
    model = umap.UMAP(
        n_components=dims,
        n_neighbors=neighbours,
        min_dist=min_dist,
        spread=spread,
        metric=metric,
        n_epochs=epochs,
        random_state=seed,
        n_jobs=1 if seed is not None else -1,
    )
    check(cancel)
    if progress is not None:
        progress(0, 100)
    out = np.empty((n, dims))
    with warnings.catch_warnings():  # "inverse_transform will be unavailable" (not used)
        warnings.filterwarnings("ignore", message=".*inverse_transform.*")
        # Identical fingerprints are "disconnected" under Jaccard; UMAP still places them.
        warnings.filterwarnings("ignore", message=".*isconnect.*")
        out[fit_rows] = model.fit_transform(take(fit_rows))
    notes = []
    rest = np.setdiff1d(np.arange(n), fit_rows)
    if len(rest):
        notes.append(f"fitted on {len(fit_rows):,} sampled compounds; the rest placed onto it")
        chunks = np.array_split(rest, max(1, len(rest) // 10_000))
        for i, chunk in enumerate(chunks):
            check(cancel)
            out[chunk] = model.transform(take(chunk))
            if progress is not None:
                progress(50 + 50 * (i + 1) // len(chunks), 100)
    if progress is not None:
        progress(100, 100)
    return _scatter(features, out, "umap", dims, notes=notes)


def tsne_map(
    features: Features,
    *,
    perplexity: float = 30.0,
    learning_rate: float | None = None,
    iterations: int = 1000,
    exaggeration: float = 12.0,
    reduce_to: int | None = None,
    init: str = "pca",
    seed: int | None = 42,
    dims: int = 2,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """t-SNE (Barnes-Hut) in *dims* (2 or 3) dimensions.

    With *reduce_to* it runs on that many PCA dimensions (Euclidean); without, on the features
    themselves: Jaccard (Tanimoto) on fingerprint bits, Euclidean on z-scored columns.
    *learning_rate* None is scikit-learn's "auto" (rows / exaggeration / 4); *init* is
    "pca" (stable layouts) or "random".
    """
    _need(features, 5)
    from sklearn.manifold import TSNE

    dims = _dims(dims)
    if init not in ("pca", "random"):
        raise ValueError(f"unknown t-SNE start {init!r}")
    x = _features_for(features, reduce_to, progress=progress, cancel=cancel)
    jaccard = features.fps is not None and not reduce_to
    check(cancel)
    n = len(features)
    perplexity = min(perplexity, (n - 1) / 3)
    metric = "euclidean"
    if jaccard:
        # The neighbours t-SNE needs (3 × perplexity), found exactly on the bits, as a sparse
        # distance graph: much faster than scikit-learn's brute-force Jaccard.
        from scipy.sparse import csr_matrix

        k = min(n - 1, int(3.0 * perplexity + 1))
        near, sims = bit_neighbours(features.fps, k, scores=True, cancel=cancel)
        distances = np.clip(1.0 - sims.astype(np.float64), 1e-9, None)  # 0 would be dropped
        # Each row also lists itself (scikit-learn asks for it, then leaves it out).
        index = np.column_stack([np.arange(n), near])
        distances = np.column_stack([np.full(n, 1e-12), distances])
        graph = csr_matrix(
            (distances.ravel(), index.ravel(), np.arange(0, n * (k + 1) + 1, k + 1)),
            shape=(n, n),
        )
        if init == "pca":  # scikit-learn's own PCA start, which it skips for a graph
            start = np.asarray(reduce(features, dims, cancel=cancel), dtype=np.float32)
            if start.shape[1] < dims:
                start = np.column_stack([start, np.zeros((n, dims - start.shape[1]))])
            init = start / max(float(np.std(start[:, 0])), 1e-12) * 1e-4
        x, metric = graph, "precomputed"
    model = TSNE(
        n_components=dims,
        perplexity=perplexity,
        learning_rate="auto" if learning_rate is None else learning_rate,
        max_iter=max(250, iterations),
        early_exaggeration=exaggeration,
        init=init,
        metric=metric,
        random_state=seed,
        n_jobs=-1,
    )
    with warnings.catch_warnings():  # sparse graph efficiency notes
        warnings.simplefilter("ignore")
        coords = model.fit_transform(x)
    return _scatter(features, coords, "tsne", dims)


def pacmap_map(
    features: Features,
    *,
    neighbours: int = 10,
    mn_ratio: float = 0.5,
    fp_ratio: float = 2.0,
    distance: str | None = None,
    reduce_to: int | None = None,
    seed: int | None = 42,
    dims: int = 2,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """PaCMAP (Wang et al. 2021, JMLR 22:201): neighbour, mid-near and further pairs, so the
    map keeps both local groups and the overall layout.

    With *reduce_to* it runs on that many PCA dimensions (fingerprints are reduced by the
    streamed PCA, so no dense bit matrix is built); without, on the features themselves.
    *mn_ratio* and *fp_ratio* are the mid-near and further pairs per neighbour pair; *distance*
    is a key of ``PACMAP_METRICS`` (None: Jaccard on fingerprint bits, else Euclidean).
    Jaccard and Hamming find the neighbour pairs exactly on the bits (see :func:`bit_neighbours`);
    the mid-near and further pairs are sampled by Hamming distance.
    """
    _need(features, 5)
    dims = _dims(dims)
    bits = features.fps is not None and not reduce_to
    if distance is None:
        distance = "jaccard" if bits else "euclidean"
    if distance not in PACMAP_METRICS:
        raise ValueError(f"unknown PaCMAP distance {distance!r}")
    if distance in BIT_METRICS and not bits:
        raise ValueError(
            f"{PACMAP_METRICS[distance]} distance needs fingerprints without PCA dimensions"
        )
    import pacmap  # imports numba and faiss: only when a map is made

    x = _features_for(features, reduce_to, progress=progress, cancel=cancel)
    check(cancel)
    n = len(features)
    n_neighbors = max(2, min(neighbours, n - 2))
    pairs = None
    if distance in BIT_METRICS:
        near = bit_neighbours(features.fps, n_neighbors, distance, cancel=cancel)
        pairs = np.column_stack([np.repeat(np.arange(n), n_neighbors), near.ravel()])
        pairs = np.ascontiguousarray(pairs, dtype=np.int32)
        check(cancel)
    model = pacmap.PaCMAP(
        n_components=dims,
        n_neighbors=n_neighbors,
        MN_ratio=mn_ratio,
        FP_ratio=fp_ratio,
        pair_neighbors=pairs,
        distance="hamming" if distance in BIT_METRICS else distance,
        apply_pca=False,  # reduced above when asked
        random_state=seed,
    )
    quiet = logging.getLogger("pacmap")  # "random state is set", "high-dimensional data"
    level = quiet.level
    quiet.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coords = model.fit_transform(x, init="pca")
    finally:
        quiet.setLevel(level)
    check(cancel)
    return _scatter(features, np.asarray(coords, dtype=np.float64), "pacmap", dims)


OPTIONS = {
    "pca": ("components", "dims", "progress", "cancel"),
    "umap": (
        "neighbours",
        "min_dist",
        "spread",
        "metric",
        "epochs",
        "seed",
        "fit_max",
        "dims",
        "progress",
        "cancel",
    ),
    "tsne": (
        "perplexity",
        "learning_rate",
        "iterations",
        "exaggeration",
        "reduce_to",
        "init",
        "seed",
        "dims",
        "progress",
        "cancel",
    ),
    "pacmap": (
        "neighbours",
        "mn_ratio",
        "fp_ratio",
        "distance",
        "reduce_to",
        "seed",
        "dims",
        "progress",
        "cancel",
    ),
}


def embed(features: Features, method: str, **options) -> Embedding:
    """A 2D or 3D (``dims=3``) map by *method*; options another method takes (e.g. perplexity
    for UMAP) are ignored, so one call can carry every setting."""
    if method not in OPTIONS:
        raise ValueError(f"unknown method {method!r}")
    used = {k: v for k, v in options.items() if k in OPTIONS[method]}
    makers = {"pca": pca_map, "umap": umap_map, "tsne": tsne_map, "pacmap": pacmap_map}
    return makers[method](features, **used)
