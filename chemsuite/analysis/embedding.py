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
#: PaCMAP distances (on z-scored columns, or on the PCA-reduced fingerprints)
PACMAP_METRICS = {"euclidean": "Euclidean", "manhattan": "Manhattan", "angular": "Angular"}
DIMENSIONS = (2, 3)
#: UMAP distances: metric name -> title. The bit metrics need fingerprints.
UMAP_METRICS = {
    "jaccard": "Jaccard (Tanimoto)",
    "dice": "Dice",
    "rogerstanimoto": "Rogers-Tanimoto",
    "euclidean": "Euclidean",
    "manhattan": "Manhattan",
    "cosine": "Cosine",
    "correlation": "Correlation",
}
BIT_METRICS = frozenset({"jaccard", "dice", "rogerstanimoto"})
REDUCED_DIMS = 50  # PCA dimensions handed to t-SNE and the SOM
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
    reduce_to: int = REDUCED_DIMS,
    init: str = "pca",
    seed: int | None = 42,
    dims: int = 2,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """t-SNE (Barnes-Hut) in *dims* (2 or 3) dimensions, on up to *reduce_to* PCA dimensions.

    *learning_rate* None is scikit-learn's "auto" (rows / exaggeration / 4); *init* is
    "pca" (stable layouts) or "random".
    """
    _need(features, 5)
    from sklearn.manifold import TSNE

    dims = _dims(dims)
    if init not in ("pca", "random"):
        raise ValueError(f"unknown t-SNE start {init!r}")
    x = reduce(features, max(2, reduce_to), progress=progress, cancel=cancel)
    check(cancel)
    perplexity = min(perplexity, (len(features) - 1) / 3)
    model = TSNE(
        n_components=dims,
        perplexity=perplexity,
        learning_rate="auto" if learning_rate is None else learning_rate,
        max_iter=max(250, iterations),
        early_exaggeration=exaggeration,
        init=init,
        random_state=seed,
        n_jobs=-1,
    )
    coords = model.fit_transform(np.asarray(x, dtype=np.float32))
    return _scatter(features, coords, "tsne", dims)


def pacmap_map(
    features: Features,
    *,
    neighbours: int = 10,
    mn_ratio: float = 0.5,
    fp_ratio: float = 2.0,
    distance: str = "euclidean",
    reduce_to: int = REDUCED_DIMS,
    seed: int | None = 42,
    dims: int = 2,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """PaCMAP (Wang et al. 2021, JMLR 22:201): neighbour, mid-near and further pairs, so the
    map keeps both local groups and the overall layout.

    Runs on up to *reduce_to* PCA dimensions (fingerprints are reduced by the streamed PCA, so
    no dense bit matrix is built). *mn_ratio* and *fp_ratio* are the mid-near and further pairs
    per neighbour pair; *distance* is a key of ``PACMAP_METRICS``.
    """
    _need(features, 5)
    dims = _dims(dims)
    if distance not in PACMAP_METRICS:
        raise ValueError(f"unknown PaCMAP distance {distance!r}")
    import pacmap  # imports numba and faiss: only when a map is made

    x = reduce(features, max(2, reduce_to), progress=progress, cancel=cancel)
    check(cancel)
    n = len(features)
    model = pacmap.PaCMAP(
        n_components=dims,
        n_neighbors=max(2, min(neighbours, n - 2)),
        MN_ratio=mn_ratio,
        FP_ratio=fp_ratio,
        distance=distance,
        apply_pca=False,  # already reduced above
        random_state=seed,
    )
    quiet = logging.getLogger("pacmap")  # "random state is set", "high-dimensional data"
    level = quiet.level
    quiet.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coords = model.fit_transform(np.ascontiguousarray(x, dtype=np.float32), init="pca")
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
