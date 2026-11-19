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

from .clustering import scale_columns, standardize
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
UMAP_INITS = {"spectral": "Spectral", "pca": "PCA", "random": "Random"}
#: t-SNE distances for columns, counts and PCA dimensions (bits use Jaccard)
TSNE_METRICS = {
    "euclidean": "Euclidean",
    "cosine": "Cosine",
    "manhattan": "Manhattan",
    "chebyshev": "Chebyshev",
}
TSNE_EXACT_MAX = 5_000  # exact t-SNE is quadratic in time and memory
PACMAP_VARIANTS = {"pacmap": "PaCMAP", "localmap": "LocalMAP"}
GUIDE_GROUPS = 5  # a number column guiding UMAP is cut into this many groups


@dataclass
class Embedding:
    coords: np.ndarray  # n × 2 or n × 3 (NaN for rows without features)
    method: str
    explained: tuple[float, ...] = ()  # PCA: share of variance on each axis
    notes: list[str] = field(default_factory=list)
    #: PCA with more components written: n × k scores of PC1…PCk (NaN without features)
    components: np.ndarray | None = None
    #: PCA of columns alone: d × k loadings of each kept column (see ``kept``) on PC1…PCk
    loadings: np.ndarray | None = None
    kept: np.ndarray | None = None  # per input column: kept (not constant)

    @property
    def placed(self) -> int:
        return int(np.isfinite(self.coords[:, 0]).sum())


@dataclass
class Features:
    """What a map is built from: fingerprints (bits or counts), or scaled number columns
    (complete rows), or both."""

    n: int  # input rows
    rows: np.ndarray  # input rows that have features
    fps: Fingerprints | None = None  # of the rows in *rows*
    matrix: np.ndarray | None = None  # len(rows) × d (counts are uint8)
    kept: np.ndarray | None = None  # columns alone: which input columns were kept

    @classmethod
    def from_fingerprints(cls, fps: Fingerprints) -> Features:
        rows = np.flatnonzero(fps.valid)
        return cls(len(fps), rows, fps=fps.take(rows))

    @classmethod
    def from_counts(cls, counts: np.ndarray, valid: np.ndarray) -> Features:
        """Count fingerprints (``fingerprints.compute_counts``): not bits, so measured by
        Euclidean-type distances only."""
        rows = np.flatnonzero(valid)
        return cls(len(valid), rows, matrix=np.ascontiguousarray(counts[rows]))

    @classmethod
    def from_columns(cls, columns: Sequence[np.ndarray], scaling: str = "zscore") -> Features:
        x, rows, kept = scale_columns(columns, scaling)
        n = len(columns[0]) if columns else 0
        return cls(n, rows, matrix=np.ascontiguousarray(x, dtype=np.float64), kept=kept)

    @classmethod
    def combined(
        cls,
        structures: Fingerprints | Features,
        columns: Sequence[np.ndarray],
        *,
        weight: float = 0.5,
        scaling: str = "zscore",
        dims: int = REDUCED_DIMS,
        progress: Progress | None = None,
        cancel: CancelToken | None = None,
    ) -> Features:
        """Fingerprints (bits, or counts as :class:`Features`) and number columns together,
        for rows that have both.

        The fingerprints are reduced to at most *dims* principal components (streamed, so no
        dense bit matrix), and the columns are scaled (*scaling*, see
        ``clustering.SCALINGS``). Each block is then brought to the same total variance and
        weighted: *weight* is the fingerprints' share (0.5, the default, counts both alike, so
        a 2048-bit fingerprint does not drown out three property columns, nor the reverse).
        Distances on the result are Euclidean.
        """
        if not 0.0 <= weight <= 1.0:
            raise ValueError(f"the fingerprint weight is 0 to 1, not {weight}")
        if isinstance(structures, Fingerprints):
            structures = cls.from_fingerprints(structures)
        x, column_rows = standardize(columns, scaling)
        rows = np.intersect1d(structures.rows, column_rows)
        n = structures.n
        if not len(rows) or not x.shape[1]:
            return cls(n, rows, matrix=np.empty((len(rows), 0)))
        part = structures.subset(np.searchsorted(structures.rows, rows))
        scores = reduce(part, dims, progress=progress, cancel=cancel)
        numbers = x[np.searchsorted(column_rows, rows)]

        def unit(block: np.ndarray, share: float) -> np.ndarray:
            block = np.asarray(block, dtype=np.float64)
            total = float(block.var(axis=0).sum())
            block = block / np.sqrt(total) if total > 1e-12 else block
            return block * np.sqrt(2.0 * share)  # equal shares: each block at variance 1

        blocks = [unit(scores, weight)] if weight > 0 else []
        if weight < 1:
            blocks.append(unit(numbers, 1.0 - weight))
        return cls(n, rows, matrix=np.ascontiguousarray(np.column_stack(blocks)))

    def subset(self, index: np.ndarray) -> Features:
        """The features at positions *index* (into :attr:`rows`)."""
        if self.fps is not None:
            return Features(self.n, self.rows[index], fps=self.fps.take(index))
        return Features(self.n, self.rows[index], matrix=self.matrix[index], kept=self.kept)

    def positions(self, rows: np.ndarray | None) -> np.ndarray | None:
        """Positions (into :attr:`rows`) of the input *rows* that have features; None for
        None, or when they are every row with features (nothing is left to place)."""
        if rows is None:
            return None
        found = np.flatnonzero(np.isin(self.rows, rows))
        return None if len(found) == len(self.rows) else found

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
    fit: np.ndarray | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(scores, explained variance ratios)`` for the rows that have features; with *fit*
    (positions into ``features.rows``), the axes come from those rows alone and every row is
    projected onto them."""
    scores, ratio, _vectors = pca_axes(features, dims, fit=fit, progress=progress, cancel=cancel)
    return scores, ratio


def pca_axes(
    features: Features,
    dims: int = 2,
    *,
    fit: np.ndarray | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(scores, explained variance ratios, d × dims axes)``: see :func:`pca`."""
    basis = features if fit is None else features.subset(fit)
    dims = max(1, min(dims, len(basis) - 1, features.dims))
    if features.fps is None and features.dims <= 200:
        from sklearn.decomposition import PCA

        model = PCA(n_components=dims, svd_solver="full")
        model.fit(np.asarray(basis.matrix, dtype=np.float64))
        scores = model.transform(np.asarray(features.matrix, dtype=np.float64))
        return scores, model.explained_variance_ratio_, model.components_.T
    # Exact PCA from the covariance, summed tile by tile: d × d (2048² for most fingerprints)
    # is small, and the sums are matrix products, so this streams 200k rows quickly.
    tiles = basis.tiles()
    every = features.tiles()
    d = features.dims
    gram = np.zeros((d, d), dtype=np.float64)
    total = np.zeros(d, dtype=np.float64)
    steps = len(tiles) + len(every) + 1
    for i, (start, stop) in enumerate(tiles):
        check(cancel)
        block = basis.dense(start, stop)
        gram += block.T @ block
        total += block.sum(axis=0)
        if progress is not None:
            progress(i + 1, steps)
    n = len(basis)
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
    scores = np.empty((len(features), dims), dtype=np.float64)
    shift = mean @ vectors
    for i, (start, stop) in enumerate(every):
        check(cancel)
        scores[start:stop] = features.dense(start, stop, np.float64) @ vectors - shift
        if progress is not None:
            progress(len(tiles) + 2 + i, steps)
    return scores, values / variance, vectors


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
        x = features.dense(0, None)
    return np.ascontiguousarray(x, dtype=np.float32)


def reduce(
    features: Features,
    dims: int = REDUCED_DIMS,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> np.ndarray:
    """At most *dims* columns per row: PCA scores, or the features themselves if fewer."""
    if features.fps is None and features.dims <= dims:
        return features.matrix
    return pca(features, dims, progress=progress, cancel=cancel)[0]


def _fit_note(features: Features, fit: np.ndarray | None) -> list[str]:
    if fit is None:
        return []
    return [f"fitted on {len(fit):,} of {len(features):,} compounds; the rest placed onto it"]


def _need_fit(fit: np.ndarray | None, minimum: int) -> None:
    if fit is not None and len(fit) < minimum:
        raise ValueError(f"A map needs at least {minimum} compounds to fit on.")


def pca_map(
    features: Features,
    *,
    components: tuple[int, ...] = (1, 2),
    dims: int = 2,
    write: int = 0,
    fit: np.ndarray | None = None,
    **kw,
) -> Embedding:
    """PCA scores on two components (1-based, e.g. ``(1, 3)`` for PC1 against PC3); a 3D map
    (*dims* 3) is always PC1, PC2 and PC3.

    *write* (more than the components drawn) also returns the scores of PC1…PC*write* in
    ``components``; a PCA of columns alone returns their ``loadings`` on those components.
    *fit* (input rows) fits the axes on those rows and places the rest onto them.
    """
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
    fit = features.positions(fit)
    _need_fit(fit, 3)
    last = max(components[-1], write)
    scores, ratio, axes = pca_axes(features, last, fit=fit, **kw)
    if scores.shape[1] < last:  # fewer components than asked: the missing axes are flat
        pad = last - scores.shape[1]
        scores = np.column_stack([scores, np.zeros((len(scores), pad))])
        ratio = np.concatenate([np.asarray(ratio, dtype=float), np.zeros(pad)])
        axes = np.column_stack([axes, np.zeros((len(axes), pad))])
    picked = [c - 1 for c in components]
    explained = tuple(float(ratio[i]) for i in picked)
    parts = [f"PC{c} {e:.0%}" for c, e in zip(components, explained, strict=True)]
    notes = [", ".join(parts[:-1]) + f" and {parts[-1]} of the variance"]
    notes += _fit_note(features, fit)
    extra = {}
    if write > len(components):  # PC1 up to *write*, and always up to the ones drawn
        every = np.full((features.n, last), np.nan)
        every[features.rows] = scores[:, :last]
        extra["components"] = every
    if features.kept is not None:  # columns alone: what each column contributes
        extra["loadings"] = axes[:, :last]
        extra["kept"] = features.kept
    return _scatter(
        features, scores[:, picked], "pca", dims, explained=explained, notes=notes, **extra
    )


def guide_labels(values: Sequence, groups: int = GUIDE_GROUPS) -> np.ndarray:
    """Labels that guide a supervised UMAP, one per row: numbers are cut into *groups*
    groups of about equal size (by quantile), text gives one label per distinct value. A row
    with no value (NaN, None or "") is -1, which UMAP leaves unlabeled."""
    if isinstance(values, np.ndarray) and values.dtype.kind in "fiu":
        x = values.astype(np.float64)
        out = np.full(len(x), -1, dtype=np.int64)
        ok = np.isfinite(x)
        if ok.any():
            edges = np.unique(np.quantile(x[ok], np.linspace(0, 1, max(2, groups) + 1)[1:-1]))
            out[ok] = np.searchsorted(edges, x[ok], side="right")
        return out
    texts = ["" if v is None else str(v) for v in values]
    names = {t: i for i, t in enumerate(sorted({t for t in texts if t}))}
    return np.array([names.get(t, -1) for t in texts], dtype=np.int64)


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
    fit: np.ndarray | None = None,
    guide: np.ndarray | None = None,
    guide_weight: float = 0.5,
    densmap: bool = False,
    density_weight: float = 2.0,
    init: str = "spectral",
    local_connectivity: float = 1.0,
    repulsion: float = 1.0,
    negative_samples: int = 5,
    learning_rate: float = 1.0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """UMAP in *dims* (2 or 3) dimensions. The distance (*metric*) defaults to Jaccard on
    fingerprint bits and Euclidean otherwise; a fingerprint metric needs fingerprint bits
    (see ``UMAP_METRICS``).

    The map is fitted on *fit* (input rows; None: all) and the rest are placed onto it; above
    *fit_max* rows it is fitted on a random sample of those. A *seed* makes the map
    repeatable, at the cost of running on one core. *epochs* None lets UMAP choose (200 for
    large sets, 500 for small).

    *guide* (a label per input row, -1 for none; see :func:`guide_labels`) makes the map
    supervised, pulling rows with the same label together by *guide_weight* (0 to 1).
    *densmap* keeps how densely packed each region is (weighed by *density_weight*); it
    cannot place rows, so it fits every row and takes no *fit* or sample. *init* is
    ``UMAP_INITS``; the rest are UMAP's own settings.
    """
    _need(features, 5)
    if init not in UMAP_INITS:
        raise ValueError(f"unknown UMAP start {init!r}")
    import umap  # slow to import (numba): only when a map is made

    n = len(features)
    neighbours = max(2, min(neighbours, n - 1))
    metric = metric or ("jaccard" if features.fps is not None else "euclidean")
    if metric in BIT_METRICS and features.fps is None:
        raise ValueError(f"{UMAP_METRICS[metric]} distance needs fingerprint bits")
    if min_dist > spread:
        raise ValueError("Min. distance must not be larger than the spread")
    pool = features.positions(fit)
    _need_fit(pool, 5)
    if densmap and pool is not None:
        raise ValueError("densMAP cannot place rows onto a map: fit it on every row")
    if densmap and n > fit_max:
        raise ValueError(f"densMAP fits every row: map at most {fit_max:,} rows")
    pool = np.arange(n) if pool is None else pool
    rng = np.random.default_rng(seed if seed is not None else 0)
    sample = len(pool) > fit_max
    fit_rows = np.sort(rng.choice(pool, fit_max, replace=False)) if sample else pool

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
        target_weight=guide_weight,
        densmap=densmap,
        dens_lambda=density_weight,
        init=init,
        local_connectivity=local_connectivity,
        repulsion_strength=repulsion,
        negative_sample_rate=negative_samples,
        learning_rate=learning_rate,
    )
    labels = None if guide is None else np.asarray(guide, dtype=np.int64)[features.rows]
    check(cancel)
    if progress is not None:
        progress(0, 100)
    out = np.empty((n, dims))
    with warnings.catch_warnings():  # "inverse_transform will be unavailable" (not used)
        warnings.filterwarnings("ignore", message=".*inverse_transform.*")
        # Identical fingerprints are "disconnected" under Jaccard; UMAP still places them.
        warnings.filterwarnings("ignore", message=".*isconnect.*")
        y = None if labels is None else labels[fit_rows]
        out[fit_rows] = model.fit_transform(take(fit_rows), y=y)
    notes = []
    rest = np.setdiff1d(np.arange(n), fit_rows)
    if len(rest):
        notes.append(f"fitted on {len(fit_rows):,} of {n:,} compounds; the rest placed onto it")
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
    metric: str = "euclidean",
    exact: bool = False,
    angle: float = 0.5,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Embedding:
    """t-SNE in *dims* (2 or 3) dimensions: Barnes-Hut, or *exact* (at most
    ``TSNE_EXACT_MAX`` rows). *angle* is Barnes-Hut's trade-off: lower is more accurate and
    slower.

    With *reduce_to* it runs on that many PCA dimensions; without, on the features
    themselves: Jaccard (Tanimoto) on fingerprint bits, *metric* (``TSNE_METRICS``) on
    anything else. *learning_rate* None is scikit-learn's "auto" (rows / exaggeration / 4);
    *init* is "pca" (stable layouts) or "random". t-SNE cannot place rows onto a map, so it
    always fits every row.
    """
    _need(features, 5)
    from sklearn.manifold import TSNE

    dims = _dims(dims)
    if init not in ("pca", "random"):
        raise ValueError(f"unknown t-SNE start {init!r}")
    if metric not in TSNE_METRICS:
        raise ValueError(f"unknown t-SNE distance {metric!r}")
    n = len(features)
    if exact and n > TSNE_EXACT_MAX:
        raise ValueError(f"Exact t-SNE takes at most {TSNE_EXACT_MAX:,} compounds.")
    jaccard = features.fps is not None and not reduce_to
    if jaccard and exact:  # the bits themselves, every pair measured
        x, metric = features.bits(np.arange(n)), "jaccard"
    else:
        x = _features_for(features, reduce_to, progress=progress, cancel=cancel)
    check(cancel)
    perplexity = min(perplexity, (n - 1) / 3)
    if jaccard and not exact:
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
        method="exact" if exact else "barnes_hut",
        angle=angle,
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
    fit: np.ndarray | None = None,
    variant: str = "pacmap",
    learning_rate: float = 1.0,
    iterations: tuple[int, int, int] = (100, 100, 250),
    init: str = "pca",
    local_threshold: float = 10.0,
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

    *variant* "localmap" is LocalMAP (Wang et al. 2025), which sharpens local groups by
    dropping neighbour pairs farther apart than *local_threshold* late in the run.
    *iterations* are the three phases' lengths; *init* is "pca" or "random". With *fit*
    (input rows) the map is fitted on those and the rest placed onto it, which bit distances
    cannot do.
    """
    _need(features, 5)
    dims = _dims(dims)
    if variant not in PACMAP_VARIANTS:
        raise ValueError(f"unknown PaCMAP variant {variant!r}")
    if init not in ("pca", "random"):
        raise ValueError(f"unknown PaCMAP start {init!r}")
    if len(iterations) != 3 or min(iterations) < 0:
        raise ValueError("PaCMAP takes three phase lengths of 0 or more")
    bits = features.fps is not None and not reduce_to
    if distance is None:
        distance = "jaccard" if bits else "euclidean"
    if distance not in PACMAP_METRICS:
        raise ValueError(f"unknown PaCMAP distance {distance!r}")
    if distance in BIT_METRICS and not bits:
        raise ValueError(
            f"{PACMAP_METRICS[distance]} distance needs fingerprint bits without PCA dimensions"
        )
    fit = features.positions(fit)
    _need_fit(fit, 5)
    if fit is not None and distance in BIT_METRICS:
        raise ValueError(
            f"PaCMAP cannot place rows by {PACMAP_METRICS[distance]} distance: choose another "
            "distance to fit on some rows"
        )
    import pacmap  # imports numba and faiss: only when a map is made

    x = _features_for(features, reduce_to, progress=progress, cancel=cancel)
    check(cancel)
    n = len(features)
    basis = x if fit is None else x[fit]
    n_neighbors = max(2, min(neighbours, len(basis) - 2))
    pairs = None
    if distance in BIT_METRICS:
        near = bit_neighbours(features.fps, n_neighbors, distance, cancel=cancel)
        pairs = np.column_stack([np.repeat(np.arange(n), n_neighbors), near.ravel()])
        pairs = np.ascontiguousarray(pairs, dtype=np.int32)
        check(cancel)
    settings = {
        "n_components": dims,
        "n_neighbors": n_neighbors,
        "MN_ratio": mn_ratio,
        "FP_ratio": fp_ratio,
        "pair_neighbors": pairs,
        "distance": "hamming" if distance in BIT_METRICS else distance,
        "apply_pca": False,  # reduced above when asked
        "random_state": seed,
        "lr": learning_rate,
        "num_iters": tuple(int(i) for i in iterations),
        "save_tree": fit is not None,
    }
    if variant == "localmap":
        model = pacmap.LocalMAP(low_dist_thres=local_threshold, **settings)
    else:
        model = pacmap.PaCMAP(**settings)
    quiet = logging.getLogger("pacmap")  # "random state is set", "high-dimensional data"
    level = quiet.level
    quiet.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fitted = model.fit_transform(basis, init=init)
            if fit is None:
                coords = np.asarray(fitted, dtype=np.float64)
            else:
                coords = np.empty((n, dims))
                coords[fit] = fitted
                rest = np.setdiff1d(np.arange(n), fit)
                check(cancel)
                coords[rest] = model.transform(x[rest], basis=basis)
    finally:
        quiet.setLevel(level)
    check(cancel)
    return _scatter(features, coords, "pacmap", dims, notes=_fit_note(features, fit))


OPTIONS = {
    "pca": ("components", "dims", "write", "fit", "progress", "cancel"),
    "umap": (
        "neighbours",
        "min_dist",
        "spread",
        "metric",
        "epochs",
        "seed",
        "fit_max",
        "dims",
        "fit",
        "guide",
        "guide_weight",
        "densmap",
        "density_weight",
        "init",
        "local_connectivity",
        "repulsion",
        "negative_samples",
        "learning_rate",
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
        "metric",
        "exact",
        "angle",
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
        "fit",
        "variant",
        "learning_rate",
        "iterations",
        "init",
        "local_threshold",
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
