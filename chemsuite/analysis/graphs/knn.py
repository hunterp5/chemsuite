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

"""The *k* most similar compounds of every compound, by Tanimoto on packed fingerprints.

Up to :data:`EXACT_MAX` rows the neighbours are exact (``embedding.bit_neighbours``, whose time
grows with the square of the rows). Larger sets use MinHash and a forest of sorted prefixes,
after the LSH Forest of TMAP (Probst, D. and Reymond, J.-L. (2020). "Visualization of very
large high-dimensional data sets as minimum spanning trees." *J. Cheminf.* 12, 12; Bawa, M.,
Condie, T. and Ganesan, P. (2005). "LSH Forest: self-tuning indexes for similarity search."
*WWW '05*, 651–660):

- **MinHash.** Each of *perms* random orderings of the fingerprint's bits gives a row the
  rank of its first set bit. Two rows agree on one ordering with probability equal to their
  Tanimoto (Jaccard) similarity.
- **Prefix forest.** The orderings are split into *trees* bands. In each band the rows are
  sorted by their ranks, in order, so rows sharing a long prefix sit side by side, as leaves
  of the same branch of a prefix tree would. The rows within *window* places of a row, in any
  band, are its candidates.
- **Exact rerank.** Each candidate's Tanimoto is computed from the packed bits, and the best
  *k* kept, so the similarities are those every other tool in chemsuite reports.

MHFP6 (:mod:`.mhfp`) is a MinHash already: its values go straight into the forest at every
size, and candidates are ranked by the share of values two rows agree on.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from chemsuite.analysis.embedding import bit_neighbours
from chemsuite.analysis.fingerprints import Fingerprints, from_counts
from chemsuite.core.jobs import CancelToken, check

from .mhfp import MinHashes, agreement

Progress = Callable[[int, int], None]

#: more valid rows than this: approximate neighbours (exact ones take ~7 s at 10,000 rows)
EXACT_MAX = 10_000
PERMS = 128
TREES = 16
WINDOW = 10
REFINE = 1
CHUNK_ROWS = 1024  # MinHash rows at once: one numpy call holds the GIL, ~25 ms at this size
RERANK_CELLS = 200_000  # candidate pairs reranked at once (× 256 bytes of fingerprint words)
SEED = 42


def knn(
    fps: Fingerprints | MinHashes,
    k: int,
    *,
    exact_max: int = EXACT_MAX,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """``(index, similarity)``, each ``n × k``: every row's *k* most similar rows, most similar
    first. Invalid rows, and the places a row has no neighbour for, hold -1 and NaN. MHFP
    rows always go through the forest (all their pairs would take far longer than bits')."""
    n = len(fps)
    index = np.full((n, k), -1, dtype=np.int64)
    sim = np.full((n, k), np.nan)
    ok = np.flatnonzero(fps.valid)
    if len(ok) < 2 or k < 1:
        return index, sim
    sub = fps.take(ok)
    kk = min(k, len(ok) - 1)
    if isinstance(sub, Fingerprints) and len(ok) <= exact_max:
        found, scores = bit_neighbours(
            sub, kk, "jaccard", scores=True, progress=progress, cancel=cancel
        )
    else:
        found, scores = lsh_neighbours(sub, kk, seed=seed, progress=progress, cancel=cancel)
    has = found >= 0
    index[ok, :kk] = np.where(has, ok[np.maximum(found, 0)], -1)
    sim[ok, :kk] = np.where(has, scores, np.nan)
    return index, sim


def knn_edges(index: np.ndarray, sim: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The neighbour lists as undirected edges ``(a, b, similarity)`` with ``a < b``, each
    pair once."""
    n, k = index.shape
    rows = np.repeat(np.arange(n), k)
    cols = index.ravel()
    sims = sim.ravel()
    keep = cols >= 0
    rows, cols, sims = rows[keep], cols[keep], sims[keep]
    a, b = np.minimum(rows, cols), np.maximum(rows, cols)
    keep = a != b
    a, b, sims = a[keep], b[keep], sims[keep]
    _, first = np.unique(a * n + b, return_index=True)
    return a[first], b[first], sims[first]


# -- MinHash and the prefix forest ---------------------------------------------------------------
def minhash(
    fps: Fingerprints, perms: int = PERMS, *, seed: int = SEED, cancel: CancelToken | None = None
) -> np.ndarray:
    """``n × perms`` ranks (``uint16``): per ordering of the bits, each row's first set bit.
    A row with no bits set ranks past every bit in every ordering."""
    bits = fps.kind.bits
    rng = np.random.default_rng(seed)
    order = np.stack([rng.permutation(bits) for _ in range(perms)]).astype(np.uint16)
    out = np.full((len(fps), perms), bits, dtype=np.uint16)
    for start in range(0, len(fps), CHUNK_ROWS):
        check(cancel)
        stop = min(len(fps), start + CHUNK_ROWS)
        rows, cols = np.nonzero(np.unpackbits(fps.packed[start:stop], axis=1, count=bits))
        if not len(rows):
            continue
        starts = np.flatnonzero(np.r_[True, rows[1:] != rows[:-1]])
        out[start + rows[starts]] = np.minimum.reduceat(order[:, cols], starts, axis=1).T
    return out


def _candidates(
    hashes: np.ndarray, trees: int, window: int, cancel: CancelToken | None
) -> np.ndarray:
    """``n × (trees · 2 · window)`` candidate rows (-1 past either end of a band)."""
    n, perms = hashes.shape
    width = perms // trees
    out = np.full((n, trees, 2 * window), -1, dtype=np.int32)  # rows fit: half the memory
    offsets = np.r_[np.arange(-window, 0), np.arange(1, window + 1)]
    for t in range(trees):
        check(cancel)
        band = hashes[:, t * width : (t + 1) * width]
        order = np.lexsort(band.T[::-1])  # by the band's first rank, then its second, …
        place = np.empty(n, dtype=np.int64)
        place[order] = np.arange(n)
        at = place[:, None] + offsets[None, :]
        inside = (at >= 0) & (at < n)
        out[:, t, :] = np.where(inside, order[np.clip(at, 0, n - 1)], -1)
    return out.reshape(n, -1)


def _scorer(fps: Fingerprints | MinHashes) -> tuple[Callable, int]:
    """``(score(start, stop, rows) -> similarities, bytes per row)``: exact Tanimoto for bits,
    the share of agreeing values for MHFP."""
    if isinstance(fps, MinHashes):
        hashes = fps.hashes

        def agree(start: int, stop: int, rows: np.ndarray) -> np.ndarray:
            return agreement(hashes[start:stop, None, :], hashes[rows])

        return agree, hashes.shape[1] * hashes.itemsize
    words, counts = fps.words, fps.counts

    def tanimoto(start: int, stop: int, rows: np.ndarray) -> np.ndarray:
        common = np.bitwise_count(words[start:stop, None, :] & words[rows]).sum(axis=2)
        return from_counts(common, counts[start:stop, None], counts[rows], "tanimoto")

    return tanimoto, words.shape[1] * words.itemsize


def _rerank(
    fps: Fingerprints | MinHashes,
    candidates: np.ndarray,
    k: int,
    progress: Progress | None,
    cancel: CancelToken | None,
) -> tuple[np.ndarray, np.ndarray]:
    """The best *k* of each row's *candidates* (``n × c``, -1 for none) by exact Tanimoto (or,
    for MHFP, MinHash agreement): ``(index, similarity)``, most similar first, -1 and NaN where
    there are fewer."""
    n = len(fps)
    index = np.full((n, k), -1, dtype=np.int64)
    sim = np.full((n, k), np.nan)
    score_of, row_bytes = _scorer(fps)
    pairs = max(1, RERANK_CELLS * 256 // max(row_bytes, 1))  # the same memory as 2048 bits
    step = max(1, pairs // max(candidates.shape[1], 1))
    chunks = range(0, n, step)
    for done, start in enumerate(chunks):
        check(cancel)
        stop = min(n, start + step)
        cand = np.sort(candidates[start:stop], axis=1)
        repeat = np.zeros(cand.shape, dtype=bool)
        repeat[:, 1:] = cand[:, 1:] == cand[:, :-1]
        cand[repeat | (cand == np.arange(start, stop)[:, None])] = -1
        score = score_of(start, stop, np.maximum(cand, 0))
        score[cand < 0] = -1.0
        take = min(k, score.shape[1])
        top = np.argpartition(-score, take - 1, axis=1)[:, :take]
        top_score = np.take_along_axis(score, top, axis=1)
        rank = np.argsort(-top_score, axis=1, kind="stable")
        top = np.take_along_axis(top, rank, axis=1)
        top_score = np.take_along_axis(top_score, rank, axis=1)
        found = np.take_along_axis(cand, top, axis=1)
        missing = top_score < 0
        index[start:stop, :take] = np.where(missing, -1, found)
        sim[start:stop, :take] = np.where(missing, np.nan, top_score)
        if progress is not None:
            progress(done + 1, len(chunks))
    return index, sim


def lsh_neighbours(
    fps: Fingerprints | MinHashes,
    k: int,
    *,
    perms: int = PERMS,
    trees: int = TREES,
    window: int = WINDOW,
    refine: int = REFINE,
    seed: int = SEED,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Approximate :func:`knn` for valid fingerprints: ``(index, similarity)``, ``n × k``,
    most similar first, -1 and NaN where fewer than *k* candidates were found.

    After the forest, *refine* rounds look among the neighbours of each row's neighbours (a
    neighbour of a neighbour is often a neighbour: Dong, W., Moses, C. and Li, K. (2011).
    "Efficient k-nearest neighbor graph construction for generic similarity measures."
    *WWW '11*, 577–586)."""
    hashes = fps.hashes if isinstance(fps, MinHashes) else minhash(fps, perms, seed=seed,
                                                                   cancel=cancel)  # fmt: skip
    candidates = _candidates(hashes, trees, window, cancel)
    del hashes
    index, sim = _rerank(fps, candidates, k, None, cancel)
    for _ in range(refine):
        n = len(fps)
        hops = index[np.maximum(index, 0)].reshape(n, -1)
        hops[np.repeat(index < 0, k, axis=1)] = -1
        index, sim = _rerank(fps, np.concatenate([index, hops], axis=1), k, None, cancel)
    if progress is not None:
        progress(1, 1)
    return index, sim
