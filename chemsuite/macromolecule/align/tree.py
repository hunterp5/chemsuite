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
"""The guide tree of a multiple alignment and the sequence weights it gives.

Distances come from all-pairs alignments when there are few sequences (identity, with Kimura's
correction for multiple substitutions, as ClustalW: Thompson, Higgins & Gibson, Nucleic Acids
Res. 1994, 22, 4673) and from shared k-mers on a compressed alphabet when there are many
(Edgar, Nucleic Acids Res. 2004, 32, 1792). The tree is UPGMA; each sequence's weight is
ClustalW's: the branch lengths from it to the root, each shared among the leaves below it, so
a cluster of near-copies counts about as much as one distinct sequence.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import linkage

#: above this many sequences, distances come from k-mers rather than alignments
ALIGNED_DISTANCES = 120
#: Dayhoff's six groups of amino acids (AGPST, C, DENQ, FWY, HKR, ILMV), for k-mer counting
_DAYHOFF = {**dict.fromkeys("AGPST", 0), "C": 1, **dict.fromkeys("DENQBZ", 2),
            **dict.fromkeys("FWY", 3), **dict.fromkeys("HKR", 4), **dict.fromkeys("ILMV", 5)}  # fmt: skip
_K = 4


@dataclass(frozen=True)
class GuideTree:
    #: scipy's linkage: row k joins clusters ``merges[k, 0]`` and ``merges[k, 1]`` into ``n + k``
    merges: np.ndarray
    weights: np.ndarray  # one per sequence, summing to the number of sequences


def guide_tree(sequences: Sequence[str], pair_distance: Callable[[int, int], float]) -> GuideTree:
    """UPGMA over ``pair_distance`` (few sequences) or k-mer distances (many)."""
    n = len(sequences)
    if n < 2:
        return GuideTree(np.zeros((0, 4)), np.ones(n))
    if n <= ALIGNED_DISTANCES:
        condensed = np.array([pair_distance(i, j) for i in range(n) for j in range(i + 1, n)])
    else:
        condensed = kmer_distances(sequences)
    merges = linkage(condensed, method="average")
    return GuideTree(merges, clustal_weights(merges, n))


def kimura(identity: float) -> float:
    """Kimura's protein distance from the share of identical aligned residues."""
    d = 1.0 - identity
    inner = 1.0 - d - 0.2 * d * d
    return -np.log(inner) if inner > 0.01 else 10.0


def kmer_distances(sequences: Sequence[str]) -> np.ndarray:
    """Condensed distances: 1 - the share of 4-mers (Dayhoff alphabet) two sequences share."""
    bins = 6**_K
    counts = np.zeros((len(sequences), bins), np.float32)
    lengths = np.zeros(len(sequences))
    for s, seq in enumerate(sequences):
        code = np.array([_DAYHOFF.get(ch, 5) for ch in seq.upper()], np.int64)
        if len(code) < _K:
            continue
        index = np.zeros(len(code) - _K + 1, np.int64)
        for k in range(_K):
            index = index * 6 + code[k : len(code) - _K + 1 + k]
        counts[s] = np.bincount(index, minlength=bins)
        lengths[s] = len(index)
    out = []
    for i in range(len(sequences) - 1):
        shared = np.minimum(counts[i], counts[i + 1 :]).sum(axis=1)
        shortest = np.maximum(np.minimum(lengths[i], lengths[i + 1 :]), 1.0)
        out.append(1.0 - shared / shortest)
    return np.concatenate(out)


def clustal_weights(merges: np.ndarray, n: int) -> np.ndarray:
    """Each leaf's share of the branch lengths above it (UPGMA heights are half distances)."""
    height = np.zeros(2 * n - 1)
    size = np.ones(2 * n - 1)
    parent = np.full(2 * n - 1, -1)
    for k, (a, b, dist, count) in enumerate(merges):
        node = n + k
        height[node] = dist / 2.0
        size[node] = count
        parent[int(a)] = parent[int(b)] = node
    weights = np.zeros(n)
    for leaf in range(n):
        node = leaf
        while parent[node] >= 0:
            up = parent[node]
            weights[leaf] += (height[up] - height[node]) / size[node]
            node = up
    if weights.sum() <= 0:
        return np.ones(n)
    weights = np.maximum(weights, 1e-3 * weights.max())
    return weights * (n / weights.sum())


def leaves(merges: np.ndarray, n: int, node: int) -> list[int]:
    """The sequences under a node of the tree, in order."""
    if node < n:
        return [node]
    a, b = int(merges[node - n, 0]), int(merges[node - n, 1])
    return leaves(merges, n, a) + leaves(merges, n, b)
