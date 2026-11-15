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
"""Neighbours, spanning forests, layouts and TMAP."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

from chemsuite.analysis.fingerprints import KINDS, Fingerprints, compute_fingerprints, pack
from chemsuite.analysis.graphs.knn import knn, knn_edges, lsh_neighbours, minhash
from chemsuite.analysis.graphs.layout import force, layout, multilevel
from chemsuite.analysis.graphs.layout import pack as pack_layout
from chemsuite.analysis.graphs.tmap import tmap
from chemsuite.analysis.graphs.trees import components, spanning_forest
from chemsuite.core.jobs import Cancelled, CancelToken

KIND = KINDS["morgan2"]


def clustered(n: int, seed: int = 0, clusters: int = 20) -> Fingerprints:
    """Fingerprints in *clusters* families: each a random centre with a few bits flipped."""
    rng = np.random.default_rng(seed)
    centres = rng.random((clusters, KIND.bits)) < 0.03
    bits = centres[rng.integers(0, clusters, n)].copy()
    flips = rng.random(bits.shape) < 0.01
    bits ^= flips
    return Fingerprints(KIND, pack(bits, KIND), np.ones(n, bool))


def brute(fps: Fingerprints) -> np.ndarray:
    words = fps.words
    common = np.bitwise_count(words[:, None, :] & words[None, :, :]).sum(axis=2)
    union = fps.counts[:, None] + fps.counts[None, :] - common
    sim = np.where(union > 0, common / np.maximum(union, 1), 0.0)
    np.fill_diagonal(sim, -1.0)
    return sim


# -- neighbours -----------------------------------------------------------------------------------
def test_exact_neighbours_match_brute_force() -> None:
    fps = clustered(300)
    index, sim = knn(fps, 5)
    expected = np.sort(brute(fps), axis=1)[:, ::-1][:, :5]
    assert np.allclose(sim, expected)
    assert np.allclose(brute(fps)[np.arange(300)[:, None], index], sim)


def test_invalid_rows_have_no_neighbours_and_are_never_one() -> None:
    fps = clustered(50)
    valid = fps.valid.copy()
    valid[[3, 7]] = False
    fps = Fingerprints(KIND, fps.packed, valid)
    index, sim = knn(fps, 4)
    assert (index[[3, 7]] == -1).all() and np.isnan(sim[[3, 7]]).all()
    assert not np.isin(index, [3, 7]).any()


def test_approximate_neighbours_find_most_of_the_exact_ones() -> None:
    fps = clustered(3000, clusters=60)
    _, approx = lsh_neighbours(fps, 10)
    exact = np.sort(brute(fps), axis=1)[:, ::-1][:, :10]
    recall = np.mean(np.nan_to_num(approx, nan=-1) >= exact[:, -1:] - 1e-9)
    assert recall >= 0.85


def test_knn_switches_to_the_approximate_search_above_the_limit() -> None:
    fps = clustered(400)
    index, sim = knn(fps, 5, exact_max=100)
    assert (index >= 0).all() and np.all(np.diff(sim, axis=1) <= 1e-12)


def test_minhash_agreement_estimates_tanimoto() -> None:
    fps = clustered(40, clusters=4)
    hashes = minhash(fps, 512)
    agree = (hashes[:, None, :] == hashes[None, :, :]).mean(axis=2)
    sim = brute(fps)
    off = ~np.eye(40, dtype=bool)
    assert np.abs(agree[off] - sim[off]).mean() < 0.05


def test_knn_edges_are_each_pair_once() -> None:
    index = np.array([[1, 2], [0, -1], [0, 1]])
    sim = np.array([[0.9, 0.5], [0.9, np.nan], [0.5, 0.4]])
    a, b, s = knn_edges(index, sim)
    assert list(zip(a.tolist(), b.tolist(), s.tolist(), strict=True)) == [
        (0, 1, 0.9),
        (0, 2, 0.5),
        (1, 2, 0.4),
    ]


def test_cancelling_the_search_stops_it() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        lsh_neighbours(clustered(100), 5, cancel=token)


# -- trees ----------------------------------------------------------------------------------------
def test_spanning_forest_matches_scipy_and_keeps_zero_distances() -> None:
    rng = np.random.default_rng(1)
    n = 60
    a, b = rng.integers(0, n, 300), rng.integers(0, n, 300)
    keep = a != b
    a, b = np.minimum(a[keep], b[keep]), np.maximum(a[keep], b[keep])
    _, first = np.unique(a * n + b, return_index=True)
    a, b = a[first], b[first]
    w = rng.random(len(a)).round(1)  # ties and zeros
    ta, tb, tw = spanning_forest(n, a, b, w)
    labels = components(n, a, b)
    assert len(ta) == n - (labels.max() + 1)
    reference = minimum_spanning_tree(coo_matrix((w + 1.0, (a, b)), shape=(n, n)).tocsr())
    assert tw.sum() == pytest.approx(reference.sum() - len(ta))


def test_components_are_numbered_largest_first() -> None:
    labels = components(7, np.array([0, 4, 5]), np.array([6, 5, 6]))
    assert labels.tolist() == [0, 1, 2, 3, 0, 0, 0]


# -- layouts --------------------------------------------------------------------------------------
def test_edges_end_near_their_length_and_nodes_apart() -> None:
    rng = np.random.default_rng(2)
    n = 500
    parent = np.array([rng.integers(0, i) for i in range(1, n)])
    a, b = parent, np.arange(1, n)
    pos = multilevel(n, a, b)
    length = np.hypot(*(pos[a] - pos[b]).T)
    assert 0.3 < np.median(length) < 3.0
    gaps = np.hypot(*(pos[:, None] - pos[None, :]).transpose(2, 0, 1))
    np.fill_diagonal(gaps, np.inf)
    assert np.median(gaps.min(axis=1)) > 0.2


def test_layouts_are_the_same_for_the_same_seed() -> None:
    a, b = np.arange(0, 199), np.arange(1, 200)
    assert np.array_equal(multilevel(200, a, b, seed=3), multilevel(200, a, b, seed=3))


def test_force_keeps_a_pair_at_its_edge_length() -> None:
    pos = force(np.array([[0.0, 0.0], [5.0, 0.0]]), np.array([0]), np.array([1]), k=1.0)
    assert np.hypot(*(pos[0] - pos[1])) == pytest.approx(1.0, abs=0.05)


def test_packed_components_do_not_overlap() -> None:
    rng = np.random.default_rng(4)
    labels = np.repeat(np.arange(6), [40, 20, 10, 5, 1, 1])
    pos = pack_layout(rng.normal(size=(len(labels), 2)) * 3, labels)
    boxes = [(pos[labels == c].min(axis=0), pos[labels == c].max(axis=0)) for c in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            (lo1, hi1), (lo2, hi2) = boxes[i], boxes[j]
            assert (hi1 < lo2).any() or (hi2 < lo1).any()


def test_layout_places_every_component() -> None:
    a, b = np.array([0, 1, 3]), np.array([1, 2, 4])
    labels = components(6, a, b)
    pos = layout(6, a, b, labels)
    assert np.isfinite(pos).all() and len(np.unique(pos.round(6), axis=0)) == 6


# -- TMAP -----------------------------------------------------------------------------------------
def test_tmap_of_molecules() -> None:
    smiles = ["c1ccccc1O", "c1ccccc1N", "c1ccccc1C", "CCO", "CCCO", "CCCCO", "not a smiles"]
    tree = tmap(compute_fingerprints(smiles), k=3)
    assert tree.placed.tolist() == [True] * 6 + [False]
    assert tree.component[6] == -1
    assert len(tree.a) == 6 - (tree.component.max() + 1)
    assert ((tree.similarity > 0) & (tree.similarity <= 1)).all()
