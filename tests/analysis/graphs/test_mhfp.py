"""MHFP6 values and the neighbours TMAP finds with them."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMHFPFingerprint

from chemsuite.analysis.graphs.knn import knn
from chemsuite.analysis.graphs.mhfp import RADIUS, SEED, MinHashes, agreement, compute_mhfp
from chemsuite.analysis.graphs.tmap import tmap
from chemsuite.core.jobs import Cancelled, CancelToken

PAIRS = [
    ("c1ccccc1O", "c1ccccc1N"),
    ("CC(=O)Nc1ccc(O)cc1", "CC(=O)Nc1ccc(OC)cc1"),
    ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "CC(C)Cc1ccc(CC(=O)O)cc1"),
    ("CCO", "c1ccc2ccccc2c1"),
]


def shingles(smiles: str) -> set[bytes]:
    encoder = rdMHFPFingerprint.MHFPEncoder(256, SEED)
    return set(encoder.CreateShinglingFromMol(Chem.MolFromSmiles(smiles), radius=RADIUS))


def test_agreement_estimates_the_jaccard_of_the_shingles() -> None:
    values = compute_mhfp([s for pair in PAIRS for s in pair])
    for i, (a, b) in enumerate(PAIRS):
        sa, sb = shingles(a), shingles(b)
        jaccard = len(sa & sb) / len(sa | sb)
        estimate = float(agreement(values.hashes[2 * i], values.hashes[2 * i + 1]))
        assert estimate == pytest.approx(jaccard, abs=0.1), (a, b)  # 3σ at 256 values


def test_blank_and_invalid_cells_are_not_valid_and_values_repeat() -> None:
    cells = ["c1ccccc1O", "", None, "not a smiles", "CCO"]
    first, second = compute_mhfp(cells), compute_mhfp(cells, chunk_rows=2)
    assert first.valid.tolist() == [True, False, False, False, True]
    assert np.array_equal(first.hashes, second.hashes)
    assert first.hashes.dtype == np.uint32 and first.hashes.shape == (5, 256)
    assert not first.hashes[1:4].any()


def test_cancelling_stops_the_encoding() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        compute_mhfp(["CCO"] * 10, cancel=token)


def clustered(n: int, perms: int = 256, clusters: int = 30, seed: int = 0) -> MinHashes:
    """MinHash rows in families: each a random centre with some values redrawn."""
    rng = np.random.default_rng(seed)
    centres = rng.integers(0, 2**32, (clusters, perms), dtype=np.uint32)
    rows = centres[rng.integers(0, clusters, n)]
    redrawn = rng.random(rows.shape) < rng.uniform(0.1, 0.6, (n, 1))
    rows = np.where(redrawn, rng.integers(0, 2**32, rows.shape, dtype=np.uint32), rows)
    return MinHashes(rows, np.ones(n, bool))


def brute(values: MinHashes) -> np.ndarray:
    sim = agreement(values.hashes[:, None, :], values.hashes[None, :, :])
    np.fill_diagonal(sim, -1.0)
    return sim


def test_mhfp_neighbours_find_most_of_the_exact_ones() -> None:
    values = clustered(2000)
    index, sim = knn(values, 10)
    exact = np.sort(brute(values), axis=1)[:, ::-1][:, :10]
    recall = np.mean(sim >= exact[:, -1:] - 1e-9)  # by similarity: ties are not misses
    assert recall >= 0.85
    found = brute(values)[np.arange(2000)[:, None], index]
    assert np.allclose(found, sim)  # the similarities reported are the agreements


def test_a_few_mhfp_rows_get_their_exact_neighbours() -> None:
    values = clustered(12, clusters=3)
    values = MinHashes(values.hashes, np.r_[np.ones(11, bool), False])
    index, sim = knn(values, 3)
    expected = np.sort(brute(values.take(np.arange(11))), axis=1)[:, ::-1][:, :3]
    assert np.allclose(sim[:11], expected)
    assert (index[11] == -1).all() and np.isnan(sim[11]).all()


def test_tmap_of_molecules_with_mhfp() -> None:
    smiles = ["c1ccccc1O", "c1ccccc1N", "c1ccccc1C", "CCO", "CCCO", "CCCCO", "not a smiles"]
    tree = tmap(compute_mhfp(smiles), k=3)
    assert tree.placed.tolist() == [True] * 6 + [False]
    assert len(tree.a) == 6 - (tree.component.max() + 1)
    assert ((tree.similarity >= 0) & (tree.similarity <= 1)).all()
    assert np.array_equal(tmap(compute_mhfp(smiles), k=3).coords, tree.coords, equal_nan=True)
