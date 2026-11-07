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

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.analysis.embedding import Features, embed, pca
from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.analysis.som import default_size, som


def _blobs(per: int = 40, seed: int = 0):
    """Three well-separated groups in 8 dimensions, plus one row with a blank."""
    rng = np.random.default_rng(seed)
    centres = np.array([[0.0] * 8, [8.0] * 8, [0.0, 8.0] * 4])
    x = np.vstack([c + rng.normal(0, 0.5, (per, 8)) for c in centres])
    labels = np.repeat([0, 1, 2], per)
    columns = [np.append(x[:, j], np.nan if j == 0 else 1.0) for j in range(8)]
    return columns, labels


def _separated(coords: np.ndarray, labels: np.ndarray) -> bool:
    """Every point is nearer its own group's centre than any other group's."""
    centres = np.array([coords[labels == g].mean(axis=0) for g in range(3)])
    d = ((coords[:, None, :] - centres[None]) ** 2).sum(axis=2)
    return bool((d.argmin(axis=1) == labels).all())


@pytest.mark.parametrize("method", ["pca", "umap", "tsne"])
def test_maps_separate_groups(method) -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    result = embed(features, method, neighbours=10, perplexity=15)
    coords = result.coords
    assert coords.shape == (len(labels) + 1, 2)
    assert np.isnan(coords[-1]).all() and result.placed == len(labels)  # the blank row
    assert _separated(coords[:-1], labels)


def test_pca_explained_and_fingerprints() -> None:
    columns, _ = _blobs()
    result = embed(Features.from_columns(columns), "pca")
    assert result.explained[0] > result.explained[1] > 0
    smiles = ["CCO", "CCO", "CCCO", "c1ccccc1", "Cc1ccccc1", "c1ccncc1", "", "CCCCO"]
    features = Features.from_fingerprints(compute_fingerprints(smiles))
    scores, ratio = pca(features, 2)
    assert scores.shape == (7, 2) and ratio.sum() <= 1.0
    coords = embed(features, "pca").coords
    assert np.allclose(coords[0], coords[1]) and np.isnan(coords[6]).all()


def test_umap_on_fingerprints_with_sampled_fit() -> None:
    from chemsuite.analysis.embedding import umap_map

    smiles = [f"{'C' * k}O" for k in range(1, 13)] + [f"c1ccccc1{'C' * k}" for k in range(12)]
    features = Features.from_fingerprints(compute_fingerprints(smiles))
    result = umap_map(features, neighbours=5, fit_max=16)
    assert result.placed == len(smiles) and result.notes
    labels = np.repeat([0, 1], 12)
    centres = np.array([result.coords[labels == g].mean(axis=0) for g in (0, 1)])
    spread = max(
        np.linalg.norm(result.coords[labels == g] - centres[g], axis=1).mean() for g in (0, 1)
    )
    assert np.linalg.norm(centres[0] - centres[1]) > spread  # alcohols apart from aromatics


def test_self_organising_map() -> None:
    columns, labels = _blobs()
    result = som(Features.from_columns(columns), 6, 6, epochs=15)
    assert result.unit[-1] == -1 and result.counts().sum() == len(labels)
    groups = [set(result.unit[:-1][labels == g].tolist()) for g in range(3)]
    assert not (groups[0] & groups[1]) and not (groups[0] & groups[2])
    r, c = result.cell()
    assert np.isnan(r[-1]) and np.nanmin(r) >= 1 and np.nanmax(c) <= 6
    rows, columns = default_size(3000)  # about 5·√3000 ≈ 274 cells
    assert 220 <= rows * columns <= 330 and rows <= columns


def test_fingerprint_pca_matches_sklearn() -> None:
    from sklearn.decomposition import PCA

    smiles = [
        "CCO",
        "CCCO",
        "CCCCO",
        "c1ccccc1",
        "Cc1ccccc1",
        "c1ccncc1",
        "CC(=O)O",
        "CCN",
        "C1CCCCC1",
    ]
    features = Features.from_fingerprints(compute_fingerprints(smiles))
    scores, ratio = pca(features, 3)
    reference = PCA(n_components=3, svd_solver="full")
    expected = reference.fit_transform(features.dense(dtype=np.float64))
    assert np.allclose(ratio, reference.explained_variance_ratio_, atol=1e-9)
    for k in range(3):  # equal up to the sign of each axis
        assert np.allclose(np.abs(scores[:, k]), np.abs(expected[:, k]), atol=1e-6)


def test_pca_on_other_components() -> None:
    columns, _ = _blobs()
    features = Features.from_columns(columns)
    full, ratio = pca(features, 3)
    result = embed(features, "pca", components=(1, 3))
    assert np.allclose(result.coords[:-1], full[:, [0, 2]])
    assert result.explained == pytest.approx((ratio[0], ratio[2]))
    assert "PC1" in result.notes[0] and "PC3" in result.notes[0]
    with pytest.raises(ValueError, match="increasing"):
        embed(features, "pca", components=(2, 1))


def test_umap_and_tsne_settings() -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    result = embed(
        features, "umap", neighbours=10, metric="manhattan", spread=2.0, min_dist=0.5, epochs=100
    )
    assert _separated(result.coords[:-1], labels)
    with pytest.raises(ValueError, match="needs fingerprints"):
        embed(features, "umap", metric="dice")
    with pytest.raises(ValueError, match="spread"):
        embed(features, "umap", min_dist=0.9, spread=0.5)
    result = embed(
        features,
        "tsne",
        perplexity=15,
        learning_rate=100.0,
        iterations=500,
        exaggeration=8.0,
        reduce_to=4,
        init="random",
    )
    assert _separated(result.coords[:-1], labels)
