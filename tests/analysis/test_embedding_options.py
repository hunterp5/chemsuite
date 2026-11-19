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
"""The settings a map is made with: scaling, fingerprint variants and weight, fitting on some
rows, and each method's own options."""

from __future__ import annotations

import math

import numpy as np
import pytest

from chemsuite.analysis import fingerprints as fp
from chemsuite.analysis.clustering import standardize
from chemsuite.analysis.embedding import Features, embed, guide_labels
from chemsuite.analysis.som import HEX_PITCH, grid_distances2, som

SMILES = [f"{'C' * k}O" for k in range(1, 13)] + [f"c1ccccc1{'C' * k}" for k in range(12)]


def _blobs(per: int = 30, seed: int = 0):
    """Two well-separated groups in 6 dimensions."""
    rng = np.random.default_rng(seed)
    x = np.vstack([rng.normal(0, 0.5, (per, 6)), rng.normal(6, 0.5, (per, 6))])
    return [x[:, j] for j in range(6)], np.repeat([0, 1], per)


def _apart(coords: np.ndarray, labels: np.ndarray) -> float:
    """Distance between the two groups' centers over the mean spread within them."""
    a, b = coords[labels == 0], coords[labels == 1]
    spread = np.mean([np.linalg.norm(g - g.mean(axis=0), axis=1).mean() for g in (a, b)])
    return float(np.linalg.norm(a.mean(axis=0) - b.mean(axis=0)) / spread)


# -- column scaling and the fingerprint weight -------------------------------------------------
def test_column_scalings() -> None:
    column = np.array([1.0, 2.0, 3.0, 4.0, 100.0, np.nan])
    constant = np.full(6, 7.0)
    x, rows = standardize([column, constant])
    assert rows.tolist() == [0, 1, 2, 3, 4] and x.shape == (5, 1)  # constant dropped
    assert x.mean() == pytest.approx(0) and x.std() == pytest.approx(1)
    x, _ = standardize([column], "minmax")
    assert x.min() == 0 and x.max() == 1
    x, _ = standardize([column], "robust")
    assert np.median(x) == 0 and x[1, 0] - x[0, 0] == pytest.approx(0.5)  # IQR is 2
    x, _ = standardize([column], "none")
    assert x[:, 0].tolist() == [1, 2, 3, 4, 100]
    with pytest.raises(ValueError, match="scaling"):
        standardize([column], "log")


def test_fingerprint_weight_shares_the_variance() -> None:
    fps = fp.compute_fingerprints(SMILES)
    logp = np.linspace(0, 5, len(SMILES))
    mw = np.linspace(50, 300, len(SMILES)) ** 1.5

    def shares(weight: float) -> tuple[float, float]:
        m = Features.combined(fps, [logp, mw], weight=weight, dims=4).matrix
        return float(m[:, :-2].var(axis=0).sum()), float(m[:, -2:].var(axis=0).sum())

    assert shares(0.5) == pytest.approx((1.0, 1.0))
    bits, numbers = shares(0.8)
    assert bits / numbers == pytest.approx(4.0)
    assert Features.combined(fps, [logp, mw], weight=1.0, dims=4).matrix.shape[1] == 4
    assert Features.combined(fps, [logp, mw], weight=0.0, dims=4).matrix.shape[1] == 2
    with pytest.raises(ValueError, match="weight"):
        Features.combined(fps, [logp], weight=1.5)


# -- fingerprint size, radius and counts -------------------------------------------------------
def test_fingerprint_variants() -> None:
    assert fp.variant("morgan2") == "morgan2" and fp.variant("morgan2", 2048, 2) == "morgan2"
    key = fp.variant("morgan2", 1024, 3)
    kind = fp.kind_of(key)
    assert key == "morgan2@b1024,r3" and kind.bits == 1024 and kind.column == "ECFP6 1024"
    fps = fp.compute_fingerprints(["CCO", "c1ccccc1"], key)
    assert fps.packed.shape == (2, 128) and fps.valid.all()
    # a bigger radius sets at least the bits a smaller one does, on the same molecule
    small = fp.compute_fingerprints(["CCCCCCO"], fp.variant("morgan2", radius=1)).counts[0]
    large = fp.compute_fingerprints(["CCCCCCO"], fp.variant("morgan2", radius=4)).counts[0]
    assert large > small
    with pytest.raises(ValueError):
        fp.variant("maccs", 1024)
    with pytest.raises(ValueError):
        fp.variant("rdkit", radius=3)
    with pytest.raises(ValueError):
        fp.variant("morgan2", 1000)
    with pytest.raises(ValueError):
        fp.parse("morgan2@x9")


def test_count_fingerprints_map_by_counts() -> None:
    key = fp.variant("morgan2", 512, counts=True)
    counts, valid = fp.compute_counts([*SMILES, "not a smiles"], key)
    assert counts.dtype == np.uint8 and counts.shape == (25, 512)
    assert valid[:-1].all() and not valid[-1]
    assert counts[11].max() >= 2  # a long chain repeats its CH2 environment
    features = Features.from_counts(counts, valid)
    result = embed(features, "pca")
    assert result.placed == 24 and np.isnan(result.coords[-1, 0])
    labels = np.repeat([0, 1], 12)
    assert _apart(result.coords[:-1], labels) > 1.0
    both = Features.combined(features, [np.arange(25.0)], dims=3)
    assert both.matrix.shape == (24, 4)
    with pytest.raises(ValueError, match="bits"):
        embed(features, "umap", metric="jaccard")


# -- fitting on some rows ----------------------------------------------------------------------
def test_pca_fitted_on_some_rows_projects_the_rest() -> None:
    from sklearn.decomposition import PCA

    columns, _labels = _blobs()
    x = np.column_stack(columns)
    fit = np.arange(0, 60, 3)
    result = embed(Features.from_columns(columns), "pca", fit=fit)
    z = (x - x.mean(axis=0)) / x.std(axis=0)
    expected = PCA(2, svd_solver="full").fit(z[fit]).transform(z)
    assert np.allclose(np.abs(result.coords), np.abs(expected))
    assert "fitted on 20 of 60" in result.notes[-1]
    # fingerprints: the streamed PCA does the same
    fps = Features.from_fingerprints(fp.compute_fingerprints(SMILES))
    alone = embed(fps.subset(np.arange(0, 24, 2)), "pca").coords[0::2]
    placed = embed(fps, "pca", fit=np.arange(0, 24, 2)).coords[0::2]
    assert np.allclose(alone, placed)


def test_pca_writes_more_components_and_loadings() -> None:
    columns, _labels = _blobs()
    constant = np.ones(60)
    features = Features.from_columns([*columns, constant])
    result = embed(features, "pca", write=4, components=(1, 2))
    assert result.coords.shape == (60, 2) and result.components.shape == (60, 4)
    assert np.allclose(result.components[:, :2], result.coords)
    assert result.kept.tolist() == [True] * 6 + [False]  # the constant column is dropped
    assert result.loadings.shape == (6, 4)
    assert np.allclose((result.loadings**2).sum(axis=0), 1.0)  # unit axes
    fps = Features.from_fingerprints(fp.compute_fingerprints(SMILES))
    assert embed(fps, "pca").loadings is None  # bit loadings mean nothing


def test_umap_fitted_on_some_rows_and_guided() -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    result = embed(features, "umap", fit=np.arange(0, 60, 2), neighbours=5, seed=1)
    assert result.placed == 60 and "fitted on 30 of 60" in result.notes[0]
    assert _apart(result.coords, labels) > 2.0
    guide = guide_labels(np.where(labels == 1, 9.0, 1.0))
    guided = embed(features, "umap", guide=guide, guide_weight=0.9, neighbours=5, seed=1)
    assert _apart(guided.coords, labels) > 2.0
    other = embed(features, "umap", init="random", repulsion=2.0, negative_samples=3,
                  local_connectivity=2.0, learning_rate=0.5, neighbours=5, seed=1)  # fmt: skip
    assert other.placed == 60
    dense = embed(features, "umap", densmap=True, density_weight=1.0, neighbours=5, seed=1)
    assert dense.placed == 60
    with pytest.raises(ValueError, match="densMAP"):
        embed(features, "umap", densmap=True, fit=np.arange(30))
    with pytest.raises(ValueError, match="densMAP"):
        embed(features, "umap", densmap=True, fit_max=20)
    with pytest.raises(ValueError, match="start"):
        embed(features, "umap", init="tsne")


def test_guide_labels() -> None:
    values = np.array([1.0, 2.0, 3.0, 4.0, np.nan, 5.0, 6.0, 7.0, 8.0])
    labels = guide_labels(values, groups=4)
    assert labels[4] == -1
    assert sorted(set(labels.tolist()) - {-1}) == [0, 1, 2, 3]
    assert (np.diff(labels[values == values]) >= 0).all()  # in order of value
    assert guide_labels(["b", "a", None, "", "b"]).tolist() == [1, 0, -1, -1, 1]


def test_tsne_options() -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    for options in ({"metric": "cosine"}, {"exact": True}, {"angle": 0.8, "metric": "manhattan"}):
        result = embed(features, "tsne", perplexity=10, **options)
        assert _apart(result.coords, labels) > 2.0
    bits = Features.from_fingerprints(fp.compute_fingerprints(SMILES))
    assert embed(bits, "tsne", exact=True, perplexity=5).placed == 24
    with pytest.raises(ValueError, match="distance"):
        embed(features, "tsne", metric="jaccard")
    many = Features(6000, np.arange(6000), matrix=np.zeros((6000, 2)))
    with pytest.raises(ValueError, match="Exact"):
        embed(many, "tsne", exact=True)


def test_pacmap_options_and_fit() -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    result = embed(features, "pacmap", fit=np.arange(0, 60, 2), neighbours=5)
    assert result.placed == 60 and "fitted on 30 of 60" in result.notes[0]
    assert _apart(result.coords, labels) > 2.0
    settings = {"learning_rate": 0.5, "iterations": (50, 50, 100), "init": "random"}
    local = embed(
        features, "pacmap", variant="localmap", local_threshold=5.0, neighbours=5, **settings
    )
    assert _apart(local.coords, labels) > 2.0
    bits = Features.from_fingerprints(fp.compute_fingerprints(SMILES))
    with pytest.raises(ValueError, match="cannot place"):
        embed(bits, "pacmap", fit=np.arange(12))
    assert embed(bits, "pacmap", fit=np.arange(12), distance="euclidean").placed == 24
    with pytest.raises(ValueError, match="variant"):
        embed(features, "pacmap", variant="trimap")


# -- self-organizing maps ----------------------------------------------------------------------
def test_som_grid_distances() -> None:
    square = grid_distances2(3, 4)
    assert square[0, 1] == 1 and square[0, 5] == 2 and square[0, 3] == 9
    torus = grid_distances2(3, 4, wrap=True)
    assert torus[0, 3] == 1 and torus[0, 8] == 1  # round each edge
    hexes = grid_distances2(4, 4, hexagonal=True)
    # cell (0, 0) touches (0, 1) and (1, 0); odd rows sit half a cell right
    assert hexes[0, 1] == pytest.approx(1) and hexes[0, 4] == pytest.approx(1)
    assert hexes[1, 4] == pytest.approx(1) and hexes[0, 5] == pytest.approx(HEX_PITCH**2 + 2.25)
    hex_torus = grid_distances2(4, 4, hexagonal=True, wrap=True)
    assert math.sqrt(hex_torus[0, 12]) == pytest.approx(math.sqrt(HEX_PITCH**2 + 0.25))


def test_som_options() -> None:
    columns, labels = _blobs()
    features = Features.from_columns(columns)
    for options in (
        {"hexagonal": True},
        {"wrap": True},
        {"neighborhood": "bubble", "start_radius": 2.0, "end_radius": 0.5},
        {"init": "random", "seed": 3},
        {"fit": np.arange(0, 60, 2)},
    ):
        result = som(features, 4, 4, epochs=10, **options)
        assert result.counts().sum() == 60
        groups = [set(result.unit[labels == g].tolist()) for g in (0, 1)]
        assert not groups[0] & groups[1], options
    assert som(features, 4, 4, hexagonal=True, wrap=True).hexagonal
    again = som(features, 4, 4, init="random", seed=3)
    assert (again.unit == som(features, 4, 4, init="random", seed=3).unit).all()
    with pytest.raises(ValueError, match="neighborhood"):
        som(features, 4, 4, neighborhood="mexican")


def test_a_wide_bubble_never_collapses_the_map() -> None:
    from chemsuite.analysis.som import som

    rng = np.random.default_rng(0)
    x = np.vstack([rng.normal(-2, 1, (10, 2)), rng.normal(2, 1, (11, 2))])
    features = Features.from_columns([x[:, 0], x[:, 1]], "zscore")
    for wrap in (True, False):
        result = som(features, 3, 3, hexagonal=True, wrap=wrap, neighborhood="bubble",
                     init="random", start_radius=5.0)  # fmt: skip
        assert len(np.unique(result.unit)) > 2  # the cells spread out, not one mean
