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
"""Generative Topographic Mapping: EM, placing rows, landscapes and the stored map."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.analysis import gtm as g
from chemsuite.analysis.embedding import Features
from chemsuite.analysis.fingerprints import compute_fingerprints
from chemsuite.core.jobs import Cancelled, CancelToken


def _blobs(per: int = 40, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    centers = np.array([[0.0] * 6, [8.0] * 6, [0.0, 8.0] * 3])
    x = np.vstack([c + rng.normal(0, 0.5, (per, 6)) for c in centers])
    return x, np.repeat([0, 1, 2], per)


def _features(x: np.ndarray) -> Features:
    return Features.from_columns([x[:, j] for j in range(x.shape[1])], "none")


def test_defaults() -> None:
    assert [g.default_side(n) for n in (50, 299, 300, 2999, 3000, 10**6)] == [
        10, 10, 16, 16, 25, 25,
    ]  # fmt: skip
    assert [g.default_basis(s) for s in (10, 16, 25, 40)] == [3, 4, 6, 10]


@pytest.mark.parametrize("regularization", [0.0, 0.1, 1.0])
def test_em_never_lowers_the_penalized_likelihood(regularization) -> None:
    x, _ = _blobs()
    _w, _b, history = g.train(x, 8, 3, regularization=regularization, iterations=40, tolerance=0.0)
    steps = np.diff(history)
    assert (steps >= -1e-9).all(), steps.min()
    assert history[-1] > history[0]


def test_chunks_change_nothing() -> None:
    x, _ = _blobs()
    whole = g.train(x, 6, 3, iterations=10, tolerance=0.0, chunk=10_000)
    parts = g.train(x, 6, 3, iterations=10, tolerance=0.0, chunk=7)
    np.testing.assert_allclose(whole[0], parts[0], rtol=1e-8, atol=1e-8)
    assert whole[1] == pytest.approx(parts[1])
    a = g.place(x, whole[0], whole[1], 6, 3, chunk=10_000)
    b = g.place(x, whole[0], whole[1], 6, 3, chunk=5)
    for left, right in zip(a, b, strict=True):
        np.testing.assert_allclose(left, right, rtol=1e-6, atol=1e-6)


def test_groups_land_apart_and_responsibilities_sum_to_one() -> None:
    x, labels = _blobs()
    gmap = g.gtm(_features(x), 10)
    assert len(gmap) == len(x) and gmap.side == 10 and gmap.per_side == 3
    centers = np.array([gmap.position[labels == k].mean(axis=0) for k in range(3)])
    d = ((gmap.position[:, None, :] - centers[None]) ** 2).sum(axis=2)
    assert (d.argmin(axis=1) == labels).all()
    rows = np.arange(len(gmap))
    owner, _node, weight = g._expand(gmap, rows)
    np.testing.assert_allclose(np.bincount(owner, weight), 1.0, rtol=1e-5)
    assert (np.abs(gmap.position) <= 1.0 + 1e-9).all()
    assert np.isfinite(gmap.loglik).all()


def test_a_curve_maps_in_order() -> None:
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(0, 3, 150))
    x = np.column_stack([np.cos(t), np.sin(t), 0.3 * t]) * 4 + rng.normal(0, 0.05, (150, 3))
    gmap = g.gtm(_features(x), 12, 4)
    # one latent direction follows the curve: the positions are monotone in t along it
    best = max(abs(np.corrcoef(t, gmap.position[:, j])[0, 1]) for j in range(2))
    assert best > 0.95


def test_sparse_responsibilities_keep_almost_all_the_mass() -> None:
    x, _ = _blobs()
    w, beta, _ = g.train(x, 10, 3, iterations=30)
    *_, indptr, node, weight = g.place(x, w, beta, 10, 3)
    log_r, _ = g._log_resp(x, g.basis(g.latent_grid(10), 3) @ w, beta)
    full = np.exp(log_r)
    for i in range(len(x)):
        kept = node[indptr[i] : indptr[i + 1]].astype(int)
        assert full[i, kept].sum() >= g.KEEP_MASS - 1e-9 or len(kept) == g.KEEP_MAX


def test_landscapes() -> None:
    x, labels = _blobs()
    gmap = g.gtm(_features(x), 10)
    rows = np.arange(len(gmap))
    constant, mass = g.landscape(gmap, rows, np.full(len(rows), 7.0))
    reached = mass > 0
    np.testing.assert_allclose(constant[reached], 7.0)
    assert np.isnan(constant[~reached]).all()
    np.testing.assert_allclose(mass.sum(), len(rows), rtol=1e-5)
    np.testing.assert_allclose(g.density(gmap, rows), mass)
    best, purity, _ = g.class_landscape(gmap, rows, labels, 3)
    for k in range(3):  # each group's own nodes are that group's class
        modes = gmap.mode[labels == k]
        assert (best[modes] == k).all() and (purity[modes] > 0.9).all()
    # skipped rows and NaN values are left out
    half, _ = g.landscape(gmap, np.where(labels == 0, rows, -1), np.where(labels == 0, 1.0, 0.0))
    assert np.nanmax(half) == pytest.approx(1.0)
    values = np.where(labels == 0, np.nan, 1.0)
    _mean, without = g.landscape(gmap, rows, values)
    assert without.sum() == pytest.approx((labels != 0).sum(), rel=1e-5)


def test_fit_rows_train_and_every_row_is_placed() -> None:
    x, labels = _blobs()
    gmap = g.gtm(_features(x), 8, fit=np.flatnonzero(labels < 2))
    assert len(gmap) == len(x)
    # the unseen group is explained worse than the groups trained on
    assert gmap.loglik[labels == 2].mean() < gmap.loglik[labels < 2].mean()


def test_fingerprints_map() -> None:
    smiles = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccccc1C", "c1ccccc1CC", "C1CCCCC1", ""]
    fps = compute_fingerprints(smiles, "morgan2@b512")  # small: a fast PCA
    gmap = g.gtm(Features.from_fingerprints(fps), 4, 2)
    assert gmap.oids.tolist() == list(range(7))  # the blank has no features


def test_round_trip_and_lookup() -> None:
    x, _ = _blobs(per=10)
    gmap = g.gtm(_features(x), 5, 2).with_oids(np.arange(30) * 3 + 100, {"fingerprint": "x"})
    back = g.unpack(gmap.to_bytes())
    for name in ("weights", "oids", "position", "mode", "loglik", "indptr", "node", "weight"):
        np.testing.assert_array_equal(getattr(back, name), getattr(gmap, name))
    assert (back.side, back.per_side, back.beta) == (gmap.side, gmap.per_side, gmap.beta)
    assert back.recipe == {"fingerprint": "x"} and back.history == gmap.history
    assert gmap.nbytes > 0
    assert gmap.lookup(np.array([100, 103, 101, 187, 5])).tolist() == [0, 1, -1, 29, -1]


def test_bad_bytes_and_settings() -> None:
    with pytest.raises(g.GtmError):
        g.unpack(b"not a map")
    x, _ = _blobs(per=5)
    with pytest.raises(g.GtmError):
        g.gtm(_features(x), 1)
    with pytest.raises(g.GtmError):
        g.gtm(_features(x[:1]))


def test_cancel_stops_training() -> None:
    x, _ = _blobs()
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        g.train(x, 6, 3, cancel=token)
