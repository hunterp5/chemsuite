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
from rdkit import Chem

from chemsuite.analysis.fingerprints import compute_fingerprints, query_fingerprint, similarity
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.similarity import model as sm
from chemsuite.similarity.fusion import consensus, percentile_ranks, unit
from chemsuite.similarity.run import score_cells

CELLS = [
    "CC(=O)Nc1ccc(O)cc1", "CC(=O)Nc1ccc(OC)cc1", "CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1O",
    "CN1CCC[C@H]1c1cccnc1", "c1ccccc1", None, "not a smiles",
]  # fmt: skip
PARACETAMOL = sm.Reference(Chem.MolFromSmiles(CELLS[0]), "paracetamol")
ASPIRIN = sm.Reference(Chem.MolFromSmiles(CELLS[2]), "aspirin")


def _single(smiles: str, key: str = "morgan2", metric: str = "tanimoto") -> np.ndarray:
    fps = compute_fingerprints(CELLS, key)
    return similarity(fps, query_fingerprint(Chem.MolFromSmiles(smiles), key), metric)


def test_settings_are_checked() -> None:
    with pytest.raises(ValueError):
        sm.Method(kind="nope")
    with pytest.raises(ValueError):
        sm.Method(metric="nope")
    with pytest.raises(ValueError):
        sm.Method(metric="tversky", alpha=0.0, beta=0.0)
    with pytest.raises(ValueError):
        sm.Search(())
    with pytest.raises(ValueError):
        sm.Search((sm.Method(), sm.Method()))  # the same column twice
    assert sm.Method().column == "ECFP4 Tanimoto"
    assert sm.Method(fingerprint="pharm2d").column == "Pharmacophore 2D Tanimoto"
    assert sm.Method(metric="tversky").column == "ECFP4 Tversky"
    assert sm.Method(kind=sm.SHAPE).expensive and not sm.Method().expensive


def test_one_reference_matches_plain_similarity() -> None:
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(),)))
    expected = _single(CELLS[0])
    got = scores.columns["ECFP4 Tanimoto"]
    assert np.allclose(got[:6], expected[:6]) and got[0] == pytest.approx(1.0)
    assert np.isnan(got[6:]).all()
    assert scores.nearest.tolist() == [0] * 6 + [-1, -1]
    assert scores.consensus is None and "morgan2" in scores.fingerprints


def test_several_references_fuse_best_or_mean() -> None:
    refs = [PARACETAMOL, ASPIRIN]
    each = np.column_stack([_single(CELLS[0]), _single(CELLS[2])])
    best = score_cells(CELLS, refs, sm.Search((sm.Method(),), fuse="max"))
    mean = score_cells(CELLS, refs, sm.Search((sm.Method(),), fuse="mean"))
    assert np.allclose(best.columns["ECFP4 Tanimoto"][:6], each[:6].max(axis=1))
    assert np.allclose(mean.columns["ECFP4 Tanimoto"][:6], each[:6].mean(axis=1))
    assert best.nearest[:6].tolist() == each[:6].argmax(axis=1).tolist()


def test_a_reference_that_is_a_row_skips_itself() -> None:
    refs = [sm.Reference(Chem.MolFromSmiles(CELLS[0]), "row 1", row=0),
            sm.Reference(Chem.MolFromSmiles(CELLS[2]), "row 3", row=2)]  # fmt: skip
    scores = score_cells(CELLS, refs, sm.Search((sm.Method(),)))
    values = scores.columns["ECFP4 Tanimoto"]
    assert values[0] == pytest.approx(_single(CELLS[2])[0])  # its nearest other reference
    assert scores.nearest[0] == 1 and scores.nearest[2] == 0


def test_tversky_and_pharmacophore_methods() -> None:
    method = sm.Method(metric="tversky", alpha=1.0, beta=0.0)
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search((method,)))
    assert scores.columns["ECFP4 Tversky"][0] == pytest.approx(1.0)
    pharm = sm.Method(fingerprint="pharm2d")
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search((pharm,)))
    assert np.allclose(scores.columns[pharm.column][:6], _single(CELLS[0], "pharm2d")[:6])


def test_consensus_over_methods() -> None:
    methods = (sm.Method(), sm.Method(fingerprint="atompair"))
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search(methods, consensus="mean"))
    a, b = (scores.columns[m.column] for m in methods)
    assert np.allclose(scores.consensus[:6], (a[:6] + b[:6]) / 2)
    assert np.isnan(scores.consensus[6:]).all()
    ranked = score_cells(CELLS, [PARACETAMOL], sm.Search(methods, consensus="rank"))
    assert ranked.consensus[0] == pytest.approx(1.0)  # top by both


def test_fusion_pieces() -> None:
    ranks = percentile_ranks(np.array([0.2, np.nan, 0.5, 0.5, 0.9]))
    assert np.allclose(ranks[[0, 2, 3, 4]], [0.0, 0.5, 0.5, 1.0]) and np.isnan(ranks[1])
    assert unit(sm.Method(kind=sm.COLOR), np.array([1.6]))[0] == pytest.approx(0.8)
    assert unit(sm.Method(kind=sm.FIELD), np.array([-0.3]))[0] == 0.0
    methods = [sm.Method(), sm.Method(kind=sm.SCAFFOLD)]
    flat = consensus(methods, [np.array([0.2, 0.6]), np.array([1.0, 1.0])], "normalized")
    assert np.allclose(flat, [0.0, 1.0])  # a method that gives every row the same is left out


def test_cancel_keeps_the_finished_rows(monkeypatch) -> None:
    from chemsuite.analysis import fingerprints as module

    monkeypatch.setattr(module, "TILE_ROWS", 3)
    token = CancelToken()
    seen = []

    def progress(done: int, total: int) -> None:
        seen.append(done)
        if len(seen) == 1:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(),)), progress=progress,
                    cancel=token)  # fmt: skip
    partial = stopped.value.partial
    assert len(partial) == 3 and partial.columns["ECFP4 Tanimoto"][0] == pytest.approx(1.0)


def test_unbuilt_methods_say_so() -> None:
    with pytest.raises(ValueError, match="not available"):
        score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(kind=sm.SHAPE),)))
    with pytest.raises(ValueError, match="reference"):
        score_cells(CELLS, [], sm.Search((sm.Method(),)))


def test_mcs_tanimoto() -> None:
    from chemsuite.similarity.substructure import mcs_tanimoto

    scores = score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(kind=sm.MCS),), prefilter=0))
    values = scores.columns["MCS Tanimoto"]
    assert values[0] == pytest.approx(1.0)  # itself
    paracetamol = Chem.MolFromSmiles(CELLS[0])
    expected = mcs_tanimoto(Chem.MolFromSmiles(CELLS[1]), paracetamol, 2)
    assert values[1] == pytest.approx(expected) and 0.5 < values[1] < 1.0
    assert values[1] > values[4]  # phenacetin shares more than nicotine does
    assert np.isnan(values[6:]).all() and scores.nearest[:6].tolist() == [0] * 6


def test_scaffold_match() -> None:
    cells = ["c1ccc2ccccc2c1", "Cc1ccc2ccccc2c1", "c1ccccc1", "CCO", None]
    naphthalene = sm.Reference(Chem.MolFromSmiles(cells[0]), "naphthalene")
    for kind in (sm.SCAFFOLD, sm.GENERIC_SCAFFOLD):
        method = sm.Method(kind=kind)
        values = score_cells(cells, [naphthalene], sm.Search((method,))).columns[method.column]
        assert values[:4].tolist() == [1.0, 1.0, 0.0, 0.0], kind  # the methyl is a side chain
        assert np.isnan(values[4]), kind
    pyridine = sm.Reference(Chem.MolFromSmiles("c1ccncc1"), "pyridine")
    exact, generic = (sm.Method(kind=k) for k in (sm.SCAFFOLD, sm.GENERIC_SCAFFOLD))
    rings = ["c1ccccc1", "c1ccncc1"]
    strict = score_cells(rings, [pyridine], sm.Search((exact,))).columns[exact.column]
    assert strict.tolist() == [0.0, 1.0]
    loose = score_cells(rings, [pyridine], sm.Search((generic,))).columns[generic.column]
    assert loose.tolist() == [1.0, 1.0]  # generic: every atom a carbon, every bond single


def test_the_prefilter_limits_the_expensive_method() -> None:
    method = sm.Method(kind=sm.MCS)
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(), method), prefilter=2))
    mcs = scores.columns["MCS Tanimoto"]
    ecfp = scores.columns["ECFP4 Tanimoto"]
    assert np.isfinite(mcs).sum() == 2  # only the two the fingerprint ranked best
    best = np.argsort(-np.nan_to_num(ecfp, nan=-np.inf))[:2]
    assert np.isfinite(mcs[best]).all()
    assert np.isfinite(ecfp).sum() == 6  # the cheap method still scores every row
    whole = score_cells(CELLS, [PARACETAMOL], sm.Search((method,), prefilter=0))
    assert np.isfinite(whole.columns["MCS Tanimoto"]).sum() == 6
    assert np.isfinite(mcs).sum() < np.isfinite(whole.columns["MCS Tanimoto"]).sum()
