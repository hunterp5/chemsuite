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


def _cancel_in_the_mcs(cells, references, search, monkeypatch, after: int = 2):
    """Run *search* and cancel part way through its MCS stage; the Cancelled it raises."""
    from chemsuite.similarity import substructure

    monkeypatch.setitem(substructure.CHUNK_ROWS, sm.MCS, 1)
    token = CancelToken()
    calls = []

    def progress(done: int, total: int) -> None:
        calls.append(done)
        if len(calls) == after:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        score_cells(cells, references, search, progress=progress, cancel=token)
    return stopped.value


def test_a_cancelled_expensive_method_keeps_its_scores_on_their_own_rows(monkeypatch) -> None:
    # the two rows the fingerprint ranks best are the last two, not the first
    cells = ["c1ccccc1", "CCO", "CCCCN", "CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1O"]
    aspirin = sm.Reference(Chem.MolFromSmiles(cells[3]), "aspirin")
    search = sm.Search((sm.Method(), sm.Method(kind=sm.MCS)), prefilter=2)
    whole = score_cells(cells, [aspirin], search)
    assert np.flatnonzero(np.isfinite(whole.columns["MCS Tanimoto"])).tolist() == [3, 4]
    partial = _cancel_in_the_mcs(cells, [aspirin], search, monkeypatch).partial
    assert len(partial) == len(cells)  # every row of the finished cheap method is kept
    assert np.isfinite(partial.columns["ECFP4 Tanimoto"]).all()
    mcs = partial.columns["MCS Tanimoto"]
    scored = np.flatnonzero(np.isfinite(mcs))
    assert 0 < len(scored) < 2 and scored.tolist() == [3]  # the best row, where it belongs
    assert mcs[3] == pytest.approx(whole.columns["MCS Tanimoto"][3])


def test_a_cancelled_expensive_method_keeps_every_row_without_a_prefilter(monkeypatch) -> None:
    search = sm.Search((sm.Method(), sm.Method(kind=sm.MCS)))  # the default prefilter: 500 > 8
    whole = score_cells(CELLS, [PARACETAMOL], search)
    partial = _cancel_in_the_mcs(CELLS, [PARACETAMOL], search, monkeypatch).partial
    assert len(partial) == len(CELLS)  # the finished fingerprint rows are all kept
    ecfp = partial.columns["ECFP4 Tanimoto"]
    assert np.isfinite(ecfp).sum() == 6
    assert np.allclose(ecfp, whole.columns["ECFP4 Tanimoto"], equal_nan=True)
    mcs = partial.columns["MCS Tanimoto"]
    scored = np.flatnonzero(np.isfinite(mcs))
    assert 0 < len(scored) < 6  # blank where it did not reach, right where it did
    assert np.allclose(mcs[scored], whole.columns["MCS Tanimoto"][scored])


def test_a_cancelled_cheap_method_keeps_its_first_rows(monkeypatch) -> None:
    from chemsuite.analysis import fingerprints as module

    monkeypatch.setattr(module, "TILE_ROWS", 3)
    token = CancelToken()
    calls = []

    def progress(done: int, total: int) -> None:
        calls.append(done)
        if len(calls) == 1:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        score_cells(CELLS, [PARACETAMOL], sm.Search((sm.Method(), sm.Method(kind=sm.MCS)),),
                    progress=progress, cancel=token)  # fmt: skip
    partial = stopped.value.partial
    assert len(partial) == 3 and "MCS Tanimoto" not in partial.columns
    assert partial.consensus is None  # the other method never ran


def test_the_consensus_needs_every_method_to_have_scored_the_row() -> None:
    methods = (sm.Method(), sm.Method(kind=sm.MCS))
    scores = score_cells(CELLS, [PARACETAMOL], sm.Search(methods, prefilter=2, consensus="rank"))
    mcs = scores.columns["MCS Tanimoto"]
    assert np.isfinite(scores.consensus).sum() == 2  # not the rows the MCS skipped
    assert np.isfinite(scores.consensus).tolist() == np.isfinite(mcs).tolist()
    top = np.nanargmax(scores.consensus)  # no skipped row outranks a scored one
    assert np.isfinite(mcs[top])


def test_nearest_falls_back_when_the_primary_method_skipped_the_row() -> None:
    methods = (sm.Method(kind=sm.MCS), sm.Method())  # the expensive one is primary
    scores = score_cells(CELLS, [PARACETAMOL, ASPIRIN], sm.Search(methods, prefilter=2))
    mcs = scores.columns["MCS Tanimoto"]
    assert np.isfinite(mcs).sum() == 2
    assert (scores.nearest[:6] >= 0).all()  # every scored row still names a reference
    skipped = np.flatnonzero(~np.isfinite(mcs[:6]))
    fingerprint_nearest = score_cells(CELLS, [PARACETAMOL, ASPIRIN], sm.Search((sm.Method(),)))
    assert scores.nearest[skipped].tolist() == fingerprint_nearest.nearest[skipped].tolist()


def test_the_methods_run_in_worker_processes() -> None:
    from chemsuite.core.jobs.parallel import process_pool

    methods = (sm.Method(), sm.Method(kind=sm.MCS), sm.Method(kind=sm.SCAFFOLD))
    search = sm.Search(methods, prefilter=0)
    inline = score_cells(CELLS, [PARACETAMOL], search)
    with process_pool(2) as pool:
        pooled = score_cells(CELLS, [PARACETAMOL], search, executor=pool)
    for method in methods:  # the worker entry points pickle and give the same answers
        got, want = pooled.columns[method.column], inline.columns[method.column]
        assert np.allclose(got, want, equal_nan=True), method.column


def _in_3d(smiles: str, name: str) -> sm.Reference:
    from chemsuite.conformers.etkdg import generate

    return sm.Reference(generate(Chem.MolFromSmiles(smiles)).mol, name)


def test_a_flat_reference_is_turned_away_by_the_3d_methods() -> None:
    search = sm.Search((sm.Method(kind=sm.SHAPE),))
    with pytest.raises(ValueError, match="no 3D coordinates"):
        score_cells(CELLS, [PARACETAMOL], search)  # from SMILES: flat
    with pytest.raises(ValueError, match="paracetamol"):  # the message names it
        score_cells(CELLS, [PARACETAMOL], search)
    sm.check_references([PARACETAMOL], sm.Search((sm.Method(),)))  # 2D methods do not mind
    sm.check_references([_in_3d(CELLS[0], "paracetamol")], search)


def test_no_reference_is_turned_away() -> None:
    with pytest.raises(ValueError, match="reference"):
        score_cells(CELLS, [], sm.Search((sm.Method(),)))


@pytest.mark.parametrize("kind", [sm.SHAPE, sm.COLOR, sm.FIELD])
def test_3d_methods_score_a_molecule_against_itself_highest(kind) -> None:
    cells = [CELLS[0], CELLS[1], "C1C2CC3CC1CC(C2)C3"]  # paracetamol, phenacetin, adamantane
    reference = _in_3d(CELLS[0], "paracetamol")
    method = sm.Method(kind=kind, conformers=3)
    scores = score_cells(cells, [reference], sm.Search((method,), prefilter=0))
    values = scores.columns[method.column]
    assert np.isfinite(values[:3]).all()
    assert values[0] == values[:3].max()  # itself fits best
    assert values[0] > values[2]  # and adamantane fits worst
    assert scores.nearest[:3].tolist() == [0, 0, 0]


def test_3d_keeps_the_reference_each_row_fits_best() -> None:
    # rigid molecules, so the fit does not turn on which conformers were generated
    cells = ["CC(=O)Nc1ccc(O)cc1", "C1C2CC3CC1CC(C2)C3"]  # paracetamol, adamantane
    refs = [_in_3d(cells[0], "paracetamol"), _in_3d(cells[1], "adamantane")]
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    scores = score_cells(cells, refs, sm.Search((method,), prefilter=0))
    assert scores.nearest[:2].tolist() == [0, 1]  # each row fits its own shape best
    values = scores.columns[method.column]
    assert (values[:2] > 0.5).all()


def test_3d_keeps_each_row_s_pose_on_its_nearest_reference() -> None:
    from chemsuite.conformers import energies
    from chemsuite.core.mols import mol_from_blob
    from chemsuite.superpose.model import reference_id_of

    cells = ["CC(=O)Nc1ccc(O)cc1", "C1C2CC3CC1CC(C2)C3", None]  # paracetamol, adamantane
    refs = [_in_3d(cells[0], "paracetamol"), _in_3d(cells[1], "adamantane")]
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    search = sm.Search((sm.Method(), method), prefilter=0)
    assert score_cells(cells, refs, search).poses is None  # only when asked for
    scores = score_cells(cells, refs, search, pose_ids=["ref-a", "ref-b"])
    assert scores.poses[2] is None  # a blank row has no pose
    for row in (0, 1):
        blob, smiles, index, delta_e = scores.poses[row]
        mol = mol_from_blob(blob)
        assert mol.GetNumConformers() == 1 and smiles == Chem.CanonSmiles(cells[row])
        # the row's own conformers were generated, so the pose knows its energy above the lowest
        assert delta_e >= 0 and delta_e == pytest.approx(energies.delta_of(mol))
        assert index == row  # each fits its own shape best, as the 3D score says
        assert reference_id_of(mol) == ["ref-a", "ref-b"][row]
    with pytest.raises(ValueError, match="one pose id per reference"):
        score_cells(cells, refs, search, pose_ids=["ref-a"])


def test_3d_poses_follow_the_prefilter() -> None:
    cells = ["CCCCCCCCCC", "CC(=O)Nc1ccc(O)cc1", "CCO"]
    reference = _in_3d(cells[1], "paracetamol")
    search = sm.Search((sm.Method(), sm.Method(kind=sm.SHAPE, conformers=3)), prefilter=1)
    scores = score_cells(cells, [reference], search, pose_ids=["ref"])
    assert [p is not None for p in scores.poses] == [False, True, False]


def test_3d_fuses_over_the_references_as_asked() -> None:
    cells = ["CC(=O)Nc1ccc(O)cc1", "C1C2CC3CC1CC(C2)C3"]  # paracetamol, adamantane
    refs = [_in_3d(cells[0], "paracetamol"), _in_3d(cells[1], "adamantane")]
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    best = score_cells(cells, refs, sm.Search((method,), prefilter=0, fuse="max"))
    mean = score_cells(cells, refs, sm.Search((method,), prefilter=0, fuse="mean"))
    each = [
        score_cells(cells, [r], sm.Search((method,), prefilter=0)).columns[method.column]
        for r in refs
    ]
    pairs = np.column_stack(each)
    assert np.allclose(best.columns[method.column], pairs.max(axis=1))
    assert np.allclose(mean.columns[method.column], pairs.mean(axis=1))
    assert (mean.columns[method.column] < best.columns[method.column]).all()


def test_3d_uses_a_row_s_own_stored_conformers() -> None:
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.superpose.model import Settings
    from chemsuite.superpose.run import prepare_reference, superpose_one

    reference = _in_3d("CC(=O)Nc1ccc(O)cc1", "paracetamol")
    stored = _in_3d("CC(=O)Nc1ccc(OC)cc1", "phenacetin").mol  # a row that already has 3D
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    scores = score_cells([mol_to_blob(stored)], [reference], sm.Search((method,), prefilter=0))
    settings = Settings(method="shape", conformers=3)
    _moved, fit = superpose_one(stored, prepare_reference(reference.mol), settings)
    # the stored conformers are used as they are, so the score is exactly Superpose's
    assert scores.columns[method.column][0] == pytest.approx(fit.scores["Shape Tanimoto"])


def test_a_3d_reference_that_is_a_row_skips_itself() -> None:
    cells = ["CC(=O)Nc1ccc(O)cc1", "CCCCCCCCCC"]
    refs = [sm.Reference(_in_3d(cells[0], "").mol, "paracetamol", row=0)]
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    scores = score_cells(cells, refs, sm.Search((method,), prefilter=0))
    values = scores.columns[method.column]
    assert np.isnan(values[0]) and scores.nearest[0] == -1  # its only reference was itself
    assert np.isfinite(values[1])


def test_3d_runs_in_worker_processes_and_the_prefilter_applies() -> None:
    from chemsuite.core.jobs.parallel import process_pool

    cells = [CELLS[0], CELLS[1], "CCCCCCCCCC", "CCO"]
    reference = _in_3d(CELLS[0], "paracetamol")
    method = sm.Method(kind=sm.SHAPE, conformers=3)
    search = sm.Search((sm.Method(), method), prefilter=2)
    with process_pool(2) as pool:
        scores = score_cells(cells, [reference], search, executor=pool)
    shape = scores.columns[method.column]
    assert np.isfinite(shape).sum() == 2  # only the fingerprint's two best rows
    assert np.isfinite(shape[0])  # paracetamol itself among them
    assert np.isfinite(scores.columns["ECFP4 Tanimoto"]).sum() == 4


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
