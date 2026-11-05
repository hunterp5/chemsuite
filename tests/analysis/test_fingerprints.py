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
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

from chemsuite.analysis import clustering as cl
from chemsuite.analysis.diversity import maxmin_fingerprints, maxmin_numeric
from chemsuite.analysis.fingerprints import (
    KINDS,
    cache_for,
    compute_fingerprints,
    max_similarity,
    neighbours,
    query_fingerprint,
    similarity,
)
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_to_blob
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn

SMILES = [
    "CCO", "CCCO", "CCCCO", "CCCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "c1ccncc1",
    "Cc1ccncc1", "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "O=C(O)c1ccccc1O", "NCCc1ccc(O)c(O)c1",
    "CN1CCC[C@H]1c1cccnc1", "C1CCCCC1", "C1CCNCC1", "OC1CCCCC1", "Clc1ccccc1",
]  # fmt: skip
CELLS = [*SMILES, None, "not a smiles"]


@pytest.fixture(scope="module")
def fps():
    return compute_fingerprints(CELLS, "morgan2")


def _rdkit(smiles):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]


def test_fingerprints_match_rdkit(fps) -> None:
    assert fps.valid.tolist() == [True] * 20 + [False, False]
    ref = _rdkit(SMILES)
    sims = similarity(fps, fps.take(np.array([9])))
    expected = DataStructs.BulkTanimotoSimilarity(ref[9], ref)
    assert np.allclose(sims[:20], expected) and np.isnan(sims[20:]).all()
    dice = similarity(fps, fps.take(np.array([9])), "dice")
    assert np.allclose(dice[:20], DataStructs.BulkDiceSimilarity(ref[9], ref))
    cosine = similarity(fps, fps.take(np.array([9])), "cosine")
    assert np.allclose(cosine[:20], DataStructs.BulkCosineSimilarity(ref[9], ref))


def test_every_kind_and_blobs_and_3d() -> None:
    blob = mol_to_blob(Chem.MolFromSmiles("CCOc1ccccc1"))
    for key in KINDS:
        out = compute_fingerprints([blob, "CCOc1ccccc1"], key)
        assert out.valid.all() and (out.packed[0] == out.packed[1]).all(), key
        assert out.counts[0] > 0 and out.packed.shape[1] % 8 == 0, key
    maccs = compute_fingerprints(["CCOc1ccccc1"], "maccs")
    assert (
        maccs.counts[0] == MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles("CCOc1ccccc1")).GetNumOnBits()
    )
    from chemsuite.conformers.etkdg import generate

    three_d = mol_to_blob(generate(Chem.MolFromSmiles("CCOc1ccccc1")).mol)
    out = compute_fingerprints([three_d, "CCOc1ccccc1"])
    assert (out.packed[0] == out.packed[1]).all()  # explicit hydrogens do not matter


def test_query_and_parallel(fps) -> None:
    query = query_fingerprint(Chem.MolFromSmiles("CCCO"))
    assert similarity(fps, query)[1] == pytest.approx(1.0)
    from chemsuite.core.jobs.parallel import process_pool

    with process_pool(2) as pool:
        parallel = compute_fingerprints(CELLS, "morgan2", executor=pool, chunk_rows=5)
    assert (parallel.packed == fps.packed).all()


def test_max_similarity_and_neighbours(fps) -> None:
    refs = np.array([0, 4])
    best, which = max_similarity(
        fps, fps.take(refs), exclude_self=np.array([0, -1, -1, -1, 1] + [-1] * 17)
    )
    singles = np.column_stack([similarity(fps, fps.take(np.array([r]))) for r in refs])
    singles[0, 0] = singles[4, 1] = -1  # each ref ignores itself
    assert np.allclose(best[:20], singles[:20].max(axis=1))
    assert which[1] == 0 and which[5] == 1 and which[20] == -1
    indptr, indices = neighbours(fps, 0.3)
    full = np.column_stack([similarity(fps, fps.take(np.array([i]))) for i in range(20)])
    for i in range(20):
        expected = sorted(j for j in range(20) if j != i and full[i, j] >= 0.3)
        assert indices[indptr[i] : indptr[i + 1]].tolist() == expected
    assert indptr[21] == indptr[22] == indptr[20]  # invalid rows have none


def test_tiles_give_the_same_answer(fps, monkeypatch) -> None:
    from chemsuite.analysis import fingerprints as module

    whole = neighbours(fps, 0.2)
    monkeypatch.setattr(module, "TILE_ROWS", 3)
    tiled = neighbours(fps, 0.2)
    assert all((a == b).all() for a, b in zip(whole, tiled, strict=True))
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        neighbours(fps, 0.2, cancel=token)


def test_butina_matches_rdkit(fps) -> None:
    clusters = cl.butina(fps, 0.4)
    ref = _rdkit(SMILES)
    distances = []
    for i in range(1, 20):
        distances.extend(DataStructs.BulkTanimotoSimilarity(ref[i], ref[:i], returnDistance=True))
    expected = Butina.ClusterData(distances, 20, 0.6, isDistData=True)
    assert clusters.count == len(expected)
    assert sorted(np.bincount(clusters.labels[:20]).tolist()) == sorted(len(c) for c in expected)
    assert clusters.labels[20] == -1 and clusters.centres.sum() == clusters.count
    sizes = np.bincount(clusters.labels[clusters.labels >= 0])
    assert (np.diff(sizes) <= 0).all()  # numbered largest first


def test_kmeans_and_hierarchical(fps) -> None:
    km = cl.kmeans_fingerprints(fps, 4, seed=1)
    assert km.count == 4 and km.centres.sum() == 4 and (km.labels[20:] == -1).all()
    hier = cl.hierarchical_fingerprints(fps, 3, "complete")
    assert hier.count == 3 and hier.centres.sum() == 3
    with pytest.raises(ValueError, match="Ward"):
        cl.hierarchical_fingerprints(fps, 3, "ward")
    with pytest.raises(ValueError):
        cl.kmeans_fingerprints(fps, 50)


def test_numeric_clustering() -> None:
    rng = np.random.default_rng(0)
    a = np.r_[rng.normal(0, 0.1, 30), rng.normal(5, 0.1, 20), np.nan]
    b = np.r_[rng.normal(0, 0.1, 30), rng.normal(5, 0.1, 20), 1.0]
    for clusters in (cl.kmeans_numeric([a, b], 2), cl.hierarchical_numeric([a, b], 2, "ward")):
        assert (
            clusters.labels[:30].tolist() == [0] * 30
            and clusters.labels[30:50].tolist() == [1] * 20
        )
        assert clusters.labels[50] == -1 and clusters.sizes()[0] == 30
    with pytest.raises(ValueError, match="do not vary"):
        cl.kmeans_numeric([np.ones(10)], 2)
    with pytest.raises(ValueError, match="limited"):
        cl.hierarchical_numeric([np.arange(cl.MAX_HIERARCHICAL + 1.0)], 2)


def test_maxmin(fps) -> None:
    picks = maxmin_fingerprints(fps, 5, seed=3)
    assert len(picks) == 5 and len(set(picks.rows.tolist())) == 5
    assert not set(picks.rows.tolist()) & {20, 21}
    assert (np.diff(picks.distances[1:]) <= 1e-12).all()  # each pick no farther than the last
    seeded = maxmin_fingerprints(fps, 3, seeds=[0, 4])
    assert not set(seeded.rows.tolist()) & {0, 4}
    dupes = compute_fingerprints(["CCO", "CCO", "CCO"])
    assert len(maxmin_fingerprints(dupes, 3)) == 1  # duplicates add nothing
    x = np.array([0.0, 0.1, 10.0, 5.0, np.nan])
    numeric = maxmin_numeric([x], 3, seeds=[0])
    assert numeric.rows.tolist() == [2, 3, 1]


def test_cache_follows_the_table() -> None:
    table = Table()
    table.add_rows({"S": StructureColumn.from_texts(["CCO", "CCN", "CCC"])}, undoable=False)
    cache = cache_for(table)
    assert cache_for(table) is cache
    oids = table.oids.copy()
    fps = compute_fingerprints(["CCO", "CCN", "CCC"])
    cache.put("S", "morgan2", oids, fps, cache.generation("S"))
    got = cache.get("S", "morgan2", oids[::-1])
    assert (got.packed == fps.packed[::-1]).all()
    assert cache.get("S", "morgan2", np.array([999])) is None
    stale = cache.generation("S")
    table.set_cells("S", oids[:1], ["c1ccccc1"])
    assert cache.get("S", "morgan2", oids) is None
    cache.put("S", "morgan2", oids, fps, stale)  # computed before the edit: ignored
    assert cache.get("S", "morgan2", oids) is None
