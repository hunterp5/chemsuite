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

from chemsuite.core.jobs import CancelToken
from chemsuite.similarity import model as sm
from chemsuite.similarity.run import score_cells
from chemsuite.similarity.screen import screen_file

PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
#: a small library: two close to paracetamol, the rest not, and one line that is not a structure
LIBRARY = [
    ("CCO", "ethanol"),
    (PARACETAMOL, "paracetamol"),
    ("c1ccccc1", "benzene"),
    ("CC(=O)Nc1ccc(OC)cc1", "phenacetin"),
    ("CCCCCCCC", "octane"),
    ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
    ("not a structure", "junk"),
]


@pytest.fixture
def library(tmp_path):
    path = tmp_path / "library.smi"
    path.write_text("\n".join(f"{s} {n}" for s, n in LIBRARY) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def sdf_library(tmp_path):
    path = tmp_path / "library.sdf"
    writer = Chem.SDWriter(str(path))
    for smiles, name in LIBRARY[:-1]:  # an SD file holds no unparsable record
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", name)
        mol.SetProp("Vendor", f"V-{name}")
        writer.write(mol)
    writer.close()
    return path


def _reference(smiles: str = PARACETAMOL, name: str = "paracetamol") -> sm.Reference:
    return sm.Reference(Chem.MolFromSmiles(smiles), name)


def _column(table, name: str) -> list[str]:
    return table.column(name).displays(np.arange(len(table)))


def test_the_hits_come_back_best_first_with_their_scores(library) -> None:
    search = sm.Search((sm.Method(),))
    hits = screen_file(library, [_reference()], search, top_k=3)
    assert hits.scanned == 7 and hits.readable == 6  # the junk line is counted but unreadable
    assert len(hits) == 3 and not hits.cancelled
    assert _column(hits.table, "Name")[:2] == ["paracetamol", "phenacetin"]
    scores = [float(s) for s in _column(hits.table, "ECFP4 Tanimoto")]
    assert scores[0] == pytest.approx(1.0)
    assert scores == sorted(scores, reverse=True)
    assert "Nearest Reference" not in hits.table.column_names  # only one reference to name


def test_top_k_keeps_exactly_the_best_rows(library) -> None:
    search = sm.Search((sm.Method(),))
    every = screen_file(library, [_reference()], search, top_k=0)
    assert len(every) == 6  # every readable row
    ranked = _column(every.table, "Name")
    for k in (1, 2, 4):
        hits = screen_file(library, [_reference()], search, top_k=k)
        assert _column(hits.table, "Name") == ranked[:k], k


def test_the_threshold_drops_the_rest(library) -> None:
    search = sm.Search((sm.Method(),))
    hits = screen_file(library, [_reference()], search, top_k=0, threshold=0.5)
    assert _column(hits.table, "Name") == ["paracetamol", "phenacetin"]
    assert all(float(s) >= 0.5 for s in _column(hits.table, "ECFP4 Tanimoto"))
    assert not len(screen_file(library, [_reference()], search, threshold=1.01).table)


def test_the_file_s_own_columns_come_too(sdf_library) -> None:
    hits = screen_file(sdf_library, [_reference()], sm.Search((sm.Method(),)), top_k=2)
    assert hits.scanned == 6 and len(hits) == 2
    assert "Vendor" in hits.table.column_names
    assert _column(hits.table, "Vendor")[0] == "V-paracetamol"
    assert hits.table.kind("Structure") == "structure"


def test_every_method_scores_the_survivors(sdf_library) -> None:
    methods = (sm.Method(), sm.Method(kind=sm.MCS), sm.Method(kind=sm.SCAFFOLD))
    hits = screen_file(sdf_library, [_reference()], sm.Search(methods, prefilter=0), top_k=3)
    for method in methods:
        assert method.column in hits.table.column_names, method.column
    assert "Consensus" in hits.table.column_names
    # the same rows scored through the table search give the same numbers
    cells = [s for s, _n in LIBRARY if s != "not a structure"]
    table = score_cells(cells, [_reference()], sm.Search(methods, prefilter=0))
    best = _column(hits.table, "MCS Tanimoto")[0]
    assert float(best) == pytest.approx(table.columns["MCS Tanimoto"].max(), abs=5e-4)


def test_several_references_name_the_nearest(library) -> None:
    refs = [_reference(), _reference("CCCCCCCC", "octane")]
    hits = screen_file(library, refs, sm.Search((sm.Method(),)), top_k=0)
    names, nearest = _column(hits.table, "Name"), _column(hits.table, "Nearest Reference")
    named = dict(zip(names, nearest, strict=True))
    assert named["paracetamol"] == "paracetamol" and named["octane"] == "octane"


def test_cancelling_keeps_the_hits_found_so_far(library, monkeypatch) -> None:
    from chemsuite.similarity import screen as module

    # the readers' batches hold 20,000 rows, so cut the file into batches of two to cancel
    # part way through one of a size a test can write
    real = module.read_batches

    def in_twos(path, **kwargs):
        for batch in real(path, **kwargs):
            rows = max(len(values) for values in batch.values())
            for start in range(0, rows, 2):
                yield module.slice_batch(batch, np.arange(start, min(rows, start + 2)))

    monkeypatch.setattr(module, "read_batches", in_twos)
    token = CancelToken()
    batches = []

    def progress(done: int, total: int) -> None:
        batches.append(done)

    hits = screen_file(library, [_reference()], sm.Search((sm.Method(),)), top_k=0,
                       progress=progress, cancel=token)  # fmt: skip
    assert not hits.cancelled and hits.scanned == 7  # nothing cancelled: the whole file
    token.cancel()
    stopped = screen_file(library, [_reference()], sm.Search((sm.Method(),)), top_k=0,
                          cancel=token)  # fmt: skip
    assert stopped.cancelled and stopped.scanned == 2  # it stopped after the first batch
    assert len(stopped) == stopped.readable == 2  # and those rows are scored hits
    assert "ECFP4 Tanimoto" in stopped.table.column_names
    assert _column(stopped.table, "Name") == ["paracetamol", "ethanol"]  # best of the two


def test_a_token_cancelled_before_the_start_stops_at_the_first_batch(library) -> None:
    # a batch is scored whole, so an already-cancelled screen still returns that first batch's
    # hits (this small file is one batch) rather than nothing
    token = CancelToken()
    token.cancel()
    hits = screen_file(library, [_reference()], sm.Search((sm.Method(),)), cancel=token)
    assert hits.cancelled and len(hits) == hits.readable
    assert _column(hits.table, "Name")[0] == "paracetamol"


def test_a_file_without_structures_says_so(tmp_path) -> None:
    path = tmp_path / "numbers.csv"
    path.write_text("ID,Value\na,1\nb,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no structure column"):
        screen_file(path, [_reference()], sm.Search((sm.Method(),)))


def test_a_flat_reference_is_turned_away_before_the_file_is_read(library) -> None:
    search = sm.Search((sm.Method(kind=sm.SHAPE),))
    with pytest.raises(ValueError, match="no 3D coordinates"):
        screen_file(library, [_reference()], search)
