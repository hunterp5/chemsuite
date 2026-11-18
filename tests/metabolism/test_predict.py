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

"""The batch driver: one row in, one entry out, blanks and failures included."""

from __future__ import annotations

import math

import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_to_blob
from chemsuite.metabolism import predict
from chemsuite.metabolism.predict import predict_metabolism

SMILES = ("CCO", "CC(=O)Nc1ccc(O)cc1", "Cc1ccccc1")


def blobs(*smiles: str) -> list[bytes]:
    return [mol_to_blob(Chem.MolFromSmiles(text)) for text in smiles]


def test_every_row_gets_an_entry_for_every_isoform_asked_for():
    results = predict_metabolism(blobs(*SMILES), ("3A4", "2D6"))
    assert len(results) == 3
    assert set(results.top_score) == {"3A4", "2D6"}
    assert results.failed == 0
    for isoform in ("3A4", "2D6"):
        assert len(results.top_site[isoform]) == 3
        assert all(0.0 <= score <= 1.0 for score in results.top_score[isoform])
        # the top site is named by element and index, and is the first ranked atom
        for row, text in enumerate(results.top_site[isoform]):
            assert text[0].isalpha() and text[1:].isdigit()
            assert int(text[1:]) == results.ranked_atoms[isoform][row][0]


def test_the_sites_text_lists_the_best_three_with_their_scores():
    results = predict_metabolism(blobs("CC(=O)Nc1ccc(O)cc1"), ("3A4",))
    text = results.sites["3A4"][0]
    assert text.count(",") == predict.RANKED - 1
    assert text.startswith(results.top_site["3A4"][0] + " ")
    assert text.split(",")[0].split()[1] == f"{results.top_score['3A4'][0]:.2f}"


def test_a_blank_cell_is_blank_and_not_a_failure():
    results = predict_metabolism([None, b"", *blobs("CCO")], ("3A4",))
    assert len(results) == 3
    assert results.failed == 0
    assert results.top_site["3A4"][:2] == ["", ""]
    assert math.isnan(results.top_score["3A4"][0])
    assert results.ranked_atoms["3A4"][0] == []
    assert results.top_site["3A4"][2] != ""


def test_an_unparsable_cell_is_counted_as_a_failure():
    results = predict_metabolism(["not a molecule", *blobs("CCO")], ("3A4",))
    assert results.failed == 1
    assert len(results) == 2
    assert results.top_site["3A4"][0] == ""


def test_smiles_text_is_parsed_like_a_structure():
    from_text = predict_metabolism(["CCO"], ("3A4",))
    from_blob = predict_metabolism(blobs("CCO"), ("3A4",))
    assert from_text.top_site["3A4"] == from_blob.top_site["3A4"]
    assert from_text.top_score["3A4"] == from_blob.top_score["3A4"]


def test_progress_is_reported_per_chunk_and_reaches_the_total():
    seen: list[tuple[int, int]] = []
    rows = blobs(*SMILES) * 20  # 60 rows, more than one chunk
    predict_metabolism(rows, ("3A4",), progress=lambda done, total: seen.append((done, total)))
    assert seen[-1] == (len(rows), len(rows))
    assert [done for done, _ in seen] == sorted(done for done, _ in seen)
    assert len(seen) == math.ceil(len(rows) / predict.CHUNK_ROWS)


def test_a_cancel_keeps_the_rows_already_finished():
    cancel = CancelToken()
    rows = blobs(*SMILES) * 20

    def stop(done: int, _total: int) -> None:
        if done >= predict.CHUNK_ROWS:
            cancel.cancel()

    with pytest.raises(Cancelled) as raised:
        predict_metabolism(rows, ("3A4",), progress=stop, cancel=cancel)
    partial = raised.value.partial
    assert partial is not None
    assert len(partial) == predict.CHUNK_ROWS


def test_a_cancel_before_the_first_chunk_keeps_nothing():
    cancel = CancelToken()
    cancel.cancel()
    with pytest.raises(Cancelled) as raised:
        predict_metabolism(blobs("CCO"), ("3A4",), cancel=cancel)
    assert raised.value.partial is None
