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

"""BRICS and RECAP fragmentation."""

from __future__ import annotations

import math

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.sar.fragments import (
    fragment,
    fragment_rows,
    fragment_stats,
    plain_dummies,
)

PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
ACETANILIDE = "CC(=O)Nc1ccccc1"


def cut(smiles: str, method: str = "brics", **kwargs) -> list[str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return fragment(mol, method, **kwargs)


def test_brics_cuts_the_amide_and_labels_the_ends() -> None:
    # the amide C-N bond is BRICS bond type 1-16
    assert cut(PARACETAMOL) == ["[16*]c1ccc(O)cc1", "[1*]C(C)=O"]


def test_recap_cuts_the_amide_and_writes_plain_attachment_points() -> None:
    found = cut(PARACETAMOL, "recap")
    assert sorted(found) == ["*C(C)=O", "*Nc1ccc(O)cc1"]
    assert all("*" in piece for piece in found)


def test_fragments_come_back_largest_first() -> None:
    found = cut(PARACETAMOL)
    sizes = [Chem.MolFromSmiles(s, sanitize=False).GetNumAtoms() for s in found]
    assert sizes == sorted(sizes, reverse=True)


def test_plain_drops_the_brics_labels() -> None:
    assert cut(PARACETAMOL, plain=True) == ["[*]c1ccc(O)cc1", "[*]C(C)=O"]
    assert plain_dummies("[16*]c1ccccc1.[1*]C") == "[*]c1ccccc1.[*]C"


def test_min_size_leaves_out_the_small_pieces() -> None:
    assert "[1*]C(C)=O" in cut(PARACETAMOL, min_size=3)
    assert "[1*]C(C)=O" not in cut(PARACETAMOL, min_size=4)


def test_a_molecule_with_no_cuttable_bond_is_its_own_recap_fragment() -> None:
    assert cut("c1ccccc1", "recap") == ["c1ccccc1"]


def test_a_molecule_too_small_to_keep_gives_nothing() -> None:
    assert cut("CC", "recap") == []


def test_an_unknown_method_is_refused() -> None:
    with pytest.raises(ValueError, match="method must be one of"):
        cut(PARACETAMOL, "murcko")


def test_fragment_rows_keeps_row_order_and_counts_failures() -> None:
    found = fragment_rows([PARACETAMOL, ACETANILIDE, None, "not a structure"])
    assert len(found) == 4
    assert found.rows[0] == ["[16*]c1ccc(O)cc1", "[1*]C(C)=O"]
    assert found.rows[2] == [] and found.rows[3] == []
    assert found.failed == 1  # the blank row is not a failure
    assert found.widest() == 2
    assert found.fragmented() == 2


def test_fragment_rows_reports_progress() -> None:
    seen: list[tuple[int, int]] = []
    fragment_rows([PARACETAMOL] * 25, progress=lambda done, total: seen.append((done, total)))
    assert seen and seen[-1] == (25, 25)


def test_cancelling_keeps_the_rows_already_done() -> None:
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done >= 20:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        fragment_rows([PARACETAMOL] * 100, progress=progress, cancel=token)
    partial = stopped.value.partial
    assert partial is not None and 0 < len(partial) < 100
    assert partial.rows[0] == ["[16*]c1ccc(O)cc1", "[1*]C(C)=O"]


def test_fragment_stats_counts_compounds_and_averages_activity() -> None:
    found = fragment_rows([PARACETAMOL, ACETANILIDE])
    stats = fragment_stats(found, np.array([5.0, 7.0]))
    shared = next(s for s in stats if s.smiles == "[1*]C(C)=O")
    assert shared.count == 2 and shared.rows == (0, 1)
    assert shared.mean == pytest.approx(6.0)
    assert shared.sd == pytest.approx(math.sqrt(2.0))
    assert stats[0] is shared  # commonest first


def test_fragment_stats_ignores_missing_activity() -> None:
    found = fragment_rows([PARACETAMOL, ACETANILIDE])
    stats = fragment_stats(found, np.array([5.0, np.nan]))
    shared = next(s for s in stats if s.smiles == "[1*]C(C)=O")
    assert shared.count == 2
    assert shared.mean == pytest.approx(5.0) and shared.sd == 0.0


def test_fragment_stats_without_activity_has_no_numbers() -> None:
    stats = fragment_stats(fragment_rows([PARACETAMOL]))
    assert math.isnan(stats[0].mean) and math.isnan(stats[0].sd)


def test_a_fragment_repeated_in_one_row_counts_that_row_once() -> None:
    found = fragment_rows(["CC(=O)Nc1ccc(NC(C)=O)cc1"])  # two identical acetamides
    stats = fragment_stats(found)
    acetyl = next(s for s in stats if s.smiles == "[1*]C(C)=O")
    assert acetyl.count == 1 and acetyl.rows == (0,)
