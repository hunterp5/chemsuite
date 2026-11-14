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

"""R-groups around each molecule's own most central ring system."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.sar.central_ring import (
    central_ring_rows,
    decompose_central_ring,
    most_central_ring_atoms,
)

PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"


def split(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return decompose_central_ring(mol)


def core_atoms(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return most_central_ring_atoms(mol)


def test_a_ring_with_two_substituents_is_cut_at_both() -> None:
    found = split(PARACETAMOL)
    assert found is not None
    assert found.core == "c1cc([*:2])ccc1[*:1]"
    assert set(found.groups) == {"*NC(C)=O", "*O"}


def test_substituents_are_labelled_as_r_group_decomposition_labels_them() -> None:
    found = split("Cc1ccc(C(=O)Nc2ccc(Cl)cc2)cc1")
    assert found is not None
    assert "*C" in found.groups  # a methyl, not [*:1]C
    assert "*C(=O)Nc1ccc(Cl)cc1" in found.groups


def test_biphenyl_keeps_one_ring_as_the_core() -> None:
    found = split("c1ccc(-c2ccccc2)cc1")
    assert found is not None
    assert found.core == "c1ccc([*:1])cc1"
    assert found.groups == ("*c1ccccc1",)  # the other ring is a substituent


def test_a_fused_system_stays_whole() -> None:
    found = split("c1ccc2[nH]ccc2c1")  # indole
    assert found is not None
    assert found.core == "c1ccc2[nH]ccc2c1"
    assert found.groups == ()
    assert len(core_atoms("c1ccc2[nH]ccc2c1")) == 9


def test_a_spiro_system_stays_whole() -> None:
    atoms = core_atoms("C1CCC2(CC1)CCCCC2")
    assert atoms is not None and len(atoms) == 11


def test_the_middle_ring_of_three_wins() -> None:
    found = split("O=C(Nc1ccccc1)c1ccc(N2CCOCC2)cc1")
    assert found is not None
    assert found.core == "c1cc([*:2])ccc1[*:1]"
    assert set(found.groups) == {"*C(=O)Nc1ccccc1", "*N1CCOCC1"}


def test_a_molecule_with_no_ring_gives_nothing() -> None:
    assert split("CCCC") is None
    assert core_atoms("CCCC") is None


def test_a_ring_with_nothing_on_it_has_no_substituents() -> None:
    found = split("c1ccccc1")
    assert found is not None
    assert found.core == "c1ccccc1" and found.groups == ()


def test_an_empty_molecule_gives_nothing() -> None:
    assert decompose_central_ring(None) is None
    assert decompose_central_ring(Chem.Mol()) is None


def test_the_same_scaffold_numbers_its_positions_the_same_way() -> None:
    first = split("Clc1ccc(C)cc1")
    second = split("Brc1ccc(C)cc1")
    assert first is not None and second is not None
    assert first.core == second.core
    # the methyl sits at the same position in both
    assert first.groups.index("*C") == second.groups.index("*C")


def test_rows_keep_their_order_and_blanks() -> None:
    found = central_ring_rows([PARACETAMOL, "Cc1ccc(C(=O)N)cc1", None, "not a structure"])
    assert len(found) == 4
    assert found.cores[0] == "c1cc([*:2])ccc1[*:1]"
    assert found.cores[2] == "" and found.cores[3] == ""
    assert found.failed == 1  # the blank is not a failure
    assert found.positions() == ("R1", "R2")
    assert found.labels["R1"][2] == "" and found.labels["R2"][3] == ""
    assert found.decomposed() == 2


def test_rows_with_different_position_counts_line_up() -> None:
    found = central_ring_rows(["c1ccccc1", PARACETAMOL])
    assert found.positions() == ("R1", "R2")
    assert found.labels["R1"] == ["", "*NC(C)=O"]
    assert len(found.labels["R2"]) == 2


def test_rows_report_progress() -> None:
    seen: list[tuple[int, int]] = []
    central_ring_rows([PARACETAMOL] * 25, progress=lambda done, total: seen.append((done, total)))
    assert seen and seen[-1] == (25, 25)


def test_cancelling_keeps_the_rows_already_done() -> None:
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done >= 20:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        central_ring_rows([PARACETAMOL] * 100, progress=progress, cancel=token)
    partial = stopped.value.partial
    assert partial is not None and 0 < len(partial) < 100
    assert partial.cores[0] == "c1cc([*:2])ccc1[*:1]"
    assert all(len(column) == len(partial) for column in partial.labels.values())
