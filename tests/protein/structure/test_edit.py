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
"""Edits to a structure keep every atom's uid and refuse what cannot be."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.protein.structure import edit
from chemsuite.protein.structure.model import ADDITIVE, ION, LIGAND, WATER, Structure


def residue(s: Structure, label: str) -> int:
    (found,) = [r for r in range(s.n_residues) if s.residue_label(r) == label]
    return found


def atoms_of(s: Structure, label: str) -> np.ndarray:
    return s.residue_of == residue(s, label)


def test_deleting_a_water_keeps_the_others_and_every_uid(egfr) -> None:
    waters = np.flatnonzero(egfr.res_kind == WATER)
    gone = atoms_of(egfr, "A/HOH 1103")
    edited = edit.delete(egfr, gone)
    assert edited.n_residues == egfr.n_residues - 1
    assert int((edited.res_kind == WATER).sum()) == len(waters) - 1
    assert "A/HOH 1103" not in {edited.residue_label(r) for r in range(edited.n_residues)}
    assert edited.uid.tolist() == egfr.uid[~gone].tolist()
    assert np.array_equal(edited.coords[0], egfr.coords[0][~gone])


def test_an_atom_goes_with_its_alternate_locations(egfr) -> None:
    lys = residue(egfr, "A/LYS 713")
    span = egfr.atoms_of(lys)
    shown = np.zeros(egfr.n_atoms, bool)
    names = egfr.name[span]
    shown[span.start + int(np.flatnonzero(names == "CB")[0])] = True  # CB (A) only
    edited = edit.delete(egfr, shown)
    after = edited.atoms_of(residue(edited, "A/LYS 713"))
    assert "CB" not in edited.name[after].tolist()
    assert egfr.n_atoms - edited.n_atoms == int((names == "CB").sum()) == 2


def test_deleting_everything_is_refused(egfr) -> None:
    with pytest.raises(edit.EditError):
        edit.delete(egfr, np.ones(egfr.n_atoms, bool))
    assert edit.delete(egfr, np.zeros(egfr.n_atoms, bool)) is egfr


def test_element_and_charge(egfr) -> None:
    first = np.zeros(egfr.n_atoms, bool)
    first[0] = True
    made = edit.set_charge(edit.set_element(egfr, first, "cl"), first, -1)
    assert made.element[0] == 17 and made.charge[0] == -1
    assert np.array_equal(made.element[1:], egfr.element[1:])
    assert made.uid.tolist() == egfr.uid.tolist()
    with pytest.raises(edit.EditError):
        edit.set_element(egfr, first, "Xx")
    with pytest.raises(edit.EditError):
        edit.set_charge(egfr, first, 9)


def test_a_renamed_residue_takes_the_kind_of_its_new_name(egfr) -> None:
    water = residue(egfr, "A/HOH 1103")
    additive = int(np.flatnonzero(egfr.res_kind == ADDITIVE)[0])
    edited = edit.rename_residues(egfr, [water], "lig")
    assert edited.res_name[water] == "LIG" and edited.res_kind[water] == LIGAND
    back = edit.rename_residues(edited, [water], "HOH")
    assert back.res_kind[water] == WATER
    assert edit.rename_residues(egfr, [additive], "LIG").res_kind[additive] == LIGAND
    assert edit.rename_residues(egfr, [additive], "EDO").res_kind[additive] == ADDITIVE
    lys = residue(egfr, "A/LYS 713")
    assert edit.rename_residues(egfr, [lys], "ALY").res_kind[lys] == egfr.res_kind[lys]
    with pytest.raises(edit.EditError):
        edit.rename_residues(egfr, [water], "TOOLONG")


def test_an_ion_by_its_new_name(egfr) -> None:
    ion = int(np.flatnonzero(egfr.res_kind == ION)[0])
    edited = edit.rename_residues(egfr, [ion], "BR")
    assert edited.res_kind[ion] == ION


def test_rename_a_chain_unless_taken(kinase) -> None:
    edited = edit.rename_chain(kinase, 0, "X")
    assert edited.chain_name.tolist() == ["X", *kinase.chain_name.tolist()[1:]]
    with pytest.raises(edit.EditError):
        edit.rename_chain(kinase, 0, str(kinase.chain_name[1]))
    with pytest.raises(edit.EditError):
        edit.rename_chain(kinase, 0, "")


def test_renumber_refuses_a_number_already_taken(egfr) -> None:
    waters = np.flatnonzero(egfr.res_kind == WATER)[:3]
    edited = edit.renumber(egfr, waters, 2001)
    assert edited.res_seq[waters].tolist() == [2001, 2002, 2003]
    assert edited.uid.tolist() == egfr.uid.tolist()
    with pytest.raises(edit.EditError, match="already has that number"):
        edit.renumber(egfr, waters, int(egfr.res_seq[0]))


def test_renumber_one_chain_at_a_time(kinase) -> None:
    ends = [kinase.residues_of(c).start for c in range(2)]
    with pytest.raises(edit.EditError):
        edit.renumber(kinase, ends, 5000)


def test_added_atoms_go_into_their_residues_with_new_uids(crambin) -> None:
    from chemsuite.protein.structure import store

    first, third = 0, 2
    n = crambin.n_atoms
    n_first = crambin.atoms_of(first).stop
    ca = n_first - len(crambin.name[crambin.atoms_of(first)]) + 1  # residue 0's CA
    assert crambin.name[ca] == "CA"
    xyz = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    edited = edit.add_atoms(
        crambin,
        residues=np.array([third, first, first]),
        names=np.array(["X3", "X1", "X2"]),
        elements=np.array([1, 1, 6]),
        xyz=xyz,
        bonds=np.array([[ca, n + 1], [n + 1, n + 2]]),
        order=np.array([1, 2]),
    )
    assert edited.n_atoms == n + 3
    assert edited.next_uid == crambin.next_uid + 3
    # old atoms keep their uids and order; the new ones end their residues, in the order given
    assert edited.uid[np.isin(edited.uid, crambin.uid)].tolist() == crambin.uid.tolist()
    span = edited.atoms_of(first)
    assert edited.name[span][-2:].tolist() == ["X1", "X2"]
    assert edited.name[edited.atoms_of(third)][-1] == "X3"
    added = edited.indices(crambin.next_uid + np.arange(3))
    assert edited.name[added].tolist() == ["X3", "X1", "X2"]
    assert np.allclose(edited.coords[0][added], xyz)
    assert edited.element[added].tolist() == [1, 1, 6]
    # every old bond still joins the same uids, and the new ones are in
    old_pairs = {tuple(sorted(crambin.uid[b])) for b in crambin.bonds}
    new_pairs = {tuple(sorted(edited.uid[b])) for b in edited.bonds}
    x1, x2 = (int(u) for u in crambin.next_uid + np.array([1, 2]))
    assert new_pairs == old_pairs | {tuple(sorted((int(crambin.uid[ca]), x1))), (x1, x2)}
    assert (edited.bonds[:, 0] < edited.bonds[:, 1]).all()
    assert edited.n_residues == crambin.n_residues
    assert store.unpack(store.pack(edited)).n_atoms == n + 3


def test_added_atoms_need_a_place_in_every_model(nmr) -> None:
    one = np.zeros((1, 3))
    grown = edit.add_atoms(nmr, np.array([0]), np.array(["X"]), np.array([1]), one)
    assert grown.coords.shape == (2, nmr.n_atoms + 1, 3)
    with pytest.raises(edit.EditError):
        edit.add_atoms(nmr, np.array([0]), np.array(["X"]), np.array([1]), np.zeros((3, 1, 3)))
