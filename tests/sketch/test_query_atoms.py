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

"""R-groups, wildcards and atom lists: labels, molfile round trips and RDKit's reading."""

from __future__ import annotations

from dataclasses import replace

import pytest
from rdkit import Chem

from chemsuite.core.mols import display_text, has_query_atoms, parse_structure
from chemsuite.sketch import abbrev, convert, edits, render, valence
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import (
    ATOM_LIST,
    RGROUP,
    Atom,
    Bond,
    Order,
    Sketch,
    Stereo,
    next_rgroup,
    query_fields,
)


def chain(*atoms: Atom) -> Sketch:
    sketch = Sketch()
    ids = []
    for i, atom in enumerate(atoms):
        aid = sketch.new_id()
        sketch.atoms[aid] = replace(atom, x=float(i), y=0.0)
        ids.append(aid)
    for a, b in zip(ids, ids[1:], strict=False):
        sketch.bonds[sketch.new_id()] = Bond(a, b, Order.SINGLE, Stereo.NONE)
    return sketch


@pytest.mark.parametrize(
    ("text", "label"),
    [
        ("R", "R1"),
        ("R12", "R12"),
        ("A", "A"),
        ("QH", "QH"),
        ("*", "*"),
        ("[c, n,O]", "[C,N,O]"),
        ("![F,Cl]", "![F,Cl]"),
    ],
)
def test_query_labels(text, label):
    fields = query_fields(text)
    assert fields is not None
    assert Atom(0, 0, **fields).label == label


@pytest.mark.parametrize("text", ["R0", "R100", "[C,Zz]", "[]", "B", "Cl", "Rx"])
def test_not_query_labels(text):
    assert query_fields(text) is None


def test_parse_label_prefers_abbreviations_then_queries_then_elements():
    assert abbrev.parse_label("R2")[0] == "query"
    assert abbrev.parse_label("X")[1]["symbol"] == "X"
    assert abbrev.parse_label("Ac")[0] == "abbrev"  # acetyl, not a query or actinium
    assert abbrev.parse_label("N")[0] == "atom"


def test_molfile_round_trip_keeps_every_query_kind():
    atoms = [
        Atom(0, 0, symbol=RGROUP, rgroup=3),
        Atom(0, 0, symbol=ATOM_LIST, atom_list=("C", "N", "O")),
        Atom(0, 0, symbol="A"),
        Atom(0, 0, symbol=ATOM_LIST, atom_list=("F", "Cl"), list_not=True),
        Atom(0, 0, symbol="Q"),
        Atom(0, 0, symbol="X"),
        Atom(0, 0, symbol="M"),
        Atom(0, 0, symbol="AH"),
        Atom(0, 0, symbol="XH"),
        Atom(0, 0, symbol="*"),
    ]
    sketch = chain(*atoms)
    labels = [a.label for a in atoms]
    for block in (convert.drawn_molblock(sketch), convert.to_molblock(sketch)):
        back = convert.from_text(block)
        assert [back.atoms[i].label for i in convert.atom_order(back)] == labels
    exported = convert.to_molblock(sketch)
    assert "M  RGP  1   1   3" in exported
    assert "M  ALS   4  2 T F   Cl  " in exported  # CTfile layout (RDKit misaligns NOT lists)
    assert "V  " not in exported


def test_rdkit_reads_queries_as_the_ctfile_spec_defines():
    sketch = chain(Atom(0, 0, symbol="Q"), Atom(0, 0, symbol=ATOM_LIST, atom_list=("N", "O")))
    pattern = convert.to_rdkit(sketch)
    assert Chem.MolFromSmiles("NO").HasSubstructMatch(pattern)  # Q = N, list = O
    assert not Chem.MolFromSmiles("CO").HasSubstructMatch(pattern)  # Q excludes carbon
    assert Chem.MolFromSmiles("CCN").HasSubstructMatch(
        convert.to_rdkit(chain(Atom(0, 0, symbol="A"), Atom(0, 0, symbol="N")))
    )


def test_rgroups_follow_rdkit_smiles_convention():
    sketch = chain(Atom(0, 0, symbol=RGROUP, rgroup=2), Atom(0, 0), Atom(0, 0, symbol="O"))
    assert convert.to_smiles(sketch) == "[2*]CO"
    for text in ("[2*]CO", "[*:2]CO"):
        first = convert.from_text(text).atoms
        assert sorted(a.label for a in first.values()) == ["C", "O", "R2"]
    assert [a.label for a in convert.from_text("*CO").atoms.values()][0] == "*"


def test_queries_are_valid_and_labelled_without_hydrogens():
    sketch = chain(Atom(0, 0, symbol=RGROUP, rgroup=1), Atom(0, 0), Atom(0, 0, symbol="A"))
    report = valence.check(sketch)
    assert report.valid
    first, _c, last = convert.atom_order(sketch)
    assert render.shows_label(sketch, first)
    assert render.label_runs(sketch, first, report.hydrogens[first])[0] == [("R1", "normal")]
    assert render.label_runs(sketch, last, report.hydrogens[last])[0] == [("A", "normal")]


def test_set_element_switches_between_queries_and_elements():
    history = History(chain(Atom(0, 0), Atom(0, 0)))
    aid = convert.atom_order(history.sketch)[0]

    def apply(label):
        change = history.change(label)
        edits.set_element(change, aid, label)
        history.commit(change)

    apply("R1")
    assert history.sketch.atoms[aid].label == "R1"
    apply("[N,O]")
    assert history.sketch.atoms[aid].atom_list == ("N", "O")
    apply("N")
    atom = history.sketch.atoms[aid]
    assert (atom.symbol, atom.rgroup, atom.atom_list) == ("N", 0, ())


def test_next_rgroup_fills_gaps():
    sketch = chain(Atom(0, 0, symbol=RGROUP, rgroup=1), Atom(0, 0, symbol=RGROUP, rgroup=3))
    assert next_rgroup(sketch) == 2
    first = convert.atom_order(sketch)[0]
    assert next_rgroup(sketch, skip=first) == 1


def test_query_cells_show_smarts_and_keep_the_query():
    block = convert.to_molblock(chain(Atom(0, 0, symbol="Q"), Atom(0, 0)))
    blob, text = parse_structure(block)
    mol = Chem.Mol(blob)
    assert has_query_atoms(mol)
    assert text == display_text(mol) == "[!#6&!#1]-[#6]"
    rgroup = convert.to_molblock(chain(Atom(0, 0, symbol=RGROUP, rgroup=1), Atom(0, 0)))
    assert parse_structure(rgroup)[1] == "[1*]C"
