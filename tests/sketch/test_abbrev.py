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

"""Abbreviations (SKETCHER_PLAN §4.5; docs/sketcher/abbreviation-review.md).

Each entry's structure is checked against the formula written by hand from its IUPAC
substituent name, then placed in a molecule, round-tripped through a molfile and drawn.
"""

from __future__ import annotations

import itertools
import math
import re
from collections import Counter

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemsuite.sketch import abbrev, render, valence
from chemsuite.sketch.convert import from_text, to_molblock, to_smiles
from chemsuite.sketch.edits import History
from tests.sketch.helpers import canonical

ENTRIES = [pytest.param(e, id=e.label) for e in abbrev.TABLE]


def _counts(formula: str) -> Counter:
    out: Counter = Counter()
    for symbol, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out[symbol] += int(n or 1)
    return out


@pytest.mark.parametrize("entry", ENTRIES)
def test_structure_matches_the_formula_its_name_gives(entry):
    capped = Chem.MolFromSmiles(entry.smiles.replace("*", "[H]", 1))
    assert capped is not None
    counts = _counts(rdMolDescriptors.CalcMolFormula(capped))
    counts["H"] -= 1  # the hydrogen standing in for the attachment
    counts = +counts
    assert counts == _counts(entry.formula), entry.name


@pytest.mark.parametrize("entry", ENTRIES)
def test_one_attachment_point_on_the_first_atom(entry):
    mol = Chem.MolFromSmiles(entry.smiles)
    star = mol.GetAtomWithIdx(0)
    assert star.GetAtomicNum() == 0 and star.GetDegree() == 1
    assert sum(a.GetAtomicNum() == 0 for a in mol.GetAtoms()) == 1


def test_labels_and_aliases_are_unique():
    keys = [k.lower() for e in abbrev.TABLE for k in (e.label, *e.aliases)]
    assert len(keys) == len(set(keys))


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        ("iBu", "*CC(C)C"),  # 2-methylpropyl; RDKit's default list has sec-butyl here
        ("sBu", "*C(C)CC"),
        ("Mes", "*c1c(C)cc(C)cc1C"),  # Ketcher's fg.sdf has 3-methylbenzyl here
        ("Tol", "*c1ccc(C)cc1"),
        ("THP", "*C1CCCCO1"),  # attached at C2, next to the ring oxygen
        ("Bpin", "*B1OC(C)(C)C(C)(C)O1"),
    ],
)
def test_isomers_stated_by_the_names(label, expected):
    assert canonical(abbrev.find(label).smiles) == canonical(expected)


@pytest.mark.parametrize(
    ("text", "kind", "value"),
    [
        ("Ac", "abbrev", "Ac"),  # acetyl, not actinium
        ("Pr", "abbrev", "Pr"),  # propyl, not praseodymium
        ("Ts", "abbrev", "Ts"),  # tosyl, not tennessine
        ("Co", "atom", "Co"),  # cobalt
        ("co2me", "abbrev", "CO2Me"),
        ("COOH", "abbrev", "CO2H"),
        ("TBDMS", "abbrev", "TBS"),
        ("Z", "abbrev", "Cbz"),
        ("Cl", "atom", "Cl"),
        ("OH", "atom", "O"),
        ("NH2", "atom", "N"),
        ("D", "atom", "H"),
    ],
)
def test_parse_label(text, kind, value):
    parsed = abbrev.parse_label(text)
    assert parsed[0] == kind
    assert (parsed[1].label if kind == "abbrev" else parsed[1]) == value


def test_parse_label_isotope_and_charge():
    assert abbrev.parse_label("13C") == ("atom", "C", 13, 0)
    assert abbrev.parse_label("N+") == ("atom", "N", 0, 1)
    assert abbrev.parse_label("O-") == ("atom", "O", 0, -1)
    assert abbrev.parse_label("Fe2+") == ("atom", "Fe", 0, 2)
    assert abbrev.parse_label("Xyz") is None


def _on_benzene(label: str) -> History:
    history = History(from_text("c1ccccc1C"))
    methyl = next(a for a in history.sketch.atoms if history.sketch.degree(a) == 1)
    change = history.change("label")
    assert abbrev.set_label(change, methyl, label) == ""
    history.commit(change)
    return history


@pytest.mark.parametrize("entry", ENTRIES)
def test_placed_on_a_ring_gives_the_right_molecule_and_layout(entry):
    history = _on_benzene(entry.label)
    sketch = history.sketch
    assert to_smiles(sketch) == canonical(entry.smiles.replace("*", "c1ccccc1", 1))
    assert valence.check(sketch).valid
    pts = [a.pos for a in sketch.atoms.values()]
    assert min(math.dist(p, q) for p, q in itertools.combinations(pts, 2)) > 0.5
    (group,) = sketch.groups.values()
    assert group.label == entry.label and not group.expanded
    # The group is the label atom and everything new; the ring is not in it.
    assert len(group.atoms) == Chem.MolFromSmiles(entry.smiles).GetNumAtoms() - 1


def test_label_is_one_undo_step_and_refuses_inner_atoms():
    history = _on_benzene("Boc")
    history.undo()
    assert to_smiles(history.sketch) == "Cc1ccccc1"
    ring_atom = next(a for a in history.sketch.atoms if history.sketch.degree(a) == 3)
    change = history.change("x")
    assert "one bond" in abbrev.set_label(change, ring_atom, "Ph")


@pytest.mark.parametrize("label", ["CO2Me", "Boc", "Ph", "NO2", "OTf", "TBS"])
def test_molfile_keeps_abbreviations(label):
    history = _on_benzene(label)
    block = to_molblock(history.sketch)
    assert "M  STY  1   1 SUP" in block and f"M  SMT   1 {label}" in block
    assert "M  SAP" in block and "M  SBL" in block
    back = from_text(block)
    (group,) = back.groups.values()
    assert group.label == label and not group.expanded
    assert to_smiles(back) == to_smiles(history.sketch)
    # RDKit reads the same molecule; the abbreviation never changes the chemistry.
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(block)) == to_smiles(history.sketch)


def test_expanded_state_survives_a_molfile():
    history = _on_benzene("Boc")
    gid = next(iter(history.sketch.groups))
    change = history.change("expand")
    abbrev.set_expanded(change, gid, True)
    history.commit(change)
    block = to_molblock(history.sketch)
    assert "M  SDS EXP" in block
    assert next(iter(from_text(block).groups.values())).expanded


def test_contracted_group_draws_as_its_label():
    sketch = _on_benzene("CO2Me").sketch
    scene = render.build(sketch, valence.check(sketch))
    (group,) = sketch.groups.values()
    hidden = group.atoms - {group.attach}
    assert not hidden & set(scene.labels)
    assert all(
        sketch.bonds[line.bond].begin not in hidden and sketch.bonds[line.bond].end not in hidden
        for line in scene.lines
    )
    label = scene.labels[group.attach]
    text = "".join(t for t, _ in label.runs)
    neighbour = next(n for n in sketch.neighbours(group.attach) if n not in group.atoms)
    from_right = sketch.atoms[neighbour].x > sketch.atoms[group.attach].x + 0.1
    assert text == ("MeO2C" if from_right else "CO2Me")
    assert ("2", "sub") in label.runs


def test_label_mirrors_when_the_bond_comes_from_the_right():
    sketch = _on_benzene("CO2Me").sketch
    (group,) = sketch.groups.values()
    anchor = sketch.atoms[group.attach]
    for aid, atom in list(sketch.atoms.items()):  # mirror the drawing left-right
        sketch.atoms[aid] = atom.moved(-atom.x, atom.y)
    scene = render.build(sketch, valence.check(sketch))
    texts = {"".join(t for t, _ in scene.labels[group.attach].runs)}
    assert texts <= {"CO2Me", "MeO2C"} and anchor is not None
    neighbour = next(n for n in sketch.neighbours(group.attach) if n not in group.atoms)
    from_right = sketch.atoms[neighbour].x > sketch.atoms[group.attach].x
    assert texts == {"MeO2C" if from_right else "CO2Me"}


def test_erasing_an_inner_atom_dissolves_the_group():
    from chemsuite.sketch import edits

    history = _on_benzene("Boc")
    (group,) = history.sketch.groups.values()
    inner = next(iter(group.atoms - {group.attach}))
    change = history.change("erase")
    edits.erase(change, atoms=[inner])
    assert not history.sketch.groups
