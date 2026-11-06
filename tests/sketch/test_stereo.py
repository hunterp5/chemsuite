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

"""CIP labels on drawings and the edits that set R/S and E/Z (SKETCHER_PLAN §4.4)."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.sketch import convert, stereo, valence
from chemsuite.sketch.edits import History


def labelled(smiles: str):
    sketch = convert.from_text(smiles)
    found = stereo.labels(sketch)
    mol = convert.to_rdkit(sketch)
    order = convert.atom_order(sketch)
    by_symbol = {}
    for aid, code in found.atoms.items():
        atom = mol.GetAtomWithIdx(order.index(aid))
        by_symbol.setdefault(atom.GetSymbol(), []).append(code)
    return sketch, found, by_symbol


@pytest.mark.parametrize(
    ("name", "smiles", "atoms", "bonds"),
    [
        ("L-alanine", "C[C@@H](C(=O)O)N", ["S"], []),
        ("D-alanine", "C[C@H](C(=O)O)N", ["R"], []),
        ("L-threonine (2S,3R)", "C[C@H]([C@@H](C(=O)O)N)O", ["R", "S"], []),
        ("L-allo-threonine (2S,3S)", "C[C@@H]([C@@H](C(=O)O)N)O", ["S", "S"], []),
        ("(R)-glyceraldehyde", "O=C[C@H](O)CO", ["R"], []),
        ("fumaric acid", "OC(=O)/C=C/C(=O)O", [], ["E"]),
        ("maleic acid", r"OC(=O)/C=C\C(=O)O", [], ["Z"]),
        ("(Z)-2-bromo-2-butene", "C/C=C(/C)Br", [], ["Z"]),
        ("(E)-2-bromo-2-butene", r"C/C=C(\C)Br", [], ["E"]),
    ],
)
def test_reference_labels(name, smiles, atoms, bonds):
    sketch = convert.from_text(smiles)
    found = stereo.labels(sketch)
    assert sorted(found.atoms.values()) == atoms, name
    assert sorted(found.bonds.values()) == bonds, name
    # the drawing still means the same molecule
    assert convert.to_smiles(sketch) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_undefined_centres_and_bonds_are_reported():
    found = stereo.labels(convert.from_text("CC(N)C=CC"))
    assert len(found.undefined_atoms) == 1 and not found.atoms
    assert found.bonds  # 2D geometry always defines a drawn double bond's E/Z


def apply(sketch, fn):
    history = History(sketch.copy())
    change = history.change("x")
    problem = fn(change)
    history.commit(change)
    return history.sketch, problem


def test_invert_a_centre_swaps_its_wedges():
    sketch = convert.from_text("C[C@@H](C(=O)O)N")  # S
    centre = next(iter(stereo.labels(sketch).atoms))
    done, problem = apply(sketch, lambda ch: stereo.set_centre(ch, centre, "R"))
    assert problem == ""
    assert stereo.labels(done).atoms == {centre: "R"}
    assert convert.to_smiles(done) == Chem.MolToSmiles(Chem.MolFromSmiles("C[C@H](C(=O)O)N"))


@pytest.mark.parametrize("want", ["R", "S"])
def test_define_an_undefined_centre(want):
    sketch = convert.from_text("CC(N)C(=O)O")
    centre = next(iter(stereo.labels(sketch).undefined_atoms))
    done, problem = apply(sketch, lambda ch: stereo.set_centre(ch, centre, want))
    assert problem == ""
    assert stereo.labels(done).atoms == {centre: want}
    assert sum(1 for b in done.bonds.values() if b.stereo.name in ("WEDGE", "HASH")) == 1


def test_not_a_centre():
    sketch = convert.from_text("CCO")
    first = convert.atom_order(sketch)[0]
    _, problem = apply(sketch, lambda ch: stereo.set_centre(ch, first, "R"))
    assert "not a stereocentre" in problem


def test_switch_e_and_z_keeps_centres():
    smiles = "C/C=C/[C@H](O)CC"
    sketch = convert.from_text(smiles)
    before = stereo.labels(sketch)
    bond = next(iter(before.bonds))
    done, problem = apply(sketch, lambda ch: stereo.set_double_bond(ch, bond, "Z"))
    assert problem == ""
    after = stereo.labels(done)
    assert after.bonds == {bond: "Z"} and after.atoms == before.atoms
    assert valence.check(done).valid


def test_ring_double_bonds_cannot_switch():
    sketch = convert.from_text("C1CCC=CCCC1")  # cyclooctene: E/Z exists, but in a ring
    found = stereo.labels(sketch)
    bond = next(iter(found.bonds))
    _, problem = apply(sketch, lambda ch: stereo.set_double_bond(ch, bond, stereo.opposite(
        found.bonds[bond])))  # fmt: skip
    assert "ring" in problem


def test_show_and_hide_hydrogens_round_trip():
    sketch = convert.from_text("C[C@@H](C(=O)O)N")
    report = valence.check(sketch)
    done, added = apply(sketch, lambda ch: stereo.show_hydrogens(ch, list(sketch.atoms),
                                                              report.hydrogens))  # fmt: skip
    assert added == sum(report.hydrogens.values())
    assert valence.check(done).valid
    assert convert.to_smiles(done) == convert.to_smiles(sketch)
    assert sorted(stereo.labels(done).atoms.values()) == sorted(
        stereo.labels(sketch).atoms.values()
    )
    back, (removed, kept) = apply(done, lambda ch: stereo.hide_hydrogens(ch, list(done.atoms)))
    assert removed == added and kept == 0
    assert convert.to_smiles(back) == convert.to_smiles(sketch)


def test_hydrogen_on_a_wedge_is_kept():
    sketch = convert.from_text("CC(N)C(=O)O")
    centre = next(iter(stereo.labels(sketch).undefined_atoms))
    report = valence.check(sketch)
    drawn, _ = apply(sketch, lambda ch: stereo.show_hydrogens(ch, [centre], report.hydrogens))
    wedged, problem = apply(drawn, lambda ch: stereo.set_centre(ch, centre, "S"))
    assert problem == ""
    h_bond = next(b for b in wedged.bonds.values() if b.stereo.name != "NONE")
    assert wedged.atoms[h_bond.end].symbol == "H"  # the wedge went on the hydrogen
    done, (removed, kept) = apply(wedged, lambda ch: stereo.hide_hydrogens(ch, [centre]))
    assert (removed, kept) == (0, 1)
    assert stereo.labels(done).atoms == {centre: "S"}
