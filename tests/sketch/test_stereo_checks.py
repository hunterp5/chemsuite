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

"""Stereo beyond R/S labels (SKETCHER_PLAN §4.4): wedge checks, "either" double bonds,
enhanced stereo, and hard CIP cases from the CIP Validation Suite (Hanson et al. 2018).

The whole suite is ``scripts/review/cip_review.py``; docs/sketcher/stereo-review.md holds its
results (the sketcher gives exactly RDKit's labels on all 266 tetrahedral and E/Z records).
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.sketch import convert, stereo
from chemsuite.sketch.edits import History, set_enhanced
from chemsuite.sketch.model import Atom, Bond, Order, Sketch, Stereo

W, H, N = Stereo.WEDGE, Stereo.HASH, Stereo.NONE


def _drawing(atoms, bonds) -> Sketch:
    sketch = Sketch()
    for symbol, x, y in atoms:
        sketch.atoms[sketch.new_id()] = Atom(x, y, symbol=symbol)
    for a, b, wedge in bonds:
        sketch.bonds[sketch.new_id()] = Bond(a, b, Order.SINGLE, wedge)
    return sketch


FOUR = [("C", 0, 0), ("N", 0, 1), ("O", 0, -1), ("F", -1, 0), ("Cl", 1, 0)]
THREE_GOOD = [("C", 0, 0), ("N", 0, 1), ("O", -0.87, -0.5), ("F", 0.87, -0.5)]
THREE_T = [("C", 0, 0), ("N", -1, 0), ("O", 1, 0), ("F", 0, 1)]  # two plain bonds in a line


def test_sound_wedges_raise_nothing():
    assert stereo.wedge_problems(_drawing(THREE_GOOD, [(1, 2, W), (1, 3, N), (1, 4, N)])) == {}
    assert stereo.wedge_problems(convert.from_text("C[C@H](N)C(=O)O")) == {}
    assert stereo.wedge_problems(convert.from_text(
        "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O")) == {}  # fmt: skip


def test_opposite_wedge_and_hash_contradict():
    problems = stereo.wedge_problems(_drawing(FOUR, [(1, 2, W), (1, 3, H), (1, 4, N), (1, 5, N)]))
    assert set(problems.values()) == {stereo.CONTRADICTORY} and len(problems) == 2


def test_wedge_between_two_collinear_plain_bonds_is_ambiguous():
    problems = stereo.wedge_problems(_drawing(THREE_T, [(1, 2, N), (1, 3, N), (1, 4, W)]))
    assert list(problems.values()) == [stereo.AMBIGUOUS]


def test_wedge_at_a_non_centre_and_a_backwards_wedge():
    propane = _drawing([("C", 0, 0), ("C", 1, 0), ("C", 2, 0)], [(1, 2, W), (2, 3, N)])
    assert list(stereo.wedge_problems(propane).values()) == [stereo.NO_CENTRE]
    backwards = _drawing(
        [("C", -1, 0), ("C", 0, 0), ("N", 0, 1), ("O", 0.87, -0.5)],
        [(1, 2, W), (2, 3, N), (2, 4, N)],
    )  # narrow end on the methyl, wide end on the stereocentre
    assert list(stereo.wedge_problems(backwards).values()) == [stereo.WRONG_END]


def test_wavy_bonds_are_never_problems():
    sketch = _drawing(THREE_T, [(1, 2, N), (1, 3, N), (1, 4, Stereo.EITHER)])
    assert stereo.wedge_problems(sketch) == {}


# -- "either" double bonds ---------------------------------------------------------------------


def _double(sketch):
    return next(b for b, v in sketch.bonds.items() if v.order is Order.DOUBLE)


def test_either_double_bond_round_trips_as_ctfile_stereo_3():
    history = History(convert.from_text("C/C=C/Cl"))
    bid = _double(history.sketch)
    change = history.change("Unknown")
    assert stereo.set_double_bond(change, bid, stereo.UNKNOWN) == ""
    history.commit(change)
    found = stereo.labels(history.sketch)
    assert bid in found.unknown_bonds and bid not in found.bonds
    block = convert.to_molblock(history.sketch)
    assert any(line.split()[2:4] == ["2", "3"] for line in block.splitlines()[4:])
    back = convert.from_text(block)
    assert back.bonds[_double(back)].stereo is Stereo.EITHER
    assert convert.to_smiles(back) == "CC=CCl"  # no E/Z in the SMILES


def test_either_double_bond_back_to_e_or_z():
    history = History(convert.from_text("C/C=C/Cl"))
    bid = _double(history.sketch)
    for want in (stereo.UNKNOWN, "Z"):
        change = history.change(want)
        assert stereo.set_double_bond(change, bid, want) == ""
        history.commit(change)
    assert convert.to_smiles(history.sketch) == r"C/C=C\Cl"


def test_wavy_bond_beside_a_double_bond_is_not_doubled_on_import():
    mol = Chem.MolFromSmiles("CC=CCl")
    from rdkit.Chem import rdDepictor

    rdDepictor.Compute2DCoords(mol)
    mol.GetBondWithIdx(0).SetBondDir(Chem.BondDir.UNKNOWN)
    block = Chem.MolToMolBlock(mol)
    sketch = convert.from_text(block)
    assert sketch.bonds[_double(sketch)].stereo is Stereo.NONE


# -- enhanced stereo ---------------------------------------------------------------------------


def test_enhanced_stereo_groups_round_trip_through_v3000():
    sketch = convert.from_text("C[C@H](N)C[C@@H](O)C")
    centres = sorted(stereo.labels(sketch).atoms)
    history = History(sketch)
    change = history.change("Groups")
    set_enhanced(change, [centres[0]], "and1")
    set_enhanced(change, [centres[1]], "or1")
    history.commit(change)
    mol = convert.to_rdkit(history.sketch)
    kinds = sorted(str(g.GetGroupType()) for g in mol.GetStereoGroups())
    assert kinds == ["STEREO_AND", "STEREO_OR"]
    block = convert.to_molblock(history.sketch)
    assert "V3000" in block and "MDLV30/STERAC1" in block and "MDLV30/STEREL1" in block
    back = convert.from_text(block)
    assert sorted(a.enhanced for a in back.atoms.values() if a.enhanced) == ["and1", "or1"]
    cx = Chem.MolToCXSmiles(convert.to_rdkit(back))
    assert "&1:" in cx and "o1:" in cx


def test_cxsmiles_enhanced_stereo_reads_in():
    sketch = convert.from_text("C[C@H](N)C[C@@H](O)C |&1:1,o2:4|")
    assert sorted(a.enhanced for a in sketch.atoms.values() if a.enhanced) == ["and1", "or2"]


def test_enhanced_group_on_a_non_centre_is_dropped_from_the_molecule():
    sketch = convert.from_text("CCO")
    history = History(sketch)
    change = history.change("Group")
    set_enhanced(change, list(sketch.atoms), "and1")
    history.commit(change)
    assert not convert.to_rdkit(history.sketch).GetStereoGroups()
    assert "V3000" not in convert.to_molblock(history.sketch)


def test_v3000_export_keeps_explicit_valences():
    sketch = convert.from_text("F[P-](F)(F)(F)(F)F.C[C@H](N)C(=O)O")
    centre = next(iter(stereo.labels(sketch).atoms))
    history = History(sketch)
    change = history.change("Group")
    set_enhanced(change, [centre], "and1")
    history.commit(change)
    block = convert.to_molblock(history.sketch)
    assert "V3000" in block
    phosphorus = next(ln for ln in block.splitlines() if " P " in ln and ln.startswith("M  V30"))
    assert "VAL=6" in phosphorus
    assert all(len(ln) <= 80 for ln in block.splitlines())


# -- hard CIP cases (CIP Validation Suite, Hanson et al. 2018) ----------------------------------

#: (suite id, SMILES, suite labels as "<1-based SMILES atom><label>") verbatim from
#: compounds.smi of github.com/CIPValidationSuite/ValidationSuite (Hanson et al., J. Chem. Inf.
#: Model. 2018, 58, 1755). They need sequence rules 3 (Z > E), 4b (like/unlike) and 5
#: (pseudo-asymmetry, r/s), here checked through the draw -> molfile -> RDKit path.
HARD = [
    ("VS192", r"C/C=C(/CCC[C@@H](CCC/C(=C\C)/Cl)O)\Cl", "2Z 3Z 7R 11E 12E"),
    ("VS194", r"C/C=C\CCC[C@@H](CCC/C=C/C)O", "2Z 3Z 7R 11E 12E"),
    ("VS196", "Cl[C@H]1[C@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl", "2R 3S 4S 5R 6S 7R"),
    ("VS199", "O=C([C@H]([C@H]([C@@H]([C@@H]([C@H](C(O)=O)Cl)Cl)Cl)Cl)Cl)O", "3R 4R 5R 6R 7S"),
    ("VS213", "O[C@H]1CC[C@@](CC1)(Cl)C", "2r 5r"),
    ("VS219", "O[C@H]1C[C@@H](C1)O", "2r 4r"),
    ("VS221", "O[C@@H]1CC[C@H](CC1)O", "2r 5r"),
    ("VS222", "O[C@H]1CC[C@H](CC1)O", "2s 5s"),
    ("VS227", "O=C([C@H]([C@@H]([C@H](C(O)=O)Cl)Cl)Cl)O", "3R 4r 5S"),
    ("VS228", "O=C([C@H]([C@H]([C@H](C(O)=O)Cl)Cl)Cl)O", "3R 4s 5S"),
    ("VS230", "S[C@H]1[C@H](C[C@H](C1)S)S", "2R 3S 5r"),
    ("VS236", "C[C@@H]([C@H]([C@@H](C)S)S)S", "2S 3r 4R"),
    ("VS252", "CO[C@@H]1[C@@H]([C@@H]([C@H]1OC)OC)OC", "3S 4s 5R 6r"),
    ("VS253", "Cl[C@@H]1[C@@H]([C@@H]([C@H]([C@H]([C@H]1Cl)Cl)Cl)Cl)Cl", "2R 3s 4S 5S 6s 7R"),
    ("VS254", "Cl[C@H]1[C@@H]([C@H]([C@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl", "2S 3r 4R 5S 6r 7R"),
]


@pytest.mark.parametrize("smiles, labels", [h[1:] for h in HARD], ids=[h[0] for h in HARD])
def test_hard_cip_cases_through_the_drawing(smiles, labels):
    want = {int(item[:-1]) - 1: item[-1] for item in labels.split()}
    sketch = convert.from_text(smiles)
    assert len(sketch.atoms) == Chem.MolFromSmiles(smiles).GetNumAtoms()  # same numbering
    index = {aid: i for i, aid in enumerate(convert.atom_order(sketch))}
    found = stereo.labels(sketch)
    got = {index[aid]: code for aid, code in found.atoms.items()}
    for bid, code in found.bonds.items():
        bond = sketch.bonds[bid]
        got[index[bond.begin]] = got[index[bond.end]] = code
    assert got == want
