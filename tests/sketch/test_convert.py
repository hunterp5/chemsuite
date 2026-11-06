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

"""Round trips: SMILES or molfile → sketch → molfile → RDKit gives the same molecule.

The drug set covers stereocentres, E/Z double bonds, charges, isotopes, radicals, salts and
heteroaromatics.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.sketch import edits
from chemsuite.sketch.convert import (
    SketchError,
    bond_length_spread,
    drawn_molblock,
    from_text,
    to_mol,
    to_molblock,
    to_smiles,
)
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Order, Sketch, Stereo
from tests.sketch.helpers import canonical

DRUGS = [
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1",  # (S)-ibuprofen
    "CN1CCC[C@H]1c1cccnc1",  # nicotine
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",  # caffeine
    "CC(C)NC[C@@H](O)COc1cccc2ccccc12",  # (S)-propranolol
    "C[C@@H](N)C(=O)O",  # L-alanine
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",  # L-tryptophan
    "CC/C(=C(/CC)c1ccccc1)c1ccc(OCCN(C)C)cc1",  # tamoxifen (Z)
    "C/C(=C\\CO)CCC=C(C)C",  # geraniol (E)
    "OC(=O)/C=C\\C(=O)O",  # maleic acid
    "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",  # penicillin G
    "C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",  # a steroid
    "[O-][N+](=O)c1ccc(Cl)cc1",  # nitro, charge-separated
    "C[N+](C)(C)C.[Cl-]",  # salt
    "[Na+].[O-]C(=O)c1ccccc1",
    "F[P-](F)(F)(F)(F)F.CCCC[N+](CCCC)(CCCC)CCCC",
    "[2H]C([2H])([2H])Oc1ccccc1",  # isotopes
    "[13CH3]O",
    "[CH2]CC",  # radical
    "C#CCO",
    "C=C=C",
    "O=S(=O)(N)c1ccc(N)cc1",  # sulfanilamide
    "Clc1ccc2c(c1)C(c1ccccc1F)=NCC(=O)N2",  # a benzodiazepine
    "c1ccc2[nH]ccc2c1",  # indole
    "CC1=CC(=O)C=CC1=O",  # quinone
]


@pytest.mark.parametrize("smiles", DRUGS)
def test_smiles_round_trip(smiles):
    sketch = from_text(smiles)
    assert to_smiles(sketch) == canonical(smiles)


@pytest.mark.parametrize("smiles", DRUGS)
def test_molfile_round_trip_keeps_molecule_and_wedges(smiles):
    sketch = from_text(smiles)
    block = to_molblock(sketch)
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(block)) == canonical(smiles)
    again = from_text(block)
    assert to_smiles(again) == canonical(smiles)
    wedges = sorted((b.stereo.value, b.order.value) for b in sketch.bonds.values())
    assert sorted((b.stereo.value, b.order.value) for b in again.bonds.values()) == wedges


@pytest.mark.parametrize("smiles", DRUGS)
def test_imported_structures_have_unit_bonds_and_kekule_rings(smiles):
    sketch = from_text(smiles)
    assert bond_length_spread(sketch) < 0.15
    assert all(b.order is not Order.AROMATIC for b in sketch.bonds.values())


def test_drawn_wedge_survives_export_as_drawn():
    """A hash drawn on C-F is written on C-F, not moved by RDKit's own wedging."""
    sketch = from_text("C[C@@H](F)Cl")
    centre = next(a for a, atom in sketch.atoms.items() if sketch.degree(a) == 3)
    f = next(a for a, atom in sketch.atoms.items() if atom.symbol == "F")
    history = History(sketch)
    change = history.change("hash")
    for bid, bond in list(sketch.bonds.items()):
        sketch.bonds[bid] = bond.__class__(bond.begin, bond.end, bond.order, Stereo.NONE)
    bid = sketch.bond_between(centre, f)
    if sketch.bonds[bid].begin != centre:
        sketch.bonds[bid] = sketch.bonds[bid].flipped()
    edits.click_bond(change, bid, Order.SINGLE, Stereo.HASH)
    block = to_molblock(sketch)
    mol = Chem.MolFromMolBlock(block)
    Chem.ReapplyMolBlockWedging(mol)
    fl = next(
        b
        for b in mol.GetBonds()
        if "F" in (b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol())
    )
    assert fl.GetBondDir() == Chem.BondDir.BEGINDASH
    assert sum(b.GetBondDir() != Chem.BondDir.NONE for b in mol.GetBonds()) == 1


def test_wedge_direction_sets_configuration():
    """Flipping which way C-F is hashed in the drawing inverts R/S."""
    sketch = from_text("C[C@@H](F)Cl")
    before = to_smiles(sketch)
    bid, bond = next((i, b) for i, b in sketch.bonds.items() if b.stereo is not Stereo.NONE)
    other = Stereo.HASH if bond.stereo is Stereo.WEDGE else Stereo.WEDGE
    sketch.bonds[bid] = bond.__class__(bond.begin, bond.end, bond.order, other)
    assert to_smiles(sketch) == canonical("C[C@H](F)Cl") != before


def test_drawn_double_bond_geometry_sets_e_or_z():
    history = History()
    change = history.change("draw")
    atoms = edits.chain(change, (0.0, 0.0), 0.0, 3)  # zig-zag: trans
    bid = history.sketch.bond_between(atoms[1], atoms[2])
    edits.set_bond_order(change, bid, Order.DOUBLE)
    edits.set_element(change, atoms[0], "F")
    edits.set_element(change, atoms[3], "Cl")
    assert to_smiles(history.sketch) == canonical("F/C=C/Cl")


def test_invalid_drawing_raises_but_still_copies_as_drawn():
    sketch = from_text("CC(C)(C)C")
    centre = next(a for a in sketch.atoms if sketch.degree(a) == 4)
    change = History(sketch).change("x")
    edits.bond_from_atom(change, centre)
    with pytest.raises(SketchError):
        to_mol(sketch)
    assert to_molblock(sketch) == drawn_molblock(sketch)


def test_text_that_is_not_a_structure_gives_none():
    assert from_text("not a molecule") is None
    assert to_smiles(Sketch()) == ""


def test_charges_and_isotopes_over_eight_atoms_split_lines():
    smiles = ".".join(["[Na+]"] * 10)
    block = drawn_molblock(from_text(smiles))
    assert sum(line.startswith("M  CHG") for line in block.splitlines()) == 2
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(block)) == canonical(smiles)


def test_chemdraw_style_molfile_is_scaled_to_unit_bonds():
    mol = Chem.MolFromSmiles("c1ccccc1O")
    from rdkit.Chem import rdDepictor

    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x * 0.55, p.y * 0.55, 0.0))  # ChemDraw's ~0.825 bonds
    sketch = from_text(Chem.MolToMolBlock(mol))
    assert bond_length_spread(sketch) < 0.05
