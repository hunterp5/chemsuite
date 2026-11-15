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
"""The PDBQT files Vina reads and writes (Meeko)."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from chemsuite.protein.dock import pdbqt
from chemsuite.protein.structure.model import AMINO


def test_the_receptor_is_the_polymer_and_what_meeko_can_type(pocket) -> None:
    receptor = pocket.receptor
    lines = [line for line in receptor.pdbqt.splitlines() if line.startswith("ATOM")]
    assert len(lines) > 500
    assert all(line.split()[-1] for line in lines)  # every atom carries an AutoDock type
    # the ligand and the crystallisation additives are not part of a receptor
    assert "80U" in receptor.dropped
    assert not any(line[17:20].strip() == "80U" for line in lines)
    assert all(name == name.strip().upper() for name in receptor.dropped)


def test_the_receptor_can_be_narrowed_to_a_part_of_the_protein(pocket) -> None:
    s = pocket.structure
    near = np.isin(s.atom_kind, [AMINO]) & (s.atom_chain == 0)
    narrowed = pdbqt.receptor(s, keep=near)
    assert 0 < len(narrowed.pdbqt.splitlines()) <= len(pocket.receptor.pdbqt.splitlines())
    assert "80U" not in narrowed.dropped  # it was not asked for in the first place


def test_a_receptor_meeko_knows_nothing_of_is_an_error(pocket) -> None:
    s = pocket.structure
    ligand_only = s.residue_of == pocket.residue
    with pytest.raises(pdbqt.PdbqtError, match="receptor"):
        pdbqt.receptor(s, keep=ligand_only)


def test_meeko_never_has_to_download_a_template(pocket) -> None:
    """Every residue that goes into a receptor is one Meeko already has (it fetches otherwise)."""
    known = pdbqt.known_residues()
    assert {"ALA", "HIS", "HOH", "ZN"} <= known
    assert "80U" not in known and "NO3" not in known
    assert all(str(name).strip().upper() in known for name in pocket.structure.res_name[
        np.isin(pocket.structure.res_kind, pdbqt.RECEPTOR_KINDS)
        & ~np.isin(pocket.structure.res_name, list(pocket.receptor.dropped))
    ])  # fmt: skip


def test_a_ligand_survives_the_round_trip_through_pdbqt(pocket) -> None:
    text = pdbqt.ligand(pocket.mol)
    assert text.startswith("REMARK") or "ROOT" in text
    assert "TORSDOF" in text
    (back,) = pdbqt.poses(text)
    heavy, grown = Chem.RemoveHs(pocket.mol), Chem.RemoveHs(back)
    assert grown.GetNumAtoms() == heavy.GetNumAtoms()
    assert Chem.MolToSmiles(grown) == Chem.MolToSmiles(heavy)  # bond orders and all
    # Meeko writes the atoms in its own order (the torsion tree's root first), so the pose is
    # compared where it stands, not atom for atom
    assert rdMolAlign.CalcRMS(grown, heavy) < 0.01


def test_a_ligand_with_no_pose_cannot_be_docked() -> None:
    flat = Chem.MolFromSmiles("CCO")
    with pytest.raises(pdbqt.PdbqtError, match="3D"):
        pdbqt.ligand(flat)


def test_hydrogens_are_added_when_a_pose_has_none() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    text = pdbqt.ligand(Chem.RemoveHs(mol))
    assert "HD" in text  # the hydroxyl's polar hydrogen is there to donate


def test_output_that_is_not_a_pose_is_an_error() -> None:
    with pytest.raises(pdbqt.PdbqtError):
        pdbqt.poses("not a PDBQT file at all\n")
