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
"""Poses: a molecule's own coordinates as a residue, merged into a protein."""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.protein.chem import interactions, ligands, poses
from chemsuite.protein.structure.model import AROMATIC, LIGAND

from .test_ligands import only_ligand


def moved(mol: Chem.Mol, by: float) -> Chem.Mol:
    out = Chem.Mol(mol)
    conformer = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (p.x + by, p.y, p.z))
    return out


def test_a_molecule_becomes_one_ligand_residue(kinase) -> None:
    mol = ligands.to_mol(kinase, only_ligand(kinase)).mol
    pose = poses.from_mol(mol)
    assert pose.n_atoms == mol.GetNumAtoms() and pose.n_residues == 1 == pose.n_chains
    assert pose.res_kind[0] == LIGAND and pose.residue_label(0) == "pose/LIG 1"
    assert np.allclose(pose.xyz(), mol.GetConformer().GetPositions(), atol=1e-4)
    assert (pose.order == AROMATIC).sum() == 10 and len(pose.bonds) == mol.GetNumBonds()
    assert len(set(pose.name.tolist())) == pose.n_atoms  # C1, C2…, N1…: each its own name


def test_only_3d_molecules_are_poses() -> None:
    flat = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(flat)
    assert poses.pose_of(flat) is None and poses.pose_of(None) is None
    solid = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(solid, randomSeed=1)
    assert poses.pose_of(solid) is solid


def test_a_merged_pose_is_the_last_chain_in_every_model(nmr, kinase) -> None:
    mol = ligands.to_mol(kinase, only_ligand(kinase)).mol
    pose = poses.from_mol(mol)
    both = poses.merged(nmr, pose)
    assert both.n_atoms == nmr.n_atoms + pose.n_atoms and both.n_models == nmr.n_models == 2
    assert both.uid[nmr.n_atoms] == nmr.next_uid and both.next_uid == nmr.next_uid + pose.n_atoms
    assert str(both.chain_name[-1]) == "pose" and both.residue_label(both.n_residues - 1)
    for model in range(2):
        assert np.allclose(both.xyz(model)[nmr.n_atoms :], pose.xyz())
        assert np.array_equal(both.xyz(model)[: nmr.n_atoms], nmr.xyz(model))
    assert both.bonds[len(nmr.bonds) :].min() >= nmr.n_atoms


def test_a_row_read_from_the_protein_is_the_same_pose(kinase) -> None:
    residue = only_ligand(kinase)
    mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(ligands.to_mol(kinase, residue).mol))
    assert poses.same_pose(kinase, residue, mol)
    assert not poses.same_pose(kinase, residue, moved(mol, 0.5))
    assert not poses.same_pose(kinase, residue - 1, mol)  # another residue


def test_a_pose_stands_in_for_the_ligand_it_overlaps(kinase) -> None:
    residue = only_ligand(kinase)
    mol = moved(ligands.to_mol(kinase, residue).mol, 0.8)
    replaced = poses.overlapping(kinase, mol.GetConformer().GetPositions())
    assert np.array_equal(np.flatnonzero(replaced), np.flatnonzero(kinase.residue_of == residue))
    far = moved(mol, 40.0)
    assert not poses.overlapping(kinase, far.GetConformer().GetPositions()).any()


def test_a_merged_pose_finds_its_interactions(kinase) -> None:
    residue = only_ligand(kinase)
    mol = ligands.to_mol(kinase, residue).mol
    keep = ~poses.overlapping(kinase, mol.GetConformer().GetPositions())
    both = poses.merged(kinase.subset(keep), poses.from_mol(mol))
    mine = np.zeros(both.n_atoms, bool)
    mine[-mol.GetNumAtoms() :] = True
    found = interactions.find(both, mine)
    labels = {both.residue_label(int(both.residue_of[b])) for b in found.b}
    assert "E/LYS 72" in labels and not mine[found.b].any()
