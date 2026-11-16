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
"""Poses: a molecule's own 3D coordinates placed in a protein (a docked pose, a ligand moved in
the sketcher) rather than one of the protein's residues.

:func:`from_mol` makes a pose a one-residue structure; :func:`merged` puts it in the protein
(uids after the protein's), so selections, pockets and interactions work on it as on any
ligand. :func:`same_pose` says whether a molecule is just a residue of the protein as it sits
in the file (a ligand row read from the protein), and :func:`overlapping` which of the
protein's ligands a pose stands in for.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from ..structure.model import AROMATIC, COIL, LIGAND, Structure

POSE_NAME = "LIG"
POSE_CHAIN = "pose"
#: Å: a molecule this close to a residue's atoms, atom for atom, is that residue
SAME = 0.1
#: Å between heavy atoms: a ligand of the protein this close to the pose is the one it replaces
OVERLAP = 1.5


def pose_of(mol: Chem.Mol | None) -> Chem.Mol | None:
    """``mol`` when it has 3D coordinates (a pose), else None."""
    if mol is None or not mol.GetNumAtoms() or not mol.GetNumConformers():
        return None
    conformer = mol.GetConformer()
    if not conformer.Is3D():
        return None
    return mol


def from_mol(mol: Chem.Mol, name: str = POSE_NAME) -> Structure:
    """``mol`` (with any hydrogens it has) as a structure of one ligand residue."""
    n = mol.GetNumAtoms()
    atoms = list(mol.GetAtoms())
    seen: dict[str, int] = {}
    names = []
    for atom in atoms:
        symbol = atom.GetSymbol().upper()
        seen[symbol] = seen.get(symbol, 0) + 1
        names.append(f"{symbol}{seen[symbol]}"[:4])
    pairs, orders = [], []
    for bond in mol.GetBonds():
        a, b = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        pairs.append((a, b))
        orders.append(AROMATIC if bond.GetIsAromatic() else int(bond.GetBondTypeAsDouble()))
    xyz = mol.GetConformer().GetPositions().astype(np.float32)
    return Structure(
        uid=np.arange(n, dtype=np.int64),
        coords=xyz.reshape(1, n, 3),
        element=np.array([a.GetAtomicNum() for a in atoms], dtype=np.uint8),
        name=np.array(names, dtype="<U4"),
        altloc=np.full(n, "", dtype="<U1"),
        occupancy=np.ones(n, np.float32),
        b=np.zeros(n, np.float32),
        charge=np.array([a.GetFormalCharge() for a in atoms], dtype=np.int8),
        res_start=np.array([0, n], np.int32),
        res_name=np.array([name], dtype="<U5"),
        res_seq=np.ones(1, np.int32),
        icode=np.array([""], dtype="<U1"),
        res_kind=np.array([LIGAND], np.uint8),
        ss=np.array([COIL], np.uint8),
        seqres=np.full(1, -1, np.int32),
        label_asym=np.array([""], dtype="<U4"),
        entity=np.full(1, -1, np.int16),
        chain_start=np.array([0, 1], np.int32),
        chain_name=np.array([POSE_CHAIN], dtype="<U4"),
        entity_id=np.zeros(0, dtype="<U8"),
        entity_seq=np.zeros(0, dtype="<U5"),
        entity_seq_start=np.zeros(1, np.int32),
        bonds=np.array(pairs, np.int32).reshape(-1, 2),
        order=np.array(orders, np.int8),
        next_uid=n,
    )


def merged(protein: Structure, pose: Structure) -> Structure:
    """``protein`` with ``pose``'s atoms after its own, as its last chain; the pose's atoms
    take uids from the protein's ``next_uid`` on and its coordinates in every model."""
    n, r = protein.n_atoms, protein.n_residues
    pose_xyz = np.broadcast_to(pose.coords[:1], (protein.n_models, pose.n_atoms, 3))

    return Structure(
        uid=np.concatenate([protein.uid, pose.uid + protein.next_uid]),
        coords=np.concatenate([protein.coords, pose_xyz], axis=1),
        element=np.concatenate([protein.element, pose.element]),
        name=np.concatenate([protein.name, pose.name]),
        altloc=np.concatenate([protein.altloc, pose.altloc]),
        occupancy=np.concatenate([protein.occupancy, pose.occupancy]),
        b=np.concatenate([protein.b, pose.b]),
        charge=np.concatenate([protein.charge, pose.charge]),
        res_start=np.concatenate([protein.res_start, pose.res_start[1:] + n]).astype(np.int32),
        res_name=np.concatenate([protein.res_name, pose.res_name]),
        res_seq=np.concatenate([protein.res_seq, pose.res_seq]),
        icode=np.concatenate([protein.icode, pose.icode]),
        res_kind=np.concatenate([protein.res_kind, pose.res_kind]),
        ss=np.concatenate([protein.ss, pose.ss]),
        seqres=np.concatenate([protein.seqres, pose.seqres]),
        label_asym=np.concatenate([protein.label_asym, pose.label_asym]),
        entity=np.concatenate([protein.entity, pose.entity]),
        chain_start=np.concatenate([protein.chain_start, pose.chain_start[1:] + r]).astype(
            np.int32
        ),
        chain_name=np.concatenate([protein.chain_name, pose.chain_name]),
        entity_id=protein.entity_id,
        entity_seq=protein.entity_seq,
        entity_seq_start=protein.entity_seq_start,
        bonds=np.concatenate([protein.bonds, pose.bonds + n]).astype(np.int32).reshape(-1, 2),
        order=np.concatenate([protein.order, pose.order]),
        next_uid=protein.next_uid + pose.next_uid,
        meta=protein.meta,
    )


def overlapping(
    protein: Structure, xyz: np.ndarray, active: np.ndarray | None = None, model: int = 0
) -> np.ndarray:
    """The atoms of the protein's ligands that come within :data:`OVERLAP` of the points
    ``xyz`` (a pose's heavy atoms): the ligands the pose stands in for."""
    if active is None:
        active = protein.active()
    ligand = (protein.res_kind[protein.residue_of] == LIGAND) & active & (protein.element != 1)
    out = np.zeros(protein.n_atoms, bool)
    candidates = np.flatnonzero(ligand)
    if not len(candidates) or not len(xyz):
        return out
    distance, _ = cKDTree(xyz).query(protein.xyz(model)[candidates], distance_upper_bound=OVERLAP)
    residues = np.unique(protein.residue_of[candidates[np.isfinite(distance)]])
    return np.isin(protein.residue_of, residues) & active


def same_pose(
    structure: Structure,
    residue: int,
    mol: Chem.Mol,
    active: np.ndarray | None = None,
    model: int = 0,
) -> bool:
    """Whether ``mol``'s heavy atoms are ``residue``'s, atom for atom (same element) within
    :data:`SAME` Å."""
    if active is None:
        active = structure.active()
    span = structure.atoms_of(residue)
    atoms = np.arange(span.start, span.stop)
    atoms = atoms[active[atoms] & (structure.element[atoms] != 1)]
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    if len(heavy) != len(atoms) or not len(atoms):
        return False
    theirs = structure.xyz(model)[atoms].astype(np.float64)
    mine = mol.GetConformer().GetPositions()[heavy]
    elements = np.array([mol.GetAtomWithIdx(i).GetAtomicNum() for i in heavy])
    distance = np.linalg.norm(mine[:, None] - theirs[None], axis=2)
    distance[elements[:, None] != structure.element[atoms][None]] = np.inf
    return bool(distance.min(axis=1).max() <= SAME and distance.min(axis=0).max() <= SAME)
