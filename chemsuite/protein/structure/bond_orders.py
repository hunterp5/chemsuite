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
"""Bond orders of a ligand from its heavy-atom geometry alone: the last resort.

Crystal structures rarely have hydrogens, and RDKit's ``DetermineBondOrders`` needs them, so
this follows the classic geometric recipe (as in Baber & Hodgkin, 1992, and Open Babel):

1. each atom's hybridisation from its bond angles, or for a terminal atom its bond length;
2. short terminal C=O, C=N and C=S bonds become double;
3. flat five- and six-membered rings of sp2 atoms become aromatic, and RDKit's kekulisation
   places their double bonds (an aromatic nitrogen is given an H when that is what it takes);
4. the remaining sp2 pairs take double bonds, shortest first;
5. phosphorus and sulfur get their P=O and S=O.

The result is a guess: callers mark it so the user can check it.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, rdBase

_SP = 1
_SP2 = 2
_SP3 = 3
#: usual valence of the atoms that take double bonds here
_VALENCE = {6: 4, 7: 3, 8: 2, 16: 2}
#: shorter than this is a double bond to a terminal atom: (element, element) -> Å
_TERMINAL_DOUBLE = {(6, 8): 1.30, (6, 7): 1.32, (6, 16): 1.72, (7, 8): 1.28}
_TRIPLE = {(6, 6): 1.25, (6, 7): 1.20}
_FLAT = 0.15  # Å from the ring's plane


def guess(mol: Chem.Mol) -> Chem.Mol:
    """``mol`` (heavy atoms, single bonds, one conformer) with bond orders; unchanged when
    no consistent assignment is found."""
    work = Chem.RWMol(mol)
    xyz = work.GetConformer().GetPositions()
    hybrid = _hybridisation(work, xyz)
    free = {a.GetIdx(): _VALENCE.get(a.GetAtomicNum(), 0) - a.GetDegree() for a in work.GetAtoms()}
    _triples(work, xyz, free)
    _terminal_doubles(work, xyz, hybrid, free)
    aromatic = _aromatic_rings(work, xyz)
    _chain_doubles(work, xyz, hybrid, free, aromatic)
    _oxo(work, xyz)
    result = _sanitised(work, aromatic)
    return result if result is not None else mol


def _length(xyz: np.ndarray, bond: Chem.Bond) -> float:
    return float(np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()]))


def _pair(bond: Chem.Bond) -> tuple[int, int]:
    a, b = bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()
    return (a, b) if a <= b else (b, a)


def _hybridisation(mol: Chem.RWMol, xyz: np.ndarray) -> dict[int, int]:
    out = {}
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(neighbours) < 2:
            out[i] = _SP3  # decided by its bond length below
            continue
        vectors = xyz[neighbours] - xyz[i]
        vectors /= np.linalg.norm(vectors, axis=1, keepdims=True)
        angles = [np.degrees(np.arccos(np.clip(vectors[a] @ vectors[b], -1, 1)))
                  for a in range(len(vectors)) for b in range(a + 1, len(vectors))]  # fmt: skip
        mean = float(np.mean(angles))
        out[i] = _SP if mean > 155 else _SP2 if mean > 115 else _SP3
    return out


def _set(mol: Chem.RWMol, bond: Chem.Bond, order: Chem.BondType, free: dict[int, int]) -> None:
    step = 2 if order == Chem.BondType.TRIPLE else 1
    bond.SetBondType(order)
    free[bond.GetBeginAtomIdx()] -= step
    free[bond.GetEndAtomIdx()] -= step


def _triples(mol: Chem.RWMol, xyz: np.ndarray, free: dict[int, int]) -> None:
    for bond in mol.GetBonds():
        limit = _TRIPLE.get(_pair(bond))
        ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if limit and _length(xyz, bond) < limit and all(free[i] >= 2 for i in ends):
            _set(mol, bond, Chem.BondType.TRIPLE, free)


def _terminal_doubles(
    mol: Chem.RWMol, xyz: np.ndarray, hybrid: dict[int, int], free: dict[int, int]
) -> None:
    candidates = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if min(a.GetDegree(), b.GetDegree()) != 1 or bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        centre = a if a.GetDegree() > 1 else b
        limit = _TERMINAL_DOUBLE.get(_pair(bond))
        if limit and _length(xyz, bond) < limit and hybrid[centre.GetIdx()] != _SP3:
            candidates.append((_length(xyz, bond), bond.GetIdx()))
    for _, index in sorted(candidates):
        bond = mol.GetBondWithIdx(index)
        if free[bond.GetBeginAtomIdx()] > 0 and free[bond.GetEndAtomIdx()] > 0:
            _set(mol, bond, Chem.BondType.DOUBLE, free)


def _aromatic_rings(mol: Chem.RWMol, xyz: np.ndarray) -> set[int]:
    atoms: set[int] = set()
    for ring in Chem.GetSymmSSSR(mol):
        ring = tuple(ring)
        if len(ring) not in (5, 6):
            continue
        points = xyz[list(ring)]
        centre = points.mean(axis=0)
        normal = np.linalg.svd(points - centre)[2][-1]
        if np.abs((points - centre) @ normal).max() >= _FLAT:
            continue
        # a flat ring is aromatic when its atoms are sp2, which shows as their substituents
        # lying in its plane (crystal angles are too noisy, and near 108 in a 5-ring)
        if not _substituents_flat(mol, xyz, ring, centre, normal):
            continue
        atoms.update(ring)
    for bond in mol.GetBonds():
        ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if all(i in atoms for i in ends) and bond.IsInRing():
            bond.SetBondType(Chem.BondType.AROMATIC)
            bond.SetIsAromatic(True)
    for i in atoms:
        mol.GetAtomWithIdx(i).SetIsAromatic(True)
    return atoms


def _substituents_flat(
    mol: Chem.RWMol, xyz: np.ndarray, ring: tuple[int, ...], centre: np.ndarray, normal: np.ndarray
) -> bool:
    for i in ring:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetDegree() > 3:
            return False
        for n in atom.GetNeighbors():
            j = n.GetIdx()
            if j not in ring and abs(float((xyz[j] - centre) @ normal)) > 4 * _FLAT:
                return False
    return True


def _chain_doubles(
    mol: Chem.RWMol,
    xyz: np.ndarray,
    hybrid: dict[int, int],
    free: dict[int, int],
    aromatic: set[int],
) -> None:
    candidates = []
    for bond in mol.GetBonds():
        ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if bond.GetBondType() != Chem.BondType.SINGLE or any(i in aromatic for i in ends):
            continue
        if all(hybrid[i] == _SP2 for i in ends):
            candidates.append((_length(xyz, bond), bond.GetIdx()))
    for _, index in sorted(candidates):
        bond = mol.GetBondWithIdx(index)
        if free[bond.GetBeginAtomIdx()] > 0 and free[bond.GetEndAtomIdx()] > 0:
            _set(mol, bond, Chem.BondType.DOUBLE, free)


def _oxo(mol: Chem.RWMol, xyz: np.ndarray) -> None:
    """One P=O on each four-coordinate phosphorus, two S=O on a four-coordinate sulfur,
    on its shortest bonds to terminal oxygens."""
    for atom in mol.GetAtoms():
        doubles = {15: 1, 16: 2}.get(atom.GetAtomicNum(), 0)
        if not doubles or atom.GetDegree() != 4:
            continue
        terminal = [mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
                    for n in atom.GetNeighbors()
                    if n.GetAtomicNum() == 8 and n.GetDegree() == 1]  # fmt: skip
        for bond in sorted(terminal, key=lambda b: _length(xyz, b))[:doubles]:
            bond.SetBondType(Chem.BondType.DOUBLE)


def _sanitised(mol: Chem.RWMol, aromatic: set[int]) -> Chem.Mol | None:
    nitrogens = [i for i in aromatic if mol.GetAtomWithIdx(i).GetAtomicNum() == 7
                 and mol.GetAtomWithIdx(i).GetDegree() == 2]  # fmt: skip
    for extra in [None, *nitrogens]:
        attempt = Chem.RWMol(mol)
        if extra is not None:
            attempt.GetAtomWithIdx(extra).SetNumExplicitHs(1)
        for atom in attempt.GetAtoms():
            atom.SetNoImplicit(False)
        try:
            with rdBase.BlockLogs():
                Chem.SanitizeMol(attempt)
        except Chem.MolSanitizeException:
            continue
        return attempt.GetMol()
    return None
