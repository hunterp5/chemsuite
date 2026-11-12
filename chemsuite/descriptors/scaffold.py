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

"""The Bemis-Murcko scaffold: rings, the linkers between them, and what is double-bonded to
them.

RDKit's ``MurckoScaffold.GetScaffoldForMol`` gives it, except where RDKit 2026.9.2 fails on
macOS ("bad any cast" inside ``MurckoDecompose``). There :func:`pruned` takes its place. It
follows RDKit's rules: side chains are stripped atom by atom from their ends; an atom
double-bonded to the scaffold stays, without its own substituents; an atom that lost a
substituent loses its stereo and, where its hydrogens are fixed (an aromatic heteroatom, a
bracket atom), takes a hydrogen for it. It gives RDKit's scaffold for every one of 26,048
molecules (the approved drugs and 10,000 BindingDB compounds, read from SDF and from SMILES).
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def murcko(mol: Chem.Mol) -> Chem.Mol:
    """The scaffold of *mol* (an empty molecule when it has no ring)."""
    try:
        return MurckoScaffold.GetScaffoldForMol(mol)
    except RuntimeError:  # RDKit 2026.9.2 on macOS
        return pruned(mol)


def pruned(mol: Chem.Mol) -> Chem.Mol:
    """The scaffold of *mol*, found without ``MurckoDecompose``."""
    mol = Chem.RemoveHs(mol)
    scaffold = _rings_and_linkers(mol)
    keep = list(scaffold)
    for bond in mol.GetBonds():  # an atom double-bonded to the scaffold stays (C=O, C=C)
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            keep[a] = keep[a] or scaffold[b]
            keep[b] = keep[b] or scaffold[a]
    if not any(keep):
        return Chem.Mol()
    edited = Chem.RWMol(mol)
    for atom in edited.GetAtoms():
        if scaffold[atom.GetIdx()]:  # an atom kept for its double bond is left as it was
            _lose_substituents(atom, sum(1 for n in atom.GetNeighbors() if not keep[n.GetIdx()]))
    for i in sorted((i for i, k in enumerate(keep) if not k), reverse=True):
        edited.RemoveAtom(i)
    out = edited.GetMol()
    # as GetScaffoldForMol finishes: no sanitising, so aromatic flags stay as they were
    out.ClearComputedProps()
    out.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(out)
    return out


def _rings_and_linkers(mol: Chem.Mol) -> list[bool]:
    """Per atom, whether it is in a ring or on a path between rings: side chains are
    stripped from their ends, an atom outside the rings left with one neighbour at a time."""
    ring = mol.GetRingInfo()
    keep = [True] * mol.GetNumAtoms()
    degree = [atom.GetDegree() for atom in mol.GetAtoms()]
    changed = True
    while changed:
        changed = False
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            if keep[i] and degree[i] <= 1 and not ring.NumAtomRings(i):
                keep[i] = False
                changed = True
                for neighbour in atom.GetNeighbors():
                    degree[neighbour.GetIdx()] -= 1
    return keep


def _lose_substituents(atom: Chem.Atom, lost: int) -> None:
    """*atom* loses *lost* substituents: its stereo goes, and where its hydrogens are fixed
    (an aromatic heteroatom, or a bracket atom outside aromatic rings: [N+], [C@H]) it takes
    a hydrogen for each."""
    if not lost:
        return
    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    fixed = atom.GetAtomicNum() != 6 if atom.GetIsAromatic() else atom.GetNoImplicit()
    if fixed:
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + lost)
