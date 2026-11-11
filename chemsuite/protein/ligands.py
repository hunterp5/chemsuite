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
"""Ligands as RDKit molecules in their bound pose, with real bond orders.

A structure file gives a ligand's atoms and positions but not its bond orders. They come from
the first of these that fits the atoms (``AssignBondOrdersFromTemplate``):

1. a SMILES the caller has: one the user gave, or the SMILES of the ligand's table row;
2. the Chemical Component Dictionary's SMILES for the residue name (the bundled templates
   now; the online dictionary later, when the network is allowed);
3. the geometry (:mod:`.bond_orders`), marked as guessed.

Heavy atoms only: hydrogens are placed when a tool needs them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from . import bond_orders
from .bonds import template_smiles
from .model import LIGAND, Structure

SMILES = "smiles"
TEMPLATE = "template"
GEOMETRY = "geometry"


@dataclass(frozen=True)
class Ligand:
    mol: Chem.Mol
    residue: int  # residue index in the structure
    source: str  # where the bond orders came from: SMILES, TEMPLATE or GEOMETRY

    @property
    def guessed(self) -> bool:
        return self.source == GEOMETRY


def ligand_residues(structure: Structure) -> np.ndarray:
    return np.flatnonzero(structure.res_kind == LIGAND)


def to_mol(
    structure: Structure,
    residue: int,
    smiles: str | None = None,
    active: np.ndarray | None = None,
    model: int = 0,
) -> Ligand:
    """Residue ``residue`` as an RDKit molecule (heavy atoms) in its pose."""
    if active is None:
        active = structure.active()
    span = structure.atoms_of(residue)
    atoms = [i for i in range(span.start, span.stop)
             if active[i] and structure.element[i] != 1]  # fmt: skip
    bare = _connectivity(structure, residue, atoms, model)
    name = str(structure.res_name[residue])
    for source, text in ((SMILES, smiles), (TEMPLATE, template_smiles(name))):
        if text:
            fitted = _from_template(bare, text)
            if fitted is not None:
                return Ligand(fitted, residue, source)
    return Ligand(bond_orders.guess(bare), residue, GEOMETRY)


def _connectivity(structure: Structure, residue: int, atoms: list[int], model: int) -> Chem.Mol:
    mol = Chem.RWMol()
    index = {}
    chain = str(structure.chain_name[structure.chain_of[residue]])
    for i in atoms:
        atom = Chem.Atom(int(structure.element[i]))
        atom.SetFormalCharge(int(structure.charge[i]))
        atom.SetNoImplicit(True)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f" {structure.name[i]:<3}"[:4])
        info.SetResidueName(str(structure.res_name[residue]))
        info.SetResidueNumber(int(structure.res_seq[residue]))
        info.SetInsertionCode(structure.icode[residue] or " ")
        info.SetChainId(chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        atom.SetIntProp("uid", int(structure.uid[i]))
        index[i] = mol.AddAtom(atom)
    for a, b in structure.bonds:
        if int(a) in index and int(b) in index:
            mol.AddBond(index[int(a)], index[int(b)], Chem.BondType.SINGLE)
    conformer = Chem.Conformer(len(atoms))
    xyz = structure.xyz(model)
    for i, atom in enumerate(atoms):
        conformer.SetAtomPosition(i, xyz[atom].astype(float).tolist())
    mol.AddConformer(conformer, assignId=True)
    out = mol.GetMol()
    out.SetProp("_Name", structure.residue_label(residue))
    return out


def _from_template(bare: Chem.Mol, smiles: str) -> Chem.Mol | None:
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        return None
    template = Chem.RemoveHs(template)
    if template.GetNumAtoms() != bare.GetNumAtoms():
        return None
    try:
        with rdBase.BlockLogs():
            fitted = AllChem.AssignBondOrdersFromTemplate(template, bare)
    except (ValueError, Chem.MolSanitizeException):
        return None
    fitted = Chem.RWMol(fitted)
    for atom in fitted.GetAtoms():
        atom.SetNoImplicit(False)  # the template's hydrogens become implicit ones
        atom.SetNumRadicalElectrons(0)
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(fitted)
    except Chem.MolSanitizeException:
        return None
    out = fitted.GetMol()
    Chem.AssignStereochemistryFrom3D(out)
    return out


# -- ligands as table rows ----------------------------------------------------------------------

#: the columns a protein's ligands fill: which protein, which residue, the ligand in its pose
PROTEIN_COLUMN = "Protein"
LIGAND_ID_COLUMN = "Ligand ID"
LIGAND_COLUMN = "Ligand"
BOND_ORDERS_COLUMN = "Bond Orders"
SOURCE_TITLES = {SMILES: "SMILES", TEMPLATE: "CCD template", GEOMETRY: "Guessed from geometry"}


def ligand_rows(
    structure: Structure, protein: str, residues: list[int] | None = None
) -> dict[str, list[str]]:
    """One row per ligand (or per residue of ``residues``): the protein's asset id, the residue
    (as :meth:`Structure.residue_label` writes it), a molfile in the ligand's pose and where its
    bond orders came from. A ligand RDKit cannot read is left out."""
    rows: dict[str, list[str]] = {PROTEIN_COLUMN: [], LIGAND_ID_COLUMN: [], LIGAND_COLUMN: [],
                                  BOND_ORDERS_COLUMN: []}  # fmt: skip
    if residues is None:
        residues = ligand_residues(structure).tolist()
    for residue in residues:
        try:
            ligand = to_mol(structure, residue)
            block = Chem.MolToMolBlock(ligand.mol)
        except (ValueError, RuntimeError):
            continue
        rows[PROTEIN_COLUMN].append(protein)
        rows[LIGAND_ID_COLUMN].append(structure.residue_label(residue))
        rows[LIGAND_COLUMN].append(block)
        rows[BOND_ORDERS_COLUMN].append(SOURCE_TITLES[ligand.source])
    return rows


def residue_named(structure: Structure, label: str) -> int:
    """The residue :meth:`Structure.residue_label` calls ``label`` (-1 when there is none)."""
    for residue in ligand_residues(structure).tolist():
        if structure.residue_label(residue) == label:
            return int(residue)
    for residue in range(structure.n_residues):
        if structure.residue_label(residue) == label:
            return residue
    return -1
