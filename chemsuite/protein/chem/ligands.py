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
the first of these that fits the atoms:

1. a SMILES the caller has: one the user gave, or the SMILES of the ligand's table row;
2. the file's own definition of the component (its ``_chem_comp_bond`` table);
3. the Chemical Component Dictionary online
   (:class:`~chemsuite.protein.structure.components.CcdCache`), when the caller allows a
   download;
4. the templates bundled with chemsuite (common cofactors and additives);
5. the geometry (:mod:`chemsuite.protein.structure.bond_orders`), marked as guessed.

A component fits by atom name when the residue has exactly its heavy atoms, which also takes the
bonds from the dictionary; otherwise by matching the graph (``AssignBondOrdersFromTemplate``).

Heavy atoms only: hydrogens are placed when a tool needs them. Ligands become table rows
through the bridge (:mod:`chemsuite.tools.protein.rows`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from chemsuite.core.jobs import CancelToken

from ..structure import bond_orders, components
from ..structure.bonds import template_smiles
from ..structure.components import CcdCache, Component
from ..structure.model import LIGAND, Structure

SMILES = "smiles"
FILE = "file"
CCD = "ccd"
TEMPLATE = "template"
GEOMETRY = "geometry"
#: cofactors: a protein's ligand rows leave them out when it has another ligand
COFACTORS = frozenset(["ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "GTP", "GDP", "GNP", "GSP",
                       "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "HEM", "HEC",
                       "COA", "ACO", "PLP", "TPP", "BTN", "UDP", "CTP", "UTP", "TTP"])  # fmt: skip


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
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> Ligand:
    """Residue ``residue`` as an RDKit molecule (heavy atoms) in its pose. ``ccd`` allows
    asking the online dictionary."""
    if active is None:
        active = structure.active()
    span = structure.atoms_of(residue)
    atoms = [i for i in range(span.start, span.stop)
             if active[i] and structure.element[i] != 1]  # fmt: skip
    bare = _connectivity(structure, residue, atoms, model)
    name = str(structure.res_name[residue])
    if smiles:
        fitted = _from_smiles(bare, smiles)
        if fitted is not None:
            return Ligand(fitted, residue, SMILES)
    own = components.from_meta(structure.meta.get(components.META_KEY, "")).get(name)
    if own is not None:
        fitted = _from_component(bare, own)
        if fitted is not None:
            return Ligand(fitted, residue, FILE)
    if ccd is not None:
        online = ccd.get(name, cancel)
        if online is not None:
            fitted = _from_component(bare, online)
            if fitted is not None:
                return Ligand(fitted, residue, CCD)
    fitted = _from_smiles(bare, template_smiles(name)) if template_smiles(name) else None
    if fitted is not None:
        return Ligand(fitted, residue, TEMPLATE)
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


def _from_smiles(bare: Chem.Mol, smiles: str) -> Chem.Mol | None:
    with rdBase.BlockLogs():
        template = Chem.MolFromSmiles(smiles)
    return None if template is None else _from_template(bare, Chem.RemoveHs(template))


def _from_component(bare: Chem.Mol, component: Component) -> Chem.Mol | None:
    names = [atom.GetMonomerInfo().GetName().strip() for atom in bare.GetAtoms()]
    defined = {name: (element, charge) for name, element, charge in component.atoms}
    absent = defined.keys() - set(names)
    if (len(set(names)) == len(names) and absent <= component.leaving
            and all(defined.get(n, ("", 0))[0] == atom.GetSymbol().upper()
                    for n, atom in zip(names, bare.GetAtoms(), strict=True))):  # fmt: skip
        fitted = _by_name(bare, component, names, defined)
        if fitted is not None:
            return fitted
    template = component.to_mol()
    return None if template is None else _from_template(bare, template)


def _by_name(
    bare: Chem.Mol, component: Component, names: list[str], defined: dict[str, tuple[str, int]]
) -> Chem.Mol | None:
    """The dictionary's bonds, orders and charges on the residue's atoms, matched by name."""
    mol = Chem.RWMol(bare)
    for bond in list(mol.GetBonds()):
        mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    index = {name: i for i, name in enumerate(names)}
    for a, b, order, _aromatic in component.bonds:
        if a in index and b in index:  # not a leaving atom
            mol.AddBond(index[a], index[b], components.bond_type(order))
    for name, i in index.items():
        mol.GetAtomWithIdx(i).SetFormalCharge(defined[name][1])
    return _finish(mol)


def _from_template(bare: Chem.Mol, template: Chem.Mol) -> Chem.Mol | None:
    if template.GetNumAtoms() != bare.GetNumAtoms():
        return None
    try:
        with rdBase.BlockLogs():
            fitted = AllChem.AssignBondOrdersFromTemplate(template, bare)
    except (ValueError, Chem.MolSanitizeException):
        return None
    return _finish(Chem.RWMol(fitted))


def _finish(fitted: Chem.RWMol) -> Chem.Mol | None:
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


# -- which ligands ---------------------------------------------------------------------------------

#: a protein asset's meta: the ligand IDs whose bond orders are a guess, one per line
GUESSED = "guessed"


def row_residues(structure: Structure) -> list[int]:
    """The ligands that become rows: all of them, less the cofactors when there is another."""
    found = ligand_residues(structure).tolist()
    others = [r for r in found if str(structure.res_name[r]) not in COFACTORS]
    return others or found


def residue_named(structure: Structure, label: str) -> int:
    """The residue :meth:`Structure.residue_label` calls ``label`` (-1 when there is none)."""
    for residue in ligand_residues(structure).tolist():
        if structure.residue_label(residue) == label:
            return int(residue)
    for residue in range(structure.n_residues):
        if structure.residue_label(residue) == label:
            return residue
    return -1
