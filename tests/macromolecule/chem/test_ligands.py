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
"""Ligands as RDKit molecules with bond orders, in their bound pose."""

from __future__ import annotations

import numpy as np
from rdkit import Chem

from chemsuite.macromolecule.chem import ligands
from chemsuite.macromolecule.structure import bond_orders
from chemsuite.macromolecule.structure.bonds import template_smiles

from ..conftest import only_ligand, without_file_components


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def flat(mol: Chem.Mol) -> str:
    copy = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copy)
    return Chem.MolToSmiles(copy)


def carbon_centres(mol: Chem.Mol) -> list[str]:
    centres = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    return sorted(label for i, label in centres if mol.GetAtomWithIdx(i).GetSymbol() == "C")


def test_a_cofactor_takes_its_bond_orders_and_stereo_from_the_dictionary(kinase) -> None:
    template = Chem.MolFromSmiles(template_smiles("ATP"))
    for structure, source in ((kinase, ligands.FILE),
                              (without_file_components(kinase), ligands.TEMPLATE)):  # fmt: skip
        ligand = ligands.to_mol(structure, only_ligand(structure))
        assert ligand.source == source and not ligand.guessed
        assert flat(ligand.mol) == flat(template)
        # the ribose's stereocentres come from the pose (a phosphate's are only resonance)
        assert carbon_centres(ligand.mol) == carbon_centres(template)


def test_a_given_smiles_comes_first(lysozyme) -> None:
    smiles = "CCCc1ccccc1O"
    ligand = ligands.to_mol(lysozyme, only_ligand(lysozyme), smiles=smiles)
    assert ligand.source == ligands.SMILES
    assert Chem.MolToSmiles(ligand.mol) == canonical(smiles)


def test_a_smiles_that_does_not_fit_is_passed_over(lysozyme) -> None:
    ligand = ligands.to_mol(lysozyme, only_ligand(lysozyme), smiles="c1ccccc1")
    assert ligand.source == ligands.FILE
    bare = without_file_components(lysozyme)
    ligand = ligands.to_mol(bare, only_ligand(bare), smiles="c1ccccc1")
    assert ligand.source == ligands.GEOMETRY and ligand.guessed


def test_without_a_template_the_geometry_gives_the_bond_orders(lysozyme, thrombin) -> None:
    lysozyme = without_file_components(lysozyme)
    jz4 = ligands.to_mol(lysozyme, only_ligand(lysozyme))
    assert jz4.guessed
    assert Chem.MolToSmiles(jz4.mol) == canonical("CCCc1ccccc1O")  # 2-propylphenol
    ppack = ligands.to_mol(thrombin, only_ligand(thrombin))
    # the file's own 0G6 has the chlorine that left on binding His57, and does not flag it
    assert ppack.guessed
    assert flat(ppack.mol) == canonical("CC(O)C(CCCN=C(N)N)NC(=O)C1CCCN1C(=O)C(N)Cc1ccccc1")


def test_the_geometry_alone_recovers_atp(kinase) -> None:
    residue = only_ligand(kinase)
    bare = Chem.RemoveHs(ligands.to_mol(kinase, residue).mol)
    single = Chem.RWMol(bare)
    for bond in single.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in single.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
    guessed = bond_orders.guess(single.GetMol())
    assert flat(guessed) == flat(Chem.MolFromSmiles(template_smiles("ATP")))


def test_the_pose_atoms_and_uids_are_kept(kinase) -> None:
    residue = only_ligand(kinase)
    mol = ligands.to_mol(kinase, residue).mol
    span = kinase.atoms_of(residue)
    heavy = [i for i in range(span.start, span.stop) if kinase.element[i] != 1]
    assert mol.GetNumAtoms() == len(heavy)
    uids = [atom.GetIntProp("uid") for atom in mol.GetAtoms()]
    assert uids == kinase.uid[heavy].tolist()
    assert np.allclose(mol.GetConformer().GetPositions(), kinase.coords[0][heavy], atol=1e-4)
    info = mol.GetAtomWithIdx(0).GetMonomerInfo()
    assert (info.GetResidueName(), info.GetChainId(), info.GetResidueNumber()) == ("ATP", "E", 355)
    assert mol.GetProp("_Name") == "E/ATP 355"


def test_the_prepared_form_has_the_structures_own_hydrogens_and_charges(kinase) -> None:
    from chemsuite.macromolecule.prepare import pipeline

    s = pipeline.prepare(kinase, pipeline.Options(waters=pipeline.NONE)).structure
    residue = only_ligand(s)
    mol = ligands.prepared_form(s, residue).mol
    span = s.atoms_of(residue)
    atoms = list(range(span.start, span.stop))
    assert mol.GetNumAtoms() == len(atoms)  # every hydrogen is an atom, none added
    by_uid = {atom.GetIntProp("uid"): atom for atom in mol.GetAtoms()}
    assert set(by_uid) == set(s.uid[atoms].tolist())
    positions = mol.GetConformer().GetPositions()
    for i in atoms:
        atom = by_uid[int(s.uid[i])]
        assert atom.GetAtomicNum() == s.element[i] and atom.GetFormalCharge() == s.charge[i]
        assert np.allclose(positions[atom.GetIdx()], s.xyz()[i], atol=1e-4)
        assert atom.GetMonomerInfo().GetName().strip() == s.name[i]
    assert sum(a.GetNumImplicitHs() for a in mol.GetAtoms()) == 0
    assert Chem.GetFormalCharge(mol) == int(s.charge[atoms].sum()) < 0  # ATP's phosphates


def test_a_ligand_bonded_to_the_protein_has_no_form_of_its_own(egfr) -> None:
    import pytest

    residue = [r for r in ligands.ligand_residues(egfr).tolist() if egfr.res_name[r] == "80U"][0]
    with pytest.raises(ligands.FormError, match="CYS 797"):
        ligands.prepared_form(egfr, residue)


def test_a_hydrogen_short_is_an_error_not_a_radical(kinase) -> None:
    import pytest

    from chemsuite.macromolecule.prepare import pipeline
    from chemsuite.macromolecule.structure import edit

    s = pipeline.prepare(kinase, pipeline.Options(waters=pipeline.NONE)).structure
    residue = only_ligand(s)
    span = s.atoms_of(residue)
    one_h = [i for i in range(span.start, span.stop) if s.element[i] == 1][0]
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[one_h] = True
    with pytest.raises(ligands.FormError, match="missing a hydrogen"):
        ligands.prepared_form(edit.delete(s, gone), residue)
