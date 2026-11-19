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
"""PDBQT: the files AutoDock Vina reads and writes, through Meeko.

This is the only module that imports Meeko. :func:`receptor` turns a prepared structure into
the rigid receptor, :func:`ligand` a molecule with 3D coordinates into a flexible ligand, and
:func:`poses` reads Vina's output back as molecules whose bond orders survive (Meeko writes the
ligand's SMILES into the file and reads it back from there).

Meeko types a receptor residue by residue from templates it bundles, and **downloads the
component from the wwPDB for any residue name it does not know**. chemsuite never lets it:
:func:`receptor` keeps only the residues :func:`known_residues` lists, and names the rest in
:attr:`Receptor.dropped`, so the caller can say what was left out. Ligands, additives and
anything unrecognised are no part of a receptor anyway — a ligand is what is being docked.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

import numpy as np
from rdkit import Chem

from ..chem import poses as pose_tools
from ..structure import io
from ..structure.model import AMINO, GLYCAN, ION, NUCLEIC, WATER, Structure

#: the residue kinds a receptor is built from; a ligand is docked, not part of the receptor
RECEPTOR_KINDS = (AMINO, NUCLEIC, WATER, ION, GLYCAN)
#: what Meeko raises when it will not have the input (it has no exceptions of its own for
#: molecules: a ligand with no hydrogens or no conformer is a ``ValueError``, a PDBQT it cannot
#: read a ``RuntimeError``), plus what a molecule that breaks its own assumptions reaches us as
MEEKO_ERRORS = (ValueError, KeyError, IndexError, TypeError, AttributeError, RuntimeError)


class PdbqtError(Exception):
    """A structure or molecule Meeko could not turn into PDBQT, or an output it could not read."""


@dataclass(frozen=True)
class Receptor:
    """The rigid receptor as Vina reads it."""

    pdbqt: str
    #: residue names left out, because Meeko has no template for them
    dropped: tuple[str, ...] = ()


@cache
def known_residues() -> frozenset[str]:
    """Every residue name Meeko can type from the templates it bundles (nothing downloaded)."""
    from meeko import ResidueChemTemplates

    templates = ResidueChemTemplates.create_from_defaults()
    names = {str(name).split("_")[0] for name in templates.residue_templates}
    return frozenset(names | {str(name) for name in templates.ambiguous})


def receptor(s: Structure, keep: np.ndarray | None = None, model: int = 0) -> Receptor:
    """*s* (prepared, with hydrogens) as the rigid receptor.

    *keep* is an atom mask narrowing what goes in; of it, the residue kinds in
    :data:`RECEPTOR_KINDS` whose names :func:`known_residues` knows are kept.
    """
    from meeko import PDBQTWriterLegacy, Polymer, PolymerCreationError, ResidueChemTemplates

    known = known_residues()
    names = [str(n).strip().upper() for n in s.res_name]
    typed = np.isin(s.res_kind, RECEPTOR_KINDS) & np.array([n in known for n in names], bool)
    wanted = np.ones(s.n_residues, bool) if keep is None else _residues_in(s, keep)
    chosen = typed & wanted
    if not chosen.any():
        raise PdbqtError("There is nothing in this protein that Vina can use as a receptor.")
    dropped = {name for name, out in zip(names, wanted & ~typed, strict=True) if out}
    atoms = chosen[s.residue_of]
    if keep is not None:
        atoms = atoms & np.asarray(keep, bool)
    one = s.with_coords(s.xyz(model)) if s.n_models > 1 else s
    try:
        polymer = Polymer.from_pdb_string(
            io.to_text(one, "pdb", keep=atoms), ResidueChemTemplates.create_from_defaults()
        )
        rigid, _flexible = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    except (PolymerCreationError, ValueError, KeyError, RuntimeError) as err:
        raise PdbqtError(f"Meeko could not type this receptor: {err}") from err
    if not rigid.strip():
        raise PdbqtError("Meeko wrote an empty receptor.")
    missed = set(polymer.monomers) - set(polymer.get_valid_monomers())
    return Receptor(rigid, tuple(sorted(dropped | {str(m) for m in missed})))


def _residues_in(s: Structure, keep: np.ndarray) -> np.ndarray:
    """The residues with at least one atom in the mask *keep*."""
    chosen = np.zeros(s.n_residues, bool)
    chosen[np.unique(s.residue_of[np.asarray(keep, bool)])] = True
    return chosen


def ligand(mol: Chem.Mol) -> str:
    """*mol*, with 3D coordinates, as a flexible ligand (Meeko finds its rotatable bonds)."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    posed = pose_tools.pose_of(mol)
    if posed is None:
        raise PdbqtError("A ligand needs 3D coordinates before it can be docked.")
    if not any(atom.GetAtomicNum() == 1 for atom in posed.GetAtoms()):
        posed = Chem.AddHs(posed, addCoords=True)
    try:
        setups = MoleculePreparation().prepare(posed)
    except MEEKO_ERRORS as err:
        raise PdbqtError(f"Meeko could not prepare this ligand: {err}") from err
    if not setups:
        raise PdbqtError("Meeko made nothing of this ligand.")
    text, ok, why = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise PdbqtError(f"Meeko could not write this ligand: {str(why).strip()}")
    return text


def poses(text: str) -> list[Chem.Mol]:
    """Vina's output PDBQT as one molecule per pose, bond orders and hydrogens restored, and
    stereo read from the pose's own geometry (a ligand docked with a center or double bond left
    undefined is still one isomer in each pose, as a pose read from a file would say)."""
    from meeko import PDBQTMolecule, RDKitMolCreate

    try:
        parsed = PDBQTMolecule(text, skip_typing=True)
        mols = [m for m in RDKitMolCreate.from_pdbqt_mol(parsed) if m is not None]
    except MEEKO_ERRORS as err:
        raise PdbqtError(f"Could not read the docked poses: {err}") from err
    out = []
    for mol in mols:
        for conformer in mol.GetConformers():
            one = Chem.Mol(mol)
            one.RemoveAllConformers()
            one.AddConformer(Chem.Conformer(conformer), assignId=True)
            Chem.AssignStereochemistryFrom3D(one)
            out.append(one)
    if not out:
        raise PdbqtError("The docking output held no poses.")
    return out
