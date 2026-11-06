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

"""The sketch to and from RDKit, molfiles and SMILES.

RDKit always reads a sketch the way it reads any molfile: the sketch is written as a V2000
block exactly as drawn (no hydrogen counts, wedges as drawn) and parsed with
``Chem.MolFromMolBlock``. A drawing is therefore valid here exactly when the table, search and
descriptors can use it (docs/sketcher/valence-review.md), and double-bond geometry and wedges
are interpreted by RDKit's molfile reader, not by code of our own.
"""

from __future__ import annotations

import math
from statistics import median

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor

from chemsuite.core.mols import mol_from_text
from chemsuite.sketch import ctfile
from chemsuite.sketch.model import Atom, Bond, Order, Sketch, Stereo

#: Molfile coordinates per sketch bond length. RDKit draws bonds 1.5 long; files keep that.
FILE_BOND_LENGTH = 1.5
MAX_V2000_ATOMS = 999


class SketchError(ValueError):
    """The drawing is not a molecule RDKit accepts."""


_RADICAL_CODES = {1: 2, 2: 3}  # unpaired electrons -> molfile RAD (doublet, triplet)
_STEREO_CODES = {Stereo.NONE: 0, Stereo.WEDGE: 1, Stereo.HASH: 6, Stereo.EITHER: 4}
_BOND_DIRS = {
    Chem.BondDir.BEGINWEDGE: Stereo.WEDGE,
    Chem.BondDir.BEGINDASH: Stereo.HASH,
    Chem.BondDir.UNKNOWN: Stereo.EITHER,
}


def atom_order(sketch: Sketch) -> list[int]:
    """Sketch atom ids in molfile order (atom *i* in RDKit is ``atom_order(sketch)[i]``)."""
    return sorted(sketch.atoms)


def drawn_molblock(sketch: Sketch) -> str:
    """The sketch as a V2000 molfile, exactly as drawn: no hydrogen counts or valences.

    Raises SketchError for more than 999 atoms or bonds, the V2000 limit.
    """
    order = atom_order(sketch)
    index = {aid: i + 1 for i, aid in enumerate(order)}
    bonds = sorted(sketch.bonds.items())
    if len(order) > MAX_V2000_ATOMS or len(bonds) > MAX_V2000_ATOMS:
        raise SketchError("The drawing has more than 999 atoms or bonds.")
    lines = ["", "  chemsuite sketch", ""]
    lines.append(f"{len(order):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000")
    charges, isotopes, radicals = [], [], []
    for aid in order:
        atom = sketch.atoms[aid]
        x, y = atom.x * FILE_BOND_LENGTH, atom.y * FILE_BOND_LENGTH
        lines.append(
            f"{x:10.4f}{y:10.4f}{0:10.4f} {atom.symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
        if atom.charge:
            charges.append((index[aid], atom.charge))
        if atom.isotope:
            isotopes.append((index[aid], atom.isotope))
        if atom.radical:
            radicals.append((index[aid], _RADICAL_CODES.get(atom.radical, 3)))
    for _, bond in bonds:
        stereo = _STEREO_CODES[bond.stereo] if bond.order is Order.SINGLE else 0
        lines.append(f"{index[bond.begin]:3d}{index[bond.end]:3d}{bond.order.value:3d}{stereo:3d}")
    for tag, values in (("CHG", charges), ("ISO", isotopes), ("RAD", radicals)):
        for start in range(0, len(values), 8):  # at most eight entries per line
            chunk = values[start : start + 8]
            body = "".join(f" {i:3d} {v:3d}" for i, v in chunk)
            lines.append(f"M  {tag}{len(chunk):3d}{body}")
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def to_rdkit(sketch: Sketch, *, sanitize: bool = True) -> Chem.Mol:
    """The sketch read by RDKit, drawn hydrogens kept, atoms in ``atom_order``.

    Raises SketchError when RDKit rejects it (sanitisation failed).
    """
    with rdBase.BlockLogs():
        mol = Chem.MolFromMolBlock(drawn_molblock(sketch), sanitize=sanitize, removeHs=False)
    if mol is None:
        raise SketchError("RDKit cannot read this drawing.")
    return mol


def to_mol(sketch: Sketch) -> Chem.Mol:
    """The molecule as the table would store it: sanitised, drawn hydrogens removed."""
    mol = to_rdkit(sketch)
    with rdBase.BlockLogs():
        return Chem.RemoveHs(mol, updateExplicitCount=True)


def to_smiles(sketch: Sketch) -> str:
    """Canonical isomeric SMILES; raises SketchError for an invalid drawing."""
    if sketch.is_empty():
        return ""
    return Chem.MolToSmiles(to_mol(sketch))


def to_molblock(sketch: Sketch) -> str:
    """A molfile for other programs, with the wedges as drawn.

    Bonds are Kekulé. Any atom whose hydrogen count a CTfile reader would work out differently
    from RDKit carries an explicit valence (for example 6 on the phosphorus of PF₆⁻), so other
    programs read the same structure. A drawing RDKit rejects is written exactly as drawn, so a
    copy never fails.
    """
    try:
        mol = to_rdkit(sketch)
    except SketchError:
        return drawn_molblock(sketch)
    Chem.ReapplyMolBlockWedging(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    if mol.GetNumAtoms() > MAX_V2000_ATOMS or mol.GetNumBonds() > MAX_V2000_ATOMS:
        return Chem.MolToMolBlock(mol, kekulize=False, forceV3000=True)
    lines = Chem.MolToMolBlock(mol, kekulize=False).split("\n")
    for atom in mol.GetAtoms():
        valence = _ctfile_valence(atom)
        if valence is not None:
            line = lines[4 + atom.GetIdx()]
            lines[4 + atom.GetIdx()] = f"{line[:48]}{valence:3d}{line[51:]}"
    return "\n".join(lines)


def _ctfile_valence(atom: Chem.Atom) -> int | None:
    """The molfile valence field this atom needs, or None if CTfile readers agree with RDKit.

    The field holds the total valence (bond orders plus hydrogens); 15 means zero.
    """
    bonds = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
    hydrogens = atom.GetTotalNumHs()
    symbol = atom.GetSymbol()
    if symbol in ctfile.ELEMENTS:
        theirs = ctfile.default_hydrogens(
            symbol, atom.GetFormalCharge(), bonds, atom.GetNumRadicalElectrons()
        )
        if theirs == hydrogens:
            return None
    total = bonds + hydrogens
    return total if total else 15


def from_mol(mol: Chem.Mol) -> Sketch:
    """A sketch of a molecule: its 2D coordinates (computed if it has none), Kekulé bonds, and
    wedges as in its molfile, or chosen by RDKit."""
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0 or mol.GetConformer().Is3D():
        rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    if any(b.HasProp("_MolFileBondStereo") for b in mol.GetBonds()):
        Chem.ReapplyMolBlockWedging(mol)
    else:
        Chem.WedgeMolBonds(mol, conf)
    try:
        with rdBase.BlockLogs():
            Chem.Kekulize(mol, clearAromaticFlags=True)
    except (ValueError, RuntimeError):
        pass  # keep aromatic bonds that cannot be made alternating
    lengths = [
        (
            conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())
        ).Length()
        for b in mol.GetBonds()
    ]
    scale = median(lengths) if lengths and median(lengths) > 1e-6 else FILE_BOND_LENGTH
    sketch = Sketch()
    ids: list[int] = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        aid = sketch.new_id()
        sketch.atoms[aid] = Atom(
            p.x / scale,
            p.y / scale,
            symbol=atom.GetSymbol(),
            charge=atom.GetFormalCharge(),
            isotope=atom.GetIsotope(),
            radical=atom.GetNumRadicalElectrons(),
        )
        ids.append(aid)
    for bond in mol.GetBonds():
        kind = bond.GetBondType()
        order = {
            Chem.BondType.DOUBLE: Order.DOUBLE,
            Chem.BondType.TRIPLE: Order.TRIPLE,
            Chem.BondType.AROMATIC: Order.AROMATIC,
        }.get(kind, Order.SINGLE)
        stereo = _BOND_DIRS.get(bond.GetBondDir(), Stereo.NONE)
        if bond.HasProp("_MolFileBondStereo") and bond.GetIntProp("_MolFileBondStereo") == 4:
            stereo = Stereo.EITHER
        sketch.bonds[sketch.new_id()] = Bond(
            ids[bond.GetBeginAtomIdx()], ids[bond.GetEndAtomIdx()], order, stereo
        )
    return sketch


def from_text(text: str) -> Sketch | None:
    """A sketch from SMILES, a molfile or an InChI; None when the text is not a structure."""
    mol = mol_from_text(text)
    return None if mol is None else from_mol(mol)


def placed(sketch: Sketch, at: tuple[float, float]) -> Sketch:
    """A copy of the sketch moved so the centre of its box is at ``at``."""
    x0, y0, x1, y1 = sketch.bounds()
    dx, dy = at[0] - (x0 + x1) / 2, at[1] - (y0 + y1) / 2
    out = sketch.copy()
    out.atoms = {aid: a.moved(a.x + dx, a.y + dy) for aid, a in sketch.atoms.items()}
    return out


def bond_length_spread(sketch: Sketch) -> float:
    """Largest relative deviation of a bond length from 1 (a layout check for tests)."""
    worst = 0.0
    for bond in sketch.bonds.values():
        a, b = sketch.atoms[bond.begin], sketch.atoms[bond.end]
        worst = max(worst, abs(math.hypot(a.x - b.x, a.y - b.y) - 1.0))
    return worst
