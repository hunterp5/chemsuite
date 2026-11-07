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

"""2D structure pictures for table cells, tooltips and exports (Qt-free).

A molecule keeps the 2D coordinates it was drawn with (e.g. from ChemDraw via SDF) so chemists
see their own orientation. Molecules with no coordinates, or with 3D coordinates, get a fresh
2D layout (CoordGen, which lays out macrocycles and bridged rings far better than the default).
"""

from __future__ import annotations

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

#: Bond length in device pixels at scale 1: keeps small molecules from being blown up so a
#: column of compounds reads at one consistent size; large ones still shrink to fit.
_BOND_PX = 22.0


def depiction_mol(mol: Chem.Mol) -> Chem.Mol:
    """*mol* itself if it has 2D coordinates, else a copy with a computed 2D layout.

    3D structures (conformers) are drawn without their explicit hydrogens, as chemists draw.
    """
    if mol.GetNumConformers() and not mol.GetConformer().Is3D():
        return mol
    with rdBase.BlockLogs():
        copy = Chem.RemoveHs(mol) if mol.GetNumConformers() else Chem.Mol(mol)
    copy.RemoveAllConformers()
    with rdBase.BlockLogs():
        rdDepictor.Compute2DCoords(copy)
    return copy


def fitted_size(
    mol: Chem.Mol, max_width: int, max_height: int, *, min_width: int = 56, min_height: int = 40
) -> tuple[int, int]:
    """The box (in pixels at scale 1) that holds *mol* drawn at the normal bond length, no
    larger than ``max_width x max_height``: a small molecule gets a small box instead of a
    large empty one, and a large one fills the box (shrinking only as much as it must)."""
    shown = depiction_mol(mol)
    if not shown.GetNumConformers() or shown.GetNumAtoms() == 0:
        return max_width, max_height
    positions = shown.GetConformer().GetPositions()[:, :2]
    lengths = [
        float(((positions[b.GetBeginAtomIdx()] - positions[b.GetEndAtomIdx()]) ** 2).sum()) ** 0.5
        for b in shown.GetBonds()
    ]
    bond = sum(lengths) / len(lengths) if lengths else 1.5
    extent = (positions.max(axis=0) - positions.min(axis=0)) / max(bond, 1e-6) * _BOND_PX
    margin = 1.6 * _BOND_PX  # room for atom labels and padding on both sides
    width = min(max_width, max(min_width, round(extent[0] + margin)))
    height = min(max_height, max(min_height, round(extent[1] + margin)))
    return int(width), int(height)


def render_png(
    mol: Chem.Mol,
    width: int,
    height: int,
    *,
    scale: float = 1.0,
    dark: bool = False,
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
) -> bytes:
    """PNG bytes of *mol* centred in a transparent ``width x height`` pixel box.

    *scale* is the device pixel ratio (2 on Retina/HiDPI screens), so lines and labels stay the
    same physical size on every display. Highlight indices refer to :func:`depiction_mol`'s
    atoms and bonds (the same as *mol*'s unless *mol* is 3D).
    """
    drawer = rdMolDraw2D.MolDraw2DCairo(max(int(width), 1), max(int(height), 1))
    options = drawer.drawOptions()
    if dark:
        rdMolDraw2D.SetDarkMode(options)
    options.clearBackground = False
    options.padding = 0.06
    options.fixedBondLength = _BOND_PX * scale
    options.bondLineWidth = max(1.0, 1.4 * scale)
    options.minFontSize = int(9 * scale)
    options.maxFontSize = int(15 * scale)
    options.addStereoAnnotation = False
    with rdBase.BlockLogs():
        try:
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer,
                depiction_mol(mol),
                highlightAtoms=highlight_atoms or [],
                highlightBonds=highlight_bonds or [],
            )
        except (RuntimeError, ValueError):  # e.g. a structure that cannot be kekulized
            drawer = rdMolDraw2D.MolDraw2DCairo(max(int(width), 1), max(int(height), 1))
            drawer.drawOptions().clearBackground = False
            drawer.drawOptions().prepareMolsBeforeDrawing = False
            drawer.DrawMolecule(depiction_mol(mol))
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def render_reaction_png(
    rxn: rdChemReactions.ChemicalReaction,
    width: int,
    height: int,
    *,
    scale: float = 1.0,
    dark: bool = False,
) -> bytes | None:
    """PNG bytes of *rxn* (reactants, arrow with agents above it, products) in a transparent
    ``width x height`` box, drawn as :func:`render_png` draws molecules; None if RDKit cannot
    draw it. Atom map numbers are left out: in a cell they are clutter, and the sketcher shows
    them."""
    copy = rdChemReactions.ChemicalReaction(rxn)
    for template in [*copy.GetReactants(), *copy.GetAgents(), *copy.GetProducts()]:
        for atom in template.GetAtoms():
            atom.SetAtomMapNum(0)
    drawer = rdMolDraw2D.MolDraw2DCairo(max(int(width), 1), max(int(height), 1))
    options = drawer.drawOptions()
    if dark:
        rdMolDraw2D.SetDarkMode(options)
    options.clearBackground = False
    options.padding = 0.04
    options.fixedBondLength = _BOND_PX * scale
    options.bondLineWidth = max(1.0, 1.4 * scale)
    options.minFontSize = int(9 * scale)
    options.maxFontSize = int(15 * scale)
    options.addStereoAnnotation = False
    with rdBase.BlockLogs():
        try:
            drawer.DrawReaction(copy)
        except (RuntimeError, ValueError):
            return None
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def alignment_core(text: str) -> Chem.Mol | None:
    """The core that drawings are aligned to, from SMILES (or SMARTS), with 2D coordinates.

    None when *text* reads as neither.
    """
    text = text.strip()
    if not text:
        return None
    with rdBase.BlockLogs():
        core = Chem.MolFromSmiles(text) or Chem.MolFromSmarts(text)
        if core is None:
            return None
        core.UpdatePropertyCache(strict=False)
        try:
            rdDepictor.Compute2DCoords(core)
        except (RuntimeError, ValueError):
            return None
    return core


def aligned_mol(mol: Chem.Mol, core: Chem.Mol) -> Chem.Mol:
    """A copy of *mol* laid out so the atoms matching *core* sit where they do in *core*.

    Molecules without the core keep their usual drawing, so a column can be aligned even when
    only some rows share it. The molecule itself (and its saved coordinates) is not changed.
    """
    shown = depiction_mol(mol)
    with rdBase.BlockLogs():
        if not shown.HasSubstructMatch(core):
            return shown
        copy = Chem.Mol(shown)
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(copy, core, acceptFailure=True)
        except (RuntimeError, ValueError):
            return shown
    return copy
