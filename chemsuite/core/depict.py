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

from collections.abc import Mapping, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdDepictor.SetPreferCoordGen(True)

#: Bond length in device pixels at scale 1: keeps small molecules from being blown up so a
#: column of compounds reads at one consistent size; large ones still shrink to fit.
_BOND_PX = 22.0
#: Property carrying an atom's index in the molecule as given, while hydrogens are removed.
_SOURCE = "_depictSourceAtom"


def depiction_mol(mol: Chem.Mol) -> Chem.Mol:
    """*mol* itself if it has 2D coordinates, else a copy with a computed 2D layout.

    3D structures (conformers) are drawn without their explicit hydrogens, as chemists draw.
    Each kept atom of such a copy remembers its index in *mol* (see :func:`drawn_atoms`), since
    removing hydrogens renumbers the rest.
    """
    if mol.GetNumConformers() and not mol.GetConformer().Is3D():
        return mol
    with rdBase.BlockLogs():
        if mol.GetNumConformers():
            copy = Chem.Mol(mol)
            for atom in copy.GetAtoms():
                atom.SetIntProp(_SOURCE, atom.GetIdx())
            copy = Chem.RemoveHs(copy)
        else:
            copy = Chem.Mol(mol)
    copy.RemoveAllConformers()
    with rdBase.BlockLogs():
        rdDepictor.Compute2DCoords(copy)
    return copy


def drawn_atoms(shown: Chem.Mol, atoms: Sequence[int]) -> list[int]:
    """*atoms*, which index the molecule as given, as indices into *shown*.

    The two numberings differ only when :func:`depiction_mol` removed hydrogens; an atom that
    is not drawn (a hydrogen, or an index past the end) is left out.
    """
    where = _drawn_index(shown)
    if where is None:
        return [atom for atom in atoms if atom < shown.GetNumAtoms()]
    return [where[atom] for atom in atoms if atom in where]


def _drawn_index(shown: Chem.Mol) -> dict[int, int] | None:
    """Where each atom of the molecule as given ended up in *shown*; None when unchanged."""
    first = shown.GetAtomWithIdx(0) if shown.GetNumAtoms() else None
    if first is None or not first.HasProp(_SOURCE):
        return None
    return {atom.GetIntProp(_SOURCE): atom.GetIdx() for atom in shown.GetAtoms()}


def drawn_highlight(
    shown: Chem.Mol, atoms: Sequence[tuple[int, float | None]], dark: bool = False
) -> tuple[list[int], dict[int, tuple[float, float, float]]]:
    """The atoms of *shown* to highlight and their colors, from ``(atom, score)`` pairs that
    index the molecule as given (best first; see :func:`graded_colors` for the shading)."""
    where = _drawn_index(shown)
    kept: list[tuple[int, float | None]] = []
    for atom, score in atoms:
        if where is None:
            if atom < shown.GetNumAtoms():
                kept.append((atom, score))
        elif atom in where:
            kept.append((where[atom], score))
    return [atom for atom, _ in kept], graded_colors(kept, dark)


#: The strongest highlight: a warm orange that reads on white and on the dark theme's base.
_HIGHLIGHT = (0.93, 0.44, 0.18)
#: What a faded highlight blends toward, light theme and dark.
_PAPER, _INK = (1.0, 1.0, 1.0), (0.12, 0.12, 0.12)
#: How far each weaker rank fades toward the background, and as far as fading goes.
_FADE_STEP, _FADE_MOST = 0.34, 0.68


def graded_colors(
    atoms: Sequence[tuple[int, float | None]], dark: bool = False
) -> dict[int, tuple[float, float, float]]:
    """A highlight color per atom of *atoms* (``(index, score)``, best first): the best site in
    full strength, each weaker score fainter.

    Atoms of equal score share a color — symmetry-equivalent atoms score identically, and
    shading them differently would invent a ranking the model did not make. Atoms with no score
    at all are all drawn in full strength.
    """
    background = _INK if dark else _PAPER
    colors: dict[int, tuple[float, float, float]] = {}
    rank = 0
    last: float | None = atoms[0][1] if atoms else None
    for atom, score in atoms:
        if score != last:
            rank += 1
        last = score
        colors[atom] = _faded(background, min(_FADE_MOST, _FADE_STEP * rank))
    return colors


def _faded(background: tuple[float, float, float], fade: float) -> tuple[float, float, float]:
    """The highlight color blended *fade* of the way toward *background*."""
    red, green, blue = _HIGHLIGHT
    paper_red, paper_green, paper_blue = background
    keep = 1.0 - fade
    return (
        red * keep + paper_red * fade,
        green * keep + paper_green * fade,
        blue * keep + paper_blue * fade,
    )


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
    highlight_colors: Mapping[int, tuple[float, float, float]] | None = None,
    attachments: bool = False,
) -> bytes:
    """PNG bytes of *mol* centred in a transparent ``width x height`` pixel box.

    *scale* is the device pixel ratio (2 on Retina/HiDPI screens), so lines and labels stay the
    same physical size on every display. Highlight indices refer to :func:`depiction_mol`'s
    atoms and bonds (the same as *mol*'s unless *mol* is 3D); *highlight_colors* gives an atom
    a color of its own (see :func:`graded_colors`) instead of the one color for all. With
    *attachments*, dummy atoms (``*``) are drawn as attachment points (a wavy line) and methyl
    groups are written out, as for a substituent.
    """
    drawer = rdMolDraw2D.MolDraw2DCairo(max(int(width), 1), max(int(height), 1))
    options = drawer.drawOptions()
    if dark:
        rdMolDraw2D.SetDarkMode(options)
    options.clearBackground = False
    options.dummiesAreAttachments = attachments
    options.explicitMethyl = attachments  # a substituent's lone carbon reads as CH3, not a stub
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
                highlightAtomColors=dict(highlight_colors or {}),
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
