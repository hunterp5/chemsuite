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
from rdkit.Chem import rdDepictor
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


def render_png(
    mol: Chem.Mol,
    width: int,
    height: int,
    *,
    scale: float = 1.0,
    dark: bool = False,
    highlight_atoms: list[int] | None = None,
) -> bytes:
    """PNG bytes of *mol* centred in a transparent ``width x height`` pixel box.

    *scale* is the device pixel ratio (2 on Retina/HiDPI screens), so lines and labels stay the
    same physical size on every display.
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
                drawer, depiction_mol(mol), highlightAtoms=highlight_atoms or []
            )
        except (RuntimeError, ValueError):  # e.g. a structure that cannot be kekulized
            drawer = rdMolDraw2D.MolDraw2DCairo(max(int(width), 1), max(int(height), 1))
            drawer.drawOptions().clearBackground = False
            drawer.drawOptions().prepareMolsBeforeDrawing = False
            drawer.DrawMolecule(depiction_mol(mol))
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
