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
"""The Browser's 3D molecules: quick conformers for 2D rows, overlay colors, and the
reference a superposed row was fitted to (Chemistry ▸ 3D ▸ Superpose)."""

from __future__ import annotations

from rdkit import Chem

from chemsuite.molview.geometry import POLAR_H, THIN, molecule_part

#: carbon colors for overlays (distinct on white and on dark backgrounds)
PALETTE = [
    (0.55, 0.55, 0.55), (0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17),
    (0.58, 0.4, 0.74), (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.74, 0.74, 0.13),
    (0.09, 0.75, 0.81), (0.84, 0.15, 0.16),
]  # fmt: skip
#: the reference's carbons: the overlay's second color, so it stands apart from the row's gray
REFERENCE_CARBON = PALETTE[1]


def is_3d(mol: Chem.Mol | None) -> bool:
    return mol is not None and mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def make_3d(mols: dict[int, Chem.Mol]) -> dict[int, Chem.Mol | None]:
    """A quick conformer for each structure that has none (None when embedding fails)."""
    from chemsuite.conformers.etkdg import ConformerSettings, generate

    settings = ConformerSettings(count=1, energy_window=None, timeout=10)
    out: dict[int, Chem.Mol | None] = {}
    for oid, mol in mols.items():
        try:
            out[oid] = generate(mol, settings).mol
        except Exception:  # noqa: BLE001 - one bad structure must not stop the rest
            out[oid] = None
    return out


def reference_part(reference: Chem.Mol | None):
    """The reference drawn as thin sticks with its own carbon color (None: nothing to draw)."""
    if not is_3d(reference):
        return None
    return molecule_part(reference, style=THIN, hydrogens=POLAR_H, carbon=REFERENCE_CARBON)
