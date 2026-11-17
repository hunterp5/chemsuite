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

from collections.abc import Callable, Sequence
from dataclasses import dataclass

from rdkit import Chem

from chemsuite.conformers import energies
from chemsuite.molview.align import aligned_conformers
from chemsuite.molview.geometry import NO_H, POLAR_H, THIN, WIRE, molecule_part

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
            made = generate(mol, settings).mol
            if made is not None:
                energies.clear(made)  # one quick conformer is no ensemble to compare with
            out[oid] = made
        except Exception:  # noqa: BLE001 - one bad structure must not stop the rest
            out[oid] = None
    return out


@dataclass(frozen=True)
class Look:
    """How the Browser draws its 3D molecules. *conformer* is the one shown of each row, unless
    *overlay* draws every conformer, or several are *chosen* (positions, 0 first); *rows*
    draws several rows at once. Any of these makes the drawing thin sticks without the
    hydrogens on carbon, as 3Dmol overlays read best."""

    style: str
    hydrogens: str
    by_color: bool
    align: bool
    overlay: bool
    rows: bool
    conformer: int
    reference: bool
    limit: int
    chosen: tuple[int, ...] = ()


def scene_parts(
    molecules: Sequence[tuple[int, Chem.Mol]],
    look: Look,
    reference_of: Callable[[Chem.Mol], Chem.Mol | None],
) -> tuple[list, list[int]]:
    """The parts to draw and the row each belongs to (-1 for the reference)."""
    style, hydrogens = look.style, look.hydrogens
    several = len(look.chosen) > 1 and not look.overlay
    if look.overlay or look.rows or several:
        style = WIRE if style == WIRE else THIN
        hydrogens = NO_H if hydrogens == NO_H else POLAR_H
    parts, rows = [], []
    for oid, mol in molecules:
        if look.align and mol.GetNumConformers() > 1:
            mol = aligned_conformers(mol)
        ids = [c.GetId() for c in mol.GetConformers()]
        if look.overlay:
            chosen = ids
        elif several:
            chosen = [ids[i] for i in look.chosen if i < len(ids)] or ids[:1]
        else:
            chosen = [ids[min(look.conformer, len(ids) - 1)]]
        for conf_id in chosen:
            if len(parts) >= look.limit:
                break
            carbon = PALETTE[len(parts) % len(PALETTE)] if look.by_color else None
            parts.append(molecule_part(mol, conf_id, style=style, hydrogens=hydrogens,
                                       carbon=carbon))  # fmt: skip
            rows.append(oid)
    if look.reference and molecules and len(parts) < look.limit:
        reference = reference_part(reference_of(molecules[0][1]))
        if reference is not None:
            parts.append(reference)
            rows.append(-1)  # not a row
    return parts, rows


def reference_part(reference: Chem.Mol | None):
    """The reference drawn as thin sticks with its own carbon color (None: nothing to draw)."""
    if not is_3d(reference):
        return None
    return molecule_part(reference, style=THIN, hydrogens=POLAR_H, carbon=REFERENCE_CARBON)
