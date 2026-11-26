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

"""A part of a protein (a chain, a ligand, its waters…) included in calculations or left out,
deleted or duplicated: each one undoable step of the document.

A part is the atoms a selection expression picks in the whole structure, alternate locations
and excluded atoms included, so an excluded part can still be included, deleted or copied.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from chemsuite.core.assets import PutAsset
from chemsuite.macromolecule.structure import edit as edits
from chemsuite.macromolecule.structure import parts
from chemsuite.macromolecule.structure.select import select

from .session import WorkbenchSession


def atoms_of(session: WorkbenchSession, protein: str, expression: str) -> np.ndarray:
    """The atoms of *protein* that *expression* picks (none when there is no such protein)."""
    shown = session.display(protein)
    if shown is None:
        return np.zeros(0, dtype=bool)
    return select(shown.structure, expression, model=shown.model)


def included_share(session: WorkbenchSession, protein: str, atoms: np.ndarray) -> float:
    """How much of *atoms* is included: 1 all, 0 none."""
    shown = session.display(protein)
    if shown is None or not atoms.any():
        return 1.0
    return 1.0 - float(shown.excluded[atoms].mean())


def set_included(
    session: WorkbenchSession, protein: str, atoms: np.ndarray, on: bool, what: str
) -> bool:
    """Include *atoms* of *protein* in calculations (*on*) or leave them out; False when
    nothing changed."""
    asset, shown = session.store.assets.get(protein), session.display(protein)
    if asset is None or shown is None or not atoms.any():
        return False
    s = shown.structure
    before = parts.excluded(asset.meta, s)
    after = before & ~atoms if on else before | atoms
    if np.array_equal(before, after):
        return False
    meta = parts.with_excluded(asset.meta, s, after)
    label = f"{'Include' if on else 'Exclude'} {what}"
    session.store.edit(label, [PutAsset(protein, asset, replace(asset, meta=meta))])
    return True


def delete(session: WorkbenchSession, protein: str, atoms: np.ndarray, what: str) -> int:
    """Delete *atoms* of *protein*; how many went. Raises EditError."""
    shown = session.display(protein)
    if shown is None or not atoms.any():
        return 0
    made = edits.delete(shown.structure, atoms)
    session.edit(f"Delete {what}", {protein: made})
    return shown.structure.n_atoms - made.n_atoms


def duplicate(
    session: WorkbenchSession, protein: str, atoms: np.ndarray, what: str, chain: bool
) -> int:
    """Copy the residues *atoms* belong to into *protein*: each in its own chain, numbered on,
    or (*chain*) together as a new chain. How many residues were added. Raises EditError."""
    shown = session.display(protein)
    if shown is None or not atoms.any():
        return 0
    s = shown.structure
    residues = np.unique(s.residue_of[atoms])
    name = parts.free_chain_name(s) if chain else ""
    session.edit(f"Duplicate {what}", {protein: parts.copy_residues(s, residues, name)})
    return len(residues)
