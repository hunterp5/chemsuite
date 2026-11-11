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
"""Interaction profiles of ligand rows: for each row, how many interactions of each kind its
ligand makes with its protein, and with which residues (``GLU121 HB; MET120 HB HY``).

A row's ligand is its own pose when its 3D molecule is not just the residue its ligand ID names
(:mod:`.poses`), else that residue. These columns are what a protein–ligand interaction
fingerprint across many poses builds on.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken

from ..structure.model import WATER, Structure
from . import interactions as contacts
from . import poses
from .ligands import residue_named

#: one count column per kind of interaction, then the residues touched
COUNT_COLUMNS = {
    contacts.HBOND: "H-Bonds",
    contacts.SALT_BRIDGE: "Salt Bridges",
    contacts.PI_STACK: "π-Stacking",
    contacts.CATION_PI: "Cation–π",
    contacts.HYDROPHOBIC: "Hydrophobic Contacts",
    contacts.HALOGEN: "Halogen Bonds",
    contacts.METAL: "Metal Contacts",
    contacts.WATER_BRIDGE: "Water Bridges",
}
CONTACTS_COLUMN = "Contacts"
#: how the Contacts column names each kind
CODES = {
    contacts.HBOND: "HB",
    contacts.SALT_BRIDGE: "SB",
    contacts.PI_STACK: "PI",
    contacts.CATION_PI: "CP",
    contacts.HYDROPHOBIC: "HY",
    contacts.HALOGEN: "X",
    contacts.METAL: "M",
    contacts.WATER_BRIDGE: "W",
}


@dataclass(frozen=True)
class Row:
    """What a row names: its protein's structure, its ligand ID, its molecule and SMILES."""

    structure: Structure
    label: str = ""
    mol: Chem.Mol | None = None
    smiles: str = ""


@dataclass(frozen=True)
class Profile:
    counts: tuple[int, ...]  # in the order of COUNT_COLUMNS
    contacts: str


def query(row: Row) -> tuple[Structure, np.ndarray, int] | None:
    """The structure to look in, the ligand's atoms and its residue: the row's own pose
    (merged into the protein) or the residue it names; None when it names neither."""
    s = row.structure
    residue = residue_named(s, row.label) if row.label else -1
    pose = poses.pose_of(row.mol)
    if pose is not None and (residue < 0 or not poses.same_pose(s, residue, pose)):
        placed = poses.from_mol(pose)
        heavy = placed.xyz()[placed.element != 1]
        both = poses.merged(s.subset(~poses.overlapping(s, heavy)), placed)
        mask = np.zeros(both.n_atoms, bool)
        mask[both.n_atoms - placed.n_atoms :] = True
        return both, mask, both.n_residues - 1
    if residue < 0:
        return None
    mask = np.zeros(s.n_atoms, bool)
    mask[s.atoms_of(residue)] = True
    return s, mask, residue


def profile(row: Row) -> Profile | None:
    found = query(row)
    if found is None:
        return None
    s, mask, residue = found
    hits = contacts.find(s, mask, smiles={residue: row.smiles} if row.smiles else None)
    counts = hits.counts()
    return Profile(tuple(counts.get(k, 0) for k in COUNT_COLUMNS), describe(s, hits))


def describe(s: Structure, hits: contacts.Interactions) -> str:
    """The residues touched and how, in sequence order: ``GLU121 HB; MET120 HB HY``. Waters
    are left out (they show as water bridges to the residue beyond)."""
    kinds: dict[int, set[int]] = {}
    for kind, partner in zip(hits.kind.tolist(), hits.b.tolist(), strict=True):
        residue = int(s.residue_of[partner])
        if s.res_kind[residue] != WATER:
            kinds.setdefault(residue, set()).add(kind)
    parts = []
    for residue in sorted(kinds):
        name = f"{s.res_name[residue]}{s.res_seq[residue]}{s.icode[residue]}"
        codes = " ".join(CODES[k] for k in COUNT_COLUMNS if k in kinds[residue])
        parts.append(f"{name} {codes}")
    return "; ".join(parts)


def for_rows(
    rows: Sequence[Row | None],
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> list[Profile | None]:
    """The profile of each row (None for a row naming no ligand). Cancelling raises
    ``Cancelled`` carrying the profiles of the rows finished, in order."""
    out: list[Profile | None] = []
    for row in rows:
        if cancel is not None and cancel.cancelled:
            raise Cancelled(out)
        out.append(profile(row) if row is not None else None)
        if progress is not None:
            progress(len(out), len(rows))
    return out
