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
"""Docking as table rows: the work behind Protein ▸ Dock… and Rescore Poses….

The docking itself is :mod:`chemsuite.protein.dock` (Qt-free, table-free). This module is the
bridge: it reads the ligands out of a structure column, gives each one 3D coordinates if the
column holds none, runs the receptor and the ligands past AutoDock Vina in a background job,
and writes each pose as a row linked to its protein — so a pose follows into its pocket like
any ligand row, and plots, filters and SAR work on poses as on anything else.

Cancelling keeps the ligands that finished: the engine raises ``Cancelled(partial)`` and the
dialogs pass ``keep_partial=True`` (docking is per-row work).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

import numpy as np
from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.conformers.etkdg import ConformerSettings, generate
from chemsuite.core.jobs import Cancelled
from chemsuite.core.mols import mol_from_blob, mol_from_text
from chemsuite.core.wording import plural
from chemsuite.protein.chem import poses as pose_tools
from chemsuite.protein.dock import pdbqt, run, vina
from chemsuite.protein.dock.box import Box
from chemsuite.protein.structure.model import Structure

from .rows import LIGAND_COLUMN, PROTEIN_COLUMN

POSE_COLUMN = "Pose"
AFFINITY_COLUMN = "Affinity (kcal/mol)"
RMSD_LB_COLUMN = "RMSD l.b."
RMSD_UB_COLUMN = "RMSD u.b."
PARENT_COLUMN = "Docked From"
#: the columns docking writes, in the order they are added
POSE_COLUMNS = (
    PROTEIN_COLUMN,
    LIGAND_COLUMN,
    POSE_COLUMN,
    AFFINITY_COLUMN,
    RMSD_LB_COLUMN,
    RMSD_UB_COLUMN,
    PARENT_COLUMN,
)
HYDROGEN = 1
#: one conformer, for a ligand whose column holds no 3D coordinates
ONE_CONFORMER = ConformerSettings(count=1)


def has_hydrogens(s: Structure) -> bool:
    """Has this protein been prepared? Vina types its atoms by the hydrogens on them."""
    return bool((s.element == HYDROGEN).any())


def ligand_mols(cells: Sequence[bytes | str | None]) -> list[Chem.Mol | None]:
    """The molecules of a structure column's cells, in order (None for a blank or unreadable
    cell)."""
    out: list[Chem.Mol | None] = []
    for cell in cells:
        if isinstance(cell, bytes):
            out.append(mol_from_blob(cell) if cell else None)
        elif isinstance(cell, str):
            out.append(mol_from_text(cell))
        else:
            out.append(None)
    return out


def with_coordinates(mol: Chem.Mol | None) -> Chem.Mol | None:
    """*mol* in 3D: as it is when it has a pose, else one conformer of it."""
    if mol is None:
        return None
    if pose_tools.pose_of(mol) is not None:
        return mol
    made = generate(mol, ONE_CONFORMER)
    return made.mol


def pose_rows(
    protein: str, results: Sequence[run.Docked], names: Sequence[str]
) -> dict[str, list[str]]:
    """Every pose of *results* as a row: the protein it is in, the pose, and its scores."""
    rows: dict[str, list[str]] = {name: [] for name in POSE_COLUMNS}
    for docked in results:
        parent = names[docked.index] if docked.index < len(names) else ""
        for pose in docked.poses:
            rows[PROTEIN_COLUMN].append(protein)
            rows[LIGAND_COLUMN].append(Chem.MolToMolBlock(pose.mol))
            rows[POSE_COLUMN].append(str(pose.rank))
            rows[AFFINITY_COLUMN].append(f"{pose.affinity:.2f}")
            rows[RMSD_LB_COLUMN].append(f"{pose.rmsd_lb:.2f}")
            rows[RMSD_UB_COLUMN].append(f"{pose.rmsd_ub:.2f}")
            rows[PARENT_COLUMN].append(parent)
    return rows


def dock_rows(
    ctx: ToolContext,
    protein: str,
    structure: Structure,
    column: str,
    oids: np.ndarray,
    box: Box | None,
    options: vina.Options,
    rescore: bool = False,
):
    """Dock (or rescore) the ligands of *column* at *oids* against *protein*, as a job.

    Returns the job, or None when there is nothing to dock.
    """
    cells = ctx.structure_inputs(column, oids)
    names = row_names(ctx, oids)
    if not cells:
        return None
    label = "Rescoring poses…" if rescore else f"Docking into {protein}…"
    cache = ctx.cache_folder("vina")

    def work(job):
        receptor = pdbqt.receptor(structure)
        made = [with_coordinates(mol) for mol in ligand_mols(cells)]
        # a blank cell, or a structure RDKit could not read, is skipped rather than docked
        usable = [(i, mol) for i, mol in enumerate(made) if mol is not None]
        skipped = len(made) - len(usable)
        job.report(0, len(usable))
        many = run.rescore_many if rescore else run.dock_many
        try:
            results = many(receptor.pdbqt, [mol for _i, mol in usable], box, options,
                           progress=job.report, cancel=job.cancel_token, cache=cache)  # fmt: skip
        except Cancelled as stop:  # the ligands that finished are still worth their rows
            kept = _for_rows(cancelled_partial(stop), usable)
            raise Cancelled((receptor, kept, skipped)) from stop
        return receptor, _for_rows(results, usable), skipped

    def done(result) -> None:
        if result is None:
            return
        receptor, results, skipped = result
        rows = pose_rows(protein, results, names)
        poses = len(rows[PROTEIN_COLUMN])
        if poses:
            step = f"Rescore poses of {protein}" if rescore else f"Dock into {protein}"
            ctx.add_rows(rows, step, structures=(LIGAND_COLUMN,))
        failed = sum(1 for d in results if d.error)
        message = f"{plural(poses, 'pose')} from {plural(len(results) - failed, 'ligand')}"
        if failed:
            message += f"; {failed} could not be docked"
        if skipped:
            message += f"; {plural(skipped, 'row')} had no structure"
        if receptor.dropped:
            left = ", ".join(receptor.dropped)
            message += f". Left out of the receptor: {left}"
        ctx.show_message(message)

    return ctx.run(label, work, done, keep_partial=True)


def _for_rows(
    results: Sequence[run.Docked], usable: Sequence[tuple[int, Chem.Mol]]
) -> list[run.Docked]:
    """The results with each one's ``index`` back on the row it came from (blank rows were
    not sent to Vina, so the engine counted only the ligands it was given)."""
    return [replace(d, index=usable[d.index][0]) for d in results if d.index < len(usable)]


def row_names(ctx: ToolContext, oids: np.ndarray) -> list[str]:
    """What each row is called: its name column's cell, else its row number."""
    table = ctx.table
    name_column = next((n for n in ("Name", "ID", "Compound") if n in table), None)
    if name_column is None:
        return [str(int(oid)) for oid in oids]
    return [str(table.cell(name_column, int(oid)) or int(oid)) for oid in oids]


def docking_problem(structure: Structure | None) -> str:
    """Why this protein cannot be docked into, or ``""``."""
    if structure is None:
        return "There is no such protein in this document"
    if not has_hydrogens(structure):
        return "This protein has no hydrogens: prepare it first (Protein ▸ Prepare…)"
    return ""


def cancelled_partial(err: Cancelled) -> list[run.Docked]:
    """The ligands a cancelled run finished."""
    partial = err.partial
    return list(partial) if isinstance(partial, list) else []
