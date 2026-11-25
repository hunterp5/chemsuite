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
"""Docking from the table's ligands into the viewer's Poses panel: Tools ▸ Dock… and
Rescore Poses….

The docking itself is :mod:`chemsuite.macromolecule.dock` (Qt-free, table-free). This module is the
bridge: it reads the ligands out of a structure column, gives each one 3D coordinates if the
column holds none, runs the receptor and the ligands past AutoDock Vina in a background job,
and keeps the poses with the protein (:mod:`chemsuite.macromolecule.workspace.poses`), shown in the
**Poses** panel beside the pocket they were docked into.

Poses reach the table when the chemist asks (:func:`poses_to_table`, the panel's Copy to
Table): each becomes a row linked to its protein, so it follows into its pocket like any
ligand row and plots, filters and SAR work on it as on anything else.

Cancelling keeps the ligands that finished: the engine raises ``Cancelled(partial)`` and the
dialogs pass ``keep_partial=True`` (docking is per-row work).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import replace

import numpy as np
from rdkit import Chem

from chemsuite.app.document.tools import ToolContext
from chemsuite.conformers.etkdg import ConformerSettings, generate
from chemsuite.core.jobs import Cancelled
from chemsuite.core.mols import mol_from_blob, mol_from_text
from chemsuite.core.wording import plural
from chemsuite.macromolecule.chem import poses as pose_tools
from chemsuite.macromolecule.dock import export, pdbqt, run, vina
from chemsuite.macromolecule.dock.box import Box
from chemsuite.macromolecule.structure.model import Structure
from chemsuite.macromolecule.workspace.poses import DOCKED, RESCORED, DockedPoses, Pose

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


def docked_poses(
    protein: str,
    results: Sequence[run.Docked],
    names: Sequence[str],
    *,
    rescore: bool = False,
    skipped: int = 0,
    scoring: str = "",
) -> DockedPoses:
    """What a run produced, in the form the document keeps."""
    poses = [
        Pose(
            names[docked.index] if docked.index < len(names) else "",
            pose.rank,
            pose.affinity,
            pose.rmsd_lb,
            pose.rmsd_ub,
            Chem.MolToMolBlock(pose.mol),
        )
        for docked in results
        for pose in docked.poses
    ]
    failed = sum(1 for d in results if d.error)
    kind = RESCORED if rescore else DOCKED
    return DockedPoses(protein, tuple(poses), kind, failed, skipped, scoring)


def pose_rows(protein: str, poses: Sequence[Pose]) -> dict[str, list[str]]:
    """Every pose as a row: the protein it is in, the pose, and its scores."""
    rows: dict[str, list[str]] = {name: [] for name in POSE_COLUMNS}
    for pose in poses:
        rows[PROTEIN_COLUMN].append(protein)
        rows[LIGAND_COLUMN].append(pose.molblock)
        rows[POSE_COLUMN].append(str(pose.rank))
        rows[AFFINITY_COLUMN].append(f"{pose.affinity:.2f}")
        rows[RMSD_LB_COLUMN].append(f"{pose.rmsd_lb:.2f}")
        rows[RMSD_UB_COLUMN].append(f"{pose.rmsd_ub:.2f}")
        rows[PARENT_COLUMN].append(pose.ligand)
    return rows


def poses_to_table(ctx: ToolContext, protein: str, poses: Sequence[Pose]) -> int:
    """Copy *poses* into the table, one row each. Returns how many were copied."""
    if not poses:
        return 0
    rows = pose_rows(protein, poses)
    ctx.add_rows(rows, f"Copy {plural(len(poses), 'pose')} of {protein}",
                 structures=(LIGAND_COLUMN,))  # fmt: skip
    ctx.show_message(f"{plural(len(poses), 'pose')} copied to the table")
    return len(poses)


def single_conformers(mol: Chem.Mol) -> list[Chem.Mol]:
    """*mol*'s 3D conformers, one molecule each (empty when it has none)."""
    posed = pose_tools.pose_of(mol)
    if posed is None:
        return []
    out = []
    for conformer in posed.GetConformers():
        if not conformer.Is3D():
            continue
        one = Chem.Mol(posed)
        one.RemoveAllConformers()
        one.AddConformer(Chem.Conformer(conformer), assignId=True)
        out.append(one)
    return out


def ligand_inputs(
    mols: Sequence[Chem.Mol | None], names: Sequence[str], oids: np.ndarray, conformers: int
) -> tuple[list[export.LigandInput], int]:
    """The ligand files a save should write, and how many rows had no structure.

    A row keeps the conformers it already carries (the Conformers tool's, or a pose); a row
    with none, or with a flat drawing, gets *conformers* made for it. Vina searches the
    rotatable bonds itself, so more than one is worth sending only where the geometry it holds
    fixed differs — the pucker of a saturated ring.
    """
    settings = ConformerSettings(count=max(1, int(conformers)))
    out: list[export.LigandInput] = []
    skipped = 0
    for i, mol in enumerate(mols):
        if mol is None:
            skipped += 1
            continue
        name = names[i] if i < len(names) else ""
        row = int(oids[i]) if i < len(oids) else -1
        made = single_conformers(mol)
        if not made:
            made = single_conformers(generate(mol, settings).mol) if settings.count else []
        if not made:
            skipped += 1
            continue
        for k, one in enumerate(made, 1):
            out.append(export.LigandInput(one, name, row, k, len(made)))
    return out, skipped


def save_input(
    ctx: ToolContext,
    protein: str,
    structure: Structure | None,
    column: str,
    oids: np.ndarray,
    box: Box,
    options: vina.Options,
    folder: str,
    *,
    receptor: bool = True,
    ligands: bool = True,
    conformers: int = 1,
):
    """Write this docking run into *folder* to run at the command line, as a job.

    *receptor* and *ligands* say how much is written: the config alone (the search box and the
    conditions), the receptor beside it, or the ligands too with the manifest and the run
    scripts (:mod:`chemsuite.macromolecule.dock.export`).
    """
    cells = ctx.structure_inputs(column, oids) if ligands else []
    names = row_names(ctx, oids) if ligands else []

    def work(job):
        receptor_pdbqt = None
        dropped: tuple[str, ...] = ()
        if receptor and structure is not None:
            made = pdbqt.receptor(structure)
            receptor_pdbqt, dropped = made.pdbqt, made.dropped
        inputs, skipped = ([], 0)
        if ligands:
            inputs, skipped = ligand_inputs(ligand_mols(cells), names, oids, conformers)
        saved = export.save(folder, box, options, receptor_pdbqt, inputs,
                            cancel=job.cancel_token, progress=job.report)  # fmt: skip
        return saved, skipped, dropped

    def done(result) -> None:
        if result is None:
            return
        saved, skipped, dropped = result
        if saved.ligands:
            message = (f"{plural(saved.ligands, 'ligand file')} from "
                       f"{plural(saved.rows, 'row')} saved to {saved.folder}")  # fmt: skip
        else:
            message = f"Saved to {saved.folder}"
        if saved.failed:
            message += f"; {plural(len(saved.failed), 'ligand')} could not be written"
        if skipped:
            message += f"; {plural(skipped, 'row')} had no structure"
        if dropped:
            message += f". Left out of the receptor: {', '.join(dropped)}"
        ctx.show_message(message)

    title = f"Saving docking files for {protein}…" if protein else "Saving docking files…"
    return ctx.run(title, work, done)


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
    if not cells:
        return None
    return _dock_job(ctx, protein, structure, lambda: ligand_mols(cells),
                     row_names(ctx, oids), box, options, rescore)  # fmt: skip


def rescore_kept(
    ctx: ToolContext,
    protein: str,
    structure: Structure,
    poses: Sequence[Pose],
    options: vina.Options,
):
    """Score the poses already kept with *protein* where they stand, as a job.

    Returns the job, or None when there are no poses.
    """
    if not poses:
        return None
    return _dock_job(ctx, protein, structure, lambda: [mol_from_text(p.molblock) for p in poses],
                     [pose.ligand for pose in poses], None, options, True)  # fmt: skip


def _dock_job(
    ctx: ToolContext,
    protein: str,
    structure: Structure,
    ligands: Callable[[], list[Chem.Mol | None]],
    names: Sequence[str],
    box: Box | None,
    options: vina.Options,
    rescore: bool,
):
    """The background run behind docking and rescoring: *ligands* is read in the job (reading
    and embedding structures is work), and the poses are kept with the protein."""
    label = "Rescoring poses…" if rescore else f"Docking into {protein}…"
    cache = ctx.cache_folder("vina")

    def work(job):
        receptor = pdbqt.receptor(structure)
        made = [with_coordinates(mol) for mol in ligands()]
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
        run_poses = docked_poses(protein, results, names, rescore=rescore, skipped=skipped,
                                 scoring=options.scoring)  # fmt: skip
        poses = len(run_poses.poses)
        if poses:
            show_poses(ctx, run_poses)
        failed = run_poses.failed
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
        return "This protein has no hydrogens: prepare it first (Tools ▸ Prepare…)"
    return ""


def cancelled_partial(err: Cancelled) -> list[run.Docked]:
    """The ligands a cancelled run finished."""
    partial = err.partial
    return list(partial) if isinstance(partial, list) else []


# -- the poses panel ----------------------------------------------------------------------------
def show_poses(ctx: ToolContext, run_poses: DockedPoses):
    """Keep *run_poses* with the document, beside its protein's earlier runs, and show it in
    the Poses panel."""
    from .tool import macromolecule_workspace

    kept = macromolecule_workspace(ctx).add_docked(run_poses)
    return show_runs(ctx, kept.protein, kept.number)


def show_runs(ctx: ToolContext, protein: str, number: int | None = None):
    """Open the Poses panel on *protein*'s runs, run *number* shown (None: the newest), or
    fill the open one. It docks in the viewer's window, beside the pocket. None when the
    protein has no runs."""
    from chemsuite.macromolecule.ui.poses import PosesPanel

    from .tool import TOOL_ID, macromolecule_workspace

    runs = macromolecule_workspace(ctx).runs(protein)
    if not runs:
        return None
    panel = ctx.find_panel(PosesPanel.tool_id)
    if isinstance(panel, PosesPanel):
        panel.fill(runs, number)
        return panel
    # the protein is read from the panel when asked, not captured now: a later run into
    # another protein fills this same panel
    panel = PosesPanel(
        runs,
        show=lambda pose: show_pose(ctx, panel.run.protein, pose),
        to_table=lambda poses: poses_to_table(ctx, panel.run.protein, poses),
        remove=lambda run: remove_run(ctx, panel, run),
    )
    panel.fill(runs, number)
    ctx.add_panel(panel, panel.title(), beside=ctx.find_panel(TOOL_ID))
    return panel


def remove_run(ctx: ToolContext, panel, run: DockedPoses) -> bool:
    """Forget *run* (Remove Run): the panel shows the protein's next newest, or closes when
    none are left."""
    from .tool import macromolecule_workspace

    if not macromolecule_workspace(ctx).drop_docked(run.protein, run.number):
        return False
    if show_runs(ctx, run.protein) is None:
        panel.close_panel()
        dock = panel.parent()
        if dock is not None:
            dock.close()
    ctx.show_message(f"Removed run {run.number} of {run.protein}")
    return True


def open_poses(ctx: ToolContext, state: dict):
    """Reopen the Poses panel with the document (the runs were saved with the protein
    workspace), on the run it showed."""
    protein = state.get("protein")
    if not isinstance(protein, str):
        return None
    number = state.get("run")
    return show_runs(ctx, protein, number if isinstance(number, int) else None)


def show_pose(ctx: ToolContext, protein: str, pose: Pose) -> bool:
    """Draw *pose* in *protein*'s pocket (a pose was chosen in the panel)."""
    from .tool import macromolecule_workspace

    mol = mol_from_text(pose.molblock)
    if mol is None:
        return False
    return macromolecule_workspace(ctx).session().show_pose(protein, mol)
