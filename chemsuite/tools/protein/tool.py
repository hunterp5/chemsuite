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
"""Protein ▸ Protein Viewer: a workspace of its own for proteins, in a window of its own.

Proteins are opened and fetched from the viewer's File menu (this module hands it those
commands). A protein is stored in the document as an asset (``protein``, see
:mod:`chemsuite.protein`). Loading one does not touch the table: its ligands become rows, linked
to it by the ``Protein`` column, only when asked (Ligand ▸ Ligands to Table, Add Ligand to
Table, Interactions to Table); the viewer then shows a row's ligand in its pocket when that row
is current. Reading, fetching and working out ligands' bond orders run in the background.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

from chemsuite.app.document.tools import Tool, ToolContext
from chemsuite.app.shell.jobs import QUICK
from chemsuite.core.table import Asset, PutAsset
from chemsuite.protein import fetch as pdb
from chemsuite.protein import io, ligands, profile, store
from chemsuite.protein.components import META_KEY, CcdCache
from chemsuite.protein.model import Structure

TOOL_ID = "protein"
KIND = "protein"
MENU = ("Protein",)
#: the viewer's window: this share of the main window, centred over it
WINDOW_SHARE = 0.85
FILE_FILTER = "Protein files (*.pdb *.ent *.cif *.mmcif *.pdb.gz *.ent.gz *.cif.gz)"


def open_viewer(ctx: ToolContext, state: dict | None = None):
    """Bring the open Protein Viewer forward, or open one (with *state* when reopening)."""
    from chemsuite.app.protein.panel import (
        EXTRACT,
        FETCH,
        LIGAND,
        OPEN,
        PROFILE,
        REMOVE,
        SHOW_TABLE,
        ProteinPanel,
    )

    if state is None:
        panel = ctx.find_panel(TOOL_ID)
        if panel is not None:
            return panel
    commands = {
        OPEN: lambda: open_protein(ctx),
        FETCH: lambda: fetch_from_pdb(ctx),
        EXTRACT: lambda: extract_ligands(ctx),
        LIGAND: lambda protein, residue: add_ligand(ctx, protein, residue),
        PROFILE: lambda: interactions_to_table(ctx),
        REMOVE: lambda protein: remove_protein(ctx, protein),
        SHOW_TABLE: lambda: show_table(ctx),
    }

    def run(label, work, on_done):  # surfaces and interactions: never hold up closing
        return ctx.run(label, work, on_done, lane=QUICK, background=True)

    panel = ProteinPanel(
        ctx.table,
        current=ctx.current_oid,
        watch_current=ctx.watch_current,
        commands=commands,
        run=run,
    )
    if state is not None:
        panel.set_state(state)
    window = ctx.add_window(panel, panel.title(), menu_bar=panel.menu_bar,
                            on_close=panel.close_panel)  # fmt: skip
    if state is None:
        centre_over(window, ctx.parent)
    panel.setFocus()
    return panel


def centre_over(window, parent) -> None:
    """Size *window* to most of *parent*'s size, centred over it."""
    frame = parent.geometry()
    width, height = int(frame.width() * WINDOW_SHARE), int(frame.height() * WINDOW_SHARE)
    window.setGeometry(frame.x() + (frame.width() - width) // 2,
                       frame.y() + (frame.height() - height) // 2, width, height)  # fmt: skip


def show_table(ctx: ToolContext) -> None:
    """Bring the document's table window forward (Window ▸ Show Table in the viewer)."""
    window = ctx.parent
    if window.isMinimized():
        window.showNormal()
    window.raise_()
    window.activateWindow()


def dialog_parent(ctx: ToolContext):
    """The viewer's window when it is open (dialogs opened from it appear over it), else the
    document's."""
    panel = ctx.find_panel(TOOL_ID)
    return panel.window() if panel is not None else ctx.parent


def unique_id(ctx: ToolContext, name: str) -> str:
    """An asset id for a protein called ``name``: the name, or ``name-2``… when taken."""
    taken = set(ctx.table.assets)
    if name not in taken:
        return name
    n = 2
    while f"{name}-{n}" in taken:
        n += 1
    return f"{name}-{n}"


def ccd(ctx: ToolContext, offline: bool = False) -> CcdCache:
    """The Chemical Component Dictionary files kept on disk (``offline``: no downloads, for
    work on the GUI thread)."""
    return CcdCache(ctx.cache_folder("ccd"), offline=offline)


def ligand_work(ctx: ToolContext):
    """For a background job: the rows of a structure's ligands, the CCD asked when needed."""
    cache = ccd(ctx)

    def rows(structure: Structure, job) -> dict[str, list[str]]:
        return ligands.ligand_rows(structure, "", ccd=cache, cancel=job.cancel_token)

    return rows


def add_protein(
    ctx: ToolContext,
    structure: Structure,
    name: str,
    source: str,
    rows: dict[str, list[str]] | None = None,
) -> str:
    """Add the protein to the document (one undo step), show it, and return its id.

    The table is not touched: *rows* (from :func:`ligand_work`, worked out in the background)
    only say which ligands' bond orders are a guess, for the ⚠ in the proteins list."""
    protein = unique_id(ctx, name)
    if rows is None:
        rows = ligands.ligand_rows(structure, protein, ccd=ccd(ctx, offline=True))
    meta = {k: str(v) for k, v in structure.meta.items() if v and k != META_KEY}
    meta["source"] = source
    if guessed := ligands.guessed(rows):
        meta[ligands.GUESSED] = "\n".join(guessed)
    ctx.edit(f"Add protein {protein}",
             [PutAsset(protein, None, Asset(protein, KIND, protein, structure, meta))])  # fmt: skip
    count = len(rows[ligands.PROTEIN_COLUMN])
    ctx.show_message(f"Added {protein} ({structure.n_atoms:,} atoms, {count} ligand(s))")
    open_viewer(ctx)
    return protein


def open_protein(ctx: ToolContext, paths: list[str] | None = None):
    """Read protein files in the background and add each one."""
    if paths is None:
        from PySide6.QtWidgets import QFileDialog

        paths, _ = QFileDialog.getOpenFileNames(dialog_parent(ctx), "Open Protein", "", FILE_FILTER)
    files = io.readable(paths)
    if not files:
        return None
    rows_of = ligand_work(ctx)

    def work(job):
        out = []
        for path in files:
            try:
                structure = io.read(path)
            except (io.ProteinFileError, OSError) as err:
                out.append((path, None, None, str(err)))
                continue
            out.append((path, structure, rows_of(structure, job), ""))
        return out

    def done(read) -> None:
        problems = []
        for path, structure, rows, error in read:
            if structure is None:
                problems.append(f"{Path(path).name}: {error}")
            else:
                name = pdb.name_of(structure, _stem(path))
                add_protein(ctx, structure, name, Path(path).name, rows)
        if problems:
            ctx.show_message("Could not read " + "; ".join(problems), 10_000)

    return ctx.run("Reading proteins…", work, done)


def _stem(path: Path) -> str:
    name = Path(path).name
    for suffix in (".gz", ".cif", ".mmcif", ".pdb", ".ent"):
        name = name.removesuffix(suffix)
    return name


def fetch_from_pdb(ctx: ToolContext):
    from .dialogs import FetchDialog

    return ctx.show_dialog(FetchDialog(ctx))


def extract_ligands(ctx: ToolContext, then: Callable[[], object] | None = None):
    """Ligand ▸ Ligands to Table: add rows for the ligands of every protein that has none in
    the table yet. A background job (the CCD may be asked for bond orders); *then* runs once
    the rows are in. Returns the job, or None when there is nothing to add."""
    table = ctx.table
    linked = set()
    if ligands.PROTEIN_COLUMN in table:
        linked = {table.cell(ligands.PROTEIN_COLUMN, int(o)) for o in table.oids}
    pending = [(asset.id, store.as_structure(asset.payload)) for asset in table.assets.values()
               if asset.kind == KIND and asset.id not in linked]  # fmt: skip
    if not pending:
        ctx.show_message("No new ligands to add")
        return None
    cache = ccd(ctx)

    def work(job):
        return [(protein, ligands.ligand_rows(structure, protein, ccd=cache,
                                              cancel=job.cancel_token))
                for protein, structure in pending]  # fmt: skip

    def done(found) -> None:
        added = 0
        for protein, rows in found:
            if protein in ctx.table.assets and rows[ligands.PROTEIN_COLUMN]:
                ctx.add_rows(rows, f"Ligands of {protein} to table",
                             structures=(ligands.LIGAND_COLUMN,))  # fmt: skip
                added += len(rows[ligands.PROTEIN_COLUMN])
        ctx.show_message(f"Added {added} ligand row(s)" if added else "No ligands to add")
        if then is not None and added:
            then()

    return ctx.run("Reading ligands…", work, done)


def add_ligand(ctx: ToolContext, protein: str, residue: int) -> int:
    """Add a row for one ligand (a cofactor, say), unless the table has it already."""
    asset = ctx.table.assets.get(protein)
    if asset is None or asset.kind != KIND:
        return 0
    structure = store.as_structure(asset.payload)
    label = structure.residue_label(residue)
    table = ctx.table
    if ligands.PROTEIN_COLUMN in table and ligands.LIGAND_ID_COLUMN in table:
        for oid in table.oids:
            o = int(oid)
            if (table.cell(ligands.PROTEIN_COLUMN, o) == protein
                    and table.cell(ligands.LIGAND_ID_COLUMN, o) == label):  # fmt: skip
                ctx.show_message(f"{label} of {protein} is in the table already")
                return 0
    rows = ligands.ligand_rows(structure, protein, [residue], ccd=ccd(ctx, offline=True))
    if not rows[ligands.PROTEIN_COLUMN]:
        ctx.show_message(f"{label} could not be read as a molecule")
        return 0
    ctx.add_rows(rows, f"Add ligand {label}", structures=(ligands.LIGAND_COLUMN,))
    ctx.show_message(f"Added {label} of {protein}")
    return 1


KEEP_ROWS = "keep"
DELETE_ROWS = "delete"


def linked_rows(ctx: ToolContext, protein: str) -> np.ndarray:
    """The rows whose Protein is ``protein``."""
    table = ctx.table
    if ligands.PROTEIN_COLUMN not in table:
        return np.zeros(0, np.int64)
    values = table.read([ligands.PROTEIN_COLUMN])[ligands.PROTEIN_COLUMN]
    return table.oids[np.array([str(v) == protein for v in values.tolist()], bool)]


def remove_protein(ctx: ToolContext, protein: str, rows: str | None = None) -> bool:
    """Remove the protein, and its ligand rows unless ``rows`` is KEEP_ROWS (asked when None
    and it has rows), in one undo step. False when cancelled or there is no such protein."""
    asset = ctx.table.assets.get(protein)
    if asset is None or asset.kind != KIND:
        return False
    linked = linked_rows(ctx, protein)
    if rows is None and len(linked):
        rows = _ask_rows(ctx, asset.name, len(linked))
        if rows is None:
            return False
    deleting = rows != KEEP_ROWS
    ctx.remove_rows(linked if deleting else [], f"Remove protein {protein}",
                    extra_ops=[PutAsset(protein, asset, None)])  # fmt: skip
    gone = f" and {len(linked)} ligand row(s)" if deleting and len(linked) else ""
    ctx.show_message(f"Removed {protein}{gone}")
    return True


def _ask_rows(ctx: ToolContext, name: str, count: int) -> str | None:
    """Delete or keep the protein's ligand rows (None: cancel)."""
    from PySide6.QtWidgets import QMessageBox

    question = f"Remove {name} and its {count} ligand row(s)?"
    box = QMessageBox(
        QMessageBox.Icon.Question, "Remove Protein", question, parent=dialog_parent(ctx)
    )
    delete = box.addButton("Delete Rows", QMessageBox.ButtonRole.DestructiveRole)
    keep = box.addButton("Keep Rows", QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(delete)
    box.exec()
    clicked = box.clickedButton()
    return DELETE_ROWS if clicked is delete else KEEP_ROWS if clicked is keep else None


def interactions_to_table(ctx: ToolContext):
    """Count each ligand row's interactions with its protein, by kind, and name the residues
    (the selected rows, or every row): a queued job that keeps the rows finished if cancelled."""
    table = ctx.table
    if ligands.PROTEIN_COLUMN not in table:  # no ligand rows yet: asking for this adds them
        return extract_ligands(ctx, then=lambda: interactions_to_table(ctx))
    structures: dict[str, Structure] = {}
    oids, names = [], []
    for oid in ctx.target_oids(selected_only=len(ctx.selection) > 0).tolist():
        protein = str(table.cell(ligands.PROTEIN_COLUMN, oid) or "")
        asset = table.assets.get(protein)
        if asset is not None and asset.kind == KIND:
            if protein not in structures:
                structures[protein] = store.as_structure(asset.payload)
            oids.append(int(oid))
            names.append(protein)
    if not oids:
        ctx.show_message("No rows name a protein in this document")
        return None
    labels, smiles, inputs = _row_ligands(ctx, oids)

    def work(job):
        each = zip(names, labels, inputs, smiles, strict=True)
        rows = [profile.Row(structures[p], label, _mol(given), text)
                for p, label, given, text in each]  # fmt: skip
        return profile.for_rows(rows, progress=job.report, cancel=job.cancel_token)

    def done(found) -> None:
        done_oids = np.array(oids[: len(found)], np.int64)
        results = {
            title: (done_oids, np.array([p.counts[k] if p else np.nan for p in found], float))
            for k, title in enumerate(profile.COUNT_COLUMNS.values())
        }
        results[profile.CONTACTS_COLUMN] = (done_oids, [p.contacts if p else "" for p in found])
        ctx.apply_results(results, label="Interactions to table", overwrite=True)
        ctx.show_message(f"Interactions of {len(found):,} of {len(oids):,} ligand row(s)")

    return ctx.run(f"Finding interactions of {len(oids):,} ligands…", work, done,
                   keep_partial=True)  # fmt: skip


def _row_ligands(ctx: ToolContext, oids: list[int]) -> tuple[list[str], list[str], list]:
    """The rows' ligand IDs, their ligands' SMILES, and the ligands as a job takes them."""
    table = ctx.table
    labels = [""] * len(oids)
    if ligands.LIGAND_ID_COLUMN in table:
        labels = [str(table.cell(ligands.LIGAND_ID_COLUMN, o) or "") for o in oids]
    if ligands.LIGAND_COLUMN not in table or table.kind(ligands.LIGAND_COLUMN) != "structure":
        return labels, [""] * len(oids), [None] * len(oids)
    smiles = [str(table.cell(ligands.LIGAND_COLUMN, o) or "") for o in oids]
    return labels, smiles, ctx.structure_inputs(ligands.LIGAND_COLUMN, oids)


def _mol(given):
    """A molecule from :meth:`ToolContext.structure_inputs` (a blob, unparsed text or None)."""
    from chemsuite.core.mols import mol_from_blob, mol_from_text

    if isinstance(given, bytes):
        return mol_from_blob(given)
    return mol_from_text(given) if given else None


TOOLS = [
    Tool(
        id=TOOL_ID,
        title="Protein Viewer",
        menu=MENU,
        open=open_viewer,
        restore=lambda ctx, state: open_viewer(ctx, state),
        requires=("gemmi",),
        shortcut="Ctrl+Shift+R",
        tip="A workspace for proteins: open or fetch structures, cartoon and surface, "
        "selections, interactions and ligand pockets",
    ),
]
