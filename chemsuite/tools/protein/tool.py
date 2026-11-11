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
"""Protein ▸ Protein Viewer: a workspace of its own for proteins, opened as its own window.

Proteins are opened, fetched and turned into ligand rows from the viewer's File and Pocket
menus (this module hands it those commands). A protein is stored in the document as an asset
(``protein``, see :mod:`chemsuite.protein`); its ligands become rows linked to it by the
``Protein`` column, so the viewer shows a row's ligand in its pocket when that row is current.
Adding a protein and its ligand rows is one undo step. Reading and fetching run in the
background.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from chemsuite.app.jobs import QUICK
from chemsuite.app.tools import Tool, ToolContext
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
        ProteinPanel,
    )

    if state is None:
        panel = ctx.find_panel(TOOL_ID)
        if panel is not None:
            panel.setFocus()
            return panel
    commands = {
        OPEN: lambda: open_protein(ctx),
        FETCH: lambda: fetch_from_pdb(ctx),
        EXTRACT: lambda: extract_ligands(ctx),
        LIGAND: lambda protein, residue: add_ligand(ctx, protein, residue),
        PROFILE: lambda: interactions_to_table(ctx),
        REMOVE: lambda protein: remove_protein(ctx, protein),
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
    dock = ctx.add_panel(panel, panel.title())
    if state is None:
        float_over(dock, ctx.parent)
    panel.setFocus()
    return panel


def float_over(dock, window) -> None:
    """Undock *dock* as its own window, most of *window*'s size, centred over it."""
    dock.setFloating(True)
    frame = window.geometry()
    width, height = int(frame.width() * WINDOW_SHARE), int(frame.height() * WINDOW_SHARE)
    dock.setGeometry(frame.x() + (frame.width() - width) // 2,
                     frame.y() + (frame.height() - height) // 2, width, height)  # fmt: skip


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
    """Add the protein and its ligands (as rows, from :func:`ligand_work` when given) in one
    undo step, show it, and return its id."""
    protein = unique_id(ctx, name)
    if rows is None:
        rows = ligands.ligand_rows(structure, protein, ccd=ccd(ctx, offline=True))
    rows[ligands.PROTEIN_COLUMN] = [protein] * len(rows[ligands.PROTEIN_COLUMN])
    meta = {k: str(v) for k, v in structure.meta.items() if v and k != META_KEY}
    meta["source"] = source
    if guessed := ligands.guessed(rows):
        meta[ligands.GUESSED] = "\n".join(guessed)
    op = PutAsset(protein, None, Asset(protein, KIND, protein, structure, meta))
    label = f"Add protein {protein}"
    if rows[ligands.PROTEIN_COLUMN]:
        ctx.add_rows(rows, label, structures=(ligands.LIGAND_COLUMN,), extra_ops=[op])
    else:
        ctx.edit(label, [op])
    count = len(rows[ligands.PROTEIN_COLUMN])
    ctx.show_message(f"Added {protein} ({structure.n_atoms:,} atoms, {count} ligand(s))")
    open_viewer(ctx)
    return protein


def open_protein(ctx: ToolContext, paths: list[str] | None = None):
    """Read protein files in the background and add each one."""
    if paths is None:
        from PySide6.QtWidgets import QFileDialog

        paths, _ = QFileDialog.getOpenFileNames(ctx.parent, "Open Protein", "", FILE_FILTER)
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


def extract_ligands(ctx: ToolContext) -> int:
    """Add rows for the ligands of every protein that has none in the table yet."""
    table = ctx.table
    linked = set()
    if ligands.PROTEIN_COLUMN in table:
        linked = {table.cell(ligands.PROTEIN_COLUMN, int(o)) for o in table.oids}
    added = 0
    for asset in list(table.assets.values()):
        if asset.kind != KIND or asset.id in linked:
            continue
        rows = ligands.ligand_rows(store.as_structure(asset.payload), asset.id,
                                   ccd=ccd(ctx, offline=True))  # fmt: skip
        if rows[ligands.PROTEIN_COLUMN]:
            ctx.add_rows(rows, f"Extract ligands of {asset.id}",
                         structures=(ligands.LIGAND_COLUMN,))  # fmt: skip
            added += len(rows[ligands.PROTEIN_COLUMN])
    ctx.show_message(f"Added {added} ligand row(s)" if added else "No new ligands to add")
    return added


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
    box = QMessageBox(QMessageBox.Icon.Question, "Remove Protein", question, parent=ctx.parent)
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
    if ligands.PROTEIN_COLUMN not in table:
        ctx.show_message("No ligand rows: add a protein first")
        return None
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
