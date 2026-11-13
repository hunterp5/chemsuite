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

How the proteins look, named selections, focus, labels, measurements and the view live in the
document's protein workspace (:class:`~chemsuite.protein.workspace.workspace.ProteinWorkspace`,
a document part saved under ``protein_workspace``), not in the window: they are saved with the
document whether or not the viewer is open, and the window attaches to the workspace's session
when it opens. Documents saved before kept them with the window; they move into the workspace
when the window is restored.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from chemsuite.app.document.tools import Tool, ToolContext
from chemsuite.app.shell.jobs import QUICK
from chemsuite.core.table import Asset, PutAsset
from chemsuite.core.wording import plural
from chemsuite.protein.chem import ligands
from chemsuite.protein.structure import fetch as pdb
from chemsuite.protein.structure import io
from chemsuite.protein.structure.components import META_KEY, CcdCache
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.workspace.workspace import ProteinWorkspace

from .follow import TableRow
from .rows import PROTEIN_COLUMN, ligand_rows
from .rows import guessed as guessed_ligands

TOOL_ID = "protein"
#: the document part holding the protein workspace
WORKSPACE_KEY = "protein_workspace"
KIND = "protein"
MENU = ("Protein",)
#: the viewer's window: this share of the main window, centred over it
WINDOW_SHARE = 0.85
FILE_FILTER = "Protein files (*.pdb *.ent *.cif *.mmcif *.pdb.gz *.ent.gz *.cif.gz)"


def open_viewer(ctx: ToolContext, state: dict | None = None):
    """Bring the open Protein Viewer forward, or open one (with *state* when reopening)."""
    from chemsuite.protein.ui.panel import (  # the workbench loads when first opened
        EXTRACT,
        FETCH,
        LIGAND,
        OPEN,
        PROFILE,
        REMOVE,
        SHOW_TABLE,
        ProteinPanel,
    )

    from .to_table import add_ligand, extract_ligands, interactions_to_table

    if state is None:
        panel = ctx.find_panel(TOOL_ID)
        if panel is not None:
            return panel
    workbench = protein_workspace(ctx)
    if state is not None and "session" in state:  # saved before the workspace kept it
        workbench.adopt(state)
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
        session=workbench.session(),
        current=ctx.current_oid,
        watch_current=ctx.watch_current,
        commands=commands,
        run=run,
        linked=TableRow(ctx),
    )
    panel.set_view_state(workbench.view)
    if state is not None:
        panel.set_state({k: v for k, v in state.items() if k in ("columns", "sequence")})
    workbench.view_source = panel.view_state

    def closed() -> None:
        if workbench.view_source == panel.view_state:
            workbench.set_view(panel.view_state())
            workbench.view_source = None
        panel.close_panel()

    window = ctx.add_window(panel, panel.title(), menu_bar=panel.menu_bar, on_close=closed)
    if state is None:
        centre_over(window, ctx.parent)
    panel.setFocus()
    return panel


def protein_workspace(ctx: ToolContext) -> ProteinWorkspace:
    """The document's protein workspace (made, from what the document saved, when first
    asked for)."""
    rows = TableRow(ctx)

    def make(saved) -> ProteinWorkspace:
        return ProteinWorkspace(ctx.table, saved, relabel=rows.relabel, count_named=rows.naming)

    return ctx.document_part(WORKSPACE_KEY, make)


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
        return ligand_rows(structure, "", ccd=cache, cancel=job.cancel_token)

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
        rows = ligand_rows(structure, protein, ccd=ccd(ctx, offline=True))
    meta = {k: str(v) for k, v in structure.meta.items() if v and k != META_KEY}
    meta["source"] = source
    if guessed := guessed_ligands(rows):
        meta[ligands.GUESSED] = "\n".join(guessed)
    ctx.edit(f"Add protein {protein}",
             [PutAsset(protein, None, Asset(protein, KIND, protein, structure, meta))])  # fmt: skip
    count = len(rows[PROTEIN_COLUMN])
    ctx.show_message(f"Added {protein} ({structure.n_atoms:,} atoms, {plural(count, 'ligand')})")
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
            ctx.show_message("Could not read " + "; ".join(problems))

    return ctx.run("Reading proteins…", work, done)


def _stem(path: Path) -> str:
    name = Path(path).name
    for suffix in (".gz", ".cif", ".mmcif", ".pdb", ".ent"):
        name = name.removesuffix(suffix)
    return name


def fetch_from_pdb(ctx: ToolContext):
    from .dialogs import FetchDialog

    return ctx.show_dialog(FetchDialog(ctx))


KEEP_ROWS = "keep"
DELETE_ROWS = "delete"


def linked_rows(ctx: ToolContext, protein: str) -> np.ndarray:
    """The rows whose Protein is ``protein``."""
    table = ctx.table
    if PROTEIN_COLUMN not in table:
        return np.zeros(0, np.int64)
    values = table.read([PROTEIN_COLUMN])[PROTEIN_COLUMN]
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
