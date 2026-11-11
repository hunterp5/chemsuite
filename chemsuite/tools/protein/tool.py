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

from chemsuite.app.jobs import QUICK
from chemsuite.app.tools import Tool, ToolContext
from chemsuite.core.table import Asset, PutAsset
from chemsuite.protein import fetch as pdb
from chemsuite.protein import io, ligands, store
from chemsuite.protein.model import Structure

TOOL_ID = "protein"
KIND = "protein"
MENU = ("Protein",)
#: the viewer's window: this share of the main window, centred over it
WINDOW_SHARE = 0.85
FILE_FILTER = "Protein files (*.pdb *.ent *.cif *.mmcif *.pdb.gz *.ent.gz *.cif.gz)"


def open_viewer(ctx: ToolContext, state: dict | None = None):
    """Bring the open Protein Viewer forward, or open one (with *state* when reopening)."""
    from chemsuite.app.protein.panel import EXTRACT, FETCH, OPEN, ProteinPanel

    if state is None:
        panel = ctx.find_panel(TOOL_ID)
        if panel is not None:
            panel.setFocus()
            return panel
    commands = {
        OPEN: lambda: open_protein(ctx),
        FETCH: lambda: fetch_from_pdb(ctx),
        EXTRACT: lambda: extract_ligands(ctx),
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


def add_protein(ctx: ToolContext, structure: Structure, name: str, source: str) -> str:
    """Add the protein and its ligands (as rows) in one undo step, show it, and return its id."""
    protein = unique_id(ctx, name)
    meta = {k: str(v) for k, v in structure.meta.items() if v}
    meta["source"] = source
    op = PutAsset(protein, None, Asset(protein, KIND, protein, structure, meta))
    rows = ligands.ligand_rows(structure, protein)
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

    def work(_job):
        out = []
        for path in files:
            try:
                out.append((path, io.read(path), ""))
            except (io.ProteinFileError, OSError) as err:
                out.append((path, None, str(err)))
        return out

    def done(read) -> None:
        problems = []
        for path, structure, error in read:
            if structure is None:
                problems.append(f"{Path(path).name}: {error}")
            else:
                add_protein(ctx, structure, pdb.name_of(structure, _stem(path)), Path(path).name)
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
        rows = ligands.ligand_rows(store.as_structure(asset.payload), asset.id)
        if rows[ligands.PROTEIN_COLUMN]:
            ctx.add_rows(rows, f"Extract ligands of {asset.id}",
                         structures=(ligands.LIGAND_COLUMN,))  # fmt: skip
            added += len(rows[ligands.PROTEIN_COLUMN])
    ctx.show_message(f"Added {added} ligand row(s)" if added else "No new ligands to add")
    return added


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
