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
"""The viewer's workspace as a whole, and each structure in it, on the document's side:

- File ▸ Open… (structure files and session files), Save Session…, Clear Session…;
- a structure's Rename…, Duplicate and Delete… (the list's right-click menu).

Each change to the document's structures is one undo step. Session files themselves are
:mod:`chemsuite.macromolecule.workspace.file`.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.assets import PutAsset
from chemsuite.macromolecule.workspace import file as session_file
from chemsuite.macromolecule.workspace.session import DERIVED_FROM

from .tool import (
    DELETE_ROWS,
    FILE_FILTER,
    KEEP_ROWS,
    KIND,
    TOOL_ID,
    dialog_parent,
    linked_rows,
    macromolecule_workspace,
    open_protein,
    open_viewer,
)

SESSION_FILTER = "Sessions (*.cspx)"
#: how a session file joins structures already open
REPLACE, ADD = "replace", "add"


def structures(ctx: ToolContext) -> dict:
    """The document's structures, by id, in order."""
    return {id: a for id, a in ctx.table.assets.items() if a.kind == KIND}


# -- File ▸ Open… --------------------------------------------------------------------------------
def open_files(ctx: ToolContext, paths: list[str] | None = None):
    """Open structure files and session files (asked for when *paths* is None)."""
    if paths is None:
        from PySide6.QtWidgets import QFileDialog

        both = f"Structures and sessions ({_patterns(FILE_FILTER)} *.cspx)"
        paths, _ = QFileDialog.getOpenFileNames(
            dialog_parent(ctx), "Open", "", f"{both};;{FILE_FILTER};;{SESSION_FILTER}"
        )
    sessions = [p for p in paths if session_file.is_session(p)]
    others = [p for p in paths if not session_file.is_session(p)]
    for path in sessions:
        open_session(ctx, path)
    return open_protein(ctx, others) if others else None


def _patterns(dialog_filter: str) -> str:
    return dialog_filter[dialog_filter.index("(") + 1 : dialog_filter.rindex(")")]


def open_session(ctx: ToolContext, path: str, how: str | None = None):
    """Read the session file at *path* in the background, then add its structures and how
    they looked: in place of the structures open (*how* REPLACE), beside them (ADD), or as
    asked when there are some (None)."""

    def work(job) -> session_file.Session | str:
        try:
            return session_file.read(path)
        except (session_file.SessionFileError, OSError) as error:
            return str(error)  # a file that is not one: said, not a crash box

    def done(found: session_file.Session | str) -> None:
        name = Path(path).name
        if isinstance(found, str):
            ctx.show_message(f"Could not open {name}: {found}")
            return
        chosen = how
        if chosen is None:
            chosen = _ask_how(ctx, name) if structures(ctx) else REPLACE
        if chosen is None:
            return
        take_session(ctx, found, chosen, name)

    return ctx.run(f"Opening {Path(path).name}…", work, done)


def take_session(ctx: ToolContext, found: session_file.Session, how: str, name: str) -> None:
    """Add a session's structures (one undo step) and its state to the workspace."""
    workspace = macromolecule_workspace(ctx)
    ops = []
    if how == REPLACE:
        ops = [PutAsset(id, asset, None) for id, asset in structures(ctx).items()]
        workspace.clear()
        taken = {id for id in ctx.table.assets if id not in structures(ctx)}
    else:
        taken = set(ctx.table.assets)

    def unique(wanted: str) -> str:
        found_id, n = wanted, 2
        while found_id in taken:
            found_id, n = f"{wanted}-{n}", n + 1
        return found_id

    incoming = [a for a in found.assets if a.kind == KIND]
    ids = session_file.new_ids([a.id for a in incoming], unique)
    assets = session_file.remap_assets(incoming, ids)
    ops += [PutAsset(a.id, None, a) for a in assets]
    ctx.edit(f"Open {name}", ops)
    if found.state is not None:
        workspace.take(session_file.remap(found.state, ids), replace=how == REPLACE)
    panel = open_viewer(ctx)
    if how == REPLACE and workspace.view:
        panel.set_view_state(workspace.view)
    count = len(assets)
    ctx.show_message(f"Opened {name}: {count} structure{'s' if count != 1 else ''}")


def _ask_how(ctx: ToolContext, name: str) -> str | None:
    """Replace the structures open, or add the session's beside them (None: cancel)."""
    from PySide6.QtWidgets import QMessageBox

    box = QMessageBox(
        QMessageBox.Icon.Question,
        "Open Session",
        f"Open {name} in place of the structures here, or beside them?",
        parent=dialog_parent(ctx),
    )
    replacing = box.addButton("Replace", QMessageBox.ButtonRole.DestructiveRole)
    adding = box.addButton("Add", QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(adding)
    box.exec()
    clicked = box.clickedButton()
    return REPLACE if clicked is replacing else ADD if clicked is adding else None


# -- File ▸ Save Session… ------------------------------------------------------------------------
def save_session(ctx: ToolContext, path: str | None = None):
    """Write every structure and the workspace's state to a session file, in the background.
    Returns the job, or None when there is nothing to save or the dialog was cancelled."""
    assets = list(structures(ctx).values())
    if not assets:
        return None
    if path is None:
        from PySide6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getSaveFileName(dialog_parent(ctx), "Save Session",
                                              "session.cspx", SESSION_FILTER)  # fmt: skip
    if not path:
        return None
    if not session_file.is_session(path):
        path += session_file.SUFFIX
    state = macromolecule_workspace(ctx).capture()

    def work(job):
        session_file.save(path, assets, state, job.cancel_token)
        return path

    def done(saved: str) -> None:
        ctx.show_message(f"Saved session {Path(saved).name}")

    return ctx.run(f"Saving {Path(path).name}…", work, done)


# -- File ▸ Clear Session… -----------------------------------------------------------------------
def clear_session(ctx: ToolContext, rows: str | None = None) -> bool:
    """Remove every structure, and their ligand rows unless *rows* is KEEP_ROWS (asked once
    when None and there are some), in one undo step; the focus, named selections, labels and
    measurements go too. False when cancelled or there is nothing to clear."""
    here = structures(ctx)
    if not here:
        return False
    linked = np.unique(np.concatenate([linked_rows(ctx, id) for id in here]))
    if rows is None and len(linked):
        rows = _ask_rows(ctx, len(here), len(linked))
        if rows is None:
            return False
    deleting = rows != KEEP_ROWS
    macromolecule_workspace(ctx).clear()
    removing = [PutAsset(id, asset, None) for id, asset in here.items()]
    ctx.remove_rows(linked if deleting else [], "Clear session", extra_ops=removing)
    gone = f" and {len(linked)} ligand row(s)" if deleting and len(linked) else ""
    ctx.show_message(f"Cleared {len(here)} structure(s){gone}")
    return True


def _ask_rows(ctx: ToolContext, count: int, rows: int) -> str | None:
    from PySide6.QtWidgets import QMessageBox

    box = QMessageBox(
        QMessageBox.Icon.Question,
        "Clear Session",
        f"Remove {count} structure(s) and their {rows} ligand row(s)?",
        parent=dialog_parent(ctx),
    )
    delete = box.addButton("Delete Rows", QMessageBox.ButtonRole.DestructiveRole)
    keep = box.addButton("Keep Rows", QMessageBox.ButtonRole.AcceptRole)
    box.addButton(QMessageBox.StandardButton.Cancel)
    box.setDefaultButton(delete)
    box.exec()
    clicked = box.clickedButton()
    return DELETE_ROWS if clicked is delete else KEEP_ROWS if clicked is keep else None


# -- a structure's Rename… and Duplicate ---------------------------------------------------------
def rename_structure(ctx: ToolContext, id: str, name: str | None = None) -> bool:
    """Give structure *id* a new name (asked for when None), one undo step. The id, which the
    table's rows link to, stays."""
    asset = structures(ctx).get(id)
    if asset is None:
        return False
    if name is None:
        from PySide6.QtWidgets import QInputDialog

        name, ok = QInputDialog.getText(dialog_parent(ctx), "Rename", "Name", text=asset.name)
        if not ok:
            return False
    name = name.strip()
    if not name or name == asset.name:
        return False
    ctx.edit(f"Rename {asset.name}", [PutAsset(id, asset, replace(asset, name=name))])
    return True


def duplicate_structure(ctx: ToolContext, id: str) -> str:
    """Add a copy of structure *id* that looks as it does (one undo step), with its
    preparation record. Returns the copy's id ("" when there is no such structure)."""
    asset = structures(ctx).get(id)
    if asset is None:
        return ""
    new = ctx.table.unique_asset_id(f"{asset.name}-copy")
    workspace = macromolecule_workspace(ctx)
    meta = {k: v for k, v in asset.meta.items() if k != DERIVED_FROM}
    workspace.session().copy_look(id, new)
    ctx.edit(f"Duplicate {asset.name}",
             [PutAsset(new, None, replace(asset, id=new, name=new, meta=meta))])  # fmt: skip
    record = workspace.prepared(id)
    if record is not None:
        workspace.set_prepared(replace(record, protein=new))
    return new


__all__ = [
    "ADD",
    "REPLACE",
    "TOOL_ID",
    "clear_session",
    "duplicate_structure",
    "open_files",
    "open_session",
    "rename_structure",
    "save_session",
]
