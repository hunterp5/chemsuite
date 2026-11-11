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
"""Chemistry ▸ Structure Viewer and Data ▸ Browser: the table's structures one row at a time,
in 2D or 3D, in a panel."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext
from chemsuite.app.shell.jobs import QUICK

TOOL_ID = "viewer"


def open_viewer(
    ctx: ToolContext, mode: str | None, state: dict | None = None, at: int | None = None
):
    """Bring the open viewer forward in *mode* (None: as it is; a new one opens in 2D), or
    open one (with *state*, when reopening a document); *at* first makes that row current."""
    from chemsuite.app.viewer.panel import StructureViewer  # imported on first use

    if at is not None:
        ctx.set_current_oid(at)
    if state is None:
        panel = ctx.find_panel(TOOL_ID)
        if panel is not None:
            if mode is not None:
                panel.set_mode(mode)
            panel.setFocus()
            return panel

    def run(label, work, on_done):  # 3D made for the view: never holds up closing
        return ctx.run(label, work, on_done, lane=QUICK, background=True)

    panel = StructureViewer(
        ctx.table,
        ctx.selection,
        current=ctx.current_oid,
        run=run,
        move_to=ctx.set_current_oid,
        watch_current=ctx.watch_current,
        mode=mode or "2d",
    )
    if state is not None:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    panel.setFocus()
    return panel


def _restore(ctx: ToolContext, state: dict):
    return open_viewer(ctx, state.get("mode", "3d"), state)


TOOL = Tool(
    id=TOOL_ID,
    title="Structure Viewer",
    menu=("Chemistry",),
    open=lambda ctx: open_viewer(ctx, "3d"),
    restore=_restore,
    aliases=("molview",),  # the 3D Viewer, before it drew in 2D too
    requires=("rdkit",),
    needs_structures=True,
    shortcut="Ctrl+Shift+3",
    tip="Browse the table's structures one row at a time, in 2D or 3D: step through or "
    "overlay conformers, measure, show column values over the structure.",
    # Data ▸ Browser opens it in 2D; a cell's Browse brings it to that row, as it is
    open_at=lambda ctx, oid: open_viewer(ctx, "2d" if oid is None else None, at=oid),
)
