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
"""Chemistry ▸ 3D Viewer: the selected rows' structures in 3D, in a panel."""

from __future__ import annotations

from chemsuite.app.jobs import QUICK
from chemsuite.app.tools import Tool, ToolContext


def _open(ctx: ToolContext, state: dict | None = None):
    from chemsuite.app.molview.panel import MolViewPanel  # imported on first use

    def run(label, work, on_done):
        return ctx.run(label, work, on_done, lane=QUICK)

    panel = MolViewPanel(
        ctx.table,
        ctx.selection,
        current=ctx.current_oid,
        run=run,
        move_to=ctx.set_current_oid,
        watch_current=ctx.watch_current,
    )
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOL = Tool(
    id="molview",
    title="3D Viewer",
    menu=("Chemistry",),
    open=_open,
    restore=_open,
    requires=("rdkit",),
    needs_structures=True,
    shortcut="Ctrl+Shift+3",
    tip="Browse the table's structures in 3D: step through or overlay conformers, measure.",
)
