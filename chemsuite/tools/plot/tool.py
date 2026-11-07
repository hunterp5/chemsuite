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
"""Registration of the plot tools."""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _open(ctx: ToolContext, kind: str, state: dict | None = None):
    from chemsuite.app.plots.panel import PlotPanel

    def image(oid: int, size):
        column = ctx.default_structure_column()
        return ctx.structure_image(column, oid, size, fit=True) if column else None

    panel = PlotPanel(ctx.table, ctx.selection, image, kind=kind)
    panel.tool_id = "scatter"  # any plot tool restores any plot kind
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


def _restore(ctx: ToolContext, state: dict):
    return _open(ctx, str(state.get("kind", "scatter")), state)


TOOLS = [
    Tool(
        id="scatter",
        title="New Scatter Plot",
        menu=("Plot",),
        open=lambda ctx: _open(ctx, "scatter"),
        restore=_restore,
        shortcut="Ctrl+Shift+P",
        tip="Plot one number column against another; select points to select rows.",
    ),
    Tool(
        id="histogram",
        title="New Histogram",
        menu=("Plot",),
        open=lambda ctx: _open(ctx, "histogram"),
        tip="Count rows per value range (numbers) or per value (text).",
    ),
    Tool(
        id="box",
        title="New Box Plot",
        menu=("Plot",),
        open=lambda ctx: _open(ctx, "box"),
        tip="Compare the spread of a number column across groups (e.g. pIC50 by series).",
    ),
]
