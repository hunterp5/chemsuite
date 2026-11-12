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

from chemsuite.app.document.tools import Tool, ToolContext
from chemsuite.app.shell.jobs import QUICK

#: rows of structures whose cores are worth the worker processes
CORES_POOL_MIN_ROWS = 500


def _open(ctx: ToolContext, kind: str, state: dict | None = None):
    from chemsuite.app.plots.panel import PlotPanel

    def image(oid: int, size):
        column = ctx.default_structure_column()
        return ctx.structure_image(column, oid, size, fit=True) if column else None

    def run(label, work, on_done):  # finding cores ahead: never holds up closing
        return ctx.run(label, work, on_done, lane=QUICK, background=True)

    panel = PlotPanel(
        ctx.table,
        ctx.selection,
        image,
        kind=kind,
        run=run,
        workers=lambda rows: ctx.workers_for(rows, CORES_POOL_MIN_ROWS),
    )
    panel.tool_id = "scatter"  # any plot tool restores any plot kind
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


def _restore(ctx: ToolContext, state: dict):
    return _open(ctx, str(state.get("kind", "scatter")), state)


def _restore_3d(ctx: ToolContext, state: dict):
    """Documents saved while 3D Scatter was a panel of its own."""
    return _open(ctx, "scatter3d", {"kind": "scatter3d", **state})


TOOLS = [
    Tool(
        id="scatter",
        title="Scatter",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "scatter"),
        restore=_restore,
        shortcut="Ctrl+Shift+P",
        tip="Plot one number column against another; select points to select rows.",
    ),
    Tool(
        id="line",
        title="Line",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "line"),
        tip="Points joined in X order, one line per category (e.g. a time course per series).",
    ),
    Tool(
        id="bar",
        title="Bar Chart",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "bar"),
        tip="The mean, median, sum or count of a column per group, with error bars.",
    ),
    Tool(
        id="histogram",
        title="Histogram",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "histogram"),
        tip="Count rows per value range (numbers) or per value (text).",
    ),
    Tool(
        id="box",
        title="Box Plot",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "box"),
        tip="Compare the spread of a number column across groups (e.g. pIC50 by series).",
    ),
    Tool(
        id="density",
        title="Density",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "density"),
        tip="Where tens of thousands of points crowd together: a 2D histogram of X and Y.",
    ),
    Tool(
        id="radar",
        title="Radar",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "radar"),
        tip="Several properties at once on spokes: the selected compounds against the mean.",
    ),
    Tool(
        id="scatter3d",
        title="3D Scatter",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "scatter3d"),
        restore=_restore_3d,
        tip="Three number columns in a cube you can turn, coloured and sized by others.",
    ),
    Tool(
        id="landscape",
        title="3D Landscape",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "landscape"),
        tip="An activity landscape: a surface of one column (pIC50) over two others, such as "
        "a chemical-space map. Peaks next to valleys are activity cliffs.",
        separator_after=True,
    ),
    Tool(
        id="constellation",
        title="Constellation Plot",
        menu=("Data", "Plot"),
        open=lambda ctx: _open(ctx, "constellation"),
        needs_structures=True,
        tip="The cores analog series share, placed in chemical space: sized by their "
        "compounds, coloured by activity (Naveja and Medina-Franco, 2019).",
    ),
]
