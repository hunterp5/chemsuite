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

"""Chemistry ▸ Retrosynthesis, and the Routes pane it fills (docs/RETRO_PLAN.md)."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

ROUTES_TOOL_ID = "retro_routes"


def _search(ctx: ToolContext):
    from .dialog import RetrosynthesisDialog  # imported on first use: keeps startup fast

    return ctx.show_dialog(RetrosynthesisDialog(ctx))


def _routes(ctx: ToolContext, state: dict | None = None):
    """Open (or reopen) the Routes pane. Saved routes come back through *state*."""
    from .routes_panel import RoutesPanel

    panel = RoutesPanel(ctx)
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOLS = (
    Tool(
        id="retrosynthesis",
        title="Retrosynthesis…",
        menu=("Chemistry", "Reaction"),
        open=_search,
        requires=("rdkit",),
        needs_structures=True,
        tip="Plan synthetic routes to each structure with a local AiZynthFinder install, "
        "and show the routes in the Routes pane.",
    ),
    Tool(
        id=ROUTES_TOOL_ID,
        title="Routes",
        menu=("Chemistry",),
        open=_routes,
        restore=_routes,
        requires=("rdkit",),
        tip="The retrosynthetic routes found in this document.",
    ),
)
