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
"""Analysis ▸ Chemical Space ▸ Chemical Space Map…, Self-Organizing Map…"""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext

MENU = ("Analysis", "Chemical Space")


def _map(ctx: ToolContext):
    from .map_dialog import SpaceMapDialog

    dialog = SpaceMapDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _som(ctx: ToolContext):
    from .som_ui import SomDialog

    dialog = SomDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _som_panel(ctx: ToolContext, state: dict | None = None):
    from .som_ui import SomPanel

    panel = SomPanel(ctx)
    if state:
        panel.set_state(state)
    ctx.add_panel(panel, panel.title())
    return panel


TOOLS = (
    Tool(
        id="space_map",
        title="Space Map…",
        menu=MENU,
        open=_map,
        requires=("sklearn",),
        tip="Place compounds on a 2D map (PCA, UMAP or t-SNE) by fingerprints or properties.",
    ),
    Tool(
        id="som",
        title="Self-Organizing Map…",
        menu=MENU,
        open=_som,
        restore=_som_panel,
        tip="Bin compounds onto a grid where similar compounds share or neighbour a cell.",
    ),
)
