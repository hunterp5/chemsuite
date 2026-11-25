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
"""Data ▸ Dimensionality Reduction ▸ PCA…, UMAP…, t-SNE…, PaCMAP… | Self-Organizing Map…"""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

MAP_MENU = ("Data", "Dimensionality Reduction")
#: the methods in the submenu, in order: title (``embedding.METHODS``), the packages each
#: needs beyond scikit-learn, and its tooltip
MAPS = {
    "pca": ("PCA", (), "Fast and exact; distances along the axes mean something, but two "
            "components often cannot separate every group."),
    "umap": ("UMAP", ("umap",), "Keeps neighborhoods together, so series form islands; "
             "distances between islands mean little. Above 50,000 rows it is fitted on a "
             "sample."),
    "tsne": ("t-SNE", (), "Like UMAP, good at local groups; slower, and only for the rows "
             "given."),
    "pacmap": ("PaCMAP", ("pacmap",), "Keeps local groups and the overall layout together, so "
               "distances between islands mean more than in UMAP or t-SNE; fast on large "
               "sets."),
}  # fmt: skip


def _map(method: str):
    def open_map(ctx: ToolContext):
        from .map_dialog import SpaceMapDialog

        dialog = SpaceMapDialog(ctx, method)
        ctx.show_dialog(dialog)
        return dialog

    return open_map


def _loadings(ctx: ToolContext, state: dict):
    from .loadings import restore

    return restore(ctx, state)


def _som(ctx: ToolContext):
    from .som_ui import SomDialog

    dialog = SomDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _gtm(ctx: ToolContext):
    from .gtm_ui import GtmDialog

    dialog = GtmDialog(ctx)
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
    *(
        Tool(
            id=f"space_map_{method}",
            title=f"{title}…",
            menu=MAP_MENU,
            open=_map(method),
            requires=("sklearn", *packages),
            tip=f"Map compounds to 2D or 3D by fingerprints, properties or both. {tip}",
            separator_after=method == "pacmap",  # the Self-Organizing Map follows
        )
        for method, (title, packages, tip) in MAPS.items()
    ),
    Tool(
        id="pca_loadings",
        title="PCA Loadings",
        menu=(),  # a PCA map opens it; this entry reopens it with the document
        open=lambda ctx: None,
        restore=_loadings,
        tip="What each number column contributes to each principal component.",
    ),
    Tool(
        id="som",
        title="Self-Organizing Map…",
        menu=MAP_MENU,
        open=_som,
        restore=_som_panel,
        tip="Bin compounds onto a grid where similar compounds share or neighbor a cell.",
    ),
    Tool(
        id="gtm_map",
        title="Generative Topographic Map…",
        menu=MAP_MENU,
        open=_gtm,
        tip="A smooth map where each compound is spread over the nodes it fits: landscapes "
        "of activity or class, and how well each compound fits the map (Bishop et al., 1998).",
    ),
)
