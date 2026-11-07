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
"""Analysis ▸ Similarity Search…, Diverse Subset…, Cluster…"""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _dialog(name: str):
    def open_dialog(ctx: ToolContext):
        from . import cluster, diverse, search

        cls = {
            "search": search.SimilarityDialog,
            "diverse": diverse.DiverseSubsetDialog,
            "cluster": cluster.ClusterDialog,
        }[name]
        dialog = cls(ctx)
        ctx.show_dialog(dialog)
        return dialog

    return open_dialog


TOOLS = (
    Tool(
        id="similarity",
        title="Similarity Search…",
        menu=("Analysis",),
        open=_dialog("search"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Similarity of every row to a structure, or to its nearest selected compound.",
    ),
    Tool(
        id="diverse_subset",
        title="Diverse Subset…",
        menu=("Analysis",),
        open=_dialog("diverse"),
        requires=("rdkit",),
        tip="Pick compounds that cover the chemistry (MaxMin), optionally adding to a selection.",
    ),
    Tool(
        id="cluster",
        title="Cluster…",
        menu=("Analysis",),
        open=_dialog("cluster"),
        requires=("rdkit", "sklearn", "scipy"),
        tip="Group compounds by Butina, k-means or hierarchical clustering.",
        separator_after=True,
    ),
)
