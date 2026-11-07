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

"""Chemistry ▸ Reactions ▸ Reaction Search, Split Reaction and Enumerate Products."""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _search(ctx: ToolContext):
    from .search import ReactionSearchDialog  # imported on first use: keeps startup fast

    return ctx.show_dialog(ReactionSearchDialog(ctx))


def _split(ctx: ToolContext):
    from .split import SplitReactionDialog

    return ctx.show_dialog(SplitReactionDialog(ctx))


def _enumerate(ctx: ToolContext):
    from .enumerate import EnumerateDialog

    return ctx.show_dialog(EnumerateDialog(ctx))


TOOLS = (
    Tool(
        id="reaction-search",
        title="Reaction Search…",
        menu=("Chemistry", "Reactions"),
        open=_search,
        requires=("rdkit",),
        needs_reactions=True,
        tip="Find reactions that contain a drawn reaction, or a structure as a reactant, "
        "product or agent.",
    ),
    Tool(
        id="split-reaction",
        title="Split Reaction…",
        menu=("Chemistry", "Reactions"),
        open=_split,
        requires=("rdkit",),
        needs_reactions=True,
        tip="The reactants, agents and products of each reaction as structure columns.",
    ),
    Tool(
        id="enumerate-products",
        title="Enumerate Products…",
        menu=("Chemistry", "Reactions"),
        open=_enumerate,
        requires=("rdkit",),
        needs_structures=True,
        tip="Apply a reaction template to building blocks: row by row, or every combination "
        "(a virtual library).",
    ),
)
