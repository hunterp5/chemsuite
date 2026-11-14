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

"""Chemistry ▸ Decompose (Fragments, Central Ring) and Enumerate (R-Groups on a Core)."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

MENU = ("Chemistry", "Decompose")
BUILD = ("Chemistry", "Enumerate")


def _fragments(ctx: ToolContext):
    from .decompose import FragmentsDialog  # imported on first use: keeps startup fast

    return ctx.show_dialog(FragmentsDialog(ctx))


def _central_ring(ctx: ToolContext):
    from .central_ring import CentralRingDialog

    return ctx.show_dialog(CentralRingDialog(ctx))


def _rgroup_library(ctx: ToolContext):
    from .rgroups_library import RGroupLibraryDialog

    return ctx.show_dialog(RGroupLibraryDialog(ctx))


TOOLS = (
    Tool(
        id="fragments",
        title="Fragments…",
        menu=MENU,
        open=_fragments,
        requires=("rdkit",),
        needs_structures=True,
        tip="Cut each compound at its BRICS or RECAP bonds, into fragment columns and a "
        "panel of the distinct fragments.",
    ),
    Tool(
        id="central-ring",
        title="Central Ring…",
        menu=MENU,
        open=_central_ring,
        requires=("rdkit",),
        needs_structures=True,
        tip="R-groups around each compound's own most central ring system, with no core to enter.",
    ),
    Tool(
        id="rgroup-library",
        title="R-Groups on a Core…",
        menu=BUILD,
        open=_rgroup_library,
        requires=("rdkit",),
        tip="Every combination of the substituents you choose on a numbered core, as a "
        "virtual library.",
    ),
)
