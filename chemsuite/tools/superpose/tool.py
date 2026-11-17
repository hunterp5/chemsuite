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
"""Chemistry ▸ 3D ▸ Superpose…"""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

MENU = ("Chemistry", "3D")


def _open(ctx: ToolContext):
    from .dialog import SuperposeDialog  # imported on first use: keeps startup fast

    dialog = SuperposeDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


TOOLS = (
    Tool(
        id="superpose",
        title="Superpose…",
        menu=MENU,
        open=_open,
        requires=("rdkit",),
        needs_structures=True,
        tip="Superpose molecules onto a reference (a row or a protein's ligand) by a common "
        "substructure, by atoms, by shape, by shape and features, or by field; or each row's "
        "conformers onto its first.",
    ),
)
