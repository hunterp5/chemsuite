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
"""Plot ▸ MedChem ▸ BOILED-Egg…, Golden Triangle…, LipE…"""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _dialog(kind: str):
    def open_dialog(ctx: ToolContext):
        from .dialog import MedChemPlotDialog

        dialog = MedChemPlotDialog(ctx, kind)
        dialog.open()
        return dialog

    return open_dialog


TOOLS = (
    Tool(
        id="boiled_egg",
        title="BOILED-Egg…",
        menu=("Plot", "MedChem"),
        open=_dialog("boiled_egg"),
        requires=("rdkit",),
        needs_structures=True,
        tip="WLOGP against TPSA: predicted gut absorption (white) and brain penetration (yolk).",
    ),
    Tool(
        id="golden_triangle",
        title="Golden Triangle…",
        menu=("Plot", "MedChem"),
        open=_dialog("golden_triangle"),
        requires=("rdkit",),
        needs_structures=True,
        tip="MW against logP, with the region favoured for permeability and clearance.",
    ),
    Tool(
        id="lipe",
        title="Lipophilic Efficiency (LipE)…",
        menu=("Plot", "MedChem"),
        open=_dialog("lipe"),
        requires=("rdkit",),
        needs_structures=True,
        tip="Potency against logP with lines of equal LipE (pActivity − logP).",
    ),
)
