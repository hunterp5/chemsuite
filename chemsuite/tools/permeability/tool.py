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

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext


def _open(ctx: ToolContext):
    from .dialog import PermeabilityDialog  # imported on first use: keeps startup fast

    dialog = PermeabilityDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


TOOL = Tool(
    id="predict_permeability",
    title="Permeability…",
    menu=("Chemistry", "Calculate"),
    open=_open,
    requires=("rdkit",),
    needs_structures=True,
    tip="Predict Caco-2 permeability, MDR1-MDCK efflux and PAMPA, each with an interval.",
)
