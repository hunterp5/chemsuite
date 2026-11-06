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

"""Chemistry ▸ Standardize Structures."""

from __future__ import annotations

from chemsuite.app.tools import Tool, ToolContext


def _open(ctx: ToolContext):
    from .dialog import StandardizeDialog  # imported on first use: keeps startup fast

    dialog = StandardizeDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


TOOL = Tool(
    id="standardize",
    title="Standardize Structures…",
    menu=("Chemistry",),
    open=_open,
    requires=("rdkit",),
    needs_structures=True,
    tip="Prepare structures for docking: clean up, strip salts, neutralise, and optionally "
    "pick the canonical tautomer and protonate at a pH.",
)
