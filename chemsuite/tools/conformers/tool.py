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
"""Chemistry ▸ Generate Conformers… and Export Conformers to SDF…"""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext


def _open_generate(ctx: ToolContext):
    from .dialog import ConformerDialog  # imported on first use: keeps startup fast

    dialog = ConformerDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


def _open_export(ctx: ToolContext):
    from .dialog import ExportConformersDialog

    dialog = ExportConformersDialog(ctx)
    ctx.show_dialog(dialog)
    return dialog


TOOLS = (
    Tool(
        id="conformers",
        title="Generate Conformers…",
        menu=("Chemistry",),
        open=_open_generate,
        requires=("rdkit",),
        needs_structures=True,
        tip="3D conformers with ETKDG and MMFF94, stored lowest energy first in a new column.",
        separator_after=True,
    ),
    Tool(
        id="export_conformers",
        title="Conformers (SDF)…",
        menu=("File", "Export"),
        open=_open_export,
        requires=("rdkit",),
        needs_structures=True,
        tip="Write every conformer as its own SD record, with the row's other columns.",
    ),
)
