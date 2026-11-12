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

"""Tools ▸ Sketcher: a blank sketcher; Add to Table appends the drawing as a new row."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext


def _open(ctx: ToolContext):
    from chemsuite.app.sketcher import SketcherDialog

    dialog = SketcherDialog("", "Sketcher", ctx.parent, ok_text="Add to Table")

    def accepted() -> None:
        if dialog.molblock and ctx.add_compound(dialog.molblock) is not None:
            ctx.show_message("Added the drawing as a new row")

    dialog.accepted.connect(accepted)
    ctx.show_dialog(dialog)
    return dialog


TOOL = Tool(
    id="sketcher",
    title="Sketcher…",
    menu=("Tools",),
    open=_open,
    requires=("rdkit",),
    tip="Draw a structure and add it to the table as a new row.",
    separator_after=True,
)
