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
"""Tools ▸ Notepad: each use opens a new note in a pane of its own; the notes are saved with
the document (a document part), so typing marks it as needing a save."""

from __future__ import annotations

from chemsuite.app.document.tools import Tool, ToolContext

from .notes import KEY, Notes

TOOL_ID = "notepad"


def open_notepad(ctx: ToolContext, note_id: str | None = None):
    """A pane for *note_id* (made if the document lacks it). Without one: a note the document
    has but shows nowhere (its pane was lost), else a new note."""
    from .panel import NotepadPanel  # imported on first use

    notes = ctx.document_part(KEY, Notes)
    if note_id is None:
        note_id = next((n for n in notes.ids if not notes.is_open(n)), None) or notes.new()
    elif note_id not in notes:
        notes.new(note_id)
    panel = NotepadPanel(notes, note_id)
    ctx.add_panel(panel, panel.title())
    panel.editor.setFocus()
    return panel


TOOL = Tool(
    id=TOOL_ID,
    title="Notepad",
    menu=("Tools",),
    open=lambda ctx: open_notepad(ctx),
    restore=lambda ctx, state: open_notepad(ctx, str(state.get("note") or "") or None),
    tip="Notes with fonts, colors, bullets and numbering, kept with the document; each "
    "opens in a pane of its own.",
)
