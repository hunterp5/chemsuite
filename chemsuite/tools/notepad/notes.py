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
"""The document's notes: a document part, saved whether or not a note's pane is the one in
front. Each note is a title and rich text (HTML)."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

#: the document part's key (a saved name: never change it)
KEY = "notepad"
TITLE = "Notepad"


class Notes:
    """Every note of one document. An open pane registers where its live text is
    (:meth:`attach`), so a save takes what is typed now; typing only bumps the revision."""

    def __init__(self, saved: Any = None) -> None:
        self._notes: dict[str, dict[str, str]] = {}
        found = saved.get("notes") if isinstance(saved, dict) else None
        if isinstance(found, dict):
            for note_id, note in found.items():
                if isinstance(note_id, str) and isinstance(note, dict):
                    title, html = note.get("title"), note.get("html")
                    self._notes[note_id] = {
                        "title": title if isinstance(title, str) and title else TITLE,
                        "html": html if isinstance(html, str) else "",
                    }
        self._live: dict[str, Callable[[], str]] = {}
        self._revision = 0

    # -- the document part ---------------------------------------------------------------------
    @property
    def revision(self) -> int:
        return self._revision

    def capture(self) -> Any:
        if not self._notes:
            return {}
        notes = {
            note_id: {"title": note["title"], "html": self.html(note_id)}
            for note_id, note in self._notes.items()
        }
        return {"notes": notes}

    def close(self) -> None:
        self._live.clear()

    # -- notes ---------------------------------------------------------------------------------
    def __contains__(self, note_id: object) -> bool:
        return note_id in self._notes

    @property
    def ids(self) -> list[str]:
        return list(self._notes)

    def new(self, note_id: str | None = None) -> str:
        """A new, empty note (under *note_id* if given) titled Notepad, Notepad 2, …"""
        if note_id is None:
            number = len(self._notes) + 1
            while f"note-{number}" in self._notes:
                number += 1
            note_id = f"note-{number}"
        titles = {note["title"] for note in self._notes.values()}
        title, number = TITLE, 1
        while title in titles:
            number += 1
            title = f"{TITLE} {number}"
        self._notes[note_id] = {"title": title, "html": ""}
        self._revision += 1
        return note_id

    def title(self, note_id: str) -> str:
        return self._notes[note_id]["title"]

    def set_title(self, note_id: str, title: str | None) -> None:
        """Rename a note; empty or None names it Notepad again."""
        title = (title or "").strip() or TITLE
        if self._notes[note_id]["title"] != title:
            self._notes[note_id]["title"] = title
            self._revision += 1

    def html(self, note_id: str) -> str:
        live = self._live.get(note_id)
        return live() if live is not None else self._notes[note_id]["html"]

    def attach(self, note_id: str, live: Callable[[], str]) -> None:
        """A pane shows *note_id*: its text is read from *live* from now on."""
        self._live[note_id] = live

    def detach(self, note_id: str) -> None:
        """The pane closed: keep the text it had."""
        live = self._live.pop(note_id, None)
        if live is not None and note_id in self._notes:
            self._notes[note_id]["html"] = live()

    def is_open(self, note_id: str) -> bool:
        return note_id in self._live

    def touch(self) -> None:
        """A note's text changed (the document needs saving)."""
        self._revision += 1

    def remove(self, note_id: str) -> None:
        self._live.pop(note_id, None)
        if self._notes.pop(note_id, None) is not None:
            self._revision += 1
