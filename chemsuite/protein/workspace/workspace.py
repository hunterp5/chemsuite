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
"""A document's protein workspace: its session and how its view was left, saved with it.

There is one per document, whether or not the viewer's window is open:

- Opening a document keeps the saved state as it is (no protein is decoded or drawn) until the
  viewer first asks for the :meth:`ProteinWorkspace.session`; a document saved before then
  saves that state back untouched.
- Closing the viewer keeps the session, so the window opens again as it was left.
- :attr:`ProteinWorkspace.revision` changes with whatever is worth saving, so the document
  knows it needs saving.

The saved form (:meth:`capture`) is versioned: ``{"version": 1, "session": …, "view": …}``,
``view`` holding the background, depth cue and camera. A state from a newer chemsuite is kept
as it was and saved back unchanged rather than half read.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from typing import Any

from chemsuite.core.assets import AssetHost

from .session import WorkbenchSession

log = logging.getLogger(__name__)

VERSION = 1
#: the view settings a workspace keeps for its window (camera moves are not changes to save)
VIEW_KEYS = ("background", "depth_cue", "camera")
SAVED_VIEW_KEYS = ("background", "depth_cue")


class ProteinWorkspace:
    def __init__(
        self,
        store: AssetHost,
        saved: Any = None,
        *,
        deferred: bool = True,
        relabel: Callable[[str, Mapping[str, str]], Sequence[object]] | None = None,
        count_named: Callable[[str, Any], int] | None = None,
    ) -> None:
        self.store = store
        self._deferred = deferred
        self._relabel = relabel
        self._count_named = count_named
        self._session: WorkbenchSession | None = None
        self._pending: dict = {}  # the session's saved state, until the session is made
        #: background, depth cue and camera, as the viewer left them
        self.view: dict = {}
        #: while the viewer is open: its view settings now
        self.view_source: Callable[[], dict] | None = None
        self._newer: Any = None  # a state from a newer chemsuite, kept as it was
        self._changes = 0  # changes not counted by the session (the view, an adopted state)
        #: nothing was saved for it (an older document may hand over its viewer's state)
        self.fresh = saved is None
        self._load(saved)

    def _load(self, saved: Any) -> None:
        if saved is None:
            return
        if not isinstance(saved, dict) or not isinstance(saved.get("version"), int):
            log.warning("ignored a damaged protein workspace")
            return
        if saved["version"] > VERSION:
            log.warning("kept a protein workspace from a newer chemsuite as it was")
            self._newer = saved
            return
        session = saved.get("session")
        self._pending = dict(session) if isinstance(session, dict) else {}
        view = saved.get("view")
        self.view = {k: view[k] for k in VIEW_KEYS if isinstance(view, dict) and k in view}

    # -- the session --------------------------------------------------------------------------
    @property
    def started(self) -> bool:
        return self._session is not None

    def session(self) -> WorkbenchSession:
        """The workspace's session, made (and its saved state put back) when first asked for."""
        if self._session is None:
            session = WorkbenchSession(
                self.store,
                deferred=self._deferred,
                relabel=self._relabel,
                count_named=self._count_named,
            )
            if self._pending:
                session.set_state(self._pending)
            session.revision = 0
            self._session = session
            self._pending = {}
        return self._session

    def adopt(self, viewer_state: Mapping) -> bool:
        """Take the state an older document saved with its viewer window (its ``session``,
        background, depth cue and camera); only when the document saved no workspace."""
        if not self.fresh:
            return False
        session = viewer_state.get("session")
        if isinstance(session, dict):
            if self._session is not None:
                self._session.set_state(session)
            else:
                self._pending = dict(session)
        self.view = {k: viewer_state[k] for k in VIEW_KEYS if k in viewer_state}
        self.fresh = False
        return True

    def set_view(self, view: Mapping) -> None:
        """The viewer's view settings now (when it closes)."""
        new = {k: view[k] for k in VIEW_KEYS if k in view}
        if any(new.get(k) != self.view.get(k) for k in SAVED_VIEW_KEYS):
            self._changes += 1
        self.view = new

    # -- saving -------------------------------------------------------------------------------
    @property
    def revision(self) -> int:
        """Changes since the document opened that are worth saving."""
        session = self._session.revision if self._session is not None else 0
        live = self._live_view()
        moved = any(live.get(k) != self.view.get(k) for k in SAVED_VIEW_KEYS)
        return self._changes + session + (1 if moved else 0)

    def _live_view(self) -> dict:
        if self.view_source is None:
            return self.view
        return {k: v for k, v in self.view_source().items() if k in VIEW_KEYS}

    def capture(self) -> Any:
        """What the document saves (None: nothing worth saving)."""
        if self._newer is not None:
            return self._newer
        session = self._session.state() if self._session is not None else self._pending
        view = self._live_view()
        if not session and not view:
            return None
        return {"version": VERSION, "session": session, "view": view}

    def close(self) -> None:
        """The document is closing."""
        self.view_source = None
        if self._session is not None:
            self._session.close()
