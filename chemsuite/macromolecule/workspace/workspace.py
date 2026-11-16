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
  viewer first asks for the :meth:`MacromoleculeWorkspace.session`; a document saved before then
  saves that state back untouched.
- Closing the viewer keeps the session, so the window opens again as it was left.
- :attr:`MacromoleculeWorkspace.revision` changes with whatever is worth saving, so the document
  knows it needs saving.

The saved form (:meth:`capture`) is versioned: ``{"version": 1, "session": …, "view": …,
"reports": …, "poses": …, "reviews": …, "prepared": …}``, ``view`` holding the background,
depth cue and camera, ``reports`` each protein's last preparation report (:mod:`.reports`),
``poses`` each protein's docking runs (:mod:`.poses`), ``reviews`` each preparation still
being reviewed (:mod:`.reviews`) and ``prepared`` what each preparation left besides the
structure (:mod:`.prepared`; the last four were added later, so each is optional and an older
chemsuite simply ignores it). A state from a newer chemsuite is kept as it was and saved back
unchanged rather than half read.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import replace
from typing import Any

from chemsuite.core.assets import AssetHost

from .poses import DockedPoses, runs_from_json
from .prepared import PreparedRecord
from .reports import PrepareReport
from .reviews import Review
from .session import WorkbenchSession

log = logging.getLogger(__name__)

VERSION = 1
#: the view settings a workspace keeps for its window (camera moves are not changes to save)
VIEW_KEYS = ("background", "depth_cue", "camera")
SAVED_VIEW_KEYS = ("background", "depth_cue")


class MacromoleculeWorkspace:
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
        #: each protein's last preparation report, by protein id
        self.reports: dict[str, PrepareReport] = {}
        #: each protein's docking runs, oldest first, by protein id
        self.poses: dict[str, list[DockedPoses]] = {}
        #: each preparation being reviewed and not yet applied, by protein id
        self.reviews: dict[str, Review] = {}
        #: what each protein's preparation (and minimisation) left beside it, by protein id
        self.records: dict[str, PreparedRecord] = {}
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
        for protein, report in dict(saved.get("reports") or {}).items():
            made = PrepareReport.from_json(report)
            if made is not None and isinstance(protein, str):
                self.reports[protein] = made
        for protein, runs in dict(saved.get("poses") or {}).items():
            made = runs_from_json(runs)
            if made and isinstance(protein, str):
                self.poses[protein] = made
        for protein, pending in dict(saved.get("reviews") or {}).items():
            review = Review.from_json(pending)
            if review is not None and isinstance(protein, str):
                self.reviews[protein] = review
        for protein, record in dict(saved.get("prepared") or {}).items():
            made = PreparedRecord.from_json(record)
            if made is not None and isinstance(protein, str):
                self.records[protein] = made

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

    # -- preparation reports ------------------------------------------------------------------
    def report(self, protein: str) -> PrepareReport | None:
        """What the last preparation of *protein* did, if it was prepared."""
        return self.reports.get(protein)

    def set_report(self, report: PrepareReport) -> None:
        """Keep what a preparation did (replacing that protein's last report)."""
        self.reports[report.protein] = report
        self._changes += 1

    def drop_report(self, protein: str) -> bool:
        """Forget *protein*'s report (its preparation was undone, or the protein went)."""
        if self.reports.pop(protein, None) is None:
            return False
        self._changes += 1
        return True

    # -- what a preparation left besides the structure ----------------------------------------
    def prepared(self, protein: str) -> PreparedRecord | None:
        """The record of *protein*'s last preparation (or minimisation), if it has one. Check
        :meth:`PreparedRecord.matches` before using it: the protein may have been edited."""
        return self.records.get(protein)

    def set_prepared(self, record: PreparedRecord) -> None:
        self.records[record.protein] = record
        self._changes += 1

    def drop_prepared(self, protein: str) -> bool:
        if self.records.pop(protein, None) is None:
            return False
        self._changes += 1
        return True

    # -- docked poses -------------------------------------------------------------------------
    def runs(self, protein: str) -> tuple[DockedPoses, ...]:
        """*protein*'s docking runs, oldest first."""
        return tuple(self.poses.get(protein, ()))

    def docked(self, protein: str, number: int | None = None) -> DockedPoses | None:
        """*protein*'s run numbered *number* (None: its latest), if there is one."""
        runs = self.poses.get(protein, [])
        if number is None:
            return runs[-1] if runs else None
        return next((run for run in runs if run.number == number), None)

    def add_docked(self, run: DockedPoses) -> DockedPoses:
        """Keep a docking run beside *run.protein*'s earlier ones; returns it numbered."""
        runs = self.poses.setdefault(run.protein, [])
        run = replace(run, number=max((kept.number for kept in runs), default=0) + 1)
        runs.append(run)
        self._changes += 1
        return run

    def drop_docked(self, protein: str, number: int | None = None) -> bool:
        """Forget *protein*'s run numbered *number* (None: all of its runs)."""
        runs = self.poses.get(protein, [])
        kept = [run for run in runs if number is not None and run.number != number]
        if len(kept) == len(runs):
            return False
        if kept:
            self.poses[protein] = kept
        else:
            del self.poses[protein]
        self._changes += 1
        return True

    # -- preparations being reviewed ----------------------------------------------------------
    def review(self, protein: str) -> Review | None:
        """The review of *protein*'s preparation still open, if there is one."""
        return self.reviews.get(protein)

    def set_review(self, review: Review) -> None:
        """Keep what has been decided so far (it is saved with the document)."""
        if self.reviews.get(review.protein) == review:
            return
        self.reviews[review.protein] = review
        self._changes += 1

    def drop_review(self, protein: str) -> bool:
        """Forget *protein*'s review (applied, or abandoned)."""
        if self.reviews.pop(protein, None) is None:
            return False
        self._changes += 1
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
        here = self.store.assets
        kept = {id: report for id, report in self.reports.items() if id in here}
        docked = {id: run for id, run in self.poses.items() if id in here}
        pending = {id: review for id, review in self.reviews.items() if id in here}
        records = {id: record for id, record in self.records.items() if id in here}
        if not (session or view or kept or docked or pending or records):
            return None
        saved = {"version": VERSION, "session": session, "view": view}
        if kept:
            saved["reports"] = {id: report.to_json() for id, report in kept.items()}
        if docked:
            saved["poses"] = {id: [run.to_json() for run in runs] for id, runs in docked.items()}
        if pending:
            saved["reviews"] = {id: review.to_json() for id, review in pending.items()}
        if records:
            saved["prepared"] = {id: record.to_json() for id, record in records.items()}
        return saved

    def close(self) -> None:
        """The document is closing."""
        self.view_source = None
        if self._session is not None:
            self._session.close()
