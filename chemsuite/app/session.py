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
"""The working session saved with a document: everything about the window beside the table.

A saved document reopens exactly as it was left. The table file itself holds the data, column
order and kinds, styles and the sort; :data:`PARTS` holds the rest, one entry per piece of
window state, each saved under its own key in the document's extras:

* ``filters`` – filter cards and filter sets
* ``advanced_search`` – the advanced query, its recent searches, and whether it is applied
* ``panels`` – tool panels (plots, statistics, maps, SAR results…), via each tool's ``restore``
* ``view`` – selected rows, the current cell and the scroll position
* ``layout`` – which panes are open and where, splitter sizes, column widths, the search bar

**Adding a feature with state of its own?** Add a :class:`Part` here (or, for a tool panel,
give the panel ``tool_id`` and ``state()`` and its tool a ``restore``). Parts are restored in
order, so ``layout`` comes last: it arranges panes the earlier parts created. A part's
``restore`` gets only what its ``capture`` returned, and must skip damaged values rather than
raise. Row lists are document data, saved by the workspace itself.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtCore import QItemSelectionModel
from PySide6.QtWidgets import QAbstractItemView

from . import dock_layout
from .advanced_search import EXTRAS_KEY as ADVANCED_KEY
from .filters.panel import EXTRAS_KEY as FILTERS_KEY

if TYPE_CHECKING:
    from .main_window import MainWindow

log = logging.getLogger(__name__)

PANELS_KEY = "panels"
VIEW_KEY = "view"
LAYOUT_KEY = dock_layout.LAYOUT_KEY
SELECTION_ORIGIN = "session"


@dataclass(frozen=True)
class Part:
    key: str
    #: the window's state as JSON-friendly data; empty (None, {}, []) is not saved
    capture: Callable[[MainWindow], Any]
    #: put saved data back; called only when the document has this key
    restore: Callable[[MainWindow, Any], None]


# -- the view -----------------------------------------------------------------------------------
def capture_view(window: MainWindow) -> dict[str, Any]:
    view, model = window.view, window.workspace.model
    out: dict[str, Any] = {}
    selected = window.workspace.selection.as_array()
    if len(selected):
        out["selected"] = selected.tolist()
    current = view.currentIndex()
    if current.isValid():
        out["current"] = [model.oid_at(current.row()), model.column_name(current.column())]
    top = view.indexAt(view.viewport().rect().topLeft())
    if top.isValid() and top.row() > 0:
        out["top"] = model.oid_at(top.row())
    left = view.horizontalScrollBar().value()
    if left:
        out["left"] = left
    return out


def restore_view(window: MainWindow, saved: Any) -> None:
    if not isinstance(saved, dict):
        return
    view, model, table = window.view, window.workspace.model, window.workspace.table
    current = saved.get("current")
    if isinstance(current, list) and len(current) == 2 and current[1] in table:
        rows = model.rows_of_oids([current[0]]) if isinstance(current[0], int) else []
        if len(rows) and rows[0] >= 0:
            index = model.index(int(rows[0]), model.column_index(current[1]))
            view.selectionModel().setCurrentIndex(index, QItemSelectionModel.SelectionFlag.NoUpdate)
    scroll_view(window, saved)
    # last: moving the current cell may select its row
    selected = saved.get("selected")
    if isinstance(selected, list):
        oids = np.asarray([o for o in selected if isinstance(o, int)], dtype=np.int64)
        oids = oids[np.isin(oids, table.oids)]
        if len(oids):
            window.workspace.selection.set(oids.tolist(), SELECTION_ORIGIN)


def scroll_view(window: MainWindow, saved: Any) -> None:
    """Scroll back to the saved position (again after the first layout pass, which may move
    it while the window settles)."""
    if not isinstance(saved, dict):
        return
    view, model = window.view, window.workspace.model
    top = saved.get("top")
    if isinstance(top, int):
        rows = model.rows_of_oids([top])
        if len(rows) and rows[0] >= 0:
            view.scrollTo(model.index(int(rows[0]), 0), QAbstractItemView.ScrollHint.PositionAtTop)
    left = saved.get("left")
    if isinstance(left, int):
        view.horizontalScrollBar().setValue(left)


# -- the parts, in restore order ----------------------------------------------------------------
def _restore_filters(window: MainWindow, saved: Any) -> None:
    if isinstance(saved, dict):
        window.filters.load_state(saved)


def _restore_search(window: MainWindow, saved: Any) -> None:
    if isinstance(saved, dict):
        window.advanced_search.load_state(saved)


def _restore_panels(window: MainWindow, saved: Any) -> None:
    if isinstance(saved, list):
        window._restore_panels(saved)


PARTS: tuple[Part, ...] = (
    Part(FILTERS_KEY, lambda w: w.filters.state(), _restore_filters),
    Part(ADVANCED_KEY, lambda w: w.advanced_search.state(), _restore_search),
    # names the panels Panel0… as the layout refers to them, so it is captured before it
    Part(PANELS_KEY, lambda w: w.panel_states(), _restore_panels),
    Part(VIEW_KEY, capture_view, restore_view),
    Part(LAYOUT_KEY, dock_layout.capture, lambda w, saved: dock_layout.restore(w, saved)),
)
KEYS = frozenset(part.key for part in PARTS)


def capture(window: MainWindow) -> dict[str, Any]:
    """The window's session, keyed as the document's extras."""
    out: dict[str, Any] = {}
    for part in PARTS:
        try:
            data = part.capture(window)
        except Exception:  # one broken part must not stop the document saving
            log.exception("could not save the %s of the session", part.key)
            continue
        if data:
            out[part.key] = data
    return out


def settle(window: MainWindow, extras: dict[str, Any]) -> None:
    """After the window's first layout pass: apply the saved splitter sizes and scroll
    position once more (that pass shares the space out afresh)."""
    if LAYOUT_KEY in extras:
        dock_layout.restore_docks(window, extras[LAYOUT_KEY])
    if VIEW_KEY in extras:
        scroll_view(window, extras[VIEW_KEY])


def restore(window: MainWindow, extras: dict[str, Any]) -> None:
    """Put a saved session back, part by part; a damaged part is skipped."""
    for part in PARTS:
        if part.key not in extras:
            continue
        try:
            part.restore(window, extras[part.key])
        except Exception:
            log.exception("could not restore the %s of the session", part.key)
