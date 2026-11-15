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

"""Bring rows selected elsewhere (a plot) into view in the table."""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QAbstractItemView, QTableView

from chemsuite.core.selection import Selection, SelectionListener


def follow_selection(view: QTableView, selection: Selection, apply: SelectionListener):
    """Subscribe *apply* (the view's own listener) to *selection*; when a plot changes the
    selection, the view then scrolls to the rows it added. Returns the unsubscribe function."""
    before = selection.oids

    def listener(oids: frozenset[int], origin: str) -> None:
        nonlocal before
        added, before = oids - before, oids
        apply(oids, origin)
        if added and origin.startswith("plot:"):
            snap_to(view, rows_for(view.model(), added))

    return selection.subscribe(listener)


def rows_for(model, oids: frozenset[int]) -> np.ndarray:
    """The display rows of *oids* that the table shows, ascending."""
    if not oids:
        return np.empty(0, dtype=np.int64)
    rows = model.rows_of_oids(np.fromiter(oids, dtype=np.int64, count=len(oids)))
    return np.sort(rows[rows >= 0])


def snap_to(view: QTableView, rows: np.ndarray) -> None:
    """Scroll *view* to the first of *rows* (display rows, ascending) unless one is in view.

    Rows the user can already see leave the table where it is, and the table never scrolls
    sideways.
    """
    if not len(rows):
        return
    viewport = view.viewport().rect()
    top = view.rowAt(viewport.top())
    if top >= 0:
        bottom = view.rowAt(viewport.bottom())
        bottom = bottom if bottom >= 0 else view.model().rowCount() - 1
        if np.any((rows >= top) & (rows <= bottom)):
            return
    sideways = view.horizontalScrollBar().value()
    view.scrollTo(
        view.model().index(int(rows[0]), 0), QAbstractItemView.ScrollHint.PositionAtCenter
    )
    view.horizontalScrollBar().setValue(sideways)
