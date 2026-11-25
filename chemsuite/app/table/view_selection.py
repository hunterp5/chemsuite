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

"""Working with a table view's selection without Qt's per-row or per-range slow paths.

A selection of scattered rows is one range per row (50,000 for every other row of 100,000).
Some of Qt's selection calls are fine with that; others take seconds or minutes:

* ``select(…, ClearAndSelect)`` over an existing selection compares the old ranges with the
  new ones pairwise: 26 s to swap one 50,000-range selection for another. Clearing first and
  then selecting into the empty model is linear (:func:`replace_selection`).
* ``isColumnSelected`` asks the model about every row of the column, through Python: about a
  second per column on 100,000 rows. :func:`full_columns` reads the ranges instead.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QItemSelection, QItemSelectionModel
from PySide6.QtWidgets import QHeaderView


def replace_selection(model: QItemSelectionModel, selection: QItemSelection) -> None:
    """Make *selection* the whole selection of *model* (the current index stays)."""
    if len(model.selection()):
        model.clearSelection()
    if len(selection):
        model.select(selection, QItemSelectionModel.SelectionFlag.Select)


def full_columns(selection: QItemSelection, rows: int) -> list[int]:
    """The columns that *selection* covers from the first row to the last (*rows* rows),
    whether by one range or several."""
    if not rows or not len(selection):
        return []
    spans = np.array([(r.left(), r.right(), r.top(), r.bottom()) for r in selection])
    left, right, top, bottom = spans.T
    out = []
    for col in range(int(left.min()), int(right.max()) + 1):
        hit = (left <= col) & (right >= col)
        if not hit.any():
            continue
        order = np.argsort(top[hit], kind="stable")
        starts, ends = top[hit][order], bottom[hit][order]
        reach = np.maximum.accumulate(ends)  # the furthest row covered so far
        gaps = starts[1:] > reach[:-1] + 1
        if starts[0] == 0 and reach[-1] >= rows - 1 and not gaps.any():
            out.append(col)
    return out


def mirror_columns(header: QHeaderView, columns: list[int]) -> None:
    """Give a column header's own selection model (``quiet_column_header``) the whole columns
    *columns*, and nothing else: Qt asks it whether a Ctrl+click on the header adds that column
    or takes it away."""
    model = header.model()
    last = model.rowCount() - 1
    selection = QItemSelection()
    for col in columns:
        selection.select(model.index(0, col), model.index(last, col))
    mirror = header.selectionModel()
    if mirror.selection() != selection:  # (a model reset empties it)
        replace_selection(mirror, selection)
