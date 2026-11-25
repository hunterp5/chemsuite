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

"""A table view's selection read and replaced without Qt's slow per-row paths."""

from __future__ import annotations

from PySide6.QtCore import QItemSelection, QItemSelectionModel
from PySide6.QtGui import QStandardItemModel

from chemsuite.app.table.view_selection import full_columns, replace_selection
from tests.gui import gui

pytestmark = gui


def _model(rows: int = 10, columns: int = 4) -> QStandardItemModel:
    return QStandardItemModel(rows, columns)


def _ranges(model, *boxes) -> QItemSelection:
    selection = QItemSelection()
    for top, left, bottom, right in boxes:
        selection.select(model.index(top, left), model.index(bottom, right))
    return selection


def test_full_columns_need_every_row_in_one_range_or_several() -> None:
    model = _model()
    assert full_columns(_ranges(model, (0, 1, 9, 1)), 10) == [1]
    assert full_columns(_ranges(model, (0, 0, 4, 2), (5, 2, 9, 3)), 10) == [2]  # split rows
    assert full_columns(_ranges(model, (0, 0, 4, 0), (6, 0, 9, 0)), 10) == []  # row 5 missing
    assert full_columns(_ranges(model, (0, 0, 6, 0), (3, 0, 9, 0)), 10) == [0]  # overlapping
    assert full_columns(_ranges(model, (1, 0, 9, 3)), 10) == []  # row 0 missing
    assert full_columns(QItemSelection(), 10) == []


def test_replace_selection_selects_exactly_the_new_ranges() -> None:
    model = _model()
    chooser = QItemSelectionModel(model)
    chooser.setCurrentIndex(model.index(2, 2), QItemSelectionModel.SelectionFlag.NoUpdate)
    chooser.select(
        _ranges(model, (0, 0, 0, 3), (2, 0, 2, 3)), QItemSelectionModel.SelectionFlag.Select
    )
    replace_selection(chooser, _ranges(model, (1, 0, 1, 3), (3, 1, 5, 2)))
    selected = {(i.row(), i.column()) for i in chooser.selectedIndexes()}
    assert selected == {(1, c) for c in range(4)} | {(r, c) for r in (3, 4, 5) for c in (1, 2)}
    assert chooser.currentIndex() == model.index(2, 2)  # the current cell stays
    replace_selection(chooser, QItemSelection())
    assert not chooser.hasSelection()
