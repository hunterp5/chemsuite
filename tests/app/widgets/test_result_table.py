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
"""The analysis panels' results table: a model, sorted in place, rows keeping their data."""

from __future__ import annotations

import math

from PySide6.QtCore import Qt

from chemsuite.app.widgets.result_table import ResultTable


def test_sorting_keeps_each_row_with_its_data_and_blanks_last(qtbot) -> None:
    table = ResultTable(["Name", "Value"])
    qtbot.addWidget(table)
    table.fill([["b", 2.0], ["a", math.nan], ["c", 10.0], ["d", 1]], ["B", "A", "C", "D"])
    assert table.rowCount() == 4 and table.row_data(0) == "B"  # the given order at first
    table.sortByColumn(1, Qt.SortOrder.DescendingOrder)
    assert [table.row_data(r) for r in range(4)] == ["C", "B", "D", "A"]  # 10 > 2 > 1, blank
    table.sortByColumn(1, Qt.SortOrder.AscendingOrder)
    assert [table.row_data(r) for r in range(4)] == ["D", "B", "C", "A"]  # numbers, not text
    table.sortByColumn(0, Qt.SortOrder.AscendingOrder)
    assert [table.row_data(r) for r in range(4)] == ["A", "B", "C", "D"]
    table.selectRow(2)
    assert table.selected_data() == ["C"]
    assert table.copy_text().splitlines() == ["Name\tValue", "a\t", "b\t2.0", "c\t10.0", "d\t1"]


def test_a_large_fill_is_quick(qtbot) -> None:
    import time

    table = ResultTable(["From", "To", "Pairs"])
    qtbot.addWidget(table)
    rows = [[f"*C{i}", f"*N{i}", i] for i in range(100_000)]
    start = time.perf_counter()
    table.fill(rows, rows)
    assert table.rowCount() == 100_000
    assert time.perf_counter() - start < 0.5  # no per-cell items: cells format when drawn


def test_drawing_cells_sort_and_copy_by_text_and_draw_once_per_size(qtbot) -> None:
    from rdkit import Chem

    from chemsuite.app.widgets.result_table import DRAWING_ROLE, Drawing
    from chemsuite.core.depict import render_png

    def drawing(smiles: str) -> Drawing:
        mol = Chem.MolFromSmiles(smiles)
        return Drawing(
            smiles,
            lambda w, h, scale, dark: render_png(mol, w, h, scale=scale, dark=dark),
            caption=smiles,
        )

    table = ResultTable(["Structure", "Value"], row_height=80, drawing_width=120)
    qtbot.addWidget(table)
    table.resize(400, 300)
    table.fill([[drawing("CCO"), 2], [drawing("CC"), 1], [drawing("c1ccccc1"), 3]], [0, 1, 2])
    table.sortByColumn(0, Qt.SortOrder.AscendingOrder)
    assert [table.row_data(r) for r in range(3)] == [2, 1, 0]  # c1ccccc1, CC, CCO
    assert table.copy_text().splitlines()[1:] == ["c1ccccc1\t3", "CC\t1", "CCO\t2"]
    index = table.result_model.index(0, 0)
    assert isinstance(index.data(DRAWING_ROLE), Drawing)
    assert index.data(Qt.ItemDataRole.ToolTipRole) == "c1ccccc1"
    assert table.columnWidth(0) >= 120
    table.show()
    qtbot.waitExposed(table)
    table.grab()
    drawn = table.drawings.drawn
    assert drawn == 3  # each visible structure once
    table.grab()
    assert table.drawings.drawn == drawn  # repainting reuses the pictures
