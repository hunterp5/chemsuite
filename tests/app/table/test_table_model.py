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

from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import QPersistentModelIndex, Qt

from chemsuite.app.table.table_model import KIND_ROLE, OID_ROLE, STATUS_ROLE, TableModel
from chemsuite.core.mols import parse_structures
from chemsuite.core.table import RangeFilter, Table
from chemsuite.core.table.columns import StructureColumn

DISPLAY = Qt.ItemDataRole.DisplayRole


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "Structure": StructureColumn.from_texts(["CCO", "C1CC", "", "c1ccccc1"]),
            "ID": ["CPD-10", "CPD-2", "CPD-1", "CPD-3"],
            "pIC50": ["7.50", "6.2", "", "8"],
        },
        undoable=False,
    )
    return t


@pytest.fixture
def model(qapp, qtmodeltester, table: Table) -> TableModel:
    m = TableModel(table)
    qtmodeltester.check(m, force_py=True)  # Qt's model tester stays attached for every change
    return m


def column(model: TableModel, name: str) -> list[str]:
    col = model.column_index(name)
    return [model.data(model.index(r, col), DISPLAY) for r in range(model.rowCount())]


def test_shape_and_roles(model: TableModel) -> None:
    assert (model.rowCount(), model.columnCount()) == (4, 3)
    assert model.headerData(1, Qt.Orientation.Horizontal, DISPLAY) == "ID"
    assert model.headerData(0, Qt.Orientation.Vertical, DISPLAY) == "1"
    index = model.index(0, 2)
    assert model.data(index, DISPLAY) == "7.50"
    assert model.data(index, Qt.ItemDataRole.TextAlignmentRole) & Qt.AlignmentFlag.AlignRight
    assert model.data(index, OID_ROLE) == 1
    assert model.data(index, KIND_ROLE) == "number"


def test_structure_status(model: TableModel, table: Table) -> None:
    statuses = [model.data(model.index(r, 0), STATUS_ROLE) for r in range(4)]
    assert statuses == ["pending", "pending", "blank", "pending"]
    table.fill_structures(
        "Structure", table.oids, parse_structures(["CCO", "C1CC", "", "c1ccccc1"])
    )
    statuses = [model.data(model.index(r, 0), STATUS_ROLE) for r in range(4)]
    assert statuses == ["ok", "invalid", "blank", "ok"]
    assert model.data(model.index(1, 0), Qt.ItemDataRole.ToolTipRole).startswith("Invalid")


def test_editing_goes_through_the_table_and_is_undoable(model: TableModel, table: Table) -> None:
    changed: list[tuple[int, int]] = []
    model.dataChanged.connect(lambda tl, br: changed.append((tl.row(), br.row())))
    assert model.setData(model.index(1, 1), "NEW")
    assert table.cell("ID", 2) == "NEW"
    assert changed == [(1, 1)]
    assert not model.setData(model.index(1, 1), "NEW")  # unchanged text is not an edit
    table.undo()
    assert column(model, "ID")[1] == "CPD-2"


def test_sorting_keeps_persistent_indexes_on_the_same_row(model: TableModel, table: Table) -> None:
    tracked = QPersistentModelIndex(model.index(0, 1))  # CPD-10
    model.sort(1, Qt.SortOrder.AscendingOrder)
    assert column(model, "ID") == ["CPD-1", "CPD-2", "CPD-3", "CPD-10"]
    assert tracked.row() == 3 and tracked.data(DISPLAY) == "CPD-10"
    model.sort(0)  # structure columns cannot be sorted: nothing happens
    assert column(model, "ID") == ["CPD-1", "CPD-2", "CPD-3", "CPD-10"]


def test_filters_rows_and_columns_reset_the_model(model: TableModel, table: Table) -> None:
    resets: list[int] = []
    model.modelReset.connect(lambda: resets.append(model.rowCount()))
    table.set_filter("p", RangeFilter("pIC50", low=7))
    assert column(model, "ID") == ["CPD-10", "CPD-3"]
    table.add_rows({"ID": ["CPD-99"], "pIC50": ["9"]})
    table.remove_rows([1])
    table.add_column("Notes")
    table.remove_column("Notes")
    assert resets == [2, 3, 2, 2, 2]
    assert column(model, "ID") == ["CPD-3", "CPD-99"]


def test_type_change_refreshes_without_reset(model: TableModel, table: Table) -> None:
    resets: list[int] = []
    model.modelReset.connect(lambda: resets.append(1))
    table.set_cells("pIC50", [1], [">10"])  # number column becomes text
    assert resets == []
    assert model.data(model.index(0, 2), KIND_ROLE) == "text"


def test_rows_of_oids(model: TableModel, table: Table) -> None:
    table.set_sort([("ID", False)])
    rows = model.rows_of_oids(np.array([1, 3, 99]))
    assert rows.tolist() == [0, 3, -1]
    assert model.rows_of_oids([]).tolist() == []


def test_big_model_stays_fast(qapp, table: Table) -> None:
    import time

    big = Table()
    big.add_rows({"x": np.random.default_rng(0).random(200_000).tolist()}, undoable=False)
    m = TableModel(big)
    start = time.perf_counter()
    big.set_sort([("x", True)])  # layout change with 200k rows
    big.set_filter("f", RangeFilter("x", low=0.5))  # reset
    assert time.perf_counter() - start < 0.5
    assert m.rowCount() == len(big.view())
