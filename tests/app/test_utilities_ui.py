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

"""Table (utilities): the Calculator, random numbers and rows."""

from __future__ import annotations

import pytest

from .conftest import gui

pytestmark = gui


@pytest.fixture
def data(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["A", "B", "C", "D"],
            "SMILES": ["CCO", "c1ccccc1O", "", "CC(=O)O"],
            "IC50 (nM)": ["10", "1000", "5", ""],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.convert_column("IC50 (nM)", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def cells(table, name):
    return [table.cell(name, int(o)) for o in table.oids]


def finish(qtbot, dialog):
    job = dialog.job
    qtbot.waitUntil(lambda: not job.running, timeout=10_000)
    qtbot.wait(20)


def test_calculator_menu_and_live_check(data):
    action = data.tool_action("calculator")
    assert action.text() == "Calculator…"
    dialog = data.open_tool("calculator")
    assert not dialog.primary.isEnabled()
    dialog.edit.setPlainText("9 - log10([IC50 (nM)]")
    assert "not closed" in dialog.status.text()
    dialog.edit.setPlainText("9 - log10([IC50 (nM)])")
    assert dialog.primary.isEnabled()
    assert dialog.preview.text().startswith("Number: 8  ·  6  ·  8.301029996  ·  (blank)")
    dialog.reject()


def test_calculator_writes_a_column_in_one_undo_step(data, qtbot):
    table = data.workspace.table
    dialog = data.open_tool("calculator")
    dialog.edit.setPlainText("9 - log10([IC50 (nM)])")
    dialog.name.setText("pIC50")
    dialog.accept()
    finish(qtbot, dialog)
    assert table.kind("pIC50") == "number"
    assert cells(table, "pIC50")[:2] == ["8", "6"] and cells(table, "pIC50")[3] == ""
    assert table.history.undo_label == "Calculate pIC50"
    table.undo()
    assert "pIC50" not in table


def test_calculator_chemistry_and_text(data, qtbot):
    table = data.workspace.table
    dialog = data.open_tool("calculator")
    dialog.insert("Formula()")
    assert dialog.edit.toPlainText() == "Formula()"
    dialog.edit.setPlainText('if(substructure(SMILES, "c1ccccc1") > 0, "aryl", Formula(SMILES))')
    dialog.name.setText("Kind")
    dialog.accept()
    finish(qtbot, dialog)
    assert cells(table, "Kind") == ["C2H6O", "aryl", "", "C2H4O2"]


def test_calculator_selected_rows_only(data, qtbot):
    table = data.workspace.table
    data.workspace.selection.set(table.oids[:2], "test")
    dialog = data.open_tool("calculator")
    dialog.rows.set_value("selected")
    dialog.edit.setPlainText("[IC50 (nM)] * 2")
    assert dialog.verb() == "Calculate 2 Rows"
    dialog.accept()
    finish(qtbot, dialog)
    assert cells(table, "Result") == ["20", "2000", "", ""]


def test_calculator_column_menu_lists_columns(data):
    dialog = data.open_tool("calculator")
    dialog.column_menu.aboutToShow.emit()
    titles = [a.text() for a in dialog.column_menu.actions()]
    assert titles[0] == "ID  (text)" and "SMILES  (structure)" in titles
    dialog.column_menu.actions()[2].trigger()
    assert dialog.edit.toPlainText() == "[IC50 (nM)]"
    dialog.reject()


def test_random_numbers_fill_a_column_reproducibly(data):
    table = data.workspace.table
    dialog = data.open_tool("random_numbers")
    dialog.distribution.setCurrentIndex(dialog.distribution.findData("integer"))
    dialog.low.setValue(1)
    dialog.high.setValue(6)
    assert not dialog.method.is_row_visible(dialog.mean)
    dialog.seed.setText("7")
    dialog.accept()
    first = cells(table, "Random")
    assert all(v in {"1", "2", "3", "4", "5", "6"} for v in first)
    dialog = data.open_tool("random_numbers")
    dialog.distribution.setCurrentIndex(dialog.distribution.findData("integer"))
    dialog.low.setValue(1)
    dialog.high.setValue(6)
    dialog.seed.setText("7")
    dialog.accept()
    assert cells(table, "Random (2)") == first


def test_random_numbers_check_their_settings(data):
    dialog = data.open_tool("random_numbers")
    dialog.low.setValue(5)
    dialog.high.setValue(1)
    assert "Low must not be above high" in dialog.status.text()
    dialog.distribution.setCurrentIndex(dialog.distribution.findData("lognormal"))
    assert dialog.method.is_row_visible(dialog.mean) and not dialog.method.is_row_visible(
        dialog.low
    )
    dialog.mean.setValue(0)
    assert "median" in dialog.status.text()
    dialog.reject()


def test_random_rows_select_or_save_a_list(data):
    table = data.workspace.table
    selection = data.workspace.selection
    dialog = data.open_tool("random_rows")
    dialog.unit.set_value("percent")
    dialog.amount.setValue(50)
    assert dialog.verb() == "Pick 2 Rows"
    dialog.accept()
    assert len(selection) == 2 and set(selection.as_array()) <= set(table.oids)
    dialog = data.open_tool("random_rows")
    dialog.amount.setValue(3)
    dialog.action.set_value("list")
    dialog.list_name.setText("Test set")
    dialog.accept()
    assert data.workspace.lists.count("Test set") == 3
    dialog = data.open_tool("random_rows")
    dialog.action.set_value("list")
    dialog.list_name.setText("Test set")
    assert "already exists" in dialog.status.text()
    dialog.reject()
