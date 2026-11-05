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

import pytest

from .conftest import gui

pytestmark = gui


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "SMILES": ["CCO", "c1ccccc1", "not smiles"]}, undoable=False
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_registry_builds_the_chemistry_menu(window) -> None:
    titles = [a.text() for a in window.menuBar().actions()]
    assert titles.index("&Chemistry") == titles.index("&Data") + 1
    action = window.tool_action("descriptors")
    assert action.text() == "Calculate Descriptors…"
    assert not action.isEnabled()  # the fixture table has no structures


def test_tool_enabled_once_there_are_structures(chem) -> None:
    assert chem.tool_action("descriptors").isEnabled()


def test_missing_requirement_disables_tool() -> None:
    from chemsuite.app.tools import Tool

    tool = Tool("x", "X", ("Chemistry",), open=lambda ctx: None, requires=("no_such_module_q",))
    assert tool.missing() == ["no_such_module_q"]


def test_context_skips_rows_deleted_while_running(chem) -> None:
    from chemsuite.app.tools import ToolContext

    ctx = ToolContext(chem)
    table = ctx.table
    oids = table.oids.copy()
    table.remove_rows(oids[1:2])
    names = ctx.apply_results({"N": (oids, [1, 2, 3])}, label="t")
    assert names == {"N": "N"}
    assert [table.cell("N", o) for o in table.oids] == ["1", "3"]
    assert ctx.default_structure_column() == "SMILES"


def test_descriptor_dialog_end_to_end(chem, qtbot) -> None:
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(["MolWt", "FORMULA"])
    assert dialog.summary.text() == "2 descriptor(s) for 3 row(s)"
    dialog.accept()
    table = chem.workspace.table
    qtbot.waitUntil(lambda: "MW" in table and "Formula" in table, timeout=20_000)
    assert [table.cell("MW", o) for o in table.oids] == ["46.07", "78.11", ""]
    assert table.kind("MW") == "number"
    assert table.history.undo_label == "Calculate 2 descriptor(s)"
    # Running again updates the same columns (only the selected row here).
    table.set_cells("MW", table.oids[:1], ["0"])
    chem.workspace.selection.set([int(table.oids[0])], "test")
    dialog = chem.open_tool("descriptors")
    assert dialog.chosen_keys() == ["MolWt", "FORMULA"]  # remembered
    dialog.selected_only.setChecked(True)
    dialog.accept()
    qtbot.waitUntil(lambda: table.cell("MW", int(table.oids[0])) == "46.07", timeout=20_000)
    assert "MW (2)" not in table


def test_descriptor_dialog_new_columns_when_not_overwriting(chem, qtbot) -> None:
    table = chem.workspace.table
    table.add_column("TPSA", ["x"] * 3)
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(["TPSA"])
    dialog.overwrite.setChecked(False)
    dialog.accept()
    qtbot.waitUntil(lambda: "TPSA (2)" in table, timeout=20_000)
    assert table.cell("TPSA", int(table.oids[0])) == "x"


def test_descriptor_dialog_filter_and_presets(chem) -> None:
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(())
    ok = dialog.buttons.button(dialog.buttons.StandardButton.Ok)
    assert not ok.isEnabled()
    dialog.search.setText("inchi")
    visible = [i.text(0) for i in dialog._items() if not i.isHidden()]
    assert visible == ["InChIKey", "InChI"]
    dialog.set_checked(["QED"])
    assert dialog.chosen_keys() == ["QED"] and ok.isEnabled()
    dialog.reject()


def test_tool_menu_shortcut_opens_dialog(chem) -> None:
    from chemsuite.tools.descriptors.dialog import DescriptorDialog

    chem.tool_action("descriptors").trigger()
    dialog = chem.findChild(DescriptorDialog)
    assert dialog is not None and dialog.structure.currentText() == "SMILES"
    dialog.reject()
