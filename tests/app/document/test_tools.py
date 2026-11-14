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

from tests.gui import gui

pytestmark = gui


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["A", "B", "C"], "SMILES": ["CCO", "c1ccccc1", "not smiles"]}, undoable=False
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_menu_bar_order(window) -> None:
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert titles[:6] == ["File", "Edit", "Tools", "Chemistry", "Data", "Protein"]
    assert titles[-1] == "Help"
    if window.menuBar().isNativeMenuBar():  # macOS: Layout and Log in a Window menu
        assert titles[6:-1] == ["Window"]
    else:  # at the right end of the menu bar
        assert titles[6:-1] == []
        assert window.layout_button.text() == "Layout"
        assert window.log_button.text() == "Log"
    data = window.menuBar().actions()[4].menu()
    groups: list[list[str]] = [[]]
    for a in data.actions():
        if a.isSeparator():
            groups.append([])
        else:
            groups[-1].append(a.text().replace("&", ""))
    # every group in order, the Browser in a section of its own at the bottom
    assert groups == [
        ["Table", "Lists"],
        ["Scoring", "Similarity"],
        ["Dimensionality Reduction…", "Self-Organizing Map…"],
        ["Statistics…", "SAR", "Plot"],
        ["Filters", "Search"],
        ["Browser"],
    ]
    assert data.actions()[-1].shortcut().toString() == "Ctrl+Shift+B"


def test_tools_fill_slots_among_menu_commands(window) -> None:
    def visible(menu) -> list[str]:
        return [a.text().replace("&", "") for a in menu.actions() if a.isVisible()]

    menus = {a.text().replace("&", ""): a.menu() for a in window.menuBar().actions()}
    search = next(a.menu() for a in menus["Data"].actions() if a.text() == "&Search")
    assert visible(search)[4:7] == [
        "Structure Search…",
        "Reaction Search…",
        "Clear Highlighting",
    ]

    def submenu(menu, title):
        return next(a.menu() for a in menu.actions() if a.text().replace("&", "") == title)

    assert [t for t in visible(menus["Tools"]) if t] == ["Utilities", "Sketcher…"]
    utilities = submenu(menus["Tools"], "Utilities")
    assert [t for t in visible(utilities) if t] == ["Calculator…", "Random"]
    assert visible(submenu(utilities, "Random")) == ["Random Numbers…", "Random Rows…"]
    table = [t for t in visible(submenu(menus["Data"], "Table")) if t]
    assert table == ["Add Compound…", "Add Rows…", "Add Column…", "Split Column…",
                     "Join Columns…"]  # fmt: skip
    # two modules fill Similarity, in MENU_ORDER's order
    assert visible(submenu(menus["Data"], "Similarity")) == [
        "Cluster…",
        "Bin Values…",
        "Similarity Search…",
        "Diverse Subset…",
    ]
    assert visible(submenu(menus["Data"], "Scoring")) == ["MPO Score…", "Selectivity…"]
    decompose = submenu(menus["Chemistry"], "Decompose")
    assert visible(decompose) == ["R-Group Decomposition…", "Fragments…", "Central Ring…"]
    enumerate_ = submenu(menus["Chemistry"], "Enumerate")
    assert visible(enumerate_) == [
        "Reaction Products…",
        "R-Groups on a Core…",
        "Recombine Fragments…",
    ]
    assert not {"Structure Viewer", "Browser"} & set(visible(menus["Chemistry"]))
    calculate = submenu(menus["Chemistry"], "Calculate")
    assert visible(calculate) == ["Descriptors…", "pKa…"]
    reaction = submenu(menus["Chemistry"], "Reaction")
    assert visible(reaction) == ["Split Reactions…", "Retrosynthesis…"]
    assert "Preferences…" in visible(menus["Edit"])
    assert "Status Bar" in visible(menus["Edit"])
    assert "Sketcher…" not in visible(menus["Chemistry"])
    edit = visible(menus["Edit"])
    assert "Delete Selection" in edit
    assert not {"Copy Structure As", "Clear Cells", "Edit Structure…", "Lists"} & set(edit)
    assert "View" not in [a.text() for a in menus["Data"].actions()]


def test_registry_builds_the_chemistry_menu(window) -> None:
    action = window.tool_action("descriptors")
    assert action.text() == "Descriptors…"
    assert not action.isEnabled()  # the fixture table has no structures


def test_tool_enabled_once_there_are_structures(chem) -> None:
    assert chem.tool_action("descriptors").isEnabled()


def test_missing_requirement_disables_tool() -> None:
    from chemsuite.app.document.tools import Tool

    tool = Tool("x", "X", ("Chemistry",), open=lambda ctx: None, requires=("no_such_module_q",))
    assert tool.missing() == ["no_such_module_q"]


def test_context_skips_rows_deleted_while_running(chem) -> None:
    from chemsuite.app.document.tools import ToolContext

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
    assert len(dialog.chosen_keys()) == 2
    dialog.accept()
    table = chem.workspace.table
    qtbot.waitUntil(lambda: "MW" in table and "Formula" in table, timeout=20_000)
    assert [table.cell("MW", o) for o in table.oids] == ["46.07", "78.11", ""]
    assert table.kind("MW") == "number"
    assert table.history.undo_label == "Calculate 2 descriptors"
    # Running again keeps the earlier columns and adds new ones (only the selected row here).
    table.set_cells("MW", table.oids[:1], ["0"])
    chem.workspace.selection.set([int(table.oids[0])], "test")
    dialog = chem.open_tool("descriptors")
    assert dialog.chosen_keys() == ["MolWt", "FORMULA"]  # remembered
    dialog.rows.set_value("selected")
    assert dialog.primary.text() == "Calculate 2 Descriptors for 1 Rows"
    dialog.accept()
    qtbot.waitUntil(lambda: "MW (2)" in table and "Formula (2)" in table, timeout=20_000)
    first = int(table.oids[0])
    assert table.cell("MW", first) == "0" and table.cell("MW (2)", first) == "46.07"
    assert table.cell("MW (2)", int(table.oids[1])) == ""


def test_descriptor_dialog_never_replaces_a_column(chem, qtbot) -> None:
    table = chem.workspace.table
    table.add_column("TPSA", ["x"] * 3)
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(["TPSA"])
    dialog.accept()
    qtbot.waitUntil(lambda: "TPSA (2)" in table, timeout=20_000)
    assert table.cell("TPSA", int(table.oids[0])) == "x"


def test_descriptor_dialog_fingerprints(chem, qtbot) -> None:
    from chemsuite.analysis.fingerprints import mol_bits
    from chemsuite.core.mols import mol_from_text

    table = chem.workspace.table
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(["MolWt"])
    dialog.set_fingerprints(["morgan2", "maccs"])
    assert dialog.primary.text() == "Calculate 1 Descriptor and 2 Fingerprints for 3 Rows"
    dialog.set_checked(())
    assert dialog.primary.isEnabled()
    assert dialog.primary.text() == "Calculate 2 Fingerprints for 3 Rows"
    dialog.accept()
    qtbot.waitUntil(lambda: "ECFP4" in table and "MACCS Keys" in table, timeout=20_000)
    assert table.kind("ECFP4") == "text"
    first, last = int(table.oids[0]), int(table.oids[-1])
    smiles = table.cell("SMILES", first)
    bits = mol_bits(mol_from_text(smiles), "morgan2")
    stored = bytes.fromhex(table.cell("ECFP4", first))
    assert len(stored) == 256
    assert (np.unpackbits(np.frombuffer(stored, np.uint8)).astype(bool) == bits).all()
    assert len(table.cell("MACCS Keys", first)) == 2 * 21  # 167 bits in 21 bytes
    assert table.cell("ECFP4", last) == ""  # no structure, no fingerprint
    dialog = chem.open_tool("descriptors")
    assert dialog.chosen_fingerprints() == ["morgan2", "maccs"]  # remembered
    dialog.reject()


def test_descriptor_dialog_filter_and_presets(chem) -> None:
    dialog = chem.open_tool("descriptors")
    dialog.set_checked(())
    dialog.set_fingerprints(())
    ok = dialog.primary
    assert (
        not ok.isEnabled()
        and dialog.status.text() == "Tick at least one descriptor or fingerprint."
    )
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


def test_work_that_holds_the_gil_runs_in_a_worker_process(window) -> None:
    import time
    from types import SimpleNamespace

    from chemsuite.app.document.tools import ToolContext
    from chemsuite.core.jobs import Cancelled, CancelToken

    ctx = ToolContext(window)
    job = SimpleNamespace(cancel_token=CancelToken())
    assert ctx.in_process(job, pow, 2, 10) == 1024
    job.cancel_token.cancel()  # a cancel stops the wait at once, the call left to finish
    started = time.monotonic()
    with pytest.raises(Cancelled):
        ctx.in_process(job, time.sleep, 5)
    assert time.monotonic() - started < 2
