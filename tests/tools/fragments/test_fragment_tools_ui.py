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

"""Chemistry ▸ Decompose ▸ Fragments."""

from __future__ import annotations

import pytest

from tests.gui import gui

pytestmark = gui

SMILES = [
    "CC(=O)Nc1ccc(O)cc1",  # paracetamol
    "CC(=O)Nc1ccccc1",  # acetanilide: shares the acetyl
    "c1ccccc1",  # nothing to cut
    "",
]
ACTIVITY = [5.0, 7.0, 6.0, None]


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i}" for i in range(len(SMILES))],
            "SMILES": SMILES,
            "pIC50": ACTIVITY,
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_fragments_writes_columns_and_a_count(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("fragments")
    assert dialog.primary.text() == "Fragment 4 Rows"
    dialog.panel.setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: "BRICS 1" in table.column_names, timeout=20_000)
    para, anilide, benzene, blank = (int(o) for o in table.oids)
    assert table.cell("BRICS 1", para) == "[16*]c1ccc(O)cc1"  # largest first
    assert table.cell("BRICS 2", para) == "[1*]C(C)=O"
    assert table.cell("BRICS Fragments", para) == "2"
    assert table.cell("BRICS 1", benzene) == "c1ccccc1"  # nothing to cut: itself
    assert table.cell("BRICS 1", blank) == ""
    assert table.cell("BRICS Fragments", anilide) == "2"


def test_recap_writes_its_own_columns(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("fragments")
    dialog.method.set_value("recap")
    dialog.panel.setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: "RECAP 1" in table.column_names, timeout=20_000)
    para = int(table.oids[0])
    assert table.cell("RECAP 1", para) == "*Nc1ccc(O)cc1"
    assert "BRICS 1" not in table.column_names


def test_plain_attachment_points_drop_the_labels(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("fragments")
    dialog.plain.setChecked(True)
    dialog.panel.setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: "BRICS 1" in table.column_names, timeout=20_000)
    assert table.cell("BRICS 1", int(table.oids[0])) == "[*]c1ccc(O)cc1"


def test_smallest_piece_leaves_out_the_little_fragments(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("fragments")
    dialog.min_size.setValue(5)
    dialog.panel.setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: "BRICS Fragments" in table.column_names, timeout=20_000)
    para = int(table.oids[0])
    assert table.cell("BRICS Fragments", para) == "1"  # the acetyl is too small now
    assert table.cell("BRICS 1", para) == "[16*]c1ccc(O)cc1"


def test_the_panel_lists_distinct_fragments_with_activity(chem, qtbot) -> None:
    dialog = chem.open_tool("fragments")
    dialog.columns.setChecked(False)
    dialog.activity.setCurrentText("pIC50")
    dialog.run()
    qtbot.waitUntil(lambda: chem.tool_panels.find("fragments") is not None, timeout=20_000)
    panel = chem.tool_panels.find("fragments")
    assert panel.title() == "BRICS Fragments"
    acetyl = next(s for s in panel.stats if s.smiles == "[1*]C(C)=O")
    assert acetyl.count == 2  # paracetamol and acetanilide
    assert acetyl.mean == pytest.approx(6.0)
    assert panel.table.horizontalHeaderItem(2).text() == "Mean pIC50"


def test_selecting_a_fragment_selects_its_compounds(chem, qtbot) -> None:
    dialog = chem.open_tool("fragments")
    dialog.columns.setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: chem.tool_panels.find("fragments") is not None, timeout=20_000)
    panel = chem.tool_panels.find("fragments")
    row = next(
        r for r in range(panel.table.rowCount()) if panel.table.row_data(r).smiles == "[1*]C(C)=O"
    )
    panel.table.selectRow(row)
    chosen = panel.select_rows()
    para, anilide = (int(o) for o in chem.workspace.table.oids[:2])
    assert chosen == [para, anilide]
    assert set(chem.workspace.selection.oids) == {para, anilide}


def test_choosing_neither_output_is_refused(chem) -> None:
    dialog = chem.open_tool("fragments")
    dialog.columns.setChecked(False)
    dialog.panel.setChecked(False)
    assert "Choose fragment columns" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_central_ring_writes_a_core_and_r_columns(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("central-ring")
    assert dialog.primary.text() == "Decompose 4 Rows"
    dialog.run()
    qtbot.waitUntil(lambda: "Core" in table.column_names, timeout=20_000)
    para, anilide, benzene, blank = (int(o) for o in table.oids)
    assert table.cell("Core", para) == "c1cc([*:2])ccc1[*:1]"
    assert {table.cell("R1", para), table.cell("R2", para)} == {"*NC(C)=O", "*O"}
    # acetanilide has one substituent on the same ring
    assert table.cell("Core", anilide) == "c1ccc([*:1])cc1"
    assert table.cell("R1", anilide) == "*NC(C)=O"
    assert table.cell("R2", anilide) == ""
    assert table.cell("Core", benzene) == "c1ccccc1"  # a ring with nothing on it
    assert table.cell("Core", blank) == ""


def test_central_ring_needs_no_core_to_be_entered(chem) -> None:
    dialog = chem.open_tool("central-ring")
    assert not hasattr(dialog, "core")
    assert dialog.primary.isEnabled()
