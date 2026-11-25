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
"""Data ▸ SAR ▸ Scaffold Tree: the dialog, the tree panel and reopening it."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QModelIndex, Qt

from chemsuite.app.widgets.result_table import DRAWING_ROLE
from tests.gui import gui

pytestmark = gui

COMPOUNDS = [  # (SMILES, pIC50)
    ("Cc1ccc2[nH]ccc2c1", "6.0"),  # indole
    ("Clc1ccc2[nH]ccc2c1", "7.0"),  # indole
    ("c1ccc2[nH]c(-c3ccccc3)cc2c1", "8.0"),  # 2-phenylindole: under indole
    ("Cc1cc[nH]c1", "5.0"),  # pyrrole itself
    ("c1ccc2ncccc2c1", "4.0"),  # quinoline: under pyridine
    ("CCO", "3.0"),  # no ring
]


@pytest.fixture
def doc(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"c{i}" for i in range(len(COMPOUNDS))],
            "SMILES": [s for s, _ in COMPOUNDS],
            "pIC50": [a for _, a in COMPOUNDS],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _build(doc, qtbot):
    dialog = doc.open_tool("scaffold_tree")
    assert dialog.structure.currentText() == "SMILES"
    dialog.accept()
    qtbot.waitUntil(lambda: not dialog.job.running, timeout=30_000)
    return dialog.panel


def _rows(panel, parent=QModelIndex()):  # noqa: B008
    proxy = panel.proxy
    out = []
    for row in range(proxy.rowCount(parent)):
        index = proxy.index(row, 0, parent)
        count = proxy.index(row, 1, parent).data()
        out.append((index.data(), count, index))
    return out


def test_the_tree_panel(doc, qtbot) -> None:
    panel = _build(doc, qtbot)
    pyrrole, pyridine = _rows(panel)  # largest branch first
    assert (pyrrole[0], pyrrole[1]) == ("c1cc[nH]c1", "4")
    assert (pyridine[0], pyridine[1]) == ("c1ccncc1", "1")
    (indole,) = _rows(panel, pyrrole[2])
    assert indole[1] == "3"
    assert [r[1] for r in _rows(panel, indole[2])] == ["1"]  # 2-phenylindole
    tip = indole[2].data(Qt.ItemDataRole.ToolTipRole)
    assert "3 compounds in this branch, 2 with this scaffold" in tip
    assert "rule 8: the ring with the fewest heteroatoms goes" in tip
    drawing = indole[2].data(DRAWING_ROLE)
    assert drawing is not None and indole[2].data(DRAWING_ROLE) is drawing  # cached once

    # the activity is guessed and summarized over each branch, tinted
    assert panel.show_column.currentText() == "pIC50"
    summary = panel.proxy.index(indole[2].row(), 3, pyrrole[2])
    assert summary.data() == "7.00"
    assert summary.data(Qt.ItemDataRole.BackgroundRole) is not None
    panel.statistic.setCurrentIndex(panel.statistic.findData("max"))
    assert panel.proxy.index(0, 3).data() == "8.00"  # the pyrrole branch's best
    panel.show_column.setCurrentText("None")
    assert panel.view.isColumnHidden(3)

    # clicking a node selects its branch (indexes fetched again: the summary re-sorted)
    (indole,) = _rows(panel, _rows(panel)[0][2])
    oids = doc.workspace.table.oids
    assert sorted(panel.select_branch(indole[2])) == sorted(oids[:3].tolist())
    assert doc.workspace.selection.oids == frozenset(oids[:3].tolist())

    # small branches hide; a parent stays while a child passes
    panel.least.setValue(2)
    assert [r[0] for r in _rows(panel)] == ["c1cc[nH]c1"]
    assert len(_rows(panel, _rows(panel)[0][2])) == 1


def test_the_panel_reopens_with_the_document(doc, qtbot) -> None:
    from chemsuite.tools.sar.common import RerunPanel

    panel = _build(doc, qtbot)
    pyrrole = _rows(panel)[0][2]
    panel.view.expand(pyrrole)
    panel.statistic.setCurrentIndex(panel.statistic.findData("median"))
    panel.higher.setChecked(False)
    saved = doc.panel_states()
    assert [s["tool"] for s in saved] == ["scaffold_tree"]
    assert saved[0]["state"]["expanded"] == ["c1cc[nH]c1"]
    for dock in list(doc.tool_panels.docks):
        dock.close()
    qtbot.waitUntil(lambda: not doc.tool_panels.docks)
    doc.restore_panels(saved)
    (holder,) = doc.panels
    assert isinstance(holder, RerunPanel)
    qtbot.waitUntil(lambda: holder.panel is not None, timeout=30_000)
    assert holder.state() == saved[0]["state"]
    again = holder.panel
    assert again.view.isExpanded(_rows(again)[0][2]) and not again.higher.isChecked()
