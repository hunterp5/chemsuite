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
"""Chemistry ▸ Calculate ▸ Structural Alerts…"""

from __future__ import annotations

import pytest

from tests.gui import gui

pytestmark = gui

CATECHOL = "Oc1ccccc1O"  # PAINS catechol A, Brenk catechol
SMILES = ["CC(=O)Oc1ccccc1C(=O)O", CATECHOL, ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": ["aspirin", "catechol", "blank"], "SMILES": SMILES}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _cells(table, name):
    return [table.cell(name, int(oid)) for oid in table.oids]


def test_matching_writes_a_count_and_a_names_column(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("structural_alerts")
    assert dialog.chosen_keys() == ["PAINS", "BRENK", "NIH"]  # the common set, to start
    dialog.set_checked(["PAINS"])
    assert dialog.verb() == "Match Alert Set for 3 Rows"
    assert not dialog.status.text().strip()
    dialog.accept()
    qtbot.waitUntil(lambda: "PAINS Alerts" in table, timeout=20_000)
    assert _cells(table, "PAINS Alerts") == ["0", "1", ""]
    assert _cells(table, "PAINS Alert Names") == ["", "catechol A", ""]
    assert "Brenk Alerts" not in table  # only the ticked set


def test_count_only_leaves_the_names_out(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("structural_alerts")
    dialog.set_checked(["BRENK"])
    dialog.name_them.setChecked(False)
    dialog.accept()
    qtbot.waitUntil(lambda: "Brenk Alerts" in table, timeout=20_000)
    assert _cells(table, "Brenk Alerts") == ["1", "1", ""]  # aspirin: phenol ester
    assert "Brenk Alert Names" not in table


def test_no_set_ticked_blocks_running(chem) -> None:
    dialog = chem.open_tool("structural_alerts")
    dialog.set_checked(())
    assert dialog.status.text() == "Tick at least one filter set."
    assert not dialog.primary.isEnabled()
    dialog.reject()
