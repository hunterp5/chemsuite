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

"""Chemistry ▸ Standardize Structures."""

from __future__ import annotations

import pytest

from .conftest import gui

pytestmark = gui

SMILES = ["CN1CCCCC1.Cl", "CC(=O)[O-].[Na+]", "CC(N)C(=O)O", "C1CC", ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": [f"C{i}" for i in range(len(SMILES))], "SMILES": SMILES}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_standardize_writes_a_new_column_and_notes(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("standardize")
    assert dialog.primary.text() == "Standardize 5 Rows"
    assert "SMILES (Standardized)" in dialog.status.text()
    for box in dialog.steps.values():
        box.setChecked(False)
    assert "Choose at least one step" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.steps["strip_salts"].setChecked(True)
    dialog.steps["neutralise"].setChecked(True)
    dialog.protonate.setChecked(True)
    assert dialog.ph.isEnabled() and "protonate at pH 7.4" in dialog.status.text()
    dialog.run()
    qtbot.waitUntil(lambda: "SMILES (Standardized)" in table.column_names, timeout=20_000)
    amine, acetate, alanine, bad, blank = (int(o) for o in table.oids)
    out = "SMILES (Standardized)"
    assert table.kind(out) == "structure"
    assert table.cell(out, amine) == "C[NH+]1CCCCC1"
    assert table.cell(out, acetate) == "CC(=O)[O-]"
    assert table.cell(out, alanine) == "CC([NH3+])C(=O)[O-]"
    assert table.cell(out, bad) == "" and table.cell(out, blank) == ""
    assert "undefined stereocentre" in table.cell("Standardization Notes", alanine)
    assert table.cell("Standardization Notes", bad) == "not a valid structure"
    table.undo()
    assert out not in table.column_names


def test_replacing_keeps_structures_that_fail(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("standardize")
    dialog.target.set_value("source")
    dialog.notes.setChecked(False)
    assert dialog.status.text().endswith("→ SMILES")
    dialog.run()
    first = int(table.oids[0])
    qtbot.waitUntil(lambda: table.cell("SMILES", first) == "CN1CCCCC1", timeout=20_000)
    assert table.cell("SMILES", int(table.oids[3])) == "C1CC"  # invalid: kept as it was
    assert "Standardization Notes" not in table.column_names
