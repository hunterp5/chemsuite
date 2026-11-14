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

from tests.gui import gui

pytestmark = gui

SMILES = ["CN1CCCCC1.Cl", "CC(=O)[O-].[Na+]", "CC(N)C(=O)O", "C1CC", ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
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
    assert not dialog.footer.isVisibleTo(dialog)  # no explanatory text when nothing is wrong
    assert set(dialog.target.buttons) == {"source", "new"}
    assert dialog.target.buttons["source"].text() == "Source"
    assert dialog.target.buttons["new"].text() == "New Column"
    for box in dialog.steps.values():
        box.setChecked(False)
    assert "Choose at least one step" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.steps["strip_salts"].setChecked(True)
    dialog.steps["neutralise"].setChecked(True)
    dialog.protonate.setChecked(True)
    assert dialog.ph.isEnabled() and dialog.primary.isEnabled()
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
    dialog.run()
    first = int(table.oids[0])
    qtbot.waitUntil(lambda: table.cell("SMILES", first) == "CN1CCCCC1", timeout=20_000)
    assert table.cell("SMILES", int(table.oids[3])) == "C1CC"  # invalid: kept as it was
    assert "Standardization Notes" in table.column_names  # the notes are always written


def test_running_again_keeps_the_earlier_results(chem, qtbot) -> None:
    table = chem.workspace.table
    for _ in range(2):
        chem.open_tool("standardize").run()
        qtbot.waitUntil(lambda: not chem._runner.active_jobs, timeout=20_000)
    assert "SMILES (Standardized)" in table.column_names
    assert "SMILES (Standardized) (2)" in table.column_names
    assert "Standardization Notes (2)" in table.column_names


def test_stripped_salts_go_to_their_own_column(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("standardize")
    for box in dialog.steps.values():
        box.setChecked(False)
    dialog.steps["strip_salts"].setChecked(True)
    dialog.run()
    qtbot.waitUntil(lambda: "Salts" in table.column_names, timeout=20_000)
    amine, acetate, alanine, _bad, _blank = (int(o) for o in table.oids)
    assert table.cell("Salts", amine) == "Cl"  # the hydrochloride's counter-ion is kept
    assert table.cell("Salts", acetate) == "[Na+]"
    assert table.cell("Salts", alanine) == ""  # nothing was stripped


def test_no_salts_column_when_nothing_was_stripped(chem, qtbot) -> None:
    table = chem.workspace.table
    table.set_cells("SMILES", table.oids[:2], ["CCN", "CCO"])
    dialog = chem.open_tool("standardize")
    for box in dialog.steps.values():
        box.setChecked(False)
    dialog.steps["strip_salts"].setChecked(True)
    dialog.rows.set_value("all")
    dialog.run()
    qtbot.waitUntil(lambda: "Standardization Notes" in table.column_names, timeout=20_000)
    assert "Salts" not in table.column_names


def test_form_salt_adds_a_counterion(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("standardize")
    for box in dialog.steps.values():
        box.setChecked(False)
    dialog.steps["strip_salts"].setChecked(True)
    dialog.steps["neutralise"].setChecked(True)
    assert not dialog.counterion.isEnabled()
    dialog.form_salt.setChecked(True)
    assert dialog.counterion.isEnabled() and dialog.salt_count.isEnabled()
    dialog.counterion.setCurrentText("Mesylate")
    assert dialog.options().salt == "mesylate"
    dialog.run()
    qtbot.waitUntil(lambda: "SMILES (Standardized)" in table.column_names, timeout=20_000)
    amine, _acetate, alanine, _bad, _blank = (int(o) for o in table.oids)
    # the hydrochloride is stripped and the mesylate put on in its place
    assert table.cell("SMILES (Standardized)", amine) == "CS(=O)(=O)[O-].C[NH+]1CCCCC1"
    assert "mesylate salt formed" in table.cell("Standardization Notes", amine)
    assert table.cell("Salts", amine) == "Cl"
    # alanine's acid is the strongest site, so a mesylate has nowhere to go but its amine
    assert "[NH3+]" in table.cell("SMILES (Standardized)", alanine)


def test_form_salt_notes_a_molecule_with_no_site(chem, qtbot) -> None:
    table = chem.workspace.table
    table.set_cells("SMILES", table.oids[:1], ["CCO"])
    dialog = chem.open_tool("standardize")
    dialog.form_salt.setChecked(True)
    dialog.counterion.setCurrentText("Hydrochloride (HCl)")
    dialog.run()
    qtbot.waitUntil(lambda: "Standardization Notes" in table.column_names, timeout=20_000)
    first = int(table.oids[0])
    assert "no site for a hydrochloride salt" in table.cell("Standardization Notes", first)
    assert table.cell("SMILES (Standardized)", first) == "CCO"
