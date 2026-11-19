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

"""Chemistry ▸ Calculate ▸ Sites of Metabolism: the dialog and the columns it writes."""

from __future__ import annotations

import pytest

from tests.gui import gui

pytestmark = gui
#: testosterone (3A4 attacks C6), paracetamol, toluene, and a blank row
SMILES = [
    "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",
    "CC(=O)Nc1ccc(O)cc1",
    "Cc1ccccc1",
    "",
]


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


def test_the_dialog_offers_four_isoforms_and_folds_the_other_five(chem) -> None:
    dialog = chem.open_tool("predict_som")
    assert dialog.primary.text() == "Find Sites for 4 Rows"
    assert not dialog.status.text().strip()
    assert dialog.chosen_isoforms() == ["1A2", "2C9", "2D6", "3A4"]  # the model's own order
    assert set(dialog.isoforms) == {
        "1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4",
    }  # fmt: skip
    assert dialog.chosen() == ["top", "sites"]


def test_a_run_needs_an_isoform_and_a_column(chem) -> None:
    dialog = chem.open_tool("predict_som")
    for box in dialog.isoforms.values():
        box.setChecked(False)
    assert "at least one isoform" in dialog.status.text()
    assert not dialog.primary.isEnabled()
    dialog.isoforms["3A4"].setChecked(True)
    assert not dialog.status.text().strip()
    for box in dialog.outputs.values():
        box.setChecked(False)
    assert "at least one kind of column" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_it_writes_the_ranked_sites_for_each_chosen_isoform(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_som")
    for isoform, box in dialog.isoforms.items():
        box.setChecked(isoform in ("3A4", "2D6"))
    dialog.outputs["domain"].setChecked(True)
    dialog.run()
    qtbot.waitUntil(lambda: "Sites (3A4)" in table.column_names, timeout=120_000)
    assert [name for name in table.column_names if "(3A4)" in name] == [
        "Top Site (3A4)",
        "Score (3A4)",
        "Sites (3A4)",
        "Applicability (3A4)",
    ]
    assert "Top Site (2D6)" in table.column_names
    assert "Top Site (2C9)" not in table.column_names  # an isoform not asked for

    testosterone, paracetamol, toluene, blank = (int(oid) for oid in table.oids)
    # 3A4's 6-beta hydroxylation of testosterone: the allylic carbon, atom 7
    assert table.cell("Top Site (3A4)", testosterone) == "C7"
    assert float(table.cell("Score (3A4)", testosterone)) > 0.5
    # paracetamol's amide nitrogen, the route to NAPQI
    assert table.cell("Top Site (3A4)", paracetamol) == "N3"
    # the ranked text names the top site first, with its score
    sites = table.cell("Sites (3A4)", toluene)
    assert sites.startswith(table.cell("Top Site (3A4)", toluene) + " ")
    assert sites.count(",") == 2
    # every atom of these three is inside the training set
    assert float(table.cell("Applicability (3A4)", paracetamol)) == pytest.approx(1.0)
    # and a blank row stays blank rather than failing
    assert table.cell("Top Site (3A4)", blank) == ""
    assert table.cell("Sites (3A4)", blank) == ""


def test_every_run_writes_new_columns_with_one_suffix(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_som")
    for isoform, box in dialog.isoforms.items():
        box.setChecked(isoform == "3A4")
    dialog.run()
    qtbot.waitUntil(lambda: "Sites (3A4)" in table.column_names, timeout=120_000)

    dialog = chem.open_tool("predict_som")
    for isoform, box in dialog.isoforms.items():
        box.setChecked(isoform == "3A4")
    dialog.run()
    qtbot.waitUntil(lambda: "Sites (3A4) (2)" in table.column_names, timeout=120_000)
    assert "Top Site (3A4) (2)" in table.column_names
    assert "Score (3A4) (2)" in table.column_names


def test_a_run_marks_one_isoform_sites_on_the_structures(chem, qtbot) -> None:
    from chemsuite.app.table.highlights import sources
    from chemsuite.tools.metabolism.dialog import shown_isoform

    table = chem.workspace.table
    dialog = chem.open_tool("predict_som")
    for isoform, box in dialog.isoforms.items():
        box.setChecked(isoform in ("1A2", "3A4"))
    dialog.run()
    qtbot.waitUntil(lambda: "Sites (3A4)" in table.column_names, timeout=120_000)

    # the ranked column of the most prescribed isoform chosen, not the first in model order
    assert shown_isoform(["1A2", "3A4"]) == "3A4"
    assert sources(chem.view) == {"SMILES": "Sites (3A4)"}
    assert "marked on the structures" in chem.statusBar().currentMessage()

    # with the ranked column turned off, the top site alone is marked
    chem.clear_structure_highlights()
    dialog = chem.open_tool("predict_som")
    for isoform, box in dialog.isoforms.items():
        box.setChecked(isoform == "2D6")
    dialog.outputs["sites"].setChecked(False)
    dialog.run()
    qtbot.waitUntil(lambda: "Top Site (2D6)" in table.column_names, timeout=120_000)
    assert sources(chem.view) == {"SMILES": "Top Site (2D6)"}
