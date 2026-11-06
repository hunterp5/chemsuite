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

SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "c1ccncc1", "C1CCCCC1", ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i}" for i in range(len(SMILES))],
            "SMILES": SMILES,
            "MW": [f"{10 * i}" for i in range(len(SMILES))],
            "logP": [f"{(i % 3) * 2.0}" for i in range(len(SMILES))],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    table.convert_column("MW", "number")
    table.convert_column("logP", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _cells(table, name):
    return [table.cell(name, int(o)) for o in table.oids]


def test_similarity_to_a_structure(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.view.setCurrentIndex(chem.workspace.model.index(1, 0))
    dialog = chem.open_tool("similarity")
    assert dialog.query.text() == "CCCO"  # the current row
    dialog.select.setChecked(True)
    dialog.threshold.setValue(0.99)
    job = dialog.calculate()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    values = _cells(table, "Similarity")
    assert values[1] == "1.000" and values[-1] == ""
    assert chem.workspace.selection.as_array().tolist() == [int(table.oids[1])]
    assert table.sort_keys[0] == ("Similarity", False)
    dialog.query.setText("not smiles")
    assert not dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()
    dialog.reject()


def test_similarity_to_the_selection(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.workspace.selection.set(table.oids[[0, 3]].tolist(), "test")
    dialog = chem.open_tool("similarity")
    dialog.to_selection.setChecked(True)
    assert dialog.output.text() == "Max Similarity to Selection"
    dialog.sort.setChecked(False)
    job = dialog.calculate()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    nearest = _cells(table, "Nearest Selected")
    assert nearest[1] == "C0" and nearest[4] == "C3"
    assert nearest[0] == "C3"  # a selected row is not its own nearest neighbour
    # The second run reuses cached fingerprints.
    from chemsuite.analysis.fingerprints import cache_for

    assert cache_for(table).get("SMILES", "morgan2", table.oids) is not None
    dialog.reject()


def test_diverse_subset(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("diverse_subset")
    dialog.count.setValue(3)
    job = dialog.pick()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    ranks = _cells(table, "Diversity Rank")
    assert sorted(r for r in ranks if r) == ["1", "2", "3"]
    assert len(chem.workspace.selection) == 3
    # Numbers, adding to the current selection.
    chem.workspace.selection.set([int(table.oids[0])], "test")
    dialog = chem.open_tool("diverse_subset")
    dialog.features.use_numbers.setChecked(True)
    assert "Tick" in dialog.status.text()
    dialog.features.set_number_columns(["MW"])
    dialog.seeded.setChecked(True)
    dialog.count.setValue(1)
    job = dialog.pick()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    assert _cells(table, "Diversity Rank")[8] == "1"  # MW 80: farthest from MW 0


def test_cluster(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("cluster")
    dialog.cutoff.setValue(0.3)
    job = dialog.cluster()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    clusters = _cells(table, "Cluster")
    assert clusters[-1] == "" and clusters[0] != ""
    assert _cells(table, "Cluster Centre").count("yes") == len({c for c in clusters if c})
    assert table.history.undo_label == "Cluster"
    # K-means on numbers; Butina is not offered for numbers.
    dialog = chem.open_tool("cluster")
    dialog.features.use_numbers.setChecked(True)
    dialog.features.set_number_columns(["MW", "logP"])
    assert "Butina works on fingerprints" in dialog.status.text()
    dialog.method.setCurrentIndex(dialog.method.findData("kmeans"))
    dialog.k.setValue(2)
    job = dialog.cluster()
    qtbot.waitUntil(
        lambda: table.cell("Cluster Size", int(table.oids[0])) != "" and not job.running,
        timeout=20_000,
    )
    sizes = [int(s) for s in _cells(table, "Cluster Size")]
    assert sum(1 for s in sizes if s) == 9 and set(sizes) <= {1, 2, 3, 4, 5, 6, 7, 8}


def test_hierarchical_limit(chem, monkeypatch) -> None:
    from chemsuite.analysis import clustering

    monkeypatch.setattr(clustering, "MAX_HIERARCHICAL", 3)
    dialog = chem.open_tool("cluster")
    dialog.method.setCurrentIndex(dialog.method.findData("hierarchical"))
    assert "up to 3 rows" in dialog.status.text()
    assert not dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()
    dialog.reject()


def test_cluster_dialog_uses_the_dialog_kit(chem, qtbot) -> None:
    """The kit's shared behaviour, through the first dialog built on it."""
    table = chem.workspace.table
    dialog = chem.open_tool("cluster")
    # One label column: the feature rows live in the dialog's Input section.
    assert dialog.features.form is None
    # Help is a tooltip on the field and its label, not text in the dialog.
    assert "tighter clusters" in dialog.cutoff.toolTip()
    # The dialog does not block the window, so the row counts follow the selection live.
    assert not dialog.isModal()
    rows = dialog.rows
    assert rows.buttons["all"].text() == "All (9)"
    rows.set_value("selected")  # nothing selected yet: the footer says so
    assert "No rows are selected" in dialog.status.text() and not dialog.primary.isEnabled()
    chem.workspace.selection.set([int(o) for o in table.oids[:4]], "test")
    assert rows.buttons["selected"].text() == "Selected (4)"
    assert dialog.primary.text() == "Cluster 4 Rows" and dialog.primary.isEnabled()
    # The footer says what will happen, then turns into the problem and blocks the button.
    assert "Butina at 0.60" in dialog.status.text() and dialog.primary.isEnabled()
    assert not dialog.warning_mark.isVisibleTo(dialog)
    dialog.features.use_numbers.setChecked(True)
    assert dialog.features.problem() == dialog.status.text()  # no columns ticked
    assert dialog.warning_mark.isVisibleTo(dialog) and not dialog.primary.isEnabled()
    dialog.features.use_fingerprints.setChecked(True)
    # Method settings come and go with the method.
    dialog.method.setCurrentIndex(dialog.method.findData("kmeans"))
    assert dialog._method_section.is_row_visible(dialog.k)
    assert not dialog._method_section.is_row_visible(dialog.cutoff)
    # "Add a new column" keeps an earlier result.
    dialog.method.setCurrentIndex(dialog.method.findData("butina"))
    job = dialog.cluster()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    dialog = chem.open_tool("cluster")
    dialog.if_exists.set_value("new")
    job = dialog.cluster()
    qtbot.waitUntil(lambda: "Cluster (2)" in table.column_names, timeout=20_000)
    assert "Cluster" in table.column_names


def test_older_dialogs_follow_the_selection_and_close_with_the_document(chem) -> None:
    """Dialogs do not block the table, so their selected-row counts stay current."""
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = chem.workspace.table
    dialog = chem.open_tool("diverse_subset")
    assert not dialog.isModal() and not dialog.seeded.isEnabled()
    selected = dialog.rows.findData("selected")
    chem.workspace.selection.set([int(o) for o in table.oids[:3]], "test")
    assert dialog.seeded.isEnabled() and "(3)" in dialog.seeded.text()
    assert dialog.rows.itemText(selected) == "Selected rows (3)"
    dialog.seeded.setChecked(True)
    chem.workspace.selection.set([], "test")
    assert not dialog.seeded.isEnabled() and not dialog.seeded.isChecked()
    # A new document closes the dialog: it belongs to the old table.
    chem._install(Workspace(Table()))
    assert not dialog.isVisible()
