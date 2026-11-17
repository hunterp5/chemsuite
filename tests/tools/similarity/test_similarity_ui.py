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

from tests.gui import gui

pytestmark = gui

SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "c1ccncc1", "C1CCCCC1", ""]


@pytest.fixture(autouse=True)
def fresh_settings(monkeypatch):
    """Similarity Search remembers the last run's settings; each test starts without them."""
    from chemsuite.tools.similarity import search

    monkeypatch.setattr(search, "_last", None)


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
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


def _run(dialog, qtbot):
    assert dialog.problem() == "", dialog.problem()
    dialog.accept()
    qtbot.waitUntil(lambda: not dialog.job.running, timeout=30_000)


def test_similarity_to_a_structure(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.view.setCurrentIndex(chem.workspace.model.index(1, 0))
    dialog = chem.open_tool("similarity")
    assert [e.name for e in dialog.refs.entries()] == ["C1"]  # the current row
    assert dialog.primary.text() == "Compare 9 Rows"
    assert not dialog.refs.section.is_row_visible(dialog.refs.fuse)  # one reference
    dialog.select.setChecked(True)
    dialog.threshold.setValue(0.99)
    _run(dialog, qtbot)
    values = _cells(table, "ECFP4 Tanimoto")
    assert values[1] == "1.000" and values[-1] == ""  # the row itself scores as itself
    assert "Nearest Reference" not in table and "Consensus" not in table
    assert chem.workspace.selection.as_array().tolist() == [int(table.oids[1])]
    assert table.sort_keys[0] == ("ECFP4 Tanimoto", False)
    assert table.history.undo_label == "Similarity search"


def test_references_are_added_typed_and_removed(chem) -> None:
    dialog = chem.open_tool("similarity")
    refs = dialog.refs
    assert not refs.entries() and "Add a reference" in dialog.status.text()
    refs.smiles.setText("not smiles")
    assert not refs.add_button.isEnabled()
    refs.smiles.setText("c1ccncc1")
    refs.add_button.click()
    refs.smiles.setText("CCO")
    refs.smiles.returnPressed.emit()
    assert [e.name for e in refs.entries()] == ["c1ccncc1", "CCO"] and not refs.smiles.text()
    assert refs.section.is_row_visible(refs.fuse)  # several: Best or Mean
    refs.list.setCurrentRow(0)
    refs.remove.click()
    assert [e.name for e in refs.entries()] == ["CCO"] and dialog.primary.isEnabled()
    dialog.reject()


def test_similarity_to_the_selected_rows(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.workspace.selection.set(table.oids[[0, 3]].tolist(), "test")
    dialog = chem.open_tool("similarity")
    assert dialog.refs.source.buttons["selection"].text() == "Selected Rows (2)"
    dialog.refs.source.set_value("selection")
    assert not dialog.refs.section.is_row_visible(dialog.refs.list)
    dialog.sort.setChecked(False)
    _run(dialog, qtbot)
    nearest = _cells(table, "Nearest Reference")
    assert nearest[1] == "C0" and nearest[4] == "C3"
    assert nearest[0] == "C3"  # a selected row is not its own nearest neighbor
    assert _cells(table, "ECFP4 Tanimoto")[0] != "1.000"
    # the fingerprints are cached for the next search
    from chemsuite.analysis.fingerprints import cache_for

    assert cache_for(table).get("SMILES", "morgan2", table.oids) is not None
    chem.workspace.selection.set([], "test")
    assert "No rows are selected" in dialog.refs.problem()


def test_several_methods_give_a_consensus(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.view.setCurrentIndex(chem.workspace.model.index(4, 0))  # toluene
    dialog = chem.open_tool("similarity")
    methods = dialog.methods
    assert not methods.section.is_row_visible(methods.consensus)
    methods.use_pharmacophore.setChecked(True)
    methods.use_substructure.setChecked(True)
    methods.substructure.setCurrentIndex(methods.substructure.findData("scaffold"))
    assert methods.section.is_row_visible(methods.consensus)
    _run(dialog, qtbot)
    for name in ("ECFP4 Tanimoto", "Pharmacophore 2D Tanimoto", "Scaffold Match", "Consensus"):
        assert name in table
    assert _cells(table, "Scaffold Match")[3] == "1.000"  # benzene: the same scaffold
    assert table.sort_keys[0] == ("Consensus", False)
    # the next dialog remembers the methods
    again = chem.open_tool("similarity")
    assert again.methods.use_substructure.isChecked() and again.methods.count() == 3
    again.reject()


def test_3d_keeps_the_aligned_poses_with_their_reference(window, qtbot) -> None:
    from chemsuite.superpose.model import REFERENCE_KIND, reference_id_of
    from tests.tools.superpose.test_superpose_ui import _install

    chem = _install(window)
    table = chem.workspace.table
    chem.view.setCurrentIndex(chem.workspace.model.index(0, 0))
    dialog = chem.open_tool("similarity")
    assert not dialog.output.is_row_visible(dialog.poses)
    dialog.methods.use_3d.setChecked(True)
    assert "no 3D coordinates" in dialog.status.text()  # the SMILES column is flat
    dialog.structure.setCurrentText("Mol 3D")
    dialog.refs.list.setCurrentRow(0)
    dialog.refs.remove.click()
    dialog.refs.add_current_row()
    assert dialog.output.is_row_visible(dialog.poses)
    dialog.poses.setChecked(True)
    _run(dialog, qtbot)
    assert "Shape Tanimoto" in table and "Mol 3D superposed" in table
    assert table.history.undo_label == "Similarity search"
    poses = table.column("Mol 3D superposed")
    assert poses.mol(1) is not None and poses.mol(1).GetNumConformers() == 1
    found = {reference_id_of(poses.mol(i)) for i in range(len(table)) if poses.mol(i)}
    assert len(found) == 1  # one reference, kept once
    assert "Pose ΔE (kcal/mol)" not in table  # molfile conformers carry no energies
    asset = table.assets[found.pop()]
    assert asset.kind == REFERENCE_KIND and asset.name == "A"
    assert asset.meta["source"] == "row" and asset.meta["method"] == "shape"
    table.undo()  # the asset goes with the column
    assert "Mol 3D superposed" not in table and not table.assets


def test_screening_a_file_opens_the_hits(chem, qtbot, tmp_path) -> None:
    path = tmp_path / "library.smi"
    path.write_text("CCO one\nCCCO two\nc1ccccc1 three\nCCCCO four\n", encoding="utf-8")
    opened = []
    chem.open_table = lambda table, message="": opened.append((table, message))
    dialog = chem.open_tool("similarity")
    dialog.where.set_value("file")
    assert dialog.status.text() == "Choose a library file to screen."
    dialog.path.setText(str(tmp_path / "missing.smi"))
    assert dialog.status.text() == "The file was not found."
    dialog.path.setText(str(path))
    dialog.refs.smiles.setText("CCO")
    dialog.refs.add_typed()
    dialog.top.setValue(2)
    assert dialog.primary.text() == "Screen File"
    _run(dialog, qtbot)
    qtbot.waitUntil(lambda: bool(opened), timeout=10_000)
    hits, message = opened[0]
    assert len(hits) == 2 and "ECFP4 Tanimoto" in hits.column_names
    assert message == "2 hits from 4 compounds of library.smi"
    names = [hits.cell(n, int(o)) for o in hits.oids for n in hits.column_names
             if hits.kind(n) != "structure" and n != "ECFP4 Tanimoto"]  # fmt: skip
    assert names[0] == "one"  # itself first


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
    assert _cells(table, "Diversity Rank (2)")[8] == "1"  # MW 80: farthest from MW 0
    assert sorted(r for r in _cells(table, "Diversity Rank") if r) == ["1", "2", "3"]


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
        lambda: "Cluster Size (2)" in table and not job.running,
        timeout=20_000,
    )
    sizes = [int(s) for s in _cells(table, "Cluster Size (2)")]
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
    assert not dialog.status.text().strip() and dialog.primary.isEnabled()
    assert not dialog.warning_mark.isVisibleTo(dialog)
    dialog.features.use_numbers.setChecked(True)
    assert dialog.features.problem() == dialog.status.text()  # no columns ticked
    assert dialog.warning_mark.isVisibleTo(dialog) and not dialog.primary.isEnabled()
    dialog.features.use_fingerprints.setChecked(True)
    # Method settings come and go with the method.
    dialog.method.setCurrentIndex(dialog.method.findData("kmeans"))
    assert dialog._method_section.is_row_visible(dialog.k)
    assert not dialog._method_section.is_row_visible(dialog.cutoff)
    # A rerun keeps the earlier result and adds a new column.
    dialog.method.setCurrentIndex(dialog.method.findData("butina"))
    job = dialog.cluster()
    qtbot.waitUntil(lambda: not job.running, timeout=20_000)
    dialog = chem.open_tool("cluster")
    job = dialog.cluster()
    qtbot.waitUntil(lambda: "Cluster (2)" in table.column_names, timeout=20_000)
    assert "Cluster" in table.column_names


def test_older_dialogs_follow_the_selection_and_close_with_the_document(chem) -> None:
    """Dialogs do not block the table, so their selected-row counts stay current."""
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = chem.workspace.table
    dialog = chem.open_tool("diverse_subset")
    assert not dialog.isModal() and not dialog.seeded.isEnabled()
    chem.workspace.selection.set([int(o) for o in table.oids[:3]], "test")
    assert dialog.seeded.isEnabled() and "(3)" in dialog.seeded.text()
    assert dialog.rows.buttons["selected"].text() == "Selected (3)"
    dialog.seeded.setChecked(True)
    chem.workspace.selection.set([], "test")
    assert not dialog.seeded.isEnabled() and not dialog.seeded.isChecked()
    # A new document closes the dialog: it belongs to the old table.
    chem._install(Workspace(Table()))
    assert not dialog.isVisible()
