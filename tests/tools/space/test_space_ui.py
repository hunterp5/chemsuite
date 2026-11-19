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

ALCOHOLS = [f"{'C' * k}O" for k in range(1, 11)]
AROMATICS = [f"c1ccccc1{'C' * k}" for k in range(10)]


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    rng = np.random.default_rng(0)
    smiles = ALCOHOLS + AROMATICS + [""]
    group = [0] * 10 + [1] * 10 + [0]
    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i}" for i in range(len(smiles))],
            "SMILES": smiles,
            "pIC50": [f"{5 + g + rng.normal(0, 0.1):.2f}" for g in group],
            "A": [f"{8 * g + rng.normal(0, 0.5):.2f}" for g in group],
            "B": [f"{-8 * g + rng.normal(0, 0.5):.2f}" for g in group],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    for name in ("pIC50", "A", "B"):
        table.convert_column(name, "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _numbers(table, name):
    from chemsuite.plotting.data import numbers

    return numbers(table, name, table.oids)


def _run(qtbot, dialog):
    dialog.accept()
    qtbot.waitUntil(lambda: not dialog.job.running, timeout=120_000)


def test_pca_map_of_number_columns(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_pca")
    assert dialog.windowTitle() == "PCA"
    assert dialog.features.compare is None  # fingerprints and columns are chosen freely
    dialog.features.set_fingerprint(None)
    assert dialog.problem().startswith("Choose a fingerprint")  # nothing to map yet
    dialog.features.set_number_columns(["A", "B"])
    _run(qtbot, dialog)
    x = _numbers(table, "PCA 1")
    first, second = x[:10], x[10:20]
    assert np.sign(first.mean()) != np.sign(second.mean())  # the two groups split on PC1
    assert dialog.plot.spec().x == "PCA 1" and dialog.plot.spec().color == "pIC50"


def test_umap_map_of_fingerprints(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_umap")
    dialog.settings.neighbours.setValue(5)
    _run(qtbot, dialog)
    x, y = _numbers(table, "UMAP 1"), _numbers(table, "UMAP 2")
    assert np.isnan(x[-1]) and np.isfinite(x[:-1]).all()  # the blank structure is left out
    points = np.column_stack([x[:-1], y[:-1]])
    centres = [points[:10].mean(axis=0), points[10:].mean(axis=0)]
    spread = max(np.linalg.norm(points[:10] - centres[0], axis=1).mean(),
                 np.linalg.norm(points[10:] - centres[1], axis=1).mean())  # fmt: skip
    assert np.linalg.norm(centres[0] - centres[1]) > spread


def test_self_organising_map(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("som")
    dialog.features.use_numbers.setChecked(True)
    dialog.features.set_number_columns(["A", "B"])
    dialog.auto.setChecked(False)
    dialog.grid_rows.setValue(3)
    dialog.grid_columns.setValue(4)
    _run(qtbot, dialog)
    panel = dialog.panel
    rows = _numbers(table, "SOM Row")
    assert np.nanmax(rows) <= 3 and np.nanmax(_numbers(table, "SOM Column")) <= 4
    assert panel.grid.rowCount() == 3 and panel.grid.columnCount() == 4
    assert panel.value.currentText() == "pIC50"
    # Pick the cell holding C0 and check that clicking it selects C0 among its compounds.
    r = int(rows[0]) - 1
    c = int(_numbers(table, "SOM Column")[0]) - 1
    chosen = panel.select_cell(r, c)
    assert int(table.oids[0]) in chosen
    ids = {table.cell("ID", o) for o in chosen}
    assert all(int(i[1:]) < 10 or i == "C20" for i in ids)  # only first-group compounds
    panel.value.setCurrentText("Count")
    assert panel.grid.text(r, c) == str(len(chosen))
    assert not panel.grid.hexagonal  # square cells by default
    saved = chem.panel_states()
    assert [s["tool"] for s in saved] == ["som"] and saved[0]["state"]["value"] == "Count"
    chem.restore_panels(saved)
    assert [d.widget().tool_id for d in chem.tool_panels.docks].count("som") == 2


def test_map_options_follow_the_method(chem) -> None:
    dialog = chem.open_tool("space_map_tsne")
    more, settings = dialog._more_section, dialog.settings
    assert more.is_open()  # the optional settings start open
    assert more.is_row_visible(settings.iterations)
    assert not dialog._inputs.is_row_visible(dialog.fit)  # t-SNE cannot place rows
    settings.exact.setChecked(True)
    assert not more.is_row_visible(settings.angle)
    dialog = chem.open_tool("space_map_pca")
    assert dialog._method_section.is_row_visible(dialog.settings.components)
    assert dialog._inputs.is_row_visible(dialog.fit)
    dialog = chem.open_tool("space_map_umap")
    settings = dialog.settings
    dialog.features.set_fingerprint(None)
    dialog.features.set_number_columns(["A", "B"])
    settings.distance.setCurrentIndex(settings.distance.findData("manhattan"))
    settings.seed.setValue(-1)
    options = dialog.options()
    assert options["metric"] == "manhattan" and options["seed"] is None
    assert options["epochs"] is None and options["guide"] is None and options["dims"] == 2


def test_pca_map_on_other_components(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_pca")
    dialog.settings.components.setCurrentText("PC2 and PC3")
    _run(qtbot, dialog)
    assert "PCA 1" in table.column_names  # the columns are named by axis, not component


def test_pca_writes_components_and_shows_loadings(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_pca")
    dialog.features.set_fingerprint(None)
    dialog.features.set_number_columns(["pIC50", "A", "B"])
    dialog.settings.write.setValue(3)
    dialog.settings.components.setCurrentText("PC1 and PC3")
    dialog.settings.loadings.setChecked(True)
    dialog.size.setCurrentText("A")
    _run(qtbot, dialog)
    assert {"PC1", "PC2", "PC3"} <= set(table.column_names) and "PCA 1" not in table
    spec = dialog.plot.spec()
    assert (spec.x, spec.y, spec.size) == ("PC1", "PC3", "A")
    panel = dialog.loadings_panel
    assert panel.table.rowCount() == 3 and panel.columns == ["pIC50", "A", "B"]
    saved = [s for s in chem.panel_states() if s["tool"] == "pca_loadings"]
    assert len(saved) == 1 and len(saved[0]["state"]["loadings"]) == 3
    chem.restore_panels(saved)
    docks = [d.widget().tool_id for d in chem.tool_panels.docks]
    assert docks.count("pca_loadings") == 2


def test_pacmap_map_in_3d(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_pacmap")
    settings = dialog.settings
    assert settings.neighbours.value() == 10  # PaCMAP's own default
    assert dialog._more_section.is_row_visible(settings.mn_ratio)
    assert not dialog._more_section.is_row_visible(settings.local_threshold)
    settings.variant.setCurrentIndex(settings.variant.findData("localmap"))
    assert dialog._more_section.is_row_visible(settings.local_threshold)
    dialog.dims.set_value("3")
    settings.neighbours.setValue(5)
    assert dialog.options()["dims"] == 3 and dialog.options()["variant"] == "localmap"
    _run(qtbot, dialog)
    x = _numbers(table, "PaCMAP 3")
    assert np.isnan(x[-1]) and np.isfinite(x[:-1]).all()
    spec = dialog.plot.spec()
    assert (spec.kind, spec.x, spec.y, spec.z) == ("scatter3d", "PaCMAP 1", "PaCMAP 2", "PaCMAP 3")
    assert chem.open_tool("space_map_umap").settings.neighbours.value() == 15  # UMAP's own
    dialog2 = chem.open_tool("space_map_pca")
    dialog2.dims.set_value("3")
    assert not dialog2._method_section.is_row_visible(dialog2.settings.components)  # PC1-3


def test_maps_never_overwrite_and_combine_features(chem, qtbot) -> None:
    from chemsuite.tools.space.map_dialog import map_columns

    assert map_columns(["ID"], "UMAP", 2) == ["UMAP 1", "UMAP 2"]
    assert map_columns(["UMAP 2"], "UMAP", 2) == ["UMAP 1 (1)", "UMAP 2 (1)"]
    assert map_columns(["UMAP 1", "UMAP 1 (1)"], "UMAP", 3) == [
        "UMAP 1 (2)",
        "UMAP 2 (2)",
        "UMAP 3 (2)",
    ]
    table = chem.workspace.table
    for _run_number in range(2):
        dialog = chem.open_tool("space_map_pca")
        assert dialog.features.fingerprints() and not dialog.features.both()
        dialog.features.set_number_columns(["A", "B"])
        assert dialog.features.both() and not dialog.problem()
        _run(qtbot, dialog)
    assert {"PCA 1", "PCA 2", "PCA 1 (1)", "PCA 2 (1)"} <= set(table.column_names)
    assert np.allclose(_numbers(table, "PCA 1"), _numbers(table, "PCA 1 (1)"), equal_nan=True)
    x = _numbers(table, "PCA 1")
    assert np.isnan(x[-1])  # the blank structure has no fingerprint
    assert np.sign(x[:10].mean()) != np.sign(x[10:20].mean())  # both blocks split the groups
    assert dialog.plot.spec().color == "pIC50"  # opens colored by the activity column


def test_map_defaults_and_bit_distances(chem) -> None:
    umap = chem.open_tool("space_map_umap")
    pacmap = chem.open_tool("space_map_pacmap")
    tsne = chem.open_tool("space_map_tsne")
    assert pacmap.settings.reduce_to.value() == 0 and pacmap.settings.reduce_to.text() == "Off"
    assert umap.settings.distance.currentData() == "jaccard"
    assert pacmap.settings.distance.currentData() == "jaccard"
    assert tsne.settings.init.currentData() == "random"  # t-SNE starts from a random layout
    for dialog in (umap, pacmap):  # fingerprints and a column: bit distances do not fit
        dialog.features.set_number_columns(["A"])
        assert dialog.features.both() and dialog.settings.distance.currentData() == "euclidean"
        jaccard = dialog.settings.distance.findData("jaccard")
        assert not dialog.settings.distance.model().item(jaccard).isEnabled()
        dialog.features.set_number_columns([])  # fingerprints alone again: back to Jaccard
        assert dialog.settings.distance.currentData() == "jaccard"
    pacmap.settings.distance.setCurrentIndex(pacmap.settings.distance.findData("hamming"))
    pacmap.settings.reduce_to.setValue(20)  # PCA dimensions: no bit distances for PaCMAP
    assert pacmap.settings.distance.currentData() == "euclidean"
    pacmap.settings.reduce_to.setValue(0)
    assert pacmap.settings.distance.currentData() == "hamming" and not pacmap.problem()
    umap.features.count.setChecked(True)  # counts are not bits
    assert umap.settings.distance.currentData() == "euclidean"
    umap.features.set_fingerprint(None)
    umap.features.set_number_columns(["A", "B"])
    section = umap.features.section
    assert not section.is_row_visible(umap.features.structure)  # no fingerprints
    assert section.is_row_visible(umap.features.kind)
    assert section.is_row_visible(umap.features.columns)


def test_feature_options_follow_the_fingerprint(chem) -> None:
    dialog = chem.open_tool("space_map_umap")
    features, more = dialog.features, dialog._more_section
    assert features.fingerprint_key() == "morgan2" and features.radius.value() == 2
    assert more.is_row_visible(features.radius) and not more.is_row_visible(features.weight)
    features.set_fingerprint("morgan3")
    assert features.radius.value() == 3  # each kind starts from its own radius
    features.radius.setValue(4)
    features.size.setCurrentIndex(features.size.findData(1024))
    assert features.fingerprint_key() == "morgan3@b1024,r4"
    features.set_fingerprint("maccs")
    assert not more.is_row_visible(features.size) and not more.is_row_visible(features.count)
    assert features.fingerprint_key() == "maccs"
    features.set_fingerprint("morgan2")
    features.set_number_columns(["A"])
    assert more.is_row_visible(features.weight) and more.is_row_visible(features.scaling)


def test_umap_fit_guide_and_counts(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.workspace.lists.add("Reference", table.oids[:12].tolist())
    dialog = chem.open_tool("space_map_umap")
    settings = dialog.settings
    dialog.fit.setCurrentIndex(dialog.fit.findData("list:Reference"))
    settings.neighbours.setValue(4)
    settings.densmap.setChecked(True)
    assert "densMAP" in dialog.problem()  # densMAP cannot place rows
    settings.densmap.setChecked(False)
    settings.guide.setCurrentText("pIC50")
    assert dialog._method_section.is_row_visible(settings.groups)
    dialog.features.count.setChecked(True)
    dialog.open_plot.setChecked(False)
    assert not dialog._output.is_row_visible(dialog.color)
    options = dialog.options()
    assert sorted(set(options["guide"].tolist())) == list(range(5))
    _run(qtbot, dialog)
    x = _numbers(table, "UMAP 1")
    assert np.isfinite(x[:-1]).all() and dialog.plot is None  # columns only
    assert len(chem.tool_panels.docks) == 0


def test_dialogs_remember_their_settings(chem, qtbot) -> None:
    dialog = chem.open_tool("space_map_tsne")
    dialog.settings.perplexity.setValue(12)
    dialog.settings.exact.setChecked(True)
    dialog.features.set_fingerprint("fcfp4")
    dialog.features.radius.setValue(3)
    dialog.features.set_number_columns(["A", "B"])  # the table's own: not remembered
    _run(qtbot, dialog)
    again = chem.open_tool("space_map_tsne")
    assert again.settings.perplexity.value() == 12 and again.settings.exact.isChecked()
    assert again.features.fingerprint_key() == "fcfp4@r3"
    assert again.features.number_columns() == []
    assert chem.open_tool("space_map_umap").features.kind.currentData() == "morgan2"
    again.defaults_button.click()
    assert again.settings.perplexity.value() == 30 and not again.settings.exact.isChecked()
    assert again.features.fingerprint_key() == "morgan2"


def test_som_hexagons_and_options(chem, qtbot) -> None:
    table = chem.workspace.table
    chem.workspace.lists.add("Reference", table.oids[:15].tolist())
    dialog = chem.open_tool("som")
    assert dialog._more_section.is_open()
    dialog.features.use_numbers.setChecked(True)
    dialog.features.set_number_columns(["A", "B"])
    assert not dialog._more_section.is_row_visible(dialog.features.radius)  # no fingerprint
    dialog.auto.setChecked(False)
    dialog.grid_rows.setValue(3)
    dialog.grid_columns.setValue(3)
    dialog.cells.set_value("hexagonal")
    dialog.wrap.setChecked(True)
    dialog.neighborhood.setCurrentIndex(dialog.neighborhood.findData("bubble"))
    dialog.init.setCurrentIndex(dialog.init.findData("random"))
    dialog.start_radius.setValue(0.5)
    assert "start radius" in dialog.problem()
    dialog.start_radius.setValue(2.0)
    dialog.fit.setCurrentIndex(dialog.fit.findData("list:Reference"))
    options = dialog.options()
    assert options["hexagonal"] and options["wrap"] and options["neighborhood"] == "bubble"
    assert len(options["fit"]) == 15 and options["start_radius"] == 2.0
    _run(qtbot, dialog)
    panel = dialog.panel
    assert panel.grid.hexagonal and panel.state()["hexagonal"] and panel.state()["wrap"]
    assert (panel.grid.rowCount(), panel.grid.columnCount()) == (3, 3)  # the trained grid
    assert np.isfinite(_numbers(table, "SOM Row")[:-1]).all()  # unfitted rows are placed too
    # Clicking the middle of a drawn hexagon picks that cell.
    panel.grid.resize(300, 300)
    point = panel.grid.center(1, 1)
    assert panel.grid.cell_at(point) == (1, 1)
    assert panel.grid.shape(1, 1).count() == 6
    saved = chem.panel_states()
    chem.restore_panels(saved)
    restored = [d.widget() for d in chem.tool_panels.docks if d.widget().tool_id == "som"]
    assert all(w.grid.hexagonal for w in restored)
    again = chem.open_tool("som")  # remembered
    assert again.cells.value() == "hexagonal" and again.wrap.isChecked()


def test_som_panels_of_older_documents_are_square(chem) -> None:
    chem.restore_panels([{"tool": "som", "state": {"value": "Count"}}])  # no "hexagonal"
    panel = chem.tool_panels.docks[-1].widget()
    assert panel.tool_id == "som" and not panel.grid.hexagonal
    assert panel.state()["hexagonal"] is False and panel.state()["wrap"] is False


def test_pca_writes_up_to_the_plotted_components_in_3d_with_size(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map_pca")
    dialog.features.set_fingerprint(None)
    dialog.features.set_number_columns(["pIC50", "A", "B"])
    dialog.settings.components.setCurrentText("PC3 and PC4")
    dialog.settings.write.setValue(3)  # fewer than the plotted PC4: PC1 to PC4 are written
    _run(qtbot, dialog)
    assert {"PC1", "PC2", "PC3", "PC4"} <= set(table.column_names)
    assert (dialog.plot.spec().x, dialog.plot.spec().y) == ("PC3", "PC4")
    dialog = chem.open_tool("space_map_umap")
    dialog.dims.set_value("3")
    dialog.settings.neighbours.setValue(5)
    dialog.size.setCurrentText("A")
    _run(qtbot, dialog)
    spec = dialog.plot.spec()
    assert spec.kind == "scatter3d" and spec.z == "UMAP 3"


def test_a_remembered_fingerprint_stays_off_without_structures(chem, qtbot, window) -> None:
    dialog = chem.open_tool("space_map_pca")
    _run(qtbot, dialog)  # remembers a fingerprint
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"A": ["1", "2", "3", "4"], "B": ["4", "1", "3", "2"]}, undoable=False)
    for name in ("A", "B"):
        table.convert_column(name, "number")
    window._install(Workspace(table))
    again = window.open_tool("space_map_pca")
    assert again.features.kind.currentData() is None and again.features.fingerprint_key() is None
