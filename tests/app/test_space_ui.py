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

from .conftest import gui

pytestmark = gui

ALCOHOLS = [f"{'C' * k}O" for k in range(1, 11)]
AROMATICS = [f"c1ccccc1{'C' * k}" for k in range(10)]


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
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
    dialog = chem.open_tool("space_map")
    dialog.method.setCurrentIndex(dialog.method.findData("pca"))
    assert dialog.windowTitle() == "Dimensionality Reduction"
    dialog.features.use_numbers.setChecked(True)
    assert dialog.problem()  # no columns ticked yet
    dialog.features.set_number_columns(["A", "B"])
    _run(qtbot, dialog)
    x = _numbers(table, "PCA 1")
    first, second = x[:10], x[10:20]
    assert np.sign(first.mean()) != np.sign(second.mean())  # the two groups split on PC1
    assert dialog.plot.spec().x == "PCA 1" and dialog.plot.spec().color == "pIC50"


def test_umap_map_of_fingerprints(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map")
    dialog.method.setCurrentIndex(dialog.method.findData("umap"))
    dialog.neighbours.setValue(5)
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
    assert panel.grid.item(r, c).text() == str(len(chosen))
    saved = chem.panel_states()
    assert [s["tool"] for s in saved] == ["som"] and saved[0]["state"]["value"] == "Count"
    chem._restore_panels(saved)
    assert [d.widget().tool_id for d in chem._panels].count("som") == 2


def test_map_options_follow_the_method(chem) -> None:
    dialog = chem.open_tool("space_map")
    more = dialog._more_section
    assert not more.is_open()  # optional settings start folded
    more.set_open(True)
    dialog.method.setCurrentIndex(dialog.method.findData("tsne"))
    assert more.is_row_visible(dialog.iterations) and not more.is_row_visible(dialog.spread)
    dialog.method.setCurrentIndex(dialog.method.findData("pca"))
    assert not more.isVisibleTo(dialog)  # PCA has no optional settings
    assert dialog._method_section.is_row_visible(dialog.components)
    dialog.method.setCurrentIndex(dialog.method.findData("umap"))
    dialog.features.use_numbers.setChecked(True)
    dialog.features.set_number_columns(["A", "B"])
    dialog.distance.setCurrentIndex(dialog.distance.findData("dice"))
    assert "needs fingerprints" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.distance.setCurrentIndex(dialog.distance.findData("manhattan"))
    dialog.seed.setValue(-1)
    options = dialog.options()
    assert options["metric"] == "manhattan" and options["seed"] is None
    assert options["epochs"] is None and options["learning_rate"] is None


def test_pca_map_on_other_components(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map")
    dialog.method.setCurrentIndex(dialog.method.findData("pca"))
    dialog.components.setCurrentText("PC2 and PC3")
    _run(qtbot, dialog)
    assert "PCA 1" in table.column_names  # the columns are named by axis, not component


def test_pacmap_map_in_3d(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("space_map")
    dialog.method.setCurrentIndex(dialog.method.findData("pacmap"))
    assert dialog.neighbours.value() == 10  # PaCMAP's own default
    assert dialog._method_section.is_row_visible(dialog.neighbours)
    assert dialog._more_section.is_row_visible(dialog.mn_ratio)
    dialog.dims.set_value("3")
    dialog.neighbours.setValue(5)
    assert dialog.options()["dims"] == 3
    _run(qtbot, dialog)
    x = _numbers(table, "PaCMAP 3")
    assert np.isnan(x[-1]) and np.isfinite(x[:-1]).all()
    spec = dialog.plot.spec()
    assert (spec.kind, spec.x, spec.y, spec.z) == ("scatter3d", "PaCMAP 1", "PaCMAP 2", "PaCMAP 3")
    dialog2 = chem.open_tool("space_map")
    dialog2.method.setCurrentIndex(dialog2.method.findData("pacmap"))
    dialog2.method.setCurrentIndex(dialog2.method.findData("umap"))
    assert dialog2.neighbours.value() == 15  # back to UMAP's default
    dialog2.dims.set_value("3")
    dialog2.method.setCurrentIndex(dialog2.method.findData("pca"))
    assert not dialog2._method_section.is_row_visible(dialog2.components)  # 3D: PC1-3


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
        dialog = chem.open_tool("space_map")
        dialog.method.setCurrentIndex(dialog.method.findData("pca"))
        dialog.features.compare.set_value("both")
        assert dialog.features.both() and dialog.problem() == "Tick at least one number column."
        dialog.features.set_number_columns(["A", "B"])
        assert not dialog.problem()
        _run(qtbot, dialog)
    assert {"PCA 1", "PCA 2", "PCA 1 (1)", "PCA 2 (1)"} <= set(table.column_names)
    assert np.allclose(_numbers(table, "PCA 1"), _numbers(table, "PCA 1 (1)"), equal_nan=True)
    x = _numbers(table, "PCA 1")
    assert np.isnan(x[-1])  # the blank structure has no fingerprint
    assert np.sign(x[:10].mean()) != np.sign(x[10:20].mean())  # both blocks split the groups
    assert dialog.plot.spec().color == "pIC50"  # opens coloured by the activity column
