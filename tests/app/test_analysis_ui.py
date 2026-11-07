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

import math

import numpy as np
import pytest

from .conftest import gui

pytestmark = gui


@pytest.fixture
def data(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    rng = np.random.default_rng(3)
    n = 40
    pic50 = rng.normal(7, 1, n).round(2)
    table = Table()
    table.add_rows(
        {
            "ID": [f"C{i}" for i in range(n)],
            "pIC50": [f"{v}" for v in pic50],
            "hERG pIC50": [f"{v}" for v in (pic50 - rng.uniform(0, 3, n)).round(2)],
            "MW": [f"{v:.1f}" for v in rng.uniform(250, 550, n)],
            "cLogP": [f"{v:.2f}" for v in rng.uniform(0, 6, n)],
        },
        undoable=False,
    )
    table.set_cells("pIC50", table.oids[:1], ["30"])  # an outlier
    table.set_cells("MW", table.oids[1:2], [""])  # a gap
    for name in ("pIC50", "hERG pIC50", "MW", "cLogP"):
        table.convert_column(name, "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _cell(table, name, i):
    return table.cell(name, int(table.oids[i]))


def test_menus(data) -> None:
    for tool in ("statistics", "mpo", "selectivity", "binning"):
        assert data.tool_action(tool).isEnabled()


def test_statistics_panel_tabs(data, qtbot) -> None:
    panel = data.open_tool("statistics")
    qtbot.waitUntil(lambda: panel.summary.grid.rowCount() == 4)
    grid = panel.summary.grid
    assert grid.verticalHeaderItem(2).text() == "MW"
    assert grid.item(2, 0).text() == "39" and grid.item(2, 1).text() == "1"
    # Scope follows the selection.
    table = data.workspace.table
    panel.scope.setCurrentIndex(panel.scope.findData("selected"))
    data.workspace.selection.set(table.oids[:5].tolist(), "test")
    qtbot.waitUntil(lambda: grid.item(0, 0).text() == "5")
    assert panel.rows_label.text() == "5 row(s)"

    panel.scope.setCurrentIndex(panel.scope.findData("all"))
    panel.tabs.setCurrentWidget(panel.correlation)
    qtbot.waitUntil(lambda: panel.correlation.grid.rowCount() == 4)
    hit = panel.correlation.grid.item(0, 1)
    assert -1 <= float(hit.text()) <= 1 and "n = 40" in hit.toolTip()
    panel.correlation.method.setCurrentIndex(1)  # Spearman
    assert panel.correlation.result.r[0, 0] == pytest.approx(1)
    plots_before = len(data.panels)
    panel.correlation._scatter(0, 1)
    assert len(data.panels) == plots_before + 1
    assert data.panels[-1].spec().x == "hERG pIC50"

    panel.tabs.setCurrentWidget(panel.outliers)
    outliers = panel.outliers
    outliers.column.setCurrentText("pIC50")
    assert outliers.label_column.currentText() == "ID"
    assert outliers.grid.item(0, 0).text() == "C0" and outliers.grid.item(0, 1).text() == "30"
    outliers.select_rows()
    assert int(table.oids[0]) in data.workspace.selection.as_array()

    panel.tabs.setCurrentWidget(panel.fit)
    panel.fit.x.setCurrentText("pIC50")
    panel.fit.y.setCurrentText("hERG pIC50")
    assert "R²" in panel.fit.report.toPlainText()

    panel.tabs.setCurrentWidget(panel.tests)
    panel.tests.test.setCurrentIndex(panel.tests.test.findData("shapiro"))
    assert "W =" in panel.tests.report.toPlainText()
    state = panel.state()
    assert state["tab"] == 4 and state["scope"] == "all"


def test_statistics_panel_follows_edits_and_restores(data, qtbot) -> None:
    panel = data.open_tool("statistics")
    table = data.workspace.table
    table.remove_column("cLogP")
    qtbot.waitUntil(lambda: panel.summary.grid.rowCount() == 3)
    assert panel.outliers.column.findText("cLogP") < 0
    states = data.panel_states()
    entry = next(s for s in states if s["tool"] == "statistics")
    data._restore_panels([entry])
    assert sum(1 for p in data.panels if getattr(p, "tool_id", "") == "statistics") == 2


def test_big_correlations_wait_for_the_button(data, qtbot, monkeypatch) -> None:
    from chemsuite.tools.statistics import panel as module

    monkeypatch.setattr(module, "AUTO_CORRELATION_CELLS", 10)
    panel = data.open_tool("statistics")
    panel.tabs.setCurrentWidget(panel.correlation)
    qtbot.waitUntil(lambda: "Calculate" in panel.correlation.note.text())
    assert panel.correlation.grid.rowCount() == 0
    job = panel.correlation.calculate()
    qtbot.waitUntil(lambda: panel.correlation.grid.rowCount() == 4, timeout=10_000)
    assert not job.running


def test_mpo_dialog(data) -> None:
    from chemsuite.analysis.mpo import Desirability, MpoProfile

    dialog = data.open_tool("mpo")
    row = dialog.add_property()
    assert dialog._widget(row, 0).currentText() == "pIC50"
    profile = MpoProfile(
        (
            Desirability("pIC50", "maximize", 6, 8),
            Desirability("MW", "minimize", 350, 500, weight=2),
        ),
        "arithmetic",
    )
    dialog.set_profile(profile)
    assert dialog.profile() == profile
    dialog.per_property.setChecked(True)
    dialog.score(dialog.profile())
    table = data.workspace.table
    assert {"MPO Score", "MPO Score: pIC50", "MPO Score: MW"} <= set(table.column_names)
    assert _cell(table, "MPO Score", 1) == ""  # MW missing
    assert _cell(table, "MPO Score: pIC50", 0) == "1.000"
    # Bad input disables Score.
    dialog._widget(0, 4).setText("abc")
    dialog.set_profile(MpoProfile((Desirability("nope"),)))
    assert "Not number columns" in dialog.status.text()
    dialog.reject()


def test_mpo_profile_files(data, tmp_path, monkeypatch) -> None:
    from PySide6.QtWidgets import QFileDialog

    path = tmp_path / "cns.json"
    dialog = data.open_tool("mpo")
    dialog.add_property()
    dialog.add_property()
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
    dialog.save_dialog()
    saved = dialog.profile()
    dialog.grid.setRowCount(0)
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(path), ""))
    dialog.load_dialog()
    assert dialog.profile() == saved
    dialog.reject()


def test_binning_dialog(data) -> None:
    dialog = data.open_tool("binning")
    dialog.column.setCurrentText("cLogP")
    dialog.method.setCurrentIndex(dialog.method.findData("custom"))
    assert "Enter bin edges" in dialog.status.text()
    dialog.edges.setText("1, 3, 5")
    assert dialog.preview.rowCount() == 4
    assert dialog.preview.item(0, 0).text() == "< 1"
    dialog.names.setText("A, B")
    assert "4 names" in dialog.status.text()
    dialog.names.setText("Very low, Low, Medium, High")
    total = sum(int(dialog.preview.item(r, 2).text()) for r in range(4))
    assert total == 40
    dialog.apply()
    table = data.workspace.table
    value = float(_cell(table, "cLogP", 0))
    expected = (
        "Very low" if value < 1 else "Low" if value < 3 else "Medium" if value < 5 else "High"
    )
    assert _cell(table, "cLogP Class", 0) == expected
    assert table.history.undo_label == "Bin cLogP"
    dialog.reject()


def test_selectivity_dialog(data) -> None:
    dialog = data.open_tool("selectivity")
    assert dialog.primary.currentText() == "pIC50"
    assert dialog.chosen_off_targets() == ["hERG pIC50"]
    dialog.window.setValue(1.0)
    dialog.cutoff.setValue(7.0)
    dialog.apply()
    table = data.workspace.table
    for i in range(1, 40):
        primary = float(_cell(table, "pIC50", i))
        sel = float(_cell(table, "logSel vs hERG pIC50", i))
        assert sel == pytest.approx(primary - float(_cell(table, "hERG pIC50", i)), abs=0.011)
        expected = "yes" if primary >= 7 and sel < 1 else "no"
        assert _cell(table, "Narrow window", i) == expected or math.isclose(sel, 1, abs_tol=0.01)
    assert "Min logSel" not in table  # only one off-target
    dialog.set_off_targets(["hERG pIC50", "MW"])
    assert not dialog.status.text().strip()
    dialog.reject()
