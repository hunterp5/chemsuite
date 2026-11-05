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
from PySide6.QtCore import Qt

from .conftest import gui

pytestmark = gui

NONE = Qt.KeyboardModifier.NoModifier
SHIFT = Qt.KeyboardModifier.ShiftModifier
CTRL = Qt.KeyboardModifier.ControlModifier


@pytest.fixture
def data(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["A", "B", "C", "D", "E"],
            "MW": ["100", "200", "300", "400", ""],
            "pIC50": ["5", "6", "7", "8", "9"],
            "Series": ["x", "y", "x", "z", "x"],
        },
        undoable=False,
    )
    window._install(Workspace(table))
    return window


def ids(window, oids) -> list[str]:
    table = window.workspace.table
    return sorted(table.cell("ID", int(o)) for o in oids)


def test_scatter_opens_in_a_dock_with_defaults(data) -> None:
    panel = data.open_tool("scatter")
    assert panel in data.panels
    assert panel.spec().x == "MW" and panel.spec().y == "pIC50"
    assert panel.status.text().startswith("4 points · 1 rows without both values")
    assert data.findChild(type(data._panels[0])).windowTitle() == "pIC50 vs MW"


def test_box_lasso_and_modifiers_select_rows(data) -> None:
    panel = data.open_tool("scatter")
    selection = data.workspace.selection
    panel.view_box.boxSelected.emit(150, 0, 350, 10, NONE)
    assert ids(data, selection.oids) == ["B", "C"]
    panel.view_box.boxSelected.emit(350, 0, 450, 10, SHIFT)
    assert ids(data, selection.oids) == ["B", "C", "D"]
    panel.view_box.boxSelected.emit(250, 0, 450, 10, CTRL)  # toggles C and D
    assert ids(data, selection.oids) == ["B"]
    square = np.array([[50, 4], [250, 4], [250, 6.5], [50, 6.5]])
    panel.view_box.lassoSelected.emit(square, NONE)
    assert ids(data, selection.oids) == ["A", "B"]


def test_selection_elsewhere_highlights_points(data) -> None:
    panel = data.open_tool("scatter")
    table = data.workspace.table
    data.workspace.selection.set(table.oids[[0, 3]].tolist(), "table")
    overlay = panel.renderer._overlay
    assert sorted(overlay.data["x"].tolist()) == [100, 400]
    assert panel.renderer._layers[0].opacity() < 1


def test_click_selects_nearest_point_and_empty_click_clears(data) -> None:
    panel = data.open_tool("scatter")
    selection = data.workspace.selection
    oids = panel.renderer.hit_point(300, 7, 1.0, 0.01)
    assert ids(data, oids) == ["C"]
    panel._select(oids, NONE)
    assert ids(data, selection.oids) == ["C"]
    panel._select(panel.renderer.hit_point(1000, 100, 1.0, 0.01), NONE)
    assert not selection.oids


def test_plot_follows_filters_and_edits(data, qtbot) -> None:
    from chemsuite.core.table.filters import RangeFilter

    panel = data.open_tool("scatter")
    table = data.workspace.table
    table.set_filter("t", RangeFilter("pIC50", 6, 8))
    qtbot.waitUntil(lambda: len(panel.renderer.oids) == 3)
    table.set_cells("MW", [int(table.oids[4])], ["500"])
    table.set_filter("t", None)
    qtbot.waitUntil(lambda: len(panel.renderer.oids) == 5)


def test_color_by_category_and_number(data) -> None:
    panel = data.open_tool("scatter")
    panel.set_spec(color="Series")
    names = [layer.name() for layer in panel.renderer._layers]
    assert names == ["x (2)", "y (1)", "z (1)"]  # E has no MW, so x counts 2 plotted points
    panel.set_spec(color="pIC50", size="pIC50", fit=True)
    assert len(panel.renderer._layers) == 4  # four distinct colour steps
    assert "R² = 1.000" in panel.status.text()


def test_histogram_numbers_and_text(data) -> None:
    panel = data.open_tool("histogram")
    panel.set_spec(x="pIC50", bins=5)
    assert panel.renderer.counts.tolist() == [1, 1, 1, 1, 1]
    rows = panel.renderer.hit_box(4.9, 0, 6.1, 1)
    assert ids(data, rows) == ["A", "B"]
    panel.set_spec(x="Series")
    assert panel.renderer.labels == ["x", "y", "z"]
    hover = panel.renderer.hover(data.workspace.table, 0, 1, 1, 1)
    assert hover.title == "x" and hover.lines == ("3 rows",)
    assert ids(data, panel.renderer.hit_point(0, 1, 1, 1)) == ["A", "C", "E"]


def test_hover_card_shows_values(data) -> None:
    panel = data.open_tool("scatter")
    hover = panel.renderer.hover(data.workspace.table, 200, 6, 1.0, 0.01)
    assert hover.lines == ("MW: 200", "pIC50: 6")
    panel.show_hover(hover)
    assert panel.card.text.text().splitlines()[0] == "B"
    assert not panel.card.image.isVisible()  # no structure column


def test_removed_column_and_closing(data, qtbot) -> None:
    panel = data.open_tool("scatter")
    table = data.workspace.table
    table.remove_column("MW")
    panel.rebuild()
    assert panel.spec().x == "pIC50"  # the next number column
    data._panels[0].close()
    qtbot.waitUntil(lambda: not data._panels)
    table.set_cells("pIC50", [int(table.oids[0])], ["1"])  # a closed panel no longer listens


def test_export_png(data, tmp_path) -> None:
    panel = data.open_tool("scatter")
    target = tmp_path / "plot.png"
    panel.export(target)
    assert target.stat().st_size > 1000


def test_plot_menu(data) -> None:
    titles = [a.text() for a in data.menuBar().actions()]
    assert "&Plot" in titles
    assert data.tool_action("scatter").isEnabled()


def test_hover_card_draws_the_structure(window, qtbot) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"SMILES": ["CCO", "c1ccccc1"], "MW": ["46", "78"], "LogP": ["0", "2"]}, undoable=False
    )
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    panel = window.open_tool("scatter")
    hover = panel.renderer.hover(table, 78, 2, 1.0, 0.01)
    panel.show_hover(hover)
    assert panel.card.image.isVisible()
    assert not panel.card.image.pixmap().isNull()
    panel.card.hide()


def test_box_plot_groups_and_selection(data) -> None:
    panel = data.open_tool("box")
    panel.set_spec(x="Series", y="pIC50")
    renderer = panel.renderer
    assert renderer.labels == ["x", "y", "z"]
    assert [b.n for b in renderer.stats] == [3, 1, 1]
    assert renderer.stats[0].median == 7.0
    assert panel.title() == "pIC50 by Series"
    assert panel._labels[panel.x].text() == "Group by:"
    # A box drawn around group "x" selects its three rows.
    panel.view_box.boxSelected.emit(-0.5, 0, 0.5, 10, NONE)
    assert ids(data, data.workspace.selection.oids) == ["A", "C", "E"]
    hover = renderer.hover(data.workspace.table, 0.3, 7.0, 0.001, 0.001)
    assert hover is not None and hover.title == "x" and "median 7" in hover.lines[0]


def test_box_plot_without_groups(data) -> None:
    panel = data.open_tool("box")
    panel.set_spec(x=None, y="MW")
    assert panel.renderer.labels == ["MW"]
    assert "4 values in 1 group · 1 blank" in panel.status.text()


def test_plots_are_saved_and_reopened_with_the_document(data, qtbot, tmp_path) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    data.open_tool("scatter").set_spec(x="MW", y="pIC50", color="Series", fit=True)
    data.open_tool("histogram").set_spec(x="Series")
    target = tmp_path / "plots.csd"
    assert data._save_to(target)
    data._install(Workspace(Table()))
    assert not data.panels
    data.open_path(target)
    qtbot.waitUntil(lambda: len(data.panels) == 2, timeout=10_000)
    scatter, histogram = data.panels
    assert scatter.state() == {
        "kind": "scatter",
        "x": "MW",
        "y": "pIC50",
        "color": "Series",
        "size": None,
        "bins": None,
        "fit": True,
        "curve": None,
        "guide": None,
    }
    assert histogram.spec().x == "Series" and histogram.spec().kind == "histogram"


def test_damaged_saved_panels_are_skipped(data) -> None:
    data._restore_panels([{"tool": "nope"}, "junk", {"tool": "scatter", "state": {"kind": 3}}])
    assert len(data.panels) == 1  # the scatter entry, with its bad kind ignored
    assert data.panels[0].spec().kind == "scatter"


def test_histogram_curve_and_series_fits(data) -> None:
    panel = data.open_tool("histogram")
    panel.set_spec(x="pIC50", curve="normal")
    assert not panel.curve.isHidden() and panel.fit.isHidden()
    assert "Normal fit: μ = 7, σ = 1.581" in panel.status.text()
    assert panel.renderer._curve is not None
    panel.set_spec(curve="lognormal")
    assert "Log-normal fit" in panel.status.text()
    assert panel.state()["curve"] == "lognormal"
    panel.set_spec(x="Series")  # a bar chart: no curve
    assert panel.renderer._curve is None
    scatter = data.open_tool("scatter")
    scatter.set_spec(x="MW", y="pIC50", fit=True)
    assert len(scatter.renderer._fit_lines) == 1 and "R² = 1.000" in scatter.status.text()
    table = data.workspace.table
    table.add_rows(
        {"MW": ["150", "250", "350"], "pIC50": ["9", "8", "7"], "Series": ["y", "y", "y"]}
    )
    scatter.set_spec(color="Series")
    # x has two plotted points (no fit); y has four, so it gets its own line.
    assert len(scatter.renderer._fit_lines) == 1
    assert "y: R² =" in scatter.status.text()
