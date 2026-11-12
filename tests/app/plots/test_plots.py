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

from tests.gui import gui

pytestmark = gui


@pytest.fixture(autouse=True)
def gpu(monkeypatch):
    """The offscreen platform has no GPU drawing, so a 3D view would report it and the plot
    fall back to 2D; the tests act as if drawing works (the 3D scene is built all the same)
    and call the real handler when testing the fallback."""
    from chemsuite.app.plots.panel import PlotPanel

    real = PlotPanel._3d_failed
    monkeypatch.setattr(PlotPanel, "_3d_failed", lambda self: None)
    return real


NONE = Qt.KeyboardModifier.NoModifier
SHIFT = Qt.KeyboardModifier.ShiftModifier
CTRL = Qt.KeyboardModifier.ControlModifier


@pytest.fixture
def data(window):
    from chemsuite.app.document.main_window import Workspace
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
    assert panel.summary.text().startswith("4 points · 1 rows without both values")
    assert data.findChild(type(data.tool_panels.docks[0])).windowTitle() == "pIC50 vs MW"


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
    assert "R² = 1.000" in panel.summary.text()


def test_colour_maps_ranges_and_point_sizes(data) -> None:
    from chemsuite.plotting.style import scale_colors

    panel = data.open_tool("scatter")
    form, controls = panel.options.form, panel.controls
    panel.set_spec(color="Series")
    assert not form.isRowVisible(controls.colors) and not form.isRowVisible(controls.size_range)
    panel.set_spec(
        color="pIC50",
        cmap="turbo",
        reverse=True,
        color_range=[6, 8],
        size="pIC50",
        size_range=[None, 7],
        point_sizes=[2, 30],
    )
    assert form.isRowVisible(controls.colors) and form.isRowVisible(controls.color_range)
    assert form.isRowVisible(controls.size_range) and form.isRowVisible(controls.point_sizes)
    spec = panel.spec()
    assert (spec.cmap, spec.reverse, spec.color_range) == ("turbo", True, (6.0, 8.0))
    assert spec.size_range == (None, 7.0) and spec.point_sizes == (2.0, 30.0)
    # A (pIC50 5) lies below the range: with B (6) it takes the first colour of reversed Turbo
    first = scale_colors(24, "turbo", reverse=True)[0]
    layer = next(layer for layer in panel.renderer._layers if len(layer.data) == 2)
    assert layer.opts["brush"].color().getRgb()[:3] == first
    assert sorted(panel.renderer.sizes) == [2.0, 16.0, 30.0, 30.0]
    assert controls.labels[controls.color_range].isVisibleTo(panel.options)

    again = data.open_tool("scatter")
    again.set_state(panel.state())
    assert again.spec() == panel.spec()
    again.set_state({"kind": "scatter"})  # an older document: what it lacks stays
    assert again.spec().cmap == "turbo"


def test_typed_ranges_that_are_not_numbers_go_back_to_auto(data, qtbot) -> None:
    panel = data.open_tool("scatter")
    panel.set_spec(color="pIC50")
    edit = panel.controls.color_range
    with qtbot.waitSignal(panel.controls.changed):
        edit.low.setText("6")
        edit.low.editingFinished.emit()
    assert panel.spec().color_range == (6.0, None)
    edit.high.setText("lots")
    edit.high.editingFinished.emit()
    assert edit.high.text() == "" and panel.spec().color_range == (6.0, None)


def test_density_and_landscape_take_a_colour_map(data) -> None:
    panel = data.open_tool("scatter")
    panel.set_spec(kind="density", cmap="greys")
    assert panel.options.form.isRowVisible(panel.controls.colors)
    assert not panel.options.form.isRowVisible(panel.controls.color_range)


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
    data.tool_panels.docks[0].close()
    qtbot.waitUntil(lambda: not data.tool_panels.docks)
    table.set_cells("pIC50", [int(table.oids[0])], ["1"])  # a closed panel no longer listens


def test_export_png(data, tmp_path) -> None:
    panel = data.open_tool("scatter")
    target = tmp_path / "plot.png"
    panel.export(target)
    assert target.stat().st_size > 1000


def test_plot_menu(data) -> None:
    menus = {a.text(): a.menu() for a in data.menuBar().actions()}
    assert "&Plot" in [a.text() for a in menus["&Data"].actions()]
    assert data.tool_action("scatter").isEnabled()


def test_hover_card_draws_the_structure(window, qtbot) -> None:
    from chemsuite.app.document.main_window import Workspace
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
    assert panel.controls.labels[panel.x].text() == "Group by:"
    # A box drawn around group "x" selects its three rows.
    panel.view_box.boxSelected.emit(-0.5, 0, 0.5, 10, NONE)
    assert ids(data, data.workspace.selection.oids) == ["A", "C", "E"]
    hover = renderer.hover(data.workspace.table, 0.3, 7.0, 0.001, 0.001)
    assert hover is not None and hover.title == "x" and "median 7" in hover.lines[0]


def test_box_plot_without_groups(data) -> None:
    panel = data.open_tool("box")
    panel.set_spec(x=None, y="MW")
    assert panel.renderer.labels == ["MW"]
    assert "4 values in 1 group · 1 blank" in panel.summary.text()


def test_plots_are_saved_and_reopened_with_the_document(data, qtbot, tmp_path) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.app.plots.cards import CardOptions
    from chemsuite.core.table import Table

    first = data.open_tool("scatter")
    first.set_spec(x="MW", y="pIC50", color="Series", fit=True)
    cards = CardOptions(when="both", values=False, columns=("ID",), size="large")
    first.set_card_options(cards)
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
        "aggregate": "mean",
        "errors": False,
        "columns": ["MW", "pIC50"],
        "z": None,
        "drop": False,
        "cmap": "viridis",
        "reverse": False,
        "color_range": [None, None],
        "size_range": [None, None],
        "point_sizes": [4, 18],
        "camera": None,
        "title": None,
        "cards": cards.state(),
    }
    assert histogram.spec().x == "Series" and histogram.spec().kind == "histogram"


def test_damaged_saved_panels_are_skipped(data) -> None:
    data.restore_panels([{"tool": "nope"}, "junk", {"tool": "scatter", "state": {"kind": 3}}])
    assert len(data.panels) == 1  # the scatter entry, with its bad kind ignored
    assert data.panels[0].spec().kind == "scatter"


def test_histogram_curve_and_series_fits(data) -> None:
    panel = data.open_tool("histogram")
    panel.set_spec(x="pIC50", curve="normal")
    assert not panel.curve.isHidden() and panel.fit.isHidden()
    assert "Normal fit: μ = 7, σ = 1.581" in panel.summary.text()
    assert panel.renderer._curve is not None
    panel.set_spec(curve="lognormal")
    assert "Log-normal fit" in panel.summary.text()
    assert panel.state()["curve"] == "lognormal"
    panel.set_spec(x="Series")  # a bar chart: no curve
    assert panel.renderer._curve is None
    scatter = data.open_tool("scatter")
    scatter.set_spec(x="MW", y="pIC50", fit=True)
    assert len(scatter.renderer._fit_lines) == 1 and "R² = 1.000" in scatter.summary.text()
    table = data.workspace.table
    table.add_rows(
        {"MW": ["150", "250", "350"], "pIC50": ["9", "8", "7"], "Series": ["y", "y", "y"]}
    )
    scatter.set_spec(color="Series")
    # x has two plotted points (no fit); y has four, so it gets its own line.
    assert len(scatter.renderer._fit_lines) == 1
    assert "y: R² =" in scatter.summary.text()


def test_options_open_from_the_dock_settings_button(data, qtbot) -> None:
    from PySide6.QtCore import Qt

    data.show()
    panel = data.open_tool("scatter")
    dock = data.tool_panels.docks[-1]
    # Beside the table, not under it.
    assert data.dockWidgetArea(dock) == Qt.DockWidgetArea.RightDockWidgetArea
    bar = dock.title_bar
    assert bar.title.text() == "pIC50 vs MW"
    assert panel.kind.window() is panel.options  # the options are not in the panel itself
    assert not panel.options.isVisible()
    bar.settings.click()
    qtbot.waitUntil(panel.options.isVisible)
    assert not panel.summary.isVisibleTo(panel.options)  # no explanatory text in the options
    panel.options.hide()
    # Histogram options replace the scatter ones.
    panel.set_spec(kind="histogram")
    assert not panel.options.form.isRowVisible(panel.color)
    assert panel.options.form.isRowVisible(panel.bins)
    bar.float_button.click()
    assert dock.isFloating()
    bar.float_button.click()
    assert not dock.isFloating()
    # A second plot joins the first as a tab.
    data.open_tool("histogram")
    assert data.tabifiedDockWidgets(dock) == [data.tool_panels.docks[-1]]
    bar.close_button.click()
    qtbot.waitUntil(lambda: panel not in data.panels)


def test_panels_without_options_have_no_settings_button(data) -> None:
    from PySide6.QtWidgets import QLabel

    dock = data.add_panel(QLabel("plain"), "Plain")
    assert not dock.title_bar.settings.isVisibleTo(dock.title_bar)


def test_rename_a_plot_from_its_header(data, qtbot, tmp_path) -> None:
    from PySide6.QtCore import Qt

    data.show()
    panel = data.open_tool("scatter")
    dock = data.tool_panels.docks[-1]
    title = dock.title_bar.title
    assert panel.toolbar.parent() is not None and dock.isAncestorOf(panel.toolbar)  # in header
    assert title.isReadOnly()
    title.start_editing()
    assert title.editing
    qtbot.keyClicks(title, "Potency map")
    qtbot.keyClick(title, Qt.Key.Key_Return)
    assert not title.editing
    assert panel.title() == "Potency map" and dock.windowTitle() == "Potency map"
    panel.set_spec(x="pIC50", y="MW")  # a named plot keeps its name
    assert title.text() == "Potency map"
    # Escape puts the old title back.
    title.start_editing()
    qtbot.keyClicks(title, "Oops")
    qtbot.keyClick(title, Qt.Key.Key_Escape)
    assert title.text() == "Potency map" and panel.title() == "Potency map"
    # Saved with the document.
    assert panel.state()["title"] == "Potency map"
    target = tmp_path / "named.csd"
    assert data._save_to(target)
    # An empty title goes back to the automatic one.
    title.start_editing()
    title.clear()
    qtbot.keyClick(title, Qt.Key.Key_Return)
    assert panel.title() == "MW vs pIC50" and panel.custom_title is None


def test_hover_card_options(data) -> None:
    from chemsuite.app.plots.cards import CardOptions

    panel = data.open_tool("scatter")
    hover = panel.renderer.hover(data.workspace.table, 200, 6, 1.0, 0.01)
    panel.set_card_options(CardOptions(label=False, columns=("Series",)))
    panel.show_hover(hover)
    assert panel.card.text.text().splitlines() == ["MW: 200", "pIC50: 6", "Series: y"]
    panel.set_card_options(CardOptions(values=False))
    panel.show_hover(hover)
    assert panel.card.text.text() == "B"
    panel.card.hide()
    panel.card_controls.card_when.setCurrentIndex(panel.card_controls.card_when.findData("off"))
    assert panel.card_options.when == "off" and not panel.card_options.hovering
    state = panel.state()["cards"]
    assert state["when"] == "off" and CardOptions.from_state(state) == panel.card_options
    assert CardOptions.from_state({"when": "nonsense", "size": 3}) == CardOptions()


def test_selected_points_keep_their_cards(data, qtbot) -> None:
    data.show()
    qtbot.waitExposed(data)
    panel = data.open_tool("scatter")
    table = data.workspace.table
    first, second = (int(o) for o in table.oids[:2])
    data.workspace.selection.set([first, second], "test")
    assert not panel.pinned.cards  # "On Hover" does not pin
    panel.card_controls.card_when.setCurrentIndex(panel.card_controls.card_when.findData("both"))
    assert set(panel.pinned.cards) == {first, second}
    card = panel.pinned.cards[first][0]
    assert card.text.text().splitlines()[0] == "A"
    qtbot.waitUntil(card.isVisible)
    before = card.pos()
    panel.plot.setXRange(0, 1000, padding=0)  # zooming out moves the card with its point
    qtbot.waitUntil(lambda: card.pos() != before)
    data.workspace.selection.set([first], "test")
    assert set(panel.pinned.cards) == {first}
    data.workspace.selection.set([], "test")
    assert not panel.pinned.cards


def test_pinned_cards_never_overlap() -> None:
    from PySide6.QtCore import QPoint, QRect, QSize

    from chemsuite.app.plots.cards import _free_spot

    bounds, size, taken = QRect(0, 0, 800, 600), QSize(150, 120), []
    for x in (300, 310, 320, 330, 340):  # five points almost on top of each other
        rect = _free_spot(QPoint(x, 300), size, bounds, taken)
        assert bounds.contains(rect) and not any(rect.intersects(t) for t in taken)
        taken.append(rect)


def test_pinned_cards_keep_clear_of_the_other_points() -> None:
    from PySide6.QtCore import QPoint, QRect, QSize

    from chemsuite.app.plots.cards import POINT_CLEARANCE, _free_spot

    bounds, size, taken = QRect(0, 0, 800, 600), QSize(110, 80), []
    points = [QPoint(400 + dx, 300 + dy) for dx in (-30, 0, 30) for dy in (-30, 0, 30)]
    for point in points:  # a tight cluster of nine selected points
        others = [p for p in points if p != point]
        rect = _free_spot(point, size, bounds, taken, others)
        assert bounds.contains(rect) and not any(rect.intersects(t) for t in taken)
        c = POINT_CLEARANCE
        assert not any(rect.intersects(QRect(p.x() - c, p.y() - c, 2 * c, 2 * c)) for p in points)
        taken.append(rect)


def test_card_drawings_fit_the_molecule() -> None:
    from rdkit import Chem

    from chemsuite.app.plots.cards import SIZES
    from chemsuite.core.depict import fitted_size

    largest = SIZES["medium"][1]
    small = fitted_size(Chem.MolFromSmiles("CCO"), largest.width(), largest.height())
    assert small[0] < 0.6 * largest.width() and small[1] < 0.6 * largest.height()
    big = Chem.MolFromSmiles("CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O")
    assert fitted_size(big, largest.width(), largest.height()) == (
        largest.width(),
        largest.height(),
    )


def test_many_pinned_cards_get_smaller_drawings() -> None:
    from PySide6.QtCore import QSize

    from chemsuite.app.plots.cards import SIZES, TINY, pinned_size

    large, plot = SIZES["large"][1], QSize(900, 700)
    assert pinned_size(large, 2, plot) == large
    shrunk = pinned_size(large, 8, plot)
    assert shrunk is not None and shrunk.width() < large.width()
    assert pinned_size(large, 12, QSize(300, 300)) is None  # text only
    assert pinned_size(SIZES["small"][1], 12, QSize(700, 500)) in (SIZES["small"][1], TINY)


def test_bar_chart_sums_up_groups(data) -> None:
    panel = data.open_tool("bar")
    panel.set_spec(x="Series", y="pIC50", aggregate="mean", errors=True)
    renderer = panel.renderer
    assert renderer.labels == ["x", "y", "z"]
    assert renderer.counts.tolist() == [7.0, 6.0, 8.0]  # x: mean of 5, 7, 9
    assert renderer.spread[0] == pytest.approx(2.0) and renderer._errors is not None
    assert panel.title() == "Mean pIC50 by Series"
    assert ids(data, renderer.hit_point(0, 3, 1, 1)) == ["A", "C", "E"]
    hover = renderer.hover(data.workspace.table, 0, 3, 1, 1)
    assert hover.title == "x" and hover.lines[0] == "Mean: 7 ± 2 (SD)"
    panel.set_spec(y=None)  # (none): a row count
    assert renderer.counts.tolist() == [3.0, 1.0, 1.0] and panel.title() == "Rows by Series"
    data.workspace.selection.set([int(data.workspace.table.oids[1])], "test")
    assert renderer._overlay is not None


def test_line_plot_joins_points_in_x_order(data) -> None:
    panel = data.open_tool("line")
    panel.set_spec(x="MW", y="pIC50", color="Series")
    lines = panel.renderer._lines
    assert len(lines) == 1  # only series x has two or more points
    xs, ys = lines[0].getData()
    assert xs.tolist() == [100.0, 300.0] and ys.tolist() == [5.0, 7.0]
    assert ids(data, panel.renderer.hit_box(0, 0, 250, 10)) == ["A", "B"]


def test_density_selects_the_rows_beneath(data) -> None:
    panel = data.open_tool("density")
    panel.set_spec(x="MW", y="pIC50", bins=2)
    renderer = panel.renderer
    assert renderer.grid.counts.sum() == 4  # E has no MW
    assert ids(data, renderer.hit_point(120, 5.2, 1, 1)) == ["A", "B"]
    assert renderer.hover(data.workspace.table, 120, 5.2, 1, 1).title == "2 rows"
    assert ids(data, renderer.hit_box(250, 0, 450, 10)) == ["C", "D"]
    assert panel.title() == "Density: pIC50 vs MW"


def test_radar_draws_selected_rows(data) -> None:
    table = data.workspace.table
    table.add_column("LogP", ["1", "2", "3", "4", "5"])
    panel = data.open_tool("radar")
    panel.set_spec(columns=["MW", "pIC50"])
    assert "at least three" in panel.summary.text()
    panel.set_spec(columns=["MW", "pIC50", "LogP"])
    renderer = panel.renderer
    assert renderer.columns == ["MW", "pIC50", "LogP"]
    first = int(table.oids[3])
    data.workspace.selection.set([first], "test")
    assert set(renderer._voids.tolist()) == {-1, first}
    # D is the highest on every axis: its first corner is at the top of the first spoke
    corner = (renderer._voids == first) & (renderer._vaxis == 0)
    assert renderer._vy[corner][0] == pytest.approx(1.0)
    assert renderer.hit_point(0, 1, 0.01, 0.01).tolist() == [first]
    assert panel.state()["columns"] == ["MW", "pIC50", "LogP"]
    panel.set_spec(kind="scatter")  # leaving the radar brings the axes back
    assert panel.plot.getAxis("left").isVisible()


@pytest.fixture
def cube(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": ["A", "B", "C", "D"],
            "MW": ["100", "200", "300", ""],
            "LogP": ["1", "2", "3", "4"],
            "pIC50": ["5", "6", "7", "8"],
            "Series": ["x", "y", "x", "y"],
        },
        undoable=False,
    )
    window._install(Workspace(table))
    window.resize(900, 600)
    window.show()
    return window


def test_3d_scatter_selects_what_you_see(cube, qtbot) -> None:
    from PySide6.QtCore import QPoint, QRect

    panel = cube.open_tool("scatter3d")
    spec = panel.spec()
    assert (spec.kind, spec.x, spec.y, spec.z) == ("scatter3d", "MW", "LogP", "pIC50")
    assert spec.color is None
    assert panel.title() == "3D: MW, LogP, pIC50"
    assert "1 rows without" in panel.summary.text()
    scene = panel.scene3d
    assert panel.stack.currentWidget() is panel.view3d
    assert scene.pos.min() == -1 and scene.pos.max() == 1  # each axis fills the cube
    qtbot.waitUntil(lambda: panel.view3d.width() > 50)
    xy, front = scene.screen_points()
    assert front.all() and len(xy) == 3
    first = QPoint(round(xy[0, 0]), round(xy[0, 1]))
    assert scene.hit_index(first) == 0
    panel.view3d.clicked.emit(first, NONE)
    assert ids(cube, cube.workspace.selection.oids) == ["A"]
    everything = QRect(0, 0, panel.view3d.width(), panel.view3d.height())
    panel.view3d.boxed.emit(everything, NONE)
    assert ids(cube, cube.workspace.selection.oids) == ["A", "B", "C"]
    around_first = xy[0] + np.array([[-5, -5], [5, -5], [5, 5], [-5, 5]])
    panel.view3d.lassoed.emit(around_first, NONE)
    assert ids(cube, cube.workspace.selection.oids) == ["A"]
    card = panel.card_content(int(cube.workspace.table.oids[1]))
    assert card.title == "B" and card.lines == ["MW: 200", "LogP: 2", "pIC50: 6"]


def test_3d_scatter_is_saved_with_its_camera(cube) -> None:
    panel = cube.open_tool("scatter3d")
    panel.set_spec(color="Series", size="pIC50", drop=True)
    scene = panel.scene3d
    assert scene.colors.shape == (3, 4)
    assert not np.allclose(scene.colors[0], scene.colors[1])  # x and y differ
    assert len(set(np.round(scene.sizes, 3))) == 3  # sized by pIC50
    panel.set_camera({"distance": 6.0, "azimuth": 30.0, "elevation": 10.0})
    state = panel.state()
    assert state["color"] == "Series" and state["camera"]["azimuth"] == 30.0 and state["drop"]
    cube.restore_panels([{"tool": "scatter", "state": state}])
    reopened = cube.panels[-1]
    assert reopened.state() == state
    reopened.set_state({"x": "gone", "camera": "junk"})  # damaged values are skipped
    assert reopened.spec().x == "MW"
    # a document saved while 3D Scatter was a panel of its own
    cube.restore_panels([{"tool": "scatter3d", "state": {"x": "LogP", "y": "MW", "z": "pIC50"}}])
    old = cube.panels[-1]
    assert old.spec().kind == "scatter3d" and old.spec().x == "LogP"


def test_switching_between_2d_and_3d_keeps_the_pane(cube) -> None:
    panel = cube.open_tool("scatter")
    docks = len(cube.panels)
    selection = cube.workspace.selection
    selection.set([int(cube.workspace.table.oids[1])], "test")
    move = panel.mode_buttons["zoom"]
    assert move.text() == "Pan and Zoom"
    panel.set_mode("box")
    panel.set_spec(kind="scatter3d")
    assert len(cube.panels) == docks and panel.stack.currentWidget() is panel.view3d
    assert move.text() == "Turn"  # the same toolbar, its first button now turns the cube
    assert panel.view3d.mode == "box"  # the chosen mode carries over
    assert panel.spec().x == "MW" and panel.spec().y == "LogP" and panel.spec().z == "pIC50"
    assert panel.scene3d.selected_mask(selection.oids).sum() == 1
    assert panel.renderer is None and not panel.plot.listDataItems()  # nothing 2D left behind
    panel.set_spec(kind="landscape")
    assert panel.title() == "Landscape: pIC50 over MW, LogP"
    assert panel.options.form.isRowVisible(panel.bins) and not panel.options.form.isRowVisible(
        panel.color
    )
    panel.set_spec(kind="histogram")
    assert panel.stack.currentWidget() is panel.plot_widget and move.text() == "Pan and Zoom"
    assert panel.spec().z is None and panel.state()["z"] is None


def test_landscape_draws_a_surface(cube) -> None:
    from chemsuite.core.table import Table

    rng = np.random.default_rng(0)
    table = cube.workspace.table
    n = 60
    table.add_rows(
        {
            "ID": [f"N{i}" for i in range(n)],
            "MW": [f"{v:.2f}" for v in rng.uniform(100, 300, n)],
            "LogP": [f"{v:.2f}" for v in rng.uniform(1, 3, n)],
            "pIC50": [f"{v:.2f}" for v in rng.uniform(5, 8, n)],
        },
        undoable=False,
    )
    assert isinstance(table, Table)
    panel = cube.open_tool("landscape")
    panel.set_spec(bins=12)
    assert "surface on a 12×12 grid" in panel.summary.text()
    assert panel.scene3d.has_surface and panel.view3d.surface_triangles > 0
    assert panel.state()["bins"] == 12 and panel.state()["kind"] == "landscape"


def test_radar_fits_its_pane_once_laid_out(data, qtbot) -> None:
    """A radar drawn while its pane is still tiny used to keep that tiny view."""
    table = data.workspace.table
    table.add_column("LogP", ["1", "2", "3", "4", "5"])
    panel = data.open_tool("radar")
    panel.set_spec(columns=["MW", "pIC50", "LogP"])
    panel.plot_widget.resize(60, 40)
    panel.reset_view()
    qtbot.wait(20)
    panel.plot_widget.resize(800, 600)
    qtbot.wait(20)
    (x0, x1), (y0, y1) = panel.view_box.viewRange()
    assert x1 - x0 < 6 and y1 - y0 < 6  # the radar (±1.6 × ±1.25) fills the pane


def test_scatter_item_places_points_like_pyqtgraph(qtbot) -> None:
    """The quick path for points of one style gives each point the atlas entry pyqtgraph
    would; points with their own sizes take pyqtgraph's own path."""
    import numpy as np
    import pyqtgraph as pg

    from chemsuite.app.plots.scatter_item import ScatterItem

    x, y = np.arange(500.0), np.arange(500.0)
    sizes = np.where(np.arange(500) % 2, 5.0, 9.0)
    for kwargs in ({"size": 7.0}, {"size": sizes}):
        ours = ScatterItem(x=x, y=y, pen=None, brush=pg.mkBrush(10, 20, 200), **kwargs)
        stock = pg.ScatterPlotItem(x=x, y=y, pen=None, brush=pg.mkBrush(10, 20, 200), **kwargs)
        rects = ours.data["sourceRect"]
        assert (rects["w"] > 0).all()
        assert len({tuple(r) for r in rects.tolist()}) == len(
            {tuple(r) for r in stock.data["sourceRect"].tolist()}
        )
        ours.setData(x=x[:10], y=y[:10], size=11.0)  # the selection overlay's redraw
        assert (ours.data["sourceRect"]["w"] > 0).all()


def test_3d_kinds_say_why_when_the_graphics_cannot_draw(cube, monkeypatch) -> None:
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", graphics.NO_GPU)
    panel = cube.open_tool("scatter3d")
    item = panel.kind.model().item(panel.kind.findData("scatter3d"))
    assert not item.isEnabled() and item.toolTip() == f"⚠ {graphics.NO_GPU}"
    assert panel.view3d is None  # drawn flat instead of failing


def test_3d_that_fails_when_drawn_falls_back_to_2d(cube, gpu, monkeypatch) -> None:
    """Graphics that pass the start-up check but cannot draw (found only when drawing, as on
    Windows): the plot goes flat and the 3D kinds say why, here and in later plots."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    panel = cube.open_tool("scatter3d")
    gpu(panel)
    assert panel.spec().kind == "scatter" and not panel.is_3d()
    item = panel.kind.model().item(panel.kind.findData("landscape"))
    assert not item.isEnabled() and item.toolTip() == f"⚠ {graphics.NO_GPU}"
    assert graphics.problem() == graphics.NO_GPU
