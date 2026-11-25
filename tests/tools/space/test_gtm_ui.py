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
"""The Generative Topographic Map tool and its GTM Landscape plot."""

from __future__ import annotations

import numpy as np

from tests.gui import gui
from tests.tools.space.test_space_ui import _numbers, _run, chem  # noqa: F401 (a fixture)

pytestmark = gui

COLUMNS = ("GTM X", "GTM Y", "GTM Log-likelihood")


def _train(chem, qtbot, side: int = 5):  # noqa: F811
    dialog = chem.open_tool("gtm_map")
    dialog.features.set_fingerprint(None)
    dialog.features.use_numbers.setChecked(True)
    dialog.features.set_number_columns(["A", "B"])
    dialog.auto.setChecked(False)
    dialog.side.setValue(side)
    dialog.basis.setValue(3)
    _run(qtbot, dialog)
    return dialog


def test_training_writes_columns_and_the_map_and_opens_the_plot(chem, qtbot) -> None:  # noqa: F811
    from chemsuite.analysis.gtm import ASSET_KIND
    from chemsuite.app.plots.gtm import gtm_maps

    table = chem.workspace.table
    dialog = _train(chem, qtbot)
    assert set(COLUMNS) <= set(table.column_names)
    x = _numbers(table, "GTM X")
    assert np.isfinite(x).all() and (np.abs(x) <= 1 + 1e-6).all()
    maps = gtm_maps(table)
    assert len(maps) == 1 and table.assets[maps[0][0]].kind == ASSET_KIND
    plot = dialog.plot
    spec = plot.spec()
    assert (spec.kind, spec.gtm_map, spec.gtm_value) == ("gtm", maps[0][0], "pIC50")
    assert plot.title() == "GTM: pIC50"
    renderer = plot.renderer
    assert len(renderer.oids) == len(table.oids) and renderer._image is not None
    assert "5 × 5 nodes" in plot.summary.text()
    # the two groups (A high / B low against the reverse) sit apart on the map
    group = _numbers(table, "A") > 4
    a, b = renderer.x[group].mean(), renderer.x[~group].mean()
    c, d = renderer.y[group].mean(), renderer.y[~group].mean()
    assert np.hypot(a - b, c - d) > 0.5

    table.undo()  # one step takes the columns and the map away
    assert not set(COLUMNS) & set(table.column_names) and not gtm_maps(table)
    plot.rebuild()
    assert plot.summary.text().startswith("No map")
    table.redo()
    plot.rebuild()
    assert len(plot.renderer.oids) == len(table.oids)


def test_landscapes_select_and_hover(chem, qtbot) -> None:  # noqa: F811
    from PySide6.QtCore import Qt

    dialog = _train(chem, qtbot)
    plot = dialog.plot
    table = chem.workspace.table
    for value, note in ((None, ""), ("pIC50", "landscape of pIC50"), ("ID", "classes of ID")):
        plot.set_spec(gtm_value=value)
        assert plot.renderer._image is not None
        assert (note in plot.summary.text()) if note else ("landscape" not in plot.summary.text())
    plot.set_spec(gtm_value=None)
    assert plot.title() == "GTM: Density"
    renderer = plot.renderer
    x, y = renderer.x[0], renderer.y[0]
    plot._select(renderer.hit_point(x, y, 1e-3, 1e-3), Qt.KeyboardModifier.NoModifier)
    assert int(renderer.oids[0]) in chem.workspace.selection.oids
    hover = renderer.hover(table, x, y, 1e-3, 1e-3)
    assert hover is not None and hover.lines[-1].startswith("Log-likelihood: ")


def test_a_second_map_can_be_chosen(chem, qtbot) -> None:  # noqa: F811
    from chemsuite.app.plots.gtm import gtm_maps

    first = _train(chem, qtbot, 4).plot
    second = _train(chem, qtbot, 6).plot
    maps = gtm_maps(chem.workspace.table)
    assert len(maps) == 2 and second.spec().gtm_map == maps[1][0]
    assert "6 × 6 nodes" in second.summary.text()
    first.set_spec(gtm_map=maps[0][0])
    assert "4 × 4 nodes" in first.summary.text()


def test_the_map_and_plot_come_back_when_the_document_reopens(chem, qtbot, tmp_path) -> None:  # noqa: F811
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.app.plots.gtm import gtm_maps
    from chemsuite.core.table import Table

    dialog = _train(chem, qtbot)
    before = dialog.map
    path = tmp_path / "with_a_map.csd"
    assert chem._save_to(path)
    chem._install(Workspace(Table()))
    chem.open_path(path)
    qtbot.waitUntil(lambda: chem.workspace.path == path, timeout=20_000)
    qtbot.wait(20)
    table = chem.workspace.table
    maps = gtm_maps(table)
    assert len(maps) == 1
    plots = [d.widget() for d in chem.tool_panels.docks if d.widget().tool_id == "scatter"]
    gtm_plots = [p for p in plots if p.spec().kind == "gtm"]
    assert len(gtm_plots) == 1, "the plot reopens with the document"
    renderer = gtm_plots[0].renderer
    assert renderer.map is not None and renderer.map is not before  # decoded from the file
    np.testing.assert_allclose(renderer.map.position, before.position)
    assert gtm_plots[0].spec().gtm_value == "pIC50"


def test_automatic_nodes_follow_the_rows(chem) -> None:  # noqa: F811
    from chemsuite.analysis.gtm import default_basis, default_side

    dialog = chem.open_tool("gtm_map")
    assert dialog.auto.isChecked() and not dialog.side.isEnabled()
    side = default_side(len(dialog.rows.oids()))
    assert dialog.side.value() == side and dialog.basis.value() == default_basis(side)
