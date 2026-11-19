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
"""TMAP and Chemical Space Network plots: built in the background, drawn, selected and saved."""

from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import Qt

from tests.gui import gui

pytestmark = gui

ROWS = [
    ("A", "c1ccccc1C(=O)NC", "5"),
    ("B", "c1ccccc1C(=O)NCC", "6"),
    ("C", "c1ccccc1C(=O)NCCC", "7"),
    ("D", "c1ccccc1C(=O)NCCCC", ""),
    ("E", "CCCCCCCCO", "8"),
    ("F", "CCCCCCCCCO", "9"),
    ("G", "CCCCCCCCCCO", "9"),
    ("H", "c1ccncc1", "4"),
]


def open_plot(window, qtbot, kind: str):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    names, smiles, pic50 = zip(*ROWS, strict=True)
    table.add_rows({"ID": names, "SMILES": smiles, "pIC50": pic50}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    window._install(Workspace(table))
    panel = window.open_tool(kind)
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    return panel


def ids(window, oids) -> list[str]:
    table = window.workspace.table
    return sorted(table.cell("ID", int(o)) for o in oids)


@pytest.fixture
def network(window, qtbot):
    panel = open_plot(window, qtbot, "network")
    panel.set_spec(threshold=0.4)
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    return window, panel


def test_a_network_joins_alike_compounds(network) -> None:
    _window, panel = network
    graph = panel.renderer.data.graph
    assert len(graph.a) > 0 and (graph.similarity >= 0.4).all()
    assert len(panel.renderer.oids) == 8  # singletons are drawn too
    assert "edges" in panel.summary.text() and "singletons" in panel.summary.text()
    assert panel.title() == "Network: SMILES"
    form, controls = panel.options.form, panel.controls
    assert form.isRowVisible(controls.threshold) and not form.isRowVisible(controls.neighbours)
    assert not form.isRowVisible(controls.fingerprint)


def test_clicking_a_node_selects_its_row_and_hover_shows_it(network) -> None:
    window, panel = network
    renderer = panel.renderer
    i = int(np.flatnonzero(renderer.degree > 0)[0])
    x, y = renderer.x[i], renderer.y[i]
    panel._select(renderer.hit_point(x, y, 1e-3, 1e-3), Qt.KeyboardModifier.NoModifier)
    assert window.workspace.selection.oids == {int(renderer.oids[i])}
    hover = renderer.hover(window.workspace.table, x, y, 1e-3, 1e-3)
    assert hover.oid == int(renderer.oids[i]) and hover.lines[-1].startswith("Neighbours: ")


def test_box_and_lasso_select_rows(network) -> None:
    window, panel = network
    renderer = panel.renderer
    lo, hi = renderer.x.min() - 1, renderer.x.max() + 1
    bottom, top = renderer.y.min() - 1, renderer.y.max() + 1
    every = ids(window, renderer.oids)
    assert ids(window, renderer.hit_box(lo, bottom, hi, top)) == every
    square = np.array([[lo, bottom], [hi, bottom], [hi, top], [lo, top]])
    assert ids(window, renderer.hit_lasso(square)) == every


def test_hiding_singletons_and_colouring_edges(network, qtbot) -> None:
    _window, panel = network
    panel.set_spec(singletons=False, color="pIC50", edge_color="difference")
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    renderer = panel.renderer
    assert (renderer.degree > 0).all()
    assert len(renderer._edges) >= 1
    assert panel.title() == "Network: pIC50"


def test_matched_pair_edges(network, qtbot) -> None:
    window, panel = network
    panel.set_spec(edges="mmp")
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    graph = panel.renderer.data.graph
    pairs = {
        tuple(ids(window, panel.renderer.data.oids[[a, b]]))
        for a, b in zip(graph.a, graph.b, strict=True)
    }
    assert ("A", "B") in pairs and np.isnan(graph.similarity).all()
    assert not panel.options.form.isRowVisible(panel.controls.threshold)


def test_too_many_edges_says_so(network, qtbot, monkeypatch) -> None:
    import chemsuite.app.plots.graph as module
    from chemsuite.analysis.fingerprints import TooManyPairs

    def crowded(*_args, **_kwargs):
        raise TooManyPairs(60_000)

    monkeypatch.setattr(module, "similarity_edges", crowded)
    _window, panel = network
    panel.set_spec(threshold=0.3)
    qtbot.waitUntil(lambda: "raise the threshold" in panel.summary.text(), timeout=60_000)
    assert panel.summary.text().startswith("⚠ Over 60,000 pairs at 0.30")


def test_tmap_draws_a_tree(window, qtbot) -> None:
    panel = open_plot(window, qtbot, "tmap")
    graph = panel.renderer.data.graph
    trees = graph.component.max() + 1
    assert len(graph.a) == 8 - trees
    assert panel.summary.text().startswith("8 compounds · ")
    assert panel.options.form.isRowVisible(panel.controls.neighbours)
    assert panel.options.form.isRowVisible(panel.controls.fingerprint)
    assert not panel.options.form.isRowVisible(panel.controls.threshold)


def test_tmap_on_mhfp6_is_saved_and_reopens(window, qtbot) -> None:
    panel = open_plot(window, qtbot, "tmap")
    before = panel.renderer.data
    panel.set_spec(fingerprint="mhfp6")
    qtbot.waitUntil(
        lambda: panel.renderer.data is not None and panel.renderer.data is not before,
        timeout=60_000,
    )
    graph = panel.renderer.data.graph
    assert len(graph.a) == 8 - (graph.component.max() + 1)
    state = panel.state()
    assert state["fingerprint"] == "mhfp6"
    again = window.open_tool("scatter")
    again.set_state(state)
    qtbot.waitUntil(lambda: again.renderer.data is not None, timeout=60_000)
    assert again.spec() == panel.spec()
    again.set_state({**state, "fingerprint": "no such"})  # damaged: back to ECFP4
    assert again.spec().fingerprint == "ecfp4"


def test_the_options_are_saved_with_the_plot(network, qtbot) -> None:
    window, panel = network
    panel.set_spec(edges="similarity", threshold=0.55, edge_color="similarity", singletons=False)
    state = panel.state()
    assert state["kind"] == "network" and state["threshold"] == 0.55
    again = window.open_tool("scatter")
    again.set_state(state)
    qtbot.waitUntil(lambda: again.renderer.data is not None, timeout=60_000)
    assert again.spec() == panel.spec()


def test_a_filter_finds_the_graph_again(network, qtbot) -> None:
    from chemsuite.core.table.filters import RangeFilter

    window, panel = network
    before = panel.renderer.data
    window.workspace.table.set_filter("t", RangeFilter("pIC50", 5, 9))
    qtbot.waitUntil(
        lambda: panel.renderer.data is not None and panel.renderer.data is not before,
        timeout=60_000,
    )
    assert len(panel.renderer.data.oids) == 6  # D has no pIC50, H is 4


def test_tmap_keeps_its_mhfp6_values(window, qtbot, monkeypatch) -> None:
    import chemsuite.app.plots.graph as graph

    panel = open_plot(window, qtbot, "tmap")
    panel.set_spec(fingerprint="mhfp6")
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    calls = []
    real = graph.compute_mhfp
    monkeypatch.setattr(
        graph,
        "compute_mhfp",
        lambda cells, *a, **kw: calls.append(len(cells)) or real(cells, *a, **kw),
    )
    panel.set_spec(neighbours=4)
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    assert calls == []  # the first build kept them


def test_tmap_keeps_its_fingerprints(window, qtbot, monkeypatch) -> None:
    import chemsuite.analysis.fingerprints as fingerprints

    panel = open_plot(window, qtbot, "tmap")
    calls = []
    real = fingerprints.compute_fingerprints
    monkeypatch.setattr(
        fingerprints,
        "compute_fingerprints",
        lambda cells, *a, **kw: calls.append(len(cells)) or real(cells, *a, **kw),
    )
    panel.set_spec(neighbours=4)
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    assert calls == []  # the first build stored them
