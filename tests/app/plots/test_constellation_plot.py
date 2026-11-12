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
"""The Constellation Plot: cores found in the background, drawn, selected and saved."""

from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import Qt

from tests.gui import gui

pytestmark = gui

SCAFFOLD = "Cn1c(-c2ccc(C(=O){})cc2)nc2ccccc21"  # four amides share its acyl core
OTHER = "O=C({})c1ccc2ccccc2c1"  # three naphthamides share theirs
ROWS = [
    ("A", SCAFFOLD.format("NC"), "5"),
    ("B", SCAFFOLD.format("NCC"), "6"),
    ("C", SCAFFOLD.format("NC3CC3"), "7"),
    ("D", SCAFFOLD.format("N(C)C"), ""),
    ("E", OTHER.format("NC"), "8"),
    ("F", OTHER.format("NCC"), "9"),
    ("G", OTHER.format("NCCC"), "9"),
    ("H", "CCO", "4"),
]


@pytest.fixture
def plotted(window, qtbot):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    names, smiles, pic50 = zip(*ROWS, strict=True)
    table.add_rows({"ID": names, "SMILES": smiles, "pIC50": pic50}, undoable=False)
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    panel = window.open_tool("constellation")
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    return window, panel


def ids(window, oids) -> list[str]:
    table = window.workspace.table
    return sorted(table.cell("ID", int(o)) for o in oids)


def test_cores_are_found_in_the_background_and_drawn(plotted) -> None:
    _window, panel = plotted
    con = panel.renderer.data.con
    assert sorted(con.n.tolist()) == [3, 4] and len(con.series_size) == 2
    assert panel.summary.text().startswith("2 cores · 2 analog series · 7 of 8 compounds")
    assert panel.title() == "Constellation: SMILES"
    assert panel.options.form.isRowVisible(panel.controls.structure)
    assert not panel.options.form.isRowVisible(panel.controls.core_summary)


def test_colour_summarises_each_cores_compounds(plotted) -> None:
    _window, panel = plotted
    panel.set_spec(color="pIC50", summary="median")
    renderer = panel.renderer
    four = int(np.flatnonzero(renderer.data.con.n == 4)[0])
    assert renderer.values[four] == 6.0  # D has no pIC50
    assert panel.title() == "Constellation: Median pIC50"
    assert panel.options.form.isRowVisible(panel.controls.core_summary)


def test_clicking_a_core_selects_its_compounds(plotted) -> None:
    window, panel = plotted
    renderer = panel.renderer
    three = int(np.flatnonzero(renderer.data.con.n == 3)[0])
    x, y = renderer.cx[three], renderer.cy[three]
    panel._select(renderer.hit_point(x, y, 1e-3, 1e-3), Qt.KeyboardModifier.NoModifier)
    assert ids(window, window.workspace.selection.oids) == ["E", "F", "G"]
    assert len(renderer._overlay.data) == 1  # the core holding them is ringed
    hover = renderer.hover(window.workspace.table, x, y, 1e-3, 1e-3)
    assert hover.lines[0].endswith("3 compounds") and hover.smiles == renderer.data.con.cores[three]
    panel.show_hover(hover)
    assert panel.card.isVisible()
    panel.card.hide()


def test_a_filter_regroups_without_finding_cores_again(plotted, qtbot, monkeypatch) -> None:
    import chemsuite.app.plots.constellation as module

    window, panel = plotted
    seen = []
    real = module.find_cores
    monkeypatch.setattr(
        module, "find_cores", lambda cells, **kw: seen.append(cells) or real(cells, **kw)
    )
    panel.set_spec(min_compounds=4)
    qtbot.waitUntil(lambda: panel.renderer.data is not None and len(panel.renderer.data.con) == 1)
    assert seen == [[]]  # every structure's cores were known


def test_the_options_are_saved_with_the_plot(plotted, qtbot) -> None:
    window, panel = plotted
    panel.set_spec(color="pIC50", summary="max", layout="pca", core_percent=50, labels=False)
    state = panel.state()
    assert state["kind"] == "constellation" and state["structure"] == "SMILES"
    again = window.open_tool("scatter")
    again.set_state(state)
    qtbot.waitUntil(lambda: again.renderer.data is not None, timeout=60_000)
    assert again.spec() == panel.spec()


def test_no_shared_core_says_so(window, qtbot) -> None:
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"SMILES": ["CCO", "c1ccccc1"]}, undoable=False)
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    panel = window.open_tool("constellation")
    qtbot.waitUntil(lambda: panel.renderer.data is not None, timeout=60_000)
    assert panel.summary.text().startswith("No core is shared by 2 or more of the 2 compounds")
