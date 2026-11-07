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

WAX = "C" * 26 + "C(=O)O"  # cLogP ≈ 9.7: too lipophilic for the egg
SMILES = ["Cn1c(=O)c2c(ncn2C)n(C)c1=O", "CC(=O)Oc1ccccc1C(=O)O", WAX, ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {"ID": ["caffeine", "aspirin", "wax", "blank"], "SMILES": SMILES,
         "pIC50": ["6", "7", "8", "5"]},
        undoable=False,
    )  # fmt: skip
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _cells(table, name):
    return [table.cell(name, int(o)) for o in table.oids]


def test_boiled_egg_plot(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("boiled_egg")
    assert not dialog.status.text().strip()
    dialog.accept()
    qtbot.waitUntil(lambda: dialog.plot is not None, timeout=20_000)
    assert _cells(table, "BOILED-Egg") == ["HIA (white)", "BBB (yolk)", "Outside", ""]
    plot = dialog.plot
    assert plot.spec().guide == "boiled_egg" and plot.spec().color == "BOILED-Egg"
    assert plot.spec().x == "TPSA (N,O,S,P)" and plot.spec().y == "cLogP"
    assert len(plot.renderer._guide_items) == 2
    assert plot.state()["guide"] == "boiled_egg"
    # Descriptors now exist, so a second plot needs no job.
    second = chem.open_tool("golden_triangle")
    assert not second.status.text().strip()
    second.accept()
    qtbot.waitUntil(lambda: second.plot is not None, timeout=20_000)
    assert _cells(table, "Golden Triangle")[:3] == ["Outside"] * 3  # MW < 200; too lipophilic


def test_lipe_plot(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("lipe")
    assert dialog.activity.currentText() == "pIC50"
    dialog.accept()
    qtbot.waitUntil(lambda: dialog.plot is not None, timeout=20_000)
    lipe = _cells(table, "LipE")
    logp = _cells(table, "cLogP")
    assert float(lipe[1]) == pytest.approx(7 - float(logp[1]), abs=0.011)
    assert dialog.plot.spec().y == "pIC50" and dialog.plot.spec().guide == "lipe"
    # Switching the guide off removes the lines.
    dialog.plot.set_spec(guide=None)
    assert dialog.plot.renderer._guide_items == []
