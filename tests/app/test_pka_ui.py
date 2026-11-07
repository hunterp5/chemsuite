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

"""Chemistry ▸ Predict pKa."""

from __future__ import annotations

import pytest

from .conftest import gui

pytestmark = gui

SMILES = ["CC(=O)O", "C1CCNCC1", "NCC(=O)O", "c1ccccc1", ""]


@pytest.fixture
def chem(window):
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": [f"C{i}" for i in range(len(SMILES))], "SMILES": SMILES}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def test_predict_pka_writes_the_chosen_columns(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_pka")
    assert dialog.primary.text() == "Predict pKa for 5 Rows"
    assert not dialog.status.text().strip()
    dialog.ph.setValue(7.0)
    dialog.outputs["species"].setChecked(True)
    assert not dialog.status.text().strip()
    for box in dialog.outputs.values():
        box.setChecked(False)
    assert "Choose at least one" in dialog.status.text() and not dialog.primary.isEnabled()
    dialog.outputs["pka"].setChecked(True)
    dialog.outputs["species"].setChecked(True)
    dialog.outputs["charge"].setChecked(True)
    dialog.run()
    qtbot.waitUntil(lambda: "Major Species @ pH 7" in table.column_names, timeout=20_000)
    first, piperidine, glycine, benzene, blank = (int(o) for o in table.oids)
    assert float(table.cell("pKa (most acidic)", first)) == pytest.approx(4.47, abs=0.01)
    assert table.cell("pKa (most basic)", benzene) == ""
    assert float(table.cell("Charge @ pH 7", piperidine)) == pytest.approx(1.0, abs=0.01)
    assert table.cell("Major Species @ pH 7", glycine) == "[NH3+]CC(=O)[O-]"
    assert table.kind("Major Species @ pH 7") == "structure"
    assert table.cell("Major Species @ pH 7", blank) == ""


def test_every_run_writes_new_columns_with_one_suffix(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_pka")
    assert not hasattr(dialog, "if_exists")
    for key, box in dialog.outputs.items():
        box.setChecked(key in ("pka", "logd"))
    dialog.run()
    qtbot.waitUntil(lambda: "cLogD 7.4" in table.column_names, timeout=20_000)
    first = int(table.oids[0])
    # an all-blank column from the first run is still not reused...
    table.set_cells("pKa (most acidic)", table.oids, [""] * len(table.oids))
    # ...and a free name still takes the run's shared suffix
    table.remove_column("pKa (most basic)")
    dialog = chem.open_tool("predict_pka")
    for key, box in dialog.outputs.items():
        box.setChecked(key in ("pka", "logd"))
    dialog.run()
    qtbot.waitUntil(lambda: "cLogD 7.4 (2)" in table.column_names, timeout=20_000)
    assert "pKa (most acidic) (2)" in table.column_names
    assert "pKa (most basic) (2)" in table.column_names
    assert "pKa (most basic)" not in table.column_names
    assert table.cell("pKa (most acidic)", first) == ""  # the first run's column is untouched
    assert float(table.cell("pKa (most acidic) (2)", first)) == pytest.approx(4.47, abs=0.01)
