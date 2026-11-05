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

import pytest

from .conftest import gui

pytestmark = gui

R1 = {"H": "", "Cl": "Cl", "C": "C", "OC": "OC"}
R2 = {"H": "", "F": "F", "N": "N"}
EFFECT1 = {"H": 0.0, "Cl": 1.0, "C": 0.5, "OC": -0.3}
EFFECT2 = {"H": 0.0, "F": 0.7, "N": -0.4}
CORE = "[*:1]c1ccc(C(N)=O)cc1[*:2]"


@pytest.fixture
def series(window):
    """Benzamides with additive activity; Cl/F is not made. Plus one compound off-series."""
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    smiles, ids, activity = [], [], []
    for a, sa in R1.items():
        for b, sb in R2.items():
            if (a, b) == ("Cl", "F"):
                continue
            smiles.append(f"c1cc({sa})c({sb})cc1C(=O)N".replace("()", ""))
            ids.append(f"{a}/{b}")
            activity.append(f"{5 + EFFECT1[a] + EFFECT2[b]:.2f}")
    smiles.append("CCO")
    ids.append("ethanol")
    activity.append("4.0")
    table = Table()
    table.add_rows({"ID": ids, "SMILES": smiles, "pIC50": activity}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.convert_column("pIC50", "number")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


def _cells(table, name):
    return [table.cell(name, int(o)) for o in table.oids]


def _run(qtbot, dialog):
    dialog.accept()
    qtbot.waitUntil(lambda: not dialog.job.running, timeout=30_000)


def test_rgroups_then_matrix(series, qtbot) -> None:
    table = series.workspace.table
    dialog = series.open_tool("rgroups")
    assert not dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()
    dialog.core.core.setText(CORE)
    _run(qtbot, dialog)
    r1 = _cells(table, "R1")
    assert set(r1[:-1]) == {"H", "*Cl", "*C", "*OC"} and r1[-1] == ""
    assert _cells(table, "Core Match")[-1] == ""

    panel = series.open_tool("sar_matrix")
    assert panel.rows_axis.currentText() == "R1" and panel.columns_axis.currentText() == "R2"
    assert panel.activity.currentText() == "pIC50"
    m = panel.matrix
    assert m.row_labels[0] == "H" and len(m.row_labels) == 4 and len(m.column_labels) == 3
    r, c = m.row_labels.index("*Cl"), m.column_labels.index("*F")
    gap = panel.grid.item(r, c)
    assert gap.font().italic() and math.isclose(float(gap.text()), 6.7, abs_tol=0.05)
    h = m.row_labels.index("H")
    selected = panel.select_cell(h, 0)
    assert [table.cell("ID", o) for o in selected] == ["H/H"]
    assert "*Cl" in panel.copy_text()
    state = panel.state()
    assert state["rows"] == "R1"
    panel.set_state({**state, "aggregate": "max", "predict": False})
    assert panel.grid.item(r, c) is None or panel.grid.item(r, c).text() == ""
    # The matrix follows the table: deleting a row empties its cell.
    table.remove_rows([table.oids[0]])
    panel.rebuild()
    assert panel.grid.item(h, 0) is None or panel.grid.item(h, 0).text() == ""
    panel.close_panel()


def test_suggest_core_from_selection(series) -> None:
    table = series.workspace.table
    series.workspace.selection.set(table.oids[:-1].tolist(), "test")
    dialog = series.open_tool("rgroups")
    smarts = dialog.core.suggest()
    assert smarts and dialog.core.problem() == ""
    dialog.reject()


def test_free_wilson(series, qtbot) -> None:
    table = series.workspace.table
    dialog = series.open_tool("free_wilson")
    dialog.core.core.setText(CORE)
    assert dialog.activity.currentText() == "pIC50"
    _run(qtbot, dialog)
    panel = dialog.panel
    assert panel.model.fit.r2 > 0.99
    predicted = _cells(table, "FW Predicted")
    assert predicted[-1] == "" and math.isclose(float(predicted[0]), 5.0, abs_tol=0.01)
    best = panel.unmade[0]
    assert best[0] == ("*Cl", "*F") and math.isclose(best[1], 6.7, abs_tol=0.01)
    before = len(table)
    new = panel.add_to_table([best])
    assert len(new) == 1 and len(table) == before + 1
    assert table.cell("Source", int(new[0])) == "Free-Wilson (not made)"
    assert table.cell("SMILES", int(new[0])) == "NC(=O)c1ccc(Cl)c(F)c1"
    panel.terms.selectRow(next(i for i, t in enumerate(panel.model.terms) if t.label == "*OC"))
    carriers = panel.select_carriers()
    assert sorted(table.cell("ID", o) for o in carriers) == ["OC/F", "OC/H", "OC/N"]


def test_hansch(series, qtbot, tmp_path) -> None:
    table = series.workspace.table
    dialog = series.open_tool("hansch")
    dialog.core.core.setText(CORE)
    path = tmp_path / "constants.csv"
    dialog.save_template(str(path))
    assert "smiles,pi" in path.read_text()
    assert dialog.problem() == ""  # a header-only file is chosen; loading reports it
    dialog.accept()
    assert "no substituents" in dialog.status.text()
    with path.open("a") as handle:  # made-up π values that track the R1 effects
        handle.write("*Cl,1.0\n*C,0.5\n*OC,-0.3\n*F,0.7\n*N,-0.4\n")
    dialog.electronic.setCurrentIndex(dialog.electronic.findData(""))
    _run(qtbot, dialog)
    model = dialog.panel.model
    assert model.terms == ["pi"] and model.fit.r2 > 0.99
    assert math.isclose(model.fit.coefficients[0], 1.0, abs_tol=1e-6)
    assert _cells(table, "Hansch Predicted")[-1] == ""


def test_mmp(series, qtbot) -> None:
    table = series.workspace.table
    dialog = series.open_tool("mmp")
    assert dialog.activity.currentText() == "pIC50"
    dialog.min_pairs.setValue(1)
    _run(qtbot, dialog)
    panel = dialog.panel
    found = {(t.from_group, t.to_group): t for t in panel.found}
    transform = found[("[*:1][H]", "Cl[*:1]")]  # H → Cl, with R2 = H and N (not F)
    assert transform.count == 2 and math.isclose(transform.mean, 1.0, abs_tol=1e-6)
    row = next(
        r for r in range(panel.transforms.rowCount()) if panel.transforms.row_data(r) is transform
    )
    panel.transforms.selectRow(row)
    assert panel.pairs.rowCount() == 2
    ids = sorted(table.cell("ID", o) for o in panel.select_rows())
    assert ids == ["Cl/H", "Cl/N", "H/H", "H/N"]
    assert "Pairs" in panel.transforms.copy_text()


def test_sali(series, qtbot) -> None:
    table = series.workspace.table
    dialog = series.open_tool("sali")
    dialog.min_similarity.setValue(0.3)
    _run(qtbot, dialog)
    sali = _cells(table, "Max SALI")
    partners = _cells(table, "Cliff Partner")
    assert all(s for s in sali[:-1]) and sali[-1] == "" and partners[-1] == ""
    panel = dialog.panel
    assert panel.table.rowCount() > 0
    top = panel.table.row_data(0)
    panel.table.selectRow(0)
    oids = panel.select_pairs()
    assert len(oids) == 2 and sorted(oids) == sorted(table.oids[list(top)].tolist())


def test_sar_panels_restore(series, qtbot) -> None:
    dialog = series.open_tool("rgroups")
    dialog.core.core.setText(CORE)
    _run(qtbot, dialog)
    series.open_tool("sar_matrix").aggregate.setCurrentIndex(2)
    saved = series.panel_states()
    assert [s["tool"] for s in saved] == ["sar_matrix"]
    series._restore_panels(saved)
    restored = [d.widget() for d in series._panels if d.widget().tool_id == "sar_matrix"]
    assert len(restored) == 2 and restored[-1].aggregate.currentIndex() == 2
