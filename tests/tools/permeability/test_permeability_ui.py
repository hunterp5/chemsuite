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

"""Chemistry ▸ Calculate ▸ Permeability: the dialog and the columns it writes."""

from __future__ import annotations

import pytest

from chemsuite.permeability import ENDPOINTS
from chemsuite.permeability import models as shipped
from tests.gui import gui

pytestmark = gui

#: propranolol (highly permeable), atenolol (poorly permeable), verapamil (a P-gp substrate),
#: and a blank row
SMILES = [
    "CC(C)NCC(O)COc1cccc2ccccc12",
    "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
    "",
]

installed = pytest.mark.skipif(
    not all(shipped.installed(e) for e in ENDPOINTS),
    reason="the permeability models are not built in this checkout",
)


@pytest.fixture
def chem(window):
    from chemsuite.app.document.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows({"ID": [f"C{i}" for i in range(len(SMILES))], "SMILES": SMILES}, undoable=False)
    table.convert_column("SMILES", "structure")
    table.history.mark_clean()
    window._install(Workspace(table))
    return window


@installed
def test_the_dialog_offers_every_endpoint_and_checks_them(chem) -> None:
    dialog = chem.open_tool("predict_permeability")
    assert dialog.primary.text() == "Predict for 4 Rows"
    # PAMPA was measured and held back; see chemsuite/permeability/__init__.py
    assert set(dialog.endpoints) == {"mdck_er", "caco2"}
    assert [e.id for e in dialog.chosen_endpoints()] == ["mdck_er", "caco2"]
    assert dialog.chosen() == ["prediction", "interval"]


@installed
def test_a_run_needs_an_endpoint(chem) -> None:
    dialog = chem.open_tool("predict_permeability")
    for box in dialog.endpoints.values():
        box.setChecked(False)
    assert "at least one endpoint" in dialog.status.text()
    assert not dialog.primary.isEnabled()


def test_the_dialog_says_so_when_no_model_is_installed(chem, monkeypatch) -> None:
    monkeypatch.setattr(shipped, "installed", lambda _endpoint: False)
    dialog = chem.open_tool("predict_permeability")
    assert "No permeability model is installed" in dialog.status.text()
    assert not dialog.primary.isEnabled()
    # nothing to choose from, and nothing offered as checked
    assert not dialog.chosen_endpoints()


@installed
def test_the_dialog_shows_no_helper_text(chem) -> None:
    """The kit's rule: the only text a dialog shows is a problem that blocks running."""
    dialog = chem.open_tool("predict_permeability")
    assert not dialog.status.text().strip()


@installed
def test_it_writes_a_prediction_and_an_interval_per_endpoint(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_permeability")
    dialog.outputs["domain"].setChecked(True)
    dialog.run()
    qtbot.waitUntil(lambda: "Caco-2 Papp Predicted" in table.column_names, timeout=120_000)

    assert [name for name in table.column_names if name.startswith("Caco-2 Papp")] == [
        "Caco-2 Papp Predicted",
        "Caco-2 Papp Low",
        "Caco-2 Papp High",
        "Caco-2 Papp Distance",
        "Caco-2 Papp In Domain",
    ]
    assert "MDCK ER Predicted" in table.column_names

    propranolol, atenolol, _verapamil, blank = (int(oid) for oid in table.oids)
    # the interval brackets the prediction, and a blank row gets nothing
    for oid in (propranolol, atenolol):
        low = float(table.cell("Caco-2 Papp Low", oid))
        high = float(table.cell("Caco-2 Papp High", oid))
        assert low <= float(table.cell("Caco-2 Papp Predicted", oid)) <= high
    assert table.cell("Caco-2 Papp Predicted", blank) == ""

    # the ordering the assay is known for: propranolol crosses a Caco-2 monolayer, atenolol
    # does not. This is the model's sanity check, not its evaluation.
    assert float(table.cell("Caco-2 Papp Predicted", propranolol)) > float(
        table.cell("Caco-2 Papp Predicted", atenolol)
    )


@installed
def test_an_endpoint_not_asked_for_writes_nothing(chem, qtbot) -> None:
    table = chem.workspace.table
    dialog = chem.open_tool("predict_permeability")
    for endpoint_id, box in dialog.endpoints.items():
        box.setChecked(endpoint_id == "caco2")
    dialog.run()
    qtbot.waitUntil(lambda: "Caco-2 Papp Predicted" in table.column_names, timeout=120_000)

    assert not [name for name in table.column_names if name.startswith("MDCK ER")]


def test_a_classification_endpoint_would_write_a_label_and_a_class_set() -> None:
    """No classification endpoint ships today — PAMPA was measured and held back — but the
    code path is live for the day one does, so the column contract is pinned here rather
    than left to be discovered later."""
    from chemsuite.tools.qsar.columns import names_for

    assert names_for("PAMPA", ["prediction", "interval"], "classification") == [
        "PAMPA Predicted",
        "PAMPA Confidence",
        "PAMPA Possible",
    ]
    assert names_for("Caco-2 Papp", ["prediction", "interval"], "regression") == [
        "Caco-2 Papp Predicted",
        "Caco-2 Papp Low",
        "Caco-2 Papp High",
    ]
