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
"""Minimising in the workbench: Prepare's Minimise and Protein ▸ Minimise…, each result kept
as a new protein, run here on the fake OpenMM worker."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from chemsuite.app.document.tools import ToolContext
from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure import edit, store
from chemsuite.macromolecule.ui.report import TOOL_ID as REPORT_ID
from chemsuite.macromolecule.ui.report import PrepareReportPanel
from chemsuite.simulate import runner
from chemsuite.simulate.config import Interpreter
from chemsuite.simulate.spec import FF19SB, VACUUM, MinimiseSettings
from chemsuite.tools.macromolecule import minimise, prepare_dialog, tool
from chemsuite.tools.macromolecule.minimise_dialog import MinimiseDialog
from chemsuite.tools.macromolecule.minimise_options import MinimiseFields
from tests.gui import gui
from tests.simulate.test_runner import FAKE

pytestmark = gui
DATA = Path(__file__).resolve().parents[2] / "data" / "protein"


@pytest.fixture
def fake(monkeypatch):
    """OpenMM 'set up', its worker the fake one."""
    worker = runner.Simulator(Interpreter(sys.executable))
    worker.command = [sys.executable, "-I", "-c", FAKE]
    worker.env = None
    monkeypatch.setattr(minimise, "ready", lambda: "")
    monkeypatch.setattr(minimise, "simulator", lambda: worker)
    yield worker
    worker.close()


@pytest.fixture
def opened(window, qtbot, monkeypatch):
    """Crambin opened from a file (small, and prepared in a second)."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    job = tool.open_protein(ToolContext(window), [str(DATA / "1crn.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    return window


def wait(qtbot, job) -> None:
    assert job is not None
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)


def prepare(window, qtbot, minimising: MinimiseSettings | None = None) -> None:
    ctx = ToolContext(window)
    wait(qtbot, prepare_dialog.prepare_protein(ctx, "1CRN", pipeline.Options(), minimising))


def shown(window, protein: str) -> bool:
    display = tool.macromolecule_workspace(ToolContext(window)).session().display(protein)
    assert display is not None
    return display.visible


def test_minimise_is_a_tool_and_in_the_viewers_tools_menu(opened) -> None:
    assert tool.MINIMISE_ID in [t.id for t in tool.TOOLS]
    panel = ToolContext(opened).find_panel("protein")
    menu = next(a.menu() for a in panel.menu_bar.actions() if a.text() == "Tools")
    menu.aboutToShow.emit()
    action = next(a for a in menu.actions() if a.text() == "Minimise…")
    assert action.isEnabled()
    action.trigger()
    dialog = next(d for d in panel.window().findChildren(MinimiseDialog))
    assert dialog.protein.currentData() == "1CRN"
    dialog.close()


def test_without_openmm_minimise_cannot_be_ticked_or_run(opened, monkeypatch) -> None:
    monkeypatch.setattr(minimise, "ready", lambda: "OpenMM is not set up")
    ctx = ToolContext(opened)
    dialog = prepare_dialog.PrepareDialog(ctx)
    assert not dialog.minimise.isEnabled() and not dialog.minimise.isChecked()
    assert dialog.minimise_settings() is None
    assert minimise.problem(ctx, "1CRN", MinimiseSettings()) == "OpenMM is not set up"
    dialog.close()


def test_minimise_needs_a_full_repair(opened, fake) -> None:
    dialog = prepare_dialog.PrepareDialog(ToolContext(opened))
    assert dialog.minimise.isEnabled()
    dialog.minimise.setChecked(True)
    assert dialog.minimise_settings() == minimise.last_settings()
    dialog.fields.caps.setChecked(False)
    assert not dialog.minimise.isEnabled() and dialog.minimise_settings() is None
    assert "Cap Ends" in dialog.minimise.toolTip()
    dialog.close()


def test_the_fields_give_the_settings_they_were_given(qtbot) -> None:
    from chemsuite.app.widgets.dialogs import ToolDialog

    dialog = ToolDialog(None, "Test")
    qtbot.addWidget(dialog)
    settings = MinimiseSettings(FF19SB, VACUUM, 2.0, 0.45)
    fields = MinimiseFields(dialog.section("Method"), settings)
    assert fields.settings() == settings
    assert MinimiseFields(dialog.section("Other")).settings() == MinimiseSettings()


def test_only_a_prepared_protein_is_minimised(opened, fake) -> None:
    ctx = ToolContext(opened)
    assert minimise.problem(ctx, "1CRN", MinimiseSettings()) == "Prepare 1CRN first"
    dialog = MinimiseDialog(ctx)
    assert dialog.problem() == "Prepare 1CRN first"
    assert minimise.minimise_protein(ctx, "1CRN", MinimiseSettings()) is None
    dialog.close()


def test_prepare_and_minimise_adds_two_proteins_in_one_step(opened, qtbot, fake) -> None:
    table = opened.workspace.table
    before = store.as_structure(table.assets["1CRN"].payload)
    prepare(opened, qtbot, MinimiseSettings())
    assert list(table.assets) == ["1CRN", "1CRN-prep", "1CRN-min"]
    assert store.as_structure(table.assets["1CRN"].payload).n_atoms == before.n_atoms
    assert minimise.stage_of(ToolContext(opened), "1CRN-min") == minimise.MINIMISED
    assert table.assets["1CRN-min"].meta[minimise.DERIVED_FROM] == "1CRN-prep"
    # the newest copy is shown in place of the protein it came from
    assert not shown(opened, "1CRN") and not shown(opened, "1CRN-prep")
    assert shown(opened, "1CRN-min")
    workspace = tool.macromolecule_workspace(ToolContext(opened))
    record = workspace.prepared("1CRN-min")
    assert record is not None and record.minimised is not None
    assert record.minimised.source == "1CRN-prep"
    assert workspace.prepared("1CRN-prep") is not None
    panel = ToolContext(opened).find_panel(REPORT_ID)
    assert isinstance(panel, PrepareReportPanel) and panel.report.protein == "1CRN-min"
    assert any(row.text.startswith("Heavy atoms moved") for row in panel.report.rows)
    assert table.history.undo_label == "Prepare and minimise 1CRN"
    table.undo()
    assert list(table.assets) == ["1CRN"] and shown(opened, "1CRN")


def test_a_minimisation_that_cannot_run_still_adds_the_prepared_protein(
    opened, qtbot, fake, monkeypatch
) -> None:
    def broken(*_a, **_k):
        raise runner.SimulateError("ff14SB has no template for A/TYR 29\nmore")

    monkeypatch.setattr(minimise, "relax", broken)
    prepare(opened, qtbot, MinimiseSettings())
    assert list(opened.workspace.table.assets) == ["1CRN", "1CRN-prep"]
    report = tool.macromolecule_workspace(ToolContext(opened)).report("1CRN-prep")
    assert report is not None
    first = report.rows[0]
    assert first.warn and first.text == "Not minimised: ff14SB has no template for A/TYR 29"


def test_minimise_adds_the_minimised_protein_beside_the_prepared_one(opened, qtbot, fake) -> None:
    ctx = ToolContext(opened)
    prepare(opened, qtbot)
    assert minimise.problem(ctx, "1CRN", MinimiseSettings()) == ""
    wait(qtbot, minimise.minimise_protein(ctx, "1CRN", MinimiseSettings()))
    table = opened.workspace.table
    assert list(table.assets) == ["1CRN", "1CRN-min"]
    assert table.history.undo_label == "Minimise 1CRN"
    s = store.as_structure(table.assets["1CRN-min"].payload)
    workspace = tool.macromolecule_workspace(ctx)
    record = workspace.prepared("1CRN-min")
    assert record is not None and record.matches(s) and record.minimised is not None
    report = workspace.report("1CRN-min")
    assert report is not None and report.rows[0].text.startswith("Energy")
    # never minimised twice: the restraints would hold it to the first result, not the deposit
    assert minimise.problem(ctx, "1CRN-min", MinimiseSettings()) == (
        "1CRN-min is minimised already: minimise 1CRN instead"
    )
    dialog = prepare_dialog.PrepareDialog(ctx)
    dialog.protein.setCurrentIndex(dialog.protein.findData("1CRN-min"))
    assert dialog.problem() == "1CRN-min is minimised already: prepare 1CRN instead"
    dialog.close()
    # a second minimisation of the prepared protein is a new copy
    wait(qtbot, minimise.minimise_protein(ctx, "1CRN", MinimiseSettings()))
    assert list(table.assets) == ["1CRN", "1CRN-min", "1CRN-min-2"]


def test_an_edit_after_preparing_asks_for_a_new_preparation(opened, qtbot, fake) -> None:
    ctx = ToolContext(opened)
    prepare(opened, qtbot)
    s = prepare_dialog.structure_of(ctx, "1CRN")
    assert s is not None
    tool.macromolecule_workspace(ctx).session().edit(
        "Cut", {"1CRN": edit.delete(s, s.name == "OXT")}
    )
    assert minimise.problem(ctx, "1CRN", MinimiseSettings()) == (
        "1CRN was edited after it was prepared: prepare it again"
    )
