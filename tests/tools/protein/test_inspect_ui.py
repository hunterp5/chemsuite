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
"""Protein ▸ Inspect…: a preparation surveyed, overruled line by line and applied."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemsuite.app.document.tools import ToolContext
from chemsuite.protein.prepare import inspect as inspection
from chemsuite.protein.prepare import pipeline, sites
from chemsuite.protein.structure import edit, store
from chemsuite.protein.ui.inspector import TOOL_ID as INSPECT_ID
from chemsuite.protein.ui.inspector import InspectorPanel
from chemsuite.protein.ui.report import PrepareReportPanel
from chemsuite.tools.protein import inspect_dialog, tool
from tests.gui import gui

pytestmark = gui
DATA = Path(__file__).resolve().parents[2] / "data" / "protein"
HYDROGEN = 1


@pytest.fixture
def opened(window, qtbot, monkeypatch):
    """Crambin, missing an arginine's tip so there is a repair to leave alone, in the viewer."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    ctx = ToolContext(window)
    job = tool.open_protein(ctx, [str(DATA / "1crn.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    s = store.as_structure(window.workspace.table.assets["1CRN"].payload)
    r = next(i for i in range(s.n_residues) if s.residue_label(i) == "A/ARG 10")
    span = s.atoms_of(r)
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[[i for i in range(span.start, span.stop) if s.name[i] in {"CZ", "NH1", "NH2"}]] = True
    tool.protein_workspace(ctx).session().edit("Cut", {"1CRN": edit.delete(s, gone)})
    tool.open_viewer(ctx)
    return window


def surveyed(window, qtbot, options: pipeline.Options | None = None) -> InspectorPanel:
    ctx = ToolContext(window)
    job = inspect_dialog.inspect_protein(ctx, "1CRN", options or pipeline.Options())
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)
    panel = ctx.find_panel(INSPECT_ID)
    assert isinstance(panel, InspectorPanel)
    return panel


def structure(window):
    return store.as_structure(window.workspace.table.assets["1CRN"].payload)


def test_inspect_is_in_the_viewers_tools_menu(opened) -> None:
    assert InspectorPanel.tool_id == tool.INSPECT_ID
    panel = ToolContext(opened).find_panel(tool.TOOL_ID)
    menu = next(a.menu() for a in panel.menu_bar.actions() if a.text() == "Tools")
    menu.aboutToShow.emit()
    action = next(a for a in menu.actions() if a.text() == "Inspect…")
    assert action.isEnabled()
    action.trigger()
    dialog = next(d for d in panel.window().findChildren(inspect_dialog.InspectDialog))
    assert dialog.verb() == "Inspect 1CRN" and dialog.windowTitle() == "Inspect Protein"
    dialog.close()


def test_surveying_commits_nothing_and_docks_in_the_viewer(opened, qtbot) -> None:
    table = opened.workspace.table
    before = structure(opened)
    undo = table.history.undo_label
    panel = surveyed(opened, qtbot)
    viewer = ToolContext(opened).find_panel(tool.TOOL_ID)
    assert panel.window() is viewer.window()
    after = structure(opened)
    assert after.n_atoms == before.n_atoms and table.history.undo_label == undo
    kinds = {line.data(0, 0x0100).kind for line in panel.lines()}  # LINE role
    assert inspection.MISSING in kinds and sites.AMIDE in kinds


def test_an_unticked_repair_and_an_overruled_flip_are_applied_as_one_step(opened, qtbot) -> None:
    table = opened.workspace.table
    panel = surveyed(opened, qtbot)
    missing = next(i for i in panel.survey.issues if i.kind == inspection.MISSING)
    assert panel.set_fixed(missing.id, False)
    amide = next(c for c in panel.survey.choices if c.kind == sites.AMIDE)
    other = next(s.label for k, s in enumerate(amide.states) if k != amide.chosen)
    assert panel.set_state(amide.id, other)
    assert panel.decisions() == ((missing.id,), {amide.id: other})
    assert "2 by hand" in panel.title()
    # what was decided is kept with the document before it is applied
    saved = tool.protein_workspace(ToolContext(opened)).review("1CRN")
    assert saved is not None and saved.skipped == (missing.id,)
    assert saved.overrides == {amide.id: other}

    ctx = ToolContext(opened)
    job = inspect_dialog.apply_review(ctx, panel, *panel.decisions())
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)
    s = structure(opened)
    assert (s.element == HYDROGEN).any()
    r = next(i for i in range(s.n_residues) if s.residue_label(i) == "A/ARG 10")
    span = s.atoms_of(r)
    assert "CZ" not in {str(s.name[i]) for i in range(span.start, span.stop)}  # left alone
    assert table.history.undo_label == "Prepare 1CRN (reviewed)"
    # the review becomes the report, led by what was done by hand
    report = ctx.find_panel(PrepareReportPanel.tool_id)
    assert isinstance(report, PrepareReportPanel)
    texts = [row.text for row in report.report.rows]
    assert texts[0].startswith("Left as in the file, by hand")
    assert any(t.startswith(f"Chosen by hand: Amide flip {other}") for t in texts)
    assert tool.protein_workspace(ctx).review("1CRN") is None
    qtbot.waitUntil(lambda: ctx.find_panel(INSPECT_ID) is None, timeout=5_000)  # its work done
    table.undo()  # one step takes it all back
    assert not (structure(opened).element == HYDROGEN).any()


def test_a_resurvey_keeps_what_was_decided(opened, qtbot) -> None:
    panel = surveyed(opened, qtbot)
    amide = next(c for c in panel.survey.choices if c.kind == sites.AMIDE)
    other = next(s.label for k, s in enumerate(amide.states) if k != amide.chosen)
    panel.set_state(amide.id, other)
    panel.resurvey()
    qtbot.waitUntil(lambda: panel.apply_button.isEnabled(), timeout=120_000)
    held = panel.survey.choice(amide.id)
    assert held is not None and held.state.label == other  # the solve held it
    # and it is still the chemist's decision, recorded as theirs
    assert panel.decisions()[1] == {amide.id: other}
    assert panel.line(amide.id).font(2).bold()


def test_only_warnings_hides_the_rest(opened, qtbot) -> None:
    panel = surveyed(opened, qtbot)
    panel.only_warnings.setChecked(True)
    shown = [line for line in panel.lines() if not line.isHidden()]
    assert all(line.text(0) == "⚠" for line in shown)
    panel.only_warnings.setChecked(False)
    assert not any(line.isHidden() for line in panel.lines())


def test_reset_puts_every_line_back(opened, qtbot) -> None:
    panel = surveyed(opened, qtbot)
    missing = next(i for i in panel.survey.issues if i.kind == inspection.MISSING)
    panel.set_fixed(missing.id, False)
    panel.reset()
    assert panel.decisions() == ((), {})


def test_a_half_finished_review_comes_back_with_the_document(opened, qtbot, tmp_path) -> None:
    panel = surveyed(opened, qtbot)
    missing = next(i for i in panel.survey.issues if i.kind == inspection.MISSING)
    panel.set_fixed(missing.id, False)
    path = tmp_path / "reviewing.csd"
    assert opened._save_to(path)
    assert INSPECT_ID in [entry.get("tool") for entry in opened.panel_states()]
    opened.open_path(path)
    qtbot.waitUntil(lambda: opened.workspace.path == path, timeout=10_000)
    ctx = ToolContext(opened)

    def back() -> bool:
        found = ctx.find_panel(INSPECT_ID)
        return isinstance(found, InspectorPanel) and found.apply_button.isEnabled()

    qtbot.waitUntil(back, timeout=120_000)
    again = ctx.find_panel(INSPECT_ID)
    assert again.decisions() == ((missing.id,), {})


def _overruled_and_resurveyed(window, qtbot):
    panel = surveyed(window, qtbot)
    amide = next(c for c in panel.survey.choices if c.kind == sites.AMIDE)
    proposed = amide.state.label
    other = next(s.label for k, s in enumerate(amide.states) if k != amide.chosen)
    panel.set_state(amide.id, other)
    panel.resurvey()
    qtbot.waitUntil(lambda: panel.apply_button.isEnabled(), timeout=120_000)
    assert panel.survey.choice(amide.id).state.label == other  # shown held
    return panel, amide, proposed, other


def test_after_a_resurvey_the_report_still_says_what_was_proposed(opened, qtbot) -> None:
    """The survey on screen already holds the choice; the report's "proposed" must come from
    the survey with nothing held, or it records the chemist's choice as the program's."""
    panel, amide, proposed, other = _overruled_and_resurveyed(opened, qtbot)
    assert proposed != other
    ctx = ToolContext(opened)
    job = inspect_dialog.apply_review(ctx, panel, *panel.decisions())
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)
    report = ctx.find_panel(PrepareReportPanel.tool_id)
    lines = [row.text for row in report.report.rows if row.text.startswith("Chosen by hand")]
    assert lines == [f"Chosen by hand: Amide flip {other} (proposed: {proposed})"]


def test_reset_after_a_resurvey_applies_exactly_what_prepare_would(opened, qtbot) -> None:
    """Reset must drop a held choice for good, not leave it riding on the shown survey."""
    before = structure(opened)
    panel, _amide, _proposed, _other = _overruled_and_resurveyed(opened, qtbot)
    panel.reset()
    assert panel.decisions() == ((), {})
    ctx = ToolContext(opened)
    job = inspect_dialog.apply_review(ctx, panel, *panel.decisions())
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)
    plain = pipeline.prepare(before, pipeline.Options())
    applied = structure(opened)
    assert applied.n_atoms == plain.structure.n_atoms
    assert np.array_equal(applied.uid, plain.structure.uid)
    assert np.allclose(applied.coords, plain.structure.coords)
    report = ctx.find_panel(PrepareReportPanel.tool_id)
    assert not any("by hand" in row.text for row in report.report.rows)
