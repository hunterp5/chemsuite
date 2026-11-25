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
"""Tools ▸ Prepare…: the dialog, the one undoable edit, and the report panel."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemsuite.app.document.tools import ToolContext
from chemsuite.macromolecule.prepare import pipeline, sites
from chemsuite.macromolecule.structure import edit, store
from chemsuite.macromolecule.ui.panel import MacromoleculePanel
from chemsuite.macromolecule.ui.report import TOOL_ID as REPORT_ID
from chemsuite.macromolecule.ui.report import PrepareReportPanel
from chemsuite.tools.macromolecule import prepare_dialog, prepare_options, tool
from tests.gui import gui

pytestmark = gui
DATA = Path(__file__).resolve().parents[2] / "data" / "protein"
HYDROGEN = 1


@pytest.fixture
def opened(window, qtbot, monkeypatch):
    """Crambin opened from a file (small, and prepared in a second)."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    job = tool.open_protein(ToolContext(window), [str(DATA / "1crn.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    return window


def hydrogens(window, protein: str = "1CRN") -> int:
    structure = store.as_structure(window.workspace.table.assets[protein].payload)
    return int((structure.element == HYDROGEN).sum())


def prepared(window, qtbot, options: pipeline.Options | None = None) -> PrepareReportPanel:
    """Prepare the protein and wait for the report panel."""
    ctx = ToolContext(window)
    job = prepare_dialog.prepare_protein(ctx, "1CRN", options or pipeline.Options())
    qtbot.waitUntil(lambda: not job.running, timeout=120_000)
    panel = ctx.find_panel(REPORT_ID)
    assert isinstance(panel, PrepareReportPanel)
    return panel


def test_prepare_is_in_the_viewers_tools_menu(opened) -> None:
    assert PrepareReportPanel.tool_id == tool.PREPARE_ID  # the tool reopens the report panel
    panel = ToolContext(opened).find_panel("protein")
    menu = next(a.menu() for a in panel.menu_bar.actions() if a.text() == "Tools")
    menu.aboutToShow.emit()
    action = next(a for a in menu.actions() if a.text() == "Prepare…")
    assert action.isEnabled()
    action.trigger()
    # the dialog belongs to the viewer's window, which it was opened from
    dialog = next(d for d in panel.window().findChildren(prepare_dialog.PrepareDialog))
    assert dialog.protein.currentData() == "1CRN"
    dialog.close()


def test_the_dialog_needs_a_protein(window) -> None:
    dialog = prepare_dialog.PrepareDialog(ToolContext(window))
    assert dialog.problem() == "There are no proteins in this document"
    assert not dialog.primary.isEnabled() and dialog.run() is None


def test_the_dialog_offers_the_shown_protein_and_its_defaults(opened) -> None:
    dialog = prepare_dialog.PrepareDialog(ToolContext(opened))
    assert dialog.protein.currentData() == "1CRN" and dialog.problem() == ""
    assert dialog.verb() == "Prepare 1CRN"
    options = dialog.options()
    assert options.ph == 7.0 and options.waters == pipeline.ALL
    assert options.repair and options.optimise
    assert options.titrate == sites.OFF  # off by default: less accurate on average (P2)
    assert options == pipeline.Options()  # every finer setting at the preparation's default


def test_every_setting_goes_through_the_dialog_and_the_document(opened) -> None:
    changed = pipeline.Options(
        ph=6.5, waters=pipeline.NEAR_LIGANDS, repair=False, optimise=False,
        titrate=sites.POCKET, ligand_states=False,
        orient_waters=False, side_chains=False, fill_gaps=False, loop_max=4, caps=False,
        altloc=pipeline.FIRST, water_reach=7.5, water_b_max=40.0, his="HIE", flips=False,
        ambiguous=2.0,
    )  # fmt: skip
    assert set(changed.changed()) == set(prepare_options.NAMES)  # not one left out
    dialog = prepare_dialog.PrepareDialog(ToolContext(opened))
    dialog.fields.set_options(changed)
    assert dialog.options() == changed
    saved = prepare_options.settings_of(changed)
    assert prepare_options.options_from(saved) == changed
    assert prepare_options.options_from({}) == pipeline.Options()  # an older document
    junk = {"his": "HXX", "water_reach": "far", "caps": "yes", "water_b_max": None}
    assert prepare_options.options_from(junk) == pipeline.Options()


def test_settings_under_another_follow_it(opened) -> None:
    fields = prepare_dialog.PrepareDialog(ToolContext(opened)).fields
    assert not fields.water_reach.isEnabled()  # Waters ▸ All
    fields.waters.set_value(pipeline.NEAR_LIGANDS)
    assert fields.water_reach.isEnabled()
    fields.repair.setChecked(False)
    assert not fields.caps.isEnabled() and not fields.side_chains.isEnabled()
    fields.optimise.setChecked(False)
    assert not fields.his.isEnabled() and not fields.flips.isEnabled()


def test_the_report_says_which_settings_were_changed(opened, qtbot) -> None:
    panel = prepared(opened, qtbot, pipeline.Options(caps=False, his="HID"))
    first = panel.report.rows[0]
    assert first.label == prepare_dialog.SETTINGS
    assert first.text == "Cap Ends off, His HID"
    default = prepared(opened, qtbot)
    assert all(row.label != prepare_dialog.SETTINGS for row in default.report.rows)


def test_preparing_adds_hydrogens_in_one_undoable_step(opened, qtbot) -> None:
    table = opened.workspace.table
    assert hydrogens(opened) == 0
    panel = prepared(opened, qtbot)
    assert hydrogens(opened) > 300  # 46 residues' worth
    assert panel.report.protein == "1CRN" and panel.report.rows
    assert table.history.undo_label == "Prepare 1CRN"
    table.undo()
    assert hydrogens(opened) == 0 and list(table.assets) == ["1CRN"]


def test_preparing_keeps_its_record_and_a_second_run_keeps_what_the_first_built(
    opened, qtbot
) -> None:
    ctx = ToolContext(opened)
    workspace = tool.macromolecule_workspace(ctx)
    s = prepare_dialog.structure_of(ctx, "1CRN")
    tip = (s.name == "CZ") & (s.res_name[s.residue_of] == "ARG")  # Arg 10's, to be rebuilt
    workspace.session().edit("Cut", {"1CRN": edit.delete(s, tip)})
    prepared(opened, qtbot, pipeline.Options(ph=6.0))
    first = workspace.prepared("1CRN")
    s = prepare_dialog.structure_of(ctx, "1CRN")
    assert first is not None and first.matches(s) and first.ph == 6.0
    assert first.settings["ph"] == 6.0 and first.built
    prepared(opened, qtbot)
    again = workspace.prepared("1CRN")
    assert again is not None and first.built <= again.built  # found already there: still built


def test_the_report_docks_in_the_viewers_window(opened, qtbot) -> None:
    """A result worked out in the viewer shows there, not beside the table in another window."""
    from chemsuite.app.document.panel_window import PanelWindow

    viewer = ToolContext(opened).find_panel(tool.TOOL_ID)
    assert isinstance(viewer.window(), PanelWindow)
    panel = prepared(opened, qtbot)
    assert panel.window() is viewer.window()
    assert panel not in [dock.widget() for dock in opened.tool_panels.docks]
    # and it is saved with the document, as a panel beside the table would be
    assert REPORT_ID in [entry.get("tool") for entry in opened.panel_states()]


def test_the_report_docks_beside_the_table_with_the_viewer_closed(opened, qtbot) -> None:
    viewer = ToolContext(opened).find_panel(tool.TOOL_ID)
    viewer.window().close()
    qtbot.waitUntil(lambda: ToolContext(opened).find_panel(tool.TOOL_ID) is None, timeout=5_000)
    panel = prepared(opened, qtbot)
    assert panel.window() is opened  # never out of reach


def test_the_report_lines_select_their_residue_in_the_workbench(opened, qtbot) -> None:
    panel = prepared(opened, qtbot)
    session = tool.macromolecule_workspace(ToolContext(opened)).session()
    rows = [i for i in range(panel.rows.topLevelItemCount())
            if panel.report.rows[i].residue >= 0]  # fmt: skip
    assert rows, "a report of nothing in particular"
    item = panel.rows.topLevelItem(rows[0])
    panel.rows.setCurrentItem(item)
    chosen = panel.report.rows[rows[0]].residue
    shown = session.display("1CRN")
    assert session.selected_count() == int((shown.structure.residue_of == chosen).sum())
    # a line about two residues (a disulfide, a chain break) names both and selects the first
    assert item.text(1).startswith(shown.structure.residue_label(chosen))


@pytest.fixture
def framed(opened, monkeypatch) -> list[np.ndarray]:
    """What the viewer is asked to zoom onto (the view itself cannot draw offscreen; its
    camera is tested in tests/macromolecule/ui)."""
    calls: list[np.ndarray] = []
    monkeypatch.setattr(MacromoleculePanel, "frame_atoms", lambda self, xyz: calls.append(xyz))
    return calls


def test_a_report_line_zooms_to_its_residue(opened, qtbot, framed) -> None:
    panel = prepared(opened, qtbot)
    session = tool.macromolecule_workspace(ToolContext(opened)).session()
    rows = [i for i, row in enumerate(panel.report.rows) if row.residue >= 0]
    panel.rows.setCurrentItem(panel.rows.topLevelItem(rows[-1]))
    shown = session.display("1CRN")
    atoms = shown.xyz[shown.structure.residue_of == panel.report.rows[rows[-1]].residue]
    assert len(framed) == 1 and np.array_equal(framed[0], atoms)


def test_a_line_naming_two_residues_frames_both(opened, framed) -> None:
    ctx = ToolContext(opened)
    session = tool.macromolecule_workspace(ctx).session()
    s = session.display("1CRN").structure
    first, last = 0, s.n_residues - 1
    while s.res_kind[last] != s.res_kind[first]:
        last -= 1
    line = f"{s.residue_label(first)}{prepare_dialog.SPAN}{s.residue_label(last)}"
    prepare_dialog.select(ctx, "1CRN", first, line)
    both = session.display("1CRN").xyz[(s.residue_of == first) | (s.residue_of == last)]
    assert len(framed) == 1 and np.array_equal(framed[0], both)
    assert session.selected_count() == int((s.residue_of == first).sum())  # the first only


def test_a_hidden_protein_is_selected_but_not_framed(opened, framed) -> None:
    ctx = ToolContext(opened)
    session = tool.macromolecule_workspace(ctx).session()
    session.set_visible("1CRN", False)
    prepare_dialog.select(ctx, "1CRN", 3)
    assert framed == []
    assert session.display("1CRN").selected.any()


def test_the_warned_lines_are_marked_and_counted(opened, qtbot) -> None:
    panel = prepared(opened, qtbot)
    warned = [i for i, row in enumerate(panel.report.rows) if row.warn]
    for i in warned:
        assert panel.rows.topLevelItem(i).text(0) == "⚠"
    assert f"net charge {panel.report.net_charge:+d}" in panel.title()
    assert f"pH {panel.report.ph:g}" in panel.title()
    assert panel.copy_text().startswith("Residue\tWhat was done")


def test_the_report_comes_back_with_the_document(opened, qtbot) -> None:
    panel = prepared(opened, qtbot)
    rows = list(panel.report.rows)
    workspace = tool.macromolecule_workspace(ToolContext(opened))
    assert workspace.capture()["reports"]["1CRN"] == panel.report.to_json()
    entry = next(e for e in opened.panel_states() if e["tool"] == tool.PREPARE_ID)
    assert entry["state"] == {"protein": "1CRN"}
    opened.tool_panels.close_all()
    assert ToolContext(opened).find_panel(REPORT_ID) is None
    opened.restore_panels([entry])
    again = ToolContext(opened).find_panel(REPORT_ID)
    assert isinstance(again, PrepareReportPanel) and list(again.report.rows) == rows


def test_preparing_again_fills_the_open_panel(opened, qtbot) -> None:
    panel = prepared(opened, qtbot)
    again = prepared(opened, qtbot, pipeline.Options(waters=pipeline.NONE))
    assert again is panel  # one report panel, not two
    assert len([e for e in opened.panel_states() if e["tool"] == tool.PREPARE_ID]) == 1


def test_the_longest_loop_is_off_at_one_follows_repair_and_is_kept_in_range(opened) -> None:
    fields = prepare_dialog.PrepareDialog(ToolContext(opened)).fields
    assert fields.loop_max.value() == pipeline.Options().loop_max == 8
    fields.loop_max.setValue(prepare_options.LOOP_OFF)
    assert fields.loop_max.text() == "Off" and fields.options().loop_max == 1
    fields.repair.setChecked(False)
    assert not fields.loop_max.isEnabled()
    fields.repair.setChecked(True)
    assert fields.loop_max.isEnabled()
    read = prepare_options.options_from
    assert read({"loop_max": 99}).loop_max == prepare_options.LOOP_MOST == 12
    assert read({"loop_max": 0}).loop_max == 1
    assert read({"loop_max": True}).loop_max == 8  # not understood: the default
    assert prepare_options.describe(pipeline.Options(loop_max=1)) == "Longest Loop off"
    assert prepare_options.describe(pipeline.Options(loop_max=4)) == "Longest Loop 4 residues"
