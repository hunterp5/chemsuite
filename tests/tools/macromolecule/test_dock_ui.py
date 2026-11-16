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
"""Tools ▸ Dock… and Rescore Poses…: the dialogs, the box, and the rows a run writes.

Vina itself is stood in for here (``tests/macromolecule/dock`` runs the real binary); these tests
are about the bridge: what the dialog asks for, and what lands in the table.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.app.document.tools import ToolContext
from chemsuite.macromolecule.dock import run as dock_run
from chemsuite.macromolecule.dock.vina import Pose
from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure import store
from chemsuite.tools.macromolecule import dock, dock_dialog, tool
from chemsuite.tools.macromolecule.rows import LIGAND_COLUMN, PROTEIN_COLUMN
from tests.gui import gui

pytestmark = gui
DATA = Path(__file__).resolve().parents[2] / "data" / "protein"


@pytest.fixture
def opened(window, qtbot, monkeypatch):
    """Crambin in the document, prepared, with one ligand row to dock."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    ctx = ToolContext(window)
    job = tool.open_protein(ctx, [str(DATA / "1crn.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    asset = window.workspace.table.assets["1CRN"]
    prepared = pipeline.prepare(store.as_structure(asset.payload), pipeline.Options())
    tool.macromolecule_workspace(ctx).session().edit("Prepare", {"1CRN": prepared.structure})
    ctx.add_rows({"Name": ["ethanol"], "Ligand": ["CCO"]}, "Add", structures=("Ligand",))
    return window


def _mol() -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=5)
    return mol


@pytest.fixture
def no_vina(monkeypatch):
    """Vina stood in for: two poses a ligand, scored by rank."""

    def fake(receptor, mol, box, options=None, cancel=None, cache=None):
        return [Pose(mol, 1, -8.5, 0.0, 0.0), Pose(mol, 2, -7.25, 1.5, 2.5)]

    monkeypatch.setattr(dock_run, "dock_one", fake)
    monkeypatch.setattr(dock_run, "score_only", lambda *a, **k: -6.75)


def test_the_dialog_asks_for_a_prepared_protein(window) -> None:
    dialog = dock_dialog.DockDialog(ToolContext(window))
    assert dialog.problem() == "There are no proteins in this document"
    dialog.close()


def test_an_unprepared_protein_is_sent_to_prepare_first(window, qtbot, monkeypatch) -> None:
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    ctx = ToolContext(window)
    job = tool.open_protein(ctx, [str(DATA / "1crn.cif.gz")])
    qtbot.waitUntil(lambda: not job.running, timeout=30_000)
    ctx.add_rows({"Ligand": ["CCO"]}, "Add", structures=("Ligand",))
    dialog = dock_dialog.DockDialog(ctx)
    assert "prepare" in dialog.problem().lower()
    dialog.close()


def test_the_box_comes_from_the_atoms_it_is_asked_for(opened) -> None:
    dialog = dock_dialog.DockDialog(ToolContext(opened))
    dialog.source.set_value(dock_dialog.WHOLE)
    dialog._update()
    whole = dialog.box()
    assert whole is not None and whole.volume > 0
    structure = store.as_structure(opened.workspace.table.assets["1CRN"].payload)
    assert whole.holds(structure.xyz())
    assert [field.value() for field in dialog.centre] == pytest.approx(
        [round(v, 2) for v in whole.centre], abs=0.01
    )
    # nothing is focused and nothing selected, so those sources have no box to give
    dialog.source.set_value(dock_dialog.FOCUSED)
    dialog._update()
    assert dialog.box() is None
    assert "box" in dialog.problem()
    # typed numbers are taken as they stand
    dialog.source.set_value(dock_dialog.TYPED)
    dialog._update()
    for field, value in zip(dialog.centre, (1.0, 2.0, 3.0), strict=True):
        field.setValue(value)
    for field in dialog.size:
        field.setValue(15.0)
    typed = dialog.box()
    assert typed is not None
    assert typed.centre == pytest.approx((1.0, 2.0, 3.0))
    assert typed.size == pytest.approx((15.0, 15.0, 15.0))
    dialog.close()


def test_every_pose_is_kept_in_the_viewer_and_copied_to_a_row_on_request(
    opened, qtbot, no_vina
) -> None:
    from chemsuite.macromolecule.ui.poses import PosesPanel

    ctx = ToolContext(opened)
    dialog = dock_dialog.DockDialog(ctx)
    dialog.structure.setCurrentText("Ligand")
    dialog.source.set_value(dock_dialog.WHOLE)
    dialog._update()
    assert dialog.problem() == ""
    assert dialog.verb() == "Dock 4 Rows"  # the three rows without a ligand are skipped
    job = dialog.run()
    qtbot.waitUntil(lambda: not job.running, timeout=60_000)
    dialog.close()
    table = opened.workspace.table
    panel = ctx.find_panel(PosesPanel.tool_id)
    assert isinstance(panel, PosesPanel)
    assert [pose.rank for pose in panel.run.poses] == [1, 2]
    assert [pose.affinity for pose in panel.run.poses] == [-8.5, -7.25]
    assert {pose.ligand for pose in panel.run.poses} == {"ethanol"}
    assert dock.POSE_COLUMN not in table  # nothing written until asked for
    panel.rows.selectAll()
    assert len(panel.send_to_table()) == 2
    poses = [
        oid for oid in table.oids.tolist() if str(table.cell(dock.POSE_COLUMN, oid) or "") != ""
    ]
    assert len(poses) == 2  # the two poses of the one ligand row
    assert [str(table.cell(dock.POSE_COLUMN, o)) for o in poses] == ["1", "2"]
    assert [float(table.cell(dock.AFFINITY_COLUMN, o)) for o in poses] == [-8.5, -7.25]
    assert all(str(table.cell(PROTEIN_COLUMN, o)) == "1CRN" for o in poses)
    assert all(str(table.cell(dock.PARENT_COLUMN, o)) == "ethanol" for o in poses)
    # the pose is a 3D structure, so the row follows into the pocket like any ligand row
    column = table.column(LIGAND_COLUMN)
    mol = column.obj(int(table.slots_of(np.array([poses[0]]))[0]))
    assert mol is not None and mol.GetNumConformers() == 1


def test_rescoring_keeps_one_pose_a_row_in_the_viewer(opened, qtbot, no_vina) -> None:
    """A rescore is work done in the pocket, so its scores stay in the viewer's panel and
    the table is left alone until they are asked for."""
    from chemsuite.macromolecule.ui.poses import PosesPanel

    ctx = ToolContext(opened)
    ctx.add_rows(
        {"Name": ["pose"], "Ligand": [Chem.MolToMolBlock(_mol())]},
        "Add pose",
        structures=("Ligand",),
    )
    table = opened.workspace.table
    before = len(table.oids)
    dialog = dock_dialog.RescoreDialog(ctx)
    assert dialog.problem() == ""
    job = dialog.run()
    qtbot.waitUntil(lambda: not job.running, timeout=60_000)
    dialog.close()
    panel = ctx.find_panel(PosesPanel.tool_id)
    assert isinstance(panel, PosesPanel)
    run_poses = panel.run
    assert len(run_poses.poses) == 2  # the 2D row got a conformer, the pose taken as it stands
    assert {pose.affinity for pose in run_poses.poses} == {-6.75}
    assert dock.AFFINITY_COLUMN not in table and len(table.oids) == before
    # Copy to Table is what puts them there
    panel.rows.selectAll()
    assert len(panel.send_to_table()) == 2
    scored = [
        oid for oid in table.oids.tolist() if str(table.cell(dock.AFFINITY_COLUMN, oid) or "") != ""
    ]
    assert len(scored) == 2
    assert {float(table.cell(dock.AFFINITY_COLUMN, o)) for o in scored} == {-6.75}


def test_a_ligand_with_no_3d_gets_a_conformer() -> None:
    flat = Chem.MolFromSmiles("CCO")
    assert flat.GetNumConformers() == 0
    grown = dock.with_coordinates(flat)
    assert grown is not None and grown.GetConformer().Is3D()
    assert dock.with_coordinates(None) is None


def test_the_poses_of_a_run_name_their_source() -> None:
    mol = _mol()
    results = [
        dock_run.Docked(0, (Pose(mol, 1, -9.0, 0.0, 0.0),)),
        dock_run.Docked(1, error="Meeko could not prepare this ligand"),
    ]
    run_poses = dock.docked_poses("1CRN", results, ["first", "second"])
    assert run_poses.protein == "1CRN" and run_poses.failed == 1
    assert [pose.ligand for pose in run_poses.poses] == ["first"]  # the failed one keeps none
    assert run_poses.best == -9.0
    assert "V2000" in run_poses.poses[0].molblock or "V3000" in run_poses.poses[0].molblock
    rows = dock.pose_rows("1CRN", run_poses.poses)
    assert rows[PROTEIN_COLUMN] == ["1CRN"]
    assert rows[dock.PARENT_COLUMN] == ["first"]
    assert rows[dock.AFFINITY_COLUMN] == ["-9.00"]
    assert "V2000" in rows[LIGAND_COLUMN][0] or "V3000" in rows[LIGAND_COLUMN][0]


def _dock(window, qtbot):
    """Dock the ligand row and return the Poses panel."""
    from chemsuite.macromolecule.ui.poses import PosesPanel

    ctx = ToolContext(window)
    dialog = dock_dialog.DockDialog(ctx)
    dialog.structure.setCurrentText("Ligand")
    dialog.source.set_value(dock_dialog.WHOLE)
    dialog._update()
    job = dialog.run()
    qtbot.waitUntil(lambda: not job.running, timeout=60_000)
    dialog.close()
    panel = ctx.find_panel(PosesPanel.tool_id)
    assert isinstance(panel, PosesPanel)
    return panel


def test_the_poses_panel_docks_in_the_viewers_window(opened, qtbot, no_vina) -> None:
    from chemsuite.app.document.panel_window import PanelWindow

    viewer = tool.open_viewer(ToolContext(opened))
    panel = _dock(opened, qtbot)
    assert isinstance(viewer.window(), PanelWindow)
    assert panel.window() is viewer.window()  # beside the pocket, not beside the table
    assert panel not in [dock.widget() for dock in opened.tool_panels.docks]


def test_choosing_a_pose_shows_it_in_the_pocket(opened, qtbot, no_vina) -> None:
    tool.open_viewer(ToolContext(opened))
    panel = _dock(opened, qtbot)
    session = tool.macromolecule_workspace(ToolContext(opened)).session()
    assert session.pose is None
    top = panel.rows.topLevelItem(0)
    panel.rows.setCurrentItem(top.child(0))
    assert len(panel.chosen()) == 1 and panel.chosen()[0].rank == 1
    # the pose is drawn, in 1CRN's pocket
    assert session.pose is not None and session.pose.name == "1CRN"


def test_the_poses_come_back_with_the_document(opened, qtbot, no_vina, tmp_path) -> None:
    from chemsuite.macromolecule.ui.poses import PosesPanel

    panel = _dock(opened, qtbot)
    kept = panel.run
    path = tmp_path / "docked.csd"
    assert opened._save_to(path)
    assert PosesPanel.tool_id in [entry.get("tool") for entry in opened.panel_states()]
    opened.open_path(path)
    qtbot.waitUntil(lambda: opened.workspace.path == path, timeout=10_000)
    qtbot.wait(20)
    ctx = ToolContext(opened)
    assert tool.macromolecule_workspace(ctx).docked("1CRN") == kept
    again = ctx.find_panel(PosesPanel.tool_id)
    assert isinstance(again, PosesPanel) and again.run == kept and kept.number == 1


def test_rescore_offers_the_docked_poses_first(opened, qtbot, no_vina) -> None:
    """A protein that has been docked into is rescored from its own poses, not the table."""
    ctx = ToolContext(opened)
    dialog = dock_dialog.RescoreDialog(ctx)
    assert dialog.source.value() == dock_dialog.TABLE  # nothing docked yet
    assert dialog.problem() == "" and dialog.verb().endswith("Rows")
    dialog.close()

    panel = _dock(opened, qtbot)
    dialog = dock_dialog.RescoreDialog(ctx)
    assert dialog.source.value() == dock_dialog.KEPT
    assert dialog.verb() == "Rescore 2 Poses"
    job = dialog.run()
    qtbot.waitUntil(lambda: not job.running, timeout=60_000)
    dialog.close()
    # the rescore is a run of its own, beside the dock it scored; the table is untouched
    assert {pose.affinity for pose in panel.run.poses} == {-6.75}
    assert panel.run.kind == "rescored" and panel.run.scoring == "vina"
    assert panel.run.label == "Run 2 · Rescored (vina)"
    runs = tool.macromolecule_workspace(ctx).runs("1CRN")
    assert [run.kind for run in runs] == ["docked", "rescored"]
    assert [pose.affinity for pose in runs[0].poses] == [-8.5, -7.25]  # the dock kept as it was
    assert panel.choice.count() == 2 and dock.AFFINITY_COLUMN not in opened.workspace.table


def test_rescore_scores_the_run_shown_in_the_poses_panel(opened, qtbot, no_vina) -> None:
    ctx = ToolContext(opened)
    panel = _dock(opened, qtbot)
    first = panel.run
    _dock(opened, qtbot)  # a second run, now shown
    panel.choice.setCurrentIndex(panel.choice.findData(first.number))
    assert panel.run == first
    dialog = dock_dialog.RescoreDialog(ctx)
    assert dialog.run_choice.currentData() == first.number and dialog.run_choice.count() == 2
    dialog.close()


def test_a_run_can_be_removed_and_the_last_one_closes_the_panel(opened, qtbot, no_vina) -> None:
    from chemsuite.macromolecule.ui.poses import PosesPanel

    ctx = ToolContext(opened)
    panel = _dock(opened, qtbot)
    _dock(opened, qtbot)
    workspace = tool.macromolecule_workspace(ctx)
    assert [run.number for run in workspace.runs("1CRN")] == [1, 2]
    assert panel.run.number == 2 and "2 runs" in panel.title()
    assert panel.remove_run().number == 2
    assert [run.number for run in workspace.runs("1CRN")] == [1]
    assert panel.run.number == 1 and panel.choice.count() == 1
    panel.remove_run()
    assert workspace.runs("1CRN") == ()
    qtbot.waitUntil(lambda: ctx.find_panel(PosesPanel.tool_id) is None, timeout=5_000)


def test_rescore_says_so_when_a_protein_has_no_docked_poses(opened) -> None:
    dialog = dock_dialog.RescoreDialog(ToolContext(opened))
    dialog.source.set_value(dock_dialog.KEPT)
    assert dialog.problem() == (
        "This protein has no docked poses: dock into it first (Tools ▸ Dock…)"
    )
    dialog.close()


def test_a_later_run_into_another_protein_copies_with_that_protein(opened, qtbot) -> None:
    """The panel is filled again by a run into another protein, and what it copies to the
    table must say which protein the poses are in now, not which it was first opened for."""
    from chemsuite.macromolecule.workspace.poses import DockedPoses
    from chemsuite.macromolecule.workspace.poses import Pose as KeptPose

    ctx = ToolContext(opened)
    block = Chem.MolToMolBlock(_mol())
    first = DockedPoses("1CRN", (KeptPose("a", 1, -7.0, 0.0, 0.0, block),))
    second = DockedPoses("2XYZ", (KeptPose("b", 1, -8.0, 0.0, 0.0, block),))
    panel = dock.show_poses(ctx, first)
    assert dock.show_poses(ctx, second) is panel  # the same panel, filled again
    panel.rows.selectAll()
    assert len(panel.send_to_table()) == 1
    table = opened.workspace.table
    copied = [oid for oid in table.oids.tolist()
              if str(table.cell(dock.AFFINITY_COLUMN, oid) or "") == "-8.00"]  # fmt: skip
    assert [str(table.cell(PROTEIN_COLUMN, oid)) for oid in copied] == ["2XYZ"]
