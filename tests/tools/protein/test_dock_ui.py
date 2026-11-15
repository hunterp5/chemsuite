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
"""Protein ▸ Dock… and Rescore Poses…: the dialogs, the box, and the rows a run writes.

Vina itself is stood in for here (``tests/protein/dock`` runs the real binary); these tests
are about the bridge: what the dialog asks for, and what lands in the table.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.app.document.tools import ToolContext
from chemsuite.protein.dock import run as dock_run
from chemsuite.protein.dock.vina import Pose
from chemsuite.protein.prepare import pipeline
from chemsuite.protein.structure import store
from chemsuite.tools.protein import dock, dock_dialog, tool
from chemsuite.tools.protein.rows import LIGAND_COLUMN, PROTEIN_COLUMN
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
    tool.protein_workspace(ctx).session().edit("Prepare", {"1CRN": prepared.structure})
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


def test_every_pose_becomes_a_row_linked_to_its_protein(opened, qtbot, no_vina) -> None:
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


def test_rescoring_writes_one_pose_a_row(opened, qtbot, no_vina) -> None:
    ctx = ToolContext(opened)
    ctx.add_rows(
        {"Name": ["pose"], "Ligand": [Chem.MolToMolBlock(_mol())]},
        "Add pose",
        structures=("Ligand",),
    )
    dialog = dock_dialog.RescoreDialog(ctx)
    assert dialog.problem() == ""
    job = dialog.run()
    qtbot.waitUntil(lambda: not job.running, timeout=60_000)
    dialog.close()
    table = opened.workspace.table
    scored = [
        oid for oid in table.oids.tolist() if str(table.cell(dock.AFFINITY_COLUMN, oid) or "") != ""
    ]
    assert len(scored) == 2  # the 2D row got a conformer, the pose was taken as it stands
    assert {float(table.cell(dock.AFFINITY_COLUMN, o)) for o in scored} == {-6.75}


def test_a_ligand_with_no_3d_gets_a_conformer() -> None:
    flat = Chem.MolFromSmiles("CCO")
    assert flat.GetNumConformers() == 0
    grown = dock.with_coordinates(flat)
    assert grown is not None and grown.GetConformer().Is3D()
    assert dock.with_coordinates(None) is None


def test_the_pose_rows_of_a_run_name_their_source() -> None:
    mol = _mol()
    results = [
        dock_run.Docked(0, (Pose(mol, 1, -9.0, 0.0, 0.0),)),
        dock_run.Docked(1, error="Meeko could not prepare this ligand"),
    ]
    rows = dock.pose_rows("1CRN", results, ["first", "second"])
    assert rows[PROTEIN_COLUMN] == ["1CRN"]
    assert rows[dock.PARENT_COLUMN] == ["first"]  # the failed ligand writes no row
    assert rows[dock.AFFINITY_COLUMN] == ["-9.00"]
    assert "V2000" in rows[LIGAND_COLUMN][0] or "V3000" in rows[LIGAND_COLUMN][0]
