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
"""The 3D Viewer panel: what it draws for the selection, conformers, measuring, saving."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from .conftest import gui

pytestmark = gui

ROOT = Path(__file__).parents[2]


def _conformers(smiles: str, count: int) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=7)
    return mol


@pytest.fixture
def chem(window):
    """Rows A (butanol, 3 conformers), B (benzene, 1 conformer) and C (piperidine, 2D only)."""
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.mols import mol_to_blob
    from chemsuite.core.table import Table
    from chemsuite.core.table.columns import StructureColumn

    table = Table()
    smiles = ["CCCCO", "c1ccccc1", "C1CCNCC1"]
    table.add_rows({"ID": ["A", "B", "C"], "SMILES": smiles}, undoable=False)
    table.convert_column("SMILES", "structure")
    mols = [_conformers("CCCCO", 3), _conformers("c1ccccc1", 1)]
    column = StructureColumn.from_parsed([mol_to_blob(m) for m in mols] + [b""], [*smiles[:2], ""])
    table.add_column("SMILES 3D", column)
    table.history.mark_clean()
    window._install(Workspace(table))
    window.resize(1000, 700)
    window.show()
    return window


ATOMS = {1: 15, 2: 12, 3: 17}  # each row's atoms, with every hydrogen


def _current(chem, oid: int) -> None:
    """Move the table's cursor to *oid*, leaving the selection alone."""
    from PySide6.QtCore import QItemSelectionModel

    model = chem.workspace.model
    chem.view.selectionModel().setCurrentIndex(
        model.index(int(model.rows_of_oids([oid])[0]), 0),
        QItemSelectionModel.SelectionFlag.NoUpdate,
    )


def _open(chem, qtbot, oid: int = 1):
    """The viewer opened with *oid* the table's current row."""
    _current(chem, oid)
    panel = chem.open_tool("molview")
    _shows(panel, qtbot, oid)
    return panel


def _shows(panel, qtbot, oid: int) -> None:
    qtbot.waitUntil(
        lambda: panel.oid == oid and len(panel.renderer.scene) == ATOMS[oid], timeout=30_000
    )


def _current_oid(chem) -> int:
    return chem.workspace.model.oid_at(chem.view.currentIndex().row())


def test_opens_on_the_current_row_in_the_3d_column(chem, qtbot) -> None:
    panel = _open(chem, qtbot, 2)
    assert panel in chem.panels
    assert panel.column_name() == "SMILES 3D"
    assert panel.row_counter.text() == "Row 2 of 3"


def test_browsing_moves_the_table_and_follows_it(chem, qtbot) -> None:
    panel = _open(chem, qtbot, 1)
    panel.go_by(1)
    assert panel.oid == 2 and _current_oid(chem) == 2
    assert chem.workspace.selection.oids == set()  # browsing never selects
    assert panel.previous_row.isEnabled() and panel.next_row.isEnabled()
    _current(chem, 1)  # moved in the table
    assert panel.oid == 1 and not panel.previous_row.isEnabled()


def test_browse_selected_and_visible(chem, qtbot) -> None:
    from chemsuite.core.table.filters import OidFilter

    table = chem.workspace.table
    panel = _open(chem, qtbot, 1)
    assert panel.browsed.tolist() == [1, 2, 3]  # Browse Visible, nothing hidden
    chem.workspace.selection.set([1, 2], "test")
    panel.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1, 2], timeout=5_000)
    table.set_filter("test", OidFilter([1, 3]))  # row 2 hidden
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1], timeout=5_000)
    panel.browse_visible.setChecked(False)  # the hidden selected row too
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1, 2], timeout=5_000)
    panel.browse_selected.setChecked(False)  # every row
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1, 2, 3], timeout=5_000)
    panel.go(1)
    _shows(panel, qtbot, 2)  # a hidden row can be viewed; the table's cursor stays
    assert _current_oid(chem) == 1
    # In the table's sort order.
    table.set_sort([("ID", False)])
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [3, 2, 1], timeout=5_000)
    assert panel.oid == 2 and panel.row_counter.text() == "Row 2 of 3"


def test_deselecting_the_row_shown_moves_on(chem, qtbot) -> None:
    from PySide6.QtCore import Qt

    chem.workspace.selection.set([1, 2], "test")
    panel = _open(chem, qtbot, 1)
    panel.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1, 2], timeout=5_000)
    panel._clicked(0, Qt.KeyboardModifier.ControlModifier)  # Ctrl-click: deselect row 1
    assert chem.workspace.selection.oids == {2}
    _shows(panel, qtbot, 2)
    assert panel.row_counter.text() == "Row 1 of 1"


def test_conformer_stepper_and_overlay(chem, qtbot) -> None:
    panel = _open(chem, qtbot, 1)
    first = panel.renderer.scene.positions.copy()
    assert panel.counter.text() == "1 / 3"
    assert not panel.previous.isEnabled() and panel.next.isEnabled()
    panel.step(1)
    assert panel.counter.text() == "2 / 3"
    assert (panel.renderer.scene.positions != first).any()
    panel.step(10)
    assert panel.counter.text() == "3 / 3" and not panel.next.isEnabled()
    panel.go_by(1)  # another row starts at its first conformer
    _shows(panel, qtbot, 2)
    assert panel.conformer == 0 and not panel.next.isVisible()  # one conformer: no stepper
    panel.go_by(-1)
    _shows(panel, qtbot, 1)
    panel.overlay_action.setChecked(True)
    qtbot.waitUntil(lambda: panel.counter.text() == "3 Conformers", timeout=5_000)
    # every conformer, hydrogens only on the oxygen: 5 heavy atoms + 1 H, three times
    assert len(panel.renderer.scene) == 18
    assert not panel.previous.isEnabled() and not panel.next.isEnabled()


def test_keys_step_rows_and_conformers(chem, qtbot) -> None:
    from PySide6.QtCore import Qt

    panel = _open(chem, qtbot, 1)
    qtbot.keyClick(panel, Qt.Key.Key_Down)
    assert panel.conformer == 1
    qtbot.keyClick(panel, Qt.Key.Key_End)
    assert panel.oid == 3
    qtbot.keyClick(panel, Qt.Key.Key_Left)
    assert panel.oid == 2
    qtbot.keyClick(panel, Qt.Key.Key_Home)
    assert panel.oid == 1 and panel.conformer == 0


def test_measuring_marks_atoms(chem, qtbot) -> None:
    from rdkit.Chem import rdMolTransforms

    panel = _open(chem, qtbot, 1)
    panel.set_measuring(True)
    panel._clicked(0, None)
    panel._clicked(1, None)
    distance = rdMolTransforms.GetBondLength(_conformers("CCCCO", 3).GetConformer(0), 0, 1)
    assert panel.measurement() == f"{distance:.2f} Å"
    panel._clicked(2, None)
    assert panel.measurement().endswith("°")
    panel.set_measuring(False)
    assert panel.measurement() == ""


def test_a_2d_row_is_made_3d_in_the_background(chem, qtbot) -> None:
    panel = _open(chem, qtbot, 1)
    panel.column.setCurrentText("SMILES")  # 2D only
    panel.go(2)
    _shows(panel, qtbot, 3)
    spread = panel.renderer.scene.positions.std(axis=0)
    assert spread.min() > 0.3  # a chair, not the flat 2D drawing
    _assert_fitted(panel)  # fitted when its 3D structure arrived, not while it was empty


def _assert_fitted(panel) -> None:
    renderer = panel.renderer
    scene, camera = renderer.scene, renderer.camera
    width, height = renderer.width(), renderer.height()
    xy, depth = camera.project(scene.positions, width, height)
    size = scene.radii * np.array([camera.pixels_per_unit(d, height) for d in depth])
    lo, hi = (xy - size[:, None]).min(axis=0), (xy + size[:, None]).max(axis=0)
    assert np.all(np.abs((lo + hi) / 2 - [width / 2, height / 2]) < 3)
    assert max((hi[0] - lo[0]) / width, (hi[1] - lo[1]) / height) > 0.9


def test_each_row_browsed_to_is_centred_and_fitted(chem, qtbot) -> None:
    panel = _open(chem, qtbot, 1)
    panel.renderer.camera.turn(120, 60)  # the user turned and zoomed
    panel.renderer.camera.zoom(3.0)
    panel.go_by(1)
    _shows(panel, qtbot, 2)
    _assert_fitted(panel)
    panel.step(0)
    panel.window().resize(1300, 700)  # wider: still fitted
    qtbot.wait(20)
    _assert_fitted(panel)


def test_state_round_trip(chem, qtbot) -> None:
    from chemsuite.app.molview.panel import BY_ROW
    from chemsuite.molview.geometry import SPACEFILL

    panel = _open(chem, qtbot, 1)
    panel.style_box.setCurrentIndex(panel.style_box.findData(SPACEFILL))
    panel.colour.setCurrentIndex(panel.colour.findData(BY_ROW))
    panel.browse_visible.setChecked(False)
    panel.step(2)
    panel.renderer.camera.turn(40, 10)
    panel.renderer.camera.zoom(1.5)
    state = json.loads(json.dumps(panel.state()))  # as saved in the document
    _current(chem, 3)
    again = chem.open_tool("molview")
    again.set_state(state)
    assert again.state() == state
    assert _current_oid(chem) == 1  # the saved row is the table's current row again


def test_state_survives_save_and_reopen(chem, qtbot, tmp_path) -> None:
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    panel = _open(chem, qtbot, 1)
    panel.step(1)
    panel.go_by(0)
    saved = panel.state()
    path = tmp_path / "viewer.csd"
    assert chem._save_to(path)
    chem._install(Workspace(Table()))
    chem.open_path(path)
    qtbot.waitUntil(lambda: chem.workspace.path == path, timeout=10_000)
    reopened = [p for p in chem.panels if getattr(p, "tool_id", None) == "molview"]
    assert len(reopened) == 1
    assert reopened[0].state()["conformer"] == 1
    assert reopened[0].state()["row"] == 1
    assert reopened[0].state()["camera"] == saved["camera"]


def test_shaders_are_built_from_the_current_sources() -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import build_shaders
    finally:
        sys.path.pop(0)
    manifest = json.loads(build_shaders.MANIFEST.read_text())
    assert manifest == build_shaders.hashes(), "shaders changed: run scripts/build_shaders.py"
    for source in build_shaders.sources():
        assert source.with_name(source.name + ".qsb").exists()


RENDER = """
import sys
from chemsuite.app.molview.renderer import MolRenderer, prepare_opengl
prepare_opengl()
from PySide6.QtWidgets import QApplication, QHBoxLayout, QWidget
from rdkit import Chem
from rdkit.Chem import AllChem
app = QApplication([])
import pyqtgraph.opengl as gl
from chemsuite.molview.geometry import SPACEFILL, combine, molecule_part
mol = Chem.AddHs(Chem.MolFromSmiles("O=C=O"))
AllChem.EmbedMolecule(mol, randomSeed=1)
window = QWidget()
layout = QHBoxLayout(window)
layout.addWidget(gl.GLViewWidget())  # a 3D plot in the same window, drawn first
window.resize(400, 200)
window.show()
app.processEvents()
view = MolRenderer()
failed = []
view.failed.connect(lambda: failed.append(1))
layout.addWidget(view)
view.set_scene(combine([molecule_part(mol, style=SPACEFILL, hydrogens="all")]), fit=True)
for _ in range(20):
    app.processEvents()
image = view.grabFramebuffer()
if image.isNull():
    print("NO-RHI")
    sys.exit(0)
print("FAILED" if failed else "OK")
xy, _ = view.camera.project(view.scene.positions, image.width(), image.height())
for (x, y), number in zip(xy, view.scene.atom_number):
    c = image.pixelColor(int(x), int(y))
    print(int(number), c.red(), c.green(), c.blue())
"""


@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_atoms_are_drawn_in_their_colours() -> None:
    """Draws with the GPU beside a 3D plot in one window, in a subprocess: the offscreen
    platform the other tests use has no GPU drawing."""
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER], capture_output=True, text=True, env=env, cwd=ROOT,
        timeout=120,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        pytest.skip("no GPU drawing on this machine")
    assert "OK" in out.stdout.split(), out.stdout + out.stderr
    rows = [line.split() for line in out.stdout.splitlines() if line[:1].isdigit()]
    assert len(rows) == 3
    for number, r, g, b in ((int(v) for v in row) for row in rows):
        if number == 8:  # oxygen: red
            assert r > 120 and r > 2 * g and r > 2 * b, (r, g, b)
        else:  # carbon: grey
            assert abs(r - g) < 30 and abs(g - b) < 30 and r > 40, (r, g, b)
