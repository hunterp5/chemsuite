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
from PySide6.QtCore import Qt
from rdkit import Chem
from rdkit.Chem import AllChem

from .conftest import gui

pytestmark = gui

ROOT = Path(__file__).parents[2]


@pytest.fixture(autouse=True)
def gpu(monkeypatch):
    """The offscreen platform has no GPU drawing, so the 3D view would report it and the
    viewer fall back to 2D; the tests act as if drawing works (the 3D scene is built all the
    same) and call the real handler when testing the fallback."""
    from chemsuite.app.viewer.panel import StructureViewer

    real = StructureViewer._failed
    monkeypatch.setattr(StructureViewer, "_failed", lambda self, *_: None)
    return real


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
    """The viewer opened in 3D (Chemistry ▸ Structure Viewer) with *oid* the table's current
    row."""
    _current(chem, oid)
    panel = chem.open_tool("viewer")
    assert panel.mode == "3d"
    _shows(panel, qtbot, oid)
    return panel


def _rows_text(panel) -> str:
    return f"{panel.position.value()} {panel.count_label.text()}"


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
    assert _rows_text(panel) == "2 of 3"


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
    before = chem.view.currentIndex()
    panel.go(1)
    _shows(panel, qtbot, 2)  # a hidden row can be viewed; the table's cursor stays
    assert chem.view.currentIndex() == before
    # In the table's sort order.
    table.set_sort([("ID", False)])
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [3, 2, 1], timeout=5_000)
    assert panel.oid == 2 and _rows_text(panel) == "2 of 3"


def test_deselecting_the_row_shown_moves_on(chem, qtbot) -> None:
    from PySide6.QtCore import Qt

    chem.workspace.selection.set([1, 2], "test")
    panel = _open(chem, qtbot, 1)
    panel.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: panel.browsed.tolist() == [1, 2], timeout=5_000)
    panel._clicked(0, Qt.KeyboardModifier.ControlModifier)  # Ctrl-click: deselect row 1
    assert chem.workspace.selection.oids == {2}
    _shows(panel, qtbot, 2)
    assert _rows_text(panel) == "1 of 1"


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
    from chemsuite.app.tools import ToolContext
    from chemsuite.app.viewer.panel import BY_ROW
    from chemsuite.molview.geometry import SPACEFILL
    from chemsuite.tools.viewer.tool import open_viewer

    panel = _open(chem, qtbot, 1)
    panel.style_box.setCurrentIndex(panel.style_box.findData(SPACEFILL))
    panel.colour.setCurrentIndex(panel.colour.findData(BY_ROW))
    panel.browse_visible.setChecked(False)
    panel.step(2)
    panel.renderer.camera.turn(40, 10)
    panel.renderer.camera.zoom(1.5)
    state = json.loads(json.dumps(panel.state()))  # as saved in the document
    _current(chem, 3)
    again = open_viewer(ToolContext(chem), "2d", state)  # a second viewer, as on reopening
    assert again is not panel
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
    reopened = [p for p in chem.panels if getattr(p, "tool_id", None) == "viewer"]
    assert len(reopened) == 1
    assert reopened[0].state()["conformer"] == 1
    assert reopened[0].state()["row"] == 1
    assert reopened[0].state()["camera"] == saved["camera"]
    assert reopened[0].mode == "3d"


# -- 2D (Data ▸ Browser) ------------------------------------------------------------------------


@pytest.fixture
def browser(window, qtbot):
    """Data ▸ Browser on rows A–D (D's SMILES is invalid), with a number column."""
    from chemsuite.app.main_window import Workspace
    from chemsuite.core.table import Table

    table = Table()
    table.add_rows(
        {
            "ID": list("ABCD"),
            "SMILES": ["CCO", "c1ccccc1", "CC(=O)O", "not a smiles"],
            "pIC50": [5.1, 6.25, 7.0, 4.5],
        },
        undoable=False,
    )
    table.convert_column("SMILES", "structure")
    window._install(Workspace(table))
    window.show()
    qtbot.waitExposed(window)
    return window.show_browser()


def ids(window, oids) -> list[str]:
    table = window.workspace.table
    return [table.cell("ID", int(o)) for o in oids]


def shown_id(viewer) -> str:
    return viewer.table.cell("ID", viewer.oid)


def test_browser_is_the_viewer_in_2d(window, browser) -> None:
    assert browser.mode == "2d" and browser in window.panels
    assert browser.renderer is None  # 2D never makes the 3D view (or touches the GPU)
    assert not browser.depiction.pixmap().isNull()
    browser.set_mode("3d")
    assert window.show_browser() is browser and browser.mode == "2d"  # brought forward


def test_a_cells_browse_opens_the_browser_at_its_row(window, browser) -> None:
    from .test_column_actions import trigger

    model = window.workspace.model
    browser.set_mode("3d")
    trigger(window.cell_menu(model.index(2, 0)), "Browse")
    assert shown_id(browser) == model.index(2, 0).data() and browser.mode == "3d"  # as it was
    window.panels[0].close()


def test_3d_made_after_the_viewer_closed_is_dropped(window, browser, monkeypatch) -> None:
    import shiboken6

    browser.set_mode("3d")
    key = (browser.column_name(), browser.oid)
    finish = []
    monkeypatch.setattr(browser, "_run", lambda _label, _work, done: finish.append(done))
    browser._made.pop(key, None)
    browser._make_3d(key[0], {key[1]: None})
    shiboken6.delete(browser)  # closed while its Make 3D structures job ran
    finish[0]({})


def test_a_cells_browse_opens_a_new_browser_in_2d(chem) -> None:
    from .test_column_actions import trigger

    model = chem.workspace.model
    trigger(chem.cell_menu(model.index(1, 0)), "Browse")
    viewer = chem.panels[-1]
    assert viewer.mode == "2d" and viewer.table.cell("ID", viewer.oid) == "B"


def test_browses_selected_rows_in_display_order(window, browser, qtbot) -> None:
    table = window.workspace.table
    window.workspace.selection.set(table.oids[[3, 0, 2]].tolist(), "test")
    table.set_sort([("ID", False)])
    browser.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: len(browser.browsed) == 3)
    assert ids(window, browser.browsed) == ["D", "C", "A"]
    browser.go(0)
    assert shown_id(browser) == "D"
    browser.go_by(1)
    assert shown_id(browser) == "C"
    # The table cursor follows, without touching the selection.
    assert window.workspace.model.oid_at(window.view.currentIndex().row()) == browser.oid
    assert len(window.workspace.selection) == 3


def test_keys_and_messages(window, browser, qtbot) -> None:
    qtbot.keyClick(browser, Qt.Key.Key_End)
    assert shown_id(browser) == "D"
    assert browser.depiction.text() == "No structure"  # the invalid SMILES
    qtbot.keyClick(browser, Qt.Key.Key_Up)  # rows, in 2D
    assert shown_id(browser) == "C"
    assert not browser.depiction.pixmap().isNull()
    qtbot.keyClick(browser, Qt.Key.Key_Space)
    assert browser.oid in window.workspace.selection and browser.select_box.isChecked()
    qtbot.keyClick(browser, Qt.Key.Key_Home)
    assert shown_id(browser) == "A" and not browser.previous_row.isEnabled()


def test_jump_box_and_follow_table_cursor(window, browser) -> None:
    browser.position.setValue(3)
    assert shown_id(browser) == "C"
    window.view.setCurrentIndex(window.workspace.model.index(1, 0))
    assert shown_id(browser) == "B"


def test_deselecting_while_browsing_selected_stays_in_place(window, browser, qtbot) -> None:
    window.workspace.selection.set(window.workspace.table.oids[:3].tolist(), "test")
    browser.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: len(browser.browsed) == 3)
    browser.go(1)
    browser.set_selected(False)  # B leaves the selected rows
    qtbot.waitUntil(lambda: len(browser.browsed) == 2)
    assert shown_id(browser) == "C"
    assert browser.count_label.text() == "of 2"


def test_empty_scope(window, browser, qtbot) -> None:
    browser.browse_selected.setChecked(True)
    qtbot.waitUntil(lambda: browser.oid is None)
    assert browser.depiction.text() == "No rows selected"
    assert not browser.select_box.isEnabled()


def test_caption_shows_the_chosen_columns(window, browser) -> None:
    table = window.workspace.table
    assert browser.caption_names == ["ID"]  # an ID-like column by default
    assert browser.caption.lines == [("ID", "A")]
    browser.set_caption(["ID", "pIC50"])
    browser.go(1)
    assert browser.caption.lines == [("ID", "B"), ("pIC50", table.cell("pIC50", browser.oid))]
    assert browser.caption_button.text() == "ID, pIC50"
    browser._fill_caption_menu()
    names = [a.text() for a in browser.caption_menu.actions() if a.text()]
    assert names == ["None", "ID", "pIC50"]  # no structure column
    browser._toggle_caption("ID", False)
    assert browser.caption_names == ["pIC50"]
    browser.set_caption([])
    assert browser.caption.lines == [] and browser.caption_button.text() == "None"


def test_the_drawing_fits_at_once_when_the_view_is_resized(browser, qtbot) -> None:
    depiction = browser.depiction
    depiction.resize(600, 500)
    depiction.show_mol(depiction.mol)
    depiction.resize(300, 250)  # smaller: fitted before the sharp redraw
    pixmap = depiction.pixmap()
    dpr = depiction.devicePixelRatioF()
    assert pixmap.width() <= 300 * dpr + 1 and pixmap.height() <= 250 * dpr + 1


def test_export_includes_the_caption(window, browser, tmp_path) -> None:
    from PySide6.QtGui import QImage

    browser.set_caption([])
    plain = browser.picture()
    browser.set_caption(["ID", "pIC50"])
    captioned = browser.picture()
    assert not plain.isNull() and plain.size() == captioned.size()
    assert plain != captioned
    assert browser.export(tmp_path / "view")
    assert not QImage(str(tmp_path / "view.png")).isNull()


# -- one viewer, two drawings ---------------------------------------------------------------------


def test_switching_keeps_the_row_column_rows_and_caption(chem, qtbot) -> None:
    _current(chem, 2)
    viewer = chem.show_browser()
    viewer.column.setCurrentText("SMILES 3D")
    viewer.browse_visible.setChecked(False)
    viewer.set_caption(["ID"])
    viewer.refresh()
    before = (viewer.oid, viewer.column_name(), viewer.browsed.tolist(), viewer.caption.lines)
    viewer.set_mode("3d")
    assert viewer.pages.currentIndex() == 1
    _shows(viewer, qtbot, 2)
    assert (
        viewer.oid,
        viewer.column_name(),
        viewer.browsed.tolist(),
        viewer.caption.lines,
    ) == before
    camera = viewer.renderer.camera.state()
    viewer.renderer.camera.turn(30, 10)  # the user turned it
    turned = viewer.renderer.camera.state()
    assert turned != camera
    viewer.set_mode("2d")
    assert viewer.pages.currentIndex() == 0
    assert (
        viewer.oid,
        viewer.column_name(),
        viewer.browsed.tolist(),
        viewer.caption.lines,
    ) == before
    viewer.set_mode("3d")
    assert viewer.renderer.camera.state() == turned  # back as it was, not refitted


def test_3d_controls_only_in_3d(chem) -> None:
    viewer = chem.show_browser()
    form = viewer.options.form
    three_d = [viewer.overlay_action, viewer.measure_action, viewer.reset_action]
    assert not any(a.isVisible() for a in three_d)
    assert not form.isRowVisible(viewer.style_box) and form.isRowVisible(viewer.column)
    assert form.isRowVisible(viewer.caption_button)
    viewer.set_mode("3d")
    assert all(a.isVisible() for a in three_d)
    assert form.isRowVisible(viewer.style_box)
    viewer.set_measuring(True)
    viewer.set_mode("2d")
    assert not viewer.measuring and not viewer.measure_action.isChecked()


def test_3d_is_made_ahead_for_the_rows_around(chem, qtbot) -> None:
    _current(chem, 1)
    viewer = chem.show_browser()
    viewer.column.setCurrentText("SMILES")  # 2D only
    viewer.refresh()
    viewer.make_ahead()
    assert viewer._ahead == {("SMILES", 1)}  # the row shown first, one row at a time
    viewer.go(1)
    viewer.make_ahead()  # nothing more is queued meanwhile
    assert viewer._waiting == {("SMILES", 1)}
    keys = {("SMILES", oid) for oid in (1, 2, 3)}  # then the rows either side, in turn
    qtbot.waitUntil(lambda: keys <= viewer._made.keys(), timeout=30_000)
    viewer.go(0)
    viewer.set_mode("3d")
    assert len(viewer.renderer.scene) > 0  # ready at once, no wait for a conformer


def test_without_3d_graphics_the_viewer_stays_2d(chem, qtbot, gpu) -> None:
    viewer = _open(chem, qtbot, 1)
    gpu(viewer)  # the real handler: the graphics cannot draw
    assert viewer.mode == "2d" and viewer.gpu_failed
    button = viewer.mode_buttons["3d"]
    assert not button.isEnabled() and button.toolTip().startswith("⚠")
    viewer.set_mode("3d")
    assert viewer.mode == "2d"
    assert not viewer.depiction.pixmap().isNull()
    viewer.make_ahead()
    assert not viewer._ahead  # no 3D made for a view that cannot show it


def test_a_saved_3d_viewer_reopens_in_3d(chem) -> None:
    chem._restore_panels([{"tool": "molview", "state": {"row": 2, "column": "SMILES 3D"}}])
    viewer = chem.panels[-1]
    assert viewer.tool_id == "viewer" and viewer.mode == "3d" and viewer.oid == 2
    assert chem.panel_states()[-1]["tool"] == "viewer"


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
import numpy as np
from PySide6.QtWidgets import QApplication, QHBoxLayout, QWidget
from rdkit import Chem
from rdkit.Chem import AllChem
app = QApplication([])
from chemsuite.app.molview.renderer import MolRenderer
from chemsuite.app.plots.view3d import View3D
from chemsuite.molview.geometry import SPACEFILL, combine, molecule_part
mol = Chem.AddHs(Chem.MolFromSmiles("O=C=O"))
AllChem.EmbedMolecule(mol, randomSeed=1)
window = QWidget()
layout = QHBoxLayout(window)
plot = View3D()  # a 3D plot in the same window: one graphics API for both
spots = np.array([[0.0, 0.0, 0.0], [0.6, -0.6, 0.0]])
plot.set_points(spots, np.array([[0.9, 0.1, 0.1, 1.0], [0.1, 0.2, 0.9, 1.0]]), 24.0)
plot.set_lines(np.array([[-1, -1, -1], [1, 1, 1]], float), (0.2, 0.2, 0.2, 1.0))
view = MolRenderer()
failed = []
view.failed.connect(lambda: failed.append(1))
plot.failed.connect(lambda: failed.append(1))
layout.addWidget(plot)
layout.addWidget(view)
window.resize(500, 250)
view.set_scene(combine([molecule_part(mol, style=SPACEFILL, hydrogens="all")]), fit=True)
try:
    window.show()
    for _ in range(20):
        app.processEvents()
except Exception as error:  # graphics too old to draw at all
    print("NO-RHI", type(error).__name__)
    sys.exit(0)
image = view.grabFramebuffer()
picture = plot.grabFramebuffer()
if image.isNull() or picture.isNull():
    print("NO-RHI")
    sys.exit(0)
print("FAILED" if failed else "OK", view.rhi().backendName())
xy, _ = view.camera.project(view.scene.positions, image.width(), image.height())
for (x, y), number in zip(xy, view.scene.atom_number):
    c = image.pixelColor(int(x), int(y))
    print(int(number), c.red(), c.green(), c.blue())
xy, _ = plot.camera.project(spots, picture.width(), picture.height())
for (x, y) in xy:
    c = picture.pixelColor(int(x), int(y))
    print("PLOT", c.red(), c.green(), c.blue())
"""


@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_molecules_and_3d_plots_draw_side_by_side() -> None:
    """Draws a molecule and a 3D plot with the GPU in one window, in a subprocess: the
    offscreen platform the other tests use has no GPU drawing."""
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER], capture_output=True, text=True, env=env, cwd=ROOT,
        timeout=120,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        if os.environ.get("CHEMSUITE_REQUIRE_GPU"):  # CI's Linux 3D step: must draw
            pytest.fail("no GPU drawing: " + out.stderr)
        pytest.skip("no GPU drawing on this machine")
    assert "OK" in out.stdout.split(), out.stdout + out.stderr
    rows = [line.split() for line in out.stdout.splitlines() if line[:1].isdigit()]
    assert len(rows) == 3
    for number, r, g, b in ((int(v) for v in row) for row in rows):
        if number == 8:  # oxygen: red
            assert r > 120 and r > 2 * g and r > 2 * b, (r, g, b)
        else:  # carbon: grey
            assert abs(r - g) < 30 and abs(g - b) < 30 and r > 40, (r, g, b)
    (red, blue) = [
        tuple(int(v) for v in line.split()[1:]) for line in out.stdout.splitlines()
        if line.startswith("PLOT")
    ]  # fmt: skip
    assert red[0] > 180 and red[1] < 80 and red[2] < 80, red  # the points, in their colours
    assert blue[2] > 180 and blue[0] < 80, blue


def test_3d_known_unavailable_at_start_up_is_never_tried(chem, qtbot, gpu, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.viewer.panel import StructureViewer

    monkeypatch.setattr(StructureViewer, "_failed", gpu)
    monkeypatch.setattr(graphics, "_problem", graphics.TURNED_OFF)
    _current(chem, 1)
    viewer = chem.open_tool("viewer")
    assert viewer.mode == "2d" and viewer.renderer is None
    assert viewer.mode_buttons["3d"].toolTip() == f"⚠ {graphics.TURNED_OFF}"
