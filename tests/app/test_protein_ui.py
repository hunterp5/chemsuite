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
"""The Protein Viewer panel: session, tree, sequence strip, selection bar, menus, state."""

from __future__ import annotations

from functools import cache
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QObject, QPoint, Qt, Signal

from chemsuite.core.table import Asset, PutAsset, Table
from chemsuite.protein import io, scene, store
from chemsuite.protein.model import Structure

DATA = Path(__file__).resolve().parents[1] / "data" / "protein"


@cache
def structure(name: str) -> Structure:
    return io.read(DATA / name)


def add_protein(table: Table, name: str, file: str, packed: bool = False) -> None:
    s = structure(file)
    payload = s.to_bytes() if packed else s
    table.edit(f"Open {name}", [PutAsset(name, None, Asset(name, "protein", name, payload, {}))])


@pytest.fixture
def table() -> Table:
    t = Table()
    add_protein(t, "1ATP", "1atp.cif.gz")
    return t


@pytest.fixture
def panel(qtbot, table, monkeypatch):
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    widget = ProteinPanel(table)
    qtbot.addWidget(widget)
    widget.resize(900, 600)
    return widget


# the session follows the document


def test_the_session_follows_the_documents_proteins(panel, table) -> None:
    session = panel.session
    assert [d.id for d in session.displays] == ["1ATP"]
    add_protein(table, "3HTB", "3htb.cif.gz", packed=True)  # stored bytes are decoded
    assert [d.id for d in session.displays] == ["1ATP", "3HTB"]
    assert panel.tree.topLevelItemCount() == 2
    table.undo()
    assert [d.id for d in session.displays] == ["1ATP"]
    table.redo()
    assert session.display("3HTB").structure.n_atoms == structure("3htb.cif.gz").n_atoms


def test_hiding_a_protein_takes_it_out_of_the_scene(panel, table) -> None:
    add_protein(table, "3HTB", "3htb.cif.gz")
    session = panel.session
    both = session.scene.atom_offsets[-1]
    item = panel.tree.topLevelItem(1)
    item.setCheckState(0, Qt.CheckState.Unchecked)
    assert not session.display("3HTB").visible
    assert session.scene.atom_offsets[-1] == both - structure("3htb.cif.gz").n_atoms
    assert len(panel.sequence.strip.rows) == 2  # 1ATP's chains E and I only


# selecting


def test_the_selection_bar_selects_and_marks_mistakes(panel) -> None:
    bar, session = panel.selection_bar, panel.session
    bar.field.setText("ligand")
    assert bar.apply()
    assert session.selected_count() == 31
    assert bar.count.text() == "31 atoms, 1 residues"
    bar.field.setText("protein and")
    assert not bar.apply()
    assert bar.field.toolTip().startswith("⚠")
    assert bar.field.cursorPosition() == len("protein and")
    assert session.selected_count() == 31  # a mistake changes nothing
    bar.field.setText("")
    bar.apply()
    assert session.selected_count() == 0


def test_saved_selections_survive_by_uid(panel) -> None:
    session = panel.session
    session.select_text("ligand")
    session.save_selection("atp")
    session.clear_selection()
    assert session.select_text("$atp") == 31
    with pytest.raises(Exception, match="not a name"):
        session.save_selection("two words")


def test_a_selection_only_changes_the_look(panel, qtbot) -> None:
    session = panel.session
    before = session.scene
    with qtbot.waitSignal(session.lookChanged), qtbot.assertNotEmitted(session.sceneChanged):
        session.select_text("chain I")
    after = session.scene
    assert after.geometry is before.geometry and after.cartoon is before.cartoon
    selected = after.sphere_look[:, 4] == scene.SELECTED
    chain_i = session.display("1ATP").structure.atom_chain[after.geometry.atoms] == 1
    assert np.array_equal(selected, chain_i)


def test_tree_entries_select_what_they_name(panel) -> None:
    top = panel.tree.topLevelItem(0)
    labels = [top.child(i).text(0) for i in range(top.childCount())]
    assert labels[:3] == ["Chain E (336)", "Chain I (20)", "ATP E355"]
    assert "Ions (2)" in labels and "Waters (103)" in labels
    panel.tree.select_item(top.child(2))
    s = structure("1atp.cif.gz")
    assert panel.session.display("1ATP").selected.sum() == 31
    panel.tree.select_item(top.child(1))
    selected = panel.session.display("1ATP").selected
    assert np.array_equal(np.unique(s.atom_chain[selected]), [1])


def test_clicking_letters_selects_residues(panel, qtbot) -> None:
    strip = panel.sequence.strip
    row = strip.rows[0]
    first = int(np.flatnonzero(row.sequence.modelled)[0])

    def point(position: int) -> QPoint:
        return QPoint(strip.label_width + position * strip.cell + 2, strip.line + 2)

    qtbot.mouseClick(strip, Qt.MouseButton.LeftButton, pos=point(first))
    shown = panel.session.display("1ATP")
    s = shown.structure
    residues = np.unique(s.residue_of[shown.selected])
    assert residues.tolist() == [int(row.sequence.residue_at[first])]
    qtbot.mouseClick(strip, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.ShiftModifier,
                     pos=point(first + 4))  # fmt: skip
    assert len(np.unique(s.residue_of[shown.selected])) == 5
    # a missing residue (the unmodelled start) selects nothing new
    missing = int(np.flatnonzero(~row.sequence.modelled)[0])
    qtbot.mouseClick(strip, Qt.MouseButton.LeftButton, pos=point(missing))
    assert not shown.selected.any()


def test_clicks_in_the_view_select_residues_atoms_and_chains(panel) -> None:
    session = panel.session
    s = structure("1atp.cif.gz")
    atom = int(np.flatnonzero(s.name == "CA")[20])
    panel._clicked(atom, Qt.KeyboardModifier.NoModifier)
    assert session.display("1ATP").selected.sum() == len(
        range(*s.atoms_of(int(s.residue_of[atom])).indices(s.n_atoms))
    )
    panel._clicked(atom, Qt.KeyboardModifier.AltModifier)
    assert session.display("1ATP").selected.sum() == 1
    panel._double_clicked(atom)
    assert session.display("1ATP").selected.sum() == int((s.atom_chain == 0).sum())
    panel._clicked(-1, Qt.KeyboardModifier.NoModifier)
    assert session.selected_count() == 0


# showing


def test_show_hide_and_colour_act_on_the_selection(panel, qtbot) -> None:
    session = panel.session
    session.select_text("chain I")
    with qtbot.waitSignal(session.sceneChanged):
        session.show_as(scene.SPACEFILL)
    shown = session.display("1ATP")
    chain_i = shown.structure.atom_chain == 1
    assert (shown.styles[chain_i & shown.active] == scene.SPACEFILL).all()
    session.show_cartoon(False)
    assert not shown.cartoon[shown.structure.chain_of == 1].any()
    session.select_text("ligand")
    session.show_surface(True)
    assert len(session.scene.surface) > 0
    assert set(np.unique(session.scene.surface.owner)) <= set(np.flatnonzero(
        shown.structure.atom_kind == 2).tolist())  # fmt: skip
    session.colour((1.0, 0.0, 1.0))
    assert np.allclose(shown.colours[shown.structure.atom_kind == 2, :3], (1, 0, 1))


def test_the_context_menu_offers_what_was_clicked(panel) -> None:
    from chemsuite.app.protein.menus import context_menu

    s = structure("1atp.cif.gz")
    atom = int(np.flatnonzero(s.atom_kind == 2)[0])
    menu = context_menu(panel.session, atom, panel, focus=panel.frame_selection)
    titles = [a.text() for a in menu.actions()]
    assert titles[0] == s.atom_label(atom)
    assert "Select Residue" in titles and "Select Chain" in titles
    next(a for a in menu.actions() if a.text() == "Select Residues Within 5 Å").trigger()
    assert panel.session.selected_count() > 31
    empty = context_menu(panel.session, -1, panel, focus=panel.frame_selection)
    assert [a.text() for a in empty.actions()][:3] == ["Show", "Hide", "Colour"]


def test_alternate_locations_can_be_switched(qtbot, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    t = Table()
    add_protein(t, "3HTB", "3htb.cif.gz")
    panel = ProteinPanel(t)
    qtbot.addWidget(panel)
    shown = panel.session.display("3HTB")
    residue = next(iter(shown.structure.default_altlocs))
    other = "B" if shown.structure.default_altlocs[residue] == "A" else "A"
    span = shown.structure.atoms_of(residue)
    panel.session.choose_altloc("3HTB", residue, other)
    visible = shown.styles[span] != scene.HIDDEN
    assert set(shown.structure.altloc[span][visible]) <= {"", other}


# saving and the GPU fallback


def test_the_panel_reopens_as_it_was(panel, table, qtbot, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    session = panel.session
    session.select_text("chain I")
    session.show_as(scene.STICK)
    session.colour_by("chain", everything=True)
    session.save_selection("peptide")
    panel.view.camera.zoom(0.5)
    state = panel.state()
    monkeypatch.setattr(graphics, "_problem", None)
    again = ProteinPanel(table)
    qtbot.addWidget(again)
    again.set_state(state)
    old, new = session.display("1ATP"), again.session.display("1ATP")
    assert np.array_equal(old.styles, new.styles)
    assert np.array_equal(old.cartoon, new.cartoon)
    assert new.scheme == "chain"
    assert again.session.named == session.named
    assert again.view.camera.distance == pytest.approx(panel.view.camera.distance)


def test_the_models_bar_shows_only_while_a_shown_protein_has_several(panel, table) -> None:
    bar = panel.models_bar
    assert not bar.isVisibleTo(panel)  # 1ATP has one model
    add_protein(table, "1D3Z", "1d3z_2models.cif.gz")
    assert bar.isVisibleTo(panel)
    assert bar.count.text() == "of 2" and bar.position.value() == 1
    assert not bar._protein_action.isVisible()  # one protein to step: no list
    panel.session.set_visible("1D3Z", False)
    assert not bar.isVisibleTo(panel)


def test_stepping_models_moves_the_atoms_and_keeps_the_selection(panel, table, qtbot) -> None:
    add_protein(table, "1D3Z", "1d3z_2models.cif.gz")
    session, bar = panel.session, panel.models_bar
    nmr = session.display("1D3Z")
    session.select_text("resi 10-20")
    selected = session.selected_count()
    assert not bar.previous.isEnabled() and bar.next.isEnabled()
    bar.next.trigger()
    assert nmr.model == 1 and bar.position.value() == 2 and bar.slider.value() == 1
    drawn = session.scene
    d = drawn.displays.index(nmr)
    start, stop = drawn.atom_offsets[d], drawn.atom_offsets[d + 1]
    own = (drawn.pick_atoms >= start) & (drawn.pick_atoms < stop)
    local = drawn.pick_atoms[own] - start
    assert own.any() and np.allclose(drawn.pick_centres[own], nmr.structure.xyz(1)[local])
    assert session.selected_count() == selected
    assert bar.previous.isEnabled() and not bar.next.isEnabled()
    bar.position.setValue(1)  # typed: back to the first
    assert nmr.model == 0
    panel.models_bar.slider.setValue(1)
    assert nmr.model == 1
    state = panel.state()
    again = type(panel)(table)
    qtbot.addWidget(again)
    again.set_state(state)
    assert again.session.display("1D3Z").model == 1
    assert again.models_bar.position.value() == 2


def test_playing_loops_and_stops_when_the_protein_goes(panel, table) -> None:
    add_protein(table, "1D3Z", "1d3z_2models.cif.gz")
    bar, nmr = panel.models_bar, panel.session.display("1D3Z")
    bar.play.setChecked(True)
    assert bar.playing and bar.play.text() == "Pause"
    bar._tick()
    bar._tick()
    assert nmr.model == 0  # 0 → 1 → back to 0
    table.undo()  # the protein is gone
    assert not bar.playing and not bar.play.isChecked()


def test_the_list_picks_which_protein_the_bar_steps(panel, table) -> None:
    add_protein(table, "1D3Z", "1d3z_2models.cif.gz")
    add_protein(table, "1D3Z-2", "1d3z_2models.cif.gz")
    bar = panel.models_bar
    assert bar._protein_action.isVisible() and bar.protein.count() == 2
    bar.protein.setCurrentIndex(1)
    bar.next.trigger()
    assert panel.session.display("1D3Z-2").model == 1
    assert panel.session.display("1D3Z").model == 0


class _Job(QObject):
    done = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.cancelled = False

    def cancel(self) -> None:
        self.cancelled = True


class _Runner:
    """Holds background work until a test finishes it, in any order."""

    def __init__(self) -> None:
        self.jobs: list[tuple[str, object, object, _Job]] = []

    def __call__(self, label, work, on_done) -> _Job:
        job = _Job()
        self.jobs.append((label, work, on_done, job))
        return job

    def finish(self, index: int) -> None:
        _label, work, on_done, job = self.jobs[index]
        on_done(work(job))
        job.done.emit()


@pytest.fixture
def deferring(qtbot, table, monkeypatch):
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    runner = _Runner()
    widget = ProteinPanel(table, run=runner)
    qtbot.addWidget(widget)
    return widget, runner


def test_a_surface_is_built_in_the_background_and_a_stale_one_dropped(deferring) -> None:
    panel, runner = deferring
    session = panel.session
    session.select_text("chain E")
    session.show_surface(True)
    assert not len(session.scene.surface.triangles)  # drawn without waiting
    assert [label for label, *_ in runner.jobs] == ["Building a protein surface"]
    session.built()  # drawn again: asking for the version being built starts nothing
    assert len(runner.jobs) == 1
    session.select_text("chain I")
    session.show_surface(True)  # changed before the first arrived
    assert len(runner.jobs) == 2 and runner.jobs[0][3].cancelled
    runner.finish(1)
    newest = session.scene.surface
    assert len(newest.triangles)
    runner.finish(0)  # the stale one arrives last: dropped
    assert session.scene.surface is newest  # nothing was redrawn
    assert not panel.session.builder.busy()


def test_interactions_are_found_in_the_background(deferring) -> None:
    panel, runner = deferring
    session = panel.session
    session.show_interactions(True)
    assert not len(session.scene.lines.triangles)
    runner.finish(0)
    assert len(session.scene.lines.triangles)


def test_work_for_a_protein_that_went_is_ignored(deferring, table) -> None:
    panel, runner = deferring
    panel.session.select_text("all")
    panel.session.show_surface(True)
    assert len(runner.jobs) == 1
    table.undo()  # the protein is gone before its surface arrives
    runner.finish(0)
    assert panel.session.displays == [] and not panel.session.builder.busy()


def test_closing_the_viewer_cancels_its_work(deferring) -> None:
    panel, runner = deferring
    panel.session.select_text("all")
    panel.session.show_surface(True)
    assert len(runner.jobs) == 1
    panel.close()
    assert runner.jobs[0][3].cancelled


def test_state_for_a_protein_not_loaded_yet_waits_for_it(qtbot, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel
    from chemsuite.protein.display import runs

    monkeypatch.setattr(graphics, "_problem", None)
    t = Table()
    panel = ProteinPanel(t)
    qtbot.addWidget(panel)
    n = structure("1crn.cif.gz").n_atoms
    panel.set_state({"session": {"proteins": {"1CRN": {"styles": runs(np.full(n, 4))}}}})
    add_protein(t, "1CRN", "1crn.cif.gz")
    assert (panel.session.display("1CRN").styles == scene.SPACEFILL).all()


def test_bad_state_is_ignored(panel) -> None:
    panel.set_state({"session": {"proteins": {"1ATP": {"styles": [1, 2, 3], "scheme": "x"}}},
                     "columns": "junk", "camera": {"distance": "far"}})  # fmt: skip
    shown = panel.session.display("1ATP")
    assert len(shown.styles) == shown.structure.n_atoms and shown.scheme == "element"


def test_without_3d_the_view_says_why_and_the_rest_works(qtbot, table, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", graphics.TURNED_OFF)
    panel = ProteinPanel(table)
    qtbot.addWidget(panel)
    assert panel.view is None
    assert panel.views.currentWidget().text() == f"⚠ {graphics.TURNED_OFF}"
    panel.selection_bar.field.setText("ligand")
    assert panel.selection_bar.apply()
    panel.frame_selection()
    panel.reset_view()
    assert "camera" not in panel.state()


def test_a_view_that_fails_when_drawn_is_replaced(panel, monkeypatch) -> None:
    from chemsuite.app import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    panel.view.failed.emit()
    assert panel.view is None
    assert panel.views.currentWidget().text().startswith("⚠")


def test_stored_proteins_round_trip_through_a_document(tmp_path) -> None:
    from chemsuite.core.document import load_document, save_document

    t = Table()
    add_protein(t, "1ATP", "1atp.cif.gz")  # held as a Structure, written as bytes
    path = tmp_path / "protein.csd"
    save_document(t, path)
    back, _extras = load_document(path)
    restored = store.unpack(back.assets["1ATP"].data())
    assert np.array_equal(restored.uid, structure("1atp.cif.gz").uid)


# pocket work: interactions, focus, measurements


def _atp_atom(name: str) -> int:
    s = structure("1atp.cif.gz")
    return int(np.flatnonzero((s.res_name[s.residue_of] == "ATP") & (s.name == name))[0])


def test_interactions_of_the_ligands_or_the_selection(panel) -> None:
    session = panel.session
    panel.interactions_action.trigger()  # nothing selected: every ligand
    assert panel.interactions_action.isChecked()
    shown = session.display("1ATP")
    s = shown.structure
    assert np.array_equal(shown.query, s.res_kind[s.residue_of] == 2)
    assert len(session.scene.lines) > 0
    assert "interactions" in panel.view.mesh_names
    session.select_text("chain I")
    session.show_interactions()
    assert np.array_equal(shown.query, s.atom_chain == 1)
    panel.interactions_action.trigger()
    assert not panel.interactions_action.isChecked() and not shown.query.any()
    assert len(session.scene.lines) == 0


def test_hovering_an_atom_names_its_interactions(panel) -> None:
    panel.session.show_interactions()
    panel._hovered("1ATP", _atp_atom("N6"))
    assert "Hydrogen Bond" in panel.status.text() and "GLU 121/O" in panel.status.text()


def test_focus_shows_the_pocket_and_fades_the_rest(panel, qtbot) -> None:
    session = panel.session
    s = structure("1atp.cif.gz")
    atp = int(s.residue_of[_atp_atom("N6")])
    with qtbot.waitSignal(session.focusChanged):
        session.focus("1ATP", atp)
    shown = session.display("1ATP")
    pocket = ~shown.ghost
    assert pocket[s.residue_of == atp].all()
    site = np.unique(s.residue_of[pocket & (s.res_kind[s.residue_of] == 0)])
    assert {"LYS", "GLU", "VAL"} <= set(s.res_name[site].tolist()) and len(site) < 40
    assert (
        shown.styles[pocket & shown.active & (s.res_kind[s.residue_of] == 0)] == scene.STICK
    ).all()
    assert shown.query[s.residue_of == atp].all() and len(session.scene.lines) > 0
    ions = s.res_kind[s.residue_of] == 4
    assert pocket[ions].all()  # both Mn bind the phosphates
    far_waters = (s.res_kind[s.residue_of] == 3) & ~pocket
    assert far_waters.any() and (shown.styles[far_waters] == scene.HIDDEN).all()
    near, far = panel.view.slab
    assert 0 < far - near < 60
    panel.view.camera.zoom(0.5)
    panel.view.camera_changed()
    assert panel.view.slab != (near, far)  # the slab follows the camera
    session.clear_focus()
    assert not shown.ghost.any() and panel.view.slab[1] > 1e8


def test_the_focus_menu_lists_the_ligands(panel) -> None:
    from PySide6.QtWidgets import QMenu

    menu = QMenu()
    panel.fill_focus(menu)
    titles = [a.text() for a in menu.actions() if a.text()]
    assert titles == ["E/ATP 355", "Clear Focus"]
    menu.actions()[0].trigger()
    assert panel.session.focused_residue() is not None
    panel.fill_focus(menu)
    assert menu.actions()[0].isChecked() and menu.actions()[-1].isEnabled()


def test_measuring_takes_the_clicks(panel) -> None:
    session = panel.session
    a, b, c = _atp_atom("PG"), _atp_atom("PB"), _atp_atom("PA")
    session.measure(2)
    panel._clicked(a, Qt.KeyboardModifier.NoModifier)
    assert session.picked and session.selected_count() == 0  # not a selection
    panel._clicked(b, Qt.KeyboardModifier.NoModifier)
    assert len(session.measurements) == 1 and not session.picked
    (points, text), = session.measured()  # fmt: skip
    assert text.endswith(" Å") and 2.5 < float(text.split()[0]) < 3.5
    assert panel.overlay.labels()[0][1] == text
    session.measure(3)
    for atom in (a, b, c):
        panel._clicked(atom, Qt.KeyboardModifier.NoModifier)
    assert session.measured()[1][1].endswith("°")
    panel._escape()
    assert session.measuring == 0
    panel._clicked(a, Qt.KeyboardModifier.NoModifier)
    assert session.selected_count() > 0  # selecting again
    session.clear_measurements()
    assert session.measured() == []


def test_pocket_work_reopens_as_it_was(panel, table, qtbot, monkeypatch) -> None:
    from chemsuite.app import graphics
    from chemsuite.app.protein.panel import ProteinPanel

    session = panel.session
    atp = int(structure("1atp.cif.gz").residue_of[_atp_atom("N6")])
    session.focus("1ATP", atp)
    session.measure(2)
    panel._clicked(_atp_atom("PG"), Qt.KeyboardModifier.NoModifier)
    panel._clicked(_atp_atom("PB"), Qt.KeyboardModifier.NoModifier)
    state = panel.state()
    monkeypatch.setattr(graphics, "_problem", None)
    again = ProteinPanel(table)
    qtbot.addWidget(again)
    again.resize(900, 600)
    again.set_state(state)
    old, new = session.display("1ATP"), again.session.display("1ATP")
    assert np.array_equal(old.ghost, new.ghost) and np.array_equal(old.query, new.query)
    assert again.session.measured()[0][1] == session.measured()[0][1]
    assert again.session.focused == session.focused
    assert again.view.slab[1] < 1e8


def test_bad_pocket_state_is_ignored(panel) -> None:
    measurements = [[["1ATP", 1]], "junk", [["1ATP", "y"]]]
    panel.set_state({"session": {"focused": ["1ATP", "x"], "measurements": measurements}})
    assert panel.session.focused is None
    assert panel.session.measurements == [(("1ATP", 1),)]
