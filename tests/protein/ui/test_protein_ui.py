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
from chemsuite.protein.look import scene
from chemsuite.protein.structure import io, store
from chemsuite.protein.structure.model import Structure
from chemsuite.protein.structure.select import SelectionError, select
from chemsuite.protein.ui.sequence_view import POLYMER

DATA = Path(__file__).resolve().parents[2] / "data" / "protein"


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


def menu_of(panel, title: str):
    """The viewer's menu-bar menu *title*, filled as when it opens."""
    menu = next(a.menu() for a in panel.menu_bar.actions() if a.text() == title)
    menu.aboutToShow.emit()
    return menu


def action_in(menu, text: str):
    return next(a for a in menu.actions() if a.text() == text)


@pytest.fixture
def panel(qtbot, table, monkeypatch):
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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


def test_a_reopened_protein_keeps_its_display_when_another_arrives(qtbot) -> None:
    """A reopened document holds stored bytes: adding a protein (or any asset undo) must not
    decode it again and lose its selection and focus."""
    from chemsuite.protein.ui.panel import ProteinPanel

    table = Table()
    add_protein(table, "1ATP", "1atp.cif.gz", packed=True)
    widget = ProteinPanel(table)
    qtbot.addWidget(widget)
    session = widget.session
    shown = session.display("1ATP")
    session.select_text("resn ATP")
    add_protein(table, "3HTB", "3htb.cif.gz", packed=True)
    assert session.display("1ATP") is shown and shown.selected.any()
    table.undo()
    assert session.display("1ATP") is shown


def test_hiding_a_protein_takes_it_out_of_the_scene(panel, table) -> None:
    add_protein(table, "3HTB", "3htb.cif.gz")
    session = panel.session
    both = session.scene.atom_offsets[-1]
    item = panel.tree.topLevelItem(1)
    item.setCheckState(0, Qt.CheckState.Unchecked)
    # the tree is not rebuilt under Qt's feet: the ticked item lives on (a crash on Linux)
    assert panel.tree.topLevelItem(1) is item
    assert not session.display("3HTB").visible
    assert session.scene.atom_offsets[-1] == both - structure("3htb.cif.gz").n_atoms
    chains = [r.label for r in panel.sequence.strip.rows if r.kind == POLYMER]
    assert chains == ["1ATP E", "1ATP I"]  # 1ATP's chains only


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


def test_clicks_in_the_view_select_atoms_residues_and_chains(panel) -> None:
    session = panel.session
    shown = session.display("1ATP")
    s = shown.structure
    atom = int(np.flatnonzero(s.name == "CA")[20])
    residue = (s.residue_of == s.residue_of[atom]) & shown.active
    # under a cartoon no atom is drawn: a click there takes the residue
    assert shown.styles[atom] == scene.HIDDEN
    panel._clicked(atom, Qt.KeyboardModifier.NoModifier)
    assert np.array_equal(shown.selected, residue)
    session.show_as(scene.STICK)  # the residue drawn: now a click takes one atom
    panel._clicked(atom, Qt.KeyboardModifier.NoModifier)
    assert np.flatnonzero(shown.selected).tolist() == [atom]
    panel._clicked(atom + 1, Qt.KeyboardModifier.ShiftModifier)
    assert np.flatnonzero(shown.selected).tolist() == [atom, atom + 1]
    panel._double_clicked(atom)
    assert np.array_equal(shown.selected, residue)
    panel._clicked(atom, Qt.KeyboardModifier.AltModifier)
    assert np.array_equal(shown.selected, (s.atom_chain == 0) & shown.active)
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


def test_hydrogens_follow_their_atom_and_the_choice(qtbot, monkeypatch) -> None:
    from PySide6.QtWidgets import QMenu

    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.menus import fill_hydrogens
    from chemsuite.protein.ui.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    table = Table()
    add_protein(table, "1D3Z", "1d3z_2models.cif.gz")  # NMR: it has hydrogens
    panel = ProteinPanel(table)
    qtbot.addWidget(panel)
    session = panel.session
    nmr = session.display("1D3Z")
    s = nmr.structure
    session.select_text("resi 10-20")
    session.show_as(scene.STICK)

    def drawn(text: str) -> int:
        atoms = select(s, text)
        return int(np.isin(session.scene.geometry.atoms, np.flatnonzero(atoms)).sum())

    polar, nonpolar = "resi 10-20 and hydrogen and name H", "resi 10-20 and name HA"
    assert drawn(polar) > 0 and drawn(nonpolar) == 0  # the default: polar ones
    menu = QMenu()
    fill_hydrogens(menu, session)
    assert [a.text() for a in menu.actions() if a.isChecked()] == ["Polar Hydrogens"]
    session.show_hydrogens("all")
    assert drawn(nonpolar) == len(np.flatnonzero(select(s, nonpolar)))
    session.show_hydrogens("none")
    assert drawn("hydrogen") == 0
    assert nmr.styles[select(s, polar)].max() == scene.STICK  # styles are never rewritten


def test_surface_transparency_changes_only_the_look(panel, qtbot) -> None:
    from PySide6.QtWidgets import QMenu

    from chemsuite.protein.ui.menus import fill_transparency

    session = panel.session
    session.select_text("chain I")
    session.show_surface(True)
    assert session.scene.surface_look[:, 3] == pytest.approx(0.6)
    with qtbot.assertNotEmitted(session.sceneChanged):
        session.set_transparency(0.8)
    assert session.scene.surface_look[:, 3] == pytest.approx(0.2)
    menu = QMenu()
    fill_transparency(menu, session)
    assert [a.text() for a in menu.actions()] == ["Opaque", "20%", "40%", "60%", "80%"]
    assert [a.text() for a in menu.actions() if a.isChecked()] == ["80%"]
    again = session.display("1ATP").state()
    assert again["transparency"] == 0.8 and again["hydrogens"] == "polar"


def test_label_layout_goes_front_to_back_and_skips_overlaps() -> None:
    from chemsuite.protein.ui import labels

    boxes = np.array([[0, 0, 50, 10], [10, 5, 50, 10], [100, 0, 50, 10], [200, 0, 50, 10]])
    depth = np.array([5.0, 1.0, 3.0, 40.0])
    assert labels.layout(boxes, depth).tolist() == [1, 2, 3]  # 0 hides behind 1
    assert labels.layout(boxes, depth, near=2.0, far=30.0).tolist() == [2, 0]  # in the slab
    assert labels.layout(boxes, -depth).tolist() == []  # behind the camera
    many = np.array([[i * 60.0, 0, 50, 10] for i in range(400)])
    assert len(labels.layout(many, np.arange(400.0))) == labels.MAX_LABELS


def test_residues_and_atoms_are_labelled_and_kept_by_uid(panel, table, qtbot) -> None:
    from chemsuite.protein.ui.session import ATOM_LABEL, RESIDUE_LABEL

    session = panel.session
    s = structure("1atp.cif.gz")
    session.select_text("chain E and resi 72")
    with qtbot.waitSignal(session.labelsChanged):
        session.label(RESIDUE_LABEL)
    positions, texts = session.labelled()
    ca = np.flatnonzero(select(s, "chain E and resi 72 and name CA"))[0]
    assert texts == ["LYS 72"] and np.allclose(positions[0], s.xyz(0)[ca])
    session.select_text("resn ATP and name N6")
    session.label(ATOM_LABEL)
    session.select_text("resn ATP")
    session.label(RESIDUE_LABEL)  # a ligand's label sits in its middle
    positions, texts = session.labelled()
    assert texts == ["LYS 72", "N6", "ATP 355"]
    atp = select(s, "resn ATP") & s.active()
    assert np.allclose(positions[2], s.xyz(0)[atp].mean(axis=0), atol=1e-4)
    session.focus("1ATP", int(s.residue_of[np.flatnonzero(atp)[0]]))
    assert "LYS 72" in session.labelled()[1]  # in the pocket: not faded
    session.select_text("chain E and resi 200")
    session.label(RESIDUE_LABEL)
    assert "GLU 200" not in session.labelled()[1]  # faded by the focus
    state = panel.state()
    again = type(panel)(table)
    qtbot.addWidget(again)
    again.set_state(state)
    assert again.session.labels == session.labels
    session.clear_labels()
    assert session.labelled()[1] == []
    again.set_state({"session": {"labels": [["1ATP", "x", "atom"], ["1ATP", 3, "bond"], 7]}})
    assert again.session.labels == {}


def test_the_overlay_draws_labels_without_overlaps(panel) -> None:
    from chemsuite.protein.ui.session import ATOM_LABEL

    panel.view.resize(800, 600)
    panel.overlay.resize(800, 600)
    panel.view.fit()
    panel.session.select_text("resn ATP")
    panel.session.label(ATOM_LABEL)
    boxes, texts = panel.overlay.placed_labels()
    assert 0 < len(texts) <= 31
    for i in range(len(boxes)):
        x, y, w, h = boxes[i]
        others = np.delete(boxes, i, axis=0)
        assert not (
            (others[:, 0] < x + w) & (x < others[:, 0] + others[:, 2])
            & (others[:, 1] < y + h) & (y < others[:, 1] + others[:, 3])
        ).any()  # fmt: skip
    panel.overlay.grab()  # paints without error


def test_the_context_menu_offers_what_was_clicked(panel) -> None:
    from chemsuite.protein.ui.menus import context_menu

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


def test_the_context_menu_labels_measures_and_adds_a_ligand(panel) -> None:
    from chemsuite.protein.ui.menus import context_menu

    s = structure("1atp.cif.gz")
    atom = int(np.flatnonzero(s.atom_kind == 2)[0])
    added = []
    menu = context_menu(panel.session, atom, panel, focus=panel.frame_selection,
                        add_ligand=lambda p, r: added.append((p, r)))  # fmt: skip
    actions = {a.text(): a for a in menu.actions()}
    actions["Add Ligand to Table"].trigger()
    assert added == [("1ATP", int(s.residue_of[atom]))]
    labels = {a.text(): a for a in actions["Label"].menu().actions()}
    labels["Residue"].trigger()
    assert panel.session.labelled()[1] == ["ATP 355"]
    labels["Atom"].trigger()
    assert panel.session.labelled()[1] == ["ATP 355", str(s.name[atom])]
    measuring = {a.text(): a for a in actions["Measure From Here"].menu().actions()}
    measuring["Angle"].trigger()
    assert panel.session.measuring == 3 and panel.session.picked == [("1ATP", int(s.uid[atom]))]
    protein_atom = int(np.flatnonzero(s.atom_kind == 0)[0])
    plain = context_menu(panel.session, protein_atom, panel, focus=panel.frame_selection,
                         add_ligand=lambda p, r: None)  # fmt: skip
    assert "Add Ligand to Table" not in [a.text() for a in plain.actions()]


def test_around_adds_what_lies_near_the_selection(panel) -> None:
    from chemsuite.protein.ui.selection_bar import AROUND

    bar, session = panel.selection_bar, panel.session
    around = action_in(menu_of(panel, "Select"), "Around")
    assert not around.isEnabled()  # nothing selected
    bar.pick("resn ATP")
    assert bar.field.text() == "resn ATP"
    atp = session.selected_count()
    around = action_in(menu_of(panel, "Select"), "Around").menu()
    titles = [a.text() for a in around.actions() if a.text()]
    assert len(titles) == 2 * len(AROUND) and titles[1] == "Residues Within 5 Å"
    next(a for a in around.actions() if a.text() == "Residues Within 5 Å").trigger()
    shown = session.display("1ATP")
    s = shown.structure
    grown = session.selected_count()
    assert grown > atp and bar.field.text() == ""
    for residue in np.unique(s.residue_of[shown.selected]):  # whole residues
        span = s.atoms_of(int(residue))
        assert shown.selected[span][shown.active[span]].all()
    with pytest.raises(SelectionError):
        session.save_selection("selection")


def test_alternate_locations_can_be_switched(qtbot, monkeypatch) -> None:
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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
    from chemsuite.app.shell import graphics
    from chemsuite.protein.look.display import runs
    from chemsuite.protein.ui.panel import ProteinPanel

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
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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
    from chemsuite.app.shell import graphics

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
    action_in(menu_of(panel, "Ligand"), "Interactions").trigger()  # nothing selected: ligands
    assert action_in(menu_of(panel, "Ligand"), "Interactions").isChecked()
    shown = session.display("1ATP")
    s = shown.structure
    assert np.array_equal(shown.query, s.res_kind[s.residue_of] == 2)
    assert len(session.scene.lines) > 0
    assert "interactions" in panel.view.mesh_names
    session.select_text("chain I")
    session.show_interactions()
    assert np.array_equal(shown.query, s.atom_chain == 1)
    action_in(menu_of(panel, "Ligand"), "Interactions").trigger()
    assert not action_in(menu_of(panel, "Ligand"), "Interactions").isChecked()
    assert not shown.query.any()
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
    for residue in np.unique(s.residue_of[pocket]):  # whole residues, never half faded
        span = s.atoms_of(int(residue))
        assert pocket[span][shown.active[span]].all(), s.residue_label(int(residue))
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
    fog_start, fog_end = panel.view.fog
    assert near < fog_start < fog_end == pytest.approx(far)  # the back half of the pocket fades
    panel.view.camera.zoom(0.5)
    panel.view.camera_changed()
    assert panel.view.slab != (near, far)  # the slab follows the camera
    session.clear_focus()
    assert not shown.ghost.any() and panel.view.slab[1] > 1e8 and panel.view.fog is None


def test_focus_draws_a_skin_over_the_pocket(panel) -> None:
    session = panel.session
    s = structure("1atp.cif.gz")
    atp = int(s.residue_of[_atp_atom("N6")])
    session.focus("1ATP", atp)
    shown = session.display("1ATP")
    skin = session.scene.pocket
    assert len(skin.triangles)
    owners = np.unique(skin.owner)
    assert (s.res_kind[s.residue_of[owners]] == 0).all()  # the protein's side of the pocket
    assert not shown.ghost[owners].any() and (s.element[owners] != 1).all()
    assert session.scene.pocket_look[:, 3] == pytest.approx(1 - shown.transparency)
    session.show_pocket_surface(False)
    assert not len(session.scene.pocket.triangles)
    session.show_pocket_surface(True)
    assert len(session.scene.pocket.triangles)
    session.clear_focus()
    assert not len(session.scene.pocket.triangles)


def test_a_saved_focus_brings_its_pocket_skin_back(panel, table, qtbot, monkeypatch) -> None:
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

    s = structure("1atp.cif.gz")
    panel.session.focus("1ATP", int(s.residue_of[_atp_atom("N6")]))
    state = panel.state()
    monkeypatch.setattr(graphics, "_problem", None)
    again = ProteinPanel(table)
    qtbot.addWidget(again)
    again.set_state(state)
    assert len(again.session.scene.pocket.triangles)
    panel.session.show_pocket_surface(False)
    later = ProteinPanel(table)
    qtbot.addWidget(later)
    later.set_state(panel.state())
    assert not later.session.pocket_surface and not len(later.session.scene.pocket.triangles)


def test_the_binding_site_picks_whole_residues(panel) -> None:
    from chemsuite.protein.ui.selection_bar import PICKS

    session = panel.session
    for title in ("Binding Site", "Ligands and Binding Site"):
        session.select_text(dict(PICKS)[title])
        shown = session.display("1ATP")
        s = shown.structure
        for residue in np.unique(s.residue_of[shown.selected & (s.res_kind[s.residue_of] == 0)]):
            span = s.atoms_of(int(residue))
            assert shown.selected[span][shown.active[span]].all(), (title, s.residue_label(residue))


def test_the_focus_menu_lists_the_ligands(panel) -> None:
    ligand = menu_of(panel, "Ligand")
    focus = action_in(ligand, "Focus").menu()
    assert [a.text() for a in focus.actions()] == ["E/ATP 355"]
    assert not action_in(ligand, "Clear Focus").isEnabled()
    focus.actions()[0].trigger()
    assert panel.session.focused_residue() is not None
    ligand = menu_of(panel, "Ligand")
    assert action_in(ligand, "Focus").menu().actions()[0].isChecked()
    assert action_in(ligand, "Clear Focus").isEnabled()


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
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

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


def test_every_part_built_in_the_background_has_a_job_title() -> None:
    from chemsuite.protein.look import display
    from chemsuite.protein.ui import builds

    assert set(display.DEFERRED) <= set(builds.TITLES)


# the polish pass: selection feedback, keys, light and dark


def test_the_tree_shows_what_is_selected(panel) -> None:
    from PySide6.QtCore import Qt

    tree = panel.tree
    top = tree.topLevelItem(0)
    chain_e, chain_i = top.child(0), top.child(1)
    panel.session.select_text("chain I and polymer")
    assert chain_i.background(0).style() != Qt.BrushStyle.NoBrush  # all of it: highlighted
    assert chain_e.background(0).style() == Qt.BrushStyle.NoBrush and not chain_e.font(0).bold()
    assert top.font(0).bold()  # some of the protein
    panel.session.clear_selection()
    assert chain_i.background(0).style() == Qt.BrushStyle.NoBrush and not top.font(0).bold()


def test_the_viewers_own_keys_while_focus_is_inside(panel, qtbot) -> None:
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QGuiApplication

    from chemsuite.app.shell.windows import ToolWindow

    window = ToolWindow(panel, "Protein Viewer", name=lambda: "", menu_bar=panel.menu_bar)
    with qtbot.waitActive(window):  # a shortcut needs its window active
        window.show()
        window.activateWindow()
    panel.tree.setFocus()
    qtbot.keyClick(panel.tree, Qt.Key.Key_A, Qt.KeyboardModifier.ControlModifier)
    shown = panel.session.display("1ATP")
    assert shown.selected[shown.active].all()
    panel.session.select_text("ligand")
    qtbot.keyClick(panel.tree, Qt.Key.Key_C, Qt.KeyboardModifier.ControlModifier)
    assert QGuiApplication.clipboard().text().count("HETATM") == 31
    framed: list[int] = []
    panel.frame_selection = lambda: framed.append(1)
    qtbot.keyClick(panel.tree, Qt.Key.Key_F)
    assert framed == [1]
    field = panel.selection_bar.field  # a text field keeps its own ⌘A and its letters
    assert field.text() == "all"  # Select ▸ All shows what it selected, as the other picks do
    field.clear()
    field.setFocus()
    qtbot.keyClicks(field, "f")
    assert field.text() == "f" and framed == [1]
    panel.setParent(None)  # the fixture's panel outlives the window


def test_the_sequence_strip_by_keyboard(panel, qtbot) -> None:
    from PySide6.QtCore import Qt

    strip = panel.sequence.strip
    shown = panel.session.display("1ATP")
    s = shown.structure
    modelled = np.flatnonzero(strip.rows[0].sequence.residue_at >= 0)
    at = int(modelled[0]) + 3  # the first residues of the sequence are not in the model
    strip.setFocus()
    for _ in range(at):
        qtbot.keyClick(strip, Qt.Key.Key_Right)
    assert strip._cursor == (0, at) and not shown.selected.any()
    qtbot.keyClick(strip, Qt.Key.Key_Return)
    first = np.unique(s.residue_of[shown.selected])
    assert len(first) == 1 and first[0] == strip.residue_at(0, at)
    for _ in range(2):
        qtbot.keyClick(strip, Qt.Key.Key_Right, Qt.KeyboardModifier.ShiftModifier)
    assert len(np.unique(s.residue_of[shown.selected])) == 3
    qtbot.keyClick(strip, Qt.Key.Key_Down)
    assert strip._cursor[0] == 1
    qtbot.keyClick(strip, Qt.Key.Key_End)
    assert strip._cursor == (1, len(strip.rows[1].letters) - 1)


def test_the_view_follows_light_and_dark(panel) -> None:
    from PySide6.QtGui import QColor, QPalette

    if panel.view is None:
        return
    dark = QPalette(panel.palette())
    dark.setColor(QPalette.ColorRole.Base, QColor(20, 20, 24))
    panel.setPalette(dark)
    assert panel.view.background == QColor(20, 20, 24)


def test_closing_the_viewer_window_lets_go_of_the_document(qtbot, table, monkeypatch) -> None:
    """Closing the window must unsubscribe the session: a deleted session left listening
    raised on the next table change, and the listeners after it never heard of the change."""
    from PySide6.QtCore import QCoreApplication, QEvent

    from chemsuite.app.shell import graphics
    from chemsuite.app.shell.windows import ToolWindow
    from chemsuite.protein.ui.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    stopped = []
    panel = ProteinPanel(table, watch_current=lambda follow: lambda: stopped.append(True))
    window = ToolWindow(panel, "Protein Viewer", name=lambda: "Kinases",
                        menu_bar=panel.menu_bar, on_close=panel.close_panel)  # fmt: skip
    assert window.windowTitle() == "Protein Viewer — Kinases"
    assert window.menuBar() is panel.menu_bar
    window.show()
    window.close()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    assert stopped == [True]
    heard = []
    table.subscribe(heard.append)
    add_protein(table, "2ATP", "1atp.cif.gz")
    assert heard


def test_the_menu_bar_holds_every_command(panel) -> None:
    from chemsuite.protein.ui.menubar import MENUS

    assert [a.text() for a in panel.menu_bar.actions()] == list(MENUS)
    display = [a.text() for a in menu_of(panel, "Display").actions() if a.text()]
    assert display == ["Show", "Hide", "Colour", "Label", "Hydrogens", "Surface Transparency"]
    select = menu_of(panel, "Select")
    action_in(select, "All").trigger()
    assert panel.session.selected_count() == panel.session.display("1ATP").active.sum()
    assert action_in(menu_of(panel, "Edit"), "Copy").isEnabled()


# figures: export size, background, projection, depth cue


def test_export_image_sizes_keep_the_views_shape(panel, qtbot) -> None:
    from chemsuite.protein.ui.export import CENTIMETRES, INCHES, MAX_PIXELS, ExportImageDialog

    panel.view.resize(400, 300)
    dialog = ExportImageDialog(panel)
    qtbot.addWidget(dialog)
    assert dialog.pixel_size() == (800, 600)  # twice the view
    dialog.width.setValue(1000)
    assert dialog.height.value() == 750 and dialog.primary.text() == "Export 1,000 × 750 px"
    dialog.height.setValue(300)
    assert dialog.width.value() == 400
    dialog.units.set_value(INCHES)  # 400 px at 300 dpi
    assert dialog.width.value() == pytest.approx(400 / 300, abs=0.01)
    dialog.width.setValue(2.0)
    assert dialog.pixel_size() == (600, 450)
    dialog.units.set_value(CENTIMETRES)
    assert dialog.width.value() == pytest.approx(5.08, abs=0.01)
    dialog.dpi.setValue(600)
    assert dialog.pixel_size() == (1200, 900)
    dialog.units.set_value("px")
    dialog.width.setValue(MAX_PIXELS + 1)
    assert "At most" in dialog.status.text() and not dialog.primary.isEnabled()


def test_the_view_menu_sets_background_projection_and_depth_cue(panel) -> None:
    from PySide6.QtGui import QColor

    view = menu_of(panel, "View")
    background = action_in(view, "Background").menu()
    action_in(background, "Black").trigger()
    assert panel.view.background == QColor(0, 0, 0) and panel.background == "black"
    action_in(view, "Orthographic").trigger()
    assert panel.view.camera.orthographic
    action_in(menu_of(panel, "View"), "Depth Cue").trigger()
    assert not panel.view.depth_cue
    view = menu_of(panel, "View")
    assert action_in(view, "Orthographic").isChecked()
    assert not action_in(view, "Depth Cue").isChecked()
    state = panel.state()
    panel.set_background("theme")
    panel.set_orthographic(False)
    panel.set_depth_cue(True)
    assert panel.view.background == panel.palette().base().color()
    panel.set_state(state)
    assert panel.background == "black" and panel.view.camera.orthographic
    assert not panel.view.depth_cue
