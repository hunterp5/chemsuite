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
from PySide6.QtCore import QPoint, Qt

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
