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
"""Editing in the Protein Viewer, and its Sequence window: every edit one undoable step."""

from __future__ import annotations

import numpy as np
import pytest
from PySide6.QtCore import QPoint, Qt

from chemsuite.core.table import Table
from chemsuite.protein.look import scene
from chemsuite.protein.structure.model import WATER
from chemsuite.protein.ui.sequence_view import POLYMER

from .test_protein_ui import action_in, add_protein, menu_of


@pytest.fixture
def table() -> Table:
    t = Table()
    add_protein(t, "5GNK", "5gnk_pocket.cif.gz")
    return t


@pytest.fixture
def panel(qtbot, table, monkeypatch):
    from chemsuite.app.shell import graphics
    from chemsuite.protein.ui.panel import ProteinPanel

    monkeypatch.setattr(graphics, "_problem", None)
    widget = ProteinPanel(table)
    qtbot.addWidget(widget)
    widget.resize(900, 600)
    return widget


def waters_row(panel) -> int:
    rows = panel.sequence.strip.rows
    return next(r for r, row in enumerate(rows) if row.kind == WATER)


def test_the_sequence_has_its_own_window_from_the_last_menu(panel, qtbot) -> None:
    titles = [a.text() for a in panel.menu_bar.actions()]
    assert titles[-1] == "Sequence"
    assert panel._sequence_window is None  # made when first asked for
    action_in(menu_of(panel, "Sequence"), "Sequence Viewer").trigger()
    window = panel.sequence_window
    assert window.isVisible() and window.isWindow()
    assert "Sequence" not in [a.text() for a in menu_of(panel, "View").actions()]
    rows = window.strip.rows
    assert [r.label for r in rows if r.kind == POLYMER] == ["5GNK A"]
    assert [r.label for r in rows if r.kind != POLYMER] == [
        "5GNK Ligands", "5GNK Additives", "5GNK Ions", "5GNK Waters"]  # fmt: skip
    waters = rows[waters_row(panel)]
    assert len(waters.letters) == 50 and set(waters.letters) == {"HOH"}
    # the strip is as wide as its longest row from the start (no row cut short)
    assert window.strip.width() == window.strip.sizeHint().width() > window.width()
    action_in(menu_of(panel, "Sequence").actions()[-1].menu(), "Waters").trigger()
    assert all(r.kind != WATER for r in window.strip.rows)
    assert panel.sequence_hidden_kinds() == {WATER}


def test_one_water_is_picked_out_and_deleted_and_undone(panel, table, qtbot) -> None:
    session = panel.session
    shown = session.display("5GNK")
    s = shown.structure
    strip = panel.sequence.strip
    row = waters_row(panel)
    first, second = (int(r) for r in strip.rows[row].sequence.residue_at[:2])
    # the look of the waters left must survive the edit
    shown.styles[s.res_kind[s.residue_of] == WATER] = scene.SPACEFILL
    for k, modifiers in ((0, Qt.KeyboardModifier.NoModifier),
                         (1, Qt.KeyboardModifier.ControlModifier)):  # fmt: skip
        qtbot.mouseClick(strip, Qt.MouseButton.LeftButton, modifiers,
                         strip.box(row, k).center() + QPoint(0, strip.line // 2))  # fmt: skip
    chosen = np.unique(s.residue_of[shown.selected]).tolist()
    assert chosen == [first, second]
    gone = [s.residue_label(first), s.residue_label(second)]
    assert panel.delete_selection() == 2  # a water is one oxygen here
    after = session.display("5GNK")
    labels = {after.structure.residue_label(r) for r in range(after.structure.n_residues)}
    assert not labels & set(gone)
    assert int((after.structure.res_kind == WATER).sum()) == 48
    left = after.structure.res_kind[after.structure.residue_of] == WATER
    assert (after.styles[left] == scene.SPACEFILL).all()
    assert len(panel.sequence.strip.rows[waters_row(panel)].letters) == 48
    assert table.history.undo_label == "Delete 2 atoms from 5GNK"
    assert panel.undo() == "Delete 2 atoms from 5GNK"
    back = session.display("5GNK").structure
    assert back.uid.tolist() == s.uid.tolist()
    assert len(panel.sequence.strip.rows[waters_row(panel)].letters) == 50


def test_a_row_label_selects_the_whole_row(panel, qtbot) -> None:
    strip = panel.sequence.strip
    row = waters_row(panel)
    qtbot.mouseClick(strip, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(4, row * strip.row_height + strip.line + 2))  # fmt: skip
    shown = panel.session.display("5GNK")
    s = shown.structure
    assert np.array_equal(shown.selected, (s.res_kind[s.residue_of] == WATER) & shown.active)


def test_element_and_charge_from_the_edit_menu(panel) -> None:
    session = panel.session
    session.select_text("resn HOH and resi 1103")
    edit = menu_of(panel, "Edit")
    action_in(action_in(edit, "Element").menu(), "S").trigger()
    action_in(action_in(menu_of(panel, "Edit"), "Charge").menu(), "-2").trigger()
    s = session.display("5GNK").structure
    oxygen = s.residue_of == np.flatnonzero(s.res_seq == 1103)[0]
    assert (s.element[oxygen] == 16).all() and (s.charge[oxygen] == -2).all()
    undo = action_in(menu_of(panel, "Edit"), "Undo Set the charge of 1 atom to -2")
    undo.trigger()
    assert (session.display("5GNK").structure.charge[oxygen] == 0).all()


def test_rename_and_renumber_dialogs(panel, qtbot) -> None:
    session = panel.session
    session.select_text("resn HOH and resi 1103")
    dialog = panel.rename_residue()
    qtbot.addWidget(dialog)
    dialog.field.setText("dod")
    dialog.accept()
    s = session.display("5GNK").structure
    assert "A/DOD 1103" in {s.residue_label(r) for r in range(s.n_residues)}
    session.select_text("resn DOD")
    dialog = panel.renumber()
    qtbot.addWidget(dialog)
    dialog.field.setValue(int(s.res_seq[0]))  # taken by the first residue
    dialog.accept()
    assert dialog.isVisible() and "already has that number" in dialog.status.text()
    dialog.field.setValue(5001)
    dialog.accept()
    assert not dialog.isVisible()
    s = session.display("5GNK").structure
    assert "A/DOD 5001" in {s.residue_label(r) for r in range(s.n_residues)}
    session.select_text("all")
    dialog = panel.rename_chain()
    qtbot.addWidget(dialog)
    dialog.field.setText("B")
    dialog.accept()
    assert session.display("5GNK").structure.chain_name.tolist() == ["B"]


def test_delete_residue_from_the_right_click_menu(panel) -> None:
    from chemsuite.protein.ui import menus

    shown = panel.session.display("5GNK")
    s = shown.structure
    lys = int(np.flatnonzero(s.res_name == "LYS")[0])
    atom = s.atoms_of(lys).start
    menus._delete_here(panel, "5GNK", atom, residue=True)
    after = panel.session.display("5GNK").structure
    assert after.n_residues == s.n_residues - 1  # its alternate locations too
    assert "Deleted" in panel.status.text()


def test_deleting_everything_is_refused(panel) -> None:
    panel.session.select_text("all")
    assert panel.delete_selection() == 0
    assert "Remove Protein" in panel.status.text()


def test_the_sequence_window_is_saved_and_restored(panel, table, qtbot) -> None:
    from chemsuite.protein.ui.panel import ProteinPanel

    panel.show_sequence()
    panel.hide_sequence_kind(WATER, True)
    state = panel.state()
    assert state["sequence"]["open"] and state["sequence"]["hidden"] == [WATER]
    again = ProteinPanel(table)
    qtbot.addWidget(again)
    again.set_state(state)
    assert again.sequence_shown()
    assert all(r.kind != WATER for r in again.sequence.strip.rows)
    old = {**state, "rows": [700, 120]}  # a document saved with the strip under the view
    old.pop("sequence")
    third = ProteinPanel(table)
    qtbot.addWidget(third)
    third.set_state(old)
    assert not third.sequence_shown()
