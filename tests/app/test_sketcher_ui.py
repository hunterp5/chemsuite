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

"""The sketcher driven by synthetic mouse and key events (interaction-spec §1–3)."""

from __future__ import annotations

import math

import pytest
from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QGuiApplication, QKeyEvent, QMouseEvent
from rdkit import Chem

from chemsuite.app.sketcher import SketcherDialog
from chemsuite.app.sketcher import tools as sk_tools
from chemsuite.sketch.convert import from_text, to_smiles
from chemsuite.sketch.model import Order

from .conftest import gui

pytestmark = gui

NO_MODS = Qt.KeyboardModifier.NoModifier
LEFT = Qt.MouseButton.LeftButton
ALANINE = Chem.MolToSmiles(Chem.MolFromSmiles("C[C@@H](N)C(=O)O"))  # L-alanine


@pytest.fixture
def dialog(qapp):
    dlg = SketcherDialog("")
    dlg.resize(800, 600)
    dlg.show()
    qapp.processEvents()
    yield dlg
    dlg.close()
    dlg.deleteLater()


def _send(canvas, kind, point, buttons=LEFT, mods=NO_MODS, button=LEFT):
    screen = canvas.to_screen(point)
    event = QMouseEvent(kind, screen, canvas.mapToGlobal(screen), button, buttons, mods)
    QGuiApplication.sendEvent(canvas, event)


def click(canvas, point, mods=NO_MODS):
    _send(canvas, QEvent.Type.MouseButtonPress, point, mods=mods)
    _send(canvas, QEvent.Type.MouseButtonRelease, point, buttons=Qt.MouseButton.NoButton, mods=mods)


def drag(canvas, a, b, mods=NO_MODS, steps=6):
    _send(canvas, QEvent.Type.MouseButtonPress, a, mods=mods)
    for i in range(1, steps + 1):
        t = i / steps
        p = (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)
        _send(canvas, QEvent.Type.MouseMove, p, mods=mods, button=Qt.MouseButton.NoButton)
    _send(canvas, QEvent.Type.MouseButtonRelease, b, buttons=Qt.MouseButton.NoButton, mods=mods)


def hover(canvas, point):
    _send(
        canvas,
        QEvent.Type.MouseMove,
        point,
        buttons=Qt.MouseButton.NoButton,
        button=Qt.MouseButton.NoButton,
    )


def key(canvas, text, k=None, mods=NO_MODS):
    k = k if k is not None else Qt.Key(ord(text.upper())) if text else Qt.Key.Key_unknown
    QGuiApplication.sendEvent(canvas, QKeyEvent(QEvent.Type.KeyPress, k, mods, text))


def _angle(canvas, bid):
    bond = canvas.sketch.bonds[bid]
    a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 360


def _atom(canvas, symbol):
    return next(aid for aid, a in canvas.sketch.atoms.items() if a.symbol == symbol)


def test_bond_tool_click_empty_then_atoms_and_bond(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["single"].click()
    assert isinstance(canvas.tool, sk_tools.BondTool)
    click(canvas, (0.0, 0.0))
    (bid,) = canvas.sketch.bonds
    assert _angle(canvas, bid) == pytest.approx(30)  # Ketcher: single bonds rise 30 degrees
    end = canvas.sketch.bonds[bid].end
    click(canvas, canvas.sketch.atoms[end].pos)
    click(canvas, canvas.sketch.atoms[max(canvas.sketch.atoms)].pos)
    assert to_smiles(canvas.sketch) == "CCCC"
    # Clicking a bond with the single-bond tool cycles single -> double -> triple.
    mid = canvas.sketch.atoms[canvas.sketch.bonds[bid].begin].pos
    mid = (
        (mid[0] + canvas.sketch.atoms[end].pos[0]) / 2,
        (mid[1] + canvas.sketch.atoms[end].pos[1]) / 2,
    )
    click(canvas, mid)
    assert canvas.sketch.bonds[bid].order is Order.DOUBLE
    click(canvas, mid)
    assert canvas.sketch.bonds[bid].order is Order.TRIPLE
    assert "C4H" in dialog.info.text()


def test_dragged_bond_snaps_to_15_degrees_and_ctrl_frees_it(dialog):
    canvas = dialog.canvas
    canvas.set_tool(sk_tools.BondTool())
    tip = (math.cos(math.radians(50)), math.sin(math.radians(50)))
    drag(canvas, (0.0, 0.0), tip)
    (bid,) = canvas.sketch.bonds
    assert _angle(canvas, bid) == pytest.approx(45, abs=0.01)
    assert not canvas.history.redo_name
    canvas.history.undo()
    assert canvas.sketch.is_empty()
    drag(canvas, (0.0, 0.0), tip, mods=Qt.KeyboardModifier.ControlModifier)
    (bid,) = canvas.sketch.bonds
    assert _angle(canvas, bid) == pytest.approx(50, abs=0.5)


def test_drag_from_atom_onto_atom_closes_a_ring(dialog):
    canvas = dialog.canvas
    canvas.set_tool(sk_tools.ChainTool())
    drag(canvas, (0.0, 0.0), (4.4, 0.0))  # 5 bonds
    assert len(canvas.sketch.bonds) == 5
    canvas.set_tool(sk_tools.BondTool())
    ids = sorted(canvas.sketch.atoms)
    first, last = canvas.sketch.atoms[ids[0]].pos, canvas.sketch.atoms[ids[-1]].pos
    drag(canvas, first, last)
    assert len(canvas.sketch.atoms) == 6 and len(canvas.sketch.bonds) == 6
    assert to_smiles(canvas.sketch) == "C1CCCCC1"


def test_chain_adds_one_bond_per_bond_length(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["chain"].click()
    drag(canvas, (0.0, 0.0), (3.4, 0.2))
    assert len(canvas.sketch.bonds) == 4
    assert canvas.history.undo_name == "Chain"
    canvas.history.undo()
    assert canvas.sketch.is_empty()


def test_hover_keys_edit_what_is_under_the_pointer(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CCC"))
    end = min(canvas.sketch.atoms)
    hover(canvas, canvas.sketch.atoms[end].pos)
    key(canvas, "o")
    assert to_smiles(canvas.sketch) == "CCO"
    o = _atom(canvas, "O")
    key(canvas, "=", Qt.Key.Key_Equal)
    assert canvas.sketch.atoms[o].charge == 1
    key(canvas, "-", Qt.Key.Key_Minus)
    key(canvas, "-", Qt.Key.Key_Minus)
    assert to_smiles(canvas.sketch) == "CC[O-]"
    key(canvas, "s")
    assert to_smiles(canvas.sketch) == "CC[S-]"
    # Over a bond, digits set the order.
    bid = next(iter(canvas.sketch.bonds))
    bond = canvas.sketch.bonds[bid]
    a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
    hover(canvas, ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2))
    key(canvas, "2")
    assert canvas.sketch.bonds[bid].order is Order.DOUBLE
    # Over an atom, a digit adds a bond.
    hover(canvas, canvas.sketch.atoms[_atom(canvas, "S")].pos)
    before = len(canvas.sketch.bonds)
    key(canvas, "1")
    assert len(canvas.sketch.bonds) == before + 1


def test_keys_without_hover_pick_tools_and_the_toolbar_follows(dialog):
    canvas = dialog.canvas
    hover(canvas, (50.0, 50.0))
    key(canvas, "n")
    assert isinstance(canvas.tool, sk_tools.AtomTool) and canvas.tool.label == "N"
    assert dialog.toolbar.buttons["atom:N"].isChecked()
    key(canvas, "2")
    assert dialog.toolbar.buttons["double"].isChecked()
    key(canvas, "", Qt.Key.Key_Escape)
    assert dialog.toolbar.buttons["select"].isChecked()


def test_atom_tool_sets_elements_and_places_atoms(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["atom:N"].click()
    click(canvas, (0.0, 0.0))
    assert to_smiles(canvas.sketch) == "N"
    drag(canvas, (0.0, 0.0), (1.0, 0.0))
    dialog.toolbar.buttons["atom:O"].click()
    click(canvas, (0.0, 0.0))
    assert to_smiles(canvas.sketch) == "NO"


def test_erase_stroke_is_one_undo_step(dialog):
    canvas = dialog.canvas
    canvas.set_tool(sk_tools.ChainTool())
    drag(canvas, (0.0, 0.0), (3.6, 0.0))
    dialog.toolbar.buttons["erase"].click()
    pts = [canvas.sketch.atoms[a].pos for a in sorted(canvas.sketch.atoms)]
    drag(canvas, pts[-1], pts[-2], steps=3)
    assert len(canvas.sketch.bonds) < 4
    canvas.history.undo()
    assert len(canvas.sketch.bonds) == 4


def test_select_box_and_move_merges_onto_atoms_then_undo(dialog):
    canvas = dialog.canvas
    canvas.set_tool(sk_tools.ChainTool())
    drag(canvas, (0.0, 0.0), (4.4, 0.0))  # six atoms, five bonds
    canvas.set_tool(sk_tools.SelectTool())
    ids = sorted(canvas.sketch.atoms)
    last = canvas.sketch.atoms[ids[-1]].pos
    first = canvas.sketch.atoms[ids[0]].pos
    click(canvas, last)
    assert canvas.selected_atoms == {ids[-1]}
    drag(canvas, last, first)
    assert to_smiles(canvas.sketch) == "C1CCCC1"
    key(canvas, "z", Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    assert to_smiles(canvas.sketch) == "CCCCCC"
    # A box selects what it encloses; Delete erases the selection.
    drag(canvas, (-0.5, -2.0), (2.0, 2.0))
    assert canvas.selected_atoms and len(canvas.selected_atoms) < len(canvas.sketch.atoms)
    key(canvas, "", Qt.Key.Key_Delete)
    assert len(canvas.sketch.atoms) < 6


def test_copy_and_paste_round_trip(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("C[C@@H](N)C(=O)O"))
    assert canvas.copy()
    mime = QGuiApplication.clipboard().mimeData()
    assert "chemical/x-mdl-molfile" in mime.formats()
    block = bytes(mime.data("chemical/x-mdl-molfile")).decode()
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(block)) == ALANINE
    canvas.set_sketch(from_text("O"))
    assert canvas.paste()
    assert to_smiles(canvas.sketch) == Chem.MolToSmiles(Chem.MolFromSmiles(f"{ALANINE}.O"))
    canvas.history.undo()
    assert to_smiles(canvas.sketch) == "O"


def test_valence_problem_disables_ok_and_explains(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CC(C)(C)C"))
    ok = dialog.buttons.button(dialog.buttons.StandardButton.Ok)
    assert ok.isEnabled()
    centre = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 4)
    hover(canvas, canvas.sketch.atoms[centre].pos)
    key(canvas, "1")
    assert not ok.isEnabled()
    assert "valence problem" in dialog.info.text()
    assert "Too many bonds" in canvas.report.problems[centre]
    key(canvas, "z", Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    assert ok.isEnabled()


def test_escape_cancels_a_drag_preview(dialog):
    canvas = dialog.canvas
    canvas.set_tool(sk_tools.ChainTool())
    _send(canvas, QEvent.Type.MouseButtonPress, (0.0, 0.0))
    for x in (0.5, 1.5, 2.5):
        _send(canvas, QEvent.Type.MouseMove, (x, 0.0), button=Qt.MouseButton.NoButton)
    assert len(canvas.sketch.bonds) == 3
    canvas.cancel_preview()
    assert canvas.sketch.is_empty() and not canvas.history.can_undo


def test_ok_returns_a_molfile_and_the_structure_editor_uses_it(qapp, monkeypatch):
    from chemsuite.app.structure_editor import StructureEditor

    def fake_exec(self):
        self.canvas.history.replace(from_text("c1ccccc1O"))
        self.accept()
        return True

    monkeypatch.setattr(SketcherDialog, "exec", fake_exec)
    editor = StructureEditor("CCO")
    editor.draw()
    assert editor.current.kind == "molfile"
    assert Chem.MolToSmiles(editor.current.mol) == "Oc1ccccc1"
    editor.close()


def test_sketcher_opens_on_the_editor_structure(qapp):
    dlg = SketcherDialog("CC(=O)Oc1ccccc1C(=O)O")
    assert to_smiles(dlg.canvas.sketch) == "CC(=O)Oc1ccccc1C(=O)O"
    assert not dlg.canvas.history.can_undo  # opening is not an edit
    dlg.accept()
    assert "M  END" in dlg.molblock
    dlg.deleteLater()


def test_dragging_a_selected_bond_moves_the_whole_selection(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CCCC"))
    canvas.select_atoms(canvas.sketch.atoms)
    bid = next(iter(canvas.sketch.bonds))
    bond = canvas.sketch.bonds[bid]
    a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
    mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    before = {aid: atom.pos for aid, atom in canvas.sketch.atoms.items()}
    drag(canvas, mid, (mid[0] + 3.0, mid[1]))
    for aid, (x, y) in before.items():
        assert canvas.sketch.atoms[aid].pos == pytest.approx((x + 3.0, y))
    assert canvas.history.undo_name == "Move"


def test_ring_tool_places_fuses_and_spiro(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["template:Benzene"].click()
    assert isinstance(canvas.tool, sk_tools.TemplateTool)
    click(canvas, (0.0, 0.0))
    assert to_smiles(canvas.sketch) == "c1ccccc1"
    bid = min(canvas.sketch.bonds)
    bond = canvas.sketch.bonds[bid]
    a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
    click(canvas, ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2))
    assert to_smiles(canvas.sketch) == Chem.CanonSmiles("c1ccc2ccccc2c1")
    assert canvas.history.undo_name == "Fuse Benzene"


def test_ring_drag_from_atom_past_one_bond_attaches_by_a_bond(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CC"))
    atom = max(canvas.sketch.atoms, key=lambda a: canvas.sketch.atoms[a].x)
    p = canvas.sketch.atoms[atom].pos
    dialog.toolbar.buttons["template:Cyclohexane"].click()
    drag(canvas, p, (p[0] + 2.0, p[1]))
    assert to_smiles(canvas.sketch) == Chem.CanonSmiles("CCC1CCCCC1")


def test_t_cycles_the_rings(dialog):
    canvas = dialog.canvas
    hover(canvas, (50.0, 50.0))
    key(canvas, "t")
    assert dialog.toolbar.buttons["template:Benzene"].isChecked()
    key(canvas, "t")
    assert dialog.toolbar.buttons["template:Cyclopropane"].isChecked()
