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
from chemsuite.app.sketcher import canvas as sk_canvas
from chemsuite.app.sketcher import tools as sk_tools
from chemsuite.sketch import convert
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


def test_sketcher_opens_on_a_structure(qapp):
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


def _label_editor(canvas):
    from PySide6.QtWidgets import QLineEdit

    found = [e for e in canvas.findChildren(QLineEdit, "sketch-label-editor") if e.isVisible()]
    return found[-1] if found else None


def _type_label(qapp, canvas, atom, text):
    canvas.edit_label(atom)
    editor = _label_editor(canvas)
    assert editor is not None and editor.isVisible()
    editor.setText(text)
    editor.returnPressed.emit()
    qapp.processEvents()


def test_double_click_types_an_abbreviation_and_it_draws_contracted(dialog, qapp):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("c1ccccc1C"))
    methyl = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 1)
    canvas.set_tool(sk_tools.SelectTool())
    _send(canvas, QEvent.Type.MouseButtonDblClick, canvas.sketch.atoms[methyl].pos)
    editor = _label_editor(canvas)
    assert editor is not None and editor.text() == "C"
    editor.setText("CO2Me")
    editor.returnPressed.emit()
    qapp.processEvents()
    assert to_smiles(canvas.sketch) == Chem.CanonSmiles("COC(=O)c1ccccc1")
    assert canvas.history.undo_name == "Label CO2Me"
    (group,) = canvas.sketch.groups.values()
    # Hidden atoms cannot be picked; the label is picked through its attachment atom.
    for aid in group.atoms - {group.attach}:
        assert canvas.hit_at(canvas.sketch.atoms[aid].pos).id != aid
    assert canvas.hit_at(canvas.sketch.atoms[group.attach].pos).id == group.attach


def test_selecting_or_erasing_a_label_takes_the_whole_group(dialog, qapp):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("c1ccccc1C"))
    methyl = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 1)
    _type_label(qapp, canvas, methyl, "Boc")
    (group,) = canvas.sketch.groups.values()
    canvas.set_tool(sk_tools.SelectTool())
    click(canvas, canvas.sketch.atoms[group.attach].pos)
    assert group.atoms <= canvas.selected_atoms
    hover(canvas, canvas.sketch.atoms[group.attach].pos)
    key(canvas, "", Qt.Key.Key_Delete)
    assert to_smiles(canvas.sketch) == "c1ccccc1"


def test_element_key_on_a_label_replaces_the_abbreviation(dialog, qapp):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("c1ccccc1C"))
    methyl = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 1)
    _type_label(qapp, canvas, methyl, "OTf")
    hover(canvas, canvas.sketch.atoms[methyl].pos)
    key(canvas, "f")
    assert to_smiles(canvas.sketch) == Chem.CanonSmiles("Fc1ccccc1")
    assert not canvas.sketch.groups


def test_bad_label_explains_and_changes_nothing(dialog, qapp):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CC(C)C"))
    centre = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 3)
    _type_label(qapp, canvas, centre, "Ph")
    assert "one bond" in dialog.info.text()
    _type_label(qapp, canvas, centre, "Qq")
    assert "not an element" in dialog.info.text()
    assert to_smiles(canvas.sketch) == "CC(C)C"
    assert canvas.history.undo_name == "Replace drawing"  # the failed labels left no edit


def test_enter_over_an_atom_opens_the_label_editor_and_sets_elements(dialog, qapp):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CC"))
    atom = min(canvas.sketch.atoms)
    hover(canvas, canvas.sketch.atoms[atom].pos)
    key(canvas, "\r", Qt.Key.Key_Return)
    editor = _label_editor(canvas)
    assert editor is not None
    editor.setText("NH2")
    editor.returnPressed.emit()
    assert to_smiles(canvas.sketch) == "CN"


def test_expand_and_contract_keep_the_molecule(dialog, qapp):
    from chemsuite.sketch import abbrev

    canvas = dialog.canvas
    canvas.set_sketch(from_text("c1ccccc1C"))
    methyl = next(a for a in canvas.sketch.atoms if canvas.sketch.degree(a) == 1)
    _type_label(qapp, canvas, methyl, "tBu")
    gid = next(iter(canvas.sketch.groups))
    before = to_smiles(canvas.sketch)
    canvas.edit("Expand", lambda ch: abbrev.set_expanded(ch, gid, True))
    assert not canvas.sketch.hidden_atoms() and to_smiles(canvas.sketch) == before
    canvas.edit("Ungroup", lambda ch: abbrev.ungroup(ch, gid))
    assert not canvas.sketch.groups and to_smiles(canvas.sketch) == before


def test_r_tool_numbers_each_new_rgroup_and_menu_offers_queries(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["atom:R"].click()
    click(canvas, (0.0, 0.0))
    click(canvas, (3.0, 0.0))
    assert sorted(a.label for a in canvas.sketch.atoms.values()) == ["R1", "R2"]
    first = next(aid for aid, a in canvas.sketch.atoms.items() if a.label == "R1")
    click(canvas, (0.0, 0.0))  # clicking an R-group again keeps its number
    assert canvas.sketch.atoms[first].label == "R1"
    dialog.toolbar.buttons["atom:A"].click()
    click(canvas, (0.0, 0.0))
    assert canvas.sketch.atoms[first].label == "A"
    assert " A " in convert.to_molblock(canvas.sketch)


def test_label_editor_accepts_atom_lists(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["atom:N"].click()
    click(canvas, (0.0, 0.0))
    aid = _atom(canvas, "N")
    canvas.edit_label(aid, "![F,Cl]")
    canvas._apply_label(aid, "![F,Cl]")
    assert canvas.sketch.atoms[aid].label == "![F,Cl]"


def _menu(canvas, kind, id_):
    from chemsuite.app.sketcher.canvas import Hit

    return canvas.context_menu(Hit(kind, id_))


def _action(menu, *path):
    for title in path[:-1]:
        menu = next(a.menu() for a in menu.actions() if a.text() == title)
    return next(a for a in menu.actions() if a.text() == path[-1])


def test_stereo_labels_and_menus(qapp):
    from chemsuite.sketch import stereo

    dialog = SketcherDialog("C/C=C/C[C@H](N)O")  # (R), (E)
    dialog.resize(800, 600)
    dialog.show()
    qapp.processEvents()
    canvas = dialog.canvas
    centre, code = next(iter(canvas.stereo.atoms.items()))
    bond = next(iter(canvas.stereo.bonds))
    texts = sorted(t for t, _ in canvas.stereo_label_positions())
    assert texts == sorted([f"({code})", "(E)"])
    flip = stereo.opposite(code)
    _action(_menu(canvas, "atom", centre), "Stereocentre (R/S)", flip).trigger()
    assert canvas.stereo.atoms[centre] == flip
    assert canvas.history.undo_name == f"Make {flip}"
    _action(_menu(canvas, "bond", bond), "Double Bond (E/Z)", "Z").trigger()
    assert canvas.stereo.bonds[bond] == "Z" and canvas.stereo.atoms[centre] == flip
    _action(_menu(canvas, "bond", bond), "Double Bond (E/Z)", "Show Label").trigger()
    assert "(Z)" not in [t for t, _ in canvas.stereo_label_positions()]
    _action(_menu(canvas, "", 0), "Show Stereo Labels (R/S, E/Z)").trigger()
    assert canvas.stereo_label_positions() == []
    _action(_menu(canvas, "", 0), "Show Stereo Labels (R/S, E/Z)").trigger()
    assert len(canvas.stereo_label_positions()) == 2
    dialog.close()


def test_explicit_hydrogens_from_the_atom_menu(qapp):
    dialog = SketcherDialog("CC(N)O")
    dialog.show()
    qapp.processEvents()
    canvas = dialog.canvas
    oxygen = _atom(canvas, "O")
    _action(_menu(canvas, "atom", oxygen), "Show Explicit Hydrogens").trigger()
    hydrogens = [a for a in canvas.sketch.atoms.values() if a.symbol == "H"]
    assert len(hydrogens) == 1 and to_smiles(canvas.sketch) == "CC(N)O"
    hide = _action(_menu(canvas, "atom", oxygen), "Hide Explicit Hydrogens")
    assert hide.isEnabled()
    hide.trigger()
    assert not any(a.symbol == "H" for a in canvas.sketch.atoms.values())
    _action(_menu(canvas, "", 0), "Show All Hydrogens").trigger()
    assert sum(a.symbol == "H" for a in canvas.sketch.atoms.values()) == 7  # CH3, CH, NH2, OH
    dialog.close()


def test_lasso_selects_the_atoms_it_encircles(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CCCCCC"))
    dialog.toolbar.buttons["lasso"].click()
    assert canvas.tool.lasso and dialog.toolbar.buttons["lasso"].isChecked()
    xs = sorted(a.x for a in canvas.sketch.atoms.values())
    cut = (xs[2] + xs[3]) / 2  # loop round the three left-hand atoms
    x0, y0, _x1, y1 = canvas.sketch.bounds()
    loop = [(x0 - 0.5, y0 - 0.5), (cut, y0 - 0.5), (cut, y1 + 0.5), (x0 - 0.5, y1 + 0.5)]
    _send(canvas, QEvent.Type.MouseButtonPress, loop[0])
    for point in [*loop[1:], loop[0]]:
        _send(canvas, QEvent.Type.MouseMove, point, button=Qt.MouseButton.NoButton)
    _send(canvas, QEvent.Type.MouseButtonRelease, loop[0], buttons=Qt.MouseButton.NoButton)
    assert {canvas.sketch.atoms[a].x for a in canvas.selected_atoms} == set(xs[:3])
    assert len(canvas.selected_bonds) == 2
    assert sk_tools.inside_polygon((0.5, 0.5), [(0, 0), (1, 0), (1, 1), (0, 1)])
    assert not sk_tools.inside_polygon((1.5, 0.5), [(0, 0), (1, 0), (1, 1), (0, 1)])


def test_rotate_handle_turns_the_selection_in_15_degree_steps(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CC.O"))
    carbons = [a for a, atom in canvas.sketch.atoms.items() if atom.symbol == "C"]
    oxygen = _atom(canvas, "O")
    still = canvas.sketch.atoms[oxygen]
    canvas.set_tool(sk_tools.SelectTool())
    canvas.select_atoms(carbons)
    before = _angle(canvas, next(iter(canvas.sketch.bonds)))
    at, centre = sk_tools.SelectTool.handle(canvas)
    start = canvas.to_model(at)
    radius = math.dist(start, centre)
    turn = math.radians(50)  # snaps to 45
    base = math.atan2(start[1] - centre[1], start[0] - centre[0])
    end = (centre[0] + radius * math.cos(base + turn), centre[1] + radius * math.sin(base + turn))
    drag(canvas, start, end)
    after = _angle(canvas, next(iter(canvas.sketch.bonds)))
    assert (after - before) % 180 == pytest.approx(45, abs=0.01)
    assert canvas.history.undo_name == "Rotate"
    assert canvas.sketch.atoms[oxygen] == still
    assert set(carbons) <= canvas.selected_atoms


def test_clean_shortcut_tidies_the_drawing_and_undoes_in_one_step(dialog):
    canvas = dialog.canvas
    sketch = from_text("C[C@H](N)C(=O)O")
    for aid, atom in list(sketch.atoms.items()):
        sketch.atoms[aid] = atom.moved(atom.x * 1.3 + 0.2 * (aid % 2), atom.y * 0.8)
    canvas.set_sketch(sketch)
    smiles = to_smiles(canvas.sketch)
    assert convert.bond_length_spread(canvas.sketch) > 0.1
    key(canvas, "K", Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier
        | Qt.KeyboardModifier.ShiftModifier)  # fmt: skip
    assert canvas.history.undo_name == "Clean Structure"
    assert convert.bond_length_spread(canvas.sketch) < 0.02
    assert to_smiles(canvas.sketch) == smiles
    canvas.history.undo()
    assert convert.bond_length_spread(canvas.sketch) > 0.1


def test_salts_menu_adds_a_solvent_beside_the_drawing(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("c1ccccc1N"))
    menu = canvas.context_menu(sk_canvas.NOTHING)
    _action(menu, "Add Salt or Solvent", "Polar aprotic", "N,N-dimethylformamide (DMF)").trigger()
    assert to_smiles(canvas.sketch) == Chem.CanonSmiles("Nc1ccccc1.CN(C)C=O")
    assert len(canvas.selected_atoms) == 5
    assert canvas.history.undo_name == "Add N,N-dimethylformamide"
    texts = [a.text() for a in canvas.context_menu(sk_canvas.NOTHING).actions()]
    assert "Clean Selection\t" + sk_canvas._shortcut("Ctrl+Shift+K") in texts


def test_either_double_bond_from_the_menu_draws_crossed(dialog):
    from chemsuite.sketch import stereo

    canvas = dialog.canvas
    canvas.set_sketch(from_text("C/C=C/C"))
    bond = next(iter(canvas.stereo.bonds))
    _action(_menu(canvas, "bond", bond), "Double Bond (E/Z)", "Unknown (crossed bond)").trigger()
    assert bond in canvas.stereo.unknown_bonds and not canvas.stereo.bonds
    assert canvas.history.undo_name == "Make E/Z unknown"
    lines = [ln for ln in canvas.scene.lines if ln.bond == bond]
    assert len(lines) == 2
    (a1, b1), (a2, b2) = (lines[0].a, lines[0].b), (lines[1].a, lines[1].b)
    from chemsuite.sketch.layout import _crosses

    assert _crosses(a1, b1, a2, b2)  # the two lines of the bond cross
    menu = _menu(canvas, "bond", bond)
    assert _action(menu, "Double Bond (E/Z)", "Unknown (crossed bond)").isChecked()
    _action(menu, "Double Bond (E/Z)", "E").trigger()
    assert canvas.stereo.bonds[bond] == "E"
    assert stereo.UNKNOWN == "either"


def test_bad_wedge_is_marked_and_explained_without_blocking_ok(dialog):
    from dataclasses import replace

    from chemsuite.sketch import stereo
    from chemsuite.sketch.model import Stereo

    canvas = dialog.canvas
    sketch = from_text("CCC")
    first = min(sketch.bonds)
    sketch.bonds[first] = replace(sketch.bonds[first], stereo=Stereo.WEDGE)
    canvas.set_sketch(sketch)
    assert canvas.wedge_problems == {first: stereo.NO_CENTRE}
    assert dialog.info.isVisible() and "wedge" in dialog.info.text()
    assert dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()


def test_enhanced_stereo_menu_labels_and_round_trip(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("C[C@H](N)C[C@@H](O)C"))
    first, second = sorted(canvas.stereo.atoms)
    and_group = "AND Group (racemic: this and its mirror image)"
    _action(_menu(canvas, "atom", first), "Enhanced Stereo", and_group, "New Group (&1)").trigger()
    assert canvas.sketch.atoms[first].enhanced == "and1"
    assert canvas.history.undo_name == "Enhanced stereo: and1"
    canvas.select_atoms([second])
    menu = canvas.popup_enhanced()  # Alt+E on the selection
    _action(menu, "OR Group (one of this or its mirror image, not known which)",
            "New Group (or1)").trigger()  # fmt: skip
    menu.close()
    assert canvas.sketch.atoms[second].enhanced == "or1"
    texts = [t for t, _ in canvas.stereo_label_positions()]
    assert any(t.endswith(" &1") for t in texts) and any(t.endswith(" or1") for t in texts)
    dialog.accept()
    assert "MDLV30/STERAC1" in dialog.molblock and "MDLV30/STEREL1" in dialog.molblock


def test_bond_query_menu_and_zero_key(dialog):
    from chemsuite.sketch.model import BOND_QUERIES, TOPOLOGY

    canvas = dialog.canvas
    canvas.set_sketch(from_text("CCN"))
    first, second = sorted(canvas.sketch.bonds)
    menu = _menu(canvas, "bond", first)
    _action(menu, "Bond Query", BOND_QUERIES[Order.SINGLE_OR_DOUBLE]).trigger()
    assert canvas.sketch.bonds[first].order is Order.SINGLE_OR_DOUBLE
    _action(_menu(canvas, "bond", first), "Bond Query", TOPOLOGY[1]).trigger()
    assert canvas.sketch.bonds[first].topology == 1
    assert [t for _, t in canvas.scene.tags] == ["S/D Rng"]
    _action(_menu(canvas, "bond", first), "Bond Query", "Back to Single").trigger()
    assert canvas.sketch.bonds[first].order is Order.SINGLE
    bond = canvas.sketch.bonds[second]
    a, b = canvas.sketch.atoms[bond.begin].pos, canvas.sketch.atoms[bond.end].pos
    hover(canvas, ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2))
    key(canvas, "0", Qt.Key.Key_0)
    assert canvas.sketch.bonds[second].order is Order.ANY
    assert dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()


def test_attachment_point_menu(dialog):
    canvas = dialog.canvas
    canvas.set_sketch(from_text("CCC"))
    end = min(canvas.sketch.atoms)
    _action(_menu(canvas, "atom", end), "Attachment Point", "First").trigger()
    assert canvas.sketch.atoms[end].attach == 1
    assert canvas.history.undo_name == "Attachment point: first"
    assert len(canvas.scene.stubs) == 1
    dialog.accept()
    assert "M  APO" in dialog.molblock


def test_drawn_query_reaches_structure_search_as_smarts(qapp, monkeypatch):
    from dataclasses import replace

    from chemsuite.app.sketcher import SketcherDialog
    from chemsuite.app.structure_search import QueryInput
    from chemsuite.sketch.model import ATOM_LIST

    sketch = from_text("Clc1ccccc1")
    halogen = convert.atom_order(sketch)[0]
    sketch.atoms[halogen] = replace(sketch.atoms[halogen], symbol=ATOM_LIST, atom_list=("Cl", "Br"))
    block = convert.to_molblock(sketch)
    seen: list[str] = []
    real_init = SketcherDialog.__init__

    def init(self, text="", *args, **kwargs):
        seen.append(text)
        real_init(self, block, *args, **kwargs)

    monkeypatch.setattr(SketcherDialog, "__init__", init)

    def accept(self):
        self.molblock = block
        return True

    monkeypatch.setattr(SketcherDialog, "exec", accept)
    field = QueryInput()
    field.draw()
    assert field.text() == "[#17,#35]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    field.draw()  # drawing again starts from the molfile, not the SMARTS
    assert Chem.MolFromMolBlock(seen[1]) is not None and "M  END" in seen[1]


# -- reactions (REACTIONS_PLAN R2) ----------------------------------------------------------------


def _paste_text(canvas, text):
    from PySide6.QtCore import QMimeData

    mime = QMimeData()
    mime.setText(text)
    QGuiApplication.clipboard().setMimeData(mime)
    assert canvas.paste()


def test_draw_a_reaction_with_the_arrow_plus_and_map_tools(dialog):
    from chemsuite.sketch import reaction

    canvas = dialog.canvas
    canvas.set_sketch(convert.placed(from_text("CO"), (0, 0)), "Start")
    dialog.toolbar.buttons["atom:C"].click()  # a second reactant, then a plus between them
    click(canvas, (3.0, 0.0))
    dialog.toolbar.buttons["rxn-plus"].click()
    click(canvas, (1.8, 0.0))
    assert len(canvas.sketch.pluses) == 1
    dialog.toolbar.buttons["arrow"].click()
    drag(canvas, (4.5, 0.0), (7.5, 0.2))  # snaps to horizontal
    (arrow,) = canvas.sketch.arrows.values()
    assert arrow.y1 == pytest.approx(arrow.y0) and arrow.x1 > arrow.x0
    dialog.toolbar.buttons["atom:O"].click()
    click(canvas, (9.0, 0.0))
    found = reaction.components(canvas.sketch)
    assert len(found.reactants) == 2 and len(found.products) == 1
    oxygen = next(a for a, at in canvas.sketch.atoms.items() if at.symbol == "O" and at.x < 3)
    product = next(a for a, at in canvas.sketch.atoms.items() if at.x > 8)
    dialog.toolbar.buttons["map"].click()
    drag(canvas, canvas.sketch.atoms[oxygen].pos, canvas.sketch.atoms[product].pos)
    assert canvas.sketch.atoms[oxygen].mapping == canvas.sketch.atoms[product].mapping == 1
    assert len(canvas.scene.maps) == 2 and canvas.scene.strokes
    assert canvas.history.undo_name == "Map atoms"


def test_arrows_select_move_and_erase(dialog):
    canvas = dialog.canvas
    dialog.toolbar.buttons["arrow"].click()
    click(canvas, (0.0, 0.0))
    (rid,) = canvas.sketch.arrows
    before = canvas.sketch.arrows[rid]
    dialog.toolbar.buttons["select"].click()
    drag(canvas, (0.0, 0.0), (0.0, 2.0))  # grab the arrow's middle and move it up
    moved = canvas.sketch.arrows[rid]
    assert moved.y0 == pytest.approx(before.y0 + 2.0) and rid in canvas.selected_objects
    key(canvas, "", Qt.Key.Key_Delete)
    assert not canvas.sketch.arrows


def test_arrow_menu_changes_its_kind(dialog):
    from chemsuite.sketch.model import ArrowKind

    canvas = dialog.canvas
    dialog.toolbar.buttons["arrow"].click()
    click(canvas, (0.0, 0.0))
    (rid,) = canvas.sketch.arrows
    menu = canvas.context_menu(sk_canvas.Hit("arrow", rid))
    arrows = next(a.menu() for a in menu.actions() if a.text() == "Arrow")
    next(a for a in arrows.actions() if a.text() == "Equilibrium Arrow").trigger()
    assert canvas.sketch.arrows[rid].kind is ArrowKind.EQUILIBRIUM
    menu = canvas.context_menu(sk_canvas.Hit("arrow", rid))
    next(a for a in menu.actions() if a.text() == "Reverse Arrow").trigger()
    arrow = canvas.sketch.arrows[rid]
    assert arrow.x0 > arrow.x1


def test_reaction_copy_and_paste(dialog):
    canvas = dialog.canvas
    text = "[CH3:1][OH:2].[CH3:3][Cl:4]>CN(C)C=O>[CH3:1][O:2][CH3:3]"
    _paste_text(canvas, text)
    assert len(canvas.sketch.arrows) == 1 and len(canvas.sketch.pluses) == 1
    assert canvas.selected_objects  # the pasted arrow and plus come selected
    canvas.clear_selection()
    assert canvas.copy()
    mime = QGuiApplication.clipboard().mimeData()
    assert bytes(mime.data("chemical/x-mdl-rxnfile").data()).startswith(b"$RXN")
    from rdkit.Chem import rdChemReactions

    want = rdChemReactions.ReactionToSmiles(
        rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
    )
    assert mime.text() == want


def test_dialog_accepts_reactions_only_when_asked(qapp):
    text = "CCO>>CC=O"
    plain = SketcherDialog(text)
    plain.accept()
    assert plain.result() == 0 and "not a reaction" in plain.info.text()
    allowed = SketcherDialog(text, reactions=True)
    allowed.accept()
    assert allowed.rxn.startswith("$RXN") and allowed.molblock == ""
    for dlg in (plain, allowed):
        dlg.deleteLater()


def test_template_library_picks_a_drawing_tool(dialog, qapp) -> None:
    from chemsuite.app.sketcher.library import TemplateLibrary

    opened: list = []
    dialog.canvas.library_requested.disconnect()
    dialog.canvas.library_requested.connect(lambda: opened.append(dialog.open_library()))
    key(dialog.canvas, "T", Qt.Key.Key_T, Qt.KeyboardModifier.ShiftModifier)  # Shift+T
    assert opened and isinstance(opened[0], TemplateLibrary)
    window = opened[0]
    assert len(window.names()) == window.tiles.count() >= 60
    window.search.setText("morph")
    assert window.names() == ["Morpholine", "Thiomorpholine"]
    window.search.clear()
    window.groups.setCurrentRow(window.groups.count() - 1)  # the last group: bridged and spiro
    assert "Adamantane" in window.names() and "Pyridine" not in window.names()
    window.tiles.setCurrentRow(window.names().index("Adamantane"))
    window.accept()
    assert isinstance(dialog.canvas.tool, sk_tools.TemplateTool)
    assert dialog.canvas.tool.template.name == "Adamantane"
    assert dialog.toolbar.buttons["library"].isChecked()  # the ring bar shows it in use
    click(dialog.canvas, (0.0, 0.0))
    assert Chem.MolToSmiles(convert.to_rdkit(dialog.canvas.sketch)) == Chem.MolToSmiles(
        Chem.MolFromSmiles("C1C2CC3CC1CC(C2)C3")
    )
