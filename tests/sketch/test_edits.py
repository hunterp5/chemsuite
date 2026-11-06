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

"""Undoable edits and the editing rules taken from Ketcher at 2b031cc2."""

from __future__ import annotations

import math

import pytest

from chemsuite.sketch import edits
from chemsuite.sketch.convert import from_text, to_smiles
from chemsuite.sketch.edits import History
from chemsuite.sketch.model import Atom, Order, Sketch, Stereo
from tests.sketch.helpers import canonical


def _snapshot(sketch: Sketch):
    return (dict(sketch.atoms), dict(sketch.bonds))


def _history(smiles: str) -> History:
    return History(from_text(smiles))


def _bond(history: History, a_symbol: str, b_symbol: str) -> int:
    for bid, bond in history.sketch.bonds.items():
        pair = {history.sketch.atoms[bond.begin].symbol, history.sketch.atoms[bond.end].symbol}
        if pair == {a_symbol, b_symbol}:
            return bid
    raise LookupError


# -- undo and redo ---------------------------------------------------------------------------


def test_every_edit_undoes_exactly_and_redoes_exactly():
    history = _history("CC(=O)Oc1ccccc1C(=O)O")
    states = [_snapshot(history.sketch)]
    steps = [
        lambda ch: edits.bond_from_atom(ch, next(iter(ch.sketch.atoms))),
        lambda ch: edits.set_element(ch, next(iter(ch.sketch.atoms)), "N"),
        lambda ch: edits.add_charge(ch, next(iter(ch.sketch.atoms)), +1),
        lambda ch: edits.chain(ch, (10.0, 0.0), 0.0, 3),
        lambda ch: edits.erase(ch, bonds=[next(iter(ch.sketch.bonds))]),
        lambda ch: edits.flip(ch, list(ch.sketch.atoms), horizontal=True),
        lambda ch: edits.rotate(ch, list(ch.sketch.atoms), (0, 0), 0.3),
    ]
    for step in steps:
        change = history.change("step")
        step(change)
        assert history.commit(change) is not None
        states.append(_snapshot(history.sketch))
    for state in reversed(states[:-1]):
        assert history.undo()
        assert _snapshot(history.sketch) == state
    assert not history.undo()
    for state in states[1:]:
        assert history.redo()
        assert _snapshot(history.sketch) == state
    assert not history.redo()


def test_a_new_edit_clears_redo_and_empty_changes_are_dropped():
    history = _history("CC")
    change = history.change("noop")
    assert history.commit(change) is None and not history.can_undo
    change = history.change("N")
    edits.set_element(change, next(iter(history.sketch.atoms)), "N")
    history.commit(change)
    history.undo()
    assert history.can_redo
    change = history.change("O")
    edits.set_element(change, next(iter(history.sketch.atoms)), "O")
    history.commit(change)
    assert not history.can_redo and history.undo_name == "O"


def test_discard_rolls_back_a_preview():
    history = _history("CC")
    before = _snapshot(history.sketch)
    change = history.change("drag")
    edits.chain(change, next(iter(history.sketch.atoms)), 0.0, 5)
    history.discard(change)
    assert _snapshot(history.sketch) == before and not history.can_undo


def test_listeners_hear_commits_undo_and_redo():
    history = _history("C")
    heard = []
    history.listeners.append(lambda: heard.append(1))
    change = history.change("N")
    edits.set_element(change, next(iter(history.sketch.atoms)), "N")
    history.commit(change)
    history.undo()
    history.redo()
    assert len(heard) == 3


def test_replace_swaps_the_drawing_in_one_undoable_step():
    history = _history("CCO")
    history.replace(from_text("c1ccccc1"), "Paste")
    assert to_smiles(history.sketch) == "c1ccccc1"
    history.undo()
    assert to_smiles(history.sketch) == "CCO"


# -- bondChangingAction ----------------------------------------------------------------------


def test_single_bond_tool_cycles_plain_bonds():
    history = _history("CC")
    (bid,) = history.sketch.bonds
    seen = []
    for _ in range(4):
        change = history.change("click")
        edits.click_bond(change, bid, Order.SINGLE)
        history.commit(change)
        seen.append(history.sketch.bonds[bid].order)
    assert seen == [Order.DOUBLE, Order.TRIPLE, Order.SINGLE, Order.DOUBLE]


def test_other_bond_tools_set_their_type():
    history = _history("CC")
    (bid,) = history.sketch.bonds
    for order in (Order.TRIPLE, Order.TRIPLE, Order.DOUBLE):
        change = history.change("click")
        edits.click_bond(change, bid, order)
        history.commit(change)
        assert history.sketch.bonds[bid].order is order  # no cycling for these tools


def test_same_wedge_tool_flips_the_wedge():
    history = _history("CC")
    (bid,) = history.sketch.bonds
    first = history.sketch.bonds[bid]
    change = history.change("wedge")
    edits.click_bond(change, bid, Order.SINGLE, Stereo.WEDGE)
    assert history.sketch.bonds[bid].stereo is Stereo.WEDGE
    assert history.sketch.bonds[bid].begin == first.begin
    edits.click_bond(change, bid, Order.SINGLE, Stereo.WEDGE)
    assert history.sketch.bonds[bid].begin == first.end  # narrow end moved
    edits.click_bond(change, bid, Order.SINGLE, Stereo.HASH)
    assert history.sketch.bonds[bid].stereo is Stereo.HASH
    edits.click_bond(change, bid, Order.SINGLE)  # single tool on a wedge: plain single
    assert history.sketch.bonds[bid].stereo is Stereo.NONE
    assert history.sketch.bonds[bid].order is Order.SINGLE


def test_hover_key_sets_the_order_and_drops_a_wedge():
    history = _history("CC")
    (bid,) = history.sketch.bonds
    change = history.change("x")
    edits.click_bond(change, bid, Order.SINGLE, Stereo.WEDGE)
    edits.set_bond_order(change, bid, Order.DOUBLE)
    assert history.sketch.bonds[bid].order is Order.DOUBLE
    assert history.sketch.bonds[bid].stereo is Stereo.NONE


# -- atoms -----------------------------------------------------------------------------------


def test_set_element_keeps_charge_and_handles_deuterium():
    history = _history("C[NH3+]")
    n = next(a for a, atom in history.sketch.atoms.items() if atom.symbol == "N")
    change = history.change("P")
    edits.set_element(change, n, "P")
    history.commit(change)
    assert to_smiles(history.sketch) == canonical("C[PH3+]")
    change = history.change("D")
    c = next(a for a, atom in history.sketch.atoms.items() if atom.symbol == "C")
    edits.set_element(change, c, "D")
    history.commit(change)
    assert history.sketch.atoms[c].symbol == "H" and history.sketch.atoms[c].isotope == 2


def test_charge_keys_step_the_charge():
    history = _history("CN")
    n = next(a for a, atom in history.sketch.atoms.items() if atom.symbol == "N")
    change = history.change("+")
    edits.add_charge(change, n, +1)
    edits.add_charge(change, n, +1)
    edits.add_charge(change, n, -1)
    history.commit(change)
    assert to_smiles(history.sketch) == "C[NH3+]"


def test_atom_on_canvas():
    history = History()
    change = history.change("atom")
    edits.atom_on_canvas(change, (1.0, 2.0), "Cl")
    history.commit(change)
    assert to_smiles(history.sketch) == "Cl"


# -- erase -----------------------------------------------------------------------------------


def test_erasing_a_bond_removes_atoms_it_leaves_alone():
    """Ketcher ``fromBondDeletion``: an end atom left with no bonds goes too."""
    history = _history("CCO")
    change = history.change("erase")
    edits.erase(change, bonds=[_bond(history, "C", "O")])
    history.commit(change)
    assert to_smiles(history.sketch) == "CC"


def test_erasing_the_only_bond_removes_both_atoms():
    history = _history("CC")
    change = history.change("erase")
    edits.erase(change, bonds=list(history.sketch.bonds))
    assert history.sketch.is_empty()


def test_erasing_an_atom_keeps_its_neighbours():
    history = _history("CC(C)O")
    o = next(a for a, atom in history.sketch.atoms.items() if atom.symbol == "O")
    change = history.change("erase")
    edits.erase(change, atoms=[o])
    assert to_smiles(history.sketch) == "CCC"


# -- move and merge --------------------------------------------------------------------------


def test_moving_an_atom_onto_another_merges_them():
    """Ketcher ``atomMerge``: dragging an end atom onto another atom closes a ring."""
    history = History()
    change = history.change("chain")
    atoms = edits.chain(change, (0.0, 0.0), 0.0, 5)
    history.commit(change)
    first, last = history.sketch.atoms[atoms[0]], history.sketch.atoms[atoms[-1]]
    change = history.change("move")
    merged = edits.move(change, [atoms[-1]], first.x - last.x, first.y - last.y)
    history.commit(change)
    assert merged == {atoms[-1]: atoms[0]}
    assert to_smiles(history.sketch) == "C1CCCC1"
    history.undo()
    assert to_smiles(history.sketch) == "CCCCCC"


def test_merge_drops_duplicate_bonds():
    sketch = Sketch()
    history = History(sketch)
    change = history.change("x")
    a = change.add_atom(Atom(0, 0))
    b = change.add_atom(Atom(1, 0))
    c = change.add_atom(Atom(1, 0.05))
    change.add_bond(a, b)
    change.add_bond(a, c, Order.DOUBLE)
    change.merge_atom(c, b)
    assert len(sketch.atoms) == 2 and len(sketch.bonds) == 1


def test_rotate_and_flip_keep_bond_lengths():
    history = _history("CC(C)CC(=O)O")
    change = history.change("x")
    ids = list(history.sketch.atoms)
    edits.rotate(change, ids, (0.0, 0.0), math.radians(37))
    edits.flip(change, ids, horizontal=False)
    for bond in history.sketch.bonds.values():
        a, b = history.sketch.atoms[bond.begin], history.sketch.atoms[bond.end]
        assert math.dist(a.pos, b.pos) == pytest.approx(1.0, abs=0.02)
