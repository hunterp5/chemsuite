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

"""Every change to a sketch, as an undoable edit.

An edit is a list of steps; each step replaces one atom or bond value (``None`` meaning absent)
and remembers what was there, so undo replays the steps backwards with old and new swapped.

Build an edit with ``Change``, which applies each step as it is recorded so later steps see
the result, then hand it to ``History.commit``. The functions below the classes are the edits
the tools make, each following the Ketcher rule it cites.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field, replace

from chemsuite.sketch import geometry
from chemsuite.sketch.model import Atom, Bond, Order, Sketch, Stereo

Point = tuple[float, float]


@dataclass(frozen=True)
class Step:
    kind: str  # "atom" or "bond"
    id: int
    old: Atom | Bond | None
    new: Atom | Bond | None

    def apply(self, sketch: Sketch, forward: bool = True) -> None:
        value = self.new if forward else self.old
        store = sketch.atoms if self.kind == "atom" else sketch.bonds
        if value is None:
            store.pop(self.id, None)
        else:
            store[self.id] = value


@dataclass
class Edit:
    name: str
    steps: list[Step] = field(default_factory=list)

    def undo(self, sketch: Sketch) -> None:
        for step in reversed(self.steps):
            step.apply(sketch, forward=False)

    def redo(self, sketch: Sketch) -> None:
        for step in self.steps:
            step.apply(sketch)


class Change:
    """Records an edit while applying it to the sketch."""

    def __init__(self, sketch: Sketch, name: str):
        self.sketch = sketch
        self.edit = Edit(name)

    def _step(self, kind: str, id_: int, new: Atom | Bond | None) -> None:
        store = self.sketch.atoms if kind == "atom" else self.sketch.bonds
        old = store.get(id_)
        if old == new:
            return
        step = Step(kind, id_, old, new)
        step.apply(self.sketch)
        self.edit.steps.append(step)

    def add_atom(self, atom: Atom) -> int:
        aid = self.sketch.new_id()
        self._step("atom", aid, atom)
        return aid

    def set_atom(self, aid: int, atom: Atom) -> None:
        self._step("atom", aid, atom)

    def remove_atom(self, aid: int) -> None:
        for bid, _ in list(self.sketch.bonds_of(aid)):
            self._step("bond", bid, None)
        self._step("atom", aid, None)

    def add_bond(
        self, a: int, b: int, order: Order = Order.SINGLE, stereo: Stereo = Stereo.NONE
    ) -> int:
        """A bond from ``a`` to ``b``; an existing bond between them is changed instead."""
        existing = self.sketch.bond_between(a, b)
        if existing is not None:
            self._step("bond", existing, Bond(a, b, order, stereo))
            return existing
        bid = self.sketch.new_id()
        self._step("bond", bid, Bond(a, b, order, stereo))
        return bid

    def set_bond(self, bid: int, bond: Bond) -> None:
        self._step("bond", bid, bond)

    def remove_bond(self, bid: int) -> None:
        self._step("bond", bid, None)

    def move_atom(self, aid: int, x: float, y: float) -> None:
        self.set_atom(aid, self.sketch.atoms[aid].moved(x, y))

    def merge_atom(self, source: int, target: int) -> None:
        """Move ``source``'s bonds onto ``target`` and delete ``source`` (Ketcher ``atomMerge``).

        A bond that would join ``target`` to itself, or duplicate one it already has, is
        dropped; the target's own element and charge are kept.
        """
        for bid, bond in list(self.sketch.bonds_of(source)):
            other = bond.other(source)
            self.remove_bond(bid)
            if other == target or self.sketch.bond_between(other, target) is not None:
                continue
            ends = (target, other) if bond.begin == source else (other, target)
            self._step("bond", self.sketch.new_id(), replace(bond, begin=ends[0], end=ends[1]))
        self._step("atom", source, None)


class History:
    """The sketch with its undo and redo stacks."""

    def __init__(self, sketch: Sketch | None = None):
        self.sketch = sketch if sketch is not None else Sketch()
        self._undo: list[Edit] = []
        self._redo: list[Edit] = []
        self.listeners: list[Callable[[], None]] = []

    def change(self, name: str) -> Change:
        return Change(self.sketch, name)

    def commit(self, change: Change) -> Edit | None:
        """Keep a finished change for undo. Changes that changed nothing are dropped."""
        if not change.edit.steps:
            return None
        self._undo.append(change.edit)
        self._redo.clear()
        self._notify()
        return change.edit

    def discard(self, change: Change) -> None:
        """Roll back a change that was being previewed (for example, a cancelled drag)."""
        change.edit.undo(self.sketch)

    def undo(self) -> bool:
        if not self._undo:
            return False
        edit = self._undo.pop()
        edit.undo(self.sketch)
        self._redo.append(edit)
        self._notify()
        return True

    def redo(self) -> bool:
        if not self._redo:
            return False
        edit = self._redo.pop()
        edit.redo(self.sketch)
        self._undo.append(edit)
        self._notify()
        return True

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    @property
    def undo_name(self) -> str:
        return self._undo[-1].name if self._undo else ""

    @property
    def redo_name(self) -> str:
        return self._redo[-1].name if self._redo else ""

    def replace(self, sketch: Sketch, name: str = "Replace drawing") -> None:
        """Swap in a whole new drawing (paste over, open), as one undoable edit."""
        change = self.change(name)
        for aid in list(self.sketch.atoms):
            change.remove_atom(aid)
        mapping = {aid: change.add_atom(atom) for aid, atom in sketch.atoms.items()}
        for bond in sketch.bonds.values():
            change.add_bond(mapping[bond.begin], mapping[bond.end], bond.order, bond.stereo)
        self.commit(change)

    def _notify(self) -> None:
        for listener in list(self.listeners):
            listener()


# -- the edits the tools make ----------------------------------------------------------------


def bond_from_atom(
    change: Change, atom: int, order: Order = Order.SINGLE, stereo: Stereo = Stereo.NONE
) -> int:
    """A new bond from ``atom`` in the direction Ketcher's ``atomForNewBond`` picks.

    If the far end lands on an existing atom (within 0.1 bond lengths), the bond joins it.
    """
    end, existing = geometry.new_bond_target(change.sketch, atom, order)
    if existing is None:
        existing = change.add_atom(Atom(*end))
    return change.add_bond(atom, existing, order, stereo)


def bond_to_point(
    change: Change,
    atom: int,
    point: Point,
    order: Order = Order.SINGLE,
    stereo: Stereo = Stereo.NONE,
    snap_to: float = 0.0,
) -> int:
    """A bond from ``atom`` to the atom at ``point`` (within ``snap_to``), or to a new atom."""
    target = change.sketch.atom_at(*point, within=snap_to, exclude={atom}) if snap_to else None
    if target is None:
        target = change.add_atom(Atom(*point))
    return change.add_bond(atom, target, order, stereo)


def bond_on_canvas(
    change: Change, point: Point, order: Order = Order.SINGLE, stereo: Stereo = Stereo.NONE
) -> int:
    """A lone bond clicked onto empty canvas (Ketcher ``tool/bond.ts``)."""
    begin = change.add_atom(Atom(*point))
    end = change.add_atom(Atom(*geometry.along(point, geometry.first_bond_angle(order))))
    return change.add_bond(begin, end, order, stereo)


PLAIN_CYCLE = (Order.SINGLE, Order.DOUBLE, Order.TRIPLE)


def click_bond(change: Change, bid: int, order: Order, stereo: Stereo = Stereo.NONE) -> None:
    """Clicking a bond with a bond tool (Ketcher ``bondChangingAction``).

    - The plain single-bond tool on a plain bond steps single → double → triple → single.
    - A wedge tool on a bond that already has that wedge flips its direction.
    - Otherwise the bond becomes the tool's type.
    """
    bond = change.sketch.bonds[bid]
    if stereo is not Stereo.NONE and bond.order is order and bond.stereo is stereo:
        change.set_bond(bid, bond.flipped())
        return
    if order is Order.SINGLE and stereo is Stereo.NONE and bond.is_plain:
        order = PLAIN_CYCLE[(PLAIN_CYCLE.index(bond.order) + 1) % len(PLAIN_CYCLE)]
    change.set_bond(bid, replace(bond, order=order, stereo=stereo))


def set_bond_order(change: Change, bid: int, order: Order) -> None:
    """A hover key over a bond: that order, without a wedge."""
    change.set_bond(bid, replace(change.sketch.bonds[bid], order=order, stereo=Stereo.NONE))


#: Labels that are an element with an isotope.
ISOTOPE_LABELS = {"D": ("H", 2), "T": ("H", 3)}


def set_element(change: Change, aid: int, label: str) -> None:
    """Set an atom's element. Charge and radical stay; a different element clears the isotope."""
    symbol, isotope = ISOTOPE_LABELS.get(label, (label, 0))
    atom = change.sketch.atoms[aid]
    if symbol != atom.symbol or label in ISOTOPE_LABELS:
        change.set_atom(aid, replace(atom, symbol=symbol, isotope=isotope))


def atom_on_canvas(change: Change, point: Point, label: str = "C") -> int:
    symbol, isotope = ISOTOPE_LABELS.get(label, (label, 0))
    return change.add_atom(Atom(*point, symbol=symbol, isotope=isotope))


def add_charge(change: Change, aid: int, delta: int) -> None:
    atom = change.sketch.atoms[aid]
    change.set_atom(aid, replace(atom, charge=atom.charge + delta))


def chain(change: Change, start: int | Point, angle: float, bonds: int) -> list[int]:
    """A zig-zag carbon chain from an atom or a point (Ketcher ``fromChain``).

    Returns the chain's atoms, starting atom first.
    """
    first = start if isinstance(start, int) else change.add_atom(Atom(*start))
    origin = change.sketch.atoms[first].pos
    atoms = [first]
    for point in geometry.chain_points(origin, angle, bonds):
        nxt = change.add_atom(Atom(*point))
        change.add_bond(atoms[-1], nxt)
        atoms.append(nxt)
    return atoms


def erase(change: Change, atoms: Iterable[int] = (), bonds: Iterable[int] = ()) -> None:
    """Delete atoms (with their bonds) and bonds.

    Deleting a bond also deletes an end atom left with no bonds (Ketcher ``fromBondDeletion``);
    deleting an atom leaves its neighbours in place.
    """
    for aid in list(atoms):
        if aid in change.sketch.atoms:
            change.remove_atom(aid)
    for bid in list(bonds):
        bond = change.sketch.bonds.get(bid)
        if bond is None:
            continue
        change.remove_bond(bid)
        for end in (bond.begin, bond.end):
            if end in change.sketch.atoms and change.sketch.degree(end) == 0:
                change.remove_atom(end)


def move(
    change: Change, atoms: Iterable[int], dx: float, dy: float, merge: bool = True
) -> dict[int, int]:
    """Move atoms; with ``merge``, a moved atom that lands on a still one joins it.

    Returns the merges made, moved atom -> atom it joined (Ketcher ``atomMerge``).
    """
    moving = set(atoms)
    for aid in moving:
        atom = change.sketch.atoms[aid]
        change.move_atom(aid, atom.x + dx, atom.y + dy)
    merged: dict[int, int] = {}
    if merge:
        for aid in sorted(moving):
            atom = change.sketch.atoms[aid]
            target = change.sketch.atom_at(
                atom.x, atom.y, within=geometry.MERGE_DISTANCE * 2, exclude=moving
            )
            if target is not None:
                change.merge_atom(aid, target)
                merged[aid] = target
    return merged


def rotate(change: Change, atoms: Iterable[int], centre: Point, angle: float) -> None:
    c, s = math.cos(angle), math.sin(angle)
    for aid in atoms:
        atom = change.sketch.atoms[aid]
        x, y = atom.x - centre[0], atom.y - centre[1]
        change.move_atom(aid, centre[0] + x * c - y * s, centre[1] + x * s + y * c)


def flip(change: Change, atoms: Iterable[int], horizontal: bool) -> None:
    """Mirror atoms about the centre of their box (Ketcher Alt+H, Alt+V)."""
    ids = list(atoms)
    if not ids:
        return
    xs = [change.sketch.atoms[a].x for a in ids]
    ys = [change.sketch.atoms[a].y for a in ids]
    cx, cy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    for aid in ids:
        atom = change.sketch.atoms[aid]
        if horizontal:
            change.move_atom(aid, 2 * cx - atom.x, atom.y)
        else:
            change.move_atom(aid, atom.x, 2 * cy - atom.y)
