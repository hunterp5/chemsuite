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

"""The drawing: atoms and bonds with 2D coordinates.

Coordinates are in bond lengths with y pointing **up**, as in a molfile; the canvas flips y
when it draws. Atoms and bonds are immutable values kept under stable integer ids, so an edit
records old and new values and undo is exact (see ``edits.py``).

The model holds states RDKit would reject or rewrite, such as a five-bonded carbon in the
middle of an edit. RDKit judges the chemistry through ``valence.py`` and ``convert.py``.
"""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass, field, replace
from enum import Enum


class Order(Enum):
    """Bond order as drawn. ``AROMATIC`` comes only from files or the aromatic bond key."""

    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4


class Stereo(Enum):
    """Wedge style of a single bond, pointing from its begin atom (the narrow end)."""

    NONE = 0
    WEDGE = 1  # solid: towards the viewer
    HASH = 6  # hashed: away from the viewer
    EITHER = 4  # wavy: unknown configuration


@dataclass(frozen=True)
class Atom:
    x: float
    y: float
    symbol: str = "C"
    charge: int = 0
    isotope: int = 0  # mass number; 0 means natural abundance
    radical: int = 0  # unpaired electrons

    @property
    def pos(self) -> tuple[float, float]:
        return (self.x, self.y)

    def moved(self, x: float, y: float) -> Atom:
        return replace(self, x=x, y=y)


@dataclass(frozen=True)
class Bond:
    begin: int
    end: int
    order: Order = Order.SINGLE
    stereo: Stereo = Stereo.NONE

    def other(self, atom: int) -> int:
        return self.end if atom == self.begin else self.begin

    def flipped(self) -> Bond:
        return replace(self, begin=self.end, end=self.begin)

    @property
    def is_plain(self) -> bool:
        """A single, double or triple bond without a wedge."""
        return self.order is not Order.AROMATIC and self.stereo is Stereo.NONE


@dataclass(frozen=True)
class Group:
    """An abbreviation (a superatom S-group): its atoms stay real atoms in the drawing.

    ``attach`` is the group atom bonded to the rest of the molecule (or any group atom for a
    free-standing group). Contracted, the canvas draws ``label`` at ``attach`` in place of
    the atoms; molfiles carry it as a ``SUP`` S-group.
    """

    label: str
    atoms: frozenset[int]
    attach: int
    expanded: bool = False


@dataclass
class Sketch:
    atoms: dict[int, Atom] = field(default_factory=dict)
    bonds: dict[int, Bond] = field(default_factory=dict)
    next_id: int = 1
    groups: dict[int, Group] = field(default_factory=dict)

    # -- queries ------------------------------------------------------------------------------

    def is_empty(self) -> bool:
        return not self.atoms

    def bonds_of(self, atom: int) -> Iterator[tuple[int, Bond]]:
        for bid, bond in self.bonds.items():
            if atom in (bond.begin, bond.end):
                yield bid, bond

    def neighbours(self, atom: int) -> list[int]:
        return [bond.other(atom) for _, bond in self.bonds_of(atom)]

    def degree(self, atom: int) -> int:
        return sum(1 for _ in self.bonds_of(atom))

    def bond_between(self, a: int, b: int) -> int | None:
        for bid, bond in self.bonds.items():
            if {bond.begin, bond.end} == {a, b}:
                return bid
        return None

    def atom_at(
        self, x: float, y: float, within: float, exclude: set[int] | frozenset[int] = frozenset()
    ) -> int | None:
        """The atom closest to (x, y) within the given distance, if any."""
        best, best_d = None, within
        for aid, atom in self.atoms.items():
            if aid in exclude:
                continue
            d = math.hypot(atom.x - x, atom.y - y)
            if d <= best_d:
                best, best_d = aid, d
        return best

    def bond_at(self, x: float, y: float, within: float) -> int | None:
        """The bond whose segment passes closest to (x, y) within the given distance."""
        best, best_d = None, within
        for bid, bond in self.bonds.items():
            a, b = self.atoms[bond.begin], self.atoms[bond.end]
            d = _segment_distance((x, y), a.pos, b.pos)
            if d <= best_d:
                best, best_d = bid, d
        return best

    def fragment(self, atom: int) -> set[int]:
        """Every atom connected to ``atom``, including itself."""
        seen, todo = {atom}, [atom]
        while todo:
            for n in self.neighbours(todo.pop()):
                if n not in seen:
                    seen.add(n)
                    todo.append(n)
        return seen

    def bounds(self) -> tuple[float, float, float, float]:
        """(min x, min y, max x, max y) of the atoms; zeros for an empty sketch."""
        if not self.atoms:
            return (0.0, 0.0, 0.0, 0.0)
        xs = [a.x for a in self.atoms.values()]
        ys = [a.y for a in self.atoms.values()]
        return (min(xs), min(ys), max(xs), max(ys))

    def copy(self) -> Sketch:
        return Sketch(dict(self.atoms), dict(self.bonds), self.next_id, dict(self.groups))

    def group_of(self, atom: int) -> tuple[int, Group] | None:
        """The abbreviation an atom belongs to, if any."""
        for gid, group in self.groups.items():
            if atom in group.atoms:
                return gid, group
        return None

    def hidden_atoms(self) -> set[int]:
        """Atoms drawn as part of a contracted abbreviation's label (all but its anchor)."""
        out: set[int] = set()
        for group in self.groups.values():
            if not group.expanded:
                out |= group.atoms - {group.attach}
        return out

    # -- raw mutation (edits.py records these so they can be undone) --------------------------

    def new_id(self) -> int:
        nid = self.next_id
        self.next_id += 1
        return nid


def _segment_distance(p: tuple[float, float], a: tuple[float, float], b: tuple[float, float]):
    ax, ay = a
    dx, dy = b[0] - ax, b[1] - ay
    length2 = dx * dx + dy * dy
    t = 0.0 if length2 == 0 else max(0.0, min(1.0, ((p[0] - ax) * dx + (p[1] - ay) * dy) / length2))
    return math.hypot(p[0] - (ax + t * dx), p[1] - (ay + t * dy))
