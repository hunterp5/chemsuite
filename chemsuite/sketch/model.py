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
import re
from collections.abc import Iterator
from dataclasses import dataclass, field, replace
from enum import Enum


class Order(Enum):
    """Bond order as drawn; the value is the CTfile bond type.

    ``AROMATIC`` comes only from files or the aromatic bond key. 5–8 are CTfile query bonds,
    which match more than one kind of bond in a structure search.
    """

    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4
    SINGLE_OR_DOUBLE = 5
    SINGLE_OR_AROMATIC = 6
    DOUBLE_OR_AROMATIC = 7
    ANY = 8

    @property
    def is_query(self) -> bool:
        return self.value >= 5


#: Bond query names as menus show them.
BOND_QUERIES: dict[Order, str] = {
    Order.ANY: "Any bond",
    Order.SINGLE_OR_DOUBLE: "Single or double",
    Order.SINGLE_OR_AROMATIC: "Single or aromatic",
    Order.DOUBLE_OR_AROMATIC: "Double or aromatic",
    Order.AROMATIC: "Aromatic",
}
#: CTfile bond topology: 0 either, 1 ring, 2 chain.
TOPOLOGY = {0: "Ring or chain", 1: "Ring only", 2: "Chain only"}
#: CTfile attachment points (``M  APO``): 1 first, 2 second, 3 both.
ATTACH_POINTS = {0: "None", 1: "First", 2: "Second", 3: "First and second"}


class Stereo(Enum):
    """Wedge style of a single bond, pointing from its begin atom (the narrow end)."""

    NONE = 0
    WEDGE = 1  # solid: towards the viewer
    HASH = 6  # hashed: away from the viewer
    EITHER = 4  # wavy single bond, or crossed double bond: unknown configuration


#: CTfile query atoms (BIOVIA CTfile formats, atom block "atom symbol"), as RDKit reads them.
QUERY_ATOMS: dict[str, str] = {
    "A": "Any atom except hydrogen",
    "Q": "Any heteroatom (not carbon or hydrogen)",
    "X": "Any halogen",
    "M": "Any metal",
    "AH": "Any atom, including hydrogen",
    "QH": "Any atom except carbon",
    "XH": "Any halogen or hydrogen",
    "MH": "Any metal or hydrogen",
    "*": "Any atom",
}
RGROUP = "R#"  # an R-group; ``Atom.rgroup`` is its number
ATOM_LIST = "L"  # one of (or none of) ``Atom.atom_list``


MAX_RGROUP = 99
_RGROUP_TEXT = re.compile(r"R(\d{0,2})")
_LIST_TEXT = re.compile(r"(!?)\[([A-Za-z]{1,2}(?:\s*,\s*[A-Za-z]{1,2})*)\]")


def query_fields(text: str) -> dict | None:
    """The Atom fields for a typed query label, or None when *text* is not one.

    ``R`` / ``R1``-``R99`` (an R-group; plain ``R`` means R1), the CTfile wildcards in
    QUERY_ATOMS (``A``, ``Q``, ``X``, ``M``, their ``H`` forms and ``*``), and atom lists
    ``[C,N,O]`` or ``![F,Cl]`` ("none of these"). List members must be real elements.
    """
    text = text.strip()
    if text in QUERY_ATOMS:
        return {"symbol": text, "rgroup": 0, "atom_list": (), "list_not": False, "isotope": 0}
    match = _RGROUP_TEXT.fullmatch(text)
    if match:
        number = int(match.group(1) or 1)
        if not 1 <= number <= MAX_RGROUP:
            return None
        return {"symbol": RGROUP, "rgroup": number, "atom_list": (), "list_not": False,
                "isotope": 0}  # fmt: skip
    match = _LIST_TEXT.fullmatch(text)
    if match:
        from rdkit import Chem, rdBase

        table = Chem.GetPeriodicTable()
        elements: list[str] = []
        for item in match.group(2).split(","):
            name = item.strip()
            name = name[:1].upper() + name[1:].lower()
            try:
                with rdBase.BlockLogs():
                    if table.GetAtomicNumber(name) <= 0:
                        return None
            except (RuntimeError, ValueError):
                return None
            if name not in elements:
                elements.append(name)
        if not 1 <= len(elements) <= 16:  # CTfile allows up to 16 list members
            return None
        return {"symbol": ATOM_LIST, "rgroup": 0, "atom_list": tuple(elements),
                "list_not": bool(match.group(1)), "isotope": 0}  # fmt: skip
    return None


_ENHANCED = re.compile(r"(abs)|(and|or)([1-9]\d?)")


def parse_enhanced(text: str) -> tuple[str, int] | None:
    """("abs" | "and" | "or", group number) for an ``Atom.enhanced`` value; None for ""."""
    match = _ENHANCED.fullmatch(text)
    if not match:
        return None
    if match.group(1):
        return ("abs", 0)
    return (match.group(2), int(match.group(3)))


def enhanced_text(value: str) -> str:
    """How a group is shown beside its centre: abs, &1, or1."""
    parsed = parse_enhanced(value)
    if parsed is None:
        return ""
    kind, n = parsed
    return {"abs": "abs", "and": f"&{n}", "or": f"or{n}"}[kind]


def next_enhanced(sketch: Sketch, kind: str) -> int:
    """The lowest unused group number of *kind* ("and" or "or")."""
    used = {
        p[1] for a in sketch.atoms.values() if (p := parse_enhanced(a.enhanced)) and p[0] == kind
    }
    return next(n for n in range(1, 100) if n not in used)


def next_rgroup(sketch: Sketch, skip: int | None = None) -> int:
    """The lowest R-group number not used in *sketch* (ignoring atom *skip*)."""
    used = {a.rgroup for aid, a in sketch.atoms.items() if a.symbol == RGROUP and aid != skip}
    return next(n for n in range(1, MAX_RGROUP + 2) if n not in used)


#: The fields that make a plain element atom (no R-group number or list).
ELEMENT_FIELDS = {"rgroup": 0, "atom_list": (), "list_not": False}


@dataclass(frozen=True)
class Atom:
    x: float
    y: float
    symbol: str = "C"  # an element, a QUERY_ATOMS key, RGROUP or ATOM_LIST
    charge: int = 0
    isotope: int = 0  # mass number; 0 means natural abundance
    radical: int = 0  # unpaired electrons
    rgroup: int = 0  # R-group number when symbol is RGROUP
    atom_list: tuple[str, ...] = ()  # elements when symbol is ATOM_LIST
    list_not: bool = False  # the list means "none of these"
    enhanced: str = ""  # enhanced stereo group: "abs", "and1", "or2" … ("" = absolute)
    attach: int = 0  # R-group member attachment point: 1 first, 2 second, 3 both (M  APO)
    mapping: int = 0  # reaction atom-atom map number (CTfile aamap); 0 = unmapped

    @property
    def pos(self) -> tuple[float, float]:
        return (self.x, self.y)

    @property
    def is_query(self) -> bool:
        """An R-group, wildcard or atom list rather than one element."""
        return self.symbol in QUERY_ATOMS or self.symbol in (RGROUP, ATOM_LIST)

    @property
    def label(self) -> str:
        """The text drawn for the atom's symbol: R1, A, [C,N,O], ![F,Cl], or the element."""
        if self.symbol == RGROUP:
            return f"R{self.rgroup}" if self.rgroup else "R"
        if self.symbol == ATOM_LIST:
            return ("!" if self.list_not else "") + "[" + ",".join(self.atom_list) + "]"
        return self.symbol

    def moved(self, x: float, y: float) -> Atom:
        return replace(self, x=x, y=y)


@dataclass(frozen=True)
class Bond:
    begin: int
    end: int
    order: Order = Order.SINGLE
    stereo: Stereo = Stereo.NONE
    topology: int = 0  # query: 0 ring or chain, 1 ring only, 2 chain only (CTfile)

    def other(self, atom: int) -> int:
        return self.end if atom == self.begin else self.begin

    def flipped(self) -> Bond:
        return replace(self, begin=self.end, end=self.begin)

    @property
    def is_plain(self) -> bool:
        """A single, double or triple bond without a wedge."""
        return self.order in (Order.SINGLE, Order.DOUBLE, Order.TRIPLE) and (
            self.stereo is Stereo.NONE and not self.topology
        )

    @property
    def is_query(self) -> bool:
        """A query bond: a CTfile query type or a ring/chain restriction."""
        return self.order.is_query or bool(self.topology)


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


class ArrowKind(Enum):
    """Reaction arrow styles (a subset of Ketcher's ``RxnArrowMode``, ``rxnArrow.ts``)."""

    PLAIN = "plain"  # reactants -> products
    EQUILIBRIUM = "equilibrium"  # two half-headed arrows
    RETROSYNTHETIC = "retrosynthetic"  # open double-shafted arrow: target => precursors
    FAILED = "failed"  # a plain arrow crossed out: the reaction did not work


ARROW_TITLES = {
    ArrowKind.PLAIN: "Reaction Arrow",
    ArrowKind.EQUILIBRIUM: "Equilibrium Arrow",
    ArrowKind.RETROSYNTHETIC: "Retrosynthetic Arrow",
    ArrowKind.FAILED: "Failed Reaction Arrow",
}


@dataclass(frozen=True)
class Arrow:
    """A reaction arrow from tail (x0, y0) to head (x1, y1)."""

    x0: float
    y0: float
    x1: float
    y1: float
    kind: ArrowKind = ArrowKind.PLAIN

    @property
    def tail(self) -> tuple[float, float]:
        return (self.x0, self.y0)

    @property
    def head(self) -> tuple[float, float]:
        return (self.x1, self.y1)

    @property
    def centre(self) -> tuple[float, float]:
        return ((self.x0 + self.x1) / 2, (self.y0 + self.y1) / 2)

    def moved(self, dx: float, dy: float) -> Arrow:
        return replace(self, x0=self.x0 + dx, y0=self.y0 + dy, x1=self.x1 + dx, y1=self.y1 + dy)


@dataclass(frozen=True)
class Plus:
    """A reaction plus sign at (x, y)."""

    x: float
    y: float

    @property
    def pos(self) -> tuple[float, float]:
        return (self.x, self.y)

    def moved(self, dx: float, dy: float) -> Plus:
        return Plus(self.x + dx, self.y + dy)


@dataclass
class Sketch:
    atoms: dict[int, Atom] = field(default_factory=dict)
    bonds: dict[int, Bond] = field(default_factory=dict)
    next_id: int = 1
    groups: dict[int, Group] = field(default_factory=dict)
    arrows: dict[int, Arrow] = field(default_factory=dict)
    pluses: dict[int, Plus] = field(default_factory=dict)

    # -- queries ------------------------------------------------------------------------------

    def is_empty(self) -> bool:
        """Nothing drawn: no atoms, arrows or plus signs."""
        return not self.atoms and not self.arrows and not self.pluses

    @property
    def is_reaction(self) -> bool:
        """The drawing has a reaction arrow (see ``chemsuite.sketch.reaction``)."""
        return bool(self.arrows)

    def arrow_at(self, x: float, y: float, within: float) -> int | None:
        best, best_d = None, within
        for rid, arrow in self.arrows.items():
            d = _segment_distance((x, y), arrow.tail, arrow.head)
            if d <= best_d:
                best, best_d = rid, d
        return best

    def plus_at(self, x: float, y: float, within: float) -> int | None:
        best, best_d = None, within
        for pid, plus in self.pluses.items():
            d = math.hypot(plus.x - x, plus.y - y)
            if d <= best_d:
                best, best_d = pid, d
        return best

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
        """(min x, min y, max x, max y) of atoms, arrows and plus signs; zeros when empty."""
        xs = [a.x for a in self.atoms.values()] + [p.x for p in self.pluses.values()]
        ys = [a.y for a in self.atoms.values()] + [p.y for p in self.pluses.values()]
        for arrow in self.arrows.values():
            xs += [arrow.x0, arrow.x1]
            ys += [arrow.y0, arrow.y1]
        if not xs:
            return (0.0, 0.0, 0.0, 0.0)
        return (min(xs), min(ys), max(xs), max(ys))

    def copy(self) -> Sketch:
        return Sketch(
            dict(self.atoms),
            dict(self.bonds),
            self.next_id,
            dict(self.groups),
            dict(self.arrows),
            dict(self.pluses),
        )

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
