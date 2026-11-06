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

"""What to draw for a sketch, in model units, so the canvas only paints.

Conventions follow common structure-drawing practice (ACS style, as RDKit's
``SetACS1996Mode`` also draws): carbon is unlabelled unless it is alone or carries a charge,
isotope or radical; a heteroatom shows its hydrogens; a ring double bond has its second line
inside the ring; a chain double bond with substituents on one side has its second line on that
side, and otherwise the two lines are centred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from chemsuite.sketch import geometry
from chemsuite.sketch.model import Group, Order, Sketch, Stereo
from chemsuite.sketch.valence import Report

Point = tuple[float, float]

#: Second line of a double bond: offset and how far it is pulled in at each end (bond lengths).
DOUBLE_OFFSET = 0.18
INNER_TRIM = 0.14
TRIPLE_OFFSET = 0.15
WEDGE_WIDTH = 0.2
HASH_STEP = 0.11
#: Radius around a labelled atom that bonds stop short of.
LABEL_GAP = 0.24


@dataclass(frozen=True)
class Line:
    a: Point
    b: Point
    bond: int
    dashed: bool = False


@dataclass(frozen=True)
class Wedge:
    """A filled triangle: ``tip`` at the stereocentre, ``base`` the two far corners."""

    tip: Point
    base: tuple[Point, Point]
    bond: int


@dataclass(frozen=True)
class Wavy:
    points: tuple[Point, ...]
    bond: int


@dataclass(frozen=True)
class Label:
    """An atom label as runs of (text, "normal" | "sub" | "sup").

    ``centre`` is the index of the run that sits on the atom (the element symbol); other runs
    extend to its left or right as they come.
    """

    atom: int
    pos: Point
    runs: tuple[tuple[str, str], ...]
    centre: int
    symbol: str
    problem: str = ""


@dataclass
class Scene:
    lines: list[Line] = field(default_factory=list)
    wedges: list[Wedge] = field(default_factory=list)
    wavy: list[Wavy] = field(default_factory=list)
    labels: dict[int, Label] = field(default_factory=dict)
    problems: dict[int, str] = field(default_factory=dict)
    #: small text beside query bonds and attachment points: (position, text)
    tags: list[tuple[Point, str]] = field(default_factory=list)
    #: attachment-point stubs: wavy lines from the atom
    stubs: list[Wavy] = field(default_factory=list)


#: Tags on query bonds (Ketcher and ChemDraw label these bonds in the same spirit).
QUERY_TAGS = {
    Order.SINGLE_OR_DOUBLE: "S/D",
    Order.SINGLE_OR_AROMATIC: "S/A",
    Order.DOUBLE_OR_AROMATIC: "D/A",
}
TOPOLOGY_TAGS = {1: "Rng", 2: "Chn"}


def shows_label(sketch: Sketch, aid: int) -> bool:
    atom = sketch.atoms[aid]
    if atom.symbol != "C" or atom.is_query:
        return True
    return bool(atom.charge or atom.isotope or atom.radical) or sketch.degree(aid) == 0


def _charge_text(charge: int) -> str:
    if not charge:
        return ""
    sign = "+" if charge > 0 else "−"
    return f"{abs(charge) if abs(charge) > 1 else ''}{sign}"


def label_runs(sketch: Sketch, aid: int, hydrogens: int) -> tuple[list[tuple[str, str]], int]:
    """The runs of an atom's label and the index of its element symbol."""
    atom = sketch.atoms[aid]
    symbol = atom.label
    runs: list[tuple[str, str]] = []
    if atom.is_query:
        pass  # R1, A, [C,N,O]: no isotope; the label is the whole symbol
    elif atom.symbol == "H" and atom.isotope in (2, 3):
        symbol = "D" if atom.isotope == 2 else "T"
    elif atom.isotope:
        runs.append((str(atom.isotope), "sup"))
    centre = len(runs)
    runs.append((symbol, "normal"))
    h_runs: list[tuple[str, str]] = []
    if hydrogens:
        h_runs.append(("H", "normal"))
        if hydrogens > 1:
            h_runs.append((str(hydrogens), "sub"))
    tail = []
    if atom.radical:
        tail.append(("•" * atom.radical, "sup"))
    if atom.charge:
        tail.append((_charge_text(atom.charge), "sup"))
    if h_runs and _hydrogens_left(sketch, aid):
        runs = h_runs + runs
        centre += len(h_runs)
        runs += tail
    else:
        runs += h_runs + tail
    return runs, centre


def _hydrogens_left(sketch: Sketch, aid: int) -> bool:
    """Hydrogens go on the side away from the bonds: left when every bond points right."""
    here = sketch.atoms[aid]
    xs = [sketch.atoms[n].x - here.x for n in sketch.neighbours(aid)]
    return bool(xs) and all(x > 0.1 for x in xs)


def group_runs(sketch: Sketch, group: Group) -> tuple[list[tuple[str, str]], int]:
    """A contracted abbreviation's label, written from its attachment atom outwards.

    The first letter sits on the attachment atom (CO₂Me); when the bond comes in from the
    right the mirrored form is used and its last letter sits there (MeO₂C).
    """
    from chemsuite.sketch.abbrev import display, find

    entry = find(group.label)
    here = sketch.atoms[group.attach]
    outside = [n for n in sketch.neighbours(group.attach) if n not in group.atoms]
    from_right = bool(outside) and all(sketch.atoms[n].x > here.x + 0.1 for n in outside)
    text = entry.left_label if (entry and from_right) else group.label
    runs = display(text)
    if from_right:
        last, kind = runs[-1]
        runs[-1:] = [(last[:-1], kind), (last[-1], kind)] if len(last) > 1 else [(last, kind)]
        runs = [r for r in runs if r[0]]
        return runs, len(runs) - 1
    first, kind = runs[0]
    runs[:1] = [(first[0], kind), (first[1:], kind)] if len(first) > 1 else [(first, kind)]
    return [r for r in runs if r[0]], 0


def build(sketch: Sketch, report: Report) -> Scene:
    scene = Scene(problems=dict(report.problems))
    hidden = sketch.hidden_atoms()
    labelled = {aid for aid in sketch.atoms if shows_label(sketch, aid)}
    labelled |= set(report.problems)
    labelled -= hidden
    contracted = {g.attach: g for g in sketch.groups.values() if not g.expanded}
    for aid, group in contracted.items():
        problem = next((report.problems[a] for a in group.atoms if a in report.problems), "")
        runs, centre = group_runs(sketch, group)
        scene.labels[aid] = Label(aid, sketch.atoms[aid].pos, tuple(runs), centre, "C", problem)
        labelled.add(aid)
    for aid in labelled - set(contracted):
        runs, centre = label_runs(sketch, aid, report.hydrogens.get(aid, 0))
        scene.labels[aid] = Label(
            aid,
            sketch.atoms[aid].pos,
            tuple(runs),
            centre,
            sketch.atoms[aid].symbol,
            report.problems.get(aid, ""),
        )
    rings = geometry.ring_centres(sketch)
    for bid, bond in sketch.bonds.items():
        if bond.begin in hidden or bond.end in hidden:
            continue
        a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
        if math.dist(a, b) < 1e-6:
            continue
        p, q = _trim(a, b, bond.begin in labelled, bond.end in labelled)
        tag = " ".join(t for t in (QUERY_TAGS.get(bond.order, ""),
                                   TOPOLOGY_TAGS.get(bond.topology, "")) if t)  # fmt: skip
        if tag:
            off = _offset(a, b, -0.3)
            scene.tags.append((((off[0][0] + off[1][0]) / 2, (off[0][1] + off[1][1]) / 2), tag))
        if bond.order.is_query:  # any or either-of: a dashed line, the tag says which
            scene.lines.append(Line(p, q, bid, dashed=True))
            continue
        if bond.stereo is Stereo.WEDGE and bond.order is Order.SINGLE:
            scene.wedges.append(_wedge(p, q, bid))
        elif bond.stereo is Stereo.HASH and bond.order is Order.SINGLE:
            scene.lines.extend(_hashes(p, q, bid))
        elif bond.stereo is Stereo.EITHER and bond.order is Order.SINGLE:
            scene.wavy.append(Wavy(_wave(p, q), bid))
        elif bond.order is Order.SINGLE:
            scene.lines.append(Line(p, q, bid))
        elif bond.order is Order.TRIPLE:
            scene.lines.append(Line(p, q, bid))
            for side in (1, -1):
                scene.lines.append(Line(*_offset(p, q, side * TRIPLE_OFFSET), bid))
        elif bond.order is Order.DOUBLE and bond.stereo is Stereo.EITHER:
            # unknown E/Z: the two lines cross (CTfile bond stereo 3, as ChemDraw draws it)
            a1, b1 = _offset(p, q, DOUBLE_OFFSET / 2)
            a2, b2 = _offset(p, q, -DOUBLE_OFFSET / 2)
            scene.lines.extend((Line(a1, b2, bid), Line(a2, b1, bid)))
        else:
            side = _double_side(sketch, bid, rings)
            dashed = bond.order is Order.AROMATIC
            if side == 0:
                for s in (1, -1):
                    scene.lines.append(Line(*_offset(p, q, s * DOUBLE_OFFSET / 2), bid))
            else:
                scene.lines.append(Line(p, q, bid))
                inner = _offset(p, q, side * DOUBLE_OFFSET)
                trim_a = 0.0 if bond.begin in labelled else INNER_TRIM
                trim_b = 0.0 if bond.end in labelled else INNER_TRIM
                scene.lines.append(Line(*_shorten(*inner, trim_a, trim_b), bid, dashed=dashed))
    for aid, atom in sketch.atoms.items():
        if atom.attach and aid not in hidden:
            _attachment(sketch, aid, scene)
    return scene


def _attachment(sketch: Sketch, aid: int, scene: Scene) -> None:
    """A wavy stub (and its number) for each attachment point of an R-group member."""
    here = sketch.atoms[aid].pos
    taken = [geometry.angle_of(here, sketch.atoms[n].pos) for n in sketch.neighbours(aid)]
    for number in (1, 2):
        if not sketch.atoms[aid].attach & number:
            continue
        angle = _free_direction(taken)
        taken.append(angle)
        start = geometry.along(here, angle, LABEL_GAP)
        end = geometry.along(here, angle, 0.7)
        scene.stubs.append(Wavy(_wave(start, end, waves=2), -aid))
        scene.tags.append((geometry.along(here, angle, 0.9), str(number)))


def _free_direction(taken: list[float]) -> float:
    if not taken:
        return 0.0
    if len(taken) == 1:
        return taken[0] + math.pi
    ordered = sorted(t % (2 * math.pi) for t in taken)
    gaps = [((ordered[(i + 1) % len(ordered)] - a) % (2 * math.pi) or 2 * math.pi, a)
            for i, a in enumerate(ordered)]  # fmt: skip
    width, start = max(gaps)
    return start + width / 2


def _trim(a: Point, b: Point, cut_a: bool, cut_b: bool) -> tuple[Point, Point]:
    length = math.dist(a, b)
    ca = LABEL_GAP if cut_a else 0.0
    cb = LABEL_GAP if cut_b else 0.0
    if ca + cb >= length * 0.9:
        return a, b
    return _shorten(a, b, ca, cb)


def _shorten(a: Point, b: Point, from_a: float, from_b: float) -> tuple[Point, Point]:
    length = math.dist(a, b)
    ux, uy = (b[0] - a[0]) / length, (b[1] - a[1]) / length
    return (a[0] + ux * from_a, a[1] + uy * from_a), (b[0] - ux * from_b, b[1] - uy * from_b)


def _offset(a: Point, b: Point, d: float) -> tuple[Point, Point]:
    """The segment a-b moved ``d`` to its left (positive) or right."""
    length = math.dist(a, b)
    nx, ny = -(b[1] - a[1]) / length, (b[0] - a[0]) / length
    return (a[0] + nx * d, a[1] + ny * d), (b[0] + nx * d, b[1] + ny * d)


def _double_side(sketch: Sketch, bid: int, rings: dict[int, Point]) -> int:
    """+1 to draw the second line on the left of begin→end, -1 on the right, 0 centred."""
    bond = sketch.bonds[bid]
    a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
    if bid in rings:
        return 1 if geometry.side_of(a, b, rings[bid]) > 0 else -1
    left = right = 0
    for end in (bond.begin, bond.end):
        for n in sketch.neighbours(end):
            if n in (bond.begin, bond.end):
                continue
            s = geometry.side_of(a, b, sketch.atoms[n].pos)
            if s > 1e-6:
                left += 1
            elif s < -1e-6:
                right += 1
    ends = (sketch.degree(bond.begin), sketch.degree(bond.end))
    if left == right or (1 in ends and left and right):
        return 0
    return 1 if left > right else -1


def _wedge(tip: Point, end: Point, bid: int) -> Wedge:
    base = _offset(tip, end, WEDGE_WIDTH / 2)[1], _offset(tip, end, -WEDGE_WIDTH / 2)[1]
    return Wedge(tip, base, bid)


def _hashes(tip: Point, end: Point, bid: int) -> list[Line]:
    length = math.dist(tip, end)
    count = max(3, int(length / HASH_STEP))
    lines = []
    for i in range(1, count + 1):
        t = i / count
        centre = (tip[0] + (end[0] - tip[0]) * t, tip[1] + (end[1] - tip[1]) * t)
        half = WEDGE_WIDTH / 2 * t
        a, b = _offset(tip, end, half), _offset(tip, end, -half)
        da = (a[0][0] - tip[0], a[0][1] - tip[1])
        db = (b[0][0] - tip[0], b[0][1] - tip[1])
        lines.append(
            Line(
                (centre[0] + da[0], centre[1] + da[1]), (centre[0] + db[0], centre[1] + db[1]), bid
            )
        )
    return lines


def _wave(a: Point, b: Point, waves: int = 4) -> tuple[Point, ...]:
    length = math.dist(a, b)
    ux, uy = (b[0] - a[0]) / length, (b[1] - a[1]) / length
    points = []
    steps = waves * 8
    for i in range(steps + 1):
        t = i / steps
        w = 0.06 * math.sin(t * waves * 2 * math.pi)
        points.append((a[0] + ux * length * t - uy * w, a[1] + uy * length * t + ux * w))
    return tuple(points)
