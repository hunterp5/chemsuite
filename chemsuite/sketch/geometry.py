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

"""Where new atoms go (SKETCHER_PLAN §4.2).

The rules follow Ketcher at commit 2b031cc2 (Apache-2.0). Ketcher stores y pointing down and
this model y pointing up, so each rule is stated here by what it does on the page, which does
not depend on the axis direction, and tested that way in tests/sketch/test_geometry.py.

- Angle snapping: ``fracAngle`` in ``ketcher-core/src/application/editor/shared/utils.ts``.
- New bond from an atom: ``atomForNewBond`` in
  ``ketcher-core/src/application/editor/actions/utils.ts``.
- Chains: ``fromChain`` in ``ketcher-core/src/application/editor/actions/chain.ts`` and the
  chain tool in ``ketcher-react/src/script/editor/tool/chain.ts``.
"""

from __future__ import annotations

import math

from chemsuite.sketch.model import Order, Sketch

#: Dragged bonds snap to multiples of 15 degrees (Ketcher's ``FRAC = Math.PI / 12``).
FRAC = math.pi / 12
#: Atoms closer than this (in bond lengths) are the same place: a new bond that would end there
#: joins the existing atom instead (``findClosestAtom(..., 0.1)``).
MERGE_DISTANCE = 0.1
#: A first bond drawn on empty canvas: single bonds rise at 30 degrees, others are horizontal.
FIRST_SINGLE_ANGLE = math.pi / 6

Point = tuple[float, float]


def snap_angle(angle: float) -> float:
    """The nearest multiple of 15 degrees."""
    return round(angle / FRAC) * FRAC


def angle_of(a: Point, b: Point) -> float:
    return math.atan2(b[1] - a[1], b[0] - a[0])


def along(origin: Point, angle: float, length: float = 1.0) -> Point:
    return (origin[0] + length * math.cos(angle), origin[1] + length * math.sin(angle))


def first_bond_angle(order: Order) -> float:
    """Direction of a bond clicked onto empty canvas (Ketcher ``tool/bond.ts``)."""
    return FIRST_SINGLE_ANGLE if order is Order.SINGLE else 0.0


def _wants_straight(previous: Order, new: Order) -> bool:
    """Ketcher draws these pairs collinear: an alkyne, or cumulated double bonds (allenes)."""
    if previous is new and new in (Order.DOUBLE, Order.TRIPLE):
        return True
    return {previous, new} == {Order.SINGLE, Order.TRIPLE}


def new_bond_angle(sketch: Sketch, atom: int, order: Order = Order.SINGLE) -> float:
    """The direction for a new bond from ``atom`` (Ketcher ``atomForNewBond``).

    - No neighbours: horizontal.
    - One neighbour: 120 degrees from it, on the side away from that neighbour's own other
      bonds, so chains zig-zag; or straight on (180 degrees) for an alkyne or an allene.
    - Two or more: bisect the largest gap between the existing bonds.
    """
    here = sketch.atoms[atom].pos
    neighbours = [n for n in sketch.neighbours(atom) if _dist(sketch.atoms[n].pos, here) >= 0.1]
    if not neighbours:
        return 0.0
    if len(neighbours) == 1:
        (n,) = neighbours
        back = angle_of(here, sketch.atoms[n].pos)  # towards the neighbour
        previous = sketch.bonds[sketch.bond_between(atom, n)].order
        if _wants_straight(previous, order):
            return _norm(back + math.pi)
        forward = _norm(back + math.pi)  # neighbour -> atom
        # Ketcher's zig-zag check: if every other bond of the neighbour lies on one side of the
        # line neighbour -> atom, the new bond goes to the other side (a trans, zig-zag chain).
        sides = [
            _norm(angle_of(sketch.atoms[n].pos, sketch.atoms[m].pos) - forward)
            for m in sketch.neighbours(n)
            if m != atom and _dist(sketch.atoms[m].pos, sketch.atoms[n].pos) >= 0.1
        ]
        if sides and all(s >= 0.99 * math.pi for s in sides):
            return _norm(forward + math.pi / 3)  # others are clockwise, so turn anticlockwise
        return _norm(forward - math.pi / 3)
    angles = sorted(angle_of(here, sketch.atoms[n].pos) for n in neighbours)
    best_gap, best_start = -1.0, angles[0]
    for i, start in enumerate(angles):
        gap = _norm(angles[(i + 1) % len(angles)] - start) or 2 * math.pi
        if gap > best_gap:
            best_gap, best_start = gap, start
    return _norm(best_start + best_gap / 2)


def new_bond_target(
    sketch: Sketch, atom: int, order: Order = Order.SINGLE
) -> tuple[Point, int | None]:
    """Where a new bond from ``atom`` ends, and the existing atom there if one is that close."""
    end = along(sketch.atoms[atom].pos, new_bond_angle(sketch, atom, order))
    return end, sketch.atom_at(*end, within=MERGE_DISTANCE, exclude={atom})


def chain_length(drag: float) -> int:
    """Bonds in a chain dragged ``drag`` bond lengths: one per bond length begun."""
    return max(1, math.ceil(drag - 1e-9))


def chain_points(start: Point, angle: float, bonds: int) -> list[Point]:
    """Atoms of a zig-zag chain of ``bonds`` bonds from ``start`` along ``angle``.

    Ketcher ``fromChain``: atom *i* sits at (cos 30° · (i + 1), 0 or sin 30°) in the frame of
    the drag, alternating, so bonds are one bond length and 120 degrees apart. The first new
    atom is on the clockwise side of the drag line, as it is on Ketcher's page.
    """
    dx, dy = math.cos(math.pi / 6), math.sin(math.pi / 6)
    c, s = math.cos(angle), math.sin(angle)
    points = []
    for i in range(bonds):
        u, v = dx * (i + 1), (0.0 if i % 2 else -dy)
        points.append((start[0] + u * c - v * s, start[1] + u * s + v * c))
    return points


def side_of(a: Point, b: Point, p: Point) -> float:
    """Positive when ``p`` is left of the line a→b, negative when right."""
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])


def ring_centres(sketch: Sketch) -> dict[int, Point]:
    """For each ring bond, the centre of the smallest ring through it (up to 8 atoms)."""
    out: dict[int, Point] = {}
    for bid, bond in sketch.bonds.items():
        path = shortest_path(sketch, bond.begin, bond.end, skip=bid, limit=7)
        if path:
            pts = [sketch.atoms[a].pos for a in path]
            out[bid] = (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))
    return out


def shortest_path(sketch: Sketch, start: int, goal: int, skip: int, limit: int) -> list[int]:
    """Atoms on the shortest path start→goal not using bond ``skip`` (at most ``limit`` bonds)."""
    adjacency: dict[int, list[int]] = {}
    for bid, bond in sketch.bonds.items():
        if bid == skip:
            continue
        adjacency.setdefault(bond.begin, []).append(bond.end)
        adjacency.setdefault(bond.end, []).append(bond.begin)
    previous = {start: start}
    frontier = [start]
    for _ in range(limit):
        nxt = []
        for atom in frontier:
            for n in adjacency.get(atom, ()):
                if n not in previous:
                    previous[n] = atom
                    if n == goal:
                        path = [goal]
                        while path[-1] != start:
                            path.append(previous[path[-1]])
                        return path
                    nxt.append(n)
        frontier = nxt
    return []


def _dist(a: Point, b: Point) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _norm(angle: float) -> float:
    """The angle in [0, 2π)."""
    return angle % (2 * math.pi)
