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

"""Where new atoms go (SKETCHER_PLAN §4.2), rule by rule from Ketcher at 2b031cc2.

Angles are checked to 0.1 degree. Expectations are stated as they look on the page (y up).
"""

from __future__ import annotations

import itertools
import math

import pytest

from chemsuite.sketch import geometry
from chemsuite.sketch.edits import History, bond_from_atom, bond_on_canvas, chain
from chemsuite.sketch.model import Atom, Order, Sketch
from tests.sketch.helpers import degrees


def _star(*angles_deg: float) -> tuple[Sketch, int]:
    """A centre atom with one bond at each given angle."""
    sketch = Sketch()
    centre = sketch.new_id()
    sketch.atoms[centre] = Atom(0.0, 0.0)
    history = History(sketch)
    change = history.change("star")
    for a in angles_deg:
        end = change.add_atom(Atom(math.cos(math.radians(a)), math.sin(math.radians(a))))
        change.add_bond(centre, end)
    history.commit(change)
    return sketch, centre


def _min_separation(sketch: Sketch) -> float:
    pts = [a.pos for a in sketch.atoms.values()]
    return min(math.dist(p, q) for p, q in itertools.combinations(pts, 2))


# -- fracAngle -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "snapped"), [(0, 0), (7, 0), (8, 15), (22, 15), (23, 30), (-8, -15), (181, 180)]
)
def test_angles_snap_to_15_degrees(raw, snapped):
    assert math.degrees(geometry.snap_angle(math.radians(raw))) == pytest.approx(snapped)


# -- bond tool on empty canvas ---------------------------------------------------------------


def test_first_single_bond_rises_30_degrees_others_are_horizontal():
    assert degrees(geometry.first_bond_angle(Order.SINGLE)) == pytest.approx(30)
    assert degrees(geometry.first_bond_angle(Order.DOUBLE)) == pytest.approx(0)
    assert degrees(geometry.first_bond_angle(Order.TRIPLE)) == pytest.approx(0)


# -- atomForNewBond --------------------------------------------------------------------------


def test_lone_atom_grows_horizontally():
    sketch = Sketch()
    aid = sketch.new_id()
    sketch.atoms[aid] = Atom(0, 0)
    assert degrees(geometry.new_bond_angle(sketch, aid)) == pytest.approx(0)


def test_one_neighbour_gives_120_degrees():
    sketch, centre = _star(180)
    angle = degrees(geometry.new_bond_angle(sketch, centre))
    assert angle == pytest.approx(300, abs=0.1) or angle == pytest.approx(60, abs=0.1)


@pytest.mark.parametrize(
    ("existing", "expected"),
    [
        ((0, 120), 240),  # the 240 degree gap is bisected
        ((90, 210, 330), 30),  # three bonds 120 apart: first of the equal gaps, bisected
        ((0, 90), 225),
        ((0, 90, 180), 270),
        ((10, 20, 30), 200),
    ],
)
def test_largest_gap_is_bisected(existing, expected):
    sketch, centre = _star(*existing)
    assert degrees(geometry.new_bond_angle(sketch, centre)) == pytest.approx(expected, abs=0.1)


def test_chain_built_by_clicking_zig_zags():
    """Each click on the newest atom continues a trans zig-zag, never folding back."""
    history = History()
    change = history.change("draw")
    bid = bond_on_canvas(change, (0, 0))
    tip = history.sketch.bonds[bid].end
    tips = [history.sketch.bonds[bid].begin, tip]
    for _ in range(6):
        nb = bond_from_atom(change, tip)
        tip = history.sketch.bonds[nb].end
        tips.append(tip)
    pts = [history.sketch.atoms[t].pos for t in tips]
    for p, q, r in zip(pts, pts[1:], pts[2:], strict=False):
        v1 = (q[0] - p[0], q[1] - p[1])
        v2 = (r[0] - q[0], r[1] - q[1])
        turn = math.degrees(
            math.atan2(v1[0] * v2[1] - v1[1] * v2[0], v1[0] * v2[0] + v1[1] * v2[1])
        )
        assert abs(turn) == pytest.approx(60, abs=0.1)  # 120 degree bond angle
    turns = []
    for p, q, r in zip(pts, pts[1:], pts[2:], strict=False):
        cross = (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0])
        turns.append(cross > 0)
    assert all(a != b for a, b in itertools.pairwise(turns))  # alternating: zig-zag
    assert _min_separation(history.sketch) > 0.99
    # A first single bond at +30 then zig-zag: every other atom on the same baseline.
    assert {round(p[1], 6) for p in pts[::2]} == {0.0}


@pytest.mark.parametrize(
    ("first", "second"),
    [
        (Order.TRIPLE, Order.SINGLE),  # alkyne: C-C#C-C is linear
        (Order.SINGLE, Order.TRIPLE),
        (Order.DOUBLE, Order.DOUBLE),  # allene: C=C=C is linear
        (Order.TRIPLE, Order.TRIPLE),
    ],
)
def test_alkynes_and_allenes_continue_straight(first, second):
    history = History()
    change = history.change("draw")
    bid = bond_on_canvas(change, (0, 0), first)
    end = history.sketch.bonds[bid].end
    angle = geometry.new_bond_angle(history.sketch, end, second)
    assert degrees(angle) == pytest.approx(degrees(geometry.first_bond_angle(first)), abs=0.1)


def test_double_after_single_is_not_straight():
    history = History()
    change = history.change("draw")
    bid = bond_on_canvas(change, (0, 0), Order.SINGLE)
    end = history.sketch.bonds[bid].end
    turn = degrees(geometry.new_bond_angle(history.sketch, end, Order.DOUBLE)) - 30
    assert abs(((turn + 180) % 360) - 180) == pytest.approx(60, abs=0.1)


def test_new_bond_joins_an_atom_already_there():
    """Ketcher reuses an atom within 0.1 bond lengths of the new position."""
    sketch, centre = _star(180)
    angle = geometry.new_bond_angle(sketch, centre)
    target = geometry.along((0, 0), angle)
    other = sketch.new_id()
    sketch.atoms[other] = Atom(target[0] + 0.05, target[1])
    _, existing = geometry.new_bond_target(sketch, centre)
    assert existing == other
    history = History(sketch)
    change = history.change("bond")
    bid = bond_from_atom(change, centre)
    assert other in (sketch.bonds[bid].begin, sketch.bonds[bid].end)
    assert len(sketch.atoms) == 3


def test_atom_beyond_merge_distance_is_not_reused():
    sketch, centre = _star(180)
    target = geometry.along((0, 0), geometry.new_bond_angle(sketch, centre))
    other = sketch.new_id()
    sketch.atoms[other] = Atom(target[0] + 0.15, target[1])
    assert geometry.new_bond_target(sketch, centre)[1] is None


def test_filling_a_carbon_gives_a_tetrahedral_style_cross():
    sketch, centre = _star(180)
    history = History(sketch)
    change = history.change("fill")
    for _ in range(3):
        bond_from_atom(change, centre)
    angles = sorted(
        degrees(geometry.angle_of((0, 0), sketch.atoms[n].pos)) for n in sketch.neighbours(centre)
    )
    gaps = [(b - a) % 360 for a, b in itertools.pairwise([*angles, angles[0] + 360])]
    assert min(gaps) >= 60 - 0.1  # no two bonds closer than 60 degrees
    assert _min_separation(sketch) > 0.99


# -- chains ----------------------------------------------------------------------------------


@pytest.mark.parametrize(("drag", "bonds"), [(0.2, 1), (1.0, 1), (1.01, 2), (3.5, 4), (6, 6)])
def test_chain_adds_one_bond_per_bond_length_begun(drag, bonds):
    assert geometry.chain_length(drag) == bonds


def test_chain_points_are_a_regular_zig_zag_along_the_drag():
    for drag_deg in (0, 30, 45, 90, 200):
        a = math.radians(drag_deg)
        pts = [(0.0, 0.0), *geometry.chain_points((0.0, 0.0), a, 5)]
        for p, q in itertools.pairwise(pts):
            assert math.dist(p, q) == pytest.approx(1.0)
        for p, q, r in zip(pts, pts[1:], pts[2:], strict=False):
            v1 = (p[0] - q[0], p[1] - q[1])
            v2 = (r[0] - q[0], r[1] - q[1])
            cos = (v1[0] * v2[0] + v1[1] * v2[1]) / (math.hypot(*v1) * math.hypot(*v2))
            assert math.degrees(math.acos(cos)) == pytest.approx(120, abs=0.1)
        # Every second atom lies on the drag line.
        for p in pts[2::2]:
            assert math.sin(geometry.angle_of(pts[0], p) - a) == pytest.approx(0, abs=1e-9)


def test_chain_first_atom_is_clockwise_of_the_drag_as_in_ketcher():
    (p,) = geometry.chain_points((0.0, 0.0), 0.0, 1)
    assert p[1] < 0  # dragging right, the first atom sits below the line


def test_chain_edit_from_an_atom_and_from_canvas():
    history = History()
    change = history.change("chain")
    atoms = chain(change, (0.0, 0.0), 0.0, 4)
    assert len(atoms) == 5 and len(history.sketch.bonds) == 4
    more = chain(change, atoms[-1], math.radians(90), 2)
    assert more[0] == atoms[-1] and len(history.sketch.atoms) == 7
