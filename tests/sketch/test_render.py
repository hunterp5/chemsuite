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

"""Drawing conventions: labels, hydrogens, double-bond placement and wedges."""

from __future__ import annotations

import math

from chemsuite.sketch import render, valence
from chemsuite.sketch.convert import from_text
from chemsuite.sketch.model import Atom, Bond, Order, Sketch, Stereo


def _scene(smiles: str):
    sketch = from_text(smiles)
    return sketch, render.build(sketch, valence.check(sketch))


def _text(label: render.Label) -> str:
    return "".join(text for text, _ in label.runs)


def _labels(smiles: str) -> list[str]:
    return sorted(_text(label) for label in _scene(smiles)[1].labels.values())


def test_carbon_is_unlabelled_unless_alone_or_decorated():
    assert _labels("CCC") == []
    assert _labels("C") == ["CH4"]
    assert _labels("[13CH3]C")[0] in ("13CH3", "H313C")  # isotope stays on its element
    assert _labels("C[CH2+]") == ["CH2+"]
    assert _labels("C[CH2]") == ["CH2•"]


def test_heteroatoms_show_hydrogens_with_subscripts():
    sketch, scene = _scene("NCCO")
    by_symbol = {label.symbol: label for label in scene.labels.values()}
    assert by_symbol["N"].runs[-1] == ("2", "sub")
    assert _text(by_symbol["N"]) in ("NH2", "H2N")
    assert "H" in _text(by_symbol["O"])


def test_hydrogens_go_left_when_every_bond_points_right():
    sketch = Sketch()
    o, c = sketch.new_id(), sketch.new_id()
    sketch.atoms[o] = Atom(0.0, 0.0, "O")
    sketch.atoms[c] = Atom(1.0, 0.0)
    sketch.bonds[sketch.new_id()] = Bond(o, c)
    scene = render.build(sketch, valence.check(sketch))
    label = scene.labels[o]
    assert _text(label) == "HO"
    assert label.runs[label.centre] == ("O", "normal")


def test_charges_and_deuterium():
    assert _labels("C[N+](C)(C)C") == ["N+"]
    assert _labels("CC(=O)[O-]") == ["O", "O−"]
    assert _labels("[2H]C(Cl)(Cl)Cl") == ["Cl", "Cl", "Cl", "D"]
    assert _labels("[O-2]") == ["O2−"]


def _second_line(scene, bid):
    return [line for line in scene.lines if line.bond == bid]


def test_ring_double_bonds_are_drawn_inside_the_ring():
    sketch, scene = _scene("c1ccccc1")
    xs = [a.x for a in sketch.atoms.values()]
    ys = [a.y for a in sketch.atoms.values()]
    centre = (sum(xs) / len(xs), sum(ys) / len(ys))
    doubles = [bid for bid, b in sketch.bonds.items() if b.order is Order.DOUBLE]
    assert len(doubles) == 3
    for bid in doubles:
        bond = sketch.bonds[bid]
        outer, inner = sorted(
            _second_line(scene, bid),
            key=lambda line: (
                -math.dist(((line.a[0] + line.b[0]) / 2, (line.a[1] + line.b[1]) / 2), centre)
            ),
        )
        assert math.dist(inner.a, inner.b) < math.dist(outer.a, outer.b)  # pulled in at the ends
        a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
        assert math.dist(outer.a, a) < 1e-9 and math.dist(outer.b, b) < 1e-9


def test_symmetric_double_bonds_are_centred():
    """Acetone's C=O and ethylene: two equal lines either side of the bond axis."""
    for smiles in ("CC(C)=O", "C=C"):
        sketch, scene = _scene(smiles)
        (bid,) = [b for b, bond in sketch.bonds.items() if bond.order is Order.DOUBLE]
        lines = _second_line(scene, bid)
        assert len(lines) == 2
        lengths = [math.dist(line.a, line.b) for line in lines]
        assert abs(lengths[0] - lengths[1]) < 1e-9


def test_terminal_double_bond_with_one_substituent_is_offset_towards_it():
    sketch, scene = _scene("C=CC")
    (bid,) = [b for b, bond in sketch.bonds.items() if bond.order is Order.DOUBLE]
    lengths = sorted(math.dist(line.a, line.b) for line in _second_line(scene, bid))
    assert lengths[0] < lengths[1]


def test_aromatic_bonds_have_a_dashed_inner_line():
    sketch = Sketch()
    ids = []
    for i in range(6):
        t = math.pi / 3 * i
        aid = sketch.new_id()
        sketch.atoms[aid] = Atom(math.cos(t), math.sin(t))
        ids.append(aid)
    for i in range(6):
        sketch.bonds[sketch.new_id()] = Bond(ids[i], ids[(i + 1) % 6], Order.AROMATIC)
    scene = render.build(sketch, valence.check(sketch))
    assert sum(line.dashed for line in scene.lines) == 6


def test_wedge_tip_sits_at_the_begin_atom_and_widens():
    sketch = Sketch()
    a, b = sketch.new_id(), sketch.new_id()
    sketch.atoms[a] = Atom(0.0, 0.0)
    sketch.atoms[b] = Atom(1.0, 0.0)
    sketch.bonds[sketch.new_id()] = Bond(a, b, Order.SINGLE, Stereo.WEDGE)
    scene = render.build(sketch, valence.check(sketch))
    (wedge,) = scene.wedges
    assert wedge.tip == (0.0, 0.0)
    assert math.dist(*wedge.base) == render.WEDGE_WIDTH
    sketch.bonds = {9: Bond(a, b, Order.SINGLE, Stereo.HASH)}
    scene = render.build(sketch, valence.check(sketch))
    widths = [math.dist(line.a, line.b) for line in scene.lines]
    assert widths == sorted(widths) and len(widths) >= 3


def test_bonds_stop_short_of_labels():
    sketch, scene = _scene("CO")
    (line,) = scene.lines
    o = next(a for a in sketch.atoms.values() if a.symbol == "O")
    assert min(math.dist(line.a, o.pos), math.dist(line.b, o.pos)) >= render.LABEL_GAP - 1e-9


def test_problem_atoms_are_labelled_with_the_reason():
    sketch, scene = _scene("CC(C)(C)C")
    centre = next(a for a in sketch.atoms if sketch.degree(a) == 4)
    from chemsuite.sketch import edits

    edits.bond_from_atom(edits.Change(sketch, "x"), centre)
    scene = render.build(sketch, valence.check(sketch))
    assert scene.labels[centre].problem.startswith("Too many bonds")
