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

"""Rings and templates: placing them on empty canvas, on an atom, or fused onto a bond.

The placement rules follow Ketcher at 2b031cc2 (Apache-2.0):
``ketcher-core/src/application/editor/actions/template.ts`` (``fromTemplateOnCanvas``,
``fromTemplateOnAtom``, ``fromTemplateOnBond``) and the template tool
``ketcher-react/src/script/editor/tool/template.ts``.

- **On an atom** the template's attachment atom *becomes* that atom (a spiro or ring-end
  attachment) and the template is turned so its centre lies in the direction a new bond would
  take (``atomForNewBond``). Dragging more than one bond length away instead attaches the
  template through a new bond (a substituent).
- **On a bond** the template's attachment bond is laid onto that bond and the template goes on
  the side away from the rest of the molecule: away from the ring the bond is in, else away
  from more of its neighbours. A drag picks the side the pointer is on.
- New atoms landing within 0.1 bond lengths of existing atoms merge with them, so fusing onto
  a ring shares the right atoms; a shared bond keeps the higher of the two orders.

One rule differs from Ketcher's, which hard-codes Kekulé patterns for fusing benzene (issue
#1855): after an aromatic template is fused, the whole aromatic ring system it joined is
Kekulé-assigned again by RDKit, so naphthalene, anthracene or phenanthrene come out valid
whichever bond was clicked.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rdkit import Chem, rdBase

from chemsuite.sketch import convert, geometry
from chemsuite.sketch.edits import Change
from chemsuite.sketch.model import Atom, Bond, Order, Sketch

Point = tuple[float, float]
MERGE = geometry.MERGE_DISTANCE


@dataclass(frozen=True)
class Template:
    """A drawing to place: atoms (x, y, element), bonds (i, j, order), attachment atom and
    bond (indexes), and whether its ring is aromatic (drawn Kekulé, re-assigned on fusing)."""

    name: str
    atoms: tuple[tuple[float, float, str], ...]
    bonds: tuple[tuple[int, int, Order], ...]
    attach_atom: int = 0
    attach_bond: int = 0
    aromatic: bool = False

    @property
    def centre(self) -> Point:
        n = len(self.atoms)
        return (sum(a[0] for a in self.atoms) / n, sum(a[1] for a in self.atoms) / n)

    @property
    def angle0(self) -> float:
        """Direction from the attachment atom to the centre (Ketcher's ``angle0``)."""
        a = self.atoms[self.attach_atom]
        return geometry.angle_of((a[0], a[1]), self.centre)


def ring(size: int, aromatic: bool = False, name: str = "") -> Template:
    """A regular ring of unit bonds with a vertex at the top (vertical sides for six)."""
    radius = 1 / (2 * math.sin(math.pi / size))
    atoms = tuple(
        (
            radius * math.cos(math.pi / 2 + 2 * math.pi * i / size),
            radius * math.sin(math.pi / 2 + 2 * math.pi * i / size),
            "C",
        )
        for i in range(size)
    )
    bonds = tuple(
        (
            i,
            (i + 1) % size,
            Order.DOUBLE if aromatic and i % 2 == 0 else Order.SINGLE,
        )
        for i in range(size)
    )
    return Template(name or f"ring{size}", atoms, bonds, aromatic=aromatic)


BENZENE = ring(6, aromatic=True, name="Benzene")
#: The ring bar, as in Ketcher's default template bar.
RING_BAR = (
    BENZENE,
    ring(3, name="Cyclopropane"),
    ring(4, name="Cyclobutane"),
    ring(5, name="Cyclopentane"),
    ring(6, name="Cyclohexane"),
    ring(7, name="Cycloheptane"),
    ring(8, name="Cyclooctane"),
)


def _transform(points, origin: Point, angle: float, scale: float, at: Point) -> list[Point]:
    c, s = math.cos(angle), math.sin(angle)
    out = []
    for x, y in points:
        dx, dy = (x - origin[0]) * scale, (y - origin[1]) * scale
        out.append((at[0] + dx * c - dy * s, at[1] + dx * s + dy * c))
    return out


def _place(change: Change, template: Template, positions: list[Point], fixed: dict[int, int]):
    """Add the template's atoms at ``positions`` (reusing ``fixed`` and any atom within 0.1)
    and its bonds; returns (template index -> atom id, ids of the bonds it touched)."""
    sketch = change.sketch
    mapping = dict(fixed)
    for i, (x, y) in enumerate(positions):
        if i in mapping:
            continue
        existing = sketch.atom_at(x, y, within=MERGE, exclude=set(mapping.values()))
        if existing is not None:
            mapping[i] = existing
        else:
            mapping[i] = change.add_atom(Atom(x, y, template.atoms[i][2]))
    touched = []
    for i, j, order in template.bonds:
        a, b = mapping[i], mapping[j]
        if a == b:
            continue
        bid = sketch.bond_between(a, b)
        if bid is None:
            touched.append(change.add_bond(a, b, order))
        else:
            old = sketch.bonds[bid]
            if order.value > old.order.value and old.order is not Order.AROMATIC:
                change.set_bond(bid, Bond(old.begin, old.end, order, old.stereo))
            touched.append(bid)
    return mapping, touched


def on_canvas(change: Change, template: Template, centre: Point, angle: float = 0.0) -> list[int]:
    """The template centred on ``centre``, turned by ``angle``. Returns its atom ids."""
    points = [(a[0], a[1]) for a in template.atoms]
    positions = _transform(points, template.centre, angle, 1.0, centre)
    mapping, _ = _place(change, template, positions, {})
    return [mapping[i] for i in range(len(template.atoms))]


def on_atom(
    change: Change,
    template: Template,
    atom: int,
    angle: float | None = None,
    extra_bond: bool = False,
) -> list[int]:
    """Attach at an atom: sharing it, or through a new bond with ``extra_bond``."""
    sketch = change.sketch
    if angle is None:
        angle = geometry.new_bond_angle(sketch, atom)
    anchor = atom
    if extra_bond:
        end = geometry.along(sketch.atoms[atom].pos, angle)
        existing = sketch.atom_at(*end, within=MERGE, exclude={atom})
        anchor = existing if existing is not None else change.add_atom(Atom(*end))
        change.add_bond(atom, anchor)
    points = [(a[0], a[1]) for a in template.atoms]
    origin = points[template.attach_atom]
    positions = _transform(points, origin, angle - template.angle0, 1.0, sketch.atoms[anchor].pos)
    aromatic_before = _aromatic_bonds(sketch) if template.aromatic else set()
    mapping, touched = _place(change, template, positions, {template.attach_atom: anchor})
    if template.aromatic:
        _kekulize_around(change, touched, aromatic_before)
    return [mapping[i] for i in range(len(template.atoms))]


def fuse_side(sketch: Sketch, bond_id: int) -> int:
    """+1 to fuse on the left of begin→end, -1 on the right: away from the molecule."""
    bond = sketch.bonds[bond_id]
    a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
    centre = geometry.ring_centres(sketch).get(bond_id)
    if centre is not None:
        return -1 if geometry.side_of(a, b, centre) > 0 else 1
    left = right = 0
    for end in (bond.begin, bond.end):
        for n in sketch.neighbours(end):
            if n in (bond.begin, bond.end):
                continue
            s = geometry.side_of(a, b, sketch.atoms[n].pos)
            left += s > 1e-6
            right += s < -1e-6
    return -1 if left > right else 1


def on_bond(change: Change, template: Template, bond_id: int, side: int | None = None) -> list[int]:
    """Fuse the template onto a bond, on ``side`` (+1 left of begin→end, -1 right)."""
    sketch = change.sketch
    side = side if side is not None else fuse_side(sketch, bond_id)
    bond = sketch.bonds[bond_id]
    ti, tj, _ = template.bonds[template.attach_bond]
    points = [(a[0], a[1]) for a in template.atoms]
    a, b = sketch.atoms[bond.begin].pos, sketch.atoms[bond.end].pos
    # Lay template bond ti->tj onto begin->end; if its centre lands on the wrong side, lay
    # tj->ti instead (Ketcher's shouldFlip).
    for first, second in ((ti, tj), (tj, ti)):
        t_angle = geometry.angle_of(points[first], points[second])
        scale = math.dist(a, b) / math.dist(points[first], points[second])
        delta = geometry.angle_of(a, b) - t_angle
        positions = _transform(points, points[first], delta, scale, a)
        centre = (
            sum(p[0] for p in positions) / len(positions),
            sum(p[1] for p in positions) / len(positions),
        )
        if (geometry.side_of(a, b, centre) > 0) == (side > 0):
            break
    fixed = {first: bond.begin, second: bond.end}
    aromatic_before = _aromatic_bonds(sketch) if template.aromatic else set()
    mapping, touched = _place(change, template, positions, fixed)
    if template.aromatic:
        _kekulize_around(change, touched, aromatic_before)
    return [mapping[i] for i in range(len(template.atoms))]


def _kekulize_around(change: Change, new_bonds: list[int], aromatic_before: set[int]) -> None:
    """Give the aromatic ring system that the new ring joined a valid Kekulé structure.

    The new bonds and every bond RDKit already finds aromatic in the ring system they touch
    are made aromatic, RDKit Kekulé-assigns the molecule, and those bonds take its orders. If
    RDKit cannot (the drawing is invalid elsewhere), nothing changes.
    """
    sketch = change.sketch
    index = {aid: i for i, aid in enumerate(convert.atom_order(sketch))}
    system = set(new_bonds)
    aromatic_now = {bid for bid in aromatic_before if bid in sketch.bonds}
    grown = True
    while grown:
        atoms = {x for bid in system for x in (sketch.bonds[bid].begin, sketch.bonds[bid].end)}
        extra = {
            bid
            for bid in aromatic_now - system
            if sketch.bonds[bid].begin in atoms or sketch.bonds[bid].end in atoms
        }
        grown = bool(extra)
        system |= extra
    trial = sketch.copy()
    trial.bonds = dict(sketch.bonds)
    for bid in system:
        old = trial.bonds[bid]
        trial.bonds[bid] = Bond(old.begin, old.end, Order.AROMATIC, old.stereo)
    try:
        mol = convert.to_rdkit(trial)
        with rdBase.BlockLogs():
            Chem.Kekulize(mol, clearAromaticFlags=True)
    except (convert.SketchError, ValueError, RuntimeError):
        return
    for bid in system:
        old = sketch.bonds[bid]
        rd = mol.GetBondBetweenAtoms(index[old.begin], index[old.end])
        order = Order.DOUBLE if rd.GetBondType() == Chem.BondType.DOUBLE else Order.SINGLE
        if order is not old.order:
            change.set_bond(bid, Bond(old.begin, old.end, order, old.stereo))


def _aromatic_bonds(sketch: Sketch) -> set[int]:
    """Bonds RDKit perceives as aromatic (none if the drawing is invalid)."""
    index = {aid: i for i, aid in enumerate(convert.atom_order(sketch))}
    try:
        mol = convert.to_rdkit(sketch)
    except convert.SketchError:
        return set()
    out = set()
    for bid, bond in sketch.bonds.items():
        rd = mol.GetBondBetweenAtoms(index[bond.begin], index[bond.end])
        if rd is not None and rd.GetIsAromatic():
            out.add(bid)
    return out
