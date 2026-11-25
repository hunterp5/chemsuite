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
"""What a protein view draws, as GPU-ready rows, in two parts that change at different rates.

- :class:`Geometry`: one sphere per shown atom (``x y z radius``) and one cylinder per bond
  between shown atoms (``x0 y0 z0 radius x1 y1 z1 0``); in sticks and ball and stick, a double
  or triple bond is two or three thinner ones side by side in the plane of a neighbouring atom,
  and aromatic rings are drawn in a Kekulé form (:func:`drawn_orders`), as the table's 3D view
  draws them (:mod:`chemsuite.molview.geometry`). Rebuilt when atoms move or the styles change.
- The *look*: per sphere ``r g b a flags`` and per cylinder each end's colour and flags
  (:func:`sphere_look`, :func:`cylinder_look`). Rebuilt on a colour or selection change, which
  then uploads only these rows; hovering rewrites one atom's flags.

Each atom has a style (:data:`STYLES`); hidden atoms are simply left out. ``flags`` add up
:data:`SELECTED`, :data:`HOVER` and :data:`GHOST` (the shaders decode them with float maths:
older GLSL has no bit operations).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemsuite.molview.elements import vdw_radii
from chemsuite.molview.geometry import ALL_H, NO_H, POLAR_H

from ..structure.model import ADDITIVE, AMINO, AROMATIC, GLYCAN, ION, LIGAND, NUCLEIC, Structure

HIDDEN = 0
WIRE = 1
STICK = 2
BALL_AND_STICK = 3
SPACEFILL = 4
STYLES = {WIRE: "Wire", STICK: "Stick", BALL_AND_STICK: "Ball and Stick", SPACEFILL: "Spacefill"}

SELECTED = 1.0
HOVER = 2.0
GHOST = 4.0

SPHERE_FLOATS = 4
CYLINDER_FLOATS = 8
SPHERE_LOOK_FLOATS = 5
CYLINDER_LOOK_FLOATS = 10

#: bond radius per style (0: no bonds), and atom radius (a scale of the van der Waals radius
#: when the second number is set, else the first)
_BOND = np.array([0.0, 0.06, 0.15, 0.1, 0.0], np.float32)
_ATOM_FIXED = np.array([0.0, 0.06, 0.15, 0.0, 0.0], np.float32)
_ATOM_SCALE = np.array([0.0, 0.0, 0.0, 0.25, 1.0], np.float32)
#: a double or triple bond's cylinders: their radius as a share of the bond's, and their offsets
#: across it in bond radii (as :mod:`chemsuite.molview.geometry` draws them)
_DOUBLE = (0.6, (-0.9, 0.9))
_TRIPLE = (0.45, (-1.25, 0.0, 1.25))
#: styles whose bonds show their order (a wire stays one line)
_ORDERED = (STICK, BALL_AND_STICK)
#: an aromatic system larger than this keeps its bonds single (no Kekulé form is searched)
_KEKULE_MAX = 40
#: an ion shown as a stick or wire would vanish: it gets a ball
_ION_SCALE = 0.4
#: pick radius at least this, so thin wires can be clicked
MIN_PICK_RADIUS = 0.5


@dataclass(frozen=True)
class Geometry:
    atoms: np.ndarray  # int64 (s,): the atom each sphere draws
    spheres: np.ndarray  # float32 (s, 4)
    bonds: np.ndarray  # int64 (b,): the structure bond each cylinder draws
    ends: np.ndarray  # int64 (b, 2): its two atoms
    cylinders: np.ndarray  # float32 (b, 8)

    @property
    def pick_radii(self) -> np.ndarray:
        return np.maximum(self.spheres[:, 3], MIN_PICK_RADIUS)

    @classmethod
    def empty(cls) -> Geometry:
        return cls(
            np.zeros(0, np.int64),
            np.zeros((0, SPHERE_FLOATS), np.float32),
            np.zeros(0, np.int64),
            np.zeros((0, 2), np.int64),
            np.zeros((0, CYLINDER_FLOATS), np.float32),
        )


#: residues this close to a ligand are the binding site, shown as sticks with a cartoon
SITE = 5.0


def default_styles(
    structure: Structure, active: np.ndarray | None = None, cartoon: bool = False
) -> np.ndarray:
    """A first look at a complex: the ligands as sticks, ions as balls, additives and
    glycans as sticks; the protein as wire, or with ``cartoon`` hidden but for the binding
    site's residues as sticks; waters and non-polar hydrogens hidden; of alternate locations
    only the shown one."""
    kind = structure.atom_kind
    styles = np.full(structure.n_atoms, HIDDEN, np.uint8)
    polymer = np.isin(kind, (AMINO, NUCLEIC))
    if cartoon:
        styles[polymer & binding_site(structure)] = STICK
    else:
        styles[polymer] = WIRE
    styles[kind == LIGAND] = STICK
    styles[kind == ION] = BALL_AND_STICK
    styles[np.isin(kind, (ADDITIVE, GLYCAN))] = STICK
    styles[nonpolar_hydrogens(structure)] = HIDDEN
    if active is None:
        active = structure.active()
    styles[~active] = HIDDEN
    return styles


def binding_site(structure: Structure, distance: float = SITE) -> np.ndarray:
    """Atoms of the residues with an atom within ``distance`` of a ligand."""
    from ..structure.select import within

    ligand = structure.atom_kind == LIGAND
    near = within(structure.xyz(0), ligand, distance)
    residues = np.zeros(structure.n_residues, bool)
    residues[structure.residue_of[near]] = True
    return residues[structure.residue_of]


def nonpolar_hydrogens(structure: Structure) -> np.ndarray:
    """Hydrogens bonded to carbon."""
    hydrogen = structure.element == 1
    out = np.zeros(structure.n_atoms, bool)
    if not hydrogen.any() or not len(structure.bonds):
        return out
    a, b = structure.bonds[:, 0], structure.bonds[:, 1]
    elements = structure.element
    out[a[hydrogen[a] & (elements[b] == 6)]] = True
    out[b[hydrogen[b] & (elements[a] == 6)]] = True
    return out


def hydrogen_parents(structure: Structure) -> np.ndarray:
    """The heavy atom each hydrogen is bonded to (-1 for none, and for other atoms)."""
    out = np.full(structure.n_atoms, -1, np.int64)
    if not len(structure.bonds):
        return out
    a, b = structure.bonds[:, 0], structure.bonds[:, 1]
    hydrogen = structure.element == 1
    first, second = hydrogen[a] & ~hydrogen[b], hydrogen[b] & ~hydrogen[a]
    out[a[first]] = b[first]
    out[b[second]] = a[second]
    return out


def shown_styles(
    structure: Structure, styles: np.ndarray, active: np.ndarray, hydrogens: str = POLAR_H
) -> np.ndarray:
    """The style each atom is drawn in. A hydrogen takes its heavy atom's, so it follows
    whatever that atom is shown as; then ``hydrogens`` (all, polar, none: the Structure
    Viewer's choices) hides none, the non-polar, or all of them. Of alternate locations only
    the shown one is drawn."""
    out = np.array(styles, np.uint8)
    hydrogen = np.flatnonzero(structure.element == 1)
    if len(hydrogen):
        parent = hydrogen_parents(structure)[hydrogen]
        bonded = parent >= 0
        out[hydrogen[bonded]] = styles[parent[bonded]]
        if hydrogens == NO_H:
            out[hydrogen] = HIDDEN
        elif hydrogens != ALL_H:
            out[nonpolar_hydrogens(structure)] = HIDDEN
    out[~active] = HIDDEN
    return out


def geometry(structure: Structure, styles: np.ndarray, model: int = 0) -> Geometry:
    xyz = structure.xyz(model)
    atoms = np.flatnonzero(styles != HIDDEN)
    style = styles[atoms]
    radius = np.where(
        _ATOM_SCALE[style] > 0,
        vdw_radii(structure.element[atoms]) * _ATOM_SCALE[style],
        _ATOM_FIXED[style],
    )
    ion = structure.atom_kind[atoms] == ION
    radius[ion] = np.maximum(radius[ion], vdw_radii(structure.element[atoms[ion]]) * _ION_SCALE)
    spheres = np.empty((len(atoms), SPHERE_FLOATS), np.float32)
    spheres[:, :3] = xyz[atoms]
    spheres[:, 3] = radius
    i, j = structure.bonds[:, 0], structure.bonds[:, 1]
    bond_radius = np.minimum(_BOND[styles[i]], _BOND[styles[j]])
    shown = np.flatnonzero((styles[i] != HIDDEN) & (styles[j] != HIDDEN) & (bond_radius > 0))
    order = np.ones(len(shown), np.int8)
    ordered = np.isin(styles[i[shown]], _ORDERED) & np.isin(styles[j[shown]], _ORDERED)
    if ordered.any():
        order[ordered] = drawn_orders(structure)[shown[ordered]]
    drawn, step, scale = [shown], [np.zeros(len(shown))], [np.ones(len(shown))]
    for n, (share, offsets) in ((2, _DOUBLE), (3, _TRIPLE)):
        these = shown[order == n]
        drawn[0] = drawn[0][np.isin(drawn[0], these, invert=True)]
        step[0], scale[0] = np.zeros(len(drawn[0])), np.ones(len(drawn[0]))
        for offset in offsets:
            drawn.append(these)
            step.append(np.full(len(these), offset))
            scale.append(np.full(len(these), share))
    bonds = np.concatenate(drawn).astype(np.int64)
    steps, radius = np.concatenate(step), bond_radius[bonds] * np.concatenate(scale)
    keep = np.argsort(bonds, kind="stable")  # each bond's cylinders together
    bonds, steps, radius = bonds[keep], steps[keep], radius[keep]
    shift = np.zeros((len(bonds), 3))
    multi = steps != 0
    if multi.any():
        across = _sides(structure, xyz, bonds[multi])
        shift[multi] = across * (steps[multi] * bond_radius[bonds[multi]])[:, None]
    cylinders = np.zeros((len(bonds), CYLINDER_FLOATS), np.float32)
    cylinders[:, :3] = xyz[i[bonds]] + shift
    cylinders[:, 3] = radius
    cylinders[:, 4:7] = xyz[j[bonds]] + shift
    ends = np.stack([i[bonds], j[bonds]], axis=1).astype(np.int64)
    return Geometry(atoms, spheres, bonds, ends, cylinders)


def _sides(structure: Structure, xyz: np.ndarray, bonds: np.ndarray) -> np.ndarray:
    """A unit vector across each of *bonds* (structure bond indices), in the plane of a
    neighbouring heavy atom when there is one, else any perpendicular."""
    s = structure
    heavy = s.element != 1
    neighbours: dict[int, list[int]] = {int(a): [] for a in np.unique(s.bonds[bonds])}
    for a, b in s.bonds.tolist():
        if a in neighbours and heavy[b]:
            neighbours[a].append(b)
        if b in neighbours and heavy[a]:
            neighbours[b].append(a)
    out = np.zeros((len(bonds), 3))
    for row, (a, b) in enumerate(s.bonds[bonds].tolist()):
        axis = xyz[b] - xyz[a]
        axis = axis / max(float(np.linalg.norm(axis)), 1e-6)
        out[row] = _across(xyz, axis, neighbours, a, b)
    return out


def _across(xyz, axis, neighbours, a: int, b: int) -> np.ndarray:
    for here, other in ((a, b), (b, a)):
        for k in neighbours[here]:
            if k != other:
                v = xyz[k] - xyz[here]
                v = v - axis * float(v @ axis)
                norm = float(np.linalg.norm(v))
                if norm > 1e-3:
                    return v / norm
    guess = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = np.cross(axis, guess)
    return v / np.linalg.norm(v)


def drawn_orders(structure: Structure) -> np.ndarray:
    """Each bond's order as drawn (int8: 1, 2 or 3), aromatic bonds in a Kekulé form. The
    double bonds go on the ring atoms that can take one: not an O, S or Se, an atom already
    double-bonded outside the ring (a pyridone's C=O), or an N with a hydrogen or a third
    neighbour (as in pyrrole), unless it is charged. Where an odd system leaves an atom
    without one, it is an N if it can be. A system over :data:`_KEKULE_MAX` atoms stays single."""
    s = structure
    order = np.where(np.isin(s.order, (2, 3)), s.order, 1).astype(np.int8)
    aromatic = np.flatnonzero(s.order == AROMATIC)
    if not len(aromatic):
        return order
    ring: dict[int, list[tuple[int, int]]] = {}  # atom: (neighbour, aromatic bond)
    for k, (a, b) in zip(aromatic.tolist(), s.bonds[aromatic].tolist(), strict=True):
        ring.setdefault(a, []).append((b, k))
        ring.setdefault(b, []).append((a, k))
    a, b = s.bonds[:, 0], s.bonds[:, 1]
    degree = np.bincount(s.bonds.ravel(), minlength=s.n_atoms)
    hydrogen = s.element == 1
    hydrogens = np.bincount(a[hydrogen[b]], minlength=s.n_atoms) + np.bincount(
        b[hydrogen[a]], minlength=s.n_atoms
    )
    doubled = np.zeros(s.n_atoms, bool)
    multiple = np.isin(s.order, (2, 3))
    doubled[a[multiple]] = doubled[b[multiple]] = True

    def takes(atom: int) -> bool:
        element = int(s.element[atom])
        if element in (8, 16, 34) or doubled[atom]:
            return False
        if element == 7 and s.charge[atom] <= 0:
            return bool(hydrogens[atom] == 0 and degree[atom] < 3)
        return True

    seen: set[int] = set()
    for start in ring:
        if start in seen:
            continue
        component, todo = [], [start]
        seen.add(start)
        while todo:
            atom = todo.pop()
            component.append(atom)
            for other, _k in ring[atom]:
                if other not in seen:
                    seen.add(other)
                    todo.append(other)
        if len(component) <= _KEKULE_MAX:
            order[_kekule([x for x in component if takes(x)], ring, s.element)] = 2
    return order


def _kekule(atoms: list[int], ring, element: np.ndarray) -> list[int]:
    """The aromatic bonds between *atoms* (those of one ring system that take a double bond)
    that are double in its best Kekulé form: the most double bonds, then the most of the atoms
    left without one nitrogens. An exhaustive search, pruned; ring systems are small."""
    allowed = set(atoms)
    atoms = sorted(atoms, key=lambda x: int(element[x]) == 7)  # carbons matched first
    best: tuple[int, int, list[int]] = (-1, -1, [])
    used: set[int] = set()
    chosen: list[int] = []
    skipped: list[int] = []

    def search(at: int) -> None:
        nonlocal best
        while at < len(atoms) and atoms[at] in used:
            at += 1
        if at == len(atoms):
            left = sum(1 for x in skipped if int(element[x]) == 7)
            if (len(chosen), left) > best[:2]:
                best = (len(chosen), left, list(chosen))
            return
        if len(chosen) + (len(atoms) - at) // 2 < best[0]:
            return  # cannot match as many as the best found
        atom = atoms[at]
        used.add(atom)
        for other, k in ring[atom]:
            if other in allowed and other not in used:
                used.add(other)
                chosen.append(k)
                search(at + 1)
                chosen.pop()
                used.discard(other)
        skipped.append(atom)  # this atom left without a double bond
        search(at + 1)
        skipped.pop()
        used.discard(atom)

    search(0)
    return best[2]


def sphere_look(geometry: Geometry, colours: np.ndarray, flags: np.ndarray) -> np.ndarray:
    look = np.empty((len(geometry.atoms), SPHERE_LOOK_FLOATS), np.float32)
    look[:, :4] = colours[geometry.atoms]
    look[:, 4] = flags[geometry.atoms]
    return look


def cylinder_look(geometry: Geometry, colours: np.ndarray, flags: np.ndarray) -> np.ndarray:
    ends = geometry.ends
    look = np.empty((len(ends), CYLINDER_LOOK_FLOATS), np.float32)
    look[:, 0:4] = colours[ends[:, 0]]
    look[:, 4:8] = colours[ends[:, 1]]
    look[:, 8] = flags[ends[:, 0]]
    look[:, 9] = flags[ends[:, 1]]
    return look


def flags(n_atoms: int, selected: np.ndarray | None = None, ghost: np.ndarray | None = None):
    out = np.zeros(n_atoms, np.float32)
    if selected is not None:
        out[selected] += SELECTED
    if ghost is not None:
        out[ghost] += GHOST
    return out
