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
  between shown atoms (``x0 y0 z0 radius x1 y1 z1 0``). Rebuilt when atoms move or the styles
  change.
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

from .model import ADDITIVE, AMINO, GLYCAN, ION, LIGAND, NUCLEIC, Structure

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
    """A first look at a complex: the ligands as ball and stick, ions as balls, additives and
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
    styles[kind == LIGAND] = BALL_AND_STICK
    styles[kind == ION] = BALL_AND_STICK
    styles[np.isin(kind, (ADDITIVE, GLYCAN))] = STICK
    styles[nonpolar_hydrogens(structure)] = HIDDEN
    if active is None:
        active = structure.active()
    styles[~active] = HIDDEN
    return styles


def binding_site(structure: Structure, distance: float = SITE) -> np.ndarray:
    """Atoms of the residues with an atom within ``distance`` of a ligand."""
    from .select import within

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
    cylinders = np.zeros((len(shown), CYLINDER_FLOATS), np.float32)
    cylinders[:, :3] = xyz[i[shown]]
    cylinders[:, 3] = bond_radius[shown]
    cylinders[:, 4:7] = xyz[j[shown]]
    ends = np.stack([i[shown], j[shown]], axis=1).astype(np.int64)
    return Geometry(atoms, spheres, shown.astype(np.int64), ends, cylinders)


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
