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
"""Atoms and bonds as GPU-ready rows: one sphere per atom, one cylinder per bond (or two or
three for a double or triple bond, side by side in the plane of a neighbouring atom, as
3Dmol.js draws them).

Sphere rows are ``x y z radius r g b a`` and cylinder rows ``x0 y0 z0 radius x1 y1 z1 0``
followed by the colour of each half (``r g b a`` twice), all float32: the layouts the
renderer's shaders read.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from .elements import colours, vdw_radii

SPHERE_FLOATS = 8
CYLINDER_FLOATS = 16

BALL_AND_STICK = "ball_and_stick"
STICK = "stick"
SPACEFILL = "spacefill"
WIRE = "wire"
STYLES = {BALL_AND_STICK: "Ball and Stick", STICK: "Stick", SPACEFILL: "Spacefill", WIRE: "Wire"}
THIN = "thin"  # overlays: many conformers on top of each other stay readable

ALL_H = "all"
POLAR_H = "polar"
NO_H = "none"
HYDROGENS = {ALL_H: "All Hydrogens", POLAR_H: "Polar Hydrogens", NO_H: "No Hydrogens"}


@dataclass(frozen=True)
class Look:
    """How atoms and bonds are sized for a style."""

    atom_scale: float  # times the van der Waals radius; 0: use atom_radius
    atom_radius: float
    bond_radius: float  # 0: no bonds


LOOKS = {
    BALL_AND_STICK: Look(atom_scale=0.25, atom_radius=0.0, bond_radius=0.1),
    STICK: Look(atom_scale=0.0, atom_radius=0.15, bond_radius=0.15),
    SPACEFILL: Look(atom_scale=1.0, atom_radius=0.0, bond_radius=0.0),
    WIRE: Look(atom_scale=0.0, atom_radius=0.045, bond_radius=0.045),
    THIN: Look(atom_scale=0.0, atom_radius=0.09, bond_radius=0.09),
}

#: a thinner cylinder per extra bond, and how far apart they sit (in bond radii)
_DOUBLE = (0.6, (-0.9, 0.9))
_TRIPLE = (0.45, (-1.25, 0.0, 1.25))


@dataclass
class Scene:
    """Everything one view draws, plus what picking and labels need to know per atom."""

    spheres: np.ndarray  # (n, 8) float32
    cylinders: np.ndarray  # (m, 16) float32
    atom_molecule: np.ndarray  # (n,) which molecule of the scene an atom belongs to
    atom_index: np.ndarray  # (n,) the atom's index in its molecule
    atom_number: np.ndarray  # (n,) atomic number

    @property
    def positions(self) -> np.ndarray:
        return self.spheres[:, :3]

    @property
    def radii(self) -> np.ndarray:
        return self.spheres[:, 3]

    def __len__(self) -> int:
        return len(self.spheres)

    @classmethod
    def empty(cls) -> Scene:
        return combine([])


@dataclass
class _Part:
    spheres: np.ndarray
    cylinders: np.ndarray
    atom_index: np.ndarray
    atom_number: np.ndarray


def molecule_part(
    mol: Chem.Mol,
    conf_id: int = -1,
    *,
    style: str = BALL_AND_STICK,
    hydrogens: str = ALL_H,
    carbon: tuple[float, float, float] | None = None,
) -> _Part:
    """One conformer of *mol* drawn in *style*; *carbon* recolours carbon atoms."""
    look = LOOKS[style]
    if hydrogens != NO_H and any(a.GetTotalNumHs() for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)  # hydrogens placed on the conformer
    conf = mol.GetConformer(conf_id)
    xyz = np.array(conf.GetPositions(), np.float32).reshape(-1, 3)
    numbers = np.array([a.GetAtomicNum() for a in mol.GetAtoms()], np.int64)
    keep = _kept(mol, numbers, hydrogens)
    rgb = colours(numbers)
    if carbon is not None:
        rgb = rgb.copy()
        rgb[numbers == 6] = carbon
    radius = (
        vdw_radii(numbers) * look.atom_scale
        if look.atom_scale
        else np.full(len(numbers), look.atom_radius, np.float32)
    )
    spheres = np.zeros((len(numbers), SPHERE_FLOATS), np.float32)
    spheres[:, :3] = xyz
    spheres[:, 3] = radius
    spheres[:, 4:7] = rgb
    spheres[:, 7] = 1.0
    cylinders = (
        _bonds(mol, xyz, rgb, keep, look.bond_radius) if look.bond_radius else _no_cylinders()
    )
    index = np.flatnonzero(keep)
    return _Part(spheres[index], cylinders, index, numbers[index])


def combine(parts: Sequence[_Part]) -> Scene:
    """One scene from several molecules (or conformers), in order."""
    if not parts:
        return Scene(
            np.zeros((0, SPHERE_FLOATS), np.float32),
            _no_cylinders(),
            np.zeros(0, np.int64),
            np.zeros(0, np.int64),
            np.zeros(0, np.int64),
        )
    return Scene(
        np.concatenate([p.spheres for p in parts]),
        np.concatenate([p.cylinders for p in parts]),
        np.concatenate([np.full(len(p.atom_index), i, np.int64) for i, p in enumerate(parts)]),
        np.concatenate([p.atom_index for p in parts]),
        np.concatenate([p.atom_number for p in parts]),
    )


def _kept(mol: Chem.Mol, numbers: np.ndarray, hydrogens: str) -> np.ndarray:
    keep = numbers != 1
    if hydrogens == ALL_H:
        keep[:] = True
    elif hydrogens == POLAR_H:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1 and any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors()):
                keep[atom.GetIdx()] = True
    return keep


def _no_cylinders() -> np.ndarray:
    return np.zeros((0, CYLINDER_FLOATS), np.float32)


def _orders(mol: Chem.Mol) -> list[int]:
    """Bond orders with aromatic rings shown as alternating single and double bonds."""
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Exception:  # noqa: BLE001 - an unkekulizable ring is drawn with single bonds
        kekule = mol
    out = []
    for bond in kekule.GetBonds():
        order = bond.GetBondTypeAsDouble()
        out.append(3 if order == 3 else 2 if order == 2 else 1)
    return out


def _bonds(mol, xyz, rgb, keep, radius) -> np.ndarray:
    rows = []
    for bond, order in zip(mol.GetBonds(), _orders(mol), strict=True):
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if not (keep[i] and keep[j]):
            continue
        a, b = xyz[i], xyz[j]
        if order == 1:
            rows.append(_cylinder(a, b, radius, rgb[i], rgb[j]))
            continue
        scale, offsets = _DOUBLE if order == 2 else _TRIPLE
        side = _side(mol, xyz, keep, i, j) * radius
        for k in offsets:
            shift = side * k
            rows.append(_cylinder(a + shift, b + shift, radius * scale, rgb[i], rgb[j]))
    return np.array(rows, np.float32).reshape(-1, CYLINDER_FLOATS)


def _cylinder(a, b, radius, ca, cb) -> list[float]:
    return [*a, radius, *b, 0.0, *ca, 1.0, *cb, 1.0]


def _side(mol, xyz, keep, i, j) -> np.ndarray:
    """A unit vector across bond i-j, in the plane of a neighbouring atom when there is one."""
    axis = xyz[j] - xyz[i]
    axis = axis / max(float(np.linalg.norm(axis)), 1e-6)
    for here, other in ((i, j), (j, i)):
        for n in mol.GetAtomWithIdx(here).GetNeighbors():
            k = n.GetIdx()
            if k != other and keep[k]:
                v = xyz[k] - xyz[here]
                v = v - axis * float(v @ axis)
                norm = float(np.linalg.norm(v))
                if norm > 1e-3:
                    return v / norm
    guess = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    v = np.cross(axis, guess)
    return v / np.linalg.norm(v)
