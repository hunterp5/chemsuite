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
"""The search box: where docking is allowed to put a ligand.

A :class:`Box` is a centre and a size in ångströms, as Vina takes them. It comes from the
atoms of a bound ligand, from a residue, or from whatever is selected, each padded so a ligand
of about that size has room to move and to grow a little
(:data:`PADDING`). :meth:`Box.holds` says whether a pose stayed inside it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..structure.model import Structure

#: Å added on each side of the atoms the box is built from
PADDING = 4.0
#: Å: no side of a box is smaller than this (a tiny ligand still needs room to turn)
SMALLEST = 10.0
#: Å: Vina's own limit on a side, and far more than a pocket needs
LARGEST = 100.0


@dataclass(frozen=True)
class Box:
    centre: tuple[float, float, float]
    size: tuple[float, float, float]

    @property
    def volume(self) -> float:
        return float(np.prod(self.size))

    def holds(self, xyz: np.ndarray) -> bool:
        """Is every point of *xyz* (n, 3) inside the box?"""
        points = np.asarray(xyz, dtype=float).reshape(-1, 3)
        if not len(points):
            return True
        half = np.asarray(self.size, dtype=float) / 2
        return bool((np.abs(points - np.asarray(self.centre, dtype=float)) <= half).all())

    def arguments(self) -> list[str]:
        """The box as Vina's command line takes it."""
        keys = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        values = (*self.centre, *self.size)
        return [f"--{k}={v:.3f}" for k, v in zip(keys, values, strict=True)]


def from_points(xyz: np.ndarray, padding: float = PADDING) -> Box:
    """A box around the points *xyz* (n, 3), *padding* ångströms clear on every side."""
    points = np.asarray(xyz, dtype=float).reshape(-1, 3)
    if not len(points):
        raise ValueError("a box needs at least one atom")
    low, high = points.min(axis=0), points.max(axis=0)
    centre = (low + high) / 2
    size = np.clip(high - low + 2 * padding, SMALLEST, LARGEST)
    cx, cy, cz = (float(v) for v in centre)
    sx, sy, sz = (float(v) for v in size)
    return Box((cx, cy, cz), (sx, sy, sz))


def from_atoms(s: Structure, atoms: np.ndarray, padding: float = PADDING, model: int = 0) -> Box:
    """A box around the structure's atoms *atoms* (a mask or indices)."""
    chosen = np.asarray(atoms)
    xyz = s.xyz(model)
    picked = xyz[chosen] if chosen.dtype == bool else xyz[chosen.astype(np.int64)]
    if not len(picked):
        raise ValueError("a box needs at least one atom")
    return from_points(picked, padding)


def from_residue(s: Structure, residue: int, padding: float = PADDING, model: int = 0) -> Box:
    """A box around one residue: the bound ligand whose pocket is being docked into."""
    span = s.atoms_of(residue)
    return from_points(s.xyz(model)[span], padding)
