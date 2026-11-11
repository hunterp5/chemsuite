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
"""Triangle meshes the protein view draws (cartoons, surfaces), and their look.

Every vertex belongs to an *owner*: a residue for a cartoon, an atom for a surface. Colours
and flags are given per owner and gathered per vertex (:func:`look`), so a selection or colour
change rebuilds no mesh.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

VERTEX_FLOATS = 6  # x y z nx ny nz
LOOK_FLOATS = 5  # r g b a flags


@dataclass(frozen=True)
class Mesh:
    vertices: np.ndarray  # float32 (v, 6): position, unit normal
    triangles: np.ndarray  # uint32 (t, 3), counter-clockwise seen from outside
    owner: np.ndarray  # int64 (v,): residue or atom index of each vertex

    @classmethod
    def empty(cls) -> Mesh:
        return cls(
            np.zeros((0, VERTEX_FLOATS), np.float32),
            np.zeros((0, 3), np.uint32),
            np.zeros(0, np.int64),
        )

    @property
    def positions(self) -> np.ndarray:
        return self.vertices[:, :3]

    @property
    def normals(self) -> np.ndarray:
        return self.vertices[:, 3:]

    def __len__(self) -> int:
        return len(self.vertices)


def join(meshes: list[Mesh]) -> Mesh:
    meshes = [m for m in meshes if len(m)]
    if not meshes:
        return Mesh.empty()
    offsets = np.cumsum([0] + [len(m) for m in meshes[:-1]])
    return Mesh(
        np.concatenate([m.vertices for m in meshes]),
        np.concatenate([m.triangles + np.uint32(o) for m, o in zip(meshes, offsets, strict=True)]),
        np.concatenate([m.owner for m in meshes]),
    )


def look(mesh: Mesh, colours: np.ndarray, flags: np.ndarray) -> np.ndarray:
    """Per vertex ``r g b a flags`` from per-owner colours (k, 4) and flags (k,)."""
    out = np.empty((len(mesh), LOOK_FLOATS), np.float32)
    out[:, :4] = colours[mesh.owner]
    out[:, 4] = flags[mesh.owner]
    return out


def open_edges(mesh: Mesh) -> int:
    """How many edges are used by other than exactly two triangles (0: watertight)."""
    t = mesh.triangles.astype(np.int64)
    edges = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return int((counts != 2).sum())
