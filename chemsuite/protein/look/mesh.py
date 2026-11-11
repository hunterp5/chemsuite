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


#: the most vertices a chunk can index with 16-bit indices
SHORT_INDICES = 65535


@dataclass(frozen=True)
class Chunk:
    """Part of a mesh small enough for 16-bit indices: which of the mesh's vertices it uses,
    and its triangles numbered within those."""

    vertex_map: np.ndarray  # int64 (k,), k <= the limit: the mesh's vertex for each
    triangles: np.ndarray  # uint16 (t, 3)


def chunks(mesh: Mesh, limit: int = SHORT_INDICES) -> list[Chunk]:
    """The mesh cut into chunks of at most ``limit`` vertices each, for GPUs without 32-bit
    indices. Together the chunks hold every triangle once, though not in the mesh's order
    (the view tests depth, so the order of triangles does not change the picture).

    Meshes are built piece by piece (a cartoon's segments, a surface's grid cells), so a
    triangle's vertices are close in number. Triangles are put in buckets by their lowest
    vertex, ``limit // 2`` vertex numbers to a bucket: a triangle spanning fewer than that
    reaches at most into the next bucket, so each bucket's vertices fit in one chunk. The rare
    triangle that spans more (a cap closing a long tube) goes with others like it, ``limit //
    3`` triangles to a chunk. Any mesh can so be drawn.
    """
    tri = np.asarray(mesh.triangles, np.int64)
    if not len(tri):
        return []
    width = max(limit // 2, 1)
    lo, hi = tri.min(axis=1), tri.max(axis=1)
    narrow = hi - lo < width
    out = []
    near = tri[narrow]
    bucket = lo[narrow] // width
    order = np.argsort(bucket, kind="stable")
    near, bucket = near[order], bucket[order]
    breaks = np.flatnonzero(np.diff(bucket)) + 1
    for part in np.split(near, breaks):
        if len(part):
            out.append(_chunk(part))
    wide = tri[~narrow]
    step = max(limit // 3, 1)
    for first in range(0, len(wide), step):
        out.append(_chunk(wide[first : first + step]))
    return out


def _chunk(triangles: np.ndarray) -> Chunk:
    vertices, numbered = np.unique(triangles, return_inverse=True)
    return Chunk(vertices, numbered.reshape(-1, 3).astype(np.uint16))


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


DASH = 0.25  # Å of tube, then a gap
GAP = 0.2
DASH_RADIUS = 0.05
_DASH_SIDES = 6


def dashes(starts: np.ndarray, ends: np.ndarray, radius: float = DASH_RADIUS) -> Mesh:
    """Dashed lines from ``starts`` to ``ends`` (k, 3) as thin open tubes; a vertex's owner is
    its line's index. The dashes are centred, so a line looks the same from either end."""
    starts = np.asarray(starts, np.float64).reshape(-1, 3)
    ends = np.asarray(ends, np.float64).reshape(-1, 3)
    if not len(starts):
        return Mesh.empty()
    axis = ends - starts
    length = np.linalg.norm(axis, axis=1)
    count = np.maximum(1, np.floor((length + GAP) / (DASH + GAP))).astype(np.int64)
    line = np.repeat(np.arange(len(starts)), count)
    k = np.arange(count.sum()) - np.repeat(np.cumsum(count) - count, count)
    used = count * DASH + (count - 1) * GAP
    first = (length - np.minimum(used, length)) / 2
    unit = axis / np.maximum(length, 1e-9)[:, None]
    begin = first[line] + k * (DASH + GAP)
    stop = np.minimum(begin + DASH, length[line])
    a = starts[line] + unit[line] * begin[:, None]
    b = starts[line] + unit[line] * stop[:, None]
    # a ring of sides around each end, on two axes perpendicular to the line
    d = unit[line]
    helper = np.where(np.abs(d[:, :1]) < 0.9, [[1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0]])
    u = np.cross(d, helper)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v = np.cross(d, u)
    angles = np.linspace(0.0, 2 * np.pi, _DASH_SIDES, endpoint=False)
    ring = np.cos(angles)[None, :, None] * u[:, None] + np.sin(angles)[None, :, None] * v[:, None]
    n = len(a)
    positions = np.concatenate([a[:, None] + radius * ring, b[:, None] + radius * ring], axis=1)
    normals = np.concatenate([ring, ring], axis=1)
    vertices = np.concatenate([positions, normals], axis=2).reshape(-1, VERTEX_FLOATS)
    s = _DASH_SIDES
    j = np.arange(s)
    nxt = (j + 1) % s
    quad = np.stack([np.stack([j, nxt, s + nxt], 1), np.stack([j, s + nxt, s + j], 1)], 1)
    triangles = (quad.reshape(-1, 3)[None] + (2 * s * np.arange(n))[:, None, None]).reshape(-1, 3)
    return Mesh(
        vertices.astype(np.float32),
        triangles.astype(np.uint32),
        np.repeat(line, 2 * s).astype(np.int64),
    )
