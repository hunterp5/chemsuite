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
"""A molecular surface: a Gaussian density of the atoms, contoured with surface nets.

Each atom adds a Gaussian whose level set at :data:`LEVEL` lies at its van der Waals radius
plus :data:`PADDING` (so neighbours blend into a smooth skin, as in VMD's QuickSurf). Atoms
are splatted onto a grid and blurred with one Gaussian per radius (scipy), which is fast for
large proteins. Surface nets put one vertex in each grid cell the surface crosses, at the mean
of the crossings on its edges, and one quad across every grid edge that crosses; the grid is
padded, so the surface is closed. Normals are the density's gradient; each vertex's owner is
its nearest atom, so the surface takes the atoms' colours and selection.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.spatial import cKDTree

from chemsuite.molview.elements import vdw_radii

from .mesh import VERTEX_FLOATS, Mesh
from .model import Structure

LEVEL = 0.3
PADDING = 0.3
#: grid spacing (Å), widened for big boxes so the grid stays under MAX_VOXELS
SPACING = 0.6
MAX_VOXELS = 12_000_000
#: a Gaussian's sigma for an atom of radius r: its LEVEL contour at r
_SIGMA_PER_RADIUS = 1.0 / np.sqrt(2.0 * np.log(1.0 / LEVEL))


def build(
    structure: Structure, atoms: np.ndarray, model: int = 0, spacing: float = SPACING
) -> Mesh:
    """The surface of ``atoms`` (indices), owned by those atoms."""
    atoms = np.asarray(atoms, np.int64)
    if not len(atoms):
        return Mesh.empty()
    centres = structure.xyz(model)[atoms].astype(np.float64)
    radii = vdw_radii(structure.element[atoms]).astype(np.float64) + PADDING
    density, origin, step = _density(centres, radii, spacing)
    vertices, quads = _surface_nets(density, LEVEL)
    if not len(vertices):
        return Mesh.empty()
    gradient = np.stack(
        [map_coordinates(g, vertices.T, order=1, mode="nearest") for g in np.gradient(density)],
        axis=1,
    )
    normals = -gradient / np.maximum(np.linalg.norm(gradient, axis=1, keepdims=True), 1e-9)
    positions = origin + vertices * step
    _, nearest = cKDTree(centres).query(positions)
    out = np.empty((len(positions), VERTEX_FLOATS), np.float32)
    out[:, :3] = positions
    out[:, 3:] = normals
    triangles = np.concatenate([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    return Mesh(out, triangles.astype(np.uint32), atoms[nearest])


def _density(
    centres: np.ndarray, radii: np.ndarray, spacing: float
) -> tuple[np.ndarray, np.ndarray, float]:
    reach = radii.max() * 2.0 + 2 * spacing
    low = centres.min(axis=0) - reach
    high = centres.max(axis=0) + reach
    size = high - low
    step = max(spacing, float(np.cbrt(np.prod(size) / MAX_VOXELS)))
    shape = np.ceil(size / step).astype(int) + 1
    density = np.zeros(shape)
    position = (centres - low) / step
    # one blur per radius, to 0.2 Å (carbon, nitrogen and oxygen: two or three blurs)
    rounded = np.round(radii / 0.2) * 0.2
    for radius in np.unique(rounded):
        group = rounded == radius
        splat = np.zeros(shape)
        _trilinear(splat, position[group])
        sigma = radius * _SIGMA_PER_RADIUS / step
        # a blurred unit splat peaks at 1 / (2 pi sigma^2)^(3/2): scale it back to 1
        density += gaussian_filter(splat, sigma, mode="constant", truncate=3.0) * (
            (2 * np.pi * sigma * sigma) ** 1.5
        )
    return density, low, step


def _trilinear(grid: np.ndarray, points: np.ndarray) -> None:
    base = np.floor(points).astype(int)
    frac = points - base
    for corner in range(8):
        offset = np.array([(corner >> 2) & 1, (corner >> 1) & 1, corner & 1])
        weight = np.prod(np.where(offset, frac, 1 - frac), axis=1)
        index = base + offset
        np.add.at(grid, (index[:, 0], index[:, 1], index[:, 2]), weight)


def _surface_nets(density: np.ndarray, level: float) -> tuple[np.ndarray, np.ndarray]:
    """Vertices (in grid units) and quads (counter-clockwise from outside)."""
    inside = density > level
    shape = np.array(density.shape)
    cells = shape - 1
    cell_index = np.full(cells, -1, np.int64)
    total = np.zeros((*cells, 3))
    count = np.zeros(cells)
    crossings = []
    for axis in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[axis], b[axis] = slice(None, -1), slice(1, None)
        changed = inside[tuple(a)] != inside[tuple(b)]
        points = np.argwhere(changed)
        d0, d1 = density[tuple(a)][changed], density[tuple(b)][changed]
        t = (level - d0) / (d1 - d0)
        where = points.astype(float)
        where[:, axis] += t
        # the cells around this edge: the other two axes at offset 0 or -1
        others = [k for k in range(3) if k != axis]
        for du in (0, -1):
            for dv in (0, -1):
                cell = points.copy()
                cell[:, others[0]] += du
                cell[:, others[1]] += dv
                ok = ((cell >= 0) & (cell < cells)).all(axis=1)
                np.add.at(total, tuple(cell[ok].T), where[ok])
                np.add.at(count, tuple(cell[ok].T), 1)
        crossings.append((axis, points, inside[tuple(a)][changed]))
    active = count > 0
    cell_index[active] = np.arange(int(active.sum()))
    vertices = total[active] / count[active][:, None]
    quads = []
    for axis, points, first_inside in crossings:
        u, v = [k for k in range(3) if k != axis]
        corners = []
        for du, dv in ((0, 0), (-1, 0), (-1, -1), (0, -1)):
            cell = points.copy()
            cell[:, u] += du
            cell[:, v] += dv
            corners.append(cell)
        ok = np.all([((c >= 0) & (c < cells)).all(axis=1) for c in corners], axis=0)
        ids = np.stack([cell_index[tuple(c[ok].T)] for c in corners], axis=1)
        # wind so the normal points from inside to outside
        flip = ~(first_inside[ok] ^ (axis == 1))
        ids[flip] = ids[flip][:, ::-1]
        quads.append(ids)
    return vertices, np.concatenate(quads) if quads else np.zeros((0, 4), np.int64)
