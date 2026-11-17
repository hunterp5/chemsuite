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
"""Molecular surfaces, contoured with surface nets: four kinds.

- **SES** (solvent excluded, Connolly's molecular surface; the default): where a probe sphere
  rolling over the atoms cannot reach. On a grid, the solvent-accessible inside is found, and
  the SES is the set of points at least the probe radius from the solvent outside it (the
  Euclidean distance transform; EDTSurf's method, Xu & Zhang, PLoS ONE 2009).
- **SAS** (solvent accessible): the atoms' radii grown by the probe's, where its center can go.
- **vdW**: the union of the atoms' van der Waals spheres.
- **Blob** (Gaussian): each atom adds a Gaussian whose level set at :data:`LEVEL` lies at its
  van der Waals radius plus :data:`PADDING`, so neighbours blend into a smooth skin, as in
  VMD's QuickSurf. Atoms are splatted onto a grid and blurred with one Gaussian per radius.

The vdW and SAS kinds are signed distance fields, ``min(|x - c| - r)`` over the atoms: per
radius, the distance transform's feature map finds each grid point's nearest atom, and the
distance to that atom's true center is exact, so there is no loop over atoms.

Surface nets put one vertex in each grid cell the surface crosses, at the mean of the
crossings on its edges, and one quad across every grid edge that crosses; the grid is padded,
so the surface is closed. *Smoothing* (Taubin's lambda-mu filter, which does not shrink) takes
the grid's terraces off the distance-field kinds. Each vertex's owner is its nearest atom, so
the surface takes the atoms' colors and selection.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates
from scipy.spatial import cKDTree

from chemsuite.molview.elements import vdw_radii

from ..structure.model import Structure
from .mesh import VERTEX_FLOATS, Mesh

SES, SAS, VDW, BLOB = "ses", "sas", "vdw", "blob"
KINDS = {
    SES: "Molecular (SES)",
    SAS: "Solvent Accessible (SAS)",
    VDW: "Van der Waals",
    BLOB: "Gaussian",
}
#: grid spacing (Å) of each quality for the distance-field kinds; the blob's is 1.2 times it
QUALITIES = {"coarse": 0.8, "normal": 0.5, "fine": 0.35}
#: Taubin passes of each smoothing setting
SMOOTHING = {"off": 0, "low": 4, "high": 16}
PROBE = 1.4  # water
PROBES = (1.0, 1.2, 1.4, 1.6, 2.0)

LEVEL = 0.3
PADDING = 0.3
#: grid spacing (Å) of the blob, widened for big boxes so the grid stays under MAX_VOXELS
SPACING = 0.6
MAX_VOXELS = 12_000_000
#: the distance-field kinds hold a few grids at once (the transform's feature map among
#: them), so they stay under fewer voxels
MAX_FIELD_VOXELS = 8_000_000
#: a Gaussian's sigma for an atom of radius r: its LEVEL contour at r
_SIGMA_PER_RADIUS = 1.0 / np.sqrt(2.0 * np.log(1.0 / LEVEL))
#: Taubin's weights: a smoothing step, then an inflating one slightly larger
_LAMBDA, _MU = 0.5, -0.53


@dataclass(frozen=True)
class SurfaceOptions:
    """How a surface is made: its *kind*, the *probe* radius (Å; SES and SAS), the grid's
    *quality* and the *smoothing* (keys of :data:`QUALITIES` and :data:`SMOOTHING`)."""

    kind: str = SES
    probe: float = PROBE
    quality: str = "normal"
    smoothing: str = "low"

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown surface kind: {self.kind!r}")
        if not 0.5 <= self.probe <= 3.0:
            raise ValueError("the probe radius must be between 0.5 and 3 Å")
        if self.quality not in QUALITIES:
            raise ValueError(f"unknown surface quality: {self.quality!r}")
        if self.smoothing not in SMOOTHING:
            raise ValueError(f"unknown smoothing: {self.smoothing!r}")

    @property
    def spacing(self) -> float:
        step = QUALITIES[self.quality]
        return step * SPACING / QUALITIES["normal"] if self.kind == BLOB else step

    def state(self) -> dict:
        return {
            "kind": self.kind,
            "probe": self.probe,
            "quality": self.quality,
            "smoothing": self.smoothing,
        }

    @classmethod
    def from_state(cls, state: object) -> SurfaceOptions:
        """Options read back from :meth:`state`; anything missing or wrong is a Blob, as every
        surface was before there were kinds."""
        if not isinstance(state, dict):
            return BLOB_OPTIONS
        try:
            return cls(
                kind=str(state.get("kind", BLOB)),
                probe=float(state.get("probe", PROBE)),
                quality=str(state.get("quality", "normal")),
                smoothing=str(state.get("smoothing", "low")),
            )
        except (TypeError, ValueError):
            return BLOB_OPTIONS


DEFAULT_OPTIONS = SurfaceOptions()
BLOB_OPTIONS = SurfaceOptions(kind=BLOB)


def build(
    structure: Structure,
    atoms: np.ndarray,
    model: int = 0,
    options: SurfaceOptions = BLOB_OPTIONS,
) -> Mesh:
    """The surface of ``atoms`` (indices), owned by those atoms."""
    atoms = np.asarray(atoms, np.int64)
    if not len(atoms):
        return Mesh.empty()
    centres = structure.xyz(model)[atoms].astype(np.float64)
    vdw = vdw_radii(structure.element[atoms]).astype(np.float64)
    if options.kind == BLOB:
        density, origin, step = _density(centres, vdw + PADDING, options.spacing)
        vertices, quads = _surface_nets(density, LEVEL)
        if not len(vertices):
            return Mesh.empty()
        gradient = np.stack(
            [map_coordinates(g, vertices.T, order=1, mode="nearest") for g in np.gradient(density)],
            axis=1,
        )
        normals = -gradient / np.maximum(np.linalg.norm(gradient, axis=1, keepdims=True), 1e-9)
        positions = origin + vertices * step
        triangles = _triangles(quads)
    else:
        field, level, origin, step = _field(centres, vdw, options)
        vertices, quads = _surface_nets(field, level)
        if not len(vertices):
            return Mesh.empty()
        triangles = _triangles(quads)
        positions = smooth(origin + vertices * step, triangles, SMOOTHING[options.smoothing])
        normals = vertex_normals(positions, triangles)
    _, nearest = cKDTree(centres).query(positions)
    out = np.empty((len(positions), VERTEX_FLOATS), np.float32)
    out[:, :3] = positions
    out[:, 3:] = normals
    return Mesh(out, triangles, atoms[nearest])


def _triangles(quads: np.ndarray) -> np.ndarray:
    return np.concatenate([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]]).astype(np.uint32)


def _grid(
    centres: np.ndarray, reach: float, spacing: float, limit: int
) -> tuple[np.ndarray, float, np.ndarray]:
    """``(origin, step, shape)`` of a grid over *centres* with *reach* to spare, at *spacing*
    unless that would pass *limit* voxels."""
    low = centres.min(axis=0) - reach
    size = centres.max(axis=0) + reach - low
    step = max(spacing, float(np.cbrt(np.prod(size) / limit)))
    shape = np.ceil(size / step).astype(int) + 1
    return low, step, shape


def _field(
    centres: np.ndarray, vdw: np.ndarray, options: SurfaceOptions
) -> tuple[np.ndarray, float, np.ndarray, float]:
    """``(field, level, origin, step)``: a grid whose values are above *level* inside the
    surface (the negated signed distance for vdW and SAS; for SES, the distance to the
    solvent)."""
    grown = vdw + (0.0 if options.kind == VDW else options.probe)
    reach = float(grown.max()) + 3 * options.spacing
    origin, step, shape = _grid(centres, reach, options.spacing, MAX_FIELD_VOXELS)
    signed = _signed_distance(centres, grown, origin, step, shape)
    if options.kind != SES:
        return -signed, 0.0, origin, step
    inside = signed < 0
    # each inside point's distance to the solvent: to its nearest outside grid point, less
    # how far that point lies beyond the SAS (known exactly), so the SES does not sit half a
    # grid step out. The SES is where this passes the probe radius; outside, the (negative)
    # distance to the SAS keeps the field going.
    distance, index = distance_transform_edt(inside, sampling=step, return_indices=True)
    beyond = signed[index[0], index[1], index[2]]
    del index
    depth = (distance - np.where(inside, beyond, 0.0)).astype(np.float32)
    return np.where(inside, depth, -signed), options.probe, origin, step


def _signed_distance(
    centres: np.ndarray, radii: np.ndarray, origin: np.ndarray, step: float, shape: np.ndarray
) -> np.ndarray:
    """``min(|x - c| - r)`` over the atoms at every grid point (negative inside a sphere).

    Each radius's atoms are worked on in the box they span plus a margin: beyond it they are
    more than two grid steps away, so they cannot be the nearest surface there, and a rare
    radius (sulfur, a metal) costs only its own small box.
    """
    size = np.array([int(n) for n in shape])
    out = np.full(tuple(size), np.inf, np.float32)
    voxel = np.clip(np.round((centres - origin) / step).astype(np.int64), 0, size - 1)
    rounded = np.round(radii / 0.05) * 0.05
    for radius in np.unique(rounded):
        members = np.flatnonzero(rounded == radius)
        margin = int(np.ceil(radius / step)) + 3
        low = np.maximum(voxel[members].min(axis=0) - margin, 0)
        high = np.minimum(voxel[members].max(axis=0) + margin + 1, size)
        box = tuple(int(h - lo) for lo, h in zip(low, high, strict=True))
        local = voxel[members] - low
        free = np.ones(box, bool)
        owner = np.zeros(box, np.int32)
        free[tuple(local.T)] = False
        owner[tuple(local.T)] = members
        index = distance_transform_edt(free, return_distances=False, return_indices=True)
        nearest = owner[index[0], index[1], index[2]]
        del index
        squared = np.zeros(box, np.float32)
        for k in range(3):
            shape_k = [1, 1, 1]
            shape_k[k] = -1
            grid = (origin[k] + step * np.arange(low[k], high[k])).astype(np.float32)
            squared += (grid.reshape(shape_k) - centres[nearest, k].astype(np.float32)) ** 2
        region = tuple(slice(int(lo), int(h)) for lo, h in zip(low, high, strict=True))
        np.minimum(out[region], np.sqrt(squared) - np.float32(radius), out=out[region])
    return out


def smooth(positions: np.ndarray, triangles: np.ndarray, passes: int) -> np.ndarray:
    """*positions* moved by Taubin's lambda-mu smoothing over the mesh's edges, *passes* times
    (a shrinking step and an inflating one each, so the surface keeps its size)."""
    if not passes or not len(triangles):
        return positions
    from scipy.sparse import coo_matrix

    t = triangles.astype(np.int64)
    n = len(positions)
    a = np.concatenate([t[:, 0], t[:, 1], t[:, 2]])
    b = np.concatenate([t[:, 1], t[:, 2], t[:, 0]])
    # each edge is in two triangles, once each way round: a symmetric 0/1 adjacency
    adjacency = coo_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    adjacency = adjacency.tocsr()
    adjacency.data[:] = 1.0  # duplicates summed on conversion: back to one per neighbor
    degree = np.maximum(np.asarray(adjacency.sum(axis=1)).ravel(), 1.0)[:, None]
    out = positions.astype(np.float64).copy()
    for _ in range(passes):
        for weight in (_LAMBDA, _MU):
            out += weight * (adjacency @ out / degree - out)
    return out


def vertex_normals(positions: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Unit normals per vertex: the area-weighted sum of its triangles' outward normals."""
    t = triangles.astype(np.int64)
    p = positions.astype(np.float64)
    faces = np.cross(p[t[:, 1]] - p[t[:, 0]], p[t[:, 2]] - p[t[:, 0]])
    normals = np.zeros_like(p)
    for k in range(3):
        np.add.at(normals, t[:, k], faces)
    return normals / np.maximum(np.linalg.norm(normals, axis=1, keepdims=True), 1e-12)


def _density(
    centres: np.ndarray, radii: np.ndarray, spacing: float
) -> tuple[np.ndarray, np.ndarray, float]:
    reach = radii.max() * 2.0 + 2 * spacing
    low, step, shape = _grid(centres, reach, spacing, MAX_VOXELS)
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
