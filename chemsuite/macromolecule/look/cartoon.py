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
"""The cartoon: helices as ribbons, strands as arrows, the rest as a tube.

Each unbroken stretch of a chain is a Catmull-Rom spline through its CA atoms (P for nucleic
acids). Along it a profile is swept: a superellipse whose half-width, half-thickness and
squareness change smoothly with the secondary structure (round tube, flat ellipse for a helix,
square slab for a strand, which widens and narrows to an arrow on its last residue). All
profiles have the same number of points, so rings join without seams; both ends of a stretch
are capped, so the mesh is closed. The profile's flat side follows the peptide planes (CA to
O, flipped where a strand's planes alternate), so ribbons do not twist.

All of it is numpy, per stretch. Every vertex's owner is its residue.
"""

from __future__ import annotations

import numpy as np

from ..structure.model import AMINO, HELIX, NUCLEIC, STRAND, Structure
from .mesh import VERTEX_FLOATS, Mesh, join

#: spline samples per residue, and points around the profile
SAMPLES = 6
POINTS = 12
#: (half-width, half-thickness, squareness) of each profile
TUBE = (0.25, 0.25, 2.0)
HELIX_PROFILE = (0.8, 0.16, 2.0)
STRAND_PROFILE = (0.85, 0.22, 6.0)
ARROW_WIDTH = 1.6
ARROW_TIP = 0.08
NUCLEIC_TUBE = (0.45, 0.45, 2.0)
#: CA-CA (or P-P) longer than this breaks the chain
_BREAK = {AMINO: 4.3, NUCLEIC: 8.5}
_TRACE = {AMINO: ("CA",), NUCLEIC: ("P", "C4'")}
_GUIDE = {AMINO: "O", NUCLEIC: "C1'"}


def trace(structure: Structure, active: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Residues the cartoon passes through and their trace atom (CA, or P / C4')."""
    if active is None:
        active = structure.active()
    residues: list[int] = []
    atoms: list[int] = []
    for kind, names in _TRACE.items():
        for name in names:
            found = np.flatnonzero(
                active & (structure.name == name) & (structure.atom_kind == kind)
            )
            have = set(residues)
            for atom in found:
                residue = int(structure.residue_of[atom])
                if residue not in have:
                    residues.append(residue)
                    atoms.append(int(atom))
                    have.add(residue)
    order = np.argsort(residues)
    return np.array(residues, np.int64)[order], np.array(atoms, np.int64)[order]


def build(
    structure: Structure,
    residues: np.ndarray | None = None,
    model: int = 0,
    active: np.ndarray | None = None,
) -> Mesh:
    """The cartoon of ``residues`` (all polymer residues by default)."""
    if active is None:
        active = structure.active()
    traced, atoms = trace(structure, active)
    if residues is not None:
        keep = np.isin(traced, residues)
        traced, atoms = traced[keep], atoms[keep]
    if len(traced) < 2:
        return Mesh.empty()
    xyz = structure.xyz(model).astype(np.float64)
    guides = _guides(structure, traced, atoms, active, xyz)
    meshes = []
    for segment in _segments(structure, traced, xyz[atoms]):
        if len(segment) < 2:
            continue
        kind = int(structure.res_kind[traced[segment[0]]])
        meshes.append(
            _sweep(xyz[atoms[segment]], guides[segment], traced[segment], structure.ss, kind)
        )
    return join(meshes)


def _guides(
    structure: Structure,
    residues: np.ndarray,
    atoms: np.ndarray,
    active: np.ndarray,
    xyz: np.ndarray,
) -> np.ndarray:
    """Per residue, a direction across the backbone (CA -> O), or zeros when unknown."""
    out = np.zeros((len(residues), 3))
    for kind, name in _GUIDE.items():
        found = np.flatnonzero(active & (structure.name == name) & (structure.atom_kind == kind))
        by_residue = dict(zip(structure.residue_of[found].tolist(), found.tolist(), strict=True))
        for i, residue in enumerate(residues.tolist()):
            atom = by_residue.get(residue)
            if atom is not None and int(structure.res_kind[residue]) == kind:
                out[i] = xyz[atom] - xyz[atoms[i]]
    return out


def _segments(structure: Structure, residues: np.ndarray, points: np.ndarray) -> list[np.ndarray]:
    """Index runs into ``residues`` that are one unbroken stretch of one chain."""
    if not len(residues):
        return []
    chain = structure.chain_of[residues]
    kind = structure.res_kind[residues]
    gap = np.linalg.norm(np.diff(points, axis=0), axis=1)
    limit = np.array([_BREAK.get(int(k), 4.3) for k in kind[1:]])
    breaks = (chain[1:] != chain[:-1]) | (kind[1:] != kind[:-1]) | (gap > limit)
    starts = np.concatenate([[0], np.flatnonzero(breaks) + 1, [len(residues)]])
    return [np.arange(a, b) for a, b in zip(starts[:-1], starts[1:], strict=True)]


def _sweep(
    points: np.ndarray, guides: np.ndarray, residues: np.ndarray, ss: np.ndarray, kind: int
) -> Mesh:
    m = len(points)
    centre, tangent, u = _spline(points)
    owner_index = np.minimum(np.floor(u + 0.5).astype(np.int64), m - 1)
    side = _frames(points, guides, tangent, u)
    binormal = np.cross(tangent, side)
    width, thickness, square = _profiles(ss[residues], kind, u)
    ring = _rings(centre, side, binormal, width, thickness, square)
    n = len(centre)
    triangles = _tube_triangles(n)
    # caps: a centre vertex at each end, fanned to its ring
    first, last = n * POINTS, n * POINTS + 1
    k = np.arange(POINTS)
    start_cap = np.stack([np.full(POINTS, first), (k + 1) % POINTS, k], axis=1)
    end_cap = np.stack([np.full(POINTS, last), (n - 1) * POINTS + k,
                        (n - 1) * POINTS + (k + 1) % POINTS], axis=1)  # fmt: skip
    caps = np.array([[*centre[0], *(-tangent[0])], [*centre[-1], *tangent[-1]]], dtype=np.float32)
    vertices = np.concatenate([ring, caps])
    owner = np.concatenate([np.repeat(residues[owner_index], POINTS), residues[[0, -1]]])
    return Mesh(
        vertices.astype(np.float32),
        np.concatenate([triangles, start_cap, end_cap]).astype(np.uint32),
        owner.astype(np.int64),
    )


def _spline(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Catmull-Rom samples through ``points``: positions, unit tangents and the residue
    parameter (residue i at u = i)."""
    m = len(points)
    padded = np.concatenate([[2 * points[0] - points[1]], points, [2 * points[-1] - points[-2]]])
    t = np.arange(SAMPLES) / SAMPLES
    i = np.repeat(np.arange(m - 1), SAMPLES)
    t = np.tile(t, m - 1)
    i = np.concatenate([i, [m - 2]])
    t = np.concatenate([t, [1.0]])
    p0, p1, p2, p3 = (padded[i + k] for k in range(4))
    t1 = t[:, None]
    t2, t3 = t1 * t1, t1 * t1 * t1
    position = 0.5 * (
        2 * p1
        + (p2 - p0) * t1
        + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2
        + (3 * p1 - p0 - 3 * p2 + p3) * t3
    )
    velocity = 0.5 * (
        (p2 - p0) + 2 * (2 * p0 - 5 * p1 + 4 * p2 - p3) * t1 + 3 * (3 * p1 - p0 - 3 * p2 + p3) * t2
    )
    tangent = velocity / np.maximum(np.linalg.norm(velocity, axis=1, keepdims=True), 1e-9)
    return position, tangent, i + t


def _frames(
    points: np.ndarray, guides: np.ndarray, tangent: np.ndarray, u: np.ndarray
) -> np.ndarray:
    """A unit vector across the backbone at every sample, perpendicular to the tangent."""
    guides = guides.copy()
    known = np.linalg.norm(guides, axis=1) > 1e-6
    if known.sum() < 2:
        return _transported(tangent)
    # carry missing guides over from neighbours, then keep successive guides on one side
    index = np.where(known, np.arange(len(guides)), -1)
    index = np.maximum.accumulate(index)
    first = int(np.argmax(known))
    index[index < 0] = first
    guides = guides[index]
    for j in range(1, len(guides)):
        if guides[j] @ guides[j - 1] < 0:
            guides[j] = -guides[j]
    lo = np.floor(u).astype(np.int64).clip(0, len(guides) - 1)
    hi = np.minimum(lo + 1, len(guides) - 1)
    f = (u - lo)[:, None]
    g = guides[lo] * (1 - f) + guides[hi] * f
    g -= (g * tangent).sum(axis=1, keepdims=True) * tangent
    length = np.linalg.norm(g, axis=1, keepdims=True)
    bad = length[:, 0] < 1e-6
    side = g / np.maximum(length, 1e-9)
    if bad.any():
        side[bad] = _transported(tangent)[bad]
    return side


def _transported(tangent: np.ndarray) -> np.ndarray:
    """Sides by parallel transport, for traces without guide atoms."""
    side = np.zeros_like(tangent)
    helper = np.array([0.0, 0.0, 1.0]) if abs(tangent[0, 2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    s = np.cross(tangent[0], helper)
    side[0] = s / np.linalg.norm(s)
    for j in range(1, len(tangent)):
        s = side[j - 1] - (side[j - 1] @ tangent[j]) * tangent[j]
        norm = np.linalg.norm(s)
        side[j] = s / norm if norm > 1e-9 else side[j - 1]
    return side


def _profiles(
    ss: np.ndarray, kind: int, u: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Half-width, half-thickness and squareness at every sample."""
    if kind == NUCLEIC:
        a, b, n = NUCLEIC_TUBE
        return np.full(len(u), a), np.full(len(u), b), np.full(len(u), n)
    table = np.array([TUBE, HELIX_PROFILE, STRAND_PROFILE])
    per_residue = table[np.where(np.isin(ss, (HELIX, STRAND)), ss, 0)]
    lo = np.floor(u).astype(np.int64).clip(0, len(ss) - 1)
    hi = np.minimum(lo + 1, len(ss) - 1)
    f = u - lo
    f = (f * f * (3 - 2 * f))[:, None]  # smooth step between residues
    values = per_residue[lo] * (1 - f) + per_residue[hi] * f
    width, thickness, square = values[:, 0], values[:, 1], values[:, 2]
    # arrows: over the last residue of each strand, wide to a point
    strand = ss == STRAND
    ends = np.flatnonzero(strand & ~np.append(strand[1:], False))
    for end in ends:
        if end == 0 or not strand[end - 1]:
            continue
        over = (u >= end - 1) & (u <= end)
        width[over] = ARROW_WIDTH + (ARROW_TIP - ARROW_WIDTH) * (u[over] - (end - 1))
        thickness[over] = STRAND_PROFILE[1]
        square[over] = STRAND_PROFILE[2]
        after = (u > end) & (u < end + 1)  # from the tip, no wider than the tube
        width[after] = np.minimum(width[after], TUBE[0])
    return width, thickness, square


def _rings(
    centre: np.ndarray,
    side: np.ndarray,
    binormal: np.ndarray,
    width: np.ndarray,
    thickness: np.ndarray,
    square: np.ndarray,
) -> np.ndarray:
    angle = 2 * np.pi * np.arange(POINTS) / POINTS
    c, s = np.cos(angle)[None, :], np.sin(angle)[None, :]
    e = (2.0 / square)[:, None]
    x = width[:, None] * np.sign(c) * np.abs(c) ** e
    y = thickness[:, None] * np.sign(s) * np.abs(s) ** e
    nx = np.sign(c) * np.abs(c) ** (2 - e) / width[:, None]
    ny = np.sign(s) * np.abs(s) ** (2 - e) / thickness[:, None]
    position = (
        centre[:, None, :] + x[..., None] * side[:, None, :] + y[..., None] * binormal[:, None, :]
    )
    normal = nx[..., None] * side[:, None, :] + ny[..., None] * binormal[:, None, :]
    normal /= np.maximum(np.linalg.norm(normal, axis=2, keepdims=True), 1e-9)
    out = np.empty((len(centre), POINTS, VERTEX_FLOATS))
    out[..., :3] = position
    out[..., 3:] = normal
    return out.reshape(-1, VERTEX_FLOATS)


def _tube_triangles(rings: int) -> np.ndarray:
    j = np.arange(rings - 1)[:, None]
    k = np.arange(POINTS)[None, :]
    a = j * POINTS + k
    b = j * POINTS + (k + 1) % POINTS
    c = (j + 1) * POINTS + k
    d = (j + 1) * POINTS + (k + 1) % POINTS
    first = np.stack([a, b, d], axis=-1).reshape(-1, 3)
    second = np.stack([a, d, c], axis=-1).reshape(-1, 3)
    return np.concatenate([first, second])
