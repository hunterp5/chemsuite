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
"""Cartoons and surfaces: closed, facing out, and every vertex owned."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.protein import cartoon, mesh, surface
from chemsuite.protein.model import Structure
from chemsuite.protein.select import select

from .conftest import load


def edge_counts(m: mesh.Mesh) -> np.ndarray:
    t = m.triangles.astype(np.int64)
    edges = np.sort(np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1)
    return np.unique(edges, axis=0, return_counts=True)[1]


def volume(m: mesh.Mesh) -> float:
    v = m.positions.astype(np.float64)
    t = m.triangles.astype(np.int64)
    return float(np.einsum("ij,ij->i", v[t[:, 0]], np.cross(v[t[:, 1]], v[t[:, 2]])).sum() / 6)


def normals_agree(m: mesh.Mesh) -> float:
    v = m.positions.astype(np.float64)
    t = m.triangles.astype(np.int64)
    face = np.cross(v[t[:, 1]] - v[t[:, 0]], v[t[:, 2]] - v[t[:, 0]])
    return float(((face * m.normals[t[:, 0]]).sum(axis=1) > 0).mean())


@pytest.mark.parametrize(
    "name", ["1crn.cif.gz", "1atp.cif.gz", "1ppb.cif.gz", "1d3z_2models.cif.gz"]
)
def test_cartoons_are_closed_face_out_and_belong_to_residues(name: str) -> None:
    s = load(name)
    m = cartoon.build(s)
    assert len(m) and m.triangles.max() < len(m)
    assert mesh.open_edges(m) == 0  # every edge shared by exactly two triangles
    assert volume(m) > 0  # wound counter-clockwise from outside
    assert normals_agree(m) > 0.95
    residues, _ = cartoon.trace(s)
    assert np.isin(m.owner, residues).all()
    assert set(np.unique(m.owner)) == set(residues.tolist())


def test_a_chain_break_splits_the_cartoon(kinase: Structure) -> None:
    whole = cartoon.build(kinase)
    # two chains (E and I), each closed with a cap at both ends: one centre vertex per cap
    caps = len(whole) % cartoon.POINTS
    assert caps == 2 * 2
    residues, atoms = cartoon.trace(kinase)
    gaps = np.linalg.norm(np.diff(kinase.coords[0][atoms], axis=0), axis=1)
    assert (gaps > 4.3).sum() == 1  # only the step from chain E to chain I


def test_a_cartoon_of_some_residues_only(crambin: Structure) -> None:
    first = np.arange(10)
    m = cartoon.build(crambin, residues=first)
    assert set(np.unique(m.owner)) == set(first.tolist())
    assert mesh.open_edges(m) == 0
    assert len(cartoon.build(crambin, residues=np.array([3]))) == 0  # one residue: no tube


def test_helices_are_wider_than_the_coil(crambin: Structure) -> None:
    m = cartoon.build(crambin)
    residues, atoms = cartoon.trace(crambin)
    ca = dict(zip(residues.tolist(), crambin.coords[0][atoms], strict=True))

    def reach(residue: int) -> float:
        rows = m.owner == residue
        return float(np.linalg.norm(m.positions[rows] - ca[residue], axis=1).max())

    helix = [r for r in residues if crambin.ss[r] == 1 and crambin.ss[r - 1] == 1]
    coil = [r for r in residues if crambin.ss[r] == 0 and crambin.ss[r - 1] == 0]
    assert np.median([reach(r) for r in helix]) > np.median([reach(r) for r in coil])


def test_surfaces_are_closed_face_out_and_belong_to_atoms(kinase: Structure) -> None:
    atoms = np.flatnonzero(select(kinase, "ligand"))
    m = surface.build(kinase, atoms)
    counts = edge_counts(m)
    assert (counts % 2 == 0).all()  # closed (an edge where the skin touches itself: four)
    assert (counts == 2).mean() > 0.99
    assert volume(m) > 0 and normals_agree(m) > 0.99
    assert np.isin(m.owner, atoms).all()
    # the skin lies a little beyond the padded van der Waals radii (neighbours add up)
    from scipy.spatial import cKDTree

    distance, _ = cKDTree(kinase.coords[0][atoms]).query(m.positions)
    assert 1.7 < np.median(distance) < 2.8


def test_a_lone_atom_has_a_sphere_at_its_padded_radius(kinase: Structure) -> None:
    from chemsuite.molview.elements import vdw_radii

    (ion,) = np.flatnonzero(select(kinase, "ion"))[:1]
    m = surface.build(kinase, np.array([ion]))
    distance = np.linalg.norm(m.positions - kinase.coords[0][ion], axis=1)
    expected = float(vdw_radii(kinase.element[[ion]])[0]) + surface.PADDING
    assert np.allclose(distance, expected, atol=0.12)


def test_mesh_look_gathers_owner_colours_and_flags(crambin: Structure) -> None:
    m = cartoon.build(crambin)
    colours = np.random.default_rng(1).random((crambin.n_residues, 4)).astype(np.float32)
    flags = np.zeros(crambin.n_residues, np.float32)
    flags[5] = 1.0
    look = mesh.look(m, colours, flags)
    assert np.allclose(look[:, :4], colours[m.owner])
    assert np.array_equal(look[:, 4] == 1.0, m.owner == 5)


def test_joined_meshes_keep_their_triangles(crambin: Structure) -> None:
    a = cartoon.build(crambin, residues=np.arange(10))
    b = cartoon.build(crambin, residues=np.arange(20, 30))
    both = mesh.join([a, b, mesh.Mesh.empty()])
    assert len(both) == len(a) + len(b)
    assert mesh.open_edges(both) == 0
    assert np.array_equal(both.triangles[len(a.triangles) :], b.triangles + len(a))
