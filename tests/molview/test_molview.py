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
"""3D scenes without a window: geometry, camera, picking, measuring and alignment."""

from __future__ import annotations

import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from chemsuite.molview.align import aligned_conformers
from chemsuite.molview.camera import Camera
from chemsuite.molview.geometry import (
    ALL_H,
    BALL_AND_STICK,
    NO_H,
    POLAR_H,
    SPACEFILL,
    STICK,
    combine,
    molecule_part,
)
from chemsuite.molview.measure import angle, describe, dihedral, distance
from chemsuite.molview.pick import pick_sphere


def embedded(smiles: str, count: int = 1) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=7)
    return mol


@pytest.mark.parametrize(
    "hydrogens, atoms", [(ALL_H, 9), (POLAR_H, 4), (NO_H, 3)]
)  # ethanol: C C O + 5 H on carbon + 1 on oxygen
def test_hydrogens_all_polar_or_none(hydrogens, atoms) -> None:
    part = molecule_part(embedded("CCO"), hydrogens=hydrogens)
    assert len(part.spheres) == atoms
    assert len(part.cylinders) == atoms - 1  # a tree: one bond fewer than atoms


def test_hydrogens_are_placed_when_the_input_has_none() -> None:
    mol = Chem.RemoveHs(embedded("CCO"))
    assert len(molecule_part(mol, hydrogens=ALL_H).spheres) == 9
    assert len(molecule_part(mol, hydrogens=NO_H).spheres) == 3


def test_double_and_triple_bonds_are_parallel_cylinders() -> None:
    part = molecule_part(embedded("C=CC#N"), hydrogens=NO_H, style=STICK)
    cyl = part.cylinders
    assert len(cyl) == 2 + 1 + 3
    axes = cyl[:, 4:7] - cyl[:, 0:3]
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    double, triple = axes[:2], axes[3:]
    assert abs(double[0] @ double[1]) == pytest.approx(1.0, abs=1e-5)
    assert np.allclose(np.abs(triple @ triple[0]), 1.0, atol=1e-5)
    assert np.all(cyl[:2, 3] < cyl[2, 3])  # thinner than a single bond


def test_aromatic_rings_alternate() -> None:
    part = molecule_part(embedded("c1ccccc1"), hydrogens=NO_H)
    assert len(part.cylinders) == 3 + 3 * 2


def test_bonds_take_each_atoms_colour_and_carbon_can_be_recoloured() -> None:
    part = molecule_part(embedded("CO"), hydrogens=NO_H, carbon=(0.0, 0.5, 1.0))
    (bond,) = part.cylinders
    assert tuple(bond[8:11]) == pytest.approx((0.0, 0.5, 1.0))
    assert tuple(bond[12:15]) == pytest.approx((1.0, 13 / 255, 13 / 255))


def test_spacefill_uses_van_der_waals_radii_and_no_bonds() -> None:
    part = molecule_part(embedded("O"), style=SPACEFILL)
    assert len(part.cylinders) == 0
    assert part.spheres[part.atom_number == 8][0, 3] == pytest.approx(1.55, abs=0.01)


def test_combine_numbers_molecules_in_order() -> None:
    mol = embedded("CCO", 2)
    scene = combine([molecule_part(mol, 0), molecule_part(mol, 1, style=BALL_AND_STICK)])
    assert len(scene) == 18
    assert scene.atom_molecule.tolist() == [0] * 9 + [1] * 9
    assert len(combine([])) == 0


def test_camera_fit_sees_every_atom() -> None:
    mol = embedded("c1ccc2ccccc2c1CCCCCCN")
    scene = combine([molecule_part(mol)])
    camera = Camera()
    camera.turn(40, -25)
    camera.fit(scene.positions, scene.radii, aspect=1.5)
    xy, depth = camera.project(scene.positions, 600, 400)
    assert np.all(depth > 0)
    assert np.all((xy >= 0) & (xy <= [600, 400]))
    near, far = camera.clip_range()
    assert np.all((depth > near) & (depth < far))


def drawn_box(camera: Camera, scene, width: int, height: int):
    """The atoms' outline on screen: (offset of its centre from the view's centre, the share
    of the view it fills in its tighter direction)."""
    xy, depth = camera.project(scene.positions, width, height)
    size = scene.radii * np.array([camera.pixels_per_unit(d, height) for d in depth])
    lo, hi = (xy - size[:, None]).min(axis=0), (xy + size[:, None]).max(axis=0)
    centre = (lo + hi) / 2 - [width / 2, height / 2]
    return centre, max((hi[0] - lo[0]) / width, (hi[1] - lo[1]) / height)


@pytest.mark.parametrize("size", [(520, 480), (900, 300), (300, 800)])
@pytest.mark.parametrize("smiles", ["CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CCCCCCCCCCCCCCCCCC(=O)O"])
def test_camera_fit_centres_and_fills_the_view(smiles, size) -> None:
    from chemsuite.molview.camera import MARGIN

    scene = combine([molecule_part(embedded(smiles), style=BALL_AND_STICK)])
    rng = np.random.default_rng(len(smiles) + size[0])
    for _ in range(3):  # from any direction
        q, _r = np.linalg.qr(rng.normal(size=(3, 3)))
        camera = Camera(rotation=q * np.sign(np.linalg.det(q)))
        camera.fit(scene.positions, scene.radii, aspect=size[0] / size[1])
        centre, fill = drawn_box(camera, scene, *size)
        assert np.all(np.abs(centre) < 3), centre  # pixels
        assert fill == pytest.approx(1 / MARGIN, abs=0.01)


def test_orient_lays_the_long_axis_along_the_longer_side() -> None:
    scene = combine([molecule_part(embedded("CCCCCCCCCCCCCCCC"))])
    for aspect, along in ((2.0, 0), (0.5, 1)):
        camera = Camera()
        camera.turn(70, 40)
        camera.orient(scene.positions, aspect)
        assert np.allclose(camera.rotation @ camera.rotation.T, np.eye(3), atol=1e-9)
        assert np.linalg.det(camera.rotation) == pytest.approx(1.0)
        view = (scene.positions - scene.positions.mean(axis=0)) @ camera.rotation.T
        spread = view.std(axis=0)
        assert spread[along] == spread.max() and spread[2] == spread.min()


def test_camera_ray_and_project_agree() -> None:
    camera = Camera()
    camera.turn(120, 30)
    camera.fit(np.array([[0.0, 0, 0], [3, 1, -2]]))
    origin, direction = camera.ray(150, 260, 400, 300)
    point = origin + direction * camera.distance
    xy, _ = camera.project(point[None], 400, 300)
    assert xy[0] == pytest.approx([150, 260], abs=1e-6)


def test_an_orthographic_camera_keeps_the_framing_and_its_rays_agree() -> None:
    camera = Camera()
    camera.turn(40, -25)
    camera.fit(np.array([[0.0, 0, 0], [6, 2, -3], [-2, 5, 1]]))
    centre_xy, _ = camera.project(camera.centre[None], 400, 300)
    near = camera.centre + camera.rotation.T @ np.array([2.0, 1.0, 0.0])  # on the centre plane
    before, _ = camera.project(near[None], 400, 300)
    camera.orthographic = True
    after, _ = camera.project(near[None], 400, 300)
    assert after[0] == pytest.approx(before[0])  # the centre's plane keeps its scale
    deeper = near + camera.rotation.T @ np.array([0.0, 0.0, -4.0])  # further from the eye
    assert camera.project(deeper[None], 400, 300)[0][0] == pytest.approx(after[0])  # parallel
    origin, direction = camera.ray(150, 260, 400, 300)
    assert direction == pytest.approx(camera.rotation.T @ np.array([0, 0, -1.0]))
    point = origin + direction * 7.5
    xy, depth = camera.project(point[None], 400, 300)
    assert xy[0] == pytest.approx([150, 260], abs=1e-6) and depth[0] == pytest.approx(7.5)
    clip = camera.projection(400 / 300) @ np.append(camera.view()[:3] @ np.append(near, 1), 1)
    assert (clip[:2] / clip[3]) == pytest.approx((after[0] / [200, 150] - 1) * [1, -1])
    assert centre_xy[0] == pytest.approx([200, 150], abs=1e-6)
    assert Camera.from_state(camera.state()).orthographic


def test_camera_state_round_trips() -> None:
    camera = Camera()
    camera.turn(10, 20)
    camera.pan(5, 5, 300)
    camera.zoom(0.5)
    again = Camera.from_state(camera.state())
    assert np.allclose(again.view(), camera.view())
    assert Camera.from_state({"rotation": "junk"}).distance == Camera().distance


def test_pick_finds_the_nearest_sphere_on_the_ray() -> None:
    centres = np.array([[0.0, 0, 0], [0, 0, -5], [3, 0, 0]])
    radii = np.array([1.0, 1.0, 1.0])
    origin, down = np.array([0.0, 0, 10]), np.array([0.0, 0, -1])
    assert pick_sphere(origin, down, centres, radii) == 0
    assert pick_sphere(origin, down, centres[1:], radii[1:]) == 0  # the one behind
    assert pick_sphere(origin, np.array([0.0, 1, 0]), centres, radii) is None


def test_measurements_match_rdkit() -> None:
    mol = embedded("CCCC")
    conf = mol.GetConformer()
    p = conf.GetPositions()
    assert distance(p[0], p[1]) == pytest.approx(rdMolTransforms.GetBondLength(conf, 0, 1))
    assert angle(p[0], p[1], p[2]) == pytest.approx(rdMolTransforms.GetAngleDeg(conf, 0, 1, 2))
    assert dihedral(*p[:4]) == pytest.approx(rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3))
    assert describe([p[0], p[1]]).endswith("Å")
    assert describe([p[0]]) == ""


def test_aligned_conformers_sit_on_the_first() -> None:
    mol = embedded("CCOc1ccccc1C(=O)N", 4)
    for conf in list(mol.GetConformers())[1:]:  # scatter them
        rdMolTransforms.TransformConformer(conf, _moved(conf.GetId()))
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]

    def spread(m):
        ref = m.GetConformer(0).GetPositions()[heavy]
        return max(
            np.sqrt(((c.GetPositions()[heavy] - ref) ** 2).sum(1).mean())
            for c in list(m.GetConformers())[1:]
        )

    aligned = aligned_conformers(mol)
    assert spread(aligned) < spread(mol)
    assert spread(aligned) < 2.0
    assert mol.GetConformer(1).GetPositions()[0] is not None  # the input is left alone


def _moved(seed: int) -> np.ndarray:
    m = np.eye(4)
    c, s = math.cos(seed), math.sin(seed)
    m[:3, :3] = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    m[:3, 3] = [10 * seed, -5, 3]
    return m
