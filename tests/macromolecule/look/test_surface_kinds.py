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
"""The surface kinds (SES, SAS, vdW), their options, and the mesh helpers for the styles."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.look import mesh, surface
from chemsuite.macromolecule.look.surface import SurfaceOptions
from chemsuite.macromolecule.structure.model import Structure
from chemsuite.macromolecule.structure.select import select
from chemsuite.molview.elements import vdw_radii

from .test_mesh import edge_counts, normals_agree, volume

KINDS = (surface.SES, surface.SAS, surface.VDW)


def _area(m: mesh.Mesh) -> float:
    p = m.positions.astype(np.float64)
    t = m.triangles.astype(np.int64)
    return float(0.5 * np.linalg.norm(np.cross(p[t[:, 1]] - p[t[:, 0]],
                                               p[t[:, 2]] - p[t[:, 0]]), axis=1).sum())  # fmt: skip


@pytest.mark.parametrize("kind", KINDS)
def test_every_kind_is_closed_faces_out_and_belongs_to_its_atoms(kinase, kind) -> None:
    atoms = np.flatnonzero(select(kinase, "ligand"))
    m = surface.build(kinase, atoms, 0, SurfaceOptions(kind=kind))
    counts = edge_counts(m)
    assert (counts % 2 == 0).all()  # closed (an edge where the skin touches itself: four)
    assert (counts == 2).mean() > 0.99
    assert volume(m) > 0 and normals_agree(m) > 0.99
    assert np.isin(m.owner, atoms).all()


# measured at Normal quality (0.5 Å) with Low smoothing: the mean radius comes within 0.03 Å
# and the area within 4% (vdW), 1.5% (SAS) and 2% (SES)
@pytest.mark.parametrize(("kind", "probe", "area"), [(surface.VDW, 0.0, 0.05),
                                                     (surface.SAS, 1.4, 0.02),
                                                     (surface.SES, 0.0, 0.03)])  # fmt: skip
def test_a_lone_atom_is_a_sphere_of_the_right_size(kinase, kind, probe, area) -> None:
    (ion,) = np.flatnonzero(select(kinase, "ion"))[:1]
    radius = float(vdw_radii(kinase.element[[ion]])[0]) + probe
    m = surface.build(kinase, np.array([ion]), 0, SurfaceOptions(kind=kind))
    distance = np.linalg.norm(m.positions - kinase.coords[0][ion], axis=1)
    assert distance.mean() == pytest.approx(radius, abs=0.05)
    assert np.abs(distance - radius).max() < 0.15
    assert _area(m) == pytest.approx(4 * np.pi * radius**2, rel=area)


def test_the_ses_lies_between_the_atoms_and_the_probe_s_reach(lysozyme: Structure) -> None:
    heavy = np.flatnonzero((lysozyme.element != 1) & select(lysozyme, "protein"))
    volumes = {kind: volume(surface.build(lysozyme, heavy, 0, SurfaceOptions(kind=kind)))
               for kind in KINDS}  # fmt: skip
    assert volumes[surface.VDW] < volumes[surface.SES] < volumes[surface.SAS]


def _two_atoms(gap: float) -> Structure:
    from chemsuite.macromolecule.structure.io import parse

    lines = [
        f"HETATM{i + 1:5d}  C{i + 1}  UNL A   1    {x:8.3f}   0.000   0.000  1.00  0.00           C"
        for i, x in enumerate((0.0, gap))
    ]
    return parse("\n".join([*lines, "END"]) + "\n", "pdb")


def test_the_ses_fills_the_cleft_between_two_atoms() -> None:
    pair = _two_atoms(3.6)
    atoms = np.arange(pair.n_atoms)

    def waist(kind: str) -> float:
        m = surface.build(pair, atoms, 0, SurfaceOptions(kind=kind, smoothing="off"))
        middle = np.abs(m.positions[:, 0] - 1.8) < 0.3
        return float(np.linalg.norm(m.positions[middle][:, 1:], axis=1).min())

    # the vdW spheres (1.7 Å) pinch to about 0.8 Å between them; a 1.4 Å probe cannot get in
    assert waist(surface.VDW) < 1.1
    assert waist(surface.SES) > waist(surface.VDW) + 0.5


def test_smoothing_keeps_the_surface_closed_and_close(kinase: Structure) -> None:
    atoms = np.flatnonzero(select(kinase, "ligand"))
    rough = surface.build(kinase, atoms, 0, SurfaceOptions(smoothing="off"))
    smooth = surface.build(kinase, atoms, 0, SurfaceOptions(smoothing="high"))
    assert np.array_equal(rough.triangles, smooth.triangles)  # the same mesh, moved
    moved = np.linalg.norm(smooth.positions - rough.positions, axis=1)
    assert moved.mean() > 0 and moved.max() < 0.3
    assert (edge_counts(smooth) % 2 == 0).all()


def test_options_are_checked_and_saved() -> None:
    with pytest.raises(ValueError):
        SurfaceOptions(kind="nope")
    with pytest.raises(ValueError):
        SurfaceOptions(probe=0.1)
    options = SurfaceOptions(kind=surface.SAS, probe=1.6, quality="fine", smoothing="off")
    assert SurfaceOptions.from_state(options.state()) == options
    # a surface saved before there were kinds, or a damaged one, is the old Gaussian
    assert SurfaceOptions.from_state(None) == surface.BLOB_OPTIONS
    assert SurfaceOptions.from_state({"kind": "nope"}) == surface.BLOB_OPTIONS
    assert SurfaceOptions().kind == surface.SES  # new surfaces are molecular surfaces
    assert surface.BLOB_OPTIONS.spacing == surface.SPACING


def test_edges_dots_and_carving(kinase: Structure) -> None:
    atoms = np.flatnonzero(select(kinase, "ligand"))
    m = surface.build(kinase, atoms, 0, SurfaceOptions())
    lines = mesh.edges(m)
    assert len(lines) == len(np.unique(np.sort(lines, axis=1), axis=0))
    assert len(lines) == pytest.approx(1.5 * len(m.triangles), rel=0.01)  # closed: 3 per 2
    dotted = mesh.dots(m)
    assert 0 < len(dotted) // 6 < len(m) and len(dotted.triangles) == 8 * (len(dotted) // 6)
    assert np.isin(dotted.owner, m.owner).all()
    centre = m.positions.mean(axis=0)
    carved = mesh.carve(m, centre[None], 6.0)
    assert 0 < len(carved.triangles) < len(m.triangles)
    assert (np.linalg.norm(carved.positions - centre, axis=1) <= 6.0).all()
    assert int(carved.triangles.max()) < len(carved)
    assert len(mesh.carve(m, np.zeros((0, 3)), 6.0)) == 0
