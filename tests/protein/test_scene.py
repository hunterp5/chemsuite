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
"""What the protein view draws: styles, geometry, look, colours and picking."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.molview.camera import Camera
from chemsuite.molview.pick import pick_sphere
from chemsuite.protein import colours, scene
from chemsuite.protein.model import ION, LIGAND, WATER, Structure
from chemsuite.protein.pick import ScreenPicker
from chemsuite.protein.select import select


def test_default_styles(kinase: Structure, nmr: Structure) -> None:
    styles = scene.default_styles(kinase)
    assert (styles[kinase.atom_kind == WATER] == scene.HIDDEN).all()
    assert (styles[kinase.atom_kind == LIGAND] == scene.BALL_AND_STICK).all()
    assert (styles[select(kinase, "protein")] == scene.WIRE).all()
    nmr_styles = scene.default_styles(nmr)
    polar = select(nmr, "hydrogen and name H")  # backbone amide hydrogens stay
    assert (nmr_styles[polar] == scene.WIRE).all()
    assert (nmr_styles[select(nmr, "hydrogen and name HA")] == scene.HIDDEN).all()


def test_only_the_shown_alternate_location_is_drawn(lysozyme: Structure) -> None:
    styles = scene.default_styles(lysozyme)
    assert not (styles[~lysozyme.active()] != scene.HIDDEN).any()


def test_geometry_draws_shown_atoms_and_bonds_between_them(kinase: Structure) -> None:
    styles = scene.default_styles(kinase)
    geometry = scene.geometry(kinase, styles)
    assert np.array_equal(geometry.atoms, np.flatnonzero(styles != scene.HIDDEN))
    assert np.allclose(geometry.spheres[:, :3], kinase.coords[0][geometry.atoms])
    ends = kinase.bonds[geometry.bonds]
    assert np.array_equal(ends, geometry.ends)
    assert (styles[ends] != scene.HIDDEN).all()
    assert np.allclose(geometry.cylinders[:, :3], kinase.coords[0][ends[:, 0]])
    assert np.allclose(geometry.cylinders[:, 4:7], kinase.coords[0][ends[:, 1]])
    # ball and stick: bonds thinner than the balls; ions are visible balls
    ligand = np.isin(geometry.atoms, np.flatnonzero(kinase.atom_kind == LIGAND))
    assert (geometry.spheres[ligand, 3] > 0.25).all()
    ion = np.isin(geometry.atoms, np.flatnonzero(kinase.atom_kind == ION))
    assert ion.sum() == 2 and (geometry.spheres[ion, 3] > 0.5).all()


def test_spacefill_atoms_have_no_bonds_and_hidden_atoms_nothing(crambin) -> None:
    styles = np.full(crambin.n_atoms, scene.SPACEFILL, np.uint8)
    styles[:10] = scene.HIDDEN
    geometry = scene.geometry(crambin, styles)
    assert len(geometry.cylinders) == 0
    assert len(geometry.spheres) == crambin.n_atoms - 10
    assert scene.geometry(crambin, np.zeros(crambin.n_atoms, np.uint8)).spheres.shape == (0, 4)


def test_look_rows_follow_the_geometry(kinase: Structure) -> None:
    geometry = scene.geometry(kinase, scene.default_styles(kinase))
    rgba = colours.colour(kinase)
    selected = select(kinase, "ligand")
    flags = scene.flags(kinase.n_atoms, selected=selected, ghost=select(kinase, "chain I"))
    spheres = scene.sphere_look(geometry, rgba, flags)
    assert spheres.shape == (len(geometry.atoms), scene.SPHERE_LOOK_FLOATS)
    assert np.allclose(spheres[:, :4], rgba[geometry.atoms])
    assert set(np.unique(spheres[:, 4])) == {0.0, scene.SELECTED, scene.GHOST}
    cylinders = scene.cylinder_look(geometry, rgba, flags)
    assert np.allclose(cylinders[:, 0:4], rgba[geometry.ends[:, 0]])
    assert np.allclose(cylinders[:, 9], flags[geometry.ends[:, 1]])


@pytest.mark.parametrize("scheme", list(colours.SCHEMES))
def test_every_scheme_colours_every_atom(kinase: Structure, scheme: str) -> None:
    rgba = colours.colour(kinase, scheme)
    assert rgba.shape == (kinase.n_atoms, 4) and rgba.dtype == np.float32
    assert (rgba >= 0).all() and (rgba <= 1).all()


def test_ligand_carbons_stand_out_from_the_protein(kinase: Structure) -> None:
    rgba = colours.colour(kinase)
    ligand_carbon = rgba[select(kinase, "ligand and elem C")]
    protein_carbon = rgba[select(kinase, "protein and elem C")]
    assert np.allclose(ligand_carbon[:, :3], colours.LIGAND_CARBON)
    assert (ligand_carbon[:, 1] > protein_carbon[:, 1].max()).all()
    oxygen = rgba[select(kinase, "elem O")]
    assert (oxygen[:, 0] > 0.9).all() and (oxygen[:, 1] < 0.2).all()


def test_the_screen_picker_agrees_with_ray_testing_every_sphere(kinase: Structure) -> None:
    geometry = scene.geometry(kinase, scene.default_styles(kinase))
    centres, radii = geometry.spheres[:, :3].astype(float), geometry.pick_radii.astype(float)
    camera = Camera()
    camera.orient(centres, 1.25)
    camera.fit(centres, radii, 1.25)
    width, height = 500, 400
    picker = ScreenPicker(camera, centres, radii, width, height)
    rng = np.random.default_rng(3)
    hits = 0
    for x, y in rng.uniform((0, 0), (width, height), size=(300, 2)):
        origin, direction = camera.ray(x, y, width, height)
        expected = pick_sphere(origin, direction, centres, radii)
        assert picker.pick(x, y) == (-1 if expected is None else expected)
        hits += expected is not None
    assert hits > 50


def test_the_picker_ignores_what_the_slab_cuts_away(crambin: Structure) -> None:
    geometry = scene.geometry(crambin, np.full(crambin.n_atoms, scene.SPACEFILL, np.uint8))
    centres, radii = geometry.spheres[:, :3].astype(float), geometry.pick_radii.astype(float)
    camera = Camera()
    camera.fit(centres, radii, 1.0)
    whole = ScreenPicker(camera, centres, radii, 300, 300)
    front = whole.pick(150, 150)
    _, depth = camera.project(centres[front : front + 1], 300, 300)
    cut = ScreenPicker(camera, centres, radii, 300, 300, slab=(float(depth[0]) + 1.0, 1e9))
    assert cut.pick(150, 150) not in (front, -1)


def test_another_model_moves_the_drawing_and_keeps_the_look(nmr: Structure) -> None:
    from chemsuite.protein.display import Display, compose

    shown = Display("1D3Z", "1D3Z", nmr)
    shown.show(np.ones(nmr.n_atoms, bool), scene.STICK)
    shown.selected[:40] = True
    before = compose([shown])
    assert len(before.geometry.atoms)
    assert shown.set_model(1)
    after = compose([shown])
    atoms = after.geometry.atoms
    assert np.allclose(after.geometry.spheres[:, :3], nmr.xyz(1)[atoms])
    assert not np.allclose(before.geometry.spheres[:, :3], after.geometry.spheres[:, :3])
    assert not np.allclose(before.cartoon.vertices, after.cartoon.vertices)
    assert np.array_equal(before.sphere_look, after.sphere_look)  # colours and selection stay
    assert np.allclose(after.pick_centres, nmr.xyz(1)[after.pick_atoms])
    assert not shown.set_model(1) and not shown.set_model(2) and not shown.set_model(-1)


def test_the_model_shown_is_saved_and_a_bad_one_is_ignored(nmr: Structure) -> None:
    from chemsuite.protein.display import Display

    shown = Display("1D3Z", "1D3Z", nmr)
    shown.set_model(1)
    again = Display("1D3Z", "1D3Z", nmr)
    again.set_state(shown.state())
    assert again.model == 1
    for bad in (7, -1, "1", None):
        again.set_state({**shown.state(), "model": bad})
        assert again.model == 0


def test_a_deferring_display_hands_out_slow_work_and_keeps_the_last_result(kinase) -> None:
    from chemsuite.protein.display import Display

    shown = Display("1ATP", "1ATP", kinase, deferred=True)
    assert shown.pending() == []  # the cartoon of 350 residues is quick: built inline
    assert len(shown.cartoon_mesh().triangles)
    shown.show_surface(np.ones(kinase.n_atoms, bool))
    assert not len(shown.surface_mesh().triangles)  # nothing built yet
    (first,) = shown.pending()
    assert first.part == "surface"
    built = first()
    assert shown.accept("surface", first.version, built)
    assert shown.surface_mesh() is built and shown.pending() == []
    shown.show_surface(select(kinase, "chain I"), on=False)
    assert shown.surface_mesh() is built  # the last one stays until the new one arrives
    (second,) = shown.pending()
    assert not shown.accept("surface", first.version, built)  # stale: dropped
    assert shown.accept("surface", second.version, second())
    shown.show_surface(np.ones(kinase.n_atoms, bool), on=False)
    assert shown.pending() == [] and not len(shown.surface_mesh().triangles)  # empty: inline


def test_deferred_work_reads_copies_taken_when_it_was_handed_out(kinase) -> None:
    from chemsuite.protein.display import Display

    shown = Display("1ATP", "1ATP", kinase, deferred=True)
    shown.show_interactions(kinase.atom_kind == LIGAND)
    (build,) = shown.pending()
    query = build.args[1].copy()
    shown.query[:] = False
    assert np.array_equal(build.args[1], query) and len(build())


def test_without_deferring_everything_is_built_at_once(kinase) -> None:
    from chemsuite.protein.display import Display

    shown = Display("1ATP", "1ATP", kinase)
    shown.show_surface(np.ones(kinase.n_atoms, bool))
    assert shown.pending() == [] and len(shown.surface_mesh().triangles)


def test_hydrogens_and_transparency_are_saved_and_bad_ones_ignored(nmr: Structure) -> None:
    from chemsuite.protein.display import SURFACE_TRANSPARENCY, Display

    shown = Display("1D3Z", "1D3Z", nmr)
    shown.show_hydrogens("all")
    shown.transparency = 0.2
    again = Display("1D3Z", "1D3Z", nmr)
    again.set_state(shown.state())
    assert again.hydrogens == "all" and again.transparency == 0.2
    again.set_state({**shown.state(), "hydrogens": "some", "transparency": 1.5})
    assert again.hydrogens == "polar" and again.transparency == SURFACE_TRANSPARENCY


def test_a_hydrogen_is_drawn_as_its_atom_is(nmr: Structure) -> None:
    styles = np.full(nmr.n_atoms, scene.HIDDEN, np.uint8)
    amide = select(nmr, "resi 5 and name N")
    styles[amide] = scene.SPACEFILL
    shown = scene.shown_styles(nmr, styles, nmr.active(), "polar")
    assert (shown[select(nmr, "resi 5 and name H")] == scene.SPACEFILL).all()
    assert (shown[select(nmr, "resi 5 and name HA")] == scene.HIDDEN).all()
