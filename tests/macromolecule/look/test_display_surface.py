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
"""A display's surface settings: kind, color, transparency, carving; saved, and kept by edits."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.look import colours, surface
from chemsuite.macromolecule.look import display as displays
from chemsuite.macromolecule.look.display import FLAT, Display
from chemsuite.macromolecule.look.surface import SurfaceOptions
from chemsuite.macromolecule.structure.model import Structure
from chemsuite.macromolecule.structure.select import select

COARSE = SurfaceOptions(quality="coarse")


def _surfaced(kinase: Structure, options: SurfaceOptions = COARSE) -> Display:
    shown = Display("1ATP", "1ATP", kinase, surface_options=options)
    shown.show_surface(select(kinase, "chain I"))
    return shown


def test_new_surfaces_are_molecular_and_old_documents_open_as_gaussian(kinase) -> None:
    shown = Display("1ATP", "1ATP", kinase)
    assert shown.surface_options.kind == surface.SES
    saved = shown.state()
    del saved["surface_options"]  # as a document saved before surfaces had kinds
    again = Display("1ATP", "1ATP", kinase)
    again.set_state(saved)
    assert again.surface_options == surface.BLOB_OPTIONS


def test_surface_settings_are_saved_and_bad_ones_ignored(kinase) -> None:
    shown = _surfaced(kinase, SurfaceOptions(kind=surface.SAS, probe=1.6, smoothing="high"))
    shown.colour_surface(FLAT, (0.2, 0.4, 0.6))
    shown.pocket_transparency = 0.6
    shown.carve_surface(kinase.coords[0][:3], 8.0)
    again = Display("1ATP", "1ATP", kinase)
    again.set_state(shown.state())
    assert again.surface_options == shown.surface_options
    assert again.surface_scheme == FLAT and again.surface_rgb == pytest.approx((0.2, 0.4, 0.6))
    assert again.pocket_transparency == 0.6 and again.carve_distance == 8.0
    assert np.allclose(again.carve_points, kinase.coords[0][:3], atol=1e-3)
    bad = {
        **shown.state(),
        "surface_rgb": [2, 0, 0],
        "carve": "x",
        "carve_distance": -1,
        "pocket_transparency": True,
        "surface_options": {"kind": "ses", "probe": "big"},
    }
    again = Display("1ATP", "1ATP", kinase)
    again.set_state(bad)
    assert again.surface_rgb == displays.FLAT_RGB and again.carve_points is None
    assert again.carve_distance == displays.CARVE_DISTANCE
    assert again.pocket_transparency == again.transparency
    assert again.surface_options == surface.BLOB_OPTIONS


def test_an_edit_keeps_the_surface_settings(kinase) -> None:
    shown = _surfaced(kinase, SurfaceOptions(kind=surface.VDW, quality="coarse"))
    shown.colour_surface(FLAT, (1.0, 0.0, 0.0))
    shown.transparency, shown.pocket_transparency = 0.2, 0.8
    shown.carve_surface(kinase.coords[0][:5], 3.0)
    first = int(kinase.residue_of[np.flatnonzero(select(kinase, "chain E"))[0]])
    edited = Display("1ATP", "1ATP", kinase.subset(kinase.residue_of != first))
    edited.carry(shown)
    for name in ("surface_options", "surface_scheme", "surface_rgb", "transparency",
                 "pocket_transparency", "carve_distance"):  # fmt: skip
        assert getattr(edited, name) == getattr(shown, name), name
    assert np.array_equal(edited.carve_points, shown.carve_points)
    assert edited.surface.sum() == shown.surface.sum()


def test_carving_cuts_the_built_surface_without_rebuilding_it(kinase) -> None:
    shown = _surfaced(kinase)
    whole = shown.surface_mesh()
    centre = kinase.coords[0][select(kinase, "chain I and resi 10")]
    shown.carve_surface(centre, 5.0)
    carved = shown.surface_mesh()
    assert 0 < len(carved.triangles) < len(whole.triangles)
    assert shown.surface_mesh() is carved  # kept while nothing changes
    assert shown._part("surface") is whole  # the whole surface was not rebuilt
    shown.carve_surface(None)
    assert shown.surface_mesh() is whole


def test_options_rebuild_the_surface_and_the_pocket(kinase) -> None:
    shown = _surfaced(kinase)
    before = shown.surface_mesh()
    shown.set_surface_options(COARSE)  # the same: nothing to rebuild
    assert shown.surface_mesh() is before
    shown.set_surface_options(SurfaceOptions(kind=surface.SAS, quality="coarse"))
    after = shown.surface_mesh()
    assert after is not before and len(after.triangles)


def test_surface_and_pocket_colors_are_kept_apart(kinase) -> None:
    shown = _surfaced(kinase)
    shown.transparency, shown.pocket_transparency = 0.2, 0.6
    surface_rgba = shown.surface_colours()
    pocket_rgba = shown.surface_colours(pocket=True)
    assert surface_rgba[:, 3] == pytest.approx(0.8) and pocket_rgba[:, 3] == pytest.approx(0.4)
    # both kept: asking again does not color every atom over
    assert shown.surface_colours() is surface_rgba
    assert shown.surface_colours(pocket=True) is pocket_rgba
    shown.colour_surface(FLAT, (0.0, 0.5, 1.0))
    flat = shown.surface_colours()
    assert flat[:, :3] == pytest.approx(np.tile([0.0, 0.5, 1.0], (kinase.n_atoms, 1)))
    shown.colour_surface(colours.CHAIN)
    by_chain = colours.colour(kinase, colours.CHAIN)[:, :3]
    assert shown.surface_colours()[:, :3] == pytest.approx(by_chain)
    shown.colour_surface("nope")
    assert shown.surface_scheme == colours.CHAIN


def test_styles_route_surfaces_to_filled_and_wire_meshes(kinase) -> None:
    from chemsuite.macromolecule.look import mesh

    shown = _surfaced(kinase)
    other = _surfaced(kinase)  # a second copy of the protein, drawn as a wire frame
    other.set_surface_style(displays.MESH)
    scene = displays.compose([shown, other])
    whole = shown.surface_mesh()
    assert len(scene.surface) == len(whole) and len(scene.surface_wire) == len(whole)
    assert (scene.surface_wire.owner >= kinase.n_atoms).all()  # the second protein's atoms
    assert len(scene.surface_wire_look) == len(scene.surface_wire)
    shown.set_surface_style(displays.DOTS)
    dotted = displays.compose([shown, other])
    assert len(dotted.surface) == len(mesh.dots(whole))
    assert shown.drawn_mesh()[0] is shown.drawn_mesh()[0]  # the dots are kept
    shown.selected[:] = True
    again = displays.relook(dotted)
    assert (again.surface_look[:, 4] > 0).all() and (again.surface_wire_look[:, 4] == 0).all()
    shown.set_surface_style("nope")
    assert shown.surface_style == displays.DOTS


def test_styles_are_saved_and_kept_by_edits(kinase) -> None:
    shown = _surfaced(kinase)
    shown.set_surface_style(displays.MESH)
    shown.set_surface_style(displays.DOTS, pocket=True)
    again = Display("1ATP", "1ATP", kinase)
    again.set_state(shown.state())
    assert (again.surface_style, again.pocket_style) == (displays.MESH, displays.DOTS)
    again.set_state({**shown.state(), "surface_style": 3})
    assert again.surface_style == displays.SOLID
    edited = Display("1ATP", "1ATP", kinase.subset(np.arange(kinase.n_atoms) > 0))
    edited.carry(shown)
    assert (edited.surface_style, edited.pocket_style) == (displays.MESH, displays.DOTS)
