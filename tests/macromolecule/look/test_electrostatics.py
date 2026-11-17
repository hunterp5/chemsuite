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
"""Coloring surfaces by electrostatics: the charges, the potential, and the display's part."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.look import display as displays
from chemsuite.macromolecule.look import electrostatics as es
from chemsuite.macromolecule.look import surface
from chemsuite.macromolecule.look.display import ELECTROSTATICS, Display
from chemsuite.macromolecule.structure.model import AMINO, Structure
from chemsuite.macromolecule.structure.select import select

COARSE = surface.SurfaceOptions(quality="coarse")
NET = {"LYS": 1, "ARG": 1, "ASP": -1, "GLU": -1, "HIS": 0, "SER": 0, "CYS": 0, "GLY": 0}


def _residue_charges(structure: Structure, charges: np.ndarray) -> np.ndarray:
    return np.add.reduceat(charges, structure.res_start[:-1])


def test_every_hydrogen_of_the_table_has_its_heavy_atom() -> None:
    for residue, names in es.table().items():
        heavy = [n for n in names if not n.startswith("H")]
        for hydrogen, parent in es._parents(residue).items():
            assert parent in heavy, (residue, hydrogen)
    assert es._parents("ARG")["HH21"] == "NH2" and es._parents("ASN")["HD21"] == "ND2"
    assert es._parents("NALA")["H3"] == "N" and es._parents("DA")["H61"] == "N6"
    assert es._parents("RA")["HO'2"] == "O2'" and es._parents("DA5")["H5T"] == "O5'"


def test_an_unprepared_protein_gets_united_atom_charges_of_whole_residues(kinase) -> None:
    charges = es.charges(kinase)
    assert not charges[kinase.element == 1].any()  # no hydrogens to carry any
    per_residue = _residue_charges(kinase, charges)
    amino = np.flatnonzero(kinase.res_kind == AMINO)
    chains = kinase.chain_of[amino]
    ends = set(amino[np.r_[True, chains[1:] != chains[:-1]]]) | set(
        amino[np.r_[chains[1:] != chains[:-1], True]]
    )
    for r in amino:
        name = str(kinase.res_name[r])
        if name in NET and r not in ends:
            assert per_residue[r] == pytest.approx(NET[name], abs=1e-3), (r, name)
    # a chain's first residue carries the N-terminal amine's +1, its last the carboxylate's -1
    first = int(amino[0])
    name = str(kinase.res_name[first])
    if name in NET:
        assert per_residue[first] == pytest.approx(NET[name] + 1, abs=1e-3)
    assert charges[select(kinase, "water")].sum() == 0
    assert charges[select(kinase, "ion")] == pytest.approx([2.0, 2.0])  # Mn2+, from the table


def test_prepared_and_unprepared_histidines_differ(nmr: Structure) -> None:
    """The NMR model has its hydrogens: its His has HD1 only (HID). Without them, it is the
    default HIE, united atom; either way neutral."""
    prepared = es.charges(nmr)
    bare = nmr.subset(nmr.element != 1)
    unprepared = es.charges(bare)
    (his,) = [r for r in range(nmr.n_residues) if nmr.res_name[r] == "HIS"]

    def nd1(structure: Structure, charges: np.ndarray) -> float:
        span = structure.atoms_of(his)
        return float(charges[span][structure.name[span] == "ND1"][0])

    assert nd1(nmr, prepared) == pytest.approx(-0.3811)  # HID's, its HD1 its own atom
    assert nd1(bare, unprepared) == pytest.approx(-0.5432)  # HIE's, no hydrogen on it
    for structure, charges in ((nmr, prepared), (bare, unprepared)):
        assert _residue_charges(structure, charges)[his] == pytest.approx(0.0, abs=1e-3)


def test_a_ligand_gets_gasteiger_charges_that_add_up(kinase) -> None:
    charges = es.charges(kinase)
    atp = select(kinase, "ligand")
    assert np.abs(charges[atp]).max() > 0.1  # not all zero
    assert charges[atp].sum() == pytest.approx(0.0, abs=1e-3)  # its template is neutral


def test_the_potential_follows_coulomb_with_a_4r_dielectric() -> None:
    centres = np.zeros((1, 3))
    points = np.array([[5.0, 0, 0], [0, 10.0, 0], [0, 0, 20.0], [0.5, 0, 0]])
    values = es.potential(points, centres, np.array([1.0]))
    expected = es.COULOMB / (4 * np.array([5.0, 10.0]) ** 2)
    assert values[:2] == pytest.approx(expected, rel=1e-5)
    assert values[2] == 0  # beyond the 15 Å cutoff
    assert values[3] == pytest.approx(es.COULOMB / 4)  # r held at 1 Å
    pair = es.potential(points[:1], np.array([[0, 0, 0], [10.0, 0, 0]]), np.array([1.0, -1.0]))
    assert pair[0] == pytest.approx(0.0, abs=1e-5)  # halfway between opposite charges
    rgb = es.colour(np.array([-20.0, 0.0, 5.0, 20.0]))
    assert rgb[0] == pytest.approx(es.RED) and rgb[1] == pytest.approx(es.WHITE)
    assert rgb[3] == pytest.approx(es.BLUE) and rgb[2, 2] > rgb[2, 0]  # halfway to blue


def test_lysines_make_their_surface_blue_and_acids_red(kinase) -> None:
    protein = np.flatnonzero(select(kinase, "protein") & (kinase.element != 1))
    mesh = surface.build(kinase, protein, 0, COARSE)
    values = es.surface_potential(kinase, mesh).values
    tree_points = mesh.positions

    def near(names: str) -> list[float]:
        out = []
        for atom in np.flatnonzero(select(kinase, names)):
            close = np.linalg.norm(tree_points - kinase.coords[0][atom], axis=1) < 4.0
            if close.sum() > 5:  # an exposed one
                out.append(float(values[close].mean()))
        return out

    basic, acidic = (
        near("resn LYS and name NZ"),
        near("(resn ASP and name OD1) or (resn GLU and name OE1)"),
    )
    assert len(basic) > 5 and len(acidic) > 5
    assert np.median(basic) > 2 and np.median(acidic) < -2


def test_a_display_hands_out_the_potential_once_its_surface_is_built(kinase) -> None:
    shown = Display("1ATP", "1ATP", kinase, deferred=True, surface_options=COARSE)
    shown.show_surface(select(kinase, "chain I"))
    shown.colour_surface(ELECTROSTATICS)
    assert [b.part for b in shown.pending()] == ["surface"]  # the potential waits for it
    (build,) = shown.pending()
    assert shown.accept("surface", build.version, build())
    (potential,) = shown.pending()
    assert potential.part == "potential"
    white = shown.drawn_potential()
    assert len(white) == len(shown.surface_mesh()) and not white.any()  # 0 until it arrives
    assert shown.accept("potential", potential.version, potential())
    assert shown.drawn_potential().any() and shown.pending() == []
    shown.show_surface(select(kinase, "chain I"), on=False)  # a new surface: a new potential
    assert not shown.accept("potential", potential.version, potential())


def test_electrostatic_colors_follow_carving_dots_and_a_selection(kinase) -> None:
    shown = Display("1ATP", "1ATP", kinase, surface_options=COARSE)
    shown.show_surface(select(kinase, "chain I"))
    shown.colour_surface(ELECTROSTATICS)
    whole = shown.drawn_potential()
    scene = displays.compose([shown])
    assert np.isfinite(scene.surface_potential).all()
    assert scene.surface_look[:, :3] == pytest.approx(es.colour(whole))
    shown.selected[:] = True
    again = displays.relook(scene)  # a selection keeps the potential's colors
    assert again.surface_look[:, :3] == pytest.approx(scene.surface_look[:, :3])
    assert (again.surface_look[:, 4] > 0).all()
    centre = kinase.coords[0][select(kinase, "chain I and resi 10")]
    shown.carve_surface(centre, 5.0)
    _, _, used = shown._drawn(False)
    assert shown.drawn_potential() == pytest.approx(whole[used])
    shown.set_surface_style(displays.DOTS)
    mesh, _, source = shown._drawn(False)
    assert len(shown.drawn_potential()) == len(mesh)
    assert shown.drawn_potential() == pytest.approx(whole[source])
    shown.colour_surface(displays.FLAT)
    assert shown.drawn_potential() is None
    assert np.isnan(displays.compose([shown]).surface_potential).all()
