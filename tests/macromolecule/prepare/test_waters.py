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

"""Pocket waters: orientations, and hydrogens written for the waters near ligands only."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial import cKDTree

from chemsuite.macromolecule.prepare import pipeline, sites, waters
from chemsuite.macromolecule.structure.model import LIGAND, WATER
from tests.macromolecule.conftest import load
from tests.macromolecule.prepare.test_ligand_sites import _protonated


def test_each_orientation_is_a_water() -> None:
    p = _protonated("5mo0_pocket.cif.gz")
    found = waters.find(p)
    assert found and all(site.kind == sites.WATER and not site.titrating for site in found)
    for site in found:
        assert 2 <= len(site.states) <= waters.MAX_STATES
        for st in site.states:
            o, h1, h2 = st.xyz[0]
            assert np.linalg.norm(h1 - o) == pytest.approx(waters.OH)
            assert np.linalg.norm(h2 - o) == pytest.approx(waters.OH)
            cos = np.dot(h1 - o, h2 - o) / waters.OH**2
            assert np.degrees(np.arccos(cos)) == pytest.approx(104.52, abs=0.01)


def test_a_water_points_a_hydrogen_at_each_pair_of_neighbours() -> None:
    a, b = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    h1, h2 = waters._orientations([a, b])[0]
    bisector = (a + b) / np.linalg.norm(a + b)
    assert np.dot(h1, bisector) == pytest.approx(np.dot(h2, bisector))
    assert np.cross(h1, h2) @ np.cross(a, b) != 0  # in their plane, one each side


def test_waters_near_ligands_are_sites_and_every_water_gets_hydrogens() -> None:
    s = load("5mo0_pocket.cif.gz")
    s = s.subset(s.element != 1)
    done = pipeline.prepare(s, pipeline.Options(ph=7.5))
    out = done.structure
    ligand = out.xyz()[(out.atom_kind == LIGAND) & (out.element != 1)]
    tree = cKDTree(ligand)
    solved = {c.site.residue for c in done.choices if c.site.kind == sites.WATER}
    far = 0
    for r in np.flatnonzero(out.res_kind == WATER).tolist():
        span = out.atoms_of(r)
        oxygen = [i for i in range(span.start, span.stop) if out.element[i] == 8][0]
        near = bool(tree.query_ball_point(out.xyz()[oxygen], waters.REACH))
        assert (r in solved) == near
        assert int((out.element[span] == 1).sum()) == 2
        far += not near
    assert far  # the pocket has waters away from the ligand too
    assert any(f.label == "Waters" for f in done.report)
    off = pipeline.prepare(s, pipeline.Options(ph=7.5, orient_waters=False))
    o = off.structure
    assert not any(c.site.kind == sites.WATER for c in off.choices)
    assert (o.element[o.atom_kind == WATER] == 1).sum() == 2 * (o.res_kind == WATER).sum()


def test_far_waters_take_tip3p_geometry_and_point_at_their_neighbours() -> None:
    s = load("5mo0_pocket.cif.gz")
    s = s.subset((s.element != 1) & (s.atom_kind != LIGAND))  # no ligand: every water is far
    out = waters.complete(s)
    x = out.xyz()
    assert set(s.uid.tolist()) <= set(out.uid.tolist())
    polar = np.flatnonzero(np.isin(out.element, [7, 8]))
    tree = cKDTree(x[polar])
    bonded = 0
    hydrogens = 0
    for r in np.flatnonzero(out.res_kind == WATER).tolist():
        span = out.atoms_of(r)
        names = out.name[span].tolist()
        assert sorted(names) == ["H1", "H2", "O"]
        o, h1, h2 = (x[span][names.index(n)] for n in ("O", "H1", "H2"))
        assert np.linalg.norm(h1 - o) == pytest.approx(waters.OH, abs=1e-3)
        angle = np.degrees(np.arccos(np.dot(h1 - o, h2 - o) / waters.OH**2))
        assert angle == pytest.approx(104.52, abs=0.1)
        for h in (h1, h2):
            hydrogens += 1
            near = [polar[k] for k in tree.query_ball_point(h, waters.BOND)]
            bonded += any(out.residue_of[j] != r for j in near)
    assert bonded / hydrogens > 0.5  # most hydrogens reach an N or O
    assert waters.complete(out) is out  # nothing left to do
