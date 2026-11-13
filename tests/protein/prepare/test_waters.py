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

from chemsuite.protein.prepare import pipeline, sites, waters
from chemsuite.protein.structure.model import LIGAND, WATER
from tests.protein.conftest import load
from tests.protein.prepare.test_ligand_sites import _protonated


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


def test_only_waters_near_ligands_get_hydrogens() -> None:
    s = load("5mo0_pocket.cif.gz")
    s = s.subset(s.element != 1)
    done = pipeline.prepare(s, pipeline.Options(ph=7.5))
    out = done.structure
    ligand = out.xyz()[(out.atom_kind == LIGAND) & (out.element != 1)]
    tree = cKDTree(ligand)
    for r in np.flatnonzero(out.res_kind == WATER).tolist():
        span = out.atoms_of(r)
        oxygen = [i for i in range(span.start, span.stop) if out.element[i] == 8][0]
        near = bool(tree.query_ball_point(out.xyz()[oxygen], waters.REACH))
        n_h = int((out.element[span] == 1).sum())
        assert n_h == (2 if near else 0)
    assert any(f.label == "Waters" for f in done.report)
    off = pipeline.prepare(s, pipeline.Options(ph=7.5, orient_waters=False))
    assert not (off.structure.element[off.structure.atom_kind == WATER] == 1).any()
