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
"""The search box."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.macromolecule.dock import box as box_tools
from chemsuite.macromolecule.dock.box import LARGEST, PADDING, SMALLEST, Box


def test_a_box_pads_the_atoms_it_is_built_from() -> None:
    points = np.array([[0.0, 0.0, 0.0], [20.0, 4.0, 0.0]])
    box = box_tools.from_points(points, padding=4.0)
    assert box.centre == pytest.approx((10.0, 2.0, 0.0))
    # 20 Å, 4 Å and 0 Å wide, each padded by 4 Å on both sides; the flat side clamped
    assert box.size == pytest.approx((28.0, 12.0, SMALLEST))
    assert box.holds(points) and not box.holds(np.array([[30.0, 0.0, 0.0]]))


def test_a_box_is_never_smaller_than_a_ligand_can_turn_in_nor_bigger_than_vina_takes() -> None:
    one = box_tools.from_points(np.zeros((1, 3)))
    assert one.size == pytest.approx((SMALLEST, SMALLEST, SMALLEST))
    wide = box_tools.from_points(np.array([[-500.0, 0, 0], [500.0, 0, 0]]))
    assert wide.size[0] == pytest.approx(LARGEST)
    with pytest.raises(ValueError, match="at least one atom"):
        box_tools.from_points(np.zeros((0, 3)))


def test_the_box_of_a_bound_ligand_holds_it(pocket) -> None:
    box = box_tools.from_residue(pocket.structure, pocket.residue)
    atoms = pocket.structure.xyz()[pocket.structure.atoms_of(pocket.residue)]
    assert box.holds(atoms)
    span = pocket.structure.atoms_of(pocket.residue)
    mask = np.zeros(pocket.structure.n_atoms, bool)
    mask[span] = True
    indices = np.arange(span.start, span.stop)
    # a mask and a list of indices mean the same thing
    assert box_tools.from_atoms(pocket.structure, mask).size == pytest.approx(box.size)
    assert box_tools.from_atoms(pocket.structure, indices).centre == pytest.approx(box.centre)
    assert box_tools.from_atoms(pocket.structure, mask, padding=0.0).volume < box.volume


def test_vinas_command_line_takes_the_centre_and_the_size() -> None:
    arguments = Box((1.0, 2.0, 3.0), (10.0, 11.0, 12.0)).arguments()
    assert arguments == [
        "--center_x=1.000",
        "--center_y=2.000",
        "--center_z=3.000",
        "--size_x=10.000",
        "--size_y=11.000",
        "--size_z=12.000",
    ]
    assert PADDING > 0
