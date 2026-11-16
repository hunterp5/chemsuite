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

"""Prepared structures, made once for the simulation tests."""

from __future__ import annotations

from functools import cache

import pytest

from chemsuite.protein.prepare import pipeline
from tests.protein.conftest import load


@cache
def _prepared(name: str, waters: str) -> pipeline.Prepared:
    return pipeline.prepare(load(name), pipeline.Options(waters=waters))


@pytest.fixture
def crambin() -> pipeline.Prepared:
    """Crambin prepared: 46 residues, three disulfides, no ligand."""
    return _prepared("1crn.cif.gz", pipeline.NONE)


@pytest.fixture
def lysozyme() -> pipeline.Prepared:
    """T4 lysozyme L99A with its ligand, prepared: a His, a C-terminal NME cap."""
    return _prepared("3htb.cif.gz", pipeline.NONE)


@pytest.fixture
def kinase() -> pipeline.Prepared:
    """PKA with ATP and two Mn, prepared: phosphorylated residues no force field here takes."""
    return _prepared("1atp.cif.gz", pipeline.NONE)


@pytest.fixture
def egfr() -> pipeline.Prepared:
    """The EGFR pocket of 5GNK, prepared, its waters kept: a covalent inhibitor, caps."""
    return _prepared("5gnk_pocket.cif.gz", pipeline.ALL)
