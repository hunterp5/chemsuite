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
"""Protein test files (``tests/data/protein``, made by ``scripts/build_protein_fixtures.py``)."""

from __future__ import annotations

from functools import cache
from pathlib import Path

import pytest

from chemsuite.protein import io
from chemsuite.protein.model import Structure

DATA = Path(__file__).resolve().parents[1] / "data" / "protein"


@cache
def load(name: str) -> Structure:
    return io.read(DATA / name)


@pytest.fixture
def crambin() -> Structure:
    return load("1crn.cif.gz")


@pytest.fixture
def lysozyme() -> Structure:
    """T4 lysozyme L99A/M102Q with 2-propylphenol (JZ4), waters and alternate locations."""
    return load("3htb.cif.gz")


@pytest.fixture
def kinase() -> Structure:
    """Protein kinase A with ATP, two Mn and an inhibitor peptide (chains E and I)."""
    return load("1atp.cif.gz")


@pytest.fixture
def thrombin() -> Structure:
    """Insertion codes, and PPACK (0G6) bound covalently to His57 and Ser195."""
    return load("1ppb.cif.gz")


@pytest.fixture
def nmr() -> Structure:
    """Ubiquitin, two NMR models, with hydrogens."""
    return load("1d3z_2models.cif.gz")


@pytest.fixture
def egfr() -> Structure:
    """EGFR T790M cut to the pocket of LXX-6-34 (80U): covalent on Cys797, a halogen bond."""
    return load("5gnk_pocket.cif.gz")
