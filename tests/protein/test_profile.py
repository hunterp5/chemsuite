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
"""Interaction profiles of ligand rows: counts by kind and the residues touched."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.protein.chem import interactions, ligands, profile

from .test_poses import moved


def row_of(structure, smiles: str = "") -> profile.Row:
    rows = ligands.ligand_rows(structure, "P")
    mol = Chem.MolFromMolBlock(rows[ligands.LIGAND_COLUMN][0])
    return profile.Row(structure, rows[ligands.LIGAND_ID_COLUMN][0], mol, smiles)


def count(found: profile.Profile, kind: int) -> int:
    return found.counts[list(profile.COUNT_COLUMNS).index(kind)]


def test_a_kinase_inhibitor_row(egfr) -> None:
    found = profile.profile(row_of(egfr))
    parts = found.contacts.split("; ")
    assert "MET793 HB" in parts and "LEU788 X" in parts and "GLN791 HB" in parts
    assert count(found, interactions.HALOGEN) == 1 and count(found, interactions.HBOND) >= 3
    assert "HOH" not in found.contacts  # waters show as the bridges they make


def test_atp_with_its_lysine_and_metals(kinase) -> None:
    found = profile.profile(row_of(kinase))
    assert "LYS72 SB" in found.contacts.split("; ")
    assert "MN" in found.contacts and count(found, interactions.METAL) == 4
    assert len(found.counts) == len(profile.COUNT_COLUMNS)


def test_a_moved_pose_is_profiled_where_it_is(kinase) -> None:
    row = row_of(kinase)
    far = profile.Row(kinase, row.label, moved(row.mol, 30.0))  # out of the protein
    assert profile.profile(far).counts == (0,) * len(profile.COUNT_COLUMNS)
    nameless = profile.Row(kinase, "", row.mol)  # no ligand ID: its pose alone
    assert profile.profile(nameless) == profile.profile(row)
    assert profile.profile(profile.Row(kinase)) is None  # names nothing


def test_cancelling_keeps_the_rows_finished(kinase, egfr) -> None:
    rows = [row_of(egfr), None, row_of(kinase)]
    cancel = CancelToken()
    seen: list[tuple[int, int]] = []

    def progress(done: int, total: int) -> None:
        seen.append((done, total))
        if done == 2:
            cancel.cancel()

    with pytest.raises(Cancelled) as stopped:
        profile.for_rows(rows, progress, cancel)
    assert len(stopped.value.partial) == 2 and stopped.value.partial[1] is None
    assert seen == [(1, 3), (2, 3)]
