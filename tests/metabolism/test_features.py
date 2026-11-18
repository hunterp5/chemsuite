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

"""The 43 columns the sites-of-metabolism model reads."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.metabolism import ISOFORMS
from chemsuite.metabolism import features as f


def test_the_column_layout_is_the_models_input_contract():
    # The exported model indexes columns by position, so the names, their count and the
    # order of the isoform block are all part of the archive's format.
    assert f.N_ATOM_FEATURES == 34
    assert f.N_FEATURES == 43
    assert f.FEATURE_NAMES[:10] == (
        "el_C", "el_N", "el_O", "el_S", "el_P", "el_F", "el_Cl", "el_Br", "el_I",
        "el_other",
    )  # fmt: skip
    assert f.FEATURE_NAMES[34:] == tuple(f"isoform_{name}" for name in ISOFORMS)


def test_a_block_has_one_row_per_heavy_atom():
    block = f.atom_block(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    assert block.shape == (11, f.N_ATOM_FEATURES)
    assert np.isfinite(block).all()


def test_the_columns_say_what_they_claim():
    # toluene: the methyl is benzylic, terminal, sp3, with three hydrogens
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    block = f.atom_block(mol)
    column = {name: i for i, name in enumerate(f.ATOM_FEATURE_NAMES)}
    assert block[0, column["is_benzylic_or_allylic_CH"]] == 1.0
    assert block[0, column["is_terminal"]] == 1.0
    assert block[0, column["hyb_sp3"]] == 1.0
    assert block[0, column["n_hydrogens"]] == 3.0
    assert block[0, column["aromatic_neighbours"]] == 1.0
    # and the ring carbon it hangs off is aromatic, in one six-ring, with no hydrogen
    assert block[1, column["is_aromatic"]] == 1.0
    assert block[1, column["min_ring_size"]] == 6.0
    assert block[1, column["n_rings_at_atom"]] == 1.0
    assert block[1, column["n_hydrogens"]] == 0.0


def test_an_unlisted_element_lands_in_el_other():
    block = f.atom_block(Chem.MolFromSmiles("C[Si](C)(C)C"))
    column = f.ATOM_FEATURE_NAMES.index("el_other")
    assert block[1, column] == 1.0
    assert block[0, column] == 0.0


def test_non_finite_descriptors_become_zero():
    assert f._safe(float("nan")) == 0.0
    assert f._safe(float("inf")) == 0.0
    assert f._safe(-2.5) == -2.5


def test_the_design_matrix_runs_isoform_major():
    block = f.atom_block(Chem.MolFromSmiles("CCO"))
    matrix = f.design_matrix(block, ("3A4", "2D6"))
    assert matrix.shape == (6, f.N_FEATURES)
    # each isoform's copy holds the same atom block
    assert np.array_equal(matrix[:3, :34], block)
    assert np.array_equal(matrix[3:, :34], block)
    # and its own one-hot, the same on every row of the copy
    assert np.array_equal(matrix[0, 34:], f.isoform_one_hot("3A4"))
    assert np.array_equal(matrix[2, 34:], f.isoform_one_hot("3A4"))
    assert np.array_equal(matrix[3, 34:], f.isoform_one_hot("2D6"))


def test_an_unknown_isoform_is_refused_rather_than_scored_as_none():
    with pytest.raises(ValueError, match="unknown isoform"):
        f.isoform_one_hot("2C18")


def test_every_atom_has_environments_to_measure_the_domain_with():
    mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    found = f.atom_environments(mol)
    assert len(found) == mol.GetNumAtoms()
    assert all(pairs for pairs in found)
    radii = {radius for pairs in found for radius, _ in pairs}
    assert radii <= set(f.DOMAIN_RADII)
    # the two symmetry-equivalent pairs of ring CH share their environments exactly
    assert found[5] == found[10]
    assert found[6] == found[9]
