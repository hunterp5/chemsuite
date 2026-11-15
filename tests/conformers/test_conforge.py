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

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from chemsuite.conformers import conforge
from chemsuite.conformers.etkdg import ConformerSettings, generate, generate_chunk
from chemsuite.core.mols import mol_from_blob

needs_cdpkit = pytest.mark.skipif(not conforge.available(), reason="cdpkit is not installed")

IBUPROFEN = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
# stereocentres and an E double bond, to check they survive the trip through CDPKit
STEREO = r"C[C@H](N)[C@@H](O)/C=C/c1ccccc1"
MACROCYCLE = "C1CCCCCCCCCCC(=O)NCCCCCC1"


def _settings(count: int = 10, energy_window: float | None = 10.0) -> ConformerSettings:
    return ConformerSettings(count=count, energy_window=energy_window, method="conforge")


def test_unknown_method_is_refused() -> None:
    with pytest.raises(ValueError):
        ConformerSettings(method="omega")


@needs_cdpkit
def test_conformers_keep_atom_order_and_stereo() -> None:
    mol = Chem.MolFromSmiles(STEREO)
    result = generate(mol, _settings(count=20))
    assert result.mol is not None and result.count > 1
    assert result.force_field == "MMFF94s"
    expected = Chem.AddHs(mol)
    assert [a.GetSymbol() for a in result.mol.GetAtoms()] == [
        a.GetSymbol() for a in expected.GetAtoms()
    ]
    for conformer in result.mol.GetConformers():
        posed = Chem.Mol(result.mol, confId=conformer.GetId())
        Chem.AssignStereochemistryFrom3D(posed, confId=conformer.GetId())
        assert Chem.MolToInchi(posed) == Chem.MolToInchi(mol)
    assert mol.GetNumConformers() == 0  # the input is left alone


@needs_cdpkit
def test_energies_sorted_and_settings_honoured() -> None:
    result = generate(Chem.MolFromSmiles(IBUPROFEN), _settings(count=8, energy_window=4.0))
    assert 1 < result.count <= 8
    assert result.energies == sorted(result.energies)
    assert result.energies[-1] - result.energies[0] <= 4.0 + 1e-6
    assert [c.GetId() for c in result.mol.GetConformers()] == list(range(result.count))
    heavy = Chem.RemoveHs(result.mol)
    ids = [c.GetId() for c in heavy.GetConformers()]
    rms = [rdMolAlign.GetBestRMS(heavy, heavy, ids[0], i) for i in ids[1:]]
    assert min(rms) >= 0.5 - 0.05  # prune_rms, measured CONFORGE's way within a tolerance


@needs_cdpkit
def test_macrocycle_gets_conformers() -> None:
    result = generate(Chem.MolFromSmiles(MACROCYCLE), _settings(count=10))
    assert result.count > 1
    assert np.ptp(result.mol.GetConformer().GetPositions(), axis=0).min() > 0.5  # 3D


@needs_cdpkit
def test_failure_is_an_error_not_an_exception() -> None:
    mol = Chem.MolFromSmiles("[Xe]")  # no MMFF94 parameters
    result = generate(mol, _settings())
    assert result.mol is None and result.error


@needs_cdpkit
def test_batch_path_runs_conforge() -> None:
    rows = generate_chunk((_settings(count=5), [IBUPROFEN, None]))
    blob, display, count, span = rows[0]
    assert count >= 1 and mol_from_blob(blob).GetNumConformers() == count
    assert np.isnan(rows[1][2])


def test_missing_cdpkit_is_reported(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(conforge, "available", lambda: False)
    result = generate(Chem.MolFromSmiles(IBUPROFEN), _settings())
    assert result.mol is None and "cdpkit" in result.error
