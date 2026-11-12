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

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.conformers.etkdg import (
    ConformerSettings,
    generate,
    generate_chunk,
    generate_conformers,
)
from chemsuite.core.depict import depiction_mol
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_from_blob, mol_to_blob
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn
from chemsuite.descriptors.compute import compute_chunk
from chemsuite.io.sdf import write_sdf

FLEXIBLE = "CCCCCCCCO"


def test_generate_sorts_by_energy_within_window() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=12, energy_window=3))
    assert result.force_field == "MMFF94"
    assert 1 < result.count <= 12
    assert result.energies == sorted(result.energies)
    assert result.energies[-1] - result.energies[0] <= 3
    assert [c.GetId() for c in result.mol.GetConformers()] == list(range(result.count))
    assert result.mol.GetConformer().Is3D()
    assert any(a.GetAtomicNum() == 1 for a in result.mol.GetAtoms())  # explicit Hs kept


def test_generate_is_reproducible_and_leaves_input_alone() -> None:
    mol = Chem.MolFromSmiles(FLEXIBLE)
    a = generate(mol, ConformerSettings(count=4))
    b = generate(mol, ConformerSettings(count=4))
    assert np.allclose(a.mol.GetConformer().GetPositions(), b.mol.GetConformer().GetPositions())
    assert mol.GetNumConformers() == 0 and mol.GetNumAtoms() == 9


def test_without_optimisation_keeps_embedding_order() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=3, optimize=None))
    assert result.force_field == "" and np.isnan(result.energies).all()


def test_uff_fallback_and_rigid_molecules() -> None:
    assert generate(Chem.MolFromSmiles("C[Fe]")).force_field == "UFF"
    assert generate(Chem.MolFromSmiles("c1ccccc1"), ConformerSettings(count=10)).count == 1


def test_settings_validate() -> None:
    with pytest.raises(ValueError):
        ConformerSettings(count=0)
    with pytest.raises(ValueError):
        ConformerSettings(optimize="amber")


def test_chunk_handles_blank_invalid_and_text() -> None:
    blob = mol_to_blob(Chem.MolFromSmiles("CCO"))
    rows = generate_chunk((ConformerSettings(count=2), [blob, "CCN", None, "nonsense"]))
    assert [r[1] for r in rows] == ["CCO", "CCN", "", ""]
    assert mol_from_blob(rows[0][0]).GetNumConformers() >= 1
    assert rows[2][0] == b"" and np.isnan(rows[2][2])


def test_batch_progress_and_cancel() -> None:
    seen = []
    batch = generate_conformers(
        ["CCO", "CCC", "CCN"],
        ConformerSettings(count=2),
        progress=lambda d, t: seen.append(d),
        chunk_rows=1,
    )
    assert seen == [1, 2, 3] and batch.failed == 0 and (batch.counts >= 1).all()
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        generate_conformers(["CCO"], cancel=token)


def test_cancel_keeps_the_conformers_already_made() -> None:
    token = CancelToken()

    def progress(done: int, total: int) -> None:
        if done == 2:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        generate_conformers(["CCO", "CCN", "CCC", "CCCl"], ConformerSettings(count=2),
                            progress=progress, cancel=token)  # fmt: skip
    batch = stopped.value.partial
    assert len(batch.blobs) == len(batch.displays) == len(batch.counts) == 2


def test_conformers_survive_the_structure_column() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=5))
    column = StructureColumn.from_parsed([mol_to_blob(result.mol)], [FLEXIBLE])
    assert column.mol(0).GetNumConformers() == result.count


def test_3d_structures_draw_without_hydrogens() -> None:
    mol = generate(Chem.MolFromSmiles("CCO")).mol
    drawn = depiction_mol(mol)
    assert drawn.GetNumAtoms() == 3 and not drawn.GetConformer().Is3D()


def test_shape_descriptors_need_3d_and_2d_ones_ignore_hydrogens() -> None:
    mol = generate(Chem.MolFromSmiles("OCCc1ccccc1")).mol
    keys = ("NPR1", "PBF", "SMILES", "HeavyAtomCount")
    out = compute_chunk((keys, [mol_to_blob(mol), "OCCc1ccccc1"]))
    assert 0 < out["NPR1"][0] < 1 and np.isnan(out["NPR1"][1])
    assert out["SMILES"] == ["OCCc1ccccc1", "OCCc1ccccc1"]
    assert out["HeavyAtomCount"] == [9.0, 9.0]


def test_sdf_writes_each_conformer(tmp_path: Path) -> None:
    mol = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=5)).mol
    count = mol.GetNumConformers()
    assert count > 1
    table = Table()
    table.add_rows(
        {
            "Structure 3D": StructureColumn.from_parsed([mol_to_blob(mol), b""], [FLEXIBLE, ""]),
            "ID": ["A", "B"],
        },
        undoable=False,
    )
    path = tmp_path / "out.sdf"
    write_sdf(table, path, all_conformers=True)
    records = list(Chem.SDMolSupplier(str(path), removeHs=False))
    assert len(records) == count + 1  # every conformer, plus the blank row once
    assert [r.GetProp("Conformer") for r in records[:count]] == [str(i + 1) for i in range(count)]
    first = mol.GetConformer(0).GetPositions()
    assert np.allclose(records[0].GetConformer().GetPositions(), first, atol=1e-3)
    assert not np.allclose(records[1].GetConformer().GetPositions(), first, atol=1e-3)
    assert {r.GetProp("ID") for r in records[:count]} == {"A"}
    write_sdf(table, path)  # default: one record per row
    assert len(list(Chem.SDMolSupplier(str(path)))) == 2
