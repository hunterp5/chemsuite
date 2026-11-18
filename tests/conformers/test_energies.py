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

import math

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.conformers import conforge
from chemsuite.conformers import energies as en
from chemsuite.conformers.etkdg import ConformerSettings, generate, generate_chunk
from chemsuite.core.mols import mol_from_blob, mol_to_blob

FLEXIBLE = "CCCCCCCCO"


def _embedded(count: int) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=7)
    return mol


def test_boltzmann_populations() -> None:
    assert en.boltzmann([]) == []
    pops = en.boltzmann([0.0, 0.0])
    assert pops == pytest.approx([0.5, 0.5])
    pops = en.boltzmann([0.0, 1.0, 3.0])
    assert sum(pops) == pytest.approx(1.0)
    assert pops[0] > pops[1] > pops[2]
    # 1 kcal/mol at 298.15 K is a factor of about 5.4
    assert pops[0] / pops[1] == pytest.approx(math.exp(1 / (en.GAS_CONSTANT * 298.15)))


def test_annotate_and_read_back() -> None:
    mol = _embedded(3)
    en.annotate(mol, [-10.0, -8.5, -7.0], "MMFF94")
    rows = en.rows(mol)
    assert [r.number for r in rows] == [1, 2, 3]
    assert [r.delta for r in rows] == pytest.approx([0.0, 1.5, 3.0])
    assert [r.energy for r in rows] == pytest.approx([-10.0, -8.5, -7.0])
    assert sum(r.population for r in rows) == pytest.approx(1.0)
    assert en.force_field(mol) == "MMFF94"
    assert en.delta_of(mol, 1) == pytest.approx(1.5)
    en.clear(mol)
    assert en.rows(mol) == [] and en.force_field(mol) == "" and math.isnan(en.delta_of(mol))


def test_nothing_is_stored_without_a_force_field_or_with_a_missing_energy() -> None:
    mol = _embedded(2)
    en.annotate(mol, [1.0, 2.0], "")
    en.annotate(mol, [1.0, math.nan], "UFF")
    en.annotate(mol, [1.0], "UFF")  # one energy for two conformers
    assert en.rows(mol) == []


def test_blob_keeps_conformer_energies_but_not_molecule_fields() -> None:
    mol = _embedded(2)
    mol.SetProp("IC50", "12")  # a data field read from an SD file
    en.annotate(mol, [0.0, 2.0], "MMFF94")
    back = mol_from_blob(mol_to_blob(mol))
    assert [r.delta for r in en.rows(back)] == pytest.approx([0.0, 2.0])
    assert not back.HasProp("IC50")
    assert mol.HasProp("IC50")  # the caller's molecule is left alone
    plain = Chem.MolFromSmiles("CCO")
    plain.SetProp("IC50", "12")
    assert not mol_from_blob(mol_to_blob(plain)).HasProp("IC50")


def test_etkdg_energies_follow_the_kept_conformers() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=20, energy_window=3))
    rows = en.rows(result.mol)
    assert len(rows) == result.count > 1
    assert [r.energy for r in rows] == pytest.approx(result.energies)
    assert rows[0].delta == 0.0 and all(r.delta <= 3 for r in rows)
    assert en.force_field(result.mol) == "MMFF94"


def test_etkdg_without_optimization_stores_no_energies() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=3, optimize=None))
    assert result.count and en.rows(result.mol) == []


def test_chunk_blobs_carry_the_energies() -> None:
    [(blob, _display, count, span)] = generate_chunk((ConformerSettings(count=6), [FLEXIBLE]))
    rows = en.rows(mol_from_blob(blob))
    assert len(rows) == count
    assert rows[-1].delta == pytest.approx(span)


@pytest.mark.skipif(not conforge.available(), reason="cdpkit is not installed")
def test_conforge_energies_follow_its_sort() -> None:
    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=8, method="conforge"))
    rows = en.rows(result.mol)
    assert [r.energy for r in rows] == pytest.approx(result.energies)
    assert en.force_field(result.mol) == conforge.FORCE_FIELD


def test_energies_reopen_with_the_document(tmp_path) -> None:
    from chemsuite.core.document import load_document, save_document
    from chemsuite.core.table import Table
    from chemsuite.core.table.columns import StructureColumn

    result = generate(Chem.MolFromSmiles(FLEXIBLE), ConformerSettings(count=5))
    table = Table()
    column = StructureColumn.from_parsed([mol_to_blob(result.mol)], [FLEXIBLE])
    table.add_rows({"Mol 3D": column}, undoable=False)
    save_document(table, tmp_path / "doc.chemsuite")
    reopened, _extras = load_document(tmp_path / "doc.chemsuite")
    rows = en.rows(reopened.column("Mol 3D").mol(0))
    assert [r.energy for r in rows] == pytest.approx(result.energies)


def test_a_structure_changed_by_standardize_drops_its_energies() -> None:
    from chemsuite.analysis.standardize import Options, standardize_chunk

    cells = [
        mol_to_blob(generate(Chem.MolFromSmiles(s), ConformerSettings(count=3)).mol)
        for s in ("CCCC(=O)[O-]", "CCCC(=O)O")  # neutralized, unchanged
    ]
    result = standardize_chunk((Options(), cells))
    changed, kept = (mol_from_blob(b) for b in result.blobs)
    assert Chem.MolToSmiles(Chem.RemoveHs(changed)) == "CCCC(=O)O"
    assert changed.GetNumConformers() and en.rows(changed) == []
    assert en.rows(kept)  # nothing changed: the energies still hold


def test_a_charged_major_species_drops_its_energies() -> None:
    from chemsuite.ionisation import profile
    from chemsuite.ionisation.predict import predict_chunk
    from tests.ionisation.test_profile import SmilesStub

    blob = mol_to_blob(generate(Chem.MolFromSmiles("CCCC(=O)O"), ConformerSettings(count=3)).mol)
    profile.use(SmilesStub())  # a model in place of the weights
    try:
        result = predict_chunk((7.4, [blob]))
    finally:
        profile.use(None)
    species = mol_from_blob(result.species_blobs[0])
    assert Chem.MolToSmiles(Chem.RemoveHs(species)) == "CCCC(=O)[O-]"
    assert en.rows(species) == []
