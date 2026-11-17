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
from rdkit.Chem import AllChem, rdMolAlign
from scipy.spatial.transform import Rotation

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.mols import mol_from_blob, mol_to_blob
from chemsuite.core.table.columns import StructureColumn
from chemsuite.superpose import field, kabsch
from chemsuite.superpose.model import (
    ANCHORED,
    COLOR,
    COMBO,
    CORE,
    CORE_ATOMS,
    CORE_RMSD,
    FIELD,
    FIELD_SIMILARITY,
    MAX_RMSD,
    MCS,
    MEAN_RMSD,
    METHODS,
    O3A,
    SHAPE,
    Settings,
    blob_with_reference,
    reference_id_of,
)
from chemsuite.superpose.run import (
    align_conformers_rows,
    prepare_reference,
    superpose_one,
    superpose_rows,
)

SULFONAMIDE = "Cc1ccc(cc1)S(=O)(=O)Nc1ncccn1"
ATP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"
PROTEINS = Path(__file__).resolve().parents[1] / "data" / "protein"


def embedded(smiles: str, seed: int = 7, count: int = 1) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMultipleConfs(mol, count, randomSeed=seed)
    return mol


def moved(mol: Chem.Mol, seed: int = 0, shift=(3.0, -2.0, 5.0)) -> Chem.Mol:
    """A copy of *mol* turned and shifted at random (all conformers alike)."""
    out = Chem.Mol(mol)
    rotation = Rotation.random(random_state=seed).as_matrix()
    for conf in out.GetConformers():
        xyz = conf.GetPositions() @ rotation.T + np.asarray(shift)
        for i, (x, y, z) in enumerate(xyz):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    return out


def position(mol: Chem.Mol, symbol: str) -> np.ndarray:
    index = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == symbol)
    return mol.GetConformer().GetPositions()[index]


@pytest.mark.parametrize("method", METHODS)
def test_each_method_puts_a_moved_copy_back(method: str) -> None:
    mol = embedded(SULFONAMIDE)
    out, fit = superpose_one(moved(mol), prepare_reference(mol),
                             Settings(method=method, core="c1ccccc1"))  # fmt: skip
    assert kabsch.rmsd(out.GetConformer().GetPositions(), mol.GetConformer().GetPositions()) < 0.1
    if method in (SHAPE, COLOR, FIELD):
        assert min(v for k, v in fit.scores.items()) > 0.99
    if method == COLOR:
        assert fit.scores[COMBO] == pytest.approx(2.0, abs=0.01)


@pytest.mark.parametrize("method", sorted(ANCHORED))
def test_an_anchor_holds_its_atoms_on_the_reference_s(method: str) -> None:
    """Shape alone leaves 4-ethylaniline's nitrogen 0.2 Å off 4-methylaniline's."""
    reference = prepare_reference(embedded("Nc1ccc(C)cc1", seed=1))
    probe = moved(embedded("Nc1ccc(CC)cc1"))
    settings = Settings(method=method, anchor="Nc1ccccc1")
    out, fit = superpose_one(probe, reference, settings)
    assert np.linalg.norm(position(out, "N") - position(reference, "N")) < 0.1
    assert fit.scores[CORE_RMSD] < 0.1
    assert settings.scores[-1] == CORE_RMSD


def test_rows_without_the_anchor_are_fitted_without_it() -> None:
    reference = embedded(SULFONAMIDE)
    cells = [mol_to_blob(moved(reference)), "c1ccccc1C"]
    settings = Settings(method=SHAPE, anchor="c1ncccn1")
    batch = superpose_rows(cells, reference, "r", settings)
    assert all(batch.blobs) and batch.scores["Shape Tanimoto"][1] > 0
    assert batch.scores[CORE_RMSD][0] < 0.1 and np.isnan(batch.scores[CORE_RMSD][1])


def test_an_anchor_does_not_apply_to_the_substructure_methods() -> None:
    assert Settings(method=MCS, anchor="c1ccccc1").scores == (CORE_RMSD, CORE_ATOMS)


def test_mcs_fits_the_shared_ring_of_an_analog() -> None:
    reference = prepare_reference(embedded("CCc1ccccc1"))
    _out, fit = superpose_one(moved(embedded("Cc1ccccc1", seed=3)), reference, Settings(method=MCS))
    assert fit.scores[CORE_ATOMS] == 7  # the ring and the first carbon of the chain
    assert fit.scores[CORE_RMSD] < 0.1


def test_a_core_that_is_not_there_gives_no_fit() -> None:
    mol = embedded("CCO")
    with pytest.raises(ValueError):
        superpose_one(mol, prepare_reference(mol), Settings(method=CORE, core="c1ccccc1"))


def test_the_best_of_the_stored_conformers_is_kept() -> None:
    mol = embedded("CCCCCC(=O)NCCc1ccccc1", seed=11, count=5)
    reference = prepare_reference(Chem.Mol(mol, confId=3))
    _out, fit = superpose_one(moved(mol), reference, Settings(method=SHAPE))
    assert fit.conf_id == 3


def test_a_2d_row_gets_conformers() -> None:
    reference = prepare_reference(embedded(SULFONAMIDE))
    out, fit = superpose_one(Chem.MolFromSmiles(SULFONAMIDE), reference,
                             Settings(method=COLOR, conformers=10))  # fmt: skip
    assert out.GetNumConformers() == 1 and out.GetConformer().Is3D()
    assert fit.scores[COMBO] > 1.2


@pytest.mark.parametrize("method", [O3A, FIELD, COLOR])
def test_a_reference_without_hydrogens_works(method: str) -> None:
    """A protein's ligand comes as heavy atoms only; MMFF types and charges need hydrogens."""
    mol = embedded(SULFONAMIDE)
    reference = prepare_reference(Chem.RemoveHs(mol))
    out, _fit = superpose_one(moved(mol), reference, Settings(method=method))
    assert kabsch.rmsd(out.GetConformer().GetPositions(), mol.GetConformer().GetPositions()) < 0.3


def test_a_flat_reference_is_refused() -> None:
    with pytest.raises(ValueError, match="3D"):
        prepare_reference(Chem.MolFromMolBlock(Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))))


@pytest.mark.parametrize(("weight", "lands"), [(0.0, False), (0.5, True)])
def test_the_field_puts_the_ring_nitrogen_on_the_reference_s(weight: float, lands: bool) -> None:
    """The methyl of 4-methylpyridine misleads shape alone; the field turns the ring so its N
    lies on pyridine's."""
    reference = prepare_reference(embedded("c1ccncc1", seed=1))
    for seed in range(3):
        out, fit = superpose_one(moved(embedded("Cc1ccncc1", seed=2), seed=seed), reference,
                                 Settings(method=FIELD, field_weight=weight))  # fmt: skip
        distance = float(np.linalg.norm(position(out, "N") - position(reference, "N")))
        assert (distance < 0.5) == lands
        if lands:
            assert fit.scores[FIELD_SIMILARITY] > 0.8


def test_the_field_overlap_approximates_one_over_r() -> None:
    r = np.linspace(1.0, 8.0, 200)
    approx = (field.POTENTIAL_A[:, None] * np.exp(-field.POTENTIAL_B[:, None] * r**2)).sum(0)
    assert np.abs(approx * r - 1.0).max() < 0.035


def test_superposes_onto_a_protein_s_ligand() -> None:
    from chemsuite.macromolecule.chem import ligands
    from chemsuite.macromolecule.structure import io

    structure = io.read(PROTEINS / "1atp.cif.gz")
    residue = next(int(r) for r in ligands.ligand_residues(structure)
                   if str(structure.res_name[r]) == "ATP")  # fmt: skip
    crystal = ligands.to_mol(structure, residue, smiles=ATP).mol
    out, fit = superpose_one(Chem.MolFromSmiles(ATP), prepare_reference(crystal),
                             Settings(method=COLOR, conformers=30))  # fmt: skip
    assert fit.scores[COMBO] > 0.9  # ATP is very flexible: 1.07 measured with 30 conformers
    assert rdMolAlign.CalcRMS(Chem.RemoveHs(out), crystal) < 2.0  # in place, not refitted


def test_rows_name_their_reference_and_survive_the_structure_column() -> None:
    reference = embedded(SULFONAMIDE)
    cells = [mol_to_blob(moved(reference)), "c1ccccc1C", "", "not a smiles"]
    batch = superpose_rows(cells, reference, "superpose-1", Settings(method=SHAPE))
    assert [bool(b) for b in batch.blobs] == [True, True, False, False]
    assert batch.failed == 2 and batch.conformers[0] == 1
    assert batch.scores["Shape Tanimoto"][0] > 0.99 and np.isnan(batch.scores["Shape Tanimoto"][2])
    column = StructureColumn.from_parsed(batch.blobs, batch.displays)
    assert reference_id_of(column.mol(0)) == "superpose-1"
    assert reference_id_of(mol_from_blob(mol_to_blob(Chem.MolFromSmiles("C")))) is None


def test_the_pose_keeps_its_conformer_s_energy_above_the_ensemble_s_lowest() -> None:
    from chemsuite.conformers import energies
    from chemsuite.conformers.etkdg import ConformerSettings, generate

    reference = embedded(SULFONAMIDE)
    ensemble = generate(Chem.MolFromSmiles(SULFONAMIDE), ConformerSettings(count=8)).mol
    cells = [mol_to_blob(ensemble), "c1ccccc1C", mol_to_blob(moved(reference)), ""]
    batch = superpose_rows(cells, reference, "r", Settings(method=SHAPE, conformers=5))
    kept = int(batch.conformers[0]) - 1
    assert batch.delta_e[0] == pytest.approx(energies.delta_of(ensemble, kept))
    assert energies.delta_of(mol_from_blob(batch.blobs[0])) == pytest.approx(batch.delta_e[0])
    assert np.isfinite(batch.delta_e[1])  # a 2D row's own conformers carry energies
    assert np.isnan(batch.delta_e[2]) and np.isnan(batch.delta_e[3])  # none stored, blank


def test_the_reference_link_survives_a_blob_round_trip() -> None:
    blob = blob_with_reference(embedded("CCO"), "superpose-7")
    assert reference_id_of(mol_from_blob(blob)) == "superpose-7"


def test_cancel_keeps_the_rows_already_superposed() -> None:
    reference = embedded(SULFONAMIDE)
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done == 2:
            token.cancel()

    cells = [mol_to_blob(moved(reference, seed=s)) for s in range(4)]
    with pytest.raises(Cancelled) as stopped:
        superpose_rows(cells, reference, "r", Settings(method=SHAPE), progress=progress,
                       cancel=token, chunk_rows=1)  # fmt: skip
    batch = stopped.value.partial
    assert len(batch.blobs) == len(batch.displays) == len(batch.scores["Shape Tanimoto"]) == 2


def test_each_row_s_conformers_go_onto_its_first() -> None:
    mol = embedded("CCCCCC(=O)NCCc1ccccc1", seed=5, count=4)
    scattered = Chem.Mol(mol)
    for k, conf in enumerate(scattered.GetConformers()):
        xyz = conf.GetPositions() + k * 10.0
        for i, (x, y, z) in enumerate(xyz):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    batch = align_conformers_rows([mol_to_blob(scattered), "CCO"])
    out = mol_from_blob(batch.blobs[0])
    first = out.GetConformer(0).GetPositions().mean(axis=0)
    for conf in out.GetConformers():
        assert np.linalg.norm(conf.GetPositions().mean(axis=0) - first) < 1.5
    assert 0 < batch.scores[MEAN_RMSD][0] <= batch.scores[MAX_RMSD][0] < 5
    assert batch.blobs[1] == b""  # no 3D: nothing to align
