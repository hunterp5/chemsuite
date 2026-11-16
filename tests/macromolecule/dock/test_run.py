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
"""Many ligands: order, progress, failures and cancel.

The binary itself is tested in ``test_vina.py``; here it is stood in for, so these tests are
about what :mod:`chemsuite.macromolecule.dock.run` does with its results.
"""

from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.macromolecule.dock import run, vina
from chemsuite.macromolecule.dock.box import Box
from chemsuite.macromolecule.dock.vina import Pose

BOX = Box((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))


def _mol(smiles: str = "CCO", seed: int = 1) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


@pytest.fixture
def stub(monkeypatch) -> list[Chem.Mol]:
    """Vina stood in for: each ligand gets one pose whose affinity says which ligand it was."""
    order: list[Chem.Mol] = []

    def fake(receptor, mol, box, options=None, cancel=None, cache=None):
        order.append(mol)
        if mol.GetNumAtoms() < 4:
            raise vina.VinaError("too small to dock")
        return [Pose(mol, 1, -float(mol.GetNumAtoms())), Pose(mol, 2, -1.0)]

    monkeypatch.setattr(run, "dock_one", fake)
    return order


def test_every_ligand_comes_back_in_the_order_it_was_given(stub) -> None:
    mols = [_mol("CCO"), _mol("c1ccccc1"), _mol("CCCCO")]
    seen: list[tuple[int, int]] = []
    results = run.dock_many("receptor", mols, BOX, progress=lambda k, n: seen.append((k, n)))
    assert [d.index for d in results] == [0, 1, 2]
    assert [len(d.poses) for d in results] == [2, 2, 2]
    assert [round(-d.best.affinity) for d in results] == [m.GetNumAtoms() for m in mols]
    assert seen[-1] == (3, 3) and len(seen) == 3


def test_a_ligand_vina_will_not_dock_is_an_error_beside_the_others(stub) -> None:
    results = run.dock_many("receptor", [_mol("CCO"), Chem.AddHs(Chem.MolFromSmiles("O"))], BOX)
    assert results[0].poses and not results[0].error
    assert results[1].error == "too small to dock" and results[1].best is None


def test_a_ligand_with_no_3d_is_an_error_not_a_crash() -> None:
    (result,) = run.dock_many("receptor", [Chem.MolFromSmiles("CCO")], BOX)
    assert "3D" in result.error


def test_nothing_to_dock_is_nothing_to_do(stub) -> None:
    assert run.dock_many("receptor", [], BOX) == []


def test_a_cancel_keeps_the_ligands_that_finished(monkeypatch) -> None:
    cancel = CancelToken()

    def fake(receptor, mol, box, options=None, cancel=None, cache=None):
        if mol.GetNumAtoms() > 9 and cancel is not None:  # benzene, the second, cancels the run
            cancel.cancel()
        return [Pose(mol, 1, -1.0)]

    monkeypatch.setattr(run, "dock_one", fake)
    mols = [_mol("CCO"), _mol("c1ccccc1"), _mol("CCO"), _mol("CCO")]
    with pytest.raises(Cancelled) as stopped:
        run.dock_many("receptor", mols, BOX, cancel=cancel, workers=1)
    kept = stopped.value.partial
    # what finished is kept, in order, and the rest was never started
    assert [d.index for d in kept] == list(range(len(kept)))
    assert 0 < len(kept) < len(mols)


def test_rescoring_gives_one_pose_a_ligand(monkeypatch) -> None:
    def fake(receptor, mol, box=None, options=None, cancel=None, cache=None):
        return -7.5

    monkeypatch.setattr(run, "score_only", fake)
    results = run.rescore_many("receptor", [_mol("CCO"), _mol("CCCO")])
    assert [d.best.affinity for d in results] == [-7.5, -7.5]
    assert all(len(d.poses) == 1 and d.best.rank == 1 for d in results)
