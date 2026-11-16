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
"""The bundled AutoDock Vina binary."""

from __future__ import annotations

import json

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from chemsuite.core.jobs import CancelToken
from chemsuite.macromolecule.dock import pdbqt, vina


def _output(affinity: float = -9.191) -> str:
    """An output PDBQT as Vina writes one: Meeko's ligand with Vina's result line added."""
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=3)
    result = f"REMARK VINA RESULT:   {affinity:8.3f}      0.000      1.234"
    return "\n".join(["MODEL 1", result, *pdbqt.ligand(mol).splitlines(), "ENDMDL"]) + "\n"


def test_the_binary_ships_with_its_licence_and_its_hashes() -> None:
    licence = vina.BIN / "LICENSE-AutoDock-Vina.txt"
    assert "Apache License" in licence.read_text(encoding="utf-8")
    assert "AutoDock Vina" in (vina.BIN / "NOTICE").read_text(encoding="utf-8")
    pinned = json.loads((vina.BIN / "checksums.json").read_text(encoding="utf-8"))
    assert set(pinned["files"]) == set(vina.BINARIES.values())


def test_the_bundled_binary_runs_and_says_what_it_is() -> None:
    assert vina.executable().exists()
    assert vina.version().startswith("AutoDock Vina v")


def test_an_unknown_scoring_function_is_refused(pocket) -> None:
    with pytest.raises(vina.VinaError, match="scoring"):
        vina.dock_one("", pocket.mol, pocket.box, vina.Options(scoring="glide"))


def test_a_cancelled_run_stops_the_binary(pocket) -> None:
    from chemsuite.core.jobs import Cancelled

    cancel = CancelToken()
    cancel.cancel()
    with pytest.raises(Cancelled):
        vina.dock_one(
            pocket.receptor.pdbqt, pocket.mol, pocket.box, vina.Options(exhaustiveness=32), cancel
        )


def test_the_affinities_come_from_the_output() -> None:
    (pose,) = vina._poses(_output(-9.191))
    assert (pose.rank, pose.affinity, pose.rmsd_lb, pose.rmsd_ub) == (1, -9.191, 0.0, 1.234)
    assert Chem.MolToSmiles(Chem.RemoveHs(pose.mol)) == "CCO"


def test_an_output_whose_scores_are_missing_is_an_error() -> None:
    silent = _output().replace("REMARK VINA RESULT", "REMARK NOTHING")
    with pytest.raises(vina.VinaError, match="scores"):
        vina._poses(silent)


def test_docking_a_ligand_back_into_its_own_pocket_finds_its_pose(pocket) -> None:
    """The real binary, at the lowest exhaustiveness there is: poses, ranked, inside the box."""
    poses = vina.dock_one(
        pocket.receptor.pdbqt, pocket.mol, pocket.box, vina.Options(exhaustiveness=1, poses=3)
    )
    assert poses
    assert [p.rank for p in poses] == list(range(1, len(poses) + 1))
    assert poses == sorted(poses, key=lambda p: p.affinity)
    for pose in poses:
        heavy = Chem.RemoveHs(pose.mol)
        assert pocket.box.holds(heavy.GetConformer().GetPositions())
        assert heavy.GetNumAtoms() == Chem.RemoveHs(pocket.mol).GetNumAtoms()
    assert poses[0].affinity < -5  # a 35-atom inhibitor in its own pocket


def test_a_pose_can_be_scored_where_it_stands(pocket) -> None:
    crystal = vina.score_only(pocket.receptor.pdbqt, pocket.mol)
    assert -20 < crystal < 0
    pushed = Chem.Mol(pocket.mol)
    conformer = pushed.GetConformer()
    for i in range(pushed.GetNumAtoms()):
        x, y, z = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (x + 30.0, y, z))
    assert vina.score_only(pocket.receptor.pdbqt, pushed) > crystal  # nowhere near the protein


@pytest.mark.bench
def test_the_top_pose_is_the_crystal_pose(pocket) -> None:
    """Redocking 80U into 5GNK at Vina's own default: within 2 Å of where the crystal has it."""
    poses = vina.dock_one(pocket.receptor.pdbqt, pocket.mol, pocket.box)
    rmsd = rdMolAlign.CalcRMS(Chem.RemoveHs(poses[0].mol), Chem.RemoveHs(pocket.mol))
    assert rmsd < 2.0
