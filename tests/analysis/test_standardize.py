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

"""Standardize Structures: each step, the defaults, stereo, failures and cancelling."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.standardize import DOCKING, Options, standardize, standardize_cells
from chemsuite.core.jobs import Cancelled, CancelToken

FULL = Options(remove_isotopes=True, tautomer=True, protonate=True)


def run(smiles: str, options: Options = DOCKING) -> tuple[str, list[str]]:
    outcome = standardize(Chem.MolFromSmiles(smiles), options)
    return Chem.MolToSmiles(outcome.mol), outcome.notes


@pytest.mark.parametrize(
    "smiles, want, note",
    [
        ("CCN(CC)CC.Cl", "CCN(CC)CC", "salts stripped"),  # hydrochloride
        ("CC(=O)[O-].[Na+]", "CC(=O)O", "neutralised"),  # sodium salt
        ("CC(=O)Oc1ccccc1C(=O)O.O", "CC(=O)Oc1ccccc1C(=O)O", "salts stripped"),  # hydrate
        ("[Na]OC(=O)c1ccccc1", "O=C(O)c1ccccc1", "cleaned up"),  # metal drawn bonded
        ("C[n+]1ccccc1[O-]", "Cn1ccccc1=O", "cleaned up"),  # charge-separated pyridone
        ("[NH3+]CCC(=O)[O-]", "NCCC(=O)O", "neutralised"),  # zwitterion
    ],
)
def test_default_steps(smiles, want, note):
    got, notes = run(smiles)
    assert got == Chem.CanonSmiles(want)
    assert note in notes


def test_quaternary_nitrogen_stays_charged():
    got, _ = run("CC[N+](C)(C)C.[Cl-]")
    assert got == "CC[N+](C)(C)C"


def test_optional_steps_are_off_by_default():
    assert run("Oc1ccccn1")[0] == "Oc1ccccn1"
    assert run("CN1CCCC1")[0] == "CN1CCCC1"
    assert run("[2H]C(Cl)Cl")[0] == "[2H]C(Cl)Cl"


def test_tautomer_and_isotopes():
    assert run("Oc1ccccn1", Options(tautomer=True)) == ("O=c1cccc[nH]1", ["tautomer changed"])
    assert run("[2H]C(Cl)Cl", Options(remove_isotopes=True))[0] == "ClCCl"


def test_protonation_at_ph_uses_molgpka():
    amine, notes = run("CN1CCC[C@H]1c1cccnc1", Options(protonate=True))  # nicotine
    assert amine == "C[NH+]1CCC[C@H]1c1cccnc1" and "net charge +1" in notes
    assert run("CC(=O)O", Options(protonate=True, ph=2.0))[0] == "CC(=O)O"
    assert run("CC(=O)O", Options(protonate=True, ph=7.4))[0] == "CC(=O)[O-]"


def test_protonation_never_makes_a_carbanion():
    ticarcillin = "CC1(C)SC2C(NC(=O)[C@H](C(=O)O)c3ccsc3)C(=O)N2C1C(=O)O"
    got, _ = run(ticarcillin, Options(protonate=True))
    mol = Chem.MolFromSmiles(got)
    assert all(a.GetFormalCharge() >= 0 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    assert got.count("[O-]") == 2


def test_protonation_keeps_counter_ions_when_salts_are_kept():
    got, _ = run("CN1CCCC1.Cl", Options(strip_salts=False, neutralise=False, protonate=True))
    assert set(got.split(".")) == {"C[NH+]1CCCC1", "Cl"}


def test_drawn_stereo_is_kept_by_every_step():
    for smiles in ["C[C@H](N)C(=O)O", "C/C=C/C(=O)O", "C[C@@H]1CCC(=O)C1"]:
        got, notes = run(smiles, FULL)
        mol = Chem.MolFromSmiles(got)
        assert not any("undefined" in n for n in notes), (smiles, notes)
        assert (
            Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
            == (
                Chem.FindMolChiralCenters(Chem.MolFromSmiles(smiles), useLegacyImplementation=False)
            )
            or "/" in smiles
        )


def test_undefined_stereo_is_reported():
    assert run("CC(N)C(=O)O")[1] == ["1 undefined stereocentre"]
    assert run("CC=CC")[1] == ["1 undefined double bond"]
    assert run("CC(O)C(N)CC")[1] == ["2 undefined stereocentres"]
    assert run("CN=[N+]=[N-]")[1] == []  # an azide is linear: no E/Z


def test_cells_failures_blanks_and_counts():
    results = standardize_cells(["CCN.Cl", None, "C1CC", "CCO", b""])
    assert results.texts == ["CCN", "", "", "CCO", ""]
    assert results.notes == [
        "salts stripped",
        "",
        "not a valid structure",
        "",
        "not a valid structure",
    ]
    assert results.changed == 1 and results.failed == 2


def test_cancel_keeps_finished_rows():
    token = CancelToken()
    done: list[int] = []

    def progress(n, total):
        done.append(n)
        token.cancel()

    with pytest.raises(Cancelled) as caught:
        standardize_cells(["CCO"] * 100, progress=progress, cancel=token)
    partial = caught.value.partial
    assert partial is not None and 0 < len(partial) < 100
    assert partial.texts == ["CCO"] * len(partial)
