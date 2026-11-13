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

"""MolGpKa port, pH-dependent properties and the batch predictor."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.ionisation import molgpka
from chemsuite.ionisation.molgpka import ACID, BASE, predict_sites, prepare
from chemsuite.ionisation.predict import predict_pka
from chemsuite.ionisation.profile import clogd, profile

GOLDEN = json.loads((Path(__file__).parent / "molgpka_golden.json").read_text())


def _upstream_prepared(smiles: str) -> Chem.Mol:
    """Upstream's preparation (no largest-fragment step), to compare like with like."""
    mol = rdMolStandardize.Uncharger().uncharge(Chem.MolFromSmiles(smiles))
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))


@pytest.mark.parametrize(("smiles", "expected"), GOLDEN, ids=[g[0] for g in GOLDEN])
def test_port_matches_upstream_molgpka(smiles, expected) -> None:
    """Golden values from upstream MolGpKa's own featurisation and the original .pth weights
    run with PyTorch Geometric: same sites, same pKa to within float rounding."""
    got = sorted([s.atom, s.kind, s.pka] for s in predict_sites(_upstream_prepared(smiles)))
    assert [(a, k) for a, k, _ in got] == [(a, k) for a, k, _ in expected]
    for (_, _, pka), (_, _, want) in zip(got, expected, strict=True):
        assert pka == pytest.approx(want, abs=1e-3)


def test_predictions_are_close_to_experiment_for_common_groups() -> None:
    def strongest(smiles: str, kind: str) -> float:
        values = [
            s.pka for s in predict_sites(prepare(Chem.MolFromSmiles(smiles))) if s.kind == kind
        ]
        return min(values) if kind == ACID else max(values)

    assert strongest("CC(=O)O", ACID) == pytest.approx(4.76, abs=0.5)  # acetic acid
    assert strongest("Oc1ccccc1", ACID) == pytest.approx(10.0, abs=0.5)  # phenol
    assert strongest("c1ccncc1", BASE) == pytest.approx(5.2, abs=0.5)  # pyridine
    assert strongest("C1CCNCC1", BASE) == pytest.approx(11.1, abs=0.7)  # piperidine


def test_an_atom_can_be_both_a_base_and_an_acid() -> None:
    """Aniline's nitrogen: a weak base, and through its hydrogens a very weak acid."""
    sites = predict_sites(prepare(Chem.MolFromSmiles("Nc1ccccc1")))
    assert [(s.atom, s.kind) for s in sites] == [(0, BASE), (0, ACID)]


def test_sites_on_a_molecule_as_it_is_keep_its_numbering() -> None:
    """``predict_on`` numbers sites by the molecule's own atoms, not the re-read SMILES's."""
    mol = Chem.MolFromSmiles("OC(=O)CCc1ccncc1")  # the acid's O first, the pyridine N late
    found = {(s.atom, s.kind): s.pka for s in molgpka.predict_on(mol)}
    canonical = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    expected = {s.kind: s.pka for s in predict_sites(canonical)}
    assert set(found) == {(0, ACID), (8, BASE)}
    assert mol.GetAtomWithIdx(8).GetSymbol() == "N"
    assert found[0, ACID] == pytest.approx(expected[ACID])
    assert found[8, BASE] == pytest.approx(expected[BASE])


def test_salts_and_charges_are_prepared_to_the_neutral_parent() -> None:
    salt = prepare(Chem.MolFromSmiles("[Na+].CC(=O)[O-]"))
    assert Chem.MolToSmiles(salt) == "CC(=O)O"
    assert clogd(Chem.MolFromSmiles("[Na+].CC(=O)[O-]")) == pytest.approx(
        clogd(Chem.MolFromSmiles("CC(=O)O"))
    )


def test_ph_dependent_properties() -> None:
    acid = profile(Chem.MolFromSmiles("CC(=O)O"))
    pka = acid.most_acidic
    assert acid.most_basic is None
    assert acid.charge(pka) == pytest.approx(-0.5)  # half ionised at its pKa
    assert acid.fraction_neutral(pka - 3) == pytest.approx(1 / 1.001, rel=1e-6)
    assert acid.logd(1.0, pka + 2) == pytest.approx(1.0 + math.log10(1 / 101))
    assert Chem.MolToSmiles(acid.major_species(7.4)) == "CC(=O)[O-]"
    assert Chem.MolToSmiles(acid.major_species(2.0)) == "CC(=O)O"

    glycine = profile(Chem.MolFromSmiles("NCC(=O)O"))
    assert Chem.MolToSmiles(glycine.major_species(7.4)) == "[NH3+]CC(=O)[O-]"
    assert glycine.charge(7.4) == pytest.approx(0.0, abs=0.05)
    assert glycine.fraction_neutral(7.4) < 0.01  # a zwitterion counts as ionised

    quaternary = profile(Chem.MolFromSmiles("C[N+](C)(C)C"))
    assert quaternary.charge(7.4) == pytest.approx(1.0)


def test_profile_runs_the_model_once_per_molecule(monkeypatch) -> None:
    calls = []
    real = molgpka.predict_sites
    monkeypatch.setattr(
        "chemsuite.ionisation.profile.predict_sites", lambda m: calls.append(1) or real(m)
    )
    mol = Chem.MolFromSmiles("CCN")
    profile(mol), profile(mol), clogd(mol)
    assert len(calls) == 1


def test_predict_pka_rows_blanks_and_failures() -> None:
    cells = ["CC(=O)O", None, "not a smiles", "c1ccccc1"]
    results = predict_pka(cells, 7.4)
    assert len(results) == 4
    assert results.most_acidic[0] == pytest.approx(4.47, abs=0.01)
    assert math.isnan(results.most_acidic[1]) and results.species_texts[1] == ""
    assert results.failed == 1  # the bad SMILES; the blank row is not a failure
    assert results.sites[3] == "" and results.neutral[3] == pytest.approx(100.0)
    assert results.species_texts[0] == "CC(=O)[O-]"
    assert results.sites[0] == "4.47 (acid)"


def test_predict_pka_keeps_finished_rows_when_cancelled(monkeypatch) -> None:
    monkeypatch.setattr("chemsuite.ionisation.predict.INLINE_CHUNK_ROWS", 2)
    token = CancelToken()
    seen = []

    def progress(done: int, total: int) -> None:
        seen.append(done)
        if done >= 2:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        predict_pka(["CCO", "CCN", "CC(=O)O", "c1ccncc1"], progress=progress, cancel=token)
    assert len(stopped.value.partial) == 2
