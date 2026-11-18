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

"""A molecule's pKa profile from microstate free energies: here a stub whose sites are
independent with textbook pKa values (:class:`SmilesStub`); with the shipped weights, a few
measured values too."""

from __future__ import annotations

import math

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.ionisation import profile
from chemsuite.ionisation.profile import ACID, BASE
from chemsuite.ionisation.unipka import ranges

LN10 = float(np.log(10.0))
MEAN = 6.0


class SmilesStub:
    """Independent sites: an amine's N⁺–H a base of pKa 10, an aromatic N⁺–H one of 5, a
    carboxylate O⁻ an acid of 4, any other O⁻ one of 16; any other charge costs 50 kT. A
    proton on a site of pKa p costs ``−ln10·(p − μ)``, so each crosses over at its pKa."""

    mean = MEAN

    def __init__(self) -> None:
        self.asked: list[str] = []

    def free_energy(self, smiles: str) -> float:
        self.asked.append(smiles)
        mol = Chem.MolFromSmiles(smiles)
        g = 0.0
        for a in mol.GetAtoms():
            q = a.GetFormalCharge()
            if q > 0 and a.GetAtomicNum() == 7 and a.GetTotalNumHs():
                g -= LN10 * ((5.0 if a.GetIsAromatic() else 10.0) - MEAN)
            elif q < 0 and a.GetAtomicNum() == 8:
                carboxyl = any(n.GetAtomicNum() == 6 and _has_carbonyl(n) for n in a.GetNeighbors())
                g += LN10 * ((4.0 if carboxyl else 16.0) - MEAN)
            elif q != 0 and not (a.GetAtomicNum() == 7 and a.GetDegree() == 4):
                g += 50.0
        return g


def _has_carbonyl(carbon: Chem.Atom) -> bool:
    double = Chem.BondType.DOUBLE
    return any(b.GetBondType() == double and b.GetOtherAtom(carbon).GetAtomicNum() == 8
               for b in carbon.GetBonds())  # fmt: skip


@pytest.fixture
def stub():
    model = SmilesStub()
    profile.use(model)
    yield model
    profile.use(None)


def _profile(smiles: str, ph: float = 7.4) -> profile.Profile:
    found = profile.profile(Chem.MolFromSmiles(smiles), ph)
    assert found is not None
    return found


def test_acetic_acid_is_an_acid_of_its_pka(stub) -> None:
    found = _profile("CC(=O)O")
    assert [(round(s.pka, 2), s.kind) for s in found.steps] == [(4.0, ACID)]
    assert found.most_acidic == pytest.approx(4.0) and found.most_basic is None
    assert found.charge(4.0) == pytest.approx(-0.5)
    assert Chem.MolToSmiles(found.major_species(7.4)) == "CC(=O)[O-]"


def test_glycine_is_marked_by_the_atom_that_changes(stub) -> None:
    """The +1 → 0 step loses the carboxyl's proton (an acid), 0 → −1 the ammonium's (a base,
    though the molecule is then a zwitterion)."""
    found = _profile("NCC(=O)O")
    assert [(round(s.pka, 1), s.kind) for s in found.steps] == [(4.0, ACID), (10.0, BASE)]
    assert found.charge(7.0) == pytest.approx(0.0, abs=1e-3)
    assert found.fraction_neutral(7.0) < 1e-3  # the zwitterion counts as ionized
    assert Chem.MolToSmiles(found.major_species(7.0)) == "[NH3+]CC(=O)[O-]"


def test_a_permanent_charge_is_not_an_ionized_site(stub) -> None:
    found = _profile("C[N+](C)(C)CC(=O)O")
    assert found.charge(2.0) == pytest.approx(1.0, abs=0.01)
    assert found.fraction_neutral(2.0) == pytest.approx(1.0, abs=0.01)
    assert found.logd(1.5, 2.0) == pytest.approx(1.5, abs=0.01)


def test_the_range_widens_to_the_charges_populated_at_the_ph(stub) -> None:
    """Propane-1,2,3-tricarboxylic acid: the trianion at pH 7 lies past the first range."""
    found = _profile("OC(=O)CC(C(=O)O)CC(=O)O", 7.0)
    assert found.charge(7.0) == pytest.approx(-3.0, abs=0.01)
    assert len([s for s in found.steps if s.kind == ACID]) == 3
    low = _profile("OC(=O)CC(C(=O)O)CC(=O)O", 1.0)
    assert min(m.charge for m in low.ensemble.states) == -1  # not needed at pH 1


def test_oxonium_microstates_are_left_out(stub) -> None:
    found = _profile("CC(=O)O")
    assert not any("[OH+]" in m.smiles or "[OH2+]" in m.smiles for m in found.ensemble.states)


def test_too_many_microstates_give_no_profile(stub, monkeypatch) -> None:
    monkeypatch.setattr(ranges, "MAX_MICROSTATES", 1)
    assert profile.profile(Chem.MolFromSmiles("NCC(=O)O")) is None


def test_a_salt_gives_its_parents_values(stub) -> None:
    salt = _profile("CC(=O)[O-].[Na+]")
    assert Chem.MolToSmiles(salt.mol) == "CC(=O)O"


def test_without_another_model_the_shipped_weights_run() -> None:
    from chemsuite.ionisation.unipka import weights

    profile.use(None)
    assert profile._current() is weights.predictor()


def test_clogd_corrects_for_the_ionized_share(stub) -> None:
    mol = Chem.MolFromSmiles("CC(=O)O")
    logp = profile.clogd(mol, 2.0)
    assert profile.clogd(mol, 6.0) == pytest.approx(logp - 2.0, abs=0.01)


# -- with the shipped weights ---------------------------------------------------------------------
@pytest.mark.parametrize(
    ("smiles", "kind", "measured"),
    [
        ("CC(=O)O", ACID, 4.76),
        ("CCN", BASE, 10.7),
        ("Oc1ccccc1", ACID, 10.0),
        ("c1ccncc1", BASE, 5.2),
    ],  # fmt: skip
)
def test_the_model_is_near_measured_values(smiles, kind, measured) -> None:
    profile.use(None)
    found = _profile(smiles)
    value = found.most_acidic if kind == ACID else found.most_basic
    assert value is not None and math.isclose(value, measured, abs_tol=1.0)
