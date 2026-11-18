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

"""Ligand forms from a microstate model (chemsuite pKa), here a stub whose free energies make
every base a pKa of :data:`BASE_PKA` and every acid one of :data:`ACID_PKA`; with the real
weights (``samples/unipka/weights/e20_f0.npz``) the named cases run too."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.ionisation.unipka import microstates
from chemsuite.macromolecule.prepare import energy, ligand_sites, pipeline, populations, sites
from chemsuite.macromolecule.prepare.energy import Tables, Weights
from chemsuite.macromolecule.prepare.ligand_pka import ModelStates
from chemsuite.macromolecule.prepare.sites import KT, RT_LN10
from tests.macromolecule.conftest import load
from tests.macromolecule.prepare.test_ligand_sites import BENZAMIDINE_PH, _forms, _protonated

LN10 = float(np.log(10.0))
MEAN = 6.0
BASE_PKA, ACID_PKA = 10.0, 4.0
UNLIKELY = 50.0  # kT
WEIGHTS = Path(__file__).parents[3] / "samples" / "unipka" / "weights" / "e20_f0.npz"


class Stub:
    """Each N⁺ a base of pKa :data:`BASE_PKA`, each O⁻ or S⁻ an acid of :data:`ACID_PKA`
    (any other charge :data:`UNLIKELY`), independent of one another: ``E(pH) = G + ln10·n·
    (pH − μ)`` crosses over at the pKa when one proton costs ``−ln10·(pKa − μ)``."""

    name = "stub"
    mean = MEAN

    def __init__(self, base_pka: float = BASE_PKA, acid_pka: float = ACID_PKA) -> None:
        self.base_pka, self.acid_pka = base_pka, acid_pka
        self.asked: list[str] = []

    def microstates(self, mol: Chem.Mol, lower: int, upper: int):
        return microstates.ensemble(mol, lower, upper)

    def free_energy(self, mol: Chem.Mol) -> float:
        self.asked.append(Chem.MolToSmiles(mol))
        g = 0.0
        for a in mol.GetAtoms():
            q, n = a.GetFormalCharge(), a.GetAtomicNum()
            if q > 0:
                g += -LN10 * (self.base_pka - MEAN) if n == 7 else UNLIKELY
            elif q < 0:
                g += LN10 * (self.acid_pka - MEAN) if n in (8, 16) else UNLIKELY
        return g


def _at(form: ligand_sites.Form, ph: float) -> float:
    return form.intrinsic + RT_LN10 * form.protons * ph


@pytest.fixture(scope="module")
def trypsin():
    return _protonated("5mo0_pocket.cif.gz")


@pytest.fixture(scope="module")
def benzamidine(trypsin) -> sites.Site:
    found, notes = ligand_sites.find(trypsin, BENZAMIDINE_PH, pka=Stub())
    assert not notes
    (site,) = found
    return site


def test_benzamidine_takes_the_models_forms(benzamidine) -> None:
    forms = _forms(benzamidine)
    assert set(forms) == {ligand_sites.AS_DRAWN, "+H N1"}
    for st in benzamidine.states:
        assert st.within == 0.0 and st.attract == 0.0 and st.tautomer_cost == 0.0


def test_alone_in_water_a_ligand_titrates_at_the_models_pka(benzamidine) -> None:
    single = energy.priors(benzamidine, BENZAMIDINE_PH)
    table = Tables([benzamidine], [single], {}, BENZAMIDINE_PH, [energy.protons(benzamidine)])
    found = populations.titration(table, Weights(ccd=0.0, tautomer=0.0)).pka[0]
    assert found == pytest.approx(BASE_PKA, abs=0.05)


def test_a_forms_free_energy_is_the_models() -> None:
    """``intrinsic = kT·G − RT ln10·n·μ``, n the protons added to the neutral parent."""
    neutral = Chem.MolFromSmiles("NC(=N)c1ccccc1")
    stub = Stub()
    forms = ligand_sites._model_forms(neutral, set(), 7.0, neutral, stub)
    charged = next(f for f in forms if sum(f.charges) == 1)
    assert charged.protons == 1
    expected = KT * stub.free_energy(charged.mol) - RT_LN10 * MEAN
    assert charged.intrinsic == pytest.approx(expected)
    assert forms[0].intrinsic == pytest.approx(0.0)  # the drawing (neutral) comes first


def test_the_charge_range_widens_while_the_best_form_is_at_its_end() -> None:
    """Propane-1,2,3-tricarboxylic acid at pH 7: three acids of pKa 4 make the trianion
    best, which the first range (charges −1 to +1) does not reach."""
    acid = Chem.MolFromSmiles("OC(=O)CC(C(=O)O)CC(=O)O")
    forms = ligand_sites._model_forms(acid, set(), 7.0, acid, Stub())
    best = min(forms, key=lambda f: _at(f, 7.0))
    assert sum(best.charges) == -3


def test_frozen_atoms_keep_their_drawing() -> None:
    """Glycine bonded through its N: no form changes the N."""
    glycine = Chem.MolFromSmiles("NCC(=O)O")
    forms = ligand_sites._model_forms(glycine, {0}, 7.0, glycine, Stub())
    assert len(forms) > 1
    assert all(f.hydrogens[0] == 2 and f.charges[0] == 0 for f in forms)


def test_tautomers_are_scored_by_the_model() -> None:
    """2-Hydroxypyridine, neutral at pH 7 (no acid or base near it): its pyridone tautomer is
    a form with the model's own free energy."""
    mol = Chem.MolFromSmiles("Oc1ccccn1")
    stub = Stub(base_pka=-10.0, acid_pka=20.0)
    forms = ligand_sites._model_forms(mol, set(), 7.0, mol, stub)
    pyridone = [f for f in forms if f.hydrogens[0] == 0 and sum(f.charges) == 0]
    assert pyridone
    assert pyridone[0].tautomer_cost == 0.0
    assert Chem.MolToSmiles(pyridone[0].mol) in stub.asked


def test_too_many_microstates_fall_back_to_molgpka(trypsin, monkeypatch) -> None:
    monkeypatch.setattr(ligand_sites, "MAX_MICROSTATES", 1)
    found, notes = ligand_sites.find(trypsin, BENZAMIDINE_PH, pka=Stub())
    assert [n.text for n in notes] == ["Ligand states by MolGpKa: 4 microstates"]
    (site,) = found
    molgpka_site = ligand_sites.find(trypsin, BENZAMIDINE_PH)[0][0]
    assert [st.intrinsic for st in site.states] == [st.intrinsic for st in molgpka_site.states]


def _prepared(name: str, ph: float, pka) -> pipeline.Prepared:
    s = load(name)
    s = s.subset(s.element != 1)  # neutron entries: their H and D atoms off
    return pipeline.prepare(s, pipeline.Options(ph=ph, waters=pipeline.NONE), ligand_model=pka)


def test_the_report_names_the_model() -> None:
    for model, name in ((Stub(), "stub"), (None, "MolGpKa")):
        done = _prepared("5mo0_pocket.cif.gz", BENZAMIDINE_PH, model)
        lines = [f.text for f in done.report if f.label == "Ligands"]
        assert lines == [f"Protonation states by {name}"]


# -- with the real weights ------------------------------------------------------------------------
needs_weights = pytest.mark.skipif(not WEIGHTS.is_file(), reason="no chemsuite pKa weights")


@pytest.fixture(scope="module")
def model() -> ModelStates:
    from chemsuite.ionisation.unipka import weights

    return ModelStates(weights.predictor(WEIGHTS))


@needs_weights
def test_benzamidine_in_trypsin_is_the_amidinium(model) -> None:
    done = _prepared("5mo0_pocket.cif.gz", BENZAMIDINE_PH, model)
    s = done.structure
    (r,) = [r for r in range(s.n_residues) if str(s.res_name[r]) == "BEN"]
    assert done.states[r] == "+H N1"


@needs_weights
def test_acetazolamide_binds_zinc_as_the_sulfonamide_anion(model) -> None:
    """Carbonic anhydrase II with acetazolamide (4G0C, neutron): N1 keeps one D, on the Zn."""
    done = _prepared("4g0c_pocket.cif.gz", 7.8, model)
    s = done.structure
    (r,) = [r for r in range(s.n_residues) if str(s.res_name[r]) == "AZM"]
    assert done.states[r] == "−H N1"


@needs_weights
@pytest.mark.parametrize(
    ("smiles", "charge"),
    [
        ("OC(=O)CC(O)(C(=O)O)CC(=O)O", -3),  # citrate
        ("CC(=O)O", -1),  # acetic acid
        ("NCCCC[C@H](N)C(=O)O", 1),  # lysine
    ],
)
def test_the_model_gives_the_textbook_charge_at_ph_7(model, smiles, charge) -> None:
    mol = Chem.MolFromSmiles(smiles)
    forms = ligand_sites._model_forms(mol, set(), 7.0, mol, model)
    best = min(forms, key=lambda f: _at(f, 7.0))
    assert sum(best.charges) == charge
