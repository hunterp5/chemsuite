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

"""Ligand sites: protomers and tautomers from MolGpKa and RDKit, written back."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemsuite.ionisation import molgpka
from chemsuite.protein.prepare import (
    energy,
    hydrogens,
    ligand_sites,
    pipeline,
    populations,
    repair,
    sites,
    solve,
    templates,
)
from chemsuite.protein.prepare.energy import Tables, Weights
from tests.protein.conftest import load
from tests.protein.prepare.test_hydrogens import _on

BENZAMIDINE_PH = 7.5


def _protonated(name: str) -> hydrogens.Protonated:
    s = load(name)
    s = s.subset(s.element != 1)  # neutron entries: their H and D atoms off
    library = templates.Library()
    return hydrogens.add_hydrogens(repair.build_missing(repair.one_location(s), library), library)


@pytest.fixture(scope="module")
def trypsin() -> hydrogens.Protonated:
    """Trypsin with benzamidine (5MO0), cut to the pocket."""
    return _protonated("5mo0_pocket.cif.gz")


@pytest.fixture(scope="module")
def benzamidine(trypsin) -> sites.Site:
    found, notes = ligand_sites.find(trypsin, BENZAMIDINE_PH)
    assert not notes
    (site,) = [x for x in found if str(trypsin.structure.res_name[x.residue]) == "BEN"]
    return site


def _forms(site: sites.Site) -> dict[str, list[sites.State]]:
    out: dict[str, list[sites.State]] = {}
    for st in site.states:
        out.setdefault(st.form, []).append(st)
    return out


def test_benzamidine_is_drawn_neutral_and_may_take_a_proton(benzamidine) -> None:
    forms = _forms(benzamidine)
    assert set(forms) == {ligand_sites.AS_DRAWN, "+H N1"}
    assert benzamidine.titrating and benzamidine.kind == sites.LIGAND
    charged = forms["+H N1"][0]
    assert charged.protons - forms[ligand_sites.AS_DRAWN][0].protons == 1
    assert charged.unlike == 1 and forms[ligand_sites.AS_DRAWN][0].unlike == 0
    assert sum(1 for p in charged.parents if p == "N1") == 2
    formal = dict(zip(charged.names, charged.formal, strict=True))
    partial = dict(zip(charged.names, charged.partial, strict=True))
    assert formal["N1"] == 1 and formal["N2"] == 0
    assert partial["N1"] == pytest.approx(0.5) and partial["N2"] == pytest.approx(0.5)


def test_alone_in_water_a_ligand_titrates_at_its_molgpka_pka(trypsin, benzamidine) -> None:
    """With nothing around it and no prior for its drawing, the amidine's protonated share
    crosses ½ at MolGpKa's pKa."""
    single = energy.priors(benzamidine, BENZAMIDINE_PH)
    table = Tables([benzamidine], [single], {}, BENZAMIDINE_PH, [energy.protons(benzamidine)])
    found = populations.titration(table, Weights(ccd=0.0, tautomer=0.0)).pka[0]
    mol = Chem.MolFromSmiles("NC(=N)c1ccccc1")
    expected = max(s.pka for s in molgpka.predict_sites(mol) if s.kind == molgpka.BASE)
    assert found == pytest.approx(expected, abs=0.05)


def test_apply_writes_the_ligand_form(trypsin, benzamidine) -> None:
    k = next(i for i, st in enumerate(benzamidine.states) if st.form == "+H N1")
    done = solve.apply(trypsin, [benzamidine], [k])
    s, r = done.structure, benzamidine.residue
    assert len(_on(s, r, "N1")) == 2 and len(_on(s, r, "N2")) == 2
    span = s.atoms_of(r)
    n1 = next(i for i in range(span.start, span.stop) if s.name[i] == "N1")
    assert int(s.charge[n1]) == 1
    assert done.states[r] == "+H N1"


def test_rotors_on_a_ligand_sites_atoms_are_left_to_it(trypsin) -> None:
    found = pipeline._sites(trypsin, pipeline.Options(ph=BENZAMIDINE_PH), None, frozenset(), [],
                            None)  # fmt: skip
    ligand = np.concatenate([x.moving for x in found if x.kind == sites.LIGAND])
    for site in found:
        if site.kind == sites.ROTOR:
            assert not np.isin(site.moving, ligand).any()


def test_bases_whose_lone_pair_is_not_free_are_left_out() -> None:
    """9-Methyladenine: N6 (an amine conjugated with ring bases) and N9 (a pyrrole N) take no
    proton; the ring N may."""
    adenine = Chem.MolFromSmiles("Cn1cnc2c(N)ncnc21")  # N9 is atom 1, N6 atom 6
    found = molgpka.predict_on(adenine)
    assert {1, 6} <= {s.atom for s in found if s.kind == molgpka.BASE}
    bases = {s.atom for s in ligand_sites._plausible(adenine, found) if s.kind == molgpka.BASE}
    assert bases and not bases & {1, 6}


def test_alcohols_are_not_acids() -> None:
    """NADP's ribose hydroxyls: MolGpKa calls them acids (pKa 2–7); only the phosphates are."""
    nadp = Chem.MolFromSmiles(
        "NC(=O)c1ccc[n+]([C@@H]2O[C@H](CO[P@@](=O)([O-])O[P@](=O)(O)OC[C@H]3O[C@@H]"
        "(n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1"
    )
    found = molgpka.predict_on(nadp)
    on_carbon = {s.atom for s in found if s.kind == molgpka.ACID
                 and nadp.GetAtomWithIdx(s.atom).GetNeighbors()[0].GetSymbol() == "C"}  # fmt: skip
    assert on_carbon  # the model's mistake, which the filter is for
    acids = [
        s.atom
        for s in ligand_sites._plausible(nadp, found)
        if s.kind == molgpka.ACID and nadp.GetAtomWithIdx(s.atom).GetSymbol() == "O"
    ]
    assert acids and all(nadp.GetAtomWithIdx(k).GetNeighbors()[0].GetSymbol() == "P"
                         for k in acids)  # fmt: skip


def test_like_charges_on_one_ligand_cost() -> None:
    """Methyl phosphate: the dianion pays for its two charges, the monoanion does not."""
    mol = Chem.AddHs(Chem.MolFromSmiles("COP(=O)(O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    xyz = mol.GetConformer().GetPositions()
    neutral = Chem.RemoveHs(mol)
    pose = xyz[: neutral.GetNumAtoms()]
    forms, _fixed, _drawn = ligand_sites._protomers(neutral, set(), 7.0, pose)
    by_charge = {sum(f.charges): f for f in forms}
    assert by_charge[-2].within > 2.0
    assert by_charge[-1].within == 0.0


def test_a_phenol_near_its_pka_may_lose_its_proton() -> None:
    """T4 lysozyme's 2-propylphenol at pH 9: the phenol and the phenolate are its forms."""
    s = _protonated("3htb.cif.gz")
    found, _notes = ligand_sites.find(s, 9.0)
    (site,) = found
    assert str(s.structure.res_name[site.residue]) == "JZ4"
    assert set(_forms(site)) == {ligand_sites.AS_DRAWN, "−H OAB"}  # OAB: its O
