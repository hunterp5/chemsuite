"""Ligand sites: forms from a microstate model (here :class:`Stub`, standing in for chemsuite
pKa), written back."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.ionisation.unipka import microstates
from chemsuite.macromolecule.prepare import (
    hydrogens,
    ligand_sites,
    pipeline,
    repair,
    sites,
    solve,
    templates,
)
from tests.macromolecule.conftest import load
from tests.macromolecule.prepare.test_hydrogens import _on

BENZAMIDINE_PH = 7.5
LN10 = float(np.log(10.0))
MEAN = 6.0
BASE_PKA, ACID_PKA = 10.0, 4.0
UNLIKELY = 50.0  # kT


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
    found, notes = ligand_sites.find(trypsin, BENZAMIDINE_PH, pka=Stub())
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
                            None, Stub())  # fmt: skip
    ligand = np.concatenate([x.moving for x in found if x.kind == sites.LIGAND])
    for site in found:
        if site.kind == sites.ROTOR:
            assert not np.isin(site.moving, ligand).any()


def test_a_phenol_near_its_pka_may_lose_its_proton() -> None:
    """T4 lysozyme's 2-propylphenol at pH 9, an acid of pKa 9.5: the phenol and the
    phenolate are its forms."""
    s = _protonated("3htb.cif.gz")
    found, _notes = ligand_sites.find(s, 9.0, pka=Stub(acid_pka=9.5))
    (site,) = found
    assert str(s.structure.res_name[site.residue]) == "JZ4"
    assert set(_forms(site)) == {ligand_sites.AS_DRAWN, "−H OAB"}  # OAB: its O
