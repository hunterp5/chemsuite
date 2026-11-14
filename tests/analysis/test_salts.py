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

"""Salt forms: what a strip removes, and making a salt of a parent."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.salts import COUNTERIONS, by_id, form_salt, removed_fragments


def salt(smiles: str, ion_id: str, count: int = 1) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    ion = by_id(ion_id)
    assert ion is not None, ion_id
    made = form_salt(mol, ion, count)
    return Chem.MolToSmiles(made) if made is not None else None


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_every_counterion_parses_and_carries_its_charge() -> None:
    for ion in COUNTERIONS:
        mol = Chem.MolFromSmiles(ion.smiles)
        assert mol is not None, ion.id
        assert Chem.GetFormalCharge(mol) == ion.charge, ion.id
        assert ion.acid == (ion.charge < 0), ion.id


def test_counterion_ids_are_unique() -> None:
    ids = [ion.id for ion in COUNTERIONS]
    assert len(ids) == len(set(ids))


def test_label_shows_the_abbreviation_when_there_is_one() -> None:
    assert by_id("hydrochloride").label == "Hydrochloride (HCl)"
    assert by_id("mesylate").label == "Mesylate"


@pytest.mark.parametrize(
    ("parent", "ion_id", "expected"),
    [
        ("CCN", "hydrochloride", "CC[NH3+].[Cl-]"),
        ("CCNCC", "mesylate", "CC[NH2+]CC.CS(=O)(=O)[O-]"),
        ("c1ccncc1", "hydrochloride", "c1cc[nH+]cc1.[Cl-]"),
        ("CC(=O)O", "sodium", "CC(=O)[O-].[Na+]"),
        ("Oc1ccccc1", "sodium", "[O-]c1ccccc1.[Na+]"),
        ("CCS(=O)(=O)O", "sodium", "CCS(=O)(=O)[O-].[Na+]"),
    ],
)
def test_simple_salts(parent: str, ion_id: str, expected: str) -> None:
    assert salt(parent, ion_id) == canonical(expected)


def test_a_doubly_charged_ion_repeats_the_parent() -> None:
    # ethylamine sulfate is the 2:1 salt
    assert salt("CCN", "sulfate") == canonical("CC[NH3+].CC[NH3+].[O-]S(=O)(=O)[O-]")
    assert salt("CC(=O)O", "calcium") == canonical("CC(=O)[O-].CC(=O)[O-].[Ca+2]")


def test_a_diamine_balances_a_sulfate_on_its_own() -> None:
    assert salt("NCCCCN", "sulfate") == canonical("[NH3+]CCCC[NH3+].[O-]S(=O)(=O)[O-]")


def test_two_counterions_give_the_di_salt() -> None:
    assert salt("NCCCCN", "hydrochloride", 2) == canonical("[NH3+]CCCC[NH3+].[Cl-].[Cl-]")


def test_the_strongest_site_is_chosen() -> None:
    # aspirin: the carboxylic acid, not the ester
    assert salt("CC(=O)Oc1ccccc1C(=O)O", "sodium") == canonical("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")
    # an amine wins over an aniline
    assert salt("NCCc1ccc(N)cc1", "hydrochloride") == canonical("[NH3+]CCc1ccc(N)cc1.[Cl-]")


def test_an_amide_nitrogen_is_not_basic() -> None:
    assert salt("CC(=O)NC", "hydrochloride") is None


def test_a_molecule_with_no_site_gives_none() -> None:
    assert salt("CCO", "hydrochloride") is None
    assert salt("CCO", "sodium") is None
    assert salt("c1ccccc1", "hydrochloride") is None


def test_charge_the_parent_already_carries_balances_the_ion() -> None:
    assert salt("C[NH3+]", "hydrochloride") == canonical("C[NH3+].[Cl-]")
    assert salt("CC(=O)[O-]", "sodium") == canonical("CC(=O)[O-].[Na+]")


def test_removed_fragments_names_what_a_strip_took() -> None:
    assert removed_fragments(Chem.MolFromSmiles("CCN.Cl"), Chem.MolFromSmiles("CCN")) == "Cl"
    assert removed_fragments(Chem.MolFromSmiles("CCN"), Chem.MolFromSmiles("CCN")) == ""


def test_removed_fragments_keeps_one_copy_of_a_repeated_fragment() -> None:
    removed = removed_fragments(
        Chem.MolFromSmiles("CC(=O)[O-].[Na+].[Na+]"), Chem.MolFromSmiles("CC(=O)[O-]")
    )
    assert removed == "[Na+].[Na+]"
