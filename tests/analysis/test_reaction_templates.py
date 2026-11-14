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

"""Every named template parses, and gives the hand-worked product for known building blocks."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.reaction_templates import (
    GROUPS,
    TEMPLATES,
    any_fits,
    by_id,
    fits,
    in_group,
    role_queries,
)
from chemsuite.analysis.reactions import apply_template, parse_template


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


def products(template_id: str, *blocks: str) -> list[str]:
    template = by_id(template_id)
    assert template is not None, template_id
    rxn = parse_template(template.smarts)
    mols = [Chem.MolFromSmiles(b) for b in blocks]
    assert all(m is not None for m in mols), blocks
    return apply_template(rxn, mols)


#: Each template, with building blocks and the product worked out by hand.
WORKED = [
    ("amide", ("CC(=O)O", "CNC"), "CC(=O)N(C)C"),
    ("amide", ("c1ccccc1C(=O)O", "NCc1ccccc1"), "O=C(NCc1ccccc1)c1ccccc1"),
    ("amide-acyl-chloride", ("CC(=O)Cl", "NCC"), "CCNC(C)=O"),
    ("sulfonamide", ("c1ccccc1S(=O)(=O)Cl", "NC1CC1"), "O=S(=O)(NC1CC1)c1ccccc1"),
    ("urea", ("CNC", "CCN=C=O"), "CCNC(=O)N(C)C"),
    ("carbamate", ("CCOC(=O)Cl", "NCC"), "CCNC(=O)OCC"),
    ("ester", ("CC(=O)O", "CCO"), "CCOC(C)=O"),
    ("suzuki", ("Brc1ccccc1", "OB(O)c1ccc(C)cc1"), "Cc1ccc(-c2ccccc2)cc1"),
    ("suzuki", ("Ic1ccncc1", "CC1(C)OB(c2ccccc2)OC1(C)C"), "c1ccc(-c2ccncc2)cc1"),
    ("buchwald", ("Brc1ccccc1", "C1CCNCC1"), "c1ccc(N2CCCCC2)cc1"),
    ("sonogashira", ("Brc1ccccc1", "C#Cc1ccccc1"), "C(#Cc1ccccc1)c1ccccc1"),
    ("chan-lam", ("OB(O)c1ccccc1", "C1CCNCC1"), "c1ccc(N2CCCCC2)cc1"),
    ("snar", ("Clc1ccncc1", "NCC"), "CCNc1ccncc1"),
    ("snar", ("Fc1ccccc1[N+](=O)[O-]", "C1CCNCC1"), "O=[N+]([O-])c1ccccc1N1CCCCC1"),
    ("reductive-amination", ("CCC=O", "C1CCNCC1"), "CCCN1CCCCC1"),
    ("reductive-amination", ("CC(=O)C", "NCc1ccccc1"), "CC(C)NCc1ccccc1"),
    ("n-alkylation", ("CCBr", "C1CCNCC1"), "CCN1CCCCC1"),
    ("ether", ("Oc1ccccc1", "CCBr"), "CCOc1ccccc1"),
    ("click", ("CCN=[N+]=[N-]", "C#Cc1ccccc1"), "CCn1cc(-c2ccccc2)nn1"),
    ("hantzsch-thiazole", ("CC(=O)CBr", "NC(=S)c1ccccc1"), "Cc1csc(-c2ccccc2)n1"),
    ("paal-knorr-pyrrole", ("CC(=O)CCC(C)=O", "NCc1ccccc1"), "Cc1ccc(C)n1Cc1ccccc1"),
    ("pyrazole", ("CC(=O)CC(C)=O", "NNc1ccccc1"), "Cc1cc(C)n(-c2ccccc2)n1"),
    ("benzimidazole", ("Nc1ccccc1N", "CC(=O)O"), "Cc1nc2ccccc2[nH]1"),
    ("oxadiazole", ("CC(=NO)N", "c1ccccc1C(=O)O"), "Cc1noc(-c2ccccc2)n1"),
    ("pictet-spengler", ("c1ccc2[nH]cc(CCN)c2c1", "CC=O"), "CC1NCCc2c1[nH]c1ccccc21"),
    ("tetrazole", ("CCC#N",), "CCc1nn[nH]n1"),
    (
        "ugi",
        ("CC(=O)O", "NCc1ccccc1", "CCC=O", "[C-]#[N+]C1CCCCC1"),
        "CCC(N(Cc1ccccc1)C(C)=O)C(=O)NC1CCCCC1",
    ),
    (
        "passerini",
        ("CC(=O)O", "CCC=O", "[C-]#[N+]C1CCCCC1"),
        "CCC(OC(C)=O)C(=O)NC1CCCCC1",
    ),
    (
        "biginelli",
        ("c1ccccc1C=O", "CC(=O)CC(=O)OCC", "NC(N)=O"),
        "CCOC(=O)C1=C(C)NC(=O)NC1c1ccccc1",
    ),
    (
        "gbb",
        ("Nc1ccccn1", "CCC=O", "[C-]#[N+]C1CCCCC1"),
        "CCc1nc2ccccn2c1NC1CCCCC1",
    ),
    ("van-leusen", ("CCC=O", "NCc1ccccc1"), "CCc1cncn1Cc1ccccc1"),
    ("boc-removal", ("CCNC(=O)OC(C)(C)C",), "CCN"),
    ("ester-hydrolysis", ("CCOC(C)=O",), "CC(=O)O"),
]


def test_every_template_parses() -> None:
    for template in TEMPLATES:
        rxn = parse_template(template.smarts)  # raises TemplateError if it cannot be applied
        assert rxn.GetNumReactantTemplates() == template.components, template.id


def test_ids_are_unique_and_groups_are_known() -> None:
    ids = [t.id for t in TEMPLATES]
    assert len(ids) == len(set(ids))
    assert {t.group for t in TEMPLATES} <= set(GROUPS)
    assert sum(len(in_group(g)) for g in GROUPS) == len(TEMPLATES)


def test_every_template_is_covered_by_a_worked_example() -> None:
    assert {t.id for t in TEMPLATES} == {case[0] for case in WORKED}


@pytest.mark.parametrize(("template_id", "blocks", "expected"), WORKED)
def test_worked_products(template_id: str, blocks: tuple[str, ...], expected: str) -> None:
    assert canonical(expected) in products(template_id, *blocks)


def test_amide_does_not_acylate_an_amide_nitrogen() -> None:
    assert products("amide", "CC(=O)O", "CC(=O)NC") == []


def test_reductive_amination_ignores_an_acid() -> None:
    assert products("reductive-amination", "CCC(=O)O", "C1CCNCC1") == []


def test_snar_needs_an_activated_ring() -> None:
    assert products("snar", "Clc1ccccc1", "NCC") == []


def test_stereocentres_away_from_the_reaction_survive() -> None:
    made = products("amide", "C[C@H](N(C)C)C(=O)O", "NCC")
    # the same arrangement around the centre as the acid, written from the amide end
    assert made == [canonical("C[C@H](N(C)C)C(=O)NCC")]


def test_fits_names_the_components_a_block_can_fill() -> None:
    acid = Chem.MolFromSmiles("CC(=O)O")
    amine = Chem.MolFromSmiles("NCC")
    assert fits("amide", acid) == [0]
    assert fits("amide", amine) == [1]
    # an amino acid can stand on either side
    assert fits("amide", Chem.MolFromSmiles("NCC(=O)O")) == [0, 1]


def test_any_fits_checks_a_whole_column() -> None:
    mols = [Chem.MolFromSmiles(s) for s in ("CCO", "CC(=O)O")]
    assert any_fits("amide", 0, mols)
    assert not any_fits("amide", 1, [Chem.MolFromSmiles("CCO")])
    assert not any_fits("amide", 5, mols)


def test_role_queries_are_cached_per_template() -> None:
    assert role_queries("amide") is role_queries("amide")
    assert len(role_queries("ugi")) == 4


def test_describe_names_the_roles_and_the_source() -> None:
    text = by_id("suzuki").describe()
    assert "Aryl halide or triflate + Boronic acid or ester" in text
    assert "Suzuki" in text or "Miyaura" in text
