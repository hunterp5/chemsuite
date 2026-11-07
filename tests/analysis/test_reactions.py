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

"""Using reactions (REACTIONS_PLAN R4), checked against hand-worked reactions."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.reactions import (
    QueryError,
    ReactionQuery,
    TemplateError,
    apply_template,
    enumerate_library,
    enumerate_rows,
    parse_reaction_query,
    parse_template,
    search_reactions,
    split_cells,
)
from chemsuite.core.jobs import Cancelled, CancelToken

AMIDE = "[C:1](=[O:2])[OH].[N;H2,H1;!$(NC=O):3]>>[C:1](=[O:2])[N:3]"
SUZUKI = "[c:1][Br,I].[c:2][B]([OH])[OH]>>[c:1]-[c:2]"
REDUCTIVE_AMINATION = "[CH1:1]=O.[N;H2,H1;!$(NC=O):2]>>[CH2:1][N:2]"
BOC_REMOVAL = "[N:1]C(=O)OC(C)(C)C>>[N:1]"

# Products worked out by hand.
HAND_CHECKED = [
    (AMIDE, ["CC(=O)O", "NCc1ccccc1"], "CC(=O)NCc1ccccc1"),  # N-benzylacetamide
    (AMIDE, ["OC(=O)c1ccccc1", "C1CCNCC1"], "O=C(c1ccccc1)N1CCCCC1"),  # secondary amine
    (SUZUKI, ["Brc1ccc(C)cc1", "OB(O)c1ccccc1"], "Cc1ccc(-c2ccccc2)cc1"),  # 4-methylbiphenyl
    (SUZUKI, ["Ic1ccncc1", "OB(O)c1ccc(F)cc1"], "Fc1ccc(-c2ccncc2)cc1"),
    (REDUCTIVE_AMINATION, ["O=Cc1ccccc1", "C1CCNCC1"], "c1ccc(CN2CCCCC2)cc1"),  # 1-benzylpiperidine
    (BOC_REMOVAL, ["CC(C)(C)OC(=O)N1CCNCC1"], "C1CNCCN1"),  # piperazine
    (BOC_REMOVAL, ["CC(C)(C)OC(=O)N[C@@H](C)C(=O)O"], "C[C@H](N)C(=O)O"),  # keeps L-alanine
]


@pytest.mark.parametrize(("template", "reactants", "expected"), HAND_CHECKED)
def test_templates_give_the_hand_checked_product(template, reactants, expected) -> None:
    products = apply_template(parse_template(template), [Chem.MolFromSmiles(s) for s in reactants])
    assert products == [Chem.CanonSmiles(expected)]


def test_a_template_that_does_not_apply_gives_nothing() -> None:
    # an amide nitrogen is excluded by the template, so no product
    assert (
        apply_template(
            parse_template(AMIDE), [Chem.MolFromSmiles("CC(=O)O"), Chem.MolFromSmiles("CC(N)=O")]
        )
        == []
    )


def test_a_building_block_with_two_sites_gives_two_products() -> None:
    products = apply_template(
        parse_template(SUZUKI),
        [Chem.MolFromSmiles("Brc1ccc(I)cc1"), Chem.MolFromSmiles("OB(O)c1ccccc1")],
    )
    assert sorted(products) == sorted(
        [Chem.CanonSmiles("Brc1ccc(-c2ccccc2)cc1"), Chem.CanonSmiles("Ic1ccc(-c2ccccc2)cc1")]
    )


@pytest.mark.parametrize(
    ("text", "reason"),
    [
        ("", "Enter"),
        ("CCO", "reactants >> products"),
        ("CC>>CO", "Map the atoms"),
        ("C>>", "at least one reactant and one product"),
        ("[C:1>>C", "Not a valid"),
    ],
)
def test_bad_templates_say_why(text: str, reason: str) -> None:
    with pytest.raises(TemplateError, match=reason):
        parse_template(text)


def test_enumerate_row_by_row() -> None:
    products = enumerate_rows(
        AMIDE,
        [["CC(=O)O", "OC(=O)c1ccccc1", None, "CCO"], ["NC", "NCC", "NC", "NC"]],
    )
    assert products == [["CNC(C)=O"], ["CCNC(=O)c1ccccc1"], [], []]
    with pytest.raises(TemplateError, match="2 reactant"):
        enumerate_rows(AMIDE, [["CC(=O)O"]])


def test_enumerate_every_combination() -> None:
    library = enumerate_library(AMIDE, [["CC(=O)O", "OC(=O)c1ccccc1"], ["NC", "NCC", "CCO"]])
    assert library.combinations == 6 and library.failed == 2  # ethanol is no amine
    assert library.products == ["CNC(C)=O", "CCNC(C)=O", "CNC(=O)c1ccccc1", "CCNC(=O)c1ccccc1"]
    assert library.blocks == [[0, 0, 1, 1], [0, 1, 0, 1]]


def test_cancelled_enumeration_keeps_finished_combinations() -> None:
    token = CancelToken()
    seen = []

    def progress(done: int, total: int) -> None:
        seen.append(done)
        token.cancel()

    blocks = [["CC(=O)O"] * 150, ["NC", "NCC"]]
    with pytest.raises(Cancelled) as stopped:
        enumerate_library(AMIDE, blocks, progress=progress, cancel=token)
    partial = stopped.value.partial
    assert partial is not None and 0 < partial.combinations < 300
    assert len(partial.products) == partial.combinations


# -- search -------------------------------------------------------------------------------------

REACTIONS = [
    "CC(=O)O.NCc1ccccc1>CN(C)C=O>CC(=O)NCc1ccccc1",  # amide coupling
    "Brc1ccccc1.OB(O)c1ccccc1>[Pd]>c1ccc(-c2ccccc2)cc1",  # Suzuki
    "O=Cc1ccccc1.C1CCNCC1>>c1ccc(CN2CCCCC2)cc1",  # reductive amination
    "CC(C)(C)OC(=O)N1CCNCC1>>C1CNCCN1",  # Boc removal
    None,
    "not a reaction",
]


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        (ReactionQuery("C(=O)O.N>>C(=O)N"), [True, False, False, False, False, False]),
        (
            ReactionQuery("[CX3](=O)[OX2H1]>>[CX3](=O)[NX3]"),
            [True, False, False, False, False, False],
        ),
        (ReactionQuery("c1ccccc1Br>>"), [False, True, False, False, False, False]),
        (ReactionQuery(">[Pd]>"), [False, True, False, False, False, False]),
        # piperidine: in the reductive amination, not in piperazine (two nitrogens)
        (ReactionQuery("C1CCNCC1", "reactant"), [False, False, True, False, False, False]),
        (ReactionQuery("C1CCNCC1", "product"), [False, False, True, False, False, False]),
        (ReactionQuery("C1CNCCN1", "product"), [False, False, False, True, False, False]),
        (ReactionQuery("CN(C)C=O", "agent"), [True, False, False, False, False, False]),
        (ReactionQuery("[Pd]", "any"), [False, True, False, False, False, False]),
        (ReactionQuery("c1ccccc1-c1ccccc1", "product"), [False, True, False, False, False, False]),
    ],
)
def test_reaction_search(query: ReactionQuery, expected: list[bool]) -> None:
    assert search_reactions(query, REACTIONS).tolist() == expected


def test_query_mistakes_are_explained() -> None:
    with pytest.raises(QueryError, match="Enter"):
        parse_reaction_query(ReactionQuery(""))
    with pytest.raises(QueryError, match="Enter a reaction"):
        parse_reaction_query(ReactionQuery("CCO"))
    with pytest.raises(QueryError, match="Reaction contains"):
        parse_reaction_query(ReactionQuery("C>>O", "product"))


# -- split --------------------------------------------------------------------------------------


def test_split_by_role_drops_atom_maps() -> None:
    mapped = "[CH3:1][OH:2].Cl>O>[CH3:1][Cl:2]"
    assert split_cells([mapped, None, "bad>>"]) == [
        {"reactant": ["CO", "Cl"], "agent": ["O"], "product": ["CCl"]},
        None,
        None,
    ]
    kept = split_cells([mapped], keep_maps=True)[0]
    assert kept["reactant"][0] == "[CH3:1][OH:2]"
