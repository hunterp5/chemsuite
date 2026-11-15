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

"""Every reagent kind parses, and gives the hand-worked fragment for a known reagent."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.reactions import parse_template
from chemsuite.analysis.reagents import (
    GROUPS,
    REAGENTS,
    Prepared,
    any_fits,
    by_id,
    in_group,
    prepare,
    prepare_cells,
    sites,
)
from chemsuite.analysis.rgroup_library import build, core_positions, parse_core, parse_group
from chemsuite.core.jobs import Cancelled, CancelToken


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


def fragments(reagent_id: str, smiles: str, *, clean: bool = True) -> tuple[str, ...]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return prepare(reagent_id, mol, clean=clean).fragments


#: Each reagent kind, with a reagent and the fragment worked out by hand.
WORKED = [
    ("acyl-chloride", "CC(=O)Cl", "*C(C)=O"),
    ("carboxylic-acid", "CC(=O)O", "*C(C)=O"),
    ("sulfonyl-chloride", "CS(=O)(=O)Cl", "*S(C)(=O)=O"),
    ("chloroformate", "CCOC(=O)Cl", "*C(=O)OCC"),
    ("isocyanate", "CCN=C=O", "*C(=O)NCC"),
    ("aryl-halide", "Ic1ccccc1", "*c1ccccc1"),
    ("boronic", "OB(O)c1ccccc1", "*c1ccccc1"),
    ("boronic", "CC1(C)OB(c2ccccc2)OC1(C)C", "*c1ccccc1"),
    ("alkyl-halide", "CCBr", "*CC"),
    ("aldehyde", "CCC=O", "*CCC"),
    ("alkyne", "C#Cc1ccccc1", "*C#Cc1ccccc1"),
    ("amine", "C1CCNCC1", "*N1CCCCC1"),
    ("alcohol", "CCO", "*OCC"),
    ("thiol", "CCS", "*SCC"),
]


@pytest.mark.parametrize(("reagent_id", "reagent", "fragment"), WORKED)
def test_worked_fragments(reagent_id: str, reagent: str, fragment: str) -> None:
    assert fragments(reagent_id, reagent) == (canonical(fragment),)


def test_every_reagent_smarts_parses() -> None:
    """Every tier of every kind is a usable one-reactant transform."""
    for reagent in REAGENTS:
        assert reagent.smarts, reagent.id
        for text in reagent.smarts:
            rxn = parse_template(text)
            assert rxn.GetNumReactantTemplates() == 1, reagent.id
            assert rxn.GetNumProductTemplates() == 1, reagent.id


def test_every_reagent_has_a_worked_example() -> None:
    assert {r.id for r in REAGENTS} == {reagent_id for reagent_id, _r, _f in WORKED}


def test_every_reagent_is_in_a_listed_group_and_describes_itself() -> None:
    for reagent in REAGENTS:
        assert reagent.group in GROUPS, reagent.id
        assert reagent.attaches in reagent.describe()
    assert sum(len(in_group(group)) for group in GROUPS) == len(REAGENTS)
    assert by_id("acyl-chloride") is not None and by_id("nonsense") is None


def test_a_reagent_of_another_kind_gives_nothing() -> None:
    assert fragments("acyl-chloride", "c1ccccc1") == ()
    assert fragments("amine", "CC(=O)NC") == ()  # an amide N is not an amine
    assert fragments("aldehyde", "CC(C)=O") == ()  # a ketone is not an aldehyde
    missed = prepare("acyl-chloride", Chem.MolFromSmiles("c1ccccc1"))
    assert missed == Prepared(reagent="c1ccccc1")  # the reagent is kept, so it can be listed


def test_the_first_matching_tier_wins() -> None:
    """Iodide before bromide before chloride, so only one leaving group is taken."""
    assert sites("aryl-halide", Chem.MolFromSmiles("Ic1ccc(Cl)cc1")) == (0, 1)
    assert sites("aryl-halide", Chem.MolFromSmiles("Brc1ccc(Cl)cc1")) == (1, 1)
    assert sites("aryl-halide", Chem.MolFromSmiles("Clc1ccccc1")) == (3, 1)
    assert fragments("aryl-halide", "Ic1ccc(Cl)cc1") == (canonical("*c1ccc(Cl)cc1"),)


def test_every_site_is_prepared_and_identical_fragments_collapse() -> None:
    """A reagent with two of the same group gives one fragment per site; a symmetric one
    gives a single fragment, and the site count is how it is spotted."""
    diacid = prepare("carboxylic-acid", Chem.MolFromSmiles("OC(=O)c1ccc(C(=O)O)cc1"))
    assert diacid.sites == 2
    assert diacid.fragments == (canonical("*C(=O)c1ccc(C(=O)O)cc1"),)  # symmetric: collapsed
    mixed = prepare("carboxylic-acid", Chem.MolFromSmiles("OC(=O)CCCC(=O)Nc1ccccc1"))
    assert mixed.sites == 1
    unsymmetric = prepare("alcohol", Chem.MolFromSmiles("OCCCc1ccc(O)cc1"))
    assert unsymmetric.sites == 2
    assert len(unsymmetric.fragments) == 2  # the alkyl and the phenol ends differ
    assert all(f.count("*") == 1 for f in unsymmetric.fragments)


def test_a_salt_is_cleaned_up_first() -> None:
    assert fragments("amine", "CC[NH3+].[Cl-]") == (canonical("*NCC"),)
    assert fragments("amine", "CC[NH3+].[Cl-]", clean=False) == ()
    assert any_fits("amine", [Chem.MolFromSmiles("CC[NH3+].[Cl-]")])
    assert not any_fits("amine", [Chem.MolFromSmiles("CC[NH3+].[Cl-]")], clean=False)


def test_stereo_away_from_the_site_is_kept() -> None:
    assert fragments("acyl-chloride", "C[C@H](c1ccccc1)C(=O)Cl") == (
        canonical("*C(=O)[C@H](C)c1ccccc1"),
    )


def test_a_prepared_fragment_joins_a_core() -> None:
    """The whole point: an acyl chloride and an N-H core give the amide."""
    (acyl,) = fragments("acyl-chloride", "CC(=O)Cl")
    core_text = "C1CCN([*:1])CC1"
    group = parse_group(acyl, 1)
    assert group is not None
    product = build(parse_core(core_text), [group], core_positions(core_text))
    assert product == canonical("CC(=O)N1CCCCC1")


def test_prepare_cells_reports_rows_and_counts() -> None:
    rows = prepare_cells("acyl-chloride", ["CC(=O)Cl", "c1ccccc1", None, "CCC(=O)Cl"])
    assert len(rows) == 4
    assert [r.fragment for r in rows] == [canonical("*C(C)=O"), "", "", canonical("*C(=O)CC")]
    assert rows.fragments == 2
    assert rows.missing == 2  # the benzene and the blank


def test_prepare_cells_keeps_finished_rows_on_cancel() -> None:
    token = CancelToken()
    cells = ["CC(=O)Cl"] * 60
    seen: list[int] = []

    def progress(done: int, _total: int) -> None:
        seen.append(done)
        token.cancel()

    with pytest.raises(Cancelled) as caught:
        prepare_cells("acyl-chloride", cells, progress=progress, cancel=token)
    partial = caught.value.partial
    assert partial is not None and 0 < len(partial) < len(cells)
    assert all(r.fragment for r in partial.rows)


def test_an_unknown_reagent_kind_is_refused() -> None:
    with pytest.raises(ValueError, match="unknown reagent kind"):
        prepare_cells("nonsense", ["CC(=O)Cl"])
