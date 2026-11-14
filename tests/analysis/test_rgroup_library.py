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

"""A library from a core and a list of substituents per position."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.rgroup_library import (
    CoreError,
    core_positions,
    enumerate_rgroups,
    parse_core,
    parse_group,
)
from chemsuite.core.jobs import Cancelled, CancelToken

CORE = "c1cc([*:2])ccc1[*:1]"
ONE = "c1ccc([*:1])cc1"


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_core_positions_are_found_in_order() -> None:
    assert core_positions(CORE) == (1, 2)
    assert core_positions(ONE) == (1,)
    assert core_positions("c1ccccc1") == ()


def test_a_core_without_attachment_points_is_refused() -> None:
    with pytest.raises(CoreError, match="no numbered attachment point"):
        parse_core("c1ccccc1")


def test_a_blank_or_invalid_core_is_refused() -> None:
    with pytest.raises(CoreError, match="Enter a core"):
        parse_core("   ")
    with pytest.raises(CoreError, match="Not a valid core"):
        parse_core("c1ccc")


def test_a_substituent_is_numbered_for_its_position() -> None:
    group = parse_group("*Cl", 2)
    assert group is not None
    dummy = next(a for a in group.GetAtoms() if a.GetAtomicNum() == 0)
    assert dummy.GetAtomMapNum() == 2


def test_hydrogen_means_no_substituent() -> None:
    assert parse_group("H", 1) is None
    assert parse_group("h", 1) is None
    assert parse_group("", 1) is None


def test_a_substituent_with_two_attachments_is_refused() -> None:
    assert parse_group("*CC*", 1) is None


def test_every_combination_is_built() -> None:
    library = enumerate_rgroups(CORE, [["*Cl", "*C"], ["*OC", "H"]])
    assert library.combinations == 4
    assert library.failed == 0
    assert sorted(library.products) == sorted(
        canonical(s) for s in ("COc1ccc(Cl)cc1", "Clc1ccccc1", "COc1ccc(C)cc1", "Cc1ccccc1")
    )


def test_hydrogen_leaves_the_position_unsubstituted() -> None:
    library = enumerate_rgroups(ONE, [["H"]])
    assert library.products == [canonical("c1ccccc1")]


def test_each_product_records_the_substituents_it_came_from() -> None:
    library = enumerate_rgroups(CORE, [["*Cl", "*C"], ["*OC"]])
    where = zip(library.blocks[0], library.blocks[1], strict=True)
    pairs = dict(zip(library.products, where, strict=True))
    assert pairs[canonical("COc1ccc(Cl)cc1")] == (0, 0)
    assert pairs[canonical("COc1ccc(C)cc1")] == (1, 0)


def test_a_substituent_that_cannot_be_attached_is_counted_as_a_failure() -> None:
    library = enumerate_rgroups(ONE, [["*Cl", "not a structure"]])
    assert library.combinations == 2
    assert library.failed == 1
    assert library.products == [canonical("Clc1ccccc1")]


def test_the_wrong_number_of_substituent_lists_is_refused() -> None:
    with pytest.raises(CoreError, match="2 attachment point"):
        enumerate_rgroups(CORE, [["*Cl"]])


def test_a_position_with_no_substituent_is_refused() -> None:
    with pytest.raises(CoreError, match="at least one substituent"):
        enumerate_rgroups(CORE, [["*Cl"], []])


def test_progress_is_reported() -> None:
    seen: list[tuple[int, int]] = []
    enumerate_rgroups(
        CORE,
        [[f"*C{'C' * i}" for i in range(20)], ["*OC", "H"]],
        progress=lambda done, total: seen.append((done, total)),
    )
    assert seen and seen[-1] == (40, 40)


def test_cancelling_keeps_the_products_already_made() -> None:
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done >= 200:
            token.cancel()

    groups = [[f"*C{'C' * i}" for i in range(40)], [f"*OC{'C' * i}" for i in range(40)]]
    with pytest.raises(Cancelled) as stopped:
        enumerate_rgroups(CORE, groups, progress=progress, cancel=token)
    partial = stopped.value.partial
    assert partial is not None and 0 < partial.combinations < 1600
    assert partial.products


def test_a_decomposition_feeds_straight_back_in() -> None:
    """The core and labels Central Ring writes are the ones this takes."""
    from chemsuite.sar.central_ring import decompose_central_ring

    split = decompose_central_ring(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    assert split is not None
    library = enumerate_rgroups(split.core, [[g] for g in split.groups])
    assert library.products == [canonical("CC(=O)Nc1ccc(O)cc1")]


def test_an_unreadable_substituent_is_not_mistaken_for_hydrogen() -> None:
    """A substituent that cannot be parsed must fail, not quietly leave the position bare."""
    library = enumerate_rgroups(ONE, [["not a structure"]])
    assert library.products == []
    assert library.failed == 1
