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

"""Building new molecules from a pool of BRICS fragments."""

from __future__ import annotations

import pytest
from rdkit import Chem

from chemsuite.analysis.recombine import (
    PRESETS,
    Limits,
    RecombineError,
    recombine,
    usable_fragments,
)
from chemsuite.core.jobs import Cancelled, CancelToken

#: paracetamol and benzene cut by BRICS, linkers kept
POOL = ["[16*]c1ccc(O)cc1", "[1*]C(C)=O", "[5*]N[5*]", "[16*]c1ccccc1"]


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_fragments_are_joined_only_where_brics_allows() -> None:
    made = recombine(POOL, max_products=50)
    # the amide is rebuilt through the N linker
    assert canonical("CC(=O)Nc1ccc(O)cc1") in made.products
    assert canonical("c1ccc(-c2ccccc2)cc1") in made.products
    # an acetyl cannot sit straight on an aromatic carbon: that is not a BRICS bond
    assert canonical("CC(=O)c1ccccc1") not in made.products


def test_products_are_distinct_and_sorted() -> None:
    made = recombine(POOL, max_products=50)
    assert made.products == sorted(set(made.products))


def test_a_pool_of_plain_attachment_points_is_refused() -> None:
    with pytest.raises(RecombineError, match="numbered attachment points"):
        recombine(["*Cl", "*C", "*OC"])


def test_a_pool_of_fewer_than_two_fragments_is_refused() -> None:
    with pytest.raises(RecombineError, match="at least two"):
        recombine(["[16*]c1ccccc1"])
    with pytest.raises(RecombineError, match="at least two"):
        recombine(["[16*]c1ccccc1", "[16*]c1ccccc1"])  # the same one twice


def test_usable_fragments_keeps_the_labelled_distinct_ones() -> None:
    assert usable_fragments(["[16*]c1ccccc1", "*Cl", "", None, "[16*]c1ccccc1"]) == [
        "[16*]c1ccccc1"
    ]


def test_the_product_cap_is_honoured() -> None:
    made = recombine(POOL, max_products=3)
    assert len(made.products) == 3


def test_depth_limits_how_many_fragments_go_into_a_product() -> None:
    shallow = recombine(POOL, max_depth=1, max_products=50)
    deep = recombine(POOL, max_depth=4, max_products=50)
    assert len(shallow.products) <= len(deep.products)


def test_property_limits_remove_what_falls_outside() -> None:
    made = recombine(POOL, max_products=50, limits=Limits({"MW": (None, 120.0)}))
    assert made.filtered
    for smiles in made.products:
        assert Chem.Descriptors.MolWt(Chem.MolFromSmiles(smiles)) <= 120.0


def test_a_preset_sets_the_limits() -> None:
    assert PRESETS["None"].ranges == {}
    assert not PRESETS["None"]
    assert PRESETS["Rule of 5"].ranges["MW"] == (None, 500.0)
    assert PRESETS["Lead-like"].ranges["MW"] == (200.0, 350.0)
    assert bool(PRESETS["Rule of 5"])


def test_limits_check_both_ends_of_a_range() -> None:
    small = Chem.MolFromSmiles("CC")
    big = Chem.MolFromSmiles("c1ccc(-c2ccccc2)cc1")
    limits = Limits({"MW": (100.0, 200.0)})
    assert not limits.passes(small)
    assert limits.passes(big)


def test_molecules_that_already_exist_are_dropped_and_counted() -> None:
    known = canonical("CC(=O)Nc1ccc(O)cc1")
    made = recombine(POOL, max_products=50, known=[known])
    assert made.known == 1
    assert known not in made.products


def test_progress_is_reported_per_product() -> None:
    seen: list[tuple[int, int]] = []
    recombine(POOL, max_products=3, progress=lambda done, total: seen.append((done, total)))
    assert seen[-1] == (3, 3)


def test_cancelling_keeps_the_products_already_made() -> None:
    token = CancelToken()

    def progress(done: int, _total: int) -> None:
        if done >= 2:
            token.cancel()

    with pytest.raises(Cancelled) as stopped:
        recombine(POOL, max_products=50, progress=progress, cancel=token)
    partial = stopped.value.partial
    assert partial is not None and partial.products


def test_a_filtered_run_gives_up_rather_than_building_for_ever() -> None:
    impossible = Limits({"MW": (None, 1.0)})  # nothing can pass
    made = recombine(POOL, max_products=50, limits=impossible)
    assert made.products == []
    assert made.filtered
