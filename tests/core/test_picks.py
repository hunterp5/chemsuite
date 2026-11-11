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
"""Selecting rows by one column's values: empty, first occurrence, duplicates."""

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.core.table import NUMBER, Table, picks
from chemsuite.core.table.columns import STRUCTURE


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "ID": ["A", "B", "C", "D", "E", "F"],
            "Series": ["S1", "", "S2", "S1", "S2", "S1"],
            "pIC50": ["5", "6", "", "5", "7.0", "7"],
            "SMILES": ["CCO", "OCC", "c1ccccc1", "", "C1=CC=CC=C1", "not a molecule"],
        },
        undoable=False,
    )
    t.convert_column("pIC50", NUMBER)
    t.convert_column("SMILES", STRUCTURE)
    return t


def ids(table: Table, oids: np.ndarray) -> list[str]:
    values = table.column("ID").values(table.slots_of(oids))
    return values.tolist()


def every(table: Table) -> np.ndarray:
    return table.visible_oids()


def test_empty_and_not_empty(table) -> None:
    assert ids(table, picks.empty(table, "Series", every(table))) == ["B"]
    assert ids(table, picks.not_empty(table, "pIC50", every(table))) == ["A", "B", "D", "E", "F"]


def test_first_occurrences_and_duplicates_split_the_values(table) -> None:
    first = picks.first_occurrences(table, "Series", every(table))
    later = picks.duplicates(table, "Series", every(table))
    assert ids(table, first) == ["A", "C"]
    assert ids(table, later) == ["D", "E", "F"]
    repeated = picks.all_duplicated(table, "Series", every(table))
    assert ids(table, repeated) == ["A", "C", "D", "E", "F"]  # B is empty


def test_numbers_compare_by_value(table) -> None:
    assert ids(table, picks.duplicates(table, "pIC50", every(table))) == ["D", "F"]  # 7 == 7.0


def test_structures_compare_as_molecules(table) -> None:
    # CCO and OCC are ethanol; the Kekulé and aromatic benzenes are one molecule
    assert ids(table, picks.first_occurrences(table, "SMILES", every(table))) == ["A", "C", "F"]
    assert ids(table, picks.duplicates(table, "SMILES", every(table))) == ["B", "E"]


def test_the_order_given_decides_which_is_first(table) -> None:
    table.set_sort([("ID", False)])
    assert ids(table, picks.first_occurrences(table, "Series", every(table))) == ["F", "E"]
    shown = every(table)[:3]  # F, E, D: as if a filter hid the rest
    assert ids(table, picks.duplicates(table, "Series", shown)) == ["D"]
