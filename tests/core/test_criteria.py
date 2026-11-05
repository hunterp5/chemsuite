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

from __future__ import annotations

import numpy as np
import pytest

from chemsuite.core.table import Table
from chemsuite.core.table.criteria import (
    ANY_COLUMN,
    CriteriaFilter,
    CriteriaSearch,
    Criterion,
    CriterionError,
)


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "ID": ["CPD-1", "CPD-2", "cpd-3", "X-4", "CPD-5"],
            "pIC50": ["5.5", "7.2", "", "8.1", "6.0"],
            "Series": ["A", "B", "A", "", "B"],
            "SMILES": ["CCO", "c1ccccc1", "CCN", "", "c1ccncc1"],
        },
        undoable=False,
    )
    t.convert_column("pIC50", "number")
    t.convert_column("SMILES", "structure")
    return t


def rows(table: Table, mask: np.ndarray) -> list[str]:
    return [table.cell("ID", int(o)) for o in table.oids[mask]]


def test_text_operators(table) -> None:
    def hit(op, value, **kw):
        return rows(table, Criterion("ID", op, value, **kw).matches(table))

    assert hit("contains", "cpd") == ["CPD-1", "CPD-2", "cpd-3", "CPD-5"]
    assert hit("contains", "cpd", case_sensitive=True) == ["cpd-3"]
    assert hit("equals", "x-4") == ["X-4"]
    assert hit("starts", "CPD", case_sensitive=True) == ["CPD-1", "CPD-2", "CPD-5"]
    assert hit("ends", "-5") == ["CPD-5"]
    assert hit("regex", r"^CPD-[12]$") == ["CPD-1", "CPD-2"]
    assert rows(table, Criterion("Series", "blank").matches(table)) == ["X-4"]


def test_number_operators(table) -> None:
    def hit(op, value, value2=""):
        return rows(table, Criterion("pIC50", op, value, value2).matches(table))

    assert hit("gt", "6") == ["CPD-2", "X-4"]
    assert hit("ge", "6") == ["CPD-2", "X-4", "CPD-5"]
    assert hit("lt", "6") == ["CPD-1"]
    assert hit("le", "6") == ["CPD-1", "CPD-5"]
    assert hit("eq", "6.0") == ["CPD-5"]
    assert hit("ne", "6") == ["CPD-1", "CPD-2", "X-4"]  # a blank is not "≠ 6"
    assert hit("between", "8", "6") == ["CPD-2", "CPD-5"]  # bounds in either order
    assert hit("blank", "") == ["cpd-3"]


def test_negate_and_any_column(table) -> None:
    assert rows(table, Criterion(ANY_COLUMN, "equals", "b").matches(table)) == ["CPD-2", "CPD-5"]
    assert rows(table, Criterion(ANY_COLUMN, "contains", "7.2").matches(table)) == ["CPD-2"]
    # Structure cells are never text-searched.
    assert rows(table, Criterion(ANY_COLUMN, "contains", "ccccc").matches(table)) == []
    not_a = Criterion("Series", "equals", "A", negate=True)
    assert rows(table, not_a.matches(table)) == ["CPD-2", "X-4", "CPD-5"]


def test_all_and_any(table) -> None:
    active = Criterion("pIC50", "ge", "6")
    series_b = Criterion("Series", "equals", "B")
    assert rows(table, CriteriaSearch((active, series_b), "all").mask(table)) == [
        "CPD-2",
        "CPD-5",
    ]
    assert rows(table, CriteriaSearch((active, series_b), "any").mask(table)) == [
        "CPD-2",
        "X-4",
        "CPD-5",
    ]


def test_structure_conditions_use_hits_found_elsewhere(table) -> None:
    ring = Criterion("SMILES", "substructure", "c1ccccc1")
    search = CriteriaSearch((ring, Criterion("pIC50", "gt", "7")))
    with pytest.raises(CriterionError, match="has not run"):
        search.mask(table)
    hits = {0: table.oids[[1]]}
    assert rows(table, search.mask(table, hits)) == ["CPD-2"]
    assert rows(table, Criterion("SMILES", "blank").matches(table)) == ["X-4"]


def test_checks(table) -> None:
    with pytest.raises(CriterionError, match="not a number"):
        Criterion("pIC50", "gt", "high").check(table)
    with pytest.raises(CriterionError, match="regular expression"):
        Criterion("ID", "regex", "(").check(table)
    with pytest.raises(CriterionError, match="no column"):
        Criterion("Gone", "contains", "x").check(table)
    with pytest.raises(CriterionError, match="does not apply"):
        Criterion("pIC50", "contains", "5").check(table)
    with pytest.raises(CriterionError, match="at least one"):
        CriteriaSearch(()).check(table)


def test_filter_follows_renames_and_is_dropped_with_its_column(table) -> None:
    search = CriteriaSearch((Criterion("Series", "equals", "A"),))
    table.set_filter("advanced", CriteriaFilter(search))
    assert rows(table, np.isin(table.oids, table.oids[table.view()])) == ["CPD-1", "cpd-3"]
    table.rename_column("Series", "Group")
    assert table.filters["advanced"].columns == ("Group",)
    assert len(table.view()) == 2
    table.remove_columns(["Group"])
    assert "advanced" not in table.filters


def test_json_round_trip() -> None:
    search = CriteriaSearch(
        (
            Criterion("pIC50", "between", "6", "8", negate=True),
            Criterion(ANY_COLUMN, "contains", "x", case_sensitive=True),
        ),
        "any",
    )
    assert CriteriaSearch.from_json(search.to_json()) == search
    assert search.describe() == "not (pIC50 between 6 and 8) or any column contains x"
    with pytest.raises(ValueError):
        CriteriaSearch.from_json({"combine": "most", "criteria": []})
