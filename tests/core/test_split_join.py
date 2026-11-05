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

import pytest

from chemsuite.core.table import Table
from chemsuite.core.table.split_join import (
    detect_delimiter,
    join_columns,
    keep_fields,
    plan_split,
    split_column,
    split_text,
)


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "ID": ["A", "B", "C", "D"],
            "Assay": ["IC50; 5.2", "Ki; 7", "", "EC50; 6.1; flagged"],
            "Series": ["S1", "", "S2", "S3"],
        },
        undoable=False,
    )
    return t


def test_detect_delimiter() -> None:
    assert detect_delimiter(["a;b", "c;d,e"]) == ";"
    assert detect_delimiter(["a b", "c"]) == " "
    assert detect_delimiter(["abc"]) == ","


def test_split_text_quotes_and_whitespace() -> None:
    assert split_text('"a, b", c', ",") == ["a, b", "c"]
    assert split_text(" x  y ", " ") == ["x", "y"]
    assert split_text("a| b |", "|") == ["a", "b", ""]
    assert split_text("   ", ",") == []


def test_keep_fields_prefers_numbers() -> None:
    assert keep_fields(["x", "5", "12", ""], "largest") == ["12"]
    assert keep_fields(["x", "5", "12"], "smallest") == ["5"]
    assert keep_fields(["b", "a"], "smallest") == ["a"]
    assert keep_fields(["", "b", "a", ""], "first") == ["b"]
    assert keep_fields([], "last") == [""]


def test_split_column_inserts_typed_columns_after_source(table: Table) -> None:
    names = split_column(table, "Assay")  # auto-detects ";"
    assert names == ["Assay_1", "Assay_2", "Assay_3"]
    assert table.column_names == ["ID", "Assay", "Assay_1", "Assay_2", "Assay_3", "Series"]
    assert table.kind("Assay_2") == "number"
    assert [table.cell("Assay_3", o) for o in table.oids] == ["", "", "", "flagged"]
    assert table.history.undo_label == "Split Assay"
    table.undo()
    assert table.column_names == ["ID", "Assay", "Series"]


def test_split_only_some_rows_and_keep_one(table: Table) -> None:
    oids = table.oids[:2]
    names = split_column(table, "Assay", ";", prefix="P", keep="largest", oids=oids)
    assert names == ["P_largest"]
    assert [table.cell("P_largest", o) for o in table.oids] == ["5.2", "7", "", ""]


def test_split_refuses_too_many_columns(table: Table) -> None:
    table.add_column("Long", [",".join("x" * 300)] * 4)
    with pytest.raises(ValueError, match="limit"):
        split_column(table, "Long", ",")


def test_split_existing_name_gets_suffix(table: Table) -> None:
    table.add_column("Assay_1", ["x"] * 4)
    assert split_column(table, "Assay", ";")[0] == "Assay_1 (2)"


def test_plan_split_reports_delimiter(table: Table) -> None:
    plan = plan_split(table, "Assay")
    assert plan.delimiter == ";" and plan.width == 3


def test_join_columns(table: Table) -> None:
    name = join_columns(table, ["ID", "Series"], "-")
    assert name == "ID_Series"
    assert table.column_names[-1] == "ID_Series"
    assert [table.cell(name, o) for o in table.oids] == ["A-S1", "B", "C-S2", "D-S3"]
    keep = join_columns(table, ["ID", "Series"], "-", name="Both", skip_blank=False)
    assert table.cell(keep, int(table.oids[1])) == "B-"
    table.undo()
    assert "Both" not in table


def test_join_needs_two_columns(table: Table) -> None:
    with pytest.raises(ValueError):
        join_columns(table, ["ID"], ",")
