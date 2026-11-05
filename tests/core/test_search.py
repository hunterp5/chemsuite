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
from chemsuite.core.table.columns import StructureColumn
from chemsuite.core.table.filters import SearchFilter
from chemsuite.core.table.search import Query, find, replace_all, replace_cell


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "Structure": StructureColumn.from_texts(["CCO", "c1ccccc1", "CCN", ""]),
            "ID": ["CPD-1", "cpd-2", "CPD-12", "X-9"],
            "Assay": ["kinase A", "Kinase B", "", "GPCR"],
            "pIC50": ["7.50", "6.2", "8", "7.5"],
        },
        undoable=False,
    )
    return t


def cells(matches) -> list[tuple[int, str]]:
    return [matches.cell(i) for i in range(len(matches))]


def test_find_is_row_major_in_display_order_and_skips_structures(table: Table) -> None:
    assert cells(find(table, Query("c"))) == [(1, "ID"), (2, "ID"), (3, "ID"), (4, "Assay")]
    table.set_sort([("ID", False)])  # natural order, descending: X-9, CPD-12, cpd-2, CPD-1
    assert [oid for oid, _ in cells(find(table, Query("cpd")))] == [3, 2, 1]


def test_modes_and_case(table: Table) -> None:
    assert cells(find(table, Query("7.5", "whole"))) == [(4, "pIC50")]
    assert cells(find(table, Query("CPD", case_sensitive=True))) == [(1, "ID"), (3, "ID")]
    assert cells(find(table, Query(r"^CPD-\d$", "regex"))) == [(1, "ID"), (2, "ID")]
    assert len(find(table, Query(""))) == 0
    with pytest.raises(Exception):  # noqa: B017 (re.error)
        Query("[", "regex").predicate()


def test_shown_only_and_hidden_rows(table: Table) -> None:
    from chemsuite.core.table import TextFilter

    table.set_filter("f", TextFilter("Assay", "kinase"))
    assert [o for o, _ in cells(find(table, Query("CPD")))] == [1, 2]
    assert [o for o, _ in cells(find(table, Query("CPD"), shown_only=False))] == [1, 2, 3]


def test_replace_all_is_one_undo_step(table: Table) -> None:
    changed = replace_all(table, Query("kinase"), "Kinase")
    assert changed == 1  # "Kinase B" already reads the same
    assert table.cell("Assay", 1) == "Kinase A"
    assert (
        replace_all(table, Query(r"CPD-(\d+)", "regex", case_sensitive=True), r"C\1", ["ID"]) == 2
    )
    assert [table.cell("ID", o) for o in (1, 2, 3)] == ["C1", "cpd-2", "C12"]
    table.undo()
    assert table.cell("ID", 3) == "CPD-12"
    assert table.history.undo_label == "Replace “kinase”"


def test_replace_in_numbers_demotes_only_when_needed(table: Table) -> None:
    replace_all(table, Query("7.5", "whole"), "7.6")
    assert table.kind("pIC50") == "number" and table.cell("pIC50", 4) == "7.6"
    replace_all(table, Query("8", "whole"), ">8")
    assert table.kind("pIC50") == "text" and table.cell("pIC50", 3) == ">8"


def test_literal_replacement_keeps_backslashes(table: Table) -> None:
    assert replace_cell(table, 4, "ID", Query("X-9"), r"a\1b")
    assert table.cell("ID", 4) == r"a\1b"
    assert not replace_cell(table, 4, "ID", Query("zzz"), "q")


def test_search_filter_matches_any_column(table: Table) -> None:
    table.set_filter("search", SearchFilter(Query("7.5")))
    assert table.visible_oids().tolist() == [1, 4]
    table.set_filter("search", SearchFilter(Query("ccc")))  # SMILES text is not searched
    assert table.visible_oids().tolist() == []
    table.remove_column("Assay")
    assert "search" in table.filters  # not tied to a column
