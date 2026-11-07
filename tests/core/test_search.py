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


# -- field syntax -------------------------------------------------------------------------------
def _ids(table: Table, flt) -> list[str]:
    keep = flt.mask(table)
    return [table.cell("ID", int(o)) for o, k in zip(table.oids, keep, strict=True) if k]


def _terms(table: Table, text: str):
    from chemsuite.core.table.search import parse_terms

    return SearchFilter(Query(text), terms=parse_terms(text, table.column_names))


def test_parse_terms() -> None:
    from chemsuite.core.table.search import Term, parse_terms

    columns = ["ID", "pIC50", "Mol Wt", "Series"]
    assert parse_terms('id:CPD-1 pIC50>7 "kinase a" -Series:B', columns) == (
        Term("CPD-1", "ID"),
        Term("7", "pIC50", ">", 7.0),
        Term("kinase a"),
        Term("B", "Series", negate=True),
    )
    assert parse_terms('"Mol Wt"<=400', columns) == (Term("400", "Mol Wt", "<=", 400.0),)
    assert parse_terms("Series=A Series!=B", columns) == (
        Term("A", "Series", "equals"),
        Term("B", "Series", "equals", negate=True),
    )
    assert parse_terms("pIC50=7.5", columns) == (Term("7.5", "pIC50", "=", 7.5),)
    # not a column, or still being typed: plain text, or skipped
    assert parse_terms("http://x 10:30", columns) == (Term("http://x"), Term("10:30"))
    assert parse_terms("ID: pIC50>", columns) == ()
    assert parse_terms('"open phrase', columns) == (Term("open phrase"),)
    assert parse_terms('"-5"', columns) == (Term("-5"),)  # quoted: not an exclusion


def test_terms_mask(table: Table) -> None:
    assert _ids(table, _terms(table, "kinase")) == ["CPD-1", "cpd-2"]
    assert _ids(table, _terms(table, "kinase b")) == ["cpd-2"]  # every word, anywhere
    assert _ids(table, _terms(table, '"kinase b"')) == ["cpd-2"]
    assert _ids(table, _terms(table, "ID:cpd-1")) == ["CPD-1", "CPD-12"]
    assert _ids(table, _terms(table, "pIC50>=7.5")) == ["CPD-1", "CPD-12", "X-9"]
    assert _ids(table, _terms(table, "pIC50=7.5")) == ["CPD-1", "X-9"]  # 7.50 and 7.5
    assert _ids(table, _terms(table, "pIC50!=7.5")) == ["cpd-2", "CPD-12"]
    assert _ids(table, _terms(table, "pIC50>7 -Assay:gpcr")) == ["CPD-1", "CPD-12"]
    assert _ids(table, _terms(table, 'Assay="kinase a"')) == ["CPD-1"]  # whole cell, any case
    # the bar's column applies to terms without a field
    flt = _terms(table, "7.5")
    assert _ids(table, SearchFilter(flt.query, column="ID", terms=flt.terms)) == []


def test_a_term_follows_a_renamed_column(table: Table) -> None:
    flt = _terms(table, "pIC50>7")
    table.set_filter("search", flt)
    table.rename_column("pIC50", "Potency")
    renamed = table.filters["search"]
    assert renamed.terms[0].column == "Potency"
    assert _ids(table, renamed) == ["CPD-1", "CPD-12", "X-9"]


# -- find and replace: scope and groups -------------------------------------------------------
def test_find_within_cells(table: Table) -> None:
    oids = [int(o) for o in table.oids]
    within = {"ID": [oids[0], oids[2]], "Assay": [oids[1]]}
    matches = find(table, Query("cpd"), within=within)
    assert [matches.cell(i) for i in range(len(matches))] == [(oids[0], "ID"), (oids[2], "ID")]
    assert len(find(table, Query("kinase"), within=within)) == 1  # Assay only in row 2
    assert len(find(table, Query("kinase"), ["ID"], within=within)) == 0  # both limits apply


def test_regex_replacement_groups(table: Table) -> None:
    from chemsuite.core.table.search import regex_template

    assert regex_template(r"$2-${name} $$ \1") == r"\g<2>-\g<name> $ \1"
    query = Query(r"(?i)(cpd)-(?P<n>\d+)", "regex")
    assert query.replace("CPD-12", "${n}/$1") == "12/CPD"
    assert query.replace("CPD-12", r"\g<n>:\1") == "12:CPD"
    assert query.replacement_problem("$1 ${n}") == ""
    assert "invalid group reference 3" in query.replacement_problem("$3")
    assert "unknown group name 'x'" in query.replacement_problem("${x}")
    assert Query("cpd").replacement_problem("$9") == ""  # literal outside regex mode
    oids = [int(o) for o in table.oids]
    within = {"ID": [oids[2]]}
    assert replace_all(table, Query(r"CPD-(\d+)", "regex", True), "C$1", within=within) == 1
    assert [table.cell("ID", o) for o in oids] == ["CPD-1", "cpd-2", "C12", "X-9"]
