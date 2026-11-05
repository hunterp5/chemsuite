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

from chemsuite.core.table import (
    NUMBER,
    TEXT,
    AboutToChange,
    CategoryFilter,
    Changed,
    OidFilter,
    RangeFilter,
    Table,
    TextFilter,
)
from chemsuite.core.table import events as ev


def cells(table: Table, name: str) -> list[str]:
    """Visible cells of one column, in display order."""
    return [table.cell(name, int(oid)) for oid in table.visible_oids()]


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "ID": ["CPD-10", "CPD-2", "CPD-1", "cpd-3"],
            "pIC50": ["7.50", "6.2", None, "8"],
            "Series": ["A", "B", "A", ""],
        },
        undoable=False,
    )
    return t


# -- columns and types --------------------------------------------------------------------------


def test_kinds_are_inferred_and_display_is_lossless(table: Table) -> None:
    assert table.kind("ID") == TEXT
    assert table.kind("pIC50") == NUMBER
    assert cells(table, "pIC50") == ["7.50", "6.2", "", "8"]


def test_computed_values_and_python_numbers() -> None:
    t = Table()
    t.add_rows({"MW": [350.123456, 12, float("nan")]})
    assert t.kind("MW") == NUMBER
    assert cells(t, "MW") == ["350.1235", "12", ""]
    np.testing.assert_array_equal(t.read(["MW"])["MW"][:2], [350.123456, 12.0])


def test_oids_are_stable_and_lookup_rejects_unknown(table: Table) -> None:
    assert table.oids.tolist() == [1, 2, 3, 4]
    with pytest.raises(KeyError):
        table.slots_of([99])
    with pytest.raises(KeyError):
        Table().slots_of([1])
    with pytest.raises(ValueError):
        table.oids[0] = 5  # read-only view


def test_text_into_number_column_demotes_it_in_one_undoable_step(table: Table) -> None:
    table.set_cells("pIC50", [1], [">10"])
    assert table.kind("pIC50") == TEXT
    assert cells(table, "pIC50") == [">10", "6.2", "", "8"]
    assert table.undo() == "Edit cells"
    assert table.kind("pIC50") == NUMBER
    assert cells(table, "pIC50") == ["7.50", "6.2", "", "8"]
    table.redo()
    assert table.kind("pIC50") == TEXT


def test_convert_column(table: Table) -> None:
    table.set_cells("pIC50", [1], ["n/a"])
    with pytest.raises(ValueError):
        table.convert_column("pIC50", NUMBER)
    assert table.convert_column("pIC50", NUMBER, force=True) == 1
    assert cells(table, "pIC50") == ["", "6.2", "", "8"]
    table.undo()
    assert cells(table, "pIC50") == ["n/a", "6.2", "", "8"]


def test_column_add_remove_rename(table: Table) -> None:
    table.add_column("Notes", position=1)
    assert table.column_names == ["ID", "Notes", "pIC50", "Series"]
    with pytest.raises(ValueError):
        table.add_column("Notes")
    table.rename_column("Notes", "Comment")
    table.remove_column("ID")
    assert table.column_names == ["Comment", "pIC50", "Series"]
    table.undo()
    table.undo()
    table.undo()
    assert table.column_names == ["ID", "pIC50", "Series"]


def test_new_rows_fill_missing_columns_with_blanks(table: Table) -> None:
    oids = table.add_rows({"ID": ["CPD-5"], "LogP": ["2.1"]})
    assert oids.tolist() == [5]
    assert table.cell("pIC50", 5) == ""
    assert table.cell("LogP", 5) == "2.1" and table.cell("LogP", 1) == ""
    assert table.kind("LogP") == NUMBER
    table.undo()
    assert len(table) == 4 and "LogP" not in table


def test_deleted_rows_come_back_in_place_with_their_oids(table: Table) -> None:
    table.remove_rows([2, 3])
    assert table.oids.tolist() == [1, 4]
    table.undo()
    assert table.oids.tolist() == [1, 2, 3, 4]
    assert cells(table, "ID") == ["CPD-10", "CPD-2", "CPD-1", "cpd-3"]
    # oids are never reused, even after undo
    assert table.add_rows(count=1).tolist() == [5]


# -- view: sort and filter ----------------------------------------------------------------------


def test_text_sort_is_natural_and_case_insensitive(table: Table) -> None:
    table.set_sort([("ID", True)])
    assert cells(table, "ID") == ["CPD-1", "CPD-2", "cpd-3", "CPD-10"]
    table.set_sort([("ID", False)])
    assert cells(table, "ID") == ["CPD-10", "cpd-3", "CPD-2", "CPD-1"]


def test_blanks_sort_last_in_both_directions(table: Table) -> None:
    table.set_sort([("pIC50", True)])
    assert cells(table, "pIC50") == ["6.2", "7.50", "8", ""]
    table.set_sort([("pIC50", False)])
    assert cells(table, "pIC50") == ["8", "7.50", "6.2", ""]


def test_multi_column_sort_is_stable(table: Table) -> None:
    table.set_sort([("Series", True), ("pIC50", False)])
    assert cells(table, "ID") == ["CPD-10", "CPD-1", "CPD-2", "cpd-3"]


def test_sort_never_reorders_storage(table: Table) -> None:
    table.set_sort([("ID", True)])
    assert table.oids.tolist() == [1, 2, 3, 4]


def test_filters_combine_with_and(table: Table) -> None:
    table.set_filter("potency", RangeFilter("pIC50", low=7))
    assert cells(table, "ID") == ["CPD-10", "cpd-3"]
    table.set_filter("series", CategoryFilter("Series", frozenset({"A"})))
    assert cells(table, "ID") == ["CPD-10"]
    table.set_filter("series", None)
    table.set_filter("potency", RangeFilter("pIC50", low=7, include_blank=True))
    assert cells(table, "ID") == ["CPD-10", "CPD-1", "cpd-3"]


def test_text_filter_modes(table: Table) -> None:
    table.set_filter("t", TextFilter("ID", "cpd-1"))
    assert cells(table, "ID") == ["CPD-10", "CPD-1"]
    table.set_filter("t", TextFilter("ID", "cpd-1", mode="equals"))
    assert cells(table, "ID") == ["CPD-1"]
    table.set_filter("t", TextFilter("ID", "cpd-1", case_sensitive=True))
    assert cells(table, "ID") == []
    table.set_filter("t", TextFilter("ID", r"-\d$", mode="regex", negate=True))
    assert cells(table, "ID") == ["CPD-10"]
    table.set_filter("t", TextFilter("pIC50", ".5"))
    assert cells(table, "ID") == ["CPD-10"]


def test_oid_filter(table: Table) -> None:
    table.set_filter("hits", OidFilter(frozenset({4, 2})))
    assert table.visible_oids().tolist() == [2, 4]


def test_edits_to_a_filtered_column_update_the_view(table: Table) -> None:
    table.set_filter("potency", RangeFilter("pIC50", low=7))
    table.set_cells("pIC50", [2], ["9.1"])
    assert cells(table, "ID") == ["CPD-10", "CPD-2", "cpd-3"]


def test_renaming_or_removing_a_column_updates_filters_and_sort(table: Table) -> None:
    table.set_filter("potency", RangeFilter("pIC50", low=7))
    table.set_sort([("pIC50", True)])
    table.rename_column("pIC50", "pKi")
    assert table.filters["potency"].columns == ("pKi",)
    assert table.sort_keys == (("pKi", True),)
    assert cells(table, "ID") == ["CPD-10", "cpd-3"]
    table.remove_column("pKi")
    assert table.filters == {} and table.sort_keys == ()
    assert len(table.view()) == 4


# -- results from tools -------------------------------------------------------------------------


def test_apply_results_reuses_blank_columns_and_never_overwrites_data(table: Table) -> None:
    table.add_column("MW")
    names = table.apply_results({"MW": ([1, 2], [300.5, 250.25]), "pIC50": ([1], [5.0])})
    assert names == {"MW": "MW", "pIC50": "pIC50 (2)"}
    assert table.kind("MW") == NUMBER
    assert cells(table, "MW") == ["300.5", "250.25", "", ""]
    assert cells(table, "pIC50") == ["7.50", "6.2", "", "8"]
    assert cells(table, "pIC50 (2)") == ["5", "", "", ""]
    table.undo()
    assert table.column_names == ["ID", "pIC50", "Series", "MW"]
    assert table.column("MW").is_all_blank()


# -- events -------------------------------------------------------------------------------------


def test_each_operation_emits_one_bracketed_event_pair(table: Table) -> None:
    seen: list[object] = []
    unsubscribe = table.subscribe(seen.append)
    table.set_cells("Series", [1, 2], ["C", "C"])
    assert isinstance(seen[0], AboutToChange) and not seen[0].structural
    assert isinstance(seen[1], Changed) and not seen[1].structural
    (change,) = seen[1].changes
    assert isinstance(change, ev.CellsChanged) and change.oids.tolist() == [1, 2]

    seen.clear()
    table.set_sort([("Series", True)])
    table.set_cells("Series", [3], ["Z"])  # now affects the view
    assert [e.structural for e in seen] == [True, True, True, True]

    seen.clear()
    unsubscribe()
    table.remove_rows([1])
    assert seen == []


def test_failed_command_leaves_table_unchanged_and_still_notifies(table: Table) -> None:
    seen: list[object] = []
    table.subscribe(seen.append)
    before = cells(table, "ID")
    with pytest.raises(KeyError):
        table.set_cells("ID", [1, 99], ["x", "y"])
    assert cells(table, "ID") == before
    assert not table.history.can_undo


def test_history_clean_state() -> None:
    t = Table()
    t.history.mark_clean()
    assert t.history.is_clean
    t.add_rows({"A": ["1"]})
    assert not t.history.is_clean
    t.undo()
    assert t.history.is_clean
    t.redo()
    t.history.mark_clean()
    t.undo()
    t.add_rows({"B": ["x"]})  # the saved state is now unreachable
    t.undo()
    assert not t.history.is_clean


def test_set_block_is_one_undo_step(table: Table) -> None:
    table.set_block({"ID": ([1, 2], ["X", "Y"]), "pIC50": ([1], ["n/a"])}, label="Paste")
    assert cells(table, "ID")[:2] == ["X", "Y"] and table.kind("pIC50") == TEXT
    assert table.undo() == "Paste"
    assert cells(table, "ID")[:2] == ["CPD-10", "CPD-2"] and table.kind("pIC50") == NUMBER


def test_snapshot_is_independent(table: Table) -> None:
    from chemsuite.core.table.columns import StructureColumn

    table.add_column("S", StructureColumn.from_texts(["CCO", "", "", "C"]))
    snap = table.snapshot()
    table.set_cells("ID", [1], ["changed"])
    table.fill_structures("S", [1], [b"x"])
    table.remove_rows([2])
    assert snap.cell("ID", 1) == "CPD-10" and len(snap) == 4
    assert snap.column("S").pending_mask().tolist() == [True, False, False, True]
    assert snap.next_oid == table.next_oid


def test_mark_clean_at_checkpoint() -> None:
    t = Table()
    t.add_rows({"a": ["1"]})
    saved = t.history.checkpoint()
    t.set_cells("a", [1], ["2"])  # edited while the save was running
    t.history.mark_clean(saved)
    assert not t.history.is_clean
    t.undo()
    assert t.history.is_clean


def test_listeners_see_history_updated() -> None:
    t = Table()
    t.add_rows({"a": ["1"]}, undoable=False)
    seen: list[tuple[str | None, bool]] = []

    def listener(event) -> None:
        if isinstance(event, Changed):
            seen.append((t.history.undo_label, t.history.is_clean))

    t.subscribe(listener)
    t.set_cells("a", [1], ["2"], label="Edit a")
    t.undo()
    assert seen == [("Edit a", False), (None, True)]
