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

from chemsuite.core.lists import ListError, RowLists
from chemsuite.core.table import Table


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows({"ID": list("ABCDE")}, undoable=False)
    return t


def test_add_remove_and_counts(table: Table) -> None:
    lists = RowLists(table)
    oids = table.oids
    assert lists.add(" hits ", oids[[0, 2]]) == 2
    assert lists.names == ["hits"]
    assert lists.add("hits", [int(oids[2]), int(oids[3])]) == 1
    assert lists.oids("hits").tolist() == oids[[0, 2, 3]].tolist()
    assert lists.remove("hits", oids[[0, 4]]) == 1
    assert lists.count("hits") == 2
    assert lists.remove("hits", oids) == 2 and "hits" in lists  # empty lists stay


def test_revision_counts_only_real_changes(table: Table) -> None:
    lists = RowLists(table)
    calls: list[int] = []
    lists.subscribe(lambda: calls.append(lists.revision))
    lists.add("a", table.oids[:2])
    lists.add("a", table.oids[:1])  # nothing new
    lists.remove("a", table.oids[3:])  # nothing there
    lists.rename("a", "a")
    assert calls == [1]


def test_rename_and_delete(table: Table) -> None:
    lists = RowLists(table)
    lists.add("a", table.oids[:1])
    lists.add("b", table.oids[1:2])
    lists.rename("a", "c")
    assert lists.names == ["c", "b"]  # keeps its place
    with pytest.raises(ListError):
        lists.rename("c", "b")
    with pytest.raises(ListError):
        lists.add("  ", table.oids)
    lists.delete("c")
    assert lists.names == ["b"]


def test_deleted_rows_stop_counting_but_come_back_on_undo(table: Table) -> None:
    lists = RowLists(table)
    lists.add("a", table.oids[:3])
    table.remove_rows(table.oids[:2])
    assert lists.count("a") == 1
    table.undo()
    assert lists.count("a") == 3


def test_json_round_trip_drops_deleted_rows_and_bad_entries(table: Table) -> None:
    lists = RowLists(table)
    lists.add("a", table.oids[:3])
    lists.add("b", table.oids[4:])
    table.remove_rows(table.oids[:1])
    data = lists.to_json()
    assert data[0] == {"name": "a", "oids": table.oids[:2].tolist()}
    other = RowLists(table)
    other.load_json([*data, {"name": "", "oids": [1]}, {"name": "x", "oids": ["?"]}, 3])
    assert other.names == ["a", "b"]
    assert other.oids("b").tolist() == lists.oids("b").tolist()
    other.load_json(None)
    assert other.names == []
