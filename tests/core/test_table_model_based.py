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

"""Property-based tests: random edits, undo and redo checked against a plain-Python model.

The model stores every cell as its display string in dicts. If the table ever disagrees with it
after any sequence of operations (including undo/redo across type changes and deletions), the
test fails and Hypothesis shrinks the sequence to a minimal reproduction.
"""

from __future__ import annotations

import copy

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.stateful import (
    RuleBasedStateMachine,
    initialize,
    invariant,
    precondition,
    rule,
)

from chemsuite.core.table import NUMBER, Table
from chemsuite.core.table.columns import natural_key

CELL = st.one_of(
    st.none(),
    st.sampled_from(["", "1", "2.50", "-3", "007", "1e5", ">10", "A", "b", "CPD-2", "CPD-10"]),
    st.integers(-1000, 1000).map(str),
    st.decimals(-100, 100, places=2).map(str),
    st.text(alphabet="abcXYZ-_ 09", max_size=6),
)
NAMES = st.sampled_from(["A", "B", "C", "D", "E"])


class Model:
    def __init__(self) -> None:
        self.columns: list[str] = []
        self.rows: dict[int, dict[str, str]] = {}


class TableMachine(RuleBasedStateMachine):
    @initialize()
    def start(self) -> None:
        self.table = Table()
        self.model = Model()
        self.undo_states: list[Model] = []
        self.redo_states: list[Model] = []
        self.next_oid = 1

    # -- helpers ------------------------------------------------------------------------------
    def _record(self) -> None:
        self.undo_states.append(copy.deepcopy(self.model))
        self.redo_states.clear()

    def _oids(self) -> list[int]:
        return sorted(self.model.rows)

    # -- rules --------------------------------------------------------------------------------
    @rule(data=st.data(), n=st.integers(0, 4))
    def add_rows(self, data, n: int) -> None:
        names = data.draw(st.lists(NAMES, unique=True, max_size=3))
        values = {name: data.draw(st.lists(CELL, min_size=n, max_size=n)) for name in names}
        self._record()
        oids = self.table.add_rows(values, count=n)
        assert oids.tolist() == list(range(self.next_oid, self.next_oid + n))
        self.next_oid += n
        for name in names:
            if name not in self.model.columns:
                self.model.columns.append(name)
                for row in self.model.rows.values():
                    row[name] = ""
        for i, oid in enumerate(oids.tolist()):
            self.model.rows[oid] = {
                c: ("" if values.get(c, [None] * n)[i] is None else values[c][i])
                for c in self.model.columns
            }

    @precondition(lambda self: self.model.rows)
    @rule(data=st.data())
    def remove_rows(self, data) -> None:
        oids = data.draw(st.lists(st.sampled_from(self._oids()), min_size=1, unique=True))
        self._record()
        self.table.remove_rows(oids)
        for oid in oids:
            del self.model.rows[oid]

    @precondition(lambda self: self.model.rows and self.model.columns)
    @rule(data=st.data())
    def set_cells(self, data) -> None:
        name = data.draw(st.sampled_from(self.model.columns))
        oids = data.draw(st.lists(st.sampled_from(self._oids()), min_size=1, unique=True))
        values = data.draw(st.lists(CELL, min_size=len(oids), max_size=len(oids)))
        self._record()
        self.table.set_cells(name, oids, values)
        for oid, value in zip(oids, values, strict=True):
            self.model.rows[oid][name] = "" if value is None else value

    @rule(name=NAMES)
    def add_column(self, name: str) -> None:
        if name in self.model.columns:
            return
        self._record()
        self.table.add_column(name)
        self.model.columns.append(name)
        for row in self.model.rows.values():
            row[name] = ""

    @precondition(lambda self: self.model.columns)
    @rule(data=st.data())
    def remove_column(self, data) -> None:
        name = data.draw(st.sampled_from(self.model.columns))
        self._record()
        self.table.remove_column(name)
        self.model.columns.remove(name)
        for row in self.model.rows.values():
            del row[name]

    @precondition(lambda self: self.model.columns)
    @rule(data=st.data(), new=NAMES)
    def rename_column(self, data, new: str) -> None:
        old = data.draw(st.sampled_from(self.model.columns))
        if new in self.model.columns:
            return
        self._record()
        self.table.rename_column(old, new)
        self.model.columns[self.model.columns.index(old)] = new
        for row in self.model.rows.values():
            row[new] = row.pop(old)

    @rule()
    def undo(self) -> None:
        label = self.table.undo()
        assert (label is None) == (not self.undo_states)
        if label is not None:
            self.redo_states.append(self.model)
            self.model = self.undo_states.pop()

    @rule()
    def redo(self) -> None:
        label = self.table.redo()
        assert (label is None) == (not self.redo_states)
        if label is not None:
            self.undo_states.append(self.model)
            self.model = self.redo_states.pop()

    # -- the check ----------------------------------------------------------------------------
    @invariant()
    def table_matches_model(self) -> None:
        if not hasattr(self, "table"):
            return
        assert self.table.column_names == self.model.columns
        assert self.table.oids.tolist() == self._oids()
        for name in self.model.columns:
            column = self.table.column(name)
            assert len(column) == len(self.model.rows)
            got = column.displays(np.arange(len(column)))
            assert got == [self.model.rows[oid][name] for oid in self._oids()], name


TestTableAgainstModel = TableMachine.TestCase
TestTableAgainstModel.settings = settings(max_examples=150, stateful_step_count=40, deadline=None)


@given(
    st.lists(st.one_of(st.none(), st.floats(-1e6, 1e6), st.integers(-50, 50)), max_size=60),
    st.booleans(),
)
def test_number_sort_matches_python(values, ascending: bool) -> None:
    table = Table()
    table.add_rows(count=len(values))
    table.add_column("x", values, kind=NUMBER)
    table.set_sort([("x", ascending)])
    got = table.read(["x"], table.visible_oids())["x"].tolist()
    present = sorted((float(v) for v in values if v is not None), reverse=not ascending)
    expected = present + [None] * (len(values) - len(present))
    assert [None if np.isnan(g) else g for g in got] == expected


@given(
    st.lists(st.one_of(st.none(), st.text(alphabet="aAbB-0129", max_size=5)), max_size=60),
    st.booleans(),
)
def test_text_sort_matches_python(values, ascending: bool) -> None:
    table = Table()
    table.add_rows({"x": [None if v is None else "t" + v for v in values]})  # keep as text
    table.set_sort([("x", ascending)])
    got = [table.cell("x", int(o)) for o in table.visible_oids()]
    present = sorted((v for v in got if v), key=natural_key, reverse=not ascending)
    assert got == present + [""] * (len(got) - len(present))
