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

import sqlite3
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from chemsuite.core import document
from chemsuite.core.document import DocumentError, load_document, save_document
from chemsuite.core.mols import parse_structures
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn


def snapshot(table: Table) -> dict:
    slots = np.arange(len(table))
    out: dict = {
        "oids": table.oids.tolist(),
        "next_oid": table.next_oid,
        "sort": table.sort_keys,
        "kinds": [table.kind(n) for n in table.column_names],
    }
    for name in table.column_names:
        column = table.column(name)
        out[name] = column.displays(slots)
        if isinstance(column, StructureColumn):
            out[name + ":state"] = (
                column.pending_mask().tolist(),
                column.invalid_mask().tolist(),
                column.blobs(slots),
            )
    return out


@pytest.fixture
def table() -> Table:
    t = Table()
    t.add_rows(
        {
            "Structure": StructureColumn.from_texts(["CCO", "C1CC", "", "c1ccccc1"]),
            "ID": ["A", "B", "C", "multi\nline ✓"],
            "pIC50": ["7.50", "", "-3", "0.001"],
            "Blank": [None] * 4,
        }
    )
    column = t.column("Structure")
    t.fill_structures("Structure", t.oids[:2], parse_structures(["CCO", "C1CC"]))
    assert column.pending_mask().tolist() == [False, False, False, True]  # one left unparsed
    t.remove_rows([3])
    t.add_rows({"ID": ["E"]})
    t.set_sort([("pIC50", False)])
    return t


def test_round_trip_is_exact(table: Table, tmp_path: Path) -> None:
    path = tmp_path / "doc.csd"
    save_document(table, path, extras={"lists": {"hits": [1, 2]}})
    loaded, extras = load_document(path)
    assert snapshot(loaded) == snapshot(table)
    assert extras == {"lists": {"hits": [1, 2]}}
    assert loaded.next_oid == 6  # oid 3 was deleted and must never come back
    assert not loaded.history.can_undo and loaded.history.is_clean


def test_large_table_spans_several_chunks(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(document, "CHUNK_ROWS", 7)
    t = Table()
    t.add_rows({"n": [str(i) for i in range(50)], "t": [f"x{i % 3}" for i in range(50)]})
    save_document(t, tmp_path / "big.csd")
    assert snapshot(load_document(tmp_path / "big.csd")[0]) == snapshot(t)


def test_empty_table(tmp_path: Path) -> None:
    save_document(Table(), tmp_path / "e.csd")
    loaded, _ = load_document(tmp_path / "e.csd")
    assert len(loaded) == 0 and loaded.column_names == []


def test_failed_save_keeps_previous_version(table: Table, tmp_path: Path, monkeypatch) -> None:
    path = tmp_path / "doc.csd"
    save_document(table, path)
    before = path.read_bytes()

    def boom(column):
        raise OSError("disk full")

    monkeypatch.setattr(document, "_encode", boom)
    with pytest.raises(OSError):
        save_document(table, path)
    assert path.read_bytes() == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["doc.csd"]


def test_rejects_other_files(tmp_path: Path, table: Table) -> None:
    garbage = tmp_path / "x.csd"
    garbage.write_bytes(b"not sqlite at all" * 100)
    with pytest.raises(DocumentError):
        load_document(garbage)

    other = tmp_path / "other.csd"
    sqlite3.connect(other).execute("CREATE TABLE t (x)").connection.commit()
    with pytest.raises(DocumentError):
        load_document(other)

    newer = tmp_path / "newer.csd"
    save_document(table, newer)
    con = sqlite3.connect(newer)
    con.execute("UPDATE meta SET value = ? WHERE key = 'version'", (b"99",))
    con.commit()
    con.close()
    with pytest.raises(DocumentError, match="newer chemsuite"):
        load_document(newer)

    with pytest.raises(FileNotFoundError):
        load_document(tmp_path / "missing.csd")


CELL = st.one_of(
    st.none(), st.sampled_from(["", "1", "2.50", "007", "x", "CCO", "C1CC"]), st.text(max_size=5)
)


@settings(max_examples=60, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.data())
def test_random_tables_round_trip(tmp_path: Path, data) -> None:
    n = data.draw(st.integers(0, 12))
    names = data.draw(st.lists(st.text(min_size=1, max_size=4), unique=True, max_size=4))
    t = Table()
    t.add_rows({name: data.draw(st.lists(CELL, min_size=n, max_size=n)) for name in names}, count=n)
    if n:
        t.add_rows(
            {
                "S": StructureColumn.from_texts(
                    data.draw(
                        st.lists(st.sampled_from(["CCO", "C1CC", ""]), min_size=1, max_size=3)
                    )
                )
            }
        )
    path = tmp_path / "r.csd"
    save_document(t, path)
    assert snapshot(load_document(path)[0]) == snapshot(t)
