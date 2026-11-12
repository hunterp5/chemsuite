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
"""The SQLite container under every chemsuite file, and the .csd staying as it was."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from chemsuite.core import container
from chemsuite.core.assets import Asset
from chemsuite.core.container import ContainerError
from chemsuite.core.document import DocumentError, load_document, save_document

#: saved by the document code as it was before the container was shared (Phase B)
PHASE_B = Path(__file__).parent / "data" / "phase_b.csd"


def _dump(path: Path) -> dict[str, list]:
    connection = sqlite3.connect(path)
    try:
        out = {"schema": connection.execute("SELECT type, name, sql FROM sqlite_master").fetchall()}
        meta = connection.execute("SELECT key, value FROM meta ORDER BY rowid").fetchall()
        out["meta"] = [(k, v) for k, v in meta if k != "app_version"]
        for name in ("columns", "chunks", "assets"):
            out[name] = connection.execute(f"SELECT * FROM {name} ORDER BY rowid").fetchall()
        return out
    finally:
        connection.close()


def test_a_document_saved_before_the_shared_container_saves_back_the_same(tmp_path) -> None:
    table, extras = load_document(PHASE_B)
    assert list(table.assets) == ["1ABC"] and extras["lists"] == {"x": 1}
    again = tmp_path / "again.csd"
    save_document(table, again, extras=extras)
    assert _dump(again) == _dump(PHASE_B)


def test_saving_is_atomic(tmp_path) -> None:
    target = tmp_path / "file.db"
    target.write_bytes(b"before")

    def fail(connection: sqlite3.Connection) -> None:
        connection.executescript(container.META)
        raise RuntimeError("disk full")

    with pytest.raises(RuntimeError):
        container.atomic_save(target, fail)
    assert target.read_bytes() == b"before"
    assert [p.name for p in tmp_path.iterdir()] == ["file.db"]  # no temporary file left


def _write(path: Path, version: int, assets=()) -> None:
    def write(connection: sqlite3.Connection) -> None:
        connection.executescript(container.META)
        meta = {"format": b"test-format", "version": str(version).encode()}
        container.write_meta(connection, meta)
        container.write_assets(connection, assets)

    container.atomic_save(path, write)


def _read(path: Path, newest: int = 2) -> tuple[int, list[Asset]]:
    def read(connection: sqlite3.Connection) -> tuple[int, list[Asset]]:
        meta = container.read_meta(connection)
        version = container.check_version(meta, "test-format", newest, what="test file")
        return version, container.read_assets(connection)

    return container.read_container(path, read, what="test file")


def test_older_versions_open_and_newer_ones_are_refused(tmp_path) -> None:
    path = tmp_path / "f.db"
    _write(path, 1)
    assert _read(path) == (1, [])
    _write(path, 3)
    with pytest.raises(ContainerError, match="newer chemsuite"):
        _read(path)


def test_assets_round_trip(tmp_path) -> None:
    path = tmp_path / "f.db"
    assets = [Asset("a", "protein", "A", b"x" * 100, {"source": "s"}), Asset("b", "k", "B", b"")]
    _write(path, 2, assets)
    _version, read = _read(path)
    assert [(a.id, a.kind, a.name, a.data(), dict(a.meta)) for a in read] == [
        ("a", "protein", "A", b"x" * 100, {"source": "s"}),
        ("b", "k", "B", b"", {}),
    ]


def test_other_files_are_refused_with_the_callers_error(tmp_path) -> None:
    path = tmp_path / "f.csd"
    path.write_bytes(b"not sqlite at all" * 100)
    with pytest.raises(DocumentError, match="not a chemsuite document"):
        load_document(path)
    _write(path, 1)
    with pytest.raises(DocumentError, match="not a chemsuite document"):
        load_document(path)
