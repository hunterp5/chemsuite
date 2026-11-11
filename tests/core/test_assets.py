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
"""Assets stored with a document (proteins): undo, events, snapshots, and the file format."""

from __future__ import annotations

import sqlite3

import pytest

from chemsuite.core import document
from chemsuite.core.document import DocumentError, load_document, save_document
from chemsuite.core.table import Asset, AssetsChanged, Changed, PutAsset, Table


class Decoded:
    """A payload a package decoded (as the protein package keeps its structures)."""

    def __init__(self, text: str) -> None:
        self.text = text

    def to_bytes(self) -> bytes:
        return self.text.encode()

    @property
    def nbytes(self) -> int:
        return len(self.text)


def _protein(asset_id: str = "1ABC", payload: bytes | Decoded = b"atoms") -> Asset:
    return Asset(asset_id, "protein", f"{asset_id} kinase", payload, {"source": "test"})


def _table() -> Table:
    table = Table()
    table.add_rows({"ID": ["a", "b"]}, undoable=False)
    table.history.mark_clean()
    return table


def test_adding_replacing_and_removing_undo_as_one_history() -> None:
    table = _table()
    first, second = _protein(), _protein(payload=b"edited")
    table.edit("Open Protein", [PutAsset("1ABC", None, first)])
    table.edit("Edit Protein", [PutAsset("1ABC", first, second)])
    assert table.assets["1ABC"] is second and not table.history.is_clean
    assert table.undo() == "Edit Protein" and table.assets["1ABC"] is first
    assert table.undo() == "Open Protein" and "1ABC" not in table.assets
    assert table.history.is_clean
    table.redo()
    table.edit("Delete Protein", [PutAsset("1ABC", first, None)])
    assert dict(table.assets) == {}
    table.undo()
    assert table.assets["1ABC"] is first


def test_asset_changes_are_announced_without_disturbing_the_rows() -> None:
    table = _table()
    events = []
    table.subscribe(events.append)
    table.edit("Open Protein", [PutAsset("1ABC", None, _protein())])
    changed = [e for e in events if isinstance(e, Changed)]
    assert len(changed) == 1 and not changed[0].structural
    assert changed[0].changes == (AssetsChanged(("1ABC",)),)


def test_a_protein_and_its_ligand_rows_are_one_undo_step() -> None:
    table = _table()
    oids = table.add_rows(
        {"ID": ["ATP"], "Protein": ["1ABC"]},
        label="Open Protein",
        extra_ops=[PutAsset("1ABC", None, _protein())],
    )
    assert len(table) == 3 and "1ABC" in table.assets
    assert table.undo() == "Open Protein"
    assert len(table) == 2 and "1ABC" not in table.assets and "Protein" not in table
    table.redo()
    assert table.has_oid(int(oids[0])) and "1ABC" in table.assets


def test_assets_are_read_only_and_ids_stay_unique() -> None:
    table = _table()
    table.edit("Open", [PutAsset("1ABC", None, _protein())])
    with pytest.raises(TypeError):
        table.assets["X"] = _protein("X")  # type: ignore[index]
    assert table.unique_asset_id("1ABC") == "1ABC-2"
    assert table.unique_asset_id("2XYZ") == "2XYZ"


def test_a_snapshot_keeps_the_assets_of_the_moment() -> None:
    table = _table()
    table.edit("Open", [PutAsset("1ABC", None, _protein())])
    copy = table.snapshot()
    table.edit("Delete", [PutAsset("1ABC", table.assets["1ABC"], None)])
    assert "1ABC" in copy.assets and "1ABC" not in table.assets


def _version(path) -> int:
    with sqlite3.connect(path) as connection:
        (value,) = connection.execute("SELECT value FROM meta WHERE key = 'version'").fetchone()
    return int(value)


def test_documents_without_assets_stay_in_the_first_format(tmp_path) -> None:
    path = tmp_path / "plain.csd"
    save_document(_table(), path)
    assert _version(path) == 1
    table, _ = load_document(path)
    assert dict(table.assets) == {}


def test_assets_round_trip_in_order_with_their_details(tmp_path) -> None:
    table = _table()
    table.edit(
        "Open",
        [
            PutAsset("2XYZ", None, _protein("2XYZ", Decoded("decoded atoms"))),
            PutAsset("1ABC", None, _protein("1ABC", b"\x00\x01binary")),
        ],
    )
    path = tmp_path / "proteins.csd"
    save_document(table, path)
    assert _version(path) == 2
    again, _ = load_document(path)
    assert list(again.assets) == ["2XYZ", "1ABC"]
    first, second = again.assets.values()
    assert first.payload == b"decoded atoms" and second.payload == b"\x00\x01binary"
    assert (second.kind, second.name, dict(second.meta)) == (
        "protein",
        "1ABC kinase",
        {"source": "test"},
    )
    assert again.history.is_clean and not again.history.can_undo


def test_an_older_chemsuite_refuses_a_document_with_assets(tmp_path, monkeypatch) -> None:
    table = _table()
    table.edit("Open", [PutAsset("1ABC", None, _protein())])
    path = tmp_path / "proteins.csd"
    save_document(table, path)
    monkeypatch.setattr(document, "VERSION", 1)  # as a version from before assets
    with pytest.raises(DocumentError, match="newer chemsuite"):
        load_document(path)
