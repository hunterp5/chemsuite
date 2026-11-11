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

"""The chemsuite document (``.csd``): one SQLite file holding a table exactly as it was.

Everything round-trips losslessly: oids (and the next oid, so oids are never reused), column
order and kinds, each number's display format, blank vs. invalid vs. unparsed structures, the
sort order and each column's style (precision, log10 marker, colours). Columns are stored in
compressed chunks of 64k rows, so saving and opening stream instead of building one huge blob,
and a later version can read single columns lazily.

Saving writes a complete new file next to the target and swaps it into place, so a crash or a
full disk never damages the previous version. ``extras`` is a JSON object for state owned by
the application (filters, named lists, plots); the core stores it without interpreting it.

The table's assets (proteins, see :mod:`chemsuite.core.table.assets`) go in their own table,
compressed. A document with assets is format 2 and one without is still format 1, so documents
without proteins keep opening in older versions, and an older version refuses a document with
proteins (rather than open it and lose them on the next save).
"""

from __future__ import annotations

import json
import os
import sqlite3
import struct
import tempfile
import zlib
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

import chemsuite
from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.table import NUMBER, TEXT, Asset, ColumnStyle, Table
from chemsuite.core.table.columns import (
    REACTION,
    STRUCTURE,
    Column,
    NumberColumn,
    ParsedColumn,
    ReactionColumn,
    StructureColumn,
    TextColumn,
    blank_column,
    concat_columns,
    object_array,
)

FORMAT = "chemsuite-document"
VERSION = 2  # the newest format this version reads
PLAIN_VERSION = 1  # written when a document has no assets
SUFFIX = ".csd"
CHUNK_ROWS = 65_536
#: Level 1 is ~5x faster than the default and only slightly larger on chemistry data.
_ZLIB_LEVEL = 1

_SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value BLOB NOT NULL);
CREATE TABLE columns (position INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE,
                      kind TEXT NOT NULL);
CREATE TABLE chunks (position INTEGER NOT NULL, chunk INTEGER NOT NULL, data BLOB NOT NULL,
                     PRIMARY KEY (position, chunk));
"""
_ASSETS = """
CREATE TABLE assets (position INTEGER PRIMARY KEY, id TEXT NOT NULL UNIQUE, kind TEXT NOT NULL,
                     name TEXT NOT NULL, meta TEXT NOT NULL, data BLOB NOT NULL);
"""


class DocumentError(ValueError):
    """The file is not a chemsuite document, is damaged, or needs a newer chemsuite."""


# -- public API ---------------------------------------------------------------------------------


def save_document(
    table: Table,
    path: str | Path,
    *,
    extras: Mapping[str, Any] | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    path = Path(path)
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        connection = sqlite3.connect(tmp)
        try:
            _write(connection, table, dict(extras or {}), progress, cancel)
            connection.commit()
        finally:
            connection.close()
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def load_document(
    path: str | Path,
    *,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[Table, dict[str, Any]]:
    """Return ``(table, extras)``. Raises :class:`DocumentError` for unreadable files."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    connection = sqlite3.connect(path)  # read only by convention; never written here
    try:
        return _read(connection, progress, cancel)
    except sqlite3.DatabaseError as err:
        raise DocumentError(f"{path.name} is not a chemsuite document ({err})") from err
    except (KeyError, ValueError, TypeError, zlib.error, struct.error) as err:
        if isinstance(err, DocumentError):
            raise
        raise DocumentError(f"{path.name} is damaged ({err!r})") from err
    finally:
        connection.close()


# -- writing ------------------------------------------------------------------------------------


def _write(connection, table, extras, progress, cancel) -> None:
    connection.executescript(_SCHEMA)
    assets = list(table.assets.values())
    meta = {
        "format": FORMAT.encode(),
        "version": str(VERSION if assets else PLAIN_VERSION).encode(),
        "app_version": chemsuite.__version__.encode(),
        "rows": str(len(table)).encode(),
        "next_oid": str(table.next_oid).encode(),
        "oids": zlib.compress(table.oids.astype("<i8").tobytes()),
        "sort": json.dumps(table.sort_keys).encode(),
        "styles": json.dumps(
            {
                name: table.column(name).style.to_json()
                for name in table.column_names
                if not table.column(name).style.is_default
            }
        ).encode(),
        "extras": json.dumps(extras).encode(),
    }
    connection.executemany("INSERT INTO meta VALUES (?, ?)", meta.items())
    names = table.column_names
    for position, name in enumerate(names):
        check(cancel)
        column = table.column(name)
        connection.execute("INSERT INTO columns VALUES (?, ?, ?)", (position, name, column.kind))
        for chunk_index, start in enumerate(range(0, len(column), CHUNK_ROWS)):
            part = column.take(np.arange(start, min(start + CHUNK_ROWS, len(column))))
            connection.execute(
                "INSERT INTO chunks VALUES (?, ?, ?)",
                (position, chunk_index, zlib.compress(_encode(part), _ZLIB_LEVEL)),
            )
        if progress is not None:
            progress(position + 1, len(names))
    if assets:
        connection.executescript(_ASSETS)
    for position, asset in enumerate(assets):
        check(cancel)
        connection.execute(
            "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?)",
            (
                position,
                asset.id,
                asset.kind,
                asset.name,
                json.dumps(dict(asset.meta)),
                zlib.compress(asset.data(), _ZLIB_LEVEL),
            ),
        )


def _encode(column: Column) -> bytes:
    if isinstance(column, NumberColumn):
        return _pack(
            [
                column._values.astype("<f8").tobytes(),
                np.packbits(column._valid).tobytes(),
                column._decimals.tobytes(),
            ]
        )
    if isinstance(column, TextColumn):
        compact = column.compacted()
        return _pack([compact._codes.astype("<u4").tobytes(), _json(compact._dict.strings)])
    if isinstance(column, ParsedColumn):
        blobs = column._blobs.tolist()
        lengths = np.array([-1 if b is None else len(b) for b in blobs], dtype="<i4")
        return _pack(
            [lengths.tobytes(), b"".join(b for b in blobs if b), _json(column._texts.tolist())]
        )
    raise TypeError(f"cannot store a {type(column).__name__}")


def _json(values: list[str]) -> bytes:
    return json.dumps(values, ensure_ascii=False, separators=(",", ":")).encode("utf-8")


def _pack(sections: list[bytes]) -> bytes:
    header = struct.pack(f"<I{len(sections)}Q", len(sections), *map(len, sections))
    return header + b"".join(sections)


def _unpack(data: bytes) -> list[bytes]:
    (count,) = struct.unpack_from("<I", data)
    lengths = struct.unpack_from(f"<{count}Q", data, 4)
    offset = 4 + 8 * count
    out = []
    for length in lengths:
        out.append(data[offset : offset + length])
        offset += length
    if offset != len(data):
        raise DocumentError("damaged column data")
    return out


# -- reading ------------------------------------------------------------------------------------


def _read(connection, progress, cancel) -> tuple[Table, dict[str, Any]]:
    meta = dict(connection.execute("SELECT key, value FROM meta"))
    if meta.get("format") != FORMAT.encode():
        raise DocumentError("not a chemsuite document")
    version = int(meta["version"])
    if version > VERSION:
        raise DocumentError(
            f"this document was saved by a newer chemsuite (format {version}); please update"
        )
    rows = int(meta["rows"])
    oids = np.frombuffer(zlib.decompress(meta["oids"]), dtype="<i8").astype(np.int64)
    if len(oids) != rows:
        raise DocumentError("damaged document: row count does not match")
    specs = connection.execute("SELECT position, name, kind FROM columns ORDER BY position")
    specs = specs.fetchall()
    columns: dict[str, Column] = {}
    for done, (position, name, kind) in enumerate(specs, start=1):
        check(cancel)
        chunks = connection.execute(
            "SELECT data FROM chunks WHERE position = ? ORDER BY chunk", (position,)
        )
        parts = [_decode(kind, _unpack(zlib.decompress(data))) for (data,) in chunks]
        column = concat_columns(parts) if parts else blank_column(kind, 0)
        if len(column) != rows:
            raise DocumentError(f"damaged document: column {name!r} has {len(column)} rows")
        columns[name] = column
        if progress is not None:
            progress(done, len(specs))
    for name, style in json.loads(meta.get("styles", b"{}")).items():
        if name in columns:
            columns[name].style = ColumnStyle.from_json(style).for_kind(columns[name].kind)
    table = Table.from_columns(columns, oids=oids, next_oid=int(meta["next_oid"]))
    if version >= 2:
        table.set_assets(_read_assets(connection))
    sort = [(name, bool(asc)) for name, asc in json.loads(meta.get("sort", b"[]"))]
    if sort:
        table.set_sort(sort)
    table.history.clear()
    table.history.mark_clean()
    return table, json.loads(meta.get("extras", b"{}"))


def _read_assets(connection) -> list[Asset]:
    rows = connection.execute("SELECT id, kind, name, meta, data FROM assets ORDER BY position")
    assets = []
    for asset_id, kind, name, meta, data in rows:
        facts = json.loads(meta)
        if not isinstance(facts, dict):
            raise DocumentError(f"damaged document: asset {asset_id!r}")
        facts = {str(k): str(v) for k, v in facts.items()}
        assets.append(Asset(asset_id, kind, name, zlib.decompress(data), facts))
    return assets


def _decode(kind: str, sections: list[bytes]) -> Column:
    if kind == NUMBER:
        values_raw, valid_raw, decimals_raw = sections
        n = len(decimals_raw)
        values = np.frombuffer(values_raw, dtype="<f8").astype(np.float64)
        valid = np.unpackbits(np.frombuffer(valid_raw, dtype=np.uint8), count=n).astype(bool)
        return NumberColumn(values, valid, np.frombuffer(decimals_raw, dtype=np.int8).copy())
    if kind == TEXT:
        codes_raw, strings_raw = sections
        column = TextColumn.from_texts(json.loads(strings_raw))
        lookup = column._codes  # strings were stored unique and in code order
        codes = np.frombuffer(codes_raw, dtype="<u4").astype(np.uint32)
        return TextColumn(lookup[codes], column._dict)
    if kind in (STRUCTURE, REACTION):
        lengths_raw, blob_data, texts_raw = sections
        lengths = np.frombuffer(lengths_raw, dtype="<i4")
        blobs: list[bytes | None] = []
        offset = 0
        for length in lengths.tolist():
            if length < 0:
                blobs.append(None)
            else:
                blobs.append(blob_data[offset : offset + length])
                offset += length
        cls = StructureColumn if kind == STRUCTURE else ReactionColumn
        return cls(object_array(blobs), object_array(json.loads(texts_raw)))
    raise DocumentError(f"unknown column kind {kind!r}")
