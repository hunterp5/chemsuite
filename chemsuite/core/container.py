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
"""The SQLite container under every chemsuite file (``.csd``, and the protein workspace's).

A container is one SQLite file with a ``meta`` table of key → bytes, saying at least which
``format`` it is and its ``version``, and optionally an ``assets`` table (see
:mod:`chemsuite.core.assets`). Each format adds tables of its own: the document's column chunks
(:mod:`chemsuite.core.document`).

- **Saving** writes a complete new file next to the target and swaps it into place, so a crash or
  a full disk never damages the previous version.
- **Opening** refuses a file of another format, and one of a *newer* version than the reader
  knows (it would lose what it cannot read on the next save); any older version is read.
"""

from __future__ import annotations

import json
import os
import sqlite3
import struct
import tempfile
import zlib
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import TypeVar

from .assets import Asset
from .jobs import CancelToken, check

T = TypeVar("T")

META = "CREATE TABLE meta (key TEXT PRIMARY KEY, value BLOB NOT NULL);\n"
ASSETS = """
CREATE TABLE assets (position INTEGER PRIMARY KEY, id TEXT NOT NULL UNIQUE, kind TEXT NOT NULL,
                     name TEXT NOT NULL, meta TEXT NOT NULL, data BLOB NOT NULL);
"""
#: Level 1 is ~5x faster than the default and only slightly larger on chemistry data.
ZLIB_LEVEL = 1


class ContainerError(ValueError):
    """The file is not of the expected kind, is damaged, or needs a newer chemsuite."""


# -- saving --------------------------------------------------------------------------------------


def atomic_save(path: str | Path, write: Callable[[sqlite3.Connection], None]) -> None:
    """Make a new container at *path*: *write* fills a fresh database, which then replaces
    the file in one step (nothing is left behind if *write* fails or is cancelled)."""
    path = Path(path)
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        connection = sqlite3.connect(tmp)
        try:
            write(connection)
            connection.commit()
        finally:
            connection.close()
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_meta(connection: sqlite3.Connection, meta: Mapping[str, bytes]) -> None:
    """Fill the ``meta`` table, in *meta*'s order (the table must exist: :data:`META`)."""
    connection.executemany("INSERT INTO meta VALUES (?, ?)", meta.items())


def write_assets(
    connection: sqlite3.Connection,
    assets: Iterable[Asset],
    *,
    cancel: CancelToken | None = None,
) -> None:
    """Add the ``assets`` table holding *assets*, compressed (nothing when there are none)."""
    for position, asset in enumerate(assets):
        check(cancel)
        if position == 0:
            connection.executescript(ASSETS)
        connection.execute(
            "INSERT INTO assets VALUES (?, ?, ?, ?, ?, ?)",
            (
                position,
                asset.id,
                asset.kind,
                asset.name,
                json.dumps(dict(asset.meta)),
                zlib.compress(asset.data(), ZLIB_LEVEL),
            ),
        )


# -- opening -------------------------------------------------------------------------------------


def read_container(
    path: str | Path,
    read: Callable[[sqlite3.Connection], T],
    *,
    error: type[ContainerError] = ContainerError,
    what: str = "chemsuite file",
) -> T:
    """Open the container at *path* and return what *read* makes of it. SQLite and decoding
    failures become *error* ("not a …", "damaged"); *read*'s own *error* passes through."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    connection = sqlite3.connect(path)  # read only by convention; never written here
    try:
        return read(connection)
    except sqlite3.DatabaseError as err:
        raise error(f"{path.name} is not a {what} ({err})") from err
    except (KeyError, ValueError, TypeError, zlib.error, struct.error) as err:
        if isinstance(err, ContainerError):
            raise
        raise error(f"{path.name} is damaged ({err!r})") from err
    finally:
        connection.close()


def read_meta(connection: sqlite3.Connection) -> dict[str, bytes]:
    return dict(connection.execute("SELECT key, value FROM meta"))


def check_version(
    meta: Mapping[str, bytes],
    form: str,
    newest: int,
    *,
    error: type[ContainerError] = ContainerError,
    what: str = "chemsuite file",
) -> int:
    """The container's version, after checking it is *form* and not newer than *newest*."""
    if meta.get("format") != form.encode():
        raise error(f"not a {what}")
    version = int(meta["version"])
    if version > newest:
        raise error(f"this {what} was saved by a newer chemsuite (format {version}); please update")
    return version


def has_assets(connection: sqlite3.Connection) -> bool:
    found = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'assets'"
    )
    return found.fetchone() is not None


def read_assets(
    connection: sqlite3.Connection, *, error: type[ContainerError] = ContainerError
) -> list[Asset]:
    """The assets, in the order they were saved (none when the file has no ``assets`` table)."""
    if not has_assets(connection):
        return []
    rows = connection.execute("SELECT id, kind, name, meta, data FROM assets ORDER BY position")
    assets = []
    for asset_id, kind, name, meta, data in rows:
        facts = json.loads(meta)
        if not isinstance(facts, dict):
            raise error(f"damaged file: asset {asset_id!r}")
        facts = {str(k): str(v) for k, v in facts.items()}
        assets.append(Asset(asset_id, kind, name, zlib.decompress(data), facts))
    return assets
