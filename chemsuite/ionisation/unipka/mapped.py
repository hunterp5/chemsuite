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

"""The network's parameters as one memory-mapped float32 file, shared by every process.

The shipped archive is compressed float64; each worker that loads it holds its own float32
copy (about 180 MB). Written once to a cache file and mapped read-only, the parameters are
pages of that file: the operating system keeps one copy for every worker and can drop it under
memory pressure, since it is backed by the file.

The cache is named after the archive's path, size and modification time, so a changed archive
gets a new cache. It is written to a temporary file and renamed into place, the manifest last;
one process writes it while the others wait on a lock file. When it cannot be written or read,
:func:`load` returns None and the caller loads the archive as before.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

from chemsuite.core.folders import derived_cache

FORMAT = 1
#: a lock older than this is a crashed writer's (seconds); writing takes a second or two
STALE_LOCK_S = 120.0
#: how long a process waits for another to write the cache before loading the archive itself
WAIT_S = 60.0
_POLL_S = 0.1
_ITEM = np.dtype(np.float32).itemsize

Params = dict[str, dict[str, np.ndarray]]  # fold ("f0") → parameter name → array


def _stem(source: Path) -> str:
    stat = source.stat()
    where = hashlib.sha256(str(source.resolve()).encode("utf-8")).hexdigest()[:12]
    return f"{source.stem}-{where}-{stat.st_size}-{stat.st_mtime_ns}-v{FORMAT}"


def _read(data: Path, manifest: Path) -> Params | None:
    """The mapped parameters, or None if the cache is missing or incomplete."""
    try:
        index = json.loads(manifest.read_text("utf-8"))
        if index.get("format") != FORMAT or data.stat().st_size != index["total"] * _ITEM:
            return None
        flat = np.memmap(data, dtype=np.float32, mode="r", shape=(index["total"],))
    except (OSError, ValueError, KeyError, TypeError):
        return None
    out: Params = {}
    for key, offset, shape in index["arrays"]:
        fold, name = key.split("/", 1)
        count = int(np.prod(shape, dtype=np.int64))
        # a plain ndarray view of the map: numpy treats it as any read-only array
        view = np.asarray(flat[offset : offset + count]).reshape(shape)
        out.setdefault(fold, {})[name] = view
    return out


def _write(source: Path, data: Path, manifest: Path) -> None:
    """Convert the archive's fold parameters to one float32 file, then its manifest."""
    arrays, offset = [], 0
    tmp = data.with_name(f"{data.name}.{os.getpid()}.tmp")
    try:
        with np.load(source) as archive, open(tmp, "wb") as out:
            for key in archive.files:
                if "/" not in key or not key.startswith("f"):
                    continue
                values = np.ascontiguousarray(archive[key], dtype=np.float32)
                out.write(values.tobytes())
                arrays.append([key, offset, list(values.shape)])
                offset += values.size
        os.replace(tmp, data)
    finally:
        with contextlib.suppress(OSError):
            tmp.unlink()
    text = json.dumps({"format": FORMAT, "total": offset, "arrays": arrays})
    tmp = manifest.with_name(f"{manifest.name}.{os.getpid()}.tmp")
    tmp.write_text(text, "utf-8")
    os.replace(tmp, manifest)


def _take(lock: Path) -> bool:
    """Create *lock* (True), or False if another process holds it. A stale lock is taken over."""
    try:
        os.close(os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        with contextlib.suppress(OSError):
            if time.time() - lock.stat().st_mtime > STALE_LOCK_S:
                lock.unlink()
                return _take(lock)
        return False


def load(source: Path, folder: Path | None = None) -> Params | None:
    """The fold parameters of the archive *source*, mapped from the cache in *folder*
    (default: :func:`~chemsuite.core.folders.derived_cache`), written first if needed."""
    try:
        folder = folder or derived_cache()
        stem = _stem(source)
        folder.mkdir(parents=True, exist_ok=True)
    except OSError:
        return None
    data, manifest, lock = (folder / f"{stem}{ext}" for ext in (".f32", ".json", ".lock"))
    found = _read(data, manifest)
    if found is not None:
        return found
    deadline = time.monotonic() + WAIT_S
    while True:
        if _take(lock):
            try:
                found = _read(data, manifest)  # written while this process waited
                if found is None:
                    _write(source, data, manifest)
                    found = _read(data, manifest)
            except (OSError, ValueError, KeyError):
                found = None
            finally:
                with contextlib.suppress(OSError):
                    lock.unlink()
            return found
        found = _read(data, manifest)
        if found is not None or time.monotonic() > deadline:
            return found
        time.sleep(_POLL_S)
