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

"""Downloading AiZynthFinder's public data and writing a ``config.yml`` for it.

The files are the ones AiZynthFinder's own ``download_public_data`` fetches: the USPTO
expansion policy and its templates, the ring-breaker policy and its templates, the USPTO filter
policy (all on Zenodo) and the ZINC stock (figshare). The stock alone is about 1.3 GB, so:

- each file is written to ``<name>.part`` and renamed only when complete, so a finished file is
  never a truncated one, and a file already in the folder is not fetched again;
- an interrupted ``.part`` resumes with an HTTP range request where the server allows it;
- a cancel stops within one chunk and keeps the ``.part`` for next time.

Downloads go through :mod:`chemsuite.core.download` (HTTPS only, the audit log, the IT
policy's ``network_allowed``). The *opener* is injectable so tests run without a network.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from chemsuite.core.download import DownloadError, Opener, download_file
from chemsuite.core.jobs import CancelToken, Progress, check

log = logging.getLogger(__name__)

CHUNK = 1 << 20  # 1 MiB
TIMEOUT_S = 60.0
CONFIG_NAME = "config.yml"
#: first line of a config chemsuite wrote, so it may be rewritten but a user's own is not
MARKER = "# Written by chemsuite for AiZynthFinder's public data."


@dataclass(frozen=True)
class DataFile:
    name: str
    url: str


#: the same files and sources as ``aizynthfinder.tools.download_public_data`` (4.4.x)
FILES: tuple[DataFile, ...] = (
    DataFile("uspto_model.onnx", "https://zenodo.org/record/7797465/files/uspto_model.onnx"),
    DataFile(
        "uspto_templates.csv.gz",
        "https://zenodo.org/record/7341155/files/uspto_unique_templates.csv.gz",
    ),
    DataFile(
        "uspto_ringbreaker_model.onnx",
        "https://zenodo.org/record/7797465/files/uspto_ringbreaker_model.onnx",
    ),
    DataFile(
        "uspto_ringbreaker_templates.csv.gz",
        "https://zenodo.org/record/7341155/files/uspto_ringbreaker_unique_templates.csv.gz",
    ),
    DataFile(
        "uspto_filter_model.onnx",
        "https://zenodo.org/record/7797465/files/uspto_filter_model.onnx",
    ),
    DataFile("zinc_stock.hdf5", "https://ndownloader.figshare.com/files/23086469"),
)

#: rough total, for progress before the servers have said how big each file is
APPROX_BYTES = 1_466_000_000

PURPOSE = "AiZynthFinder public data"


def missing(folder: Path) -> list[DataFile]:
    """The files not yet complete in *folder*."""
    return [item for item in FILES if not (folder / item.name).is_file()]


def download_public_data(
    folder: Path | str,
    *,
    opener: Opener | None = None,
    cancel: CancelToken | None = None,
    progress: Progress | None = None,
) -> Path:
    """Fetch every missing file into *folder*, write the config, and return its path.

    *progress* is told ``(done, total)`` in MB over all files still to fetch.
    """
    folder = Path(folder).expanduser()
    try:
        folder.mkdir(parents=True, exist_ok=True)
    except OSError as err:
        raise DownloadError(f"Could not create {folder}: {err}") from err
    open_url = opener
    todo = missing(folder)
    done_bytes = 0
    total_bytes = APPROX_BYTES if todo else 0

    def report(extra: int) -> None:
        if progress is not None:
            progress((done_bytes + extra) >> 20, max(total_bytes, done_bytes + extra) >> 20)

    for item in todo:
        check(cancel)
        written = _fetch(item, folder, open_url, cancel, report)
        done_bytes += written
        report(0)
    return write_config(folder)


def _fetch(
    item: DataFile,
    folder: Path,
    open_url: Opener | None,
    cancel: CancelToken | None,
    report: Callable[[int], None],
) -> int:
    """Download one file into *folder* and return the bytes it has."""
    return download_file(item.url, folder / item.name, purpose=PURPOSE, opener=open_url,
                         cancel=cancel, report=report)  # fmt: skip


def write_config(folder: Path) -> Path:
    """Write the ``config.yml`` naming *folder*'s files, and return its path.

    A ``config.yml`` already in the folder that chemsuite did not write is left alone, and the
    config goes to ``chemsuite_config.yml`` beside it instead.
    """
    folder = Path(folder).expanduser().resolve()
    path = folder / CONFIG_NAME
    if path.is_file():
        try:
            ours = path.read_text(encoding="utf-8").startswith(MARKER)
        except OSError:
            ours = False
        if not ours:
            path = folder / "chemsuite_config.yml"

    def at(name: str) -> str:
        return json.dumps(str(folder / name))  # a double-quoted YAML scalar, any path is safe

    body = "\n".join(
        [
            MARKER,
            "expansion:",
            "  uspto:",
            f"    - {at('uspto_model.onnx')}",
            f"    - {at('uspto_templates.csv.gz')}",
            "  ringbreaker:",
            f"    - {at('uspto_ringbreaker_model.onnx')}",
            f"    - {at('uspto_ringbreaker_templates.csv.gz')}",
            "filter:",
            f"  uspto: {at('uspto_filter_model.onnx')}",
            "stock:",
            f"  zinc: {at('zinc_stock.hdf5')}",
            "",
        ]
    )
    try:
        path.write_text(body, encoding="utf-8")
    except OSError as err:
        raise DownloadError(f"Could not write {path}: {err}") from err
    return path
