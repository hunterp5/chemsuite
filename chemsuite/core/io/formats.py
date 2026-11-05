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

"""File-type detection, transparent gzip, text-encoding detection and atomic writes."""

from __future__ import annotations

import gzip
import io
import os
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import IO

from chemsuite.core.table import Column

#: One chunk of rows read from a file: column name -> values (or a ready-made Column).
Batch = Mapping[str, Sequence[object] | Column]

SDF = "sdf"
CSV = "csv"
TSV = "tsv"
SMILES = "smiles"
XLSX = "xlsx"

_SUFFIXES = {
    ".sdf": SDF,
    ".sd": SDF,
    ".mol": SDF,
    ".csv": CSV,
    ".tsv": TSV,
    ".tab": TSV,
    ".txt": TSV,
    ".smi": SMILES,
    ".smiles": SMILES,
    ".xlsx": XLSX,
}

#: Header names (case-insensitive) that mark a column of structures in tabular files.
STRUCTURE_HEADERS = frozenset(
    {"smiles", "smi", "structure", "mol", "molecule", "canonical_smiles", "isomeric_smiles"}
)


class UnsupportedFormat(ValueError):
    pass


def split_suffix(path: str | Path) -> tuple[str, bool]:
    """``(format, gzipped)`` for *path*: ``compounds.sdf.gz`` -> ``("sdf", True)``."""
    suffixes = [s.lower() for s in Path(path).suffixes]
    gzipped = bool(suffixes) and suffixes[-1] == ".gz"
    if gzipped:
        suffixes.pop()
    fmt = _SUFFIXES.get(suffixes[-1] if suffixes else "")
    if fmt is None:
        raise UnsupportedFormat(f"unsupported file type: {Path(path).name}")
    if gzipped and fmt == XLSX:
        raise UnsupportedFormat("compressed Excel files are not supported")
    return fmt, gzipped


def guess_structure_column(headers: Sequence[str]) -> str | None:
    """The first header that names a structure column, else None (no structures)."""
    for header in headers:
        name = header.strip().lower()
        if name in STRUCTURE_HEADERS or name.endswith("smiles"):
            return header
    return None


def unique_headers(headers: Sequence[str]) -> list[str]:
    """Non-empty, unique column names: blanks become ``Column N``, repeats get `` (2)``."""
    out: list[str] = []
    seen: set[str] = set()
    for i, header in enumerate(headers, start=1):
        base = header.strip() or f"Column {i}"
        name, n = base, 2
        while name in seen:
            name, n = f"{base} ({n})", n + 1
        seen.add(name)
        out.append(name)
    return out


@contextmanager
def open_binary(path: str | Path) -> Iterator[tuple[IO[bytes], IO[bytes]]]:
    """``(stream, raw)``: *stream* yields file content (decompressed for ``.gz``); *raw* is the
    file on disk, whose ``tell()`` measures progress against the file size."""
    _, gzipped = split_suffix(path)
    with open(path, "rb") as raw:
        if gzipped:
            with gzip.GzipFile(fileobj=raw) as stream:
                yield stream, raw
        else:
            yield raw, raw


def detect_encoding(sample: bytes) -> str:
    """UTF-8 (with or without BOM) when the sample decodes as UTF-8, else Windows-1252.

    Excel on Windows still writes CSV in the local code page; reading it as UTF-8 would turn
    every µ, ° or Å into an error.
    """
    if sample.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    try:
        sample.decode("utf-8")
    except UnicodeDecodeError as err:
        if err.start < len(sample) - 3:  # not just a character cut off at the sample's end
            return "cp1252"
    return "utf-8"


@contextmanager
def open_text(path: str | Path) -> Iterator[tuple[IO[str], IO[bytes]]]:
    """``(text stream, raw file)`` with the encoding detected from the first 64 KB."""
    with open_binary(path) as (stream, raw):
        buffered = io.BufferedReader(stream) if not hasattr(stream, "peek") else stream
        encoding = detect_encoding(buffered.peek(65536)[:65536])
        text = io.TextIOWrapper(buffered, encoding=encoding, errors="replace", newline="")
        try:
            yield text, raw
        finally:
            text.detach()


@contextmanager
def atomic_write(path: str | Path, mode: str = "w", **kwargs) -> Iterator[IO]:
    """Write to a temporary file next to *path* and move it into place only on success.

    A crash or an exception mid-write never leaves a truncated file where the user's file was.
    """
    path = Path(path)
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, mode, **kwargs) as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def file_size(path: str | Path) -> int:
    return os.path.getsize(path)
