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
"""Downloading an entry from the Protein Data Bank by its ID.

The mmCIF file comes from RCSB's public file server (``files.rcsb.org``) through
:mod:`chemsuite.core.download`, so it is HTTPS only, audited, and refused when IT policy
turns the network off. Only the ID is sent.

IDs are the classic four characters (``1ATP``: a digit, then three letters or digits) or the
extended form (``pdb_00001atp``).
"""

from __future__ import annotations

import gzip
import re

from chemsuite.core.download import DownloadError, Opener, fetch_bytes
from chemsuite.core.jobs import CancelToken

from . import io
from .model import Structure

URL = "https://files.rcsb.org/download/{id}.cif.gz"
_CLASSIC = re.compile(r"^[1-9][A-Za-z0-9]{3}$")
_EXTENDED = re.compile(r"^pdb_[0-9]{4}[A-Za-z0-9]{4}$", re.IGNORECASE)


def normalise(text: str) -> str | None:
    """The ID as the archive writes it (``1ATP``, ``pdb_00001atp``), or None if it is not one."""
    text = text.strip()
    if _CLASSIC.match(text):
        return text.upper()
    if _EXTENDED.match(text):
        return text.lower()
    return None


def problem(text: str) -> str | None:
    """Why ``text`` is not a PDB ID (None when it is)."""
    if not text.strip():
        return "Enter a PDB ID"
    if normalise(text) is None:
        return f"“{text.strip()}” is not a PDB ID (four characters, such as 1ATP)"
    return None


def fetch(pdb_id: str, *, opener: Opener | None = None, cancel: CancelToken | None = None):
    """Download and read entry ``pdb_id``. Raises DownloadError (a ProteinFileError when the
    file cannot be read)."""
    code = normalise(pdb_id)
    if code is None:
        raise DownloadError(problem(pdb_id) or "not a PDB ID")
    data = fetch_bytes(URL.format(id=code), purpose=f"PDB entry {code}", opener=opener,
                       cancel=cancel)  # fmt: skip
    try:
        text = gzip.decompress(data).decode("utf-8", errors="replace")
    except (OSError, EOFError) as err:
        raise DownloadError(f"{code} did not come back as a gzipped mmCIF file") from err
    structure: Structure = io.parse(text, io.MMCIF)
    return structure


def name_of(structure: Structure, fallback: str) -> str:
    """The name to show a structure under: its ID, else ``fallback``."""
    return str(structure.meta.get("id") or fallback).upper()
