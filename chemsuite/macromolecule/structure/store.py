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
"""A structure as bytes, for a document asset: its arrays, exactly.

The payload is an uncompressed ``.npz`` (the document compresses assets) of the structure's
fields, plus a format number and the metadata as JSON. Reading never unpickles: only plain
numeric and fixed-width text arrays are accepted. Packing keeps everything a file would lose:
atom uids, bond orders, every alternate location and model.
"""

from __future__ import annotations

import io
import json
import zipfile
from dataclasses import fields
from typing import Any

import numpy as np

from .model import Structure

FORMAT = 1
_FORMAT_KEY = "_format"
_META_KEY = "_meta"
_SCALARS = ("next_uid",)


class StoreError(ValueError):
    """The payload is not a structure this version can read."""


def pack(structure: Structure) -> bytes:
    arrays: dict[str, np.ndarray] = {
        _FORMAT_KEY: np.array(FORMAT, dtype=np.int32),
        _META_KEY: np.frombuffer(json.dumps(dict(structure.meta)).encode(), dtype=np.uint8),
    }
    for f in fields(structure):
        value = structure.__dict__[f.name]
        if isinstance(value, np.ndarray):
            arrays[f.name] = value
        elif f.name in _SCALARS:
            arrays[f.name] = np.array(value, dtype=np.int64)
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    return buffer.getvalue()


def unpack(data: bytes) -> Structure:
    try:
        with np.load(io.BytesIO(data), allow_pickle=False) as npz:
            arrays = {key: npz[key] for key in npz.files}
    except (ValueError, OSError, zipfile.BadZipFile) as error:
        raise StoreError(f"not a stored structure: {error}") from error
    version = int(arrays.pop(_FORMAT_KEY, np.array(0)))
    if version > FORMAT:
        raise StoreError(f"stored by a newer chemsuite (format {version}; this one reads {FORMAT})")
    if version < 1:
        raise StoreError("not a stored structure: it has no format number")
    values: dict[str, Any] = {}
    for f in fields(Structure):
        if f.name == "meta":
            values["meta"] = json.loads(arrays[_META_KEY].tobytes().decode())
        elif f.name in _SCALARS:
            values[f.name] = int(arrays[f.name])
        elif f.name not in arrays:
            raise StoreError(f"stored structure has no {f.name}")
        else:
            values[f.name] = arrays[f.name]
    return Structure(**values)


def as_structure(payload) -> Structure:
    """An asset's payload as a structure: a Structure as it is, stored bytes unpacked."""
    return payload if isinstance(payload, Structure) else unpack(payload)
