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

"""The permeability models chemsuite ships, in ``data/``, pinned by size and SHA-256.

Each file is a ``.csqm`` archive written by ``scripts/permeability/train.py``: a plain ``.npz``
read with ``allow_pickle=False``, so **opening one never runs code**. The pins below are what
``--ship`` printed when the model was built, and ``tests/permeability/test_models.py`` checks
the shipped copies against them — the arrangement
:mod:`chemsuite.ionisation.unipka.weights` uses for the pKa weights.

The models are GPL-3.0-or-later, as chemsuite is. What they were fitted on is public data whose
terms are recorded in ``data/NOTICE-chemsuite-permeability.txt`` beside the licence of
Biogen's MIT-licensed set; a ``.csqm`` carries its own training keys and values, so that NOTICE
travels inside the model as well as beside it.

One :class:`~chemsuite.qsar.model.Model` is kept per file per process: the archives are a few
megabytes and a run over a table would otherwise read one per chunk.
"""

from __future__ import annotations

import hashlib
import threading
from dataclasses import dataclass
from pathlib import Path

from chemsuite.qsar.model import Model, load

from . import Endpoint, by_id

DATA = Path(__file__).parent / "data"


@dataclass(frozen=True)
class Pin:
    """What a shipped archive must be: its size in bytes and its SHA-256."""

    size: int
    sha256: str


#: endpoint id -> the pin for its archive. Filled by ``train.py --ship``, which prints these
#: lines; a model refitted without updating them fails the test rather than shipping silently.
PINS: dict[str, Pin] = {
    "mdck_er": Pin(3714371, "57fc718bdcece446b2d5de16c109dec3740482c888984efa96cf96090aa77ff3"),
    "caco2": Pin(8750375, "726b733a92ec370c468958c4e70ab3f18f5380a2c9aa65049a914bc5b2091e7c"),
}

_lock = threading.Lock()
_models: dict[Path, Model] = {}


def path_for(endpoint: Endpoint | str) -> Path:
    """Where the archive for *endpoint* is shipped."""
    found = by_id(endpoint) if isinstance(endpoint, str) else endpoint
    return DATA / found.file


def digest(path: Path) -> str:
    """*path*'s SHA-256, read in blocks so a large file needs no more memory than a small one."""
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            sha.update(block)
    return sha.hexdigest()


def mismatch(endpoint: Endpoint | str) -> str:
    """Why the shipped archive is not the pinned one, or "" when it is.

    A test calls this for every endpoint. Loading does not, because the file ships inside the
    package: a digest per prediction run would buy nothing a release check does not.
    """
    found = by_id(endpoint) if isinstance(endpoint, str) else endpoint
    pin = PINS.get(found.id)
    path = path_for(found)
    if pin is None or not pin.sha256:
        return f"{found.id}: no pin recorded; rerun scripts/permeability/train.py --ship"
    if not path.is_file():
        return f"{found.id}: {path.name} is not shipped"
    size = path.stat().st_size
    if size != pin.size:
        return f"{found.id}: {path.name} is {size} bytes, pinned at {pin.size}"
    found_digest = digest(path)
    if found_digest != pin.sha256:
        return f"{found.id}: {path.name} hashes to {found_digest}, pinned at {pin.sha256}"
    return ""


def model(endpoint: Endpoint | str) -> Model:
    """The shipped model for *endpoint*, read once per process."""
    path = path_for(endpoint)
    with _lock:
        if path not in _models:
            if not path.is_file():
                raise FileNotFoundError(
                    f"the permeability model {path.name} is not installed at {path}"
                )
            _models[path] = load(path)
        return _models[path]


def installed(endpoint: Endpoint | str) -> bool:
    """Whether the archive for *endpoint* is there to be loaded."""
    return path_for(endpoint).is_file()
