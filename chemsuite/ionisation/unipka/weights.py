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

"""The chemsuite pKa weights the app uses: one archive, fetched once from chemsuite's GitHub
release and checked against its pinned SHA-256.

It is fold 0 of the 20-epoch training (``e20_f0.npz``, docs/UNIPKA_PLAN.md U2), exactly as
evaluated: one fold and 3 conformers plus the 2D layout (:data:`N_3D`) gave RMSE 0.81 on the
external tests, five folds 0.80 at five times the cost. The weights are GPL-3.0-or-later, as
chemsuite, fine-tuned from Uni-Mol's MIT-licensed backbone; the release carries both licenses
and a NOTICE. The download goes through
:mod:`chemsuite.core.download` (HTTPS only, audited, refused when IT policy turns the network
off; an interrupted download resumes).
"""

from __future__ import annotations

import hashlib
import threading
from collections.abc import Callable
from pathlib import Path

from chemsuite.core.download import DownloadError, Opener, download_file
from chemsuite.core.jobs import CancelToken

from . import Predictor

NAME = "chemsuite-pka-1.npz"
URL = f"https://github.com/hunterp5/chemsuite/releases/download/pka-weights-1/{NAME}"
SHA256 = "2a50506a62f671e5a0122093a6fda71dc91cd1cdca69e9bbff8ad068514c8090"
SIZE = 87_943_627
N_3D = 3  # 3D conformers per microstate, with the 2D layout as a fourth

_lock = threading.Lock()
_checked: dict[tuple[Path, float], bool] = {}
_predictors: dict[Path, Predictor] = {}


def _digest(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def cached(folder: Path) -> Path | None:
    """The weights in *folder*, when they are there and whole (checked once per process)."""
    path = folder / NAME
    if not path.is_file() or path.stat().st_size != SIZE:
        return None
    key = (path, path.stat().st_mtime)
    with _lock:
        if key not in _checked:
            _checked[key] = _digest(path) == SHA256
        return path if _checked[key] else None


def fetch(
    folder: Path,
    *,
    opener: Opener | None = None,
    cancel: CancelToken | None = None,
    report: Callable[[int], None] = lambda _bytes: None,
) -> Path:
    """The weights in *folder*, downloaded first when they are not there."""
    if (path := cached(folder)) is not None:
        return path
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / NAME
    download_file(URL, path, purpose="chemsuite pKa model", opener=opener, cancel=cancel,
                  report=report)  # fmt: skip
    if _digest(path) != SHA256:
        path.unlink()
        raise DownloadError(f"{NAME} did not match its checksum; it was removed.")
    with _lock:
        _checked[(path, path.stat().st_mtime)] = True
    return path


def predictor(path: Path) -> Predictor:
    """One :class:`Predictor` per weights file and process, so free energies are cached across
    ligands and runs."""
    with _lock:
        if path not in _predictors:
            _predictors[path] = Predictor(path, n_3d=N_3D)
        return _predictors[path]
