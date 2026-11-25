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

"""The chemsuite pKa weights the app uses, shipped with it in ``data/`` (:data:`PATH`).

It is fold 0 of the 20-epoch training (``e20_f0.npz``, docs/UNIPKA_PLAN.md U2), exactly as
evaluated: one fold and 3 conformers plus the 2D layout (:data:`N_3D`) gave RMSE 0.81 on the
external tests, five folds 0.80 at five times the cost. The weights are GPL-3.0-or-later, as
chemsuite, fine-tuned from Uni-Mol's MIT-licensed backbone (``data/LICENSE-Uni-Mol.txt`` and
``data/NOTICE-chemsuite-pka.txt``). The same file is GitHub release ``pka-weights-1``;
:data:`SHA256` pins it, and a test checks the shipped copy against it.
"""

from __future__ import annotations

import threading
from pathlib import Path

from . import Predictor

NAME = "chemsuite-pka-1.npz"
PATH = Path(__file__).parent / "data" / NAME
SHA256 = "2a50506a62f671e5a0122093a6fda71dc91cd1cdca69e9bbff8ad068514c8090"
SIZE = 87_943_627
N_3D = 3  # 3D conformers per microstate, with the 2D layout as a fourth

_lock = threading.Lock()
_predictors: dict[Path, Predictor] = {}


def predictor(path: Path = PATH) -> Predictor:
    """One :class:`Predictor` per weights file and process, so free energies are cached
    across molecules and runs. Its parameters (about 180 MB) are mapped from a cache file that
    every process shares."""
    with _lock:
        if path not in _predictors:
            _predictors[path] = Predictor(path, n_3d=N_3D)
        return _predictors[path]
