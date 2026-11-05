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
"""Turning the feature chooser's choice into :class:`Features` inside a job."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from chemsuite.analysis.embedding import Features
from chemsuite.app.tools import ToolContext
from chemsuite.tools.similarity.common import FeatureChooser, fingerprint_task


def features_task(
    ctx: ToolContext, chooser: FeatureChooser, oids: np.ndarray
) -> tuple[Callable[[object], Features], Callable[[], None]]:
    """``(build(job) -> Features, remember())``: *remember* caches fingerprints afterwards."""
    if chooser.fingerprints():
        compute, keep = fingerprint_task(
            ctx, chooser.structure.currentText(), chooser.kind.currentData(), oids
        )
        made = {}

        def build(job) -> Features:
            made["fps"] = compute(job)
            return Features.from_fingerprints(made["fps"])

        def remember() -> None:
            if "fps" in made:
                keep(made["fps"])

        return build, remember
    columns = chooser.number_values(oids)
    return (lambda _job: Features.from_columns(columns)), (lambda: None)
