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
from chemsuite.analysis.fingerprints import counts_task, parse, variant
from chemsuite.app.document.tools import ToolContext
from chemsuite.tools.similarity.common import FeatureChooser, fingerprint_task


def features_task(
    ctx: ToolContext, chooser: FeatureChooser, oids: np.ndarray
) -> tuple[Callable[[object], Features], Callable[[], None]]:
    """``(build(job) -> Features, remember())``: *remember* caches fingerprints afterwards.

    *build* takes the running job (its ``report`` and ``cancel_token``)."""
    key = chooser.fingerprint_key()
    use = chooser.uses_columns() and chooser.number_columns()
    columns = chooser.number_values(oids) if use else []
    scaling = chooser.column_scaling()
    weight = chooser.fingerprint_weight()

    def combine(structures, job) -> Features:
        if not columns:
            return structures
        return Features.combined(
            structures, columns, weight=weight, scaling=scaling, cancel=job.cancel_token
        )

    if key is None:
        return (lambda _job: Features.from_columns(columns, scaling)), (lambda: None)
    structure = chooser.structure.currentText()
    if chooser.counts():
        v = parse(key)
        counted = variant(v.base, v.bits, v.radius, counts=True)
        compute_counts = counts_task(ctx.table, structure, counted, oids, ctx.process_pool)

        def build_counts(job) -> Features:
            counts, valid = compute_counts(job.report, job.cancel_token)
            return combine(Features.from_counts(counts, valid), job)

        return build_counts, (lambda: None)
    compute, keep = fingerprint_task(ctx, structure, key, oids)
    made = {}

    def build(job) -> Features:
        made["fps"] = compute(job)
        return combine(Features.from_fingerprints(made["fps"]), job)

    def remember() -> None:
        if "fps" in made:
            keep(made["fps"])

    return build, remember
