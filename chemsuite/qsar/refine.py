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
"""Refining a model: the same recipe, fitted again on the old rows plus the new ones.

A project's data arrives in batches, and a model that cannot take the next batch is a model
that is abandoned after a month. Refining is therefore a **refit, not an update**: the recipe
is the same, the rows are the union, and everything learned — fill values, which features to
keep, the hyperparameters, the calibration — is learned again from the larger set.

That is deliberate, and it is the honest option. A model nudged by new rows keeps a feature
selection and a scaling that were chosen on the old ones, and its reported score was measured
against folds the new rows were never in; after a few rounds, nobody can say what the number
means. A refit's score is as honest as the first fit's, because it is the same procedure.

What a refit must not do is quietly report the old score. :func:`refine` returns a model whose
scores are its own and whose lineage carries one line per round, so a panel can show that the
model has been round the loop and what it was worth each time.

**Which rows win.** Identity is the InChIKey. A new row whose key the model was trained on
replaces that training value — a measurement has been repeated or corrected — and a row with
a new key joins. A compound the model remembers that is no longer in the table cannot be
refitted, because its features are not in the pool this table built; it is left out, and the
lineage says how many. This is the one case where a refit sees fewer compounds than the model
did, so it is reported rather than hidden.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

import numpy as np

from chemsuite.core.jobs import CancelToken

from .fit import Fitted, fit
from .model import Model
from .recipe import Recipe
from .search import NoModel, Progress


@dataclass(frozen=True)
class Merged:
    """The training set a refit will use: the old rows, updated and extended.

    *rows* indexes the pool the caller built, with -1 for a row that is only in the old
    model's training set, whose values come from *keys* and *y* instead.
    """

    keys: np.ndarray
    y: np.ndarray
    rows: np.ndarray
    added: int
    changed: int
    carried: int

    def __len__(self) -> int:
        return len(self.y)

    @property
    def summary(self) -> str:
        parts = []
        if self.added:
            parts.append(f"{self.added:,} new")
        if self.changed:
            parts.append(f"{self.changed:,} remeasured")
        if self.carried:
            parts.append(f"{self.carried:,} carried over")
        return ", ".join(parts) or "no change"


def merge(model: Model, keys: Sequence[str], y: np.ndarray, rows: np.ndarray) -> Merged:
    """The union of the model's training rows and the table's, by InChIKey.

    *rows* are the usable rows of the pool the caller built, and *keys* and *y* cover the
    whole pool. A row with no key — a structure RDKit could not read — cannot be matched to
    anything and is treated as new.
    """
    y = np.asarray(y, dtype=np.float64)
    rows = np.asarray(rows, dtype=np.int64)
    old_keys = [str(k) for k in model.training.keys.tolist()]
    old_y = np.asarray(model.training.y, dtype=np.float64)

    merged_keys: list[str] = list(old_keys)
    merged_y: list[float] = [float(v) for v in old_y.tolist()]
    merged_rows: list[int] = [-1] * len(old_keys)
    where = {key: i for i, key in enumerate(old_keys) if key}

    added = changed = 0
    for row in rows.tolist():
        key = str(keys[row])
        at = where.get(key) if key else None
        if at is None:
            merged_keys.append(key)
            merged_y.append(float(y[row]))
            merged_rows.append(row)
            if key:
                where[key] = len(merged_keys) - 1
            added += 1
            continue
        if merged_rows[at] < 0:
            changed += int(not _same(merged_y[at], float(y[row])))
        merged_y[at] = float(y[row])
        merged_rows[at] = row
    carried = sum(1 for row in merged_rows if row < 0)
    return Merged(
        np.asarray(merged_keys, dtype=np.str_),
        np.asarray(merged_y, dtype=np.float64),
        np.asarray(merged_rows, dtype=np.int64),
        added,
        changed,
        carried,
    )


def _same(a: float, b: float) -> bool:
    if a != a and b != b:  # both blank
        return True
    return bool(np.isclose(a, b, rtol=1e-12, atol=0.0))


def refine(
    model: Model,
    pool: np.ndarray,
    names: Sequence[str],
    keys: Sequence[str],
    y: np.ndarray,
    rows: np.ndarray,
    *,
    merged: Merged | None = None,
    recipe: Recipe | None = None,
    classes: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    order: np.ndarray | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Fitted:
    """Fit *model*'s recipe again on its own rows plus the table's, and carry the lineage over.

    Rows the model remembers but the table no longer has cannot be refitted — their features
    are not in *pool* — so they are left out of the fit, and out of the refitted model's own
    training set, with the count recorded in the lineage. Everything the table has is used.
    """
    merged = merged or merge(model, keys, y, rows)
    usable = merged.rows >= 0
    if not usable.any():
        raise NoModel("none of the model's compounds are in this table, so there is nothing to fit")
    wanted = recipe or model.recipe
    if classes is None:
        classes = model.classes
    again = fit(
        pool,
        names,
        _scattered(len(pool), merged.rows[usable], merged.y[usable]),
        merged.rows[usable],
        wanted,
        keys=keys,
        classes=classes,
        groups=groups,
        order=order,
        progress=progress,
        cancel=cancel,
    )
    note = f"refit: {merged.summary}"
    if merged.carried:
        note += f" ({merged.carried:,} not in this table, left out)"
    # The refit's own "fitted" line is replaced by one "refit" line, so the history reads as
    # one entry per round rather than two.
    kept = replace(
        again.model,
        label=model.label,
        created=model.created or again.model.created,
        lineage=model.lineage,
    ).with_lineage(note, again.scores, int(usable.sum()))
    return Fitted(kept, again.outcome, again.rows, again.held, again.measured)


def _scattered(width: int, rows: np.ndarray, values: np.ndarray) -> np.ndarray:
    """*values* put back where their rows are in the pool, 0 elsewhere (never fitted on)."""
    out = np.zeros(width, dtype=np.float64)
    out[rows] = values
    return out
