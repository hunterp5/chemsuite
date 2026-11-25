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
"""Reading the table, and running a fit as a background job.

The split between this module and the engine is the whole point of the package layout: here
is where a ``Table``, a ``ToolContext`` and a ``Job`` exist, and :mod:`chemsuite.qsar` has
never heard of any of them. This module turns rows into arrays, hands them over, and gives the
answer back to the dialog.

A fit is whole-dataset work — half a search is not half a model — so it is **not** a
``keep_partial`` job: cancelling it throws the work away, as ``ToolContext.run`` describes for
an embedding or a clustering.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.shell.jobs import Job
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.qsar import features
from chemsuite.qsar.fit import Fitted, as_groups, cluster_groups, fit, scaffold_groups
from chemsuite.qsar.model import Model, Prediction
from chemsuite.qsar.recipe import Recipe
from chemsuite.qsar.refine import Merged, merge, refine

#: below this many rows the pool is built in this process, not farmed out
PARALLEL_FROM = 200

#: what the "cluster" split clusters on; the pool's own fingerprints are a different thing
#: (they follow the recipe, and may be absent from it altogether)
CLUSTER_FINGERPRINT = "morgan2"


@dataclass(frozen=True)
class Inputs:
    """What a fit needs from the table, read on the GUI thread before the job starts."""

    cells: list
    oids: np.ndarray
    columns: dict[str, np.ndarray]
    y: np.ndarray
    classes: np.ndarray | None
    order: np.ndarray | None

    def __len__(self) -> int:
        return len(self.oids)


def read_columns(
    ctx: ToolContext, wanted: Sequence[str], oids: np.ndarray
) -> dict[str, np.ndarray]:
    """The named columns as float arrays, with all-blank for one this document does not have.

    A model carries its recipe, and a recipe can name a column another document never had.
    The pool's own fill values then stand in, which is the same thing that happens to a blank
    cell — but whether that matters depends on the feature, so
    :func:`missing_features` is what decides whether to go ahead.
    """
    blank = np.full(len(np.asarray(oids)), np.nan, dtype=np.float64)
    table = ctx.table
    return {
        name: (numbers(table, name, oids) if name in table else blank.copy()) for name in wanted
    }


def missing_features(ctx: ToolContext, model: Model) -> list[str]:
    """Table columns the model actually relies on that this document does not have.

    Only the features the fit kept are counted: a recipe offers a wide pool, and a column the
    model looked at and discarded is no loss.
    """
    table = ctx.table
    return sorted(
        {name[2:] for name in model.kept_names() if name.startswith("c:") and name[2:] not in table}
    )


def read_inputs(ctx: ToolContext, recipe: Recipe, oids: np.ndarray) -> Inputs:
    """Everything the engine needs, as arrays. Reads the table, so it runs on the GUI thread."""
    table = ctx.table
    oids = np.asarray(oids, dtype=np.int64)
    cells = ctx.flat_structure_inputs(recipe.structures, oids)
    columns = read_columns(ctx, recipe.columns, oids)
    classes: np.ndarray | None = None
    if recipe.task == "classification":
        labels = ctx.texts(recipe.target, oids)
        classes = np.array(sorted({text for text in labels if text}), dtype=np.str_)
        index = {name: i for i, name in enumerate(classes.tolist())}
        y = np.array([index.get(text, -1) for text in labels], dtype=np.float64)
        y[y < 0] = np.nan
    else:
        y = numbers(table, recipe.target, oids)
    order = None
    if recipe.split == "column" and recipe.split_column in number_columns(table):
        order = numbers(table, recipe.split_column, oids)
    return Inputs(cells, oids, columns, y, classes, order)


def groups_for(recipe: Recipe, pool: features.Pool, cells: Sequence) -> np.ndarray | None:
    """Which rows must stay on the same side of a split, or None when any row may go anywhere.

    A scaffold or a cluster that straddles the split would let the model see a near-copy of
    the compound it is being judged on, which is how an honest-looking model turns out to
    predict nothing.
    """
    if recipe.split == "scaffold":
        return as_groups(scaffold_groups(features.mols_for(cells)))
    if recipe.split == "cluster":
        from chemsuite.analysis.clustering import butina
        from chemsuite.analysis.fingerprints import compute_fingerprints

        fps = compute_fingerprints(cells, CLUSTER_FINGERPRINT)
        return cluster_groups(butina(fps, cutoff=recipe.cluster_similarity).labels)
    return None


def _fitted(inputs: Inputs, recipe: Recipe, job: Job, executor) -> Fitted:
    """The work itself, on a worker thread: build the pool, then fit."""
    pool = features.build(
        inputs.cells,
        recipe,
        columns=inputs.columns,
        executor=executor,
        progress=job.report,
        cancel=job.cancel_token,
    )
    usable = pool.usable & np.isfinite(inputs.y)
    # A row with no measurement is not fitted on, so its target never matters; filling the
    # blanks keeps a NaN out of the integer class indices.
    y = np.nan_to_num(inputs.y, nan=0.0)
    keys = features.keys_for(inputs.cells)
    return fit(
        pool.x,
        pool.names,
        y,
        np.flatnonzero(usable),
        recipe,
        keys=keys,
        classes=inputs.classes,
        groups=groups_for(recipe, pool, inputs.cells),
        order=inputs.order,
        progress=job.report,
        cancel=job.cancel_token,
    )


def fit_in_job(
    ctx: ToolContext,
    recipe: Recipe,
    oids: np.ndarray,
    *,
    chosen: Sequence[str],
    on_done: Callable[[Fitted, np.ndarray], None],
) -> Job:
    """Start a fit. *on_done* is called on the GUI thread with the result and the rows."""
    inputs = read_inputs(ctx, recipe, oids)
    title = f"Fitting a model for {recipe.target}…"

    executor = ctx.workers_for(len(inputs), PARALLEL_FROM)

    def work(job: Job) -> Fitted:
        return _fitted(inputs, recipe, job, executor)

    def done(fitted: Fitted) -> None:
        on_done(fitted, inputs.oids)

    return ctx.run(title, work, done)


def predict_in_job(
    ctx: ToolContext,
    model: Model,
    oids: np.ndarray,
    *,
    on_done: Callable[[Prediction, np.ndarray, np.ndarray], None],
) -> Job:
    """Apply *model* to *oids*. *on_done* gets the prediction, the oids and the usable mask."""
    recipe = model.recipe
    structures = recipe.structures if recipe.structures in ctx.structure_columns() else ""
    if not structures:
        structures = ctx.default_structure_column() or ""
    oids = np.asarray(oids, dtype=np.int64)
    cells = ctx.flat_structure_inputs(structures, oids)
    columns = read_columns(ctx, recipe.columns, oids)

    executor = ctx.workers_for(len(oids), PARALLEL_FROM)

    def work(job: Job) -> tuple[Prediction, np.ndarray]:
        pool = features.build(
            cells,
            recipe,
            columns=columns,
            executor=executor,
            progress=job.report,
            cancel=job.cancel_token,
        )
        return model.predict(pool.x), pool.usable

    def done(result: tuple[Prediction, np.ndarray]) -> None:
        said, usable = result
        on_done(said, oids, usable)

    return ctx.run(f"Predicting {recipe.target}…", work, done)


def refine_in_job(
    ctx: ToolContext,
    part,
    asset: str,
    model: Model,
    oids: np.ndarray,
    *,
    budget: str,
    chosen: Sequence[str],
    replace_old: bool,
) -> Job:
    """Fit *model*'s recipe again on its own compounds plus *oids*, and keep the result."""
    from dataclasses import replace as _replace

    from .columns import write_refit

    recipe = _replace(model.recipe, budget=budget)
    inputs = read_inputs(ctx, recipe, oids)
    executor = ctx.workers_for(len(inputs), PARALLEL_FROM)

    def work(job: Job) -> tuple[Fitted, Merged]:
        pool = features.build(
            inputs.cells,
            recipe,
            columns=inputs.columns,
            executor=executor,
            progress=job.report,
            cancel=job.cancel_token,
        )
        keys = features.keys_for(inputs.cells)
        y = np.nan_to_num(inputs.y, nan=0.0)
        rows = np.flatnonzero(pool.usable & np.isfinite(inputs.y))
        merged = merge(model, keys, y, rows)
        fitted = refine(
            model,
            pool.x,
            pool.names,
            keys,
            y,
            rows,
            merged=merged,
            recipe=recipe,
            classes=inputs.classes,
            groups=groups_for(recipe, pool, inputs.cells),
            order=inputs.order,
            progress=job.report,
            cancel=job.cancel_token,
        )
        return fitted, merged

    def done(result: tuple[Fitted, Merged]) -> None:
        fitted, merged = result
        write_refit(
            ctx, part, asset, fitted, inputs.oids, list(chosen), merged, replace_old=replace_old
        )

    return ctx.run(f"Refining the model for {recipe.target}…", work, done)
