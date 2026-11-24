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
"""The columns a fit or a prediction writes.

The names put the target first — "pIC50 Predicted", "pIC50 Low" — so they sort together
beside the column they came from in a wide table, and read plainly on a plot's axis.

What the numbers mean after a **fit** is the part worth being careful about: every one of them
comes from a fold that did not contain its own row. The prediction in the table is therefore
not the model's opinion of a compound it learned from — it is what the model would have said
about that compound had it never seen it. That is the only honest thing to show, and it is
also why these columns are not what you get from Apply Model… on the same rows.
"""

from __future__ import annotations

import numpy as np

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.table.columns import NumberColumn
from chemsuite.qsar.fit import Fitted
from chemsuite.qsar.metrics import TITLES
from chemsuite.qsar.model import Model, Prediction

#: how many decimals a written number keeps
DECIMALS = 3
#: the order the columns are written in
ORDER = ("prediction", "interval", "residual", "domain", "fold")


def names_for(target: str, chosen: list[str], task: str) -> list[str]:
    """The column names, target first, in the order they are written."""
    out: list[str] = []
    for key in ORDER:
        if key not in chosen:
            continue
        if key == "prediction":
            out.append(f"{target} Predicted")
            if task == "classification":
                out.append(f"{target} Confidence")
        elif key == "interval":
            out.extend(
                [f"{target} Low", f"{target} High"]
                if task == "regression"
                else [f"{target} Possible"]
            )
        elif key == "residual" and task == "regression":
            out.append(f"{target} Residual")
        elif key == "residual":
            out.append(f"{target} Correct")
        elif key == "domain":
            out.extend([f"{target} Distance", f"{target} In Domain"])
        elif key == "fold":
            out.append(f"{target} Fold")
    return out


def _numbers(values: np.ndarray) -> NumberColumn:
    return NumberColumn.from_floats(np.asarray(values, dtype=np.float64), DECIMALS)


def _yes_no(flags: np.ndarray) -> list[str]:
    return ["Yes" if flag else "No" for flag in np.asarray(flags, dtype=np.bool_).tolist()]


def values_for(
    model: Model,
    said: Prediction,
    chosen: list[str],
    *,
    measured: np.ndarray | None = None,
    fold: np.ndarray | None = None,
) -> list[object]:
    """The columns' contents, in the same order :func:`names_for` gives their names."""
    task = model.recipe.task
    classes = model.classes
    out: list[object] = []
    for key in ORDER:
        if key not in chosen:
            continue
        if key == "prediction":
            if task == "regression":
                out.append(_numbers(said.predicted))
            else:
                index = np.asarray(said.predicted, dtype=np.int64)
                labels = [] if classes is None else classes.tolist()
                out.append([labels[i] if 0 <= i < len(labels) else "" for i in index.tolist()])
                scores = said.scores
                best = (
                    np.full(len(index), np.nan)
                    if scores is None
                    else scores[np.arange(len(index)), np.clip(index, 0, scores.shape[1] - 1)]
                )
                out.append(_numbers(best))
        elif key == "interval":
            if task == "regression":
                out.append(_numbers(said.low))
                out.append(_numbers(said.high))
            else:
                out.append([] if classes is None else said.set_names(classes))
        elif key == "residual":
            if task == "regression":
                residual = (
                    np.full(len(said.predicted), np.nan)
                    if measured is None
                    else np.asarray(measured, dtype=np.float64) - said.predicted
                )
                out.append(_numbers(residual))
            else:
                index = np.asarray(said.predicted, dtype=np.int64)
                if measured is None:
                    out.append(["" for _ in index.tolist()])
                else:
                    truth = np.asarray(measured, dtype=np.float64)
                    right = np.isfinite(truth) & (truth.astype(np.int64) == index)
                    out.append(
                        [
                            "" if not np.isfinite(t) else ("Yes" if r else "No")
                            for t, r in zip(truth.tolist(), right.tolist(), strict=True)
                        ]
                    )
        elif key == "domain":
            out.append(_numbers(said.distance))
            out.append(_yes_no(said.inside))
        elif key == "fold":
            out.append(
                ["" for _ in range(len(said.predicted))]
                if fold is None
                else [str(f + 1) if f >= 0 else "" for f in np.asarray(fold).tolist()]
            )
    return out


def write_fit(ctx: ToolContext, fitted: Fitted, oids: np.ndarray, chosen: list[str]) -> None:
    """Write a fit's out-of-fold columns, keep the model with the document, and open the panel.

    The columns and the model are one undo step: undoing the fit takes the model away with its
    columns, rather than leaving a document holding a model whose results have gone.
    """
    from .export_python import remember
    from .panel import show_panel
    from .part import PART_KEY, QsarModels

    model = fitted.model
    rows = fitted.rows
    names = names_for(model.recipe.target, chosen, model.recipe.task)
    values = values_for(model, fitted.held, chosen, measured=fitted.measured, fold=fitted.fold)
    part = ctx.document_part(PART_KEY, lambda saved: QsarModels(ctx.table, saved))
    asset, op = part.put(model)
    # Kept for this session only, so "Export for Python…" has a scikit-learn object to pickle.
    remember(part, asset, fitted.outcome.model.estimator)
    fresh = ctx.table.unique_names(names)
    row_oids = np.asarray(oids, dtype=np.int64)[rows]
    written = {name: (row_oids, value) for name, value in zip(fresh, values, strict=True)}
    ctx.apply_results(written, label=f"Auto QSAR ({model.recipe.target})", extra_ops=[op])

    show_panel(ctx)
    ctx.show_message(_summary(model, len(rows)))


def _summary(model: Model, rows: int) -> str:
    """The status line: what was fitted, and how well it did on rows it never saw."""
    scores = model.scores
    key = model.objective
    shown = [k for k in ("r2", "rmse", "roc_auc", "mcc") if k in scores] or [key]
    parts = [f"{TITLES.get(k, k)} {scores[k]:.3g}" for k in shown[:2] if k in scores]
    return f"Fitted {model.candidate} on {rows:,} rows, {model.features} features" + (
        f"; {', '.join(parts)} out of fold" if parts else ""
    )


def write_refit(
    ctx: ToolContext,
    part,
    asset: str,
    fitted: Fitted,
    oids: np.ndarray,
    chosen: list[str],
    merged,
    *,
    replace_old: bool,
) -> None:
    """Write a refit's columns, store the refitted model, and say what changed.

    Replacing the old model and adding the new one is one undo step, so a refit the user
    dislikes comes off cleanly rather than leaving two models and two sets of columns.
    """
    from .export_python import remember
    from .panel import show_panel

    model = fitted.model
    rows = fitted.rows
    names = names_for(model.recipe.target, chosen, model.recipe.task)
    values = values_for(model, fitted.held, chosen, measured=fitted.measured, fold=fitted.fold)
    new_asset, put = part.put(model)
    remember(part, new_asset, fitted.outcome.model.estimator)
    ops = [put]
    if replace_old:
        drop = part.forget(asset)
        if drop is not None:
            ops.append(drop)
    fresh = ctx.table.unique_names(names)
    row_oids = np.asarray(oids, dtype=np.int64)[rows]
    written = {name: (row_oids, value) for name, value in zip(fresh, values, strict=True)}
    ctx.apply_results(written, label=f"Refine Model ({model.recipe.target})", extra_ops=ops)

    show_panel(ctx)
    ctx.show_message(f"Refitted on {len(rows):,} compounds ({merged.summary}); {_scores(model)}")


def _scores(model: Model) -> str:
    """The headline numbers, for a status line."""
    shown = [k for k in ("r2", "rmse", "roc_auc", "mcc") if k in model.scores]
    return (
        ", ".join(f"{TITLES.get(k, k)} {model.scores[k]:.3g}" for k in shown[:2]) or "fitted"
    ) + " out of fold"
