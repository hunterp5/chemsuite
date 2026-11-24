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
"""The search: which model, with which settings, over which features — and what it scores.

Two loops, and the difference between them is the whole point.

**The inner loop chooses.** Over a training fold's own folds, it samples settings for every
candidate, throws out the worse half, widens the folds, and repeats until one pipeline is
left. Filling, feature selection and scaling are fitted on each inner fold's training rows, so
nothing a candidate is judged on has touched the rows it was fitted on.

**The outer loop judges.** The inner loop runs again from scratch on every outer fold's
training rows, and its winner predicts the rows held out. Those out-of-fold predictions are
where every reported number comes from, so the numbers measure *the automatic procedure* and
not one model that got lucky — the only thing worth reporting about a tool that chooses for
you.

The leaderboard is a softer, second thing: each candidate scored the same out-of-fold way. It
says what kind of model suits the data. Its maximum is optimistic, because that winner won a
comparison on those same folds, so the headline score stays the procedure's.

The model the app keeps is the inner loop run once over all the rows. Its score is never its
own fit; it is the outer loop's.

One thing is shared rather than repeated: **the basis** — the fill values and the selected
columns — depends on the training rows and the recipe, not on the candidate or its settings.
It is computed once per fold and reused by every configuration tried there. That is a hundred
times less work and changes no number, because it is still fitted on that fold's training rows
alone.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from . import metrics, zoo
from .prepare import Prepare
from .recipe import MIN_ROWS, Recipe
from .select import Selection, select
from .splits import Fold, make_folds

#: Called with (fits done, fits expected) as the search runs.
Progress = Callable[[int, int], None]


class NoModel(ValueError):
    """The data cannot carry a model: no feature with any information, or no fold that fits."""


@dataclass(frozen=True)
class Basis:
    """What a set of training rows says about the features, before any model is chosen."""

    prepare: Prepare
    selection: Selection

    @property
    def features(self) -> int:
        return len(self.selection)


@dataclass(frozen=True)
class Pipeline:
    """A fitted model, and everything needed to get the same numbers again."""

    task: str
    candidate: str
    params: dict[str, Any]
    prepare: Prepare
    selection: Selection
    estimator: Any = field(repr=False)
    classes: np.ndarray | None = None

    @property
    def family(self) -> str:
        return zoo.candidate(self.task, self.candidate).family

    @property
    def title(self) -> str:
        return zoo.candidate(self.task, self.candidate).title


@dataclass(frozen=True)
class Trial:
    """One candidate's settings, and the inner-fold score they earned."""

    candidate: str
    params: dict[str, Any]
    score: float
    folds: int  # the inner folds it survived to be scored on


@dataclass(frozen=True)
class Standing:
    """A candidate's place on the leaderboard.

    *scores* are out-of-fold, from the outer loop. *params* are the settings the final search
    over all the rows chose for this candidate — the ones a user would want to read — which
    are not necessarily the settings every outer fold used.
    """

    candidate: str
    title: str
    params: dict[str, Any]
    scores: dict[str, float]


@dataclass(frozen=True)
class Outcome:
    """A finished fit: the model to keep, its honest score, and the predictions behind it."""

    task: str
    model: Pipeline
    #: the procedure's out-of-fold metrics: the headline numbers
    scores: dict[str, float]
    leaderboard: tuple[Standing, ...]
    #: out-of-fold predictions, row for row — NaN, or -1 for a class, where a row was held out
    #: by no fold (only a column split leaves rows out)
    predicted: np.ndarray
    #: class score columns per row, for a classification
    class_scores: np.ndarray | None
    #: which outer fold held each row out (-1 for none)
    fold: np.ndarray
    classes: np.ndarray | None
    trials: tuple[Trial, ...]
    objective: str


# -- one fold's basis, and one pipeline on it ----------------------------------------------
def basis_for(x: np.ndarray, y: np.ndarray, rows: np.ndarray, recipe: Recipe) -> Basis:
    """Fill values and selected columns from *rows* alone."""
    prepare = Prepare.fit(x[rows])
    chosen = select(
        prepare.filled(x[rows]),
        y[rows],
        task=recipe.task,
        max_features=recipe.max_features,
        correlation=recipe.correlation,
        seed=recipe.seed,
    )
    if not len(chosen):
        raise NoModel("no feature carries any information on these rows")
    return Basis(prepare, chosen)


def fit_pipeline(
    x: np.ndarray,
    y: np.ndarray,
    rows: np.ndarray,
    candidate: zoo.Candidate,
    params: dict[str, Any],
    recipe: Recipe,
    basis: Basis,
    classes: np.ndarray | None = None,
) -> Pipeline:
    """Fit *candidate* on *rows*, over *basis*'s columns. Nothing outside *rows* is read."""
    prepare = basis.prepare.with_scaling(
        basis.prepare.filled(x[rows]), basis.selection.kept, candidate.scaled
    )
    design = prepare.design(x[rows], basis.selection.kept)
    settled = candidate.settle(params, len(rows), basis.features)
    estimator = candidate.make(settled, len(rows), basis.features, recipe.seed)
    estimator.fit(design, y[rows])
    return Pipeline(
        recipe.task, candidate.key, settled, prepare, basis.selection, estimator, classes
    )


def predict(model: Pipeline, x: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """The model's prediction for *rows*: a number, or a class index."""
    design = model.prepare.design(x[rows], model.selection.kept)
    out = np.asarray(model.estimator.predict(design))
    if model.task == "regression":
        return out.reshape(len(rows)).astype(np.float64)  # PLS predicts a column
    return out.reshape(len(rows)).astype(np.int64)


def class_scores(model: Pipeline, x: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """One score column per class: a probability, or a decision value for support vectors.

    Dispatched on the family rather than by asking the estimator what it can do: an estimator
    in the zoo that cannot score its classes is a mistake to fix here, not to work around.
    """
    design = model.prepare.design(x[rows], model.selection.kept)
    if model.family != "svm":
        return np.asarray(model.estimator.predict_proba(design), dtype=np.float64)
    decision = np.asarray(model.estimator.decision_function(design), dtype=np.float64)
    if decision.ndim == 1:  # two classes: one margin, which is the second class's score
        return np.column_stack([-decision, decision])
    return decision


def score_pipeline(
    model: Pipeline, x: np.ndarray, y: np.ndarray, rows: np.ndarray
) -> dict[str, float]:
    """Every metric for *model* on rows it was not fitted on."""
    predicted = predict(model, x, rows)
    scores = class_scores(model, x, rows) if model.task == "classification" else None
    return metrics.scores_for(model.task, y[rows], predicted, scores=scores, classes=model.classes)


# -- the inner loop -------------------------------------------------------------------------
def _objective_for(task: str, classes: np.ndarray | None) -> str:
    return metrics.objective_of(task, len(classes) if classes is not None else 2)


@contextmanager
def _quietly() -> Iterator[None]:
    """Convergence is a candidate's problem: a model that will not settle scores badly, and
    need not say so a thousand times over in the log."""
    from sklearn.exceptions import ConvergenceWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        yield


def halving_rounds(inner: int) -> list[int]:
    """How many inner folds each round scores on: 1, 2, 4… up to all of them."""
    rounds: list[int] = []
    folds = 1
    while folds < inner:
        rounds.append(folds)
        folds *= 2
    rounds.append(inner)
    return rounds


def _arena(
    candidates: Sequence[zoo.Candidate], configurations: int, rng: np.random.Generator
) -> list[tuple[zoo.Candidate, dict[str, Any]]]:
    """Every candidate gets the same number of configurations, so a cheap model is not crowded
    out before it has been tried. A small space is not sampled twice."""
    drawn: list[tuple[zoo.Candidate, dict[str, Any]]] = []
    for candidate in candidates:
        seen: list[dict[str, Any]] = []
        for _ in range(configurations):
            params = candidate.draw(rng)
            if params not in seen:
                seen.append(params)
                drawn.append((candidate, params))
    return drawn


def _ordered(alive: list[int], means: np.ndarray, objective: str) -> list[int]:
    """Best first, a score that could not be measured last."""
    higher = objective in metrics.HIGHER_IS_BETTER
    return sorted(
        alive,
        key=lambda i: (
            math.isnan(means[i]),
            -means[i] if higher else means[i],
            i,  # ties by draw order, so a seed decides everything
        ),
    )


@dataclass
class _Tally:
    """The running score of every configuration in the arena."""

    total: np.ndarray
    counts: np.ndarray

    def mean(self) -> np.ndarray:
        with np.errstate(invalid="ignore"):
            return np.where(self.counts > 0, self.total / np.maximum(self.counts, 1), np.nan)


def choose(
    x: np.ndarray,
    y: np.ndarray,
    rows: np.ndarray,
    recipe: Recipe,
    *,
    candidates: Sequence[zoo.Candidate],
    classes: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    order: np.ndarray | None = None,
    seed: int | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    spent: int = 0,
    expected: int = 0,
) -> tuple[list[Trial], int]:
    """The inner loop over *rows*: sample, halve, and return every candidate's best.

    The trials come back best first. The second value is how many fits have been spent in all,
    so an outer loop can keep one progress count across its folds.
    """
    objective = _objective_for(recipe.task, classes)
    budget = recipe.budget_spec
    seed = recipe.seed if seed is None else seed
    rng = np.random.default_rng(seed)
    arena = _arena(candidates, budget.configurations, rng)
    if not arena:
        raise NoModel("there is no candidate model to try")
    inner = make_folds(
        recipe.split, budget.inner, rows, y=y, groups=groups, order=order, seed=seed,
        stratify=recipe.task == "classification",
    )  # fmt: skip
    if not inner:
        raise NoModel("these rows cannot be divided into folds")
    bases = _bases_for(x, y, inner, recipe)

    tally = _Tally(np.zeros(len(arena)), np.zeros(len(arena), dtype=np.int64))
    alive = list(range(len(arena)))
    best: dict[str, Trial] = {}
    for used in halving_rounds(len(inner)):
        spent = _run_round(
            x, y, arena, alive, used, tally, bases, inner, recipe, classes, objective,
            progress=progress, cancel=cancel, spent=spent, expected=expected,
        )  # fmt: skip
        means = tally.mean()
        for i in alive:
            candidate, params = arena[i]
            trial = Trial(candidate.key, params, float(means[i]), int(tally.counts[i]))
            kept = best.get(candidate.key)
            if kept is None or metrics.is_better(objective, trial.score, kept.score):
                best[candidate.key] = trial
        if len(alive) <= 1:
            break
        ranked = _ordered(alive, means, objective)
        alive = ranked[: max(1, len(ranked) // 2)]
    trials = [best[key] for key in best]
    trials.sort(key=lambda t: (math.isnan(t.score), -t.score if objective in
                               metrics.HIGHER_IS_BETTER else t.score))  # fmt: skip
    if not trials or all(math.isnan(t.score) for t in trials):
        raise NoModel("no candidate could be fitted on these rows")
    return trials, spent


def _bases_for(x: np.ndarray, y: np.ndarray, folds: Sequence[Fold], recipe: Recipe) -> list[Basis]:
    """One basis per fold, computed once and shared by every configuration tried on it."""
    return [basis_for(x, y, fold.train, recipe) for fold in folds]


def _run_round(
    x: np.ndarray,
    y: np.ndarray,
    arena: Sequence[tuple[zoo.Candidate, dict[str, Any]]],
    alive: Sequence[int],
    used: int,
    tally: _Tally,
    bases: Sequence[Basis],
    folds: Sequence[Fold],
    recipe: Recipe,
    classes: np.ndarray | None,
    objective: str,
    *,
    progress: Progress | None,
    cancel: CancelToken | None,
    spent: int,
    expected: int,
) -> int:
    """Score every surviving configuration on the folds it has not seen yet."""
    with _quietly():
        for i in alive:
            candidate, params = arena[i]
            while tally.counts[i] < used:
                check(cancel)
                fold = int(tally.counts[i])
                value = _fold_score(
                    x, y, folds[fold], bases[fold], candidate, params, recipe, classes, objective
                )
                tally.total[i] += value if math.isfinite(value) else metrics.worst(objective)
                tally.counts[i] += 1
                spent += 1
                if progress is not None:
                    progress(spent, max(expected, spent))
    return spent


def _fold_score(
    x: np.ndarray,
    y: np.ndarray,
    fold: Fold,
    basis: Basis,
    candidate: zoo.Candidate,
    params: dict[str, Any],
    recipe: Recipe,
    classes: np.ndarray | None,
    objective: str,
) -> float:
    """One configuration's objective on one inner fold, or NaN where it cannot be fitted."""
    try:
        model = fit_pipeline(x, y, fold.train, candidate, params, recipe, basis, classes)
        return score_pipeline(model, x, y, fold.test)[objective]
    except (ValueError, ArithmeticError, MemoryError):
        return math.nan  # too few rows, a degenerate fold, or a space this model refuses


# -- the outer loop -------------------------------------------------------------------------
def expected_fits(recipe: Recipe, candidates: int) -> int:
    """Roughly how many estimator fits a run will cost, for the progress bar."""
    budget = recipe.budget_spec
    alive = max(candidates * budget.configurations, 1)
    scored = 0
    inner = 0
    for used in halving_rounds(budget.inner):
        inner += alive * (used - scored)
        scored = used
        alive = max(1, alive // 2)
    return (budget.outer + 1) * inner + budget.outer * candidates + 1


def run(
    x: np.ndarray,
    y: np.ndarray,
    rows: np.ndarray,
    recipe: Recipe,
    *,
    classes: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    order: np.ndarray | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Outcome:
    """Fit a model on *rows*, judged by the outer folds and kept from the final search.

    *x* is the whole feature pool, row for row with *y*; *rows* are the rows that have a value.
    For a classification, *y* holds class indices and *classes* their labels.
    """
    rows = np.sort(np.asarray(rows, dtype=np.int64))
    if len(rows) < MIN_ROWS:
        raise NoModel(
            f"{len(rows)} rows with a value: {MIN_ROWS} are the fewest a model can be judged on"
        )
    if recipe.task == "classification" and classes is None:
        raise ValueError("a classification needs its class labels")
    candidates = zoo.for_task(
        recipe.task,
        most=recipe.budget_spec.candidates,
        without=recipe.without,
        classes=0 if classes is None else len(classes),
    )
    if not candidates:
        raise NoModel("no candidate model can be fitted to this target")
    objective = _objective_for(recipe.task, classes)
    expected = expected_fits(recipe, len(candidates))

    outer = make_folds(
        recipe.split, recipe.budget_spec.outer, rows, y=y, groups=groups, order=order,
        seed=recipe.seed, stratify=recipe.task == "classification",
    )  # fmt: skip
    if not outer:
        raise NoModel("these rows cannot be divided into folds")

    held = _HeldOut(len(y), recipe.task, len(classes) if classes is not None else 0)
    spent = 0
    for number, fold in enumerate(outer):
        trials, spent = choose(
            x, y, fold.train, recipe, candidates=candidates, classes=classes, groups=groups,
            order=order, seed=recipe.seed + number, progress=progress, cancel=cancel,
            spent=spent, expected=expected,
        )  # fmt: skip
        basis = basis_for(x, y, fold.train, recipe)
        for rank, trial in enumerate(trials):
            check(cancel)
            candidate = zoo.candidate(recipe.task, trial.candidate)
            with _quietly():
                model = fit_pipeline(
                    x, y, fold.train, candidate, trial.params, recipe, basis, classes
                )
                held.take(model, x, fold, number, procedure=rank == 0)
            spent += 1
            if progress is not None:
                progress(spent, max(expected, spent))

    final, spent = choose(
        x, y, rows, recipe, candidates=candidates, classes=classes, groups=groups, order=order,
        seed=recipe.seed, progress=progress, cancel=cancel, spent=spent, expected=expected,
    )  # fmt: skip
    basis = basis_for(x, y, rows, recipe)
    winner = zoo.candidate(recipe.task, final[0].candidate)
    with _quietly():
        model = fit_pipeline(x, y, rows, winner, final[0].params, recipe, basis, classes)
    return Outcome(
        task=recipe.task,
        model=model,
        scores=held.scores(y, classes),
        leaderboard=held.leaderboard(y, classes, final, recipe.task),
        predicted=held.predicted,
        class_scores=held.class_scores,
        fold=held.fold,
        classes=classes,
        trials=tuple(final),
        objective=objective,
    )


class _HeldOut:
    """Out-of-fold predictions as the outer loop fills them in: the procedure's, and each
    candidate's, so the leaderboard is measured the same honest way as the headline."""

    def __init__(self, rows: int, task: str, classes: int) -> None:
        self.task = task
        self.classes = classes
        self.predicted = np.full(rows, -1 if task == "classification" else np.nan)
        self.class_scores = np.full((rows, classes), np.nan) if task == "classification" else None
        self.fold = np.full(rows, -1, dtype=np.int64)
        self.by_candidate: dict[str, np.ndarray] = {}
        self.scores_by_candidate: dict[str, np.ndarray] = {}

    def take(
        self, model: Pipeline, x: np.ndarray, fold: Fold, number: int, *, procedure: bool
    ) -> None:
        predicted = predict(model, x, fold.test)
        scores = class_scores(model, x, fold.test) if self.task == "classification" else None
        mine = self.by_candidate.setdefault(model.candidate, self._blank())
        mine[fold.test] = predicted
        if scores is not None:
            held = self.scores_by_candidate.setdefault(
                model.candidate, np.full((len(self.predicted), self.classes), np.nan)
            )
            held[fold.test] = scores
        if procedure:
            self.predicted[fold.test] = predicted
            self.fold[fold.test] = number
            if scores is not None and self.class_scores is not None:
                self.class_scores[fold.test] = scores

    def _blank(self) -> np.ndarray:
        """Room for one candidate's out-of-fold predictions, as "not predicted"."""
        if self.task == "classification":
            return np.full(len(self.predicted), -1, dtype=np.int64)
        return np.full(len(self.predicted), np.nan)

    def tested(self) -> np.ndarray:
        return np.flatnonzero(self.fold >= 0)

    def scores(self, y: np.ndarray, classes: np.ndarray | None) -> dict[str, float]:
        rows = self.tested()
        return metrics.scores_for(
            self.task,
            y[rows],
            self.predicted[rows],
            scores=None if self.class_scores is None else self.class_scores[rows],
            classes=classes,
        )

    def leaderboard(
        self,
        y: np.ndarray,
        classes: np.ndarray | None,
        final: Sequence[Trial],
        task: str,
    ) -> tuple[Standing, ...]:
        rows = self.tested()
        objective = _objective_for(task, classes)
        standings: list[Standing] = []
        for trial in final:
            predicted = self.by_candidate.get(trial.candidate)
            if predicted is None:
                continue
            held = self.scores_by_candidate.get(trial.candidate)
            scores = metrics.scores_for(
                task,
                y[rows],
                predicted[rows],
                scores=None if held is None else held[rows],
                classes=classes,
            )
            standings.append(
                Standing(
                    trial.candidate,
                    zoo.candidate(task, trial.candidate).title,
                    trial.params,
                    scores,
                )
            )
        standings.sort(
            key=lambda s: (
                math.isnan(s.scores.get(objective, math.nan)),
                -s.scores[objective] if objective in metrics.HIGHER_IS_BETTER
                else s.scores[objective],
            )
        )  # fmt: skip
        return tuple(standings)
