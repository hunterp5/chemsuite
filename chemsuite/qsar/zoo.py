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
"""The candidates a search tries, and the range of settings each one is tried over.

Every candidate is scikit-learn's, and every one of them can be written out as flat arrays
(:mod:`chemsuite.qsar.export`), which is the condition for being in the zoo at all: a model the
app cannot apply in numpy is a model the app would have to unpickle.

A candidate declares

- its **family**, which decides how it is exported and applied;
- whether it wants **scaled** columns (a penalized, kernel or distance model does; a tree does
  not, and scaling a count fingerprint only loses resolution in the archive);
- its **rank**, the order candidates are preferred in when a budget takes only some;
- its **space**: what the search samples from, and nothing fitted inside the model itself.

Settings that are not worth searching are fixed instead, with the reason: 500 trees because
more rarely helps and always costs, early stopping because the number of boosting rounds is
better found by the model than by the search.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

# -- what a space is -----------------------------------------------------------------------


@dataclass(frozen=True)
class Choice:
    """One of these values, each as likely."""

    values: tuple[Any, ...]

    def draw(self, rng: np.random.Generator) -> Any:
        return self.values[int(rng.integers(len(self.values)))]


@dataclass(frozen=True)
class Uniform:
    low: float
    high: float

    def draw(self, rng: np.random.Generator) -> float:
        return float(rng.uniform(self.low, self.high))


@dataclass(frozen=True)
class LogUniform:
    """Spread evenly over orders of magnitude, which is how a penalty or a rate behaves."""

    low: float
    high: float

    def draw(self, rng: np.random.Generator) -> float:
        return float(np.exp(rng.uniform(np.log(self.low), np.log(self.high))))


@dataclass(frozen=True)
class IntRange:
    low: int
    high: int  # inclusive

    def draw(self, rng: np.random.Generator) -> int:
        return int(rng.integers(self.low, self.high + 1))


Space = Choice | Uniform | LogUniform | IntRange


# -- a candidate ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Candidate:
    """One kind of model, with the range it is searched over."""

    key: str
    title: str
    task: str
    family: str
    rank: int
    scaled: bool
    space: dict[str, Space]
    build: Callable[[dict[str, Any], int], Any] = field(repr=False)
    #: narrows a drawn setting to what the data allows (components, neighbors)
    limit: Callable[[dict[str, Any], int, int], dict[str, Any]] | None = field(
        repr=False, default=None
    )
    #: the most classes this candidate can be written out for, or None for any number
    most_classes: int | None = None

    def draw(self, rng: np.random.Generator) -> dict[str, Any]:
        """One configuration, sampled from the space."""
        return {name: space.draw(rng) for name, space in self.space.items()}

    def settle(self, params: dict[str, Any], rows: int, features: int) -> dict[str, Any]:
        """*params* narrowed to what this much data allows."""
        return self.limit(params, rows, features) if self.limit else params

    def make(self, params: dict[str, Any], rows: int, features: int, seed: int = 0) -> Any:
        """A fresh, unfitted estimator for *params*."""
        return self.build(self.settle(params, rows, features), seed)


def _linear(params, seed):
    from sklearn.linear_model import Ridge

    return Ridge(alpha=params["alpha"], random_state=seed)


def _elasticnet(params, seed):
    from sklearn.linear_model import ElasticNet

    return ElasticNet(
        alpha=params["alpha"], l1_ratio=params["l1_ratio"], max_iter=5000, random_state=seed
    )


def _logistic(params, seed):
    from sklearn.linear_model import LogisticRegression

    return LogisticRegression(C=params["C"], max_iter=5000, random_state=seed)


def _pls(params, _seed):
    from sklearn.cross_decomposition import PLSRegression

    return PLSRegression(n_components=params["n_components"])


def _svr(params, _seed):
    from sklearn.svm import SVR

    return SVR(C=params["C"], gamma=params["gamma"], epsilon=params["epsilon"])


def _svc(params, _seed):
    from sklearn.svm import SVC

    # Left at its default of no Platt scaling, and not even named: scikit-learn deprecated
    # the parameter in 1.9. Conformal calibration wants a decision value, and probability=True
    # would fit the whole model five more times over to get one.
    return SVC(C=params["C"], gamma=params["gamma"])


def _forest(params, seed, *, classify: bool):
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

    kind = RandomForestClassifier if classify else RandomForestRegressor
    return kind(
        n_estimators=500,  # fixed: more trees rarely help and always cost
        max_features=params["max_features"],
        min_samples_leaf=params["min_samples_leaf"],
        max_depth=params["max_depth"],
        random_state=seed,
        n_jobs=1,  # the search runs the folds in parallel, not the trees
    )


def _boosting(params, seed, *, classify: bool):
    from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

    kind = HistGradientBoostingClassifier if classify else HistGradientBoostingRegressor
    return kind(
        max_iter=500,  # with early stopping, the model finds its own number of rounds
        early_stopping=True,
        n_iter_no_change=20,
        validation_fraction=0.15,
        learning_rate=params["learning_rate"],
        max_leaf_nodes=params["max_leaf_nodes"],
        min_samples_leaf=params["min_samples_leaf"],
        l2_regularization=params["l2_regularization"],
        random_state=seed,
    )


def _knn(params, _seed, *, classify: bool):
    from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor

    kind = KNeighborsClassifier if classify else KNeighborsRegressor
    return kind(
        n_neighbors=params["n_neighbors"], weights=params["weights"], metric="euclidean", n_jobs=1
    )


def _mlp(params, seed, *, classify: bool):
    from sklearn.neural_network import MLPClassifier, MLPRegressor

    kind = MLPClassifier if classify else MLPRegressor
    return kind(
        hidden_layer_sizes=tuple(params["hidden_layer_sizes"]),
        alpha=params["alpha"],
        learning_rate_init=params["learning_rate_init"],
        max_iter=1000,
        early_stopping=True,
        n_iter_no_change=20,
        random_state=seed,
    )


def _limit_components(params: dict[str, Any], rows: int, features: int) -> dict[str, Any]:
    """PLS cannot take more components than the data has dimensions to give."""
    most = max(1, min(params["n_components"], features, rows - 1))
    return {**params, "n_components": most}


def _limit_neighbors(params: dict[str, Any], rows: int, _features: int) -> dict[str, Any]:
    """A neighborhood cannot be larger than the training fold."""
    return {**params, "n_neighbors": max(1, min(params["n_neighbors"], rows - 1))}


_FOREST_SPACE: dict[str, Space] = {
    "max_features": Choice((0.1, 0.3, "sqrt", 1.0)),
    "min_samples_leaf": Choice((1, 2, 5)),
    "max_depth": Choice((None, 10, 20)),
}
_BOOSTING_SPACE: dict[str, Space] = {
    "learning_rate": LogUniform(0.01, 0.3),
    "max_leaf_nodes": Choice((15, 31, 63)),
    "min_samples_leaf": Choice((5, 10, 20)),
    "l2_regularization": Uniform(0.0, 1.0),
}
_KNN_SPACE: dict[str, Space] = {
    "n_neighbors": IntRange(1, 15),
    "weights": Choice(("uniform", "distance")),
}
_MLP_SPACE: dict[str, Space] = {
    "hidden_layer_sizes": Choice(((64,), (128,), (256, 64))),
    "alpha": LogUniform(1e-5, 1e-1),
    "learning_rate_init": LogUniform(1e-4, 1e-2),
}

_ALL: tuple[Candidate, ...] = (
        Candidate(
            "hgb", "Gradient boosting", "regression", "boosting", 1, False, _BOOSTING_SPACE,
            lambda p, s: _boosting(p, s, classify=False),
        ),
        Candidate(
            "rf", "Random forest", "regression", "forest", 2, False, _FOREST_SPACE,
            lambda p, s: _forest(p, s, classify=False),
        ),
        Candidate(
            "svr", "Support vectors (RBF)", "regression", "svm", 3, True,
            {"C": LogUniform(1e-2, 1e3), "gamma": LogUniform(1e-4, 1e0),
             "epsilon": Uniform(0.01, 0.5)},
            _svr,
        ),
        Candidate(
            "ridge", "Ridge", "regression", "linear", 4, True,
            {"alpha": LogUniform(1e-3, 1e3)}, _linear,
        ),
        Candidate(
            "pls", "Partial least squares", "regression", "pls", 5, True,
            {"n_components": IntRange(1, 20)}, _pls, _limit_components,
        ),
        Candidate(
            "elasticnet", "Elastic net", "regression", "linear", 6, True,
            {"alpha": LogUniform(1e-4, 1e1), "l1_ratio": Uniform(0.05, 0.95)}, _elasticnet,
        ),
        Candidate(
            "knn", "Nearest neighbors", "regression", "knn", 7, True, _KNN_SPACE,
            lambda p, s: _knn(p, s, classify=False), _limit_neighbors,
        ),
        Candidate(
            "mlp", "Neural network", "regression", "mlp", 8, True, _MLP_SPACE,
            lambda p, s: _mlp(p, s, classify=False),
        ),
        Candidate(
            "hgb", "Gradient boosting", "classification", "boosting", 1, False, _BOOSTING_SPACE,
            lambda p, s: _boosting(p, s, classify=True),
        ),
        Candidate(
            "rf", "Random forest", "classification", "forest", 2, False, _FOREST_SPACE,
            lambda p, s: _forest(p, s, classify=True),
        ),
        Candidate(
            "svc", "Support vectors (RBF)", "classification", "svm", 3, True,
            {"C": LogUniform(1e-2, 1e3), "gamma": LogUniform(1e-4, 1e0)}, _svc,
            # Above two classes scikit-learn's margin is a vote between pairs, which
            # chemsuite.qsar.export will not approximate, so the search never offers it.
            most_classes=2,
        ),
        Candidate(
            "logistic", "Logistic regression", "classification", "linear", 4, True,
            {"C": LogUniform(1e-3, 1e3)}, _logistic,
        ),
        Candidate(
            "knn", "Nearest neighbors", "classification", "knn", 5, True, _KNN_SPACE,
            lambda p, s: _knn(p, s, classify=True), _limit_neighbors,
        ),
        Candidate(
            "mlp", "Neural network", "classification", "mlp", 6, True, _MLP_SPACE,
            lambda p, s: _mlp(p, s, classify=True),
        ),
)  # fmt: skip

#: Keyed by task and key together, because "hgb" is a candidate of both tasks.
CANDIDATES: dict[tuple[str, str], Candidate] = {(c.task, c.key): c for c in _ALL}


def candidate(task: str, key: str) -> Candidate:
    """The candidate *key* for *task*."""
    try:
        return CANDIDATES[(task, key)]
    except KeyError:
        raise ValueError(f"no {task} candidate called {key!r}") from None


def for_task(
    task: str, *, most: int | None = None, without: Sequence[str] = (), classes: int = 0
) -> list[Candidate]:
    """The candidates a fit will try, in the order they are preferred.

    *most* keeps only that many — what a Quick budget does — and *without* leaves some out.
    *classes* drops any candidate that cannot be written out for that many classes, so a
    three-class target is never offered a model the archive would refuse at the end of a long
    search.
    """
    mine = [
        c
        for c in _ALL
        if c.task == task
        and c.key not in without
        and (c.most_classes is None or classes <= c.most_classes)
    ]
    mine.sort(key=lambda c: c.rank)
    return mine[:most] if most else mine
