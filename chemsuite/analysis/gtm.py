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
"""Generative Topographic Mapping (Bishop, Svensén and Williams, *Neural Comput.* 1998, 10,
215–234), as used for "chemography" (Gaspar et al., *J. Chem. Inf. Model.* 2015, 55, 84).

A square grid of *k* × *k* nodes in a 2D latent space is bent into the data space by a smooth
map ``y = Φ(u) W`` (Gaussian radial basis functions plus a bias), and fitted by EM so the
data are a mixture of Gaussians around the bent nodes. Each compound then has a
*responsibility* over the nodes rather than one cell (as on a SOM):

- its mean latent position, for a scatter of compounds;
- a log-likelihood, how well the map explains it (an applicability-domain score);
- property landscapes: the responsibility-weighted mean of a column per node, or the
  dominant class of a category column.

The EM runs over the rows in chunks, so it never holds a rows × nodes matrix. The result keeps
only each row's largest responsibilities (sparse), which is all a landscape needs.
"""

from __future__ import annotations

import io
import json
import math
import zipfile
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from chemsuite.core.jobs import CancelToken, check

from .embedding import Features, reduce

Progress = Callable[[int, int], None]
#: what a map is stored as among a document's assets
ASSET_KIND = "gtm_map"
#: what a saved map's manifest says it is (raise it when the arrays change)
FORMAT = 1
MAX_SIDE = 40
CHUNK_ROWS = 8192
ITERATIONS = 100
REGULARIZATION = 0.1  # α, the precision of the Gaussian prior on W (ugtm's default)
RBF_WIDTH = 1.0  # times the spacing between basis-function centers (Netlab's convention)
TOLERANCE = 1e-4  # stop when the penalized log-likelihood per row and dimension gains less
KEEP_MASS = 0.999  # a row keeps its largest responsibilities up to this share...
KEEP_MAX = 64  # ...and no more nodes than this
_FORMAT_KEY = "format"
_META_KEY = "manifest"


class GtmError(ValueError):
    """A map that cannot be trained, or bytes that are not a map."""


def default_side(n: int) -> int:
    """Nodes per side for *n* training rows: 10 below 300, 16 below 3,000, else 25."""
    if n < 300:
        return 10
    return 16 if n < 3000 else 25


def default_basis(side: int) -> int:
    """Basis functions per side for *side* nodes per side: fewer, so the sheet stays smooth."""
    return max(3, round(side / 4))


def latent_grid(side: int) -> np.ndarray:
    """``side² × 2`` points on [−1, 1]², row-major (the first coordinate varies slowest)."""
    axis = np.linspace(-1.0, 1.0, side) if side > 1 else np.zeros(1)
    a, b = np.meshgrid(axis, axis, indexing="ij")
    return np.column_stack([a.ravel(), b.ravel()])


def basis(nodes: np.ndarray, per_side: int, width: float = RBF_WIDTH) -> np.ndarray:
    """Φ: each node's value on *per_side*² Gaussian basis functions, plus a bias column."""
    centers = latent_grid(per_side)
    spacing = 2.0 / (per_side - 1) if per_side > 1 else 2.0
    sigma = width * spacing
    d2 = ((nodes[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
    return np.column_stack([np.exp(-d2 / (2.0 * sigma * sigma)), np.ones(len(nodes))])


@dataclass(frozen=True, eq=False)
class GtmMap:
    """A trained map and where it put each mapped row.

    Row arrays run over :attr:`oids` (the rows that had features): mean latent position,
    the node with the largest responsibility, the log-likelihood per dimension, and the
    largest responsibilities in CSR form (``indptr``, ``node``, ``weight``; each row's
    weights sum to 1)."""

    side: int
    per_side: int
    width: float
    regularization: float
    weights: np.ndarray  # (per_side² + 1) × d
    beta: float
    oids: np.ndarray
    position: np.ndarray  # rows × 2
    mode: np.ndarray  # rows
    loglik: np.ndarray  # rows, per dimension
    indptr: np.ndarray  # rows + 1
    node: np.ndarray
    weight: np.ndarray
    iterations: int = 0
    history: tuple[float, ...] = ()  # penalized log-likelihood per row and dimension, per step
    recipe: dict[str, Any] = field(default_factory=dict)  # how the features were made (info)

    @property
    def nodes(self) -> np.ndarray:
        return latent_grid(self.side)

    @property
    def dims(self) -> int:
        return self.weights.shape[1]

    def __len__(self) -> int:
        return len(self.oids)

    def with_oids(self, oids: np.ndarray, recipe: dict[str, Any] | None = None) -> GtmMap:
        """The same map with its rows renamed (input positions → table oids)."""
        extra = {} if recipe is None else {"recipe": dict(recipe)}
        return replace(self, oids=np.asarray(oids, dtype=np.int64), **extra)

    def lookup(self, oids: np.ndarray) -> np.ndarray:
        """Each oid's row in the map, or −1 for oids it did not map."""
        oids = np.asarray(oids, dtype=np.int64)
        out = np.full(len(oids), -1, dtype=np.int64)
        if not len(self.oids):
            return out
        order = np.argsort(self.oids, kind="stable")
        known = self.oids[order]
        at = np.clip(np.searchsorted(known, oids), 0, len(known) - 1)
        hit = known[at] == oids
        out[hit] = order[at[hit]]
        return out

    # -- storing ---------------------------------------------------------------------------
    def to_bytes(self) -> bytes:
        """An uncompressed ``.npz``: a format number, a JSON manifest and numeric arrays,
        read back with ``allow_pickle=False``."""
        meta = {
            "side": self.side,
            "per_side": self.per_side,
            "width": self.width,
            "regularization": self.regularization,
            "beta": self.beta,
            "iterations": self.iterations,
            "history": list(self.history),
            "recipe": self.recipe,
        }
        arrays = {
            _FORMAT_KEY: np.array(FORMAT, dtype=np.int32),
            _META_KEY: np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
            "weights": self.weights,
            "oids": self.oids,
            "position": self.position,
            "mode": self.mode,
            "loglik": self.loglik,
            "indptr": self.indptr,
            "node": self.node,
            "weight": self.weight,
        }
        buffer = io.BytesIO()
        np.savez(buffer, **arrays)
        return buffer.getvalue()

    @property
    def nbytes(self) -> int:
        arrays = (self.weights, self.oids, self.position, self.mode, self.loglik)
        return sum(a.nbytes for a in (*arrays, self.indptr, self.node, self.weight))


def unpack(data: bytes) -> GtmMap:
    """A map from :meth:`GtmMap.to_bytes`. Never unpickles."""
    try:
        with np.load(io.BytesIO(data), allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files}
    except (ValueError, OSError, zipfile.BadZipFile) as error:
        raise GtmError(f"not a chemsuite map: {error}") from error
    version = int(np.ravel(arrays.pop(_FORMAT_KEY, np.array(0)))[0])
    if version < 1:
        raise GtmError("not a chemsuite map: it has no format number")
    if version > FORMAT:
        raise GtmError(
            f"saved by a newer chemsuite (map format {version}; this one reads {FORMAT})"
        )
    names = ("weights", "oids", "position", "mode", "loglik", "indptr", "node", "weight")
    missing = [name for name in (_META_KEY, *names) if name not in arrays]
    if missing:
        raise GtmError(f"the map is missing {', '.join(missing)}")
    meta = json.loads(arrays[_META_KEY].tobytes().decode())
    return GtmMap(
        side=int(meta["side"]),
        per_side=int(meta["per_side"]),
        width=float(meta["width"]),
        regularization=float(meta["regularization"]),
        weights=arrays["weights"],
        beta=float(meta["beta"]),
        oids=arrays["oids"].astype(np.int64),
        position=arrays["position"],
        mode=arrays["mode"].astype(np.int64),
        loglik=arrays["loglik"],
        indptr=arrays["indptr"].astype(np.int64),
        node=arrays["node"],
        weight=arrays["weight"],
        iterations=int(meta.get("iterations", 0)),
        history=tuple(float(v) for v in meta.get("history", ())),
        recipe=dict(meta.get("recipe", {})),
    )


# -- EM ------------------------------------------------------------------------------------
def _log_resp(x: np.ndarray, y: np.ndarray, beta: float) -> tuple[np.ndarray, np.ndarray]:
    """``(log responsibilities, log p(x))`` per row of the chunk *x* (rows × nodes, rows)."""
    d = x.shape[1]
    d2 = (x * x).sum(axis=1)[:, None] + (y * y).sum(axis=1)[None, :] - 2.0 * (x @ y.T)
    np.maximum(d2, 0.0, out=d2)
    log_p = -0.5 * beta * d2  # up to the constants below, the same for every node
    top = log_p.max(axis=1, keepdims=True)
    total = top[:, 0] + np.log(np.exp(log_p - top).sum(axis=1))
    log_r = log_p - total[:, None]
    const = 0.5 * d * math.log(beta / (2.0 * math.pi)) - math.log(y.shape[0])
    return log_r, total + const


def _initial(x: np.ndarray, nodes: np.ndarray, phi: np.ndarray, sample: int, seed: int):
    """W and β from the first principal components (Bishop et al. §2.3)."""
    mean = x.mean(axis=0)
    rng = np.random.default_rng(seed)
    part = x[rng.permutation(len(x))[:sample]] - mean
    _u, s, vt = np.linalg.svd(part, full_matrices=False)
    variance = s**2 / max(1, len(part) - 1)
    k = min(2, len(variance))
    target = np.tile(mean, (len(nodes), 1))
    for j in range(k):
        target += nodes[:, j : j + 1] * math.sqrt(variance[j]) * vt[j][None, :]
    weights = np.linalg.lstsq(phi, target, rcond=None)[0]
    y = phi @ weights
    side = int(round(math.sqrt(len(nodes))))
    if side > 1:  # half the distance between neighboring nodes, in data space
        grid = y.reshape(side, side, -1)
        step = np.concatenate(
            [np.linalg.norm(np.diff(grid, axis=0), axis=2).ravel(),
             np.linalg.norm(np.diff(grid, axis=1), axis=2).ravel()]
        )  # fmt: skip
        half = float(step.mean()) / 2.0
    else:
        half = 1.0
    third = float(variance[2]) if len(variance) > 2 else 0.0
    inverse = max(third, half * half, 1e-12)
    return weights, 1.0 / inverse


def train(
    x: np.ndarray,
    side: int,
    per_side: int,
    *,
    width: float = RBF_WIDTH,
    regularization: float = REGULARIZATION,
    iterations: int = ITERATIONS,
    tolerance: float = TOLERANCE,
    fit: np.ndarray | None = None,
    chunk: int = CHUNK_ROWS,
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, float, list[float]]:
    """EM on the rows of *x* (*fit*: those rows alone): ``(W, β, history)``, where *history*
    is the penalized log-likelihood per row and dimension before each M-step."""
    every = np.asarray(x, dtype=np.float64)
    x = every if fit is None else every[fit]
    n, d = x.shape
    if n < 2:
        raise GtmError("A map needs at least two compounds to train on.")
    nodes = latent_grid(side)
    phi = basis(nodes, per_side, width)
    weights, beta = _initial(x, nodes, phi, 20_000, seed)
    eye = np.eye(phi.shape[1])
    x2 = float((x * x).sum())
    history: list[float] = []
    for step in range(iterations):
        check(cancel)
        y = phi @ weights
        g = np.zeros(len(nodes))
        rx = np.zeros((len(nodes), d))
        loglik = 0.0
        for start in range(0, n, chunk):
            part = x[start : start + chunk]
            log_r, log_px = _log_resp(part, y, beta)
            r = np.exp(log_r)
            g += r.sum(axis=0)
            rx += r.T @ part
            loglik += float(log_px.sum())
        penalized = loglik - 0.5 * regularization * float((weights * weights).sum())
        history.append(penalized / (n * d))
        a = phi.T @ (g[:, None] * phi) + (regularization / beta) * eye
        weights = np.linalg.solve(a, phi.T @ rx)
        y = phi @ weights
        spread = x2 - 2.0 * float((y * rx).sum()) + float((g * (y * y).sum(axis=1)).sum())
        beta = (n * d) / max(spread, 1e-12)
        if progress is not None:
            progress(step + 1, iterations)
        if len(history) > 1 and abs(history[-1] - history[-2]) < tolerance:
            break
    return weights, beta, history


def place(
    x: np.ndarray,
    weights: np.ndarray,
    beta: float,
    side: int,
    per_side: int,
    width: float = RBF_WIDTH,
    *,
    chunk: int = CHUNK_ROWS,
    cancel: CancelToken | None = None,
) -> tuple[np.ndarray, ...]:
    """Every row of *x* on the map: ``(position, mode, loglik per dimension, indptr, node,
    weight)`` (see :class:`GtmMap`)."""
    x = np.asarray(x, dtype=np.float64)
    nodes = latent_grid(side)
    y = basis(nodes, per_side, width) @ weights
    n, d = x.shape
    keep = min(KEEP_MAX, len(nodes))
    position = np.empty((n, 2))
    mode = np.empty(n, dtype=np.int64)
    loglik = np.empty(n)
    counts = np.empty(n, dtype=np.int64)
    kept_nodes, kept_weights = [], []
    for start in range(0, n, chunk):
        check(cancel)
        part = x[start : start + chunk]
        stop = start + len(part)
        log_r, log_px = _log_resp(part, y, beta)
        r = np.exp(log_r)
        position[start:stop] = r @ nodes
        mode[start:stop] = r.argmax(axis=1)
        loglik[start:stop] = log_px / d
        top = np.argsort(-r, axis=1, kind="stable")[:, :keep]
        values = np.take_along_axis(r, top, axis=1)
        mass = np.cumsum(values, axis=1)
        count = np.minimum((mass < KEEP_MASS).sum(axis=1) + 1, keep)
        mask = np.arange(keep)[None, :] < count[:, None]
        values = values / np.where(mask, values, 0.0).sum(axis=1, keepdims=True)
        counts[start:stop] = count
        kept_nodes.append(top[mask])
        kept_weights.append(values[mask])
    indptr = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    node = np.concatenate(kept_nodes).astype(np.int16) if kept_nodes else np.empty(0, np.int16)
    weight = (
        np.concatenate(kept_weights).astype(np.float32) if kept_weights else np.empty(0, np.float32)
    )
    return position, mode, loglik, indptr, node, weight


def gtm(
    features: Features,
    side: int | None = None,
    per_side: int | None = None,
    *,
    width: float = RBF_WIDTH,
    regularization: float = REGULARIZATION,
    iterations: int = ITERATIONS,
    tolerance: float = TOLERANCE,
    fit: np.ndarray | None = None,
    seed: int = 0,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> GtmMap:
    """A map of *features* (reduced to at most 50 principal components, as for the SOM).
    *side* and *per_side* (None: :func:`default_side`, :func:`default_basis`) are nodes and
    basis functions per side; *fit* (input rows) trains on those rows alone, and every row is
    then placed. The map's oids are the input rows placed (see :meth:`GtmMap.with_oids`)."""
    if len(features) < 2:
        raise GtmError("A map needs at least two compounds with features.")
    fit = features.positions(fit)
    if fit is not None and len(fit) < 2:
        raise GtmError("A map needs at least two compounds to train on.")
    if side is None:
        side = default_side(len(features) if fit is None else len(fit))
    if not 2 <= side <= MAX_SIDE:
        raise GtmError(f"nodes per side are 2 to {MAX_SIDE}, not {side}")
    per_side = default_basis(side) if per_side is None else per_side
    if per_side < 2:
        raise GtmError("a map needs at least 2 basis functions per side")
    if width <= 0 or regularization < 0:
        raise GtmError("the basis width must be above 0 and the regularization not below 0")
    x = np.asarray(reduce(features, cancel=cancel), dtype=np.float64)
    weights, beta, history = train(
        x, side, per_side, width=width, regularization=regularization, iterations=iterations,
        tolerance=tolerance, fit=fit, seed=seed, progress=progress, cancel=cancel,
    )  # fmt: skip
    placed = place(x, weights, beta, side, per_side, width, cancel=cancel)
    return GtmMap(
        side, per_side, width, regularization, weights, beta,
        features.rows.astype(np.int64), *placed,
        iterations=len(history), history=tuple(history),
    )  # fmt: skip


# -- landscapes ----------------------------------------------------------------------------
def _expand(gmap: GtmMap, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(which of rows, node, weight)`` for every kept responsibility of *rows*."""
    starts = gmap.indptr[rows]
    lengths = gmap.indptr[rows + 1] - starts
    owner = np.repeat(np.arange(len(rows)), lengths)
    offset = np.arange(int(lengths.sum())) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    at = np.repeat(starts, lengths) + offset
    return owner, gmap.node[at].astype(np.int64), gmap.weight[at].astype(np.float64)


def density(gmap: GtmMap, rows: np.ndarray) -> np.ndarray:
    """Per node, the responsibility the map's *rows* (−1: skipped) put on it."""
    rows = np.asarray(rows, dtype=np.int64)
    rows = rows[rows >= 0]
    _owner, node, weight = _expand(gmap, rows)
    return np.bincount(node, weight, minlength=gmap.side**2)


def landscape(gmap: GtmMap, rows: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(value, density)`` per node: the responsibility-weighted mean of *values* (one per
    entry of *rows*, the map's rows; −1 or a NaN value is skipped), NaN where no row with a
    value reaches, and the responsibility those rows put on the node."""
    rows = np.asarray(rows, dtype=np.int64)
    values = np.asarray(values, dtype=np.float64)
    ok = (rows >= 0) & np.isfinite(values)
    rows, values = rows[ok], values[ok]
    owner, node, weight = _expand(gmap, rows)
    k = gmap.side**2
    mass = np.bincount(node, weight, minlength=k)
    total = np.bincount(node, weight * values[owner], minlength=k)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(mass > 0, total / np.where(mass > 0, mass, 1.0), np.nan)
    return mean, mass


def class_landscape(
    gmap: GtmMap, rows: np.ndarray, labels: np.ndarray, classes: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(dominant class, purity, density)`` per node for integer *labels* (0 to
    *classes* − 1, −1 for none): the class with the most responsibility (−1 where none
    reaches), its share of the node's responsibility, and that responsibility."""
    rows = np.asarray(rows, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int64)
    ok = (rows >= 0) & (labels >= 0)
    rows, labels = rows[ok], labels[ok]
    owner, node, weight = _expand(gmap, rows)
    k = gmap.side**2
    per = np.bincount(node * classes + labels[owner], weight, minlength=k * max(classes, 1))
    per = per.reshape(k, max(classes, 1))
    mass = per.sum(axis=1)
    best = np.where(mass > 0, per.argmax(axis=1), -1)
    with np.errstate(invalid="ignore", divide="ignore"):
        purity = np.where(mass > 0, per.max(axis=1) / np.where(mass > 0, mass, 1.0), 0.0)
    return best, purity, mass
