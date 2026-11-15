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

"""TMAP and Chemical Space Network plots: compounds as points joined by edges.

:class:`GraphFinder` builds the tree or network in the background
(``chemsuite.analysis.graphs`` does the work); :class:`GraphRenderer` draws it. A point is a
row, so colouring, sizing, selecting, hovering and point cards work as on a scatter plot.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg

from chemsuite.analysis.fingerprints import DEFAULT_KIND, TooManyPairs, fingerprint_task
from chemsuite.analysis.graphs.csn import Network, mmp_edges, network, similarity_edges
from chemsuite.analysis.graphs.tmap import Tree, tmap
from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import numbers
from chemsuite.plotting.style import scale_bins, scale_colors

from .finder import Finder, Outcome
from .renderers import ScatterRenderer
from .spec import Hover, PlotSpec

TMAP = "tmap"
NETWORK = "network"
GRAPH_KINDS = (TMAP, NETWORK)
#: how a network's edges join compounds (:attr:`PlotSpec.edges`)
EDGES = {"similarity": "Similarity", "mmp": "Matched Pairs"}
#: how edges are coloured (:attr:`PlotSpec.edge_color`)
EDGE_COLORS = {
    "plain": "Plain",
    "similarity": "Similarity",
    "difference": "Activity Difference",
}
EDGE = (150, 150, 150, 150)
EDGE_STEPS = 8


@dataclass(frozen=True)
class Drawn:
    """A tree or network over some rows (its indices run over :attr:`oids`)."""

    oids: np.ndarray
    graph: Tree | Network


class GraphFinder(Finder):
    """Finds the TMAP or network of the rows a plot shows."""

    waiting = "Finding neighbours…"

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._cells: list = []

    def key(self, spec: PlotSpec, oids: np.ndarray) -> tuple | None:
        column = spec.structure
        table = self.table
        if not column or column not in table:
            self.message = "Choose a structure column."
            return None
        if len(oids) < 3:
            self.message = "⚠ Fewer than three structures."
            return None
        self._cells = table.column(column).worker_inputs(table.slots_of(oids))
        if spec.kind == TMAP:
            settings: tuple = (spec.neighbours,)
        else:
            settings = (spec.edges, round(spec.threshold, 3), spec.singletons)
        return (spec.kind, column, *settings, oids.tobytes(), hash(tuple(self._cells)))

    def work(self, spec: PlotSpec, oids: np.ndarray) -> Callable[[object], Outcome]:
        oids = oids.copy()
        n = len(oids)
        executor = self.executor(n)
        if spec.kind == NETWORK and spec.edges == "mmp":
            cells = self._cells
            valid = np.array([c is not None and c != "" for c in cells], dtype=bool)

            def pairs(handle=None) -> Outcome:
                a, b = mmp_edges(cells, executor=executor, **_job_args(handle))
                net = network(
                    n, a, b, valid=valid, singletons=spec.singletons, cancel=_cancel(handle)
                )
                return Drawn(oids, net), "", None

            return pairs
        compute, remember = fingerprint_task(
            self.table, spec.structure or "", DEFAULT_KIND, oids, lambda: executor
        )

        def job(handle=None) -> Outcome:
            fps = compute(handle)
            if spec.kind == TMAP:
                return Drawn(oids, tmap(fps, spec.neighbours, **_job_args(handle))), "", None
            try:
                a, b, s = similarity_edges(fps, spec.threshold, **_job_args(handle))
            except TooManyPairs as many:
                problem = f"Over {many.count:,} pairs at {spec.threshold:.2f}: raise the threshold."
                return None, problem, lambda: remember(fps)
            net = network(
                n, a, b, s, valid=fps.valid, singletons=spec.singletons, cancel=_cancel(handle)
            )
            return Drawn(oids, net), "", lambda: remember(fps)

        return job


def _cancel(handle):
    return handle.cancel_token if handle is not None else None


def _job_args(handle) -> dict:
    if handle is None:
        return {}
    return {"progress": handle.report, "cancel": handle.cancel_token}


class GraphRenderer(ScatterRenderer):
    """Draws a :class:`Drawn` tree or network (``data``, set by the panel before :meth:`build`)."""

    def __init__(self, plot: pg.PlotItem) -> None:
        super().__init__(plot)
        self.data: Drawn | None = None
        self.waiting = GraphFinder.waiting
        self.degree = np.empty(0, dtype=np.int64)
        self._edges: list[pg.PlotDataItem] = []

    def clear(self) -> None:
        super().clear()
        for item in self._edges:
            self.plot.removeItem(item)
        self._edges = []
        for side in ("bottom", "left"):
            self.plot.showAxis(side)

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        """Draw; returns a status line."""
        self.clear()
        self.spec = spec
        for side in ("bottom", "left"):
            self.plot.hideAxis(side)  # the positions mean nothing on their own
        data = self.data
        if data is None:
            self.oids = np.empty(0, dtype=np.int64)
            self.x = self.y = self.sizes = np.empty(0)
            return self.waiting
        graph = data.graph
        placed = graph.placed
        self.oids = data.oids[placed]
        self.x, self.y = graph.coords[placed, 0], graph.coords[placed, 1]
        self.degree = graph.degree()[placed]
        self._draw_edges(table, spec, data)
        self._draw_points(table, spec)
        return self._status(spec, data)

    def _draw_edges(self, table: Table, spec: PlotSpec, data: Drawn) -> None:
        graph = data.graph
        if not len(graph.a):
            return
        coords = graph.coords
        how = spec.edge_color if spec.kind == NETWORK else "plain"
        if how == "similarity" and np.isfinite(graph.similarity).any():
            values = graph.similarity
        elif how == "difference" and spec.color and table.kind(spec.color) == NUMBER:
            activity = numbers(table, spec.color, data.oids)
            values = np.abs(activity[graph.a] - activity[graph.b])
        else:
            self._add_edges(coords, graph.a, graph.b, EDGE)
            return
        bins, _lo, _hi = scale_bins(values, EDGE_STEPS)
        colors = scale_colors(EDGE_STEPS, spec.cmap, spec.reverse)
        self._add_edges(coords, graph.a[bins < 0], graph.b[bins < 0], EDGE)
        for step in range(EDGE_STEPS):
            hit = bins == step
            if hit.any():
                self._add_edges(coords, graph.a[hit], graph.b[hit], (*colors[step], 220))

    def _add_edges(self, coords: np.ndarray, a: np.ndarray, b: np.ndarray, rgba) -> None:
        if not len(a):
            return
        x = np.column_stack([coords[a, 0], coords[b, 0]]).ravel()
        y = np.column_stack([coords[a, 1], coords[b, 1]]).ravel()
        item = pg.PlotDataItem(
            x, y, connect="pairs", pen=pg.mkPen(rgba, width=1), skipFiniteCheck=True
        )
        item.setZValue(-1)
        self.plot.addItem(item)
        self._edges.append(item)

    def _status(self, spec: PlotSpec, data: Drawn) -> str:
        graph = data.graph
        parts = graph.component[graph.component >= 0]
        sizes = np.bincount(parts) if len(parts) else np.empty(0, dtype=np.int64)
        shown = int(graph.placed.sum())
        status = [f"{shown:,} compounds"]
        if spec.kind == TMAP:
            status.append(f"{len(sizes):,} tree{'s' if len(sizes) != 1 else ''}")
        else:
            status.append(f"{len(graph.a):,} edges")
            status.append(f"{int((sizes > 1).sum()):,} clusters")
            status.append(f"{int((sizes == 1).sum()):,} singletons")
        missing = len(data.oids) - shown
        if missing and (spec.kind == TMAP or spec.singletons):
            status.append(f"{missing:,} without a structure")
        return " · ".join(status)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._nearest(x, y, sx, sy)
        if i is None:
            return None
        oid = int(self.oids[i])
        spec = self.spec
        lines = [
            f"{extra}: {table.cell(extra, oid)}"
            for extra in (spec.color, spec.size)
            if extra and extra in table
        ]
        word = "Neighbours" if spec.kind == NETWORK else "Branches"
        lines.append(f"{word}: {int(self.degree[i])}")
        return Hover(oid, "", tuple(lines))
