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
"""The Constellation Plot: shared cores of analog series, placed in chemical space.

:class:`CoreFinder` finds the cores and their map in the background (``chemsuite.sar
.constellation`` does the chemistry) and keeps each structure's cores, so a filter change only
regroups and maps them again. :class:`ConstellationRenderer` draws what it found: a circle per
core sized by its compounds and coloured by a summary of a number column, lines between cores
that share compounds, and a label per analog series.
"""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import Executor
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QObject, Signal
from PySide6.QtGui import QImage
from rdkit import Chem

from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import numbers
from chemsuite.plotting.geometry import box_mask, lasso_mask
from chemsuite.plotting.style import (
    DEFAULT_COLOR,
    GREY,
    marker_sizes,
    scale_bins,
    scale_colors,
)
from chemsuite.sar.constellation import (
    MAPS,
    SUMMARIES,
    Constellation,
    constellation,
    find_cores,
    place,
    summarise,
)

from .spec import Hover, PlotSpec

SCALE_STEPS = 24
HIGHLIGHT = (230, 70, 20)
EDGE = (150, 150, 150, 110)
MAX_LABELS = 400  # more cores than this: no series labels (they would bury the plot)
HIT_PADDING = 3.0  # pixels round a circle that still hit it
DIMMED = 0.35  # opacity of the cores while some are selected
WAITING = "Finding cores…"

#: a job runner: ``run(label, work(job), on_done(result))``
Runner = Callable[[str, Callable, Callable], object]


@dataclass
class Placed:
    """Cores found for some rows, and where they sit on the map."""

    oids: np.ndarray  # the rows; the constellation's members index into it
    con: Constellation
    coords: np.ndarray  # k × 2


class CoreFinder(QObject):
    """Finds the constellation of the rows a plot shows, in the background when it can run
    jobs (without a runner, at once)."""

    ready = Signal()

    def __init__(
        self,
        table: Table,
        run: Runner | None = None,
        workers: Callable[[int], Executor | None] | None = None,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self.table = table
        self._run = run
        self._workers = workers
        self._cores: dict[tuple[float, object], tuple[str, ...]] = {}  # (fraction, cell) → cores
        self._key: tuple | None = None
        self._result: Placed | None = None
        self.message = WAITING

    def result(self, spec: PlotSpec, oids: np.ndarray) -> Placed | None:
        """The constellation for *spec* over *oids*; None while it is being found (then
        :attr:`ready` says when) or when it cannot be, :attr:`message` saying why."""
        column = spec.structure
        table = self.table
        if not column or column not in table:
            self.message = "Choose a structure column."
            return None
        cells = table.column(column).worker_inputs(table.slots_of(oids))
        fraction = spec.core_percent / 100.0
        key = (fraction, spec.min_compounds, spec.layout, oids.tobytes(), hash(tuple(cells)))
        if key == self._key:
            return self._result
        self._key, self._result, self.message = key, None, WAITING
        cache = self._cores
        known = {i: cache[(fraction, c)] for i, c in enumerate(cells) if (fraction, c) in cache}
        missing = [i for i in range(len(cells)) if i not in known]
        executor = self._workers(len(missing)) if self._workers is not None else None
        oids = oids.copy()

        def work(job=None) -> tuple:
            cancel = job.cancel_token if job is not None else None
            found = find_cores(
                [cells[i] for i in missing],
                fraction=fraction,
                executor=executor,
                progress=job.report if job is not None else None,
                cancel=cancel,
            )
            learned = {(fraction, cells[i]): cores for i, cores in zip(missing, found, strict=True)}
            lists = {**known, **dict(zip(missing, found, strict=True))}
            con = constellation(
                [lists[i] for i in range(len(cells))],
                min_compounds=spec.min_compounds,
                cancel=cancel,
            )
            try:
                return key, learned, Placed(oids, con, place(con, spec.layout, cancel=cancel)), ""
            except ValueError as error:  # a map the cores cannot make
                return key, learned, None, str(error)

        if self._run is None:
            self._done(work())
        else:
            self._run(WAITING, work, self._done)
        return self._result

    def close(self) -> None:
        """The plot closed: a job still running finishes unheard."""
        self._key = None

    def _done(self, outcome: tuple) -> None:
        key, learned, placed, problem = outcome
        self._cores.update(learned)
        if key != self._key:
            return  # a newer request is under way
        self._result = placed
        self.message = f"⚠ {problem}" if problem else ""
        self.ready.emit()


@lru_cache(maxsize=256)
def core_image(smiles: str, width: int, height: int, ratio: float, dark: bool) -> QImage | None:
    """A drawing of a core for its hover card."""
    from chemsuite.core.depict import fitted_size, render_png

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    w, h = fitted_size(mol, width, height)
    image = QImage.fromData(
        render_png(mol, round(w * ratio), round(h * ratio), scale=ratio, dark=dark)
    )
    image.setDevicePixelRatio(ratio)
    return image


class ConstellationRenderer:
    """Draws a :class:`Placed` constellation (``data``, set by the panel before :meth:`build`)."""

    def __init__(self, plot: pg.PlotItem) -> None:
        self.plot = plot
        self.data: Placed | None = None
        self.waiting = WAITING
        self.spec = PlotSpec(kind="constellation")
        self.cx = np.empty(0)
        self.cy = np.empty(0)
        self.sizes = np.empty(0)
        self.values = np.empty(0)
        self._member_oids: list[np.ndarray] = []
        self._items: list = []
        self._points: pg.ScatterPlotItem | None = None
        self._overlay = pg.ScatterPlotItem(
            pen=pg.mkPen(HIGHLIGHT, width=2.0), brush=None, pxMode=True
        )
        self._overlay.setZValue(10)
        self._legend: pg.LegendItem | None = None

    @property
    def oids(self) -> np.ndarray:
        return self.data.oids if self.data is not None else np.empty(0, dtype=np.int64)

    def clear(self) -> None:
        for item in self._items:
            self.plot.removeItem(item)
        self._items = []
        self._points = None
        if self._overlay.scene() is not None:
            self.plot.removeItem(self._overlay)
        if self._legend is not None:
            self._legend.clear()
            if self._legend.scene() is not None:
                self._legend.scene().removeItem(self._legend)
            self.plot.legend = None
            self._legend = None
        self.cx = self.cy = self.sizes = self.values = np.empty(0)
        self._member_oids = []

    # -- drawing --------------------------------------------------------------------------------
    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        """Draw; returns a status line."""
        self.clear()
        self.spec = spec
        data = self.data
        layout = MAPS.get(spec.layout, "Map")
        self.plot.setLabel("bottom", f"{layout} 1")
        self.plot.setLabel("left", f"{layout} 2")
        if data is None:
            return self.waiting
        con = data.con
        if not len(con):
            return (
                f"No core is shared by {spec.min_compounds} or more of the "
                f"{con.compounds:,} compounds."
            )
        self.cx, self.cy = data.coords[:, 0], data.coords[:, 1]
        self._member_oids = [data.oids[m] for m in con.members]
        n = con.n.astype(np.float64)
        self.sizes = marker_sizes(n, *spec.point_sizes, *spec.size_range)
        colors = self._colors(table, spec, data)
        if len(con.edges):
            a, b = con.edges[:, 0], con.edges[:, 1]
            x = np.column_stack([self.cx[a], self.cx[b]]).ravel()
            y = np.column_stack([self.cy[a], self.cy[b]]).ravel()
            self._add(pg.PlotDataItem(x, y, connect="pairs", pen=pg.mkPen(EDGE, width=1)))
        order = np.argsort(-self.sizes, kind="stable")  # small circles drawn over large ones
        self._points = pg.ScatterPlotItem(
            x=self.cx[order],
            y=self.cy[order],
            size=self.sizes[order],
            brush=[pg.mkBrush(*colors[i], 210) for i in order],
            pen=pg.mkPen(60, 60, 60, 160, width=0.8),
            pxMode=True,
        )
        self._add(self._points)
        if spec.labels and len(con) <= MAX_LABELS:
            self._draw_labels(con)
        self.plot.addItem(self._overlay)
        shown = sum(con.series_size.values())
        return (
            f"{len(con):,} cores · {len(con.series_size):,} analog series · "
            f"{shown:,} of {con.compounds:,} compounds in a series"
        )

    def _add(self, item) -> None:
        self.plot.addItem(item)
        self._items.append(item)

    def _colors(self, table: Table, spec: PlotSpec, data: Placed) -> list[tuple[int, int, int]]:
        k = len(data.con)
        if not spec.color or spec.color not in table or table.kind(spec.color) != NUMBER:
            self.values = np.full(k, np.nan)
            return [DEFAULT_COLOR] * k
        self.values = summarise(data.con, numbers(table, spec.color, data.oids), spec.summary)
        bins, lo, hi = scale_bins(self.values, SCALE_STEPS, *spec.color_range)
        steps = scale_colors(SCALE_STEPS, spec.cmap, spec.reverse)
        self._legend = self.plot.addLegend(offset=(-10, 10))
        title = f"{SUMMARIES[spec.summary]} {spec.color}"
        for step in (SCALE_STEPS - 1, SCALE_STEPS // 2, 0):
            at = lo + (hi - lo) * (step + 0.5) / SCALE_STEPS
            swatch = pg.ScatterPlotItem(size=10, brush=pg.mkBrush(*steps[step]), pen=None)
            label = f"{title}: {at:.3g}" if step == SCALE_STEPS - 1 else f"{at:.3g}"
            self._legend.addItem(swatch, label)
        if np.isnan(self.values).any():
            self._legend.addItem(pg.ScatterPlotItem(size=10, brush=pg.mkBrush(*GREY)), "No values")
        return [steps[b] if b >= 0 else GREY for b in bins]

    def _draw_labels(self, con: Constellation) -> None:
        """Each analog series named once, beside its largest core."""
        named: set[int] = set()
        for i in np.argsort(-con.n, kind="stable"):
            series = int(con.series[i])
            if series in named:
                continue
            named.add(series)
            text = pg.TextItem(str(series), color=(70, 70, 70), anchor=(-0.15, 1.1))
            text.setPos(float(self.cx[i]), float(self.cy[i]))
            self._add(text)

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        """Ring every core that holds a selected row."""
        chosen = (
            np.fromiter(selected, dtype=np.int64) if isinstance(selected, frozenset) else selected
        )
        hit = np.array(
            [len(chosen) > 0 and bool(np.isin(m, chosen).any()) for m in self._member_oids],
            dtype=bool,
        )
        if not len(hit):
            self._overlay.setData(x=[], y=[])
            return
        self._overlay.setData(x=self.cx[hit], y=self.cy[hit], size=self.sizes[hit] + 5)
        if self._points is not None:
            self._points.setOpacity(DIMMED if hit.any() else 1.0)

    # -- hits -----------------------------------------------------------------------------------
    def _rows(self, cores: np.ndarray) -> np.ndarray:
        if not len(cores):
            return np.empty(0, dtype=np.int64)
        return np.unique(np.concatenate([self._member_oids[i] for i in cores]))

    def hit_box(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        return self._rows(np.flatnonzero(box_mask(self.cx, self.cy, x0, x1, y0, y1)))

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        return self._rows(np.flatnonzero(lasso_mask(self.cx, self.cy, polygon)))

    def _core_at(self, x: float, y: float, sx: float, sy: float) -> int | None:
        if not len(self.cx):
            return None
        distance = np.hypot((self.cx - x) / max(sx, 1e-12), (self.cy - y) / max(sy, 1e-12))
        inside = distance <= self.sizes / 2 + HIT_PADDING
        if not inside.any():
            return None
        candidates = np.flatnonzero(inside)
        return int(candidates[np.argmin(self.sizes[candidates])])  # the smaller one, on top

    def hit_point(self, x: float, y: float, sx: float, sy: float) -> np.ndarray:
        i = self._core_at(x, y, sx, sy)
        return self._rows(np.array([i])) if i is not None else np.empty(0, dtype=np.int64)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._core_at(x, y, sx, sy)
        if i is None or self.data is None:
            return None
        con, spec = self.data.con, self.spec
        n = int(con.n[i])
        lines = [f"Series {int(con.series[i])} · {n:,} compound{'s' if n != 1 else ''}"]
        if len(self.values) and np.isfinite(self.values[i]):
            lines.append(f"{SUMMARIES[spec.summary]} {spec.color}: {self.values[i]:.3g}")
        return Hover(None, "", tuple(lines), smiles=con.cores[i])
