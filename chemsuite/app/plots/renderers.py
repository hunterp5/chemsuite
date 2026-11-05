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
"""What each plot type draws, and which rows a click, box, lasso or hover hits.

A renderer draws into a pyqtgraph ``PlotItem`` from table data and answers hit tests in data
coordinates with row oids; the panel owns the controls and the shared selection.

Scatter points are drawn as a few *layers* (one per colour step or category) rather than one
brush per point: pyqtgraph then draws 200k points quickly, and categories get a legend free.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg

from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import BLANK_LABEL, OTHER_LABEL, categories, numbers, xy
from chemsuite.plotting.geometry import box_mask, jitter, lasso_mask, nearest
from chemsuite.plotting.stats import (
    bin_of,
    box_stats,
    fit_distribution,
    histogram_edges,
    linear_fit,
)
from chemsuite.plotting.style import (
    DEFAULT_COLOR,
    GREY,
    category_color,
    marker_sizes,
    scale_bins,
    scale_colors,
)

SCALE_STEPS = 24
POINT_SIZE = 7.0
HIGHLIGHT = (230, 70, 20)
DIMMED = 0.3  # opacity of unselected points while something is selected
HOVER_RADIUS = 10.0
SPECIAL = (OTHER_LABEL, BLANK_LABEL)  # categories drawn grey


@dataclass(frozen=True)
class PlotSpec:
    kind: str = "scatter"  # or "histogram"
    x: str | None = None
    y: str | None = None
    color: str | None = None
    size: str | None = None
    bins: int | None = None  # None: automatic
    fit: bool = False  # scatter: least-squares line(s)
    curve: str | None = None  # histogram: "normal" or "lognormal" density overlay


@dataclass(frozen=True)
class Hover:
    """What is under the mouse: one row (``oid``) or a group of rows (a histogram bar)."""

    oid: int | None
    title: str
    lines: tuple[str, ...]


def _brush(rgb: tuple[int, int, int], alpha: int = 200):
    return pg.mkBrush(*rgb, alpha)


class ScatterRenderer:
    def __init__(self, plot: pg.PlotItem) -> None:
        self.plot = plot
        self.oids = np.empty(0, dtype=np.int64)
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.sizes = np.empty(0)
        self._layers: list[pg.ScatterPlotItem] = []
        self._legend: pg.LegendItem | None = None
        self._overlay = pg.ScatterPlotItem(
            pen=pg.mkPen(HIGHLIGHT, width=1.6), brush=_brush(HIGHLIGHT, 160), pxMode=True
        )
        self._overlay.setZValue(10)
        self._fit_lines: list[pg.PlotDataItem] = []
        self._color_codes: np.ndarray | None = None
        self.spec = PlotSpec()
        self._color_values: np.ndarray | None = None
        self._color_labels: list[str] | None = None

    def clear(self) -> None:
        for layer in self._layers:
            self.plot.removeItem(layer)
        self._layers = []
        if self._overlay.scene() is not None:
            self.plot.removeItem(self._overlay)
        for line in self._fit_lines:
            self.plot.removeItem(line)
        self._fit_lines = []
        if self._legend is not None:
            self._legend.clear()
            if self._legend.scene() is not None:
                self._legend.scene().removeItem(self._legend)
            self.plot.legend = None
            self._legend = None

    # -- drawing --------------------------------------------------------------------------------
    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        """Draw; returns a status line."""
        self.clear()
        self.spec = spec
        if not spec.x or not spec.y:
            return "Choose number columns for X and Y."
        data = xy(table, spec.x, spec.y, oids)
        self.oids, self.x, self.y = data.oids, data.x, data.y
        n = len(self.oids)
        self.sizes = (
            marker_sizes(numbers(table, spec.size, self.oids), 4.0, 18.0)
            if spec.size
            else np.full(n, POINT_SIZE)
        )
        self._color_values = None
        self._color_labels = None
        self._color_codes = None
        groups = self._groups(table, spec)
        for mask, rgb, label in groups:
            if not mask.any():
                continue
            layer = pg.ScatterPlotItem(
                x=self.x[mask],
                y=self.y[mask],
                # one size for all points is much faster than an array of equal sizes
                size=self.sizes[mask] if spec.size else POINT_SIZE,
                pen=None,
                brush=_brush(rgb),
                pxMode=True,
                name=label,
            )
            self.plot.addItem(layer)  # a named layer joins the legend by itself
            self._layers.append(layer)
        self.plot.addItem(self._overlay)
        self.plot.setLabel("bottom", spec.x)
        self.plot.setLabel("left", spec.y)
        parts = [f"{n:,} points"]
        if data.missing:
            parts.append(f"{data.missing:,} rows without both values")
        if spec.fit:
            parts.extend(self._draw_fits())
        return " · ".join(parts)

    def _draw_fits(self) -> list[str]:
        """One line overall, or one per category when colouring by a text column."""
        if self._color_labels is None or self._color_codes is None:
            groups = [(np.ones(len(self.x), dtype=np.bool_), (120, 120, 120), None)]
        else:
            groups = [
                (self._color_codes == code, category_color(code, label, special=SPECIAL), label)
                for code, label in enumerate(self._color_labels)
            ]
        notes = []
        for mask, rgb, label in groups:
            fit = linear_fit(self.x[mask], self.y[mask])
            if fit is None:
                continue
            xs = np.array([self.x[mask].min(), self.x[mask].max()])
            line = self.plot.plot(xs, fit.slope * xs + fit.intercept, pen=pg.mkPen(rgb, width=2))
            line.setZValue(8)
            self._fit_lines.append(line)
            notes.append(
                fit.describe() if label is None else f"{label}: R² = {fit.r2:.2f} (n = {fit.n:,})"
            )
        return notes

    def _groups(self, table: Table, spec: PlotSpec):
        n = len(self.oids)
        if not spec.color or spec.color not in table:
            return [(np.ones(n, dtype=np.bool_), DEFAULT_COLOR, None)]
        self._legend = self.plot.addLegend(offset=(-10, 10))
        if table.kind(spec.color) == NUMBER:
            values = numbers(table, spec.color, self.oids)
            self._color_values = values
            bins, lo, hi = scale_bins(values, SCALE_STEPS)
            colors = scale_colors(SCALE_STEPS)
            groups = [(bins == -1, GREY, BLANK_LABEL)]
            for step in range(SCALE_STEPS):
                at = lo + (hi - lo) * (step + 0.5) / SCALE_STEPS
                show = step in (0, SCALE_STEPS // 2, SCALE_STEPS - 1)
                groups.append((bins == step, colors[step], f"{at:.3g}" if show else None))
            return groups  # the legend names only the low, middle and high steps
        cats = categories(table, spec.color, self.oids)
        self._color_labels = cats.labels
        self._color_codes = cats.codes
        return [
            (
                cats.codes == code,
                category_color(code, label, special=SPECIAL),
                f"{label} ({cats.counts[code]:,})",
            )
            for code, label in enumerate(cats.labels)
        ]

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        chosen = (
            np.fromiter(selected, dtype=np.int64) if isinstance(selected, frozenset) else selected
        )
        mask = (
            np.isin(self.oids, chosen)
            if len(chosen) and len(self.oids)
            else np.zeros(len(self.oids), bool)
        )
        size = self.sizes[mask] + 2 if self.spec.size else POINT_SIZE + 2
        self._overlay.setData(x=self.x[mask], y=self.y[mask], size=size)
        opacity = DIMMED if mask.any() else 1.0
        for layer in self._layers:
            layer.setOpacity(opacity)

    # -- hits -----------------------------------------------------------------------------------
    def hit_box(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        return self.oids[box_mask(self.x, self.y, x0, x1, y0, y1)]

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        return self.oids[lasso_mask(self.x, self.y, polygon)]

    def _nearest(self, x: float, y: float, sx: float, sy: float) -> int | None:
        return nearest(self.x, self.y, x, y, sx, sy, HOVER_RADIUS)

    def hit_point(self, x: float, y: float, sx: float, sy: float) -> np.ndarray:
        i = self._nearest(x, y, sx, sy)
        return self.oids[[i]] if i is not None else np.empty(0, dtype=np.int64)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._nearest(x, y, sx, sy)
        if i is None:
            return None
        oid = int(self.oids[i])
        spec = self.spec
        lines = [f"{spec.x}: {table.cell(spec.x, oid)}", f"{spec.y}: {table.cell(spec.y, oid)}"]
        for extra in (spec.color, spec.size):
            if extra and extra in table and extra not in (spec.x, spec.y):
                lines.append(f"{extra}: {table.cell(extra, oid)}")
        return Hover(oid, "", tuple(lines))


class HistogramRenderer:
    """Counts per bin for a number column, or per value for a text column (a bar chart)."""

    def __init__(self, plot: pg.PlotItem) -> None:
        self.plot = plot
        self.oids = np.empty(0, dtype=np.int64)
        self.bins = np.empty(0, dtype=np.int64)  # bin per oid (-1: not drawn)
        self.left = np.empty(0)
        self.right = np.empty(0)
        self.counts = np.empty(0)
        self.labels: list[str] | None = None
        self.spec = PlotSpec(kind="histogram")
        self._bars: pg.BarGraphItem | None = None
        self._overlay: pg.BarGraphItem | None = None
        self._curve: pg.PlotDataItem | None = None

    def clear(self) -> None:
        for item in (self._bars, self._overlay, self._curve):
            if item is not None:
                self.plot.removeItem(item)
        self._bars = self._overlay = self._curve = None
        self.plot.getAxis("bottom").setTicks(None)

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        self.clear()
        self.spec = spec
        if not spec.x or spec.x not in table:
            return "Choose a column."
        self.oids = oids
        if table.kind(spec.x) == NUMBER:
            values = numbers(table, spec.x, oids)
            edges = histogram_edges(values, spec.bins)
            self.bins = bin_of(values, edges)
            self.left, self.right = edges[:-1], edges[1:]
            self.labels = None
            missing = int((self.bins < 0).sum())
            fit = fit_distribution(values, spec.curve) if spec.curve else None
        else:
            cats = categories(table, spec.x, oids, limit=40)
            self.bins = cats.codes
            positions = np.arange(len(cats.labels), dtype=float)
            self.left, self.right = positions - 0.4, positions + 0.4
            self.labels = cats.labels
            self.plot.getAxis("bottom").setTicks(
                [[(i, label[:18]) for i, label in enumerate(cats.labels)]]
            )
            missing = 0
            fit = None
        nbins = len(self.left)
        self.counts = np.bincount(self.bins[self.bins >= 0], minlength=nbins).astype(float)
        self._bars = pg.BarGraphItem(
            x0=self.left,
            x1=self.right,
            height=self.counts,
            brush=_brush(DEFAULT_COLOR, 190),
            pen=pg.mkPen(255, 255, 255, 120),
        )
        self.plot.addItem(self._bars)
        self.plot.setLabel("bottom", spec.x)
        self.plot.setLabel("left", "Rows")
        unit = "bins" if self.labels is None else "groups"
        status = f"{int(self.counts.sum()):,} rows in {nbins} {unit}"
        if missing:
            status += f" · {missing:,} blank or not a number"
        if fit is not None:
            xs, ys = fit.curve(edges)
            self._curve = self.plot.plot(xs, ys, pen=pg.mkPen((40, 40, 40), width=2))
            self._curve.setZValue(12)
            status += " · " + fit.describe()
        elif spec.curve and self.labels is None:
            status += " · too few distinct values to fit"
        return status

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        chosen = (
            np.fromiter(selected, dtype=np.int64) if isinstance(selected, frozenset) else selected
        )
        if self._overlay is not None:
            self.plot.removeItem(self._overlay)
            self._overlay = None
        if not len(chosen) or not len(self.oids):
            if self._bars is not None:
                self._bars.setOpacity(1.0)
            return
        mask = np.isin(self.oids, chosen) & (self.bins >= 0)
        heights = np.bincount(self.bins[mask], minlength=len(self.left)).astype(float)
        self._overlay = pg.BarGraphItem(
            x0=self.left, x1=self.right, height=heights, brush=_brush(HIGHLIGHT, 220), pen=None
        )
        self._overlay.setZValue(10)
        self.plot.addItem(self._overlay)
        if self._bars is not None:
            self._bars.setOpacity(0.55)

    def _bins_between(self, x0: float, x1: float) -> np.ndarray:
        lo, hi = sorted((x0, x1))
        return np.flatnonzero((self.right >= lo) & (self.left <= hi))

    def _rows_in(self, bins: np.ndarray) -> np.ndarray:
        return self.oids[np.isin(self.bins, bins)] if len(bins) else np.empty(0, dtype=np.int64)

    def hit_box(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        return self._rows_in(self._bins_between(x0, x1))

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        return self._rows_in(self._bins_between(polygon[:, 0].min(), polygon[:, 0].max()))

    def _bar_at(self, x: float, y: float) -> int | None:
        hit = np.flatnonzero((self.left <= x) & (x <= self.right))
        if not len(hit):
            return None
        i = int(hit[0])
        return i if 0 <= y <= max(self.counts[i], 0.0) or self.counts[i] == 0 else None

    def hit_point(self, x: float, y: float, sx: float, sy: float) -> np.ndarray:
        i = self._bar_at(x, y)
        return self._rows_in(np.array([i])) if i is not None else np.empty(0, dtype=np.int64)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._bar_at(x, y)
        if i is None:
            return None
        if self.labels is not None:
            title = self.labels[i]
        else:
            title = f"{self.left[i]:.4g} to {self.right[i]:.4g}"
        return Hover(None, title, (f"{int(self.counts[i]):,} rows",))


class BoxRenderer(ScatterRenderer):
    """Box plots of a number column per group of a text column, over a jittered strip.

    The strip has one point per row so rows can be clicked, boxed, lassoed and hovered just
    like a scatter plot (it reuses the scatter hit tests and selection overlay).
    """

    BOX_HALF = 0.35
    JITTER = 0.5
    MAX_GROUPS = 30

    def __init__(self, plot: pg.PlotItem) -> None:
        super().__init__(plot)
        self.codes = np.empty(0, dtype=np.int64)
        self.labels: list[str] = []
        self.stats: list = []
        self._box_items: list = []

    def clear(self) -> None:
        super().clear()
        for item in self._box_items:
            self.plot.removeItem(item)
        self._box_items = []
        self.plot.getAxis("bottom").setTicks(None)

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        self.clear()
        self.spec = PlotSpec(kind="box", x=spec.x, y=spec.y)
        if not spec.y or spec.y not in table:
            return "Choose a number column."
        values = numbers(table, spec.y, oids)
        ok = np.isfinite(values)
        missing = int((~ok).sum())
        self.oids, self.y = oids[ok], values[ok]
        if spec.x and spec.x in table:
            cats = categories(table, spec.x, self.oids, limit=self.MAX_GROUPS)
            self.codes, self.labels = cats.codes, cats.labels
        else:
            self.codes, self.labels = np.zeros(len(self.oids), dtype=np.int64), [spec.y]
        self.x = self.codes + jitter(self.oids, self.JITTER)
        self.sizes = np.full(len(self.oids), POINT_SIZE - 2)
        strip = pg.ScatterPlotItem(
            x=self.x,
            y=self.y,
            size=POINT_SIZE - 2,
            pen=None,
            brush=_brush(DEFAULT_COLOR, 110),
            pxMode=True,
        )
        self.plot.addItem(strip)
        self._layers = [strip]
        self.plot.addItem(self._overlay)
        self._draw_boxes()
        self.plot.getAxis("bottom").setTicks(
            [[(i, label[:18]) for i, label in enumerate(self.labels)]]
        )
        self.plot.setLabel("bottom", spec.x or "")
        self.plot.setLabel("left", spec.y)
        groups = len(self.labels)
        status = f"{len(self.oids):,} values in {groups} group{'s' if groups != 1 else ''}"
        return status + (f" · {missing:,} blank or not a number" if missing else "")

    def _draw_boxes(self) -> None:
        half = self.BOX_HALF
        self.stats = [box_stats(self.y[self.codes == i]) for i in range(len(self.labels))]
        drawn = [(i, b) for i, b in enumerate(self.stats) if b is not None]
        if not drawn:
            return
        boxes = pg.BarGraphItem(
            x0=[i - half for i, _ in drawn],
            x1=[i + half for i, _ in drawn],
            y0=[b.q1 for _, b in drawn],
            height=[b.q3 - b.q1 for _, b in drawn],
            brush=_brush(DEFAULT_COLOR, 45),
            pen=pg.mkPen(DEFAULT_COLOR, width=1.5),
        )
        xs: list[float] = []
        ys: list[float] = []
        for i, b in drawn:  # median, two whiskers and their caps, as line segments
            for x0, y0, x1, y1 in (
                (i - half, b.median, i + half, b.median),
                (i, b.q3, i, b.high),
                (i, b.q1, i, b.low),
                (i - half / 2, b.high, i + half / 2, b.high),
                (i - half / 2, b.low, i + half / 2, b.low),
            ):
                xs += [x0, x1]
                ys += [y0, y1]
        lines = pg.PlotDataItem(xs, ys, connect="pairs", pen=pg.mkPen(DEFAULT_COLOR, width=2))
        for item in (boxes, lines):
            item.setZValue(5)
            self.plot.addItem(item)
            self._box_items.append(item)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._nearest(x, y, sx, sy)
        if i is not None:
            oid = int(self.oids[i])
            lines = [f"{self.spec.y}: {table.cell(self.spec.y, oid)}"]
            if self.spec.x:
                lines.append(f"{self.spec.x}: {self.labels[int(self.codes[i])]}")
            return Hover(oid, "", tuple(lines))
        group = int(round(x))
        if 0 <= group < len(self.stats) and abs(x - group) <= self.BOX_HALF:
            box = self.stats[group]
            if box is not None and box.low <= y <= box.high:
                return Hover(None, self.labels[group], (box.describe(),))
        return None


RENDERERS = {"scatter": ScatterRenderer, "histogram": HistogramRenderer, "box": BoxRenderer}
