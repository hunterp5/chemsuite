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
from PySide6.QtCore import QPointF
from PySide6.QtGui import QColor, QPolygonF
from PySide6.QtWidgets import QGraphicsPolygonItem

from chemsuite.analysis.medchem_space import GUIDES
from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import BLANK_LABEL, OTHER_LABEL, categories, numbers, xy
from chemsuite.plotting.geometry import box_mask, jitter, lasso_mask, nearest
from chemsuite.plotting.stats import (
    AGGREGATES,
    aggregate,
    bin_of,
    box_stats,
    density_grid,
    fit_distribution,
    histogram_edges,
    linear_fit,
    mean_by_x,
    unit_scale,
)
from chemsuite.plotting.style import (
    DEFAULT_COLOR,
    GREY,
    category_color,
    marker_sizes,
    scale_bins,
    scale_colors,
)

from .scatter_item import ScatterItem

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
    guide: str | None = None  # scatter: a key of GUIDES (BOILED-Egg, Golden Triangle, LipE)
    aggregate: str = "mean"  # bar: a key of AGGREGATES
    errors: bool = False  # bar: standard deviation whiskers (with "mean")
    columns: tuple[str, ...] = ()  # radar: the axes
    z: str | None = None  # 3D: the third axis (the landscape's height)
    drop: bool = False  # 3D scatter: lines from each point down to the floor


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
        self._overlay = ScatterItem(
            pen=pg.mkPen(HIGHLIGHT, width=1.6), brush=_brush(HIGHLIGHT, 160), pxMode=True
        )
        self._overlay.setZValue(10)
        self._fit_lines: list[pg.PlotDataItem] = []
        self._guide_items: list = []
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
        for item in self._guide_items:
            self.plot.removeItem(item)
        self._guide_items = []
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
        self._draw_guide(spec.guide)
        groups = self._groups(table, spec)
        for mask, rgb, label in groups:
            if not mask.any():
                continue
            layer = ScatterItem(
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

    def _draw_guide(self, key: str | None) -> None:
        """Regions and reference lines behind the points, in data coordinates."""
        guide = GUIDES.get(key or "")
        if guide is None:
            return
        for region in guide.regions:
            polygon = QPolygonF([QPointF(x, y) for x, y in region.polygon])
            item = QGraphicsPolygonItem(polygon)
            fill = QColor(*region.rgb, 70)
            item.setBrush(fill)
            pen = pg.mkPen(QColor(*region.rgb, 230), width=1.5)
            pen.setCosmetic(True)
            item.setPen(pen)
            item.setZValue(-10)
            item.setToolTip(region.label)
            self.plot.addItem(item)
            self._guide_items.append(item)
        for slope, intercept, label in guide.lines:
            x = np.array([-50.0, 50.0])
            line = pg.PlotDataItem(
                x,
                slope * x + intercept,
                pen=pg.mkPen((150, 150, 150), width=1, style=pg.QtCore.Qt.PenStyle.DashLine),
            )
            line.setZValue(-10)
            self.plot.addItem(line, ignoreBounds=True)
            self._guide_items.append(line)
            # Label each line where it leaves the top of the data, clear of the points.
            finite = self.y[np.isfinite(self.y)]
            top = float(finite.max()) if len(finite) else intercept
            text = pg.TextItem(label, color=(130, 130, 130), anchor=(1, 0))
            text.setPos((top - intercept) / slope, top)
            self.plot.addItem(text, ignoreBounds=True)
            self._guide_items.append(text)

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
        strip = ScatterItem(
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


def _chosen(selected: frozenset[int] | np.ndarray) -> np.ndarray:
    return np.fromiter(selected, dtype=np.int64) if isinstance(selected, frozenset) else selected


class BarRenderer(HistogramRenderer):
    """One bar per group of a column: the mean (median, sum…) of a number column, or the row
    count. Clicking, boxing or lassoing a bar selects its rows."""

    MAX_GROUPS = 40

    def __init__(self, plot: pg.PlotItem) -> None:
        super().__init__(plot)
        self.spec = PlotSpec(kind="bar")
        self.how = "count"
        self.spread = np.empty(0)
        self.members = np.empty(0, dtype=np.int64)
        self._errors: pg.ErrorBarItem | None = None

    def clear(self) -> None:
        super().clear()
        if self._errors is not None:
            self.plot.removeItem(self._errors)
            self._errors = None

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        self.clear()
        self.spec = spec
        if not spec.x or spec.x not in table:
            return "Choose a column to group by."
        self.how = spec.aggregate if spec.y and spec.y in table else "count"
        self.oids = oids
        cats = categories(table, spec.x, oids, limit=self.MAX_GROUPS)
        self.bins, self.labels = cats.codes, cats.labels
        positions = np.arange(len(cats.labels), dtype=float)
        self.left, self.right = positions - 0.4, positions + 0.4
        values = numbers(table, spec.y, oids) if self.how != "count" else None
        groups = aggregate(values, self.bins, len(cats.labels), self.how)
        self.counts = np.nan_to_num(groups.value)  # bar heights
        self.spread, self.members = groups.spread, groups.count
        self._bars = pg.BarGraphItem(
            x0=self.left,
            x1=self.right,
            height=self.counts,
            brush=_brush(DEFAULT_COLOR, 190),
            pen=pg.mkPen(255, 255, 255, 120),
        )
        self.plot.addItem(self._bars)
        if spec.errors and self.how == "mean":
            spread = np.nan_to_num(self.spread)
            self._errors = pg.ErrorBarItem(
                x=positions,
                y=self.counts,
                top=spread,
                bottom=spread,
                beam=0.25,
                pen=pg.mkPen((60, 60, 60), width=1.5),
            )
            self._errors.setZValue(12)
            self.plot.addItem(self._errors)
        self.plot.getAxis("bottom").setTicks(
            [[(i, label[:18]) for i, label in enumerate(cats.labels)]]
        )
        self.plot.setLabel("bottom", spec.x)
        what = AGGREGATES[self.how]
        self.plot.setLabel("left", "Rows" if self.how == "count" else f"{what} {spec.y}")
        status = f"{len(oids):,} rows in {len(cats.labels)} groups"
        if self.how != "count":
            blank = len(oids) - int(self.members.sum())
            if blank:
                status += f" · {blank:,} blank or not a number"
        return status

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        chosen = _chosen(selected)
        if self._overlay is not None:
            self.plot.removeItem(self._overlay)
            self._overlay = None
        if self._bars is not None:
            self._bars.setOpacity(1.0)
        if not len(chosen) or not len(self.oids):
            return
        groups = np.unique(self.bins[np.isin(self.oids, chosen) & (self.bins >= 0)])
        if not len(groups):
            return
        self._overlay = pg.BarGraphItem(
            x0=self.left[groups],
            x1=self.right[groups],
            height=self.counts[groups],
            brush=_brush(HIGHLIGHT, 200),
            pen=None,
        )
        self._overlay.setZValue(10)
        self.plot.addItem(self._overlay)
        self._bars.setOpacity(0.55)

    def _bar_at(self, x: float, y: float) -> int | None:
        hit = np.flatnonzero((self.left <= x) & (x <= self.right))
        if not len(hit):
            return None
        i = int(hit[0])
        height = float(self.counts[i])
        return i if min(0.0, height) <= y <= max(0.0, height) or height == 0 else None

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._bar_at(x, y)
        if i is None:
            return None
        rows = int((self.bins == i).sum())
        if self.how == "count":
            return Hover(None, self.labels[i], (f"{rows:,} rows",))
        value = f"{AGGREGATES[self.how]}: {self.counts[i]:.4g}"
        if self.how == "mean" and np.isfinite(self.spread[i]):
            value += f" ± {self.spread[i]:.3g} (SD)"
        return Hover(None, self.labels[i], (value, f"{int(self.members[i]):,} of {rows:,} rows"))


class LineRenderer(ScatterRenderer):
    """Points joined in X order: one line, or one per category when coloured by text. Rows
    sharing an X value (replicates, a time point) are joined through their mean."""

    def __init__(self, plot: pg.PlotItem) -> None:
        super().__init__(plot)
        self._lines: list[pg.PlotDataItem] = []

    def clear(self) -> None:
        super().clear()
        for line in self._lines:
            self.plot.removeItem(line)
        self._lines = []

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        status = super().build(table, spec, oids)
        if not len(self.oids):
            return status
        if self._color_codes is not None and self._color_labels is not None:
            groups = [
                (self._color_codes == code, category_color(code, label, special=SPECIAL))
                for code, label in enumerate(self._color_labels)
            ]
        else:
            groups = [(np.ones(len(self.oids), dtype=np.bool_), DEFAULT_COLOR)]
        for mask, rgb in groups:
            if mask.sum() < 2:
                continue
            xs, ys = mean_by_x(self.x[mask], self.y[mask])
            if len(xs) < 2:
                continue
            line = pg.PlotDataItem(xs, ys, pen=pg.mkPen(rgb, width=2.5))
            line.setZValue(8)  # over the points, as fitted lines are
            self.plot.addItem(line)
            self._lines.append(line)
        return status


class DensityRenderer:
    """Where the points are crowded: a 2D histogram of X and Y, shaded by count (log scale).
    The points are not drawn, but box, lasso and click select the rows beneath."""

    def __init__(self, plot: pg.PlotItem) -> None:
        self.plot = plot
        self.oids = np.empty(0, dtype=np.int64)
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.grid = None
        self.spec = PlotSpec(kind="density")
        self._image: pg.ImageItem | None = None
        self._overlay = ScatterItem(
            pen=pg.mkPen(HIGHLIGHT, width=1.2), brush=_brush(HIGHLIGHT, 160), pxMode=True
        )
        self._overlay.setZValue(10)
        self._bar: pg.ColorBarItem | None = None

    def clear(self) -> None:
        if self._image is not None:
            self.plot.removeItem(self._image)
            self._image = None
        if self._overlay.scene() is not None:
            self.plot.removeItem(self._overlay)

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        self.clear()
        self.spec = spec
        if not spec.x or not spec.y:
            return "Choose number columns for X and Y."
        data = xy(table, spec.x, spec.y, oids)
        self.oids, self.x, self.y = data.oids, data.x, data.y
        self.grid = density_grid(self.x, self.y, spec.bins)
        self.plot.setLabel("bottom", spec.x)
        self.plot.setLabel("left", spec.y)
        if self.grid is None:
            return "No rows with both values."
        grid = self.grid
        shade = np.where(grid.counts > 0, np.log1p(grid.counts), np.nan)
        self._image = pg.ImageItem(shade)
        self._image.setColorMap(pg.colormap.get("viridis"))
        x0, x1 = grid.x_edges[0], grid.x_edges[-1]
        y0, y1 = grid.y_edges[0], grid.y_edges[-1]
        self._image.setRect(pg.QtCore.QRectF(x0, y0, x1 - x0, y1 - y0))
        self._image.setZValue(-5)
        self.plot.addItem(self._image)
        self.plot.addItem(self._overlay)
        nx, ny = grid.counts.shape
        parts = [f"{len(self.oids):,} points on a {nx} × {ny} grid"]
        parts.append(f"densest cell {int(grid.counts.max()):,}")
        if data.missing:
            parts.append(f"{data.missing:,} rows without both values")
        return " · ".join(parts)

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        chosen = _chosen(selected)
        mask = (
            np.isin(self.oids, chosen)
            if len(chosen) and len(self.oids)
            else np.zeros(len(self.oids), bool)
        )
        self._overlay.setData(x=self.x[mask], y=self.y[mask], size=POINT_SIZE)

    def hit_box(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        return self.oids[box_mask(self.x, self.y, x0, x1, y0, y1)]

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        return self.oids[lasso_mask(self.x, self.y, polygon)]

    def _cell_rows(self, x: float, y: float) -> tuple[tuple[int, int] | None, np.ndarray]:
        cell = self.grid.cell_at(x, y) if self.grid is not None else None
        if cell is None:
            return None, np.empty(0, dtype=np.int64)
        inside = (self.grid.cell[:, 0] == cell[0]) & (self.grid.cell[:, 1] == cell[1])
        return cell, self.oids[inside]

    def hit_point(self, x: float, y: float, sx: float, sy: float) -> np.ndarray:
        return self._cell_rows(x, y)[1]

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        cell, rows = self._cell_rows(x, y)
        if cell is None or not len(rows):
            return None
        g = self.grid
        i, j = cell
        return Hover(
            None,
            f"{len(rows):,} rows",
            (
                f"{self.spec.x}: {g.x_edges[i]:.4g} to {g.x_edges[i + 1]:.4g}",
                f"{self.spec.y}: {g.y_edges[j]:.4g} to {g.y_edges[j + 1]:.4g}",
            ),
        )


class RadarRenderer:
    """Several number columns as spokes, each scaled to its own range (centre: lowest value,
    rim: highest). Shows the mean of the plotted rows, and each selected row (up to
    ``MAX_ROWS``) as its own outline."""

    MAX_ROWS = 12
    RINGS = (0.25, 0.5, 0.75, 1.0)

    def __init__(self, plot: pg.PlotItem) -> None:
        self.plot = plot
        self.oids = np.empty(0, dtype=np.int64)
        self.spec = PlotSpec(kind="radar")
        self.columns: list[str] = []
        self.scaled = np.empty((0, 0))
        self.lows = np.empty(0)
        self.highs = np.empty(0)
        self._frame: list = []  # spokes, rings and labels
        self._shapes: list = []  # the mean and the selected rows
        self._vx = np.empty(0)  # outline corners for hover and hit tests
        self._vy = np.empty(0)
        self._voids = np.empty(0, dtype=np.int64)  # row per corner (-1: the mean)
        self._vaxis = np.empty(0, dtype=np.int64)

    def clear(self) -> None:
        for item in self._frame + self._shapes:
            self.plot.removeItem(item)
        self._frame, self._shapes = [], []
        self.plot.showAxis("left")
        self.plot.showAxis("bottom")
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.getViewBox().setAspectLocked(False)

    def _angles(self) -> np.ndarray:
        n = len(self.columns)
        return np.pi / 2 - 2 * np.pi * np.arange(n) / n  # first spoke straight up, clockwise

    def _corners(self, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r = np.nan_to_num(values, nan=0.0)
        angles = self._angles()
        return r * np.cos(angles), r * np.sin(angles)

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        self.clear()
        self.spec = spec
        self.columns = [c for c in spec.columns if c in table and table.kind(c) == NUMBER]
        if len(self.columns) < 3:
            return "Choose at least three number columns for the axes."
        self.oids = oids
        self.scaled, self.lows, self.highs = unit_scale(
            [numbers(table, c, oids) for c in self.columns]
        )
        self.plot.hideAxis("left")
        self.plot.hideAxis("bottom")
        self.plot.showGrid(x=False, y=False)
        self.plot.getViewBox().setAspectLocked(True)
        self._draw_frame()
        return f"{len(oids):,} rows on {len(self.columns)} axes"

    def _draw_frame(self) -> None:
        ink = (150, 150, 150)
        angles = self._angles()
        # invisible corners: fitting the view to the data keeps room for the axis labels
        room = ScatterItem(x=[-1.6, 1.6], y=[-1.25, 1.25], size=0, pen=None, brush=None)
        self.plot.addItem(room)
        self._frame.append(room)
        for ring in self.RINGS:
            xs = np.append(ring * np.cos(angles), ring * np.cos(angles[0]))
            ys = np.append(ring * np.sin(angles), ring * np.sin(angles[0]))
            item = pg.PlotDataItem(xs, ys, pen=pg.mkPen(*ink, 90, width=1))
            self.plot.addItem(item)
            self._frame.append(item)
        for name, angle in zip(self.columns, angles, strict=True):
            spoke = pg.PlotDataItem(
                [0, np.cos(angle)], [0, np.sin(angle)], pen=pg.mkPen(*ink, 140, width=1)
            )
            self.plot.addItem(spoke)
            cx, cy = np.cos(angle), np.sin(angle)
            anchor = (0.5 - 0.5 * cx, 0.5 + 0.5 * cy)
            label = pg.TextItem(name, color=(110, 110, 110), anchor=anchor)
            label.setPos(1.08 * cx, 1.08 * cy)
            # Text is sized in pixels: in a pane not yet laid out it would span a huge range,
            # and fitting the view to it drew the radar tiny. The "room" corners keep space.
            self.plot.addItem(label, ignoreBounds=True)
            self._frame += [spoke, label]

    def _outline(self, values: np.ndarray, rgb, alpha: int, width: float, oid: int) -> None:
        xs, ys = self._corners(values)
        item = pg.PlotDataItem(
            np.append(xs, xs[0]),
            np.append(ys, ys[0]),
            pen=pg.mkPen(rgb, width=width),
            fillLevel=None,
            symbol="o",
            symbolSize=5,
            symbolPen=None,
            symbolBrush=_brush(rgb, 230),
        )
        fill = QGraphicsPolygonItem(QPolygonF([QPointF(x, y) for x, y in zip(xs, ys, strict=True)]))
        fill.setBrush(QColor(*rgb, alpha))
        fill.setPen(pg.mkPen(None))
        for shape in (fill, item):
            self.plot.addItem(shape)
            self._shapes.append(shape)
        self._vx = np.append(self._vx, xs)
        self._vy = np.append(self._vy, ys)
        self._voids = np.append(self._voids, np.full(len(xs), oid))
        self._vaxis = np.append(self._vaxis, np.arange(len(xs)))

    def show_selection(self, selected: frozenset[int] | np.ndarray) -> None:
        for item in self._shapes:
            self.plot.removeItem(item)
        self._shapes = []
        self._vx, self._vy = np.empty(0), np.empty(0)
        self._voids = np.empty(0, dtype=np.int64)
        self._vaxis = np.empty(0, dtype=np.int64)
        if len(self.columns) < 3 or not len(self.oids):
            return
        with np.errstate(all="ignore"):
            mean = np.nanmean(self.scaled, axis=0)
        self._outline(mean, (120, 120, 120), 50, 2.0, -1)
        chosen = _chosen(selected)
        rows = np.flatnonzero(np.isin(self.oids, chosen))[: self.MAX_ROWS] if len(chosen) else []
        for n, row in enumerate(rows):
            rgb = category_color(n, str(n), special=())
            self._outline(self.scaled[row], rgb, 35, 1.8, int(self.oids[row]))

    def _nearest(self, x: float, y: float, sx: float, sy: float) -> int | None:
        return nearest(self._vx, self._vy, x, y, sx, sy, HOVER_RADIUS)

    def hit_box(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        inside = box_mask(self._vx, self._vy, x0, x1, y0, y1)
        return np.unique(self._voids[inside & (self._voids >= 0)])

    def hit_lasso(self, polygon: np.ndarray) -> np.ndarray:
        inside = lasso_mask(self._vx, self._vy, polygon)
        return np.unique(self._voids[inside & (self._voids >= 0)])

    def hit_point(self, x: float, y: float, sx: float, sy: float) -> np.ndarray:
        i = self._nearest(x, y, sx, sy)
        if i is None or self._voids[i] < 0:
            return np.empty(0, dtype=np.int64)
        return self._voids[[i]]

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._nearest(x, y, sx, sy)
        if i is None:
            return None
        oid, axis = int(self._voids[i]), int(self._vaxis[i])
        name = self.columns[axis]
        if oid < 0:
            with np.errstate(all="ignore"):
                mean = np.nanmean(numbers(table, name, self.oids))
            span = f"{self.lows[axis]:.4g} to {self.highs[axis]:.4g}"
            return Hover(None, "Mean of the plotted rows", (f"{name}: {mean:.4g} ({span})",))
        lines = tuple(f"{c}: {table.cell(c, oid)}" for c in self.columns)
        return Hover(oid, "", lines)


RENDERERS = {
    "scatter": ScatterRenderer,
    "line": LineRenderer,
    "histogram": HistogramRenderer,
    "bar": BarRenderer,
    "box": BoxRenderer,
    "density": DensityRenderer,
    "radar": RadarRenderer,
}
