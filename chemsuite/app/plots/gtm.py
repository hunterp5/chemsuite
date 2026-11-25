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
"""The GTM Landscape plot: a trained Generative Topographic Map (``analysis.gtm``) as a
smooth image of a column over the map's nodes, with each compound a point at its mean
position (colored, sized, selected and hovered as on a scatter plot).

The map lives in the document as an asset (kind ``analysis.gtm.ASSET_KIND``), so the plot
needs only the table: the Generative Topographic Map tool trains it and opens this plot on it.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

from chemsuite.analysis.gtm import (
    ASSET_KIND,
    GtmError,
    GtmMap,
    class_landscape,
    landscape,
    unpack,
)
from chemsuite.core.assets import Asset
from chemsuite.core.table import Table
from chemsuite.core.table.columns import NUMBER
from chemsuite.plotting.data import categories, numbers
from chemsuite.plotting.style import category_color, scale_colors

from .renderers import SPECIAL, ScatterRenderer
from .spec import Hover, PlotSpec

GTM = "gtm"
#: :attr:`PlotSpec.gtm_value` for the landscape of where the compounds are
DENSITY = "Density"
UPSAMPLE = 8  # image pixels between neighboring nodes
LUT_STEPS = 256
#: a node holding this share of the median node's responsibility (or more) is drawn opaque;
#: emptier nodes fade out
FULL_AT = 0.5
MIN_ALPHA = 0.15
#: points with no Color column: dark, so they read over every color map
POINT_INK = (50, 50, 50)


def gtm_maps(table: Table) -> list[tuple[str, str]]:
    """``(asset id, name)`` of the document's maps, oldest first."""
    return [(a.id, a.name) for a in table.assets.values() if a.kind == ASSET_KIND]


def _upsample(grid: np.ndarray, factor: int) -> np.ndarray:
    """*grid* (side × side, or side × side × c) linearly interpolated to factor·(side−1)+1
    samples per side, through the node values."""
    side = grid.shape[0]
    if side < 2 or factor < 2:
        return grid
    fine = np.linspace(0.0, side - 1.0, factor * (side - 1) + 1)
    lower = np.minimum(np.floor(fine).astype(int), side - 2)
    frac = fine - lower
    shape = (-1,) + (1,) * (grid.ndim - 1)
    rows = grid[lower] * (1 - frac).reshape(shape) + grid[lower + 1] * frac.reshape(shape)
    shape = (1, -1) + (1,) * (grid.ndim - 2)
    return rows[:, lower] * (1 - frac).reshape(shape) + rows[:, lower + 1] * frac.reshape(shape)


def _alpha(mass: np.ndarray) -> np.ndarray:
    """Per node, how opaque: full from :data:`FULL_AT` of the median filled node, fading to
    transparent with no compounds."""
    filled = mass[mass > 0]
    full = FULL_AT * float(np.median(filled)) if len(filled) else 1.0
    alpha = np.clip(mass / max(full, 1e-12), 0.0, 1.0)
    return np.where(mass > 0, MIN_ALPHA + (1.0 - MIN_ALPHA) * alpha, 0.0)


class GtmRenderer(ScatterRenderer):
    """Draws the map named by :attr:`PlotSpec.gtm_map` (the newest map when None)."""

    def __init__(self, plot: pg.PlotItem) -> None:
        super().__init__(plot)
        self.map: GtmMap | None = None
        self.rows = np.empty(0, dtype=np.int64)  # each point's row in the map
        self._image: pg.ImageItem | None = None
        self._decoded: tuple[Asset, GtmMap] | None = None

    def clear(self) -> None:
        super().clear()
        if self._image is not None:
            self.plot.removeItem(self._image)
            self._image = None
        for side in ("bottom", "left"):
            self.plot.showAxis(side)
        self.plot.getViewBox().setAspectLocked(False)  # the next kind may not be square

    def _map(self, table: Table, wanted: str | None) -> GtmMap | None:
        maps = gtm_maps(table)
        if not maps:
            return None
        ids = [asset_id for asset_id, _name in maps]
        asset = table.assets[wanted if wanted in ids else ids[-1]]
        if self._decoded is not None and self._decoded[0] is asset:
            return self._decoded[1]
        payload = asset.payload
        try:
            gmap = payload if isinstance(payload, GtmMap) else unpack(payload)
        except GtmError:
            return None
        self._decoded = (asset, gmap)
        return gmap

    def build(self, table: Table, spec: PlotSpec, oids: np.ndarray) -> str:
        """Draw; returns a status line."""
        self.clear()
        self.spec = spec
        for side in ("bottom", "left"):
            self.plot.hideAxis(side)  # latent coordinates mean nothing on their own
        self.plot.getViewBox().setAspectLocked(True)  # the map is square
        gmap = self.map = self._map(table, spec.gtm_map)
        if gmap is None:
            self.oids = self.rows = np.empty(0, dtype=np.int64)
            self.x = self.y = self.sizes = np.empty(0)
            return (
                "No map: train one with Data ▸ Dimensionality Reduction ▸ "
                "Generative Topographic Map."
            )
        rows = gmap.lookup(oids)
        placed = rows >= 0
        self.oids, self.rows = np.asarray(oids)[placed], rows[placed]
        self.x, self.y = gmap.position[self.rows, 0], gmap.position[self.rows, 1]
        what = self._draw_landscape(table, spec, gmap)
        self._draw_points(table, spec)
        status = [f"{len(self.oids):,} compounds", f"{gmap.side} × {gmap.side} nodes"]
        if what:
            status.append(what)
        if len(self.rows):
            status.append(f"mean log-likelihood {float(gmap.loglik[self.rows].mean()):.3g}")
        missing = len(oids) - len(self.oids)
        if missing:
            status.append(f"{missing:,} not on the map")
        return " · ".join(status)

    def _draw_landscape(self, table: Table, spec: PlotSpec, gmap: GtmMap) -> str:
        side = gmap.side
        value = spec.gtm_value or DENSITY
        rgba = np.zeros((side * side, 4))
        note = ""
        if value != DENSITY and value in table and table.kind(value) == NUMBER:
            mean, mass = landscape(gmap, self.rows, numbers(table, value, self.oids))
            rgba[:, :3] = self._scale(mean, spec)
            rgba[:, 3] = np.where(np.isfinite(mean), _alpha(mass), 0.0)
            note = f"landscape of {value}"
        elif value != DENSITY and value in table:
            cats = categories(table, value, self.oids)
            best, purity, mass = class_landscape(gmap, self.rows, cats.codes, len(cats.labels))
            for code, label in enumerate(cats.labels):
                rgba[best == code, :3] = category_color(code, label, special=SPECIAL)
            rgba[:, 3] = np.where(best >= 0, _alpha(mass) * purity, 0.0)
            note = f"classes of {value}"
        else:
            mass = landscape(gmap, self.rows, np.ones(len(self.rows)))[1]
            shade = np.where(mass > 0, np.log1p(mass), np.nan)
            rgba[:, :3] = self._scale(shade, spec)
            rgba[:, 3] = _alpha(mass)
        # premultiply, so a faded node blends to the background, not to black
        grid = rgba.reshape(side, side, 4)
        grid[..., :3] *= grid[..., 3:4]
        fine = _upsample(grid, UPSAMPLE)
        alpha = fine[..., 3:4]
        rgb = np.where(alpha > 1e-6, fine[..., :3] / np.maximum(alpha, 1e-6), 0.0)
        image = np.concatenate([rgb, alpha * 255.0], axis=2).clip(0, 255).astype(np.uint8)
        self._image = pg.ImageItem(image)
        step = 2.0 / (side - 1) if side > 1 else 2.0
        pad = step / (2 * UPSAMPLE)  # half a pixel: pixel centers sit on the nodes' values
        self._image.setRect(pg.QtCore.QRectF(-1 - pad, -1 - pad, 2 + 2 * pad, 2 + 2 * pad))
        self._image.setZValue(-5)
        self.plot.addItem(self._image)
        return note

    def _groups(self, table: Table, spec: PlotSpec):
        if spec.color and spec.color in table:
            return super()._groups(table, spec)
        return [(np.ones(len(self.oids), dtype=np.bool_), POINT_INK, None)]  # over any landscape

    def _scale(self, values: np.ndarray, spec: PlotSpec) -> np.ndarray:
        """RGB per node along the plot's color map, over the landscape's own range (the
        Color range belongs to the points); NaN: black, drawn transparent."""
        lut = np.array(scale_colors(LUT_STEPS, spec.cmap, spec.reverse), dtype=np.float64)
        finite = values[np.isfinite(values)]
        lo = float(finite.min()) if len(finite) else 0.0
        hi = float(finite.max()) if len(finite) else 1.0
        span = hi - lo if hi > lo else 1.0
        with np.errstate(invalid="ignore"):
            index = np.clip((values - lo) / span * (LUT_STEPS - 1), 0, LUT_STEPS - 1)
        index = np.where(np.isfinite(index), index, 0).astype(int)
        return np.where(np.isfinite(values)[:, None], lut[index], 0.0)

    def hover(self, table: Table, x: float, y: float, sx: float, sy: float) -> Hover | None:
        i = self._nearest(x, y, sx, sy)
        if i is None or self.map is None:
            return None
        oid = int(self.oids[i])
        spec = self.spec
        lines = [
            f"{extra}: {table.cell(extra, oid)}"
            for extra in dict.fromkeys((spec.gtm_value, spec.color, spec.size))
            if extra and extra in table
        ]
        lines.append(f"Log-likelihood: {float(self.map.loglik[self.rows[i]]):.3g}")
        return Hover(oid, "", tuple(lines))
