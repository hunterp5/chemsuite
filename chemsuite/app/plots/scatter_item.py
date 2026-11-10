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
"""A pyqtgraph scatter item that places points of one shared style without a per-point loop.

pyqtgraph looks up every point's (symbol, size, pen, brush) in Python before drawing, about
2 µs a point: a quarter of a second for 120k points, on every redraw of the data. When no
point has a style of its own (the plot's layers and selection overlay never do), all points
share one atlas entry, looked up once and copied to every point.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

#: the per-point style columns and the value meaning "use the item's own"
_DEFAULTS = (("symbol", None), ("size", -1), ("pen", None), ("brush", None))


class ScatterItem(pg.ScatterPlotItem):
    def updateSpots(self, dataSet=None):  # pyqtgraph names it
        data = self.data if dataSet is None else dataSet
        opts = self.opts
        if opts["pxMode"] and opts["useCache"] and not opts["hoverable"] and len(data):
            todo = data["sourceRect"]["w"] == 0
            if todo.any() and all(
                np.equal(data[name][todo], default).all() for name, default in _DEFAULTS
            ):
                style = (opts["symbol"], opts["size"], opts["pen"], opts["brush"])
                data["sourceRect"][todo] = self.fragmentAtlas[[style]][0]
                super().updateSpots(dataSet)  # nothing left to look up: the bookkeeping
                self.invalidate()
                return
        super().updateSpots(dataSet)
