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
"""Import the heavy libraries in the background soon after start-up.

Tools import what they need when first opened, on the GUI thread: the plot stack (pyqtgraph
and parts of SciPy) took 0.4 s on a warm disk and 3 s on a cold one, a frozen window the
first time a plot, cluster or statistics panel opened. Importing them on a thread leaves the
window free (Python switches threads every few milliseconds). Only modules that create no Qt
objects when imported belong here: an object made on this thread would belong to it.
"""

from __future__ import annotations

import importlib
import logging
import threading

MODULES = (
    "scipy.optimize",
    "scipy.sparse.csgraph",
    "scipy.stats",
    "pyqtgraph",
    "chemsuite.app.plots.panel",
    "chemsuite.app.plots.graph",
    "chemsuite.app.plots.constellation",
    "chemsuite.superpose.run",
)

_log = logging.getLogger(__name__)


def start() -> threading.Thread:
    """Import :data:`MODULES` on a daemon thread; returns it (tests wait for it)."""
    thread = threading.Thread(target=_import_all, name="warm-up imports", daemon=True)
    thread.start()
    return thread


def _import_all() -> None:
    for name in MODULES:
        try:
            importlib.import_module(name)
        except ImportError:  # an optional part missing: the tool reports it when opened
            _log.debug("warm-up: could not import %s", name, exc_info=True)
