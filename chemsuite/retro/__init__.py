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

"""Retrosynthesis with a local AiZynthFinder install (docs/RETRO_PLAN.md).

AiZynthFinder needs RDKit 2023 and NumPy 1.26, so it is never imported here: it runs in its
own interpreter (``Settings ▸ External Tools``), driven by the worker in :mod:`.worker` and
kept warm between jobs by :class:`.runner.Worker`. The search is local — the models, templates
and stock come from a ``config.yml`` on this machine, and :func:`.config.config_problem`
refuses a config that would reach a model server or a remote stock database.

Nothing of AiZynthFinder's own code is vendored; the worker calls its public API.
"""

from .config import (
    SearchOptions,
    catalogue,
    config_problem,
    ready_problem,
    resolve_config,
    resolve_python,
)
from .routes import (
    OUTPUT_COLUMNS,
    compact_routes,
    leaf_stock_counts,
    route_tree,
    sort_routes,
    summary_columns,
)
from .runner import RetroError, SearchReply, Worker, search

__all__ = [
    "OUTPUT_COLUMNS",
    "RetroError",
    "SearchOptions",
    "SearchReply",
    "Worker",
    "catalogue",
    "compact_routes",
    "config_problem",
    "leaf_stock_counts",
    "ready_problem",
    "resolve_config",
    "resolve_python",
    "route_tree",
    "search",
    "sort_routes",
    "summary_columns",
]
