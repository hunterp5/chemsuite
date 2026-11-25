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

"""Per-user folders that Qt-free code (and worker processes) can find without the app."""

from __future__ import annotations

import os
import sys
from pathlib import Path

#: overrides every cache folder (the app's downloads use it too)
CACHE_ENV = "CHEMSUITE_CACHE_DIR"


def derived_cache() -> Path:
    """Where files made from shipped data are kept (they can always be made again):
    ``CHEMSUITE_CACHE_DIR/derived``, else the platform's per-user cache folder."""
    override = os.environ.get(CACHE_ENV, "").strip()
    if override:
        return Path(override) / "derived"
    home = Path.home()
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or str(home / "AppData" / "Local")
        return Path(base) / "chemsuite" / "Cache" / "derived"
    if sys.platform == "darwin":
        return home / "Library" / "Caches" / "chemsuite" / "derived"
    base = os.environ.get("XDG_CACHE_HOME") or str(home / ".cache")
    return Path(base) / "chemsuite" / "derived"
