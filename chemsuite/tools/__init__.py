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
"""Tool packages, one folder per menu feature.

Tools are listed here explicitly rather than discovered at run time, so a frozen build
(PyInstaller) sees every import. Each module defines ``TOOL`` (see :mod:`chemsuite.app.tools`).
Keep these modules light: heavy imports belong inside the function that opens the tool.
"""

TOOL_MODULES: tuple[str, ...] = (
    "chemsuite.tools.sketcher.tool",
    "chemsuite.tools.descriptors.tool",
    "chemsuite.tools.standardize.tool",
    "chemsuite.tools.reactions.tool",
    "chemsuite.tools.pka.tool",
    "chemsuite.tools.plot.tool",
    "chemsuite.tools.conformers.tool",
    "chemsuite.tools.statistics.tool",
    "chemsuite.tools.mpo.tool",
    "chemsuite.tools.selectivity.tool",
    "chemsuite.tools.binning.tool",
    "chemsuite.tools.utilities.tool",
    "chemsuite.tools.similarity.tool",
    "chemsuite.tools.sar.tool",
    "chemsuite.tools.medchem.tool",
    "chemsuite.tools.space.tool",
)
