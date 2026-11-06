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

"""The structure sketcher's chemistry: drawing model, undoable edits and rules.

Qt-free, like ``chemsuite.core``. The Qt canvas in ``chemsuite.app.sketcher`` only turns mouse
and key input into the edits defined here, so every rule is testable without a GUI. RDKit is
the chemistry oracle (valence, aromaticity, stereo); the rules for where things go on the page
follow Ketcher (Apache-2.0), cited per rule. See docs/SKETCHER_PLAN.md.
"""
