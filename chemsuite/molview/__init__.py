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
"""3D molecule scenes, without Qt: what the 3D Viewer draws and how it is looked at.

A scene is packed numpy arrays ready for the GPU: one row per atom sphere and one per bond
cylinder (:mod:`.geometry`). The camera (:mod:`.camera`), picking (:mod:`.pick`), measuring
(:mod:`.measure`) and conformer alignment (:mod:`.align`) work on plain coordinates, so all of
it is tested without a window. ``chemsuite.app.render3d`` draws a scene.
"""
