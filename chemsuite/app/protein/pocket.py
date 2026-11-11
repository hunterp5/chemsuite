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
"""The camera for a focused pocket: framed when focused, then a slab of depth around it and
fog through its back half as the camera turns, so the front of the pocket is sharp, its back
fades, and whatever lies in front of it or behind it is cut away."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import QObject

from chemsuite.app.molview.rhi_view import NO_SLAB

if TYPE_CHECKING:
    from .panel import ProteinPanel

#: Å beyond the focused pocket kept in the slab
SLAB_PAD = 2.0


class PocketCamera(QObject):
    def __init__(self, panel: ProteinPanel) -> None:
        super().__init__(panel)
        self.panel = panel
        panel.session.focusChanged.connect(self.focus_changed)

    def watch(self, view) -> None:
        view.cameraMoved.connect(self.apply)

    def focus_changed(self) -> None:
        """Frame a newly focused pocket; without a focus, draw everything again."""
        view = self.panel.view
        if view is None:
            return
        found = self.panel.session.focused_atoms()
        if found is None:
            view.set_slab(*NO_SLAB)
            view.set_fog(None)
            return
        shown, atoms = found
        positions = shown.xyz[atoms]
        view.frame(positions, np.full(len(positions), 1.5))
        self.apply()

    def apply(self) -> None:
        """Keep the slab and fog around the focused pocket where the camera now is."""
        view = self.panel.view
        found = self.panel.session.focused_atoms()
        if view is None or found is None:
            return
        shown, atoms = found
        positions = shown.xyz[atoms]
        centre = positions.mean(axis=0)
        radius = float(np.linalg.norm(positions - centre, axis=1).max()) + SLAB_PAD
        _xy, depth = view.camera.project(centre[None], view.width(), view.height())
        middle = float(depth[0])
        view.set_slab(max(middle - radius, 0.0), middle + radius)
        view.set_fog((middle, middle + radius))
