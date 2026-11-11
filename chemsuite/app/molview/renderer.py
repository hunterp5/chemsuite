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
"""Draw a :class:`~chemsuite.molview.geometry.Scene` with the GPU, through Qt's RHI.

It draws with the system's own graphics API (Metal on macOS, Direct3D 11 on Windows, OpenGL on
Linux), as the 3D plots do: every GPU widget in a window shares one API. Atoms and bonds
are *impostors*: each atom is one small square and each bond one box, and the shaders work
out where the eye's ray meets the true sphere or cylinder, so they are smooth at any zoom and
a large scene is two draw calls. Python only uploads the scene when it changes and four small
matrices per frame; nothing is redrawn until something changes.

The camera and mouse are :class:`~.rhi_view.RhiView`'s; a click without a drag picks the atom
under the pointer.
"""

from __future__ import annotations

from PySide6.QtWidgets import QWidget

from chemsuite.molview.geometry import Scene
from chemsuite.molview.pick import pick_sphere

from .rhi_view import BOX, FLOAT2, FLOAT3, FLOAT4, QUAD, RhiView


class MolRenderer(RhiView):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.scene = Scene.empty()
        self._spheres_buf = None
        self._cylinders_buf = None
        self._sizes = (0, 0)
        self._upload = True

    # -- what to draw ---------------------------------------------------------------------------
    def set_scene(self, scene: Scene, *, fit: bool = False, orient: bool = False) -> None:
        """Show *scene*; *fit* centres it and fills the view, *orient* first turns it to its
        principal axes (a new molecule)."""
        self.scene = scene
        if orient:
            self.orient()
        if fit or orient:
            self.fit()
        self._upload = True
        self.update()

    def extent(self):
        return self.scene.positions, self.scene.radii

    def atom_at(self, x: float, y: float) -> int:
        origin, direction = self.camera.ray(x, y, self.width(), self.height())
        hit = pick_sphere(origin, direction, self.scene.positions, self.scene.radii)
        return -1 if hit is None else hit

    # -- GPU ------------------------------------------------------------------------------------
    def create_pipelines(self) -> bool:
        self._sphere_pipe = self.pipeline(
            "sphere", (8, 32), [(0, 0, FLOAT2, 0), (1, 1, FLOAT4, 0), (1, 2, FLOAT4, 16)]
        )
        self._cylinder_pipe = self.pipeline(
            "cylinder",
            (12, 64),
            [
                (0, 0, FLOAT3, 0),
                (1, 1, FLOAT4, 0),
                (1, 2, FLOAT4, 16),
                (1, 3, FLOAT4, 32),
                (1, 4, FLOAT4, 48),
            ],
        )
        return self._sphere_pipe is not None and self._cylinder_pipe is not None

    def release_buffers(self) -> None:
        self._spheres_buf = self._cylinders_buf = None
        self._sizes = (0, 0)
        self._upload = True

    def upload(self, batch) -> None:
        if self._upload:
            self._spheres_buf = self.instances(self._spheres_buf, self.scene.spheres, batch)
            self._cylinders_buf = self.instances(self._cylinders_buf, self.scene.cylinders, batch)
            self._sizes = (len(self.scene.spheres), len(self.scene.cylinders))
            self._upload = False

    def draw(self, cb, viewport) -> None:
        spheres, cylinders = self._sizes
        if cylinders:
            cb.setGraphicsPipeline(self._cylinder_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(0, [(self._box, 0), (self._cylinders_buf, 0)])
            cb.draw(len(BOX), cylinders)
        if spheres:
            cb.setGraphicsPipeline(self._sphere_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(0, [(self._quad, 0), (self._spheres_buf, 0)])
            cb.draw(len(QUAD), spheres)
