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
"""Draw proteins with the GPU: atoms and bonds as impostors, in two buffers each.

The geometry (:class:`~chemsuite.protein.scene.Geometry`) is uploaded when atoms move or
styles change; the look (colour and flags per sphere and per bond end) when colours or the
selection change, and a hover rewrites only the rows of the atoms it touches. Nothing is
rebuilt on a pick: picking is a screen-space search (:mod:`chemsuite.protein.pick`) made again
only after the camera has moved.

Atom indices in and out (picks, hover) are the caller's: the ``atoms`` of the geometry.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import QWidget

from chemsuite.app.molview.rhi_view import BOX, FLOAT, FLOAT2, FLOAT3, FLOAT4, QUAD, RhiView
from chemsuite.protein.pick import ScreenPicker
from chemsuite.protein.scene import (
    CYLINDER_LOOK_FLOATS,
    HOVER,
    SPHERE_LOOK_FLOATS,
    Geometry,
)

_SPHERE_LOOK_BYTES = SPHERE_LOOK_FLOATS * 4
_CYLINDER_LOOK_BYTES = CYLINDER_LOOK_FLOATS * 4


class ProteinView(RhiView):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.geometry_ = Geometry.empty()
        self._sphere_look = np.zeros((0, SPHERE_LOOK_FLOATS), np.float32)
        self._cylinder_look = np.zeros((0, CYLINDER_LOOK_FLOATS), np.float32)
        self._hovered = np.zeros(0, np.int64)
        self._picker: ScreenPicker | None = None
        self._buffers: dict[str, object] = {}
        self._dirty: set[str] = set()
        self._rows: dict[str, list[int]] = {"sphere_look": [], "cylinder_look": []}
        self._counts = (0, 0)

    # -- what to draw -----------------------------------------------------------------------
    def set_scene(
        self,
        geometry: Geometry,
        sphere_look: np.ndarray,
        cylinder_look: np.ndarray,
        *,
        fit: bool = False,
        orient: bool = False,
    ) -> None:
        """Draw ``geometry`` with this look; ``fit`` fills the view with it, ``orient`` first
        turns it to its principal axes."""
        self.geometry_ = geometry
        self._hovered = np.zeros(0, np.int64)
        self._picker = None
        self._dirty |= {"spheres", "cylinders"}
        self.set_look(sphere_look, cylinder_look)
        if orient:
            self.orient()
        if fit or orient:
            self.fit()

    def set_look(self, sphere_look: np.ndarray, cylinder_look: np.ndarray) -> None:
        """New colours or flags for the same geometry: only the look is uploaded."""
        self._sphere_look = np.array(sphere_look, np.float32)
        self._cylinder_look = np.array(cylinder_look, np.float32)
        hovered, self._hovered = self._hovered, np.zeros(0, np.int64)
        self._dirty |= {"sphere_look", "cylinder_look"}
        self.set_hover(hovered)

    def set_hover(self, atoms: np.ndarray) -> None:
        """Light up these atoms (and their bonds' halves) as hovered; the last ones go out."""
        atoms = np.asarray(atoms, np.int64)
        self._flag(self._hovered, -HOVER)
        self._hovered = atoms
        self._flag(atoms, HOVER)
        self.update()

    def _flag(self, atoms: np.ndarray, change: float) -> None:
        if not len(atoms):
            return
        geometry = self.geometry_
        spheres = np.flatnonzero(np.isin(geometry.atoms, atoms))
        self._sphere_look[spheres, 4] += change
        self._rows["sphere_look"].extend(spheres.tolist())
        for end in (0, 1):
            cylinders = np.flatnonzero(np.isin(geometry.ends[:, end], atoms))
            self._cylinder_look[cylinders, 8 + end] += change
            self._rows["cylinder_look"].extend(cylinders.tolist())

    # -- RhiView ----------------------------------------------------------------------------
    def extent(self):
        spheres = self.geometry_.spheres
        return spheres[:, :3], spheres[:, 3]

    def camera_changed(self) -> None:
        self._picker = None
        super().camera_changed()

    def set_slab(self, near: float, far: float) -> None:
        self._picker = None
        super().set_slab(near, far)

    def atom_at(self, x: float, y: float) -> int:
        if self.dragging:
            return -1
        if self._picker is None:
            geometry = self.geometry_
            self._picker = ScreenPicker(
                self.camera,
                geometry.spheres[:, :3],
                geometry.pick_radii,
                self.width(),
                self.height(),
                self.slab,
            )
        hit = self._picker.pick(x, y)
        return -1 if hit < 0 else int(self.geometry_.atoms[hit])

    def create_pipelines(self) -> bool:
        self._sphere_pipe = self.pipeline(
            "pro_sphere",
            (8, 16, _SPHERE_LOOK_BYTES),
            [(0, 0, FLOAT2, 0), (1, 1, FLOAT4, 0), (2, 2, FLOAT4, 0), (2, 3, FLOAT, 16)],
        )
        self._cylinder_pipe = self.pipeline(
            "pro_cylinder",
            (12, 32, _CYLINDER_LOOK_BYTES),
            [
                (0, 0, FLOAT3, 0),
                (1, 1, FLOAT4, 0),
                (1, 2, FLOAT4, 16),
                (2, 3, FLOAT4, 0),
                (2, 4, FLOAT4, 16),
                (2, 5, FLOAT2, 32),
            ],
        )
        return self._sphere_pipe is not None and self._cylinder_pipe is not None

    def release_buffers(self) -> None:
        self._buffers = {}
        self._dirty = {"spheres", "cylinders", "sphere_look", "cylinder_look"}
        self._counts = (0, 0)

    def upload(self, batch) -> None:
        arrays = {
            "spheres": self.geometry_.spheres,
            "cylinders": self.geometry_.cylinders,
            "sphere_look": self._sphere_look,
            "cylinder_look": self._cylinder_look,
        }
        for name in self._dirty:
            self._buffers[name] = self.instances(self._buffers.get(name), arrays[name], batch)
            self._rows.get(name, []).clear()
        for name, rows in self._rows.items():
            buffer = self._buffers.get(name)
            if rows and buffer is not None:
                data = arrays[name]
                stride = data.shape[1] * 4
                for row in sorted(set(rows)):
                    batch.uploadStaticBuffer(buffer, row * stride, stride, data[row].tobytes())
            rows.clear()
        self._dirty.clear()
        self._counts = (len(self.geometry_.spheres), len(self.geometry_.cylinders))

    def draw(self, cb, viewport) -> None:
        spheres, cylinders = self._counts
        buffers = self._buffers
        if cylinders:
            cb.setGraphicsPipeline(self._cylinder_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(
                0, [(self._box, 0), (buffers["cylinders"], 0), (buffers["cylinder_look"], 0)]
            )
            cb.draw(len(BOX), cylinders)
        if spheres:
            cb.setGraphicsPipeline(self._sphere_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(
                0, [(self._quad, 0), (buffers["spheres"], 0), (buffers["sphere_look"], 0)]
            )
            cb.draw(len(QUAD), spheres)
