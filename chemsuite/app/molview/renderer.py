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

The mouse works as in 3Dmol.js: drag to turn, right- or middle-drag (or Shift-drag) to pan,
the wheel to zoom. A click without a drag picks the atom under the pointer.
"""

from __future__ import annotations

from importlib.resources import files

import numpy as np
from PySide6.QtCore import QPointF, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QMouseEvent,
    QRhi,
    QRhiBuffer,
    QRhiDepthStencilClearValue,
    QRhiGraphicsPipeline,
    QRhiShaderResourceBinding,
    QRhiShaderStage,
    QRhiVertexInputAttribute,
    QRhiVertexInputBinding,
    QRhiVertexInputLayout,
    QRhiViewport,
    QShader,
    QWheelEvent,
)
from PySide6.QtWidgets import QRhiWidget, QWidget

from chemsuite.molview.camera import Camera
from chemsuite.molview.geometry import Scene
from chemsuite.molview.pick import pick_sphere

#: how far the pointer may move before a press becomes a drag (pixels)
CLICK_SLOP = 4
#: wheel steps to a doubling of the distance
WHEEL_STEPS = 6.0
#: the uniform block: view, projection, light, params, background
_FRAME_FLOATS = 16 + 16 + 4 + 4 + 4
_LIGHT = (0.35, 0.55, 1.0, 0.0)
#: a quad's corners, and a box's (a 14-vertex strip covering every face)
_QUAD = np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]], np.float32)
_BOX = np.array(
    [
        [0, 1, 1], [1, 1, 1], [0, -1, 1], [1, -1, 1], [1, -1, -1], [1, 1, 1], [1, 1, -1],
        [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1], [1, -1, -1], [0, 1, -1], [1, 1, -1],
    ],
    np.float32,
)  # fmt: skip


def _shader(name: str) -> QShader:
    data = (files("chemsuite.app.molview") / "shaders" / f"{name}.qsb").read_bytes()
    return QShader.fromSerialized(data)


def _matrix(qmatrix) -> np.ndarray:
    return np.array(qmatrix.data(), float).reshape(4, 4).T  # Qt stores columns first


def _columns(m: np.ndarray) -> np.ndarray:
    """A matrix as std140 wants it: columns first."""
    return np.asarray(m, np.float32).T.ravel()


class MolRenderer(QRhiWidget):
    atomClicked = Signal(int, object)  # scene atom index, keyboard modifiers
    atomHovered = Signal(int)  # scene atom index, -1 when none
    cameraMoved = Signal()
    failed = Signal()  # the GPU could not be used

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSampleCount(4)
        self.setMouseTracking(True)
        self.setMinimumSize(120, 120)
        self.camera = Camera()
        self.scene = Scene.empty()
        self.background = QColor(255, 255, 255)
        self._rhi = None  # the QRhi the GPU resources below belong to
        self._broken = False  # the GPU cannot run the shaders: draw nothing
        self._spheres_buf = None
        self._cylinders_buf = None
        self._sizes = (0, 0)
        self._upload = True
        self._press: QPointF | None = None
        self._last: QPointF | None = None
        self._dragged = False
        self._hover = -1
        #: fit at the next frame (the size is only final then), and keep fitting on resize
        #: until the view is turned, panned or zoomed
        self._fit_pending = False
        self._auto_fit = False
        self._orient_pending = False
        self.renderFailed.connect(self.failed)

    # -- what to draw ---------------------------------------------------------------------------
    def set_scene(self, scene: Scene, *, fit: bool = False, orient: bool = False) -> None:
        """Show *scene*; *fit* centres it and fills the view, *orient* first turns it to its
        principal axes (a new molecule)."""
        self.scene = scene
        if orient:
            self.camera.orient(scene.positions, self.width() / max(self.height(), 1))
            self._orient_pending = not self.isVisible()  # again once the size is final
        if fit or orient:
            self.fit()
        self._upload = True
        self.update()

    @property
    def auto_fitting(self) -> bool:
        """True until the user turns, moves or zooms the fitted view."""
        return self._auto_fit

    def reset_view(self) -> None:
        """Turn the scene to its principal axes and fit it."""
        self.camera.orient(self.scene.positions, self.width() / max(self.height(), 1))
        self.fit()

    def fit(self) -> None:
        """Fit the scene to the view now, and again whenever the view is resized, until the
        user moves the camera."""
        self._fit_now()
        self._fit_pending = not self.isVisible()
        self._auto_fit = True

    def _fit_now(self) -> None:
        aspect = self.width() / max(self.height(), 1)
        self.camera.fit(self.scene.positions, self.scene.radii, aspect)
        self.update()
        self.cameraMoved.emit()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if self._orient_pending and self.isVisible():
            self._orient_pending = False
            self.camera.orient(self.scene.positions, self.width() / max(self.height(), 1))
        if self._auto_fit or self._fit_pending:
            self._fit_pending = False
            self._fit_now()

    def set_camera(self, camera: Camera) -> None:
        self.camera = camera
        self._auto_fit = self._fit_pending = self._orient_pending = False
        self.update()
        self.cameraMoved.emit()

    def set_background(self, colour: QColor) -> None:
        self.background = QColor(colour)
        self.update()

    # -- QRhiWidget -----------------------------------------------------------------------------
    def initialize(self, cb) -> None:
        rhi = self.rhi()
        if rhi is self._rhi:
            return
        self._rhi = rhi
        self._spheres_buf = self._cylinders_buf = None
        self._sizes = (0, 0)
        self._upload = True
        self._quad = self._static(_QUAD)
        self._box = self._static(_BOX)
        self._frame = rhi.newBuffer(
            QRhiBuffer.Type.Dynamic, QRhiBuffer.UsageFlag.UniformBuffer, _FRAME_FLOATS * 4
        )
        self._frame.create()
        self._srb = rhi.newShaderResourceBindings()
        stages = (
            QRhiShaderResourceBinding.StageFlag.VertexStage
            | QRhiShaderResourceBinding.StageFlag.FragmentStage
        )
        self._srb.setBindings([QRhiShaderResourceBinding.uniformBuffer(0, stages, self._frame)])
        self._srb.create()
        float2 = QRhiVertexInputAttribute.Format.Float2
        float3 = QRhiVertexInputAttribute.Format.Float3
        float4 = QRhiVertexInputAttribute.Format.Float4
        self._sphere_pipe = self._pipeline(
            "sphere", 8, 32, [(0, 0, float2, 0), (1, 1, float4, 0), (1, 2, float4, 16)]
        )
        self._cylinder_pipe = self._pipeline(
            "cylinder",
            12,
            64,
            [
                (0, 0, float3, 0),
                (1, 1, float4, 0),
                (1, 2, float4, 16),
                (1, 3, float4, 32),
                (1, 4, float4, 48),
            ],
        )
        self._broken = not (
            rhi.isFeatureSupported(QRhi.Feature.Instancing)
            and self._sphere_pipe is not None
            and self._cylinder_pipe is not None
        )
        if self._broken:
            QTimer.singleShot(0, self.failed.emit)  # not from inside the frame
            return
        upload = rhi.nextResourceUpdateBatch()
        upload.uploadStaticBuffer(self._quad, _QUAD.tobytes())
        upload.uploadStaticBuffer(self._box, _BOX.tobytes())
        cb.resourceUpdate(upload)

    def releaseResources(self) -> None:
        self._rhi = None
        self._spheres_buf = self._cylinders_buf = None

    def render(self, cb) -> None:
        if self._broken:
            return
        rhi = self.rhi()
        target = self.renderTarget()
        size = target.pixelSize()
        width, height = size.width(), size.height()
        batch = rhi.nextResourceUpdateBatch()
        if self._upload:
            self._spheres_buf = self._instances(self._spheres_buf, self.scene.spheres, batch)
            self._cylinders_buf = self._instances(self._cylinders_buf, self.scene.cylinders, batch)
            self._sizes = (len(self.scene.spheres), len(self.scene.cylinders))
            self._upload = False
        frame = self._frame_data(rhi, width / max(height, 1))
        batch.updateDynamicBuffer(self._frame, 0, frame.nbytes, frame.tobytes())
        cb.beginPass(target, self.background, QRhiDepthStencilClearValue(1.0, 0), batch)
        viewport = QRhiViewport(0, 0, width, height)
        spheres, cylinders = self._sizes
        if cylinders:
            cb.setGraphicsPipeline(self._cylinder_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(0, [(self._box, 0), (self._cylinders_buf, 0)])
            cb.draw(len(_BOX), cylinders)
        if spheres:
            cb.setGraphicsPipeline(self._sphere_pipe)
            cb.setViewport(viewport)
            cb.setShaderResources()
            cb.setVertexInput(0, [(self._quad, 0), (self._spheres_buf, 0)])
            cb.draw(len(_QUAD), spheres)
        cb.endPass()

    # -- GPU resources --------------------------------------------------------------------------
    def _static(self, data: np.ndarray):
        buf = self._rhi.newBuffer(
            QRhiBuffer.Type.Immutable, QRhiBuffer.UsageFlag.VertexBuffer, data.nbytes
        )
        buf.create()
        return buf

    def _instances(self, old, data: np.ndarray, batch):
        """A buffer holding *data* (a new one when it no longer fits), or None when empty."""
        if not len(data):
            return old
        if old is None or old.size() < data.nbytes:
            if old is not None:
                old.deleteLater()
            old = self._rhi.newBuffer(
                QRhiBuffer.Type.Static, QRhiBuffer.UsageFlag.VertexBuffer, max(data.nbytes, 4096)
            )
            old.create()
        batch.uploadStaticBuffer(old, 0, data.nbytes, np.ascontiguousarray(data).tobytes())
        return old

    def _pipeline(self, name: str, stride: int, instance_stride: int, attributes):
        pipe = self._rhi.newGraphicsPipeline()
        pipe.setShaderStages(
            [
                QRhiShaderStage(QRhiShaderStage.Type.Vertex, _shader(f"{name}.vert")),
                QRhiShaderStage(QRhiShaderStage.Type.Fragment, _shader(f"{name}.frag")),
            ]
        )
        layout = QRhiVertexInputLayout()
        layout.setBindings(
            [
                QRhiVertexInputBinding(stride),
                QRhiVertexInputBinding(
                    instance_stride, QRhiVertexInputBinding.Classification.PerInstance
                ),
            ]
        )
        layout.setAttributes([QRhiVertexInputAttribute(*a) for a in attributes])
        pipe.setVertexInputLayout(layout)
        pipe.setTopology(QRhiGraphicsPipeline.Topology.TriangleStrip)
        pipe.setDepthTest(True)
        pipe.setDepthWrite(True)
        pipe.setSampleCount(self.sampleCount())
        pipe.setShaderResourceBindings(self._srb)
        pipe.setRenderPassDescriptor(self.renderTarget().renderPassDescriptor())
        return pipe if pipe.create() else None

    def _frame_data(self, rhi, aspect: float) -> np.ndarray:
        camera = self.camera
        projection = _matrix(rhi.clipSpaceCorrMatrix()) @ camera.projection(aspect)
        frame = np.zeros(_FRAME_FLOATS, np.float32)
        frame[0:16] = _columns(camera.view())
        frame[16:32] = _columns(projection)
        frame[32:36] = _LIGHT
        frame[36] = 1.0 if rhi.isClipDepthZeroToOne() else 0.0
        frame[37] = camera.distance  # fog from the centre of the scene back
        frame[38] = camera.distance + camera.radius * 1.2
        frame[40:44] = self.background.getRgbF()
        return frame

    # -- mouse ----------------------------------------------------------------------------------
    def atom_at(self, x: float, y: float) -> int:
        origin, direction = self.camera.ray(x, y, self.width(), self.height())
        hit = pick_sphere(origin, direction, self.scene.positions, self.scene.radii)
        return -1 if hit is None else hit

    def mousePressEvent(self, event: QMouseEvent) -> None:
        self._press = self._last = event.position()
        self._dragged = False

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        pos = event.position()
        if self._press is None:
            hover = self.atom_at(pos.x(), pos.y())
            if hover != self._hover:
                self._hover = hover
                self.atomHovered.emit(hover)
            return
        if (pos - self._press).manhattanLength() > CLICK_SLOP:
            self._dragged = True
        delta = pos - self._last
        self._last = pos
        if not self._dragged:
            return
        buttons = event.buttons()
        panning = buttons & (Qt.MouseButton.RightButton | Qt.MouseButton.MiddleButton) or (
            event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        )
        self._auto_fit = False
        if panning:
            self.camera.pan(delta.x(), delta.y(), self.height())
        else:
            self.camera.turn(delta.x(), delta.y())
        self.update()
        self.cameraMoved.emit()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if self._press is not None and not self._dragged:
            pos = event.position()
            self.atomClicked.emit(self.atom_at(pos.x(), pos.y()), event.modifiers())
        self._press = self._last = None
        self._dragged = False

    def wheelEvent(self, event: QWheelEvent) -> None:
        steps = event.angleDelta().y() / 120.0
        if steps:
            self._auto_fit = False
            self.camera.zoom(2.0 ** (-steps / WHEEL_STEPS))
            self.update()
            self.cameraMoved.emit()

    def leaveEvent(self, event) -> None:
        super().leaveEvent(event)
        if self._hover != -1:
            self._hover = -1
            self.atomHovered.emit(-1)
