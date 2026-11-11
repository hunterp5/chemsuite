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
"""What every molecule view drawn through Qt's RHI shares: camera, mouse, frame, GPU helpers.

A view draws with the system's own graphics API (Metal on macOS, Direct3D 11 on Windows,
OpenGL on Linux). This base owns the camera and how the mouse moves it (as in 3Dmol.js: drag
to turn, right- or middle-drag or Shift-drag to pan, the wheel to zoom; a click without a drag
picks), fitting the scene to the view, the per-frame uniform block, and helpers to make
buffers and pipelines. A subclass says what it draws:

- :meth:`extent`: atom positions and radii, to fit and orient the camera;
- :meth:`atom_at`: the atom under a pixel, or -1;
- :meth:`create_pipelines`, :meth:`upload` and :meth:`draw`: its GPU work.

The uniform block every shader sees starts with ``view, proj, light, params, background``
(the 3D Viewer's shaders read only these) and goes on with ``slab`` (near and far view depth
kept, ghost fade) and ``highlight`` (the selection colour), which the protein shaders read.
"""

from __future__ import annotations

from importlib.resources import files

import numpy as np
from PySide6.QtCore import QPointF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QImage,
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
from shiboken6 import getCppPointer

from chemsuite.molview.camera import Camera

#: how far the pointer may move before a press becomes a drag (pixels)
CLICK_SLOP = 4
#: wheel steps to a doubling of the distance
WHEEL_STEPS = 6.0
#: the uniform block: view, projection, light, params, background, slab, highlight
FRAME_FLOATS = 16 + 16 + 4 + 4 + 4 + 4 + 4
_LIGHT = (0.35, 0.55, 1.0, 0.0)
#: no slab: everything from the eye to far away is drawn
NO_SLAB = (0.0, 1e9)
#: a quad's corners, and a box's (a 14-vertex strip covering every face)
QUAD = np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]], np.float32)
BOX = np.array(
    [
        [0, 1, 1], [1, 1, 1], [0, -1, 1], [1, -1, 1], [1, -1, -1], [1, 1, 1], [1, 1, -1],
        [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1], [1, -1, -1], [0, 1, -1], [1, 1, -1],
    ],
    np.float32,
)  # fmt: skip

FLOAT = QRhiVertexInputAttribute.Format.Float
FLOAT2 = QRhiVertexInputAttribute.Format.Float2
FLOAT3 = QRhiVertexInputAttribute.Format.Float3
FLOAT4 = QRhiVertexInputAttribute.Format.Float4


def shader(name: str, package: str = __package__) -> QShader:
    """A compiled shader from *package*'s ``shaders`` folder."""
    data = (files(package) / "shaders" / f"{name}.qsb").read_bytes()
    return QShader.fromSerialized(data)


def _matrix(qmatrix) -> np.ndarray:
    return np.array(qmatrix.data(), float).reshape(4, 4).T  # Qt stores columns first


def _columns(m: np.ndarray) -> np.ndarray:
    """A matrix as std140 wants it: columns first."""
    return np.asarray(m, np.float32).T.ravel()


class RhiView(QRhiWidget):
    atomClicked = Signal(int, object)  # atom index (as the subclass counts), modifiers
    atomHovered = Signal(int)  # atom index, -1 when none
    atomDoubleClicked = Signal(int)  # atom index, -1 when none
    contextRequested = Signal(int, object)  # a right-click: atom index, global QPoint
    cameraMoved = Signal()
    failed = Signal()  # the GPU could not be used

    #: the package whose ``shaders`` folder holds the shaders :meth:`pipeline` names
    shaders = __package__

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSampleCount(4)
        self.setMouseTracking(True)
        self.setMinimumSize(120, 120)
        self.camera = Camera()
        self.background = QColor(255, 255, 255)
        self.highlight = QColor(255, 196, 0)  # the selection colour
        self.slab = NO_SLAB
        self.fog: tuple[float, float] | None = None  # view depths; None: the scene's back half
        self.depth_cue = True  # fade with depth at all (fog)
        self.ghost = 0.75  # how far ghosted atoms fade into the background
        self._rhi_key = None  # the address of the QRhi the GPU resources belong to
        self._broken = False  # the GPU cannot run the shaders: draw nothing
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

    # -- for subclasses ---------------------------------------------------------------------
    def extent(self) -> tuple[np.ndarray, np.ndarray]:
        """Positions (n, 3) and radii (n,) of what is shown, to fit the camera to."""
        return np.zeros((0, 3)), np.zeros(0)

    def atom_at(self, x: float, y: float) -> int:
        return -1

    def create_pipelines(self) -> bool:
        """Make the pipelines (after :attr:`_srb` exists); False when the GPU cannot."""
        return True

    def release_buffers(self) -> None:
        """Forget GPU buffers: the QRhi they belong to is gone."""

    def upload(self, batch) -> None:
        """Add what changed since the last frame to ``batch``."""

    def draw(self, cb, viewport: QRhiViewport) -> None:
        """Record the draw calls (inside the pass)."""

    # -- camera -----------------------------------------------------------------------------
    @property
    def auto_fitting(self) -> bool:
        """True until the user turns, moves or zooms the fitted view."""
        return self._auto_fit

    def _aspect(self) -> float:
        return self.width() / max(self.height(), 1)

    def orient(self) -> None:
        """Turn to the principal axes of what is shown (again once the size is final)."""
        self.camera.orient(self.extent()[0], self._aspect())
        self._orient_pending = not self.isVisible()

    def reset_view(self) -> None:
        """Turn what is shown to its principal axes and fit it."""
        self.camera.orient(self.extent()[0], self._aspect())
        self.fit()

    def fit(self) -> None:
        """Fit what is shown to the view now, and again whenever the view is resized, until
        the user moves the camera."""
        self._fit_now()
        self._fit_pending = not self.isVisible()
        self._auto_fit = True

    def _fit_now(self) -> None:
        positions, radii = self.extent()
        self.camera.fit(positions, radii, self._aspect())
        self.camera_changed()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if self._orient_pending and self.isVisible():
            self._orient_pending = False
            self.camera.orient(self.extent()[0], self._aspect())
        if self._auto_fit or self._fit_pending:
            self._fit_pending = False
            self._fit_now()
        else:
            self.camera_changed()

    def frame(self, positions: np.ndarray, radii: np.ndarray | None = None) -> None:
        """Fit these points to the view (the user's choice: not refitted on resize)."""
        self.camera.fit(positions, radii, self._aspect())
        self.set_camera(self.camera)

    def set_camera(self, camera: Camera) -> None:
        self.camera = camera
        self._auto_fit = self._fit_pending = self._orient_pending = False
        self.camera_changed()

    def camera_changed(self) -> None:
        """The camera or the view's size changed: redraw and say so."""
        self.update()
        self.cameraMoved.emit()

    def set_background(self, colour: QColor) -> None:
        self.background = QColor(colour)
        self.update()

    def set_slab(self, near: float, far: float) -> None:
        """Draw only what lies between these view depths (Å in front of the eye)."""
        self.slab = (near, far)
        self.update()

    def set_fog(self, fog: tuple[float, float] | None) -> None:
        """Fade into the background from the first view depth to the second; None fades the
        back half of the scene (the default)."""
        self.fog = fog
        self.update()

    def set_depth_cue(self, on: bool) -> None:
        """Fade what lies further back into the background (fog), or draw everything fully."""
        self.depth_cue = on
        self.update()

    def picture(self, width: int, height: int, transparent: bool = False) -> QImage:
        """The scene drawn at *width* x *height* pixels (the view's own shape gives the same
        framing), with a see-through background when *transparent*. Null when it cannot be
        drawn (no GPU, or the view is not on screen)."""
        background = QColor(self.background)
        if transparent:
            self.background.setAlpha(0)
        self.setFixedColorBufferSize(QSize(width, height))
        try:
            image = self.grabFramebuffer()
        finally:
            self.setFixedColorBufferSize(QSize())
            self.background = background
            self.update()
        return image

    # -- QRhiWidget -------------------------------------------------------------------------
    def initialize(self, cb) -> None:
        rhi = self.rhi()
        key = getCppPointer(rhi)[0]
        if key == self._rhi_key:
            return
        self._rhi_key = key
        self.release_buffers()
        self._quad = self.static_buffer(QUAD)
        self._box = self.static_buffer(BOX)
        self._frame = rhi.newBuffer(
            QRhiBuffer.Type.Dynamic, QRhiBuffer.UsageFlag.UniformBuffer, FRAME_FLOATS * 4
        )
        self._frame.create()
        self._srb = rhi.newShaderResourceBindings()
        stages = (
            QRhiShaderResourceBinding.StageFlag.VertexStage
            | QRhiShaderResourceBinding.StageFlag.FragmentStage
        )
        self._srb.setBindings([QRhiShaderResourceBinding.uniformBuffer(0, stages, self._frame)])
        self._srb.create()
        self._broken = not (
            rhi.isFeatureSupported(QRhi.Feature.Instancing) and self.create_pipelines()
        )
        if self._broken:
            QTimer.singleShot(0, self.failed.emit)  # not from inside the frame
            return
        upload = rhi.nextResourceUpdateBatch()
        upload.uploadStaticBuffer(self._quad, QUAD.tobytes())
        upload.uploadStaticBuffer(self._box, BOX.tobytes())
        cb.resourceUpdate(upload)

    def releaseResources(self) -> None:
        self._rhi_key = None
        self.release_buffers()

    def render(self, cb) -> None:
        if self._broken:
            return
        rhi = self.rhi()
        target = self.renderTarget()
        size = target.pixelSize()
        width, height = size.width(), size.height()
        batch = rhi.nextResourceUpdateBatch()
        self.upload(batch)
        frame = self._frame_data(rhi, width / max(height, 1))
        batch.updateDynamicBuffer(self._frame, 0, frame.nbytes, frame.tobytes())
        cb.beginPass(target, self.background, QRhiDepthStencilClearValue(1.0, 0), batch)
        self.draw(cb, QRhiViewport(0, 0, width, height))
        cb.endPass()

    # -- GPU helpers ------------------------------------------------------------------------
    def static_buffer(self, data: np.ndarray):
        buf = self.rhi().newBuffer(
            QRhiBuffer.Type.Immutable, QRhiBuffer.UsageFlag.VertexBuffer, data.nbytes
        )
        buf.create()
        return buf

    def instances(self, old, data: np.ndarray, batch, usage=QRhiBuffer.UsageFlag.VertexBuffer):
        """A buffer holding ``data`` (a new one when it no longer fits); ``old`` when empty."""
        if not len(data):
            return old
        if old is None or old.size() < data.nbytes:
            if old is not None:
                old.deleteLater()
            old = self.rhi().newBuffer(QRhiBuffer.Type.Static, usage, max(data.nbytes, 4096))
            old.create()
        batch.uploadStaticBuffer(old, 0, data.nbytes, np.ascontiguousarray(data).tobytes())
        return old

    def pipeline(
        self, name: str, strides, attributes, *, topology=None, per_vertex=False, build=True
    ):
        """A pipeline for ``name.vert``/``name.frag``; ``strides`` per binding (the first per
        vertex, the rest per instance unless ``per_vertex``); ``attributes`` as (binding,
        location, format, offset). ``build=False`` returns it uncreated, to adjust first.
        """
        pipe = self.rhi().newGraphicsPipeline()
        pipe.setShaderStages(
            [
                QRhiShaderStage(QRhiShaderStage.Type.Vertex, shader(f"{name}.vert", self.shaders)),
                QRhiShaderStage(
                    QRhiShaderStage.Type.Fragment, shader(f"{name}.frag", self.shaders)
                ),
            ]
        )
        layout = QRhiVertexInputLayout()
        rest = (
            QRhiVertexInputBinding.Classification.PerVertex
            if per_vertex
            else (QRhiVertexInputBinding.Classification.PerInstance)
        )
        layout.setBindings(
            [QRhiVertexInputBinding(strides[0])]
            + [QRhiVertexInputBinding(s, rest) for s in strides[1:]]
        )
        layout.setAttributes([QRhiVertexInputAttribute(*a) for a in attributes])
        pipe.setVertexInputLayout(layout)
        pipe.setTopology(topology or QRhiGraphicsPipeline.Topology.TriangleStrip)
        pipe.setDepthTest(True)
        pipe.setDepthWrite(True)
        pipe.setSampleCount(self.sampleCount())
        pipe.setShaderResourceBindings(self._srb)
        pipe.setRenderPassDescriptor(self.renderTarget().renderPassDescriptor())
        if not build:
            return pipe
        return pipe if pipe.create() else None

    def _frame_data(self, rhi, aspect: float) -> np.ndarray:
        camera = self.camera
        projection = _matrix(rhi.clipSpaceCorrMatrix()) @ camera.projection(aspect)
        frame = np.zeros(FRAME_FLOATS, np.float32)
        frame[0:16] = _columns(camera.view())
        frame[16:32] = _columns(projection)
        frame[32:36] = _LIGHT
        frame[36] = 1.0 if rhi.isClipDepthZeroToOne() else 0.0
        frame[37] = camera.distance  # fog from the centre of the scene back
        frame[38] = camera.distance + camera.radius * 1.2
        if self.fog is not None:
            frame[37:39] = self.fog
        if not self.depth_cue:
            frame[37:39] = (1e7, 1e7 + 1.0)  # fog starts far behind everything
        frame[39] = 1.0 if camera.orthographic else 0.0
        frame[40:44] = self.background.getRgbF()
        frame[44:47] = (self.slab[0], self.slab[1], self.ghost)
        frame[48:52] = self.highlight.getRgbF()
        return frame

    # -- mouse ------------------------------------------------------------------------------
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
        self.camera_changed()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if self._press is not None and not self._dragged:
            pos = event.position()
            atom = self.atom_at(pos.x(), pos.y())
            if event.button() == Qt.MouseButton.RightButton:
                self.contextRequested.emit(atom, event.globalPosition().toPoint())
            else:
                self.atomClicked.emit(atom, event.modifiers())
        self._press = self._last = None
        self._dragged = False

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self.atomDoubleClicked.emit(self.atom_at(pos.x(), pos.y()))

    def wheelEvent(self, event: QWheelEvent) -> None:
        steps = event.angleDelta().y() / 120.0
        if steps:
            self._auto_fit = False
            self.camera.zoom(2.0 ** (-steps / WHEEL_STEPS))
            self.camera_changed()

    def leaveEvent(self, event) -> None:
        super().leaveEvent(event)
        if self._hover != -1:
            self._hover = -1
            self.atomHovered.emit(-1)

    @property
    def dragging(self) -> bool:
        return self._press is not None and self._dragged
