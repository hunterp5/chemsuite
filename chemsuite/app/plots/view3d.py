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
"""The 3D plots' view: points, lines and a surface drawn through Qt's RHI.

It draws with whatever graphics API the system uses (Direct3D 11 on Windows, falling back to
the WARP software rasteriser by itself; Metal on macOS; OpenGL on Linux), like the Structure
Viewer. Points are instanced squares the shader rounds; lines and the surface are plain
vertex lists. Axis labels and the box or lasso being drawn are painted over it
(:class:`_Overlay`), and also onto an exported picture.

The mouse: drag to turn (Box and Lasso modes draw a selection instead), right- or
middle-drag (or Shift-drag) to pan, the wheel to zoom.
"""

from __future__ import annotations

from importlib.resources import files

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRect, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QImage,
    QPainter,
    QPainterPath,
    QPen,
    QRhi,
    QRhiBuffer,
    QRhiCommandBuffer,
    QRhiDepthStencilClearValue,
    QRhiGraphicsPipeline,
    QRhiShaderResourceBinding,
    QRhiShaderStage,
    QRhiVertexInputAttribute,
    QRhiVertexInputBinding,
    QRhiVertexInputLayout,
    QRhiViewport,
    QShader,
)
from PySide6.QtWidgets import QRhiWidget, QWidget
from shiboken6 import getCppPointer

from chemsuite.plotting.orbit import OrbitCamera

from .view_box import BOX, LASSO, ZOOM

CLICK_SLOP = 4  # pixels the mouse may move and still count as a click
BACKGROUND = (255, 255, 255)
HIGHLIGHT = (230, 70, 20)
INK = (70, 70, 70)  # labels
_FRAME_FLOATS = 16 + 16 + 4 + 4  # mvp, view, viewport, light
_LIGHT = (0.3, 0.5, 1.0, 0.0)
_QUAD = np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]], np.float32)
_FLOAT2 = QRhiVertexInputAttribute.Format.Float2
_FLOAT3 = QRhiVertexInputAttribute.Format.Float3
_FLOAT4 = QRhiVertexInputAttribute.Format.Float4


def _shader(name: str) -> QShader:
    data = (files("chemsuite.app.molview") / "shaders" / f"{name}.qsb").read_bytes()
    return QShader.fromSerialized(data)


def _columns(m: np.ndarray) -> np.ndarray:
    return np.asarray(m, np.float32).T.ravel()  # std140 wants columns first


def _qmatrix(m) -> np.ndarray:
    return np.array(m.data(), float).reshape(4, 4).T


def points_array(pos: np.ndarray, colors: np.ndarray, sizes) -> np.ndarray:
    """Per-point instance data: x, y, z, diameter (points), r, g, b, a."""
    n = len(pos)
    out = np.empty((n, 8), np.float32)
    if n:
        out[:, :3] = pos
        out[:, 3] = sizes
        out[:, 4:] = colors
    return out


def lines_array(pos: np.ndarray, colors) -> np.ndarray:
    """Per-vertex line data (pairs of vertices): x, y, z, r, g, b, a."""
    n = len(pos)
    out = np.empty((n, 7), np.float32)
    if n:
        out[:, :3] = pos
        out[:, 3:] = colors
    return out


def surface_normals(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Smooth vertex normals: each vertex averages the faces around it."""
    normals = np.zeros_like(vertices, dtype=float)
    if len(faces):
        a, b, c = (vertices[faces[:, i]] for i in range(3))
        face = np.cross(b - a, c - a)
        for i in range(3):
            np.add.at(normals, faces[:, i], face)
    length = np.linalg.norm(normals, axis=1, keepdims=True)
    return np.divide(normals, length, out=np.tile([0.0, 0.0, 1.0], (len(vertices), 1)),
                     where=length > 1e-12)  # fmt: skip


class _Overlay(QWidget):
    """Axis labels and the box or lasso being drawn, painted over the GPU picture."""

    def __init__(self, view: View3D) -> None:
        super().__init__(view)
        self.view = view
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.labels: list[tuple[np.ndarray, str]] = []
        self.points: list[QPointF] = []
        self.box: QRect | None = None
        self.tracing = False

    def paintEvent(self, _event) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.paint_labels(painter, self.width(), self.height())
        if not self.tracing:
            return
        painter.setPen(QPen(QColor(*HIGHLIGHT), 1.2, Qt.PenStyle.DashLine))
        painter.setBrush(QColor(*HIGHLIGHT, 30))
        if self.box is not None:
            painter.drawRect(self.box)
        elif len(self.points) > 1:
            path = QPainterPath(self.points[0])
            for point in self.points[1:]:
                path.lineTo(point)
            path.closeSubpath()
            painter.drawPath(path)

    def paint_labels(self, painter: QPainter, width: int, height: int) -> None:
        if not self.labels:
            return
        font = QFont(self.font())
        font.setPointSizeF(max(8.0, font.pointSizeF() - 2))
        painter.setFont(font)
        painter.setPen(QColor(*INK))
        xy, front = self.view.camera.project(np.array([p for p, _ in self.labels]), width, height)
        metrics = painter.fontMetrics()
        for (x, y), ok, (_pos, text) in zip(xy, front, self.labels, strict=True):
            if ok:
                w = metrics.horizontalAdvance(text)
                painter.drawText(QPointF(x - w / 2, y + metrics.ascent() / 2), text)


class View3D(QRhiWidget):
    """Draws a 3D plot; reports clicks, boxes, lassos, hovers and camera moves (so the panel
    can select what you see and move the cards beside selected points)."""

    clicked = Signal(QPoint, object)
    boxed = Signal(QRect, object)
    lassoed = Signal(object, object)  # (n × 2 widget pixels, modifiers)
    hovered = Signal(QPoint)
    turning = Signal()  # a drag started
    moved = Signal()  # the camera moved
    failed = Signal()  # the graphics cannot draw this

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSampleCount(4)
        self.setMouseTracking(True)
        self.setMinimumSize(120, 120)
        self.mode = ZOOM
        self.camera = OrbitCamera()
        self.overlay = _Overlay(self)
        self._points = points_array(np.empty((0, 3)), np.empty((0, 4)), 0)
        self._chosen = self._points
        self._lines = lines_array(np.empty((0, 3)), np.empty((0, 4)))
        self._surface = np.empty((0, 10), np.float32)
        self._faces = np.empty(0, np.uint32)
        self._dirty = {"points", "chosen", "lines", "surface"}
        self._rhi_key = None  # the address of the QRhi the GPU resources belong to
        self._broken = False
        self._buffers: dict[str, object] = {}
        self._counts: dict[str, int] = {}
        self._press: QPoint | None = None
        self._last: QPointF | None = None
        self._dragging = False
        self.renderFailed.connect(self.failed)

    # -- what to draw ---------------------------------------------------------------------------
    def set_points(self, pos: np.ndarray, colors: np.ndarray, sizes) -> None:
        self._points = points_array(pos, colors, sizes)
        self._changed("points")

    def set_chosen(self, pos: np.ndarray, colors: np.ndarray, sizes) -> None:
        """The selected points, drawn over the others."""
        self._chosen = points_array(pos, colors, sizes)
        self._changed("chosen")

    def set_lines(self, pos: np.ndarray, colors) -> None:
        self._lines = lines_array(pos, colors)
        self._changed("lines")

    def set_surface(self, vertices: np.ndarray, colors: np.ndarray, faces: np.ndarray) -> None:
        data = np.empty((len(vertices), 10), np.float32)
        if len(vertices):
            data[:, :3] = vertices
            data[:, 3:6] = surface_normals(vertices, faces)
            data[:, 6:] = colors
        self._surface = data
        self._faces = np.ascontiguousarray(faces, np.uint32).ravel()
        self._changed("surface")

    @property
    def surface_triangles(self) -> int:
        return len(self._faces) // 3

    def set_labels(self, labels: list[tuple[np.ndarray, str]]) -> None:
        self.overlay.labels = labels
        self.overlay.update()

    def _changed(self, part: str) -> None:
        self._dirty.add(part)
        self.update()

    # -- the camera -----------------------------------------------------------------------------
    def camera_state(self) -> dict:
        return self.camera.state()

    def set_camera(self, state: object) -> None:
        self.camera = OrbitCamera.from_state(state)
        self._camera_moved()

    def reset_camera(self) -> None:
        self.camera = OrbitCamera()
        self._camera_moved()

    def _camera_moved(self) -> None:
        self.update()
        self.overlay.update()
        self.moved.emit()

    def export(self, path) -> bool:
        """Save the picture with its labels; False when the graphics cannot draw it."""
        image = self.grabFramebuffer()
        if image.isNull():
            return False
        image = image.convertToFormat(QImage.Format.Format_ARGB32)
        image.setDevicePixelRatio(self.devicePixelRatioF())
        painter = QPainter(image)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.overlay.paint_labels(painter, self.width(), self.height())
        painter.end()
        return image.save(str(path))

    # -- QRhiWidget -----------------------------------------------------------------------------
    def initialize(self, cb) -> None:
        rhi = self.rhi()
        key = getCppPointer(rhi)[0]
        if key == self._rhi_key:
            return
        self._rhi_key = key
        self._buffers = {}
        self._counts = {}
        self._dirty = {"points", "chosen", "lines", "surface"}
        self._quad = self._buffer(QRhiBuffer.Type.Immutable, _QUAD.nbytes)
        self._frame = self._buffer(
            QRhiBuffer.Type.Dynamic, _FRAME_FLOATS * 4, QRhiBuffer.UsageFlag.UniformBuffer
        )
        self._srb = rhi.newShaderResourceBindings()
        stages = (
            QRhiShaderResourceBinding.StageFlag.VertexStage
            | QRhiShaderResourceBinding.StageFlag.FragmentStage
        )
        self._srb.setBindings([QRhiShaderResourceBinding.uniformBuffer(0, stages, self._frame)])
        self._srb.create()
        topology = QRhiGraphicsPipeline.Topology
        self._pipes = {
            "points": self._pipeline(
                "plot_point",
                topology.TriangleStrip,
                [(8, False), (32, True)],
                [(0, 0, _FLOAT2, 0), (1, 1, _FLOAT4, 0), (1, 2, _FLOAT4, 16)],
                blend=True,
            ),
            "lines": self._pipeline(
                "plot_line",
                topology.Lines,
                [(28, False)],
                [(0, 0, _FLOAT3, 0), (0, 1, _FLOAT4, 12)],
                blend=True,
            ),
            "surface": self._pipeline(
                "plot_surface",
                topology.Triangles,
                [(40, False)],
                [(0, 0, _FLOAT3, 0), (0, 1, _FLOAT3, 12), (0, 2, _FLOAT4, 24)],
                blend=False,
            ),
        }
        self._broken = not rhi.isFeatureSupported(QRhi.Feature.Instancing) or None in (
            self._pipes.values()
        )
        if self._broken:
            QTimer.singleShot(0, self.failed.emit)  # not from inside the frame
            return
        batch = rhi.nextResourceUpdateBatch()
        batch.uploadStaticBuffer(self._quad, _QUAD.tobytes())
        cb.resourceUpdate(batch)

    def releaseResources(self) -> None:
        self._rhi_key = None
        self._buffers = {}

    def render(self, cb) -> None:
        if self._broken:
            return
        rhi = self.rhi()
        target = self.renderTarget()
        size = target.pixelSize()
        width, height = size.width(), size.height()
        batch = rhi.nextResourceUpdateBatch()
        data = {
            "points": self._points,
            "chosen": self._chosen,
            "lines": self._lines,
            "surface": self._surface,
        }
        for part in self._dirty:
            self._upload(part, data[part], batch)
            self._counts[part] = len(data[part])
        if "surface" in self._dirty:
            self._upload("faces", self._faces, batch, QRhiBuffer.UsageFlag.IndexBuffer)
            self._counts["faces"] = len(self._faces)
        self._dirty.clear()
        frame = self._frame_data(rhi, width, height)
        batch.updateDynamicBuffer(self._frame, 0, frame.nbytes, frame.tobytes())
        cb.beginPass(target, QColor(*BACKGROUND), QRhiDepthStencilClearValue(1.0, 0), batch)
        viewport = QRhiViewport(0, 0, width, height)

        def use(pipe: str) -> None:
            cb.setGraphicsPipeline(self._pipes[pipe])
            cb.setViewport(viewport)
            cb.setShaderResources()

        if self._counts.get("faces"):
            use("surface")
            cb.setVertexInput(
                0,
                [(self._buffers["surface"], 0)],
                self._buffers["faces"],
                0,
                QRhiCommandBuffer.IndexFormat.IndexUInt32,
            )
            cb.drawIndexed(self._counts["faces"])
        if self._counts.get("lines"):
            use("lines")
            cb.setVertexInput(0, [(self._buffers["lines"], 0)])
            cb.draw(self._counts["lines"])
        for part in ("points", "chosen"):
            if self._counts.get(part):
                use("points")
                cb.setVertexInput(0, [(self._quad, 0), (self._buffers[part], 0)])
                cb.draw(len(_QUAD), self._counts[part])
        cb.endPass()

    # -- GPU resources --------------------------------------------------------------------------
    def _buffer(self, kind, nbytes: int, usage=QRhiBuffer.UsageFlag.VertexBuffer):
        buf = self.rhi().newBuffer(kind, usage, nbytes)
        buf.create()
        return buf

    def _upload(self, part: str, data: np.ndarray, batch, usage=None) -> None:
        if not len(data):
            return
        old = self._buffers.get(part)
        if old is None or old.size() < data.nbytes:
            if old is not None:
                old.deleteLater()
            old = self._buffer(
                QRhiBuffer.Type.Static,
                max(data.nbytes, 4096),
                usage or QRhiBuffer.UsageFlag.VertexBuffer,
            )
            self._buffers[part] = old
        batch.uploadStaticBuffer(old, 0, data.nbytes, np.ascontiguousarray(data).tobytes())

    def _pipeline(self, name: str, topology, bindings, attributes, *, blend: bool):
        pipe = self.rhi().newGraphicsPipeline()
        pipe.setShaderStages(
            [
                QRhiShaderStage(QRhiShaderStage.Type.Vertex, _shader(f"{name}.vert")),
                QRhiShaderStage(QRhiShaderStage.Type.Fragment, _shader(f"{name}.frag")),
            ]
        )
        layout = QRhiVertexInputLayout()
        per_instance = QRhiVertexInputBinding.Classification.PerInstance
        layout.setBindings(
            [
                QRhiVertexInputBinding(stride, per_instance)
                if instanced
                else QRhiVertexInputBinding(stride)
                for stride, instanced in bindings
            ]
        )
        layout.setAttributes([QRhiVertexInputAttribute(*a) for a in attributes])
        pipe.setVertexInputLayout(layout)
        pipe.setTopology(topology)
        pipe.setDepthTest(True)
        pipe.setDepthWrite(not blend)  # see-through points and lines never hide each other
        if blend:
            factor = QRhiGraphicsPipeline.BlendFactor
            target = QRhiGraphicsPipeline.TargetBlend()
            target.enable = True
            target.srcColor = factor.SrcAlpha
            target.dstColor = factor.OneMinusSrcAlpha
            target.srcAlpha = factor.One
            target.dstAlpha = factor.OneMinusSrcAlpha
            pipe.setTargetBlends([target])
        pipe.setSampleCount(self.sampleCount())
        pipe.setShaderResourceBindings(self._srb)
        pipe.setRenderPassDescriptor(self.renderTarget().renderPassDescriptor())
        return pipe if pipe.create() else None

    def _frame_data(self, rhi, width: int, height: int) -> np.ndarray:
        projection = _qmatrix(rhi.clipSpaceCorrMatrix()) @ self.camera.projection(
            width / max(height, 1)
        )
        frame = np.zeros(_FRAME_FLOATS, np.float32)
        frame[0:16] = _columns(projection @ self.camera.view())
        frame[16:32] = _columns(self.camera.view())
        frame[32:36] = (width, height, self.devicePixelRatioF(), 0.0)
        frame[36:40] = _LIGHT
        return frame

    # -- mouse ----------------------------------------------------------------------------------
    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self.overlay.resize(self.size())
        self.moved.emit()

    def mousePressEvent(self, event) -> None:
        self._press = event.position().toPoint()
        self._last = event.position()
        self._dragging = False
        overlay = self.overlay
        if event.button() == Qt.MouseButton.LeftButton and self.mode in (BOX, LASSO):
            overlay.box = QRect(self._press, QSize()) if self.mode == BOX else None
            overlay.points = [QPointF(self._press)]
            overlay.tracing = True
            overlay.update()
            return
        self.turning.emit()

    def mouseMoveEvent(self, event) -> None:
        pos = event.position()
        if event.buttons() == Qt.MouseButton.NoButton or self._press is None:
            self.hovered.emit(pos.toPoint())
            return
        overlay = self.overlay
        if overlay.tracing:
            if overlay.box is not None:
                overlay.box = QRect(self._press, pos.toPoint()).normalized()
            else:
                overlay.points.append(QPointF(pos))
            overlay.update()
            return
        delta = pos - self._last
        self._last = pos
        if not self._dragging and (pos.toPoint() - self._press).manhattanLength() <= CLICK_SLOP:
            return
        self._dragging = True
        panning = event.buttons() & (Qt.MouseButton.RightButton | Qt.MouseButton.MiddleButton) or (
            event.modifiers() & Qt.KeyboardModifier.ShiftModifier
        )
        if panning:
            self.camera.pan(delta.x(), delta.y(), self.width())
        else:
            self.camera.turn(delta.x(), delta.y())
        self.update()
        self.overlay.update()

    def mouseReleaseEvent(self, event) -> None:
        pos = event.position().toPoint()
        press, self._press = self._press, None
        overlay = self.overlay
        if overlay.tracing:
            overlay.tracing = False
            overlay.update()
            box, points = overlay.box, overlay.points
            overlay.box, overlay.points = None, []
            if box is not None and (box.width() > CLICK_SLOP or box.height() > CLICK_SLOP):
                self.boxed.emit(box, event.modifiers())
                return
            if box is None and len(points) > 2:
                polygon = np.array([[p.x(), p.y()] for p in points])
                if (polygon.max(axis=0) - polygon.min(axis=0)).max() > CLICK_SLOP:
                    self.lassoed.emit(polygon, event.modifiers())
                    return
        if press is not None and (pos - press).manhattanLength() <= CLICK_SLOP:
            self.clicked.emit(pos, event.modifiers())
        if self._dragging:
            self._dragging = False
            self.moved.emit()

    def wheelEvent(self, event) -> None:
        delta = event.angleDelta()
        units = delta.x() or delta.y()
        if units:
            self.camera.zoom(units)
            self._camera_moved()
