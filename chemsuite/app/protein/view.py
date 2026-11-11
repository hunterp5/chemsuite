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

Meshes (cartoon, surface) are layers by name, drawn after the atoms; each vertex has an owner
(a residue or an atom) and a look row like the atoms'. A transparent layer draws without
writing depth, after the opaque ones. Meshes use 32-bit indices; where the GPU lacks them
(OpenGL ES 2 without the extension) each mesh is drawn in chunks small enough for 16-bit ones
(:func:`~chemsuite.protein.mesh.chunks`), its look uploaded whole when it changes.

Atom indices in and out (picks, hover) are the caller's: the ``atoms`` of the geometry, and
``set_pick_points`` adds points that pick an atom (a cartoon's CA atoms).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from PySide6.QtGui import QRhi, QRhiBuffer, QRhiCommandBuffer, QRhiGraphicsPipeline
from PySide6.QtWidgets import QWidget

from chemsuite.app.render3d.rhi_view import BOX, FLOAT, FLOAT2, FLOAT3, FLOAT4, QUAD, RhiView
from chemsuite.protein import mesh as meshes
from chemsuite.protein.pick import ScreenPicker
from chemsuite.protein.scene import (
    CYLINDER_LOOK_FLOATS,
    HOVER,
    SPHERE_LOOK_FLOATS,
    Geometry,
)

_SPHERE_LOOK_BYTES = SPHERE_LOOK_FLOATS * 4
_CYLINDER_LOOK_BYTES = CYLINDER_LOOK_FLOATS * 4
RESIDUE = "residue"
ATOM = "atom"
_MESH_MODES = ("opaque", "depth", "blend")


@dataclass(eq=False)
class _Layer:
    mesh: meshes.Mesh
    look: np.ndarray
    owner: str  # RESIDUE or ATOM
    transparent: bool
    buffers: dict = field(default_factory=dict)
    dirty: set = field(default_factory=lambda: {"vertices", "indices", "look"})
    rows: list = field(default_factory=list)
    chunks: list | None = None  # without 32-bit indices: the mesh's chunks, made when drawn
    parts: list = field(default_factory=list)  # and each chunk's buffers

    @property
    def uploaded(self) -> bool:
        return "indices" in self.buffers or bool(self.parts)


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
        self._layers: dict[str, _Layer] = {}
        self._hovered_residues = np.zeros(0, np.int64)
        self._picks = (np.zeros(0, np.int64), np.zeros((0, 3), np.float32), np.zeros(0))
        self._wide_indices = True

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

    def set_hover(self, atoms: np.ndarray, residues: np.ndarray | None = None) -> None:
        """Light up these atoms (their bonds' halves and the mesh parts they own) and these
        residues' cartoon as hovered; the last ones go out."""
        atoms = np.asarray(atoms, np.int64)
        residues = np.zeros(0, np.int64) if residues is None else np.asarray(residues, np.int64)
        self._flag(self._hovered, -HOVER)
        self._flag_layers(self._hovered, self._hovered_residues, -HOVER)
        self._hovered, self._hovered_residues = atoms, residues
        self._flag(atoms, HOVER)
        self._flag_layers(atoms, residues, HOVER)
        self.update()

    def _flag_layers(self, atoms: np.ndarray, residues: np.ndarray, change: float) -> None:
        for layer in self._layers.values():
            owners = residues if layer.owner == RESIDUE else atoms
            if not len(owners):
                continue
            rows = np.flatnonzero(np.isin(layer.mesh.owner, owners))
            layer.look[rows, 4] += change
            layer.rows.extend(rows.tolist())

    # meshes
    def set_mesh(
        self,
        name: str,
        mesh: meshes.Mesh,
        look: np.ndarray,
        *,
        owner: str = RESIDUE,
        transparent: bool = False,
    ) -> None:
        """Draw ``mesh`` as layer ``name`` (replacing one of that name)."""
        old = self._layers.pop(name, None)
        if old is not None:
            self._release_layer(old)
        if len(mesh):
            self._layers[name] = _Layer(mesh, np.array(look, np.float32), owner, transparent)
            self._flag_layers(self._hovered, self._hovered_residues, HOVER)
        self.update()

    def set_mesh_look(self, name: str, look: np.ndarray) -> None:
        layer = self._layers.get(name)
        if layer is None:
            return
        layer.look = np.array(look, np.float32)
        layer.dirty.add("look")
        layer.rows.clear()
        hovered = self._hovered_residues if layer.owner == RESIDUE else self._hovered
        rows = np.flatnonzero(np.isin(layer.mesh.owner, hovered))
        layer.look[rows, 4] += HOVER
        self.update()

    def remove_mesh(self, name: str) -> None:
        self.set_mesh(name, meshes.Mesh.empty(), np.zeros((0, meshes.LOOK_FLOATS)))

    @property
    def mesh_names(self) -> list[str]:
        return list(self._layers)

    def _release_layer(self, layer: _Layer) -> None:
        for buffers in (layer.buffers, *layer.parts):
            for buffer in buffers.values():
                buffer.deleteLater()
        layer.buffers.clear()
        layer.parts.clear()

    def set_pick_points(self, atoms: np.ndarray, centres: np.ndarray, radii: np.ndarray) -> None:
        """Points that pick an atom besides its sphere: a cartoon's trace atoms."""
        self._picks = (np.asarray(atoms, np.int64), np.asarray(centres, np.float32),
                       np.asarray(radii, np.float64))  # fmt: skip
        self._picker = None

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
        _, centres, radii = self._picks
        return (np.concatenate([spheres[:, :3], centres]),
                np.concatenate([spheres[:, 3], radii]))  # fmt: skip

    def camera_changed(self) -> None:
        self._picker = None
        super().camera_changed()

    def set_slab(self, near: float, far: float) -> None:
        self._picker = None
        super().set_slab(near, far)

    def atom_at(self, x: float, y: float) -> int:
        if self.dragging:
            return -1
        geometry = self.geometry_
        atoms, centres, radii = self._picks
        if self._picker is None:
            self._picker = ScreenPicker(
                self.camera,
                np.concatenate([geometry.spheres[:, :3], centres]),
                np.concatenate([geometry.pick_radii, radii]),
                self.width(),
                self.height(),
                self.slab,
            )
        hit = self._picker.pick(x, y)
        if hit < 0:
            return -1
        spheres = len(geometry.atoms)
        return int(geometry.atoms[hit]) if hit < spheres else int(atoms[hit - spheres])

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
        self._wide_indices = self.rhi().isFeatureSupported(QRhi.Feature.ElementIndexUint)
        self._mesh_pipes = {mode: self._mesh_pipeline(mode) for mode in _MESH_MODES}
        return (
            self._sphere_pipe is not None
            and self._cylinder_pipe is not None
            and None not in self._mesh_pipes.values()
        )

    def _mesh_pipeline(self, mode: str):
        """``opaque``; or for a transparent layer ``depth`` (its nearest depth only, no
        colour) then ``blend`` (only its nearest layer, blended), so it never speckles."""
        triangles = QRhiGraphicsPipeline.Topology.Triangles
        pipe = self.pipeline(
            "pro_mesh",
            (meshes.VERTEX_FLOATS * 4, meshes.LOOK_FLOATS * 4),
            [(0, 0, FLOAT3, 0), (0, 1, FLOAT3, 12), (1, 2, FLOAT4, 0), (1, 3, FLOAT, 16)],
            topology=triangles,
            per_vertex=True,
            build=False,
        )
        target = QRhiGraphicsPipeline.TargetBlend()
        if mode == "depth":
            target.colorWrite = QRhiGraphicsPipeline.ColorMaskComponent(0)
            pipe.setTargetBlends([target])
        elif mode == "blend":
            target.enable = True
            target.srcColor = QRhiGraphicsPipeline.BlendFactor.SrcAlpha
            target.dstColor = QRhiGraphicsPipeline.BlendFactor.OneMinusSrcAlpha
            target.srcAlpha = QRhiGraphicsPipeline.BlendFactor.One
            target.dstAlpha = QRhiGraphicsPipeline.BlendFactor.OneMinusSrcAlpha
            pipe.setTargetBlends([target])
            pipe.setDepthOp(QRhiGraphicsPipeline.CompareOp.LessOrEqual)
            pipe.setDepthWrite(False)
        return pipe if pipe.create() else None

    def release_buffers(self) -> None:
        self._buffers = {}
        self._dirty = {"spheres", "cylinders", "sphere_look", "cylinder_look"}
        self._counts = (0, 0)
        for layer in self._layers.values():
            layer.buffers.clear()
            layer.parts.clear()
            layer.dirty = {"vertices", "indices", "look"}
            layer.rows.clear()

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
        for layer in self._layers.values():
            if self._wide_indices:
                self._upload_layer(layer, batch)
            else:
                self._upload_chunks(layer, batch)

    def _upload_chunks(self, layer: _Layer, batch) -> None:
        """A layer in chunks with 16-bit indices: geometry once, the look whole on a change."""
        if layer.chunks is None:
            layer.chunks = meshes.chunks(layer.mesh)
        vertex = QRhiBuffer.UsageFlag.VertexBuffer
        index = QRhiBuffer.UsageFlag.IndexBuffer
        restyled = "look" in layer.dirty or bool(layer.rows)
        for i, chunk in enumerate(layer.chunks):
            if i == len(layer.parts):
                layer.parts.append({})
            part = layer.parts[i]
            if "indices" not in part:
                vertices = layer.mesh.vertices[chunk.vertex_map]
                part["vertices"] = self.instances(None, vertices, batch, usage=vertex)
                part["indices"] = self.instances(None, chunk.triangles, batch, usage=index)
            if restyled or "look" not in part:
                look = layer.look[chunk.vertex_map]
                part["look"] = self.instances(part.get("look"), look, batch, usage=vertex)
        layer.rows.clear()
        layer.dirty.clear()

    def _upload_layer(self, layer: _Layer, batch) -> None:
        arrays = {
            "vertices": layer.mesh.vertices,
            "indices": layer.mesh.triangles,
            "look": layer.look,
        }
        for name in layer.dirty:
            usage = (
                QRhiBuffer.UsageFlag.IndexBuffer
                if name == "indices"
                else QRhiBuffer.UsageFlag.VertexBuffer
            )
            layer.buffers[name] = self.instances(layer.buffers.get(name), arrays[name], batch,
                                                 usage=usage)  # fmt: skip
        if "look" not in layer.dirty and layer.rows and "look" in layer.buffers:
            rows = np.unique(layer.rows)
            # contiguous runs (a residue's vertices are together) in one upload each
            breaks = np.flatnonzero(np.diff(rows) != 1) + 1
            stride = meshes.LOOK_FLOATS * 4
            for run in np.split(rows, breaks):
                first, last = int(run[0]), int(run[-1]) + 1
                data = layer.look[first:last].tobytes()
                batch.uploadStaticBuffer(layer.buffers["look"], first * stride, len(data), data)
        layer.rows.clear()
        layer.dirty.clear()

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
        layers = [layer for layer in self._layers.values() if layer.uploaded]
        for layer in layers:
            if not layer.transparent:
                self._draw_mesh(cb, viewport, layer, "opaque")
        for layer in layers:
            if layer.transparent:
                self._draw_mesh(cb, viewport, layer, "depth")
                self._draw_mesh(cb, viewport, layer, "blend")

    def _draw_mesh(self, cb, viewport, layer: _Layer, mode: str) -> None:
        cb.setGraphicsPipeline(self._mesh_pipes[mode])
        cb.setViewport(viewport)
        cb.setShaderResources()
        if not layer.parts:
            self._draw_indexed(cb, layer.buffers, layer.mesh.triangles.size, wide=True)
        for part, chunk in zip(layer.parts, layer.chunks or (), strict=False):
            self._draw_indexed(cb, part, chunk.triangles.size, wide=False)

    @staticmethod
    def _draw_indexed(cb, buffers: dict, count: int, *, wide: bool) -> None:
        width = QRhiCommandBuffer.IndexFormat
        cb.setVertexInput(
            0,
            [(buffers["vertices"], 0), (buffers["look"], 0)],
            buffers["indices"],
            0,
            width.IndexUInt32 if wide else width.IndexUInt16,
        )
        cb.drawIndexed(count)
