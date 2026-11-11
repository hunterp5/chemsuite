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
"""The protein view: what it uploads, and (with a GPU) what it draws."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from chemsuite.protein import colours, io, scene
from chemsuite.protein.select import select

ROOT = Path(__file__).resolve().parents[2]
KINASE = ROOT / "tests" / "data" / "protein" / "1atp.cif.gz"


def _scene(structure, selected=None):
    geometry = scene.geometry(structure, scene.default_styles(structure))
    rgba = colours.colour(structure)
    flags = scene.flags(structure.n_atoms, selected=selected)
    return (
        geometry,
        scene.sphere_look(geometry, rgba, flags),
        scene.cylinder_look(geometry, rgba, flags),
    )


def test_hover_rewrites_only_the_rows_of_the_atoms_it_touches(qtbot) -> None:
    from chemsuite.app.protein.view import ProteinView

    structure = io.read(KINASE)
    view = ProteinView()
    qtbot.addWidget(view)
    geometry, spheres, cylinders = _scene(structure)
    view.set_scene(geometry, spheres, cylinders)
    view._rows["sphere_look"].clear()
    view._rows["cylinder_look"].clear()
    ligand = np.flatnonzero(structure.atom_kind == 2)
    view.set_hover(ligand)
    look = view._sphere_look
    hovered = np.isin(geometry.atoms, ligand)
    assert (look[hovered, 4] == scene.HOVER).all() and (look[~hovered, 4] == 0).all()
    assert sorted(set(view._rows["sphere_look"])) == np.flatnonzero(hovered).tolist()
    touched = np.isin(geometry.ends, ligand).any(axis=1)
    assert sorted(set(view._rows["cylinder_look"])) == np.flatnonzero(touched).tolist()
    view.set_hover([])
    assert (view._sphere_look[:, 4] == 0).all() and (view._cylinder_look[:, 8:] == 0).all()


def test_a_new_look_keeps_the_hover(qtbot) -> None:
    from chemsuite.app.protein.view import ProteinView

    structure = io.read(KINASE)
    view = ProteinView()
    qtbot.addWidget(view)
    geometry, spheres, cylinders = _scene(structure)
    view.set_scene(geometry, spheres, cylinders)
    view.set_hover([int(geometry.atoms[0])])
    _, spheres, cylinders = _scene(structure, selected=select(structure, "ligand"))
    view.set_look(spheres, cylinders)
    assert view._sphere_look[0, 4] == scene.HOVER
    assert set(np.unique(view._sphere_look[:, 4])) == {0.0, scene.SELECTED, scene.HOVER}


def test_mesh_layers_light_up_by_their_owners(qtbot) -> None:
    from chemsuite.app.protein.view import ATOM, ProteinView
    from chemsuite.protein import cartoon, mesh, surface

    structure = io.read(KINASE)
    view = ProteinView()
    qtbot.addWidget(view)
    view.set_scene(*_scene(structure))
    tube = cartoon.build(structure)
    view.set_mesh("cartoon", tube, mesh.look(tube, np.ones((structure.n_residues, 4)),
                                             np.zeros(structure.n_residues)))  # fmt: skip
    ligand = np.flatnonzero(structure.atom_kind == 2)
    skin = surface.build(structure, ligand)
    view.set_mesh("surface", skin, mesh.look(skin, np.ones((structure.n_atoms, 4)),
                                             np.zeros(structure.n_atoms)), owner=ATOM,
                  transparent=True)  # fmt: skip
    assert view.mesh_names == ["cartoon", "surface"]
    residue = int(tube.owner[100])
    view.set_hover(ligand[:3], [residue])
    cartoon_look = view._layers["cartoon"].look
    assert np.array_equal(cartoon_look[:, 4] == scene.HOVER, tube.owner == residue)
    surface_look = view._layers["surface"].look
    assert np.array_equal(surface_look[:, 4] == scene.HOVER, np.isin(skin.owner, ligand[:3]))
    # a new look keeps the hover; removing a layer forgets it
    view.set_mesh_look("cartoon", mesh.look(tube, np.zeros((structure.n_residues, 4)),
                                            np.zeros(structure.n_residues)))  # fmt: skip
    assert (view._layers["cartoon"].look[tube.owner == residue, 4] == scene.HOVER).all()
    view.set_hover([])
    assert (view._layers["cartoon"].look[:, 4] == 0).all()
    view.remove_mesh("surface")
    assert view.mesh_names == ["cartoon"]


def test_a_click_on_the_cartoon_picks_its_residue(qtbot) -> None:
    from chemsuite.app.protein.view import ProteinView
    from chemsuite.protein import cartoon

    structure = io.read(KINASE)
    view = ProteinView()
    qtbot.addWidget(view)
    view.resize(400, 300)
    hidden = np.zeros(structure.n_atoms, np.uint8)
    geometry = scene.geometry(structure, hidden)
    view.set_scene(geometry, np.zeros((0, 5)), np.zeros((0, 10)))
    _, atoms = cartoon.trace(structure)
    view.set_pick_points(atoms, structure.coords[0][atoms], np.full(len(atoms), 1.2))
    view.reset_view()
    xy, depth = view.camera.project(structure.coords[0][atoms], 400, 300)
    front = int(np.argmin(depth))
    assert view.atom_at(*xy[front]) == atoms[front]
    assert view.atom_at(1, 1) == -1


RENDER = """
import sys
import numpy as np
from PySide6.QtWidgets import QApplication
app = QApplication([])
from chemsuite.app.protein.view import ProteinView
from chemsuite.protein import colours, io, scene
from chemsuite.protein.select import select
s = io.read(sys.argv[1])
styles = np.where(select(s, "ligand or within 6 of ligand"), scene.default_styles(s), 0)
geometry = scene.geometry(s, styles.astype(np.uint8))
rgba = colours.colour(s)
def look(selected):
    flags = scene.flags(s.n_atoms, selected=selected)
    return scene.sphere_look(geometry, rgba, flags), scene.cylinder_look(geometry, rgba, flags)
view = ProteinView()
failed = []
view.failed.connect(lambda: failed.append(1))
view.resize(400, 400)
view.set_scene(geometry, *look(None), fit=True)
try:
    view.show()
    for _ in range(20):
        app.processEvents()
except Exception as error:
    print("NO-RHI", type(error).__name__)
    sys.exit(0)
before = view.grabFramebuffer()
if before.isNull():
    print("NO-RHI")
    sys.exit(0)
carbons = np.flatnonzero(select(s, "ligand and elem C"))
xy, depth = view.camera.project(s.coords[0][carbons], before.width(), before.height())
front = int(np.argmin(depth))
x, y = (int(v) for v in xy[front])
point = s.coords[0][carbons[front : front + 1]]
lx, ly = (int(v) for v in view.camera.project(point, view.width(), view.height())[0][0])
picked = view.atom_at(lx, ly)
from chemsuite.molview.pick import pick_sphere
origin, direction = view.camera.ray(lx, ly, view.width(), view.height())
hit = pick_sphere(origin, direction, geometry.spheres[:, :3], geometry.pick_radii)
view.set_look(*look(select(s, "ligand")))
for _ in range(10):
    app.processEvents()
after = view.grabFramebuffer()
from chemsuite.protein import cartoon, mesh
tube = cartoon.build(s)
blue = np.tile(np.array([0.1, 0.2, 0.95, 1.0], np.float32), (s.n_residues, 1))
nothing = scene.geometry(s, np.zeros(s.n_atoms, np.uint8))
view.set_scene(nothing, np.zeros((0, 5)), np.zeros((0, 10)))
view.set_mesh("cartoon", tube, mesh.look(tube, blue, np.zeros(s.n_residues)))
_, trace = cartoon.trace(s)
view.set_pick_points(trace, s.coords[0][trace], np.full(len(trace), 1.2))
view.reset_view()
for _ in range(10):
    app.processEvents()
drawn = view.grabFramebuffer()
cxy, cdepth = view.camera.project(s.coords[0][trace], drawn.width(), drawn.height())
c = drawn.pixelColor(*(int(v) for v in cxy[int(np.argmin(cdepth))]))
print("CARTOON", c.red(), c.green(), c.blue())
from PySide6.QtGui import QRhi
print("HOLDS", any(isinstance(v, QRhi) for v in vars(view).values()))
from PySide6.QtWidgets import QVBoxLayout, QWidget
window = QWidget()
QVBoxLayout(window).addWidget(view)  # another top-level, as when its dock floats
window.resize(400, 400)
window.show()
for _ in range(20):
    app.processEvents()
moved = view.grabFramebuffer()
mxy, mdepth = view.camera.project(s.coords[0][trace], moved.width(), moved.height())
c = moved.pixelColor(*(int(v) for v in mxy[int(np.argmin(mdepth))]))
print("MOVED", c.red(), c.green(), c.blue())
print("FAILED" if failed else "OK", view.rhi().backendName())
print("PICK", picked >= 0 and picked == int(geometry.atoms[hit]))
for name, image in (("BEFORE", before), ("AFTER", after)):
    c = image.pixelColor(x, y)
    print(name, c.red(), c.green(), c.blue())
"""


@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_proteins_are_drawn_picked_and_selected() -> None:
    """Draws a binding site with the GPU in a subprocess (the offscreen platform the other
    tests use cannot draw): the ligand's carbons are green, a click picks the atom in front,
    selecting the ligand tints it with the selection colour, and a cartoon draws."""
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER, str(KINASE)], capture_output=True, text=True, env=env,
        cwd=ROOT, timeout=120,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        if os.environ.get("CHEMSUITE_REQUIRE_GPU"):
            pytest.fail("no GPU drawing: " + out.stderr)
        pytest.skip("no GPU drawing on this machine")
    lines = {line.split()[0]: line.split()[1:] for line in out.stdout.splitlines()}
    assert "OK" in lines, out.stdout + out.stderr
    assert lines["PICK"] == ["True"]
    r, g, b = (int(v) for v in lines["BEFORE"])
    assert g > 120 and g > r + 40 and g > b + 40, (r, g, b)  # ligand carbon green
    r2, g2, b2 = (int(v) for v in lines["AFTER"])
    assert r2 > r + 40 and b2 < g2, (r2, g2, b2)  # tinted towards the selection's amber
    r3, g3, b3 = (int(v) for v in lines["CARTOON"])
    assert b3 > 60 and b3 > 3 * r3 and b3 > 2 * g3, (r3, g3, b3)  # the cartoon's blue
    # holding the QRhi's wrapper outlives the QRhi when the view changes window, and PySide
    # then hands that stale wrapper back for a new object at the same address (a QMenu)
    assert lines["HOLDS"] == ["False"]
    r4, g4, b4 = (int(v) for v in lines["MOVED"])
    assert b4 > 60 and b4 > 3 * r4 and b4 > 2 * g4, (r4, g4, b4)  # still drawn after the move


RENDER_FIGURE = """
import sys
import numpy as np
from PySide6.QtWidgets import QApplication
app = QApplication([])
from chemsuite.app.protein.view import ProteinView
from chemsuite.protein import colours, io, scene
from chemsuite.protein.select import select
s = io.read(sys.argv[1])
styles = np.where(select(s, "ligand or within 6 of ligand"), scene.default_styles(s), 0)
geometry = scene.geometry(s, styles.astype(np.uint8))
rgba = colours.colour(s)
flags = scene.flags(s.n_atoms)
view = ProteinView()
view.resize(400, 300)
view.set_scene(geometry, scene.sphere_look(geometry, rgba, flags),
               scene.cylinder_look(geometry, rgba, flags), fit=True)
try:
    view.show()
    for _ in range(20):
        app.processEvents()
except Exception as error:
    print("NO-RHI", type(error).__name__)
    sys.exit(0)
if view.grabFramebuffer().isNull():
    print("NO-RHI")
    sys.exit(0)
picture = view.picture(1200, 900, transparent=True)
print("SIZE", picture.width(), picture.height())
print("CORNER", picture.pixelColor(2, 2).alpha())
print("BACK", view.background.alpha())
carbons = np.flatnonzero(select(s, "ligand and elem C"))
view.camera.orthographic = True
view.camera_changed()
for _ in range(10):
    app.processEvents()
drawn = view.grabFramebuffer()
xy, depth = view.camera.project(s.coords[0][carbons], drawn.width(), drawn.height())
front = int(np.argmin(depth))
c = drawn.pixelColor(*(int(v) for v in xy[front]))
print("ORTHO", c.red(), c.green(), c.blue())
lx, ly = (int(v) for v in view.camera.project(s.coords[0][carbons[front:front + 1]],
                                               view.width(), view.height())[0][0])
print("PICK", view.atom_at(lx, ly) >= 0)
"""


@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_figures_draw_at_any_size_and_orthographic() -> None:
    """File > Export Image draws at the size asked for, with a see-through background when
    asked, and leaves the view as it was; an orthographic view still ray-traces its atoms
    where the camera projects them and picks them."""
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER_FIGURE, str(KINASE)], capture_output=True, text=True,
        env=env, cwd=ROOT, timeout=120,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        if os.environ.get("CHEMSUITE_REQUIRE_GPU"):
            pytest.fail("no GPU drawing: " + out.stderr)
        pytest.skip("no GPU drawing on this machine")
    lines = {line.split()[0]: line.split()[1:] for line in out.stdout.splitlines()}
    assert lines["SIZE"] == ["1200", "900"], out.stdout + out.stderr
    assert lines["CORNER"] == ["0"] and lines["BACK"] == ["255"]
    r, g, b = (int(v) for v in lines["ORTHO"])
    assert g > 120 and g > r + 40 and g > b + 40, (r, g, b)  # ligand carbon green
    assert lines["PICK"] == ["True"]


RENDER_SHORT = """
import sys
import numpy as np
from PySide6.QtWidgets import QApplication
app = QApplication([])
from chemsuite.app.protein import view as views
from chemsuite.protein import cartoon, io, mesh, scene
made = []
chunks = mesh.chunks
def small(m):  # many chunks, as a large mesh would make
    parts = chunks(m, 2000)
    made.append(len(parts))
    return parts
mesh.chunks = small
create = views.ProteinView.create_pipelines
def without_wide_indices(self):
    ok = create(self)
    self._wide_indices = False  # as on OpenGL ES 2 without the extension
    return ok
views.ProteinView.create_pipelines = without_wide_indices
s = io.read(sys.argv[1])
view = views.ProteinView()
view.resize(400, 400)
nothing = scene.geometry(s, np.zeros(s.n_atoms, np.uint8))
view.set_scene(nothing, np.zeros((0, 5)), np.zeros((0, 10)))
tube = cartoon.build(s)
blue = np.tile(np.array([0.1, 0.2, 0.95, 1.0], np.float32), (s.n_residues, 1))
view.set_mesh("cartoon", tube, mesh.look(tube, blue, np.zeros(s.n_residues)))
_, trace = cartoon.trace(s)
view.set_pick_points(trace, s.coords[0][trace], np.full(len(trace), 1.2))
view.reset_view()
try:
    view.show()
    for _ in range(20):
        app.processEvents()
except Exception as error:
    print("NO-RHI", type(error).__name__)
    sys.exit(0)
drawn = view.grabFramebuffer()
if drawn.isNull():
    print("NO-RHI")
    sys.exit(0)
xy, depth = view.camera.project(s.coords[0][trace], drawn.width(), drawn.height())
c = drawn.pixelColor(*(int(v) for v in xy[int(np.argmin(depth))]))
print("CHUNKS", made[0] if made else 0)
print("CARTOON", c.red(), c.green(), c.blue())
"""


@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_meshes_draw_in_chunks_without_32_bit_indices() -> None:
    """Where the GPU has only 16-bit indices, a cartoon is drawn in chunks, not skipped."""
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER_SHORT, str(KINASE)], capture_output=True, text=True,
        env=env, cwd=ROOT, timeout=120,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        if os.environ.get("CHEMSUITE_REQUIRE_GPU"):
            pytest.fail("no GPU drawing: " + out.stderr)
        pytest.skip("no GPU drawing on this machine")
    lines = {line.split()[0]: line.split()[1:] for line in out.stdout.splitlines()}
    assert int(lines["CHUNKS"][0]) > 5, out.stdout
    r, g, b = (int(v) for v in lines["CARTOON"])
    assert b > 60 and b > 3 * r and b > 2 * g, (r, g, b)  # the cartoon's blue


RENDER_FRAMES = """
import sys, time
from dataclasses import replace
import numpy as np
from PySide6.QtWidgets import QApplication
app = QApplication([])
from chemsuite.app.protein.view import ProteinView
from chemsuite.protein import display, io, poses
base = io.read(sys.argv[1])
big, k = base, 1
while big.n_atoms < 100_000:  # copies side by side, each with its own chains
    shift = np.array([(k % 8) * 60.0, (k // 8) * 60.0, 0.0], np.float32)
    big = poses.merged(big, replace(base.with_coords(base.coords + shift),
        chain_name=np.array([f"{c}{k}" for c in base.chain_name], dtype="<U4")))
    k += 1
scene = display.compose([display.Display("big", "big", big)])
view = ProteinView()
view.resize(800, 600)
view.set_scene(scene.geometry, scene.sphere_look, scene.cylinder_look, fit=True)
view.set_mesh("cartoon", scene.cartoon, scene.cartoon_look)
submitted = []
view.frameSubmitted.connect(lambda: submitted.append(time.perf_counter()))
try:
    view.show()
    for _ in range(20):
        app.processEvents()
except Exception as error:
    print("NO-RHI", type(error).__name__)
    sys.exit(0)
if view.grabFramebuffer().isNull():
    print("NO-RHI")
    sys.exit(0)
times = []
for _ in range(60):
    count, start = len(submitted), time.perf_counter()
    view.camera.turn(6.0, 2.0)
    view.update()
    while len(submitted) == count and time.perf_counter() - start < 2.0:
        app.processEvents()
    times.append((time.perf_counter() - start) * 1000)
print("ATOMS", big.n_atoms)
print("BACKEND", view.rhi().backendName())
print("MEDIAN", f"{np.median(times):.1f}")
print("P95", f"{np.percentile(times, 95):.1f}")
"""
#: ms per frame while the camera turns around 100,000 atoms, asserted on Metal (a Mac's GPU);
#: other backends (CI's software Mesa) only report it
FRAME_MS = 33.0


@pytest.mark.bench
@pytest.mark.skipif(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    reason="needs a display for the GPU",
)
def test_turning_a_large_protein_keeps_its_frame_rate() -> None:
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_PLATFORM"}
    out = subprocess.run(
        [sys.executable, "-c", RENDER_FRAMES, str(ROOT / "tests/data/protein/3htb.cif.gz")],
        capture_output=True, text=True, env=env, cwd=ROOT, timeout=300,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    if "NO-RHI" in out.stdout:
        pytest.skip("no GPU drawing on this machine")
    lines = {line.split()[0]: line.split()[1:] for line in out.stdout.splitlines()}
    median, p95 = float(lines["MEDIAN"][0]), float(lines["P95"][0])
    backend = lines["BACKEND"][0]
    print(f"\n{lines['ATOMS'][0]} atoms on {backend}: median {median} ms, 95% {p95} ms a frame")
    if backend == "Metal":
        assert median < FRAME_MS, out.stdout
