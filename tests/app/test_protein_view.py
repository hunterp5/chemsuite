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
    and selecting the ligand tints it with the selection colour."""
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
