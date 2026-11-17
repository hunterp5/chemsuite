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
"""Display ▸ Surface and View ▸ Ligand ▸ Pocket Surface: how surfaces are made and drawn.

Display ▸ Surface changes the surfaces of the proteins with a selection, else of every visible
one, and its checks show the first of those (``WorkbenchSession.surface_displays``). Carving
cuts every protein's surface around the selected atoms, whichever protein they belong to (a
ligand or a pose can be another one).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import replace

from PySide6.QtGui import QActionGroup, QColor
from PySide6.QtWidgets import QColorDialog, QMenu, QWidget

from chemsuite.macromolecule.look import colours, surface
from chemsuite.macromolecule.look import display as displays
from chemsuite.macromolecule.workspace.session import WorkbenchSession

KIND_TIPS = {
    surface.SES: "The surface a water-sized probe rolling over the atoms cannot pass "
    "(Connolly's molecular surface)",
    surface.SAS: "Where the probe's center can go: the atoms grown by the probe's radius",
    surface.VDW: "The atoms' van der Waals spheres",
    surface.BLOB: "A smooth Gaussian skin over the atoms (fast, and how older sessions drew)",
}
QUALITY_TITLES = {"coarse": "Coarse", "normal": "Normal", "fine": "Fine"}
SMOOTHING_TITLES = {"off": "Off", "low": "Low", "high": "High"}


def fill_surface(menu: QMenu, session: WorkbenchSession, parent: QWidget) -> None:
    menu.clear()
    menu.setToolTipsVisible(True)
    targets = session.surface_displays()
    menu.setEnabled(bool(targets))
    if not targets:
        return
    first = targets[0]
    options = first.surface_options

    def options_to(**changes: object) -> None:
        session.change_surfaces(
            lambda d: d.set_surface_options(replace(d.surface_options, **changes))
        )

    kind = _choices(menu.addMenu("Kind"), surface.KINDS.items(), options.kind,
                    lambda k: options_to(kind=k))  # fmt: skip
    for action in kind.actions():
        action.setToolTip(KIND_TIPS[action.data()])
    _colour(menu.addMenu("Color"), session, parent, first)
    _choices(
        menu.addMenu("Transparency"),
        _transparencies(),
        first.transparency,
        lambda v: session.change_surfaces(lambda d: _set(d, "transparency", v), rebuilt=False),
    )
    menu.addSeparator()
    carve = menu.addMenu("Carve Around Selection")
    carve.setToolTipsVisible(True)
    carve.setEnabled(session.selected_count() > 0)
    for distance in displays.CARVE_DISTANCES:
        carve.addAction(f"{distance:g} Å", lambda x=distance: _carve(session, x)).setToolTip(
            f"Show only the surface within {distance:g} Å of the selected atoms"
        )
    whole = menu.addAction(
        "Show Whole Surface",
        lambda: session.change_surfaces(lambda d: d.carve_surface(None), every=True),
    )
    whole.setEnabled(any(d.carve_points is not None for d in session.displays))
    menu.addSeparator()
    probe = _choices(menu.addMenu("Probe Radius"), ((p, f"{p:g} Å") for p in surface.PROBES),
                     options.probe, lambda p: options_to(probe=p))  # fmt: skip
    probe.setEnabled(options.kind in (surface.SES, surface.SAS))
    probe.setToolTip("The solvent probe's radius (water: 1.4 Å)")
    _choices(menu.addMenu("Quality"), QUALITY_TITLES.items(), options.quality,
             lambda q: options_to(quality=q))  # fmt: skip
    smoothing = _choices(menu.addMenu("Smoothing"), SMOOTHING_TITLES.items(),
                         options.smoothing, lambda s: options_to(smoothing=s))  # fmt: skip
    smoothing.setEnabled(options.kind != surface.BLOB)


def fill_pocket(menu: QMenu, session: WorkbenchSession) -> None:
    """View ▸ Ligand ▸ Pocket Surface: the skin Focus draws over a ligand's pocket."""
    menu.clear()
    shown = menu.addAction("Show", lambda: session.show_pocket_surface(not session.pocket_surface))
    shown.setCheckable(True)
    shown.setChecked(session.pocket_surface)
    first = session.displays[0] if session.displays else None
    current = first.pocket_transparency if first else displays.SURFACE_TRANSPARENCY
    _choices(
        menu.addMenu("Transparency"),
        _transparencies(),
        current,
        lambda v: session.change_surfaces(
            lambda d: _set(d, "pocket_transparency", v), every=True, rebuilt=False
        ),
    )


def _colour(
    menu: QMenu, session: WorkbenchSession, parent: QWidget, first: displays.Display
) -> None:
    group = QActionGroup(menu)
    for key, title in colours.SCHEMES.items():
        action = menu.addAction(
            f"By {title}",
            lambda k=key: session.change_surfaces(lambda d: d.colour_surface(k), rebuilt=False),
        )
        action.setCheckable(True)
        action.setChecked(first.surface_scheme == key)
        group.addAction(action)
    menu.addSeparator()
    single = menu.addAction("Single Color…", lambda: _pick(session, parent, first.surface_rgb))
    single.setCheckable(True)
    single.setChecked(first.surface_scheme == displays.FLAT)
    group.addAction(single)


def _pick(session: WorkbenchSession, parent: QWidget, rgb: tuple[float, float, float]) -> None:
    chosen = QColorDialog.getColor(QColor.fromRgbF(*rgb), parent, "Surface Color")
    if chosen.isValid():
        picked = (chosen.redF(), chosen.greenF(), chosen.blueF())
        session.change_surfaces(lambda d: d.colour_surface(displays.FLAT, picked), rebuilt=False)


def _carve(session: WorkbenchSession, distance: float) -> None:
    points = session.selected_points()
    session.change_surfaces(lambda d: d.carve_surface(points, distance), every=True)


def _set(shown: displays.Display, name: str, value: float) -> None:
    if name == "transparency":
        shown.transparency = value
    else:
        shown.pocket_transparency = value


def _transparencies() -> list[tuple[float, str]]:
    return [(v, f"{v:.0%}" if v else "Opaque") for v in displays.TRANSPARENCIES]


def _choices(
    menu: QMenu,
    items: Iterable[tuple[object, str]],
    current: object,
    choose: Callable[[object], None],
) -> QMenu:
    """One checkable action per ``(value, title)``, the current value checked; each action's
    data is its value."""
    menu.setToolTipsVisible(True)
    group = QActionGroup(menu)
    for value, title in items:
        action = menu.addAction(title, lambda v=value: choose(v))
        action.setData(value)
        action.setCheckable(True)
        if isinstance(value, float) and isinstance(current, float):
            same = abs(value - current) < 1e-6
        else:
            same = value == current
        action.setChecked(bool(same))
        group.addAction(action)
    return menu
