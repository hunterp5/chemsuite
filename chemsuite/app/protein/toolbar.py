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
"""The Protein Viewer's toolbar, and the Focus and Measure menus it shares with the menu bar.

Show, Hide and Colour act on the selection; Hydrogens and Transparency on every protein;
then Interactions, Focus and Measure for pocket work, and Frame and Reset View for the camera.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtGui import QActionGroup
from PySide6.QtWidgets import QMenu, QToolBar, QToolButton

from . import menus
from .session import ProteinSession

if TYPE_CHECKING:
    from .panel import ProteinPanel

MEASURES = (("Distance", 2), ("Angle", 3), ("Dihedral", 4))


def build(panel: ProteinPanel) -> QToolBar:
    """The toolbar; sets ``panel.interactions_action`` and ``panel.frame_action``."""
    session = panel.session
    bar = QToolBar()
    bar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
    menu_button(bar, "Show", lambda m: menus.fill_show(m, session))
    menu_button(bar, "Hide", lambda m: menus.fill_hide(m, session))
    menu_button(bar, "Colour", lambda m: menus.fill_colour(m, session, panel))
    menu_button(bar, "Hydrogens", lambda m: menus.fill_hydrogens(m, session),
                "Which hydrogens to draw, on every protein")  # fmt: skip
    menu_button(bar, "Transparency", lambda m: menus.fill_transparency(m, session),
                "How see-through surfaces are")  # fmt: skip
    bar.addSeparator()
    interactions = bar.addAction("Interactions")
    interactions.setCheckable(True)
    interactions.setToolTip(
        "Show the interactions of the selection, or of every ligand when nothing is selected"
    )
    interactions.triggered.connect(session.show_interactions)
    panel.interactions_action = interactions
    menu_button(bar, "Focus", lambda m: fill_focus(m, session))
    menu_button(bar, "Measure", lambda m: fill_measure(m, session))
    bar.addSeparator()
    panel.frame_action = bar.addAction("Frame", panel.frame_selection)
    panel.frame_action.setToolTip("Fit the selection to the view (double-click an entry)")
    bar.addAction("Reset View", panel.reset_view)
    return bar


def menu_button(
    bar: QToolBar, title: str, fill: Callable[[QMenu], None], tip: str = ""
) -> QToolButton:
    button = QToolButton()
    button.setText(title)
    button.setToolTip(tip)
    button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
    menu = QMenu(button)
    menu.aboutToShow.connect(lambda: fill(menu))
    button.setMenu(menu)
    bar.addWidget(button)
    return button


def fill_focus(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    several = len([d for d in session.displays if d.visible]) > 1
    for protein, residue in session.ligands():
        shown = session.display(protein)
        label = shown.structure.residue_label(residue)
        title = f"{shown.name}  {label}" if several else label
        action = menu.addAction(title, lambda p=protein, r=residue: session.focus(p, r))
        action.setCheckable(True)
        action.setChecked(session.focused_residue() == (protein, residue))
    if not menu.actions():
        menu.addAction("No Ligands").setEnabled(False)
    menu.addSeparator()
    clear = menu.addAction("Clear Focus", session.clear_focus)
    clear.setEnabled(session.focused is not None)


def fill_measure(menu: QMenu, session: ProteinSession) -> None:
    menu.clear()
    group = QActionGroup(menu)
    for title, atoms in MEASURES:
        action = menu.addAction(title, lambda n=atoms: session.measure(n))
        action.setCheckable(True)
        action.setChecked(session.measuring == atoms)
        group.addAction(action)
    menu.addSeparator()
    stop = menu.addAction("Stop Measuring", lambda: session.measure(0))
    stop.setEnabled(bool(session.measuring))
    clear = menu.addAction("Clear Measurements", session.clear_measurements)
    clear.setEnabled(bool(session.measurements))
