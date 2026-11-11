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
"""The Protein Viewer panel: layout and wiring only.

A toolbar (Show, Hide, Colour, Frame, Reset View) and the selection bar run across the top;
below, the proteins list beside the 3D view, and the sequence strip under them. Everything
reads and changes one :class:`~.session.ProteinSession`.

In the view a click selects a residue, Alt-click an atom, a double-click its chain; Shift adds
and Ctrl (Cmd) toggles; a click on nothing clears the selection; Esc clears it too. When the
computer cannot draw 3D, the view's place says why and the rest keeps working.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QLabel,
    QMenu,
    QSplitter,
    QStackedWidget,
    QToolBar,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app import graphics
from chemsuite.core.table import Table
from chemsuite.molview.camera import Camera

from . import menus
from .selection_bar import SelectionBar
from .sequence_view import SequenceView
from .session import ADD, SET, TOGGLE, ProteinSession
from .tree import ProteinTree

TOOL_ID = "protein"
CARTOON = "cartoon"
SURFACE = "surface"


class ProteinPanel(QWidget):
    def __init__(self, table: Table, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.tool_id = TOOL_ID
        self.session = ProteinSession(table, self)
        self.view = None
        self._fitted: set[str] = set()  # proteins the camera has been fitted to
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self.toolbar = QToolBar()
        self.toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        for title, fill in (
            ("Show", lambda m: menus.fill_show(m, self.session)),
            ("Hide", lambda m: menus.fill_hide(m, self.session)),
            ("Colour", lambda m: menus.fill_colour(m, self.session, self)),
        ):
            self._menu_button(title, fill)
        self.toolbar.addSeparator()
        self.frame_action = self.toolbar.addAction("Frame", self.frame_selection)
        self.frame_action.setToolTip("Fit the selection to the view (double-click an entry)")
        self.toolbar.addAction("Reset View", self.reset_view)
        self.selection_bar = SelectionBar(self.session)

        self.tree = ProteinTree(self.session)
        self.tree.focusRequested.connect(self.frame_selection)
        self.sequence = SequenceView(self.session)
        self.sequence.strip.residueHovered.connect(self._sequence_hovered)
        self.views = QStackedWidget()
        self.status = QLabel()
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self.columns = QSplitter(Qt.Orientation.Horizontal)
        self.columns.addWidget(self.tree)
        self.columns.addWidget(self.views)
        self.columns.setStretchFactor(1, 1)
        self.columns.setSizes([200, 700])
        self.rows = QSplitter(Qt.Orientation.Vertical)
        self.rows.addWidget(self.columns)
        self.rows.addWidget(self.sequence)
        self.rows.setStretchFactor(0, 1)
        self.rows.setSizes([700, self.sequence.sizeHint().height()])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.selection_bar)
        layout.addWidget(self.rows, 1)
        layout.addWidget(self.status)

        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, self.session.clear_selection)
        self.session.sceneChanged.connect(self._redraw)
        self.session.lookChanged.connect(self._relook)
        self.session.proteinsChanged.connect(self._proteins_changed)
        self.session.hovered.connect(self._hovered)
        self._make_view()
        self._redraw()

    def _menu_button(self, title: str, fill) -> None:
        button = QToolButton()
        button.setText(title)
        button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(button)
        menu.aboutToShow.connect(lambda: fill(menu))
        button.setMenu(menu)
        self.toolbar.addWidget(button)

    def title(self) -> str:
        return "Protein Viewer"

    # the 3D view, or why there is none
    def _make_view(self) -> None:
        if reason := graphics.problem():
            self._no_view(reason)
            return
        from .view import ProteinView

        self.view = ProteinView()
        self.view.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.view.set_background(self.palette().base().color())
        self.view.failed.connect(lambda: self._no_view(graphics.NO_GPU))
        self.view.atomHovered.connect(self.session.hover)
        self.view.atomClicked.connect(self._clicked)
        self.view.atomDoubleClicked.connect(self._double_clicked)
        self.view.contextRequested.connect(self._context_menu)
        self.views.addWidget(self.view)

    def _no_view(self, reason: str) -> None:
        graphics.set_problem(reason)
        label = QLabel(f"⚠ {reason}")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setWordWrap(True)
        self.views.addWidget(label)
        self.views.setCurrentWidget(label)
        if self.view is not None:
            self.view.hide()
            self.view = None

    # redrawing
    def _redraw(self) -> None:
        if self.view is None:
            return
        from .view import ATOM

        scene = self.session.scene
        view = self.view
        view.set_pick_points(
            scene.pick_atoms, scene.pick_centres, np.full(len(scene.pick_atoms), 1.2)
        )
        view.set_mesh(CARTOON, scene.cartoon, scene.cartoon_look)
        view.set_mesh(SURFACE, scene.surface, scene.surface_look, owner=ATOM, transparent=True)
        new = {d.id for d in scene.displays} - self._fitted
        view.set_scene(scene.geometry, scene.sphere_look, scene.cylinder_look)
        if new:
            if not self._fitted:
                view.orient()
            view.fit()
            self._fitted |= new

    def _relook(self) -> None:
        if self.view is None:
            return
        scene = self.session.scene
        self.view.set_look(scene.sphere_look, scene.cylinder_look)
        self.view.set_mesh_look(CARTOON, scene.cartoon_look)
        self.view.set_mesh_look(SURFACE, scene.surface_look)

    def _proteins_changed(self) -> None:
        self._fitted &= {d.id for d in self.session.displays if d.visible}

    # the camera
    def frame_selection(self) -> None:
        """Fit the selection (everything, when nothing is selected) to the view."""
        if self.view is None:
            return
        points = [d.structure.xyz(0)[d.selected] for d in self.session.displays
                  if d.visible and d.selected.any()]  # fmt: skip
        if not points:
            self.view.fit()
            return
        positions = np.concatenate(points)
        self.view.frame(positions, np.full(len(positions), 1.5))

    def reset_view(self) -> None:
        if self.view is not None:
            self.view.reset_view()

    # picking in the view
    def _residue_mask(self, atom: int, whole_chain: bool = False, single: bool = False):
        found = self.session.scene.locate(atom)
        if found is None:
            return None
        shown, index = found
        s = shown.structure
        if single:
            mask = np.zeros(s.n_atoms, bool)
            mask[index] = True
        elif whole_chain:
            mask = s.atom_chain == s.atom_chain[index]
        else:
            mask = s.residue_of == s.residue_of[index]
        return {shown.id: mask & shown.active}

    def _clicked(self, atom: int, modifiers) -> None:
        mode = SET
        if modifiers & Qt.KeyboardModifier.ShiftModifier:
            mode = ADD
        elif modifiers & Qt.KeyboardModifier.ControlModifier:
            mode = TOGGLE
        masks = self._residue_mask(atom, single=bool(modifiers & Qt.KeyboardModifier.AltModifier))
        if masks is None:
            if mode == SET:
                self.session.clear_selection()
            return
        self.session.select(masks, mode)

    def _double_clicked(self, atom: int) -> None:
        masks = self._residue_mask(atom, whole_chain=True)
        if masks is not None:
            self.session.select(masks, SET)

    def _context_menu(self, atom: int, where) -> None:
        menu = menus.context_menu(self.session, atom, self, focus=self.frame_selection)
        menu.exec(where)

    # hover: the view and the sequence light up the residue; the status line names the atom
    def _hovered(self, protein: str, atom: int) -> None:
        shown = self.session.display(protein)
        if shown is None:
            self.status.setText("")
            self._light(None, -1)
            return
        s = shown.structure
        b = f"  B {s.b[atom]:.1f}" if s.b[atom] else ""
        self.status.setText(f"{shown.name}  {s.atom_label(atom)}{b}")
        self._light(shown, int(s.residue_of[atom]))

    def _sequence_hovered(self, protein: str, residue: int) -> None:
        shown = self.session.display(protein)
        self._light(shown, residue)
        if shown is not None and residue >= 0:
            self.status.setText(f"{shown.name}  {shown.structure.residue_label(residue)}")
        else:
            self.status.setText("")

    def _light(self, shown, residue: int) -> None:
        if self.view is None:
            return
        if shown is None or residue < 0 or shown not in self.session.scene.displays:
            self.view.set_hover([])
            return
        span = shown.structure.atoms_of(residue)
        scene = self.session.scene
        atoms = scene.global_atoms(shown, np.arange(span.start, span.stop))
        self.view.set_hover(atoms, scene.global_residues(shown, np.array([residue])))

    # saving
    def state(self) -> dict:
        state = {
            "session": self.session.state(),
            "columns": self.columns.sizes(),
            "rows": self.rows.sizes(),
        }
        if self.view is not None:
            state["camera"] = self.view.camera.state()
        return state

    def set_state(self, state: dict) -> None:
        self.session.set_state(state.get("session", {}))
        for splitter, key in ((self.columns, "columns"), (self.rows, "rows")):
            sizes = state.get(key)
            if isinstance(sizes, list) and len(sizes) == splitter.count():
                splitter.setSizes([int(v) for v in sizes])
        camera = state.get("camera")
        if isinstance(camera, dict) and self.view is not None:
            self.view.set_camera(Camera.from_state(camera))
            self._fitted = {d.id for d in self.session.displays if d.visible}

    def closeEvent(self, event) -> None:
        self.session.close()
        super().closeEvent(event)
