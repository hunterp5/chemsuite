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

The viewer is a workspace of its own (it opens as its own window): a menu bar
(:mod:`.menubar`: opening and fetching proteins, saving, copying, display, pocket work), a
toolbar (:mod:`.toolbar`: Show, Hide, Colour, Hydrogens, Transparency; Interactions, Focus,
Measure; Frame, Reset View) and the
selection bar run across the top; below, the proteins list beside the 3D view (with the
models bar under it while a shown protein has several models), and the sequence strip under
them. Everything reads and changes one :class:`~.session.ProteinSession`.

In the view a click selects a residue, Alt-click an atom, a double-click its chain; Shift adds
and Ctrl (Cmd) toggles; a click on nothing clears the selection; Esc clears it too. While
measuring, clicks pick atoms instead and Esc stops. A focused ligand keeps a slab (a slice of
depth around its pocket) as the camera moves. When the computer cannot draw 3D, the view's
place says why and the rest keeps working.

The panel follows the table: when the current row names a protein (``Protein``) and one of its
ligands (``Ligand ID``), that ligand is focused, so stepping through ligand rows walks the
pockets.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app import graphics
from chemsuite.core.table import Table
from chemsuite.molview.camera import Camera
from chemsuite.protein import io as protein_io
from chemsuite.protein import ligands

from . import menubar, menus, toolbar
from .builds import Run
from .models_bar import ModelsBar
from .pocket import PocketCamera, follow_row
from .selection_bar import SelectionBar
from .sequence_view import SequenceView
from .session import ADD, SET, TOGGLE, ProteinSession
from .tree import ProteinTree

TOOL_ID = "protein"
CARTOON = "cartoon"
SURFACE = "surface"
LINES = "interactions"
POCKET = "pocket"
#: the commands a tool gives the viewer (see :mod:`.menubar`)
OPEN, FETCH, EXTRACT, LIGAND = menubar.OPEN, menubar.FETCH, menubar.EXTRACT, menubar.LIGAND
PROFILE = menubar.PROFILE


class ProteinPanel(QWidget):
    def __init__(
        self,
        table: Table,
        parent: QWidget | None = None,
        *,
        current: Callable[[], int | None] | None = None,
        watch_current: Callable[[Callable[[int | None], None]], Callable[[], None]] | None = None,
        commands: menubar.Commands | None = None,
        run: Run | None = None,
    ) -> None:
        super().__init__(parent)
        self.tool_id = TOOL_ID
        self.session = ProteinSession(table, self, run=run)
        self.view = None
        self.overlay = None
        self._view_page: QWidget | None = None
        self._problem_page: QLabel | None = None
        commands = dict(commands or {})
        self.commands = commands
        self._fitted: set[str] = set()  # proteins the camera has been fitted to
        self.pocket_camera = PocketCamera(self)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self.toolbar = toolbar.build(self)
        self.selection_bar = SelectionBar(self.session)

        self.tree = ProteinTree(self.session)
        self.tree.focusRequested.connect(self.frame_selection)
        self.sequence = SequenceView(self.session)
        self.sequence.strip.residueHovered.connect(self._sequence_hovered)
        self.views = QStackedWidget()
        self.empty = self._empty_page(commands)
        self.views.addWidget(self.empty)
        self.models_bar = ModelsBar(self.session)
        middle = QWidget()
        stack = QVBoxLayout(middle)
        stack.setContentsMargins(0, 0, 0, 0)
        stack.setSpacing(0)
        stack.addWidget(self.views, 1)
        stack.addWidget(self.models_bar)
        self.status = QLabel()
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self.columns = QSplitter(Qt.Orientation.Horizontal)
        self.columns.addWidget(self.tree)
        self.columns.addWidget(middle)
        self.columns.setStretchFactor(1, 1)
        self.columns.setSizes([200, 700])
        self.rows = QSplitter(Qt.Orientation.Vertical)
        self.rows.addWidget(self.columns)
        self.rows.addWidget(self.sequence)
        self.rows.setStretchFactor(0, 1)
        self.rows.setSizes([700, self.sequence.sizeHint().height()])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 0, 6, 6)
        self.menu_bar = menubar.build(self, commands)
        layout.setMenuBar(self.menu_bar)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.selection_bar)
        layout.addWidget(self.rows, 1)
        layout.addWidget(self.status)

        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, self._escape)
        QShortcut(QKeySequence("["), self, lambda: self.models_bar.step(-1))
        QShortcut(QKeySequence("]"), self, lambda: self.models_bar.step(1))
        self.session.sceneChanged.connect(self._redraw)
        self.session.lookChanged.connect(self._relook)
        self.session.proteinsChanged.connect(self._proteins_changed)
        self.session.hovered.connect(self._hovered)
        self._make_view()
        self._redraw()
        self._show_page()
        self._stop_watching = watch_current(self.follow_row) if watch_current else lambda: None
        if current is not None:
            self.follow_row(current())
            self.session.rowsEdited.connect(lambda oids: self._row_edited(current(), oids))

    def _escape(self) -> None:
        if self.session.measuring:
            self.session.measure(0)
        else:
            self.session.clear_selection()

    def title(self) -> str:
        return "Protein Viewer"

    def close_viewer(self) -> None:
        """Close the viewer's window (its dock)."""
        parent = self.parentWidget()
        while parent is not None and not isinstance(parent, QDockWidget):
            parent = parent.parentWidget()
        (parent or self).close()

    # what the middle shows: the view, why there is none, or ways to add a protein
    def _empty_page(self, commands: menubar.Commands) -> QWidget:
        page = QWidget()
        row = QHBoxLayout(page)
        row.addStretch(1)
        for key, title in ((menubar.OPEN, "Open Protein…"), (menubar.FETCH, "Fetch from PDB…")):
            if key in commands:
                button = QPushButton(title)
                button.clicked.connect(lambda _c=False, k=key: commands[k]())
                row.addWidget(button)
        row.addStretch(1)
        return page

    def _show_page(self) -> None:
        if not self.session.displays:
            self.views.setCurrentWidget(self.empty)
        elif self._problem_page is not None:
            self.views.setCurrentWidget(self._problem_page)
        elif self._view_page is not None:
            self.views.setCurrentWidget(self._view_page)

    # the 3D view, or why there is none
    def _make_view(self) -> None:
        if reason := graphics.problem():
            self._no_view(reason)
            return
        from .overlay import ProteinOverlay
        from .view import ProteinView

        self.view = ProteinView()
        self.view.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.view.set_background(self.palette().base().color())
        self.view.failed.connect(lambda: self._no_view(graphics.NO_GPU))
        self.view.atomHovered.connect(self.session.hover)
        self.view.atomClicked.connect(self._clicked)
        self.view.atomDoubleClicked.connect(self._double_clicked)
        self.view.contextRequested.connect(self._context_menu)
        self.pocket_camera.watch(self.view)
        self.overlay = ProteinOverlay(self.view, self.session)
        page = QWidget()
        grid = QGridLayout(page)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.addWidget(self.view, 0, 0)
        grid.addWidget(self.overlay, 0, 0)
        self.views.addWidget(page)
        self._view_page = page

    def _no_view(self, reason: str) -> None:
        graphics.set_problem(reason)
        label = QLabel(f"⚠ {reason}")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setWordWrap(True)
        self.views.addWidget(label)
        self._problem_page = label
        if self.view is not None:
            self.view.hide()
            self.overlay.hide()
            self.view = None
        self._show_page()

    # redrawing
    def _redraw(self) -> None:
        self.interactions_action.setChecked(self.session.interactions_shown)
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
        view.set_mesh(POCKET, scene.pocket, scene.pocket_look, owner=ATOM, transparent=True)
        view.set_mesh(LINES, scene.lines, scene.lines_look, owner=ATOM)
        new = {d.id for d in scene.displays if d is not self.session.pose} - self._fitted
        view.set_scene(scene.geometry, scene.sphere_look, scene.cylinder_look)
        if new:
            if not self._fitted:
                view.orient()
            view.fit()
            self._fitted |= new
        self.pocket_camera.apply()  # a focused pocket moves with the model shown

    def _relook(self) -> None:
        if self.view is None:
            return
        scene = self.session.scene
        self.view.set_look(scene.sphere_look, scene.cylinder_look)
        self.view.set_mesh_look(CARTOON, scene.cartoon_look)
        self.view.set_mesh_look(SURFACE, scene.surface_look)
        self.view.set_mesh_look(POCKET, scene.pocket_look)

    def _proteins_changed(self) -> None:
        self._fitted &= {d.id for d in self.session.displays if d.visible}
        self._show_page()

    # the camera
    def frame_selection(self) -> None:
        """Fit the selection (everything, when nothing is selected) to the view."""
        if self.view is None:
            return
        points = [d.xyz[d.selected] for d in self.session.displays
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
        if found is None or found[0] is self.session.pose:  # a pose is not selected
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
        if self.session.pick(atom):  # measuring
            return
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
        menu = menus.context_menu(
            self.session,
            atom,
            self,
            focus=self.frame_selection,
            copy_smiles=self.copy_smiles,
            copy_as=self.copy_selection,
            add_ligand=self.commands.get(LIGAND),
        )
        menu.exec(where)

    # the table's current row: its ligand, in its pocket
    def follow_row(self, oid: int | None) -> bool:
        """Focus the pose or ligand that row *oid* names, if its protein is shown here."""
        return follow_row(self.session, oid)

    def _row_edited(self, oid: int | None, oids: np.ndarray) -> None:
        if oid is not None and oid in oids:
            self.follow_row(oid)

    # copying
    def copy_selection(self, fmt: str) -> str:
        """Copy the selected atoms as PDB or mmCIF text (each protein with a selection)."""
        texts = [protein_io.to_text(d.structure, fmt, keep=d.selected)
                 for d in self.session.displays if d.visible and d.selected.any()]  # fmt: skip
        text = "".join(texts)
        if text:
            QGuiApplication.clipboard().setText(text)
        return text

    def copy_smiles(self, shown, residue: int) -> str:
        from rdkit import Chem

        smiles = Chem.MolToSmiles(ligands.to_mol(shown.structure, residue).mol)
        QGuiApplication.clipboard().setText(smiles)
        return smiles

    # hover: the view and the sequence light up the residue; the status line names the atom
    def _hovered(self, protein: str, atom: int) -> None:
        shown = self.session.shown(protein)
        if shown is None:
            self.status.setText("")
            self._light(None, -1)
            return
        s = shown.structure
        b = f"  B {s.b[atom]:.1f}" if s.b[atom] else ""
        contacts = self.session.interactions_of(protein, atom)
        more = "".join(f"  ·  {c}" for c in contacts[:2])
        self.status.setText(f"{shown.name}  {s.atom_label(atom)}{b}{more}")
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
        self._stop_watching()
        self.session.close()
        super().closeEvent(event)
