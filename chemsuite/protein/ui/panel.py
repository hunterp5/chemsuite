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

The viewer is a workspace of its own, in a window of its own (``ToolContext.add_window``) whose
menu bar is :mod:`.menubar` (files, copying, selecting, display, ligand work, the camera). The
selection field runs across the top; below, the proteins list beside the 3D view (with the
models bar under it while a shown protein has several models). The sequence has a window of
its own (the Sequence menu). Everything reads and changes one
:class:`~.session.WorkbenchSession`.

In the view a click selects an atom (on a cartoon, where no atom is drawn, its residue), a
double-click the residue, Alt-click (Option) the chain; Shift adds and Ctrl (Cmd) toggles; a
click on nothing clears the selection; Esc clears it too. While measuring, clicks pick atoms
instead and Esc stops. The selection can be deleted, its atoms given another element or
charge, and its residues or chain renamed or renumbered (Edit menu), each one undoable step
of the document. A focused ligand keeps a slab (a slice of
depth around its pocket) as the camera moves. When the computer cannot draw 3D, the view's
place says why and the rest keeps working.

The panel follows the table: when the current row names a protein (``Protein``) and one of its
ligands (``Ligand ID``), that ligand is focused, so stepping through ligand rows walks the
pockets.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

import numpy as np
from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QGuiApplication, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.shell import graphics
from chemsuite.app.widgets.panel import ToolPanel
from chemsuite.core.assets import AssetHost
from chemsuite.molview.camera import Camera
from chemsuite.protein.chem import ligands
from chemsuite.protein.look import scene as scenes
from chemsuite.protein.structure import io as protein_io
from chemsuite.protein.structure.edit import EditError
from chemsuite.protein.workspace.session import ADD, SET, TOGGLE, WorkbenchSession

from . import menubar, menus
from .builds import Builder, Run
from .host import LinkedRow
from .models_bar import ModelsBar
from .pocket import PocketCamera
from .selection_bar import SelectionBar
from .sequence_view import SequenceView, SequenceWindow
from .session import SessionSignals
from .tree import ProteinTree

TOOL_ID = "protein"
CARTOON = "cartoon"
SURFACE = "surface"
LINES = "interactions"
POCKET = "pocket"
#: the commands a tool gives the viewer (see :mod:`.menubar`)
OPEN, FETCH, EXTRACT, LIGAND = menubar.OPEN, menubar.FETCH, menubar.EXTRACT, menubar.LIGAND
PROFILE, REMOVE, SHOW_TABLE = menubar.PROFILE, menubar.REMOVE, menubar.SHOW_TABLE
THEME, BACKGROUNDS = menubar.THEME, menubar.BACKGROUNDS
#: what a click in the view selects
ATOM, RESIDUE, CHAIN = "atom", "residue", "chain"


class ProteinPanel(ToolPanel):
    def __init__(
        self,
        host: AssetHost,
        parent: QWidget | None = None,
        *,
        session: WorkbenchSession | None = None,
        current: Callable[[], int | None] | None = None,
        watch_current: Callable[[Callable[[int | None], None]], Callable[[], None]] | None = None,
        commands: menubar.Commands | None = None,
        run: Run | None = None,
        linked: LinkedRow | None = None,
    ) -> None:
        super().__init__(parent)
        self.tool_id = TOOL_ID
        self.linked = linked
        #: the panel made its session (no document keeps one for it): it closes it too
        self._owns_session = session is None
        if session is None:
            session = WorkbenchSession(
                host,
                deferred=run is not None,
                relabel=linked.relabel if linked is not None else None,
                count_named=linked.naming if linked is not None else None,
            )
        self.session = session
        #: the session's events as Qt signals, for this window's widgets
        self.signals = SessionSignals(session, self)
        #: builds the slow parts in the background while the window is open
        self.builder = Builder(session, run, self) if run is not None else None
        if self.builder is not None:
            session.schedule = self.builder.schedule
        self.view = None
        self.overlay = None
        self._view_page: QWidget | None = None
        self._problem_page: QLabel | None = None
        commands = dict(commands or {})
        self.commands = commands
        self._fitted: set[str] = set()  # proteins the camera has been fitted to
        self.background = THEME
        self.pocket_camera = PocketCamera(self)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self.selection_bar = SelectionBar(self.signals)

        self.tree = ProteinTree(self.signals)
        self.tree.focusRequested.connect(self.frame_selection)
        if REMOVE in self.commands:
            self.tree.removeRequested.connect(self.commands[REMOVE])
        self._sequence_window: SequenceWindow | None = None
        self._sequence_state: dict = {}  # saved, for the window when it is first made
        self.views = QStackedWidget()
        self.empty = self._empty_page(commands)
        self.views.addWidget(self.empty)
        self.models_bar = ModelsBar(self.signals)
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
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        self.menu_bar = menubar.build(self, commands)
        layout.addWidget(self.selection_bar)
        layout.addWidget(self.columns, 1)
        layout.addWidget(self.status)

        QShortcut(QKeySequence(Qt.Key.Key_Escape), self, self._escape)
        QShortcut(QKeySequence("["), self, lambda: self.models_bar.step(-1))
        QShortcut(QKeySequence("]"), self, lambda: self.models_bar.step(1))
        frame = QShortcut(QKeySequence(Qt.Key.Key_F), self, lambda: self.frame_selection())
        frame.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)  # text fields keep "f"
        self.signals.sceneChanged.connect(self._redraw)
        self.signals.lookChanged.connect(self._relook)
        self.signals.proteinsChanged.connect(self._proteins_changed)
        self.signals.hovered.connect(self._hovered)
        self._make_view()
        self._redraw()
        self._show_page()
        self._stop_watching = watch_current(self.follow_row) if watch_current else lambda: None
        self._stop_rows: Callable[[], None] = lambda: None
        if current is not None:
            self.follow_row(current())
            if linked is not None:
                self._stop_rows = linked.watch(lambda oids: self._row_edited(current(), oids))
            self.signals.edited.connect(lambda: self.follow_row(current()))

        QWidget.setTabOrder(self.selection_bar.field, self.tree)
        QWidget.setTabOrder(self.tree, self.views)

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange:
            self._apply_background()  # light or dark, when following the theme

    # how the view draws: background, projection, depth cue (View menu; saved)
    def set_background(self, key: str) -> None:
        self.background = key if key in BACKGROUNDS else THEME
        self._apply_background()

    def _apply_background(self) -> None:
        if self.view is not None:
            colour = BACKGROUNDS[self.background][1]
            self.view.set_background(colour or self.palette().base().color())

    @property
    def orthographic(self) -> bool:
        return self.view is not None and self.view.camera.orthographic

    def set_orthographic(self, on: bool) -> None:
        if self.view is not None:
            self.view.camera.orthographic = on
            self.view.camera_changed()

    @property
    def depth_cue(self) -> bool:
        return self.view is None or self.view.depth_cue

    def set_depth_cue(self, on: bool) -> None:
        if self.view is not None:
            self.view.set_depth_cue(on)

    def export_image(self) -> None:
        from .export import ExportImageDialog

        ExportImageDialog(self).open()

    def copy_image(self) -> bool:
        from .export import copy_image

        return copy_image(self)

    def _escape(self) -> None:
        if self.session.measuring:
            self.session.measure(0)
        else:
            self.session.clear_selection()

    def title(self) -> str:
        return "Protein Viewer"

    def close_viewer(self) -> None:
        """Close the viewer's window."""
        self.window().close()

    def toggle_full_screen(self) -> None:
        window = self.window()
        if window is self:
            return  # not in a window of its own
        window.showNormal() if window.isFullScreen() else window.showFullScreen()

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
        self.view.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._apply_background()
        self.view.failed.connect(lambda: self._no_view(graphics.NO_GPU))
        self.view.atomHovered.connect(self.session.hover)
        self.view.atomClicked.connect(self._clicked)
        self.view.atomDoubleClicked.connect(self._double_clicked)
        self.view.contextRequested.connect(self._context_menu)
        self.pocket_camera.watch(self.view)
        self.overlay = ProteinOverlay(self.view, self.signals)
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

    # picking in the view: an atom, its residue or its chain
    def _mask(self, atom: int, scope: str):
        found = self.session.scene.locate(atom)
        if found is None or found[0] is self.session.pose:  # a pose is not selected
            return None
        shown, index = found
        s = shown.structure
        if scope == ATOM and shown.styles[index] == scenes.HIDDEN:
            scope = RESIDUE  # a cartoon picked through its trace: no atom is drawn there
        if scope == ATOM:
            mask = np.zeros(s.n_atoms, bool)
            mask[index] = True
        elif scope == CHAIN:
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
        scope = CHAIN if modifiers & Qt.KeyboardModifier.AltModifier else ATOM
        masks = self._mask(atom, scope)
        if masks is None:
            if mode == SET:
                self.session.clear_selection()
            return
        self.session.select(masks, mode)

    def _double_clicked(self, atom: int) -> None:
        masks = self._mask(atom, RESIDUE)
        if masks is not None:
            self.session.select(masks, SET)

    # editing (Edit menu): each edit is one undoable step of the document
    def delete_selection(self) -> int:
        """Delete the selected atoms; how many went (the status line says why none did)."""
        named = sum(self.session.rows_naming(d.id, d.selected)
                    for d in self.session.displays if d.visible)  # fmt: skip
        try:
            count = self.session.delete_selected()
        except EditError as error:
            self.status.setText(f"⚠ {error}")
            return 0
        if count:
            said = f"Deleted {count} atom{'s' if count != 1 else ''}"
            if named:
                said += f"; {named} table row(s) name a residue that lost atoms"
            self.status.setText(said)
        return count

    def set_element(self, symbol: str) -> bool:
        return self._edit(self.session.set_element, symbol)

    def set_charge(self, charge: int) -> bool:
        return self._edit(self.session.set_charge, charge)

    def _edit(self, edit: Callable[..., None], *args) -> bool:
        try:
            edit(*args)
        except EditError as error:
            self.status.setText(f"⚠ {error}")
            return False
        return True

    def undo(self) -> str | None:
        """Undo the document's last step (an edit here, or in the table)."""
        label = self.session.store.undo()
        self.status.setText(f"Undid {label}" if label else "")
        return label

    def redo(self) -> str | None:
        label = self.session.store.redo()
        self.status.setText(f"Redid {label}" if label else "")
        return label

    def rename_residue(self):
        from . import edit_dialogs

        return edit_dialogs.rename_residue(self)

    def rename_chain(self):
        from . import edit_dialogs

        return edit_dialogs.rename_chain(self)

    def renumber(self):
        from . import edit_dialogs

        return edit_dialogs.renumber(self)

    # the Sequence window
    @property
    def sequence_window(self) -> SequenceWindow:
        """The Sequence window, made when first asked for."""
        if self._sequence_window is None:
            window = SequenceWindow(self.signals, self)
            window.strip.residueHovered.connect(self._sequence_hovered)
            for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
                QShortcut(QKeySequence(key), window, self.delete_selection)
            QShortcut(QKeySequence(QKeySequence.StandardKey.Undo), window, self.undo)
            QShortcut(QKeySequence(QKeySequence.StandardKey.Redo), window, self.redo)
            QShortcut(QKeySequence(Qt.Key.Key_Escape), window, self.session.clear_selection)
            self._sequence_window = window
            window.set_state(self._sequence_state)
        return self._sequence_window

    @property
    def sequence(self) -> SequenceView:
        return self.sequence_window.view

    def show_sequence(self, on: bool = True) -> None:
        window = self.sequence_window
        if not on:
            window.hide()
            return
        window.show()
        window.raise_()
        window.activateWindow()

    def sequence_shown(self) -> bool:
        return self._sequence_window is not None and self._sequence_window.isVisible()

    def sequence_hidden_kinds(self) -> set[int]:
        """The kind rows (waters, ions…) the Sequence window leaves out."""
        if self._sequence_window is not None:
            return set(self._sequence_window.strip.hidden_kinds)
        return {int(k) for k in self._sequence_state.get("hidden", []) if isinstance(k, int)}

    def hide_sequence_kind(self, kind: int, hidden: bool) -> None:
        self.sequence_window.set_hidden(kind, hidden)

    def _context_menu(self, atom: int, where) -> None:
        menu = menus.context_menu(
            self.session,
            atom,
            self,
            focus=self.frame_selection,
            copy_smiles=self.copy_smiles,
            copy_as=self.copy_selection,
            add_ligand=self.commands.get(LIGAND),
            panel=self,
        )
        menu.exec(where)

    # the table's current row: its ligand, in its pocket
    def follow_row(self, oid: int | None) -> bool:
        """Focus the pose or ligand that row *oid* names, if its protein is shown here."""
        if self.linked is None:
            return False
        if self.linked.follow(self.session, oid):
            return True
        missing = self.linked.missing(self.session, oid)
        if missing:
            self.status.setText(f"Protein {missing} is not in this document")
        return False

    def _row_edited(self, oid: int | None, oids: np.ndarray) -> None:
        if oid is not None and oid in oids:
            self.follow_row(oid)

    # copying
    def copy_selection(self, fmt: str) -> str:
        """Copy the selected atoms as PDB or mmCIF text (each protein with a selection)."""
        text = self.selection_text(fmt)
        if text:
            QGuiApplication.clipboard().setText(text)
        return text

    def selection_text(self, fmt: str) -> str:
        """The selected atoms of every shown protein as one PDB or mmCIF file; chains whose ID
        another protein took first are renamed, and the status line says so."""
        shown = [d for d in self.session.displays if d.visible and d.selected.any()]
        if not shown:
            return ""
        text, renamed = protein_io.to_text_many([(d.structure, d.selected) for d in shown], fmt)
        if renamed:
            said = ", ".join(f"{shown[k].name} {old} → {new}"
                             for (k, old), new in renamed.items())  # fmt: skip
            self.status.setText(f"Chains renamed: {said}")
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
            said = f"{shown.name}  {shown.structure.residue_label(residue)}"
        else:
            said = ""
        self.status.setText(said)
        if self._sequence_window is not None:
            self._sequence_window.status.setText(said)

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
    def view_state(self) -> dict:
        """How the view draws: background, depth cue and camera (the document's protein
        workspace keeps them)."""
        state = {"background": self.background, "depth_cue": self.depth_cue}
        if self.view is not None:
            state["camera"] = self.view.camera.state()  # orthographic or not, too
        return state

    def set_view_state(self, state: Mapping) -> None:
        self.set_background(str(state.get("background", THEME)))
        self.set_depth_cue(state.get("depth_cue") is not False)
        camera = state.get("camera")
        if isinstance(camera, dict) and self.view is not None:
            self.view.set_camera(Camera.from_state(camera))
            self._fitted = {d.id for d in self.session.displays if d.visible}

    def state(self) -> dict:
        """The window's layout (column sizes, the Sequence window). A panel with a session of
        its own adds the session and the view; a document's workspace saves those itself."""
        state = {
            "columns": self.columns.sizes(),
            "sequence": (
                self._sequence_window.state()
                if self._sequence_window is not None
                else self._sequence_state
            ),
        }
        if self._owns_session:
            state["session"] = self.session.state()
            state.update(self.view_state())
        return state

    def set_state(self, state: Mapping) -> None:
        session = state.get("session")
        if isinstance(session, dict):
            self.session.set_state(session)
        sizes = state.get("columns")
        if isinstance(sizes, list) and len(sizes) == self.columns.count():
            self.columns.setSizes([int(v) for v in sizes])
        sequence = state.get("sequence")
        self._sequence_state = dict(sequence) if isinstance(sequence, dict) else {}
        if self._sequence_window is not None or self._sequence_state.get("open") is True:
            self.sequence_window.set_state(self._sequence_state)
        if any(key in state for key in ("background", "depth_cue", "camera")):
            self.set_view_state(state)

    def close_panel(self) -> None:
        """Let go of the session and the document: its window calls this when closed, then
        deletes the panel, so nothing may still call into it. A session the document keeps
        lives on, for the window to attach to when it opens again."""
        self._stop_watching()
        self._stop_watching = lambda: None
        self._stop_rows()
        self._stop_rows = lambda: None
        if self._sequence_window is not None:
            self._sequence_state = self._sequence_window.state()
            self._sequence_window.hide()
        self.signals.close()
        if self.builder is not None:
            self.builder.cancel_all()
            if self.session.schedule == self.builder.schedule:
                self.session.schedule = None
        if self._owns_session:
            self.session.close()

    def closeEvent(self, event) -> None:
        self.close_panel()
        super().closeEvent(event)
