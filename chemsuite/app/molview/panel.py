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
"""The 3D Viewer panel: the selected rows' structures in 3D.

It browses the table one row at a time, in the table's order: the selected rows, the shown
rows, both (selected and shown) or, with neither, every row. It moves the table's current row
along and follows it when the table's current row moves. A row's conformers are shown one at a
time (step through them) or overlaid. Rows whose structure has no 3D coordinates get a quick
ETKDG conformer in the background. Picking two to four atoms in measure mode gives a distance,
angle or dihedral.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from PySide6.QtCore import QEvent, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QLabel,
    QSizePolicy,
    QToolBar,
    QToolTip,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

from chemsuite.core.selection import Selection
from chemsuite.core.table import Changed, Table
from chemsuite.core.table.columns import STRUCTURE
from chemsuite.molview.align import aligned_conformers
from chemsuite.molview.geometry import (
    ALL_H,
    BALL_AND_STICK,
    HYDROGENS,
    NO_H,
    POLAR_H,
    STYLES,
    THIN,
    WIRE,
    Scene,
    combine,
    molecule_part,
)

from ..glyphs import glyph
from ..panel_dock import OptionsPopup
from .overlay import MolOverlay
from .renderer import MolRenderer

MAX_PARTS = 200  # conformers drawn at once in an overlay
REFRESH_MS = 60
ELEMENT = "element"
BY_ROW = "row"
COLOURS = {ELEMENT: "Element", BY_ROW: "Carbons by Conformer"}
#: carbon colours for overlays (distinct on white and on dark backgrounds)
PALETTE = [
    (0.55, 0.55, 0.55), (0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17),
    (0.58, 0.4, 0.74), (0.55, 0.34, 0.29), (0.89, 0.47, 0.76), (0.74, 0.74, 0.13),
    (0.09, 0.75, 0.81), (0.84, 0.15, 0.16),
]  # fmt: skip


def _is_3d(mol: Chem.Mol | None) -> bool:
    return mol is not None and mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def make_3d(mols: dict[int, Chem.Mol]) -> dict[int, Chem.Mol | None]:
    """A quick conformer for each structure that has none (None when embedding fails)."""
    from chemsuite.conformers.etkdg import ConformerSettings, generate

    settings = ConformerSettings(count=1, energy_window=None, timeout=10)
    out: dict[int, Chem.Mol | None] = {}
    for oid, mol in mols.items():
        try:
            out[oid] = generate(mol, settings).mol
        except Exception:  # noqa: BLE001 - one bad structure must not stop the rest
            out[oid] = None
    return out


class MolViewPanel(QWidget):
    titleChanged = Signal(str)

    def __init__(
        self,
        table: Table,
        selection: Selection,
        *,
        current: Callable[[], int | None] = lambda: None,
        run: Callable | None = None,
        move_to: Callable[[int], bool] | None = None,
        watch_current: Callable[[Callable[[int | None], None]], Callable[[], None]] | None = None,
    ) -> None:
        """*current* gives the table's current row, *move_to* moves it and *watch_current*
        calls back when it moves (returning the function that stops that)."""
        super().__init__()
        self.table = table
        self.selection = selection
        self.tool_id = "molview"
        self._current = current
        self._run = run
        self._move_to = move_to
        self._oids = np.empty(0, dtype=np.int64)  # the browsed rows, in the table's order
        self._index = -1
        self.oid: int | None = None  # the row shown
        self._scope_stale = False  # the table, selection or browse options changed
        self._origin = f"molview:{id(self)}"
        self._made: dict[tuple[str, int], Chem.Mol | None] = {}  # conformers made here
        self._waiting: set[tuple[str, int]] = set()
        self._shown: tuple = ()  # (column, oid, drawn) last shown: a new one is fitted
        self._look: tuple = ()  # (style, hydrogens, overlay) last drawn
        self._conformers = 0
        self.conformer = 0
        self.measuring = False
        self._marked: list[int] = []

        self.renderer = MolRenderer()
        self.renderer.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.setFocusProxy(self.renderer)
        self.overlay = MolOverlay(self.renderer)
        self.problem = QLabel()
        self.problem.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.problem.setWordWrap(True)
        self.problem.hide()
        self.renderer.failed.connect(self._failed)
        self.renderer.atomHovered.connect(self._hovered)
        self.renderer.atomClicked.connect(self._clicked)
        self.renderer.cameraMoved.connect(self.overlay.update)
        stack = QGridLayout()
        stack.setContentsMargins(0, 0, 0, 0)
        stack.addWidget(self.renderer, 0, 0)
        stack.addWidget(self.overlay, 0, 0)
        stack.addWidget(self.problem, 0, 0)

        self._build_options()
        self._build_toolbar()
        frame = QFrame()
        frame.setObjectName("MolViewFrame")
        inside = QVBoxLayout(frame)
        inside.setContentsMargins(0, 0, 0, 0)
        inside.setSpacing(0)
        inside.addWidget(self.toolbar)
        inside.addLayout(stack, 1)
        inside.addWidget(self.rowbar)
        self._frame = frame
        self._header_attached = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(frame)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(REFRESH_MS)
        self._timer.timeout.connect(self.refresh)
        self._unsubscribe_table = table.subscribe(self._on_table_event)
        self._unsubscribe_selection = selection.subscribe(self._on_selection)
        self._stop_watching = watch_current(self._on_current) if watch_current else lambda: None
        self._fill_columns()
        self._restyle()
        self._rescope(restart=True)
        self.refresh()

    # -- controls -------------------------------------------------------------------------------
    def _build_options(self) -> None:
        self.column = QComboBox()
        self.column.setToolTip("The structure column to show")
        self.style_box = self._combo(STYLES, BALL_AND_STICK, "How atoms and bonds are drawn")
        self.hydrogens = self._combo(HYDROGENS, ALL_H, "Which hydrogens to show")
        self.colour = self._combo(
            COLOURS,
            ELEMENT,
            "Element colours, or each conformer's (or row's) carbons in a colour of its own",
        )
        self.browse_selected = QCheckBox("Browse Selected")
        self.browse_selected.setToolTip("Step through the selected rows only")
        self.browse_visible = QCheckBox("Browse Visible")
        self.browse_visible.setChecked(True)
        self.browse_visible.setToolTip(
            "Step through the rows the filters and search show only; with both off, every row"
        )
        self.align = QCheckBox("Align Conformers")
        self.align.setChecked(True)
        self.align.setToolTip("Fit every conformer onto the first, on the heavy atoms")
        self.options = OptionsPopup(self)
        form = self.options.form
        form.addRow("", self.browse_selected)
        form.addRow("", self.browse_visible)
        form.addRow("Structure:", self.column)
        form.addRow("Style:", self.style_box)
        form.addRow("Hydrogens:", self.hydrogens)
        form.addRow("Color:", self.colour)
        form.addRow("", self.align)
        for box in (self.column, self.style_box, self.hydrogens, self.colour):
            box.currentIndexChanged.connect(self._schedule)
        self.align.toggled.connect(self._schedule)
        self.browse_selected.toggled.connect(self._rescope_later)
        self.browse_visible.toggled.connect(self._rescope_later)

    def _combo(self, choices: dict[str, str], current: str, tip: str) -> QComboBox:
        box = QComboBox()
        for key, text in choices.items():
            box.addItem(text, key)
        box.setCurrentIndex(box.findData(current))
        box.setToolTip(tip)
        return box

    def _build_toolbar(self) -> None:
        bar = QToolBar()
        bar.setIconSize(QSize(15, 15))
        bar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        bar.setMovable(False)
        self.toolbar = bar
        self._glyphs: dict[QAction, str] = {}

        def add(
            name: str, text: str, tip: str, slot, *, checkable: bool = False, on=bar
        ) -> QAction:
            action = on.addAction(text)
            action.setToolTip(f"{text}: {tip}")
            action.setCheckable(checkable)
            (action.toggled if checkable else action.triggered).connect(slot)
            self._glyphs[action] = name
            return action

        # Rows are stepped through under the view, centred, as in the 2D browser.
        rowbar = QToolBar()
        rowbar.setIconSize(QSize(15, 15))
        rowbar.setMovable(False)
        self.rowbar = rowbar
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        rowbar.addWidget(spacer)
        self.previous_row = add(
            "previous",
            "Previous Row",
            "one back (←; Home for the first)",
            lambda: self.go_by(-1),
            on=rowbar,
        )
        self.row_counter = QLabel()
        self.row_counter.setToolTip("Row shown, of the rows browsed")
        self.row_counter.setContentsMargins(4, 0, 4, 0)
        rowbar.addWidget(self.row_counter)
        self.next_row = add(
            "next", "Next Row", "one on (→; End for the last)", lambda: self.go_by(1), on=rowbar
        )
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        rowbar.addWidget(spacer)
        self._conformer_actions = []
        self.previous = add("previous", "Previous Conformer", "one back (↑)", lambda: self.step(-1))
        self.counter = QLabel()
        self.counter.setToolTip("Conformer shown, of how many")
        self.counter.setContentsMargins(4, 0, 4, 0)
        self._conformer_actions += [self.previous, bar.addWidget(self.counter)]
        self.next = add("next", "Next Conformer", "one on (↓)", lambda: self.step(1))
        self._conformer_actions += [self.next, bar.addSeparator()]
        self.overlay_action = add(
            "overlay",
            "Overlay",
            "every conformer at once, as thin sticks",
            self._schedule,
            checkable=True,
        )
        self.measure_action = add(
            "ruler",
            "Measure",
            "click two, three or four atoms for a distance, angle or dihedral",
            self.set_measuring,
            checkable=True,
        )
        bar.addSeparator()
        add(
            "fit",
            "Reset View",
            "turn the molecule square on and fit it to the view",
            self.renderer.reset_view,
        )
        add("export", "Export…", "save the view as a PNG image", self.export_dialog)

    # -- what is shown --------------------------------------------------------------------------
    def column_name(self) -> str | None:
        name = self.column.currentText()
        return name if name in self.table and self.table.kind(name) == STRUCTURE else None

    def style(self) -> str:
        return self.style_box.currentData()

    def rows(self) -> list[int]:
        """The rows to draw: the one browsed to."""
        return [] if self.oid is None else [self.oid]

    # -- browsing -------------------------------------------------------------------------------
    @property
    def browsed(self) -> np.ndarray:
        """The rows stepped through, in the table's order."""
        return self._oids

    def browsed_oids(self) -> np.ndarray:
        if self.browse_visible.isChecked():
            order = self.table.visible_oids()
        else:
            order = self.table.ordered_oids()
        if self.browse_selected.isChecked():
            order = order[np.isin(order, self.selection.as_array())]
        return order

    def _rescope(self, *, restart: bool = False) -> None:
        """Recompute the browsed rows, staying on the same row while it is one of them.

        A row that left (deselected while browsing the selection, say) leaves the view at its
        place in the list. *restart* begins at the table's current row when it is browsed.
        """
        self._scope_stale = False
        self._oids = self.browsed_oids()
        n = len(self._oids)
        anchor = self._current() if restart else self.oid
        hit = np.flatnonzero(self._oids == anchor) if anchor is not None and n else []
        if len(hit):
            index = int(hit[0])
        elif restart or self.oid is None:
            index = 0 if n else -1
        else:
            index = min(max(self._index, 0), n - 1)
        self._set_index(index)

    def _set_index(self, index: int) -> None:
        self._index = index
        oid = int(self._oids[index]) if index >= 0 else None
        if oid != self.oid:
            self.oid = oid
            self.conformer = 0

    def go(self, index: int) -> None:
        """Show the *index*-th browsed row (clamped) and make it the table's current row."""
        if not len(self._oids):
            return
        self._set_index(max(0, min(index, len(self._oids) - 1)))
        if self._move_to is not None and self.oid is not None:
            self._move_to(self.oid)
        self.refresh()

    def go_by(self, delta: int) -> None:
        self.go(self._index + delta)

    def _on_current(self, oid: int | None) -> None:
        """The table's current row moved: show it when it is browsed."""
        if oid is None or oid == self.oid:
            return
        hit = np.flatnonzero(self._oids == oid)
        if len(hit):
            self._set_index(int(hit[0]))
            self.refresh()

    def molecules(self) -> list[tuple[int, Chem.Mol]]:
        """(oid, 3D molecule) for each row ready to draw; starts making 3D for the others."""
        name = self.column_name()
        if name is None:
            return []
        column = self.table.column(name)
        rows = self.rows()
        out, missing = [], {}
        for oid, slot in zip(rows, self.table.slots_of(rows).tolist(), strict=True):
            mol = column.mol(int(slot))
            if mol is None:
                continue
            if not _is_3d(mol):
                key = (name, oid)
                if key in self._made:
                    mol = self._made[key]
                elif key not in self._waiting:
                    missing[oid] = mol
                    continue
                else:
                    continue
            if mol is not None:
                out.append((oid, mol))
        if missing:
            self._make_3d(name, missing)
        return out

    def _make_3d(self, name: str, mols: dict[int, Chem.Mol]) -> None:
        keys = {(name, oid) for oid in mols}
        self._waiting |= keys
        if self._run is None:
            self._made_3d(name, keys, make_3d(mols))
            return
        self._run(
            "Make 3D structures",
            lambda _job: make_3d(mols),
            lambda result: self._made_3d(name, keys, result),
        )

    def _made_3d(self, name: str, keys: set, result: dict) -> None:
        self._waiting -= keys
        for oid, mol in (result or {}).items():
            self._made[(name, oid)] = mol
        self._schedule()

    def scene_parts(self, molecules: list[tuple[int, Chem.Mol]]):
        style, hydrogens = self.style(), self.hydrogens.currentData()
        by_row = self.colour.currentData() == BY_ROW
        overlay = self.overlay_action.isChecked()
        if overlay:  # as 3Dmol overlays read best: thin sticks, no hydrogens on carbon
            style = WIRE if style == WIRE else THIN
            hydrogens = NO_H if hydrogens == NO_H else POLAR_H
        counts = [m.GetNumConformers() for _, m in molecules]
        self._conformers = max(counts, default=0)
        if counts:  # nothing shown yet (a document reopening) keeps the saved conformer
            self.conformer = min(self.conformer, self._conformers - 1)
        parts, rows = [], []
        for oid, mol in molecules:
            if self.align.isChecked() and mol.GetNumConformers() > 1:
                mol = aligned_conformers(mol)
            ids = [c.GetId() for c in mol.GetConformers()]
            chosen = ids if overlay else [ids[min(self.conformer, len(ids) - 1)]]
            for conf_id in chosen:
                if len(parts) >= MAX_PARTS:
                    break
                carbon = PALETTE[len(parts) % len(PALETTE)] if by_row else None
                parts.append(
                    molecule_part(mol, conf_id, style=style, hydrogens=hydrogens, carbon=carbon)
                )
                rows.append(oid)
        return parts, rows

    def refresh(self) -> None:
        self._timer.stop()
        if self._scope_stale:
            self._rescope()
        molecules = self.molecules()
        parts, _rows = self.scene_parts(molecules)
        scene = combine(parts) if parts else Scene.empty()
        # A new molecule (another row or column, or its 3D structure just made) is turned to
        # its principal axes, centred and fitted; new drawing settings refit only a view the
        # user has not moved. Stepping through conformers keeps the view, to compare them.
        shown = (self.column_name(), self.oid, bool(parts))
        look = (self.style(), self.hydrogens.currentData(), self.overlay_action.isChecked())
        new = shown != self._shown
        restyled = look != self._look
        self._shown, self._look = shown, look
        if new:
            self._clear_marks()
        self.renderer.set_scene(scene, orient=new, fit=restyled and self.renderer.auto_fitting)
        self._update_counter()
        self.overlay.update()

    def step(self, delta: int) -> None:
        """Show another conformer (clamped to the first and last)."""
        target = max(0, min(self.conformer + delta, max(self._conformers - 1, 0)))
        if target != self.conformer:
            self.conformer = target
            self._clear_marks()
            self.refresh()

    def _update_counter(self) -> None:
        n = len(self._oids)
        self.row_counter.setText(f"Row {self._index + 1:,} of {n:,}" if n else "No Rows")
        self.previous_row.setEnabled(self._index > 0)
        self.next_row.setEnabled(0 <= self._index < n - 1)
        total = self._conformers
        overlay = self.overlay_action.isChecked()
        self.counter.setText(
            f"{total} Conformers" if overlay
            else f"{self.conformer + 1} / {total}"
        )  # fmt: skip
        for action in self._conformer_actions:  # only for a row with several
            action.setVisible(total > 1)
        stepping = total > 1 and not overlay
        self.previous.setEnabled(stepping and self.conformer > 0)
        self.next.setEnabled(stepping and self.conformer < total - 1)

    def _schedule(self, *_args) -> None:
        self._timer.start()

    def _rescope_later(self, *_args) -> None:
        self._scope_stale = True
        self._timer.start()

    # -- picking ------------------------------------------------------------------------------
    def set_measuring(self, on: bool) -> None:
        self.measuring = bool(on)
        if self.measure_action.isChecked() != self.measuring:
            self.measure_action.setChecked(self.measuring)
        self._clear_marks()

    def _clear_marks(self) -> None:
        self._marked = []
        self.overlay.set_marked([])

    def _hovered(self, atom: int) -> None:
        self.overlay.set_hover(atom)
        if atom < 0:
            QToolTip.hideText()
            return
        QToolTip.showText(self.cursor().pos(), self.describe_atom(atom), self.renderer)

    def describe_atom(self, atom: int) -> str:
        scene = self.renderer.scene
        symbol = Chem.GetPeriodicTable().GetElementSymbol(int(scene.atom_number[atom]))
        return f"{symbol}{int(scene.atom_index[atom]) + 1}"

    def _clicked(self, atom: int, modifiers) -> None:
        if self.measuring:
            if atom < 0:
                return
            marked = [*self._marked, atom] if len(self._marked) < 4 else [atom]
            self._marked = marked
            self.overlay.set_marked(marked)
            return
        oid = self.oid
        if atom < 0 or oid is None:
            return
        if modifiers & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier):
            self.selection.set(self.selection.oids ^ {oid}, self._origin)

    def measurement(self) -> str:
        return self.overlay.measurement()

    # -- table and selection --------------------------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        if event.structural:
            self._fill_columns()
        self._made = {k: v for k, v in self._made.items() if k[0] in self.table}
        self._rescope_later()

    def _on_selection(self, _oids: frozenset[int], _origin: str) -> None:
        if self.browse_selected.isChecked():
            self._rescope_later()

    def _fill_columns(self) -> None:
        names = [n for n in self.table.column_names if self.table.kind(n) == STRUCTURE]
        current = self.column.currentText()
        blocked = self.column.blockSignals(True)
        self.column.clear()
        self.column.addItems(names)
        if current in names:
            self.column.setCurrentText(current)
        else:
            three_d = next((n for n in names if "3D" in n), None)
            self.column.setCurrentText(three_d or (names[0] if names else ""))
        self.column.blockSignals(blocked)

    # -- panel ----------------------------------------------------------------------------------
    def sizeHint(self) -> QSize:
        return QSize(520, 480)

    def keyPressEvent(self, event) -> None:
        """Left and right step through the rows, Home and End go to the first and last; up and
        down step through the conformers."""
        rows = {Qt.Key.Key_Left: -1, Qt.Key.Key_Right: 1}
        conformers = {Qt.Key.Key_Up: -1, Qt.Key.Key_Down: 1}
        key = event.key()
        if key in rows:
            self.go_by(rows[key])
        elif key == Qt.Key.Key_Home:
            self.go(0)
        elif key == Qt.Key.Key_End:
            self.go(len(self._oids) - 1)
        elif key in conformers:
            self.step(conformers[key])
        else:
            super().keyPressEvent(event)

    def title(self) -> str:
        return "3D Viewer"

    def state(self) -> dict:
        return {
            "browse_selected": self.browse_selected.isChecked(),
            "browse_visible": self.browse_visible.isChecked(),
            "row": self.oid,
            "column": self.column.currentText(),
            "style": self.style(),
            "hydrogens": self.hydrogens.currentData(),
            "colour": self.colour.currentData(),
            "align": self.align.isChecked(),
            "overlay": self.overlay_action.isChecked(),
            "conformer": self.conformer,
            "camera": self.renderer.camera.state(),
        }

    def set_state(self, state: dict) -> None:
        from chemsuite.molview.camera import Camera

        if state.get("column") in [self.column.itemText(i) for i in range(self.column.count())]:
            self.column.setCurrentText(state["column"])
        for box, key in (
            (self.style_box, "style"),
            (self.hydrogens, "hydrogens"),
            (self.colour, "colour"),
        ):
            index = box.findData(state.get(key))
            if index >= 0:
                box.setCurrentIndex(index)
        self.align.setChecked(bool(state.get("align", True)))
        self.browse_selected.setChecked(bool(state.get("browse_selected", False)))
        self.browse_visible.setChecked(bool(state.get("browse_visible", True)))
        self._rescope()
        row = state.get("row")
        hit = np.flatnonzero(self._oids == row) if isinstance(row, int) else []
        if len(hit):
            self._set_index(int(hit[0]))
        self.overlay_action.setChecked(bool(state.get("overlay", False)))
        self.conformer = max(0, int(state.get("conformer", 0) or 0))
        self.refresh()
        if self._move_to is not None and self.oid is not None:
            self._move_to(self.oid)
        if isinstance(state.get("camera"), dict):
            self.renderer.set_camera(Camera.from_state(state["camera"]))

    def export_dialog(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Export View", "3d-view.png", "PNG (*.png)")
        if path:
            self.export(Path(path))

    def export(self, path: Path) -> bool:
        image = self.renderer.grabFramebuffer()
        return not image.isNull() and image.save(str(path.with_suffix(".png")))

    def attach_header(self, attached: bool) -> None:
        self._header_attached = attached
        self._restyle()

    def close_panel(self) -> None:
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self._stop_watching()
        self._timer.stop()
        self.options.hide()

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "toolbar"):
            self._restyle()

    def _restyle(self) -> None:
        palette = self.palette()
        ink = palette.text().color()
        for action, name in self._glyphs.items():
            action.setIcon(glyph(name, ink))
        base = palette.base().color()
        self.renderer.set_background(base)
        line = palette.mid().color().name()
        hover = QColor(ink)
        hover.setAlphaF(0.08)
        checked = QColor(palette.highlight().color())
        checked.setAlphaF(0.28)
        rgba = "rgba({}, {}, {}, {:.2f})".format
        if self._header_attached:
            bar = "background: transparent; border: none; padding: 0px; spacing: 2px;"
        else:
            bar = f"background: {base.name()}; border: none; border-bottom: 1px solid {line};"
            bar += " padding: 2px 4px; spacing: 2px;"
        self._frame.setStyleSheet(
            f"QFrame#MolViewFrame {{ border: none; background: {base.name()}; }}"
        )
        buttons = (
            " QToolButton { border: none; border-radius: 4px; padding: 2px; }"
            f" QToolButton:hover {{ background: {rgba(*hover.getRgb()[:3], hover.alphaF())}; }}"
            " QToolButton:checked {"
            f" background: {rgba(*checked.getRgb()[:3], checked.alphaF())}; }}"
        )
        self.toolbar.setStyleSheet(f"QToolBar {{ {bar} }}" + buttons)
        self.rowbar.setStyleSheet(
            f"QToolBar {{ background: {base.name()}; border: none;"
            f" border-top: 1px solid {line}; padding: 2px 4px; spacing: 2px; }}" + buttons
        )

    def _failed(self) -> None:
        self.renderer.hide()
        self.overlay.hide()
        self.problem.setText("⚠ 3D drawing is not available on this computer's graphics")
        self.problem.show()
