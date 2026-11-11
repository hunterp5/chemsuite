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
"""The Structure Viewer panel: the table's structures one row at a time, in 2D or in 3D.

It browses the table one row at a time, in the table's order: the selected rows, the shown
rows, both (selected and shown) or, with neither, every row. It moves the table's current row
along and follows it when the table's current row moves. The 2D / 3D switch only changes how
the row is drawn: the row, the column, the browsed rows and the caption (chosen column values
drawn over the structure) stay as they are.

In 3D, a row's conformers are shown one at a time (step through them) or overlaid, and picking
two to four atoms in measure mode gives a distance, angle or dihedral. Rows whose structure
has no 3D coordinates get a quick ETKDG conformer in the background, made ahead for the rows
around the one shown so switching to 3D or stepping on finds it ready. The 3D view is only
created the first time it is shown, so the 2D view never needs the graphics card; where 3D
drawing is not available the switch stays on 2D.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from PySide6.QtCore import QEvent, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QColor, QImage, QPainter
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QLabel,
    QMenu,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QToolBar,
    QToolButton,
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

from .. import graphics
from ..glyphs import glyph
from ..panel_dock import OptionsPopup
from .caption import CaptionOverlay, caption_columns, default_caption, paint_caption
from .depiction import Depiction2D, structure_image

TWO_D, THREE_D = "2d", "3d"
MODES = {TWO_D: "2D", THREE_D: "3D"}
MAX_PARTS = 200  # conformers drawn at once in an overlay
REFRESH_MS = 60
PREFETCH_MS = 300  # resting on a row this long makes 3D for the rows around it
PREFETCH_AROUND = 2  # rows either side made ahead
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


class StructureViewer(QWidget):
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
        mode: str = TWO_D,
    ) -> None:
        """*current* gives the table's current row, *move_to* moves it and *watch_current*
        calls back when it moves (returning the function that stops that)."""
        super().__init__()
        self.table = table
        self.selection = selection
        self.tool_id = "viewer"
        self._current = current
        self._run = run
        self._move_to = move_to
        self.mode = mode if mode in MODES else TWO_D  # (set_mode below makes it so)
        self._oids = np.empty(0, dtype=np.int64)  # the browsed rows, in the table's order
        self._index = -1
        self.oid: int | None = None  # the row shown
        self._scope_stale = False  # the table, selection or browse options changed
        self._origin = f"viewer:{id(self)}"
        self._made: dict[tuple[str, int], Chem.Mol | None] = {}  # conformers made here
        self._waiting: set[tuple[str, int]] = set()
        self._ahead: set[tuple[str, int]] = set()  # being made ahead (one batch at a time)
        self._shown: tuple = ()  # (column, oid, drawn) last shown in 3D: a new one is fitted
        self._look: tuple = ()  # (style, hydrogens, overlay) last drawn
        self._drawn_2d: tuple = ()  # (column, oid) last drawn in 2D
        self._conformers = 0
        self.conformer = 0
        self.measuring = False
        self._marked: list[int] = []
        self.caption_names: list[str] = default_caption(table)
        self.renderer = None  # the 3D view, made the first time it is shown
        self.overlay = None
        self._camera_state: dict | None = None  # a saved camera, for when 3D is made
        self._pending_camera = None  # set once its molecule is drawn
        self.gpu_failed = False

        self.depiction = Depiction2D(lambda: self._draw_2d(force=True))
        self.pages = QStackedWidget()
        self.pages.addWidget(self.depiction)
        self.caption = CaptionOverlay()
        stack = QGridLayout()
        stack.setContentsMargins(0, 0, 0, 0)
        stack.addWidget(self.pages, 0, 0)
        stack.addWidget(self.caption, 0, 0)

        self._build_options()
        self._build_toolbar()
        self._build_rowbar()
        frame = QFrame()
        frame.setObjectName("ViewerFrame")
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
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(REFRESH_MS)
        self._timer.timeout.connect(self.refresh)
        self._ahead_timer = QTimer(self)
        self._ahead_timer.setSingleShot(True)
        self._ahead_timer.setInterval(PREFETCH_MS)
        self._ahead_timer.timeout.connect(self.make_ahead)
        self._unsubscribe_table = table.subscribe(self._on_table_event)
        self._unsubscribe_selection = selection.subscribe(self._on_selection)
        self._stop_watching = watch_current(self._on_current) if watch_current else lambda: None
        self._fill_columns()
        self._restyle()
        self._rescope(restart=True)
        if reason := graphics.problem():  # known at start-up: never try the GPU
            self._failed(reason)
        self.set_mode(mode)

    # -- controls -------------------------------------------------------------------------------
    def _build_options(self) -> None:
        self.column = QComboBox()
        self.column.setToolTip("The structure column to show")
        self.caption_button = QToolButton()
        self.caption_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.caption_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.caption_button.setToolTip("Column values drawn over the structure")
        self.caption_menu = QMenu(self.caption_button)
        self.caption_menu.aboutToShow.connect(self._fill_caption_menu)
        self.caption_button.setMenu(self.caption_menu)
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
        form.addRow("Caption:", self.caption_button)
        form.addRow("Style:", self.style_box)
        form.addRow("Hydrogens:", self.hydrogens)
        form.addRow("Color:", self.colour)
        form.addRow("", self.align)
        self._three_d_rows = [self.style_box, self.hydrogens, self.colour, self.align]
        for box in (self.style_box, self.hydrogens, self.colour):
            box.currentIndexChanged.connect(self._schedule)
        self.column.currentIndexChanged.connect(self._schedule)
        self.align.toggled.connect(self._schedule)
        self.browse_selected.toggled.connect(self._rescope_later)
        self.browse_visible.toggled.connect(self._rescope_later)
        self._show_caption_names()

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
        self.mode_buttons: dict[str, QToolButton] = {}
        for key, text, tip in (
            (TWO_D, "2D", "Draw the structure flat, as chemists draw it (Ctrl+Shift+B)"),
            (THREE_D, "3D", "Draw the structure in 3D: turn, step through conformers, measure"),
        ):
            button = QToolButton()
            button.setText(text)
            button.setCheckable(True)
            button.setAutoExclusive(True)
            button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
            button.setToolTip(tip)
            button.clicked.connect(lambda _on=False, k=key: self.set_mode(k))
            bar.addWidget(button)
            self.mode_buttons[key] = button
        self._three_d_actions: list[QAction] = [bar.addSeparator()]

        def add(name: str, text: str, tip: str, slot, *, checkable: bool = False) -> QAction:
            action = bar.addAction(text)
            action.setToolTip(f"{text}: {tip}")
            action.setCheckable(checkable)
            (action.toggled if checkable else action.triggered).connect(slot)
            self._glyphs[action] = name
            return action

        self.previous = add("previous", "Previous Conformer", "one back (↑)", lambda: self.step(-1))
        self.counter = QLabel()
        self.counter.setToolTip("Conformer shown, of how many")
        self.counter.setContentsMargins(4, 0, 4, 0)
        self._conformer_actions = [self.previous, bar.addWidget(self.counter)]
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
        self.reset_action = add(
            "fit", "Reset View", "turn the molecule square on and fit it to the view", self._reset
        )
        self._three_d_actions += [
            self.overlay_action,
            self.measure_action,
            self.reset_action,
            bar.addSeparator(),
        ]
        add("export", "Export…", "save the view as a PNG image", self.export_dialog)

    def _build_rowbar(self) -> None:
        """Rows are stepped through under the view, centred, with the row's selection tick box
        at the right."""
        rowbar = QToolBar()
        rowbar.setIconSize(QSize(15, 15))
        rowbar.setMovable(False)
        self.rowbar = rowbar
        self.select_box = QCheckBox("Selected")
        self.select_box.setToolTip("Add or remove this row from the selection (Space)")
        self.select_box.clicked.connect(self.set_selected)
        balance = QWidget()  # as wide as the tick box, so the buttons stay centred
        balance.setFixedWidth(self.select_box.sizeHint().width())
        rowbar.addWidget(balance)
        rowbar.addWidget(self._spacer())

        def add(name: str, text: str, tip: str, slot) -> QAction:
            action = rowbar.addAction(text)
            action.setToolTip(f"{text}: {tip}")
            action.triggered.connect(slot)
            self._glyphs[action] = name
            return action

        self.previous_row = add(
            "previous", "Previous Row", "one back (←; Home for the first)", lambda: self.go_by(-1)
        )
        self.position = QSpinBox()
        self.position.setKeyboardTracking(False)
        self.position.setMinimum(1)
        self.position.setButtonSymbols(QSpinBox.ButtonSymbols.NoButtons)
        self.position.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.position.setToolTip("Row shown, of the rows browsed; type a number to jump to it")
        self.position.valueChanged.connect(lambda value: self.go(value - 1))
        rowbar.addWidget(self.position)
        self.count_label = QLabel()
        self.count_label.setContentsMargins(4, 0, 4, 0)
        rowbar.addWidget(self.count_label)
        self.next_row = add(
            "next", "Next Row", "one on (→; End for the last)", lambda: self.go_by(1)
        )
        rowbar.addWidget(self._spacer())
        rowbar.addWidget(self.select_box)

    @staticmethod
    def _spacer() -> QWidget:
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        return spacer

    # -- 2D and 3D ------------------------------------------------------------------------------
    def set_mode(self, mode: str) -> None:
        """Draw in 2D or 3D; everything else (row, column, rows browsed, caption) stays."""
        if mode not in MODES:
            raise ValueError(f"unknown viewer mode {mode!r}")
        if mode == THREE_D and not self.gpu_failed:
            self._ensure_3d()
        if self.gpu_failed:
            mode = TWO_D
        self.mode = mode
        self.mode_buttons[mode].setChecked(True)
        three_d = mode == THREE_D
        if not three_d and self.measuring:
            self.set_measuring(False)
        self.pages.setCurrentIndex(1 if three_d else 0)
        self.setFocusProxy(self.renderer if three_d else None)
        for action in self._three_d_actions:
            action.setVisible(three_d)
        for widget in self._three_d_rows:
            self.options.form.setRowVisible(widget, three_d)
        self.refresh()

    def _ensure_3d(self) -> None:
        if self.renderer is not None:
            return
        from chemsuite.molview.camera import Camera

        from ..molview.overlay import MolOverlay
        from ..molview.renderer import MolRenderer

        self.renderer = MolRenderer()
        self.renderer.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.overlay = MolOverlay(self.renderer)
        self.renderer.failed.connect(self._failed)
        self.renderer.atomHovered.connect(self._hovered)
        self.renderer.atomClicked.connect(self._clicked)
        self.renderer.cameraMoved.connect(self.overlay.update)
        page = QWidget()
        grid = QGridLayout(page)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.addWidget(self.renderer, 0, 0)
        grid.addWidget(self.overlay, 0, 0)
        self.pages.addWidget(page)
        self.renderer.set_background(self.palette().base().color())
        if self._camera_state is not None:
            self._pending_camera = Camera.from_state(self._camera_state)
            self._camera_state = None

    def _failed(self, reason: str = graphics.NO_GPU) -> None:
        """The graphics cannot draw 3D: stay in 2D and say why on the 3D switch."""
        self.gpu_failed = True
        self.mode_buttons[THREE_D].setEnabled(False)
        self.mode_buttons[THREE_D].setToolTip(f"⚠ {reason}")
        if self.mode == THREE_D:
            self.set_mode(TWO_D)

    # -- what is shown --------------------------------------------------------------------------
    def column_name(self) -> str | None:
        name = self.column.currentText()
        return name if name in self.table and self.table.kind(name) == STRUCTURE else None

    def style(self) -> str:
        return self.style_box.currentData()

    def rows(self) -> list[int]:
        """The rows to draw: the one browsed to."""
        return [] if self.oid is None else [self.oid]

    def caption_lines(self) -> list[tuple[str, str]]:
        """(column, value as the table shows it) for the caption, for the row shown."""
        if self.oid is None:
            return []
        return [
            (name, self.table.cell(name, self.oid))
            for name in self.caption_names
            if name in self.table
        ]

    def set_caption(self, names: list[str]) -> None:
        self.caption_names = [n for n in dict.fromkeys(names) if isinstance(n, str)]
        self._show_caption_names()
        self.caption.set_lines(self.caption_lines())

    def _show_caption_names(self) -> None:
        shown = [n for n in self.caption_names if n in self.table]
        self.caption_button.setText(", ".join(shown) if shown else "None")

    def _fill_caption_menu(self) -> None:
        menu = self.caption_menu
        menu.clear()
        none = menu.addAction("None")
        none.setCheckable(True)
        none.setChecked(not any(n in self.table for n in self.caption_names))
        none.triggered.connect(lambda: self.set_caption([]))
        menu.addSeparator()
        for name in caption_columns(self.table):
            action = menu.addAction(name)
            action.setCheckable(True)
            action.setChecked(name in self.caption_names)
            action.toggled.connect(lambda on, n=name: self._toggle_caption(n, on))

    def _toggle_caption(self, name: str, on: bool) -> None:
        names = [n for n in self.caption_names if n != name and n in self.table]
        if on:  # in table order, as in the menu
            order = caption_columns(self.table)
            names = [n for n in order if n in {*names, name}]
        self.set_caption(names)

    # -- browsing -------------------------------------------------------------------------------
    @property
    def browsed(self) -> np.ndarray:
        """The rows stepped through, in the table's order."""
        return self._oids

    @property
    def index(self) -> int:
        return self._index

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

    def set_selected(self, on: bool) -> None:
        """Add the row shown to the selection, or take it out."""
        if self.oid is None:
            return
        if on:
            self.selection.add([self.oid], self._origin)
        else:
            self.selection.remove([self.oid], self._origin)

    # -- 3D structures --------------------------------------------------------------------------
    def _structure(self, name: str, oid: int) -> Chem.Mol | None:
        return self.table.column(name).mol(int(self.table.slots_of([oid])[0]))

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

    def needs_3d(self, oids) -> dict[int, Chem.Mol]:
        """The rows among *oids* with a 2D structure and no 3D one made or being made."""
        name = self.column_name()
        if name is None:
            return {}
        out = {}
        for oid in oids:
            key = (name, int(oid))
            if key in self._made or key in self._waiting:
                continue
            mol = self._structure(name, int(oid))
            if mol is not None and not _is_3d(mol):
                out[int(oid)] = mol
        return out

    def make_ahead(self) -> None:
        """Make 3D for the rows around the one shown, one row at a time and the row shown
        first, so the 3D view finds them ready. Stepping on while one is made only changes
        which row comes next."""
        self._ahead_timer.stop()
        name = self.column_name()
        if self.gpu_failed or name is None or self._index < 0 or self._ahead:
            return
        low = max(0, self._index - PREFETCH_AROUND)
        around = self._oids[low : self._index + PREFETCH_AROUND + 1]
        # the row shown first, then outwards
        around = sorted(around.tolist(), key=lambda o: o != self.oid)
        for oid, mol in self.needs_3d(around).items():
            self._ahead = {(name, oid)}
            self._make_3d(name, {oid: mol})
            break

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
        if keys & self._ahead:
            self._ahead = set()
            self._ahead_timer.start()  # the user may have moved on meanwhile
        if self.mode == THREE_D and (name, self.oid) in keys:
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

    # -- drawing --------------------------------------------------------------------------------
    def refresh(self) -> None:
        self._timer.stop()
        if self._scope_stale:
            self._rescope()
        if self.mode == THREE_D:
            self._refresh_3d()
        else:
            self._draw_2d()
        self.caption.set_lines(self.caption_lines())
        self._show_caption_names()
        self._update_counter()
        if not self.gpu_failed:
            self._ahead_timer.start()

    def _refresh_3d(self) -> None:
        molecules = self.molecules()
        parts, _rows = self.scene_parts(molecules)
        scene = combine(parts) if parts else Scene.empty()
        # A new molecule (another row or column, or its 3D structure just made) is turned to
        # its principal axes, centred and fitted; new drawing settings refit only a view the
        # user has not moved. Stepping through conformers keeps the view, to compare them, and
        # a trip to 2D and back keeps it too.
        shown = (self.column_name(), self.oid, bool(parts))
        look = (self.style(), self.hydrogens.currentData(), self.overlay_action.isChecked())
        new = shown != self._shown
        restyled = look != self._look
        self._shown, self._look = shown, look
        if new:
            self._clear_marks()
        self.renderer.set_scene(scene, orient=new, fit=restyled and self.renderer.auto_fitting)
        pending = self._pending_camera
        if pending is not None:  # a reopened document's view, once its molecule is there
            self.renderer.set_camera(pending)
            self._pending_camera = None
        self.overlay.update()

    def _draw_2d(self, *, force: bool = False) -> None:
        name = self.column_name()
        key = (name, self.oid)
        if key == self._drawn_2d and not force:
            return
        self._drawn_2d = key
        if self.oid is None:
            self.depiction.show_message(self._empty_text())
        elif name is None:
            self.depiction.show_message("No structure column")
        elif (mol := self._structure(name, self.oid)) is None:
            self.depiction.show_message("No structure")
        else:
            self.depiction.show_mol(mol)

    def _empty_text(self) -> str:
        if not len(self.table):
            return "The table is empty"
        if self.browse_selected.isChecked():
            return "No rows selected"
        return "No rows shown" if self.browse_visible.isChecked() else "No rows"

    def step(self, delta: int) -> None:
        """Show another conformer (clamped to the first and last)."""
        if self.mode != THREE_D:
            return
        target = max(0, min(self.conformer + delta, max(self._conformers - 1, 0)))
        if target != self.conformer:
            self.conformer = target
            self._clear_marks()
            self.refresh()

    def _update_counter(self) -> None:
        n = len(self._oids)
        blocked = self.position.blockSignals(True)
        self.position.setMaximum(max(n, 1))
        self.position.setValue(self._index + 1 if self._index >= 0 else 1)
        self.position.blockSignals(blocked)
        self.position.setEnabled(n > 0)
        digits = len(f"{max(n, 1):,}")
        self.position.setFixedWidth(self.fontMetrics().horizontalAdvance("0" * digits) + 20)
        self.count_label.setText(f"of {n:,}")
        self.previous_row.setEnabled(self._index > 0)
        self.next_row.setEnabled(0 <= self._index < n - 1)
        self.select_box.setEnabled(self.oid is not None)
        self.select_box.setChecked(self.oid is not None and self.oid in self.selection)
        total = self._conformers if self.mode == THREE_D else 0
        overlay = self.overlay_action.isChecked()
        self.counter.setText(
            f"{total} Conformers" if overlay
            else f"{self.conformer + 1} / {total}"
        )  # fmt: skip
        for action in self._conformer_actions:  # only for a row with several, in 3D
            action.setVisible(total > 1)
        stepping = total > 1 and not overlay
        self.previous.setEnabled(stepping and self.conformer > 0)
        self.next.setEnabled(stepping and self.conformer < total - 1)

    def _schedule(self, *_args) -> None:
        self._timer.start()

    def _rescope_later(self, *_args) -> None:
        self._scope_stale = True
        self._timer.start()

    def _reset(self) -> None:
        if self.renderer is not None:
            self.renderer.reset_view()

    # -- picking ------------------------------------------------------------------------------
    def set_measuring(self, on: bool) -> None:
        self.measuring = bool(on) and self.mode == THREE_D
        if self.measure_action.isChecked() != self.measuring:
            self.measure_action.setChecked(self.measuring)
        self._clear_marks()

    def _clear_marks(self) -> None:
        self._marked = []
        if self.overlay is not None:
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
        return self.overlay.measurement() if self.overlay is not None else ""

    # -- table and selection --------------------------------------------------------------------
    def _on_table_event(self, event: object) -> None:
        if not isinstance(event, Changed):
            return
        if event.structural:
            self._fill_columns()
        self._made = {k: v for k, v in self._made.items() if k[0] in self.table}
        self._drawn_2d = ()  # the structure or the caption's values may have changed
        self._rescope_later()

    def _on_selection(self, _oids: frozenset[int], _origin: str) -> None:
        if self.browse_selected.isChecked():
            self._rescope_later()
        else:
            self.select_box.setChecked(self.oid is not None and self.oid in self.selection)

    def _fill_columns(self) -> None:
        names = [n for n in self.table.column_names if self.table.kind(n) == STRUCTURE]
        current = self.column.currentText()
        blocked = self.column.blockSignals(True)
        self.column.clear()
        self.column.addItems(names)
        if current in names:
            self.column.setCurrentText(current)
        elif self.mode == THREE_D:
            three_d = next((n for n in names if "3D" in n), None)
            self.column.setCurrentText(three_d or (names[0] if names else ""))
        else:
            self.column.setCurrentText(names[0] if names else "")
        self.column.blockSignals(blocked)

    # -- panel ----------------------------------------------------------------------------------
    def sizeHint(self) -> QSize:
        return QSize(520, 480)

    def keyPressEvent(self, event) -> None:
        """Left and right step through the rows, Home and End go to the first and last, Space
        ticks the row; up and down step through a 3D row's conformers (rows in 2D)."""
        rows = {Qt.Key.Key_Left: -1, Qt.Key.Key_Right: 1}
        updown = {Qt.Key.Key_Up: -1, Qt.Key.Key_Down: 1}
        key = event.key()
        if key in rows:
            self.go_by(rows[key])
        elif key == Qt.Key.Key_Home:
            self.go(0)
        elif key == Qt.Key.Key_End:
            self.go(len(self._oids) - 1)
        elif key in updown:
            if self.mode == THREE_D:
                self.step(updown[key])
            else:
                self.go_by(updown[key])
        elif key == Qt.Key.Key_Space and self.oid is not None:
            self.set_selected(self.oid not in self.selection)
        else:
            super().keyPressEvent(event)

    def title(self) -> str:
        return "Structure Viewer"

    def state(self) -> dict:
        camera = self._camera_state if self.renderer is None else self.renderer.camera.state()
        return {
            "mode": self.mode,
            "browse_selected": self.browse_selected.isChecked(),
            "browse_visible": self.browse_visible.isChecked(),
            "row": self.oid,
            "column": self.column.currentText(),
            "caption": list(self.caption_names),
            "style": self.style(),
            "hydrogens": self.hydrogens.currentData(),
            "colour": self.colour.currentData(),
            "align": self.align.isChecked(),
            "overlay": self.overlay_action.isChecked(),
            "conformer": self.conformer,
            "camera": camera,
        }

    def set_state(self, state: dict) -> None:
        """Restore :meth:`state`. A 3D Viewer panel saved before the 2D view existed has no
        mode, and reopens in 3D."""
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
        if isinstance(state.get("caption"), list):
            self.set_caption(state["caption"])
        self._rescope()
        row = state.get("row")
        hit = np.flatnonzero(self._oids == row) if isinstance(row, int) else []
        if len(hit):
            self._set_index(int(hit[0]))
        self.overlay_action.setChecked(bool(state.get("overlay", False)))
        self.conformer = max(0, int(state.get("conformer", 0) or 0))
        camera = state.get("camera") if isinstance(state.get("camera"), dict) else None
        if self.renderer is not None and camera is not None:
            self._pending_camera = Camera.from_state(camera)
        else:
            self._camera_state = camera
        mode = state.get("mode", THREE_D)
        self.set_mode(mode if mode in MODES else TWO_D)
        if self._move_to is not None and self.oid is not None:
            self._move_to(self.oid)

    def export_dialog(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Export View", "structure.png", "PNG (*.png)")
        if path:
            self.export(Path(path))

    def picture(self) -> QImage:
        """The view as shown, caption included."""
        if self.mode == THREE_D and self.renderer is not None:
            image = self.renderer.grabFramebuffer()
            image.setDevicePixelRatio(self.renderer.devicePixelRatioF())
        else:
            size = self.depiction.size()
            dpr = self.devicePixelRatioF()
            image = QImage(
                int(size.width() * dpr), int(size.height() * dpr), QImage.Format.Format_ARGB32
            )
            image.setDevicePixelRatio(dpr)
            image.fill(self.palette().base().color())
            mol = self.depiction.mol
            if mol is not None:
                dark = self.palette().base().color().lightness() < 128
                drawing = structure_image(mol, size.width(), size.height(), dpr, dark)
                painter = QPainter(image)
                painter.drawImage(0, 0, drawing)
                painter.end()
        if not image.isNull() and self.caption.lines:
            painter = QPainter(image)
            width = image.width() / image.devicePixelRatio()
            paint_caption(painter, self.caption.lines, width, self.caption.font(), self.palette())
            painter.end()
        return image

    def export(self, path: Path) -> bool:
        image = self.picture()
        return not image.isNull() and image.save(str(path.with_suffix(".png")))

    def attach_header(self, attached: bool) -> None:
        self._header_attached = attached
        self._restyle()

    def close_panel(self) -> None:
        self._unsubscribe_table()
        self._unsubscribe_selection()
        self._stop_watching()
        self._timer.stop()
        self._ahead_timer.stop()
        self.options.hide()

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and hasattr(self, "toolbar"):
            self._restyle()
            self._draw_2d(force=True)  # drawings are coloured for a light or dark background

    def _restyle(self) -> None:
        palette = self.palette()
        ink = palette.text().color()
        for action, name in self._glyphs.items():
            action.setIcon(glyph(name, ink))
        base = palette.base().color()
        if self.renderer is not None:
            self.renderer.set_background(base)
        line = palette.mid().color().name()
        faint = palette.placeholderText().color()
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
            f"QFrame#ViewerFrame {{ border: none; background: {base.name()}; }}"
            f" QLabel#ViewerDepiction {{ background: {base.name()};"
            f" color: {rgba(*faint.getRgb()[:3], faint.alphaF())}; }}"
        )
        buttons = (
            " QToolButton { border: none; border-radius: 4px; padding: 2px 4px; }"
            f" QToolButton:hover {{ background: {rgba(*hover.getRgb()[:3], hover.alphaF())}; }}"
            " QToolButton:checked {"
            f" background: {rgba(*checked.getRgb()[:3], checked.alphaF())}; }}"
        )
        self.toolbar.setStyleSheet(f"QToolBar {{ {bar} }}" + buttons)
        self.rowbar.setStyleSheet(
            f"QToolBar {{ background: {base.name()}; border: none;"
            f" border-top: 1px solid {line}; padding: 2px 4px; spacing: 2px; }}" + buttons
        )
