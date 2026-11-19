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
"""The Matched Pairs panel: transformations drawn, the pairs behind the chosen ones (laid out
alike, the change highlighted), and how their activity shifts spread."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from functools import lru_cache

import numpy as np
from PySide6.QtCore import QItemSelectionModel, QPointF, Qt, Signal
from PySide6.QtGui import QColor, QMouseEvent, QPainter, QPaintEvent, QPen
from PySide6.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem, rdBase

from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.result_table import Drawing
from chemsuite.core.depict import render_png
from chemsuite.sar.mmp import Pair, Transform, pair_depiction, readable
from chemsuite.tools.similarity.common import labels

from .common import ResultTable, cell_mol
from .matrix_panel import HIGH, LOW

NO_ACTIVITY = "(none)"
TRANSFORM_ROW = 56  # px: room for a substituent drawing
PAIR_ROW = 130  # px: a compound drawing and its name
FRAGMENT_WIDTH = 110
COMPOUND_WIDTH = 200


def panel_title(activity: str) -> str:
    return "Matched Pairs" + (f": {activity}" if activity != NO_ACTIVITY else "")


def transform_rows(found: list[Transform]) -> list[list]:
    """The transformations table's cells, one list per transformation."""
    return [
        [readable(t.from_group), readable(t.to_group), t.count, t.mean, t.sd, t.ups, t.downs]
        for t in found
    ]


@lru_cache(maxsize=2048)
def _fragment_mol(group: str) -> Chem.Mol | None:
    """A changing group as drawn: its attachment point a wavy line, map numbers left out."""
    params = Chem.SmilesParserParams()
    params.removeHs = False  # keep the H of the hydrogen group
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(group, params)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def fragment_drawing(group: str) -> Drawing:
    def draw(width: int, height: int, scale: float, dark: bool) -> bytes | None:
        mol = _fragment_mol(group)
        if mol is None:
            return None
        return render_png(mol, width, height, scale=scale, dark=dark, attachments=True)

    return Drawing(readable(group), draw)


class DeltaStrip(QWidget):
    """One dot per pair on a Δ axis, with zero and the mean marked: how consistent a
    transformation's activity shift is."""

    #: a dot was clicked: the index of its pair in :meth:`set_deltas`'s list
    picked = Signal(int)

    MARGIN = 14
    RADIUS = 4.0

    def __init__(self) -> None:
        super().__init__()
        self.setMinimumHeight(48)
        self.setMaximumHeight(48)
        self.deltas = np.zeros(0)

    def set_deltas(self, deltas: Sequence[float]) -> None:
        self.deltas = np.asarray(deltas, dtype=np.float64)
        finite = self.deltas[np.isfinite(self.deltas)]
        self.setVisible(len(finite) > 0)
        if len(finite):
            sd = f", SD {finite.std(ddof=1):.2f}" if len(finite) > 1 else ""
            pairs = f"{len(finite)} pair{'s' if len(finite) != 1 else ''}"
            self.setToolTip(f"Δ over {pairs}: mean {finite.mean():+.2f}{sd}")
        self.update()

    def limit(self) -> float:
        finite = np.abs(self.deltas[np.isfinite(self.deltas)])
        top = float(finite.max()) if len(finite) else 0.0
        return top * 1.1 if top > 0 else 1.0

    def x_of(self, value: float) -> float:
        left, right = self.MARGIN, self.width() - self.MARGIN
        return left + (value + self.limit()) / (2 * self.limit()) * (right - left)

    def points(self) -> list[tuple[int, QPointF]]:
        """``(pair index, center)`` per finite Δ, nudged up and down so equal shifts show."""
        middle = (self.height() - 14) / 2
        out = []
        for i, value in enumerate(self.deltas.tolist()):
            if math.isfinite(value):
                offset = ((i * 7) % 5 - 2) * self.RADIUS
                out.append((i, QPointF(self.x_of(value), middle + offset)))
        return out

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        palette = self.palette()
        text = palette.text().color()
        axis = self.height() - 14
        faint = QColor(text)
        faint.setAlpha(90)
        painter.setPen(QPen(faint, 1))
        painter.drawLine(QPointF(self.MARGIN, axis), QPointF(self.width() - self.MARGIN, axis))
        zero = self.x_of(0.0)
        painter.drawLine(QPointF(zero, 2), QPointF(zero, axis))
        finite = self.deltas[np.isfinite(self.deltas)]
        if len(finite):
            mean = self.x_of(float(finite.mean()))
            painter.setPen(QPen(text, 2))
            painter.drawLine(QPointF(mean, 2), QPointF(mean, axis))
        painter.setPen(Qt.PenStyle.NoPen)
        for i, center in self.points():
            value = self.deltas[i]
            painter.setBrush(HIGH if value > 0 else LOW if value < 0 else faint)
            painter.drawEllipse(center, self.RADIUS, self.RADIUS)
        painter.setPen(text)
        limit = self.limit()
        font = painter.font()
        font.setPointSizeF(font.pointSizeF() * 0.85)
        painter.setFont(font)
        box = self.rect().adjusted(self.MARGIN, axis + 1, -self.MARGIN, 0)
        painter.drawText(box, Qt.AlignmentFlag.AlignLeft, f"{-limit:+.2f}")
        painter.drawText(box, Qt.AlignmentFlag.AlignRight, f"{limit:+.2f}")
        painter.drawText(box, Qt.AlignmentFlag.AlignHCenter, "0")
        painter.end()

    def mousePressEvent(self, event: QMouseEvent) -> None:
        where = event.position()
        near = [
            (abs(c.x() - where.x()) + abs(c.y() - where.y()), i)
            for i, c in self.points()
            if abs(c.x() - where.x()) <= 2 * self.RADIUS
        ]
        if near:
            self.picked.emit(min(near)[1])


class MmpPanel(ToolPanel):
    """Transformations, and the pairs behind the chosen ones."""

    tool_id = "mmp"

    def __init__(
        self,
        ctx: ToolContext,
        oids,
        activity,
        name: str,
        found: list[Transform],
        recipe=None,
        rows: list[list] | None = None,
        cells: Sequence | None = None,
    ):
        """*rows* are :func:`transform_rows` of *found*, made in the job (naming the groups
        is slow for thousands of transformations). *cells* are the structures the pairs were
        found in, by input row (none: the pairs are listed without drawings)."""
        super().__init__()
        self.ctx = ctx
        self.recipe = dict(recipe or {})
        self.oids = np.asarray(oids)
        self.activity = np.asarray(activity)
        self.name = name
        self.found = found
        # Only the structures of listed pairs are kept, for drawing them.
        self.cells = {r: cells[r] for t in found for r in t.rows()} if cells is not None else None
        self.mols: dict[int, Chem.Mol | None] = {}
        self.names = labels(ctx, self.oids)
        self.index = {int(oid): i for i, oid in enumerate(self.oids.tolist())}
        self.by_row: dict[int, list[int]] | None = None  # input row → its transformations
        self.shown: list[int] = []  # the transformations listed (indices into found)
        texts = rows if rows is not None else transform_rows(found)
        fragments: dict[str, Drawing] = {}
        self.rows = [
            [
                fragments.setdefault(t.from_group, fragment_drawing(t.from_group)),
                fragments.setdefault(t.to_group, fragment_drawing(t.to_group)),
                *text[2:],
            ]
            for t, text in zip(found, texts, strict=True)
        ]
        self.transforms = ResultTable(
            ["From", "To", "Pairs", "Mean Δ", "SD", "Up", "Down"],
            row_height=TRANSFORM_ROW,
            drawing_width=FRAGMENT_WIDTH,
        )
        self.transforms.setToolTip(
            f"Δ = {name} of the “to” compound minus the “from” compound."
            if name != NO_ACTIVITY
            else "No activity column was chosen."
        )
        self.strip = DeltaStrip()
        self.strip.hide()
        self.strip.picked.connect(self.pick_pair)
        self.shown_pairs: list[Pair] = []
        self.pairs = ResultTable(
            ["From compound", "To compound", "From", "To", "Δ"],
            row_height=PAIR_ROW if self.cells is not None else 22,
            drawing_width=COMPOUND_WIDTH,
        )
        self.pairs.setToolTip("Double-click a pair to go to its compound in the table.")
        self.pairs.cellDoubleClicked.connect(self.go_to)
        self.transforms.itemSelectionChanged.connect(self.show_pairs)
        select = QPushButton("Select Compounds")
        select.setToolTip("Select the compounds in the chosen transformations or pairs.")
        select.clicked.connect(lambda: self.select_rows())
        self.current_only = QPushButton("Current Compound")
        self.current_only.setCheckable(True)
        self.current_only.setToolTip(
            "List only the transformations of the table's current compound."
        )
        self.current_only.setChecked(bool(self.recipe.get("current_only", False)))
        self.current_only.toggled.connect(lambda *_: self.refilter())
        copy = QPushButton("Copy Transformations")
        copy.clicked.connect(lambda: self.transforms.copy_text())
        below = QWidget()
        below_layout = QVBoxLayout(below)
        below_layout.setContentsMargins(0, 0, 0, 0)
        below_layout.addWidget(self.strip)
        below_layout.addWidget(self.pairs, 1)
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.transforms)
        split.addWidget(below)
        buttons = QHBoxLayout()
        buttons.addWidget(select)
        buttons.addWidget(self.current_only)
        buttons.addStretch(1)
        buttons.addWidget(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(split, 1)
        layout.addLayout(buttons)
        self.stop_watching: Callable[[], None] = ctx.watch_current(lambda _oid: self.on_current())
        self.refilter()
        self.transforms.sortByColumn(2, Qt.SortOrder.DescendingOrder)

    def title(self) -> str:
        return panel_title(self.name)

    # -- the transformations listed -----------------------------------------------------------
    def listed(self) -> list[int]:
        """The transformations shown: all, or those of the current compound."""
        if not self.current_only.isChecked():
            return list(range(len(self.found)))
        current = self.ctx.current_oid()
        row = self.index.get(current) if current is not None else None
        if row is None:
            return []
        if self.by_row is None:
            self.by_row = {}
            for i, t in enumerate(self.found):
                for r in t.rows():
                    self.by_row.setdefault(r, []).append(i)
        return self.by_row.get(row, [])

    def on_current(self) -> None:
        # Moving along a row, or to a compound with the same transformations, changes nothing.
        if self.current_only.isChecked() and self.listed() != self.shown:
            self.refilter()

    def chosen(self) -> set[int]:
        return {id(t) for t in self.transforms.selected_data() if t is not None}

    def refilter(self) -> None:
        """Fill the transformations table, keeping its sort and the chosen transformations
        that are still listed. The pairs are refilled only when the choice changes, so going
        to a pair's compound keeps the pair in view."""
        chosen = self.chosen()
        header = self.transforms.horizontalHeader()
        column, order = header.sortIndicatorSection(), header.sortIndicatorOrder()
        shown = self.listed()
        self.shown = shown
        self.transforms.blockSignals(True)
        self.transforms.fill([self.rows[i] for i in shown], [self.found[i] for i in shown])
        if column >= 0:
            self.transforms.sortByColumn(column, order)
        selection = self.transforms.selectionModel()
        flags = QItemSelectionModel.SelectionFlag.Select | QItemSelectionModel.SelectionFlag.Rows
        for r in range(self.transforms.rowCount()):
            if id(self.transforms.row_data(r)) in chosen:
                selection.select(self.transforms.result_model.index(r, 0), flags)
        self.transforms.blockSignals(False)
        if self.chosen() != chosen:
            self.show_pairs()

    # -- the pairs behind them ----------------------------------------------------------------
    def compound_mol(self, row: int) -> Chem.Mol | None:
        if row not in self.mols:
            self.mols[row] = cell_mol(self.cells[row]) if self.cells is not None else None
        return self.mols[row]

    def compound(self, row: int, constant: str, group: str) -> Drawing | str:
        name = self.names[row]
        if self.cells is None:
            return name

        def draw(width: int, height: int, scale: float, dark: bool) -> bytes | None:
            mol = self.compound_mol(row)
            if mol is None:
                return None
            shown, atoms, bonds = pair_depiction(mol, constant, group)
            return render_png(
                shown,
                width,
                height,
                scale=scale,
                dark=dark,
                highlight_atoms=atoms,
                highlight_bonds=bonds,
            )

        return Drawing(name, draw, caption=name)

    def show_pairs(self) -> None:
        chosen = [t for t in self.transforms.selected_data() if t is not None]
        pairs = [p for t in chosen for p in t.pairs]
        self.shown_pairs = pairs
        heads = ["From compound", "To compound", f"From {self.name}", f"To {self.name}", "Δ"]
        self.pairs.setHorizontalHeaderLabels(heads)  # first: the columns are sized to fit them
        self.pairs.fill(
            [
                [
                    self.compound(p.a, p.constant, p.from_group),
                    self.compound(p.b, p.constant, p.to_group),
                    self.activity[p.a],
                    self.activity[p.b],
                    p.delta,
                ]
                for p in pairs
            ],
            pairs,
        )
        self.strip.set_deltas([p.delta for p in pairs] if self.name != NO_ACTIVITY else [])

    def pick_pair(self, index: int) -> None:
        """Choose the pair of a dot in the Δ strip."""
        pair = self.shown_pairs[index]
        for r in range(self.pairs.rowCount()):
            if self.pairs.row_data(r) is pair:
                self.pairs.selectRow(r)
                self.pairs.scrollTo(self.pairs.result_model.index(r, 0))
                return

    def go_to(self, row: int, column: int) -> int | None:
        """Make the double-clicked pair's compound the table's current row: the "from"
        compound from its own columns, otherwise the "to" compound."""
        pair = self.pairs.row_data(row)
        if not isinstance(pair, Pair):
            return None
        oid = int(self.oids[pair.a if column in (0, 2) else pair.b])
        return oid if self.ctx.set_current_oid(oid) else None

    def select_rows(self) -> list[int]:
        pairs = [p for p in self.pairs.selected_data() if p is not None]
        if pairs:
            rows = {r for p in pairs for r in (p.a, p.b)}
        else:
            chosen = [t for t in self.transforms.selected_data() if t is not None]
            rows = set().union(*(t.rows() for t in chosen)) if chosen else set()
        oids = self.oids[sorted(rows)].tolist()
        if oids:
            self.ctx.selection.set(oids, "mmp")
        return oids

    def state(self) -> dict:
        return {**self.recipe, "current_only": self.current_only.isChecked()}

    def close_panel(self) -> None:
        self.stop_watching()
