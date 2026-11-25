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
"""Data ▸ SAR ▸ Scaffold Tree: each compound's scaffold pruned ring by ring, as a tree."""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
from PySide6.QtCore import (
    QModelIndex,
    QSize,
    QSortFilterProxyModel,
    Qt,
)
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QSpinBox,
    QTreeView,
    QVBoxLayout,
)
from rdkit import Chem

from chemsuite.analysis.selectivity import looks_like_assay
from chemsuite.app.document.tools import ToolContext, ToolPanel
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.app.widgets.result_table import (
    DRAWING_ROLE,
    Drawing,
    DrawingDelegate,
    TreeTableModel,
)
from chemsuite.core.depict import render_png
from chemsuite.core.wording import plural
from chemsuite.plotting.data import number_columns, numbers
from chemsuite.sar.matrix import AGGREGATES
from chemsuite.sar.scaffold_tree import MEANINGS, ScaffoldTree, build

from .common import RerunPanel, activity_combo, saved_rows, structure_cells
from .matrix_panel import heat

TITLE = "Scaffold Tree"
NOTHING = "None"  # Show: no summary column
POOL_MIN_ROWS = 500  # pruning is slow per scaffold: worth the worker start-up early
PICTURE = QSize(180, 110)
SORT_ROLE = Qt.ItemDataRole.UserRole + 2
SCAFFOLD, COMPOUNDS, RINGS, SUMMARY = range(4)
STATISTICS: dict[str, Callable[[np.ndarray], float]] = {
    "mean": lambda v: float(np.mean(v)),
    "median": lambda v: float(np.median(v)),
    "max": lambda v: float(np.max(v)),
    "min": lambda v: float(np.min(v)),
}
INK = QColor(0, 0, 0)  # text on a tinted cell


def _drawing(smiles: str, caption: str) -> Drawing:
    def draw(width: int, height: int, scale: float, dark: bool) -> bytes | None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return render_png(mol, width, height, scale=scale, dark=dark)

    return Drawing(smiles, draw, caption)


def summaries(tree: ScaffoldTree, values: np.ndarray, statistic: str) -> np.ndarray:
    """*statistic* of the finite *values* (per input row) over each node's subtree (NaN for
    a node with none)."""
    reduce = STATISTICS[statistic]
    out = np.full(len(tree), np.nan)
    for at, rows in enumerate(tree.members()):
        v = values[rows]
        v = v[np.isfinite(v)]
        if len(v):
            out[at] = reduce(v)
    return out


class TreeModel(TreeTableModel):
    """The scaffold tree's nodes, children by subtree count (largest first)."""

    def __init__(self, tree: ScaffoldTree) -> None:
        counts = tree.subtree_counts()
        order = sorted(range(len(tree)), key=lambda i: (-int(counts[i]), tree.smiles[i]))
        kids: dict[int, list[int]] = {-1: []}
        for i in order:
            kids.setdefault(int(tree.parent[i]), []).append(i)
        super().__init__(kids, tree.parent.tolist(), 4)
        self.tree = tree
        self.counts = counts
        self.own = tree.own_counts()
        self.drawings: dict[int, Drawing] = {}
        self.values = np.full(len(tree), np.nan)
        self.low = self.high = math.nan
        self.higher_is_better = True
        self.summary_name = ""

    def set_summary(self, values: np.ndarray, name: str, higher_is_better: bool) -> None:
        self.values, self.summary_name = values, name
        self.higher_is_better = higher_is_better
        finite = values[np.isfinite(values)]
        self.low, self.high = (
            (float(finite.min()), float(finite.max())) if len(finite) else (math.nan,) * 2
        )
        self.refresh()

    def header(self, column: int) -> str:
        return ("Scaffold", "Compounds", "Rings", self.summary_name)[column]

    def cell(self, at: int, column: int, role: int) -> object:
        tree = self.tree
        if role == Qt.ItemDataRole.ToolTipRole:
            return self.tooltip(at)
        if column == SCAFFOLD:
            if role == DRAWING_ROLE:
                drawing = self.drawings.get(at)
                if drawing is None:  # one per node: the delegate's cache knows it again
                    drawing = self.drawings[at] = _drawing(tree.smiles[at], "")
                return drawing
            if role in (Qt.ItemDataRole.DisplayRole, SORT_ROLE):
                return tree.smiles[at]
            if role == Qt.ItemDataRole.SizeHintRole:
                return PICTURE
            return None
        number = {
            COMPOUNDS: float(self.counts[at]),
            RINGS: float(tree.rings[at]),
            SUMMARY: float(self.values[at]),
        }[column]
        if role == SORT_ROLE:
            return number if math.isfinite(number) else -math.inf
        if role == Qt.ItemDataRole.DisplayRole:
            if not math.isfinite(number):
                return ""
            return f"{number:.2f}" if column == SUMMARY else f"{int(number):,}"
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if column == SUMMARY and math.isfinite(number):
            if role == Qt.ItemDataRole.BackgroundRole:
                return heat(number, self.low, self.high, self.higher_is_better)
            if role == Qt.ItemDataRole.ForegroundRole:
                return INK
        return None

    def tooltip(self, at: int) -> str:
        tree = self.tree
        lines = [
            tree.smiles[at],
            f"{plural(int(self.counts[at]), 'compound')} in this branch, "
            f"{int(self.own[at]):,} with this scaffold",
        ]
        why = tree.rule[at]
        if why:
            rule = "" if why == "only" else f" by rule {why}"
            lines.append(f"Pruned to its parent{rule}: {MEANINGS[why]}")
        return "\n".join(lines)


class CountFilter(QSortFilterProxyModel):
    """Hides branches with fewer compounds than *least*; sorts by :data:`SORT_ROLE`."""

    def __init__(self, model: TreeModel) -> None:
        super().__init__()
        self.least = 1
        self.tree_model = model
        self.setSourceModel(model)
        self.setSortRole(SORT_ROLE)
        self.setRecursiveFilteringEnabled(True)
        self.setDynamicSortFilter(False)

    def set_least(self, least: int) -> None:
        self.beginFilterChange()
        self.least = least
        self.endFilterChange(QSortFilterProxyModel.Direction.Rows)

    def filterAcceptsRow(self, row, parent) -> bool:
        model = self.tree_model
        return int(model.counts[model.node(model.index(row, 0, parent))]) >= self.least


class ScaffoldTreePanel(ToolPanel):
    tool_id = "scaffold_tree"

    def __init__(self, ctx: ToolContext, oids, tree: ScaffoldTree, recipe: dict) -> None:
        super().__init__()
        self.ctx = ctx
        self.oids = np.asarray(oids, dtype=np.int64)
        self.tree = tree
        self.recipe = dict(recipe)
        self.model = TreeModel(tree)
        self.proxy = CountFilter(self.model)
        self.view = QTreeView()
        self.view.setModel(self.proxy)
        self.view.setItemDelegateForColumn(SCAFFOLD, DrawingDelegate(PICTURE.width(), self.view))
        self.view.setUniformRowHeights(True)
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(COMPOUNDS, Qt.SortOrder.DescendingOrder)
        self.view.setSelectionMode(QTreeView.SelectionMode.SingleSelection)
        header = self.view.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(SCAFFOLD, QHeaderView.ResizeMode.Interactive)
        header.resizeSection(SCAFFOLD, PICTURE.width() + 60)
        self.view.clicked.connect(lambda index: self.select_branch(index))

        self.show_column = QComboBox()
        self.show_column.addItem(NOTHING)
        self.show_column.addItems(number_columns(ctx.table))
        self.show_column.setToolTip("A number column summarized over each branch")
        self.statistic = QComboBox()
        for key, label in AGGREGATES.items():
            self.statistic.addItem(label, key)
        self.statistic.setToolTip("How each branch's values are summarized")
        self.higher = QCheckBox("Higher is better")
        self.higher.setToolTip("Green for high values, red for low (unchecked: the reverse)")
        self.least = QSpinBox()
        self.least.setRange(1, 1_000_000)
        self.least.setPrefix("Min. compounds ")
        self.least.setToolTip("Hide branches with fewer compounds")
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Show"))
        for widget in (self.show_column, self.statistic, self.higher):
            controls.addWidget(widget)
        controls.addStretch(1)
        controls.addWidget(self.least)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(controls)
        layout.addWidget(self.view, 1)

        shown = recipe.get("show", NOTHING)
        self.show_column.setCurrentText(shown if self.show_column.findText(shown) >= 0 else NOTHING)
        at = self.statistic.findData(recipe.get("statistic", "mean"))
        self.statistic.setCurrentIndex(max(at, 0))
        self.higher.setChecked(bool(recipe.get("higher_is_better", True)))
        self.least.setValue(int(recipe.get("min_compounds", 1)))
        self.update_summary()
        self.proxy.set_least(self.least.value())
        self.expand(recipe.get("expanded", []))
        for signal in (
            self.show_column.currentTextChanged,
            self.statistic.currentIndexChanged,
            self.higher.toggled,
        ):
            signal.connect(lambda *_: self.update_summary())
        self.least.valueChanged.connect(lambda value: self.proxy.set_least(value))

    def title(self) -> str:
        return TITLE

    def update_summary(self) -> None:
        name = self.show_column.currentText()
        hidden = name == NOTHING or name not in self.ctx.table
        self.statistic.setEnabled(not hidden)
        self.higher.setEnabled(not hidden)
        if hidden:
            values = np.full(len(self.tree), np.nan)
            label = ""
        else:
            statistic = self.statistic.currentData()
            row_values = numbers(self.ctx.table, name, self.oids)
            values = summaries(self.tree, row_values, statistic)
            label = f"{AGGREGATES[statistic]} {name}"
        self.model.set_summary(values, label, self.higher.isChecked())
        self.view.setColumnHidden(SUMMARY, hidden)

    def node_at(self, index: QModelIndex) -> int:
        return self.model.node(self.proxy.mapToSource(index)) if index.isValid() else -1

    def select_branch(self, index: QModelIndex) -> list[int]:
        """Select the rows of the clicked node's branch in the table."""
        at = self.node_at(index)
        if at < 0:
            return []
        oids = self.oids[self.tree.subtree_rows(at)].tolist()
        if oids:
            self.ctx.selection.set(oids, "scaffold_tree")
        return oids

    def expanded(self) -> list[str]:
        """The SMILES of the expanded nodes."""
        out = []
        for at in range(len(self.tree)):
            if self.model.kids.get(at):
                index = self.proxy.mapFromSource(self.model.node_index(at))
                if index.isValid() and self.view.isExpanded(index):
                    out.append(self.tree.smiles[at])
        return out

    def expand(self, smiles: list[str]) -> None:
        where = {s: i for i, s in enumerate(self.tree.smiles)}
        for s in smiles:
            at = where.get(s)
            if at is not None:
                index = self.proxy.mapFromSource(self.model.node_index(at))
                if index.isValid():
                    self.view.expand(index)

    def state(self) -> dict:
        return {
            **self.recipe,
            "show": self.show_column.currentText(),
            "statistic": self.statistic.currentData(),
            "higher_is_better": self.higher.isChecked(),
            "min_compounds": self.least.value(),
            "expanded": self.expanded(),
        }

    def close_panel(self) -> None:
        pass


class ScaffoldTreeDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, TITLE, guide="scaffold-tree")
        self.ctx = ctx
        self.panel: ScaffoldTreePanel | None = None
        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx, "shown"))
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    def problem(self) -> str:
        return "" if self.structure.currentText() else "There is no structure column."

    def verb(self) -> str:
        return f"Build the Tree of {len(self.rows.oids()):,} Rows"

    def run(self):
        return self.build()

    def recipe(self) -> dict:
        guess = activity_combo(self.ctx).currentText()
        return {
            "oids": self.rows.oids().tolist(),
            "structure": self.structure.currentText(),
            "show": guess if guess and looks_like_assay(guess) else NOTHING,
        }

    def build(self):
        recipe = self.recipe()
        return _compute(self.ctx, recipe, lambda oids, tree: self.finish(oids, tree, recipe))

    def finish(self, oids, tree: ScaffoldTree, recipe: dict) -> ScaffoldTreePanel:
        ctx = self.ctx
        self.panel = ScaffoldTreePanel(ctx, oids, tree, recipe)
        ctx.add_panel(self.panel, self.panel.title())
        placed = int((tree.node >= 0).sum())
        ctx.show_message(
            f"{plural(len(tree), 'scaffold')} in {plural(len(tree.roots()), 'tree')} "
            f"from {plural(placed, 'compound')} with a ring"
        )
        return self.panel


def _compute(ctx: ToolContext, recipe: dict, done):
    """Build the tree *recipe* describes in the background; then ``done(oids, tree)``."""
    oids = np.asarray(recipe["oids"], dtype=np.int64)
    cells = structure_cells(ctx, recipe["structure"], oids)
    executor = ctx.workers_for(len(cells), POOL_MIN_ROWS)

    def work(job):
        return build(cells, executor=executor, progress=job.report, cancel=job.cancel_token)

    return ctx.run("Building the scaffold tree…", work, lambda tree: done(oids, tree))


def restore(ctx: ToolContext, state: dict) -> RerunPanel:
    """Reopen a saved panel: rebuild the tree of the rows it covered."""
    recipe = {**state, "oids": saved_rows(ctx, state, state["structure"])}
    holder = RerunPanel(ScaffoldTreePanel.tool_id, TITLE, recipe)
    ctx.add_panel(holder, holder.title())

    def done(oids, tree) -> None:
        holder.fill(ScaffoldTreePanel(ctx, oids, tree, recipe))

    holder.job = _compute(ctx, recipe, done)
    return holder
