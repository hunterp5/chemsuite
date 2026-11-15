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
"""A protein's docking runs, their poses grouped by ligand, beside the pocket.

The panel shows one protein's runs (:class:`~chemsuite.protein.workspace.poses.DockedPoses`,
kept with the document, so they come back with it), one at a time: the Run list chooses which,
newest first. Each ligand is a row with its best affinity, and each pose under it; selecting a
pose shows it in the pocket, and the panel is told who to ask for that (*show*), so it knows
nothing about the table, the tools or RDKit.

**Copy to Table** hands the chosen poses to whoever opened the panel (*to_table*): a pose is
only a table row when the chemist asks for one. **Remove Run** asks it to forget the run shown
(*remove*), which answers by filling the panel again.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
)

from chemsuite.app.widgets.panel import ToolPanel
from chemsuite.protein.workspace.poses import DockedPoses, Pose

#: the tool that reopens this panel with the document
TOOL_ID = "protein-poses"
HEADERS = ("Ligand", "Pose", "Affinity", "RMSD l.b.", "RMSD u.b.")
#: a pose is held in its row, under this role
POSE = Qt.ItemDataRole.UserRole


class PosesPanel(ToolPanel):
    """One protein's runs; of the one shown, a row per ligand, its poses under it."""

    tool_id = TOOL_ID

    def __init__(
        self,
        runs: Sequence[DockedPoses],
        show: Callable[[Pose], None] | None = None,
        to_table: Callable[[Sequence[Pose]], None] | None = None,
        remove: Callable[[DockedPoses], None] | None = None,
    ) -> None:
        super().__init__()
        self.runs: tuple[DockedPoses, ...] = ()
        self.run = runs[-1]
        self._show = show
        self._to_table = to_table
        self._remove = remove
        self.choice = QComboBox()
        self.choice.setToolTip("Every docking and rescoring run into this protein, newest first.")
        self.choice.currentIndexChanged.connect(self._run_chosen)
        self.remove = QPushButton("Remove Run")
        self.remove.setToolTip("Forget this run and its poses (rows copied to the table stay).")
        self.remove.clicked.connect(self.remove_run)
        top = QHBoxLayout()
        top.addWidget(self.choice, 1)
        top.addWidget(self.remove)
        self.rows = QTreeWidget()
        self.rows.setHeaderLabels(list(HEADERS))
        self.rows.setUniformRowHeights(True)
        self.rows.setAlternatingRowColors(True)
        self.rows.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        header = self.rows.header()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for column in range(1, len(HEADERS)):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        self.rows.itemSelectionChanged.connect(self._selected)

        self.rows.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)
        copy = QAction("Copy Scores", self.rows)
        copy.setShortcut(QKeySequence.StandardKey.Copy)
        copy.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
        copy.triggered.connect(self.copy_text)
        self.rows.addAction(copy)

        self.to_table = QPushButton("Copy to Table")
        self.to_table.setEnabled(False)
        self.to_table.clicked.connect(self.send_to_table)
        footer = QHBoxLayout()
        footer.addStretch(1)
        footer.addWidget(self.to_table)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addLayout(top)
        layout.addWidget(self.rows, 1)
        layout.addLayout(footer)
        self.fill(runs)

    # -- contents -----------------------------------------------------------------------------
    def fill(self, runs: Sequence[DockedPoses], number: int | None = None) -> None:
        """Show *runs* (one protein's, oldest first), the one numbered *number* (None: the
        newest) on screen."""
        if not runs:
            return
        self.runs = tuple(runs)
        chosen = next((run for run in runs if run.number == number), runs[-1])
        self.choice.blockSignals(True)
        self.choice.clear()
        for run in reversed(runs):
            self.choice.addItem(f"{run.label} · {len(run.poses)} poses", run.number)
        self.choice.setCurrentIndex(self.choice.findData(chosen.number))
        self.choice.blockSignals(False)
        self.show_run(chosen)

    def _run_chosen(self) -> None:
        number = self.choice.currentData()
        run = next((run for run in self.runs if run.number == number), None)
        if run is not None:
            self.show_run(run)

    def show_run(self, run: DockedPoses) -> None:
        """List *run*'s poses."""
        self.run = run
        self.rows.clear()
        for ligand in run.ligands:
            poses = [pose for pose in run.poses if pose.ligand == ligand]
            best = min(pose.affinity for pose in poses)
            parent = QTreeWidgetItem([ligand, "", f"{best:.2f}", "", ""])
            for pose in poses:
                parent.addChild(
                    QTreeWidgetItem(
                        [
                            "",
                            str(pose.rank),
                            f"{pose.affinity:.2f}",
                            f"{pose.rmsd_lb:.2f}",
                            f"{pose.rmsd_ub:.2f}",
                        ]  # fmt: skip
                    )
                )
                parent.child(parent.childCount() - 1).setData(0, POSE, pose)
            self.rows.addTopLevelItem(parent)
            parent.setExpanded(len(run.ligands) == 1)
        self.titleChanged.emit(self.title())

    def title(self) -> str:
        run = self.run
        many = f" · {len(self.runs)} runs" if len(self.runs) > 1 else ""
        return f"Poses: {run.protein}{many}"

    def chosen(self) -> list[Pose]:
        """The selected poses (a ligand's row stands for all of its)."""
        found: list[Pose] = []
        for item in self.rows.selectedItems():
            pose = item.data(0, POSE)
            if isinstance(pose, Pose):
                found.append(pose)
            else:
                found += [item.child(i).data(0, POSE) for i in range(item.childCount())]
        return [pose for pose in dict.fromkeys(found) if isinstance(pose, Pose)]

    def copy_text(self) -> str:
        """The scores as lines of text, on the clipboard."""
        lines = ["\t".join(HEADERS)]
        lines += [
            f"{pose.ligand}\t{pose.rank}\t{pose.affinity:.2f}"
            f"\t{pose.rmsd_lb:.2f}\t{pose.rmsd_ub:.2f}"
            for pose in self.run.poses
        ]
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text

    # -- the workbench ------------------------------------------------------------------------
    def _selected(self) -> None:
        chosen = self.chosen()
        self.to_table.setEnabled(bool(chosen) and self._to_table is not None)
        if chosen and self._show is not None:
            self._show(chosen[0])  # show the first: one pose sits in the pocket at a time

    def remove_run(self) -> DockedPoses | None:
        """Ask for the run shown to be forgotten (the button)."""
        if self._remove is None:
            return None
        run = self.run
        self._remove(run)
        return run

    def send_to_table(self) -> list[Pose]:
        """Hand the chosen poses to the table (the button, or the context menu)."""
        chosen = self.chosen()
        if chosen and self._to_table is not None:
            self._to_table(chosen)
        return chosen

    def state(self) -> dict:
        return {"protein": self.run.protein, "run": self.run.number}

    def close_panel(self) -> None:
        self._show = None
        self._to_table = None
        self._remove = None
