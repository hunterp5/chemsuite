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
"""The poses of a docking run, grouped by ligand, beside the pocket they were docked into.

The panel shows a :class:`~chemsuite.protein.workspace.poses.DockedPoses` (kept with the
document, so the poses come back with it). Each ligand is a row with its best affinity, and
each pose under it; selecting a pose shows it in the pocket, and the panel is told who to ask
for that (*show*), so it knows nothing about the table, the tools or RDKit.

**Copy to Table** hands the chosen poses to whoever opened the panel (*to_table*): a pose is
only a table row when the chemist asks for one.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
)

from chemsuite.app.widgets.panel import ToolPanel
from chemsuite.protein.workspace.poses import RESCORED, DockedPoses, Pose

#: the tool that reopens this panel with the document
TOOL_ID = "protein-poses"
HEADERS = ("Ligand", "Pose", "Affinity", "RMSD l.b.", "RMSD u.b.")
#: a pose is held in its row, under this role
POSE = Qt.ItemDataRole.UserRole


class PosesPanel(ToolPanel):
    """One row per ligand, its poses under it."""

    tool_id = TOOL_ID

    def __init__(
        self,
        run: DockedPoses,
        show: Callable[[Pose], None] | None = None,
        to_table: Callable[[Sequence[Pose]], None] | None = None,
    ) -> None:
        super().__init__()
        self.run = run
        self._show = show
        self._to_table = to_table
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
        layout.addWidget(self.rows, 1)
        layout.addLayout(footer)
        self.fill(run)

    # -- contents -----------------------------------------------------------------------------
    def fill(self, run: DockedPoses) -> None:
        """Show *run* (the protein was docked into again)."""
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
        what = "Rescored" if run.kind == RESCORED else "Poses"
        return f"{what}: {run.protein} · {len(run.poses)} of {len(run.ligands)} ligands"

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

    def send_to_table(self) -> list[Pose]:
        """Hand the chosen poses to the table (the button, or the context menu)."""
        chosen = self.chosen()
        if chosen and self._to_table is not None:
            self._to_table(chosen)
        return chosen

    def state(self) -> dict:
        return {"protein": self.run.protein}

    def close_panel(self) -> None:
        self._show = None
        self._to_table = None
