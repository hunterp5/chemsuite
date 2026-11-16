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
"""The preparation report: what Prepare found and decided, one line each.

The panel shows a :class:`~chemsuite.protein.workspace.reports.PrepareReport` (kept with the
document, so it comes back with it). A line marked ⚠ needs the user's eye: a near tie between
two states, a chain break, a residue with no template. Clicking a line selects its residue in
the workbench, and the panel is told who to ask for that (*select*), so it knows nothing about
the table or the tools.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QGuiApplication, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
)

from chemsuite.app.widgets.panel import ToolPanel
from chemsuite.protein.workspace.reports import PrepareReport

#: the tool that reopens this panel with the document
TOOL_ID = "protein-prepare"
WARN = "⚠"
WARN_COLOUR = "#c77700"
HEADERS = ("", "Residue", "What was done")


class PrepareReportPanel(ToolPanel):
    """One line per finding or decision of a preparation."""

    tool_id = TOOL_ID

    def __init__(
        self, report: PrepareReport, select: Callable[[str, int, str], None] | None = None
    ) -> None:
        super().__init__()
        self.report = report
        self._select = select
        self.rows = QTreeWidget()
        self.rows.setHeaderLabels(list(HEADERS))
        self.rows.setRootIsDecorated(False)
        self.rows.setUniformRowHeights(True)
        self.rows.setAlternatingRowColors(True)
        self.rows.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        header = self.rows.header()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.rows.itemSelectionChanged.connect(self._selected)

        self.rows.setContextMenuPolicy(Qt.ContextMenuPolicy.ActionsContextMenu)
        copy = QAction("Copy Report", self.rows)
        copy.setShortcut(QKeySequence.StandardKey.Copy)
        copy.setShortcutContext(Qt.ShortcutContext.WidgetShortcut)
        copy.triggered.connect(self.copy_text)
        self.rows.addAction(copy)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.rows, 1)
        self.fill(report)

    # -- contents -----------------------------------------------------------------------------
    def fill(self, report: PrepareReport) -> None:
        """Show *report* (the protein was prepared again)."""
        self.report = report
        self.rows.clear()
        for row in report.rows:
            item = QTreeWidgetItem([WARN if row.warn else "", row.label, row.text])
            item.setData(0, Qt.ItemDataRole.UserRole, row.residue)
            if row.warn:
                item.setForeground(0, Qt.GlobalColor.darkYellow)
                item.setToolTip(0, "Worth checking: the next best choice was nearly as good")
            self.rows.addTopLevelItem(item)
        self.titleChanged.emit(self.title())

    def title(self) -> str:
        report = self.report
        return f"Prepare: {report.protein} · pH {report.ph:g} · net charge {report.net_charge:+d}"

    def copy_text(self) -> str:
        """The report as lines of text, on the clipboard."""
        lines = ["\t".join(HEADERS[1:])]
        lines += [
            f"{row.label}\t{row.text}" + (f"\t{WARN}" if row.warn else "")
            for row in self.report.rows
        ]
        text = "\n".join(lines)
        QGuiApplication.clipboard().setText(text)
        return text

    # -- the workbench ------------------------------------------------------------------------
    def _selected(self) -> None:
        if self._select is None:
            return
        for item in self.rows.selectedItems():
            residue = item.data(0, Qt.ItemDataRole.UserRole)
            if isinstance(residue, int) and residue >= 0:
                self._select(self.report.protein, residue, item.text(1))
                return

    def state(self) -> dict:
        return {"protein": self.report.protein}

    def close_panel(self) -> None:
        self._select = None
