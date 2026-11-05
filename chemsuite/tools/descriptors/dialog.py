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
"""The Calculate Descriptors dialog."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QTreeWidgetItemIterator,
    QVBoxLayout,
)

from chemsuite.app.tools import ToolContext
from chemsuite.descriptors.catalog import CATALOG, COMMON, GROUPS, in_group
from chemsuite.descriptors.compute import compute_descriptors

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 2000
KEY_ROLE = Qt.ItemDataRole.UserRole

_last_keys: tuple[str, ...] = COMMON  # remembered for the session


class DescriptorDialog(QDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent)
        self.setWindowTitle("Calculate Descriptors")
        self.ctx = ctx
        self.job = None

        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter descriptors")
        self.search.setClearButtonEnabled(True)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Descriptor", "Description"])
        self.tree.setColumnWidth(0, 230)
        self.tree.setUniformRowHeights(True)
        self._fill_tree()
        common = QPushButton("Common Set")
        common.setToolTip(", ".join(CATALOG[k].title for k in COMMON))
        everything = QPushButton("All")
        none = QPushButton("None")
        presets = QHBoxLayout()
        for button in (common, everything, none):
            presets.addWidget(button)
        presets.addStretch(1)
        self.selected_only = QCheckBox(f"Selected rows only ({len(ctx.selection):,})")
        self.selected_only.setEnabled(len(ctx.selection) > 0)
        self.overwrite = QCheckBox("Update columns that already have these names")
        self.overwrite.setChecked(True)
        self.overwrite.setToolTip(
            "On: recalculating replaces the values in an existing “MW” column.\n"
            "Off: results go to new columns such as “MW (2)”."
        )
        self.summary = QLabel()
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Calculate")

        form = QFormLayout()
        form.addRow("Structures:", self.structure)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.search)
        layout.addWidget(self.tree, 1)
        layout.addLayout(presets)
        layout.addWidget(self.selected_only)
        layout.addWidget(self.overwrite)
        layout.addWidget(self.summary)
        layout.addWidget(self.buttons)
        self.resize(640, 640)

        common.clicked.connect(lambda: self.set_checked(COMMON))
        everything.clicked.connect(lambda: self.set_checked(CATALOG))
        none.clicked.connect(lambda: self.set_checked(()))
        self.search.textChanged.connect(self._filter)
        self.tree.itemChanged.connect(lambda *_: self._update_summary())
        self.selected_only.toggled.connect(lambda _on: self._update_summary())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.set_checked(_last_keys)

    # -- the list -----------------------------------------------------------------------------
    def _fill_tree(self) -> None:
        for group in GROUPS:
            parent = QTreeWidgetItem(self.tree, [group])
            parent.setFlags(
                parent.flags() | Qt.ItemFlag.ItemIsAutoTristate | Qt.ItemFlag.ItemIsUserCheckable
            )
            for spec in in_group(group):
                item = QTreeWidgetItem(parent, [spec.title, spec.help])
                item.setData(0, KEY_ROLE, spec.key)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(0, Qt.CheckState.Unchecked)
                if spec.help:
                    item.setToolTip(0, spec.help)
                    item.setToolTip(1, spec.help)
        self.tree.expandAll()

    def _items(self):
        iterator = QTreeWidgetItemIterator(self.tree)
        while iterator.value():
            item = iterator.value()
            if item.data(0, KEY_ROLE):
                yield item
            iterator += 1

    def set_checked(self, keys) -> None:
        wanted = set(keys)
        self.tree.blockSignals(True)
        for item in self._items():
            state = (
                Qt.CheckState.Checked
                if item.data(0, KEY_ROLE) in wanted
                else Qt.CheckState.Unchecked
            )
            item.setCheckState(0, state)
        self.tree.blockSignals(False)
        self._update_summary()

    def chosen_keys(self) -> list[str]:
        return [
            item.data(0, KEY_ROLE)
            for item in self._items()
            if item.checkState(0) == Qt.CheckState.Checked
        ]

    def _filter(self, text: str) -> None:
        needle = text.strip().casefold()
        for i in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(i)
            shown = 0
            for j in range(group.childCount()):
                item = group.child(j)
                match = not needle or needle in (item.text(0) + " " + item.text(1)).casefold()
                item.setHidden(not match)
                shown += match
            group.setHidden(shown == 0)

    def _rows(self) -> int:
        return len(self.ctx.selection) if self.selected_only.isChecked() else len(self.ctx.table)

    def _update_summary(self) -> None:
        count = len(self.chosen_keys())
        self.summary.setText(f"{count} descriptor(s) for {self._rows():,} row(s)")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            count > 0 and self._rows() > 0 and bool(self.structure.currentText())
        )

    # -- running ------------------------------------------------------------------------------
    def accept(self) -> None:
        keys = self.chosen_keys()
        if not keys or not self._rows():
            return
        global _last_keys
        _last_keys = tuple(keys)
        self.job = self.calculate(keys)
        super().accept()

    def calculate(self, keys: list[str]):
        ctx = self.ctx
        table = ctx.table
        name = self.structure.currentText()
        oids = ctx.target_oids(self.selected_only.isChecked())
        column = table.column(name)
        slots = table.slots_of(oids)
        # Unparsed cells go as text; workers parse them.
        cells = [
            blob if blob is not None else (text or None)
            for blob, text in zip(column.blobs(slots), column.displays(slots), strict=True)
        ]
        overwrite = self.overwrite.isChecked()
        use_pool = len(cells) >= POOL_MIN_ROWS

        def work(job):
            return compute_descriptors(
                cells,
                keys,
                executor=ctx.process_pool() if use_pool else None,
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(columns) -> None:
            # After a cancel, the columns cover only the first rows.
            done_oids = oids[: len(next(iter(columns.values())))]
            names = ctx.apply_results(
                {CATALOG[key].title: (done_oids, columns[key]) for key in keys},
                label=f"Calculate {len(keys)} descriptor(s)",
                overwrite=overwrite,
            )
            ctx.show_message(
                f"Added {len(names)} descriptor column(s) for {len(done_oids):,} row(s)"
            )

        return ctx.run(f"Calculating {len(keys)} descriptor(s)…", work, done, keep_partial=True)
