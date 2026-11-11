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

from collections.abc import Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QTreeWidgetItemIterator,
    QWidget,
)

from chemsuite.analysis.fingerprints import KINDS, as_hex, compute_fingerprints
from chemsuite.app.dialogs import FIELD_WIDTH, RowsChoice, ToolDialog
from chemsuite.app.tools import ToolContext
from chemsuite.core.jobs import Cancelled
from chemsuite.descriptors.catalog import CATALOG, COMMON, GROUPS, in_group
from chemsuite.descriptors.compute import compute_descriptors

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 2000
KEY_ROLE = Qt.ItemDataRole.UserRole

_last_keys: tuple[str, ...] = COMMON  # remembered for the session
_last_fingerprints: tuple[str, ...] = ()
FINGERPRINTS = [kind for kind in KINDS.values() if kind.listed]


class DescriptorDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(
            ctx.parent, "Calculate Descriptors", guide="calculate-descriptors-ctrlshiftd"
        )
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        chosen = self.section("Descriptors")
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find a descriptor")
        self.search.setClearButtonEnabled(True)
        chosen.row("Find", self.search)
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setUniformRowHeights(True)
        self.tree.setFixedSize(FIELD_WIDTH, 300)
        self._fill_tree()
        chosen.row(
            "Calculate",
            self.tree,
            help="Tick the descriptors to calculate; point at one for what it is.",
        )
        presets = QWidget()
        row = QHBoxLayout(presets)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        common = QPushButton("Common Set")
        common.setToolTip(", ".join(CATALOG[k].title for k in COMMON))
        everything = QPushButton("All")
        everything.setToolTip("Tick every descriptor")
        none = QPushButton("None")
        none.setToolTip("Untick every descriptor")
        for button in (common, everything, none):
            button.setAutoDefault(False)
            row.addWidget(button)
        row.addStretch(1)
        chosen.row("", presets)

        bits = self.section("Fingerprints")
        self.fingerprints = QListWidget()
        for kind in FINGERPRINTS:
            item = QListWidgetItem(kind.title)
            item.setData(KEY_ROLE, kind.key)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(
                Qt.CheckState.Checked if kind.key in _last_fingerprints else Qt.CheckState.Unchecked
            )
            if kind.help:
                item.setToolTip(kind.help)
            self.fingerprints.addItem(item)
        rows_high = self.fingerprints.sizeHintForRow(0) * len(FINGERPRINTS)
        self.fingerprints.setFixedHeight(rows_high + 2 * self.fingerprints.frameWidth() + 2)
        bits.row(
            "Calculate",
            self.fingerprints,
            help="Each ticked fingerprint becomes a text column holding its bits in "
            "hexadecimal, bit 0 first (the high bit of the first byte).",
        )

        common.clicked.connect(lambda: self.set_checked(COMMON))
        everything.clicked.connect(lambda: self.set_checked(CATALOG))
        none.clicked.connect(lambda: self.set_checked(()))
        self.search.textChanged.connect(self._filter)
        self.tree.itemChanged.connect(lambda *_: self.refresh())
        self.fingerprints.itemChanged.connect(lambda *_: self.refresh())
        self.rows.changed.connect(self.refresh)
        self.structure.currentIndexChanged.connect(lambda _i: self.refresh())
        self.set_checked(_last_keys)

    # -- the list -----------------------------------------------------------------------------
    def _fill_tree(self) -> None:
        for group in GROUPS:
            parent = QTreeWidgetItem(self.tree, [group])
            parent.setFlags(
                parent.flags() | Qt.ItemFlag.ItemIsAutoTristate | Qt.ItemFlag.ItemIsUserCheckable
            )
            for spec in in_group(group):
                item = QTreeWidgetItem(parent, [spec.title])
                item.setData(0, KEY_ROLE, spec.key)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(0, Qt.CheckState.Unchecked)
                if spec.help:
                    item.setToolTip(0, spec.help)
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
        self.refresh()

    def chosen_keys(self) -> list[str]:
        return [
            item.data(0, KEY_ROLE)
            for item in self._items()
            if item.checkState(0) == Qt.CheckState.Checked
        ]

    def chosen_fingerprints(self) -> list[str]:
        items = (self.fingerprints.item(i) for i in range(self.fingerprints.count()))
        return [item.data(KEY_ROLE) for item in items if item.checkState() == Qt.CheckState.Checked]

    def set_fingerprints(self, keys) -> None:
        wanted = set(keys)
        for i in range(self.fingerprints.count()):
            item = self.fingerprints.item(i)
            checked = item.data(KEY_ROLE) in wanted
            item.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)

    def _filter(self, text: str) -> None:
        needle = text.strip().casefold()
        for i in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(i)
            shown = 0
            for j in range(group.childCount()):
                item = group.child(j)
                match = not needle or needle in (item.text(0) + " " + item.toolTip(0)).casefold()
                item.setHidden(not match)
                shown += match
            group.setHidden(shown == 0)

    # -- what will happen ---------------------------------------------------------------------
    def _rows(self) -> int:
        return len(self.rows.oids())

    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to calculate from."
        if not self.chosen_keys() and not self.chosen_fingerprints():
            return "Tick at least one descriptor or fingerprint."
        return ""

    def verb(self) -> str:
        counts = [
            (count, one)
            for count, one in (
                (len(self.chosen_keys()), "Descriptor"),
                (len(self.chosen_fingerprints()), "Fingerprint"),
            )
            if count
        ]
        parts = [
            one if count == 1 and len(counts) == 1 else f"{count} {one}{'s' * (count > 1)}"
            for count, one in counts
        ]
        return f"Calculate {' and '.join(parts)} for {self._rows():,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        keys = self.chosen_keys()
        fingerprints = self.chosen_fingerprints()
        global _last_keys, _last_fingerprints
        _last_keys = tuple(keys)
        _last_fingerprints = tuple(fingerprints)
        return self.calculate(keys, fingerprints)

    def calculate(self, keys: list[str], fingerprints: Sequence[str] = ()):
        ctx = self.ctx
        name = self.structure.currentText()
        oids = self.rows.oids()
        cells = ctx.structure_inputs(name, oids)  # workers parse unparsed cells

        n = len(cells)
        steps = (1 if keys else 0) + len(fingerprints)

        def work(job):
            """``(rows done, {column: values})``: descriptors, then each fingerprint."""
            executor = ctx.workers_for(len(cells), POOL_MIN_ROWS)
            results: dict[str, object] = {}
            step = 0

            def report(done: int, _total: int) -> None:
                job.report(step * n + done, steps * n)

            if keys:
                try:
                    columns = compute_descriptors(
                        cells, keys, executor=executor, progress=report, cancel=job.cancel_token
                    )
                except Cancelled as stopped:  # keep the descriptors of the first rows
                    if stopped.partial is None:
                        raise
                    columns = stopped.partial
                    first = len(next(iter(columns.values())))
                    raise Cancelled((first, {CATALOG[k].title: columns[k] for k in keys})) from None
                results.update({CATALOG[key].title: columns[key] for key in keys})
                step += 1
            for key in fingerprints:
                try:
                    fps = compute_fingerprints(
                        cells, key, executor=executor, progress=report, cancel=job.cancel_token
                    )
                except Cancelled:  # a fingerprint is written for every row or not at all
                    raise Cancelled((n, results) if results else None) from None
                results[KINDS[key].column] = as_hex(fps)
                step += 1
            return n, results

        def done(result) -> None:
            finished, columns = result  # fewer rows after a cancel
            done_oids = oids[:finished]
            names = ctx.apply_results(
                {name: (done_oids, values) for name, values in columns.items()},
                label=f"Calculate {len(columns)} descriptor(s)",
            )
            ctx.show_message(f"Added {len(names)} column(s) for {len(done_oids):,} row(s)")

        what = len(keys) + len(fingerprints)
        return ctx.run(f"Calculating {what} descriptor(s)…", work, done, keep_partial=True)
