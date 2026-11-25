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
"""The Structural Alerts dialog."""

from __future__ import annotations

from collections.abc import Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QWidget,
)

from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog
from chemsuite.core.wording import plural
from chemsuite.descriptors.alerts import (
    ALERT_SETS,
    COMMON,
    compute_alerts,
    count_column,
    names_column,
)

#: Below this many rows the work runs on the job's thread: starting worker processes costs more.
POOL_MIN_ROWS = 500
KEY_ROLE = Qt.ItemDataRole.UserRole

_last_keys: tuple[str, ...] = COMMON  # remembered for the session


class AlertsDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "Structural Alerts", guide="structural-alerts")
        self.ctx = ctx

        inputs = self.section("Input")
        self.structure = QComboBox()
        self.structure.addItems(ctx.structure_columns())
        default = ctx.default_structure_column()
        if default:
            self.structure.setCurrentText(default)
        inputs.row("Structures", self.structure)
        self.rows = inputs.row("Rows", RowsChoice(ctx))

        chosen = self.section("Alert Sets")
        self.sets = QListWidget()
        for spec in ALERT_SETS.values():
            item = QListWidgetItem(spec.title)
            item.setData(KEY_ROLE, spec.key)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)
            item.setToolTip(spec.help)
            self.sets.addItem(item)
        self.sets.setFixedHeight(self.sets.sizeHintForRow(0) * 8 + 2 * self.sets.frameWidth() + 2)
        chosen.row(
            "Match",
            self.sets,
            help="Tick the filter sets to match; point at one for what it covers.",
        )
        buttons = QWidget()
        row = QHBoxLayout(buttons)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        common = QPushButton("Common Set")
        common.setToolTip(", ".join(ALERT_SETS[key].title for key in COMMON))
        everything = QPushButton("All")
        everything.setToolTip("Tick every filter set")
        none = QPushButton("None")
        none.setToolTip("Untick every filter set")
        for button in (common, everything, none):
            button.setAutoDefault(False)
            row.addWidget(button)
        row.addStretch(1)
        chosen.row("", buttons)

        output = self.section("Output")
        self.name_them = QCheckBox("Name the alerts hit")
        self.name_them.setChecked(True)
        output.row(
            "Columns",
            self.name_them,
            help="Each set becomes a count column; with this on, also a text column listing the "
            "alerts the row hits.",
        )

        common.clicked.connect(lambda: self.set_checked(COMMON))
        everything.clicked.connect(lambda: self.set_checked(ALERT_SETS))
        none.clicked.connect(lambda: self.set_checked(()))
        self.sets.itemChanged.connect(lambda *_: self.refresh())
        self.name_them.toggled.connect(lambda *_: self.refresh())
        self.rows.changed.connect(self.refresh)
        self.structure.currentIndexChanged.connect(lambda _i: self.refresh())
        self.set_checked(_last_keys)

    # -- the list -----------------------------------------------------------------------------
    def _items(self):
        return (self.sets.item(i) for i in range(self.sets.count()))

    def set_checked(self, keys) -> None:
        wanted = set(keys)
        self.sets.blockSignals(True)
        for item in self._items():
            checked = item.data(KEY_ROLE) in wanted
            item.setCheckState(Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked)
        self.sets.blockSignals(False)
        self.refresh()

    def chosen_keys(self) -> list[str]:
        return [
            item.data(KEY_ROLE)
            for item in self._items()
            if item.checkState() == Qt.CheckState.Checked
        ]

    # -- what will happen ---------------------------------------------------------------------
    def problem(self) -> str:
        if not self.structure.currentText():
            return "There is no structure column to match against."
        if not self.chosen_keys():
            return "Tick at least one filter set."
        return ""

    def verb(self) -> str:
        count = len(self.chosen_keys())
        what = "Alert Set" if count == 1 else f"{count} Alert Sets"
        return f"Match {what} for {len(self.rows.oids()):,} Rows"

    # -- running ------------------------------------------------------------------------------
    def run(self):
        keys = self.chosen_keys()
        global _last_keys
        _last_keys = tuple(keys)
        return self.match(keys, names=self.name_them.isChecked())

    def match(self, keys: Sequence[str], *, names: bool = True):
        ctx = self.ctx
        oids = self.rows.oids()
        cells = ctx.flat_structure_inputs(self.structure.currentText(), oids)  # 2D only
        wanted = [count_column(ALERT_SETS[key].title) for key in keys]
        if names:
            wanted += [names_column(ALERT_SETS[key].title) for key in keys]

        def work(job):
            return compute_alerts(
                cells,
                keys,
                executor=ctx.workers_for(len(cells), POOL_MIN_ROWS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(columns) -> None:
            kept = {name: columns[name] for name in wanted}
            finished = len(next(iter(kept.values())))  # fewer rows after a cancel
            written = ctx.apply_results(
                {name: (oids[:finished], column) for name, column in kept.items()},
                label=f"Match {plural(len(keys), 'alert set')}",
            )
            ctx.show_message(
                f"Added {plural(len(written), 'column')} for {plural(finished, 'row')}"
            )

        return ctx.run("Matching structural alerts…", work, done, keep_partial=True)
