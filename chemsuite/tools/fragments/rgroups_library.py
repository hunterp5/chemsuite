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

"""Chemistry ▸ Enumerate ▸ R-Groups on a Core: every combination of substituents."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QPlainTextEdit,
    QVBoxLayout,
)

from chemsuite.analysis.reactions import Library, combinations_count
from chemsuite.analysis.rgroup_library import (
    POOL_MIN_COMBINATIONS,
    CoreError,
    core_positions,
    enumerate_rgroups,
    is_hydrogen,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.app.widgets.reaction_text import ReactionText
from chemsuite.core.wording import plural
from chemsuite.tools.reactions.enumerate import MAX_COMBINATIONS, finish_library

TYPED = "\x00typed"  # the combo item that asks for a typed list
MAX_POSITIONS = 6


class GroupsDialog(QDialog):
    """A list of substituents, one per line."""

    def __init__(self, parent, position: str, text: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle(position)
        self.edit = QPlainTextEdit()
        self.edit.setPlaceholderText("*Cl\n*OC\nH")
        self.edit.setPlainText(text)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(self.edit)
        layout.addWidget(buttons)

    def text(self) -> str:
        return self.edit.toPlainText()


def ask_groups(parent, position: str, text: str = "") -> str:
    dialog = GroupsDialog(parent, position, text)
    return dialog.text() if dialog.exec() else ""


class RGroupLibraryDialog(ToolDialog):
    def __init__(self, ctx: ToolContext) -> None:
        super().__init__(ctx.parent, "R-Groups on a Core", guide="r-groups-on-a-core")
        self.ctx = ctx
        self.typed: dict[int, str] = {}  # position index -> typed list

        core = self.section("Core")
        self.core = core.row(
            "Core",
            ReactionText("SMILES with [*:1], [*:2] …"),
            help="The core with numbered attachment points, as R-Group Decomposition and "
            "Central Ring write them. Draw… and use the sketcher's R-group button.",
        )
        self.core.set_reactions(False)

        groups = self.section("Substituents")
        self._groups_section = groups
        self.sources: list[QComboBox] = []
        self.hydrogens: list[QCheckBox] = []
        for i in range(MAX_POSITIONS):
            combo = QComboBox()
            for name in self.ctx.table.column_names:
                combo.addItem(name, name)
            combo.addItem("Typed…", TYPED)
            combo.activated.connect(lambda _index, slot=i: self._source_chosen(slot))
            combo.currentIndexChanged.connect(lambda _index: self.refresh())
            box = QCheckBox("and H")
            box.toggled.connect(lambda *_: self.refresh())
            self.hydrogens.append(box)
            self.sources.append(combo)
            groups.row(
                f"R{i + 1}",
                inline(combo, box, stretch=0),
                help="The substituents for this position: a column (an R1 column a "
                "decomposition wrote), or a typed list. “and H” also leaves the position "
                "unsubstituted.",
            )
        self.rows = groups.row("Rows", RowsChoice(ctx))

        self.core.changed.connect(self.refresh)
        self.rows.changed.connect(lambda *_: self.refresh())
        self.refresh()

    # -- the core and its positions -------------------------------------------------------------
    def positions(self) -> tuple[int, ...]:
        return core_positions(self.core.text())

    def _slots(self) -> range:
        return range(min(len(self.positions()), MAX_POSITIONS))

    def _source_chosen(self, slot: int) -> None:
        combo = self.sources[slot]
        if combo.currentData() != TYPED:
            self.typed.pop(slot, None)
            self.refresh()
            return
        text = ask_groups(self, f"R{slot + 1}", self.typed.get(slot, ""))
        if text.strip():
            self.typed[slot] = text
            count = len(self._split(text))
            combo.setItemText(combo.count() - 1, f"Typed: {plural(count, 'group')}")
        else:
            combo.setCurrentIndex(0)
            self.typed.pop(slot, None)
        self.refresh()

    @staticmethod
    def _split(text: str) -> list[str]:
        seen: dict[str, None] = {}
        for line in text.replace(",", "\n").splitlines():
            if line.strip():
                seen.setdefault(line.strip(), None)
        return list(seen)

    def groups_for(self, slot: int) -> list[str]:
        """The substituents chosen for this position, with H last when it is wanted."""
        combo = self.sources[slot]
        chosen = combo.currentData()
        if chosen == TYPED:
            found = self._split(self.typed.get(slot, ""))
        else:
            found = self._column_values(chosen)
        if self.hydrogens[slot].isChecked() and not any(is_hydrogen(g) for g in found):
            found = [*found, "H"]
        return found

    def _column_values(self, name: str) -> list[str]:
        """The distinct labels in *name* among the chosen rows; "H" is the "and H" box's job."""
        found = self.ctx.distinct_values(name, self.rows.oids())
        return [text for text in found if not is_hydrogen(text)]

    def combinations(self) -> int:
        return combinations_count([len(self.groups_for(i)) for i in self._slots()])

    # -- ToolDialog -------------------------------------------------------------------------------
    def refresh(self) -> None:
        if hasattr(self, "_groups_section"):
            count = len(self.positions())
            for i, combo in enumerate(self.sources):
                visible = i < min(count, MAX_POSITIONS)
                self._groups_section.set_row_visible(combo.parentWidget(), visible)
        super().refresh()

    def problem(self) -> str:
        try:
            positions = self.positions()
        except CoreError as err:
            return str(err)
        text = self.core.text().strip()
        if not text:
            return "Enter a core with numbered attachment points, such as c1ccc([*:1])cc1."
        if not positions:
            return "The core has no numbered attachment point ([*:1], [*:2] …)."
        if len(positions) > MAX_POSITIONS:
            return f"Cores with more than {MAX_POSITIONS} attachment points are not supported."
        for slot in self._slots():
            if not self.groups_for(slot):
                return f"Choose substituents for R{slot + 1}."
        count = self.combinations()
        if count > MAX_COMBINATIONS:
            return (
                f"{count:,} combinations is more than {MAX_COMBINATIONS:,}: use fewer substituents."
            )
        return ""

    def verb(self) -> str:
        count = self.combinations()
        return f"Enumerate {count:,} Combinations" if count else "Enumerate"

    def run(self):
        ctx = self.ctx
        core = self.core.text().strip()
        slots = list(self._slots())
        groups = [self.groups_for(i) for i in slots]
        names = [f"R{i + 1}" for i in slots]
        total = self.combinations()

        def work(job):
            return enumerate_rgroups(
                core,
                groups,
                executor=ctx.workers_for(total, POOL_MIN_COMBINATIONS),
                progress=job.report,
                cancel=job.cancel_token,
            )

        def done(library: Library) -> None:
            finish_library(ctx, groups, library, names=names)

        return ctx.run("Enumerating a library…", work, done, keep_partial=True)
