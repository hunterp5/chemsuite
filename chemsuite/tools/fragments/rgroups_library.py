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

from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QPlainTextEdit,
    QVBoxLayout,
)

from chemsuite.analysis.reactions import Library, TemplateError, combinations_count
from chemsuite.analysis.rgroup_library import (
    POOL_MIN_COMBINATIONS,
    CoreError,
    core_positions,
    enumerate_rgroups,
    is_hydrogen,
    parse_group,
)
from chemsuite.app.document.tools import ToolContext
from chemsuite.app.widgets.dialogs import RowsChoice, ToolDialog, inline
from chemsuite.app.widgets.reaction_text import ReactionText
from chemsuite.core.jobs import Cancelled
from chemsuite.core.wording import plural
from chemsuite.tools.reactions.enumerate import (
    FILE,
    FILE_FILTER,
    MAX_COMBINATIONS,
    finish_library,
    load_blocks,
)

TYPED = "\x00typed"  # the combo item that asks for a typed list
MAX_POSITIONS = 6


def _from_file(path: Path, add_hydrogen: bool) -> list[str]:
    """The substituents in *path*: its first column that holds attachment points."""
    try:
        found = load_blocks(path, prepared)
    except TemplateError:
        found = []
    if not any(parse_group(text, 1) is not None for text in found):
        raise CoreError(
            f"Nothing in {path.name} can be used as a substituent: it needs an attachment "
            "point, so prepare the reagents first (Prepare Reagents…)."
        )
    if add_hydrogen and not any(is_hydrogen(text) for text in found):
        found = [*found, "H"]
    return found


def prepared(texts: list[str]) -> bool:
    """Whether *texts* can be used as substituents: at least one has an attachment point.

    A saved set of prepared reagents holds the reagents beside their fragments, so this is
    what picks the fragments out of such a file.
    """
    return any(parse_group(text, 1) is not None or is_hydrogen(text) for text in texts)


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
        self.files: dict[int, Path] = {}  # position index -> a file of prepared reagents

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
            combo.addItem("From a File…", FILE)
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
                "decomposition wrote), a typed list, or a file of prepared reagents "
                "(Prepare Reagents…). “and H” also leaves the position unsubstituted.",
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
        chosen = combo.currentData()
        if chosen == TYPED:
            text = ask_groups(self, f"R{slot + 1}", self.typed.get(slot, ""))
            if text.strip():
                self.files.pop(slot, None)
                self.typed[slot] = text
                self._relabel_source(slot)
            else:
                combo.setCurrentIndex(0)
                self.typed.pop(slot, None)
        elif chosen == FILE:
            path, _ = QFileDialog.getOpenFileName(self, f"R{slot + 1}", "", FILE_FILTER)
            if path:
                self.typed.pop(slot, None)
                self.set_file(slot, Path(path))
            else:
                combo.setCurrentIndex(0)
                self.files.pop(slot, None)
        else:
            self.typed.pop(slot, None)
            self.files.pop(slot, None)
        self.refresh()

    def _relabel_source(self, slot: int) -> None:
        """Name the chosen file or typed list on its combo item."""
        combo = self.sources[slot]
        text = self.typed.get(slot, "")
        count = len(self._split(text)) if text else 0
        combo.setItemText(
            combo.findData(TYPED), f"Typed: {plural(count, 'group')}" if count else "Typed…"
        )
        path = self.files.get(slot)
        combo.setItemText(combo.findData(FILE), f"File: {path.name}" if path else "From a File…")

    def set_file(self, slot: int, path: Path) -> None:
        """Take position *slot*'s substituents from the structures in *path*."""
        self.files[slot] = path
        self._relabel_source(slot)
        self.sources[slot].setCurrentIndex(self.sources[slot].findData(FILE))
        self.refresh()

    def uses_file(self, slot: int) -> bool:
        return self.sources[slot].currentData() == FILE

    @staticmethod
    def _split(text: str) -> list[str]:
        seen: dict[str, None] = {}
        for line in text.replace(",", "\n").splitlines():
            if line.strip():
                seen.setdefault(line.strip(), None)
        return list(seen)

    def groups_for(self, slot: int) -> list[str]:
        """The substituents chosen for this position, with H last when it is wanted.

        A file's substituents are not read here: they are loaded in the job, so the dialog
        never blocks on a large file and the count is not known until the run.
        """
        combo = self.sources[slot]
        chosen = combo.currentData()
        if chosen == FILE:
            return []
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

    def combinations(self) -> int | None:
        """How many products will be made, or None while a file's substituents are unknown."""
        if any(self.uses_file(i) for i in self._slots()):
            return None
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
            if not self.uses_file(slot) and not self.groups_for(slot):
                return f"Choose substituents for R{slot + 1}."
        count = self.combinations()
        if count is not None and count > MAX_COMBINATIONS:
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
        #: a file to read in the job, or the substituents already known
        sources: list[list[str] | Path] = [
            self.files[i] if self.uses_file(i) else self.groups_for(i) for i in slots
        ]
        extra_h = [self.hydrogens[i].isChecked() for i in slots]
        names = [f"R{i + 1}" for i in slots]

        def work(job):
            groups = [
                _from_file(source, add_h) if isinstance(source, Path) else source
                for source, add_h in zip(sources, extra_h, strict=True)
            ]
            count = combinations_count([len(column) for column in groups])
            if count > MAX_COMBINATIONS:
                raise CoreError(
                    f"{count:,} combinations is more than {MAX_COMBINATIONS:,}: use fewer "
                    "substituents."
                )
            try:
                return groups, enumerate_rgroups(
                    core,
                    groups,
                    executor=ctx.workers_for(count, POOL_MIN_COMBINATIONS),
                    progress=job.report,
                    cancel=job.cancel_token,
                )
            except Cancelled as stopped:  # keep the products already made
                raise Cancelled((groups, stopped.partial) if stopped.partial else None) from None

        def done(result: tuple[list[list[str]], Library]) -> None:
            groups, library = result
            finish_library(ctx, groups, library, names=names)

        return ctx.run("Enumerating a library…", work, done, keep_partial=True)
