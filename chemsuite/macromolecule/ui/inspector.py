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
"""Tools ▸ Inspect…: everything a preparation would find and decide, to tick and overrule.

The panel shows a :class:`~chemsuite.macromolecule.prepare.review.Survey` as a tree:

- **Issues**, grouped by kind. One the preparation can act on (an alternate location, a
  missing side chain) has a box to tick: unticked, it is left as the file has it. The rest (a
  chain break, a residue with no template) are there to be seen.
- **Choices**, grouped by kind: every flip, His state, protonation, rotor, ligand form and
  water orientation, the state chosen, and what the next best would cost; and each residue's
  alternate location (with its occupancies) and whether each crystal water is kept.
  Double-click the state to choose another; each is listed with its cost against the best (an
  alternate location with its occupancy). A choice made by hand is in bold.

**Re-survey** solves again with the hand-made choices held, so their effect on the
neighbouring sites shows before anything is applied. **Apply** prepares the structure as
reviewed, in one undoable step. **Reset** puts every line back as the survey proposed.

The panel is told who to ask (*select*, *resurvey*, *apply*), so it knows nothing about the
table, the tools or the document.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping

from PySide6.QtCore import QModelIndex, Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QStyledItemDelegate,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.widgets.panel import ToolPanel
from chemsuite.macromolecule.prepare import inspect as inspection
from chemsuite.macromolecule.prepare.review import KINDS, Choice, Issue, Survey

#: the tool that reopens this panel with the document
TOOL_ID = "protein-inspect"
WARN = "⚠"
HEADERS = ("", "Residue", "Found or decided", "Next best")
#: column of the state chooser, and the roles a line keeps its data under
STATE_COLUMN = 2
LINE = Qt.ItemDataRole.UserRole  # an Issue or a Choice
CHOSEN = Qt.ItemDataRole.UserRole + 1  # a choice line's state label now
#: how an issue's kind reads as a group
ISSUE_GROUPS = {
    inspection.ALTLOCS: "Alternate locations",
    inspection.MISSING: "Missing atoms",
    inspection.BREAK: "Chain breaks",
    inspection.GAP: "Gaps in the sequence",
    inspection.TRUNCATED: "Chain ends not modelled",
    inspection.FILLED: "Residues built",
    inspection.CAP: "Caps",
    inspection.UNKNOWN: "No template",
    inspection.LINK: "Covalent links",
    inspection.DISULFIDE: "Disulfides",
    inspection.METAL: "Metal sites",
}

Select = Callable[[str, int, str], None]  # protein, residue, the line's label
Decide = Callable[[tuple[str, ...], dict[str, str]], None]


class _StateDelegate(QStyledItemDelegate):
    """A choice line's state, chosen from a list of its states with their costs."""

    def __init__(self, panel: InspectorPanel) -> None:
        super().__init__(panel)
        self.panel = panel

    def createEditor(self, parent: QWidget, option, index: QModelIndex) -> QWidget | None:
        choice = self.panel.choice_at(index)
        if choice is None or len(choice.states) < 2:
            return None
        box = QComboBox(parent)
        for state in choice.states:
            if not choice.scored:  # chosen, not solved for: no costs to show
                text = f"{state.label} ({state.note})" if state.note else state.label
            else:
                cost = "best" if state.cost <= 1e-6 else f"+{state.cost:.1f}"
                text = f"{state.label} ({cost})"
            box.addItem(text, state.label)
        return box

    def setEditorData(self, editor: QWidget, index: QModelIndex) -> None:
        if isinstance(editor, QComboBox):
            found = editor.findData(index.data(CHOSEN))
            editor.setCurrentIndex(max(found, 0))

    def setModelData(self, editor: QWidget, model, index: QModelIndex) -> None:
        if isinstance(editor, QComboBox):
            self.panel.choose(index, str(editor.currentData()))


class InspectorPanel(ToolPanel):
    """Every issue and choice of a preparation, to tick and overrule before applying."""

    tool_id = TOOL_ID
    #: the chemist changed a line (to save what has been decided so far)
    decided = Signal()

    def __init__(
        self,
        found: Survey,
        select: Select | None = None,
        resurvey: Decide | None = None,
        apply: Decide | None = None,
    ) -> None:
        super().__init__()
        #: what is shown: after a re-survey, solved with the hand-made choices held
        self.survey = found
        #: the survey run with nothing held — what the preparation proposed on its own. The
        #: chemist's answers are read against it, Reset goes back to it, and the report's
        #: "proposed" comes from it (never from a survey that already held their choice).
        self.base = found
        self._select = select
        self._resurvey = resurvey
        self._apply = apply
        self.skipped: set[str] = set()
        self.overrides: dict[str, str] = {}

        self.rows = QTreeWidget()
        self.rows.setHeaderLabels(list(HEADERS))
        self.rows.setUniformRowHeights(True)
        self.rows.setAlternatingRowColors(True)
        self.rows.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.rows.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked
                                  | QAbstractItemView.EditTrigger.SelectedClicked)  # fmt: skip
        self.rows.setItemDelegateForColumn(STATE_COLUMN, _StateDelegate(self))
        header = self.rows.header()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.rows.itemSelectionChanged.connect(self._selected)
        self.rows.itemChanged.connect(self._ticked)

        self.only_warnings = QToolButton()
        self.only_warnings.setText(f"{WARN} Only")
        self.only_warnings.setCheckable(True)
        self.only_warnings.setToolTip("Show only the lines worth a look: near ties, chain "
                                      "breaks, residues with no template")  # fmt: skip
        self.only_warnings.toggled.connect(lambda *_: self._filter())
        self.reset_button = QPushButton("Reset")
        self.reset_button.setToolTip("Every line back as the survey proposed")
        self.reset_button.clicked.connect(self.reset)
        self.resurvey_button = QPushButton("Re-survey")
        self.resurvey_button.setToolTip("Solve again with the choices made by hand held, to see "
                                        "how the sites around them answer")  # fmt: skip
        self.resurvey_button.clicked.connect(self.resurvey)
        self.apply_button = QPushButton("Apply")
        self.apply_button.setDefault(True)
        self.apply_button.setToolTip("Prepare the structure as reviewed, in one undoable step")
        self.apply_button.clicked.connect(self.apply)
        footer = QHBoxLayout()
        footer.addWidget(self.only_warnings)
        footer.addStretch(1)
        footer.addWidget(self.reset_button)
        footer.addWidget(self.resurvey_button)
        footer.addWidget(self.apply_button)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.rows, 1)
        layout.addLayout(footer)
        self.fill(found)

    # -- contents -----------------------------------------------------------------------------
    def fill(
        self,
        found: Survey,
        skipped: tuple[str, ...] | set[str] = (),
        overrides: Mapping[str, str] | None = None,
        base: Survey | None = None,
    ) -> None:
        """Show *found*, with what was already decided about it (a re-survey keeps it).

        *base* is the survey with nothing held; without one, *found* is taken as it (a first
        survey)."""
        self.survey = found
        self.base = base if base is not None else found
        issue_ids = {issue.id for issue in found.issues}
        choice_ids = {choice.id for choice in found.choices}
        self.skipped = {id for id in skipped if id in issue_ids}
        self.overrides = {id: v for id, v in (overrides or {}).items() if id in choice_ids}
        self.rows.blockSignals(True)
        self.rows.clear()
        issues = self._group("Issues", None)
        for kind, title in ISSUE_GROUPS.items():
            mine = [issue for issue in found.issues if issue.kind == kind]
            if mine:
                group = self._group(title, issues)
                for issue in mine:
                    self._issue(group, issue)
                self._tally(group)
        choices = self._group("Choices", None)
        for kind, title in KINDS.items():
            mine = [choice for choice in found.choices if choice.kind == kind]
            if mine:
                group = self._group(f"{title}s" if not title.endswith("s") else title, choices)
                for choice in mine:
                    self._choice(group, choice)
                self._tally(group)
        issues.setExpanded(True)
        choices.setExpanded(True)
        for top in (issues, choices):
            for i in range(top.childCount()):  # the short groups open, rotors and waters shut
                group = top.child(i)
                group.setExpanded(group.childCount() <= 20 or self._warned(group))
        self.rows.blockSignals(False)
        self._filter()
        self.titleChanged.emit(self.title())

    def _group(self, title: str, parent: QTreeWidgetItem | None) -> QTreeWidgetItem:
        item = QTreeWidgetItem([" ", title, "", ""])
        bold = QFont(item.font(1))
        bold.setBold(True)
        item.setFont(1, bold)
        item.setFirstColumnSpanned(False)
        if parent is None:
            self.rows.addTopLevelItem(item)
        else:
            parent.addChild(item)
        return item

    def _issue(self, group: QTreeWidgetItem, issue: Issue) -> None:
        item = QTreeWidgetItem([WARN if issue.warn else "", issue.label, issue.text, ""])
        item.setData(0, LINE, issue)
        if issue.fixable:
            flags = item.flags() | Qt.ItemFlag.ItemIsUserCheckable
            item.setFlags(flags)
            fixed = issue.id not in self.skipped
            item.setCheckState(2, Qt.CheckState.Checked if fixed else Qt.CheckState.Unchecked)
            item.setToolTip(2, "Ticked: fixed when applied. Unticked: left as the file has it")
        group.addChild(item)

    def _choice(self, group: QTreeWidgetItem, choice: Choice) -> None:
        chosen = self.overrides.get(choice.id) or (choice.state.label if choice.state else "")
        nxt = "" if choice.gap == float("inf") else f"+{choice.gap:.1f}"
        text = self._choice_text(choice, chosen)
        item = QTreeWidgetItem([WARN if choice.warn else "", choice.label, text, nxt])
        item.setData(0, LINE, choice)
        item.setData(STATE_COLUMN, CHOSEN, chosen)
        if len(choice.states) > 1:
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
            item.setToolTip(STATE_COLUMN, "Double-click to choose another state")
        if choice.id in self.overrides:
            self._embolden(item, True)
        group.addChild(item)

    @staticmethod
    def _choice_text(choice: Choice, label: str) -> str:
        text = f"{choice.what}: {label}"
        if choice.pka is not None:
            text += f" · pKa {choice.pka:.1f}"
        if choice.population is not None:
            text += f", {choice.population:.0%}"
        return text

    @staticmethod
    def _embolden(item: QTreeWidgetItem, on: bool) -> None:
        font = QFont(item.font(STATE_COLUMN))
        font.setBold(on)
        item.setFont(STATE_COLUMN, font)

    @staticmethod
    def _warned(group: QTreeWidgetItem) -> bool:
        return any(group.child(i).text(0) == WARN for i in range(group.childCount()))

    def _tally(self, group: QTreeWidgetItem) -> None:
        warned = sum(1 for i in range(group.childCount()) if group.child(i).text(0) == WARN)
        group.setText(2, f"{group.childCount()}" + (f" · {warned} {WARN}" if warned else ""))

    def title(self) -> str:
        found = self.survey
        changed = len(self.skipped) + len(self.overrides)
        tail = f" · {changed} by hand" if changed else ""
        return f"Inspect: {found.protein} · pH {found.options.ph:g}{tail}"

    # -- lines --------------------------------------------------------------------------------
    def lines(self) -> list[QTreeWidgetItem]:
        """Every issue and choice line (not the groups)."""
        out = []
        for t in range(self.rows.topLevelItemCount()):
            top = self.rows.topLevelItem(t)
            for g in range(top.childCount()):
                group = top.child(g)
                out += [group.child(i) for i in range(group.childCount())]
        return out

    def line(self, id: str) -> QTreeWidgetItem | None:
        for item in self.lines():
            found = item.data(0, LINE)
            if isinstance(found, (Issue, Choice)) and found.id == id:
                return item
        return None

    def proposed(self, id: str) -> str:
        """What the preparation chose for the choice *id* on its own (from :attr:`base`)."""
        choice = self.base.choice(id)
        return choice.state.label if choice is not None and choice.state else ""

    def choice_at(self, index: QModelIndex) -> Choice | None:
        item = self.rows.itemFromIndex(index)
        found = item.data(0, LINE) if item is not None else None
        return found if isinstance(found, Choice) else None

    def choose(self, index: QModelIndex, label: str) -> None:
        """Hold the choice at *index* to the state *label* (its proposal clears the hold)."""
        item = self.rows.itemFromIndex(index)
        choice = self.choice_at(index)
        if item is None or choice is None:
            return
        self.set_state(choice.id, label)

    def set_state(self, id: str, label: str) -> bool:
        """Hold the choice *id* to the state *label*; its own proposal lifts the hold."""
        item = self.line(id)
        choice = item.data(0, LINE) if item is not None else None
        if not isinstance(choice, Choice) or all(s.label != label for s in choice.states):
            return False
        proposed = self.proposed(id)
        if label == proposed:
            self.overrides.pop(id, None)
        else:
            self.overrides[id] = label
        self.rows.blockSignals(True)
        item.setData(STATE_COLUMN, CHOSEN, label)
        item.setText(STATE_COLUMN, self._choice_text(choice, label))
        self._embolden(item, id in self.overrides)
        self.rows.blockSignals(False)
        self._changed()
        return True

    def set_fixed(self, id: str, fixed: bool) -> bool:
        """Tick (fix) or untick (leave alone) the issue *id*."""
        item = self.line(id)
        issue = item.data(0, LINE) if item is not None else None
        if not isinstance(issue, Issue) or not issue.fixable:
            return False
        item.setCheckState(2, Qt.CheckState.Checked if fixed else Qt.CheckState.Unchecked)
        return True

    def _ticked(self, item: QTreeWidgetItem, column: int) -> None:
        issue = item.data(0, LINE)
        if column != 2 or not isinstance(issue, Issue) or not issue.fixable:
            return
        if item.checkState(2) == Qt.CheckState.Checked:
            self.skipped.discard(issue.id)
        else:
            self.skipped.add(issue.id)
        self._changed()

    def _changed(self) -> None:
        self.titleChanged.emit(self.title())
        self.decided.emit()

    def _filter(self) -> None:
        only = self.only_warnings.isChecked()
        for t in range(self.rows.topLevelItemCount()):
            top = self.rows.topLevelItem(t)
            for g in range(top.childCount()):
                group = top.child(g)
                shown = 0
                for i in range(group.childCount()):
                    line = group.child(i)
                    hide = only and line.text(0) != WARN
                    line.setHidden(hide)
                    shown += not hide
                group.setHidden(shown == 0)
                if only and shown:
                    group.setExpanded(True)

    def _selected(self) -> None:
        if self._select is None:
            return
        for item in self.rows.selectedItems():
            line = item.data(0, LINE)
            if isinstance(line, (Issue, Choice)) and line.residue >= 0:
                self._select(self.survey.protein, line.residue, line.label)
                return

    # -- acting -------------------------------------------------------------------------------
    def decisions(self) -> tuple[tuple[str, ...], dict[str, str]]:
        """What has been decided: the issues left alone, and the states chosen by hand."""
        return tuple(sorted(self.skipped)), dict(self.overrides)

    def reset(self) -> None:
        """Every line back as the survey proposed: the survey with nothing held."""
        self.fill(self.base)
        self._changed()

    def resurvey(self) -> None:
        if self._resurvey is not None:
            self._resurvey(*self.decisions())

    def apply(self) -> None:
        if self._apply is not None:
            self._apply(*self.decisions())

    def set_busy(self, busy: bool) -> None:
        """A survey or an apply is running: nothing more to start meanwhile."""
        for button in (self.reset_button, self.resurvey_button, self.apply_button):
            button.setEnabled(not busy)

    def state(self) -> dict:
        return {"protein": self.survey.protein}

    def close_panel(self) -> None:
        self._select = None
        self._resurvey = None
        self._apply = None
