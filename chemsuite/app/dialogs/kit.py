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
"""The pieces every tool dialog is built from (the "sectioned form" design).

A tool dialog is a :class:`ToolDialog`:

* a **header**: the tool's name, one sentence on what it does, and a Guide link;
* **sections** in a fixed order, *Input* (structures or columns, which rows), *Method*,
  *Output* (column names, what to do if they exist), each a :class:`Section` with one label
  column, fields of one width, and help text under the field it explains;
* a **footer**: a sentence saying what will happen, which turns into the problem (with a
  warning mark) when something must be fixed first, then Cancel and the primary button,
  which names the action ("Cluster 1,204 Rows").

Shared controls keep tools alike: :class:`RowsChoice` (All · Shown · Selected, with counts) and
:class:`OutputChoice` (If it exists: Replace it · Add a new column).

A subclass builds its sections in ``__init__``, calls :meth:`ToolDialog.refresh` when anything
changes, and overrides :meth:`~ToolDialog.problem`, :meth:`~ToolDialog.summary`,
:meth:`~ToolDialog.note`, :meth:`~ToolDialog.verb` and :meth:`~ToolDialog.run`.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from PySide6.QtCore import QEvent, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFont, QPalette
from PySide6.QtWidgets import (
    QAbstractSpinBox,
    QButtonGroup,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLayout,
    QLineEdit,
    QListWidget,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

#: every field (menu, number box, text box) has this width, so the right edges line up
FIELD_WIDTH = 300
#: the label column's minimum width, the same in every section of a dialog
LABEL_WIDTH = 116
#: the user guide, when it ships beside the code (a source checkout or an unpacked install)
GUIDE = Path(__file__).resolve().parents[3] / "docs" / "USER_GUIDE.md"

_FIXED_WIDTH_FIELDS = (QComboBox, QLineEdit, QAbstractSpinBox, QListWidget)


def muted(label: QLabel) -> None:
    """Grey text that stays readable in every theme: the palette's placeholder colour."""
    label.setForegroundRole(QPalette.ColorRole.PlaceholderText)


def help_label(text: str = "") -> QLabel:
    """Small grey text explaining the field above it."""
    label = QLabel(text)
    label.setWordWrap(True)
    label.setTextFormat(Qt.TextFormat.PlainText)
    font = label.font()
    font.setPointSizeF(max(8.0, font.pointSizeF() - 1))
    label.setFont(font)
    muted(label)
    label.setFixedWidth(FIELD_WIDTH)
    return label


class Section(QWidget):
    """A titled group of fields: a small capitalised heading and rule over a form."""

    def __init__(self, title: str) -> None:
        super().__init__()
        self.title = QLabel(title.upper())
        font = self.title.font()
        font.setPointSizeF(max(8.0, font.pointSizeF() - 1.5))
        font.setWeight(QFont.Weight.DemiBold)
        font.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, 0.8)
        self.title.setFont(font)
        muted(self.title)
        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setFrameShadow(QFrame.Shadow.Plain)
        rule.setEnabled(False)
        heading = QHBoxLayout()
        heading.setSpacing(8)
        heading.addWidget(self.title)
        heading.addWidget(rule, 1)
        self.form = QFormLayout()
        self.form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.form.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        self.form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
        self.form.setHorizontalSpacing(10)
        self.form.setVerticalSpacing(7)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 6, 0, 2)
        layout.setSpacing(8)
        layout.addLayout(heading)
        layout.addLayout(self.form)
        self._rows: dict[QWidget, QWidget] = {}  # field → what the form holds (field or box)
        self._helps: dict[QWidget, QLabel] = {}

    def row(self, label: str, field: QWidget, help: str = "") -> QWidget:
        """Add *field* with *label* ("" for none) and optional help text; returns *field*."""
        caption = QLabel(f"{label}:" if label else "")
        caption.setMinimumWidth(LABEL_WIDTH)
        caption.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if isinstance(field, _FIXED_WIDTH_FIELDS):
            field.setFixedWidth(FIELD_WIDTH)
        box = QWidget()
        column = QVBoxLayout(box)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(3)
        column.addWidget(field)
        note = help_label()
        column.addWidget(note)
        self.form.addRow(caption, box)
        self._rows[field] = box
        self._helps[field] = note
        self.set_help(field, help)
        return field

    def set_help(self, field: QWidget, text: str) -> None:
        """Change (or with "" hide) the help under *field*."""
        note = self._helps[field]
        note.setText(text)
        note.setVisible(bool(text))
        self._fit(field)

    def _fit(self, field: QWidget) -> None:
        # A word-wrapped label inside a form row: give the row its real height up front, or the
        # form squeezes it (Qt sizes rows before it knows how many lines the help needs).
        note = self._helps[field]
        height = field.sizeHint().height()
        if note.isVisibleTo(self) or note.text():
            height += 3 + note.heightForWidth(FIELD_WIDTH)
        self._rows[field].setFixedHeight(max(height, field.minimumSizeHint().height()))

    def set_row_visible(self, field: QWidget, visible: bool) -> None:
        self.form.setRowVisible(self._rows[field], visible)

    def is_row_visible(self, field: QWidget) -> bool:
        return self.form.isRowVisible(self._rows[field])

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.FontChange:
            for field in self._rows:
                self._fit(field)


class Segmented(QWidget):
    """A few joined toggle buttons, one of them chosen (e.g. All · Shown · Selected)."""

    changed = Signal(str)

    def __init__(self, choices: Sequence[tuple[str, str]], current: str | None = None) -> None:
        super().__init__()
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons: dict[str, QToolButton] = {}
        count = len(choices)
        for i, (key, text) in enumerate(choices):
            button = QToolButton()
            button.setText(text)
            button.setCheckable(True)
            button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            button.setStyleSheet(_segment_style(i, count))
            self.group.addButton(button, i)
            layout.addWidget(button)
            self.buttons[key] = button
            button.toggled.connect(lambda on, k=key: on and self.changed.emit(k))
        self.set_value(current if current in self.buttons else choices[0][0])
        self.setFixedWidth(FIELD_WIDTH)

    def value(self) -> str:
        return next(k for k, b in self.buttons.items() if b.isChecked())

    def set_value(self, key: str) -> None:
        self.buttons[key].setChecked(True)

    def set_enabled(self, key: str, enabled: bool, why: str = "") -> None:
        """Enable or disable one choice; a disabled chosen one hands over to the first."""
        button = self.buttons[key]
        button.setEnabled(enabled)
        button.setToolTip("" if enabled else why)
        if not enabled and button.isChecked():
            first = next((k for k, b in self.buttons.items() if b.isEnabled()), None)
            if first is not None:
                self.set_value(first)

    def set_text(self, key: str, text: str) -> None:
        self.buttons[key].setText(text)


def _segment_style(index: int, count: int) -> str:
    left = 5 if index == 0 else 0
    right = 5 if index == count - 1 else 0
    return (
        "QToolButton { border: 1px solid palette(mid); padding: 3px 8px;"
        f" border-top-left-radius: {left}px; border-bottom-left-radius: {left}px;"
        f" border-top-right-radius: {right}px; border-bottom-right-radius: {right}px;"
        f"{' border-left: none;' if index else ''} background: palette(button); }}"
        " QToolButton:checked { background: palette(highlight);"
        " color: palette(highlighted-text); border-color: palette(highlight); }"
        " QToolButton:disabled { color: palette(mid); }"
    )


class RowsChoice(Segmented):
    """Which rows a tool works on: All · Shown · Selected, with live counts.

    *ctx* is the tool's ``ToolContext``. Selected is greyed out while nothing is selected.
    """

    def __init__(self, ctx, current: str = "all") -> None:
        super().__init__(
            (("all", "All"), ("shown", "Shown"), ("selected", "Selected")),
            current,
        )
        self.ctx = ctx
        self.update_counts()

    def update_counts(self) -> None:
        from chemsuite.plotting.data import shown_oids

        table = self.ctx.table
        counts = {
            "all": len(table),
            "shown": len(shown_oids(table)),
            "selected": len(self.ctx.selection),
        }
        for key, title in (("all", "All"), ("shown", "Shown"), ("selected", "Selected")):
            self.set_text(key, f"{title} ({counts[key]:,})")
        self.set_enabled("selected", counts["selected"] > 0, "Select rows in the table first")

    def oids(self) -> np.ndarray:
        from chemsuite.tools.similarity.common import scope_oids

        return scope_oids(self.ctx, self.value())

    def describe(self) -> str:
        return {"all": "", "shown": "shown ", "selected": "selected "}[self.value()]


class OutputChoice(Segmented):
    """If the result column already exists: replace it, or add a new column beside it."""

    def __init__(self) -> None:
        super().__init__((("replace", "Replace it"), ("new", "Add a new column")), "replace")

    def overwrite(self) -> bool:
        return self.value() == "replace"


class ToolDialog(QDialog):
    """A tool's dialog: header, sections, footer (see the module notes)."""

    def __init__(
        self,
        parent: QWidget | None,
        title: str,
        blurb: str,
        *,
        guide: str = "",
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.job = None
        self.guide_anchor = guide

        heading = QLabel(title)
        font = heading.font()
        font.setPointSizeF(font.pointSizeF() + 3)
        font.setWeight(QFont.Weight.DemiBold)
        heading.setFont(font)
        self.blurb = QLabel(blurb)
        self.blurb.setWordWrap(True)
        self.guide_link = QLabel('<a href="guide">Guide</a>')
        self.guide_link.setToolTip("Open this tool's section of the user guide")
        self.guide_link.linkActivated.connect(lambda _href: self.open_guide())
        self.guide_link.setVisible(GUIDE.exists())
        top = QHBoxLayout()
        top.addWidget(heading)
        top.addStretch(1)
        top.addWidget(self.guide_link)

        self.sections = QVBoxLayout()
        self.sections.setSpacing(6)

        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setFrameShadow(QFrame.Shadow.Plain)
        rule.setEnabled(False)
        self.warning_mark = QLabel("⚠")
        self.warning_mark.setStyleSheet("color: #c77700; font-size: 14px;")
        self.warning_mark.setAlignment(Qt.AlignmentFlag.AlignTop)
        #: the footer sentence: the problem when there is one, else what will happen
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        #: a heads-up that does not block running ("Butina on 80,000 rows may take minutes")
        self.note_label = QLabel()
        self.note_label.setWordWrap(True)
        muted(self.note_label)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.primary = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        self.primary.setStyleSheet(
            "QPushButton { background: palette(highlight); color: palette(highlighted-text);"
            " border: 1px solid palette(highlight); border-radius: 5px; padding: 4px 14px;"
            " font-weight: 600; }"
            " QPushButton:disabled { background: palette(button); color: palette(mid);"
            " border-color: palette(mid); }"
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        words = QVBoxLayout()
        words.setSpacing(2)
        words.addWidget(self.status)
        words.addWidget(self.note_label)
        footer = QHBoxLayout()
        footer.setSpacing(8)
        footer.addWidget(self.warning_mark)
        footer.addLayout(words, 1)
        footer.addWidget(self.buttons, 0, Qt.AlignmentFlag.AlignBottom)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 14)
        layout.setSpacing(6)
        layout.addLayout(top)
        layout.addWidget(self.blurb)
        layout.addSpacing(4)
        layout.addLayout(self.sections)
        layout.addSpacing(4)
        layout.addWidget(rule)
        layout.addLayout(footer)
        # The dialog is always exactly as big as its contents: rows that appear or disappear
        # (another method's settings) resize it rather than leave gaps.
        layout.setSizeConstraint(QLayout.SizeConstraint.SetFixedSize)
        self.blurb.setFixedWidth(LABEL_WIDTH + 10 + FIELD_WIDTH)

    # -- building ---------------------------------------------------------------------------
    def section(self, title: str) -> Section:
        section = Section(title)
        self.sections.addWidget(section)
        return section

    # -- what subclasses say ----------------------------------------------------------------
    def problem(self) -> str:
        """What must be fixed before running ("" when it can run)."""
        return ""

    def summary(self) -> str:
        """One sentence: what running will do."""
        return ""

    def note(self) -> str:
        """A heads-up that does not stop the run ("this may take minutes"), shown greyed."""
        return ""

    def verb(self) -> str:
        """The primary button's text: the action and its size ("Cluster 1,204 Rows")."""
        return "OK"

    def run(self):
        """Start the work; return the job (or None)."""
        return None

    # -- behaviour --------------------------------------------------------------------------
    def refresh(self) -> None:
        """Bring the footer and the primary button up to date with the fields."""
        problem = self.problem()
        self.warning_mark.setVisible(bool(problem))
        self.status.setText(problem or self.summary())
        self.status.setStyleSheet("color: #b36b00;" if problem else "")
        heads_up = "" if problem else self.note()
        self.note_label.setText(heads_up)
        self.note_label.setVisible(bool(heads_up))
        self.primary.setText(self.verb())
        self.primary.setEnabled(not problem)

    def accept(self) -> None:
        if self.problem():
            self.refresh()
            return
        self.job = self.run()
        super().accept()

    def open_guide(self) -> None:
        url = QUrl.fromLocalFile(str(GUIDE))
        if self.guide_anchor:
            url.setFragment(self.guide_anchor)
        QDesktopServices.openUrl(url)


__all__ = [
    "FIELD_WIDTH",
    "LABEL_WIDTH",
    "OutputChoice",
    "RowsChoice",
    "Section",
    "Segmented",
    "ToolDialog",
    "help_label",
]
