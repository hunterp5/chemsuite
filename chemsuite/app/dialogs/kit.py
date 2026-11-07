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

* **sections** in a fixed order, *Input* (structures or columns, which rows), *Method*,
  *Output* (column names, what to do if they exist), each a :class:`Section` with one label
  column, fields of one width, and explanations in tooltips on the field and its label;
* a **footer**: Guide (bottom left, opening the tool's section of the user guide), Cancel and
  the primary button, which names the action ("Cluster 1,204 Rows"). Above them, only when
  something must be fixed first, the problem with a warning mark. There is no header (the
  window title names the tool) and no explanatory text: no sentence describing what will
  happen, no heads-up notes. Explanations go in tooltips.

Shared controls keep tools alike: :class:`RowsChoice` (All · Shown · Selected, with counts) and
:class:`OutputChoice` (If it exists: Replace it · Add a new column).

A subclass builds its sections in ``__init__``, calls :meth:`ToolDialog.refresh` when anything
changes, and overrides :meth:`~ToolDialog.problem`, :meth:`~ToolDialog.verb` and
:meth:`~ToolDialog.run`.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path

import numpy as np
import shiboken6
from PySide6.QtCore import Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFont, QPalette
from PySide6.QtWidgets import (
    QAbstractSpinBox,
    QButtonGroup,
    QCheckBox,
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
    QPushButton,
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


class Section(QWidget):
    """A titled group of fields: a small capitalised heading and rule over a form.

    Explanations are tooltips on the field and its label, not text in the dialog, so the
    dialog stays compact.
    """

    def __init__(self, title: str, *, collapsed: bool | None = None) -> None:
        """*collapsed* True or False makes the section foldable (folded or open to start);
        a foldable section holds optional settings whose defaults suit most runs."""
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
        #: the fold button of a foldable section (None otherwise)
        self.toggle: QToolButton | None = None
        if collapsed is not None:
            self.toggle = QToolButton()
            self.toggle.setAutoRaise(True)
            self.toggle.setCheckable(True)
            self.toggle.setChecked(not collapsed)
            self.toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
            self.toggle.setStyleSheet("QToolButton { border: none; padding: 0; }")
            self.toggle.setToolTip("Show or hide these settings")
            self.toggle.toggled.connect(self.set_open)
            heading.addWidget(self.toggle)
            self.title.setCursor(Qt.CursorShape.PointingHandCursor)
            self.title.mousePressEvent = lambda _event: self.toggle.toggle()
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
        self.body = QWidget()
        self.body.setLayout(self.form)
        self.form.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.body)
        self._captions: dict[QWidget, QLabel] = {}
        if self.toggle is not None:
            self.set_open(self.toggle.isChecked())

    def set_open(self, on: bool) -> None:
        """Unfold or fold a foldable section."""
        if self.toggle is None:
            return
        if self.toggle.isChecked() != on:
            self.toggle.setChecked(on)  # calls back here
            return
        self.toggle.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        self.body.setVisible(on)

    def is_open(self) -> bool:
        return self.toggle is None or self.toggle.isChecked()

    def row(self, label: str, field: QWidget, help: str = "") -> QWidget:
        """Add *field* with *label* ("" for none); *help* becomes its tooltip. Returns *field*."""
        caption = QLabel(f"{label}:" if label else "")
        caption.setMinimumWidth(LABEL_WIDTH)
        caption.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if isinstance(field, _FIXED_WIDTH_FIELDS):
            field.setFixedWidth(FIELD_WIDTH)
        self.form.addRow(caption, field)
        self._captions[field] = caption
        self.set_help(field, help)
        return field

    def set_help(self, field: QWidget, text: str) -> None:
        """The explanation shown when hovering *field* or its label ("" for none)."""
        field.setToolTip(text)
        self._captions[field].setToolTip(text)

    def help(self, field: QWidget) -> str:
        return field.toolTip()

    def set_row_visible(self, field: QWidget, visible: bool) -> None:
        self.form.setRowVisible(field, visible)

    def is_row_visible(self, field: QWidget) -> bool:
        return self.form.isRowVisible(field)


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


def follow_rows(widget: QWidget, ctx, update: Callable[[], None]) -> None:
    """Call *update* now and whenever the selection or table changes, while *widget* lives.

    Tool dialogs do not block the window, so anything showing row counts must keep up.
    """
    unsubscribe: list[Callable[[], None]] = []

    def stop() -> None:
        while unsubscribe:
            unsubscribe.pop()()

    def changed(*_args) -> None:
        if shiboken6.isValid(widget):
            update()
        else:
            stop()

    unsubscribe += [ctx.selection.subscribe(changed), ctx.table.subscribe(changed)]
    widget.destroyed.connect(lambda *_a: stop())
    update()


def selected_only(ctx, text: str = "Selected rows only") -> QCheckBox:
    """A "Selected rows only (N)" check box whose count follows the selection.

    It is unticked and disabled while nothing is selected.
    """
    box = QCheckBox()

    def update() -> None:
        count = len(ctx.selection)
        box.setText(f"{text} ({count:,})")
        box.setEnabled(count > 0)
        if not count:
            box.setChecked(False)

    follow_rows(box, ctx, update)
    return box


class RowsChoice(Segmented):
    """Which rows a tool works on: All · Shown · Selected, with counts that stay current.

    *ctx* is the tool's ``ToolContext``. Dialogs do not block the window, so the counts follow
    the selection and filters while the dialog is open. Choosing Selected with nothing selected
    is allowed (you may be about to select rows); :meth:`problem` then says so.
    """

    def __init__(self, ctx, current: str = "all") -> None:
        super().__init__(
            (("all", "All"), ("shown", "Shown"), ("selected", "Selected")),
            current,
        )
        self.ctx = ctx
        self._unsubscribe = [
            ctx.selection.subscribe(lambda *_a: self._follow()),
            ctx.table.subscribe(lambda *_a: self._follow()),
        ]
        self.update_counts()

    def _follow(self) -> None:
        if not shiboken6.isValid(self):  # the dialog closed with its window
            self.detach()
            return
        self.update_counts()
        self.changed.emit(self.value())

    def detach(self) -> None:
        """Stop following the table and selection (the dialog is closing)."""
        for unsubscribe in self._unsubscribe:
            unsubscribe()
        self._unsubscribe = []

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

    def problem(self) -> str:
        if self.value() == "selected" and not len(self.ctx.selection):
            return "No rows are selected. Select rows in the table, or choose All or Shown."
        return ""

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
    """A tool's dialog: sections and footer (see the module notes)."""

    def __init__(self, parent: QWidget | None, title: str, *, guide: str = "") -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.job = None
        self.guide_anchor = guide

        self.sections = QVBoxLayout()
        self.sections.setSpacing(6)

        rule = QFrame()
        rule.setFrameShape(QFrame.Shape.HLine)
        rule.setFrameShadow(QFrame.Shadow.Plain)
        rule.setEnabled(False)
        self.warning_mark = QLabel("⚠")
        self.warning_mark.setStyleSheet("color: #c77700; font-size: 14px;")
        self.warning_mark.setAlignment(Qt.AlignmentFlag.AlignTop)
        #: what must be fixed before running (hidden when nothing is wrong)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        self.status.setStyleSheet("color: #b36b00;")
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
        self.guide_button = QPushButton("Guide")
        self.guide_button.setToolTip("Open this tool's section of the user guide")
        self.guide_button.setAutoDefault(False)
        self.guide_button.clicked.connect(self.open_guide)
        self.guide_button.setVisible(GUIDE.exists())
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.footer = QWidget()
        footer = QHBoxLayout(self.footer)
        footer.setContentsMargins(0, 0, 0, 0)
        footer.setSpacing(8)
        footer.addWidget(self.warning_mark)
        footer.addWidget(self.status, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 16, 18, 14)
        layout.setSpacing(6)
        layout.addLayout(self.sections)
        layout.addSpacing(4)
        layout.addWidget(rule)
        layout.addWidget(self.footer)
        bottom = QHBoxLayout()
        bottom.addWidget(self.guide_button)
        bottom.addStretch(1)
        bottom.addWidget(self.buttons)
        layout.addLayout(bottom)
        # The dialog is always exactly as big as its contents: rows that appear or disappear
        # (another method's settings) resize it rather than leave gaps.
        layout.setSizeConstraint(QLayout.SizeConstraint.SetFixedSize)
        self.status.setFixedWidth(LABEL_WIDTH + 10 + FIELD_WIDTH - 20)
        self.footer.hide()

    # -- building ---------------------------------------------------------------------------
    def section(self, title: str, *, collapsed: bool | None = None) -> Section:
        """Add a section; ``collapsed=True`` makes a folded "More Options"-style section."""
        section = Section(title, collapsed=collapsed)
        self.sections.addWidget(section)
        return section

    # -- what subclasses say ----------------------------------------------------------------
    def problem(self) -> str:
        """What must be fixed before running ("" when it can run)."""
        return ""

    def verb(self) -> str:
        """The primary button's text: the action and its size ("Cluster 1,204 Rows")."""
        return "OK"

    def run(self):
        """Start the work; return the job (or None)."""
        return None

    # -- behaviour --------------------------------------------------------------------------
    def _problem(self) -> str:
        """The kit's own checks (rows to work on), then the tool's."""
        for choice in self.findChildren(RowsChoice):
            if choice.isVisibleTo(self) and choice.problem():
                return choice.problem()
        return self.problem()

    def refresh(self) -> None:
        """Bring the footer and the primary button up to date with the fields."""
        problem = self._problem()
        self.status.setText(problem)
        self.footer.setVisible(bool(problem))
        self.primary.setText(self.verb())
        self.primary.setEnabled(not problem)

    def accept(self) -> None:
        if self._problem():
            self.refresh()
            return
        self.job = self.run()
        super().accept()

    def done(self, result: int) -> None:
        for choice in self.findChildren(RowsChoice):
            choice.detach()
        super().done(result)

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
    "follow_rows",
    "selected_only",
]
