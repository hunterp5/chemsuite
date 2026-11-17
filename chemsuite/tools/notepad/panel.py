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
"""The Notepad pane: a rich-text editor (fonts, sizes, bold, italic, underline, strike-through,
colors, highlighting, alignment, bullets, numbering, indents) over one of the document's
notes."""

from __future__ import annotations

import re
from collections.abc import Callable

from PySide6.QtCore import QEvent, QMimeData, QPointF, QRectF, QSize, Qt
from PySide6.QtGui import (
    QAction,
    QActionGroup,
    QColor,
    QFont,
    QIcon,
    QKeyEvent,
    QKeySequence,
    QPainter,
    QPen,
    QPixmap,
    QTextBlockFormat,
    QTextCharFormat,
    QTextCursor,
    QTextListFormat,
)
from PySide6.QtWidgets import (
    QColorDialog,
    QComboBox,
    QFontComboBox,
    QMessageBox,
    QTextEdit,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chemsuite.app.document.tools import ToolPanel

from .notes import Notes

TOOL_ID = "notepad"
SIZES = (8, 9, 10, 11, 12, 14, 16, 18, 20, 24, 28, 32, 36, 48, 60, 72)
_IMAGE = re.compile(r"<img\b[^>]*>", re.IGNORECASE)
_BODY = re.compile(r"<body\b[^>]*>")
_DEFAULT_FONT = re.compile(r"\s*font-(?:family|size):[^;]*;")
#: Ctrl+B, I and U format here before the window sees them (Ctrl+I is Invert Selection there)
_FORMAT_KEYS = {Qt.Key.Key_B: "bold", Qt.Key.Key_I: "italic", Qt.Key.Key_U: "underline"}


def portable(html: str) -> str:
    """*html* without the editor's default font in ``<body>`` (Segoe UI 9 pt on Windows): the
    note opens in each machine's own font. Fonts chosen for text are kept on the text."""
    return _BODY.sub(lambda m: _DEFAULT_FONT.sub("", m.group(0)), html, count=1)


class NoteEdit(QTextEdit):
    """The text: a paste keeps its formatting but not its pictures (documents stay small)."""

    format_key = None  # set by the panel: called with "bold", "italic" or "underline"

    def insertFromMimeData(self, source: QMimeData) -> None:
        if source.hasHtml():
            self.insertHtml(_IMAGE.sub("", source.html()))
        elif source.hasText():
            self.insertPlainText(source.text())

    def canInsertFromMimeData(self, source: QMimeData) -> bool:
        return source.hasHtml() or source.hasText()

    def event(self, event: QEvent) -> bool:
        if event.type() == QEvent.Type.ShortcutOverride and self._format_of(event):
            event.accept()  # ours, not the window's shortcut
            return True
        return super().event(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        name = self._format_of(event)
        if name and self.format_key is not None:
            self.format_key(name)
            return
        super().keyPressEvent(event)

    @staticmethod
    def _format_of(event) -> str | None:
        if not isinstance(event, QKeyEvent):
            return None
        if event.modifiers() != Qt.KeyboardModifier.ControlModifier:
            return None
        return _FORMAT_KEYS.get(Qt.Key(event.key()))


class NotepadPanel(ToolPanel):
    tool_id = TOOL_ID
    renamable = True

    def __init__(self, notes: Notes, note_id: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.notes = notes
        self.note_id = note_id
        self.editor = NoteEdit()
        self.editor.setAcceptRichText(True)
        self.editor.setFrameShape(QTextEdit.Shape.NoFrame)
        self.editor.document().setDocumentMargin(8)
        self.editor.format_key = lambda name: self.actions[name].trigger()
        self.editor.setHtml(notes.html(note_id))
        self.editor.document().setModified(False)
        notes.attach(note_id, lambda: portable(self.editor.toHtml()))
        self.editor.textChanged.connect(lambda: notes.touch(note_id))
        self.actions: dict[str, QAction] = {}
        self.toolbar_row = self._toolbar()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.toolbar_row)
        layout.addWidget(self.editor, 1)
        self.editor.currentCharFormatChanged.connect(self._show_format)
        self.editor.cursorPositionChanged.connect(self._show_block)
        self._show_format(self.editor.currentCharFormat())
        self._show_block()

    # -- ToolPanel -----------------------------------------------------------------------------
    def title(self) -> str:
        return self.notes.title(self.note_id)

    def state(self) -> dict:
        return {"note": self.note_id}

    def set_custom_title(self, text: str | None) -> None:
        self.notes.set_title(self.note_id, text)
        self.titleChanged.emit(self.title())

    def user_closing(self) -> bool:
        """Closing the pane deletes the note: asked first when it has text."""
        if self.editor.toPlainText().strip():
            answer = QMessageBox.question(
                self,
                "Notepad",
                f"Delete note “{self.title()}”?",
                QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            if answer != QMessageBox.StandardButton.Discard:
                return False
        self.notes.remove(self.note_id)
        return True

    def close_panel(self) -> None:
        self.notes.detach(self.note_id)  # the document is closing: the note stays

    # -- the toolbar ---------------------------------------------------------------------------
    def _toolbar(self) -> QToolBar:
        bar = QToolBar()
        bar.setIconSize(QSize(16, 16))
        bar.setMovable(False)
        bar.setFloatable(False)
        add = self._add

        add(bar, "undo", "Undo", self.editor.undo, keys=QKeySequence.StandardKey.Undo)
        add(bar, "redo", "Redo", self.editor.redo, keys=QKeySequence.StandardKey.Redo)
        self.editor.undoAvailable.connect(self.actions["undo"].setEnabled)
        self.editor.redoAvailable.connect(self.actions["redo"].setEnabled)
        self.actions["undo"].setEnabled(False)
        self.actions["redo"].setEnabled(False)
        bar.addSeparator()
        self.font_box = QFontComboBox()
        self.font_box.setToolTip("Font")
        self.font_box.setMaximumWidth(170)
        self.font_box.currentFontChanged.connect(self._set_family)
        bar.addWidget(self.font_box)
        self.size_box = QComboBox()
        self.size_box.setToolTip("Size")
        self.size_box.setEditable(True)
        self.size_box.addItems([str(s) for s in SIZES])
        self.size_box.setMinimumContentsLength(3)
        self.size_box.textActivated.connect(self._set_size)
        bar.addWidget(self.size_box)
        bar.addSeparator()
        add(bar, "bold", "Bold (Ctrl+B)", self._bold, checkable=True)
        add(bar, "italic", "Italic (Ctrl+I)", self._italic, checkable=True)
        add(bar, "underline", "Underline (Ctrl+U)", self._underline, checkable=True)
        add(bar, "strike", "Strikethrough", self._strike, checkable=True)
        bar.addSeparator()
        add(bar, "color", "Text Color…", self._color)
        add(bar, "highlight", "Highlight…", self._highlight)
        bar.addSeparator()
        group = QActionGroup(self)
        for name, title, flag in (
            ("left", "Align Left", Qt.AlignmentFlag.AlignLeft),
            ("center", "Center", Qt.AlignmentFlag.AlignHCenter),
            ("right", "Align Right", Qt.AlignmentFlag.AlignRight),
            ("justify", "Justify", Qt.AlignmentFlag.AlignJustify),
        ):
            action = add(bar, name, title, lambda f=flag: self._align(f), checkable=True)
            group.addAction(action)
        bar.addSeparator()
        add(bar, "bullets", "Bullets", lambda: self._list(QTextListFormat.Style.ListDisc),
            checkable=True)  # fmt: skip
        add(bar, "numbering", "Numbering", lambda: self._list(QTextListFormat.Style.ListDecimal),
            checkable=True)  # fmt: skip
        add(bar, "outdent", "Decrease Indent", lambda: self._indent(-1))
        add(bar, "indent", "Increase Indent", lambda: self._indent(1))
        bar.addSeparator()
        add(bar, "clear", "Clear Formatting", self._clear)
        return bar

    def _add(
        self,
        bar: QToolBar,
        name: str,
        tip: str,
        slot: Callable[[], None],
        *,
        checkable: bool = False,
        keys: QKeySequence.StandardKey | None = None,
    ) -> QAction:
        action = QAction(_icon(name, self.palette().text().color()), tip, self)
        action.setCheckable(checkable)
        action.triggered.connect(lambda _checked=False: slot())
        if keys is not None:
            action.setToolTip(f"{tip} ({QKeySequence(keys).toString(QKeySequence.NativeText)})")
        bar.addAction(action)
        self.actions[name] = action
        return action

    def changeEvent(self, event: QEvent) -> None:
        super().changeEvent(event)
        if event.type() == QEvent.Type.PaletteChange and self.actions:
            ink = self.palette().text().color()
            for name, action in self.actions.items():
                action.setIcon(_icon(name, ink))

    # -- character formats ---------------------------------------------------------------------
    def _merge(self, fmt: QTextCharFormat) -> None:
        """Apply to the selection, or to what is typed next."""
        self.editor.mergeCurrentCharFormat(fmt)
        self.editor.setFocus()

    def _bold(self) -> None:
        fmt = QTextCharFormat()
        bold = self.actions["bold"].isChecked()
        fmt.setFontWeight(QFont.Weight.Bold if bold else QFont.Weight.Normal)
        self._merge(fmt)

    def _italic(self) -> None:
        fmt = QTextCharFormat()
        fmt.setFontItalic(self.actions["italic"].isChecked())
        self._merge(fmt)

    def _underline(self) -> None:
        fmt = QTextCharFormat()
        fmt.setFontUnderline(self.actions["underline"].isChecked())
        self._merge(fmt)

    def _strike(self) -> None:
        fmt = QTextCharFormat()
        fmt.setFontStrikeOut(self.actions["strike"].isChecked())
        self._merge(fmt)

    def _set_family(self, font: QFont) -> None:
        fmt = QTextCharFormat()
        fmt.setFontFamilies([font.family()])
        self._merge(fmt)

    def _set_size(self, text: str) -> None:
        try:
            size = float(text)
        except ValueError:
            return
        if 1 <= size <= 400:
            fmt = QTextCharFormat()
            fmt.setFontPointSize(size)
            self._merge(fmt)

    def _color(self) -> None:
        current = self.editor.currentCharFormat().foreground().color()
        color = QColorDialog.getColor(current, self, "Text Color")
        if color.isValid():
            fmt = QTextCharFormat()
            fmt.setForeground(color)
            self._merge(fmt)

    def _highlight(self) -> None:
        current = self.editor.currentCharFormat().background().color()
        start = current if current.alpha() else QColor("yellow")
        color = QColorDialog.getColor(start, self, "Highlight")
        if color.isValid():
            fmt = QTextCharFormat()
            fmt.setBackground(color)
            self._merge(fmt)

    def _clear(self) -> None:
        """Plain text again: no character formats, no list, no indent, left aligned."""
        cursor = self.editor.textCursor()
        cursor.beginEditBlock()
        if cursor.hasSelection():
            cursor.setCharFormat(QTextCharFormat())
        self.editor.setCurrentCharFormat(QTextCharFormat())
        for block_cursor in self._block_cursors():
            current = block_cursor.currentList()
            if current is not None:
                current.remove(block_cursor.block())
            block_cursor.setBlockFormat(QTextBlockFormat())
        cursor.endEditBlock()
        self.editor.setFocus()

    # -- paragraphs ----------------------------------------------------------------------------
    def _block_cursors(self) -> list[QTextCursor]:
        """A cursor at each paragraph the selection touches (the cursor's own without one)."""
        cursor = self.editor.textCursor()
        document = self.editor.document()
        first = document.findBlock(cursor.selectionStart())
        last = document.findBlock(cursor.selectionEnd())
        found, block = [], first
        while block.isValid():
            found.append(QTextCursor(block))
            if block == last:
                break
            block = block.next()
        return found

    def _align(self, flag: Qt.AlignmentFlag) -> None:
        self.editor.setAlignment(flag)
        self.editor.setFocus()

    def _list(self, style: QTextListFormat.Style) -> None:
        """Make the paragraphs a list of *style*, change their list to it, or (already one)
        take them out of the list."""
        cursor = self.editor.textCursor()
        current = cursor.currentList()
        cursor.beginEditBlock()
        if current is not None and current.format().style() == style:
            for block_cursor in self._block_cursors():
                in_list = block_cursor.currentList()
                if in_list is not None:
                    in_list.remove(block_cursor.block())
                fmt = block_cursor.blockFormat()
                fmt.setIndent(0)
                block_cursor.setBlockFormat(fmt)
        elif current is not None:
            fmt = current.format()
            fmt.setStyle(style)
            current.setFormat(fmt)
        else:
            fmt = QTextListFormat()
            fmt.setStyle(style)
            fmt.setIndent(cursor.blockFormat().indent() + 1)
            cursor.createList(fmt)
        cursor.endEditBlock()
        self._show_block()
        self.editor.setFocus()

    def _indent(self, step: int) -> None:
        """Nest a list item a level deeper (or out), else indent the paragraph."""
        cursor = self.editor.textCursor()
        current = cursor.currentList()
        cursor.beginEditBlock()
        if current is not None:
            fmt = QTextListFormat(current.format())
            level = fmt.indent() + step
            if level < 1:
                for block_cursor in self._block_cursors():
                    in_list = block_cursor.currentList()
                    if in_list is not None:
                        in_list.remove(block_cursor.block())
                    block_format = block_cursor.blockFormat()
                    block_format.setIndent(0)
                    block_cursor.setBlockFormat(block_format)
            else:
                fmt.setIndent(level)
                cursor.createList(fmt)
        else:
            for block_cursor in self._block_cursors():
                block_format = block_cursor.blockFormat()
                block_format.setIndent(max(0, block_format.indent() + step))
                block_cursor.setBlockFormat(block_format)
        cursor.endEditBlock()
        self._show_block()
        self.editor.setFocus()

    # -- following the cursor ------------------------------------------------------------------
    def _show_format(self, fmt: QTextCharFormat) -> None:
        font = fmt.font()
        for name, on in (
            ("bold", font.weight() > QFont.Weight.Normal),
            ("italic", font.italic()),
            ("underline", font.underline()),
            ("strike", font.strikeOut()),
        ):
            self.actions[name].setChecked(on)
        self.font_box.blockSignals(True)
        self.font_box.setCurrentFont(font)
        self.font_box.blockSignals(False)
        size = fmt.fontPointSize() or self.editor.font().pointSizeF()
        self.size_box.setEditText(f"{size:g}")

    def _show_block(self) -> None:
        alignment = self.editor.alignment()
        if alignment & Qt.AlignmentFlag.AlignHCenter:
            name = "center"
        elif alignment & Qt.AlignmentFlag.AlignJustify:
            name = "justify"
        elif alignment & Qt.AlignmentFlag.AlignRight:
            name = "right"
        else:
            name = "left"
        self.actions[name].setChecked(True)
        current = self.editor.textCursor().currentList()
        style = current.format().style() if current is not None else None
        self.actions["bullets"].setChecked(style == QTextListFormat.Style.ListDisc)
        self.actions["numbering"].setChecked(style == QTextListFormat.Style.ListDecimal)


# -- icons -----------------------------------------------------------------------------------------
def _icon(name: str, ink: QColor) -> QIcon:
    """A 16 px toolbar icon drawn in *ink* (sharp at any screen scale)."""
    icon = QIcon()
    for scale in (1, 2):
        pixmap = QPixmap(16 * scale, 16 * scale)
        pixmap.setDevicePixelRatio(scale)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(ink, 1.4, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
        _DRAW[name](painter, ink)
        painter.end()
        icon.addPixmap(pixmap)
    return icon


def _letter(p: QPainter, text: str, *, bold=False, italic=False, size=11.0, y=13.0) -> None:
    font = QFont()
    font.setPointSizeF(size)
    font.setBold(bold)
    font.setItalic(italic)
    p.setFont(font)
    p.drawText(QRectF(0, y - 13, 16, 14), Qt.AlignmentFlag.AlignHCenter, text)


def _lines(p: QPainter, spans: list[tuple[float, float]], top: float = 3.5) -> None:
    for k, (a, b) in enumerate(spans):
        y = top + 3 * k
        p.drawLine(QPointF(a, y), QPointF(b, y))


def _bar(p: QPainter, color: QColor) -> None:
    p.fillRect(QRectF(2, 13, 12, 2.5), color)


def _arrow(p: QPainter, inward: bool) -> None:
    _lines(p, [(7, 14), (7, 14), (7, 14), (2, 14)], top=3.5)
    x0, x1 = (1.5, 5) if inward else (5, 1.5)
    p.drawLine(QPointF(x0, 7), QPointF(x1, 7))
    tip = QPointF(x1, 7)
    back = -1.8 if inward else 1.8
    p.drawLine(tip, QPointF(x1 + back, 5.2))
    p.drawLine(tip, QPointF(x1 + back, 8.8))


def _list_icon(p: QPainter, ink: QColor, numbered: bool) -> None:
    for k in range(3):
        y = 3.5 + 4.5 * k
        if numbered:
            font = QFont()
            font.setPointSizeF(5.5)
            p.setFont(font)
            p.drawText(QRectF(0, y - 3.5, 5, 7), Qt.AlignmentFlag.AlignCenter, str(k + 1))
        else:
            p.setBrush(ink)
            p.drawEllipse(QPointF(2.5, y), 1.1, 1.1)
            p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawLine(QPointF(6.5, y), QPointF(14.5, y))


def _undo_icon(p: QPainter, redo: bool) -> None:
    p.save()
    if redo:
        p.translate(16, 0)
        p.scale(-1, 1)
    p.drawArc(QRectF(3, 5, 10, 8), 90 * 16, -200 * 16)
    p.drawLine(QPointF(8, 5), QPointF(5, 3))
    p.drawLine(QPointF(8, 5), QPointF(5.5, 7.5))
    p.restore()


_DRAW: dict[str, Callable[[QPainter, QColor], None]] = {
    "undo": lambda p, ink: _undo_icon(p, False),
    "redo": lambda p, ink: _undo_icon(p, True),
    "bold": lambda p, ink: _letter(p, "B", bold=True),
    "italic": lambda p, ink: _letter(p, "I", italic=True),
    "underline": lambda p, ink: (
        _letter(p, "U", y=12),
        p.drawLine(QPointF(4, 14.5), QPointF(12, 14.5)),
    ),
    "strike": lambda p, ink: (_letter(p, "S"), p.drawLine(QPointF(3, 8), QPointF(13, 8))),
    "color": lambda p, ink: (_letter(p, "A", y=11.5, size=10), _bar(p, QColor("#d03030"))),
    "highlight": lambda p, ink: (_letter(p, "ab", y=11.5, size=8), _bar(p, QColor("#f2d230"))),
    "left": lambda p, ink: _lines(p, [(2, 14), (2, 10), (2, 14), (2, 10)]),
    "center": lambda p, ink: _lines(p, [(2, 14), (4, 12), (2, 14), (4, 12)]),
    "right": lambda p, ink: _lines(p, [(2, 14), (6, 14), (2, 14), (6, 14)]),
    "justify": lambda p, ink: _lines(p, [(2, 14), (2, 14), (2, 14), (2, 14)]),
    "bullets": lambda p, ink: _list_icon(p, ink, False),
    "numbering": lambda p, ink: _list_icon(p, ink, True),
    "indent": lambda p, ink: _arrow(p, True),
    "outdent": lambda p, ink: _arrow(p, False),
    "clear": lambda p, ink: (_letter(p, "T", size=10), p.drawLine(QPointF(3, 14), QPointF(13, 3))),
}
