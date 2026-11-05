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

"""The structure editor: type or paste SMILES, a molfile or InChI, with a live drawing.

It is deliberately text-first (no drawing canvas): chemists draw in ChemDraw or similar and
paste here. The clipboard helpers read molfiles in the formats drawing programs put there.
"""

from __future__ import annotations

from PySide6.QtCore import QMimeData, Qt, QTimer
from PySide6.QtGui import QFontDatabase, QGuiApplication, QPixmap
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

from chemsuite.core.depict import render_png
from chemsuite.core.mols import looks_like_molblock, smiles_of
from chemsuite.core.structure_input import StructureInput, molblock_of, read_structure

#: Clipboard formats that carry a molfile (ChemDraw, MarvinSketch, Biovia Draw and others).
MOLFILE_FORMATS = (
    "chemical/x-mdl-molfile",
    "chemical/x-mdl-molfile-v3000",
    "application/x-mdl-molfile",
    "MDLCT",
)
PARSE_DELAY_MS = 150


def _decode(data: bytes) -> str:
    return data.replace(b"\x00", b"").decode("utf-8", errors="replace")


def _mdlct_to_text(data: bytes) -> str:
    """ChemDraw's Windows "MDLCT" format: each line prefixed by its length byte."""
    lines, i = [], 0
    while i < len(data):
        n = data[i]
        lines.append(data[i + 1 : i + 1 + n].decode("latin-1"))
        i += 1 + n
    return "\n".join(lines)


def structure_text_from_mime(mime: QMimeData | None) -> str:
    """The structure on the clipboard: a molfile if one is offered, else the plain text."""
    if mime is None:
        return ""
    for fmt in MOLFILE_FORMATS:
        if mime.hasFormat(fmt):
            data = bytes(mime.data(fmt).data())
            text = _mdlct_to_text(data) if fmt == "MDLCT" else _decode(data)
            if looks_like_molblock(text):
                return text
    return mime.text() if mime.hasText() else ""


def clipboard_structure_text() -> str:
    return structure_text_from_mime(QGuiApplication.clipboard().mimeData())


def copy_structure(mol: Chem.Mol | None) -> bool:
    """Put a structure on the clipboard as a molfile (for drawing programs) and SMILES text."""
    if mol is None:
        return False
    mime = QMimeData()
    mime.setData("chemical/x-mdl-molfile", molblock_of(mol).encode())
    mime.setText(smiles_of(mol))
    QGuiApplication.clipboard().setMimeData(mime)
    return True


class StructureEditor(QDialog):
    """Edit one structure. :attr:`result_input` holds the accepted structure."""

    def __init__(
        self, text: str = "", title: str = "Edit Structure", parent: QWidget | None = None
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.original = text
        self.edit = QPlainTextEdit()
        self.edit.setPlaceholderText("SMILES, molfile or InChI. Paste from ChemDraw with Paste.")
        self.edit.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.edit.setTabChangesFocus(True)
        self.edit.setMaximumHeight(110)
        self.preview = QLabel()
        self.preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.preview.setMinimumSize(420, 280)
        self.preview.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.preview.setStyleSheet("background: palette(base); border: 1px solid palette(mid)")
        self.info = QLabel(" ")
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self.paste_button = QPushButton("Paste")
        self.paste_button.setToolTip("Replace with the structure on the clipboard")
        self.clean_button = QPushButton("Canonical SMILES")
        self.clean_button.setToolTip("Rewrite as canonical SMILES (drops molfile coordinates)")
        self.copy_button = QPushButton("Copy")
        self.copy_button.setToolTip("Copy as molfile and SMILES")
        self.clear_button = QPushButton("Clear")
        tools = QHBoxLayout()
        for button in (self.paste_button, self.clean_button, self.copy_button, self.clear_button):
            tools.addWidget(button)
        tools.addStretch(1)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )

        layout = QVBoxLayout(self)
        layout.addWidget(self.edit)
        layout.addLayout(tools)
        layout.addWidget(self.preview, 1)
        layout.addWidget(self.info)
        layout.addWidget(self.buttons)
        self.resize(520, 520)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(PARSE_DELAY_MS)
        self._timer.timeout.connect(self.update_preview)
        self.edit.textChanged.connect(self._timer.start)
        self.paste_button.clicked.connect(self.paste)
        self.clean_button.clicked.connect(self.canonicalize)
        self.copy_button.clicked.connect(lambda: copy_structure(self.current.mol))
        self.clear_button.clicked.connect(self.edit.clear)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        self.current: StructureInput = read_structure("")
        self.result_input: StructureInput | None = None
        self.edit.setPlainText(text)
        self.update_preview()
        self.edit.setFocus()
        self.edit.selectAll()

    @property
    def text(self) -> str:
        return self.edit.toPlainText()

    def set_text(self, text: str) -> None:
        self.edit.setPlainText(text)
        self.update_preview()

    def paste(self) -> None:
        text = clipboard_structure_text()
        if text:
            self.set_text(text)

    def canonicalize(self) -> None:
        self.update_preview()
        if self.current.mol is not None:
            self.set_text(self.current.canonical_smiles())

    def update_preview(self) -> None:
        self._timer.stop()
        self.current = parsed = read_structure(self.text)
        ok = self.buttons.button(QDialogButtonBox.StandardButton.Ok)
        ok.setEnabled(parsed.is_valid)
        self.clean_button.setEnabled(parsed.mol is not None)
        self.copy_button.setEnabled(parsed.mol is not None)
        self.preview.setPixmap(QPixmap())
        if parsed.is_blank:
            self.preview.setText("Empty: OK clears the cell")
            self.info.setText(" ")
            self.info.setStyleSheet("")
            return
        if parsed.mol is None:
            self.preview.setText("")
            self.info.setText(parsed.error)
            self.info.setStyleSheet("color: #c62828")
            return
        self.preview.setText("")
        self.preview.setPixmap(self._draw(parsed))
        kind = {"smiles": "SMILES", "molfile": "Molfile", "inchi": "InChI"}[parsed.kind]
        self.info.setText(f"{kind} · {parsed.formula} · MW {parsed.mol_weight:.2f}")
        self.info.setStyleSheet("color: palette(placeholder-text)")

    def _draw(self, parsed: StructureInput) -> QPixmap:
        dpr = self.devicePixelRatioF()
        size = self.preview.size()
        dark = self.palette().base().color().lightness() < 128
        png = render_png(
            parsed.mol,
            int(size.width() * dpr),
            int(size.height() * dpr),
            scale=dpr * 1.6,
            dark=dark,
        )
        pixmap = QPixmap()
        pixmap.loadFromData(png, "PNG")
        pixmap.setDevicePixelRatio(dpr)
        return pixmap

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        if self.current.mol is not None:
            self._timer.start()

    def accept(self) -> None:
        self.update_preview()
        if not self.current.is_valid:
            return
        self.result_input = self.current
        super().accept()
