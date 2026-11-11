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

"""Structures on the clipboard: reading the molfile formats drawing programs put there, and
copying a structure as a molfile plus SMILES text."""

from __future__ import annotations

from PySide6.QtCore import QMimeData
from PySide6.QtGui import QGuiApplication
from rdkit import Chem

from chemsuite.core.mols import looks_like_molblock, smiles_of
from chemsuite.core.structure_input import molblock_of

#: Clipboard formats that carry a molfile (ChemDraw, MarvinSketch, Biovia Draw and others).
MOLFILE_FORMATS = (
    "chemical/x-mdl-molfile",
    "chemical/x-mdl-molfile-v3000",
    "application/x-mdl-molfile",
    "MDLCT",
)


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
