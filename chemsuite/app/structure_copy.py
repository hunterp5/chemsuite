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

"""Copy structures to the clipboard in a chosen format (the structure cell's Copy As menu).

Several structures copy as one text: SMILES, InChI and InChIKey one per line, molfiles as an
SD file (``$$$$`` between records). A molfile copy also offers the ``chemical/x-mdl-molfile``
clipboard format, so a single structure pastes into ChemDraw and other drawing programs as a
drawing. An image copies the first structure as a PNG picture.
"""

from __future__ import annotations

from collections.abc import Sequence

from PySide6.QtCore import QMimeData
from PySide6.QtGui import QGuiApplication, QImage
from rdkit import Chem, rdBase

from chemsuite.core.depict import render_png
from chemsuite.core.mols import smiles_of
from chemsuite.core.structure_input import molblock_of

#: (key, menu title) in menu order
FORMATS = (
    ("molfile", "Molfile"),
    ("smiles", "SMILES"),
    ("inchi", "InChI"),
    ("inchikey", "InChIKey"),
    ("image", "Image"),
)
IMAGE_SIZE = (600, 450)


def structure_text(mol: Chem.Mol, fmt: str) -> str:
    """*mol* written as *fmt* ("" when RDKit cannot write it in that format)."""
    with rdBase.BlockLogs():
        if fmt == "molfile":
            return molblock_of(mol)
        if fmt == "smiles":
            return smiles_of(mol)
        if fmt == "inchi":
            return Chem.MolToInchi(mol) or ""
        if fmt == "inchikey":
            return Chem.MolToInchiKey(mol) or ""
    raise ValueError(f"unknown text format {fmt!r}")


def copy_structures(mols: Sequence[Chem.Mol | None], fmt: str) -> str:
    """Put *mols* on the clipboard as *fmt*; returns a status message."""
    mols = [m for m in mols if m is not None]
    title = dict(FORMATS)[fmt]
    if not mols:
        return "No valid structure to copy"
    if fmt == "image":
        image = QImage.fromData(render_png(mols[0], *IMAGE_SIZE), "PNG")
        QGuiApplication.clipboard().setImage(image)
        return "Copied structure as an image"
    texts = [structure_text(m, fmt) for m in mols]
    missing = sum(1 for t in texts if not t)
    texts = [t for t in texts if t]
    if not texts:
        return f"RDKit could not write this structure as {title}"
    mime = QMimeData()
    if fmt == "molfile":
        records = [t if t.endswith("\n") else t + "\n" for t in texts]
        mime.setText("$$$$\n".join(records) + ("$$$$\n" if len(records) > 1 else ""))
        if len(records) == 1:
            mime.setData("chemical/x-mdl-molfile", records[0].encode())
    else:
        mime.setText("\n".join(texts))
    QGuiApplication.clipboard().setMimeData(mime)
    what = "structure" if len(texts) == 1 else f"{len(texts)} structures"
    skipped = f"; {missing} could not be written" if missing else ""
    return f"Copied {what} as {title}{skipped}"
