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

"""Converting between structure text, RDKit molecules and the compact binary the table stores.

The table keeps each structure as RDKit's binary pickle (connection table plus the file's
coordinates, no molecule properties): about 200–600 bytes, immutable, cheap to send to worker
processes and fast to turn back into a molecule. RDKit's own warnings are silenced here because
a bad cell is reported through the cell's status, not through stderr.
"""

from __future__ import annotations

from rdkit import Chem, rdBase

#: Exact pickle content, independent of any global RDKit default another module might change.
_PICKLE_FLAGS = Chem.PropertyPickleOptions.NoProps


def looks_like_molblock(text: str) -> bool:
    return "\n" in text and "M  END" in text


def mol_from_text(text: str) -> Chem.Mol | None:
    """Parse SMILES, a molfile block or an InChI. None when the text is not a valid structure."""
    # A molfile's first three lines are positional (the title may be empty): only trim the end.
    text = text.rstrip() if looks_like_molblock(text) else text.strip()
    if not text:
        return None
    with rdBase.BlockLogs():
        if looks_like_molblock(text):
            return Chem.MolFromMolBlock(text)
        if text.startswith("InChI="):
            return Chem.MolFromInchi(text)
        return Chem.MolFromSmiles(text)


def mol_to_blob(mol: Chem.Mol) -> bytes:
    return mol.ToBinary(_PICKLE_FLAGS)


def mol_from_blob(blob: bytes) -> Chem.Mol:
    return Chem.Mol(blob)


def smiles_of(mol: Chem.Mol) -> str:
    """Canonical isomeric SMILES ("" if RDKit cannot write one)."""
    with rdBase.BlockLogs():
        try:
            return Chem.MolToSmiles(mol)
        except (RuntimeError, ValueError):
            return ""


def has_query_atoms(mol: Chem.Mol) -> bool:
    """Whether *mol* has wildcard or atom-list atoms (molfile A, Q, X, M, L, ...), which SMILES
    cannot write. R-groups are not counted: SMILES writes them as ``[1*]``."""
    return any(
        atom.GetAtomicNum() == 0
        and atom.HasQuery()
        and bool(atom.GetQueryType() or atom.HasProp("_MolFileAtomQuery"))
        for atom in mol.GetAtoms()
    )


def display_text(mol: Chem.Mol) -> str:
    """How a parsed molfile is shown in a cell: canonical SMILES, or SMARTS when the drawing
    has query atoms that SMILES would turn into plain ``*``."""
    if has_query_atoms(mol):
        with rdBase.BlockLogs():
            try:
                return Chem.MolToSmarts(mol)
            except (RuntimeError, ValueError):
                pass
    return smiles_of(mol)


def parse_structure(text: str) -> tuple[bytes, str]:
    """``(blob, display text)`` for one structure cell; ``(b"", text)`` when it does not parse.

    Valid SMILES and InChI keep the text exactly as written. A molfile block is displayed as its
    canonical SMILES, or SMARTS when it has query atoms (the block's coordinates and queries
    are kept in the blob).
    """
    mol = mol_from_text(text)
    if mol is None:
        return b"", text
    return mol_to_blob(mol), (display_text(mol) if looks_like_molblock(text) else text.strip())


def parse_structures(texts: list[str]) -> list[bytes]:
    """Blobs for many texts (``b""`` where invalid). Picklable entry point for worker processes."""
    return [parse_structure(text)[0] if text else b"" for text in texts]
