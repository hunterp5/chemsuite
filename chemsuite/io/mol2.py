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
"""Tripos mol2 files: parallel reading (RDKit's parser) and writing (our own SYBYL typing).

Reading follows the SD reader: the file is split into whole ``@<TRIPOS>MOLECULE`` records,
parsed in worker processes, and a record RDKit cannot parse is **kept** with its text, flagged
invalid. Each record's name becomes the ``Title`` column.

RDKit reads mol2 but does not write it, so the writer assigns SYBYL atom and bond types
itself, following the Tripos mol2 file format specification (Tripos/Certara, "Mol2 File
Format", the ``@<TRIPOS>ATOM`` and ``@<TRIPOS>BOND`` records and their type tables). It writes
every hydrogen, the molecule's coordinates (3D when the molecule has them) and Gasteiger
partial charges, so the file is ready for programs that read charges from it.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Sequence
from concurrent.futures import Executor
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, rdPartialCharges

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.table import Column, Table
from chemsuite.core.table.columns import StructureColumn

from .formats import Batch, atomic_write, file_size, open_binary
from .sdf import STRUCTURE, TITLE, SdfBlock, first_structure_column

BLOCK_BYTES = 1 << 20
MOLECULE = b"@<TRIPOS>MOLECULE"
_RECORD_START = re.compile(r"^(?=@<TRIPOS>MOLECULE)", re.MULTILINE)


# -- reading ------------------------------------------------------------------------------------


def iter_blocks(stream, block_bytes: int | None = None) -> Iterator[bytes]:
    """Split a binary mol2 stream into chunks of whole records (each chunk ends just before a
    ``@<TRIPOS>MOLECULE`` line, or at the end of the file)."""
    block_bytes = block_bytes or BLOCK_BYTES
    buffer = b""
    while chunk := stream.read(block_bytes):
        buffer += chunk
        cut = buffer.rfind(b"\n" + MOLECULE)
        if cut > 0:
            yield buffer[: cut + 1]
            buffer = buffer[cut + 1 :]
    if buffer.strip():
        yield buffer


def parse_block(data: bytes) -> SdfBlock:
    """Parse every record in *data* (worker-process entry point)."""
    from chemsuite.core.mols import mol_to_blob, smiles_of

    text = data.decode("utf-8", errors="replace").replace("\r\n", "\n")
    out = SdfBlock()
    with rdBase.BlockLogs():
        for record in _RECORD_START.split(text):
            if not record.startswith("@<TRIPOS>MOLECULE"):
                continue  # comments before the first record
            lines = record.split("\n")
            out.titles.append(lines[1].strip() if len(lines) > 1 else "")
            mol = Chem.MolFromMol2Block(record)
            if mol is None:
                out.blobs.append(b"")
                out.texts.append(record.rstrip("\n"))
            else:
                out.blobs.append(mol_to_blob(mol))
                out.texts.append(smiles_of(mol) or record.rstrip("\n"))
    return out


def read_mol2(
    path: str | Path,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Iterator[Batch]:
    """Yield table batches from a mol2 file (optionally gzipped)."""
    total = file_size(path)
    with open_binary(path) as (stream, raw):
        for block in ordered_map(parse_block, iter_blocks(stream), executor, cancel=cancel):
            check(cancel)
            if len(block):
                batch: dict[str, Sequence[object] | Column] = {
                    STRUCTURE: StructureColumn.from_parsed(block.blobs, block.texts)
                }
                if any(block.titles):
                    batch[TITLE] = block.titles
                yield batch
            if progress is not None:
                progress(raw.tell(), total)


# -- SYBYL types --------------------------------------------------------------------------------

_AMIDE_N = Chem.MolFromSmarts("[NX3;!$(N=*);$(N[CX3]=[OX1,SX1])]")
_AMIDE_C = Chem.MolFromSmarts("[CX3](=[OX1,SX1])[NX3]")
_GUANIDINIUM_C = Chem.MolFromSmarts("[CX3;$(C(=[NX3+])([NX3])[NX3]),$(C(=[NX3+])([NX3])[#6,#1])]")


def sybyl_types(mol: Chem.Mol) -> list[str]:
    """The SYBYL atom type of every atom (the molecule should carry explicit hydrogens)."""
    amide_n = {m[0] for m in mol.GetSubstructMatches(_AMIDE_N)}
    co2 = co2_oxygens(mol)
    cation_c = {m[0] for m in mol.GetSubstructMatches(_GUANIDINIUM_C)}
    return [_atom_type(atom, amide_n, co2, cation_c) for atom in mol.GetAtoms()]


def co2_oxygens(mol: Chem.Mol) -> set[int]:
    """The oxygens of carboxylates and phosphates (SYBYL O.co2): the terminal oxygens of a
    C or P that carries two or more of them, at least one charged."""
    found: set[int] = set()
    for centre in mol.GetAtoms():
        if centre.GetSymbol() not in ("C", "P"):
            continue
        terminal = [n for n in centre.GetNeighbors() if n.GetSymbol() == "O" and n.GetDegree() == 1]
        if len(terminal) >= 2 and any(n.GetFormalCharge() < 0 for n in terminal):
            found.update(n.GetIdx() for n in terminal)
    return found


def _atom_type(atom: Chem.Atom, amide_n: set[int], co2: set[int], cation_c: set[int]) -> str:
    symbol = atom.GetSymbol()
    index = atom.GetIdx()
    if symbol == "C":
        return "C.cat" if index in cation_c else _carbon(atom)
    if symbol == "N":
        return "N.am" if index in amide_n and not atom.GetIsAromatic() else _nitrogen(atom)
    if symbol == "O":
        return "O.co2" if index in co2 else _oxygen(atom)
    if symbol == "S":
        return _sulfur(atom)
    if symbol == "P":
        return "P.3"
    return symbol


def _orders(atom: Chem.Atom, order: Chem.BondType) -> int:
    return sum(b.GetBondType() == order for b in atom.GetBonds())


def _carbon(atom: Chem.Atom) -> str:
    hyb = atom.GetHybridization()
    if atom.GetIsAromatic():
        return "C.ar"
    if _orders(atom, Chem.BondType.TRIPLE) or hyb == Chem.HybridizationType.SP:
        return "C.1"
    if _orders(atom, Chem.BondType.DOUBLE) or hyb == Chem.HybridizationType.SP2:
        return "C.2"
    return "C.3"


def _nitrogen(atom: Chem.Atom) -> str:
    doubles = _orders(atom, Chem.BondType.DOUBLE)
    if atom.GetIsAromatic():
        return "N.ar"
    if _orders(atom, Chem.BondType.TRIPLE):
        return "N.1"
    if atom.GetFormalCharge() > 0 and atom.GetDegree() == 4 and not doubles:
        return "N.4"
    if doubles:  # a nitro N or an iminium is planar with three neighbours
        return "N.pl3" if atom.GetDegree() == 3 else "N.2"
    return "N.pl3" if _next_to_pi(atom) else "N.3"


def _oxygen(atom: Chem.Atom) -> str:
    if _orders(atom, Chem.BondType.DOUBLE):
        return "O.2"
    if atom.GetDegree() == 1 and atom.GetFormalCharge() < 0 and _next_to_pi(atom):
        return "O.2"
    return "O.3"


def _sulfur(atom: Chem.Atom) -> str:
    oxo = sum(
        b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(atom).GetSymbol() == "O"
        for b in atom.GetBonds()
    )
    if oxo == 1 and atom.GetDegree() == 3:
        return "S.O"
    if oxo >= 2:
        return "S.O2"
    if _orders(atom, Chem.BondType.DOUBLE) or atom.GetIsAromatic():
        return "S.2"
    return "S.3"


def _next_to_pi(atom: Chem.Atom) -> bool:
    """Bonded to an aromatic or multiply bonded atom (an aniline N, an enamine, a phenoxide)."""
    for neighbour in atom.GetNeighbors():
        if neighbour.GetIsAromatic():
            return True
        for bond in neighbour.GetBonds():
            if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                return True
    return False


def _bond_type(bond: Chem.Bond, kekule: Chem.BondType, amide: set[frozenset], co2: set[int]) -> str:
    """Aromatic rings are written in a Kekulé form (RDKit's reader cannot place the hydrogen
    of a pyrrole-type N from "ar" bonds); the bonds to O.co2 oxygens are "ar", as the format
    marks a delocalised carboxylate or phosphate."""
    ends = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    if co2.intersection(ends):
        return "ar"
    if frozenset(ends) in amide:
        return "am"
    return {Chem.BondType.DOUBLE: "2", Chem.BondType.TRIPLE: "3"}.get(kekule, "1")


# -- writing ------------------------------------------------------------------------------------


def mol2_block(mol: Chem.Mol, title: str = "", conformer: int = -1) -> str:
    """One ``@<TRIPOS>MOLECULE`` record: every hydrogen added (placed from the coordinates),
    Gasteiger charges, coordinates from *conformer* (2D ones computed when there are none)."""
    with rdBase.BlockLogs():
        mol = Chem.Mol(mol)
        if not mol.GetNumConformers():
            rdDepictor.Compute2DCoords(mol)
            conformer = -1
        if conformer >= 0:  # keep only the wanted conformer, so the hydrogens go on it
            keep = mol.GetConformer(conformer)
            single = Chem.Conformer(keep)
            mol.RemoveAllConformers()
            mol.AddConformer(single, assignId=True)
        mol = Chem.AddHs(mol, addCoords=True)
        rdPartialCharges.ComputeGasteigerCharges(mol)
    types = sybyl_types(mol)
    amide = {frozenset((c, n)) for c, _o, n in mol.GetSubstructMatches(_AMIDE_C)}
    co2 = co2_oxygens(mol)
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule)
    except (Chem.KekulizeException, RuntimeError, ValueError):
        kekule = mol  # leave the orders as they are
    xyz = mol.GetConformer().GetPositions()
    name = (title.splitlines()[0].strip() if title else "") or "*****"
    lines = [
        "@<TRIPOS>MOLECULE",
        name,
        f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0",
        "SMALL",
        "GASTEIGER",
        "",
        "@<TRIPOS>ATOM",
    ]
    counts: dict[str, int] = {}
    for atom, sybyl, (x, y, z) in zip(mol.GetAtoms(), types, xyz, strict=True):
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        charge = atom.GetDoubleProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else 0
        charge = 0.0 if not np.isfinite(charge) else charge
        lines.append(
            f"{atom.GetIdx() + 1:>7} {symbol + str(counts[symbol]):<8}"
            f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {sybyl:<6}{1:>5}  UNL1{charge:>12.4f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond in mol.GetBonds():
        lines.append(
            f"{bond.GetIdx() + 1:>6}{bond.GetBeginAtomIdx() + 1:>6}{bond.GetEndAtomIdx() + 1:>6}"
            f" {_bond_type(bond, kekule.GetBondWithIdx(bond.GetIdx()).GetBondType(), amide, co2)}"
        )
    lines += [
        "@<TRIPOS>SUBSTRUCTURE",
        "     1 UNL1        1 TEMP              0 ****  ****    0 ROOT",
    ]
    return "\n".join(lines) + "\n"


def write_mol2(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    structure: str | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Write rows as a mol2 file, one record per row; a ``Title`` column names each record.
    A cell with no molecule (blank or invalid) is skipped: mol2 has no empty record."""
    structure = structure or first_structure_column(table)
    column = table.column(structure)
    if not isinstance(column, StructureColumn):
        raise ValueError(f"{structure!r} is not a structure column")
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    titles = table.column(TITLE) if TITLE in table else None
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        for start in range(0, len(slots), 512):
            check(cancel)
            for slot in slots[start : start + 512].tolist():
                mol = column.mol(slot)
                if mol is None or not mol.GetNumAtoms():
                    continue
                title = titles.display(slot) if titles is not None else ""
                handle.write(mol2_block(mol, title))
            if progress is not None:
                progress(min(start + 512, len(slots)), len(slots))
