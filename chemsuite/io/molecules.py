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

"""Molecule files RDKit reads whole: ChemDraw (CDXML, CDX), Marvin (MRV), Maestro (MAE,
MAEGZ), XYZ and InChI.

Like the SD reader, a record RDKit cannot parse is **kept** with its text, flagged invalid, so
nothing in the user's file silently disappears. Maestro properties become columns.

Writing covers the formats RDKit writes: Maestro (with the other columns as properties), XYZ
(3D coordinates only), InChI lists, and one molecule as CDXML or MRV (RDKit merges several
molecules written to one CDXML or MRV file into a single structure).
"""

from __future__ import annotations

import io
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdDetermineBonds

from chemsuite.core.jobs import CancelToken, check
from chemsuite.core.mols import as_drawn, as_read, mol_to_blob, smiles_of, without_hydrogens
from chemsuite.core.table import Column, Table
from chemsuite.core.table.columns import NumberColumn, StructureColumn

from .formats import (
    CDX,
    CDXML,
    INCHI,
    MAE,
    MRV,
    XYZ,
    Batch,
    atomic_write,
    decode,
    open_binary,
    split_suffix,
)
from .sdf import STRUCTURE, TITLE, first_structure_column

BATCH_ROWS = 2_000


@dataclass
class _Records:
    titles: list[str] = field(default_factory=list)
    blobs: list[bytes] = field(default_factory=list)
    texts: list[str] = field(default_factory=list)
    fields: dict[str, list[object]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.titles)

    def add(
        self, mol: Chem.Mol | None, text: str, title: str = "", props: dict | None = None
    ) -> None:
        index = len(self.titles)
        self.titles.append(title)
        if mol is None or mol.GetNumAtoms() == 0:
            self.blobs.append(b"")
            self.texts.append(text if mol is None else "")
        else:
            self.blobs.append(mol_to_blob(mol))
            self.texts.append(smiles_of(without_hydrogens(mol)) or text or "?")
        for name, value in (props or {}).items():
            column = self.fields.setdefault(name, [None] * index)
            column.append(value)
        for column in self.fields.values():
            if len(column) == index:
                column.append(None)

    def batch(self) -> Batch:
        out: dict[str, Sequence[object] | Column] = {
            STRUCTURE: StructureColumn.from_parsed(self.blobs, self.texts)
        }
        if any(self.titles):
            out[TITLE] = self.titles
        out.update({name: values for name, values in self.fields.items() if name not in out})
        return out


def read_molecules(path: str | Path, *, cancel: CancelToken | None = None) -> Iterator[Batch]:
    """Yield table batches from a CDXML, CDX, MRV, MAE, XYZ or InChI file."""
    fmt, _ = split_suffix(path)
    records = _Records()
    with rdBase.BlockLogs():
        for mol, text, title, props in _guarded(_RECORDS[fmt](path)):
            check(cancel)
            records.add(mol, text, title, props)
            if len(records) >= BATCH_ROWS:
                yield records.batch()
                records = _Records()
    if len(records):
        yield records.batch()


Record = tuple[Chem.Mol | None, str, str, dict]


def _guarded(records: Iterator[Record]) -> Iterator[Record]:
    """*records*, until RDKit gives up on the file: the molecules read so far are kept and
    the rest of the file becomes one invalid row that says why."""
    while True:
        try:
            record = next(records)
        except StopIteration:
            return
        except (RuntimeError, ValueError) as err:
            reason = str(err).strip().splitlines()
            yield None, f"unreadable: {reason[0] if reason else type(err).__name__}", "", {}
            return
        yield record


def _read_text(path: str | Path) -> str:
    with open_binary(path) as (stream, _raw):
        return decode(stream.read()).replace("\r\n", "\n")


def _title(mol: Chem.Mol) -> str:
    return mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""


def _chemdraw(path: str | Path) -> Iterator[Record]:
    params = Chem.CDXMLParserParams()
    if split_suffix(path)[0] == CDX:
        params.format = Chem.CDXMLFormat.CDX
        mols = Chem.MolsFromCDXMLFile(str(path), params)  # binary: read from the path
    else:
        mols = Chem.MolsFromCDXML(_read_text(path), params)
    found = False
    for mol in mols:
        found = True
        yield mol, "", _title(mol), {}
    if not found:
        yield None, Path(path).name, "", {}


def _marvin(path: str | Path) -> Iterator[Record]:
    text = _read_text(path)
    mol = Chem.MolFromMrvBlock(text)
    yield mol, Path(path).name, _title(mol) if mol is not None else "", {}


def _maestro(path: str | Path) -> Iterator[Record]:
    with open_binary(path) as (stream, _raw):
        data = stream.read()
    supplier = Chem.MaeMolSupplier(io.BytesIO(data), removeHs=False)
    for read in supplier:
        mol = as_read(read)
        if mol is None:
            yield None, "unreadable Maestro record", "", {}
            continue
        yield mol, "", _title(mol), mol.GetPropsAsDict(includePrivate=False, includeComputed=False)


def split_xyz(text: str) -> Iterator[tuple[str, str]]:
    """``(title, block)`` for each frame of an XYZ file (atom count, comment, atom lines)."""
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            count = int(lines[i].split()[0])
        except ValueError:
            yield "", "\n".join(lines[i:]).strip()  # not a frame: keep the rest as one record
            return
        frame = lines[i : i + count + 2]
        yield (frame[1].strip() if len(frame) > 1 else ""), "\n".join(frame)
        i += count + 2


def _xyz(path: str | Path) -> Iterator[Record]:
    for title, block in split_xyz(_read_text(path)):
        mol = Chem.MolFromXYZBlock(block)
        if mol is not None:
            try:
                rdDetermineBonds.DetermineBonds(mol, charge=0)
            except (RuntimeError, ValueError):
                mol = None  # no sensible bond orders for a neutral molecule
        yield mol, block, title, {}


def _inchi(path: str | Path) -> Iterator[Record]:
    """One InChI per line, with a name before or after it; a line without one is kept."""
    for line in _read_text(path).split("\n"):
        line = line.strip()
        if not line or line.startswith(("#", "AuxInfo=")):
            continue
        tokens = line.split()
        inchi = next((t for t in tokens if t.startswith("InChI=")), None)
        if inchi is None:
            yield None, line, "", {}
            continue
        name = " ".join(t for t in tokens if t is not inchi)
        yield Chem.MolFromInchi(inchi), inchi, name, {}


_RECORDS = {CDXML: _chemdraw, CDX: _chemdraw, MRV: _marvin, MAE: _maestro, XYZ: _xyz, INCHI: _inchi}
FORMATS = frozenset(_RECORDS)


# -- writing ------------------------------------------------------------------------------------

#: Formats chemsuite writes with :func:`write_molecules`.
WRITE_FORMATS = frozenset({CDXML, MRV, MAE, XYZ, INCHI})
_ONE_MOLECULE = frozenset({CDXML, MRV})


def write_molecules(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | None = None,
    structure: str | None = None,
    columns: Sequence[str] | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Write rows as a Maestro, XYZ, InChI, CDXML or MRV file (chosen by *path*'s suffix).

    Rows without a valid structure are left out: these formats cannot hold one. A ``Title``
    column names each molecule; for Maestro the other columns (all by default) become
    properties, numbers as real-valued ones.
    """
    fmt, _ = split_suffix(path)
    if fmt not in WRITE_FORMATS:
        raise ValueError(f"chemsuite does not write {fmt.upper()} files")
    structure = structure or first_structure_column(table)
    column = table.column(structure)
    if not isinstance(column, StructureColumn):
        raise ValueError(f"{structure!r} is not a structure column")
    slots = range(len(table)) if oids is None else table.slots_of(oids).tolist()
    titles = table.column(TITLE) if TITLE in table else None
    rows = [(slot, column.mol(slot)) for slot in slots]
    rows = [(slot, mol) for slot, mol in rows if mol is not None and mol.GetNumAtoms()]
    if not rows:
        raise ValueError("no row has a valid structure to write")
    if fmt in _ONE_MOLECULE and len(rows) != 1:
        raise ValueError(
            f"a {fmt.upper()} file holds one molecule: export one row, or use an SD file"
        )
    if fmt == XYZ:
        flat = sum(1 for _, mol in rows if not _has_3d(mol))
        if flat:
            raise ValueError(
                f"XYZ needs 3D coordinates and {flat} of the rows have none: "
                "generate conformers first"
            )
    named = []
    for slot, mol in rows:
        check(cancel)
        mol = Chem.Mol(mol)  # the column's molecule is shared
        mol.SetProp("_Name", titles.display(slot) if titles is not None else "")
        named.append((slot, mol))
    if fmt == MAE:
        _write_maestro(table, path, named, structure, columns)
        return
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        if fmt == XYZ:  # every atom, hydrogens included, as quantum chemistry programs expect
            for _, mol in named:
                handle.write(Chem.MolToXYZBlock(Chem.AddHs(mol, addCoords=True)))
        elif fmt == INCHI:
            for _, mol in named:
                name = mol.GetProp("_Name")
                handle.write(Chem.MolToInchi(as_drawn(mol)) + (f"\t{name}" if name else "") + "\n")
        elif fmt == MRV:
            handle.write(Chem.MolToMrvBlock(named[0][1]))
        else:
            handle.write(Chem.MolToCDXMLBlock(named[0][1]))


def _has_3d(mol: Chem.Mol) -> bool:
    return mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def _write_maestro(
    table: Table,
    path: str | Path,
    rows: list[tuple[int, Chem.Mol]],
    structure: str,
    columns: Sequence[str] | None,
) -> None:
    names = [c for c in (columns or table.column_names) if c not in (structure, TITLE)]
    names = [n for n in names if not isinstance(table.column(n), StructureColumn)]
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        writer = Chem.MaeWriter(handle)
        for slot, mol in rows:
            for name in names:
                values = table.column(name)
                if isinstance(values, NumberColumn):
                    number = values.value(slot)
                    if number is not None:
                        mol.SetDoubleProp(name, number)
                elif text := values.display(slot):
                    mol.SetProp(name, text)
            writer.write(mol)
        writer.close()
