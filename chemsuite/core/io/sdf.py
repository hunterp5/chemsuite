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

"""SD files: parallel, lossless reading and streaming writing.

Reading splits the file into blocks of whole records (no parsing on the main process) and
parses blocks in worker processes. Unlike RDKit's supplier, a record whose molecule does not
parse is **kept**: its data fields load normally and its structure cell holds the original
molfile, flagged invalid, so nothing in the user's file silently disappears.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import looks_like_molblock, mol_to_blob, smiles_of
from chemsuite.core.table import Table
from chemsuite.core.table.columns import StructureColumn

from .formats import Batch, atomic_write, file_size, open_binary

STRUCTURE = "Structure"
TITLE = "Title"
BLOCK_BYTES = 1 << 20  # ~300 drug-like records per worker task

_RECORD_END = re.compile(r"^\$\$\$\$[^\n]*(?:\n|\Z)", re.MULTILINE)
_FIELD_NAME = re.compile(r"<([^>]*)>")
_EMPTY_MOLBLOCK = "\n  chemsuite\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END"


# -- reading ------------------------------------------------------------------------------------


@dataclass
class SdfBlock:
    """Parsed records from one block (picklable; built in a worker process)."""

    titles: list[str] = field(default_factory=list)
    blobs: list[bytes] = field(default_factory=list)
    texts: list[str] = field(default_factory=list)
    fields: dict[str, list[str]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.titles)


def iter_blocks(stream, block_bytes: int | None = None) -> Iterator[bytes]:
    """Split a binary SDF stream into chunks that each end exactly after a ``$$$$`` line."""
    block_bytes = block_bytes or BLOCK_BYTES
    buffer = b""
    while chunk := stream.read(block_bytes):
        buffer += chunk
        cut = _last_record_end(buffer)
        if cut:
            yield buffer[:cut]
            buffer = buffer[cut:]
    if buffer.strip():
        yield buffer


def _last_record_end(buffer: bytes) -> int:
    marker = buffer.rfind(b"\n$$$$")
    start = marker + 1 if marker >= 0 else (0 if buffer.startswith(b"$$$$") else -1)
    if start < 0:
        return 0
    newline = buffer.find(b"\n", start)
    return newline + 1 if newline >= 0 else 0  # wait until the $$$$ line is complete


def parse_block(data: bytes) -> SdfBlock:
    """Parse every record in *data* (worker-process entry point)."""
    text = data.decode("utf-8", errors="replace").replace("\r\n", "\n")
    out = SdfBlock()
    with rdBase.BlockLogs():
        for record in _RECORD_END.split(text):
            if record.strip():
                _parse_record(record, out)
    return out


def _parse_record(record: str, out: SdfBlock) -> None:
    end = record.find("\nM  END")
    if end >= 0:
        molblock, rest = record[: end + 7], record[end + 7 :]
    else:
        molblock, rest = record, ""
    index = len(out.titles)
    out.titles.append(molblock.split("\n", 1)[0].strip())
    mol = Chem.MolFromMolBlock(molblock) if end >= 0 else None
    if mol is not None and mol.GetNumAtoms() == 0:
        out.blobs.append(b"")  # an empty connection table is a blank structure cell
        out.texts.append("")
    elif mol is None:
        out.blobs.append(b"")
        out.texts.append(molblock.rstrip("\n"))
    else:
        out.blobs.append(mol_to_blob(mol))
        out.texts.append(smiles_of(mol) or molblock.rstrip("\n"))
    for name, value in _parse_fields(rest):
        column = out.fields.get(name)
        if column is None:
            column = out.fields[name] = [""] * index
        if len(column) == index:
            column.append(value)
        else:  # repeated field name within one record: keep the first
            continue
    for column in out.fields.values():
        if len(column) == index:
            column.append("")


def _parse_fields(text: str) -> Iterator[tuple[str, str]]:
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if not line.startswith(">"):
            continue
        match = _FIELD_NAME.search(line)
        name = (match.group(1) if match else line[1:]).strip()
        value_lines: list[str] = []
        while i < len(lines) and lines[i].strip() != "":
            value_lines.append(lines[i].rstrip())
            i += 1
        if name:
            yield name, "\n".join(value_lines)


def read_sdf(
    path: str | Path,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Iterator[Batch]:
    """Yield table batches from an SD file (optionally gzipped)."""
    total = file_size(path)
    with open_binary(path) as (stream, raw):
        for block in ordered_map(parse_block, iter_blocks(stream), executor):
            check(cancel)
            if len(block):
                yield _batch(block)
            if progress is not None:
                progress(raw.tell(), total)


def _batch(block: SdfBlock) -> Batch:
    batch: dict[str, object] = {STRUCTURE: StructureColumn.from_parsed(block.blobs, block.texts)}
    if any(block.titles):
        batch[TITLE] = block.titles
    for name, values in block.fields.items():
        batch[f"{name} (SD)" if name in (STRUCTURE, TITLE) else name] = values
    return batch  # type: ignore[return-value]


# -- writing ------------------------------------------------------------------------------------


def write_sdf(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    structure: str | None = None,
    columns: Sequence[str] | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
    all_conformers: bool = False,
) -> None:
    """Write rows as an SD file. Other columns become data fields (all of them by default).

    *structure* defaults to the first structure column. A ``Title`` column becomes each
    record's title line. Invalid structures read from a molfile are written back as read.
    With *all_conformers*, a row whose molecule has several conformers is written once per
    conformer (lowest energy first, as stored), with a ``Conformer`` field numbering them.
    """
    structure = structure or first_structure_column(table)
    column = table.column(structure)
    if not isinstance(column, StructureColumn):
        raise ValueError(f"{structure!r} is not a structure column")
    fields = [c for c in (columns or table.column_names) if c not in (structure, TITLE)]
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    titles = table.column(TITLE) if TITLE in table else None
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        for start in range(0, len(slots), 2048):
            check(cancel)
            chunk = slots[start : start + 2048]
            values = {name: table.column(name).displays(chunk) for name in fields}
            for i, slot in enumerate(chunk.tolist()):
                title = titles.display(slot) if titles is not None else ""
                mol = column.mol(slot) if all_conformers else None
                conformers = mol.GetNumConformers() if mol is not None else 0
                for index in range(max(1, conformers)):
                    conformer = index if conformers > 1 else None
                    handle.write(_molblock(column, slot, title, conformer))
                    handle.write("\n")
                    if conformer is not None:
                        handle.write(f"> <Conformer>\n{index + 1}\n\n")
                    for name in fields:
                        handle.write(f"> <{name}>\n{values[name][i]}\n\n")
                    handle.write("$$$$\n")
            if progress is not None:
                progress(min(start + 2048, len(slots)), len(slots))


def _molblock(column: StructureColumn, slot: int, title: str, conformer: int | None = None) -> str:
    """One molfile; *conformer* picks a conformer by position (default: the first)."""
    mol = column.mol(slot)
    if mol is not None:
        with rdBase.BlockLogs():
            if not mol.GetNumConformers():
                mol = Chem.Mol(mol)
                rdDepictor.Compute2DCoords(mol)
            conf_id = mol.GetConformers()[conformer].GetId() if conformer is not None else -1
            try:
                block = Chem.MolToMolBlock(mol, confId=conf_id)
            except (RuntimeError, ValueError):
                block = Chem.MolToMolBlock(mol, kekulize=False, confId=conf_id)
    else:
        text = column.display(slot)
        block = text if looks_like_molblock(text) else _EMPTY_MOLBLOCK
    _, _, rest = block.partition("\n")  # replace the title line
    return f"{title.splitlines()[0] if title else ''}\n{rest}".rstrip("\n")


def first_structure_column(table: Table) -> str:
    for name in table.column_names:
        if isinstance(table.column(name), StructureColumn):
            return name
    raise ValueError("the table has no structure column")
