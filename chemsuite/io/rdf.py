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

"""RXN and RD files: reactions with their data (BIOVIA CTfile Formats, 2020, "RXN" and "RDfile").

An RD file is a header (``$RDFILE 1``, ``$DATM``) followed by records. A reaction record is
``$RFMT``, an RXN block and ``$DTYPE name`` / ``$DATUM value`` pairs; a molecule record is
``$MFMT``, a molfile and the same pairs. ``$RIREG`` / ``$REREG`` registry numbers go in an
``ID`` column. An RXN file is a single RXN block with no data.

As with SD files, records are split on the main process and parsed in worker processes, and a
record whose reaction does not parse is kept: its data loads and its cell holds the RXN block,
flagged invalid.

A datum runs from ``$DATUM`` to the next ``$`` line. The format wraps data at 80 characters, so
a line of exactly 80 characters is joined to the next without a line break; shorter lines keep
theirs. Data are written unwrapped.
"""

from __future__ import annotations

import datetime as dt
import re
from collections.abc import Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from chemsuite.core.jobs import CancelToken, Progress, check
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_to_blob, smiles_of
from chemsuite.core.reactions import (
    kekule_reaction,
    looks_like_rxn_block,
    parse_reaction,
    reaction_from_blob,
)
from chemsuite.core.table import Table
from chemsuite.core.table.columns import ReactionColumn, StructureColumn

from .formats import Batch, atomic_write, file_size, open_binary

REACTION = "Reaction"
STRUCTURE = "Structure"
ID = "ID"
BLOCK_BYTES = 1 << 20
WRITE_CHUNK = 256  # reactions per worker task when writing
WRAP = 80

_RECORD_START = re.compile(r"^\$[RM]FMT", re.MULTILINE)
_REGISTRY = re.compile(r"\$(?:RIREG|REREG|MIREG|MEREG)\s+(\S+)")
_EMPTY_RXN = "$RXN\n\n  chemsuite\n\n  0  0"


# -- reading ------------------------------------------------------------------------------------


@dataclass
class RdfBlock:
    """Parsed records from one block (picklable; built in a worker process)."""

    reactions: list[tuple[bytes, str]] = field(default_factory=list)
    molecules: list[tuple[bytes, str]] = field(default_factory=list)
    ids: list[str] = field(default_factory=list)
    fields: dict[str, list[str]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.ids)


def iter_blocks(stream, block_bytes: int | None = None) -> Iterator[bytes]:
    """Split a binary RD file into chunks of whole records (each chunk but the first starts at a
    ``$RFMT`` or ``$MFMT`` line)."""
    block_bytes = block_bytes or BLOCK_BYTES
    buffer = b""
    while chunk := stream.read(block_bytes):
        buffer += chunk
        cut = max(buffer.rfind(b"\n$RFMT"), buffer.rfind(b"\n$MFMT"))
        if cut > 0:
            yield buffer[: cut + 1]
            buffer = buffer[cut + 1 :]
    if buffer.strip():
        yield buffer


def parse_block(data: bytes) -> RdfBlock:
    """Parse every record in *data* (worker-process entry point). Text before the first record
    (the ``$RDFILE`` header) is skipped; a block with no records at all is one RXN file."""
    text = data.decode("utf-8", errors="replace").replace("\r\n", "\n")
    out = RdfBlock()
    starts = [m.start() for m in _RECORD_START.finditer(text)]
    with rdBase.BlockLogs():
        if not starts:
            if looks_like_rxn_block(text):
                _parse_record("$RFMT\n" + text, out)
            return out
        for begin, end in zip(starts, [*starts[1:], len(text)], strict=True):
            _parse_record(text[begin:end], out)
    return out


def _parse_record(record: str, out: RdfBlock) -> None:
    head, _, rest = record.partition("\n")
    registry = _REGISTRY.search(head)
    lines = rest.split("\n")
    if not registry and lines and _REGISTRY.match(lines[0]):
        registry = _REGISTRY.match(lines.pop(0))
    body: list[str] = []
    while lines and not lines[0].startswith("$DTYPE"):
        body.append(lines.pop(0))
    block = "\n".join(body).rstrip("\n")
    index = len(out.ids)
    out.ids.append(registry.group(1) if registry else "")
    if head.startswith("$RFMT"):
        out.reactions.append(parse_reaction(block) if block.strip() else (b"", ""))
        out.molecules.append((b"", ""))
    else:
        mol = Chem.MolFromMolBlock(block) if "M  END" in block else None
        cell = (b"", block) if mol is None else (mol_to_blob(mol), smiles_of(mol) or block)
        out.molecules.append(cell if block.strip() else (b"", ""))
        out.reactions.append((b"", ""))
    for name, value in _parse_data(lines):
        column = out.fields.get(name)
        if column is None:
            column = out.fields[name] = [""] * index
        if len(column) == index:  # a repeated name within one record: keep the first
            column.append(value)
    for column in out.fields.values():
        if len(column) == index:
            column.append("")


def _parse_data(lines: list[str]) -> Iterator[tuple[str, str]]:
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if not line.startswith("$DTYPE"):
            continue
        name = line[6:].strip()
        value = ""
        if i < len(lines) and lines[i].startswith("$DATUM"):
            value = lines[i][6:].lstrip()
            previous = lines[i]
            i += 1
            while i < len(lines) and not lines[i].startswith("$"):
                joiner = "" if len(previous) == WRAP else "\n"
                value += joiner + lines[i]
                previous = lines[i]
                i += 1
        if name:
            yield name, value.rstrip("\n")


def read_rdf(
    path: str | Path,
    *,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Iterator[Batch]:
    """Yield table batches from an RD or RXN file (optionally gzipped)."""
    total = file_size(path)
    with open_binary(path) as (stream, raw):
        for block in ordered_map(parse_block, iter_blocks(stream), executor, cancel=cancel):
            check(cancel)
            if len(block):
                yield _batch(block)
            if progress is not None:
                progress(raw.tell(), total)


def _batch(block: RdfBlock) -> Batch:
    batch: dict[str, object] = {}
    if any(text for _, text in block.reactions) or not any(t for _, t in block.molecules):
        blobs, texts = zip(*block.reactions, strict=True)
        batch[REACTION] = ReactionColumn.from_parsed(list(blobs), list(texts))
    if any(text for _, text in block.molecules):
        blobs, texts = zip(*block.molecules, strict=True)
        batch[STRUCTURE] = StructureColumn.from_parsed(list(blobs), list(texts))
    if any(block.ids):
        batch[ID] = block.ids
    for name, values in block.fields.items():
        batch[f"{name} (RDF)" if name in (REACTION, STRUCTURE, ID) else name] = values
    return batch  # type: ignore[return-value]


# -- writing ------------------------------------------------------------------------------------


def write_rdf(
    table: Table,
    path: str | Path,
    *,
    oids: Sequence[int] | np.ndarray | None = None,
    reaction: str | None = None,
    columns: Sequence[str] | None = None,
    executor: Executor | None = None,
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> None:
    """Write rows as an RD file: one ``$RFMT`` record per row with the reaction as an RXN block
    and the other columns (all by default) as ``$DTYPE`` / ``$DATUM`` data.

    *reaction* defaults to the first reaction column. Reactions with no coordinates are laid
    out as the sketcher lays them out, in worker processes with *executor*. Blank data are
    left out of a record.
    """
    reaction = reaction or first_reaction_column(table)
    column = table.column(reaction)
    if not isinstance(column, ReactionColumn):
        raise ValueError(f"{reaction!r} is not a reaction column")
    fields = [c for c in (columns or table.column_names) if c != reaction]
    slots = np.arange(len(table)) if oids is None else table.slots_of(oids)
    stamp = dt.datetime.now().strftime("%m/%d/%y %H:%M")
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"$RDFILE 1\n$DATM    {stamp}\n")
        chunks = [slots[i : i + WRITE_CHUNK] for i in range(0, len(slots), WRITE_CHUNK)]
        cells = ((column.blobs(c), column.displays(c)) for c in chunks)
        done = 0
        for chunk, blocks in zip(
            chunks, ordered_map(rxn_blocks, cells, executor, cancel=cancel), strict=False
        ):
            check(cancel)
            values = {name: table.column(name).displays(chunk) for name in fields}
            for i, block in enumerate(blocks):
                handle.write(f"$RFMT\n{block}\n")
                for name in fields:
                    value = values[name][i]
                    if value:
                        handle.write(f"$DTYPE {name}\n$DATUM {value}\n")
            done += len(chunk)
            if progress is not None:
                progress(done, len(slots))


def write_rxn(table: Table, path: str | Path, oid: int, reaction: str | None = None) -> None:
    """Write one row's reaction as an RXN file."""
    reaction = reaction or first_reaction_column(table)
    column = table.column(reaction)
    if not isinstance(column, ReactionColumn):
        raise ValueError(f"{reaction!r} is not a reaction column")
    with atomic_write(path, "w", encoding="utf-8", newline="\n") as handle:
        slots = table.slots_of([oid])
        handle.write(rxn_blocks((column.blobs(slots), column.displays(slots)))[0] + "\n")


def rxn_block(rxn: rdChemReactions.ChemicalReaction) -> str:
    """An RXN file (V2000; V3000 when a component needs it) with agents kept separate, written
    with Kekulé bonds as molfiles are, laid out if it has no coordinates."""
    from chemsuite.sketch import reaction as sketch_reaction  # its layout; core.io sits above it

    templates = (*rxn.GetReactants(), *rxn.GetAgents(), *rxn.GetProducts())
    if all(t.GetNumConformers() for t in templates):
        copy = rdChemReactions.ChemicalReaction(rxn)
        with rdBase.BlockLogs():
            try:
                return rdChemReactions.ReactionToRxnBlock(
                    kekule_reaction(copy), separateAgents=True
                )
            except (RuntimeError, ValueError):
                pass
    return sketch_reaction.to_rxn_block(sketch_reaction.from_reaction(rxn))


def rxn_blocks(cells: tuple[list[bytes | None], list[str]]) -> list[str]:
    """RXN blocks for ``(blobs, texts)`` of reaction cells (worker-process entry point). A cell
    not parsed yet is parsed here; one that does not parse is written back as read if it is an
    RXN block, else as an empty reaction."""
    out: list[str] = []
    for blob, text in zip(*cells, strict=True):
        if blob is None and text:
            blob = parse_reaction(text)[0]
        if blob:
            out.append(rxn_block(reaction_from_blob(blob)).rstrip("\n"))
        else:
            out.append(text.rstrip("\n") if looks_like_rxn_block(text) else _EMPTY_RXN)
    return out


def first_reaction_column(table: Table) -> str:
    for name in table.column_names:
        if isinstance(table.column(name), ReactionColumn):
            return name
    raise ValueError("the table has no reaction column")


__all__ = [
    "first_reaction_column",
    "parse_block",
    "read_rdf",
    "rxn_block",
    "rxn_blocks",
    "write_rdf",
    "write_rxn",
]
