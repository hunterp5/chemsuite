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
"""Reading and writing PDB and mmCIF files, through gemmi.

gemmi reads the file (gzipped or not) and sorts out entities and the SEQRES numbering; this
module turns its hierarchy into a :class:`~.model.Structure` (classified, bonded and with
secondary structure) and back. Every model of an NMR file is kept when they have the same
atoms. Documents do not store proteins as mmCIF (:mod:`.store` packs the arrays), so these
are for import and export only.
"""

from __future__ import annotations

import logging
import os
import string
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import replace
from os import PathLike
from pathlib import Path

import gemmi
import numpy as np

from . import bonds, classify, components, secondary
from .model import AMINO, COIL, HELIX, LIGAND, STRAND, Structure

log = logging.getLogger(__name__)

PDB = "pdb"
MMCIF = "cif"
#: file suffixes read, after an optional ``.gz``
SUFFIXES = {".pdb": PDB, ".ent": PDB, ".cif": MMCIF, ".mmcif": MMCIF}

_LINKS = {gemmi.ConnectionType.Covale, gemmi.ConnectionType.Disulf}


class ProteinFileError(ValueError):
    """The file could not be read as a structure."""


def read(path: str | PathLike[str]) -> Structure:
    """Read a PDB or mmCIF file (``.gz`` too)."""
    doc = gemmi.cif.Document()
    try:
        st = gemmi.read_structure(str(path), save_doc=doc)
    except (RuntimeError, ValueError) as error:
        raise ProteinFileError(f"{Path(path).name}: {error}") from error
    return from_gemmi(st, doc[0] if len(doc) else None)


def parse(text: str, fmt: str) -> Structure:
    """Read a structure from text in ``fmt`` (PDB or MMCIF)."""
    kind = gemmi.CoorFormat.Pdb if fmt == PDB else gemmi.CoorFormat.Mmcif
    doc = gemmi.cif.Document()
    try:
        st = gemmi.read_structure_string(text, format=kind, save_doc=doc)
    except (RuntimeError, ValueError) as error:
        raise ProteinFileError(str(error)) from error
    return from_gemmi(st, doc[0] if len(doc) else None)


def _suffix_format(path: str | PathLike[str]) -> str | None:
    suffixes = [s.lower() for s in Path(path).suffixes if s.lower() != ".gz"]
    return SUFFIXES.get(suffixes[-1]) if suffixes else None


def from_gemmi(st: gemmi.Structure, block: gemmi.cif.Block | None = None) -> Structure:
    """``st`` as a Structure; ``block``, the mmCIF block it was read from, gives the bond
    orders of its ligands (its chemical component tables)."""
    if not len(st) or not st[0].count_atom_sites():
        raise ProteinFileError("the file has no atoms")
    st.setup_entities()
    st.assign_label_seq_id(False)
    entity_index = {e.name: i for i, e in enumerate(st.entities)}
    subchain_entity = {sub: entity_index[e.name] for e in st.entities for sub in e.subchains}
    atoms: dict[str, list] = {k: [] for k in ("name", "element", "altloc", "occ", "b", "charge")}
    xyz: list[gemmi.Position] = []
    res: dict[str, list] = {k: [] for k in ("start", "name", "seq", "icode", "kind", "seqres",
                                            "asym", "entity")}  # fmt: skip
    chain_names: list[str] = []
    chain_start = [0]
    for chain in st[0]:
        if not len(chain):
            continue
        chain_names.append(chain.name)
        for residue in chain:
            res["start"].append(len(xyz))
            heavy: list[str] = []
            for atom in residue:
                z = atom.element.atomic_number
                atoms["name"].append(atom.name)
                atoms["element"].append(z)
                atoms["altloc"].append(atom.altloc if atom.has_altloc() else "")
                atoms["occ"].append(atom.occ)
                atoms["b"].append(atom.b_iso)
                atoms["charge"].append(atom.charge)
                xyz.append(atom.pos)
                if z != 1:
                    heavy.append(atom.element.name.upper())
            res["name"].append(residue.name)
            res["seq"].append(residue.seqid.num)
            res["icode"].append(residue.seqid.icode.strip())
            res["kind"].append(classify.residue_kind(residue.name, residue.entity_type,
                                                     tuple(heavy)))  # fmt: skip
            label_seq = residue.label_seq
            res["seqres"].append(label_seq - 1 if label_seq is not None else -1)
            res["asym"].append(residue.subchain)
            res["entity"].append(subchain_entity.get(residue.subchain, -1))
        chain_start.append(len(res["name"]))
    res["start"].append(len(xyz))
    n = len(xyz)
    coords = [np.array([p.tolist() for p in xyz], dtype=np.float32).reshape(n, 3)]
    for k in range(1, len(st)):
        other = [atom.pos.tolist() for chain in st[k] for residue in chain for atom in residue]
        if len(other) == n:
            coords.append(np.array(other, dtype=np.float32))
        else:
            log.info("model %d of %s has other atoms than the first; not kept", k + 1, st.name)
    sequences = [[gemmi.Entity.first_mon(m) for m in e.full_sequence] for e in st.entities]
    structure = Structure(
        uid=np.arange(n, dtype=np.int64),
        coords=np.stack(coords),
        element=np.array(atoms["element"], dtype=np.uint8),
        name=np.array(atoms["name"], dtype="<U4"),
        altloc=np.array(atoms["altloc"], dtype="<U1"),
        occupancy=np.array(atoms["occ"], dtype=np.float32),
        b=np.array(atoms["b"], dtype=np.float32),
        charge=np.array(atoms["charge"], dtype=np.int8),
        res_start=np.array(res["start"], dtype=np.int32),
        res_name=np.array(res["name"], dtype="<U5"),
        res_seq=np.array(res["seq"], dtype=np.int32),
        icode=np.array(res["icode"], dtype="<U1"),
        res_kind=np.array(res["kind"], dtype=np.uint8),
        ss=np.zeros(len(res["name"]), dtype=np.uint8),
        seqres=np.array(res["seqres"], dtype=np.int32),
        label_asym=np.array(res["asym"], dtype="<U4"),
        entity=np.array(res["entity"], dtype=np.int16),
        chain_start=np.array(chain_start, dtype=np.int32),
        chain_name=np.array(chain_names, dtype="<U4"),
        entity_id=np.array([e.name for e in st.entities], dtype="<U8"),
        entity_seq=np.array([m for s in sequences for m in s], dtype="<U5"),
        entity_seq_start=np.cumsum([0] + [len(s) for s in sequences]).astype(np.int32),
        bonds=np.zeros((0, 2), np.int32),
        order=np.zeros(0, np.int8),
        next_uid=n,
        meta=_meta(st),
    )
    names = {str(structure.res_name[r]) for r in np.flatnonzero(structure.res_kind == LIGAND)}
    found = components.from_block(block, names) if block is not None and names else {}
    if found:
        structure = replace(
            structure, meta={**structure.meta, components.META_KEY: components.to_meta(found)}
        )
    pairs, orders = bonds.perceive(
        structure.name,
        structure.element,
        structure.altloc,
        structure.coords[0],
        structure.res_start,
        structure.res_name,
        structure.res_kind,
        structure.chain_start,
        links=_declared_links(st, structure),
        library={name: c.template_bonds for name, c in found.items()},
    )
    structure = replace(structure, bonds=pairs, order=orders)
    ss = _ss_from_records(st, structure)
    if ss is None and (structure.res_kind == AMINO).any():
        ss = secondary.assign(structure)
    return structure if ss is None else replace(structure, ss=ss)


def _meta(st: gemmi.Structure) -> dict[str, str]:
    meta = {"id": st.name}
    for key, item in (("title", "_struct.title"), ("method", "_exptl.method")):
        if item in st.info:
            meta[key] = st.info[item]
    if st.resolution > 0:
        meta["resolution"] = f"{st.resolution:.2f}"
    return meta


def _residue_lookup(structure: Structure) -> dict[tuple[str, int, str], int]:
    chains = structure.chain_name[structure.chain_of]
    return {
        (str(c), int(s), str(i)): r
        for r, (c, s, i) in enumerate(zip(chains, structure.res_seq, structure.icode, strict=True))
    }


def _address(address: gemmi.AtomAddress) -> tuple[str, int, str]:
    seqid = address.res_id.seqid
    return address.chain_name, seqid.num, seqid.icode.strip()


def _declared_links(st: gemmi.Structure, structure: Structure) -> Iterator[tuple[int, int]]:
    links = [c for c in st.connections if c.type in _LINKS]
    if not links:
        return
    residues = _residue_lookup(structure)
    for link in links:
        ends = [_atom(structure, residues, p) for p in (link.partner1, link.partner2)]
        if ends[0] >= 0 and ends[1] >= 0:
            yield ends[0], ends[1]


def _atom(
    structure: Structure, residues: dict[tuple[str, int, str], int], address: gemmi.AtomAddress
) -> int:
    residue = residues.get(_address(address))
    if residue is None:
        return -1
    span = structure.atoms_of(residue)
    altloc = address.altloc if address.altloc != "\0" else ""
    for atom in range(span.start, span.stop):
        if structure.name[atom] == address.atom_name and structure.altloc[atom] in ("", altloc):
            return atom
    return -1


def _ss_from_records(st: gemmi.Structure, structure: Structure) -> np.ndarray | None:
    spans = [(h.start, h.end, HELIX) for h in st.helices]
    spans += [(s.start, s.end, STRAND) for sheet in st.sheets for s in sheet.strands]
    if not spans:
        return None
    residues = _residue_lookup(structure)
    ss = np.full(structure.n_residues, COIL, dtype=np.uint8)
    for start, end, kind in spans:
        first, last = residues.get(_address(start)), residues.get(_address(end))
        if first is not None and last is not None and first <= last:
            ss[first : last + 1] = kind
    return ss


# writing


def to_gemmi(structure: Structure, keep: np.ndarray | None = None) -> gemmi.Structure:
    """A gemmi structure of the atoms in ``keep`` (all by default), every model."""
    if keep is not None:
        structure = structure.subset(keep)
    st = gemmi.Structure()
    st.name = structure.meta.get("id", "chemsuite")
    for k in range(structure.n_models):
        st.add_model(_model(structure, k))
    for e, entity_id in enumerate(structure.entity_id):
        entity = gemmi.Entity(str(entity_id))
        subchains = sorted({str(a) for a in structure.label_asym[structure.entity == e]})
        entity.subchains = subchains
        kinds = structure.res_kind[structure.entity == e]
        entity.entity_type = _entity_type(kinds)
        if entity.entity_type == gemmi.EntityType.Polymer:
            entity.full_sequence = [str(m) for m in structure.sequence_of(e)]
        st.entities.append(entity)
    st.setup_cell_images()
    return st


def _entity_type(kinds: np.ndarray) -> gemmi.EntityType:
    from .model import NUCLEIC, WATER

    if len(kinds) and np.isin(kinds, (AMINO, NUCLEIC)).all():
        return gemmi.EntityType.Polymer
    if len(kinds) and (kinds == WATER).all():
        return gemmi.EntityType.Water
    return gemmi.EntityType.NonPolymer


def _model(structure: Structure, k: int) -> gemmi.Model:
    from .model import NUCLEIC

    model = gemmi.Model(k + 1)
    xyz = structure.coords[k]
    for c in range(structure.n_chains):
        chain = gemmi.Chain(str(structure.chain_name[c]))
        span = structure.residues_of(c)
        for r in range(span.start, span.stop):
            residue = gemmi.Residue()
            residue.name = str(structure.res_name[r])
            residue.seqid = gemmi.SeqId(int(structure.res_seq[r]), structure.icode[r] or " ")
            residue.subchain = str(structure.label_asym[r])
            polymer = structure.res_kind[r] in (AMINO, NUCLEIC)
            residue.het_flag = "A" if polymer else "H"
            if structure.seqres[r] >= 0:
                residue.label_seq = int(structure.seqres[r]) + 1
            for i in range(*structure.atoms_of(r).indices(structure.n_atoms)):
                residue.add_atom(_atom_of(structure, xyz, i))
            chain.add_residue(residue)
        model.add_chain(chain)
    return model


def _atom_of(structure: Structure, xyz: np.ndarray, i: int) -> gemmi.Atom:
    atom = gemmi.Atom()
    atom.name = str(structure.name[i])
    atom.element = gemmi.Element(int(structure.element[i]))
    atom.altloc = structure.altloc[i] or "\0"
    atom.occ = float(structure.occupancy[i])
    atom.b_iso = float(structure.b[i])
    atom.charge = int(structure.charge[i])
    atom.pos = gemmi.Position(*(float(v) for v in xyz[i]))
    atom.serial = i + 1
    return atom


def to_text(structure: Structure, fmt: str, keep: np.ndarray | None = None) -> str:
    st = to_gemmi(structure, keep)
    if fmt == PDB:
        return st.make_pdb_string()
    return st.make_mmcif_document().as_string()


#: chain IDs a renamed chain takes, in this order (one character, so PDB files hold them)
CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def to_text_many(
    parts: Sequence[tuple[Structure, np.ndarray | None]], fmt: str
) -> tuple[str, dict[tuple[int, str], str]]:
    """One file of several structures' atoms (``keep`` masks, None for all): the first model
    of each, one after the other. A chain whose ID an earlier structure took is renamed; the
    renames come back as ``{(part, old ID): new ID}``. A single structure keeps every model."""
    if len(parts) == 1:
        return to_text(parts[0][0], fmt, parts[0][1]), {}
    st = gemmi.Structure()
    st.name = "chemsuite"
    model = gemmi.Model(1)
    renamed: dict[tuple[int, str], str] = {}
    used: set[str] = set()
    for k, (structure, keep) in enumerate(parts):
        for chain in to_gemmi(structure, keep)[0]:
            name = chain.name
            if name in used:
                name = next((c for c in CHAIN_IDS if c not in used), f"{len(used)}")
                renamed[(k, chain.name)] = name
                chain.name = name
            used.add(name)
            model.add_chain(chain)
    st.add_model(model)
    st.add_entity_types(True)
    st.assign_subchains(force=True)  # label_asym IDs: one set for the file, not one per part
    st.ensure_entities()
    st.deduplicate_entities()
    st.setup_cell_images()
    text = st.make_pdb_string() if fmt == PDB else st.make_mmcif_document().as_string()
    return text, renamed


def write(structure: Structure, path: str | PathLike[str], keep: np.ndarray | None = None) -> None:
    """Write as PDB or mmCIF, by the file's suffix (every model, every alternate location)."""
    write_text(path, to_text(structure, format_of(path), keep))


def write_text(path: str | PathLike[str], text: str) -> None:
    """Write ``text`` through a temporary file renamed over ``path`` when complete, so an
    interrupted save leaves the old file as it was."""
    target = Path(path)
    part = target.with_name(target.name + ".part")
    try:
        part.write_text(text, encoding="utf-8")
        os.replace(part, target)
    finally:
        part.unlink(missing_ok=True)


def format_of(path: str | PathLike[str]) -> str:
    """PDB or MMCIF by the file's suffix (``.gz`` aside); MMCIF for any other suffix."""
    return _suffix_format(path) or MMCIF


def readable(paths: Iterable[str | PathLike[str]]) -> list[Path]:
    """The paths that look like structure files."""
    return [Path(path) for path in paths if _suffix_format(path) is not None]
