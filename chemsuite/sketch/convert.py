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

"""The sketch to and from RDKit, molfiles and SMILES.

RDKit always reads a sketch the way it reads any molfile: the sketch is written as a V2000
block exactly as drawn (no hydrogen counts, wedges as drawn) and parsed with
``Chem.MolFromMolBlock``. A drawing is therefore valid here exactly when the table, search and
descriptors can use it (docs/sketcher/valence-review.md), and double-bond geometry and wedges
are interpreted by RDKit's molfile reader, not by code of our own.
"""

from __future__ import annotations

import math
import re
from statistics import median

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor

from chemsuite.core.mols import mol_from_text
from chemsuite.sketch import ctfile
from chemsuite.sketch.model import (
    ATOM_LIST,
    QUERY_ATOMS,
    RGROUP,
    Atom,
    Bond,
    Group,
    Order,
    Sketch,
    Stereo,
    parse_enhanced,
)

#: Molfile coordinates per sketch bond length. RDKit draws bonds 1.5 long; files keep that.
FILE_BOND_LENGTH = 1.5
MAX_V2000_ATOMS = 999


class SketchError(ValueError):
    """The drawing is not a molecule RDKit accepts."""


_RADICAL_CODES = {1: 2, 2: 3}  # unpaired electrons -> molfile RAD (doublet, triplet)
_STEREO_CODES = {Stereo.NONE: 0, Stereo.WEDGE: 1, Stereo.HASH: 6, Stereo.EITHER: 4}
_BOND_DIRS = {
    Chem.BondDir.BEGINWEDGE: Stereo.WEDGE,
    Chem.BondDir.BEGINDASH: Stereo.HASH,
    Chem.BondDir.UNKNOWN: Stereo.EITHER,
}


def atom_order(sketch: Sketch) -> list[int]:
    """Sketch atom ids in molfile order (atom *i* in RDKit is ``atom_order(sketch)[i]``)."""
    return sorted(sketch.atoms)


def drawn_molblock(sketch: Sketch) -> str:
    """The sketch as a V2000 molfile, exactly as drawn: no hydrogen counts or valences.

    Raises SketchError for more than 999 atoms or bonds, the V2000 limit.
    """
    order = atom_order(sketch)
    index = {aid: i + 1 for i, aid in enumerate(order)}
    bonds = sorted(sketch.bonds.items())
    if len(order) > MAX_V2000_ATOMS or len(bonds) > MAX_V2000_ATOMS:
        raise SketchError("The drawing has more than 999 atoms or bonds.")
    lines = ["", "  chemsuite sketch", ""]
    lines.append(f"{len(order):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000")
    charges, isotopes, radicals, rgroups, lists = [], [], [], [], []
    for aid in order:
        atom = sketch.atoms[aid]
        x, y = atom.x * FILE_BOND_LENGTH, atom.y * FILE_BOND_LENGTH
        # ccc sss hhh bbb vvv HHH rrr iii mmm nnn eee: mmm is the atom-atom map (CTfile)
        lines.append(
            f"{x:10.4f}{y:10.4f}{0:10.4f} {atom.symbol:<3} 0  0  0  0  0  0  0  0  0"
            f"{atom.mapping:3d}  0  0"
        )
        if atom.charge:
            charges.append((index[aid], atom.charge))
        if atom.isotope and not atom.is_query:
            isotopes.append((index[aid], atom.isotope))
        if atom.symbol == RGROUP:
            rgroups.append((index[aid], atom.rgroup or 1))
        if atom.symbol == ATOM_LIST:
            lists.append((index[aid], atom))
        if atom.radical:
            radicals.append((index[aid], _RADICAL_CODES.get(atom.radical, 3)))
    for _, bond in bonds:
        stereo = _STEREO_CODES[bond.stereo] if bond.order is Order.SINGLE else 0
        if bond.order is Order.DOUBLE and bond.stereo is Stereo.EITHER:
            stereo = 3  # CTfile: double bond, cis or trans (either)
        line = f"{index[bond.begin]:3d}{index[bond.end]:3d}{bond.order.value:3d}{stereo:3d}"
        if bond.topology:  # CTfile bond block: 111222tttsssxxxrrr, rrr the topology
            line += f"  0{bond.topology:3d}"
        lines.append(line)
    for tag, values in (("CHG", charges), ("ISO", isotopes), ("RAD", radicals)):
        for start in range(0, len(values), 8):  # at most eight entries per line
            chunk = values[start : start + 8]
            body = "".join(f" {i:3d} {v:3d}" for i, v in chunk)
            lines.append(f"M  {tag}{len(chunk):3d}{body}")
    attach = [(index[aid], sketch.atoms[aid].attach) for aid in order if sketch.atoms[aid].attach]
    for start in range(0, len(attach), 8):
        chunk = attach[start : start + 8]
        lines.append(f"M  APO{len(chunk):3d}" + "".join(f" {i:3d} {v:3d}" for i, v in chunk))
    for start in range(0, len(rgroups), 8):
        chunk = rgroups[start : start + 8]
        lines.append(f"M  RGP{len(chunk):3d}" + "".join(f" {i:3d} {n:3d}" for i, n in chunk))
    for i, atom in lists:  # CTfile: M  ALS aaa nn e [L1 ] [L2 ] ...  (e: T = NOT list)
        elements = atom.atom_list[:16]
        flag = "T" if atom.list_not else "F"
        lines.append(
            f"M  ALS {i:3d}{len(elements):3d} {flag} " + "".join(f"{e:<4}" for e in elements)
        )
    lines.extend(_sgroup_lines(sketch, index, {bid: i + 1 for i, (bid, _) in enumerate(bonds)}))
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def _sgroup_lines(sketch: Sketch, atom_index: dict[int, int], bond_index: dict[int, int]):
    """Abbreviations as CTfile superatom (``SUP``) S-groups: atoms, crossing bonds, label,
    and ``SDS EXP`` for those shown expanded."""
    groups = sorted(sketch.groups.items())
    out: list[str] = []
    for start in range(0, len(groups), 8):
        chunk = groups[start : start + 8]
        body = "".join(f" {start + k + 1:3d} SUP" for k in range(len(chunk)))
        out.append(f"M  STY{len(chunk):3d}{body}")
    for n, (_, group) in enumerate(groups, start=1):
        atoms = sorted(atom_index[a] for a in group.atoms)
        for start in range(0, len(atoms), 15):
            chunk = atoms[start : start + 15]
            out.append(f"M  SAL {n:3d}{len(chunk):3d}" + "".join(f" {a:3d}" for a in chunk))
        crossing = sorted(
            bond_index[bid]
            for bid, b in sketch.bonds.items()
            if (b.begin in group.atoms) != (b.end in group.atoms)
        )
        for start in range(0, len(crossing), 15):
            chunk = crossing[start : start + 15]
            out.append(f"M  SBL {n:3d}{len(chunk):3d}" + "".join(f" {b:3d}" for b in chunk))
        out.append(f"M  SMT {n:3d} {group.label}")
    expanded = [n for n, (_, g) in enumerate(groups, start=1) if g.expanded]
    for start in range(0, len(expanded), 15):
        chunk = expanded[start : start + 15]
        out.append(f"M  SDS EXP{len(chunk):3d}" + "".join(f" {n:3d}" for n in chunk))
    return out


def to_rdkit(sketch: Sketch, *, sanitize: bool = True) -> Chem.Mol:
    """The sketch read by RDKit, drawn hydrogens kept, atoms in ``atom_order``.

    Raises SketchError when RDKit rejects it (sanitisation failed).
    """
    with rdBase.BlockLogs():
        mol = Chem.MolFromMolBlock(drawn_molblock(sketch), sanitize=sanitize, removeHs=False)
    if mol is None:
        raise SketchError("RDKit cannot read this drawing.")
    if any(a.enhanced for a in sketch.atoms.values()):
        mol = _with_stereo_groups(sketch, mol)
    return mol


_GROUP_TYPES = {
    "abs": Chem.StereoGroupType.STEREO_ABSOLUTE,
    "and": Chem.StereoGroupType.STEREO_AND,
    "or": Chem.StereoGroupType.STEREO_OR,
}
_GROUP_KINDS = {v: k for k, v in _GROUP_TYPES.items()}


def _with_stereo_groups(sketch: Sketch, mol: Chem.Mol) -> Chem.Mol:
    """*mol* with the drawing's enhanced stereo groups (V3000 ``MDLV30/STEABS``, ``STERAC``,
    ``STEREL``). Only atoms RDKit reads as stereocentres can belong to a group."""
    members: dict[tuple[str, int], list[int]] = {}
    for i, aid in enumerate(atom_order(sketch)):
        parsed = parse_enhanced(sketch.atoms[aid].enhanced)
        if parsed and mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            members.setdefault(parsed, []).append(i)
    if not members:
        return mol
    rw = Chem.RWMol(mol)
    groups = [
        Chem.CreateStereoGroup(_GROUP_TYPES[kind], rw, atoms, readId=n)
        for (kind, n), atoms in sorted(members.items())
    ]
    rw.SetStereoGroups(groups)
    return rw.GetMol()


def enhanced_of(mol: Chem.Mol) -> dict[int, str]:
    """Atom index -> ``Atom.enhanced`` value for each atom in one of *mol*'s stereo groups."""
    out: dict[int, str] = {}
    counts = {"and": 0, "or": 0}
    for group in mol.GetStereoGroups():
        kind = _GROUP_KINDS.get(group.GetGroupType())
        if kind is None:
            continue
        if kind == "abs":
            value = "abs"
        else:
            counts[kind] += 1
            value = f"{kind}{group.GetReadId() or group.GetWriteId() or counts[kind]}"
        for atom in group.GetAtoms():
            out[atom.GetIdx()] = value
    return out


def to_mol(sketch: Sketch) -> Chem.Mol:
    """The molecule as the table would store it: sanitised, drawn hydrogens removed."""
    mol = to_rdkit(sketch)
    with rdBase.BlockLogs():
        return Chem.RemoveHs(mol, updateExplicitCount=True)


def to_smiles(sketch: Sketch) -> str:
    """Canonical isomeric SMILES; raises SketchError for an invalid drawing."""
    if sketch.is_empty():
        return ""
    return Chem.MolToSmiles(to_mol(sketch))


def to_molblock(sketch: Sketch) -> str:
    """A molfile for other programs, with the wedges as drawn.

    Bonds are Kekulé. Any atom whose hydrogen count a CTfile reader would work out differently
    from RDKit carries an explicit valence (for example 6 on the phosphorus of PF₆⁻), so other
    programs read the same structure. A drawing RDKit rejects is written exactly as drawn, so a
    copy never fails.
    """
    try:
        mol = to_rdkit(sketch)
    except SketchError:
        return drawn_molblock(sketch)
    Chem.ReapplyMolBlockWedging(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    _add_attach_points(mol)
    v3000 = (
        mol.GetNumAtoms() > MAX_V2000_ATOMS
        or mol.GetNumBonds() > MAX_V2000_ATOMS
        or bool(mol.GetStereoGroups())  # enhanced stereo exists only in V3000
    )
    if v3000:
        return _v3000_with_valences(mol, [sketch.atoms[a].attach for a in atom_order(sketch)])
    lines = Chem.MolToMolBlock(mol, kekulize=False).split("\n")
    for atom in mol.GetAtoms():
        valence = _ctfile_valence(atom)
        if valence is not None:
            line = lines[4 + atom.GetIdx()]
            lines[4 + atom.GetIdx()] = f"{line[:48]}{valence:3d}{line[51:]}"
    if any(a.attach for a in sketch.atoms.values()):  # RDKit's writer drops M  APO
        ours = [ln for ln in drawn_molblock(sketch).split("\n") if ln.startswith("M  APO")]
        lines = [ln for ln in lines if not ln.startswith("M  APO")]
        lines[lines.index("M  END") : lines.index("M  END")] = ours
    if any(a.symbol == ATOM_LIST for a in sketch.atoms.values()):
        # RDKit misaligns NOT lists (``M  ALS 4   2  T``); write the lists ourselves
        ours = [ln for ln in drawn_molblock(sketch).split("\n") if ln.startswith("M  ALS")]
        lines = [ln for ln in lines if not ln.startswith("M  ALS")]
        lines[lines.index("M  END") : lines.index("M  END")] = ours
    # RDKit adds "V  aaa *" value lines for dummy atoms; they carry nothing
    return "\n".join(ln for ln in lines if not (ln.startswith("V  ") and ln.endswith(" *")))


def _v3000_with_valences(mol: Chem.Mol, attach: list[int] | None = None) -> str:
    """A V3000 molfile with ``VAL=`` on atoms that need an explicit valence (as in V2000) and
    ``ATTCHPT=`` on attachment points."""
    lines = Chem.MolToMolBlock(mol, kekulize=False, forceV3000=True).split("\n")
    start = lines.index("M  V30 BEGIN ATOM") + 1
    for atom in mol.GetAtoms():
        valence = _ctfile_valence(atom)
        line = lines[start + atom.GetIdx()]
        point = attach[atom.GetIdx()] if attach else 0
        if point and " ATTCHPT=" not in line:
            line = f"{line} ATTCHPT={-1 if point == 3 else point}"  # CTfile: -1 means both
        if valence is not None and " VAL=" not in line:
            line = f"{line} VAL={-1 if valence == 15 else valence}"
        if len(line) > 80:  # CTfile V3000: at most 80 characters, "-" continues a line
            line = f"{line[:79]}-\nM  V30 {line[79:]}"
        lines[start + atom.GetIdx()] = line
    return "\n".join(lines)


def _add_attach_points(mol: Chem.Mol) -> None:
    """Give each superatom its attachment point (``SAP``): the group atom bonded outside."""
    for k in range(len(Chem.GetMolSubstanceGroups(mol))):
        group = Chem.GetMolSubstanceGroupWithIdx(mol, k)  # a reference; the list holds copies
        if group.GetProp("TYPE") != "SUP" or list(group.GetAttachPoints()):
            continue
        members = set(group.GetAtoms())
        for index in group.GetAtoms():
            outside = [
                n.GetIdx() for n in mol.GetAtomWithIdx(index).GetNeighbors()
                if n.GetIdx() not in members
            ]  # fmt: skip
            if outside:
                group.AddAttachPoint(index, outside[0], "1 ")  # CTfile ids are two wide
                break


def _ctfile_valence(atom: Chem.Atom) -> int | None:
    """The molfile valence field this atom needs, or None if CTfile readers agree with RDKit.

    The field holds the total valence (bond orders plus hydrogens); 15 means zero.
    """
    if atom.GetAtomicNum() == 0:
        return None  # query atoms and R-groups have no valence of their own
    bonds = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
    hydrogens = atom.GetTotalNumHs()
    symbol = atom.GetSymbol()
    if symbol in ctfile.ELEMENTS:
        theirs = ctfile.default_hydrogens(
            symbol, atom.GetFormalCharge(), bonds, atom.GetNumRadicalElectrons()
        )
        if theirs == hydrogens:
            return None
    total = bonds + hydrogens
    return total if total else 15


def from_mol(mol: Chem.Mol) -> Sketch:
    """A sketch of a molecule: its 2D coordinates (computed if it has none), Kekulé bonds, and
    wedges as in its molfile, or chosen by RDKit."""
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() and not mol.GetConformer().Is3D():
        return _sketch_of(mol, set())
    # A molecule without a drawing (SMILES, reaction SMILES) gives no E/Z where its text gives
    # none; the layout must not pick one. Such bonds are drawn "either" (CTfile stereo 3), as
    # RDKit's molfile writer marks them.
    unknown = _undefined_double_bonds(mol)
    first = None
    for lay_out in _LAYOUTS:
        laid = Chem.Mol(mol)
        laid.RemoveAllConformers()
        if not lay_out(laid):
            continue
        sketch = _sketch_of(laid, unknown)
        first = first or sketch
        # RDKit places wedges for the layout it is given; on a few bridged layouts they read
        # back as contradictory and the centre is lost. Then the other layout is tried.
        if _same_stereo(mol, sketch):
            return sketch
    return first if first is not None else _sketch_of(mol, unknown)


def _layout_preferred(mol: Chem.Mol) -> bool:
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    # CoordGen (preferred app-wide, core/depict.py) places a lone atom such as [H+] at NaN
    return not any(math.isnan(conf.GetAtomPosition(i).x) for i in range(mol.GetNumAtoms()))


def _layout_rdkit(mol: Chem.Mol) -> bool:
    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(False)
    try:
        rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)
    return True


_LAYOUTS = (_layout_preferred, _layout_rdkit)


def _same_stereo(mol: Chem.Mol, sketch: Sketch) -> bool:
    """Whether the drawing reads back as *mol* (isomeric SMILES); True when either cannot be
    read, so the check never blocks a drawing."""
    with rdBase.BlockLogs():
        try:
            want = Chem.MolToSmiles(Chem.RemoveHs(mol))
            return Chem.MolToSmiles(to_mol(sketch)) == want
        except (ValueError, RuntimeError, SketchError):
            return True


def _sketch_of(mol: Chem.Mol, unknown: set[int]) -> Sketch:
    """The sketch of *mol* as laid out in its 2D conformer."""
    conf = mol.GetConformer()
    if any(b.HasProp("_MolFileBondStereo") for b in mol.GetBonds()):
        Chem.ReapplyMolBlockWedging(mol)
    else:
        Chem.WedgeMolBonds(mol, conf)
    try:
        with rdBase.BlockLogs():
            Chem.Kekulize(mol, clearAromaticFlags=True)
    except (ValueError, RuntimeError):
        pass  # keep aromatic bonds that cannot be made alternating
    lengths = [
        (
            conf.GetAtomPosition(b.GetBeginAtomIdx()) - conf.GetAtomPosition(b.GetEndAtomIdx())
        ).Length()
        for b in mol.GetBonds()
    ]
    scale = median(lengths) if lengths and median(lengths) > 1e-6 else FILE_BOND_LENGTH
    sketch = Sketch()
    ids: list[int] = []
    groups = enhanced_of(mol)
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        aid = sketch.new_id()
        sketch.atoms[aid] = Atom(
            p.x / scale,
            p.y / scale,
            charge=atom.GetFormalCharge(),
            radical=atom.GetNumRadicalElectrons(),
            enhanced=groups.get(atom.GetIdx(), ""),
            attach=_attach_point(atom),
            mapping=atom.GetAtomMapNum() if atom.GetAtomicNum() else 0,
            **_atom_kind(atom),
        )
        ids.append(aid)
    for bond in mol.GetBonds():
        kind = bond.GetBondType()
        order = {
            Chem.BondType.DOUBLE: Order.DOUBLE,
            Chem.BondType.TRIPLE: Order.TRIPLE,
            Chem.BondType.AROMATIC: Order.AROMATIC,
        }.get(kind, Order.SINGLE)
        file_type = bond.GetIntProp("_MolFileBondType") if bond.HasProp("_MolFileBondType") else 0
        if 5 <= file_type <= 8:
            order = Order(file_type)
        topology = _topology(bond)
        stereo = _BOND_DIRS.get(bond.GetBondDir(), Stereo.NONE)
        if bond.HasProp("_MolFileBondStereo") and bond.GetIntProp("_MolFileBondStereo") == 4:
            stereo = Stereo.EITHER
        if order is Order.DOUBLE:  # wedge directions mean nothing on a double bond
            wavy_beside = any(
                other.GetBondDir() == Chem.BondDir.UNKNOWN
                for end in (bond.GetBeginAtom(), bond.GetEndAtom())
                for other in end.GetBonds()
            )  # a wavy single bond beside it already says "unknown" (RDKit writes both)
            either = not wavy_beside and (
                (bond.HasProp("_MolFileBondStereo") and bond.GetIntProp("_MolFileBondStereo") == 3)
                or bond.GetStereo() == Chem.BondStereo.STEREOANY
            )
            stereo = Stereo.EITHER if either or bond.GetIdx() in unknown else Stereo.NONE
        sketch.bonds[sketch.new_id()] = Bond(
            ids[bond.GetBeginAtomIdx()], ids[bond.GetEndAtomIdx()], order, stereo, topology
        )
    for group in Chem.GetMolSubstanceGroups(mol):
        if group.GetProp("TYPE") != "SUP" or not group.HasProp("LABEL"):
            continue
        members = list(group.GetAtoms())
        if not members:
            continue
        inside = set(members)
        attach = next(
            (
                i
                for i in members
                if any(n.GetIdx() not in inside for n in mol.GetAtomWithIdx(i).GetNeighbors())
            ),
            members[0],
        )
        expanded = group.HasProp("ESTATE") and group.GetProp("ESTATE") == "E"
        sketch.groups[sketch.new_id()] = Group(
            group.GetProp("LABEL"), frozenset(ids[i] for i in members), ids[attach], expanded
        )
    return sketch


def _undefined_double_bonds(mol: Chem.Mol) -> set[int]:
    """Double bonds RDKit finds could be E or Z but are given as neither. A double bond at a
    cumulated centre (allene, azide) is linear there and is left out."""
    out = set()
    probe = Chem.Mol(mol)
    for atom in probe.GetAtoms():  # atom maps make otherwise equal groups distinct (reactions)
        if atom.GetAtomMapNum():
            atom.SetIsotope(1000 + atom.GetAtomMapNum())
    with rdBase.BlockLogs():
        try:
            found = Chem.FindPotentialStereo(probe)
        except (ValueError, RuntimeError):
            return out
    for info in found:
        if (
            info.type == Chem.StereoType.Bond_Double
            and info.specified == Chem.StereoSpecified.Unspecified
        ):
            bond = mol.GetBondWithIdx(info.centeredOn)
            ends = (bond.GetBeginAtom(), bond.GetEndAtom())
            if not any(
                sum(b.GetBondType() == Chem.BondType.DOUBLE for b in end.GetBonds()) > 1
                for end in ends
            ):
                out.add(info.centeredOn)
    return out


def _attach_point(atom: Chem.Atom) -> int:
    """``Atom.attach`` from RDKit's ``molAttchpt`` (it reads "both" as -1, as V3000 writes it)."""
    if not atom.HasProp("molAttchpt"):
        return 0
    value = atom.GetIntProp("molAttchpt")
    return 3 if value == -1 else value if value in (1, 2) else 0


def _topology(bond: Chem.Bond) -> int:
    """CTfile topology of a query bond read from a molfile: 1 ring, 2 chain, 0 either."""
    if not bond.HasQuery():
        return 0
    text = bond.DescribeQuery()
    if "BondInRing 1 != val" in text:
        return 2
    if "BondInRing 1 = val" in text:
        return 1
    return 0


_ATOMIC_NUM = re.compile(r"AtomAtomicNum (\d+) = val")


def _atom_kind(atom: Chem.Atom) -> dict:
    """``symbol`` and the fields that go with it, for an element, query atom or R-group.

    Molfile query atoms come back from RDKit with their CTfile symbol as the query type (A, Q,
    X, M and their H forms), atom lists as an OR of atomic numbers (``not`` for a NOT list) and
    R-groups as dummies with ``_MolFileRLabel``. SMILES dummies follow RDKit's R-group
    convention: ``[2*]`` and ``[*:2]`` are R2, a plain ``*`` is "any atom".
    """
    if atom.GetAtomicNum() != 0:
        return {"symbol": atom.GetSymbol(), "isotope": atom.GetIsotope()}
    if atom.HasQuery():
        kind = atom.GetQueryType()
        if kind in QUERY_ATOMS:
            return {"symbol": kind}
        if atom.HasProp("_MolFileAtomQuery"):
            text = atom.DescribeQuery()
            numbers = [int(n) for n in _ATOMIC_NUM.findall(text)]
            table = Chem.GetPeriodicTable()
            if numbers and "AtomAtomicNum" in text:
                return {
                    "symbol": ATOM_LIST,
                    "atom_list": tuple(table.GetElementSymbol(n) for n in numbers),
                    "list_not": text.lstrip().startswith("not"),
                }
    if atom.HasProp("_MolFileRLabel"):
        return {"symbol": RGROUP, "rgroup": atom.GetIntProp("_MolFileRLabel")}
    number = atom.GetIsotope() or atom.GetAtomMapNum()
    if number:
        return {"symbol": RGROUP, "rgroup": number}
    return {"symbol": "*"}


def from_text(text: str) -> Sketch | None:
    """A sketch from SMILES, a molfile or an InChI; None when the text is not a structure."""
    mol = mol_from_text(text)
    return None if mol is None else from_mol(mol)


def placed(sketch: Sketch, at: tuple[float, float]) -> Sketch:
    """A copy of the sketch moved so the centre of its box is at ``at``."""
    x0, y0, x1, y1 = sketch.bounds()
    dx, dy = at[0] - (x0 + x1) / 2, at[1] - (y0 + y1) / 2
    out = sketch.copy()
    out.atoms = {aid: a.moved(a.x + dx, a.y + dy) for aid, a in sketch.atoms.items()}
    out.arrows = {rid: r.moved(dx, dy) for rid, r in sketch.arrows.items()}
    out.pluses = {pid: p.moved(dx, dy) for pid, p in sketch.pluses.items()}
    return out


def bond_length_spread(sketch: Sketch) -> float:
    """Largest relative deviation of a bond length from 1 (a layout check for tests)."""
    worst = 0.0
    for bond in sketch.bonds.values():
        a, b = sketch.atoms[bond.begin], sketch.atoms[bond.end]
        worst = max(worst, abs(math.hypot(a.x - b.x, a.y - b.y) - 1.0))
    return worst
