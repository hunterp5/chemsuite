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
"""Coloring a surface by electrostatics: Coulomb's law with a distance-dependent dielectric.

**Charges** (:func:`charges`). Amino acids and nucleotides take Amber's RESP charges (Cornell
et al., JACS 1995, kept by ff14SB) from PDB2PQR's ``AMBER.DAT`` (BSD-3; ``data/``). The
variant is the one the hydrogens show, as ChimeraX reads them: His with HD1 and HE2 is HIP,
with HD1 alone HID, else HIE; Asp with HD2 is ASH, Glu with HE2 GLH, Lys missing a HZ is LYN;
a Cys bonded to another is CYX, a prepared one without HG CYM. A chain's first amino acid
takes the N-terminal variant, its last (or one with OXT) the C-terminal one. Each hydrogen
of the variant that the structure lacks has its charge added to the heavy atom it would be
bonded to, so an unprepared protein gets united-atom charges and every residue keeps its net
charge. Ligands, glycans, additives and residues the table lacks take Gasteiger charges
(RDKit): a prepared one in the protomer Prepare chose (its formal charges and hydrogens),
else in the form its bond orders give, with its implicit hydrogens' charges on their heavy
atoms. An ion takes its formal charge (or its usual one when the file gives none); water
nothing.

**Potential** (:func:`potential`). At each surface vertex, moved 1.4 Å out along its normal
(as ChimeraX's ``coulombic`` does: at the surface itself the nearest atom would dominate),
``sum 332 q / (eps r)`` with ``eps = 4 r`` over the atoms within 15 Å, in kcal/(mol e).
Colored red at -10, white at 0 and blue at +10 (ChimeraX's default scale).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cache
from importlib import resources

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from scipy.spatial import cKDTree

from ..structure.model import AMINO, ION, NUCLEIC, WATER, Structure
from .mesh import Mesh

COULOMB = 332.0636  # kcal Å / (mol e^2)
OFFSET = 1.4  # Å out along the normal
CUTOFF = 15.0  # Å
SCALE = 10.0  # kcal/(mol e): the potential that is fully red (negative) or blue
RED = np.array([0.85, 0.12, 0.10], np.float32)
WHITE = np.array([1.0, 1.0, 1.0], np.float32)
BLUE = np.array([0.10, 0.25, 0.90], np.float32)
#: an ion's usual charge, by atomic number, when the file gives it none
ION_CHARGES = {
    3: 1, 11: 1, 19: 1, 37: 1, 55: 1,  # Li Na K Rb Cs
    12: 2, 20: 2, 25: 2, 26: 2, 27: 2, 28: 2, 29: 2, 30: 2, 38: 2, 48: 2, 56: 2, 80: 2,
    9: -1, 17: -1, 35: -1, 53: -1,  # F Cl Br I
}  # fmt: skip
#: residue names of RNA in a PDB file, and the table's
_RNA = {"A": "RA", "C": "RC", "G": "RG", "U": "RU"}
#: atom names that differ between the PDB and the table
_ALIASES = {"OP1": "O1P", "OP2": "O2P", "H5'": "H5'1", "H5''": "H5'2", "H2'": "H2'1",
            "H2''": "H2'2", "HO2'": "HO'2"}  # fmt: skip


@dataclass(frozen=True)
class Potential:
    """The potential at each vertex of *mesh* (the surface it was computed for)."""

    mesh: Mesh
    values: np.ndarray  # float32 (v,)


@cache
def table() -> dict[str, dict[str, float]]:
    """Residue name -> atom name -> charge (``data/amber_charges.json``)."""
    data = resources.files("chemsuite.macromolecule.look").joinpath("data/amber_charges.json")
    return json.loads(data.read_text(encoding="utf-8"))["charges"]


@cache
def _parents(residue: str) -> dict[str, str]:
    """For each hydrogen of a table residue, the heavy atom it is bonded to (by name)."""
    names = list(table()[residue])
    heavy = [n for n in names if not n.startswith("H")]
    out = {}
    for h in names:
        if h.startswith("H"):
            out[h] = _parent(h, heavy)
    return out


def _parent(hydrogen: str, heavy: list[str]) -> str:
    """The heavy atom a hydrogen is named after: HB2 on CB, HD21 on ND2, HG12 on CG1, HH11
    on NH1, H and H1-H3 on N, HH31 on CH3, H5'1 on C5', H61 on N6, HO'2 on O2', H5T on O5'."""
    rest = hydrogen[1:]
    if rest in ("5T", "3T"):  # a nucleic chain's terminal hydroxyl
        return f"O{rest[0]}'" if f"O{rest[0]}'" in heavy else heavy[0]
    if rest.startswith("O") and "'" in rest:  # a sugar's hydroxyl: HO'2 or HO2' on O2'
        oxygen = "O" + "".join(c for c in rest if c.isdigit()) + "'"
        if oxygen in heavy:
            return oxygen
    for suffix in (rest, rest[:-1]):
        found = [n for n in heavy if n[1:] == suffix]
        if len(found) == 1:
            return found[0]
        if found:  # a sugar's C5' over O5'; a base's amine N6 over C6; the backbone's N
            first = "C" if "'" in suffix else "N"
            preferred = [n for n in found if n.startswith(first)]
            return (preferred or found)[0]
    return heavy[0]


def charges(structure: Structure, active: np.ndarray | None = None) -> np.ndarray:
    """A partial charge per atom (float32; 0 for water and atoms not shown)."""
    if active is None:
        active = structure.active()
    out = np.zeros(structure.n_atoms, np.float32)
    kinds = structure.res_kind
    amino = np.flatnonzero(kinds == AMINO)
    chain = structure.chain_of
    first = {int(chain[r]): int(r) for r in amino[::-1]}
    last = {int(chain[r]): int(r) for r in amino}
    for residue in range(structure.n_residues):
        span = structure.atoms_of(residue)
        atoms = [i for i in range(span.start, span.stop) if active[i]]
        if not atoms:
            continue
        kind = int(kinds[residue])
        if kind == WATER:
            continue
        if kind == ION:
            for i in atoms:
                given = int(structure.charge[i])
                out[i] = given or ION_CHARGES.get(int(structure.element[i]), 0)
            continue
        variant = None
        if kind in (AMINO, NUCLEIC):
            c = int(chain[residue])
            variant = _variant(structure, residue, atoms, first.get(c) == residue,
                               last.get(c) == residue)  # fmt: skip
        if variant is not None:
            _from_table(structure, atoms, variant, out)
        else:  # a ligand, glycan or additive, or a residue the table lacks
            _gasteiger(structure, residue, active, out)
    return out


def _variant(
    structure: Structure, residue: int, atoms: list[int], n_term: bool, c_term: bool
) -> str | None:
    """The table's name for this residue, or None when it has none."""
    known = table()
    name = str(structure.res_name[residue]).upper()
    names = {str(structure.name[i]) for i in atoms}
    if structure.res_kind[residue] == NUCLEIC:
        base = _RNA.get(name, name)
        if base not in known:
            return None
        if n_term and base + "5" in known:
            return base + "5"
        return base + "3" if c_term and base + "3" in known else base
    prepared = any(structure.element[i] == 1 for i in atoms)
    if name == "HIS":
        name = "HIP" if {"HD1", "HE2"} <= names else "HID" if "HD1" in names else "HIE"
    elif name == "ASP" and "HD2" in names:
        name = "ASH"
    elif name == "GLU" and "HE2" in names:
        name = "GLH"
    elif name == "LYS" and prepared and len({"HZ1", "HZ2", "HZ3"} & names) == 2:
        name = "LYN"
    elif name == "CYS":
        if _bonded_to_sulfur(structure, atoms):
            name = "CYX"
        elif prepared and "HG" not in names:
            name = "CYM"
    if name not in known:
        return None
    c_term = c_term or "OXT" in names
    for end, flag in (("N", n_term), ("C", c_term)):
        if flag and end + name in known:
            return end + name
    return name


def _bonded_to_sulfur(structure: Structure, atoms: list[int]) -> bool:
    sulfur = [i for i in atoms if structure.name[i] == "SG"]
    if not sulfur:
        return False
    bonds = structure.bonds
    s = sulfur[0]
    partners = np.concatenate([bonds[bonds[:, 0] == s, 1], bonds[bonds[:, 1] == s, 0]])
    return bool((structure.element[partners] == 16).any())


def _from_table(structure: Structure, atoms: list[int], variant: str, out: np.ndarray) -> None:
    charges_of = table()[variant]
    by_name = {_ALIASES.get(str(structure.name[i]), str(structure.name[i])): i for i in atoms}
    for name, charge in charges_of.items():
        if name in by_name:
            out[by_name[name]] += charge
        elif name.startswith("H"):  # a hydrogen the structure lacks: on its heavy atom
            parent = _parents(variant)[name]
            if parent in by_name:
                out[by_name[parent]] += charge


def _gasteiger(structure: Structure, residue: int, active: np.ndarray, out: np.ndarray) -> None:
    from ..chem import ligands

    span = structure.atoms_of(residue)
    prepared = bool((structure.element[span] == 1).any())
    try:
        try:  # a prepared ligand in the protomer Prepare chose, its hydrogens explicit
            if not prepared:
                raise ligands.FormError("no hydrogens")
            mol = ligands.prepared_form(structure, residue).mol
        except ligands.FormError:
            mol = ligands.to_mol(structure, residue, active=active).mol
        with rdBase.BlockLogs():
            AllChem.ComputeGasteigerCharges(mol)
    except (ValueError, RuntimeError, Chem.MolSanitizeException):
        return
    uids = np.array([a.GetIntProp("uid") for a in mol.GetAtoms()], np.int64)
    index = structure.indices(uids)
    for atom, i in zip(mol.GetAtoms(), index, strict=True):
        if i < 0:
            continue
        charge = atom.GetDoubleProp("_GasteigerCharge")
        # the charge of all its implicit hydrogens together (RDKit's _GasteigerHCharge)
        hydrogens = (
            atom.GetDoubleProp("_GasteigerHCharge") if atom.HasProp("_GasteigerHCharge") else 0
        )
        total = charge + hydrogens
        if np.isfinite(total):
            out[i] = total


def potential(
    points: np.ndarray, centres: np.ndarray, charge: np.ndarray, cutoff: float = CUTOFF
) -> np.ndarray:
    """The potential (kcal/(mol e)) at *points* from *charge* at *centres*: 332 q / (4 r^2)
    over the charges within *cutoff* (r at least 1 Å), in blocks so memory stays small."""
    points = np.asarray(points, np.float64).reshape(-1, 3)
    out = np.zeros(len(points), np.float32)
    charged = np.flatnonzero(charge != 0)
    if not len(points) or not len(charged):
        return out
    atoms = cKDTree(np.asarray(centres, np.float64)[charged])
    q = np.asarray(charge, np.float64)[charged]
    block = 4096
    for start in range(0, len(points), block):
        part = cKDTree(points[start : start + block])
        pairs = atoms.sparse_distance_matrix(part, cutoff, output_type="coo_matrix")
        r = np.maximum(pairs.data, 1.0)
        terms = COULOMB * q[pairs.row] / (4.0 * r * r)
        out[start : start + part.n] = np.bincount(pairs.col, terms, minlength=part.n)
    return out


def surface_potential(
    structure: Structure, mesh: Mesh, model: int = 0, active: np.ndarray | None = None
) -> Potential:
    """The potential at each vertex of *mesh*, 1.4 Å out along its normal, from every atom
    of the structure (not only the surface's)."""
    if not len(mesh):
        return Potential(mesh, np.zeros(0, np.float32))
    points = mesh.positions + OFFSET * mesh.normals
    values = potential(points, structure.xyz(model), charges(structure, active))
    return Potential(mesh, values)


def colour(values: np.ndarray) -> np.ndarray:
    """RGB (float32 (n, 3)) for potentials: red at -SCALE, white at 0, blue at +SCALE."""
    t = np.clip(np.asarray(values, np.float32) / SCALE, -1.0, 1.0)[:, None]
    toward = np.where(t < 0, RED, BLUE)
    return WHITE + np.abs(t) * (toward - WHITE)
