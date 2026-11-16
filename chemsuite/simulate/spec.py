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

"""What is simulated, with which force field, and how it is minimised.

A :class:`SystemSpec` is the topology of a prepared complex: every atom (hydrogens included)
with its uid, element, name and formal charge; every residue with the name of the form it is
in (``HID``, ``ASH``, ``CYX``, ``ACE``, ``HOH``, ``ZN``, or a ligand's own); every bond; the
disulfides; and each small molecule (ligand or additive) as a molblock with explicit
hydrogens, in the same atom order as its residue, with its partial charges once the worker has
computed them. Positions travel separately (a spec serves every model and every frame).

The JSON form (:meth:`SystemSpec.to_json`) is columns of plain values, versioned by
:data:`SCHEMA`, so the worker reads it with nothing but the standard library::

    {"schema": 1,
     "atoms": {"uid": [...], "element": [...], "name": [...], "formal": [...]},
     "residues": {"name": [...], "chain": [...], "seq": [...], "icode": [...],
                  "kind": [...], "label": [...], "first": [...], "count": [...]},
     "bonds": [[i, j, order], ...], "disulfides": [[i, j], ...],
     "molecules": [{"residue": r, "molblock": "...", "charges": [...],
                    "charge_model": "..."}, ...]}

Atom and residue references are indices into these columns. :meth:`SystemSpec.digest` names
the topology alone (charges left out), so stored charges and minimised models can be checked
against the structure they were made for.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, replace
from typing import Any

SCHEMA = 1

# residue kinds
PROTEIN = "protein"  # an amino acid, in an Amber form (HID, ASH, CYX…; N- or C-terminal)
CAP = "cap"  # ACE, NME
WATER = "water"
ION = "ion"
MOLECULE = "molecule"  # a ligand or an additive, parameterised as a small molecule
KINDS = (PROTEIN, CAP, WATER, ION, MOLECULE)

# force fields
FF14SB, FF19SB = "ff14SB", "ff19SB"
SAGE = "openff-2.3.0"
TIP3P, OPC = "tip3p", "opc"
GBN2, EXPLICIT, VACUUM = "gbn2", "explicit", "vacuum"
#: each protein force field with the water model it was fitted with
WATERS = {FF14SB: TIP3P, FF19SB: OPC}

# what a minimisation stage lets move
HYDROGENS = "hydrogens"  # heavy atoms massless, so held where they are
SIDE_CHAINS = "side chains"  # the protein backbone's heavy atoms restrained, the rest free
ALL = "all"  # every heavy atom restrained to where it started
MOVING = (HYDROGENS, SIDE_CHAINS, ALL)
KCAL_PER_A2 = 418.4  # kJ/mol/nm² in one kcal/mol/Å²


class SpecError(ValueError):
    """Saved or received JSON that is not a spec; the message says why."""


@dataclass(frozen=True)
class Residue:
    name: str  # the form's name: HID, ASH, CYX, NALA, CGLY, ACE, HOH, ZN, the ligand's own
    chain: str
    seq: int
    icode: str
    kind: str  # one of KINDS
    label: str  # as the viewer names it, ``A/HIS 57``
    first: int  # its first atom
    count: int  # its number of atoms


@dataclass(frozen=True)
class Molecule:
    """A ligand or additive, to be parameterised as a small molecule."""

    residue: int
    molblock: str  # explicit hydrogens, the atoms in the residue's order
    charges: tuple[float, ...] = ()  # partial charges (e) once computed, else empty
    charge_model: str = ""  # what computed them (``openff-gnn-am1bcc-1.0.0``)


@dataclass(frozen=True)
class SystemSpec:
    uid: tuple[int, ...]
    element: tuple[int, ...]
    name: tuple[str, ...]
    formal: tuple[int, ...]
    residues: tuple[Residue, ...]
    bonds: tuple[tuple[int, int, int], ...]  # atom indices (i < j) and the bond order
    disulfides: tuple[tuple[int, int], ...] = ()  # the two SG atoms
    molecules: tuple[Molecule, ...] = ()
    schema: int = SCHEMA

    @property
    def n_atoms(self) -> int:
        return len(self.uid)

    @property
    def net_charge(self) -> int:
        return sum(self.formal)

    @property
    def charged(self) -> bool:
        """Whether every small molecule has its partial charges."""
        return all(m.charges for m in self.molecules)

    def residue_of(self) -> list[int]:
        """The residue of every atom."""
        out = [0] * self.n_atoms
        for r, residue in enumerate(self.residues):
            out[residue.first : residue.first + residue.count] = [r] * residue.count
        return out

    def with_charges(self, residue: int, charges: tuple[float, ...], model: str) -> SystemSpec:
        """The spec with the small molecule in *residue* given its partial *charges*."""
        found = [m for m in self.molecules if m.residue == residue]
        if not found:
            raise SpecError(f"No small molecule in residue {residue}")
        if len(charges) != self.residues[residue].count:
            raise SpecError(
                f"{self.residues[residue].label} has "
                f"{self.residues[residue].count} atoms, not {len(charges)}"
            )
        molecules = tuple(replace(m, charges=tuple(float(q) for q in charges), charge_model=model)
                          if m.residue == residue else m for m in self.molecules)  # fmt: skip
        return replace(self, molecules=molecules)

    def topology(self) -> dict[str, Any]:
        """The JSON form less the charges: what :meth:`digest` names."""
        out = self.to_json()
        out["molecules"] = [{"residue": m["residue"], "molblock": m["molblock"]}
                            for m in out["molecules"]]  # fmt: skip
        return out

    def digest(self) -> str:
        text = json.dumps(self.topology(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    def to_json(self) -> dict[str, Any]:
        columns = ("name", "chain", "seq", "icode", "kind", "label", "first", "count")
        return {
            "schema": self.schema,
            "atoms": {"uid": list(self.uid), "element": list(self.element),
                      "name": list(self.name), "formal": list(self.formal)},
            "residues": {c: [_field(r, c) for r in self.residues] for c in columns},
            "bonds": [list(b) for b in self.bonds],
            "disulfides": [list(d) for d in self.disulfides],
            "molecules": [{"residue": m.residue, "molblock": m.molblock,
                           "charges": list(m.charges), "charge_model": m.charge_model}
                          for m in self.molecules],
        }  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> SystemSpec:
        """The spec *saved* holds; :class:`SpecError` when it is not one this chemsuite reads."""
        if not isinstance(saved, dict):
            raise SpecError("A system spec is a JSON object")
        if saved.get("schema") != SCHEMA:
            raise SpecError(f"System spec schema {saved.get('schema')!r}; this reads {SCHEMA}")
        atoms = _object(saved, "atoms")
        uid = _column(atoms, "uid", int)
        n = len(uid)
        element = _column(atoms, "element", int, n)
        name = _column(atoms, "name", str, n)
        formal = _column(atoms, "formal", int, n)
        rows = _object(saved, "residues")
        names = _column(rows, "name", str)
        k = len(names)
        residues = tuple(
            Residue(*values)
            for values in zip(names, _column(rows, "chain", str, k), _column(rows, "seq", int, k),
                              _column(rows, "icode", str, k), _column(rows, "kind", str, k),
                              _column(rows, "label", str, k), _column(rows, "first", int, k),
                              _column(rows, "count", int, k), strict=True)
        )  # fmt: skip
        at = 0
        for r in residues:
            if r.kind not in KINDS or r.first != at or r.count < 1:
                raise SpecError(f"Residue {r.label} is not where its atoms are")
            at += r.count
        if at != n:
            raise SpecError("The residues do not cover the atoms")
        bonds = tuple(_tuple(b, 3, n, "bond") for b in _list(saved, "bonds"))
        disulfides = tuple(_tuple(d, 2, n, "disulfide") for d in _list(saved, "disulfides"))
        molecules = tuple(_molecule(m, residues) for m in _list(saved, "molecules"))
        return cls(tuple(uid), tuple(element), tuple(name), tuple(formal), residues,
                   bonds, disulfides, molecules)  # fmt: skip


def _field(residue: Residue, column: str) -> Any:
    return {"name": residue.name, "chain": residue.chain, "seq": residue.seq,
            "icode": residue.icode, "kind": residue.kind, "label": residue.label,
            "first": residue.first, "count": residue.count}[column]  # fmt: skip


def _object(saved: dict, key: str) -> dict:
    value = saved.get(key)
    if not isinstance(value, dict):
        raise SpecError(f"'{key}' is missing")
    return value


def _list(saved: dict, key: str) -> list:
    value = saved.get(key, [])
    if not isinstance(value, list):
        raise SpecError(f"'{key}' is not a list")
    return value


def _column(saved: dict, key: str, kind: type, length: int | None = None) -> list:
    value = saved.get(key)
    if not isinstance(value, list) or not all(_is(v, kind) for v in value):
        raise SpecError(f"'{key}' is not a list of {kind.__name__}")
    if length is not None and len(value) != length:
        raise SpecError(f"'{key}' has {len(value)} values, not {length}")
    return value


def _is(value: Any, kind: type) -> bool:
    if kind is int:
        return isinstance(value, int) and not isinstance(value, bool)
    if kind is float:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, kind)


def _tuple(value: Any, length: int, n: int, what: str) -> tuple:
    if (
        not isinstance(value, list)
        or len(value) != length
        or not all(_is(v, int) for v in value)
        or not all(0 <= v < n for v in value[:2])
    ):
        raise SpecError(f"A {what} is not {length} atom references")
    return tuple(value)


def _molecule(saved: Any, residues: tuple[Residue, ...]) -> Molecule:
    if not isinstance(saved, dict):
        raise SpecError("A molecule is a JSON object")
    residue, molblock = saved.get("residue"), saved.get("molblock")
    if not _is(residue, int) or not 0 <= residue < len(residues) or not isinstance(molblock, str):
        raise SpecError("A molecule needs its residue and its molblock")
    charges = _column(saved, "charges", float) if "charges" in saved else []
    if charges and len(charges) != residues[residue].count:
        raise SpecError(f"{residues[residue].label}: one charge per atom")
    model = saved.get("charge_model", "")
    return Molecule(residue, molblock, tuple(float(q) for q in charges),
                    model if isinstance(model, str) else "")  # fmt: skip


@dataclass(frozen=True)
class ForceFieldSpec:
    """The force fields and the solvent; the worker turns the names into its files."""

    protein: str = FF14SB
    small_molecule: str = SAGE
    water: str = TIP3P
    solvent: str = GBN2

    def __post_init__(self) -> None:
        if self.protein not in WATERS:
            raise SpecError(f"No protein force field '{self.protein}'")
        if self.water not in WATERS.values():
            raise SpecError(f"No water model '{self.water}'")
        if self.solvent not in (GBN2, EXPLICIT, VACUUM):
            raise SpecError(f"No solvent '{self.solvent}'")

    @classmethod
    def paired(cls, protein: str) -> ForceFieldSpec:
        """*protein* with the water model it was fitted with."""
        if protein not in WATERS:
            raise SpecError(f"No protein force field '{protein}'")
        return cls(protein=protein, water=WATERS[protein])

    def to_json(self) -> dict[str, str]:
        return {"protein": self.protein, "small_molecule": self.small_molecule,
                "water": self.water, "solvent": self.solvent}  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> ForceFieldSpec:
        if not isinstance(saved, dict) or not all(
            isinstance(saved.get(k), str) for k in ("protein", "small_molecule", "water", "solvent")
        ):
            raise SpecError("A force field names its protein, small-molecule, water and solvent")
        return cls(saved["protein"], saved["small_molecule"], saved["water"], saved["solvent"])


@dataclass(frozen=True)
class Stage:
    """One minimisation stage: what moves, and how hard the rest is held (kJ/mol/nm²)."""

    name: str
    moving: str  # one of MOVING
    restraint: float = 0.0  # on the heavy atoms held, by harmonic positional restraints
    ligand_restraint: float = 0.0  # on the small molecules' heavy atoms (ALL only)

    def __post_init__(self) -> None:
        if self.moving not in MOVING:
            raise SpecError(f"A stage moves {', '.join(MOVING)}, not '{self.moving}'")
        if self.restraint < 0 or self.ligand_restraint < 0:
            raise SpecError("A restraint is not negative")


#: the staged protocol of the plan (M3): hydrogens alone, then side chains and waters with the
#: backbone held, then everything held to where it started, the ligand more loosely
STAGES = (
    Stage("Hydrogens", HYDROGENS),
    Stage("Side chains and waters", SIDE_CHAINS, restraint=10 * KCAL_PER_A2),
    Stage("Complex", ALL, restraint=1 * KCAL_PER_A2, ligand_restraint=0.25 * KCAL_PER_A2),
)


@dataclass(frozen=True)
class MinimiseSpec:
    stages: tuple[Stage, ...] = field(default=STAGES)
    tolerance: float = 10.0  # kJ/mol/nm, the largest force left on any atom
    max_iterations: int = 5000  # per stage (0: until the tolerance is met)
    rmsd_limit: float = 0.3  # Å: heavy atoms further from the deposit than this are marked ⚠

    def __post_init__(self) -> None:
        if not self.stages:
            raise SpecError("A minimisation has at least one stage")
        if self.tolerance <= 0 or self.max_iterations < 0 or self.rmsd_limit <= 0:
            raise SpecError("Tolerance and RMSD limit are positive, iterations not negative")

    def to_json(self) -> dict[str, Any]:
        return {
            "stages": [{"name": s.name, "moving": s.moving, "restraint": s.restraint,
                        "ligand_restraint": s.ligand_restraint} for s in self.stages],
            "tolerance": self.tolerance, "max_iterations": self.max_iterations,
            "rmsd_limit": self.rmsd_limit,
        }  # fmt: skip

    @classmethod
    def from_json(cls, saved: Any) -> MinimiseSpec:
        if not isinstance(saved, dict):
            raise SpecError("A minimisation is a JSON object")
        stages = []
        for s in _list(saved, "stages"):
            if (not isinstance(s, dict) or not isinstance(s.get("name"), str)
                    or not isinstance(s.get("moving"), str)
                    or not _is(s.get("restraint", 0.0), float)
                    or not _is(s.get("ligand_restraint", 0.0), float)):  # fmt: skip
                raise SpecError("A stage has a name, what moves and its restraints")
            stages.append(Stage(s["name"], s["moving"], float(s.get("restraint", 0.0)),
                                float(s.get("ligand_restraint", 0.0))))  # fmt: skip
        tolerance, iterations, limit = (saved.get(k) for k in
                                        ("tolerance", "max_iterations", "rmsd_limit"))  # fmt: skip
        if not (_is(tolerance, float) and _is(iterations, int) and _is(limit, float)):
            raise SpecError("A minimisation has a tolerance, an iteration limit and an RMSD limit")
        return cls(tuple(stages), float(tolerance), int(iterations), float(limit))
