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
"""Protein preparation P0: ground truth from neutron structures, inputs, and scoring.

A neutron (or joint X-ray/neutron) structure shows where the hydrogens and deuteriums are, so
it settles what preparation has to guess: His tautomers, which Asp/Glu/Lys/Cys/Tyr carry a
proton, which way round Asn/Gln/His side chains sit, and how many hydrogens each polar atom of
a ligand carries. See docs/prepare/p0-review.md.

- ``truth(path)`` reads the deposited structure. Exchangeable sites are often modelled as an
  H/D pair whose occupancies add up to the site's; a polar atom's hydrogen count is the summed
  occupancy of the H and D atoms bonded to it, rounded (a proton counts when it is there at
  least half the time).
- ``make_input(path, out, seed)`` writes what a method gets: heavy atoms only, the first
  conformation, waters as HOH, and half the Asn/Gln/His side chains flipped at random. The
  deposited orientation is the answer, so a method that never flips scores 50%, not 100%.
- ``observe(path, truth)`` reads a method's output with hydrogens the same way, by geometry
  rather than atom names, and ``score`` compares the two.
"""

from __future__ import annotations

import random
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

HYDROGEN = {"H", "D"}
POLAR = {"N", "O", "S"}
BOND = 1.25  # Å: a hydrogen this close to a polar atom is bonded to it
#: S–H and Se–H bonds are longer: 1.34 Å at nuclear positions (the CCD), ~1.2 Å as X-ray
#: refinement restrains them, so a 1.25 Å limit read every 1.34 Å thiol as a thiolate
LONG_BOND = {"S": 1.5, "Se": 1.6}


def bond_limit(atom: gemmi.Atom) -> float:
    return LONG_BOND.get(atom.element.name, BOND)


PRESENT = 0.5  # summed occupancy at which a proton counts
WATERS = {"HOH", "DOD", "WAT", "DIS"}
#: atoms swapped by a flip: the amide O and N, the His ring's ND1/CD2 and CE1/NE2
FLIPS = {
    "ASN": (("OD1", "ND2"),),
    "GLN": (("OE1", "NE2"),),
    "HIS": (("ND1", "CD2"), ("CE1", "NE2")),
}
#: the polar atoms whose hydrogen count decides each residue's state
TITRATABLE = {
    "HIS": ("ND1", "NE2"),
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "LYS": ("NZ",),
    "CYS": ("SG",),
    "TYR": ("OH",),
}
METALS = {"ZN", "FE", "CU", "MN", "CO", "NI", "CA", "MG", "CD", "HG"}

Key = tuple[str, int, str]  # chain, sequence number, insertion code


def key_of(chain: gemmi.Chain, residue: gemmi.Residue) -> Key:
    return (chain.name, residue.seqid.num, residue.seqid.icode.strip())


@dataclass
class Truth:
    """What a structure shows. ``states`` maps a residue to its state name ("HID", "HIE",
    "HIP", "ASP"/"ASH", "LYS"/"LYN", ...); ``flips`` holds the deposited positions of the
    atoms a flip swaps; ``ligand_h`` maps (residue, atom) to the hydrogens on each polar atom
    of a ligand."""

    pdb_id: str
    ph: float | None
    states: dict[tuple[Key, str], str] = field(default_factory=dict)
    flips: dict[tuple[Key, str], dict[str, np.ndarray]] = field(default_factory=dict)
    ligand_h: dict[tuple[Key, str, str], int] = field(default_factory=dict)
    #: the same atoms' hydrogens as the CCD draws the component
    ligand_ideal: dict[tuple[Key, str, str], int] = field(default_factory=dict)


def _first(residue: gemmi.Residue) -> dict[str, gemmi.Atom]:
    """Heavy atoms by name, first conformation only."""
    atoms: dict[str, gemmi.Atom] = {}
    for atom in residue:
        if atom.element.name not in HYDROGEN and atom.name not in atoms:
            atoms[atom.name] = atom
    return atoms


def _disordered(residue: gemmi.Residue, names: Iterable[str]) -> bool:
    wanted = set(names)
    return any(a.altloc != "\0" and a.name in wanted for a in residue)


def _h_count(atom: gemmi.Atom, hydrogens: list[gemmi.Atom]) -> int:
    occupancy = sum(h.occ for h in hydrogens if h.pos.dist(atom.pos) <= bond_limit(atom))
    return int(occupancy + 1 - PRESENT)  # rounds at PRESENT, not at 0.5 of a whole


def _modelled(atom: gemmi.Atom, hydrogens: list[gemmi.Atom]) -> int:
    """How many H/D sites the depositors placed on *atom*, whatever their occupancy (an H/D
    pair at one site counts once)."""
    sites: list[gemmi.Position] = []
    for h in hydrogens:
        if h.pos.dist(atom.pos) <= bond_limit(atom) and all(h.pos.dist(p) > 0.3 for p in sites):
            sites.append(h.pos)
    return len(sites)


#: hydrogens that must all be modelled before a missing proton is believed; otherwise an
#: unmodelled hydrogen (a spinning NH3 or OH) would read as absent
NEEDS_MODELLED = {"LYS": {"NZ": 3}, "TYR": {"OH": 1}, "CYS": {"SG": 1}}


def _hydrogens(residue: gemmi.Residue) -> list[gemmi.Atom]:
    return [a for a in residue if a.element.name in HYDROGEN]


def state_name(name: str, counts: Mapping[str, int]) -> str | None:
    """A residue's state from the hydrogens on its polar atoms; None when undecided."""
    if name == "HIS":
        d1, e2 = counts["ND1"] > 0, counts["NE2"] > 0
        return {(True, False): "HID", (False, True): "HIE", (True, True): "HIP"}.get((d1, e2))
    if name in ("ASP", "GLU"):
        protonated = any(counts[a] > 0 for a in TITRATABLE[name])
        return {"ASP": "ASH", "GLU": "GLH"}[name] if protonated else name
    if name == "LYS":
        return {3: "LYS", 2: "LYN"}.get(counts["NZ"])
    if name == "CYS":
        return {1: "CYS", 0: "CYM"}.get(counts["SG"])
    if name == "TYR":
        return {1: "TYR", 0: "TYM"}.get(counts["OH"])
    return None


def _bonded_away(structure: gemmi.Structure, atom: gemmi.Atom, residue_key: Key) -> bool:
    """A Cys SG in a disulfide or on a metal: its state is not a protonation question."""
    ns = gemmi.NeighborSearch(structure[0], structure.cell, 3.0).populate()
    for mark in ns.find_atoms(atom.pos, "\0", radius=2.6):
        cra = mark.to_cra(structure[0])
        if cra.atom.element.name in HYDROGEN or cra.atom is atom:
            continue
        if cra.atom.name == "SG" and key_of(cra.chain, cra.residue) != residue_key:
            return True
        if cra.atom.element.name.upper() in METALS:
            return True
    return False


def crystal_ph(block: gemmi.cif.Block) -> float | None:
    value = block.find_value("_exptl_crystal_grow.pH")
    try:
        return float(value) if value not in (None, "?", ".") else None
    except ValueError:
        return None


def truth(path: Path, ideal: Mapping[str, Mapping[str, int]] | None = None) -> Truth:
    """*ideal* maps a ligand's name to the hydrogens each of its atoms carries in the CCD.
    A ligand atom that should carry hydrogens but has none modelled is left out (the
    depositors did not place it, which is not the same as seeing it absent); a ligand the
    CCD does not cover is left out entirely."""
    ideal = ideal or {}
    block = gemmi.cif.read(str(path)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    found = Truth(path.stem.upper(), crystal_ph(block))
    model = structure[0]
    for chain in model:
        for residue in chain:
            name = residue.name
            if name in WATERS:
                continue
            key = key_of(chain, residue)
            heavy = _first(residue)
            hydrogens = _hydrogens(residue)
            polymer = residue.entity_type == gemmi.EntityType.Polymer
            if polymer and name in TITRATABLE:
                atoms = TITRATABLE[name]
                if all(a in heavy for a in atoms) and not _disordered(residue, atoms):
                    counts = {a: _h_count(heavy[a], hydrogens) for a in atoms}
                    needs = NEEDS_MODELLED.get(name, {})
                    skip = any(_modelled(heavy[a], hydrogens) < n for a, n in needs.items())
                    skip = skip or (name == "CYS" and _bonded_away(structure, heavy["SG"], key))
                    state = None if skip else state_name(name, counts)
                    if state is not None:
                        found.states[(key, name)] = state
            if polymer and name in FLIPS:
                atoms = [a for pair in FLIPS[name] for a in pair]
                if all(a in heavy for a in atoms) and not _disordered(residue, atoms):
                    # only where the hydrogens show which atom is the N
                    n_atom = "ND2" if name == "ASN" else "NE2" if name == "GLN" else None
                    shown = n_atom is None or _h_count(heavy[n_atom], hydrogens) == 2
                    if shown:
                        found.flips[(key, name)] = {
                            a: np.array(heavy[a].pos.tolist()) for a in atoms
                        }
            if not polymer and residue.het_flag == "H" and name.upper() not in METALS:
                expected = ideal.get(name)
                if not hydrogens or expected is None:  # nothing modelled, or no CCD entry
                    continue
                for atom_name, atom in heavy.items():
                    if atom.element.name not in POLAR or atom_name not in expected:
                        continue
                    if expected[atom_name] > 0 and not _modelled(atom, hydrogens):
                        continue
                    found.ligand_h[(key, name, atom_name)] = _h_count(atom, hydrogens)
                    found.ligand_ideal[(key, name, atom_name)] = expected[atom_name]
    return found


def make_input(path: Path, out: Path, seed: int) -> dict[tuple[Key, str], bool]:
    """Write the method's input (PDB format) and return which side chains were flipped."""
    structure = gemmi.read_structure(str(path))
    structure.remove_hydrogens()  # H and D
    structure.remove_alternative_conformations()
    structure.remove_empty_chains()
    rng = random.Random(f"{seed}:{path.stem.upper()}")  # one sequence per entry, not shared
    flipped: dict[tuple[Key, str], bool] = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.name in ("DOD", "WAT"):
                residue.name = "HOH"
            if residue.name in FLIPS:
                atoms = {a.name: a for a in residue}
                if not all(a in atoms for pair in FLIPS[residue.name] for a in pair):
                    continue
                flip = rng.random() < 0.5
                flipped[(key_of(chain, residue), residue.name)] = flip
                if flip:
                    for first, second in FLIPS[residue.name]:
                        a, b = atoms[first], atoms[second]
                        a.pos, b.pos = gemmi.Position(b.pos), gemmi.Position(a.pos)
    structure.setup_entities()
    out.parent.mkdir(parents=True, exist_ok=True)
    options = gemmi.PdbWriteOptions()
    options.minimal_file = True
    structure.write_pdb(str(out), options)
    return flipped


@dataclass
class Observed:
    """A method's answer, read the way ``truth`` reads a deposit. ``None`` values mean the
    residue or atom was not found in the method's output."""

    states: dict[tuple[Key, str], str | None] = field(default_factory=dict)
    flips: dict[tuple[Key, str], bool | None] = field(default_factory=dict)  # True: correct
    ligand_h: dict[tuple[Key, str, str], int | None] = field(default_factory=dict)


def _residues(structure: gemmi.Structure) -> dict[tuple[Key, str], gemmi.Residue]:
    found: dict[tuple[Key, str], gemmi.Residue] = {}
    for chain in structure[0]:
        for residue in chain:
            found.setdefault((key_of(chain, residue), residue.name), residue)
    return found


def observe(path: Path, expected: Truth) -> Observed:
    structure = gemmi.read_structure(str(path))
    structure.remove_alternative_conformations()
    residues = _residues(structure)
    seen = Observed()
    for (key, name), _state in expected.states.items():
        residue = residues.get((key, name)) or _renamed(residues, key, name)
        if residue is None:
            seen.states[(key, name)] = None
            continue
        heavy = _first(residue)
        hydrogens = _hydrogens(residue)
        atoms = TITRATABLE[name]
        if not all(a in heavy for a in atoms):
            seen.states[(key, name)] = None
            continue
        counts = {a: _modelled(heavy[a], hydrogens) for a in atoms}
        seen.states[(key, name)] = state_name(name, counts)
    for (key, name), positions in expected.flips.items():
        residue = residues.get((key, name)) or _renamed(residues, key, name)
        if residue is None:
            seen.flips[(key, name)] = None
            continue
        heavy = _first(residue)
        first, second = FLIPS[name][0]
        if first not in heavy:
            seen.flips[(key, name)] = None
            continue
        here = np.array(heavy[first].pos.tolist())
        right = np.linalg.norm(here - positions[first])
        wrong = np.linalg.norm(here - positions[second])
        seen.flips[(key, name)] = bool(right < wrong)
    for key, name, atom_name in expected.ligand_h:
        residue = residues.get((key, name))
        if residue is None:
            seen.ligand_h[(key, name, atom_name)] = None
            continue
        heavy = _first(residue)
        atom = heavy.get(atom_name)
        seen.ligand_h[(key, name, atom_name)] = (
            None if atom is None else _modelled(atom, _hydrogens(residue))
        )
    return seen


#: names methods give states (Amber-style HID/HIE/HIP, ASH, GLH, LYN, CYM)
_ALIASES = {
    "HIS": ("HID", "HIE", "HIP", "HSD", "HSE", "HSP"),
    "ASP": ("ASH",),
    "GLU": ("GLH",),
    "LYS": ("LYN",),
    "CYS": ("CYM", "CYX"),
}


def _renamed(residues, key: Key, name: str) -> gemmi.Residue | None:
    for alias in _ALIASES.get(name, ()):
        if (key, alias) in residues:
            return residues[(key, alias)]
    return None


@dataclass
class Tally:
    right: int = 0
    total: int = 0
    missing: int = 0

    def add(self, ok: bool | None) -> None:
        if ok is None:
            self.missing += 1
        else:
            self.total += 1
            self.right += bool(ok)

    def text(self) -> str:
        if not self.total:
            return "–"
        share = 100 * self.right / self.total
        return f"{share:.0f}% ({self.right}/{self.total})"


#: the usual state of each residue at pH 7; any other is a site preparation must find
STANDARD = {"ASP", "GLU", "LYS", "CYS", "TYR"}


def score(pairs: Iterable[tuple[Truth, Observed]]) -> dict[str, Tally]:
    """Accuracy per question over many structures."""
    tallies: dict[str, Tally] = {}

    def add(label: str, ok: bool | None) -> None:
        tallies.setdefault(label, Tally()).add(ok)

    for expected, seen in pairs:
        for site, state in expected.states.items():
            got = seen.states.get(site)
            name = site[1]
            ok = None if got is None else got == state
            if name == "HIS":
                # "HIS?": neutral, tautomer not given (PROPKA)
                add("His tautomer (HID/HIE/HIP)", None if got == "HIS?" else ok)
                add(
                    "His charge (HIP or neutral)",
                    None if got is None else ((got == "HIP") == (state == "HIP")),
                )
            else:
                add(f"{name} state", ok)
                if state not in STANDARD:  # the rare ones: does the method find them?
                    add(f"{name} uncommon state found", ok)
        for site, ok in seen.flips.items():
            add(f"{site[1]} orientation", ok)
        for site, count in expected.ligand_h.items():
            got = seen.ligand_h.get(site)
            ok = None if got is None else got == count
            add("Ligand polar-atom H count", ok)
            if count != expected.ligand_ideal.get(site, count):  # bound unlike the CCD's drawing
                add("Ligand atoms protonated unlike the CCD", ok)
    return tallies
