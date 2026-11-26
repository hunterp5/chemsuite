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

"""A prepared structure as a :class:`~.spec.SystemSpec`, or the reasons it cannot be one.

Each amino acid's form is read from the hydrogens it has, not from its name, since Prepare never
renames a residue: His with HD1 and HE2 is HIP, Asp with a proton on an oxygen ASH, Cys bonded
to another Cys CYX, and so on; a chain's first and last residues take the N- and C-terminal
forms unless a cap or a neighbour is bonded on. Where Prepare recorded the state it chose
(``Prepared.states``), the two must agree, and a disagreement is an error in chemsuite
(:class:`StateMismatch`), never something to work around.

What the force fields here cannot describe is a :class:`Problem` that blocks the build (the
m0 review's template survey, ``docs/simulate/m0-review.md``): tyrosinate (TYM), a neutral
N-terminus or C-terminus, imidazolate His, neutral Arg, ASH/GLH/LYN/CYM at a chain end, a
phospho-Ser, -Thr or -Tyr with a neutral phosphate or at a chain end (Amber's phosaa templates
take the dianion and the monoanion, mid-chain), any other residue that is not one of the 20
amino acids or a cap, an amino acid bonded to another by anything but a peptide bond or a
disulfide (a crosslink), nucleic acids and glycans, an ion
without a charge, and a ligand bonded to the protein or whose hydrogens do not make a molecule.
None of these is ever replaced by something close to it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from chemsuite.core.jobs import CancelToken, check
from chemsuite.macromolecule.chem import ligands
from chemsuite.macromolecule.structure.components import CcdCache
from chemsuite.macromolecule.structure.model import (
    ADDITIVE,
    AMINO,
    GLYCAN,
    ION,
    LIGAND,
    NUCLEIC,
    WATER,
    Structure,
)

from . import spec
from .spec import Molecule, Residue, SystemSpec

STANDARD = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)
#: phospho-Ser, -Thr and -Tyr: the Amber form (phosaa14SB, phosaa19SB) by how many protons
#: the phosphate oxygens carry, the dianion then the monoanion; a neutral one has no template
PHOSPHO = {"SEP": ("SEP", "S1P"), "TPO": ("TPO", "T1P"), "PTR": ("PTR", "Y1P")}
PHOSPHATE_O = ("O1P", "O2P", "O3P")
#: caps, by the name Amber gives them (the PDB's NH2 is Amber's NHE)
CAPS = {"ACE": "ACE", "NME": "NME", "NH2": "NHE"}
#: forms with no N- or C-terminal template (gate 4 of the m0 review)
MID_CHAIN_ONLY = frozenset(["ASH", "GLH", "LYN", "CYM", "SEP", "S1P", "TPO", "T1P", "PTR", "Y1P"])
#: the forms Prepare records (``prepare.solve.RECORDED``), as this module names them
NEUTRAL_N = "N-neutral"  # ``prepare.sites.NEUTRAL_N``
PROTEIN_STATES = frozenset(["HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYM", "TYM", NEUTRAL_N])


@dataclass(frozen=True)
class Problem:
    """Something in the structure that blocks the build, on the residue *label*."""

    label: str
    text: str

    def __str__(self) -> str:
        return f"{self.label}: {self.text}"


class BuildError(ValueError):
    """The structure cannot be built; :attr:`problems` says why, residue by residue."""

    def __init__(self, problems: list[Problem]) -> None:
        self.problems = tuple(problems)
        shown = "; ".join(str(p) for p in problems[:5])
        more = f" (and {len(problems) - 5} more)" if len(problems) > 5 else ""
        super().__init__(shown + more)


class StateMismatch(BuildError):
    """The hydrogens disagree with the state Prepare recorded: a bug, not the chemist's to fix."""


def system_spec(
    s: Structure,
    states: Mapping[int, str] | None = None,
    smiles: Mapping[int, str] | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> SystemSpec:
    """*s* (prepared: one location per atom, hydrogens placed) as a spec. *states* are
    Prepare's recorded states by residue, cross-checked against the hydrogens; *smiles* a
    ligand's SMILES by residue, for its bond orders. Raises :class:`BuildError` with every
    problem found, or :class:`StateMismatch`."""
    built = _Builder(s, states or {}, smiles or {}, ccd, cancel)
    built.run()
    if built.mismatches:
        raise StateMismatch(built.mismatches)
    if built.problems:
        raise BuildError(built.problems)
    return built.spec()


def problems(
    s: Structure,
    states: Mapping[int, str] | None = None,
    smiles: Mapping[int, str] | None = None,
    ccd: CcdCache | None = None,
    cancel: CancelToken | None = None,
) -> list[Problem]:
    """What blocks :func:`system_spec` on *s* (empty when it can be built)."""
    built = _Builder(s, states or {}, smiles or {}, ccd, cancel)
    built.run()
    return [*built.mismatches, *built.problems]


def positions(s: Structure, system: SystemSpec, model: int = 0) -> np.ndarray:
    """The coordinates (Å, (n, 3)) of *system*'s atoms in *s*'s *model*, in the spec's order.
    Raises :class:`BuildError` when an atom of the spec is no longer in *s*."""
    found = s.indices(np.array(system.uid, dtype=np.int64))
    if (found < 0).any():
        raise BuildError(
            [
                Problem(
                    s.meta.get("id", "structure"),
                    "Atoms of the system are no longer in the structure",
                )
            ]
        )
    return s.xyz(model)[found].astype(np.float64)


class _Builder:
    def __init__(
        self,
        s: Structure,
        states: Mapping[int, str],
        smiles: Mapping[int, str],
        ccd: CcdCache | None,
        cancel: CancelToken | None,
    ) -> None:
        self.s = s
        self.states = states
        self.smiles = smiles
        self.ccd = ccd
        self.cancel = cancel
        self.active = s.active()
        self.problems: list[Problem] = []
        self.mismatches: list[Problem] = []
        self.order: list[int] = []  # structure atom indices, in the spec's order
        self.residues: list[Residue] = []
        self.molecules: list[Molecule] = []
        self.neighbours: list[list[int]] = [[] for _ in range(s.n_atoms)]
        for a, b in s.bonds.tolist():
            if self.active[a] and self.active[b]:
                self.neighbours[a].append(b)
                self.neighbours[b].append(a)

    # -- the residues ------------------------------------------------------------------------
    def run(self) -> None:
        s = self.s
        if (s.altloc[self.active] != "").any():
            self.problems.append(
                Problem(
                    s.meta.get("id", "structure"),
                    "Alternate locations: prepare the structure first",
                )
            )
        for r in range(s.n_residues):
            check(self.cancel)
            atoms = [i for i in range(*s.atoms_of(r).indices(s.n_atoms)) if self.active[i]]
            if not atoms:
                continue
            kind = int(s.res_kind[r])
            label = s.residue_label(r)
            if kind == AMINO:
                self.amino(r, atoms, label)
            elif kind == WATER:
                self.water(r, atoms, label)
            elif kind == ION:
                self.ion(r, atoms, label)
            elif kind in (LIGAND, ADDITIVE):
                self.molecule(r, label)
            elif kind == NUCLEIC:
                self.problems.append(Problem(label, "Nucleic acids cannot be simulated yet"))
            elif kind == GLYCAN:
                self.problems.append(Problem(label, "Glycans cannot be simulated yet"))
            else:
                self.problems.append(Problem(label, "Of no kind a force field here describes"))

    def add(self, r: int, atoms: list[int], name: str, kind: str) -> None:
        s = self.s
        self.residues.append(Residue(name, str(s.chain_name[s.chain_of[r]]), int(s.res_seq[r]),
                                     str(s.icode[r]), kind, s.residue_label(r), len(self.order),
                                     len(atoms)))  # fmt: skip
        self.order += atoms

    def hydrogens(self, i: int) -> int:
        return sum(1 for j in self.neighbours[i] if self.s.element[j] == 1)

    def named(self, atoms: list[int]) -> dict[str, int]:
        return {str(self.s.name[i]): i for i in atoms}

    def amino(self, r: int, atoms: list[int], label: str) -> None:
        s = self.s
        name = str(s.res_name[r])
        if name in CAPS:
            self.add(r, atoms, CAPS[name], spec.CAP)
            return
        if name not in STANDARD and name not in PHOSPHO:
            self.problems.append(
                Problem(
                    label,
                    f"{name} is not one of the 20 amino acids: no "
                    "force-field template (rename it, or remove it)",
                )
            )
            return
        crosslinks = self.crosslinks(r, atoms)
        if crosslinks:
            self.problems.append(Problem(label, f"Bonded to {', '.join(crosslinks)}: no "
                                         "force-field template takes this link (delete one of "
                                         "its two atoms, then prepare again)"))  # fmt: skip
            return
        at = self.named(atoms)
        form = self.form(name, at, label)
        if form is None:
            return
        recorded = self.states.get(r)
        if recorded in PROTEIN_STATES and recorded != NEUTRAL_N and recorded != form:
            self.mismatches.append(Problem(label, f"recorded as {recorded}, but its hydrogens "
                                           f"make it {form}"))  # fmt: skip
        n_term, c_term = self.termini(r, at)
        if recorded == NEUTRAL_N and not (n_term and "N" in at and self._neutral_n(name, at)):
            self.mismatches.append(Problem(label, "recorded as a neutral N-terminus, but its "
                                           "hydrogens say otherwise"))  # fmt: skip
        if (n_term or c_term) and form in MID_CHAIN_ONLY:
            self.problems.append(Problem(label, f"{form} at a chain end has no force-field "
                                         "template"))  # fmt: skip
            return
        if n_term and "N" in at and self._neutral_n(name, at):
            self.problems.append(Problem(label, "A neutral N-terminus has no force-field "
                                         "template"))  # fmt: skip
            return
        if c_term:
            if "OXT" not in at:
                self.problems.append(Problem(label, "The chain ends without an OXT or a cap"))
                return
            if self.hydrogens(at["OXT"]) or ("O" in at and self.hydrogens(at["O"])):
                self.problems.append(Problem(label, "A neutral C-terminus has no force-field "
                                             "template"))  # fmt: skip
                return
        prefix = "N" if n_term else "C" if c_term else ""
        self.add(r, atoms, prefix + form, spec.PROTEIN)

    def form(self, name: str, at: dict[str, int], label: str) -> str | None:
        """The Amber form of a standard or phosphorylated amino acid from its hydrogens (None:
        no template)."""
        h = {n: self.hydrogens(i) for n, i in at.items()}
        if name == "HIS":
            d, e = h.get("ND1", 0), h.get("NE2", 0)
            if not d and not e:
                self.problems.append(
                    Problem(
                        label,
                        "His with no proton on its ring (imidazolate) has no force-field template",
                    )
                )
                return None
            return "HIP" if d and e else "HID" if d else "HIE"
        if name in PHOSPHO:
            forms = PHOSPHO[name]
            protons = sum(h.get(o, 0) for o in PHOSPHATE_O)
            if protons >= len(forms):
                self.problems.append(Problem(label, "A neutral phosphate has no force-field "
                                             "template: prepare the original structure "
                                             "again"))  # fmt: skip
                return None
            return forms[protons]
        if name == "ASP":
            return "ASH" if h.get("OD1", 0) + h.get("OD2", 0) else "ASP"
        if name == "GLU":
            return "GLH" if h.get("OE1", 0) + h.get("OE2", 0) else "GLU"
        if name == "LYS":
            return "LYN" if h.get("NZ", 3) == 2 else "LYS"
        if name == "TYR":
            if "OH" in at and not h.get("OH", 0):
                self.problems.append(Problem(label, "Tyrosinate (TYM) has no force-field "
                                             "template"))  # fmt: skip
                return None
            return "TYR"
        if name == "ARG":
            if h.get("NE", 0) + h.get("NH1", 0) + h.get("NH2", 0) < 5:
                self.problems.append(Problem(label, "Neutral Arg has no force-field template"))
                return None
            return "ARG"
        if name == "CYS" and "SG" in at:
            sg = at["SG"]
            partner = [j for j in self.neighbours[sg] if self.s.name[j] == "SG"
                       and self.s.residue_of[j] != self.s.residue_of[sg]]  # fmt: skip
            if partner:
                return "CYX"
            return "CYS" if h.get("SG", 0) else "CYM"
        return name

    def crosslinks(self, r: int, atoms: list[int]) -> list[str]:
        """The atoms of other amino acids that residue *r* is bonded to other than by a peptide
        bond or a disulfide (a ligand's or a glycan's link is reported on it)."""
        s = self.s
        found = []
        for i in atoms:
            for j in self.neighbours[i]:
                if s.residue_of[j] == r or s.res_kind[s.residue_of[j]] != AMINO:
                    continue
                pair = {str(s.name[i]), str(s.name[j])}
                if pair != {"C", "N"} and pair != {"SG"}:
                    found.append(f"{s.atom_label(j)} (by {s.name[i]})")
        return found

    def termini(self, r: int, at: dict[str, int]) -> tuple[bool, bool]:
        """Whether the residue starts a chain (its N bonded to no other residue) and ends one
        (its C bonded to no other residue)."""
        s = self.s

        def bonded_on(name: str) -> bool:
            i = at.get(name)
            return i is not None and any(s.residue_of[j] != r and s.atom_kind[j] != ION
                                         for j in self.neighbours[i])  # fmt: skip

        return not bonded_on("N"), not bonded_on("C")

    def _neutral_n(self, name: str, at: dict[str, int]) -> bool:
        """An N-terminal N with one proton too few (two for most, one for Pro)."""
        return self.hydrogens(at["N"]) < (2 if name == "PRO" else 3)

    def water(self, r: int, atoms: list[int], label: str) -> None:
        elements = sorted(int(self.s.element[i]) for i in atoms)
        if elements != [1, 1, 8]:
            self.problems.append(Problem(label, "A water without its two hydrogens"))
            return
        self.add(r, atoms, "HOH", spec.WATER)

    def ion(self, r: int, atoms: list[int], label: str) -> None:
        if any(int(self.s.charge[i]) == 0 for i in atoms):
            self.problems.append(Problem(label, "An ion of unknown charge: set its formal charge"))
            return
        self.add(r, atoms, str(self.s.res_name[r]), spec.ION)

    def molecule(self, r: int, label: str) -> None:
        try:
            found = ligands.prepared_form(self.s, r, self.smiles.get(r), ccd=self.ccd,
                                          cancel=self.cancel)  # fmt: skip
        except ligands.FormError as error:
            self.problems.append(Problem(label, str(error)))
            return
        mol = found.mol
        if found.guessed:
            self.problems.append(Problem(label, "Its bond orders are a guess from the geometry: "
                                         "give its SMILES first"))  # fmt: skip
            return
        atoms = [self.s.index_of_uid[a.GetIntProp("uid")] for a in mol.GetAtoms()]
        self.molecules.append(Molecule(len(self.residues), Chem.MolToMolBlock(mol)))
        self.add(r, atoms, str(self.s.res_name[r]), spec.MOLECULE)

    # -- the spec ------------------------------------------------------------------------------
    def spec(self) -> SystemSpec:
        s = self.s
        where = {atom: k for k, atom in enumerate(self.order)}
        bonds, disulfides = [], []
        for (a, b), order in zip(s.bonds.tolist(), s.order.tolist(), strict=True):
            if a not in where or b not in where or ION in (s.atom_kind[a], s.atom_kind[b]):
                continue  # a metal's contacts are not bonds to a force field here
            i, j = sorted((where[a], where[b]))
            bonds.append((i, j, int(order)))
            if s.name[a] == "SG" and s.name[b] == "SG" and s.residue_of[a] != s.residue_of[b]:
                disulfides.append((i, j))
        order = self.order
        return SystemSpec(
            uid=tuple(int(u) for u in s.uid[order]),
            element=tuple(int(e) for e in s.element[order]),
            name=tuple(str(n) for n in s.name[order]),
            formal=tuple(int(q) for q in s.charge[order]),
            residues=tuple(self.residues),
            bonds=tuple(sorted(bonds)),
            disulfides=tuple(sorted(disulfides)),
            molecules=tuple(self.molecules),
        )
