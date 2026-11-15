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

"""Turning reagents into substituents: the leaving group goes, an attachment point stays.

A reagent as a catalogue lists it is not a substituent. An acyl chloride carries a chloride
that leaves when it acylates, so what joins a core is the acyl fragment. This module names
each kind of reagent and the bond it makes, and writes the fragment with one attachment point
(``*``), which :mod:`chemsuite.analysis.rgroup_library` then numbers and joins to a core:
``CC(=O)Cl`` becomes ``*C(C)=O``, and on ``C1CCN([*:1])CC1`` that gives the amide.

Each reagent kind is a one-reactant transform whose product holds a plain ``*``. The patterns
are tried in **tiers**, and the first tier that matches wins: the leaving group of an aryl
halide is iodide before bromide before triflate before chloride, so a bromo-chloroarene is
read as one site, not two.

**Every matching site is prepared.** A reagent with two of the same group gives one fragment
per site, and identical fragments collapse, so a symmetric diacid gives one. A fragment from
a bifunctional reagent still carries the other group, live: ``sites`` counts them so they can
be found and checked.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field
from functools import cache

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from chemsuite.analysis.reaction_templates import (
    ACID,
    ALKYL,
    ALKYL_LEAVING,
    AMINE,
    HARTENFELLER,
)
from chemsuite.analysis.reactions import TemplateError, apply_template, parse_template
from chemsuite.analysis.standardize import Options, standardize
from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
CHUNK_ROWS = 200  # per task in a worker process
INLINE_CHUNK_ROWS = 20  # in the job's own thread, so a cancel is noticed quickly
#: from this many rows, preparing is worth a process pool
POOL_MIN_ROWS = 1000
#: how many structures are sampled to tell whether a column holds this kind of reagent
FIT_SAMPLE = 100

#: salts and charges are cleaned up first, so a hydrochloride still matches its free base
CLEAN = Options(cleanup=True, strip_salts=True, neutralise=True)
#: the order groups are offered in
GROUPS = ("Acylation", "Arylation and Alkylation", "Nucleophile")

_FAILURES = (ValueError, RuntimeError, ArithmeticError, KeyError, IndexError)


@dataclass(frozen=True)
class Reagent:
    """A kind of reagent and the fragment it becomes.

    *smarts* holds one transform per tier, tried in order, each with exactly one reactant and
    a plain ``*`` in its product where the new bond is made. *attaches* says what leaves and
    what attaches, for the tooltip.
    """

    id: str
    name: str
    group: str
    smarts: tuple[str, ...]
    attaches: str
    source: str = ""
    note: str = ""

    def describe(self) -> str:
        """One line for a tooltip: what attaches, what it assumes, and where it comes from."""
        parts = [self.attaches]
        if self.note:
            parts.append(self.note)
        if self.source:
            parts.append(self.source)
        return ". ".join(parts)


REAGENTS: tuple[Reagent, ...] = (
    # -- acylation ------------------------------------------------------------------------------
    Reagent(
        id="acyl-chloride",
        name="Acyl chloride",
        group="Acylation",
        smarts=("[C:1](=[O:2])[Cl]>>[C:1](=[O:2])*",),
        attaches="The chloride leaves and the acyl carbon attaches, making an amide or ester",
        source=HARTENFELLER,
    ),
    Reagent(
        id="carboxylic-acid",
        name="Carboxylic acid",
        group="Acylation",
        smarts=(f"{ACID.format(n=1, o=2)}>>[C:1](=[O:2])*",),
        attaches="The hydroxyl leaves and the acyl carbon attaches, making an amide or ester",
        source=HARTENFELLER,
        note="The acid is treated as activated (HATU, EDC)",
    ),
    Reagent(
        id="sulfonyl-chloride",
        name="Sulfonyl chloride",
        group="Acylation",
        smarts=("[S:1](=[O:2])(=[O:3])[Cl]>>[S:1](=[O:2])(=[O:3])*",),
        attaches="The chloride leaves and the sulfonyl attaches, making a sulfonamide",
        source=HARTENFELLER,
    ),
    Reagent(
        id="chloroformate",
        name="Chloroformate",
        group="Acylation",
        smarts=("[C:1](=[O:2])([Cl])[O:3]>>[C:1](=[O:2])([O:3])*",),
        attaches="The chloride leaves and the carbonyl attaches, making a carbamate",
        source=HARTENFELLER,
    ),
    Reagent(
        id="isocyanate",
        name="Isocyanate",
        group="Acylation",
        smarts=("[N:1]=[C:2]=[O:3]>>[N;H1:1][C:2](=[O:3])*",),
        attaches="Nothing leaves: the carbonyl attaches and the nitrogen takes a hydrogen, "
        "making a urea",
        source=HARTENFELLER,
    ),
    # -- arylation and alkylation ---------------------------------------------------------------
    Reagent(
        id="aryl-halide",
        name="Aryl halide or triflate",
        group="Arylation and Alkylation",
        smarts=(
            "[c:1][I]>>[c:1]*",
            "[c:1][Br]>>[c:1]*",
            "[c:1][$(OS(=O)(=O)C(F)(F)F)]>>[c:1]*",
            "[c:1][Cl]>>[c:1]*",
        ),
        attaches="The halide or triflate leaves and the aromatic carbon attaches (Suzuki, "
        "Buchwald, SNAr)",
        source=HARTENFELLER,
        note="Iodide before bromide before triflate before chloride, so only one leaves",
    ),
    Reagent(
        id="boronic",
        name="Boronic acid or ester",
        group="Arylation and Alkylation",
        smarts=("[c:1][$(B(O)O),$(B1OC(C)(C)C(C)(C)O1)]>>[c:1]*",),
        attaches="The boron leaves and the aromatic carbon attaches (Suzuki, Chan-Lam)",
        source=HARTENFELLER,
    ),
    Reagent(
        id="alkyl-halide",
        name="Alkyl halide or sulfonate",
        group="Arylation and Alkylation",
        smarts=(f"{ALKYL.format(n=1)}{ALKYL_LEAVING}>>[C:1]*",),
        attaches="The halide or sulfonate leaves and the carbon attaches (N- or O-alkylation)",
        source=HARTENFELLER,
    ),
    Reagent(
        id="aldehyde",
        name="Aldehyde",
        group="Arylation and Alkylation",
        smarts=("[C;H1:1]=[O]>>[C:1]*",),
        attaches="The oxygen leaves and the carbon attaches, as in reductive amination",
        source=HARTENFELLER,
    ),
    Reagent(
        id="alkyne",
        name="Terminal alkyne",
        group="Arylation and Alkylation",
        smarts=("[CH1:1]#[C:2]>>*[C:1]#[C:2]",),
        attaches="The terminal hydrogen leaves and that carbon attaches (Sonogashira)",
        source=HARTENFELLER,
    ),
    # -- nucleophiles ---------------------------------------------------------------------------
    Reagent(
        id="amine",
        name="Amine (N-H)",
        group="Nucleophile",
        smarts=(f"{AMINE.format(n=1)}>>[N:1]*",),
        attaches="A hydrogen leaves the nitrogen and the nitrogen attaches",
        source=HARTENFELLER,
        note="Amide, sulfonamide, hydrazine, hydroxylamine and imine nitrogens are not amines",
    ),
    Reagent(
        id="alcohol",
        name="Alcohol or phenol",
        group="Nucleophile",
        smarts=("[#6;!$([#6]=[O,N,S]):1][OX2H1:2]>>[#6:1][O:2]*",),
        attaches="The hydroxyl hydrogen leaves and the oxygen attaches, making an ether or ester",
    ),
    Reagent(
        id="thiol",
        name="Thiol",
        group="Nucleophile",
        smarts=("[#6:1][SX2H1:2]>>[#6:1][S:2]*",),
        attaches="The thiol hydrogen leaves and the sulfur attaches, making a thioether",
    ),
)

_BY_ID = {reagent.id: reagent for reagent in REAGENTS}


def by_id(reagent_id: str) -> Reagent | None:
    """The reagent kind with this id, or None."""
    return _BY_ID.get(reagent_id)


def in_group(group: str) -> tuple[Reagent, ...]:
    """The reagent kinds of *group*, in their listed order."""
    return tuple(r for r in REAGENTS if r.group == group)


@cache
def _tiers(reagent_id: str) -> tuple[rdChemReactions.ChemicalReaction, ...]:
    """Each tier of this reagent kind as a reaction, in order. A tier that will not parse is
    left out; ``test_every_reagent_smarts_parses`` is what keeps that from happening."""
    reagent = _BY_ID[reagent_id]
    out = []
    for text in reagent.smarts:
        try:
            out.append(parse_template(text))
        except TemplateError:
            continue
    return tuple(out)


def sites(reagent_id: str, mol: Chem.Mol) -> tuple[int, int]:
    """``(tier, how many sites)`` for *mol*: the first tier that matches, and its match count.

    ``(-1, 0)`` means *mol* is not this kind of reagent. The count is of distinct matches of
    the tier's reactant pattern, so a symmetric diacid counts 2 even though both sites give
    the same fragment.
    """
    with rdBase.BlockLogs():
        for index, rxn in enumerate(_tiers(reagent_id)):
            try:
                found = mol.GetSubstructMatches(rxn.GetReactants()[0], uniquify=True)
            except _FAILURES:
                continue
            if found:
                return index, len(found)
    return -1, 0


def _cleaned(mol: Chem.Mol, clean: bool) -> Chem.Mol:
    """*mol* cleaned up, its salts stripped and its charges neutralised, or *mol* itself."""
    if not clean:
        return mol
    try:
        return standardize(mol, CLEAN).mol
    except _FAILURES:
        return mol


def any_fits(reagent_id: str, mols: Sequence[Chem.Mol | None], *, clean: bool = True) -> bool:
    """Whether any of *mols* is this kind of reagent, cleaned up the way :func:`prepare` would
    clean it, so a hydrochloride counts as its free base."""
    return any(mol is not None and sites(reagent_id, _cleaned(mol, clean))[1] for mol in mols)


@dataclass(frozen=True)
class Prepared:
    """What one reagent gave: a fragment per site (identical ones collapsed), and how many
    sites were found. No fragment and no site: not this kind of reagent."""

    fragments: tuple[str, ...] = ()
    sites: int = 0

    @property
    def fragment(self) -> str:
        """The first fragment, for a column that holds one value per row."""
        return self.fragments[0] if self.fragments else ""


@dataclass
class PreparedRows:
    """:class:`Prepared` per row, in the rows' order."""

    rows: list[Prepared] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[Prepared]:
        return iter(self.rows)

    def extend(self, other: PreparedRows) -> None:
        self.rows.extend(other.rows)

    @property
    def fragments(self) -> int:
        """How many fragments were prepared in all."""
        return sum(len(row.fragments) for row in self.rows)

    @property
    def missing(self) -> int:
        """How many rows are not this kind of reagent (a blank row is not counted)."""
        return sum(1 for row in self.rows if not row.sites and not row.fragments)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def prepare(reagent_id: str, mol: Chem.Mol, *, clean: bool = True) -> Prepared:
    """*mol* as substituent fragments, one per site. *mol* is never changed.

    With *clean*, the reagent is cleaned up, its salts stripped and its charges neutralised
    first, so a hydrochloride is prepared like its free base.
    """
    work = _cleaned(mol, clean)
    tier, count = sites(reagent_id, work)
    if not count:
        return Prepared()
    return Prepared(tuple(apply_template(_tiers(reagent_id)[tier], [work])), count)


def prepare_chunk(args: tuple[str, bool, Sequence[Cell]]) -> PreparedRows:
    """One chunk of cells prepared (a picklable worker entry point)."""
    reagent_id, clean, cells = args
    out = PreparedRows()
    with rdBase.BlockLogs():
        for cell in cells:
            mol = _mol(cell)
            if mol is None:
                out.rows.append(Prepared())
                continue
            try:
                out.rows.append(prepare(reagent_id, mol, clean=clean))
            except _FAILURES:
                out.rows.append(Prepared())
    return out


def prepare_cells(
    reagent_id: str,
    cells: Sequence[Cell],
    *,
    clean: bool = True,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> PreparedRows:
    """Prepare every cell.

    If *cancel* is set, raises :class:`~chemsuite.core.jobs.Cancelled` whose ``partial`` holds
    the rows finished so far (the first rows, in order), or None if none were.
    """
    if reagent_id not in _BY_ID:
        raise ValueError(f"unknown reagent kind: {reagent_id!r}")
    total = len(cells)
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    chunks = ((reagent_id, clean, list(cells[i : i + size])) for i in range(0, total, size))
    out = PreparedRows()
    try:
        for part in ordered_map(prepare_chunk, chunks, executor, cancel=cancel):
            out.extend(part)
            if progress is not None:
                progress(len(out), total)
    except Cancelled:
        raise Cancelled(out if len(out) else None) from None
    return out


__all__ = [
    "CLEAN",
    "FIT_SAMPLE",
    "GROUPS",
    "POOL_MIN_ROWS",
    "REAGENTS",
    "Prepared",
    "PreparedRows",
    "Reagent",
    "any_fits",
    "by_id",
    "in_group",
    "prepare",
    "prepare_cells",
    "prepare_chunk",
    "sites",
]
