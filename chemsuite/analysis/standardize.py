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

"""Standardizing structures before docking (or any 3D work), one molecule at a time.

The steps and their order follow the usual preparation workflow. They are the order of the
ChEMBL structure pipeline (Bento et al., J. Cheminform. 2020, 12, 51) and of RDKit's
MolStandardize (a port of MolVS), which also underlie the ligand-preparation steps of docking
suites:

1. **Clean up** (``rdMolStandardize.Cleanup``):
   - remove explicit hydrogens;
   - disconnect metals bonded to N, O or F (salts drawn as bonds);
   - normalise functional groups to one drawing (nitro, N-oxides, sulfoxides, azides…);
   - reionise, so that the strongest acid carries the charge.
2. **Strip salts and solvents:** keep the largest organic fragment.
3. **Neutralise:** add or remove hydrogens so that charges cancel where they can.
   Quaternary nitrogens stay charged.
4. **Remove isotope labels** (optional).
5. **Canonical tautomer** (optional): RDKit's canonical tautomer. It is one consistent form for
   each molecule, so duplicates collapse; it is not necessarily the most stable or the
   bioactive form. A tautomer that would lose drawn stereo is not used.
6. **Protonate at pH** (optional): each ionisable site in its majority form at the pH, from
   MolGpKa's predicted pKa values (``chemsuite.ionisation``). Each fragment is protonated on
   its own, so counter-ions that are kept stay as they are. Carbon acids are ignored: a
   carbanion is never the form to dock.

Stereo is never invented: undefined stereocentres and double bonds are counted in the notes,
so they can be enumerated or fixed before docking. 3D coordinates are a separate step
(Generate Conformers).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.core.jobs.parallel import ordered_map
from chemsuite.core.mols import mol_from_blob, mol_from_text, mol_to_blob, smiles_of

#: A structure cell: its RDKit blob, or its text when not parsed yet (None: blank).
Cell = bytes | str | None
CHUNK_ROWS = 200  # per task in a worker process
INLINE_CHUNK_ROWS = 20  # in the job's own thread, so a cancel is noticed quickly

_FAILURES = (ValueError, RuntimeError, ArithmeticError, KeyError, IndexError)


@dataclass(frozen=True)
class Options:
    """Which steps to run. The defaults are the usual preparation for docking."""

    cleanup: bool = True
    strip_salts: bool = True
    neutralise: bool = True
    remove_isotopes: bool = False
    tautomer: bool = False
    protonate: bool = False
    ph: float = 7.4

    def steps(self) -> list[str]:
        """The chosen steps, in the order they run, as short names."""
        chosen = [
            (self.cleanup, "clean up"),
            (self.strip_salts, "strip salts"),
            (self.neutralise, "neutralise"),
            (self.remove_isotopes, "remove isotopes"),
            (self.tautomer, "canonical tautomer"),
            (self.protonate, f"protonate at pH {self.ph:g}"),
        ]
        return [name for on, name in chosen if on]


#: The usual preparation for docking (the defaults).
DOCKING = Options()


@dataclass
class Outcome:
    """One molecule's result: the standardized molecule (None if it failed) and what changed."""

    mol: Chem.Mol | None
    notes: list[str] = field(default_factory=list)


@dataclass
class StandardizeResults:
    """One entry per row, in order. ``b""`` / ``""`` where a row has no result."""

    blobs: list[bytes] = field(default_factory=list)
    texts: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    changed: int = 0
    failed: int = 0

    def __len__(self) -> int:
        return len(self.blobs)

    def extend(self, other: StandardizeResults) -> None:
        self.blobs.extend(other.blobs)
        self.texts.extend(other.texts)
        self.notes.extend(other.notes)
        self.changed += other.changed
        self.failed += other.failed


def _step(mol: Chem.Mol, note: str, notes: list[str], fn: Callable[[Chem.Mol], Chem.Mol]):
    """Apply *fn*; record *note* if the molecule's SMILES changed."""
    before = smiles_of(mol)
    out = fn(mol)
    if smiles_of(out) != before:
        notes.append(note)
    return out


def _canonical_tautomer(mol: Chem.Mol) -> Chem.Mol:
    """RDKit's canonical tautomer, keeping drawn stereocentres.

    By default RDKit clears the stereo of any sp3 atom that tautomerism touches (an
    alpha-carbon can enolise). That would turn a single enantiomer into an undefined one,
    which docking must not do. So stereo is kept, and when the canonical tautomer would
    still lose any drawn stereo (a double bond that moves), the drawn tautomer is kept.
    Over the 3,024 FDA-approved structures this keeps every drawn stereo element.
    """
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    try:
        out = rdMolStandardize.TautomerEnumerator(params).Canonicalize(mol)
    except (ValueError, RuntimeError):  # RDKit cannot kekulize some candidate: keep the drawn one
        return mol
    return out if _defined_stereo(out) >= _defined_stereo(mol) else mol


def _defined_stereo(mol: Chem.Mol) -> int:
    return sum(
        info.specified == Chem.StereoSpecified.Specified for info in Chem.FindPotentialStereo(mol)
    )


def undefined_stereo(mol: Chem.Mol) -> tuple[int, int]:
    """(stereocentres, double bonds) that could carry stereo but have none given.

    A double bond to a cumulated centre (an azide's N=N+=N-, an allene's C=C=C) is linear
    there and has no E/Z, though RDKit lists it as potential stereo.
    """
    centres = bonds = 0
    for info in Chem.FindPotentialStereo(mol):
        if info.specified != Chem.StereoSpecified.Unspecified:
            continue
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            centres += 1
        elif info.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(info.centeredOn)
            ends = (bond.GetBeginAtom(), bond.GetEndAtom())
            if not any(_double_bonds(atom) > 1 for atom in ends):
                bonds += 1
    return centres, bonds


def _double_bonds(atom: Chem.Atom) -> int:
    return sum(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())


def protonate(mol: Chem.Mol, ph: float) -> Chem.Mol:
    """*mol* with each ionisable site in its majority form at *ph* (MolGpKa), fragment by
    fragment. A fragment MolGpKa cannot prepare is kept as it is."""
    from chemsuite.ionisation.molgpka import ACID
    from chemsuite.ionisation.profile import Profile, profile

    out = None
    for fragment in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True):
        found = profile(fragment) if fragment.GetNumHeavyAtoms() > 1 else None
        if found is not None:
            # MolGpKa scores some C-H as acids (an alpha-CH between carbonyls); a carbanion is
            # never the form to dock, so only heteroatoms lose a proton
            sites = tuple(
                site for site in found.sites
                if site.kind != ACID or found.mol.GetAtomWithIdx(site.atom).GetAtomicNum() != 6
            )  # fmt: skip
            found = Profile(found.mol, sites)
        part = found.major_species(ph) if found is not None else fragment
        out = part if out is None else Chem.CombineMols(out, part)
    return out if out is not None else mol


def standardize(mol: Chem.Mol, options: Options = DOCKING) -> Outcome:
    """Run the chosen steps on *mol* (see the module notes). Never changes *mol* itself."""
    notes: list[str] = []
    work = Chem.Mol(mol)
    if options.cleanup:
        work = _step(work, "cleaned up", notes, rdMolStandardize.Cleanup)
    if options.strip_salts and len(Chem.GetMolFrags(work)) > 1:
        chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        work = _step(work, "salts stripped", notes, chooser.choose)
    if options.neutralise:
        uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        work = _step(work, "neutralised", notes, uncharger.uncharge)
    if options.remove_isotopes:
        work = _step(work, "isotopes removed", notes,
                     lambda m: Chem.RemoveHs(rdMolStandardize.IsotopeParent(m)))  # fmt: skip
    if options.tautomer:
        work = _step(work, "tautomer changed", notes, _canonical_tautomer)
    centres, bonds = undefined_stereo(work)  # before protonation adds N-H centres
    if options.protonate:
        work = _step(work, f"protonated at pH {options.ph:g}", notes,
                     lambda m: protonate(m, options.ph))  # fmt: skip
        charge = Chem.GetFormalCharge(work)
        if charge:
            notes.append(f"net charge {charge:+d}")
    if centres:
        notes.append(f"{centres} undefined stereocentre{'s' if centres > 1 else ''}")
    if bonds:
        notes.append(f"{bonds} undefined double bond{'s' if bonds > 1 else ''}")
    return Outcome(work, notes)


def _mol(cell: Cell) -> Chem.Mol | None:
    if isinstance(cell, bytes):
        return mol_from_blob(cell) if cell else None
    if isinstance(cell, str):
        return mol_from_text(cell)
    return None


def standardize_chunk(args: tuple[Options, Sequence[Cell]]) -> StandardizeResults:
    """Results for one chunk of cells (a picklable worker entry point)."""
    options, cells = args
    out = StandardizeResults()
    with rdBase.BlockLogs():
        for cell in cells:
            blob, text, note = b"", "", ""
            try:
                mol = _mol(cell)
                if mol is not None:
                    outcome = standardize(mol, options)
                    blob, text = mol_to_blob(outcome.mol), smiles_of(outcome.mol)
                    note = "; ".join(outcome.notes)
                    out.changed += text != smiles_of(mol)
                elif cell is not None:
                    note = "not a valid structure"
            except _FAILURES:
                blob, text, note = b"", "", "could not be standardized"
            if cell is not None and not blob:
                out.failed += 1
            out.blobs.append(blob)
            out.texts.append(text)
            out.notes.append(note)
    return out


def standardize_cells(
    cells: Sequence[Cell],
    options: Options = DOCKING,
    *,
    executor: Executor | None = None,
    progress: Callable[[int, int], None] | None = None,
    cancel: CancelToken | None = None,
) -> StandardizeResults:
    """Standardize each cell.

    If *cancel* is set, raises :class:`Cancelled` whose ``partial`` holds the results for the
    rows finished so far (the first rows, in order), or None if none were.
    """
    total = len(cells)
    size = CHUNK_ROWS if executor is not None else INLINE_CHUNK_ROWS
    chunks = ((options, list(cells[i : i + size])) for i in range(0, total, size))
    results = StandardizeResults()
    try:
        for part in ordered_map(standardize_chunk, chunks, executor, cancel=cancel):
            results.extend(part)
            if progress is not None:
                progress(len(results), total)
    except Cancelled:
        raise Cancelled(results if len(results) else None) from None
    return results
