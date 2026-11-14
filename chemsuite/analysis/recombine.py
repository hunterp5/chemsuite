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

"""Building new molecules from a pool of BRICS fragments.

RDKit's ``BRICS.BRICSBuild`` joins fragments only where the BRICS rules allow, using the
numbered attachment labels :mod:`chemsuite.sar.fragments` keeps (``[4*]``, ``[16*]``). A
fragment whose attachment points are plain ``*`` carries no rule and is refused: relabelling
it would invent bond types the rules never sanctioned, which is the mistake this replaces.

The build is a stream, not a fixed set of combinations: it yields candidates until the
product cap is met or it runs out. Each candidate is checked against :class:`Limits` and,
optionally, against molecules that already exist, so a run reports not only what it made but
why the rest went.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem import BRICS, Crippen, Descriptors, Lipinski, rdMolDescriptors

from chemsuite.core.jobs import Cancelled, CancelToken

Progress = Callable[[int, int], None]

DEFAULT_MAX_PRODUCTS = 2000
DEFAULT_MAX_DEPTH = 3
#: candidates examined per accepted product before a filtered run gives up
CANDIDATE_MULTIPLIER = 50
MIN_CANDIDATES = 5000
_LABELLED_DUMMY = re.compile(r"\[\d+\*\]")
_FAILURES = (ValueError, RuntimeError)


class RecombineError(ValueError):
    """The fragment pool cannot be recombined."""


#: property name → how to measure it, and how it is shown
PROPERTIES: dict[str, tuple[Callable[[Chem.Mol], float], str]] = {
    "MW": (Descriptors.MolWt, "molecular weight"),
    "cLogP": (Crippen.MolLogP, "cLogP"),
    "TPSA": (rdMolDescriptors.CalcTPSA, "TPSA"),
    "HBD": (lambda m: float(Lipinski.NumHDonors(m)), "H-bond donors"),
    "HBA": (lambda m: float(Lipinski.NumHAcceptors(m)), "H-bond acceptors"),
    "RotB": (lambda m: float(Lipinski.NumRotatableBonds(m)), "rotatable bonds"),
}


@dataclass(frozen=True)
class Limits:
    """Property ranges a product must fall in. A range left as None is not checked."""

    ranges: dict[str, tuple[float | None, float | None]] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return any(low is not None or high is not None for low, high in self.ranges.values())

    def passes(self, mol: Chem.Mol) -> bool:
        for name, (low, high) in self.ranges.items():
            measure = PROPERTIES.get(name)
            if measure is None:
                continue
            try:
                value = float(measure[0](mol))
            except _FAILURES:
                return False
            if (low is not None and value < low) or (high is not None and value > high):
                return False
        return True


#: ready-made limits, by the name a dialog shows
PRESETS: dict[str, Limits] = {
    "None": Limits(),
    # Lipinski et al., Adv. Drug Deliv. Rev. 1997, 23, 3
    "Rule of 5": Limits(
        {"MW": (None, 500.0), "cLogP": (None, 5.0), "HBD": (None, 5.0), "HBA": (None, 10.0)}
    ),
    # Teague et al., Angew. Chem. Int. Ed. 1999, 38, 3743
    "Lead-like": Limits(
        {
            "MW": (200.0, 350.0),
            "cLogP": (-1.0, 3.0),
            "HBD": (None, 3.0),
            "HBA": (None, 6.0),
            "RotB": (None, 7.0),
        }
    ),
}


@dataclass
class Recombined:
    """What a build made, and what became of the rest."""

    products: list[str] = field(default_factory=list)
    examined: int = 0  # candidates the build offered
    filtered: int = 0  # candidates a property limit removed
    known: int = 0  # candidates that already exist
    stopped: bool = False  # the candidate budget ran out before the product cap

    def __len__(self) -> int:
        return len(self.products)


def usable_fragments(fragments: Iterable[str]) -> list[str]:
    """The distinct fragments of *fragments* that carry numbered BRICS attachment points."""
    seen: dict[str, None] = {}
    for text in fragments:
        stripped = (text or "").strip()
        if stripped and _LABELLED_DUMMY.search(stripped):
            seen.setdefault(stripped, None)
    return list(seen)


def recombine(
    fragments: Sequence[str],
    *,
    max_depth: int = DEFAULT_MAX_DEPTH,
    max_products: int = DEFAULT_MAX_PRODUCTS,
    limits: Limits | None = None,
    known: Iterable[str] = (),
    progress: Progress | None = None,
    cancel: CancelToken | None = None,
) -> Recombined:
    """Join *fragments* where the BRICS rules allow, into at most *max_products* molecules.

    *max_depth* is the most fragments one product may be built from. *limits* keeps only
    products inside its property ranges; *known* is canonical SMILES that already exist,
    which are dropped and counted. Cancelled: ``partial`` holds what was made.
    """
    pool = usable_fragments(fragments)
    if len(pool) < 2:
        raise RecombineError(
            "Recombining needs at least two different fragments with numbered attachment "
            "points — run Fragments with BRICS and “Plain attachment points” off."
        )
    mols = []
    with rdBase.BlockLogs():
        for smiles in pool:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mols.append(mol)
    if len(mols) < 2:
        raise RecombineError("Fewer than two of the fragments could be read.")

    limits = limits or Limits()
    seen_known = {text for text in known if text}
    cap = max(int(max_products), 1)
    budget = max(MIN_CANDIDATES, cap * CANDIDATE_MULTIPLIER) if (limits or seen_known) else None
    out = Recombined()
    seen: set[str] = set()
    with rdBase.BlockLogs():
        builder = BRICS.BRICSBuild(
            mols,
            maxDepth=max(int(max_depth), 1),
            uniquify=True,
            scrambleReagents=False,
            onlyCompleteMols=True,
        )
        for candidate in builder:
            if cancel is not None and cancel.cancelled:
                raise Cancelled(out if out.products else None)
            out.examined += 1
            if budget is not None and out.examined > budget:
                out.stopped = True
                break
            try:
                Chem.SanitizeMol(candidate)
                smiles = Chem.MolToSmiles(candidate)
            except _FAILURES:
                continue
            if not smiles or smiles in seen:
                continue
            seen.add(smiles)
            if smiles in seen_known:
                out.known += 1
                continue
            if limits and not limits.passes(candidate):
                out.filtered += 1
                continue
            out.products.append(smiles)
            if progress is not None:
                progress(len(out.products), cap)
            if len(out.products) >= cap:
                break
    out.products.sort()
    return out


__all__ = [
    "CANDIDATE_MULTIPLIER",
    "DEFAULT_MAX_DEPTH",
    "DEFAULT_MAX_PRODUCTS",
    "PRESETS",
    "PROPERTIES",
    "Limits",
    "Recombined",
    "RecombineError",
    "recombine",
    "usable_fragments",
]
