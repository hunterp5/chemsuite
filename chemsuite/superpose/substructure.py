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
"""Superpose on a common substructure: a core the user gives, or the maximum common
substructure (RDKit's FMCS) of the probe and the reference, fitted by Kabsch. Every
symmetric match of the core in the probe is tried, so a symmetric core fits the right way."""

from __future__ import annotations

from rdkit import Chem, rdBase
from rdkit.Chem import rdFMCS

from . import kabsch
from .model import CORE, CORE_ATOMS, CORE_RMSD, Fit, Settings

#: a core smaller than this gives no meaningful fit
MIN_ATOMS = 3


def parse_core(text: str) -> Chem.Mol | None:
    """A core typed as SMARTS, or as SMILES (None if it is neither, or empty)."""
    text = text.strip()
    if not text:
        return None
    with rdBase.BlockLogs():
        query = Chem.MolFromSmarts(text)
        if query is None:
            query = Chem.MolFromSmiles(text)
    return query if query is not None and query.GetNumAtoms() else None


def mcs_query(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> Chem.Mol | None:
    """The maximum common substructure of the two molecules' heavy atoms, as a query."""
    params = rdFMCS.MCSParameters()
    if not settings.mcs_compare_elements:
        params.AtomTyper = rdFMCS.AtomCompare.CompareAny
    params.BondTyper = rdFMCS.BondCompare.CompareOrder
    params.BondCompareParameters.RingMatchesRingOnly = settings.mcs_ring_matches_ring
    params.BondCompareParameters.CompleteRingsOnly = settings.mcs_complete_rings
    params.AtomCompareParameters.RingMatchesRingOnly = settings.mcs_ring_matches_ring
    params.AtomCompareParameters.CompleteRingsOnly = settings.mcs_complete_rings
    params.Timeout = settings.mcs_timeout
    with rdBase.BlockLogs():
        found = rdFMCS.FindMCS([Chem.RemoveHs(reference), Chem.RemoveHs(probe)], params)
        if found.numAtoms < MIN_ATOMS:
            return None
        return Chem.MolFromSmarts(found.smartsString)


def query_for(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> Chem.Mol | None:
    if settings.method == CORE:
        return parse_core(settings.core)
    return mcs_query(probe, reference, settings)


def align(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit]:
    """A fit per probe conformer on the core (none if the core is not in both)."""
    query = query_for(probe, reference, settings)
    if query is None:
        return []
    reference_match = reference.GetSubstructMatch(query)
    probe_matches = probe.GetSubstructMatches(query, uniquify=False,
                                              maxMatches=kabsch.MAX_MATCHES)  # fmt: skip
    if len(reference_match) < MIN_ATOMS or not probe_matches:
        return []
    target = reference.GetConformer().GetPositions()
    fits = []
    for conf in probe.GetConformers():
        transform, rmsd = kabsch.best_fit(conf.GetPositions(), target, probe_matches,
                                          reference_match)  # fmt: skip
        scores = {CORE_RMSD: rmsd, CORE_ATOMS: float(len(reference_match))}
        fits.append(Fit(conf.GetId(), transform, scores))
    return fits
