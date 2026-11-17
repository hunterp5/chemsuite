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
"""Atom-based superposition without a common core: Open3DALIGN (Tosco, Balle and Shiri,
J. Comput.-Aided Mol. Des. 2011, 25, 777), RDKit's ``GetO3AForProbeConfs``, which pairs atoms
by MMFF94 types and charges, or by Crippen logP contributions."""

from __future__ import annotations

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdMolDescriptors
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties

from . import kabsch
from .model import O3A_SCORE, RMSD, Fit, Settings


def align(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit]:
    """A fit per probe conformer (both molecules with explicit hydrogens)."""
    work = Chem.Mol(probe)  # O3A moves the conformers it aligns
    with rdBase.BlockLogs():
        mmff = settings.o3a_typing == "mmff"
        probe_props = MMFFGetMoleculeProperties(work) if mmff else None
        reference_props = MMFFGetMoleculeProperties(reference) if mmff else None
        if mmff and probe_props is not None and reference_props is not None:
            aligners = rdMolAlign.GetO3AForProbeConfs(work, reference, 1, probe_props,
                                                      reference_props, refCid=0)  # fmt: skip
        else:  # Crippen, or MMFF lacks parameters for one of them
            probe_contribs = rdMolDescriptors._CalcCrippenContribs(work)
            reference_contribs = rdMolDescriptors._CalcCrippenContribs(reference)
            aligners = rdMolAlign.GetCrippenO3AForProbeConfs(
                work, reference, 1, probe_contribs, reference_contribs, refCid=0
            )
        fits = []
        for conf, aligner in zip(work.GetConformers(), aligners, strict=True):
            before = probe.GetConformer(conf.GetId()).GetPositions()
            rmsd = aligner.Align()
            score = aligner.Score()
            transform, _ = kabsch.fit(before, conf.GetPositions())
            fits.append(Fit(conf.GetId(), transform, {O3A_SCORE: score, RMSD: rmsd}))
    return fits
