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
"""Superpose on a substructure the user gives (the anchor), scored by another method.

Each probe conformer is fitted onto the reference by the anchor's atoms (Kabsch, every
symmetric match tried), and the method then scores it in that pose: shape, shape + features
and field score the pose as it is, and choose the conformer whose other atoms fit best. O3A
cannot score a pose it did not make, so it pairs atoms with the anchor's pairs held
(Open3DALIGN's constraints). **Core RMSD** is the anchor atoms' RMSD after the fit."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign, rdMolDescriptors, rdShapeAlign
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties

from . import field, kabsch, shape
from .model import (
    COLOR,
    COLOR_TANIMOTO,
    COMBO,
    CORE_RMSD,
    FIELD,
    FIELD_SIMILARITY,
    O3A,
    O3A_SCORE,
    RMSD,
    SHAPE_TANIMOTO,
    Fit,
    Settings,
    apply_transform,
)
from .substructure import MIN_ATOMS, parse_core

#: how hard O3A holds the anchor's atom pairs
O3A_CONSTRAINT_WEIGHT = 100.0


def _pose(
    xyz: np.ndarray, target: np.ndarray, matches: Sequence[Sequence[int]], ref_match: Sequence[int]
) -> tuple[np.ndarray, Sequence[int], float]:
    """The best Kabsch fit over the symmetric *matches*: (transform, match, RMSD)."""
    found: tuple[np.ndarray, Sequence[int], float] | None = None
    goal = target[list(ref_match)]
    for match in matches[: kabsch.MAX_MATCHES]:
        transform, rmsd = kabsch.fit(xyz[list(match)], goal)
        if found is None or rmsd < found[2]:
            found = (transform, match, rmsd)
    if found is None:
        raise ValueError("no atom pairs to fit")
    return found


def _shape_scores(work: Chem.Mol, reference: Chem.Mol, conf_id: int, colors: bool) -> dict:
    options = rdShapeAlign.ShapeInputOptions()
    options.useColors = colors
    shape, color = rdShapeAlign.ScoreMol(reference, work, options, options, 0, conf_id)
    if colors:
        return {COMBO: shape + color, SHAPE_TANIMOTO: shape, COLOR_TANIMOTO: color}
    return {SHAPE_TANIMOTO: shape}


def _o3a(work: Chem.Mol, reference: Chem.Mol, conf_id: int, pairs: list[tuple[int, int]],
         settings: Settings) -> tuple[float, float]:  # fmt: skip
    """O3A of one conformer with the anchor pairs held; moves it. (score, RMSD)"""
    weights = [O3A_CONSTRAINT_WEIGHT] * len(pairs)
    probe_props = reference_props = None
    if settings.o3a_typing == "mmff":
        probe_props = MMFFGetMoleculeProperties(work)
        reference_props = MMFFGetMoleculeProperties(reference)
    if probe_props is not None and reference_props is not None:
        aligner = rdMolAlign.GetO3A(work, reference, probe_props, reference_props,
                                    prbCid=conf_id, refCid=0, constraintMap=pairs,
                                    constraintWeights=weights)  # fmt: skip
    else:  # Crippen, or MMFF lacks parameters for one of them
        aligner = rdMolAlign.GetCrippenO3A(
            work, reference, rdMolDescriptors._CalcCrippenContribs(work),
            rdMolDescriptors._CalcCrippenContribs(reference), prbCid=conf_id, refCid=0,
            constraintMap=pairs, constraintWeights=weights,
        )  # fmt: skip
    rmsd = aligner.Align()
    return aligner.Score(), rmsd


def align(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit] | None:
    """A fit per probe conformer, posed on the anchor; None if the anchor is not in the probe
    (the caller fits it without), [] if it is not in the reference."""
    query = parse_core(settings.anchor)
    if query is None:
        return None
    probe_matches = probe.GetSubstructMatches(query, uniquify=False,
                                              maxMatches=kabsch.MAX_MATCHES)  # fmt: skip
    if not probe_matches:
        return None
    ref_match = reference.GetSubstructMatch(query)
    if len(ref_match) < MIN_ATOMS:
        return []
    target = reference.GetConformer().GetPositions()
    reference_model = (
        field.Model(reference, 0, settings.charges) if settings.method == FIELD else None
    )
    work = Chem.Mol(probe)  # each conformer is moved on its own
    fits = []
    with rdBase.BlockLogs():
        for conf in probe.GetConformers():
            cid = conf.GetId()
            before = conf.GetPositions()
            transform, match, core_rmsd = _pose(before, target, probe_matches, ref_match)
            moved = work.GetConformer(cid)
            apply_transform(moved, transform)
            rank = None
            if settings.method == O3A:
                pairs = [(int(p), int(r)) for p, r in zip(match, ref_match, strict=True)]
                score, rmsd = _o3a(work, reference, cid, pairs, settings)
                transform, _ = kabsch.fit(before, moved.GetPositions())
                after = moved.GetPositions()[list(match)]
                core_rmsd = kabsch.rmsd(after, target[list(ref_match)])
                scores = {O3A_SCORE: score, RMSD: rmsd}
            elif settings.method == FIELD:
                model = field.Model(work, cid, settings.charges)
                scorer = field.Scorer(model, reference_model, settings.field_weight)
                rank = scorer.objective(*scorer.scores(model.xyz))
                _tanimoto, carbo = scorer.scores(model.xyz)
                scores = {FIELD_SIMILARITY: carbo,
                          SHAPE_TANIMOTO: shape.score(work, reference, cid)}  # fmt: skip
            else:
                scores = _shape_scores(work, reference, cid, settings.method == COLOR)
            scores[CORE_RMSD] = core_rmsd
            fits.append(Fit(cid, transform, scores, rank=rank))
    return fits
