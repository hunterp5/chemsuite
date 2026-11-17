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
"""Shape superposition: atom-centred Gaussians, the overlap volume maximised (Grant and
Pickup, J. Phys. Chem. 1995, 99, 3503; the method of ROCS), with RDKit's ``rdShapeAlign``
(PubChem's shape code). With features, pharmacophore points (donors, acceptors, rings,
charges) overlap too ("color"), and Combo is shape plus color Tanimoto (0–2), as ROCS's
TanimotoCombo. Hydrogens are left out on both sides."""

from __future__ import annotations

from rdkit import Chem, rdBase
from rdkit.Chem import rdShapeAlign

from . import kabsch
from .model import COLOR, COLOR_TANIMOTO, COMBO, SHAPE_TANIMOTO, Fit, Settings


def align(probe: Chem.Mol, reference: Chem.Mol, settings: Settings) -> list[Fit]:
    """A fit per probe conformer."""
    colors = settings.method == COLOR
    balance = 1.0 - settings.color_weight if colors else 1.0  # rdShapeAlign: 1 is shape only
    fits = []
    work = Chem.Mol(probe)  # each conformer is moved on its own
    with rdBase.BlockLogs():
        for conf in probe.GetConformers():
            moved = work.GetConformer(conf.GetId())
            shape, color = rdShapeAlign.AlignMol(reference, work, refConfId=0,
                                                 probeConfId=conf.GetId(), useColors=colors,
                                                 opt_param=balance)  # fmt: skip
            transform, _ = kabsch.fit(conf.GetPositions(), moved.GetPositions())
            if colors:
                scores = {COMBO: shape + color, SHAPE_TANIMOTO: shape, COLOR_TANIMOTO: color}
            else:
                scores = {SHAPE_TANIMOTO: shape}
            fits.append(Fit(conf.GetId(), transform, scores))
    return fits


def score(probe: Chem.Mol, reference: Chem.Mol, conf_id: int) -> float:
    """The shape Tanimoto of *probe*'s conformer where it is (no alignment)."""
    options = rdShapeAlign.ShapeInputOptions()
    options.useColors = False
    with rdBase.BlockLogs():
        shape, _color = rdShapeAlign.ScoreMol(reference, probe, options, options, 0, conf_id)
    return float(shape)
