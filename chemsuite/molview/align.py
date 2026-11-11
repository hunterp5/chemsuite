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
"""Superposing conformers, so an overlay shows how they differ rather than where they sit."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def aligned_conformers(mol: Chem.Mol) -> Chem.Mol:
    """A copy of *mol* with every conformer fitted onto the first, on the heavy atoms."""
    out = Chem.Mol(mol)
    if out.GetNumConformers() < 2:
        return out
    heavy = [a.GetIdx() for a in out.GetAtoms() if a.GetAtomicNum() > 1]
    if heavy:
        rdMolAlign.AlignMolConformers(out, atomIds=heavy)
    else:
        rdMolAlign.AlignMolConformers(out)
    return out
