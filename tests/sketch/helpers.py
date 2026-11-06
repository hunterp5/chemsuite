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

"""Shared helpers for the sketcher tests."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

from rdkit import Chem

from chemsuite.sketch.model import Atom, Bond, Order, Sketch

REVIEW = Path(__file__).resolve().parents[2] / "scripts" / "review"


def valence_review():
    """scripts/review/valence_review.py, which holds the P0 reference species."""
    name = "valence_review"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, REVIEW / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


_ORDERS = {
    Chem.BondType.SINGLE: Order.SINGLE,
    Chem.BondType.DOUBLE: Order.DOUBLE,
    Chem.BondType.TRIPLE: Order.TRIPLE,
    Chem.BondType.AROMATIC: Order.AROMATIC,
}


def drawn(smiles: str) -> tuple[Sketch, list[int]]:
    """``smiles`` as a user would draw it: heavy atoms, charges and bonds, no hydrogen counts.

    Coordinates are a circle (valence does not depend on them). Returns the sketch and the
    atom ids in SMILES order.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    sketch = Sketch()
    ids = []
    n = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        t = 2 * math.pi * atom.GetIdx() / max(n, 1)
        aid = sketch.new_id()
        sketch.atoms[aid] = Atom(
            n * math.cos(t),
            n * math.sin(t),
            symbol=atom.GetSymbol(),
            charge=atom.GetFormalCharge(),
            isotope=atom.GetIsotope(),
            radical=atom.GetNumRadicalElectrons(),
        )
        ids.append(aid)
    for bond in mol.GetBonds():
        sketch.bonds[sketch.new_id()] = Bond(
            ids[bond.GetBeginAtomIdx()], ids[bond.GetEndAtomIdx()], _ORDERS[bond.GetBondType()]
        )
    return sketch, ids


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def degrees(angle: float) -> float:
    return math.degrees(angle) % 360
