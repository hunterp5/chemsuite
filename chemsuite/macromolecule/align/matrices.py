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
"""Amino-acid substitution matrices (BLOSUM45/62/80, PAM30/70/250) and how letters index them.

The numbers are NCBI BLAST's (``data/matrices.json``, written by ``scripts/build_matrices.py``,
with each source file's SHA-256 and citation); :data:`SHA256` pins the shipped file and a test
checks it. Letters outside the 20 standard ones follow one rule: ``B``, ``Z`` and ``X`` use
NCBI's own rows, selenocysteine ``U`` scores as cysteine, pyrrolysine ``O`` as lysine, and
anything else (``J``, unknown) as ``X``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).with_name("data") / "matrices.json"
SHA256 = "e9c5f5291a4208f2229c21f9cbc2ddc14faa6eaa04a4203d51851381fc65edf6"
DEFAULT = "BLOSUM62"
NAMES = ("BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250")

#: NCBI BLAST's default gap costs per matrix, (open, extend): a gap of n residues costs
#: ``open + n * extend``, the convention gemmi and BLAST use
GAPS = {"BLOSUM45": (15.0, 2.0), "BLOSUM62": (11.0, 1.0), "BLOSUM80": (10.0, 1.0),
        "PAM30": (9.0, 1.0), "PAM70": (10.0, 1.0), "PAM250": (14.0, 2.0)}  # fmt: skip

_ALIASES = {"U": "C", "O": "K"}


@dataclass(frozen=True, eq=False)
class Matrix:
    name: str
    letters: str
    scores: np.ndarray  # float64 (24, 24), rows and columns in ``letters`` order

    def encode(self, sequence: str) -> np.ndarray:
        """Indices into :attr:`scores` for a one-letter sequence (case ignored)."""
        index = _index(self.letters)
        unknown = index["X"]
        return np.fromiter(
            (index.get(_ALIASES.get(ch, ch), unknown) for ch in sequence.upper()),
            dtype=np.int64,
            count=len(sequence),
        )

    def pair_scores(self, a: str, b: str) -> np.ndarray:
        """The (len(a), len(b)) substitution scores of every residue pair."""
        return self.scores[np.ix_(self.encode(a), self.encode(b))]


@cache
def _index(letters: str) -> dict[str, int]:
    return {ch: i for i, ch in enumerate(letters)}


@cache
def _document() -> dict:
    return json.loads(DATA.read_text(encoding="ascii"))


@cache
def matrix(name: str = DEFAULT) -> Matrix:
    document = _document()
    scores = np.asarray(document["matrices"][name]["scores"], dtype=np.float64)
    scores.setflags(write=False)
    return Matrix(name, document["letters"], scores)
