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
"""A drawing as a file (the Sketcher's File ▸ Save…): a molfile, SD file or SMILES for a
structure, an RXN file or reaction SMILES for a reaction."""

from __future__ import annotations

from chemsuite.sketch import convert, reaction
from chemsuite.sketch.model import Sketch

#: file types offered, by suffix (the first is the default)
STRUCTURE_FILES = {".mol": "MDL Molfile", ".sdf": "SD File", ".smi": "SMILES"}
REACTION_FILES = {".rxn": "MDL RXN File", ".smi": "Reaction SMILES"}


def is_reaction(sketch: Sketch) -> bool:
    return bool(sketch.is_reaction or sketch.pluses)


def file_types(sketch: Sketch) -> dict[str, str]:
    """The file types this drawing can be saved as, by suffix."""
    return REACTION_FILES if is_reaction(sketch) else STRUCTURE_FILES


def file_text(sketch: Sketch, suffix: str) -> str:
    """The drawing as the text of a *suffix* file (one of :func:`file_types`). Raises
    ``ValueError`` for another suffix, and ``SketchError`` or ``ReactionError`` for a drawing
    that cannot be written that way."""
    suffix = suffix.lower()
    if suffix not in file_types(sketch):
        raise ValueError(f"this drawing cannot be saved as {suffix or 'a file with no suffix'}")
    if suffix == ".rxn":
        return reaction.to_rxn_block(sketch)
    if suffix == ".smi":
        if is_reaction(sketch):
            return reaction.to_smiles(sketch) + "\n"
        return convert.to_smiles(sketch) + "\n"
    block = convert.to_molblock(sketch)
    return block + "$$$$\n" if suffix == ".sdf" else block
