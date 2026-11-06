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

"""Converting between reaction text, RDKit reactions and the binary a reaction column stores.

The counterpart of :mod:`chemsuite.core.mols` for :class:`~chemsuite.core.table.columns.
ReactionColumn`. A cell's text is reaction SMILES (or SMARTS, or an RXN file block); its blob is
RDKit's ``ChemicalReaction`` binary, which keeps the components, agents, atom maps and any
coordinates read from an RXN file.
"""

from __future__ import annotations

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from chemsuite.sketch import reaction as sketch_reaction
from chemsuite.sketch.reaction import looks_like_reaction


def looks_like_rxn_block(text: str) -> bool:
    return text.lstrip().startswith("$RXN")


def _templates(rxn: rdChemReactions.ChemicalReaction) -> list[Chem.Mol]:
    return [*rxn.GetReactants(), *rxn.GetAgents(), *rxn.GetProducts()]


def _valid_smiles(rxn: rdChemReactions.ChemicalReaction) -> bool:
    """Whether every component of a reaction read as SMILES is a valid molecule (RDKit reads
    reaction templates without checking valences, so ``C(C)(C)(C)(C)C>>C`` would pass)."""
    with rdBase.BlockLogs():
        for template in _templates(rxn):
            try:
                Chem.SanitizeMol(Chem.Mol(template))
            except (ValueError, RuntimeError):
                return False
    return True


def reaction_from_text(text: str) -> rdChemReactions.ChemicalReaction | None:
    """Reaction SMILES, reaction SMARTS or an RXN block as an RDKit reaction; None if the text
    is not a valid reaction."""
    text = text.strip() if not looks_like_rxn_block(text) else text.rstrip()
    rxn = _read(text)
    if rxn is None or rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates() == 0:
        return None
    if not looks_like_rxn_block(text) and _is_smiles(text) and not _valid_smiles(rxn):
        return None
    return rxn


def _read(text: str) -> rdChemReactions.ChemicalReaction | None:
    if not looks_like_reaction(text):
        return None
    return sketch_reaction.reaction_from_text(text)


def _is_smiles(text: str) -> bool:
    """Whether a one-line reaction reads as SMILES (not only as SMARTS)."""
    with rdBase.BlockLogs():
        try:
            rdChemReactions.ReactionFromSmarts(text.split()[0], useSmiles=True)
        except ValueError:
            return False
    return True


def reaction_to_blob(rxn: rdChemReactions.ChemicalReaction) -> bytes:
    return rxn.ToBinary()


def reaction_from_blob(blob: bytes) -> rdChemReactions.ChemicalReaction:
    return rdChemReactions.ChemicalReaction(blob)


def reaction_smiles(rxn: rdChemReactions.ChemicalReaction) -> str:
    """Reaction SMILES ("" if RDKit cannot write it)."""
    with rdBase.BlockLogs():
        try:
            return rdChemReactions.ReactionToSmiles(rxn)
        except (RuntimeError, ValueError):
            return ""


def rxn_block(rxn: rdChemReactions.ChemicalReaction) -> str:
    """An RXN file (V2000; V3000 when a component needs it) with agents kept separate, written
    with Kekulé bonds as molfiles are, laid out if it has no coordinates."""
    from chemsuite.sketch import reaction as sketch_reaction

    if all(t.GetNumConformers() for t in _templates(rxn)):
        copy = rdChemReactions.ChemicalReaction(rxn)
        with rdBase.BlockLogs():
            try:
                return rdChemReactions.ReactionToRxnBlock(
                    sketch_reaction._kekule(copy), separateAgents=True
                )
            except (RuntimeError, ValueError):
                pass
    return sketch_reaction.to_rxn_block(sketch_reaction.from_reaction(rxn))


def parse_reaction(text: str) -> tuple[bytes, str]:
    """``(blob, display text)`` for one reaction cell; ``(b"", text)`` when it does not parse.

    Reaction SMILES and SMARTS keep the text as written; an RXN block is shown as its reaction
    SMILES (its coordinates are kept in the blob).
    """
    rxn = reaction_from_text(text)
    if rxn is None:
        return b"", text
    shown = reaction_smiles(rxn) if looks_like_rxn_block(text) else text.strip()
    return reaction_to_blob(rxn), shown or text.strip()


def parse_reactions(texts: list[str]) -> list[bytes]:
    """Blobs for many texts (``b""`` where invalid). Picklable entry point for worker processes."""
    return [parse_reaction(text)[0] if text else b"" for text in texts]


def is_reaction_text(text: str) -> bool:
    """Whether *text* is a valid reaction (for detecting reaction columns in files)."""
    return reaction_from_text(text) is not None


__all__ = [
    "is_reaction_text",
    "looks_like_reaction",
    "looks_like_rxn_block",
    "parse_reaction",
    "parse_reactions",
    "reaction_from_blob",
    "reaction_from_text",
    "reaction_smiles",
    "reaction_to_blob",
    "rxn_block",
]
