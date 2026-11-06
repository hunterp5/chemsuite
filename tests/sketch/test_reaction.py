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

"""Reactions in the drawing (REACTIONS_PLAN R1): components, RDKit, RXN and reaction SMILES."""

from __future__ import annotations

import pytest
from rdkit.Chem import rdChemReactions

from chemsuite.sketch import convert, reaction
from chemsuite.sketch.edits import History, map_pair, next_mapping
from chemsuite.sketch.model import Arrow, ArrowKind, Plus, Sketch

ESTER = "CC(=O)O.OCC>[H+]>CC(=O)OCC"
AMIDE = "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>CN(C)C=O>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]"
SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>[Pd]>c1ccc(-c2ccccc2)cc1"


def rdkit_smiles(text: str) -> str:
    return rdChemReactions.ReactionToSmiles(
        rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
    )


@pytest.mark.parametrize("text", [ESTER, AMIDE, SUZUKI])
def test_reaction_smiles_round_trips_through_the_drawing_and_rxn(text):
    sketch = reaction.from_text(text)
    found = reaction.components(sketch)
    assert (len(found.reactants), len(found.agents), len(found.products)) == (2, 1, 1)
    assert reaction.to_smiles(sketch) == rdkit_smiles(text)
    block = reaction.to_rxn_block(sketch)
    assert block.startswith("$RXN")
    assert reaction.to_smiles(reaction.from_text(block)) == rdkit_smiles(text)


def test_lone_ions_lay_out_with_coordgen_preferred():
    from rdkit.Chem import rdDepictor

    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(True)  # as core/depict.py sets it app-wide
    try:
        assert reaction.to_smiles(reaction.from_text(ESTER)) == rdkit_smiles(ESTER)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)


def test_atom_maps_are_kept():
    sketch = reaction.from_text(AMIDE)
    assert sorted(reaction.mapped_atoms(sketch)) == [1, 2, 3, 4, 5, 6]
    assert ":5]" in reaction.to_smiles(sketch)


def _drawn(*parts, arrow=(4.0, 0.0, 7.0, 0.0), kind=ArrowKind.PLAIN, pluses=()):
    """Molecules (SMILES, centre x, centre y) placed around one arrow."""
    sketch = Sketch()
    for smiles, x, y in parts:
        part = convert.placed(convert.from_text(smiles), (x, y))
        reaction._add(sketch, part, 0.0, 0.0)
    sketch.arrows[sketch.new_id()] = Arrow(*arrow, kind)
    for x, y in pluses:
        sketch.pluses[sketch.new_id()] = Plus(x, y)
    return sketch


def test_components_follow_plus_signs_and_the_arrow():
    # Ketcher getComponents: fragments between barriers are one component
    salt = _drawn(("CCN", 0, 0), ("Cl", 1.5, 0), ("CCO", 10, 0))
    found = reaction.components(salt)
    assert len(found.reactants) == 1 and len(found.products) == 1
    two = _drawn(("CCN", 0, 0), ("Cl", 2.5, 0), ("CCO", 10, 0), pluses=[(1.5, 0)])
    assert len(reaction.components(two).reactants) == 2


def test_agents_sit_over_the_arrow():
    sketch = _drawn(("CCN", 0, 0), ("O", 5.5, 1.2), ("CCO", 10, 0))
    found = reaction.components(sketch)
    assert len(found.agents) == 1
    assert reaction.to_smiles(sketch) == "CCN>O>CCO"


def test_arrow_direction_and_retrosynthesis():
    backwards = _drawn(("CCN", 0, 0), ("CCO", 10, 0), arrow=(7.0, 0.0, 4.0, 0.0))
    assert reaction.to_smiles(backwards) == "CCO>>CCN"
    vertical = _drawn(("CCN", 0, 6), ("CCO", 0, -6), arrow=(0.0, 2.0, 0.0, -2.0))
    assert reaction.to_smiles(vertical) == "CCN>>CCO"
    retro = _drawn(("CCN", 0, 0), ("CCO", 10, 0), kind=ArrowKind.RETROSYNTHETIC)
    assert reaction.to_smiles(retro) == "CCO>>CCN"  # target => precursors


def test_not_one_reaction():
    with pytest.raises(reaction.ReactionError):
        reaction.components(convert.from_text("CCO"))
    two = _drawn(("CCN", 0, 0), ("CCO", 10, 0))
    two.arrows[two.new_id()] = Arrow(12, 0, 14, 0)
    with pytest.raises(reaction.ReactionError):
        reaction.components(two)


def test_reaction_text_is_recognised():
    assert reaction.looks_like_reaction(ESTER)
    assert reaction.looks_like_reaction("[C:1]>>[C:1]O")
    assert not reaction.looks_like_reaction("CCO")
    assert not reaction.looks_like_reaction("[N+](=O)[O-]")
    assert reaction.reaction_from_text("C>>>>C") is None


def test_layout_puts_reactants_left_and_agents_above():
    sketch = reaction.from_text(ESTER)
    arrow = reaction.the_arrow(sketch)
    found = reaction.components(sketch)
    assert len(sketch.pluses) == 1
    for part in found.reactants:
        assert max(sketch.atoms[a].x for a in part) < arrow.x0
    for part in found.products:
        assert min(sketch.atoms[a].x for a in part) > arrow.x1
    (agent,) = found.agents
    assert min(sketch.atoms[a].y for a in agent) > arrow.y0


def test_arrows_pluses_and_maps_are_undoable_edits():
    history = History()
    change = history.change("Arrow")
    change.add_arrow(Arrow(0, 0, 3, 0))
    change.add_plus(Plus(-2, 0))
    history.commit(change)
    assert history.sketch.is_reaction and not history.sketch.is_empty()
    history.undo()
    assert history.sketch.is_empty()
    history = History(convert.from_text("CO.CO"))
    a, _, c, _ = sorted(history.sketch.atoms)
    change = history.change("Map")
    assert map_pair(change, a, c) == 1
    history.commit(change)
    assert next_mapping(history.sketch) == 2
    assert history.sketch.atoms[c].mapping == 1


def test_paste_keeps_arrows_and_query_topology():
    from chemsuite.sketch.edits import insert

    source = reaction.from_text(ESTER)
    history = History()
    change = history.change("Paste")
    insert(change, source)
    history.commit(change)
    assert len(history.sketch.arrows) == 1 and len(history.sketch.pluses) == 1
    assert reaction.to_smiles(history.sketch) == rdkit_smiles(ESTER)


def test_map_numbers_keep_clear_of_hydrogen_labels():
    """The map number of an NH2 at a chain end sits off its hydrogens (left of or above/below
    the N), not over the H2 to its right."""
    from chemsuite.sketch import render, valence

    sketch = convert.from_text("[CH3:1][NH2:2]")
    n = next(a for a, atom in sketch.atoms.items() if atom.symbol == "N")
    c = next(a for a in sketch.atoms if a != n)
    # put the NH2 on the left so its hydrogens run to the right
    sketch.atoms[n] = sketch.atoms[n].moved(0.0, 0.0)
    sketch.atoms[c] = sketch.atoms[c].moved(-1.0, 0.0)
    scene = render.build(sketch, valence.check(sketch))
    (x, y) = next(pos for aid, pos, _ in scene.maps if aid == n)
    assert not (x > 0.2 and abs(y) < 0.3)
