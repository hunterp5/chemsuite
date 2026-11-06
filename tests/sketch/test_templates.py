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

"""Rings and fusing (SKETCHER_PLAN §4.3).

Results are checked by canonical SMILES against the named ring systems, for overlapping atoms
and for unit bond lengths, and aromaticity is RDKit's perception of the result.
"""

from __future__ import annotations

import itertools
import math

import pytest
from rdkit import Chem

from chemsuite.sketch import templates
from chemsuite.sketch.convert import bond_length_spread, to_smiles
from chemsuite.sketch.edits import History
from chemsuite.sketch.templates import BENZENE, RING_BAR, ring
from tests.sketch.helpers import canonical

CYCLOHEXANE = ring(6)
CYCLOPENTANE = ring(5)


def _fresh(template=BENZENE):
    history = History()
    change = history.change("ring")
    atoms = templates.on_canvas(change, template, (0.0, 0.0))
    history.commit(change)
    return history, atoms


def _apply(history, fn):
    change = history.change("t")
    result = fn(change)
    history.commit(change)
    return result


def _no_overlaps(sketch):
    pts = [a.pos for a in sketch.atoms.values()]
    return all(math.dist(p, q) > 0.5 for p, q in itertools.combinations(pts, 2))


def _bond_where(sketch, predicate):
    for bid, bond in sketch.bonds.items():
        if predicate(bid, bond):
            return bid
    raise LookupError


def _ring_bonds_not_fused(sketch):
    """Bonds whose two atoms each have exactly two neighbours (peripheral ring bonds)."""
    return [
        bid
        for bid, b in sketch.bonds.items()
        if sketch.degree(b.begin) == 2 and sketch.degree(b.end) == 2
    ]


@pytest.mark.parametrize(
    ("template", "smiles"),
    [
        (RING_BAR[0], "c1ccccc1"),
        (RING_BAR[1], "C1CC1"),
        (RING_BAR[2], "C1CCC1"),
        (RING_BAR[3], "C1CCCC1"),
        (RING_BAR[4], "C1CCCCC1"),
        (RING_BAR[5], "C1CCCCCC1"),
        (RING_BAR[6], "C1CCCCCCC1"),
    ],
)
def test_ring_bar_on_canvas(template, smiles):
    history, atoms = _fresh(template)
    assert to_smiles(history.sketch) == canonical(smiles)
    assert bond_length_spread(history.sketch) < 1e-9
    assert len(atoms) == len(template.atoms)


def test_benzene_on_benzene_bond_gives_naphthalene_whichever_bond():
    for pick in range(6):
        history, _ = _fresh()
        bid = sorted(history.sketch.bonds)[pick]
        _apply(history, lambda ch, b=bid: templates.on_bond(ch, BENZENE, b))
        assert to_smiles(history.sketch) == canonical("c1ccc2ccccc2c1"), pick
        assert len(history.sketch.atoms) == 10
        assert _no_overlaps(history.sketch)
        assert bond_length_spread(history.sketch) < 1e-6


def test_third_benzene_gives_anthracene_or_phenanthrene_by_bond():
    """On the bond opposite the fusion: anthracene; on a bond next to it: phenanthrene."""
    seen = set()
    history, _ = _fresh()
    _apply(history, lambda ch: templates.on_bond(ch, BENZENE, min(history.sketch.bonds)))
    base = history.sketch.copy()
    for bid in _ring_bonds_not_fused(base):
        trial = History(base.copy())
        _apply(trial, lambda ch, b=bid: templates.on_bond(ch, BENZENE, b))
        assert _no_overlaps(trial.sketch)
        seen.add(to_smiles(trial.sketch))
    assert seen == {canonical("c1ccc2cc3ccccc3cc2c1"), canonical("c1ccc2c(c1)ccc1ccccc12")}


def test_fused_ring_goes_away_from_the_existing_ring():
    history, _ = _fresh(CYCLOHEXANE)
    bid = min(history.sketch.bonds)
    _apply(history, lambda ch: templates.on_bond(ch, CYCLOHEXANE, bid))
    assert to_smiles(history.sketch) == canonical("C1CCC2CCCCC2C1")  # decalin
    assert _no_overlaps(history.sketch)


@pytest.mark.parametrize(
    ("first", "second", "smiles"),
    [
        (BENZENE, CYCLOHEXANE, "c1ccc2c(c1)CCCC2"),  # tetralin
        (BENZENE, CYCLOPENTANE, "c1ccc2c(c1)CCC2"),  # indane
        (CYCLOPENTANE, BENZENE, "c1ccc2c(c1)CCC2"),
        (CYCLOHEXANE, BENZENE, "c1ccc2c(c1)CCCC2"),
    ],
)
def test_mixed_fusions(first, second, smiles):
    history, _ = _fresh(first)
    _apply(history, lambda ch: templates.on_bond(ch, second, min(history.sketch.bonds)))
    assert to_smiles(history.sketch) == canonical(smiles)
    assert _no_overlaps(history.sketch)


def test_indole_by_fusing_then_setting_nitrogen():
    from chemsuite.sketch import edits

    history, benzene = _fresh()
    five = _apply(
        history, lambda ch: templates.on_bond(ch, CYCLOPENTANE, min(history.sketch.bonds))
    )
    s = history.sketch
    fusion = [a for a in five if a in benzene]
    new = [a for a in five if a not in benzene]
    n, c3 = [a for a in new if set(s.neighbours(a)) & set(fusion)]
    (c2,) = [a for a in new if a not in (n, c3)]
    _apply(history, lambda ch: edits.set_element(ch, n, "N"))
    _apply(history, lambda ch: edits.set_bond_order(ch, s.bond_between(c2, c3), edits.Order.DOUBLE))
    assert to_smiles(s) == canonical("c1ccc2[nH]ccc2c1")


def test_steroid_core_from_four_fusions():
    """Gonane: three cyclohexanes and a cyclopentane fused angularly."""
    history, _ = _fresh(CYCLOHEXANE)
    s = history.sketch
    _apply(history, lambda ch: templates.on_bond(ch, CYCLOHEXANE, min(s.bonds)))  # decalin

    # Fuse on a peripheral bond next to a fusion atom (angular, as in phenanthrene).
    def next_to_fusion(sketch):
        """Peripheral bonds with an end next to a fusion atom (phenanthrene-style sites)."""
        fusion = {a for a in sketch.atoms if sketch.degree(a) == 3}
        for bid, b in sorted(sketch.bonds.items()):
            ends = (b.begin, b.end)
            if any(sketch.degree(e) != 2 for e in ends):
                continue
            if any(set(sketch.neighbours(e)) & fusion for e in ends):
                yield bid

    smiles = set()
    base = s.copy()
    for bid in next_to_fusion(base):
        trial = History(base.copy())
        _apply(trial, lambda ch, b=bid: templates.on_bond(ch, CYCLOHEXANE, b))
        for bid2 in list(next_to_fusion(trial.sketch)):
            t2 = History(trial.sketch.copy())
            _apply(t2, lambda ch, b=bid2: templates.on_bond(ch, CYCLOPENTANE, b))
            if _no_overlaps(t2.sketch):
                smiles.add(to_smiles(t2.sketch))
    assert canonical("C1CCC2C(C1)CCC1C2CCC2CCCC21") in smiles


def test_ring_on_atom_shares_the_atom():
    """Ketcher fromTemplateOnAtom: the clicked atom becomes the ring's attachment atom."""
    history, _ = _fresh(CYCLOHEXANE)
    atom = min(history.sketch.atoms)
    _apply(history, lambda ch: templates.on_atom(ch, CYCLOPENTANE, atom))
    assert to_smiles(history.sketch) == canonical("C1CCC2(CC1)CCCC2")  # spiro[4.5]decane
    assert _no_overlaps(history.sketch)


def test_ring_on_chain_end_and_as_a_substituent():
    from chemsuite.sketch import edits

    history = History()
    bid = _apply(history, lambda ch: edits.bond_on_canvas(ch, (0.0, 0.0)))
    end = history.sketch.bonds[bid].end
    trial = History(history.sketch.copy())
    _apply(trial, lambda ch: templates.on_atom(ch, BENZENE, end))
    assert to_smiles(trial.sketch) == canonical("Cc1ccccc1")  # the chain end is a ring atom
    trial = History(history.sketch.copy())
    _apply(trial, lambda ch: templates.on_atom(ch, BENZENE, end, extra_bond=True))
    assert to_smiles(trial.sketch) == canonical("CCc1ccccc1")  # through a new bond
    assert _no_overlaps(trial.sketch)


def test_ring_on_atom_points_into_the_free_direction():
    history, _ = _fresh(CYCLOHEXANE)
    atom = min(history.sketch.atoms)
    before = set(history.sketch.atoms)
    _apply(history, lambda ch: templates.on_atom(ch, BENZENE, atom, extra_bond=True))
    new = set(history.sketch.atoms) - before
    centre_new = [sum(history.sketch.atoms[a].pos[k] for a in new) / len(new) for k in (0, 1)]
    centre_old = [sum(history.sketch.atoms[a].pos[k] for a in before) / len(before) for k in (0, 1)]
    assert math.dist(centre_new, centre_old) == pytest.approx(3.0)  # ring radius, bond, radius


def test_fusing_is_one_undo_step():
    history, _ = _fresh()
    before = dict(history.sketch.atoms), dict(history.sketch.bonds)
    _apply(history, lambda ch: templates.on_bond(ch, BENZENE, min(history.sketch.bonds)))
    history.undo()
    assert (dict(history.sketch.atoms), dict(history.sketch.bonds)) == before


def test_shared_bond_keeps_the_higher_order():
    from chemsuite.sketch import edits

    history, _ = _fresh(CYCLOHEXANE)
    bid = min(history.sketch.bonds)
    _apply(history, lambda ch: edits.set_bond_order(ch, bid, edits.Order.DOUBLE))
    _apply(history, lambda ch: templates.on_bond(ch, CYCLOHEXANE, bid))
    assert history.sketch.bonds[bid].order is edits.Order.DOUBLE
    assert Chem.MolFromSmiles(to_smiles(history.sketch)) is not None
