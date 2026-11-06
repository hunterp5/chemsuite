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

"""Query bonds, attachment points and drawn queries in structure search (SKETCHER_PLAN P4).

Query bonds are CTfile bond types 5-8 and bond topology (BIOVIA CTfile Formats, bond block);
attachment points are ``M  APO``. Each query is checked by substructure searches against
reference molecules, through the same path a drawn query takes into Structure Search.
"""

from __future__ import annotations

from dataclasses import replace

import pytest
from rdkit import Chem

from chemsuite.analysis.structure_search import (
    StructureQuery,
    drawn_query_text,
    has_query_features,
    match_chunk,
    parse_query,
)
from chemsuite.sketch import convert, render, valence
from chemsuite.sketch.model import Order


def _sketch(smiles, bonds=(), atoms=()):
    """*smiles* drawn, then bond and atom edits as (index in SMILES order, field changes)."""
    sketch = convert.from_text(smiles)
    order = convert.atom_order(sketch)
    mol = Chem.MolFromSmiles(smiles)
    for (a, b), changes in bonds:
        bid = sketch.bond_between(order[a], order[b])
        sketch.bonds[bid] = replace(sketch.bonds[bid], **changes)
    for i, changes in atoms:
        sketch.atoms[order[i]] = replace(sketch.atoms[order[i]], **changes)
    assert mol.GetNumAtoms() == len(sketch.atoms)
    return sketch


def _query(sketch) -> str:
    """What Structure Search gets: the sketcher's molfile, read by the editor, made a query."""
    mol = Chem.MolFromMolBlock(convert.to_molblock(sketch))
    return drawn_query_text(mol)


def _hits(text: str, smiles: list[str]) -> list[bool]:
    cells = list(smiles)  # table cells as SMILES text, as Structure Search reads them
    parse_query(StructureQuery(text))
    return match_chunk((text, "substructure", False, cells)).tolist()


@pytest.mark.parametrize(
    "order, code",
    [(Order.SINGLE_OR_DOUBLE, 5), (Order.SINGLE_OR_AROMATIC, 6),
     (Order.DOUBLE_OR_AROMATIC, 7), (Order.ANY, 8)],
)  # fmt: skip
def test_query_bond_types_round_trip_as_ctfile_codes(order, code):
    sketch = _sketch("CCN", bonds=[((1, 2), {"order": order})])
    block = convert.to_molblock(sketch)
    assert any(ln.split()[2:4] == [str(code), "0"] for ln in block.splitlines()[4:])
    back = convert.from_text(block)
    assert sorted(b.order for b in back.bonds.values() if b.order.is_query) == [order]


@pytest.mark.parametrize("topology", [1, 2])
def test_topology_round_trips(topology):
    sketch = _sketch("CCN", bonds=[((0, 1), {"topology": topology})])
    back = convert.from_text(convert.to_molblock(sketch))
    assert sorted(b.topology for b in back.bonds.values()) == [0, topology]


def test_attachment_points_round_trip_as_apo():
    sketch = _sketch("CCC", atoms=[(0, {"attach": 1}), (2, {"attach": 3})])
    block = convert.to_molblock(sketch)
    assert any(ln.startswith("M  APO") for ln in block.splitlines())
    back = convert.from_text(block)
    assert sorted(a.attach for a in back.atoms.values()) == [0, 1, 3]
    scene = render.build(back, valence.check(back))
    assert len(scene.stubs) == 3 and sorted(t for _, t in scene.tags) == ["1", "1", "2"]


def test_query_bonds_are_drawn_dashed_with_a_tag_and_no_hydrogen_counts():
    sketch = _sketch("CCN", bonds=[((1, 2), {"order": Order.SINGLE_OR_DOUBLE, "topology": 2})])
    report = valence.check(sketch)
    assert report.valid
    nitrogen = convert.atom_order(sketch)[2]
    assert report.hydrogens[nitrogen] == 0  # its bond order is open, so is its H count
    scene = render.build(sketch, report)
    assert [t for _, t in scene.tags] == ["S/D Chn"]
    assert any(line.dashed for line in scene.lines)


# -- drawn queries in Structure Search ---------------------------------------------------------

TARGETS = ["CCN", "CC=N", "CC#N", "c1ccncc1", "C1CCNC1", "CCO"]


def test_plain_drawing_stays_smiles():
    sketch = convert.from_text("c1ccccc1C(=O)O")
    mol = Chem.MolFromMolBlock(convert.to_molblock(sketch))
    assert not has_query_features(mol)
    assert _query(sketch) == Chem.CanonSmiles("OC(=O)c1ccccc1")


@pytest.mark.parametrize(
    "order, want",
    [
        (Order.ANY, [True, True, True, True, True, False]),
        (Order.SINGLE_OR_DOUBLE, [True, True, False, False, True, False]),
        (Order.DOUBLE_OR_AROMATIC, [False, True, False, True, False, False]),
        (Order.SINGLE_OR_AROMATIC, [True, False, False, True, True, False]),
    ],
)
def test_query_bond_searches(order, want):
    sketch = _sketch("CN", bonds=[((0, 1), {"order": order})])
    assert _hits(_query(sketch), TARGETS) == want


def test_ring_and_chain_bond_searches():
    ring = _query(_sketch("CN", bonds=[((0, 1), {"topology": 1})]))
    chain = _query(_sketch("CN", bonds=[((0, 1), {"topology": 2})]))
    assert _hits(ring, TARGETS) == [False, False, False, False, True, False]
    assert _hits(chain, TARGETS) == [True, False, False, False, False, False]


def test_rgroup_and_atom_list_searches():
    from chemsuite.sketch.model import ATOM_LIST, RGROUP

    rgroup = _query(_sketch("Cc1ccccc1", atoms=[(0, {"symbol": RGROUP, "rgroup": 1})]))
    halogens = {"symbol": ATOM_LIST, "atom_list": ("Cl", "Br")}
    halo = _query(_sketch("Clc1ccccc1", atoms=[(0, halogens)]))
    targets = ["c1ccccc1", "Cc1ccccc1", "Brc1ccccc1", "Fc1ccccc1", "OCc1ccccc1"]
    assert _hits(rgroup, targets) == [False, True, True, True, True]  # R = any substituent
    assert _hits(halo, targets) == [False, False, True, False, False]


def test_query_atoms_search():
    from chemsuite.sketch.model import QUERY_ATOMS

    assert "X" in QUERY_ATOMS
    halogen = _query(_sketch("Clc1ccccc1", atoms=[(0, {"symbol": "X"})]))
    hetero = _query(_sketch("Oc1ccccc1", atoms=[(0, {"symbol": "Q"})]))
    targets = ["Fc1ccccc1", "Ic1ccccc1", "Oc1ccccc1", "Cc1ccccc1"]
    assert _hits(halogen, targets) == [True, True, False, False]
    assert _hits(hetero, targets) == [True, True, True, False]


def test_attachment_points_survive_v3000():
    from chemsuite.sketch.edits import History, set_enhanced
    from chemsuite.sketch.stereo import labels

    sketch = _sketch("C[C@H](N)CC", atoms=[(4, {"attach": 3})])
    history = History(sketch)
    change = history.change("Group")
    set_enhanced(change, list(labels(sketch).atoms), "and1")  # forces V3000
    history.commit(change)
    block = convert.to_molblock(history.sketch)
    assert "V3000" in block and "ATTCHPT=-1" in block
    back = convert.from_text(block)
    assert sorted(a.attach for a in back.atoms.values() if a.attach) == [3]


def test_smarts_with_smiles_lookalikes_is_read_as_smarts():
    parsed = parse_query(StructureQuery("[#6]~[#7]"))
    assert parsed.is_smarts
    assert parse_query(StructureQuery("CCN")).is_smarts is False
    assert parse_query(StructureQuery("C[C@H](N)O")).is_smarts is False
