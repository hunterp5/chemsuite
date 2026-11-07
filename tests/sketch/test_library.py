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

"""The template library: every template is the molecule it names, laid out cleanly, and
places like the ring bar's rings (on canvas, on an atom, fused onto a bond)."""

from __future__ import annotations

import itertools
import math

import pytest
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from chemsuite.sketch import convert, library, templates
from chemsuite.sketch.edits import History
from chemsuite.sketch.templates import BENZENE

NAMES = [e.name for e in library.TABLE]
#: two atoms drawn closer than this (in bond lengths) would overlap; bicyclo[1.1.1]pentane's
#: bridgeheads are drawn close in every 2D depiction, so it has its own limit
CLOSEST = 0.5
TIGHT = {"Bicyclo[1.1.1]pentane": 0.4}


def _placed(fn):
    history = History()
    change = history.change("t")
    result = fn(change)
    history.commit(change)
    return history, result


def test_names_are_unique_and_groups_are_known() -> None:
    assert len(NAMES) == len(set(NAMES)) >= 60
    assert set(library.GROUPS) == {e.group for e in library.TABLE}


@pytest.mark.parametrize("name", NAMES)
def test_each_template_is_the_molecule_it_names(name) -> None:
    entry = library.BY_NAME[name]
    template = library.template(name)
    history, _atoms = _placed(lambda ch: templates.on_canvas(ch, template, (0.0, 0.0)))
    mol = convert.to_rdkit(history.sketch)
    expected = Chem.MolFromSmiles(entry.smiles)
    assert Chem.MolToInchiKey(mol) == Chem.MolToInchiKey(expected)
    if entry.formula:  # worked out by hand, not by RDKit
        assert CalcMolFormula(mol) == entry.formula
    points = [(x, y) for x, y, _symbol in template.atoms]
    lengths = [math.dist(points[i], points[j]) for i, j, _order in template.bonds]
    assert min(lengths) > 0.6 and max(lengths) < 1.6  # Clean's bond length (1)
    bonded = {frozenset((i, j)) for i, j, _order in template.bonds}
    closest = min(
        math.dist(points[i], points[j])
        for i, j in itertools.combinations(range(len(points)), 2)
        if frozenset((i, j)) not in bonded
    )
    assert closest > TIGHT.get(name, CLOSEST)


def test_attachment_is_the_first_smiles_atom_and_a_ring_bond() -> None:
    template = library.template("Pyridine")
    assert template.atoms[0][2] == "C" and template.aromatic
    i, j, _order = template.bonds[template.attach_bond]
    assert 0 in (i, j)


def test_library_templates_attach_and_fuse_like_rings() -> None:
    pyridine = library.template("Pyridine")
    # dropped on an atom with a drag: a substituent through its 4-position (4-phenylpyridine)
    history, atoms = _placed(lambda ch: templates.on_canvas(ch, BENZENE, (0.0, 0.0)))
    change = history.change("substituent")
    templates.on_atom(change, pyridine, atoms[0], extra_bond=True)
    history.commit(change)
    assert Chem.MolToSmiles(convert.to_rdkit(history.sketch)) == Chem.MolToSmiles(
        Chem.MolFromSmiles("c1ccc(-c2ccncc2)cc1")
    )
    # pyrrole fuses through its C3–C4 bond (the first ring bond at its first atom): isoindole,
    # re-aromatised across the fusion (indole has its own template)
    history, atoms = _placed(lambda ch: templates.on_canvas(ch, BENZENE, (0.0, 0.0)))
    bond = next(iter(history.sketch.bonds))
    change = history.change("fuse")
    templates.on_bond(change, library.template("Pyrrole"), bond)
    history.commit(change)
    assert Chem.MolToSmiles(convert.to_rdkit(history.sketch)) == Chem.MolToSmiles(
        Chem.MolFromSmiles("c1ccc2c[nH]cc2c1")
    )


def test_search_by_name() -> None:
    assert [e.name for e in library.search("azaspiro")] == [
        "2-Azaspiro[3.3]heptane",
        "2-Oxa-6-azaspiro[3.3]heptane",
    ]
    assert library.search("") == list(library.TABLE)
