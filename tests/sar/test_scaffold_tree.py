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
"""The Scaffold Tree: removing rings, the prioritization rules and the tree."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from chemsuite.core.jobs import Cancelled, CancelToken
from chemsuite.sar import scaffold_tree as st


@pytest.fixture(params=["murcko", "pruned"], autouse=True)
def scaffolds(request, monkeypatch) -> None:
    """Every test also runs with macOS's fallback, the unsanitized ``pruned``."""
    from chemsuite.descriptors import scaffold

    if request.param == "pruned":
        monkeypatch.setattr(st, "murcko", scaffold.pruned)


def _parent(smiles: str) -> tuple[str, str] | None:
    return st.parent_of(st.scaffold_of(smiles))


def _key(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


# Each rule deciding: the smallest such scaffold among 10,000 BindingDB compounds, or a
# constructed case where none of them needed the rule.
@pytest.mark.parametrize(
    ("child", "parent", "rule"),
    [
        # 1: a three-membered heterocycle goes first; a fused epoxide leaves its alkene
        ("O=C(Cc1ccccc1)OC1CC2[NH2+]C(C1)C1OC21", "O=C(Cc1ccccc1)OC1CC2C=CC(C1)[NH2+]2", "1"),
        ("c1ccc(C2CN2)cc1", "c1ccccc1", "1"),
        # 2: a macrocycle stays while smaller rings can go
        ("O=C1OCCCCCCCCCCc2ccccc21", "O=C1C=CCCCCCCCCCCO1", "2"),
        # 3: the parent with fewer acyclic linker bonds
        ("c1cc2nonc2cc1N1CCCC1", "c1ccc2nonc2c1", "3"),
        ("c1ccc(CC2CC3CCC2C3)cc1", "C1CC2CCC1C2", "3"),
        # 4: bridged (and spiro) systems are kept: the largest |delta|
        ("c1ccc2c(c1)C1CCC2C1", "C1=CC2CCC1C2", "4"),
        # 5: bridged before spiro
        ("C1CC2(C1)CC1CCC2C1", "C1CC2CCC1C2", "5"),
        # 6: rings of 3, 5 and 6 atoms go first
        ("c1ccc(C2CCC2)cc1", "C1CCC1", "6"),
        ("O=C1CC2SCCN12", "O=C1CCN1", "6"),
        # 7 (chemsuite's reading): the aromatic system left behind stays aromatic
        ("C=c1sc2ncnn2c1=O", "c1nc[nH]n1", "7"),
        # 8: the ring with fewer heteroatoms goes
        ("c1nc2sccc2s1", "c1cscn1", "8"),
        ("c1ccc2ncccc2c1", "c1ccncc1", "8"),
        ("O=c1ccc2ccccc2[nH]1", "O=c1cccc[nH]1", "8"),
        # 9: then nitrogen is kept before oxygen before sulfur
        ("c1nn2cnnc2s1", "c1nnc[nH]1", "9a"),
        ("c1csc(-c2ccoc2)c1", "c1ccoc1", "9b"),
        ("c1csc(-c2cc[se]c2)c1", "c1ccsc1", "9c"),
        # 10: the smaller ring goes
        ("C1CCC2CC2C1", "C1CCCCC1", "10"),
        # 11: aromatic rings go before non-aromatic ones
        ("c1cc2c(cn1)CNCC2", "C1=CCNCC1", "11"),
        ("c1ccc2c(c1)CCNC2", "C1=CCNCC1", "8"),  # a fused ring leaves its double bond
        # 12: a ring whose linker ends on a ring heteroatom goes first
        ("c1ccc(-c2nncn2-c2ccccc2)cc1", "c1ccc(-c2nnc[nH]2)cc1", "12"),
        # 13: the tie-break, the smallest canonical SMILES
        ("C1=CC2CCC1C2", "C1=CCCC1", "13"),
    ],
)
def test_each_rule_decides(child, parent, rule) -> None:
    assert _parent(child) == (_key(parent), rule)


def test_ring_removal() -> None:
    # fused: the shared atoms stay; symmetric choices are one candidate
    assert _parent("c1ccc2ccccc2c1") == ("c1ccccc1", "only")
    # spiro: the spiro atom stays in the ring that is kept
    assert _parent("C1CCC2(CC1)CCNC2") == (_key("C1CCNC1"), "8")
    # an aromatic nitrogen that loses its substituent takes a hydrogen
    assert _parent("c1ccc(-n2cccc2)cc1")[0] == _key("c1cc[nH]c1")
    # a linker is cut back to the ring that stays, with its exolinker C=O
    assert _parent("O=C(c1ccccc1)C1CCNCC1")[0] == _key("C1CCNCC1")
    # a ring atom double-bonded to the ring that stays is kept as its exocyclic =CH2
    assert _parent("C1=CC(=C2CCNCC2)C=CS1") == (_key("C=C1CCNCC1"), "9a")
    assert _parent("O=C1C=CNC1=C1NC=CC1=O") == (_key("C=C1NC=CC1=O"), "only")
    # an exocyclic C=O on a ring that stays is kept
    assert _parent("O=C1CCCN1c1ccccc1")[0] == _key("O=C1CCCN1")
    # phenanthrene's middle ring cannot go (the rest would fall apart): only naphthalene
    phenanthrene = Chem.MolFromSmiles("c1ccc2c(c1)ccc1ccccc12")
    assert {c.smiles for c in st.parents(phenanthrene)} == {_key("c1ccc2ccccc2c1")}
    # a cage no ring can be removed from is a root, with more than one ring
    assert _parent("C1C2CC3CC1CC(C2)C3") is None


def test_scaffolds_and_normalizing() -> None:
    assert st.scaffold_of("CCO") is None and st.scaffold_of("") is None
    assert st.scaffold_of(None) is None and st.scaffold_of("not a molecule") is None
    # side chains go, stereo goes, exocyclic double bonds stay
    assert st.scaffold_of("C[C@H](N)c1ccc(cc1)C(=O)N1CCCC1") == _key("O=C(c1ccccc1)N1CCCC1")
    assert st.scaffold_of("Cn1ccc(=O)[nH]c1=O") == _key("O=c1cc[nH]c(=O)[nH]1")


def test_a_compound_scaffold_and_a_pruned_one_are_one_node() -> None:
    cells = ["CCc1ccc2ccc(=O)[nH]c2c1", "Cc1cc(=O)[nH]cc1", "CCO", None]
    tree = st.build(cells)
    pyridone = tree.smiles.index(_key("O=c1cccc[nH]1"))
    assert tree.node.tolist()[1] == pyridone and tree.node.tolist()[2:] == [-1, -1]
    assert tree.parent[tree.node[0]] == pyridone and tree.rule[tree.node[0]] == "8"
    assert tree.subtree_counts()[pyridone] == 2 and tree.own_counts()[pyridone] == 1


def test_the_tree() -> None:
    cells = [
        "c1ccc2[nH]c(-c3ccccc3)cc2c1",  # 2-phenylindole
        "c1ccc2[nH]ccc2c1",  # indole
        "c1ccc2[nH]ccc2c1",
        "c1ccc(-c2ccccn2)cc1",  # 2-phenylpyridine
        "C1C2CC3CC1CC(C2)C3",  # adamantane: a root (four rings in RDKit's symmetrized SSSR)
        "CCCC",
    ]
    tree = st.build(cells)
    sub, own = tree.subtree_counts(), tree.own_counts()
    assert own.sum() == 5 and sum(sub[r] for r in tree.roots()) == 5
    pyrrole = tree.smiles.index("c1cc[nH]c1")
    indole = tree.smiles.index(_key("c1ccc2[nH]ccc2c1"))
    assert tree.parent[indole] == pyrrole and tree.depth(tree.node[0]) == 2
    assert sub[pyrrole] == 3 and own[indole] == 2
    assert sorted(tree.subtree_rows(pyrrole).tolist()) == [0, 1, 2]
    assert tree.children(pyrrole) == [indole]
    roots = {tree.smiles[r]: int(tree.rings[r]) for r in tree.roots()}
    assert roots == {"c1cc[nH]c1": 1, "c1ccncc1": 1, _key("C1C2CC3CC1CC(C2)C3"): 4}
    for i in range(len(tree)):  # every chain ends at a root
        assert tree.depth(i) < len(tree)


def test_the_tree_is_the_same_with_the_scaffold_fallback(monkeypatch) -> None:
    """On macOS ``murcko`` falls back to the unsanitized ``pruned``: the same nodes."""
    from chemsuite.descriptors import scaffold

    cells = [
        "CC(=O)Nc1ccc(O)cc1",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "c1ccc2[nH]c(-c3ccccc3)cc2c1",
        "O=C(O)c1ccccc1OC(C)=O",
        "CN1CCN(Cc2ccc(C(=O)Nc3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc2)CC1",
    ]
    before = st.build(cells)
    monkeypatch.setattr(st, "murcko", scaffold.pruned)
    after = st.build(cells)
    assert after.smiles == before.smiles and after.parent.tolist() == before.parent.tolist()


def test_cancel() -> None:
    token = CancelToken()
    token.cancel()
    with pytest.raises(Cancelled):
        st.build(["c1ccccc1"] * 10, cancel=token)


def test_progress_reaches_the_end() -> None:
    seen = []
    cells = ["c1ccc2[nH]c(-c3ccccc3)cc2c1", "c1ccccc1"]  # parents over three levels
    st.build(cells, progress=lambda d, t: seen.append((d, t)))
    done = [d for d, _ in seen]
    assert seen[-1][0] == seen[-1][1] and np.all(np.diff(done) >= 0)
    assert len(set(done)) >= 4  # each level of parents moves it on, not all at once
