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
"""scripts/review/prepare_truth.py: the P0 ground truth read from neutron structures
(docs/prepare/p0-review.md). The structures here are made up: atoms far apart, each hydrogen
1 Å from the atom it belongs to, which is all the reading looks at."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import gemmi
import pytest

REVIEW = Path(__file__).resolve().parents[3] / "scripts" / "review"


def _module():
    name = "prepare_truth"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, REVIEW / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


pt = _module()

#: residue -> its atoms as (name, element, x offset, occupancy, altloc); heavy atoms 4 Å apart
RESIDUES = {
    ("HIS", 1): [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("CG", "C"),
                 ("ND1", "N"), ("CD2", "C"), ("CE1", "C"), ("NE2", "N")],
    ("ASP", 2): [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("CG", "C"),
                 ("OD1", "O"), ("OD2", "O")],
    ("ASN", 3): [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("CG", "C"),
                 ("OD1", "O"), ("ND2", "N")],
    ("LYS", 4): [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("CG", "C"),
                 ("CD", "C"), ("CE", "C"), ("NZ", "N")],
}  # fmt: skip


def write(path: Path, hydrogens: dict[tuple[int, str], list[tuple[str, str, float, str]]]):
    """A one-chain structure; *hydrogens* puts (name, element, occupancy, altloc) hydrogens
    on (residue number, atom)."""
    structure = gemmi.Structure()
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")
    x = 0.0
    for (name, number), atoms in RESIDUES.items():
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(number, " ")
        residue.entity_type = gemmi.EntityType.Polymer
        for atom_name, element in atoms:
            atom = gemmi.Atom()
            atom.name, atom.element, atom.occ = atom_name, gemmi.Element(element), 1.0
            atom.pos = gemmi.Position(x, 0, 0)
            residue.add_atom(atom)
            site = -1
            for h_name, h_element, occ, altloc in hydrogens.get((number, atom_name), []):
                site += altloc != "B"  # an H/D pair (altlocs A and B) shares one site
                h = gemmi.Atom()
                h.name, h.element, h.occ = h_name, gemmi.Element(h_element), occ
                h.altloc = altloc
                angle = 2 * math.pi * site / 3
                h.pos = gemmi.Position(x, math.cos(angle), math.sin(angle))
                residue.add_atom(h)
            x += 4.0
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    structure.make_mmcif_document().write_file(str(path))
    return path


def test_state_names() -> None:
    assert pt.state_name("HIS", {"ND1": 1, "NE2": 0}) == "HID"
    assert pt.state_name("HIS", {"ND1": 0, "NE2": 1}) == "HIE"
    assert pt.state_name("HIS", {"ND1": 1, "NE2": 1}) == "HIP"
    assert pt.state_name("HIS", {"ND1": 0, "NE2": 0}) is None  # undecided, not scored
    assert pt.state_name("ASP", {"OD1": 0, "OD2": 1}) == "ASH"
    assert pt.state_name("LYS", {"NZ": 2}) == "LYN"
    assert pt.state_name("CYS", {"SG": 0}) == "CYM"


def test_truth_reads_hd_pairs_by_summed_occupancy(tmp_path) -> None:
    path = write(tmp_path / "x.cif", {
        # an H/D pair at ND1 adding up to 0.9; NE2's deuterium refined to 0.2
        (1, "ND1"): [("HD1", "H", 0.3, "A"), ("DD1", "D", 0.6, "B")],
        (1, "NE2"): [("DE2", "D", 0.2, "\0")],
        (2, "OD2"): [("DD2", "D", 0.8, "\0")],
        (3, "ND2"): [("DD21", "D", 1.0, "\0"), ("DD22", "D", 1.0, "\0")],
        (4, "NZ"): [("DZ1", "D", 1.0, "\0"), ("DZ2", "D", 1.0, "\0"), ("DZ3", "D", 0.1, "\0")],
    })  # fmt: skip
    found = pt.truth(path)
    states = {key[1]: state for (key, _name), state in found.states.items()}
    assert states == {1: "HID", 2: "ASH", 4: "LYN"}
    assert {key[1] for key, _ in found.flips} == {1, 3}  # the amide N carries two D


def test_unmodelled_hydrogens_are_not_read_as_absent(tmp_path) -> None:
    # two D on NZ and no third site at all: a spinning NH3, not a neutral Lys
    path = write(tmp_path / "x.cif", {
        (4, "NZ"): [("DZ1", "D", 1.0, "\0"), ("DZ2", "D", 1.0, "\0")],
    })  # fmt: skip
    found = pt.truth(path)
    assert not [s for s in found.states if s[1] == "LYS"]


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_input_flips_are_scored_against_the_deposit(tmp_path, seed) -> None:
    path = write(tmp_path / "x.cif", {(3, "ND2"): [("DD21", "D", 1.0, "\0"),
                                                   ("DD22", "D", 1.0, "\0")]})  # fmt: skip
    found = pt.truth(path)
    out = tmp_path / "in.pdb"
    flipped = pt.make_input(path, out, seed)
    assert "DD21" not in out.read_text()  # heavy atoms only
    seen = pt.observe(out, found)
    for site, ok in seen.flips.items():
        assert ok is (not flipped[site])  # an untouched input is right only where not flipped
    tallies = pt.score([(found, seen)])
    assert tallies["ASN orientation"].total == 1
