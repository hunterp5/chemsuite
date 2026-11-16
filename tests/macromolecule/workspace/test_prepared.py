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
"""What a preparation leaves besides the structure: its saved form, whether it still
describes the structure, and the built atoms carried over a second preparation."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from chemsuite.core.assets import AssetStore
from chemsuite.macromolecule.structure import edit
from chemsuite.macromolecule.workspace.prepared import (
    Minimisation,
    PreparedRecord,
    atoms_digest,
    carried_built,
)
from chemsuite.macromolecule.workspace.workspace import VERSION, MacromoleculeWorkspace

from .test_session import _add, structure


def _record(s=None) -> PreparedRecord:
    s = s if s is not None else structure("1crn.cif.gz")
    return PreparedRecord.of("1CRN", s, ph=6.5, settings={"caps": False},
                             states={0: "N-neutral", 2: "CYM"}, pka={2: 8.1},
                             population={2: 0.3}, built=[int(s.uid[5]), int(s.uid[6])])  # fmt: skip


def test_a_record_names_residues_by_label_and_survives_its_saved_form() -> None:
    s = structure("1crn.cif.gz")
    record = _record(s)
    assert record.states == {s.residue_label(0): "N-neutral", s.residue_label(2): "CYM"}
    assert record.states_on(s) == {0: "N-neutral", 2: "CYM"}
    assert PreparedRecord.from_json(record.to_json()) == record
    minimised = Minimisation("1CRN-prep", {"protein": "amber14"}, {"uid": [1]}, {"rmsd": 0.2})
    full = replace(record, minimised=minimised)
    assert PreparedRecord.from_json(full.to_json()) == full


def test_damaged_records_are_ignored() -> None:
    good = _record().to_json()
    assert PreparedRecord.from_json("1CRN") is None
    assert PreparedRecord.from_json({**good, "atoms": 3}) is None
    assert PreparedRecord.from_json({**good, "states": {"A/CYS 3": 4}}) is None
    assert PreparedRecord.from_json({**good, "built": [1, "x"]}) is None
    assert PreparedRecord.from_json({**good, "minimised": {"source": 1}}) is None
    bare = PreparedRecord.from_json({"protein": "1CRN", "atoms": "abc"})
    assert bare is not None and bare.states == {} and bare.built == frozenset()


def test_an_edit_leaves_the_record_describing_another_structure() -> None:
    s = structure("1crn.cif.gz")
    record = _record(s)
    assert record.matches(s)
    moved = s.with_coords(s.coords + 1.0)
    assert record.matches(moved)  # moving atoms (minimising them) is not an edit
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[10] = True
    assert not record.matches(edit.delete(s, gone))
    assert atoms_digest(s) == record.atoms


def test_built_atoms_are_carried_over_a_second_preparation() -> None:
    s = structure("1crn.cif.gz")
    first = _record(s)
    assert carried_built([int(s.uid[9])], first, s) == {int(s.uid[i]) for i in (5, 6, 9)}
    gone = np.zeros(s.n_atoms, dtype=bool)
    gone[5] = True
    cut = edit.delete(s, gone)
    assert carried_built([], first, cut) == {int(s.uid[6])}  # a deleted atom drops out
    assert carried_built([3], None, s) == {3}


def test_the_workspace_keeps_records_with_the_document() -> None:
    store = AssetStore()
    _add(store, "1CRN", "1crn.cif.gz")
    workspace = MacromoleculeWorkspace(store, deferred=False)
    record = _record()
    workspace.set_prepared(record)
    assert workspace.prepared("1CRN") == record and workspace.revision == 1
    saved = workspace.capture()
    assert saved["prepared"] == {"1CRN": record.to_json()}
    again = MacromoleculeWorkspace(store, saved, deferred=False)
    assert again.prepared("1CRN") == record and again.revision == 0
    assert again.drop_prepared("1CRN") and not again.drop_prepared("1CRN")


def test_an_older_document_without_records_loads_and_a_gone_protein_is_not_saved() -> None:
    store = AssetStore()
    _add(store)
    saved = {"version": VERSION, "session": {"named": {}}, "view": {}}
    workspace = MacromoleculeWorkspace(store, saved, deferred=False)
    assert workspace.records == {} and workspace.capture() == saved
    workspace.set_prepared(_record())  # 1CRN is not in this store
    assert "prepared" not in workspace.capture()
