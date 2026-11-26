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

"""Parts left out of calculations (the proteins list's boxes): Prepare, Minimise and Dock run
without them, and the prepared and minimized copies carry them along untouched."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from chemsuite.app.document.tools import ToolContext
from chemsuite.core.assets import PutAsset
from chemsuite.macromolecule.prepare import pipeline
from chemsuite.macromolecule.structure import parts, store
from chemsuite.macromolecule.structure.model import WATER, Structure
from chemsuite.simulate.spec import MinimiseSettings
from chemsuite.tools.macromolecule import minimise, prepare_dialog, tool
from tests.gui import gui

from .test_minimise_ui import DATA, wait
from .test_minimise_ui import fake as fake  # a fixture

pytestmark = gui


@pytest.fixture
def lysozyme(window, qtbot, monkeypatch):
    """T4 lysozyme opened, its waters, its ligand JZ4 (which the fake worker cannot charge) and
    the PO4 and BME whose bond orders are guesses (which would stop a minimization) left
    out."""
    from chemsuite.app.shell import graphics

    monkeypatch.setattr(graphics, "_problem", None)
    ctx = ToolContext(window)
    wait(qtbot, tool.open_protein(ctx, [str(DATA / "3htb.cif.gz")]))
    asset = window.workspace.table.assets["3HTB"]
    s = store.as_structure(asset.payload)
    meta = parts.with_excluded(asset.meta, s, _left_out(s))
    ctx.edit("Exclude", [PutAsset("3HTB", asset, replace(asset, meta=meta))])
    return window


def _left_out(s: Structure) -> np.ndarray:
    return (s.atom_kind == WATER) | np.isin(s.res_name[s.residue_of], ("PO4", "BME", "JZ4"))


def _hydrogens_on(s: Structure, atoms: np.ndarray) -> int:
    inside = np.zeros(s.n_atoms, bool)
    inside[atoms] = True
    hydrogen = s.element == 1
    pairs = s.bonds.astype(np.int64)
    return int(((inside[pairs[:, 0]] & hydrogen[pairs[:, 1]])
                | (inside[pairs[:, 1]] & hydrogen[pairs[:, 0]])).sum())  # fmt: skip


def test_prepare_and_minimise_carry_the_parts_left_out_untouched(lysozyme, qtbot, fake) -> None:
    ctx = ToolContext(lysozyme)
    table = lysozyme.workspace.table
    s = store.as_structure(table.assets["3HTB"].payload)
    left_out = _left_out(s)
    assert prepare_dialog.included_of(ctx, "3HTB").n_atoms == s.n_atoms - left_out.sum()
    job = prepare_dialog.prepare_protein(ctx, "3HTB", pipeline.Options(), MinimiseSettings())
    wait(qtbot, job)
    assert list(table.assets) == ["3HTB", "3HTB-prep", "3HTB-min"]
    workspace = tool.macromolecule_workspace(ctx)
    for protein in ("3HTB-prep", "3HTB-min"):
        out = store.as_structure(table.assets[protein].payload)
        found = out.indices(s.uid[left_out])
        assert (found >= 0).all(), protein
        assert np.array_equal(out.xyz()[found], s.xyz()[left_out]), protein  # not moved
        assert _hydrogens_on(out, found) == 0, protein  # no hydrogens added
        excluded = parts.excluded(table.assets[protein].meta, out)
        assert np.array_equal(np.flatnonzero(excluded), np.sort(found)), protein  # still out
        record = workspace.prepared(protein)
        assert record is not None and record.matches(out), protein
        assert (out.element[~excluded] == 1).any()  # the rest was prepared
    report = workspace.report("3HTB-min")
    assert report is not None and not any(r.text.startswith("Not minimised") for r in report.rows)
    (row,) = [r for r in report.rows if r.label == minimise.LEFT_OUT]
    assert "A/JZ4 167" in row.text and "waters" in row.text
    prep_report = workspace.report("3HTB-prep")
    assert prep_report is not None
    assert any(r.label == minimise.LEFT_OUT for r in prep_report.rows)


def test_preparing_in_place_puts_the_parts_back_where_they_were(lysozyme, qtbot) -> None:
    ctx = ToolContext(lysozyme)
    table = lysozyme.workspace.table
    s = store.as_structure(table.assets["3HTB"].payload)
    wait(qtbot, prepare_dialog.prepare_protein(ctx, "3HTB", pipeline.Options()))
    out = store.as_structure(table.assets["3HTB"].payload)
    order = [out.residue_label(r) for r in range(out.n_residues)]
    kept = [label for label in (s.residue_label(r) for r in range(s.n_residues)) if label in order]
    assert [label for label in order if label in set(kept)] == kept  # the same order
    record = tool.macromolecule_workspace(ctx).prepared("3HTB")
    assert record is not None and record.matches(out)
    assert "edited after it was prepared" not in minimise.problem(ctx, "3HTB", MinimiseSettings())
