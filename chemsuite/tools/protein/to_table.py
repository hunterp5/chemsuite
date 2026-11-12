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

"""Ligands and their interactions as table rows: Ligands to Table, Add Ligand to Table and
Interactions to Table, from the Protein Viewer's Ligand menu."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from chemsuite.app.document.tools import ToolContext
from chemsuite.protein.chem import profile
from chemsuite.protein.structure import store
from chemsuite.protein.structure.model import Structure

from .rows import LIGAND_COLUMN, LIGAND_ID_COLUMN, PROTEIN_COLUMN, ligand_rows, mol_of
from .tool import KIND, ccd


def extract_ligands(ctx: ToolContext, then: Callable[[], object] | None = None):
    """Ligand ▸ Ligands to Table: add rows for the ligands of every protein that has none in
    the table yet. A background job (the CCD may be asked for bond orders); *then* runs once
    the rows are in. Returns the job, or None when there is nothing to add."""
    table = ctx.table
    linked = set()
    if PROTEIN_COLUMN in table:
        linked = {table.cell(PROTEIN_COLUMN, int(o)) for o in table.oids}
    pending = [(asset.id, store.as_structure(asset.payload)) for asset in table.assets.values()
               if asset.kind == KIND and asset.id not in linked]  # fmt: skip
    if not pending:
        ctx.show_message("No new ligands to add")
        return None
    cache = ccd(ctx)

    def work(job):
        return [(protein, ligand_rows(structure, protein, ccd=cache,
                                              cancel=job.cancel_token))
                for protein, structure in pending]  # fmt: skip

    def done(found) -> None:
        added = 0
        for protein, rows in found:
            if protein in ctx.table.assets and rows[PROTEIN_COLUMN]:
                ctx.add_rows(rows, f"Ligands of {protein} to table",
                             structures=(LIGAND_COLUMN,))  # fmt: skip
                added += len(rows[PROTEIN_COLUMN])
        ctx.show_message(f"Added {added} ligand row(s)" if added else "No ligands to add")
        if then is not None and added:
            then()

    return ctx.run("Reading ligands…", work, done)


def add_ligand(ctx: ToolContext, protein: str, residue: int) -> int:
    """Add a row for one ligand (a cofactor, say), unless the table has it already."""
    asset = ctx.table.assets.get(protein)
    if asset is None or asset.kind != KIND:
        return 0
    structure = store.as_structure(asset.payload)
    label = structure.residue_label(residue)
    table = ctx.table
    if PROTEIN_COLUMN in table and LIGAND_ID_COLUMN in table:
        for oid in table.oids:
            o = int(oid)
            if (table.cell(PROTEIN_COLUMN, o) == protein
                    and table.cell(LIGAND_ID_COLUMN, o) == label):  # fmt: skip
                ctx.show_message(f"{label} of {protein} is in the table already")
                return 0
    rows = ligand_rows(structure, protein, [residue], ccd=ccd(ctx, offline=True))
    if not rows[PROTEIN_COLUMN]:
        ctx.show_message(f"{label} could not be read as a molecule")
        return 0
    ctx.add_rows(rows, f"Add ligand {label}", structures=(LIGAND_COLUMN,))
    ctx.show_message(f"Added {label} of {protein}")
    return 1


def interactions_to_table(ctx: ToolContext):
    """Count each ligand row's interactions with its protein, by kind, and name the residues
    (the selected rows, or every row): a queued job that keeps the rows finished if cancelled."""
    table = ctx.table
    if PROTEIN_COLUMN not in table:  # no ligand rows yet: asking for this adds them
        return extract_ligands(ctx, then=lambda: interactions_to_table(ctx))
    structures: dict[str, Structure] = {}
    oids, names = [], []
    for oid in ctx.target_oids(selected_only=len(ctx.selection) > 0).tolist():
        protein = str(table.cell(PROTEIN_COLUMN, oid) or "")
        asset = table.assets.get(protein)
        if asset is not None and asset.kind == KIND:
            if protein not in structures:
                structures[protein] = store.as_structure(asset.payload)
            oids.append(int(oid))
            names.append(protein)
    if not oids:
        ctx.show_message("No rows name a protein in this document")
        return None
    labels, smiles, inputs = _row_ligands(ctx, oids)

    def work(job):
        each = zip(names, labels, inputs, smiles, strict=True)
        rows = [profile.Row(structures[p], label, mol_of(given), text)
                for p, label, given, text in each]  # fmt: skip
        return profile.for_rows(rows, progress=job.report, cancel=job.cancel_token)

    def done(found) -> None:
        done_oids = np.array(oids[: len(found)], np.int64)
        results = {
            title: (done_oids, np.array([p.counts[k] if p else np.nan for p in found], float))
            for k, title in enumerate(profile.COUNT_COLUMNS.values())
        }
        results[profile.CONTACTS_COLUMN] = (done_oids, [p.contacts if p else "" for p in found])
        ctx.apply_results(results, label="Interactions to table", overwrite=True)
        ctx.show_message(f"Interactions of {len(found):,} of {len(oids):,} ligand row(s)")

    return ctx.run(f"Finding interactions of {len(oids):,} ligands…", work, done,
                   keep_partial=True)  # fmt: skip


def _row_ligands(ctx: ToolContext, oids: list[int]) -> tuple[list[str], list[str], list]:
    """The rows' ligand IDs, their ligands' SMILES, and the ligands as a job takes them."""
    table = ctx.table
    labels = [""] * len(oids)
    if LIGAND_ID_COLUMN in table:
        labels = [str(table.cell(LIGAND_ID_COLUMN, o) or "") for o in oids]
    if LIGAND_COLUMN not in table or table.kind(LIGAND_COLUMN) != "structure":
        return labels, [""] * len(oids), [None] * len(oids)
    smiles = [str(table.cell(LIGAND_COLUMN, o) or "") for o in oids]
    return labels, smiles, ctx.structure_inputs(LIGAND_COLUMN, oids)
