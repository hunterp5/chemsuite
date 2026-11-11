# Handoff: state of play (updated 2026-11-11, on Windows)

The rules live in `CLAUDE.md`; this file is only what is in flight, so the next session on the
other machine can pick it up. Update it at the end of a session that leaves work open.

## In progress

- **By-hand check of Protein Viewer Phase 1 on Windows** (ROADMAP §7.0c).
  Done 2026-11-11 by script:
  - 5GNK: 80U's bond orders come from the file's own component table, and from the CCD when
    the file has none (the acrylamide's C=C intact, no ⚠); geometry alone loses the C=C, so
    the ⚠ there is right. 1ATP's ATP from the file. `samples/5GNK.cif` is a local copy.
  - Direct3D and WARP draw the same, light and dark. Fixed: ghosted (focus) and selected
    cartoons drew as yellow speckle on Direct3D, which interpolates the flag 4.0 as 3.9999.
  - Fixed (found by hand): closing the viewer window left its deleted session listening to
    the table, so the next Remove, Fetch or row change raised and the table view missed it.

  Still by hand:
  - Play 1D3Z's NMR models.
  - Focus Ligand: check the pocket skin and the fog.
  - Label residues. Step through ligand rows. Move a pose in the sketcher and confirm it
    overlays in the pocket.
  - Run Interactions to Table and cancel it part-way: the finished rows must stay.
  - Remove a protein, then undo.
  - Copy Shown Rows to a new document and confirm the protein travels with them.
  - Save, reopen, and confirm everything comes back.
  - WARP over RDP or with the GPU disabled (forced WARP is checked; the automatic fallback is
    not).
- **Known issues:**
  - When a pocket is focused, the faded ribbons in front of it can cross the ligand.
  - A kept row whose protein was removed is not yet marked in the table.
  - Save Protein As (mmCIF) drops the file's `_chem_comp_*` tables, so a ligand with no CCD
    entry reopened offline falls back to a geometry guess. Document saves keep them.

## Next

1. Graceful exit (ROADMAP §7.0).
2. Bioactive conformers (§7.0b) or protein Phase 1, Alignment (§7a), each starting in plan mode.

Still owed by the user: the ChemDraw hotkey list.

## Protein code map

- `chemsuite/protein/` (Qt-free):
  - `model.Structure`: a struct-of-arrays, keyed by atom uids.
  - `io`: reading and writing files.
  - `bonds`, `components`: bond orders from the file, the CCD cache or templates.
  - `ligands`: ligand rows.
  - `display`: background `Build`s with `pending()`/`accept()`.
  - `scenes`, `cartoons`, `surfaces`, `mesh`: what is drawn.
  - `contacts`, `profile`: interactions and the PLIF columns.
  - `poses`.
- `chemsuite/app/protein/`:
  - `session.py`: the one object that every other part talks to.
  - `panel.py`, `tree.py`, `sequence_view.py`, `toolbar.py`, `pocket.py`, `labels.py`,
    `models_bar.py`.
  - `builds.py`: runs builds as background jobs and drops stale results.
  - `view.py`: on `RhiView`.
- `chemsuite/tools/protein/`: the tool, its dialogs and the menubar.
- Tests: `tests/protein/`, `tests/app/test_protein_ui.py`, `test_protein_view.py`, and the
  benchmarks (`pytest -m bench tests/protein/test_bench_protein.py`).
