# Handoff: state of play (updated 2026-11-11, on Windows)

The rules live in `CLAUDE.md`; this file is only what is in flight, so the next session on the
other machine can pick it up. Update it at the end of a session that leaves work open.

## In progress

- **Architecture refactor** (plan approved 2026-11-11; a copy is in Claude's plan file on the
  Windows PC, summarised here). The Protein Viewer becomes a subproject: a Qt-free engine and
  workspace (`chemsuite/protein`, with its own `.cspx` session format), a Qt workbench UI
  (`chemsuite/protein/ui`), and one bridge to the table (`chemsuite/tools/protein`). The app
  around it gets a shared shell.
  - Phase A, done: protein-only documents ask to save; reopened proteins keep their look;
    custom colours saved; `ToolPanel` base and `DocumentPanels`; `JobHost`; layering tests
    (`core` imports no domain package except `core.io` → `sketch`, and `app` never imports a
    tool).
  - Phase B, next: the big mechanical move. `app/{shell,widgets,render3d,document,table}`,
    `protein/{structure,chem,look,workspace,ui}`, `core/io` → `chemsuite/io`; the full suite on
    each move commit.
  - Then Phase C (`core/{history,assets,container}.py`, the protein workspace engine and its
    document Part), Phase D (`.cspx`, the shared start-up, scoped commands, standalone
    `chemsuite --protein`), and Phase E (ToolContext helpers, one 3D stack, tool hooks and
    external workers, splitting long files). After that, viewer Phases 3–5 (presets and views,
    look undo, analysis).
  - Several tables in one window was considered and set aside (2026-11-11).
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
  - Test order: run after the viewer, sketcher and space tests (e.g. `pytest tests/core
    tests/sketch tests/app/test_structure_viewer.py … tests/app/test_sar_ui.py`), the SAR
    tests' result-table selections come up empty (`test_free_wilson`, `test_mmp`,
    `test_sali`). They also fail like this at `aab4e91`; the suite's normal order passes.

## Next

1. Graceful exit (ROADMAP §7.0).
2. Bioactive conformers (§7.0b) or protein Phase 1, Alignment (§7a), each starting in plan mode.

Still owed by the user: the ChemDraw hotkey list.

## Protein code map (until Phase B moves it)

- `chemsuite/protein/` (Qt-free):
  - `model.Structure`: a struct-of-arrays, keyed by atom uids; `store`: its bytes.
  - `io`, `fetch`: reading, writing and downloading files.
  - `bonds`, `bond_orders`, `components`, `classify`: bonds, bond orders (file, CCD cache,
    geometry), residue kinds.
  - `ligands`, `poses`, `interactions`, `profile`: ligands, poses, contacts, the PLIF columns.
  - `select`, `sequence`, `secondary`.
  - `display`: how each protein looks, with background `Build`s (`pending()`/`accept()`);
    `scene`, `cartoon`, `surface`, `mesh`, `colours`, `pick`: what is drawn.
- `chemsuite/app/protein/`:
  - `session.py`: the one object every other part talks to.
  - `panel.py` (a `ToolPanel`), `menubar.py`, `menus.py`, `export.py`, `tree.py`,
    `sequence_view.py`, `selection_bar.py`, `models_bar.py`, `pocket.py`, `labels.py`,
    `overlay.py`.
  - `builds.py`: runs builds as background jobs and drops stale results.
  - `view.py`: on `RhiView`.
- `chemsuite/tools/protein/`: the tool (window, open/fetch/remove, ligands and interactions
  to the table) and its dialogs.
- Tests: `tests/protein/`, `tests/app/test_protein_ui.py`, `test_protein_tool.py`,
  `test_protein_view.py`, and the benchmarks (`pytest -m bench
  tests/protein/test_bench_protein.py`).
