# Handoff: state of play (updated 2026-11-12, on Windows)

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
  - Phase B, done 2026-11-11: `app/{shell,widgets,render3d,document,table}`;
    `protein/{structure,chem,look,ui}` (`ui` was `app/protein`, with its own shaders); the
    bridge in `tools/protein` (`rows.py`, `follow.py` with the workbench's `LinkedRow`,
    `to_table.py`); `core/io` → `chemsuite/io`; tests mirror the code. Moves are done by
    `scripts/move_modules.py`. The launcher moved to `chemsuite.app.shell.main`: after
    pulling, run `pip install -e .` once on each machine (until then `chemsuite` will not
    start).
  - Not moved in Phase B: `protein/workspace` (Phase C creates it), and the plan's shared
    choosers (`tools/*/common` widgets into `app/widgets`), which merge code rather than move
    files (Phase E).
  - Phase C, done 2026-11-12: `core/history.py`, `core/assets.py` (`AssetStore`, `AssetHost`),
    `core/container.py`; the Qt-free `protein/workspace` session with `SessionSignals` in the
    UI; the workspace saved as a document part (`ToolContext.document_part`) with the viewer
    closed, old documents migrated. `Workspace` moved to `app/document/workspace.py`.
  - Next, Phase D (`.cspx`, the shared start-up, scoped commands, standalone
    `chemsuite --protein`), and Phase E (ToolContext helpers, one 3D stack, tool hooks and
    external workers, splitting long files). After that, viewer Phases 3–5 (presets and views,
    look undo, analysis).
  - Several tables in one window was considered and set aside (2026-11-11).
- **Plots** (2026-11-12): colour maps, colour and size ranges and point sizes in every plot's
  options (`app/plots/controls.py` holds the controls, `spec.py` the `PlotSpec`); the
  Constellation Plot (`sar/constellation.py`, the CCR method with a pruned RECAP search;
  `app/plots/constellation.py`, the renderer and the background `CoreFinder`).
- **Protein editing and the Sequence window** (2026-11-11/12, before Phase C, at the user's
  request): a click selects an atom (double-click the residue, Alt-click the chain); Edit ▸
  Delete, Element, Charge, Rename Residue/Chain, Renumber, Undo/Redo, each one table step
  (`protein/structure/edit.py`, `ProteinSession.edit`, the look carried by uid in
  `Display.carry`); a ligand's linked rows are relabelled in the same step
  (`LinkedRow.relabel`, `ToolContext.cell_ops`). The sequence strip moved into its own window
  (the Sequence menu, last), with a row each for ligands, glycans, additives, ions and waters.
  Phase C moves the session's edit code into `protein/workspace` with the rest.
- **CI** was red from RDKit 2026.9.2 (fixed 2026-11-12): cocaine's wedges (the sketcher now
  keeps the layout whose stereo survives nudges), Murcko scaffolds on macOS ("bad any cast":
  `descriptors/scaffold.py` stands in, equal to RDKit's on 26,048 molecules), a protein-tree
  use-after-free on Linux, and job callbacks reaching a closed window.
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

- **Protein preparation** (plan approved 2026-11-11, `docs/PREPARE_PLAN.md`): a coupled
  protein–ligand protonation and hydrogen model in `chemsuite/protein/prepare/`.
  - **P0 done 2026-11-12**: [`docs/prepare/p0-review.md`](docs/prepare/p0-review.md)
    has the licences, the 262-entry neutron set, the scores and the acceptance thresholds for
    P1–P3. The harness is `scripts/review/prepare_truth.py` and `prepare_review.py`, with the
    bench test `tests/protein/test_bench_prepare.py`.
  - The data live under `samples/prepare` (not committed) and are rebuilt by
    `prepare_review.py`'s steps: inputs, reduce (in WSL), propka and gemmi (in a bench venv),
    protoss and score.
  - **MOE** is optional (2026-11-12: its outputs are hard to get, so P1 goes ahead without
    them). If they ever come: Protonate3D on `samples/prepare/inputs/<ID>.pdb` at the pH in
    `samples/prepare/moe/entries.csv`, saved as `samples/prepare/out/moe/<ID>.pdb`, then
    re-run `score` and add the column to the review.
  - **P1 in progress** (plan approved 2026-11-12; Claude's plan file on the Windows PC,
    summarised here). Done: `structure/edit.add_atoms` (new atoms into residues, uids kept);
    `CcdCache.block`; `protein/prepare/templates.py` (CCD ideal geometry, internal
    coordinates and rotors for every hydrogen; 51 residues bundled in `prepare/data`, rebuilt
    by `scripts/build_prepare_templates.py`). Facts found: the CCD draws Lys, Arg and His
    charged and the termini neutral; NME's hydrogen is placed from the residue before it.
  - **`hydrogens.py` done 2026-11-12** (`add_hydrogens` → `Protonated`: structure, `states`
    such as HID/HIE/CYM, `Rotor`s by uid, `Note`s). Templates where every reference atom is
    present; otherwise VSEPR from the real neighbours, other residues included (the amide H
    in the peptide plane). His is HIE unless NE2 binds a metal or another residue; a soft
    metal (Zn, Fe, Cu…) within 2.9 Å takes a proton's place. Ligands: `to_mol` +
    `AddHs(addCoords=True)`, CCD hydrogen names where known. On ubiquitin's NMR hydrogens
    the fixed ones land 0.03 Å from the deposit. `prepare_review.py chemsuite` (~1.5 min)
    gives the no-flip baseline: orientation 50/48/51%, His tautomer 48%, Asp/Glu/Lys/Tyr/Cys
    the textbook numbers, ligand H count 96%. Fixed on the way: `Template.with_hydrogen`
    dropped `charges`; the harness read every 1.34 Å S–H as absent (Cys rows of
    p0-review.md re-scored, 356 Cys now).
  - **Flips, His states and rotors done 2026-11-12** (`sites.py`, `energy.py`, `solve.py`;
    `solve.optimise(protonated, ph)` → structure and `Choice`s with gaps). Weights fitted by
    `prepare_review.py fit` (training half by CRC32 of the ID) and set as the
    `energy.Weights` defaults. All 262: Asn 86.2, Gln 80.4, His orientation 89.6, His
    tautomer 63.9% (Reduce 85/79/88/62); the test half also beats or matches Reduce. Details
    and what is open (His charge, HID read as HIE) in
    [`docs/prepare/p1-review.md`](docs/prepare/p1-review.md); the bench test holds the
    marks. Iterating: build the tables once and pickle them (the fit takes ~3 min
    otherwise); re-run `chemsuite --redo` before `score`.
  - **P1 done 2026-11-12**: `inspect.py` (report rows), `repair.py` (one location; missing
    side chains and OXT from templates), `pipeline.prepare` (repair → waters → hydrogens →
    optimise; report with ⚠ for near-ties). The benchmark runs the pipeline: Asn 86.2, Gln
    80.4, His orientation 89.6, His tautomer 63.6%, ligand H 96.0%, all 262 entries.
  - **Next: P2** (titration: Asp/Glu/Lys/Tyr/Cys/His states at a pH, desolvation and the full
    Coulomb network, PKAD-2, PROPKA parity), then P3 (ligand states, waters, metals), then P4
    (the Prepare dialog and report panel, applied as one undoable edit). Each starts in plan
    mode.
  - **Relax after preparation** (optional minimisation) is on the roadmap (§7a 3) and gets its
    own plan after P1.

## Next

- **Saved for later: SQL import** (File ▸ Import from Database, read-only queries; SQLite
  built in, PostgreSQL/MySQL/Oracle/SQL Server through optional extras). The plan was approved
  2026-11-12 and is in [`docs/DATABASE_PLAN.md`](docs/DATABASE_PLAN.md). No code yet.

1. TMAP and Chemical Space Networks (ROADMAP §7.0d, high priority): plan them together in plan
   mode.
2. Graceful exit (ROADMAP §7.0).
3. Bioactive conformers (§7.0b) or protein Phase 1, Alignment (§7a), each starting in plan mode.

Still owed by the user: the ChemDraw hotkey list.

## Protein code map

- `chemsuite/protein/` (Qt-free outside `ui`):
  - `structure/`: `model.Structure` (a struct-of-arrays keyed by atom uids), `store` (its
    bytes), `io` and `fetch` (files and downloads), `bonds`, `bond_orders`, `components`
    (CCD cache), `classify`, `select`, `sequence`, `secondary`, and `data/`.
  - `chem/`: `ligands` (ligands as molecules; `GUESSED`), `poses`, `interactions`, `profile`
    (the PLIF columns).
  - `look/`: `display` (how each protein looks, with background `Build`s:
    `pending()`/`accept()`); `scene`, `cartoon`, `surface`, `mesh`, `colours`, `pick`: what is
    drawn.
  - `workspace/`: `session.py` (`WorkbenchSession`, Qt-free, plain `events.Event`s, reads an
    `AssetHost`) is the one object every other part talks to; `workspace.py`
    (`ProteinWorkspace`) is the document part (`protein_workspace`) that keeps it and the view.
  - `ui/`: the Qt workbench. `session.py` is `SessionSignals` (the session's events as Qt
    signals, one per window);
    `host.py` holds `LinkedRow`; `panel.py` (a `ToolPanel`), `menubar.py`, `menus.py`,
    `export.py`, `tree.py`, `sequence_view.py`, `selection_bar.py`, `models_bar.py`,
    `pocket.py`, `labels.py`, `overlay.py`; `builds.py` runs builds as background jobs and
    drops stale results; `view.py` is on `RhiView`; `shaders/` holds the `pro_*` shaders.
- `chemsuite/tools/protein/`, the bridge: `tool.py` (the tool, its window, open, fetch,
  remove), `rows.py` (the ligand-row columns, `ligand_rows`), `follow.py` (`TableRow`, the
  current row followed in the workbench), `to_table.py` (Ligands, Add Ligand and
  Interactions to Table), `dialogs.py`.
- Tests: `tests/protein/{structure,chem,look,ui}`, `tests/tools/protein/test_protein_tool.py`,
  and the benchmarks (`pytest -m bench tests/protein/test_bench_protein.py`).
