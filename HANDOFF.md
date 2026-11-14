# Handoff: state of play (updated 2026-11-14, on Windows)

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
  - **P2 done 2026-11-12** ([`docs/prepare/p2-review.md`](docs/prepare/p2-review.md)):
    titration forms in `sites` (with `titrate=True`), pH-linear tables and new terms in
    `energy`, `populations.py` (mean field, Monte Carlo as the check), `pipeline`
    `Options.titrate`, `Prepared.pka` and `population`. On PKAD-R's test half the RMSE is 1.32
    against PROPKA's 1.33; neither beats the null model (1.15). With titration on, the neutron
    His and Asp/Glu marks drop, so it is **off by default**. P1 reproduces byte for byte with it
    off. Data: `samples/prepare/pkad` (PKAD-R CSV saved by hand; `pkad_review.py fetch`,
    `propka`, `tables`, `fit`, `chemsuite`, `score`).
  - **Next:** P3 (ligand states, waters, metals), then P4 (the Prepare dialog and report
    panel, applied as one undoable edit), each in plan mode. A second titration round
    (PROPKA-style H-bond determinants, His titrated for charge only) is listed in
    p2-review.md §4.
  - **Relax after preparation** (optional minimisation) is on the roadmap (§7a 3) and gets its
    own plan after P1.

- **chemsuite pKa, U2 training** (`docs/UNIPKA_PLAN.md`; Windows PC only, WSL
  `~/bench/unipka`). The PC hard-reset (Kernel-Power 41) on 2026-11-12 with pKaHub
  preprocessing and the reference fold training at once; `resume.sh` now runs one heavy job
  at a time (log `logs/resume.log`): pKaHub lmdb by `scripts/unipka/prep_lmdb.py` (upstream's
  conformers with 180 s per row: rifampicin-like macrocycles hung upstream's script; dropped
  rows in `train.lmdb.dropped.tsv`, 22 of 22,132) → reference fold 0 on dwar-iBond, resumed
  from epoch 6 (`reference/`; folds 1–4 skipped unless fold 0 is off from the paper) →
  the 5 shipped pKaHub folds (`shipped/`, `shipped10.sh`, log `logs/shipped10.log`) at **10
  epochs** (a pKaHub epoch is 553 updates, 2.7 × dwar's; 20 epochs would be ~25–30 h).
  Reference fold 0 done (valid RMSE 0.94), converted to `samples/unipka/weights/reference_f0.npz`
  (`convert_weights.py` now stores each fold's training-split mean, `fold_means`, and
  `thermo.common_mean` puts the folds on their average). `evaluate.py` runs in WSL in
  `~/bench/unipka/evalenv` (the app's RDKit 2026.09.2, `PYTHONPATH` = the repo; on Windows it
  ran out of memory beside training). Paper's Table 1 RMSE (5 folds, after a ChEMBL pretraining
  stage we do not have): Novartis acid 1.061, base 0.653, SAMPL6 0.716, SAMPL7 0.735, SAMPL8
  0.878; MolGpKa there 1.287 / 1.064 / 0.773 / 0.98 / 1.150.
  **Done 2026-11-14:** the 5 shipped folds (valid RMSE mean 0.753, every best epoch the 10th)
  converted to `samples/unipka/weights/shipped.npz` (440 MB, not committed). External tests
  (`samples/unipka/shipped.eval.txt`): 10 conformers, 5 folds RMSE 0.83 overall, MolGpKa 1.02;
  3 conformers 0.84, 2D only 0.98. The training report and the novice's guide are Claude Docs,
  exported to `samples/*.docx`. Next: `evaluate.py --enumerated` with the cache, pick the
  app's conformer/fold count (3 conformers looks right), then U4 (protein ligand states).
- **Decomposition and enumeration** (done 2026-11-14, ROADMAP §6): Chemistry ▸ Decompose
  (R-Group Decomposition, Fragments, Central Ring), Chemistry ▸ Enumerate (Reaction
  Products with 28 named templates in `analysis/reaction_templates.py`, R-Groups on a Core,
  Recombine Fragments); Standardize keeps stripped salts in a `Salts` column and can form
  salts (`analysis/salts.py`). The Structure Viewer is now Data ▸ Browser (bottom section,
  Ctrl+Shift+B); Ctrl+Shift+3 is gone. Not yet merged to `main`: run the full suite first.

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
