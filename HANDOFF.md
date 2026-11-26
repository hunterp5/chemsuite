# Handoff: state of play (updated 2026-11-25, on Windows)

The rules live in `CLAUDE.md`; this file is only what is in flight, so the next session on the
other machine can pick it up. Update it at the end of a session that leaves work open.

## In progress

- **Permeability** (2026-11-25, Windows, committed to `dev` as `647054e`; ROADMAP §7.0d, plan
  approved in-session). **Done and shipping:** two models at Chemistry ▸ Calculate ▸
  Permeability…, MDR1-MDCK ER (0.360 MAE, R² 0.512, SVR) and Caco-2 (0.354 MAE, R² 0.655,
  random forest), both out of fold on scaffold splits, fitted by `chemsuite/qsar` and applied
  in numpy at ~13 ms a row. `chemsuite/permeability/data/*.csqm` (12 MB total) is pinned by
  size and SHA-256 in `models.PINS`, with Biogen's MIT licence and a NOTICE beside it, and
  `pyproject.toml` has the `package-data` entry they need to reach a wheel. Full suite green.
  Everything measured is in `docs/PERMEABILITY_PLAN.md`.

  **Tabled deliberately — four open items, none blocking:**
  1. **The `pka` arm has no accuracy number.** It was abandoned on measured cost, about 60
     CPU-seconds a molecule on two independent sets (MDCK ~67,000 CPU-s in 2 h, Caco-2
     ~50,000 in 1 h, neither finishing its pool), against ~13 ms a row for what ships. The
     same pass would run at predict time, so cost alone closes it as a *default* — but the
     user's call during planning was that the bake-off decides on accuracy, so this reverses
     that on different grounds and the accuracy question is genuinely still open. To settle
     it: `train.py --arm pka --endpoint caco2` (the smaller endpoint) overnight. If it clears
     +0.05 log MAE the answer is a second, slower, optional model, which is new scope.
  2. **`train.py` does not pass a `progress` callback into `features.build`**, so pool
     progress is invisible and the two runs above had to be killed blind. Fix this first, and
     run an expensive arm on the smallest endpoint first.
  3. **`certifi` is in `.venv` but not in `pyproject.toml`.** Python's HTTPS to
     `dataverse.harvard.edu` fails on this PC without it, so `scripts/permeability/fetch.py`
     needs it on a fresh venv. Whether it belongs in the dependencies is undecided.
  4. **Caco-2 is data-limited**, not model-limited: three feature pools and three different
     winning estimators give the same score. The lever is the Esaki 2024 ChEMBL curation
     (+1,661 entries, with its exclusions applied), not more tuning.

  **PAMPA was fitted and held back** (ROC-AUC 0.740 but balanced accuracy 0.570 on data that
  is 85% permeable). It is still curated and trainable, and restoring it is one entry in
  `ENDPOINTS` plus a pin; the classification column path stays covered by a test.

  `samples/` is gitignored, so the pinned snapshot and curated CSVs are local to this PC —
  `scripts/permeability/fetch.py` reproduces them by digest on the Mac.

- **Scaffold Tree** (2026-11-25, Windows, merged to `main`; ROADMAP §7.0h item 2, plan
  approved in-session). New: `sar/scaffold_tree.py` (the engine; rules after ScaffoldGraph,
  rule 7 chemsuite's own), `tools/sar/scaffold_tree_ui.py` (Data ▸ SAR ▸ Scaffold Tree…, a
  `QTreeView` panel), `TreeTableModel` in `app/widgets/result_table.py` (a tree model base:
  tools cannot call `index.column()` under the `.column(` ratchet), tests
  `tests/sar/test_scaffold_tree.py`, `tests/tools/sar/test_scaffold_tree_ui.py`. Targeted
  tests pass; the full suite passed before the merge (with GTM).
- **GTM, map + landscapes** (2026-11-25, Windows, merged to `main`; ROADMAP §7.0h item 1, plan
  approved in-session). New: `analysis/gtm.py`, `tools/space/gtm_ui.py` (Data ▸
  Dimensionality Reduction ▸ Generative Topographic Map…), `app/plots/gtm.py` (the GTM
  Landscape plot kind), tests `tests/analysis/test_gtm.py`, `tests/tools/space/test_gtm_ui.py`;
  edits to `app/plots/{spec,controls,panel}.py`, `tools/{space,plot}/tool.py`, the guide and
  ROADMAP. The full suite passed before the merge (every plot's saved
  state gained `gtm_map`/`gtm_value`, and documents gained the `gtm_map` asset kind). Found,
  not fixed:
  `np.linalg.eigh` of 2048 × 2048 takes ~30 s on this PC, which is most of a fingerprint
  SOM/PCA/GTM's time (PCA in `embedding.pca_axes`).

- **Memory and CPU deep dive** (2026-11-24, Windows, 6 physical cores / 12 threads / 32 GB;
  uncommitted). Symptom on the work PC (16 threads, 16 GB): pKa and MMP slow the whole
  machine, and Python's memory (Task Manager's default column, the private working set) stays
  high while idle. Harness in `scripts/review/` (`memlib.py`, `pka_pool_mem.py`,
  `mmp_mem.py`; `.venv` python, args `SDF N [workers]`). Done:
  - **pKa model shared** (`ionisation/unipka/mapped.py`, `core/folders.py`): the weights are
    written once as float32 to `<cache>/derived/` and memory-mapped read-only by every
    worker. 11 pKa workers: ~3 GB private → 0.57 GB (52 MB each) plus one shared 182 MB;
    predictions bit-identical, same speed.
  - **Idle workers end** (`JobRunner.end_idle_workers`): the shared pool after 60 s with no
    job running or queued, the wide worker after 300 s; the next job starts new ones.
  - **Worker count scales with the machine** (`core/jobs/machine.py`, `default_workers`):
    physical cores − 1, and no more than 25% of RAM (or 50% of free RAM) at 250 MB per
    worker. Cost measured here: 5 workers instead of 11, MMP 50k rows 66 s → 118 s.
  - Picture cache per result table capped at 32 MB as well as 400 pictures; `Pair` uses
    slots and MMP keeps one copy of each fragment SMILES (50k rows peak 371 → 317 MB).
  - Not done: the app process's BLAS threads (they reserve ~30 MB commit each, not in the
    Memory column; now follow `default_workers`, so fewer).
- **Table responsiveness at 100k rows** (2026-11-24, uncommitted; bench
  `scripts/review/table_ops_bench.py samples/bindingdb_100k.csd [--profile STEP]`, a real
  offscreen window, each step timed until idle). Fixed in `app/table/view_selection.py`:
  - Replacing a scattered selection used `ClearAndSelect`, which compares old and new ranges
    pairwise: every other row 11.7 s, inverting it 22 s → 0.4 / 0.7 s (clear, then Select).
  - Select All mirrored every column into the header's selection model, so each header
    repaint (any hover) asked about every row: 2.7 s per repaint → instant. The mirror is now
    filled only during a Ctrl+press on a column header (`TableView.header_press`).
  - `selected_column_names` (header menu) read from the ranges, not `isColumnSelected`.
  - Left: with ~50k scattered rows selected, a page of scrolling ~110 ms and a Ctrl/Shift+click
    ~0.5 s (Qt scans every range per painted cell / merge). Removing it means drawing row
    selections from the store instead of Qt's selection model: a bigger change to plan.
  - Found, not fixed: Edit ▸ Select All from the menu after a column-header click stays in
    column mode, so the store gets no rows (Delete Selection then deletes nothing).

- **Macromolecule Viewer** (2026-11-16, plan approved; Claude's plan file on the Windows PC,
  `whimsical-growing-ripple.md`). The Protein Viewer is renamed in the code and on screen:
  `chemsuite/protein` → `chemsuite/macromolecule`, `tools/protein` → `tools/macromolecule`
  (tests follow). **After pulling, run `pip install -e .` once, and delete any leftover
  `chemsuite/protein`, `tests/protein` and `tests/tools/protein` folders** (only
  `__pycache__` is left in them, and Python can import them as namespace packages). Saved
  names are unchanged (asset kind `protein`, the `protein_workspace` part, tool ids), so every
  document opens as before.
  - Done: session files (`.cspx`, `macromolecule/workspace/file.py`; Open, Save Session,
    Replace or Add with ids renamed); File ▸ Save (every shown structure in one file);
    Clear Session; the list's Rename / Duplicate / Delete; the Ligand and Window menus folded
    into View and Tools; Sequence pinned at the right (a corner `QMenuBar`).
  - Deferred: opening a `.cspx` from the table window's File ▸ Open (the app cannot call a
    tool until the tool hooks of Phase E); the standalone viewer (Phase D.2–3).
  - Undo of Clear Session brings the structures back with their looks, not the focus, saved
    selections, labels or measurements (the guide says so).
  - Found: a style sheet on a menu-bar corner `QToolButton` keeps the palette it was polished
    with, so the table window's **Layout** and **Log** buttons draw white text after a switch
    to the Fusion light theme at run time (Preferences). Not fixed yet.
  - Hierarchy deep dive (2026-11-16, against PyMOL 2 and ChimeraX 1.12 on the Windows PC):
    a chain is everything with that author chain id (polymer, ligands, ions, waters) in all
    three; the list's chain entries are polymer-only, as ChimeraX's chains; NMR models are
    PyMOL's states; sequence gaps come from the entity sequence (PyMOL adds coordinate-less
    CA atoms instead). Fixed: chain names over four characters merged (`AA-10`/`AA-11` in
    assembly files); a blank chain could not be selected (`chain ''` now, PyMOL's spelling);
    two residues modelled at one place (1EJG PRO/SER 22) were both shown and both kept by
    Prepare (now one group, one letter). Open, not fixed (the user's call):
    - Segment ids: done (2026-11-16). `Structure.chain_segment`, one per chain; a chain whose
      residues change segment is cut into one chain per run (`io._split_segments`, every
      model, before the entities); `segi`/`segment` select; labels `PROB/THR 1`; PDB writes
      them back, mmCIF does not (no standard field). Still keyed by chain name alone, as
      before: `ligands.py` ligand ids and `simulate/build.py` residues.
    - A model whose atoms differ from the first is dropped with only a log line (PyMOL and
      ChimeraX keep it).
    - Kinds differ by design: a glycan-only chain (2RH1 B) is listed under Glycans, not as a
      chain; PO4/SO4 are additives (PyMOL inorganic, ChimeraX ions or ligand); lipids
      (CLR, PLM) are ligands.
    - The Sequence Viewer shows the first residue at a shared place even when the other is
      picked (ChimeraX does the same).
- **Architecture refactor** (plan approved 2026-11-11; a copy is in Claude's plan file on the
  Windows PC, summarised here). The Protein Viewer becomes a subproject: a Qt-free engine and
  workspace (`chemsuite/macromolecule`, with its own `.cspx` session format), a Qt workbench UI
  (`chemsuite/macromolecule/ui`), and one bridge to the table (`chemsuite/tools/macromolecule`). The app
  around it gets a shared shell.
  - Phase A, done: protein-only documents ask to save; reopened proteins keep their look;
    custom colours saved; `ToolPanel` base and `DocumentPanels`; `JobHost`; layering tests
    (`core` imports no domain package except `core.io` → `sketch`, and `app` never imports a
    tool).
  - Phase B, done 2026-11-11: `app/{shell,widgets,render3d,document,table}`;
    `protein/{structure,chem,look,ui}` (`ui` was `app/protein`, with its own shaders); the
    bridge in `tools/macromolecule` (`rows.py`, `follow.py` with the workbench's `LinkedRow`,
    `to_table.py`); `core/io` → `chemsuite/io`; tests mirror the code. Moves are done by
    `scripts/move_modules.py`. The launcher moved to `chemsuite.app.shell.main`: after
    pulling, run `pip install -e .` once on each machine (until then `chemsuite` will not
    start).
  - Not moved in Phase B: `macromolecule/workspace` (Phase C creates it), and the plan's shared
    choosers (`tools/*/common` widgets into `app/widgets`), which merge code rather than move
    files (Phase E).
  - Phase C, done 2026-11-12: `core/history.py`, `core/assets.py` (`AssetStore`, `AssetHost`),
    `core/container.py`; the Qt-free `macromolecule/workspace` session with `SessionSignals` in the
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
- **Plots, 2026-11-14:** a selection made in a plot scrolls the table to its first new row
  (`app/table/snap.py`). **TMAP and Chemical Space Networks** (ROADMAP §7.0d, plan approved
  2026-11-14) are done:
  - The Qt-free `analysis/graphs` package holds:
    - `knn.py`: exact neighbours up to 10k rows; above that, MinHash over the Morgan bits, a
      prefix forest and one neighbour-of-neighbour pass, reranked by exact Tanimoto;
    - `trees.py`;
    - `layout.py`: a seeded multilevel force layout (Walshaw), with packing;
    - `tmap.py` and `csn.py`.
  - `app/plots/finder.py` is the shared background `Finder`, and `CoreFinder` now builds on it.
  - `app/plots/graph.py` holds `GraphFinder` and `GraphRenderer`; `GraphRenderer` is a
    `ScatterRenderer` whose points come from the layout.
  - Bench: `scripts/review/chemspace_bench.py` (`smiles`, `knn`, `tmap`, `csn`; it caches under
    `samples/chemspace`).
    - Measured on 49k BindingDB compounds: recall@10 0.89; TMAP in 15 s.
    - A 5k-compound network at 0.65: 1.6 s.
  - Not done from the plan:
    - `tmap-silicon` as a side-by-side reference (no bench venv yet);
    - TMAP above 49k rows (no larger set in `samples/`);
    - the redraw benchmark at 100k points.
  - Open: a few branch crossings in big TMAPs (more fine force steps may help); MHFP6 as a
    fingerprint choice.
- **Protein editing and the Sequence window** (2026-11-11/12, before Phase C, at the user's
  request): a click selects an atom (double-click the residue, Alt-click the chain); Edit ▸
  Delete, Element, Charge, Rename Residue/Chain, Renumber, Undo/Redo, each one table step
  (`macromolecule/structure/edit.py`, `ProteinSession.edit`, the look carried by uid in
  `Display.carry`); a ligand's linked rows are relabelled in the same step
  (`LinkedRow.relabel`, `ToolContext.cell_ops`). The sequence strip moved into its own window
  (the Sequence menu, last), with a row each for ligands, glycans, additives, ions and waters.
  Phase C moves the session's edit code into `macromolecule/workspace` with the rest.
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
  protein–ligand protonation and hydrogen model in `chemsuite/macromolecule/prepare/`.
  - **P0 done 2026-11-12**: [`docs/prepare/p0-review.md`](docs/prepare/p0-review.md)
    has the licences, the 262-entry neutron set, the scores and the acceptance thresholds for
    P1–P3. The harness is `scripts/review/prepare_truth.py` and `prepare_review.py`, with the
    bench test `tests/macromolecule/test_bench_prepare.py`.
  - The data live under `samples/prepare` (not committed) and are rebuilt by
    `prepare_review.py`'s steps: inputs, reduce (in WSL), propka and gemmi (in a bench venv),
    protoss and score.
  - **MOE** is optional (2026-11-12: its outputs are hard to get, so P1 goes ahead without
    them). If they ever come: Protonate3D on `samples/prepare/inputs/<ID>.pdb` at the pH in
    `samples/prepare/moe/entries.csv`, saved as `samples/prepare/out/moe/<ID>.pdb`, then
    re-run `score` and add the column to the review.
  - **P1 in progress** (plan approved 2026-11-12; Claude's plan file on the Windows PC,
    summarised here). Done: `structure/edit.add_atoms` (new atoms into residues, uids kept);
    `CcdCache.block`; `macromolecule/prepare/templates.py` (CCD ideal geometry, internal
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
  - **P3 done 2026-11-13**: ligand protomers and tautomers (`ligand_sites.py`), pocket waters
    (`waters.py`), metal coordination; the named cases in `tests/macromolecule/prepare/test_named_cases.py`.
  - **P4 done 2026-11-14** (the docking plan's Stage A): **Protein ▸ Prepare…**
    (`tools/macromolecule/prepare_dialog.py`) — protein, pH, waters, repair, optimise, and Titrate
    (Off / Pocket / All) folded away — run as a background job and applied with
    `session.edit` as one undo step; the report panel (`macromolecule/ui/report.py`), whose lines
    select their residue, kept with the document (`macromolecule/workspace/reports.py`, saved under
    the protein workspace's `reports`; older documents load unchanged). On the way: `titrate=POCKET`
    now actually titrates the pocket (`pipeline._pocket`), and rotor rows are reported only when
    the orientation makes a hydrogen bond *and* a different one was nearly as good
    (`pipeline.BONDED`, `DEGENERATE`), the free ones counted in one line as the waters are —
    crambin's report went from 14 ⚠ of 17 rows to 7 of 11.
  - **Stage B of the docking plan, 2026-11-14 (committed and merged, not yet tried by hand):** the
    bundled Vina 1.2.7 binaries with their licence, NOTICE and pinned hashes
    (`macromolecule/dock/bin`, `scripts/fetch_vina.py`); the Qt-free engine `macromolecule/dock/`
    (`box.py`, `pdbqt.py` the Meeko boundary, `vina.py` the binary, `run.py` many ligands with
    progress, cancel and finished ligands kept); the bridge `tools/macromolecule/dock.py` and
    `dock_dialog.py` (Protein ▸ Dock… and Rescore Poses…, also in the viewer's Tools menu);
    `meeko>=0.8` a core dependency; package data and `.gitattributes` for the binaries.
    Tests: `tests/macromolecule/dock/` (27, one real dock at exhaustiveness 1 ≈ 11 s, the redocking
    accuracy check behind `-m bench` ≈ 75 s) and `tests/tools/macromolecule/test_dock_ui.py` (7).
    Found on the way: **Meeko downloads a CCD component for any residue name it does not
    know**, so `pdbqt.receptor` keeps only `known_residues()` (polymer, waters, ions, glycans)
    and names the rest in `Receptor.dropped`; `--score_only` still needs a box, so
    `score_only` passes `--autobox` when given none. Redocking 80U into 5GNK puts the top pose
    1.13 Å from the crystal one.
    **Still to do:** the by-hand run in the app (plan's step 4).
  - **Results in the viewer, Inspect, hydrogen-aware H-bonds (2026-11-15, on `dev`):**
    - The viewer's window holds docks (`app/document/panel_window.py`: `PanelWindow`,
      `WindowPanels`; `ToolContext.add_panel(..., beside=)`). The Prepare report, the Inspect
      panel and the Poses panel dock there, falling back to the main window.
    - Poses stay in the viewer (`macromolecule/workspace/poses.py`, `macromolecule/ui/poses.py`): every
      run kept and numbered with its scoring function, a rescore a run of its own, Remove Run,
      Copy to Table on request. Saved as a list per protein (an older single run reads as run 1).
    - **Protein ▸ Inspect…** (`macromolecule/prepare/review.py`, `macromolecule/ui/inspector.py`,
      `tools/macromolecule/inspect_dialog.py`): `pipeline.Decisions` (repairs, waters, alternate
      locations, pinned states; residues named by label so edits elsewhere cannot move them),
      Re-survey, Apply as one edit, a report that records what was chosen by hand, the review
      saved with the document. One ⚠ rule, `pipeline.warns`, for the report and the review.
      ASH/GLH/LYN/CYM/TYM can be chosen only with titration on (documented; His can always
      be HIP).
    - H-bonds use the hydrogen where there is one (H···A ≤ 2.5 Å, D–H···A ≥ 100°), except a
      hydroxyl/thiol/ammonium whose H points elsewhere, judged on heavy atoms.
    - Ligands draw as purple sticks; ions amber; metal-contact lines light gold.
    - Not done: the PanelWindow's dock arrangement is not saved; nothing tried by hand yet.
  - **Done (protein): surface options**, 2026-11-17, [`docs/SURFACE_PLAN.md`](docs/SURFACE_PLAN.md)
    (commits fc40672 to the phase 5 docs). `look/surface.py` builds SES/SAS/vdW/Gaussian from
    `SurfaceOptions`; a `Display` holds the options, style, colors, transparencies and carve
    points (saved; older documents open as Gaussian; `carry` keeps them); `look/electrostatics.py`
    and the deferred `potential` part color By Electrostatics; `ui/surface_menu.py` fills
    Display ▸ Surface and View ▸ Ligand ▸ Pocket Surface ▸. Tried by hand on 1ATP (SES,
    Gaussian, Mesh, Dots, electrostatics, carving): all draw. SES of about 6,000 heavy atoms
    takes about 3 s (the target was 2 s; it is in the background) and its potential 5.5 s.
    Open: a whole-protein Mesh at Normal quality is dense (Coarse thins it); ligands without
    hydrogens take their template's charges, which can be neutral (ATP).
  - **Next (protein):** try Prepare, Inspect and docking by hand in the app.
  - **Titration, second round, done 2026-11-19** ([`docs/prepare/p5-review.md`](docs/prepare/p5-review.md)):
    - Judged on states with `scripts/review/prepare_states.py` (`build --jobs 6`, about
      2.5 min; `score --fast`, seconds). Its pickled cases live in
      `samples/prepare/state_cases.pkl`.
    - Titration still loses to titration off (test half 0.884 against 0.893), so it stays off.
    - Next: run Poisson–Boltzmann and ML pKa tools through that benchmark (licences checked
      first) before building either. Constant-pH MD is the MD-side option.
  - **Research for the next round, 2026-11-19** ([`docs/prepare/p6-research.md`](docs/prepare/p6-research.md)):
    licenses checked; KaML-CBtree and pKAI+ (MIT) first, MCCE4 + NextGenPB next, PypKa only if
    DelPhi's registration suits the user; pKa-ANI and KaML-ESM excluded. Proposed first step:
    a per-site-pKa hook in `prepare_states.py`, with PROPKA scored through it as the bar.
    Done the same day (p6 §5): `score --tools propka`; PROPKA scores 0.864 on the test half
    against P1's 0.893 (cases rebuilt with site keys; the old pickle is beside it as
    `state_cases.before-keys.pkl`).
    **External pKa tools done 2026-11-24 (p6 §6–8): none beats P1, so titration stays off.**
    `scripts/review/pka_tools.py` runs pKAI+ and KaML-CBtree in the WSL conda env `pka-bench`.
    Test-half objectives: P1 0.893, pKAI+ 0.874, PROPKA 0.866 (its CSVs had dropped about 5% of
    sites; fixed), KaML 0.779, or 0.851 with margins. Leaving out KaML's training proteins changes
    nothing. KaML's published `bbh_dis` has a bug (`x1 - y1`), which the runner corrects.
    Every tool loses on His. The bench tools are `prepare_states.py coverage`, `score --skip`
    and `margin --jobs N`; score with `--fast`, since the full run takes 3.5 h.
    MCCE4 + NextGenPB (step 4) is not worth starting unless the user wants PB checked anyway.
  - **Loop modeler** ([`docs/LOOPS_PLAN.md`](docs/LOOPS_PLAN.md); ROADMAP §7a.3): gaps of 2 to
    N residues closed by CCD in the Qt-free `prepare/loops.py`, relaxed by Minimize. **L0 and
    L1 done 2026-11-25:** `loops.py`, `prepare/data/ramachandran.json`
    (`scripts/build_ramachandran.py`, train half only), `tests/macromolecule/prepare/test_loops.py`
    (8), `scripts/review/loop_bench.py` (+ `loop_bench_fixer.py`, PDBFixer in the OpenMM env).
    Results in [`docs/loops/l0-review.md`](docs/loops/l0-review.md): on the 285 test cuts, ours
    + Minimize has median RMSD 0.22–2.4 Å at lengths 2–8 against PDBFixer's 1.7–4.8 Å, and it
    passes every L3 gate up to 8 (length 10 fails on one clash), so `loop_max` defaults to 8.
    **L2 done 2026-11-25, merged to `main`:**
    - `pipeline.Options.loop_max` (More Options ▸ Longest Loop, 1 = Off) and
      `Decisions.capped`;
    - `Prepared.gaps_built`, and a ⚠ "Loop of n residues built" line per loop;
    - a tick box per built gap in Inspect (`review.Issue.buildable`);
    - a per-gap memo in `loops.py`.

    Next: L3's run by hand in the app (Prepare + Minimize on an 8-residue gap). Open:
    length-1 gaps (not benchmarked); a Minimize failure on 1I1J's partial N-terminal Met; the
    menu item still reads "Minimise…" (British).
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
  exported to `samples/*.docx`.
  **2026-11-14, later:** enumerated evaluation done (results table in `docs/UNIPKA_PLAN.md`
  U2): the full Uni-pKa site templates are now bundled (failures 24 → 2 of 356; RMSE 0.82 with
  3 conformers × 5 folds, 0.84 with 1 fold; Uni-pKa's paper 0.83). pKaHub is aqueous (its
  paper), now a CLAUDE.md rule.
  **2026-11-15:** fold 0 at 20 epochs (`e20_fold0.sh`, 2.1 h) beat the shipped fold 0: 0.81
  against 0.84 (enumerated, 3 conformers + 2D), as good as the five 10-epoch folds (results in
  `docs/UNIPKA_PLAN.md` U2; archives `samples/unipka/weights/e20_f0.npz`, `s10_f0.npz`).
  Folds 1–4 at 20 epochs are training (`~/bench/unipka/e20_rest.sh`, log `logs/e20_rest.log`,
  about 8.5 h; nothing else heavy meanwhile). `~/bench/unipka/e20_after.sh` waits for it, then
  converts all five folds to `samples/unipka/weights/e20.npz` and writes
  `samples/unipka/e20.eval.txt` and `e20.enum.eval.txt` (log `logs/e20_after.log`, ends
  `E20_AFTER_DONE`). Both run detached from Claude Code. Then: pick the app's conformer/fold
  count, U4. The PC hard-reset again under training (Kernel-Power 41, bugcheck 0) on
  2026-11-15 at 19:54 in fold 1's epoch 19; both jobs restarted 2026-11-16 08:35 (rerunning
  `e20_rest.sh` resumes a fold from its `checkpoint_last.pt`).
  **Done 2026-11-17:** folds 1–4 at 20 epochs finished (11:22) and `e20_after.sh` (14:38):
  `samples/unipka/weights/e20.npz`; enumerated, 3 conformers + 2D: 5 folds 0.80, 1 fold 0.81
  (10 epochs: 0.82 / 0.84; table in `docs/UNIPKA_PLAN.md` U2). U2 is done.
  **U4 built 2026-11-17** (details and benchmark in `docs/UNIPKA_PLAN.md` U4):
  - Prepare and Inspect take ligand states from chemsuite pKa: fold 0 of e20, 3 conformers.
  - The weights are published as GitHub release `pka-weights-1` (GPL-3.0-or-later, with
    Uni-Mol's MIT license). Since U5 the same file ships in the package.
  - Ligand-H benchmark: tied with MolGpKa on the test half (94.7%).
  **U5 built 2026-11-18** (details in `docs/UNIPKA_PLAN.md` U5): Predict pKa, Standardize's
  Protonate and the pKa descriptors run on chemsuite pKa (`ionisation/profile.py`); MolGpKa is
  removed. The weights ship in the package (`ionisation/unipka/data/chemsuite-pka-1.npz`,
  88 MB, committed to git): no download, prompt or cache. About 5 s per drug-sized molecule
  (3 conformers + 2D).
- **CONFORGE conformers** (§7.0b, plan approved 2026-11-14): the default Generator in Generate
  Conformers since 2026-11-15 (`conformers/conforge.py`); `cdpkit` is a base dependency with
  platform markers (no wheel on Windows 3.11/3.13 or Intel Macs: ETKDG there).
  Benchmark: `scripts/review/conformer_bench.py` on Wang et al.'s 3,354 crystal ligands
  (`samples/conformers/wang`, fetched; ids from their SI xlsx). CONFORGE's two variants are done
  (3,354 ligands each). The ETKDG run crashed on 2026-11-14 after 339 `etkdg-raw` ligands:
  RDKit 2026.09 answers a failed random-coordinates embedding with id `-1` and no conformer
  (platensimycin, 3i8p_840), now dropped in `etkdg._embed` ("could not embed in 3D"). The
  four ETKDG variants were restarted after the fix (they resume; log
  `samples/conformers/wang/out/bench.log`); scored in `docs/conformers/bench.md` (2026-11-15).
  CONFORGE's macrocycle sampling is capped at 200 (provisional). Committed and merged.
- **Decomposition and enumeration** (done 2026-11-14, ROADMAP §6): Chemistry ▸ Decompose
  (R-Group Decomposition, Fragments, Central Ring), Chemistry ▸ Enumerate (Reaction
  Products with 28 named templates in `analysis/reaction_templates.py`, R-Groups on a Core,
  Recombine Fragments); Standardize keeps stripped salts in a `Salts` column and can form
  salts (`analysis/salts.py`). The Structure Viewer is now Data ▸ Browser (bottom section,
  Ctrl+Shift+B); Ctrl+Shift+3 is gone. Merged to `main` 2026-11-14 (full suite green).
- **Prepare Reagents, and R-groups from a file** (done 2026-11-15, at the user's request; the
  plan is in Claude's plan file on the Windows PC). Chemistry ▸ Enumerate ▸ **Prepare
  Reagents…** turns reagents into substituents: `chemsuite/analysis/reagents.py` holds 13 kinds
  (acyl chloride, carboxylic acid, sulfonyl chloride, chloroformate, isocyanate, aryl halide,
  boronic, alkyl halide, aldehyde, alkyne, amine, alcohol, thiol), each a one-reactant
  transform whose product carries a plain `*`, in **tiers** so an aryl halide loses iodide
  before bromide before triflate before chloride. Output is a new window (one row per
  fragment: `Fragment`, `Reagent`, `Sites Found`) or columns. **R-Groups on a Core** gained
  **From a File…** per position (`tools/fragments/rgroups_library.py`), reading the file in the
  job, so the count is unknown until the run.
  - **The user's choice:** *every* matching site is prepared (identical fragments collapse, so
    a symmetric diacid gives one). A fragment from a bifunctional reagent therefore carries the
    other group live, which `Sites Found` is there to expose; the design review would have
    skipped those instead.
  - Found and fixed on the way: a **delimited file's structure column is recognised by its
    header**, so a saved `Fragment` column loads as text and `load_blocks` saw no structure
    column — it now takes a predicate and looks past structure columns, which also lets it pick
    the fragments out of a set that lists reagents beside them. And the **Core field's Draw…**
    ran every drawing through `drawn_query_text` (right for a search, wrong for a core), so a
    drawn R1 lost its number; `rgroup_library.core_from_drawing` fixes it and `ReactionText`
    takes the converter from the field.
  - **Still to do:** the by-hand run in the app (prepare a reagent column, save it, build
    amides from a drawn core with it). Not yet tried outside the tests.
- **Unified Similarity Search** (done 2026-11-17, all five phases of `docs/SIMILARITY_PLAN.md`):
  the engine is `chemsuite/similarity`, the dialog `tools/similarity/search.py` with
  `references.py` and `methods.py`; Superpose's pose writing is shared in
  `tools/superpose/poses.py`. **Still to do:** the by-hand run in the app (plan §Verification:
  ECFP4 + 3D Shape against two references with Keep Aligned Poses; an SDF from `samples/`
  screened with 3D on; cancel midway). Known gap: a file screen's second stage (every method
  over the kept hits) does not stop on cancel, so a cancel with 3D on still runs 3D on them.

- **Merged to `main` 2026-11-24 (Windows): `fe1361a`.** Six commits — Auto QSAR, the ligands
  Kekule fix, the metal-site restraints, CYP phase C, this file, and the GPL headers. Two
  sessions were editing the repo at once, so `docs/ROADMAP.md` and `docs/USER_GUIDE.md` were
  split hunk by hunk across the commits they belong to. The full suite passed on Windows
  (3361 passed, 48 skipped, 2 xfailed, 15:11).
  - **CI's macOS job is red and was red before this work**, on
    `tests/metabolism/test_sites.py::test_the_shipped_model_reproduces_the_golden_scores`:
    paracetamol on 1A2, 2 of 11 atoms off by 0.0066 (0.0412 against the golden 0.0347), far
    outside the 5e-07 the test allows. Windows and Ubuntu pass, so it is a macOS arm64
    numeric divergence in the shipped GBM, not a wrong model. Whoever picks up the SoM work
    should decide whether the golden scores get a platform tolerance or the export is made
    bit-stable; until then a red CI on macOS means nothing new.
  - CI's `lint` job had also been failing for several pushes: six files had no GPL header.
    `scripts/check_gpl_headers.py --fix` writes CRLF, so run `ruff format` or convert to LF
    before committing on Windows.

## Next

- **Open items on finished features, closed 2026-11-19 on Windows** (ROADMAP §7.0b–f; plan in
  Claude's plan file on the Windows PC, `hashed-tickling-liskov.md`): kept rows whose protein
  was removed show a red Protein cell; faded structure in front of a focused ligand is cut (the
  three `pro_*.frag` shaders rebuilt with qsb from a scratch PySide6-Addons 6.12.0); MHFP6 for
  TMAP; TMAP measured to 250k on ChEMBL 37 and beside tmap-silicon
  ([docs/chemspace/bench.md](docs/chemspace/bench.md)); 2D work on a 3D column reads the
  cell's SMILES, and 2D descriptors no longer read InChI stereo from the first conformer.
  - **On the Mac:** run `tests/core/test_children.py`, `tests/core/test_parallel.py` and
    `tests/app/shell/test_exit.py` (§7.0's process-group path has never run), and check the
    focused-pocket cut on Metal (Focus a ligand, turn the view: no ribbon in front of it).
  - **TMAP's layout fixed 2026-11-24** (ROADMAP §7.0d, bench notes): repulsion reaches 16 edge
    lengths on coarse levels and 8 on the finest; crossings halved, coverage matches
    tmap-silicon; layout about 2× slower (27 s at 250k). Open: the small separate trees at
    250k sit in a thin strip beside the giant one (`pack`).
  - **Next chemical-space views** are on the roadmap (§7.0h): GTM, a scaffold tree view, PCoA;
    TriMap, PHATE and Isomap considered and not planned.
  - **Found, not fixed:** British spellings the user sees: Prepare's "Optimise" option, TMAP's
    "Neighbours" row, Settings "New Colours", the custom theme's "Colours"/"Choose Colour",
    plot tooltips with "colour/coloured", MPO's "Near a centre", Standardize's "stereocentre".
    Saved keys keep their spelling.
  - `samples/chembl/chembl_37_chemreps.txt.gz` (2,897,819 compounds, Windows PC only, not
    committed) feeds the benches above 49k; copy it to the Mac to rerun them there.

- **In progress: sites of metabolism** (ROADMAP §7.0f, plan and the whole record in
  [`docs/SOM_PLAN.md`](docs/SOM_PLAN.md); Claude's plan file on the Windows PC,
  `polished-wandering-candle.md`). **M0, M1 and M2 are done** (2026-11-18). The route went through
  ATTNSOM (reproduced, then found to be measured on a leaking split) to an independent benchmark
  that showed a plain atom-feature GBM ties it; the plan's own sections hold that history, and the
  **decision with the user on 2026-11-18 was the GBM, `PRIMARY_SOM` labels only**.
  - **`chemsuite/metabolism/` ships**: `features.py` (43 columns, ported from `som-benchmark`,
    MIT, and checked identical to it), `trees.py` (the exported trees walked in numpy),
    `sites.py` (`predict_sites`, `top_site`, ties broken by lowest atom index),
    `predict.py` (chunked, cancellable), `data/som_gbm.npz` (**250 KB**, vendored, no download,
    no consent prompt) and `data/LICENSE-SOMBENCH.txt`. Added to `QT_FREE`; 35 tests in
    `tests/metabolism/`.
  - **Training is in the app's venv now**, not WSL: `scripts/som/train_gbm.py` (measure and fit),
    `export_gbm.py` (sklearn → `.npz`, refusing to write if the numpy path disagrees by more than
    1e-9; measured 1.1e-16) and `golden.py` (`tests/metabolism/som_golden.json`). The whole run is
    **24 seconds** on the Windows PC. **It needs the curated table**, which the WSL
    `som-benchmark` clone produces (`scripts/01_download.py`, `02_curate.py`) and which is copied
    to `samples/som/instances.csv` — local test data, not committed, and the Zaretzki labels never
    ship.
  - **The numbers to quote** (molecule-grouped, `PRIMARY_SOM`, 10 seeds): MCC 0.475 ± 0.041,
    top-1 **0.592**, top-2 **0.727**, top-3 **0.813**; per isoform in the plan, 3A4 weakest
    (top-1 0.526) and 2C19 strongest (0.683). Never the instance-split figures, and no comparison
    against SMARTCyp/XenoSite/FAME — every published number in this field is on the inflated
    basis.
  - Chemistry spot-check is in the tests: testosterone → C6 (3A4's 6β-hydroxylation),
    paracetamol → the amide N, caffeine → C8, ethanol → the carbinol C for 2E1, each the
    top-ranked atom.
  - **M4 and M6 are done too** (2026-11-18): Chemistry ▸ Calculate ▸ Sites of Metabolism…
    (`chemsuite/tools/metabolism/`, the two-file Predict pKa pattern) with the four common
    isoforms on and the other five under More Options, and Top Site / Score / Sites /
    Applicability columns per isoform; the guide section is `### Sites of Metabolism`, with the
    accuracy table, the `PRIMARY_SOM` definition and why the numbers are below published ones.
  - **M5 is done (2026-11-18): the sites are marked on the drawings.** It is an app feature,
    not the tool's: a structure column marks the atoms *any* column names (`C7`, or
    `C7 0.81, N3 0.05`), through **Highlight Atoms ▸** in its header menu, each name checked
    against the molecule in that row. New: Qt-free `core/atom_names.py`; `graded_colors` /
    `drawn_highlight` / `render_png(highlight_colors=…)` in `core/depict.py`; and
    `app/table/highlights.py`, which now owns every kind of cell highlight (`Highlight`,
    `AtomSource`, `AtomHighlight`) — `table_view.py` ended 3 lines shorter than before, and
    `ColumnActions` takes the view instead of four callbacks. Saved with the window layout, so a
    reopened document still marks them; a finished run marks the most prescribed isoform it was
    asked for; `ToolContext.highlight_atoms` is a tool's only way in. Full suite run (shared
    rendering). The guide is updated in the column-menu list and under Sites of Metabolism.
  - **The M5 by-hand check is done (2026-11-23, Windows)**, by a script that drives the real
    window offscreen (not committed; it is in the scratchpad). A 2D SMILES column and a 3D
    column read from an SDF mark the same atoms in every row, light and dark, and the dark
    theme's faded marks stay legible on the dark base.
    - **The index shift is real in the 3D case and is handled.** The SDF reader keeps a 3D
      record's explicit hydrogens (`io/sdf.py`, `removeHs=False`), so testosterone is 49 atoms
      in the table and its top site is `C15`, where the 2D cell's is `C7`; the drawing leaves
      the hydrogens out and `depiction_mol`'s `_depictSourceAtom` tags carry each site to the
      atom it means. A molblock pasted into a cell goes the other way (`core/mols.py`,
      RDKit's default `removeHs=True`), so that column has no hydrogens to shift.
      `tests/core/test_depict.py` now renumbers a molecule with its hydrogens *between* the
      heavy atoms — the earlier test used `AddHs`, which appends them and so never moved an
      index; without the tags it marks a carbon instead of the nitrogen.
    - Sites columns are named for the molecules of the column the run saw, so pointing a 3D
      column's drawings at a 2D column's Sites marks almost nothing. That is the guard working
      (element and index are both checked), but element and index can coincide, so a stray
      mark is possible: worth a line in the guide only if it ever confuses anyone.
  - **The hover preview marks the atoms too (2026-11-23**, at the user's request). The large
    drawing a structure cell shows on hover came from the model's tooltip and so knew nothing
    of highlights; `CellDelegate.helpEvent` now answers for a marked structure cell and
    `preview_html` takes the highlight and the column's alignment, so the enlarged drawing is
    the cell's drawing: a structure search's match is marked too, and an aligned column's
    layout is kept.
    - **A theme bug came out of it:** `_dark()` read `QToolTip.palette()`, which Qt freezes the
      first time anything reads it — the preview itself, on the first structure hovered — so
      after a switch to the dark theme in Preferences the preview drew a white-backed
      structure on the dark tooltip. It reads the application's palette now. The test runs the
      start-up-then-switch sequence in a process of its own, since the state cannot be
      reproduced once a suite has touched the palette. This is the same family as the menu-bar
      corner `QToolButton` issue above (a run-time theme switch not reaching a widget).
    - `table_view.py` was at its line budget, so the drawing half moved to
      **`app/table/cells.py`** (`Alignment`, `StructureImages`, `CellDelegate`); the budget in
      `tests/test_architecture.py` came down from 1,139 to 976.
  - **Phase C, step C1 is done (2026-11-24, Windows): CYP inhibition has data and a baseline.**
    Three uncommitted-data scripts in `scripts/cyp/` — `fetch.py` (OpenADMET's
    `cyp-challenge-train-test` at the pinned commit `5ab9dd9`, each file checked by SHA-256,
    into the uncommitted `samples/cyp/`), `curate.py` (chemsuite's own standardizer, InChIKey
    dedupe) and `benchmark.py` (ECFP4 counts + the 104 cheap 2D descriptors, Murcko
    scaffold-grouped folds, 3 seeds, a shuffled-label null, the unblinded 750 as a check).
    `docs/SOM_PLAN.md`'s Phase C section is rewritten as a record of what was measured.
    - **The data is clean but small.** 5,655 rows gave 5,655 distinct compounds: no duplicates,
      no salts, nothing unreadable, no train/test overlap, and **no censoring** (the pIC50s are
      continuous, one row at the minimum per isoform). But 4,905 is the *union*: per isoform
      the training set is only 1,285–2,335 compounds.
    - **Scaffold-split out of fold:** 3A4 R² 0.570 / Spearman 0.753; 2C9 0.338 / 0.574;
      1A2 0.226 / 0.485; 2D6 0.096 / 0.331. The null came back at noise (R² −0.12 to −0.16), so
      nothing is leaking. Accuracy tracks `n`, not chemistry: the limit is data, not method.
    - **2D6 scores R² −0.670 on the unblinded file** — worse than the mean — because that set
      was built by hit expansion over 1A2/2C9/3A4 hits only (its assay was still in
      development), so its test mean is 3.27 against a training mean of 4.78. The external
      check is unusable for 2D6, and its training values may not share the others' assay
      footing.
    - **The recommendation written into the plan: ship 3A4 only**, gated on three things —
      the license (Apache-2.0 is Hugging Face *metadata only*; no license text, no copyright
      holder, no word on the measurements), QSAR step 4 for `apply.py` and `.csqm`, and a
      scaffold-held-out slice kept aside from the start. **Not yet accepted by the user.**
    - Still open from the plan: whether predicted pKa earns its cost (seconds per molecule
      against microseconds), measured on a subsample; and running the pool through
      `qsar.search.run` (`benchmark.py --engine`) to see what the full search adds. Neither was
      run, to avoid two CPU-heavy jobs at once.

- **Saved for later: SQL import** (File ▸ Import from Database, read-only queries; SQLite
  built in, PostgreSQL/MySQL/Oracle/SQL Server through optional extras). The plan was approved
  2026-11-12 and is in [`docs/DATABASE_PLAN.md`](docs/DATABASE_PLAN.md). No code yet.

- **Parked, high priority: minimise after Prepare** (ROADMAP §7.0g, plan approved 2026-11-15
  in [`docs/SIMULATE_PLAN.md`](docs/SIMULATE_PLAN.md)). **M0 is done** (2026-11-15): all six
  feasibility gates pass on the Windows PC, natively and in WSL; the result and the changes it
  makes to the plan are in [`docs/simulate/m0-review.md`](docs/simulate/m0-review.md).
  **M1 is done** (2026-11-15; what it changed is listed under M1 in the plan): Prepare caps cut ends (ACE/NME from SEQRES), builds a residue missing alone (the
  user's choice over crowded caps), puts hydrogens on every kept water, charges ions by
  component/element (`prepare/ions.py`, one table with `energy.py`), turns rebuilt side chains
  clear, renames the neutral N-terminus state `N-neutral` (saved reviews migrate);
  `ligands.prepared_form`; `edit.insert_residue`; the Qt-free `chemsuite/simulate`
  (`spec.py`, `build.py`, `result.py`; tests in `tests/simulate`). **M2 is done** too (the
  shared `core/jobs/lines.py` client, retro moved onto it, `simulate/worker.py`, `runner.py`,
  `config.py`, the External Tools rows, `envs/openmm-*.yml`; see M2 in the plan).
  **M3 is done** (2026-11-15: `Prepared.built`, restraints per class of atom, tightening,
  `simulate/minimise.py` measuring the result; the tuning findings are under M3 in the plan).
  **M4 is done** (2026-11-15: Minimise in Prepare and Protein ▸ Minimise…, each result a
  new protein `X-prep`/`X-min`; see M4 in the plan). **M5a is done** (2026-11-24, Windows PC:
  `scripts/review/minimise_review.py`, `tests/simulate/test_bench_minimise.py`; numbers under
  M5 in the plan): 167 of 262 neutron entries minimize, RMSD and clash score pass, neutron marks
  unchanged, surface Lys NH₃ swing in GBn2. The five Zn sulfonamides (4G0C…) now build (an
  RDKit radical on aromatic S), and every transition-metal site is held as Prepare left it
  (they collapsed onto the bare ion before). **Next: the 14 defects left** (TYM with no
  template; uncapped chain ends; Mn(III); 13HO's guessed bonds; heme refused as a raw OpenFF
  error), lowering the test's `DEFECTS` ratchet as each is fixed.
  The run's outputs are local (`samples/prepare/out/{chemsuite-prep,chemsuite-min,minimise}`),
  so the Mac needs its own run (`run --only ID,…` reruns some). On Windows `run` keeps the PC
  awake itself (it slept mid-run once); on the Mac use `caffeinate -i`. A by-hand run with the real env (Prepare + Minimise on 181L)
  is still owed. Two constraints for M4
  found in M3 are in the plan's M4 section (`built` is kept and unioned across re-prepares;
  minimise only from the prepared model). To run the real-env tests
  here: `$env:CHEMSUITE_OPENMM_PYTHON = "$HOME\miniforge3\envs\cs-sim-m0\python.exe"` then
  `pytest -m openmm` (about 5 minutes). On the Mac, before the merge: `pytest tests/retro`
  (its fake worker is a shell script, so it skips on Windows) and gate 1. The P1 neutron marks hold after M1 (Asn 86.2, Gln 80.3, His 89.6,
  tautomer 63.8%; Reduce 85.4/78.9/88.3/62.5).
  - **Open, found 2026-11-15, older than M1:** the P2 mark fails. PKAD-R's pickled tables were
    from before P3 (20 energy terms, now 26) and `pkad_review.py chemsuite` crashed on them, so
    `chemsuite.csv` was stale (1.32). Rebuilt (`tables --redo`): chemsuite 1.48 against PROPKA
    1.33 (null 1.15). The same (1.475) with M1's ion charges, side-chain turning and
    true-C-terminus OXT rule all switched off, so P3's terms or defaults moved it;
    `pkad_review.py fit` against the new terms is the likely fix. `test_chemsuite_pka_values_match_propka` (`-m bench`) fails until then.
  - `waters.complete` is pure Python, about 3 ms a water (linear): fine for crystal
    structures, about a minute for 20,000 waters. Vectorise it if large cryo-EM boxes get slow.
  State on the Windows PC:
  - Miniforge 26.9.4-2 is installed per user at `C:\Users\hunte\miniforge3` (not on PATH, not
    the registered Python). `~/.condarc` sets `solver: classic`, so pass `--solver libmamba`.
  - The spike envs (openmm 8.4.0, openmmforcefields 0.15.1, CPU pytorch, `cuda-version=12.6`):
    native at `C:\Users\hunte\miniforge3\envs\cs-sim-m0` (moved off D: 2026-11-16; it is
    Settings ▸ External ▸ OpenMM Python, and passes Test with CUDA on the RTX 4070), and WSL at
    `~/miniconda3/envs/cs-sim-m0`. Recipes are in the review.
  - Hardware: RTX 4070 (driver 560.94, CUDA 12.6 at most; OpenMM 8.6 needs a newer driver).
  - Still owed: gate 1 on the Mac (osx-arm64 env and its platforms), and the choice between
    OpenMM 8.4.0 with this driver and 8.6.1 with a newer one.

1. Conformers at scale (§7.0f): measured and fixed 2026-11-17
   ([`docs/conformers/scale.md`](docs/conformers/scale.md); data from
   `scripts/review/conformer_scale.py data`, local under `samples/conformers/scale/`). Open
   for a decision: 2D-only work on a 3D column loads every conformer (the report's two
   options). Owed: the by-hand pass on a 1k × 300 document and the same benches on the Mac.
   Conformer energies (2026-11-17) are also owed a by-hand run in the app with 3D drawing.
2. Protein alignment (§7a 1): plan approved 2026-11-26 in
   [`docs/ALIGN_PLAN.md`](docs/ALIGN_PLAN.md), phases A0–A6. A0 (licence decisions, the
   benchmark against MAFFT L-INS-i, the numpy DP kernel timing) is next. Bioactive conformers
   (§7.0b) still waits, starting in plan mode.
   (Graceful exit is done and checked on Windows; on the Mac, run its three tests and the
   by-hand checks in ROADMAP §7.0.)

Still owed by the user: the ChemDraw hotkey list.

## Protein code map

- `chemsuite/macromolecule/` (Qt-free outside `ui`):
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
    (`MacromoleculeWorkspace`) is the document part (`protein_workspace`) that keeps it and the view.
  - `ui/`: the Qt workbench. `session.py` is `SessionSignals` (the session's events as Qt
    signals, one per window);
    `host.py` holds `LinkedRow`; `panel.py` (a `ToolPanel`), `menubar.py`, `menus.py`,
    `export.py`, `tree.py`, `sequence_view.py`, `selection_bar.py`, `models_bar.py`,
    `pocket.py`, `labels.py`, `overlay.py`; `builds.py` runs builds as background jobs and
    drops stale results; `view.py` is on `RhiView`; `shaders/` holds the `pro_*` shaders.
- `chemsuite/tools/macromolecule/`, the bridge: `tool.py` (the tool, its window, open, fetch,
  remove), `rows.py` (the ligand-row columns, `ligand_rows`), `follow.py` (`TableRow`, the
  current row followed in the workbench), `to_table.py` (Ligands, Add Ligand and
  Interactions to Table), `dialogs.py`.
- Tests: `tests/macromolecule/{structure,chem,look,ui}`, `tests/tools/macromolecule/test_protein_tool.py`,
  and the benchmarks (`pytest -m bench tests/macromolecule/test_bench_protein.py`).
