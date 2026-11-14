# chemsuite roadmap after the Data menu (2026-11-05)

`PLAN.md` covered the rebuild of MolManager's core: the table, files, descriptors, plots,
conformers and the Data menu. This roadmap covers what comes next:

- the infrastructure MolManager has and chemsuite still lacks (job queue, process log);
- improvements to the tools chemsuite already has, starting with search;
- the rest of MolManager's features, ranked by value and install cost.

Decisions so far: QSAR (phase 7d) is skipped for now. Metabolite prediction (BioTransformer) is
dropped from the plan. Groovy and custom themes stay. External
tools (Gnina, AiZynthFinder…) will be integrated, so their paths are already in
Settings. Settings… sits in the app menu on macOS and the Edit menu elsewhere. Tool dialogs close on OK
and their jobs run in the background (§2). A robust drawing sketcher that works with the table
is needed (§6, tier 2).

---

## 1. Where things stand

| Area | chemsuite now | MolManager |
|---|---|---|
| Background work | `JobRunner`, 3 threads plus a process pool. Every job shows a **modal** progress dialog, so you can only wait or cancel. | A serial queue for heavy jobs, a fast lane for quick ones, and a **Processes** window (running, queued, elapsed time, cancel, clear queue). It keeps the computer awake while jobs run. |
| Logging | Python `logging` to stderr only. Errors pop up a dialog. | An in-memory session log of 10,000 entries (logging, status messages, job events), a **Log** window with level filter and text search, and an audit log for IT. |
| Search | A quick search bar (all columns; match case, whole cell, regex), Find and Replace, and filter cards (range, text, category). | A search panel with several criteria rows (column + text, partial or case options), **substructure search** with AND/OR term groups, saved with the session. |
| Settings | In progress: settings core, IT policy, themes, command IDs, Open Recent, status bar toggle. | Theme, fonts, hotkeys, tool paths, status bar. |

---

## 2. Job queue and Processes panel

**Goal:** long jobs no longer block the window. You can queue several, keep working, watch
progress, and cancel or clear them.

**Status (2026-11-05): built.** Two lanes, job records, the status-bar indicator, the Jobs
window (View ▸ Jobs), non-blocking failure reports, dropped results for replaced documents,
the close prompt and keeping the computer awake. Saving still shows a short modal progress
dialog, because the close prompt needs to know whether it worked. Done 2026-11-07: removing
ended jobs one by one (or several) in the window, and a desktop notification when a long job
ends unseen (osascript on macOS, notify-send on Linux, a tray balloon on Windows; Settings ▸
General can turn it off).

**Design**
- **Two lanes in `JobRunner`.**
  - The **queue** runs one heavy job at a time, in order. Heavy jobs are those that use the
    process pool or can take minutes: conformers, descriptors on big tables, MMP, UMAP and
    t-SNE, clustering, and later docking.
  - The **interactive lane** runs quick jobs in parallel: opening a small file, parsing,
    one-off similarity.
  - Each tool says which lane it uses, with a default by expected size (as the existing
    "warn above N rows" checks already do).
- **Job records.** Each job has:
  - an id, title, document, lane and state (queued, running, done, failed, cancelled);
  - progress, and the queued, start and finish times;
  - a one-line result ("3,374 pairs") or the error with its traceback.
- **No more modal progress.** A tool's dialog closes when you press OK, and the job shows in:
  - the **status bar**: one line. At the far left is the current job as text ("Generating
    conformers… (120/500) · 3 queued", or "Calculating UMAP…" when it has no count), then the
    latest message. The row and selection counts sit at the far right. There is no progress
    bar. Click the job text to open the panel.
  - a **Jobs panel** (View ▸ Jobs). It lists running, queued and finished jobs with elapsed
    time. You can cancel a job, remove a queued one, move a queued job up or down, clear
    finished jobs, and copy an error.
  - Results apply when the job finishes. `apply_results` already skips rows deleted in the
    meantime. Each tool's existing checks still run before a job is queued.
- **Cancelling.** Jobs stop within about a tenth of a second: work handed to worker processes
  is polled for a cancel instead of waited on, and in-thread work runs in small chunks. While
  it stops, the status bar says "Cancelling …". Tools that write per-row results keep the
  rows they finished: a cancel writes those and says "partial results kept". This covers
  conformers, descriptors, medchem-space descriptors, similarity search, R-group
  decomposition (with a numbered core) and diverse picking (the picks made so far). It is the
  rule for every new tool. Whole-dataset results (clustering, UMAP/t-SNE/PCA, SOM,
  Free-Wilson, Hansch, matched pairs, activity cliffs) mean nothing until complete, so they
  are discarded.
- **Document safety.** Closing a document cancels its queued jobs and asks about running ones.
  Quitting with jobs running asks first. Saving waits for a running job that writes to the same
  table.
- **Long jobs.**
  - The computer stays awake while queued jobs run, as in MolManager: `caffeinate` on macOS,
    `SetThreadExecutionState` on Windows, and `systemd-inhibit` where it exists on Linux.
  - When a job longer than 30 s finishes in the background, the window shows a notification.
- **Tests:** ordering, cancelling running and queued jobs, document close, quit, failure
  reporting, and the interactive lane not waiting behind the queue.

## 3. Process log

**Goal:** every job, warning and error is visible in the app. Support can ask for one log file.

**Status (2026-11-05): built.** Logbook, Log window (View ▸ Log), rotating log file with
Help ▸ Show Log Folder, policy `log_level`, crash handler with Copy Details, and the audit log
(policy `audit_log`: `true` or a file path) for start-up, open, save and export. Network and
external-tool runs will write audit records when those features arrive.

**Design**
- **`core/logbook.py`** (no Qt): a thread-safe ring buffer of 10,000 entries fed by a `logging`
  handler, plus structured job events (queued, started, finished with duration and row counts,
  failed with traceback).
- **Log panel** (View ▸ Log):
  - filter by level and by text, with colour per level;
  - follow mode, Copy, and Save As…;
  - a click on a job's entry jumps to it in the Jobs panel.
- **Log file:** a rotating file in the user's log folder (5 × 2 MB), with a Help ▸ Show Log
  Folder command. IT can set its level through the policy file.
- **Crash handler:** uncaught exceptions on the GUI thread or in jobs are logged in full. The
  dialog shows a short message with "Copy Details" instead of a raw traceback.
- **Audit log** (from `PLAN.md` §2.5): file opens, saves and exports, network use, and external
  tool runs. It is written only when the policy enables it.

## 4. Search and filtering

**Goal:** find compounds by structure, and by more than one condition, in one place.

**Status (2026-11-05): structure search (item 1) built.** Edit ▸ Structure Search… (Ctrl+Shift+M)
and a structure card in the Filters panel (the default card for a structure column). It covers
substructure from SMILES, SMARTS or a drawing, exact match with or without stereo, and similarity
above a threshold. The dialog closes on Search and the job runs in the background. Matches can be shown
alone (as a filter card), selected, or added to a list. Substructure matches are highlighted in the
table (Edit ▸ Clear Structure Highlighting). Speed: pattern fingerprints are cached per column, and
the bench (`tests/analysis/test_bench_structure_search.py`) measures 0.5–0.9 s over 200k rows
once they are cached. The first search of a column computes them (about 8 s at 200k rows). Not
yet: nothing (structure filters are saved with the other filters, see item 4).

**Status (2026-11-05): advanced search (item 2) built.** Edit ▸ Advanced Search… opens a panel
docked with Filters. Each condition picks a column or "any column", an operator that fits the column
(text, number or structure), and a value. Each row has Not and match-case options, and the rows
combine with all or any. Results can be shown alone (one table filter, ANDed with the Filters
panel), selected, or added to a list. Text and number conditions follow edits as you make them.
Structure conditions (substructure or exact) run as background searches and run again when the
structure column changes. The search, whether it is showing, and the last 10 searches are
saved with the document. Show All Rows clears it. Core: `chemsuite/core/table/criteria.py`.

**Status (2026-11-05): filters (item 4) built.** Filter cards are saved with the document,
including their settings and on/off state, and the Filters panel opens on load when there are
any. Named filter sets (Filters ▸ Sets) save the current cards under a name and switch between
sets. Data ▸ Copy Shown Rows to New Document opens the shown rows, in view order, as a new
unsaved document. Like the other view state, changing filters does not mark the document
modified.

1. **Structure search** (Edit ▸ Structure Search…, and a card in the Filters panel):
   - **Modes:**
     - *substructure*, from SMILES or SMARTS, drawn or pasted;
     - *exact match*, ignoring or respecting stereochemistry;
     - *similarity* above a threshold, which reuses the fingerprint cache.
   - **Speed:** RDKit pattern fingerprints screen out non-matches first, and the remaining
     matches run on the process pool. The target is under 2 s for a substructure search over
     200k rows.
   - **Highlighting:** matched atoms are highlighted in the table's structure drawings.
     `render_png` already takes `highlight_atoms`.
   - **What happens to the matches:** show only them, select them, or add them to a list.
2. **Advanced search** (Edit ▸ Advanced Search…, a panel):
   - **Rows of criteria.** Each row is a column (or "any column"), an operator and a value.
     - Text operators: contains, equals, starts with, regex, is blank.
     - Number operators: =, ≠, <, >, between.
     - Structure: substructure.
   - Rows combine with "all" or "any", and each row can be negated.
   - The result can show only matches, select them, or add them to a list.
   - Searches are saved with the document and kept in a recent-searches list.
3. ✅ **Quick search bar** (done 2026-11-07):
   - field syntax such as `ID:CPD-12`, `pIC50>7` or `"exact phrase"`;
   - a history dropdown, and a count of matches as you type.
4. **Filters:**
   - a **substructure card**;
   - filters **saved with the document** (a noted follow-up);
   - named filter sets ("Actives, MW < 450") you can switch between;
   - "Copy filtered rows to a new document".
5. ✅ **Find and Replace** (done 2026-11-07): searching within the selection or within chosen
   columns, and replacing with regex groups.

## 5. Improving the other existing tools

- ✅ **Table** (done 2026-11-07)
  - Conditional cell colouring: colour by value or category, like MolManager's
    `column_color.py`.
  - Per-column number formats (fixed, thousands, scientific, significant figures, percent).
  - Frozen columns (the first columns, up to the one chosen).
  - Column statistics in the header tooltip.
- ✅ **Structures**
  - ✅ Align all drawings to a common core (MolManager's "Align to Substructure"), done
    2026-11-07.
  - ✅ Larger hover previews, done 2026-11-07.
  - ✅ "Copy as SMILES / image" (cell menu ▸ Copy As).
- **Plots**
  - Colour-range controls and marker shapes.
  - A density mode (hexbin) for 200k points.
  - A radar plot for MPO profiles.
  - Several Y columns on one plot.
  - Saved plot templates.
- **Browser:** sort options, a size slider, and pinning compounds to compare them
  side by side.
- **SAR tools**
  - MMP neighbourhood map and transform ledger (from MolManager).
  - An activity-cliff network plot.
  - Free-Wilson proposals sent to a design sheet (MolManager's Design Studio idea, without its
    pKa heuristics).
- **Similarity and clustering:** a cluster browser (representative compound plus members), and
  Jarvis-Patrick and sphere exclusion (noted earlier).

## 6. MolManager features not yet in chemsuite

**Tier 1: no new dependencies, high daily value**

| Feature | Notes |
|---|---|
| Calculator / formula columns | A safe expression language over columns, with math functions and chemistry functions such as `MW(Structure)`. Columns recompute when their inputs change. |
| Potency standardisation | IC50, Ki or EC50 in any units to pIC50 / pKi, with unit detection. |
| Structure preparation | RDKit MolStandardize: salt strip, largest fragment, neutralise, canonical tautomer, add or remove Hs. "Fast Prepare" uses sensible defaults. |
| Scaffolds | Bemis–Murcko and generic scaffolds as columns, plus a scaffold tree or count view. |
| Decomposition and enumeration | **Done 2026-11-14:** Chemistry ▸ Decompose (R-Group Decomposition, Fragments by BRICS or RECAP, Central Ring) and Chemistry ▸ Enumerate (Reaction Products with 28 named templates, R-Groups on a Core, Recombine Fragments). |
| Documents | Duplicate document, send selected rows to a new document, and Notes (a notepad saved in the `.csd`). |
| Macros | Record and replay sequences of commands with their settings. The new command registry is the base for this. |
| Random utilities | Random numbers and a random molecule sample. Low priority. |

**Tier 2: external tools or network access (optional, controlled by policy)**

| Feature | Notes |
|---|---|
| Confab and CONFORGE conformers | Through the external-tool paths now in Settings. |
| Superpose conformers, pharmacophore screen | RDKit only. The Browser's 3D view (Data ▸ Browser) is done and aligns a molecule's conformers; superposing different molecules is next. |
| PubChem, ChEMBL, vendor, patent and SQL queries | Network access, off by default under policy. Plain HTTPS REST calls, not MolManager's client libraries. SQL through the standard library or an optional driver. |
| Retrosynthesis | AiZynthFinder through its own Python interpreter (path in Settings). |
| Sketcher | **Decided: needed.** A robust drawing sketcher that works with the table: edit a structure cell in place, draw a query for structure search and filters, and add new compounds. Plan: [SKETCHER_PLAN.md](SKETCHER_PLAN.md) (native, Ketcher rules as reference, verified; no MolManager iupac_* code). |

**Tier 3: heavy add-ons (a separate optional install)**

| Feature | Notes |
|---|---|
| pKa and permeability models | torch, unipka and chemprop make the install about 2 GB larger, so this is an optional "predictions" add-on. Settles the open question about the "(est.)" descriptors. pKa: MolGpKa is built into the base install (Chemistry ▸ Predict pKa); Uni-pKa stays an optional add-on, see [PKA_PLAN.md](PKA_PLAN.md). |
| Protein viewer and editor, sequence alignment, docking, preparation, MD simulation | Planned in phases: see **§7a Proteins** below. |
| QSAR, Quick QSAR, CoMSIA | Skipped for now (phase 7d). |

---

## 7. Proposed order

### 7.0 High priority: graceful exit (added 2026-11-11)

Closing the main window and accepting the exit prompt must shut chemsuite down completely: no
windows left open and no Python process, worker process or external program still running.
This must hold whatever is in flight: a job in the queue or quick lane, a file loading or
being exported or saved, structures being parsed, a download, or an external tool (the warm
AiZynthFinder worker, and later Vina, Gnina, OpenMM).

What happens today: the window cancels its jobs and closes. Then, on `aboutToQuit`,
`JobRunner.shutdown` cancels everything, waits up to 5 s for each thread pool, and shuts down
the process pool without waiting, and the retro worker is closed.

Gaps to close:
- work that cannot see its cancel token (a long RDKit call, a C extension) keeps its thread
  alive past the 5 s wait, and Python then waits for it at exit;
- the process pool's workers are not waited for or killed, so a busy one can outlive the app;
- external programs must be terminated, then killed after a grace period, including their
  child processes;
- floating panels (the Protein Viewer), dialogs and error reports must close with the last
  document window;
- a save or export cut short must leave the old file intact (write to a temporary file, then
  rename), never a half-written one;
- the "jobs still running" prompt must count every job that would be stopped, background
  work included.

Done when a test starts each kind of work, quits mid-way, and finds no child processes, no
live threads and no open windows; and when a by-hand check on macOS, Windows and Linux shows
the process gone from Activity Monitor, Task Manager or `ps`.

### 7.0b High priority: bioactive-quality conformers (added 2026-11-11)

Good drug-like conformers are essential to chemsuite: docking, 3D similarity, pharmacophores,
superposition and the protein work all start from them. Today Chemistry ▸ Generate
Conformers is RDKit's ETKDGv3 followed by MMFF94 (or UFF) optimisation
(`chemsuite/conformers/etkdg.py`).

The goal is a generator that is **measurably better than that baseline at reproducing the
bioactive conformation**: the pose of a ligand crystallised in its protein.

- **Measure first.** Build a benchmark of protein-bound ligand structures from the PDB (a
  Platinum-style or PDBbind-style set with a licence that allows it). Score the best heavy-atom
  RMSD to the crystal pose at 10, 50, 100 and 300 conformers, the fraction of ligands within
  1.0 Å and 2.0 Å, ensemble size and time per molecule. ETKDGv3 + MMFF is the number to beat.
  The benchmark runs as its own `-m bench` set, not in the regular suite.
- **Then improve.** Candidate steps, each kept only if it improves the benchmark:
  - more and better-diversified ETKDG embeddings;
  - macrocycle-aware settings;
  - correct protonation and tautomer states;
  - torsion-library (CSD-derived) knowledge;
  - energy windows and RMSD pruning tuned for coverage rather than minimum energy;
  - solvent-aware or implicit-solvent minimisation instead of gas-phase MMFF, which collapses
    polar ligands;
  - an optional ML or semi-empirical rescoring add-on (e.g. xtb) if it pays for its install
    size.

  Compare against open generators (CONFORGE, OMEGA-like settings) where licences allow.
- **Ship.** The new generator becomes the default behind the same tool, and ETKDG stays
  available. It keeps partial results on cancel and runs in the process pool. The user guide
  states the benchmark numbers.

**Process:** this work starts in plan mode. Write the plan (benchmark set, licences, metrics,
candidate methods, acceptance threshold) and get it approved before writing code.

### 7.0f High priority: hundreds of conformers per molecule without slowing the app (added 2026-11-12)

Users want hundreds of conformers per molecule. At about 1,000 compounds × 100 conformers
(100,000 conformers) the whole app slows down. This is separate from 7.0b, which is about
conformer quality; this item is about scale.

- **Measure first.** Profile generation, storage, the table (cells, sorting, filtering,
  copy), document save and load, memory, and the 3D views, at 1k × 10, 100 and 300
  conformers. Find which parts cost time and memory, and which ones block the UI thread.
- **Then decide** how conformers are held. Options to weigh in the plan: a conformer store
  outside the table cells (compact coordinate arrays per molecule, loaded lazily); cells that
  show a count and the selected conformer only; background loading for save and load; limits
  and pruning by default.
- **Process:** plan mode first. The plan states the numbers to reach (UI response, memory,
  save and load times at 1k × 100 and 1k × 300) before any decision on storage.

### 7.0d High priority: TMAP and Chemical Space Networks (added 2026-11-11)

Two more chemical-space plots, next to the Constellation Plot. Both are big new features, so
each starts in plan mode. Their graph layout and drawing code may be shared, so plan them
together.

- **TMAP** (Probst and Reymond, *J. Cheminf.* 2020, 12, 12). A tree map of large sets: MinHash
  fingerprints, an LSH Forest k-nearest-neighbour graph, its minimum spanning tree, and a tree
  layout. Every compound is a point, and the tree's branches show near neighbours. It should
  reach hundreds of thousands of rows. Check the licences of `tmap` and its layout code before
  depending on them, or write our own (the MST and layout steps are standard).
- **Chemical Space Networks** (Maggiora and Bajorath, *J. Comput.-Aided Mol. Des.* 2014, 28,
  795; Zwierzyna et al., *J. Cheminf.* 2015, 7, 47). Compounds are nodes, and an edge joins
  each pair above a similarity threshold (Tanimoto) or each matched molecular pair. A force
  layout places the nodes, colouring them by activity, and edges can be coloured by the
  activity difference (cliffs). It is meant for hundreds to a few thousand compounds.
- Both are plot kinds in the Plot panel. They select rows like any plot, use the shared colour
  maps and ranges, and compute in the background like the Constellation Plot.

### 7.0c High priority: the Protein Viewer as the other half of chemsuite (added 2026-11-11)

chemsuite has two halves. The **data table** is where a medicinal chemist works up and
visualises data and generates chemical information. The **Protein Viewer** is where they
simulate, inspect binding modes, analyse protein–ligand interaction fingerprints (PLIF) and
align sequences (MSA, and global or local pairwise) and structures. The aim is an
industry-standard modelling and simulation suite inside chemsuite, and the Protein Viewer gets
the same polish as the table and the dockable plots and tools.

- **Done 2026-11-11:** Protein is its own menu in the main menu bar (between Data and Help), no
  longer an item under Chemistry.
- **Done 2026-11-11: Phase 1 finished and the polish pass.** Bond orders from the file's own
  component tables or the CCD (cached; guessed ones marked ⚠); ligand rows built off the GUI
  thread; NMR models with a player bar (the future trajectory player); surfaces, interactions,
  pocket skins and large cartoons built in the background; 16-bit mesh chunks; hydrogens and
  surface transparency; Focus with a pocket skin and depth fog; residue and atom labels;
  Around ▸; a row's own 3D pose shown in its pocket; Interactions to Table (the PLIF column
  set); Remove Protein; rows copied with their protein; Save Protein As and safe saves;
  light/dark, the viewer's own ⌘A/⌘C/F, a keyboard sequence cursor, selection shown in the
  proteins list. Benchmarks at 100,000 atoms: an 8 MB mmCIF read in under a second, a selection
  redrawn in 13 ms, 60 fps turning on Metal.
- **Still to polish:** ghosted ribbons in front of a focused pocket can cross the ligand (cut
  them with a front slab or fade them further); marking a kept row whose protein was removed in
  the table itself.
- **Workspace:** the viewer stays its own window with its own menus. Its Protein menu in the main
  window grows as the phases land (align, dock, prepare, simulate), each opening or driving the
  viewer.

### 7.0e High priority: our own Uni-pKa (revised 2026-11-12)

Uni-pKa gives the free energy of each protonation microstate, so linked sites and tautomers
are handled together; MolGpKa treats each site on its own, and in protein preparation it
catches few ligand states (`docs/prepare/`). Its published weights cannot be used (licence),
so chemsuite trains its own on pKaHub's data from Uni-Mol's MIT backbone and runs them in numpy, with
no torch: [UNIPKA_PLAN.md](UNIPKA_PLAN.md). Ligand states in preparation use it first; Predict
pKa, standardise and the descriptors follow once it beats MolGpKa on the test sets.

### 7a. Proteins: from viewer to medicinal-chemistry workspace (revised 2026-11-11)

Proteins are stored in the document; ligands, poses and scores are table rows linked to their
protein by a `Protein` column, so plots, filters and SAR work on them. Each phase is usable on
its own and starts in plan mode.

**Dependencies, decided:** preparation is our own Qt-free code in the base install (hydrogens
from CCD templates: the benchmark found gemmi's `prepare_topology` misses entries without CCP4
monomers; PDBFixer and pdb2pqr are dropped). `propka` and `meeko` stay light pip dependencies (licences checked first).
OpenMM and openfe run in an interpreter of their own, set in Settings ▸ External Tools and
driven as a JSON-lines worker (the `retro/runner.py` pattern); there is no `chemsuite[protein]`
pip extra.

0. **Viewer and selection.** Done (see §7.0c).
1. **Alignment and correspondence** (base install: gemmi and numpy). Pairwise global
   (Needleman–Wunsch) and local (Smith–Waterman) alignment with chosen matrices and gap
   penalties; MSA through the registered MAFFT tool; conservation shading in the sequence
   strip; structural superposition on a whole chain or a pocket, with the linked ligand rows
   moving in the same undo step; predicted models (AlphaFold DB, PDB-REDO) fetched by
   accession. First because everything later needs residue correspondence: ensembles, PLIF
   across proteins, FEP mapping.
2. **Docking and PLIF.** AutoDock Vina as an external program on every platform, Gnina optional
   (Linux, Windows through WSL); the box from a ligand or a selection; Meeko for PDBQT; one
   process per ligand with progress, cancel and finished rows kept; poses as table rows (pose
   rank, Vina score, CNN score, parent) shown in the pocket by the pose display; covalent and
   flexible-side-chain docking; ensemble and consensus scoring; pocket detection and
   druggability; interaction fingerprints across poses, building on Interactions to Table. The
   receptor is prepared by `protein/prepare` (phase 3's P1 is enough for docking).
3. **Native preparation** (base install, high priority; plan approved 2026-11-11 in
   `docs/PREPARE_PLAN.md`).
   - **Progress:** P0 (benchmark) and P1 (repair, hydrogens, flips, His states, rotors;
     `pipeline.prepare`) done 2026-11-12, at Reduce's accuracy or better
     (`docs/prepare/p1-review.md`). P2 (titration: Asp/Glu/His/Lys/Tyr/Cys/N-terminus states
     at a pH, pKa values and populations) done 2026-11-12 at PROPKA's accuracy on PKAD-R, off
     by default (`docs/prepare/p2-review.md`). No UI yet (P4). Next: P3 (ligand states, waters,
     metals); a second titration round (PROPKA-style H-bond determinants) to beat the null
     model and turn titration on.
   - **Repair:** missing heavy atoms and caps from the CCD templates; short loops closed, longer
     gaps capped and marked ⚠.
   - **Hydrogens and states, protein and ligand together:** one coupled site graph covers His
     tautomers, Asn/Gln/His flips, titratable residues, OH/SH/NH₃ rotors, pocket waters, metal
     sites and the ligand's protomers and tautomers. Ligand states start from MolGpKa.
     Desolvation, Coulomb and hydrogen-bond terms act across the interface, so each side's state
     depends on the other's. It is solved exactly per cluster, with populations; near-ties are
     marked ⚠ and become receptor ensembles for docking or separate runs for MD.
   - **Measured first:** benchmarked against neutron structures, PKAD-2, Reduce, PROPKA3,
     Protoss (public PDB IDs through ProteinsPlus) and, optionally, MOE Protonate3D before it
     becomes the default. propka stays as a cross-check.
   - **Relax after preparation (optional, off by default; its own plan, after P1).** Minimise the
     prepared protein–ligand complex so new hydrogens, rebuilt side chains and the ligand settle
     into the pocket, with heavy atoms held by harmonic restraints (≤ 0.3 Å RMSD). One undo step,
     reporting RMSD and energy change. Two tiers:
     - **Quick, in the app:** the pocket (residues within ~8 Å of the ligand), the ligand and every
       hydrogen minimised with RDKit MMFF94, the rest of the protein fixed.
     - **Full, in OpenMM's own interpreter** (phase 5's worker): the whole complex with amber14
       and OpenFF or GAFF2 for the ligand, in implicit solvent, with the same restraints.

     A checkbox under More Options in the Prepare dialog, and a separate command for a structure
     already prepared. Open for its plan: whether MMFF on a capped pocket is good enough or the
     quick tier waits for OpenMM, and the benchmark (heavy-atom RMSD to the deposit and clash
     score on the neutron set).
4. **Editing** (every edit one undoable step). *Done 2026-11-12, ahead of the rest: delete,
   element and charge, rename and renumber, with the Sequence window's kind rows to pick out
   waters and ions.* Still to do: waters beyond X Å;
   move and rotate a ligand; bonds and orders; collapse alternate locations; split and merge;
   protonation states and flips by hand; mutations from a rotamer library (licence checked)
   with the fewest clashes; Edit in Sketcher with an MCS-constrained re-embed back into the
   pocket.
5. **MD** (OpenMM in its own interpreter). System setup (amber14 or CHARMM36 with GAFF2 or
   OpenFF; explicit or implicit solvent; ions); **restraints from any selection** (positional
   with a force constant, distances between picked atoms, frozen atoms), saved with the setup
   and drawn in the view; minimise, heat, equilibrate and produce, resuming from checkpoints and
   using the GPU when present. Results: the trajectory in the models bar, RMSD, RMSF per residue
   (cartoon coloured by it), ligand RMSD and interaction occupancy as rows for the line plots,
   and clustered representative frames as rows. Trajectories live beside the document.
6. **Binding free energy.** MM-GBSA/PBSA rescoring of poses and trajectories first; then
   relative binding free energy (openfe) with an atom mapper, network setup and convergence
   checks that refuse to report a number that has not converged.
7. **Optional add-ons, if the team asks:** local folding (Boltz or ESMFold: multi-GB weights and
   a GPU), ML charge models to avoid AmberTools' conda-only AM1-BCC, metadynamics.


**Done (2026-11-10): 3D plots on QRhi.** The 3D scatter and landscape plots moved off
pyqtgraph.opengl onto Qt's RHI, as the Structure Viewer's 3D view already was, and PyOpenGL
left the dependencies. Each system now draws with its own graphics API: Direct3D 11 on
Windows (with WARP, its built-in software renderer, behind it, so it draws on any PC,
virtual desktops included), Metal on macOS, and OpenGL on Linux. A start-up check greys 3D out
with the reason where the graphics cannot start.

1. **Finish Settings** (in progress): the Settings window and custom-theme editor.
2. **Jobs and log** (§2–3). Everything after this benefits from non-blocking jobs and a log.
3. **Search and filters** (§4).
4. **Packaging, pilot build** (phase 8 from `PLAN.md`). Putting a real installer in chemists'
   hands now gives feedback before the long tail.
5. **Table power tools:** formula columns, conditional colouring, structure preparation,
   potency standardisation, scaffolds and decomposition (§5, §6 tier 1).
6. **Visual and SAR improvements** (§5), and the 3D viewer (done: ligands and conformers).
7. **Optional add-ons, as the team asks:** databases, external tools, predictions, protein.

## 8. Open questions

- Which tier 2 and 3 features does the team use most? This decides the order of step 7.
- **Branch protection on `main`** (changes only through a pull request with green CI): noted
  2026-11-10, not turned on yet. Revisit with the pilot build.
- Is a **pilot build** after step 3 useful, or should packaging wait for the table power tools?
