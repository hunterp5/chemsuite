# chemsuite roadmap after the Data menu (2026-11-05)

`PLAN.md` covered the rebuild of MolManager's core: the table, files, descriptors, plots,
conformers and the Data menu. This roadmap covers what comes next:

- the infrastructure MolManager has and chemsuite still lacks (job queue, process log);
- improvements to the tools chemsuite already has, starting with search;
- the rest of MolManager's features, ranked by value and install cost.

Decisions so far: QSAR, Auto QSAR and CoMSIA (phase 7d), skipped at first, are back on the
roadmap since 2026-11-17 (§6 Tier 1). Metabolite prediction (BioTransformer: rule-based
metabolite structures from an external Java tool) is dropped from the plan; **ML site-of-metabolism
prediction, a different thing, is in since 2026-11-17** (§7.0f), with CYP inhibition after it.
Groovy and custom themes stay. External
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
  - MMP neighbourhood map and transform ledger (from MolManager; in §6 Tier 1).
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
| Scaffolds | Bemis–Murcko and generic scaffolds as columns, plus a scaffold tree or count view (the tree: §7.0h). |
| Decomposition and enumeration | **Done 2026-11-14:** Chemistry ▸ Decompose (R-Group Decomposition, Fragments by BRICS or RECAP, Central Ring) and Chemistry ▸ Enumerate (Reaction Products with 28 named templates, R-Groups on a Core, Recombine Fragments). **2026-11-15:** Prepare Reagents turns reagents into substituents (13 kinds, the leaving group replaced by an attachment point), and R-Groups on a Core takes a position's substituents From a File…, so a prepared set is reused across cores. |
| Documents | Duplicate document, send selected rows to a new document, and Notes (a notepad saved in the `.csd`). |
| Macros | Record and replay sequences of commands with their settings. The new command registry is the base for this. |
| Random utilities | Random numbers and a random molecule sample. Low priority. Random Numbers and Random Rows are done; random molecules from outside the table (ChEMBL, PubChem) are under External data in Tier 2. |
| QSAR, Auto QSAR | Added 2026-11-17. **Planned 2026-11-23, [QSAR_PLAN.md](QSAR_PLAN.md):** not MolManager's `qsar_models.py` / `auto_qsar*.py` port but chemsuite's own engine (`chemsuite/qsar/`) — a large descriptor and fingerprint pool, feature selection and model selection inside nested cross-validation, conformal prediction intervals and an applicability domain, regression and classification. Fitted with scikit-learn, applied in numpy; a `.csqm` model archive saves in the document and in a file, applies in another session and refines when rows are added. No new dependency; no chemprop and no torch (QSARtuna is the design reference only). **Built 2026-11-24:** the engine (`chemsuite/qsar/`, 18 modules) and the tool (Data ▸ QSAR ▸ Auto QSAR…, Apply Model…, Refine Model…, and the QSAR Models panel). Every zoo family is exported to flat arrays and checked against scikit-learn to 1e-9, so a saved model applies with numpy alone — a test applies one in a fresh interpreter with `sklearn` blocked at the import hook. Models are kept as document assets, so they save, undo and reopen with the document. |
| CoMSIA | Added 2026-11-17. MolManager's `comsia.py` (780 lines): 3D-QSAR fields on aligned conformers. Builds on Superpose for the alignment. After QSAR. |
| Sites of metabolism, CYP inhibition | Added 2026-11-17, see §7.0f and [SOM_PLAN.md](SOM_PLAN.md). Our own model: a score per atom per P450 isoform, then CYP inhibition on the same package. A gradient boosting machine over 34 RDKit atom descriptors and an isoform one-hot (decided 2026-11-18, after measuring it level with a graph network that cost 8.5M parameters). No external tool, no new dependency — scikit-learn is already one — and the trees are vendored in-tree at about 100 KB, so it works offline with no download. |
| Structural alerts | Added 2026-11-17. RDKit's FilterCatalog sets (PAINS, Brenk, NIH, ZINC and the rest) as columns or a filter: which alerts a row hits and how many. MolManager's `chem/structural_alerts.py` (213 lines) was marked to copy and never was. Probably a tab or group in Calculate Descriptors, as in MolManager. |
| MMP neighbourhood map and transform ledger | Added 2026-11-17. MolManager's `mmp_neighborhood_analysis.py` and its map, plot and ledger views: a compound's matched-pair neighbours around it, and every transform with its count and mean activity change. Builds on Matched Molecular Pairs. **2026-11-18:** the Matched Pairs panel draws transformations and pairs (the change highlighted), shows each transformation's Δ spread, and lists the current compound's transformations (Current Compound): a first step toward the neighborhood view. |

**Tier 2: external tools or network access (optional, controlled by policy)**

| Feature | Notes |
|---|---|
| Confab conformers | Through the Open Babel path in Settings. (CONFORGE is in: a Generator choice in Generate Conformers, from the optional `cdpkit` package.) |
| Superpose conformers, pharmacophore screen | **Superpose done (2026-11-16):** Chemistry ▸ 3D ▸ Superpose (`chemsuite/superpose`) fits rows onto a row or a protein's ligand by MCS or a core, O3A, shape, shape + features (RDKit's `rdShapeAlign`) or field (our own Gaussian ESP overlap), keeps each row's best conformer and its scores, and the Browser shows the reference. **Unified Similarity Search done (2026-11-17,** [SIMILARITY_PLAN.md](SIMILARITY_PLAN.md)**):** Chemistry ▸ Similarity ▸ Similarity Search… (`chemsuite/similarity`) scores the table or screens a library file against one or more references by fingerprints (Tversky too), the Gobbi 2D pharmacophore, MCS and scaffolds, and 3D shape, shape + features and field on the superpose engine, with a consensus score, a fingerprint prefilter for the slow methods, and the aligned 3D poses kept with their reference. **Next:** 3D pharmacophore point alignment (perhaps CDPKit's `PharmacophoreAlignment`), ErG and count fingerprints. |
| 3D pharmacophore screen | Added 2026-11-17. MolManager's Pharmacophore Screen: build a 3D pharmacophore (feature points with tolerances, excluded volumes) from a reference row, several aligned actives or a protein's pocket (§7a 2), then screen the table or a library file by point alignment over each row's conformers, keeping the fit score and the aligned pose. A method in Similarity Search beside shape and field, perhaps on CDPKit's `PharmacophoreAlignment`. |
| External data: PubChem, ChEMBL, vendors, SureChEMBL patents, Open Reaction Database, SQL | Added in detail 2026-11-17. Network access, off by default under policy. Plain HTTPS REST calls, not MolManager's client libraries. Look up names, properties and identifiers for rows; bring in ChEMBL activities for a target; vendor availability; patents that contain a structure; reactions from ORD; random molecules from ChEMBL or PubChem as a new document. SQL through the standard library or an optional driver. The user's details never go into a request. |
| Retrosynthesis | AiZynthFinder through its own Python interpreter (path in Settings). |
| Sketcher | **Decided: needed.** A robust drawing sketcher that works with the table: edit a structure cell in place, draw a query for structure search and filters, and add new compounds. Plan: [SKETCHER_PLAN.md](SKETCHER_PLAN.md) (native, Ketcher rules as reference, verified; no MolManager iupac_* code). |

**Tier 3: heavy add-ons (a separate optional install)**

| Feature | Notes |
|---|---|
| Permeability models | Predictions that belong in Chemistry ▸ Calculate beside Predict pKa. chemprop and torch would make the install about 2 GB larger, so as they stand they are an optional add-on; ported to numpy, as chemsuite pKa was, they could join the base install. Settles the open question about the "(est.)" descriptors. pKa itself is done: chemsuite pKa is in the base install (Chemistry ▸ Calculate ▸ pKa, Standardize, the descriptors and Prepare; [UNIPKA_PLAN.md](UNIPKA_PLAN.md) U5, 2026-11-18). |
| Protein viewer and editor, sequence alignment, docking, preparation, MD simulation | Planned in phases: see **§7a Proteins** below. Still missing from MolManager's macromolecule tools (added 2026-11-17): sequence alignment and MSA (§7a 1), PLIF across poses, pocket surfaces and protein pharmacophores (§7a 2), MD (§7a 5) and free energy, MolManager's FEP (§7a 6). |

---

## 7. Proposed order

### 7.0 Done on Windows, owed a Mac check: graceful exit (added 2026-11-11, built 2026-11-15)

Closing the last document window shuts chemsuite down completely, whatever is in flight:
- **The prompt** counts every job the user started, from every window, and names the
  background tasks that stop with them; background work alone stops without asking
  (`app/shell/quitting.py`). A window closed while others stay open asks only about its own
  document.
- **Windows:** once the last document window has closed, every other top-level window (Log,
  Sketcher, tool windows, error reports) closes and the event loop ends.
- **Work:** `JobRunner.shutdown` cancels every job, gives the threads 3 s in all, then ends the
  process pools' workers whatever they are doing (`parallel.stop_pool`: terminate, then kill)
  and every external program.
- **External programs** (Vina, the AiZynthFinder worker, the keep-awake helper) start through
  `core/jobs/children.start`: a process group of their own on macOS and Linux (SIGTERM, then
  SIGKILL), a kill-on-close Job Object on Windows, so whatever they start ends with them, even
  if chemsuite crashes. Programs chemsuite opens for the user (a browser, a file viewer) are
  left alone.
- **The last resort:** if a job's thread is stuck in a call that never looks at its cancel
  token, `main.run_app` logs a forced exit, writes the audit `exit` record and ends the process
  without waiting (saves are atomic, so nothing is half-written).

`tests/app/shell/test_exit.py` starts a stuck job, a busy pool worker, an external program with
a child of its own, background work, the Log and a tool window, closes the last window and
checks that the app and every process it started are gone. It passes on Windows, and on
2026-11-19 a scripted run docked eight ligands with the bundled Vina (four `vina.exe` at once):
closing the window mid-run ended the app in 0.4 s with no Vina left, and killing the app
outright ended them too (the Job Object).

Still to do, on macOS only: run `tests/core/test_children.py`, `tests/core/test_parallel.py`
and `tests/app/shell/test_exit.py` (the process-group path has not run yet); by hand, dock
several ligands, close the window mid-run and check Activity Monitor or `ps`; and check that
⌘Q from the application menu goes through the same prompt and shutdown.

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
- **Status (2026-11-14):** CONFORGE is wired in as an optional Generator (`conformers/conforge.py`,
  the `cdpkit` wheel; plan approved 2026-11-14). ETKDG stays the default until the benchmark
  (`scripts/review/conformer_bench.py`, `docs/conformers/bench.md`) is built and run. Wang et
  al. (JCIM 2023, 63, 6525) found ETKDG within 0.03–0.07 Å of OMEGA and ConfGenX at 1–250
  conformers, and force-field clean-up no help for bioactive poses, so the benchmark also scores
  ETKDG without MMFF, without the energy window, and at larger counts.
- **Bench run (2026-11-15, `docs/conformers/bench.md`):** CONFORGE is the most accurate
  (AMR-R 0.837 / 0.582 Å at 10 / 250 conformers, against 0.985 / 0.725 Å for the ETKDG
  default). The default's energy window and MMFF clean-up make ETKDG worse at 50 and 250;
  plain ETKDG gives 0.747 / 0.623 Å there. CONFORGE is now the default Generator, with
  `cdpkit` in the base install where it has wheels. Where ETKDG is the default (no cdpkit), its energy window starts off.
- **Shipped (checked 2026-11-19).** CONFORGE is the default behind the same tool, and ETKDG
  stays available. It keeps partial results on cancel and runs in the process pool. The user
  guide states the benchmark numbers.

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
- **Measured and fixed (2026-11-17,** [conformers/scale.md](conformers/scale.md)**):** at
  1k × 300 (190 conformers on average) every step the user waits on meets its target once
  three fixes were in: the molecule cache is limited by size (memory over the blobs 434 →
  238 MB), Overlay works out the bond graph once (559 → 51 ms), and exporting every
  conformer writes one molfile per molecule (114 → 37 s). **Done 2026-11-19:** 2D-only work
  on a 3D column loaded every conformer in the workers, 6–30× slower than on the 2D column;
  it now reads a cell's SMILES when the blob is over 4 KB (`StructureColumn.flat_inputs`,
  `ToolContext.flat_structure_inputs`): screen fingerprints 3.0 → 0.52 s, substructure search
  416 → 40 ms, MolWt + TPSA 1.7 → 0.14 s at 1k × 300. Parity tests found and fixed InChI
  reading undefined stereo from the first conformer. Conformers stay in the cells.

### 7.0d High priority: TMAP and Chemical Space Networks (added 2026-11-11)

**Done 2026-11-14:**
- both plots are built (Data ▸ Plot ▸ Chemical Space Network and TMAP), on chemsuite's own
  engine in `chemsuite/analysis/graphs` (no `tmap` dependency);
- the user guide has a section for each.

**Done 2026-11-19** ([chemspace/bench.md](chemspace/bench.md)):
- **MHFP6 as a fingerprint choice** for TMAP (ECFP4 stays the default; the network stays on
  ECFP4): RDKit's encoder, 256 MinHash values put straight into the forest
  (`analysis/graphs/mhfp.py`). Encoding is its cost: about a minute for 50,000 compounds.
- **Above 49k**, on ChEMBL 37: 250,000 compounds in 52 s (before the layout change below;
  about 63 s since) and 1.2 GB on ECFP4, 82 s (plus 152 s of encoding) and 1.5 GB on MHFP6;
  the forest's candidates went to int32 to get there.
- **Beside `tmap-silicon`** on the same MHFP6 values: chemsuite finds more of the true ten
  nearest (0.93 against 0.62–0.72), joins closer compounds and is 2–3× faster.

**Done 2026-11-24: the layout** ([chemspace/bench.md](chemspace/bench.md)). tmap-silicon's
(OGDF FM³) spread the tree with each branch apart; chemsuite's packed branches into clumps and
crossed them. Repulsion now reaches 16 edge lengths on the coarse levels and 8 on the finest
(it was 3), and each level starts its nodes toward their neighbors. On the 49k MHFP6 tree the
picture covers as much of the page as tmap-silicon's (0.31), crossings halved (0.059 to 0.026
per edge; tmap-silicon 0.020); layout takes about twice as long (27 s at 250,000).
Still open: at 250,000 the small separate trees are packed in a thin strip beside the giant one.

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

### 7.0h More chemical-space views: GTM, a scaffold tree, PCoA (added 2026-11-24)

Already built: PCA, UMAP, t-SNE, PaCMAP and LocalMAP (Data ▸ Dimensionality Reduction), SOM,
TMAP, Chemical Space Networks, the Constellation Plot, SALI and the 3D Landscape. Each item
below shows something none of those does. GTM and the scaffold tree are big new features, so
each starts in plan mode.

1. **Generative Topographic Mapping (GTM).** Bishop, C. M., Svensén, M. and Williams, C. K. I.
   (1998). "GTM: The generative topographic mapping." *Neural Comput.* 10, 215–234; used as
   "chemography" on ChEMBL-scale sets by Varnek's group (Gaspar et al., *J. Chem. Inf. Model.*
   2015, 55, 84). A probabilistic cousin of the
   SOM: a grid of nodes bent into the descriptor space by EM, where each compound gets a
   responsibility over the nodes rather than one cell. That gives what no current view gives:
   - smooth property and activity landscapes (the responsibility-weighted mean per node), and
     class landscapes for a category column;
   - a per-compound log-likelihood as an applicability-domain measure, which Auto QSAR could
     reuse;
   - new compounds (another library, a design set) projected onto a saved map without
     refitting, so two sets can be compared on the same map.

   chemsuite's own code in numpy (EM with an RBF basis, as in the paper); no `ugtm` or other
   dependency. Input like the SOM's: fingerprints reduced by PCA, or number columns. Saved
   maps follow the `.csqm` rule: arrays plus a JSON manifest, never a pickle. The plan should
   settle the grid size, the RBF count and width, the regularization, and how far it scales
   (EM costs rows × nodes per step; 100k rows on a 25 × 25 grid is the first target).

   **Map and landscapes built 2026-11-25** (`analysis/gtm.py`, `tools/space/gtm_ui.py`, the
   GTM Landscape plot kind in `app/plots/gtm.py`). Defaults: 10 / 16 / 25 nodes per side below
   300 / 3,000 / above, basis functions max(3, nodes/4) per side, width 1.0 × their spacing
   (Netlab), regularization 0.1 (ugtm), at most 100 EM iterations, stopping when the penalized
   log-likelihood per row and dimension gains < 1e-4. The map is a document asset (`gtm_map`,
   an `.npz` + JSON manifest, format 1) holding each row's largest responsibilities (99.9% of
   the mass, at most 64 nodes). Measured on the Windows PC (25 × 25, 50 dims): 10k rows 4 s,
   100k rows 61 s (17 iterations at 3.1 s, then 7.5 s placing), ~300 MB numpy peak; 10k
   BindingDB compounds on ECFP4, 26–73 s with the fingerprints and PCA. On those 10k (ECFP4,
   regularization 0.1) each compound keeps 9.6 nodes on average and 28% sit almost wholly
   (≥ 0.99) on one; 177 of the 625 nodes are some compound's best node. (Synthetic tight
   clusters put nearly every row on one node.)
   **Next:** projecting new compounds onto a saved map (format 2 must add the PCA mean and
   axes, the column scaling and the fingerprint weight), and the log-likelihood as Auto
   QSAR's applicability domain.
2. **Scaffold tree view** (the "scaffold tree" of §6's Scaffolds row). Schuffenhauer, A. et al.
   (2007). "The Scaffold Tree − visualization of the scaffold universe by hierarchical scaffold
   classification." *J. Chem. Inf. Model.* 47, 47–58, as Scaffold Hunter draws it. Each
   Bemis–Murcko scaffold is pruned one ring at a time by the paper's rules down to a single
   ring, and the scaffolds form a tree: which ring systems a set holds and how they nest. It is
   a hierarchy of substructures rather than of similarity, so TMAP and the maps cannot show it.
   - Each node shows its scaffold, its compound count and a summary of a chosen column (mean
     pIC50, say); clicking a node selects its compounds in the table, like any plot.
   - The pruning rules are reimplemented from the paper with tests (RDKit has no scaffold-tree
     builder; Scaffold Hunter's Java code is a reference only).
   - Builds on the scaffold columns that already exist (`descriptors/scaffold.py`).
3. **Principal coordinates (PCoA, classical MDS) on Tanimoto distance**, low priority. Gower,
   J. C. (1966). *Biometrika* 53, 325–338. A linear map whose axes come from the true
   Tanimoto distances, not from PCA on the bits, with the variance each axis explains. It
   needs the full distance matrix (n²), so it is for sets up to a few thousand compounds; a
   fourth method in Data ▸ Dimensionality Reduction beside PCA rather than a tool of its own.

**Considered and not planned:**
- **TriMap** (Amid and Warmuth, 2019): PaCMAP was designed from an analysis of TriMap and UMAP
  and matches or beats it on global and local structure in its own benchmarks (Wang et al.,
  *J. Mach. Learn. Res.* 2021, 22, 1–73), and it is already built. TriMap would add the
  `trimap` package (annoy, numba) for a map that mostly repeats PaCMAP's, with weaker local
  detail.
- **PHATE**: made for continuous trajectories in single-cell data, rarely useful on compound
  sets.
- **Isomap and LLE**: superseded by UMAP and PaCMAP.

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
- **Surface options, done 2026-11-17** ([SURFACE_PLAN.md](SURFACE_PLAN.md)): Display ▸ Surface
  gives SES (the new default; older documents open as Gaussian), SAS, van der Waals and Gaussian
  surfaces with probe, quality and smoothing; Solid, Mesh and Dots; color by any scheme, one
  color or electrostatics (Amber charges by the hydrogens present, Gasteiger for ligands, the
  potential at eps 4r); transparency per surface and for the pocket skin; carving around the
  selection. SES of about 6,000 heavy atoms takes 3 s and its potential 5 s, both in the
  background.
- **Polished 2026-11-19:** faded structure between the eye and a focused ligand is cut away
  (a ghost-only near plane at the ligand's front, following the camera), so ribbons no longer
  cross it; a kept row whose protein was removed shows its Protein cell in red with a tooltip.
- **Workspace:** the viewer stays its own window with its own menus. Its Protein menu in the main
  window grows as the phases land (align, dock, prepare, simulate), each opening or driving the
  viewer.

### 7.0g High priority, parked: minimise after Prepare on an OpenMM worker (added 2026-11-15)

Prepared complexes should be simulation-ready: Prepare's last step becomes a restrained
minimisation in OpenMM's own interpreter, and that worker and its `SystemSpec` are the base for
MD, MM-GBSA/PBSA, solvent analysis and free energy (§7a 3, 5, 6). Plan approved 2026-11-15 in
[SIMULATE_PLAN.md](SIMULATE_PLAN.md): ff14SB + OpenFF Sage 2.3 (AshGC charges, no AmberTools) +
TIP3P by default, ff19SB + OPC offered; no in-app MMFF tier; the minimised coordinates are a new
model beside the prepared one; on by default when the worker is set; a native env with an
optional WSL interpreter. M0 (the feasibility gates) passed on 2026-11-15
([simulate/m0-review.md](simulate/m0-review.md)); M1 (the Prepare fixes and
`chemsuite/simulate`'s spec, build and result) and M2 (the OpenMM worker, its runner and the
shared JSON-lines client) done the same day, as are M3 (the minimisation protocol and its
measured result) and M4 (Minimise in Prepare and Protein ▸ Minimise…, each result a new
protein). **M5, the benchmark,** splits: M5a is the minimization metrics on the
P0 sets, **done 2026-11-24** (167 of 262 neutron entries minimized; RMSD median 0.26 Å, clash
score 28 → 2 per 1000 atoms, neutron marks unchanged; 19 entries fail for defects, listed in the
plan, which are next); M5b is the plan's short explicit-solvent MD gate,
which needs solvation, ions, a thermostat and a production loop that `chemsuite/simulate` does
not have yet, so it is planned with §7a 5. Prepare's Minimise stays off by default until M5b.
A by-hand Prepare + Minimise run on 181L with the real env is still owed (see HANDOFF).
Loops built by Prepare ([LOOPS_PLAN.md](LOOPS_PLAN.md)) are relaxed here: their atoms are in
`Prepared.built`, so they move freely while the deposited atoms are held.

### 7.0e High priority: our own Uni-pKa (revised 2026-11-12)

Uni-pKa gives the free energy of each protonation microstate, so linked sites and tautomers
are handled together; MolGpKa treats each site on its own, and in protein preparation it
catches few ligand states (`docs/prepare/`). Its published weights cannot be used (licence),
so chemsuite trains its own on pKaHub's data from Uni-Mol's MIT backbone and runs them in numpy, with
no torch: [UNIPKA_PLAN.md](UNIPKA_PLAN.md). **Done 2026-11-18:** ligand states in preparation
(U4), then Predict pKa, Standardize and the descriptors (U5); MolGpKa is removed.

### 7.0f Sites of metabolism, then CYP inhibition (added 2026-11-17)

Nothing in chemsuite says anything about metabolic fate. Sites of metabolism comes first: a score
for every atom of a molecule for each of nine human P450 isoforms. FAME 3 and MetaQSAR are
academic-only, so chemsuite trains its own model on openly licensed data and ships the result.

**The model is a gradient boosting machine over 34 cheap RDKit atom descriptors plus a nine-wide
isoform one-hot** (decided 2026-11-18). The route there is worth recording: the plan began on
ATTNSOM's graph network (Kim et al., *Bioinformatics* 2026; MIT), which was trained and reproduced,
until the independent `som-benchmark` (MIT) was run under matched labels and showed a plain atom
feature GBM level with it — 0.601 against 0.618 top-1 on a ±0.047 seed spread, with the remaining
biases favoring the network. A tie on accuracy made the 8.5M-parameter numpy port, its `SAGPooling`
reimplementation and 17 MB of vendored weights impossible to justify against scikit-learn, which is
already a base dependency, and a model measured in kilobytes. The Hückel descriptors that first
motivated the switch were measured and dropped too. Full record in [SOM_PLAN.md](SOM_PLAN.md).

**The engine is in (2026-11-18).** `chemsuite/metabolism` holds the features, the trees walked in
numpy and the batch driver, with a 250 KB `som_gbm.npz` vendored in-tree — no download, no consent
prompt, and no scikit-learn at run time (the export is refused unless the numpy path reproduces
scikit-learn's probabilities, measured at 1.1e-16). Molecule-grouped over ten seeds it reaches
top-1 0.592, top-2 0.727 and top-3 0.813, which is the tie the decision assumed.

**Done (2026-11-18).** Chemistry ▸ Calculate ▸ Sites of Metabolism… writes the ranked sites, the
top site, its score and an applicability value per isoform, and the sites are marked on the
structure column's drawings — the best site strongest, weaker scores fainter. That last part is an
app feature rather than the tool's: a structure column can mark the atoms any column *names*
(`C7`, or `C7 0.81, N3 0.05`), through **Highlight Atoms ▸** in its header menu, with each name
checked against the molecule beside it. CYP inhibition (phase C) is next.

Accuracy is quoted on molecule-grouped splits against `PRIMARY_SOM` alone: counting the secondary
and tertiary annotation tiers as hits would add about 0.10 top-1 for nothing, so every figure in
the user guide names its label definition. CYP inhibition follows on the same package, trained on
OpenADMET's open pIC50 data and the Veith panel from PubChem.

**This is not the metabolite prediction that was dropped.** That item was BioTransformer, a
rule-based enumerator of metabolite *structures* run as an external Java tool. This is an ML model
that scores atoms on a graph: no external tool, no new runtime dependency, and no generated
structures — metabolite structure generation stays out of scope.

### 7a. Proteins: from viewer to medicinal-chemistry workspace (revised 2026-11-11)

Proteins are stored in the document; ligands are table rows linked to their protein by a
`Protein` column, so plots, filters and SAR work on them. **Results worked out in the Protein
Viewer stay in the viewer** (decided 2026-11-15, reversing "poses as table rows"): the
preparation report, the Inspect review and the docking runs dock beside the pocket in the
viewer's window and are kept with the document; poses reach the table only through **Copy to
Table**. Each phase is usable on its own and starts in plan mode.

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
   process per ligand with progress, cancel and finished rows kept; poses kept in the viewer's
   Poses panel, every run side by side (pose rank, Vina score, CNN score, parent), shown in the
   pocket by the pose display and copied to table rows on request; covalent and
   flexible-side-chain docking; ensemble and consensus scoring; pocket detection and
   druggability, with the pocket's surface drawn and coloured by property (MolManager's Pocket
   Surface); interaction fingerprints across poses, building on Interactions to Table;
   **protein pharmacophores** (added 2026-11-17, from MolManager): feature points and excluded
   volumes from a pocket and its ligand's interactions, or from the empty site, drawn in the
   viewer and used by the 3D pharmacophore screen (§6 Tier 2). The
   receptor is prepared by `macromolecule/prepare` (phase 3's P1 is enough for docking).
   - **Progress (2026-11-14):** Stage A is Prepare (phase 3's P4). Stage B: AutoDock Vina 1.2.7
     bundled for Windows, macOS and Linux with its licence and pinned hashes
     (`macromolecule/dock/bin`, `scripts/fetch_vina.py`); the Qt-free engine `macromolecule/dock/` (box,
     Meeko PDBQT, the Vina binary, many ligands with progress, cancel and finished ligands
     kept); **Protein ▸ Dock…** and **Rescore Poses…** (`tools/macromolecule/dock.py`,
     `dock_dialog.py`). Redocking 80U into 5GNK puts the top pose 1.13 Å from the crystal's.
     **2026-11-15:** poses moved out of the table into the viewer's Poses panel (Copy to Table
     on request); every run kept and numbered with its scoring function, a rescore a run
     beside the dock it scored, runs removable. Not yet tried by hand in the app. Next: Gnina,
     ensembles, covalent and flexible docking, pocket detection, PLIF across poses.
3. **Native preparation** (base install, high priority; plan approved 2026-11-11 in
   `docs/PREPARE_PLAN.md`).
   - **Progress:** P0 (benchmark) and P1 (repair, hydrogens, flips, His states, rotors;
     `pipeline.prepare`) done 2026-11-12, at Reduce's accuracy or better
     (`docs/prepare/p1-review.md`). P2 (titration: Asp/Glu/His/Lys/Tyr/Cys/N-terminus states
     at a pH, pKa values and populations) done 2026-11-12 at PROPKA's accuracy on PKAD-R, off
     by default (`docs/prepare/p2-review.md`). P3 (ligand protomers and tautomers, pocket
     waters, metal coordination) done 2026-11-13. **P4, the UI, done 2026-11-14: Protein ▸
     Prepare… and its report panel**, applied as one undoable edit and saved with the document.
     **2026-11-15: Protein ▸ Inspect…**, the same preparation laid out to review before it is
     applied: each repair can be left undone, each site's state, alternate location and
     crystal water overruled (`Decisions`, pinned states re-solved around), Re-survey, Apply as
     one edit, and the report records what was chosen by hand. The ⚠ rule is one function
     (`pipeline.warns`) for the report and the review. A one-click "use the other choice" on a
     report line was planned and dropped: Inspect is the route.
     **2026-11-19: titration judged on states, not pKa values**
     ([`prepare/p5-review.md`](prepare/p5-review.md)). On the neutron set, with a benchmark of
     states (`scripts/review/prepare_states.py`), titration still does not beat titration off:
     0.884 against 0.893 on the test half. That holds even with PROPKA-style determinants and a
     per-type threshold a group must pass to leave its form in water. So titration **stays off
     by default**. A titrating His now keeps P1's hard clash on its neutral states.
     **Kept options**, each to be run through that benchmark as a benchmark-only tool before
     anything is built:
     - Poisson–Boltzmann titration;
     - an ML protein-pKa model;
     - constant-pH MD for the sites marked ⚠ (in the OpenMM worker).

     Also next: the receptor export for docking (phase 2).
   - **Repair:** missing heavy atoms and caps from the CCD templates; a residue missing alone
     built (done 2026-11-15). **Loops of 2 to N residues closed** by CCD on Ramachandran-sampled
     backbones, with side chains from the templates and relaxed by Minimize (plan approved
     2026-11-19 in [LOOPS_PLAN.md](LOOPS_PLAN.md): L0 a benchmark of cut loops against
     PDBFixer, L1 the Qt-free closure, L2 Prepare, Inspect and the dialog, L3 measured before
     it is on by default; N set by L3, about 8–12). Longer gaps and missing terminal tails stay
     capped and marked ⚠.
   - **Hydrogens and states, protein and ligand together:** one coupled site graph covers His
     tautomers, Asn/Gln/His flips, titratable residues, OH/SH/NH₃ rotors, pocket waters, metal
     sites and the ligand's protomers and tautomers. Ligand states come from chemsuite pKa (U4).
     Desolvation, Coulomb and hydrogen-bond terms act across the interface, so each side's state
     depends on the other's. It is solved exactly per cluster, with populations; near-ties are
     marked ⚠ and become receptor ensembles for docking or separate runs for MD.
   - **Measured first:** benchmarked against neutron structures, PKAD-2, Reduce, PROPKA3,
     Protoss (public PDB IDs through ProteinsPlus) and, optionally, MOE Protonate3D before it
     becomes the default. propka stays as a cross-check.
   - **Relax after preparation** is now §7.0g, planned in [SIMULATE_PLAN.md](SIMULATE_PLAN.md):
     a staged, restrained minimisation of the whole complex in OpenMM's own interpreter (ff14SB +
     Sage 2.3 in GBn2), on by default when the worker is set; no in-app MMFF tier.
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
   **Metal sites need a real force-field model** (added 2026-11-24): minimization holds each
   transition metal's contacts at the distances Prepare left them and keeps its other
   neighbors from closing in, because a bare nonbonded Zn²⁺ pulls an anion and its neighbors
   in (the carbonic anhydrase sulfonamides, SIMULATE_PLAN M5a). Restraints are enough for a
   minimization but not for MD: a bonded model (MCPB.py-style parameters for each site) or a
   cationic dummy-atom model, chosen per metal, measured on the P0 metal entries.
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
