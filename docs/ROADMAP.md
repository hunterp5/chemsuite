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
- **Structure Viewer:** sort options, a size slider, and pinning compounds to compare them
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
| Decomposition and enumeration | BRICS, RECAP and most-central-ring decomposition. R-group enumeration from a core and R-group lists, and SMARTS reaction enumeration. |
| Documents | Duplicate document, send selected rows to a new document, and Notes (a notepad saved in the `.csd`). |
| Macros | Record and replay sequences of commands with their settings. The new command registry is the base for this. |
| Random utilities | Random numbers and a random molecule sample. Low priority. |

**Tier 2: external tools or network access (optional, controlled by policy)**

| Feature | Notes |
|---|---|
| Confab and CONFORGE conformers | Through the external-tool paths now in Settings. |
| Superpose conformers, pharmacophore screen | RDKit only. The Structure Viewer's 3D view (Chemistry ▸ Structure Viewer) is done and aligns a molecule's conformers; superposing different molecules is next. |
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

### 7a. Proteins (decided 2026-11-10)

A protein viewer and editor for medicinal chemists, on the native QRhi renderer (no WebEngine).
Proteins are stored in the document; ligands, poses and scores are table rows linked to their
protein by a `Protein` column, so plots, filters and SAR work on them. gemmi is in the base
install; preparation, minimisation and simulation need the optional `chemsuite[protein]` add-on
(OpenMM, PDBFixer, Meeko, propka/pdb2pqr). Each phase is usable on its own.

1. **Viewer and selection** (next). Open PDB/mmCIF or fetch by PDB ID; several proteins per
   document; cartoon, sticks and surface; select chains, residues and atoms in 3D, in the
   sequence strip or with a selection language; binding site with interactions; measurements;
   the table's current row shows its pose in the protein.
2. **Editing.** Delete and rename, protonation states and flips, mutations from a rotamer
   library, ligand edits in the sketcher placed back in the pocket; every edit undoable.
3. **Preparation and minimisation.** PDBFixer, propka protonation at a pH, OpenMM minimisation
   with restraints, in a separate process.
4. **Alignment.** Pairwise sequence alignment and superposition (whole chain or pocket); MSA
   through MAFFT.
5. **Docking and poses.** AutoDock Vina (an external program, all platforms), optional Gnina
   rescoring (Linux, Windows through WSL); poses and scores as table rows; interaction
   fingerprints.
6. **MD simulation.** OpenMM: system setup (force fields for protein and ligand, explicit or
   implicit solvent, ions), **restraints from any selection** (positional with a force
   constant, distances between picked atoms, frozen atoms), a minimise / heat / equilibrate /
   produce protocol that resumes from checkpoints and uses the GPU when present, and results:
   a trajectory player in the protein view, RMSD, RMSF per residue (cartoon coloured by it),
   ligand RMSD and interaction occupancy plots, and representative frames as table rows.
   Trajectories are kept beside the document. Later: MM/GBSA rescoring and relative binding
   free energy (openfe) as a separate add-on.


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
