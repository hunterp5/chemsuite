# chemsuite roadmap after the Data menu (2026-11-05)

`PLAN.md` covered the rebuild of MolManager's core: the table, files, descriptors, plots,
conformers and the Data menu. This roadmap covers what comes next:

- the infrastructure MolManager has and chemsuite still lacks (job queue, process log);
- improvements to the tools chemsuite already has, starting with search;
- the rest of MolManager's features, ranked by value and install cost.

Decisions so far: QSAR (phase 7d) is skipped for now. Groovy and custom themes stay. External
tools (Gnina, AiZynthFinder, BioTransformer…) will be integrated, so their paths are already in
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
dialog, because the close prompt needs to know whether it worked. Not yet: removing ended jobs
one by one in the window, and a desktop notification beyond the dock/taskbar alert.

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
  - the **status bar** ("2 running · 3 queued", with a progress bar for the current job). Click
    it to open the panel.
  - a **Jobs panel** (View ▸ Jobs). It lists running, queued and finished jobs with elapsed
    time. You can cancel a job, remove a queued one, move a queued job up or down, clear
    finished jobs, and copy an error.
  - Results apply when the job finishes. `apply_results` already skips rows deleted in the
    meantime. Each tool's existing checks still run before a job is queued.
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
3. **Quick search bar:**
   - field syntax such as `ID:CPD-12`, `pIC50>7` or `"exact phrase"`;
   - a history dropdown, and a count of matches as you type.
4. **Filters:**
   - a **substructure card**;
   - filters **saved with the document** (a noted follow-up);
   - named filter sets ("Actives, MW < 450") you can switch between;
   - "Copy filtered rows to a new document".
5. **Find and Replace:** searching within the selection or within chosen columns, and
   replacing with regex groups.

## 5. Improving the other existing tools

- **Table**
  - Conditional cell colouring: colour by value or category, like MolManager's
    `column_color.py`.
  - Per-column number formats.
  - Frozen ID and structure columns.
  - Column statistics in the header tooltip.
- **Structures**
  - Align all drawings to a common core (MolManager's "Align to Substructure").
  - Larger hover previews.
  - "Copy as SMILES / image".
- **Plots**
  - Colour-range controls and marker shapes.
  - A density mode (hexbin) for 200k points.
  - A radar plot for MPO profiles.
  - Several Y columns on one plot.
  - Saved plot templates.
- **Selection Browser:** sort and label options, a size slider, and pinning compounds to
  compare them.
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
| Superpose conformers, pharmacophore screen | RDKit only, but needs the 3D viewer to be useful. |
| PubChem, ChEMBL, vendor, patent and SQL queries | Network access, off by default under policy. Plain HTTPS REST calls, not MolManager's client libraries. SQL through the standard library or an optional driver. |
| Retrosynthesis | AiZynthFinder through its own Python interpreter (path in Settings). |
| Metabolite prediction | BioTransformer (Java .jar, path in Settings). |
| Sketcher | **Decided: needed.** A robust drawing sketcher that works with the table: edit a structure cell in place, draw a query for structure search and filters, and add new compounds. MolManager's is about 18k lines, so this is a project of its own and needs its own plan. |

**Tier 3: heavy add-ons (a separate optional install)**

| Feature | Notes |
|---|---|
| pKa and permeability models | torch, unipka and chemprop make the install about 2 GB larger, so this is an optional "predictions" add-on. Settles the open question about the "(est.)" descriptors. |
| Protein viewer, simulation, sequence/MSA, docking (Gnina), interaction fingerprints, pose browser | openmm, pdbfixer, meeko, gemmi, prolif and MAFFT. These need the 3D viewer decision: PyVista as an optional add-on, as discussed. |
| QSAR, Quick QSAR, CoMSIA | Skipped for now (phase 7d). |

---

## 7. Proposed order

1. **Finish Settings** (in progress): the Settings window and custom-theme editor.
2. **Jobs and log** (§2–3). Everything after this benefits from non-blocking jobs and a log.
3. **Search and filters** (§4).
4. **Packaging, pilot build** (phase 8 from `PLAN.md`). Putting a real installer in chemists'
   hands now gives feedback before the long tail.
5. **Table power tools:** formula columns, conditional colouring, structure preparation,
   potency standardisation, scaffolds and decomposition (§5, §6 tier 1).
6. **Visual and SAR improvements** (§5), and the 3D viewer.
7. **Optional add-ons, as the team asks:** databases, external tools, predictions, protein.

## 8. Open questions

- Which tier 2 and 3 features does the team use most? This decides the order of step 7.
- Is a **pilot build** after step 3 useful, or should packaging wait for the table power tools?
