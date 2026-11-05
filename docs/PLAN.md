# chemsuite — build plan

Status: approved 2026-11-04. Phases 0–2 (setup, table core, files and documents) done 2026-11-04. Decisions: GPLv3 licence, UMAP in the base install, lightweight structure editor, 100k–200k-row tables.

## Context

chemsuite is a desktop medicinal-chemistry workbench for a team of about 8 chemists, deployed by corporate IT
on Windows, macOS and Linux. It replaces MolManager
(`~/Applications/MolManager/mctoolkit`), a first attempt that grew to about 210k lines across about 40 features. Its
chemistry code is solid and mostly free of Qt. Its application wiring is not: about 60 mixins sit on one
window, dialogs reach into the window's private members (`app._x`) 271 times, and there are about 2,000
`getattr`/`hasattr` fallbacks and 1,088 bare `except Exception` blocks. Embedding Chromium through
QtWebEngine to draw Plotly plots pins it to Qt 5.15.2, makes a basic install 1–2 GB, and caused the failed
PySide6 port.

**Scope:** the full compound table, plotting, descriptor calculation, conformer generation, and everything
under MolManager's **Data** menu.

**Quality bar:** the best possible table and wiring, portable across operating systems, memory-efficient,
and small to install.

**Strategy:** start a new codebase. Copy in MolManager's Qt-free chemistry and analysis modules, with their
tests. Rebuild the table core, the application wiring, plotting and packaging.

---

## 1. Technology choices

| Concern | Choice | Why |
|---|---|---|
| Language | Python 3.12 | RDKit and PySide6 both support it. Without torch, NumPy 2 is fine. |
| GUI | **PySide6 (Qt 6), `PySide6-Essentials` only** | LGPL licence (easier for IT and legal than GPL PyQt5). Qt 6 is still supported. With no WebEngine, the Chromium problems behind the reverted port go away. |
| Plots | **pyqtgraph** (about 2 MB, pure Python) | Fast scatter, histograms and images with 100k+ points, drawn natively. No browser engine needed. |
| Chemistry | RDKit | Already the engine in MolManager. |
| Numerics | numpy, scipy, scikit-learn | The QSAR, Hansch, Free-Wilson and clustering code needs them. |
| Dimensionality reduction | umap-learn (**base install**; pulls in numba/llvmlite) | Required by the team. Imported only when a UMAP run starts, so it costs nothing at startup. |
| Dropped | QtWebEngine, plotly, pandas, sqlalchemy, orjson, psutil, cryptography, PubChem/ChEMBL client libraries | Size, attack surface, or replaceable with the standard library. |
| Excel | openpyxl (**base install**, decided 2026-11-04) | Excel is how most collaborators exchange data. Imported only when an .xlsx file is opened or saved. |
| Optional plugins | openbabel (Confab), cdpkit (CONFORGE) | Each tool says which packages it needs. If one is missing, its menu item is disabled and the tooltip says why. |
| Storage | SQLite from the standard library (molecule store and document files) | One file per document, works the same on every OS. |
| Packaging | PyInstaller one-folder build, Inno Setup (Windows), `.app` + DMG (macOS), AppImage (Linux) | MolManager's per-user install, SBOM, hashed lockfile and signing hooks carry over. |

**Size budget, enforced in CI:** installer ≤ 200 MB and unpacked ≤ 550 MB on Windows (UMAP's numba/llvmlite adds roughly 100–150 MB unpacked). Measure on the first
build and tighten from there. The heavy parts are RDKit, Qt Essentials, scipy plus sklearn, and numba/llvmlite. Leaving out
every Qt module that isn't used is part of the budget.

---

## 2. Architecture

```
chemsuite/
├── core/                 # Qt-free. Must never import PySide6; a CI test enforces this.
│   ├── table/            # Table, columns, view index, filters, sort, commands (undo)
│   ├── mols.py           # structure text <-> RDKit <-> stored binary
│   ├── io/               # SDF/SMILES/CSV/TSV/xlsx readers and streaming writers
│   ├── document/         # .csd document file (SQLite): save and open
│   ├── jobs/             # Job protocol, progress, cancellation, process-pool helpers
│   └── settings.py       # typed settings + IT policy file + offline switch
├── chem/ analysis/ descriptors/ conformers/ plotting/   # ported from MolManager (Qt-free)
├── app/                  # Qt shell
│   ├── main_window.py    # menus from the tool registry, docks, status bar; knows no tools
│   ├── workspace.py      # one per open document; owns table, mols, selection, undo, jobs
│   ├── table_view/       # QAbstractTableModel adapter, view, delegates, render pipeline
│   ├── plots/            # pyqtgraph plot widgets, lasso, selection linking
│   ├── panels/           # filter panel, search, find/replace, lists, browser
│   └── job_runner.py     # runs core jobs on threads/processes, delivers results on GUI thread
├── tools/                # one folder per menu feature
│   └── descriptors/      #   tool.py (registration) · dialog.py · compute.py
└── tests/
```

### 2.1 The table core (`core/table`): the most important part

The core is one Qt-free `Table` object. Everything else goes through its public methods.

- **Rows** are identified by a stable `oid` (`int64`). A row's index is only its current position. oids
  never change when the table is sorted, filtered or edited.
- **Typed columns**, all sharing one interface (`take`, `put`, `insert`, `delete`, `sort_key`,
  `match_text`, `values`). *(As built in phase 1.)*
  - `NumberColumn`: `float64` values, a validity mask, and an `int8` display code **per cell**, about
    10 bytes per row. Each cell keeps how it was written, so "1.50" stays "1.50" and "00123" stays
    text, while a column mixing `1.5` and `2.25` is still numeric. (MolManager kept such columns as
    text.) Integers are a number column whose cells use the integer code, so there is no separate
    int type.
  - `TextColumn`: dictionary-encoded, with `uint32` codes into an append-only list of unique values.
    Repeated values cost 4 bytes per row, and text filters and sort ranks run once per unique
    value. Undo snapshots share the dictionary instead of copying strings.
  - *Deferred:* a packed-UTF-8 mode for mostly-unique text, a boolean type, and structure and
    conformer columns (structure columns were built in phase 2; conformers come in phase 6). Packed UTF-8 waits until a benchmark shows it is
    needed: 200k unique IDs currently take about 20 MB.
- **Types** are inferred when a file loads and can be changed by the user. No column sits at a hard-coded
  position, unlike MolManager's 0 = ID and 1 = Structure.
- **One source of truth.** There is no SQLite mirror of cells and no `SMILES` column duplicating structures.
  A structure's SMILES is derived from the stored molecule when needed.
- **View index.** The table keeps one `int64` array mapping display positions to rows, produced as
  filter → sort.
  - Sorting is vectorised: `np.argsort(kind="stable")` on floats, sorting by code rank on dictionary-encoded
    text, and `np.lexsort` for multi-column sorts. Blank values always go last.
  - Filtering is vectorised: each filter returns a boolean mask, and the masks are combined with AND.
    Substructure filters are the exception; they run as jobs and their mask is cached for each SMARTS
    pattern.
  - The underlying data is never physically reordered.
- **Change events.** The core emits plain events such as `CellsChanged(oids, cols)`,
  `RowsInserted`, `RowsRemoved`, `ColumnsChanged` and `ViewChanged`. The Qt adapter turns these into the
  model signals Qt expects. The core itself never depends on Qt.
- **Commands give undo for free.** Every change is a command object (`SetCells`, `InsertRows`,
  `RemoveRows`, `AddColumns`, `RemoveColumn`, `RenameColumn`, `ConvertColumnType`). Each one stores only
  the data it touched: the old values of the cells it changed, or the rows it deleted (including their
  molecule blobs). A small Qt-free undo stack applies and reverts them. The UI's undo menu just calls it, so
  undo never needs a reference to the window.
- **Reading API for tools:** `table.read(cols, oids=None) -> dict[str, np.ndarray]`, `table.selected()`,
  `table.visible_oids()`. Tools never touch columns directly. This takes the place of MolManager's
  `TableReader` and the `app._table_model` reaches.
- **Writing API for tools:** `table.apply(ResultColumns(...))` is the single way results come back into the
  table. It carries over MolManager's naming rule from `on_calc_finished`: rename to `Name (2)` on a clash,
  but reuse a column whose cells are all blank or errored.

### 2.2 Structures and rendering

- **Structures are a column kind** (`StructureColumn`). *(Changed during phase 2; this replaces the
  SQLite `MolStore`.)* Each cell holds RDKit binary (about 200–600 bytes, with the file's
  coordinates) plus its display text. The same undo, copy and save paths as every other column
  therefore cover structures too, a table can have several structure columns (compound, core,
  R-group), and there is no SQLite connection to share between threads. At 200k rows a structure
  column takes about 100 MB. Text formats load structures unparsed and parse them in parallel
  afterwards. Invalid structures are kept and flagged, never dropped. A bounded cache holds decoded
  molecules.
- **Depictions** reuse MolManager's resolution-independent scene format and painter
  (`chem/depiction_*.py`, `table/depiction_painter.py`), plus `PriorityRenderPlan`, which renders the rows
  on screen first.
- **Rendering is lazy:** only rows on screen plus a margin, never the whole table. Scenes are stored in
  SQLite. Pixmaps go into **one** cache with a byte limit, keyed by `(oid, column, size)`. That replaces
  MolManager's three separate caches, one of which has no limit.

### 2.3 How the pieces connect

- **The `Workspace`** is created for each open document. It owns the table, the molecule store, the
  selection, the named lists, the undo stack and the job runner. The main window holds workspaces and docks.
  It contains no tool code.
- **Each tool is self-contained:**
  ```python
  # tools/descriptors/tool.py
  TOOL = Tool(
      id="descriptors", menu=("Chemistry", "Calculate"), title="Descriptors…",
      requires=("rdkit",), needs_selection=False,
      open=lambda ctx: DescriptorDialog(ctx).show(),
  )
  ```
  The dialog receives a narrow **`ToolContext`**. It can read the table, see the selection, submit a job,
  apply results and open a plot. It cannot reach the window or any private state, so the 271 `app._x`
  accesses in MolManager have nothing to reach.
- **The registry** builds every menu from the `TOOL` declarations. It also checks `requires`, disabling a
  tool and explaining why when a package is missing. Adding a tool never means editing the main window.
- **Jobs** are plain functions: `run(inputs, progress, cancel) -> Result`. `JobRunner` chooses where each
  one runs:
  - a thread pool for light work;
  - a spawn-context process pool for CPU-heavy work (descriptors, fingerprints, clustering, MMP,
    conformers);
  - and a single queue for jobs that must run one at a time.

  Progress is published through a lock-protected state object (from MolManager's `tool_progress.py`) that
  the GUI polls. Results are always applied on the GUI thread, and always through `table.apply`.
- **Selection** is one `SelectionStore` of oids, ported from MolManager. The table, plots and browser all
  read and write it. Each change records where it came from, so a view does not redraw for a change it made
  itself.
- **Errors:** jobs raise typed exceptions. The runner logs them and shows a single error dialog. Bare
  `except Exception` is banned outside the top-level handlers, and ruff enforces this from the first commit.

### 2.4 Plotting

- pyqtgraph `PlotWidget` subclasses: scatter (colour, size and shape mapped to columns), histogram, box,
  bar, heatmap (for the SAR matrix) and radar (for MPO).
- **Selection linking:**
  - click, box or **lasso** in a plot selects compounds in the SelectionStore;
  - selecting compounds anywhere highlights the matching points.

  pyqtgraph has no built-in lasso, so it is drawn as a path and tested with a vectorised numpy
  point-in-polygon check. This piece needs about 150 lines.
- **Hover** shows a structure tooltip using the same depiction painter as the table.
- The Qt-free `plotting/*.py` series and statistics code is ported. The parts that build Plotly figures are
  replaced by plain data `PlotSpec` objects that the widgets draw.
- **Export:** PNG and SVG come from pyqtgraph's own exporters.

### 2.5 Files, documents, settings

- **Readers:** port the gzip-aware readers in `table/table_file_formats.py`, the delimiter sniffing in
  `table/text_file_ingest.py`, and the parallel byte-range SDF parsing in `chem/ingest_records.py`. Rows
  load in chunks straight into the typed columns.
- **Writers** stream rows from the table core to the file (SDF, CSV, SMILES, xlsx). There is no full in-memory snapshot first.
- **The document file (`.csd`)** is a single SQLite file: a JSON manifest, one row per column holding its
  packed arrays, and the molecule, scene and ensemble tables, copied with `ATTACH` the way MolManager's
  `session_v4.py` does. **No encryption by default.** Disk encryption covers that, and leaving it out also
  drops the `cryptography` dependency.
- **Settings** are a typed dataclass stored in a JSON file under the per-user config directory from
  `QStandardPaths`. The **IT policy file** and **offline-by-default** network switch carry over from
  `platform_support/session_runtime.py`, along with the audit log.

### 2.6 Memory and performance targets (CI benchmarks)

| Operation | Target at 200k rows × 50 columns |
|---|---|
| Load 200k-row SDF | < 30 s, with the window responsive throughout (measured: 9 s on 8 cores) |
| Sort on a numeric column / a text column | < 50 ms / < 150 ms |
| Change a numeric range filter | < 30 ms |
| Scroll with structures shown | 60 fps once cached; placeholders while cells render |
| Memory with the table open | < 2 GB, structure images within their byte limit |
| Startup to first window | < 2 s; heavy imports (sklearn, scipy) wait until a tool needs them |

---

## 3. What is reused from MolManager

**Copied with tests, light renames:**
- `analysis/`: statistics, QSAR, Quick QSAR, CoMSIA, MPO, SALI, selectivity, MMP, Free-Wilson, Hansch, SAR
  matrix, medchem space, dimensionality reduction, linear fit.
- `descriptors/`: drop the import from `medchem_descriptors.py` to `ionization`, because it pulls in torch.
- `conformers/`: the Confab, CONFORGE and storage codecs.
- `chem/`: `rdkit_fingerprints.py`, `smarts_macropatterns.py`, `structural_alerts.py`, the depiction
  modules and `ingest_records.py`.
- `table/`: `column_split.py`, `column_join.py`, `column_binning.py`.
- `ui/search_query.py`: it has no Qt despite its location. Also `ui/row_lists.py` and `ui/selection_store.py`.

**Algorithms to separate from their Qt worker classes** (they are already module-level functions):
- `workers/chemistry_descriptors.py`
- `workers/fingerprint_similarity.py`
- `workers/diverse_subset_worker.py`
- `workers/cluster_worker.py`
- `workers/conformer_generation.py` (ETKDG)

**Rebuilt:** the table model and its mixins, the filter proxy, the SQLite cell mirror, undo commands, every
window mixin, every dialog (thin rewrites against `ToolContext`), the Plotly/WebEngine plotting layer, the
QRunnable workers and their signals, and the selection browser.

**Left out:** docking, FEP and MD simulation, protein prep, retrosynthesis, metabolite prediction, pKa and
permeability ML, the external database queries, macros, MolManager's 10k-line custom sketcher, the WSL launcher, and
session encryption. Any of these can come back later as tools without changing the core.

> **Structure editing (decided):** a lightweight editor: edit SMILES or a molfile with a live depiction
> preview and validation, plus paste from ChemDraw and other sketchers through the clipboard
> (MOL/SMILES). It is built in phase 3.

---

## 4. Phases

Every phase ends with tests passing on Windows, macOS and Linux in CI.

0. **Set up the project.** `git init`, `pyproject.toml`, the lockfile, ruff (including the ban on bare
   `except`), pytest, the CI matrix on three operating systems, and an import-boundary test that fails if
   `core/` imports Qt. Add the GPLv3 licence and the header check.
1. **Build the table core.** Typed columns, the view index, sort and filter, change events, commands and
   undo, the read and write APIs. Add property-based tests and the benchmark suite. *No Qt yet.*
2. **Structures and files.** Structure columns, the readers and streaming writers, and the `.csd` document file.
   Round-trip tests run against MolManager's `samples/`.
3. **Qt shell.** Main window, workspace, the table model adapter, the view and delegates, structure
   rendering, selection, copy and paste, undo, find and replace, search, the filter panel, lists, the
   browser, and the Data → Table operations (add row or column, split, join).
4. **Jobs and the tool registry.** Port **Descriptors** end to end as the template every later tool
   follows.
5. **Plotting.** pyqtgraph widgets, the lasso, selection linking, hover structures, the Plotter tool.
6. **Conformers.** ETKDG in the core build; Confab and CONFORGE as optional plugins.
7. **The rest of the Data menu,** in batches:
   - (a) Statistics, MPO, Binning, Selectivity Ratio;
   - (b) Fingerprint Similarity, Diverse Subset, Cluster;
   - (c) SAR Matrix, MMP, SALI, Free-Wilson, Hansch;
   - (d) QSAR, Quick QSAR, CoMSIA;
   - (e) MedChem plots, PCA, t-SNE, SOM and UMAP.
8. **Packaging and IT handover.** Per-OS installers, the size check, the SBOM, signing hooks, the policy
   file documentation, an architecture document for IT, and a user guide.

---

## 5. Decisions (2026-11-04)

1. **Licence:** GPLv3, the same as MolManager. Each source file carries a GPL header, checked in CI.
2. **Data size:** 100k–200k rows. All performance targets are set at 200k.
3. **Structure editing:** the lightweight editor (§3).
4. **UMAP and Excel (.xlsx):** in the base install. Confab and CONFORGE stay optional.

## 6. Verification

- `pytest`: unit and property tests for the core, ported domain tests, and Qt tests through `pytest-qt`
  with `QT_QPA_PLATFORM=offscreen`.
- `scripts/bench.py`: the §2.6 targets, run in CI. A build fails if it regresses by more than 20%.
- `scripts/check_size.py`: fails CI if the frozen build exceeds the §1 size budget.
- A manual smoke test on each OS for every release: open the sample SDF, sort, filter, calculate
  descriptors, lasso in a plot, generate conformers, save, reopen, undo.
