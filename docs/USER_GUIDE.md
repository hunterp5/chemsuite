# chemsuite user guide

chemsuite is a desktop workbench for compound tables. You can open a file of compounds, look at
the structures, calculate properties, search and filter, plot, and run SAR and chemical-space
analyses. All results are written back into the same table.

This guide follows the menus. Keyboard shortcuts are shown for Windows and Linux. On macOS,
read Ctrl as ⌘.

> **Keeping this guide current.** Whenever a tool is added or a feature changes, update this
> guide in the same change.

## Contents

1. [The window](#1-the-window)
2. [Files and documents](#2-files-and-documents)
3. [The table](#3-the-table)
4. [Structures](#4-structures)
5. [Selecting rows and lists](#5-selecting-rows-and-lists)
6. [Searching](#6-searching)
7. [Filters](#7-filters)
8. [Background jobs](#8-background-jobs)
9. [Plots](#9-plots)
10. [Chemistry: descriptors and conformers](#10-chemistry-descriptors-and-conformers)
11. [Data tools](#11-data-tools)
12. [Similarity](#12-similarity)
13. [SAR](#13-sar)
14. [Chemical space](#14-chemical-space)
15. [Settings](#15-settings)
16. [The log and reporting problems](#16-the-log-and-reporting-problems)
17. [Keyboard shortcuts](#17-keyboard-shortcuts)

---

## 1. The window

Each document has its own window. From top to bottom:

- **Menu bar.**
  - File, Edit, Data, Chemistry, Plot, View, Help.
  - Every command that a tool adds appears in these menus. If a command is greyed out, its
    tooltip says why: for example, there is no structure column, or an external program is
    missing.
- **Search bar.** Hidden until you ask for it with Edit ▸ Search Rows (Ctrl+F). See
  [Searching](#6-searching).
- **The table.** One row per compound. Structure columns show drawings.
- **Docked panels.**
  - Filters, Advanced Search and the Selection Browser dock at the side.
  - New plots and tool panels open to the right of the table, splitting the window left and
    right. A second plot joins the first as a tab.
  - Each panel has a header showing its title. Drag the header to move the panel. A plot's
    header also holds the plot toolbar. The buttons at the right are:

    | Button | What it does |
    |---|---|
    | Cog | Opens the panel's settings (plots only). |
    | Window with arrow | Undocks the panel into its own window, or docks it again. |
    | ✕ | Closes the panel. |

- **Status bar**, on one line:
  - At the left: what the current background job is doing, such as "Generating conformers…
    (120/500)" or "Calculating UMAP…", followed by "· 2 queued" or "· 1 failed" when that
    applies. Click it to open the Jobs window.
  - Next to it: the latest message, such as "Opened compounds.sdf" or "Cancelled:
    Generating conformers; partial results kept".
  - At the far right: the row and selection counts ("200 rows · 3 selected").

  View ▸ Status Bar hides or shows it.

---

## 2. Files and documents

### Opening

Use File ▸ Open (Ctrl+O) or File ▸ Open Recent.

| Format | Extensions | Notes |
|---|---|---|
| SD file | `.sdf`, `.sd`, `.mol` | Data fields become columns. A record whose structure does not parse is **kept**: its other fields still load, and the structure cell is marked invalid. Nothing in your file disappears silently. |
| CSV / TSV | `.csv`, `.tsv`, `.tab`, `.txt` | A column headed SMILES, Structure, Mol, Molecule, Canonical_SMILES or Isomeric_SMILES becomes a structure column. |
| SMILES list | `.smi`, `.smiles` | One SMILES per line, with an optional name. |
| Excel | `.xlsx` | The first sheet. |
| chemsuite document | `.csd` | Everything exactly as it was saved (see below). |

Any of these can also be gzip-compressed (`compounds.sdf.gz`). The text encoding is detected
automatically.

Large files open quickly. The text loads first and the table appears; structures are then parsed
in the background, using every processor core. You can scroll and work while that happens.

### Saving

- **File ▸ Save (Ctrl+S)** and **Save As (Ctrl+Shift+S)** write a chemsuite document (`.csd`).
  A document keeps everything:
  - the table: columns, number formats, sort order, and the difference between blank, invalid
    and unparsed structures;
  - named lists;
  - filters and filter sets;
  - the advanced search and recent searches;
  - open plots and panels with their settings.
- Saving is safe. The new file is written next to the old one and swapped in only once it is
  complete, so a crash or a full disk never damages your previous version.
- **File ▸ Export…** writes SDF, CSV/TSV, SMILES or Excel. **Export Selected…** writes only the
  selected rows.

### New and close

- **File ▸ New (Ctrl+N)** opens an empty document in a new window.
- **Close (Ctrl+F4)** closes the window. It asks first if there are unsaved changes or jobs
  still running.

---

## 3. The table

### Cells and editing

- Double-click a cell to edit it, or start typing.
- **Copy (Ctrl+C)** copies the selected cells as tab-separated text, which pastes cleanly into
  Excel.
- **Paste (Ctrl+V)** pastes tab-separated text starting at the current cell, as one step you can
  undo. A molfile on the clipboard (from ChemDraw and similar) pastes into a structure cell.
- **Clear Cells** blanks the selected cells.
- **Undo (Ctrl+Z)** and **Redo (Ctrl+Y)** cover every edit to the table: cells, rows, columns,
  types, and results written by tools.

### Columns

Each column has a type:

| Type | Holds | Behaviour |
|---|---|---|
| Text | Anything | Sorts alphabetically. |
| Number | Plain numbers | Sorts numerically and can be plotted. The number of decimals shown is kept per column. |
| Structure | Molecules | Drawn in the table; used by the chemistry tools. |

Right-click a column header for:

- **Sort Ascending / Descending.** **Then Sort…** adds a second (third…) sort key.
  **Clear Sort** returns rows to their original order (also Data ▸ Table ▸ Clear Sort).
- **Filter This Column.** Adds a filter card for the column (see [Filters](#7-filters)).
- **Rename…**, **Insert Column Left/Right…**, **Delete Column**. Delete acts on all selected
  columns when you right-click one of them.
- **Split Column…** and **Join Columns…**, each with a live preview.
- **Column Type ▸ Text / Number / Structure.** Converting to Number warns you if some values
  are not numbers; those cells would be cleared.

From the Data ▸ Table menu:

| Command | What it does |
|---|---|
| Add Compound… | Opens the structure editor and adds the compound as a new row. |
| Add Rows… | Adds empty rows. |
| Add Column… | Adds a column. |
| Split Column… / Join Columns… | As above. |
| Bin Column… | Groups numbers into ranges and named classes; see [Data tools](#11-data-tools). |

### Rows

- **Delete Rows (Ctrl+Backspace)** deletes the selected rows. You can undo it.
- **Data ▸ Show All Rows** clears every filter, search and "show only" at once.
- **Data ▸ Copy Shown Rows to New Document** puts exactly the rows you can see, in their shown
  order, into a new untitled document. Column types are kept.

---

## 4. Structures

- **Edit ▸ Edit Structure… (Ctrl+E)** opens the structure editor for the current structure
  cell.
  - Type or paste SMILES, a molfile or InChI. The drawing updates as you type.
  - The editor is text-based. Draw in ChemDraw (or similar), copy, and paste here.
- **Edit ▸ Copy Structure as Molfile** puts the current structure on the clipboard for drawing
  programs.
- A structure that cannot be read is shown as **invalid**. Its original text is kept, so you
  can fix it.
- Rows are taller when the table has a structure column. Change the height under Settings ▸
  Appearance.

---

## 5. Selecting rows and lists

### Selection

- Click, Shift-click and Ctrl-click rows in the table.
  - **Select All (Ctrl+A)**, **Invert Selection (Ctrl+I)**, **Clear Selection
    (Ctrl+Shift+A)**.
- The selection is shared by the table, the plots and the tools:
  - selecting points in a plot selects the rows;
  - selecting rows highlights their points.
- Many tools can run on all rows, the shown rows or just the selected rows.

### Selection Browser

View ▸ Selection Browser (Ctrl+Shift+B) steps through rows one at a time with a large drawing.

- It browses either the selected rows or every shown row, in the table's order.
- It keeps the table's current row in step with what it shows.
- Its tick box adds or removes the compound from the selection.

### Named lists

Lists are saved sets of rows, such as "actives", "to make" or "sent for assay". They are saved
with the document.

- **Data ▸ Lists ▸ Add Selected to List… (Ctrl+L)** adds the selected rows to a new or existing
  list.
- Each list has its own submenu:

  | Command | What it does |
  |---|---|
  | Select Rows | Selects the list's rows. |
  | Add to Selection | Adds them to the current selection. |
  | Show Only These Rows | Hides every other row. |

- **Stop Showing Only a List** brings the other rows back.
- **Manage Lists…** renames and deletes lists.
- Deleting a row takes it out of its lists. Undoing the delete puts it back.

---

## 6. Searching

### Quick search (Ctrl+F)

Edit ▸ Search Rows opens the search bar above the table. As you type, the table shows only the
rows that have a matching cell in any column.

- The **Options** button sets **Match Case**, **Whole Cell** or **Regular Expression**.
- Escape clears the search. A second Escape, or the ✕ button, hides the bar.
- Hiding the bar always clears its search, so you never have rows hidden by a search you cannot
  see.

### Find and Replace (Ctrl+H)

A window that steps through matching cells, with **Replace** and **Replace All**. A replacement
is one step you can undo. **Replace All includes hidden rows** controls whether rows hidden by
filters or searches are changed too.

### Advanced Search

Edit ▸ Advanced Search… opens a panel of conditions.

- Each condition picks a column (or **any column**), an operator and a value.
- Click **Not** to negate a condition.
- Operators depend on the column type:
  - **Text:** contains, equals, starts with, ends with, matches a regular expression, is
    blank. Each has a match-case option.
  - **Number:** =, ≠, <, ≤, >, ≥, between (either order), is blank. A blank cell never
    counts as "≠ 6".
  - **Structure:** contains a substructure, is exactly, is blank.
- **Match all / any** combines the conditions.
- A condition with a problem (for example, "between 6 and lots") is marked, and the search
  buttons stay disabled until it is fixed.
- Use the result in one of three ways:
  - **Show** shows only the matching rows. This works together with your filters.
  - **Select** selects the matching rows.
  - **Add to List** adds them to a named list.
- While a search is showing, text and number conditions follow your edits immediately, and
  structure conditions search again when a structure changes.
- The search, whether it is showing, and your recent searches are saved with the document.

### Structure Search (Ctrl+Shift+M)

Edit ▸ Structure Search… searches a structure column.

- **Query.** Draw (paste) or type the query as SMILES, a molfile or InChI. For substructure
  searches you can also type SMARTS, such as `[C,N]`.
- **Modes:**

  | Mode | Finds |
  |---|---|
  | Substructure | Rows that contain the query. |
  | Exact | Rows that are the query, with or without stereochemistry. |
  | Similarity | Rows at or above a similarity threshold, using a fingerprint and metric you pick. |

- **Results.** **Show** them alone (as a structure card in the Filters panel), **Select** them,
  or **Add to List**.
- **Highlighting.** Substructure matches can be highlighted in the table's drawings. Edit ▸ Clear
  Structure Highlighting turns that off.
- **Speed.** Searches are fast even on 200k rows: a fingerprint screen rules out most rows
  before any atom-by-atom matching.

---

## 7. Filters

View ▸ Filters (Ctrl+Shift+F) opens the Filters panel. Each **card** filters on one column, and
a row is shown only if it passes every enabled card.

| Card | For | Controls |
|---|---|---|
| Range | Number columns | A two-handle slider over a histogram of the values, plus min/max boxes. Choose whether blanks are kept. |
| Category | Text columns with few distinct values | Tick boxes, one per value, with counts. |
| Text | Other text | Contains / equals / starts / ends / regular expression, with **Not**. |
| Structure | Structure columns | Substructure, exact or similarity query; searches in the background after you stop typing. |

- Add cards with the panel's **Add Filter** button or a column header's **Filter This Column**.
- Each card has an on/off tick box and a remove button.
- **Filter sets.** The panel's **Sets** menu saves the current cards under a name (such as
  "Actives, MW < 450") so you can switch between groups of filters. Applying a set whose column
  has since been removed restores the cards it can.
- Filters are a way of *looking* at the table. They are not edits, so they are not undone with
  Ctrl+Z and do not mark the document as changed. They are saved with the document, though, and
  the panel opens on load when filters are hiding rows, so you can see why.

---

## 8. Background jobs

Heavy work runs in the background, so the window stays usable. This includes opening big files,
descriptors, conformers, clustering, maps and searches.

- **Two lanes.**
  - Heavy jobs wait in a **queue** and run one at a time, in order.
  - Quick jobs, like a small search, never wait behind them.
- **Status bar.** The job's progress shows at the far left of the status bar: "Generating
  conformers… (120/500)", or just "Calculating UMAP…" when there is no count. Click it to open
  the Jobs window.
- **View ▸ Jobs** lists running, queued and finished jobs with their elapsed time and result.
  From there you can:
  - cancel a job;
  - remove a queued job, or move it up or down;
  - clear finished jobs;
  - copy an error.
- **Cancelling.** A job stops within about a tenth of a second. While it stops, the status bar
  says "Cancelling …".
- **Partial results.** Tools that work row by row keep what they finished when you cancel:
  - the finished rows are written to the table;
  - the message says "partial results kept";
  - rows not reached are left as they were.

  | Partial results kept on cancel | Discarded on cancel |
  |---|---|
  | Calculate Descriptors | Clustering |
  | MedChem plots (descriptor step) | Chemical Space Map |
  | Generate Conformers | Self-Organizing Map |
  | Fingerprint Similarity | Free-Wilson |
  | R-Group Decomposition (numbered core) | Hansch |
  | Diverse Subset (the picks made so far) | Matched Molecular Pairs |
  | | Activity Cliffs |

  The tools in the right-hand column have no meaningful result until the whole run is complete.
- **Safety.**
  - Closing a document cancels its queued jobs and asks about running ones.
  - Quitting with jobs running asks first.
  - Results from a job are never written into a different document.
- **Long runs.**
  - The computer is kept awake while jobs run.
  - A job that takes more than 30 seconds and finishes while you are elsewhere shows a
    notification.

---

## 9. Plots

Open a plot from the **Plot** menu: **New Scatter Plot (Ctrl+Shift+P)**, **New Histogram** or
**New Box Plot**. It opens to the right of the table.

### Naming a plot

The plot's title in its header starts as an automatic title, such as "pIC50 vs MW" or
"Histogram: MW", and follows the columns you plot.

- **Double-click** the title to rename the plot.
- **Enter** keeps the new name. **Escape** puts the old one back.
- A named plot keeps its name when you change its columns.
- Clear the name and press Enter to go back to the automatic title.
- The name is saved with the document, and is the suggested file name when you export.

### Plot settings

Click the **cog** in the plot's header for its settings. Changes apply immediately.

| Setting | Scatter | Histogram | Box plot |
|---|---|---|---|
| Plot | Switch between scatter, histogram and box plot | | |
| X | Number column | Number or text column | **Group by** (an optional text column) |
| Y | Number column | — | **Values** (a number column) |
| Color | Number (gradient) or text (one colour per category, with a legend) | — | — |
| Size | Number column | — | — |
| Bins | — | Number of bins, or Auto | — |
| Curve | — | Overlay a Normal or Log-normal fit | — |
| Guide | BOILED-Egg, Golden Triangle or LipE lines drawn behind the points | — | — |
| Fit line | Least-squares line; one per category when coloured by text | — | — |

The settings popup also shows a summary of what was drawn:

- how many points there are, and how many rows were left out for missing values;
- for fits, the equation and R², or μ and σ for a fitted curve.

### Working with a plot

The toolbar in the plot's header has five buttons:

| Button | Name | What it does |
|---|---|---|
| Four arrows | Pan and Zoom | Drag to pan, scroll or right-drag to zoom, and click a point to select it. |
| Dashed box | Box Select | Drag a rectangle to select points. |
| Lasso loop | Lasso Select | Draw around points to select them. |
| Four corners | Reset View | Fits the plot to the data again. |
| Arrow into a tray | Export | Saves the plot as PNG or SVG. |

- The first three are modes, and the active one is highlighted. In Box and Lasso, Shift adds to
  the selection and Ctrl toggles points in or out of it.
- **Hover** over a point to see its structure and values.
- The plot shows only the rows the table shows, and redraws when the data, filters or columns
  change.
- Plots are saved with the document and reopen with it.

### MedChem plots

The Plot ▸ MedChem menu has three plots. Each one first calculates any properties it needs that
are not already columns, then plots them.

| Plot | Axes | What it shows |
|---|---|---|
| BOILED-Egg | WLOGP against TPSA | The white predicts gut absorption; the yolk, brain penetration (Daina & Zoete 2016). |
| Golden Triangle | MW against logP | The region favoured for permeability and clearance (Johnson et al. 2009). |
| Lipophilic Efficiency (LipE) | Potency (pActivity) against logP | Lines of equal LipE = pActivity − logP. |

---

## 10. Chemistry: descriptors and conformers

### Calculate Descriptors (Ctrl+Shift+D)

Pick a structure column and the properties you want. Each property becomes a number (or text)
column. **Recalculate descriptors that are already columns** controls whether existing results
are overwritten.

| Group | Descriptors |
|---|---|
| Physicochemical | MW, Exact Mass, cLogP, cLogD 7.4 (est.), LogS (ESOL), TPSA, TPSA (N,O,S,P), Molar Refractivity, Labute ASA, Fsp3 |
| Drug-likeness | QED, Ro5 Violations, Ro5 Pass, CNS MPO (est.), BBB Score (est.), AB-MPS (est.) |
| Counts | HBD, HBA, NH/OH Count, N+O Count, Rotatable Bonds, Heavy Atoms, Heteroatoms, Valence Electrons, Net Charge |
| Rings | Rings, Aromatic Rings, Aliphatic Rings, Saturated Rings, Spiro Atoms, Bridgehead Atoms |
| Atom Counts | C, N, O, S, P, F, Cl, Br, I |
| Complexity | Bertz Complexity, Balaban J, Hall-Kier Alpha, Kappa 1–3 |
| Partial Charges | Max/Min Partial Charge, Max/Min Abs Partial Charge |
| 3D Shape | PMI 1–3, NPR1, NPR2, PBF, Radius of Gyration, Asphericity, Eccentricity, Spherocity Index, Inertial Shape Factor |
| Identifiers | Canonical SMILES, InChIKey, InChI, Formula, Murcko Scaffold, Generic Scaffold |

- **"(est.)" values.** chemsuite has no pKa model. cLogD 7.4 is a rough correction of cLogP from
  an estimated basic pKa, and the scores built on it inherit that. Treat them as screening
  estimates.
- **3D shape descriptors** use a 3D structure when the cell has conformers.
- **Cancelling** keeps the rows already calculated.

### Generate Conformers

Generates 3D conformers with RDKit's ETKDG and optimises them with MMFF94 (UFF where MMFF has no
parameters).

- **Settings:** conformers per molecule, a minimum RMSD for two conformers to count as
  distinct, an energy window ("drop conformers more than … kcal/mol above the lowest"), the
  force field, and a random seed.
- **Output.** The conformers go into a new 3D structure column, lowest energy first, with
  columns **Conformers** (the count) and **Conformer Energy Range (kcal/mol)**.
  - Energies compare conformers of the same molecule only.
  - **Replace an existing 3D column of the same name** controls whether a previous result is
    overwritten.
- **Cancelling** keeps the molecules already done.

### Export Conformers to SDF

Writes every conformer as its own SD record, carrying the row's other columns. You can export
all rows or only the selected ones.

---

## 11. Data tools

### Statistics

Data ▸ Statistics… opens a panel with five tabs. It works on the shown rows, all rows, or the
selected rows, and recalculates as the data, filters or selection change.

| Tab | What it gives |
|---|---|
| Summary | N, missing, mean, SD, min, 5%, Q1, median, Q3, 95%, max, skew, kurtosis for every number column. |
| Correlation | Pearson r or Spearman ρ matrix. On big tables it runs as a background job when you ask. |
| Outliers | IQR (Tukey fences), Z-score or modified Z (median and MAD), with a threshold; label rows by a column. |
| Curve Fit | Polynomial (choose the degree), exponential, logarithmic or power fit of Y on X. |
| Tests | Welch t-test, paired t-test, Mann-Whitney U, one-sample t-test, Shapiro-Wilk normality. |

### MPO Score

Data ▸ Scoring ▸ MPO Score… combines several properties into one 0–1 score. Add one row per
property, each with a weight and one of these desirability shapes:

| Shape | Use it when |
|---|---|
| Higher is better | Larger values are always preferred. |
| Lower is better | Smaller values are always preferred. |
| Within a range, best at a target | There is an ideal value inside a range. |
| Near a centre (bell curve) | Values should cluster around a centre. |

Then choose how the desirabilities combine:

- **Geometric mean.** One unacceptable property sinks the compound.
- **Arithmetic mean.** Properties trade off against each other.

Options:

- **Score rows missing a property on the others.**
- **Also add a desirability column per property.**

### Selectivity

Data ▸ Scoring ▸ Selectivity… takes a primary assay and its off-target assays, and works out
log selectivity (1 = 10-fold).

- Say whether the values are p-values (pIC50, pKi…; higher is more potent) or concentrations
  (IC50, Ki… in one unit; lower is more potent).
- It adds a selectivity column per off-target and a **Min logSel** column.
- **Flag a narrow window** adds a **Narrow window** column for compounds below a threshold.

### Bin Column

Data ▸ Table ▸ Bin Column… groups a number column into bins.

- **Methods:** equal width, equal count, or custom edges.
- **Options:** a fixed width, from/to limits, and open-ended bins below the first edge and above
  the last.
- **Output:** a **range column** ("[5, 6)") and a **class column** with names you choose. Low /
  Medium / High is filled in for you.
- A live preview shows each bin's count.

---

## 12. Similarity

All three tools compare compounds either by **fingerprints** or by **number columns** (scaled
Euclidean distance). The fingerprint options are:

| Fingerprint | Notes |
|---|---|
| Morgan r=2 (ECFP4), 2048 bits | The usual default. |
| Morgan r=3 (ECFP6), 2048 bits | |
| Feature Morgan r=2 (FCFP4) | Atoms by role (donor, acceptor, aromatic…); finds scaffold hops. |
| RDKit paths (up to 7 bonds) | |
| Atom pairs | Good for overall shape. |
| Topological torsions | |
| MACCS keys (166) | Small and interpretable. |

Similarity is Tanimoto by default; Dice and Cosine can be chosen instead. Fingerprints are
cached per column, so a second run on the same column is much faster.

### Fingerprint Similarity

Data ▸ Similarity ▸ Fingerprint Similarity… scores every row (or the shown or selected rows)
against one of two references:

- a **structure** you type or paste; or
- the **selected compounds**. Each row gets its similarity to the nearest selected compound, and
  a column naming that compound. A selected row is never its own nearest neighbour.

Options:

- Sort the table, most similar first.
- Select rows with similarity of at least a threshold.

Cancelling keeps the rows already scored.

### Diverse Subset

Data ▸ Similarity ▸ Diverse Subset… picks *N* compounds that cover the chemistry, using MaxMin:
each pick is the compound farthest from all the earlier picks.

- **Pick from** all, shown or selected rows.
- You can start from the current selection, so the new picks complement compounds you have
  already chosen.
- The picks get a rank column, and can be selected.
- Cancelling keeps the picks made so far. Each of them is as diverse as a full run would make
  it.

### Cluster

Data ▸ Similarity ▸ Cluster… groups compounds using one of three methods:

| Method | You set |
|---|---|
| Butina | A similarity cutoff. A higher cutoff gives smaller, tighter clusters. |
| k-means | The number of clusters. |
| Hierarchical | The number of clusters and the linkage. |

- Clusters are numbered by size, so cluster 1 is the largest.
- The result is three columns: the cluster number, **… Size**, and **… Centre**, which marks
  each cluster's representative compound.
- Rows that cannot be clustered (no structure, a missing number) get no cluster.

---

## 13. SAR

### R-Group Decomposition

Data ▸ SAR ▸ R-Group Decomposition… splits a series onto a common core.

- **Core.** Enter it as SMILES or SMARTS, or let the tool suggest one.
  - Label attachment points `[*:1]`, `[*:2]`… to fix which position is R1, R2.
  - Without labels, RDKit numbers every open position it finds.
- **Output.** A column per position (R1, R2…) and a match column.
  - Substituents are labelled by SMILES with `*` at the attachment (`*Cl`, `*OC`).
  - **H** means unsubstituted.
- **Exhaustive matching** is slower; use it for symmetric cores.
- **Cancelling** keeps the rows decomposed so far, but only with a numbered core. An unnumbered
  core must see every compound at once.

### SAR Matrix

Data ▸ SAR ▸ SAR Matrix… pivots activity on two substituent positions (R1 × R2). A third
position makes pages.

- Each cell holds the compounds with that pair of substituents.
- **Estimate empty cells (Free-Wilson)** fills combinations nobody has made with a prediction,
  when the additive model can estimate both groups.

### Free-Wilson

Data ▸ SAR ▸ Free-Wilson… fits activity as a sum of one contribution per substituent.

- **Contributions** tab: each substituent's effect, both as Fujita-Ban (against a reference
  group) and as Free-Wilson (centred).
- **Best Unmade** tab: the highest-scoring combinations not yet made. Set whether higher activity
  is better, and the minimum number of occurrences a substituent needs.
- Substituents that always occur together cannot be separated by the data. They are reported,
  not given arbitrary values.

### Hansch Analysis

Data ▸ SAR ▸ Hansch Analysis… regresses activity on substituent constants:

    activity = a·π + b·π² + ρ·σ + c·(steric) + d

- **Constants.** Published constant tables are not shipped with chemsuite. Point the tool at
  your own CSV with these columns:

      smiles,pi,sigma_m,sigma_p,f,r,es,mr,upsilon,sigma_star

  `smiles` is `H` or an attachment SMILES such as `*Cl`. A blank means the constant is unknown,
  and rows needing it are left out.
- **Terms.** Choose the hydrophobic (π, π²), electronic (Hammett σ, which positions are meta)
  and steric terms.
- **Coefficients** either sum over positions (classical Hansch) or are fitted one per position.

### Matched Molecular Pairs

Data ▸ SAR ▸ Matched Molecular Pairs… finds pairs of compounds that differ by one small change
at one place. H → X changes are included.

- Set the largest change, in heavy atoms, and the minimum number of pairs a transformation needs
  to be listed.
- Transformations (A → B) are summarised over all their pairs, with the activity shift.

### Activity Cliffs (SALI)

Data ▸ SAR ▸ Activity Cliffs (SALI)… lists similar compounds with very different activity.

    SALI = |Δactivity| / (1 − similarity)

- Only pairs at least the **minimum similarity** alike are considered.
- **Skip identical fingerprints** ignores salts and stereoisomers, which would otherwise divide
  by zero.

---

## 14. Chemical space

### Chemical Space Map

Data ▸ Chemical Space ▸ Chemical Space Map… places compounds on a 2D map by fingerprints or
number columns. It writes the coordinates as two columns and opens a scatter plot.

| Method | Options | Notes |
|---|---|---|
| PCA | — | Fast and linear. |
| UMAP | Neighbours, minimum distance | Fits on up to 50,000 compounds, then places the rest. |
| t-SNE | Perplexity | Warns above 50,000 rows. |

- **Same map every run** makes UMAP repeatable, but slower: it runs on one processor core.
- **Colour by** any column.
- The job shows as "Calculating UMAP…" (or PCA / t-SNE) in the status bar.

### Self-Organizing Map

Data ▸ Chemical Space ▸ Self-Organizing Map… trains a grid of cells. Similar compounds land in
the same or neighbouring cells.

- **Size.** Set rows × columns and the number of training passes, or let the size follow the
  number of compounds.
- **Output.** Each compound's cell is written to a column. A panel shows the grid coloured by any
  column, with **Higher is better** for activity.

---

## 15. Settings

Edit ▸ Settings… (Ctrl+,) has four pages. Changes preview immediately in every window. **OK**
keeps them and **Cancel** puts everything back.

| Page | Contents |
|---|---|
| Appearance | Theme: System (follows your computer's light or dark mode), Light and Dark (fixed, whatever the computer uses), Groovy (a random colourful palette) or a custom theme. Also the application and table fonts, text alignment, structure row height, and whether to show the status bar. |
| Keyboard | Change any command's shortcut. |
| External Tools | Paths to optional programs (Open Babel, CONFORGE, Gnina, AmberTools, MAFFT, Java/BioTransformer, AiZynthFinder), each with a **Test** button. Tools that need a missing program are greyed out, and their tooltip says why. |
| General | Processor cores to use (0 = automatic), recent files to keep. |

Your IT department can fix some settings for everyone. Those are shown greyed out with a
padlock.

---

## 16. The log and reporting problems

- **View ▸ Log** shows the session's log. Filter it by level and text.
- **Help ▸ Show Log Folder** opens the folder of log files to attach to a bug report.
- If chemsuite hits an unexpected error, it shows a report you can copy. The details are also
  written to the log.

---

## 17. Keyboard shortcuts

You can change any of these in Settings ▸ Keyboard.

| Action | Shortcut |
|---|---|
| New / Open / Save / Save As | Ctrl+N / Ctrl+O / Ctrl+S / Ctrl+Shift+S |
| Close window | Ctrl+F4 |
| Undo / Redo | Ctrl+Z / Ctrl+Y |
| Copy / Paste | Ctrl+C / Ctrl+V |
| Edit Structure | Ctrl+E |
| Search Rows | Ctrl+F |
| Find and Replace | Ctrl+H |
| Structure Search | Ctrl+Shift+M |
| Select All / Invert / Clear Selection | Ctrl+A / Ctrl+I / Ctrl+Shift+A |
| Delete Rows | Ctrl+Backspace |
| Add Selected to List | Ctrl+L |
| Calculate Descriptors | Ctrl+Shift+D |
| New Scatter Plot | Ctrl+Shift+P |
| Filters | Ctrl+Shift+F |
| Selection Browser | Ctrl+Shift+B |
| Settings | Ctrl+, |
