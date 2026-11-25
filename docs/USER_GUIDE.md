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
10. [Chemistry: descriptors, pKa, standardizing, conformers and retrosynthesis](#10-chemistry-descriptors-pka-standardizing-conformers-and-retrosynthesis)
11. [Table and analysis tools](#11-table-and-analysis-tools)
12. [Similarity](#12-similarity)
13. [SAR and QSAR](#13-sar-and-qsar)
14. [Chemical space](#14-chemical-space)
15. [Settings](#15-settings)
16. [The log and reporting problems](#16-the-log-and-reporting-problems)
17. [Keyboard shortcuts](#17-keyboard-shortcuts)

---

## 1. The window

Each document has its own window. From top to bottom:

- **Menu bar.**

  | Menu | Holds |
  |---|---|
  | File | New, New from Shown Rows, Open, Append, Merge, Save, Export ▸ Table, Visible Rows, Selected Rows or Conformers (SDF), Close (clears the window to a new, empty document; with other windows open it closes this one) |
  | Edit | Undo, Copy and Paste, selecting rows (Select All, Invert, Clear, Delete Selection) |
  | Tools | Random ▸ (Random Numbers and Random Rows) and Standardize Values; Calculator and Notepad; Sketcher |
  | Chemistry | Standardize; Calculate ▸ (Descriptors, pKa) and 3D ▸ (Generate Conformers, Superpose); Similarity ▸ (Similarity Search, Diverse Subset), Reaction ▸ (Split Reactions, Retrosynthesis), Enumerate ▸ (Reaction Products, R-Groups on a Core, Recombine Fragments), Decompose ▸ (R-Group Decomposition, Fragments, Central Ring), Routes |
  | Data | Table ▸ and Lists ▸; Scoring ▸, Cluster and Bin Values; Plot ▸, Dimensionality Reduction ▸, MedChem ▸, Constellation Plot, TMAP and Chemical Space Network; SAR ▸; Search ▸ and Filters; Browser |
  | Data ▸ Table | Operations ▸ (Add Compound, Add Rows, Add Column, Split Column, Join Columns), Statistics |
  | Data ▸ Search | Search Rows, Find and Replace, Advanced Search, Structure and Reaction Search, Clear Highlighting, Show All Rows |
  | Data ▸ Dimensionality Reduction | PCA, UMAP, t-SNE, PaCMAP; Self-Organizing Map, Generative Topographic Map |
  | Data ▸ Plot | Scatter, Line, Bar Chart, Histogram, Box Plot, Density, Radar, 3D Scatter, 3D Landscape |
  | Data ▸ MedChem | BOILED-Egg, Golden Triangle, Lipophilic Efficiency (LipE) |
  | Macromolecule | A button, not a menu: opens the Macromolecule Viewer, the workspace for proteins, with its own menus (on macOS, a menu holding Macromolecule Viewer) |
  | Settings | Theme, Hotkeys, External and General (each its own dialog), Status Bar, User Guide (this guide), Log Folder, About |

  - Every command that a tool adds appears in these menus. If a command is greyed out, its
    tooltip says why: for example, there is no structure column, or an external program is
    missing.
- **Search bar.** Hidden until you ask for it with Data ▸ Search ▸ Search Rows (Ctrl+F). See
  [Searching](#6-searching).
- **The table.** One row per compound. Structure columns show drawings.
- **Docked panels.**
  - Filters sits at the right edge at a fixed width. Advanced Search opens as its own
    floating window.
  - New plots and tool panels open to the right of the table, splitting the window left and
    right. A second plot joins the first as a tab.
  - Panels sit flush against the table, divided by a thin line you can drag to resize them
    (the grab area is a few pixels wider than the line).
    Every panel (plots, Filters, Advanced Search, Browser) has the same header
    showing its title. Drag the header to move the panel. A plot's
    header also holds the plot toolbar. The buttons at the right are:

    | Button | What it does |
    |---|---|
    | Cog | Opens the panel's settings (plots only). |
    | Split square | **Arrange**: put the panel **Beside** or **Below** another panel (splitting that panel's space in two), **Tab With** another panel, dock it **Left of**, **Right of** or **Below the Table**, or **Float** it. |
    | Window with arrow | Undocks the panel into its own window, or docks it again. |
    | ✕ | Closes the panel. |

  - **Tabbed panels.** Panels that share a space show their tabs in the header, where the
    title would be, so no extra row is taken below the plot. The current panel's tab is
    raised. Click another tab to bring that panel forward. Drag the current tab (or the
    header's empty space) to move the panel. Double-click the current tab to rename a plot or a note.
    When the header is narrow, long tab titles are shortened with "…"; hover over a tab to
    see its full title.
  - Plots leave a little space between the axis titles and the panel's edges.
  - **Arranging panels.** For stacked plots beside the table, open a second plot (it joins
    the first as a tab), then choose Arrange ▸ Below ▸ *the first plot*. You can also drag a
    panel by its header onto the top, bottom or side of another panel to split its space, or
    onto its middle to tab them.
  - **Layout**, at the right end of the menu bar (**Window ▸ Layout** on macOS), lays out every
    docked panel at once. Each choice keeps the other two and applies straight away:
    - how: **As Tabs**, **Side by Side**, **Stacked** or as a **Grid**;
    - where: **Right of the Table**, **Left of the Table** or **Below the Table**;
    - how much room: **Table Focus** (the panels take a third of the window), **Balanced**
      (half) or **Pane Focus** (two thirds).

    **Dock Floating Panes** brings undocked panels back into the arrangement. **Reset Layout**
    docks every panel as tabs right of the table, at half the window. The menu's icons show
    each layout, and the last choice is saved with the document. **Larger Structures**,
    **Smaller Structures** and **Reset Structure Size**, at the bottom, size the table's
    drawings (see [Structure size](#structure-size)).
  - **The Filters panel stays put.** It always sits at the right edge of the window at a fixed
    width, whatever you do with the other panels: new panels, arrangements and moves all go
    to its left. It has no Arrange or Undock button and cannot be dragged elsewhere.

- **Status bar**, on one line:
  - At the left: what the current background job is doing, such as "Generating conformers…
    (120/500)" or "Calculating UMAP…", followed by "· 2 queued" or "· 1 failed" when that
    applies. Click it to see the jobs in the Log window.
  - Next to it: the latest message, such as "Opened compounds.sdf" or "Cancelled:
    Generating conformers; partial results kept". When a tool finishes, its message says what
    it did, past tense first, with the counts and the column it wrote ("Generated 1,234
    conformers for 100 rows in “Structure 3D”"); anything that went wrong follows a
    semicolon. It stays for about eight seconds.
  - At the far right: the row and selection counts ("200 rows · 3 selected").

  Settings ▸ Status Bar hides or shows it.

### Tool dialogs

Tool dialogs do not block the window. While one is open you can still select rows, change
filters and scroll the table, and the dialog's row counts follow along. Opening another
document closes any open tool dialogs.

Every tool dialog uses one layout, and so do Split Column, Join Columns, Precision, Align
Drawings and Structure Search. The window title names the tool, and every dialog reads the same
way, top to bottom:

- **Input.** What to work on: structures or number columns, and **Rows**.
  - Rows is **All · Shown · Selected**, with each option's count. The counts update as you
    select rows or change filters.
  - If you choose **Selected** with nothing selected, the dialog says so above its buttons. Select rows in the
    table and the dialog is ready to run.
- **Method.** The tool's settings. Settings that only apply to one method appear when you
  choose that method. Hover over a field or its label for an explanation.
- **More Options.** Optional settings whose defaults suit most runs sit in a folded section.
  Click its arrow or title to open it.
- **Output.** The name of the result column. A tool never replaces an existing column: when
  the name is taken, the result goes to a new column named "Cluster (2)", "Cluster (3)" and so
  on, so earlier results are kept.
- **Problems.** If something has to be fixed first, the problem appears above the buttons,
  marked ⚠, and the main button is greyed out. Otherwise there is no text there.
- **Main button.** It names the action and how many rows it covers, such as **Cluster 1,204
  Rows**.
- **Guide.** The button at the bottom left opens the tool's section of this guide.

---

## 2. Files and documents

### Opening

Use File ▸ Open (Ctrl+O) or File ▸ Open Recent, or drag files from your file manager onto an
empty window. Dropping several opens the first in that window and each of the others in a
window of its own; a window that already holds a table does not take dropped files.

| Format | Extensions | Notes |
|---|---|---|
| SD file or molfile | `.sdf`, `.sd`, `.mol` | A molfile (V2000 or V3000) opens as a one-row table. Data fields become columns. A record whose structure does not parse is **kept**: its other fields still load, and the structure cell is marked invalid. Nothing in your file disappears silently. |
| Mol2 file | `.mol2` | Tripos mol2, as docking and modelling programs write it: each molecule becomes a row with its 3D coordinates, and its name goes in a **Title** column. As with SD files, a record that does not parse is kept and marked invalid. A mol2 file without hydrogens may read with the wrong charges, since the format leaves them to be worked out from the hydrogens. |
| CSV / TSV | `.csv`, `.tsv`, `.tab`, `.txt` | A column headed SMILES, Structure, Mol, Molecule, Canonical_SMILES or Isomeric_SMILES becomes a structure column. A column whose first 20 filled cells are all reaction SMILES (such as USPTO's `rxn_smiles`) becomes a [reaction column](#reactions-in-the-table), even if its name ends in "smiles". |
| SMILES list | `.smi`, `.smiles` | One SMILES per line, with an optional name. |
| InChI list | `.inchi` | One InChI per line, with an optional name after a space or tab (it goes in a **Title** column). A line that does not parse is kept and marked invalid. |
| ChemDraw | `.cdxml`, `.cdx` | Each molecule in the drawing becomes a row. Text, arrows and other drawing objects are left out. |
| Marvin | `.mrv` | One molecule, as a one-row table. |
| Maestro | `.mae`, `.maegz` | Schrödinger's format: each structure becomes a row, its title goes in a **Title** column and its other properties become columns. |
| XYZ | `.xyz` | Each frame becomes a row, its comment line the **Title**. XYZ has no bonds, so they are worked out from the 3D coordinates, assuming a neutral molecule; a frame where that fails is kept and marked invalid. |
| Excel | `.xlsx` | The first sheet. Structure and reaction columns are found as for CSV. |
| RD file | `.rdf`, `.rd` | Reactions with their data: a **Reaction** column, and one column per data field (`$DTYPE`). Registry numbers go in an **ID** column. As with SD files, a record whose reaction does not parse is kept and marked invalid. Molecule records (`$MFMT`) go in a **Structure** column. |
| RXN file | `.rxn` | One reaction, as a one-row table. |
| chemsuite document | `.csd` | Everything exactly as it was saved (see below). |

**Hydrogens.** In SD, mol2 and Maestro files, a structure with 3D coordinates keeps the
hydrogens the file gives it (all of them, or only the polar ones of a docked pose), since their
positions are part of the structure and 3D descriptors need them. A 2D structure is read
without them, as drawn. The cell shows the SMILES without hydrogens either way.

Any of these except Excel and binary ChemDraw (`.cdx`) can also be gzip-compressed
(`compounds.sdf.gz`). The text encoding is detected
automatically (UTF-8 or Windows-1252). A character that does not fit the detected encoding is
read as Windows-1252, so a µ, ° or Å is kept wherever it sits in the file, in every text
format.

Large files open quickly. The text loads first and the table appears; structures are then parsed
in the background, using every processor core. You can scroll and work while that happens.

### Append and Merge

Both bring another file (any format above, or a chemsuite document) into the open table as
one step that Undo takes back. The table's own columns, filters and panels stay as they are.

- **File ▸ Append…** adds every row of the file at the bottom of the table and selects them.
  Columns are matched by name; a column the table does not have is added, and the table's
  other columns are blank in the new rows.
- **File ▸ Merge…** matches the file's rows to the table's and fills in what the file
  knows. Choose the table column to **Match** and the file column to match it **With** (a
  column named like ID, Name or Compound is chosen at first). Structures match however they
  were written (as canonical SMILES); other values match as the text shown.
  - **Values in the table:** **Replace them** (the file's value wins) or **Keep them, fill
    blanks**. A blank cell in the file never overwrites a value.
  - **Rows with no match:** **Add them as new rows**, or leave them out.
  - The file's other columns are added; a matched row gets their values. When the file has
    the same key twice, its first row is matched and the rest count as rows with no match.
  - Text merged into a number column (such as `>10`) turns that column into text.
- Proteins a chemsuite document's rows link to come too, unless the table already holds one
  with the same id.

### Saving

- **File ▸ Save (Ctrl+S)** and **Save As (Ctrl+Shift+S)** write a chemsuite document (`.csd`).
  A document keeps everything:
  - the table: columns, number formats, sort order, and the difference between blank, invalid
    and unparsed structures;
  - named lists;
  - filters and filter sets;
  - the advanced search and recent searches;
  - open plots and panels with their settings. Analysis panels (Free-Wilson, Hansch, Matched
    Pairs, Activity Cliffs) rerun their analysis on the same rows when the document opens
    (they show "Recalculating…" meanwhile); Hansch keeps its substituent constants, so the CSV
    file is not needed again;
  - the selected rows, the current cell and where the table was scrolled to;
  - the workspace layout: where each panel is docked (beside, stacked, tabbed, left, right or
    below the table, floating), the size of each pane (where you dragged the dividers
    between the table and the panes, and between panes), which panes are open (Filters,
    Advanced Search, the panels), column widths, the structure row height, and the search bar with its text and
    options.
- A document reopens exactly as you left it: a pane that was closed stays closed, even when
  its filters or search are hiding rows (the status bar still says how many rows are shown).
- Saving is safe. The new file is written next to the old one and swapped in only once it is
  complete, so a crash or a full disk never damages your previous version.
- **File ▸ Export ▸ Table…** writes SDF, mol2, CSV/TSV, text (.txt), SMILES, Excel, Maestro,
  XYZ, an InChI list, ChemDraw (CDXML), Marvin (MRV), an RD file or an RXN file.
  **File ▸ Export ▸ Visible Rows…** writes only the rows the filters and searches leave
  shown, in the order the table shows them; **File ▸ Export ▸ Selected Rows…** writes only
  the selected rows.
  - An **RD file** holds the first reaction column as RXN blocks and every other column as
    data. Reactions read from SMILES have no drawing, so they are laid out as the sketcher lays
    them out (using every processor core). Blank values are left out.
  - An **RXN file** holds one reaction: export a single selected row.
  - A **text file (.txt)** or a **SMILES file (.smi)** first asks for its **Delimiter** (tab,
    comma, semicolon or pipe; a SMILES file can also use a space) and whether to write a
    **Header** line. A text file holds every column, with a header by default. A SMILES file
    holds one **Structure** column and, after it, the **Names** column you pick (Name, Title
    or ID when the table has one, or none), with no header by default; its header is
    `SMILES` and the names column's name. A field that contains a comma, semicolon or pipe
    delimiter is quoted; after a tab or space, the rest of a SMILES line is the name. The
    choices are remembered for the next export. Opening a SMILES file reads any of these
    delimiters back.
  - CSV, TSV and Excel write reactions as reaction SMILES.
  - **Excel** gets numbers as numbers and everything else as text: a cell that starts with
    `=` stays text and never becomes a formula.
  - An **SD file** ends a data field at a blank line, so blank lines inside a cell are left
    out; a line starting `$$$$` (the end of a record) is indented by one space, and `<`, `>`
    in a column name become `(`, `)`. An RD file indents a line of a cell that starts with
    `$` in the same way.
  - A **mol2 file** holds the first structure column, one molecule per row (blank and invalid
    structures are left out), with every hydrogen, SYBYL atom types and Gasteiger partial
    charges. Molecules with 3D coordinates keep them; others get 2D ones. Only the **Title**
    column is written (as each molecule's name): mol2 has no data fields.
  - **Maestro, XYZ and InChI** files hold the first structure column, one molecule per row
    (blank and invalid structures are left out), named by the **Title** column. In a Maestro
    file every other column is kept as a property (numbers as numbers). XYZ needs 3D
    coordinates, so generate conformers first; every hydrogen is written.
  - A **ChemDraw (CDXML)** or **Marvin (MRV)** file holds one molecule: export a single
    selected row.

### New and close

- **File ▸ New (Ctrl+N)** opens an empty document in a new window.
- **Close (Ctrl+F4)** closes the window. It asks first if there are unsaved changes or jobs
  still running.

---

## 3. The table

### Cells and editing

- Double-click a cell to edit it, or start typing. Double-clicking a structure opens the
  sketcher (see [Structures](#4-structures)).
- Right-click a cell for **Browse**, Copy, Paste and Clear Cells. **Browse** opens the
  [Browser](#browser) at that row (in 2D, or as it is when it is already open). Structure cells
  offer **Edit in Sketcher…** and **Copy As ▸ Molfile / SMILES / InChI / InChIKey / Image** instead of Copy.
  Copy As takes every selected structure in that column, in the order shown: text formats one
  per line, molfiles as an SD file. A single molfile also pastes into ChemDraw as a drawing.
  Reaction cells offer **Edit in Sketcher…** and **Copy As ▸ RXN File / Reaction SMILES /
  Image** (see [Reactions in the table](#reactions-in-the-table)).
- **Copy (Ctrl+C)** copies the selected cells as tab-separated text, which pastes cleanly into
  Excel.
- **Paste (Ctrl+V)** pastes tab-separated text starting at the current cell, as one step you can
  undo. A molfile on the clipboard (from ChemDraw and similar) pastes into a structure cell, and
  an RXN file into a reaction cell.
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

**Selecting columns.** Click a column header to select the whole column; Shift-click selects a
range and Ctrl-click (⌘-click) adds or removes one. Selected columns are not a selection of
compounds, so plots and the Browser are left alone. **Delete** removes the selected
columns in one step you can undo. Sorting is in the header's right-click menu.

**Moving columns.** Drag a column header left or right to move the column (one undoable
step, saved with the document). Hold it near the table's left or right edge to scroll the
table: the closer to the edge, the faster it scrolls, and faster still if you drag past the
edge, so a column can go anywhere in a wide table. After the drop the table stays where you
scrolled to, the moved column stays in view and, if it was selected, stays selected.

Right-click a column header for:

- **Freeze Columns** keeps that column and every column left of it in place while
  the rest of the table scrolls sideways; **Unfreeze Columns** lets them scroll again. Drag the
  columns you want to keep in view (ID, Structure) to the front first. Frozen columns are saved
  with the document.
- **Align Drawings…** (Structure columns). Type, draw or take from the current row a core
  (SMILES or SMARTS); every structure that contains it is drawn with the core in the same place
  and orientation, so a series lines up down the column. Structures without the core keep
  their usual drawing. Only the drawings change, not the structures. **Stop Aligning
  Drawings** brings the usual drawings back. The core is saved with the document's layout, and
  the larger drawing you get by hovering over a cell is laid out on it too.
- **Highlight Atoms ▸** (Structure columns). Mark the atoms another column names on the
  drawings: the submenu lists the columns whose values name atoms of these structures, such as
  Top Site (3A4) or Sites (3A4) from [Sites of Metabolism](#sites-of-metabolism). The
  best-scoring atom is marked strongest and the others more faintly, in the order the column
  lists them, so a ranking reads at a glance; a column naming one atom marks that atom alone.
  **None** stops, and so does Data ▸ Search ▸ Clear Highlighting. A value naming an atom the
  molecule in that row does not have (the wrong element, or an index past the end) marks
  nothing, so a column written for other structures is ignored rather than drawn wrongly. Only
  the drawings change, and the choice is saved with the document's layout; write the column
  again and the marks move with it.
- **Select ▸** picks shown rows by this column's value, top to bottom: **Empty**, **Not Empty**,
  **First Occurrence** (the first row of each value), **Duplicates** (every row whose value
  appears in a row above it, so deleting them keeps one of each) and **All Duplicated Values**
  (every row whose value appears more than once). Structures count as the same when they are
  the same molecule, however they were written.
- **Sort ▸ Ascending / Descending.** **Then Ascending / Then Descending** add a second (third…)
  sort key; **Clear** returns rows to their original order.
- **Search.** Opens the search bar set to this column, ready to type.
- **Filter.** Adds a filter card for the column (see [Filters](#7-filters)).
- **Color…** Colors cell backgrounds: **Low → high** or **Low → middle → high**
  gradients for numbers, or **One color per value** for any column (the same value always gets
  the same color). Range ends left empty follow the data. Choose **None** to remove the
  colors.
- **Rename “column”…**, **Duplicate “column”** (a copy, with its formatting, to the right),
  **Delete “column”**, each naming the column you right-clicked. Delete acts on all selected
  columns when you right-click one of them (**Delete 3 Columns**).
- **Column Type ▸ Text / Number / Structure / Reaction.** Converting to Number warns you if
  some values are not numbers; those cells would be cleared. Converting to Structure or
  Reaction reads each cell's text; cells that do not read are kept and marked invalid.
- **Logarithmic (log₁₀)** (Number columns). Converts the values to log₁₀; choose it again to
  convert them back (10ˣ). The header shows "(log₁₀)". It needs every value above zero. One
  step you can undo.
- **Precision…** (Number columns). How every value is shown:

  | Format | 1234.5678 (or 0.4567 for Percent) with 2 digits |
  |---|---|
  | As typed | each cell as it was entered |
  | Fixed decimals | 1234.57 |
  | Thousands separator | 1,234.57 |
  | Scientific | 1.23e+03 |
  | Significant figures | 1200 |
  | Percent | 45.67% |

  The stored values are not rounded, so sorting, filters and tools still use the full numbers,
  and editing a cell starts from the number as typed. The quick search and Find match what the
  cells show.
- **Hover over a structure cell** for a larger drawing of the structure, drawn as the cell
  draws it: whatever the column marks (Highlight Atoms, or a structure search's match) marked
  on it, and laid out on the column's core if it is aligned.
- Precision, Logarithmic and colors are saved in the document and can be undone.
- **Hover over a column's header** for a summary of the column: how many cells are filled and
  blank, and for numbers the range, mean ± standard deviation and median, or for text the
  number of distinct values and the commonest one.

From Data ▸ Table ▸ Operations:

| Command | What it does |
|---|---|
| Add Compound… | Opens the sketcher; draw (or paste) a structure, and **File ▸ Add** appends it as a new row. |
| Add Rows… | Adds empty rows. |
| Add Column… | Adds a column. |
| Split Column… / Join Columns… | Each with a live preview of the first rows shown. Choose **All**, **Shown** or **Selected** rows; other rows stay blank. Split finds the delimiter itself (**Auto**) or uses one you choose, and keeps every field or only the first, last, largest or smallest; Join can skip blank cells so no delimiters are doubled. |

### Rows

- **Edit ▸ Delete Selection (Ctrl+Backspace)** deletes the selected rows. You can undo it.
- **Delete** deletes rows picked on the row numbers, in one step you can undo. On cells
  selected in the table (including Select All), Delete clears the cells instead.
- **The row numbers** on the left are kept narrow: past 999 they are shortened (1.2k, 12k,
  1.2M); point at one for the full number. Click a row number to select the row; drag down the
  numbers for several.
- **Right-click a row number** for:
  - **Delete Rows** and **Duplicate Rows** (for the selected rows; the copies are added at the
    end, or where the sort puts them, and are selected);
  - **Freeze Rows up to Here**: that row and the ones above it stay in place while the rest
    scroll up and down (with frozen columns as well, the corner stays put). It freezes the
    rows as they are shown, so re-sorting changes which rows those are. **Unfreeze Rows** lets
    them scroll again. Frozen rows are saved with the document.
- **Data ▸ Search ▸ Show All Rows** clears every filter, search and "show only" at once.
- **File ▸ New from Shown Rows** puts exactly the rows you can see, in their shown
  order, into a new untitled document. Column types are kept, and the proteins the rows name
  in their **Protein** column come with them.

---

## 4. Structures

- Structures are edited in the **sketcher** only. **Double-click** a structure, press **F2**
  or start typing in the cell, or right-click and choose **Edit in Sketcher…**. **File ▸ Done** (Ctrl/⌘+Enter) writes the drawing back as one step
  you can undo; **File ▸ Close** leaves the cell as it was. The
  drawing is stored as a molfile, which keeps your layout and wedges.
- To bring in SMILES, InChI or a molfile, paste it onto the drawing (Ctrl+V), or paste it
  straight into the cell from the table.
- A cell whose text is not a valid structure opens an empty sketcher. Its text stays in the
  cell until you draw a replacement.
- The right-click menu of a structure or reaction cell has **Copy As** (Molfile, SMILES,
  InChI, InChIKey, Image) rather than a plain Copy, so you choose what is copied.
  Copy As works on all the selected structures in the column. A molfile pastes into drawing
  programs.
- A structure that cannot be read is shown as **invalid**. Its original text is kept, so you
  can fix it.
- Rows are taller when the table has a structure column; see [Structure size](#structure-size).

### Structure size

A drawing fills its cell, so it grows with both the row height and the column width. Every
row has the same height.

- **Drag a row divider** on the row numbers (the pointer becomes a double arrow) to make every
  row taller or shorter; the divider stays under the pointer. Double-click a divider to go back
  to the Settings height. This works on frozen rows too.
- **Ctrl+scroll wheel** (⌘+scroll on macOS) over the table zooms the drawings: the rows and the
  structure and reaction columns grow or shrink together, and zooming back returns them to
  exactly where they were. Text and number columns keep their width.
- **Layout ▸ Larger Structures** (Ctrl+= or Ctrl++), **Smaller Structures** (Ctrl+-) and
  **Reset Structure Size** (Ctrl+0) do the same from the keyboard. Reset puts the rows back to
  the Settings height and leaves the column widths alone.
- Rows go from 40 to 800 pixels. The height you choose is saved with the document's layout.
  New documents start at the height under Settings ▸ Theme (110 pixels to begin with), and
  changing it there resizes every open table.
- Drag a column's edge to make that column wider or narrower.

### Reactions in the table

A reaction column holds one reaction per cell: reactants, agents (solvents, catalysts,
reagents) and products, with any atom maps. It is drawn in the cell with the arrow and plus
signs, agents above the arrow; map numbers are not shown in cells. Reaction columns come from
RD and RXN files, from reaction SMILES in CSV and Excel files, or from **Column Type ▸
Reaction** on a column of reaction SMILES.

- **Double-click** a reaction cell, or press **F2**, to draw it in the sketcher, with the map
  numbers. OK writes it back as one step you can undo, keeping your drawing for RD and RXN
  export. A reaction cell needs a reaction: OK asks for an arrow if the drawing has none.
- **Typing** in a reaction cell edits its reaction SMILES in place. A cell whose text is not a
  valid reaction opens that way too, so you can fix it.
- **Copy As ▸ RXN File** copies one reaction as an RXN file, which pastes into ChemDraw,
  Marvin and Ketcher as a drawing (several copy as an RD file). **Reaction SMILES** copies
  one per line; **Image** copies a picture of the first.
- Pasting an RXN file fills the selected reaction cells.
- The cell shows the reaction SMILES as typed or read. A reaction read from an RXN or RD file,
  or drawn in the sketcher, shows RDKit's reaction SMILES and keeps its drawing.
- Quick search and plots leave reaction columns out, as they do structure columns.

### Reaction Search

Data ▸ Search ▸ Reaction Search… finds reactions in a reaction column. Type a reaction SMILES
or SMARTS, or click **Draw…** and draw the query in the sketcher.

| Find | Matches a reaction when |
|---|---|
| Reaction contains | each reactant, agent and product you draw is a substructure of one of the reaction's own (draw agents above the arrow) |
| A reactant contains | the structure is in one of its reactants |
| A product contains | the structure is in one of its products |
| An agent contains | the structure is in one of its agents (solvents, catalysts, reagents) |
| Any component contains | the structure is anywhere in the reaction |

- **Output**: **Select them** selects the matching rows; **Add a column** marks every row Yes
  or No in a "Reaction Search" column.
- Search runs in the background and works on all, shown or selected rows.

### Split Reactions

Chemistry ▸ Reaction ▸ Split Reactions… writes the parts of each reaction as structure columns.

- **One per role**: a Reactants, an Agents and a Products column, each holding all of that
  role's components.
- **One per component**: Reactant 1, Reactant 2 … Product 1 …, in the order the reaction
  lists them.
- **Keep atom map numbers** is off by default, so the same molecule compares equal across
  reactions in structure search, duplicates and descriptors.

### Reaction Products

Chemistry ▸ Enumerate ▸ Reaction Products… applies a reaction template to building blocks.

- **Reaction**: a named reaction, or **Custom** to write your own. The named ones are grouped:
  - **Coupling**: amide (from an acid or an acyl chloride), sulfonamide, urea, carbamate,
    ester, Suzuki-Miyaura, Buchwald-Hartwig, Sonogashira, Chan-Lam, SNAr, reductive amination,
    N-alkylation, Williamson ether and the copper-catalysed azide-alkyne cycloaddition.
  - **Ring Formation**: Hantzsch thiazole, Paal-Knorr pyrrole, pyrazole, benzimidazole,
    1,2,4-oxadiazole, Pictet-Spengler and tetrazole from a nitrile.
  - **Multicomponent**: Ugi (four components), Passerini, Biginelli,
    Groebke-Blackburn-Bienaymé and Van Leusen imidazole.
  - **Deprotection**: Boc removal and ester hydrolysis.

  Each one names its components, so the building-block rows read "Aryl halide or triflate" and
  "Boronic acid or ester" rather than "Reactant 1" and "Reactant 2". Hover a reaction for what
  it assumes and where it comes from. The templates are chemoselective: an amide nitrogen is
  not treated as an amine, an SNAr needs a ring activated by a nitrogen or a nitro group, and
  a reagent that is never varied (TosMIC in a Van Leusen, azide in a tetrazole synthesis) is
  built into the product rather than asked for.
- **Template**: the reaction SMARTS itself. A named reaction fills it in; editing it makes the
  reaction **Custom**. For your own, write atom-mapped reaction SMARTS such as
  `[C:1](=[O:2])[OH].[N;H2,H1;!$(NC=O):3]>>[C:1](=[O:2])[N:3]` (amide coupling), or **Draw…**
  and map the atoms with the sketcher's Atom Map tool. Up to four reactants.
- **Each component**: a structure column of this table, **From a File…** (SDF, CSV or any
  table chemsuite opens), or **A Single Structure…** typed or drawn, to run one building block
  against every other.
- **More Options ▸ Swap Components** exchanges the first two components when the building
  blocks are the other way round.
- **Combine: Row by row** reacts reactant 1 and reactant 2 from the same row and writes the
  product to a **Product** column, with **Products Found** counting how many products the row
  gave (a building block with two reactive sites gives two).
- **Combine: Every combination** reacts every distinct building block with every other and
  opens the virtual library in a new window, one product per row with its building blocks.
  The building-block columns take the names of a named reaction's components. Up to 2,000,000
  combinations.
- Cancelling keeps the products already made.
- A component whose building blocks cannot fill its role is reported before the run, for
  example "Nothing in Amines is a carboxylic acid".
- Every named reaction is checked against hand-worked products, including that stereocentres
  away from the reaction are kept.

### Prepare Reagents

Chemistry ▸ Enumerate ▸ Prepare Reagents… turns reagents into substituents. A reagent as a
catalogue lists it is not a substituent: an acyl chloride carries a chloride that leaves when
it acylates, so what joins a core is the acyl fragment. This writes that fragment with an
attachment point, ready for [R-Groups on a Core](#r-groups-on-a-core): acetyl chloride
`CC(=O)Cl` becomes `*C(C)=O`, which on the core `C1CCN([*:1])CC1` gives the amide
`CC(=O)N1CCCCC1`.

- **Reagents** is a structure column of the open table, so open a supplier's file with File ▸
  Open first.
- **Reagent** is what kind they are, grouped by the bond each makes:
  - *Acylation*: acyl chloride, carboxylic acid, sulfonyl chloride, chloroformate, isocyanate.
  - *Arylation and Alkylation*: aryl halide or triflate, boronic acid or ester, alkyl halide or
    sulfonate, aldehyde, terminal alkyne.
  - *Nucleophile*: amine (N-H), alcohol or phenol, thiol.

  Each kind's tooltip says what leaves and what attaches. An aryl halide loses iodide before
  bromide before triflate before chloride, so a bromo-chloroarene loses only its bromide.
- **Output**:
  - **New Window** (the default) holds one row per fragment — `Fragment`, `Reagent` and
    `Sites Found` — to save as a reagent file and hand back to R-Groups on a Core.
  - **Columns** writes the first fragment of each row into the table, beside the reagent it
    came from.
- **More Options ▸ Clean up** (on) strips salts and neutralises first, so a hydrochloride is
  prepared like its free base.
- **Every matching site is prepared**, and identical fragments collapse: a symmetric diacid
  gives one fragment, an unsymmetric diol gives two. `Sites Found` counts them, which matters
  because **a fragment from a bifunctional reagent still carries the other group, live** —
  preparing terephthaloyl dichloride gives `*C(=O)c1ccc(C(=O)Cl)cc1`, with the second acyl
  chloride intact. Sort or filter on `Sites Found` to check those before building a library.
- Reagents of another kind are counted in the message rather than quietly skipped, and
  **cancelling** keeps the rows already prepared.
- When no reaction template fits the bond you want, this is the way to build the library;
  when one does, [Reaction Products](#reaction-products) needs no preparation step at all.

### R-Groups on a Core

Chemistry ▸ Enumerate ▸ R-Groups on a Core… makes every combination of the substituents you
choose on one numbered core — the other direction from
[R-Group Decomposition](#r-group-decomposition).

- **Core**: SMILES with numbered attachment points, such as `c1cc([*:2])ccc1[*:1]`. This is
  exactly what [R-Group Decomposition](#r-group-decomposition) and
  [Central Ring](#central-ring) write, so a decomposition feeds straight back in. **Draw…**
  opens the sketcher, where the **R** button places the attachment points: an R1 on a ring
  nitrogen comes back as `C1CCN([*:1])CC1`, so the substituent replaces that hydrogen. A row
  appears for each position, up to six.
- **Substituents** for each position come from a column, a typed list or a file.
  - A **column** uses its distinct values over the chosen rows: point R1 at the `R1` column a
    decomposition wrote and you get the substituents the series already uses.
  - **Typed…** takes one substituent per line (or comma separated), written as `*Cl`, `*OC`,
    with `H` for an unsubstituted position.
  - **From a File…** takes them from a file of prepared reagents — what
    [Prepare Reagents](#prepare-reagents) saves — so a curated set of acyl chlorides or amines
    can be pointed at one position and reused across cores. Any file chemsuite opens will do
    (SDF, SMILES, CSV, TSV, Excel). The column holding the attachment points is the one used,
    so a saved set's `Reagent` column is passed over in favour of its `Fragment` column,
    whichever order they are in. The file is read when you press Enumerate, so until then the
    button says just **Enumerate** and the number of combinations is unknown.
  - **and H** also leaves the position unsubstituted, without adding `H` to the list by hand.
    It applies to a file as well.
- **Output** is a virtual library in a new window: one product per row, with an `R1`, `R2`…
  column holding the substituent it came from. Up to 2,000,000 combinations.
- A substituent that cannot be read, or cannot be attached, is counted in the message rather
  than quietly skipped.
- **Cancelling** keeps the products already made.

### Recombine Fragments

Chemistry ▸ Enumerate ▸ Recombine Fragments… builds new molecules from the BRICS fragments
in the table. Two pieces are joined only where the BRICS rules allow (RDKit's `BRICSBuild`),
so an acyl goes onto an amine or linker, never straight onto an aromatic carbon.

- **Rows**: the fragment pool is every distinct fragment in the `BRICS 1`, `BRICS 2`…
  columns of these rows. Run [Fragments](#fragments) with BRICS first, with **Plain
  attachment points** off and, to keep amide and ether linkers, **Smallest piece** at 1.
  RECAP fragments carry no rule labels and cannot be recombined.
- **Most products** stops the build once that many molecules have been made (2,000 by
  default).
- **Most fragments joined** is how many pieces one product may be built from (3 by default).
- **More Options**:
  - **Drug-likeness** fills the property ranges: **Rule of 5** (Lipinski 1997: MW ≤ 500,
    cLogP ≤ 5, HBD ≤ 5, HBA ≤ 10) or **Lead-like** (Teague 1999: MW 200–350, cLogP −1 to 3,
    HBD ≤ 3, HBA ≤ 6, rotatable bonds ≤ 7).
  - **MW, cLogP, TPSA, HBD, HBA, RotB**: the least and most of each; `—` leaves that end
    unchecked.
  - **Skip compounds already in the table** (on by default) drops products that match a
    structure in the table.
- **Output** is a new window, one product per row. The message says how many candidates the
  limits and the novelty check removed. A run with tight limits gives up after examining
  50 candidates per product asked for (at least 5,000), and says so.
- **Cancelling** keeps the products already made.

**Add Filter** on a reaction column gives a reaction card: the same query and **Find** choices
as Reaction Search, applied as a filter together with the other cards.

### The sketcher

Opened by double-clicking a structure cell, F2, **Add Compound…**,
**Draw…** in the structure fields of Structure Search and Similarity Search, or
blank with **Tools ▸ Sketcher…**.

Every sketcher window has two menus instead of buttons:

- **File:** **Save…** (Ctrl/⌘+S) writes the drawing to a file: an MDL Molfile (.mol), SD
  file (.sdf) or SMILES (.smi), or for a reaction an RXN file (.rxn) or reaction SMILES.
  Then the window's own action (Ctrl/⌘+Enter): **Done** gives the drawing back to the cell or
  field it came from and closes the window; **Add** (Add Compound) adds it as a new row; and
  in the Sketcher opened from Tools, **Add to Table** adds it as a new, selected row and
  leaves the window open for the next drawing. **Close** (Ctrl+W, ⌘W) closes the window
  without giving the drawing back, as closing the window does.
- **Edit:** the main window's Edit items for the drawing: Undo, Redo, Cut, Copy, Paste,
  Select All, Invert Selection (Ctrl/⌘+I), Clear Selection and Delete Selection.

Drawing tools are along the top, elements down the left and rings along the bottom. Done,
Add and Save are available only when every atom has a valid valence; a line under the drawing
says what is wrong. Undo, Redo, Cut, Copy, Paste, Select All and Clear are on the keyboard (Ctrl/⌘+Z, Y, X,
C, V, A) and on the drawing's right-click menu.

**Tools**

| Tool | Click | Drag |
|---|---|---|
| Select | Select an atom or bond (Shift adds; an abbreviation is selected whole). Double-click an atom to type its label; double-click a bond to select its whole fragment. | From empty canvas: box-select. From a selection: move it. Atoms dropped on other atoms merge with them, so dragging a chain end onto its start closes a ring. With two or more atoms selected, a round handle appears above them: drag it to rotate the selection about its centre (15° steps; hold Ctrl (⌘) for any angle). |
| Lasso | As Select | From empty canvas: draw a loop round the atoms to select (Shift adds). Moving and rotating work as in Select. |
| Erase | Delete the atom or bond | Delete everything you pass over (one undo step) |
| Single, double, triple bond | On an atom: add a bond in the free direction (chains zig-zag; alkynes and allenes stay straight). On a bond: set that order; the single-bond tool steps single → double → triple. On empty canvas: a new bond. | From an atom or empty canvas: draw a bond. It snaps to 15° steps; hold Ctrl (⌘) for any angle. Release on an atom to bond to it. |
| Wedge, hash, wavy | As for bonds; clicking a bond that already has that wedge reverses it | As for bonds |
| Chain | Add a bond | Draw a zig-zag chain, one bond per bond length; the count shows by the pointer |
| Rings (benzene, 3- to 8-membered) | On empty canvas: place the ring. On an atom: the ring shares that atom, facing the free direction. On a bond: fuse the ring onto it, away from the rest of the molecule. Fusing benzene re-draws the alternating bonds of the whole aromatic system, so naphthalene or phenanthrene come out valid whichever bond you click. | On empty canvas: turn the ring (15° steps). From an atom: aim it; drag more than a bond length to attach it through a new bond. On a bond: choose the side. |
| Charge + / − | Change the atom's charge by one | — |
| Element (C, N, O …), R, A | On an atom: change its element (or make it an R-group or "any atom"). On empty canvas: a lone atom. Double-click an atom to type its label. | From an atom: a bond to a new atom of that element |
| Reaction arrow | A plain arrow, 2.5 bond lengths long, starting where you click | Draw the arrow from tail to head (15° steps; hold Ctrl (⌘) for any angle) |
| Reaction plus | Place a plus sign between two components | — |
| Atom map | On a mapped atom: clear its map number | From an atom to the same atom on the other side of the arrow: both get the next free map number |

**Keys over an atom or bond.** Point at an atom or bond and press a key:

| Key | Over an atom | Over a bond |
|---|---|---|
| `c` `n` `o` `s` `p` `f` `h` `i` | That element | — |
| `l` / `b` / `d` | Cl / Br / deuterium | — |
| Shift+`s` / Shift+`n` / Shift+`b` | Si / Na / B | — |
| `1` `2` `3` | Add a single, double or triple bond | Set the order |
| `4` | — | Aromatic bond |
| `0` | — | Any-bond query |
| `=` / `-` | Charge +1 / −1 | — |
| Delete, Backspace | Erase | Erase |
| Enter | Type a label (see below) | — |

With nothing under the pointer, the same keys choose the matching tool, and `t` steps through
the rings. Esc returns to Select.

**Template library.** Shift+T, or the last button of the ring bar, opens 64 ring systems a
medicinal chemist draws daily: aromatic and saturated heterocycles (pyridine to tetrazole,
piperazine, morpholine, oxetane…), fused bicycles (indole, 7-azaindole, quinazoline, purine,
tetrahydroisoquinoline…) and bridged and spiro systems (bicyclo[1.1.1]pentane, quinuclidine,
2-oxa-6-azaspiro[3.3]heptane…). Pick a group on the left or type part of a name; double-click
a drawing (or OK) to draw with it. It then places like a ring: on empty canvas, on an atom, or
fused onto a bond, with aromatic bonds re-drawn across the fused system.

**Other keys:** Ctrl+Z / Ctrl+Shift+Z or Ctrl+Y undo and redo; Ctrl+C, Ctrl+X, Ctrl+V copy,
cut and paste (molfile and SMILES, so ChemDraw and others can read it); Ctrl+Shift+M copies a
molfile as text; Ctrl+A selects everything; Alt+H and Alt+V flip the selection (or the whole
drawing) horizontally or vertically; Ctrl+Shift+K cleans the structure (below). Ctrl+scroll
zooms, scrolling or a middle-button drag pans.

**Clean Structure (Ctrl+Shift+K).** Redraws the molecule with even bond lengths and angles,
also on the background right-click menu. With atoms selected it cleans only the molecules they
belong to (**Clean Selection**). The cleaned drawing stays where it was and keeps its turn on
the page; it is never mirrored. Stereocentres and double bonds keep their configuration: your
wedges stay where they still mean the same thing, and are re-placed where they would not. If
a clean would change the stereochemistry, the drawing is left alone and a line says why.

**Salts and solvents.** Right-click the background for **Add Salt or Solvent**: 135 common
solvents, acids, bases, salts and buffers in groups (alcohols, ethers, polar aprotic, salts and
buffers …), named in full with their usual short name (DMF, THF, TFA). The one you pick is
drawn to the right of the drawing and selected, so you can move it straight away.

**Labels and abbreviations.** Double-click an atom (or point at it and press Enter) and type:

- an element, with an isotope, hydrogens or a charge if you like: `N`, `Cl`, `13C`, `NH2`,
  `O-`, `Fe2+`, `D`;
- or an abbreviation, which replaces the atom with that group and shows it as a label:
  Me, Et, Pr, iPr, cPr, Bu, iBu, sBu, tBu, Cy, Ph, Bn, Tol, Mes, Tr, PMB, Ac, Bz, Piv, Boc, Cbz, Fmoc, CHO, CO2H, CO2Me, CO2Et, CO2tBu, CONH2, CN, CF3, CCl3, OMe, OEt, OAc, OBn, OCF3, OTf, OTs, OMs, MOM, THP, SEM, NO2, NHAc, NHBoc, NMe2, NCO, NCS, SMe, SCN, Ms, Ts, Tf, SO2NH2, SO2Cl, SO3H, TMS, TES, TBS, TIPS, TBDPS, Bpin, PO3H2. Common spellings work too (COOH, TBDMS, Z, nBu …), and `Ac`, `Pr` and `Ts` mean
  acetyl, propyl and tosyl, not the elements.

An abbreviation goes on an atom with at most one bond. Its label is written outwards from the
bond, so CO₂Me on a bond from the right reads MeO₂C. The label moves, selects and erases as
one piece. Pressing an element key over it replaces the group with that atom. Right-click it
to **Expand** (show its atoms), **Contract** again, or **remove the label** and keep the atoms.
The molecule itself never depends on whether a label is contracted. Molfiles keep
abbreviations as superatom groups, which ChemDraw, Marvin and Ketcher show as labels.

**R-groups and query atoms.** For R-group cores, reaction templates and substructure
queries, an atom can stand for more than one element:

| Label | Means |
|---|---|
| `R1` … `R99` | R-group *n* (plain `R` is R1) |
| `A` / `AH` | any atom except hydrogen / any atom |
| `Q` / `QH` | any atom except carbon and hydrogen / any atom except carbon |
| `X` / `XH` | any halogen / any halogen or hydrogen |
| `M` / `MH` | any metal / any metal or hydrogen |
| `*` | any atom |
| `[C,N,O]` | one of the listed elements (up to 16) |
| `![F,Cl]` | none of the listed elements |

The **R** button under the elements places R-groups, numbering each new one with the next free
number (R1, R2, …); clicking an existing R-group keeps its number. The **A** button places
"any atom". Right-click an atom for the **R-Group** and **Query Atom** menus (with **Atom
List…** and **Back to Carbon**), or type any label above in the label editor. Query atoms are
drawn in the accent colour and carry no hydrogens of their own.

Molfiles write them the standard CTfile way (`R#` with `M  RGP`, `A`/`Q`/`X`/`M`, `L` with
`M  ALS`), so ChemDraw, Marvin and Ketcher read them; SMILES writes an R-group as `[1*]` and
`[1*]` or `[*:1]` read back as R1. In a table cell, a structure with query atoms shows as
SMARTS (SMILES has no wildcards) and keeps its queries, so it can serve as a substructure
query or R-group core.

**Bond queries.** For substructure searches a bond can match more than one kind of bond.
Right-click it for **Bond Query**, or point at it and press `0` for **Any bond**:

| Choice | Matches | Drawn |
|---|---|---|
| Any bond | any bond | dashed |
| Single or double / Single or aromatic / Double or aromatic | either kind | dashed, tagged S/D, S/A or D/A |
| Aromatic | aromatic bonds | dashed second line |
| Ring only / Chain only | a bond in a ring / not in a ring | tagged Rng or Chn |

**Back to Single** makes it a plain bond again. Atoms on a query bond show no hydrogen count,
since their bond order is open. Molfiles write these as standard CTfile query bonds.

**Attachment points.** For R-group members and substituent fragments, right-click an atom for
**Attachment Point**: **First**, **Second** or **First and second**. Each is drawn as a short
wavy stub with its number. Molfiles write them as `M  APO` (V3000: `ATTCHPT=`).

**Stereo labels.** Stereocentres show their CIP label, *(R)* or *(S)* (*(r)* / *(s)* for
pseudo-asymmetric centres), and double bonds show *(E)* or *(Z)*. The labels come from
RDKit's CIP labeller reading the drawing as any molfile reader would: wedges set
stereocentres and the drawn geometry sets double bonds. Right-click the background and untick
**Show Stereo Labels** to hide them all.

- Right-click a stereocentre for **Stereocentre (R/S)**: choose **R** or **S** to set it, or
  untick **Show Label** to hide just that label. Switching swaps the wedges and hashes that
  start at the atom; an undefined centre gets a new wedge, placed on a hydrogen or chain end
  where possible.
- Right-click a double bond for **Double Bond (E/Z)**: **E** or **Z** mirrors the smaller side
  of the bond, and wedges on that side swap so its stereocentres keep their labels. Double
  bonds in rings cannot be switched this way. **Unknown (crossed bond)** marks the geometry as
  not known: the bond is drawn crossed, has no E/Z label, and molfiles mark it "either".

**Wedge warnings.** A wedge or hash that does not say what it seems to is ringed in amber,
with the reason when you point at it and a note under the drawing. It is a warning only; OK
still works.

- **Defines nothing**: the atom at its narrow end is not a stereocentre.
- **Points the wrong way**: the stereocentre is at its wide end; the narrow end belongs there.
- **Contradict each other**: the wedges at an atom (for example a wedge and a hash on opposite
  bonds) give no configuration.
- **Ambiguous**: turning one bond slightly would change the configuration, as when the wedge
  sits between two bonds drawn in a straight line. Spread the bonds out.

**Enhanced stereo.** For mixtures and unknown absolute configuration, right-click a
stereocentre (or several selected ones) for **Enhanced Stereo**, or press Alt+E for the
selected centres:

- **Absolute (as drawn)**, the default, or **Absolute, labelled abs**;
- **AND Group**: a racemic mixture, this drawing and its mirror image, shown `&1`, `&2` …;
- **OR Group**: one of this drawing or its mirror image, not known which, shown `or1` …

Centres in the same group change together. The group shows after the R/S label, e.g.
*(R) &1*. A drawing with groups is saved as a V3000 molfile, which ChemDraw, Marvin and
Ketcher read; CXSMILES such as `C[C@H](N)O |&1:1|` pastes in with its groups.

**Explicit hydrogens.** Right-click an atom (or one of several selected atoms) for **Show
Explicit Hydrogens**, which draws its hydrogens as H atoms, or **Hide Explicit Hydrogens**,
which makes plain drawn hydrogens implicit again. The background menu has **Show All
Hydrogens** and **Hide All Hydrogens**. Stereocentres keep their R/S either way, and a
hydrogen on a wedge stays, since it defines its stereocentre. The molecule itself does not
change.

**Reactions.** Draw the reactants, an arrow and the products, with plus signs between the
molecules on each side. The sketcher sorts the drawing as Ketcher does:

- molecules between two plus signs (or a plus sign and the arrow) are one component, so a salt
  drawn without a plus stays together;
- a component left of the arrow's centre is a reactant, one to the right a product;
- a component above or below the arrow, within its length, is an agent (solvent, catalyst or
  reagent).

Right-click an arrow to change its kind (**Plain**, **Equilibrium**, **Retrosynthetic**,
**Failed**), **Reverse Arrow** or **Delete Arrow**. A retrosynthetic arrow points from the
product to the reactants, and is read that way round. Select, move, rotate and erase work on
arrows and plus signs as on atoms. A drawing has at most one arrow.

**Atom maps** pair each atom with the atom it becomes. Use the Atom map tool, or the atom
right-click menu: **Clear Map Number**, and **Clear All Mapping** on the background menu. Map
numbers are drawn in small blue type beside the atom.

Copy puts a reaction on the clipboard as an RXN file and as reaction SMILES, so ChemDraw,
Marvin and Ketcher can paste it. Paste takes either back, and reaction SMILES typed or pasted
in (`CC(=O)O.OCC>[H+]>CC(=O)OCC`) is laid out left to right, agents above the arrow. A
reaction cannot be put into a structure cell; the sketcher asks you to delete the arrow and
plus signs. Reactions belong in a [reaction column](#reactions-in-the-table).

**Valence.** An atom with an impossible valence is shown in red inside a dotted ring; point
at it to see why. Hydrogens are worked out by the same rules as the table, so a drawing the
sketcher accepts is always a valid structure there. A copied or saved molfile marks any atom
whose hydrogen count other programs would work out differently, so they read the same
structure.

---

## 5. Selecting rows and lists

### Selection

- Click, Shift-click and Ctrl-click cells, or click the row numbers for whole rows. Every row
  with a selected cell counts as selected; a selection of cells stays a selection of cells
  when the table is sorted or filtered.
  - **Select All (Ctrl+A)**, **Invert Selection (Ctrl+I)**, **Clear Selection
    (Ctrl+Shift+A)**. Clear Selection clears every kind of selection: cells, whole columns
    and rows, including selected rows a filter hides.
- The selection is shared by the table, the plots and the tools:
  - selecting points in a plot selects the rows, and the table scrolls to the first new one
    unless one of them is already in view;
  - selecting rows highlights their points.
- Many tools can run on all rows, the shown rows or just the selected rows.

### Named lists

Lists are saved sets of rows, such as "actives", "to make" or "sent for assay". They are saved
with the document.

- **Data ▸ Lists ▸ Add Selected to List… (Ctrl+L)** adds the selected rows to a new or existing
  list.
- **Data ▸ Lists ▸ Add Visible to List…** adds every row the table shows (after filters,
  searches and Show Only These Rows) to a new or existing list.
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

Data ▸ Search ▸ Search Rows opens the search bar above the table. As you type, the table shows only the
rows that match, and the count of matching rows appears beside the box.

- The **column menu** at the left of the bar chooses where to look: **All Columns** (every text
  and number column) or one column.
- **Field syntax.** Several words must all be found, anywhere in the row. Beyond that:

  | Type | Finds rows where |
  |---|---|
  | `kinase "cell line"` | some cell contains *kinase* and some cell contains the phrase *cell line* |
  | `ID:CPD-12` | the ID column contains CPD-12 (`ID:"two words"` for spaces) |
  | `pIC50>7` | pIC50 is above 7; also `>=`, `<`, `<=`, `=` and `!=` |
  | `Series=A` | Series is exactly A (any case); `Series!=A` for anything else |
  | `-Series:B` | Series does not contain B; a `-` or `!` excludes any term (`!kinase`) |
  | `>7` | with no column named, the column chosen in the column menu is above 7 (also `<=5`, `=A`, `!=A`) |
  | `ID:CPD*` | the whole ID starts with CPD; `*` stands for any text (`*ase`, `Series=A*`) |
  | `ID:*` | ID has a value; `-ID:*` finds blank IDs |
  | `kinase \| gpcr` | either term is met |
  | `Series=A\|B` | Series is A or B; a value after `\|` keeps the column and operator before it (`pIC50<5\|>8`, `Series!=A\|B` for neither) |
  | `"Mol Wt"<400` | quote a column name that has spaces |

  Column names ignore case. A `*` inside quotes is plain text (`"a*b"`), and so is a `|`. On
  All Columns, `>7` matches a number in any column, so a column such as molecular weight meets
  it in nearly every row: choose a column first. A word before a colon that is not a column name is searched for
  as text, so `10:30` or `http://` still work. Whole Cell and Regular Expression take the text
  exactly as typed. Hover over the box for a reminder.
- **Recent searches.** The clock button left of the box lists the last ten searches; pick one to run it again. A
  search goes on the list when you press Return or leave it in place for a moment. The list is
  saved with the document.
- **Substructure search.** Choose a structure column and type a substructure as SMILES or
  SMARTS (for example `c1ccccc1` or `[OH]a`). The table shows the rows that contain it, and
  the matching atoms are highlighted in the drawings. The search runs in the background
  shortly after you stop typing. A query that cannot be read turns the box red; hover over it
  to see why.
- The three toggles right of the box are **Match Case** (Aa), **Whole Cell** (ab) and
  **Regular Expression** (.*), for text searches only. Whole Cell and Regular Expression turn
  each other off.
- If the chosen column is removed, the bar goes back to All Columns. A renamed column stays
  chosen. The column is saved with the document, with the search text and options.
- Escape clears the search. A second Escape, or the ✕ button, hides the bar.
- Hiding the bar always clears its search, so you never have rows hidden by a search you cannot
  see.

### Find and Replace (Ctrl+H)

A window that steps through matching cells, with **Previous**, **Next**, **Replace** and
**Replace All**; the match count ("3 of 12") shows beside the buttons. A replacement is one
step you can undo. **Hidden rows too** makes Replace All also change rows hidden by filters or
searches.

- **In** chooses where to look: **All columns**, **Selected columns** (the columns of the
  selected cells), **Selected cells**, or **Chosen columns** (tick them under **Columns…**).
  The selection is read when you choose the scope or open the window, so stepping from match
  to match, which moves the cursor, does not shrink it.
- **Groups.** With **Regular expression**, the replacement can reuse what the groups matched:
  `$1` or `\1` for the first group, `${name}` or `\g<name>` for a named one, and `$$` for a
  dollar sign. Finding `(\w+)-(\d+)` and replacing with `$2-$1` turns `CPD-12` into
  `12-CPD`. A reference to a group the expression does not have is reported, and nothing is
  replaced.

### Advanced Search

Data ▸ Search ▸ Advanced Search… opens a window of conditions. It is always its own floating window, so
it never takes room from the table or the Filters panel. Move and resize it as you like; it
reopens where you left it.

- Each condition picks a column (or **any column**), an operator and a value.
- Click **Not** to negate a condition.
- Operators depend on the column type:
  - **Text:** contains, equals, starts with, ends with, matches a regular expression, is
    blank. Each has a match-case option.
  - **Number:** =, ≠, <, ≤, >, ≥, between (either order), is blank. A blank cell never
    counts as "≠ 6".
  - **Structure:** contains a substructure, is exactly, is blank.
- **Match all / any** combines the conditions: all means every condition must match (and),
  any means one is enough (or).
- **Groups.** **Add Group** adds a group of conditions with its own **Match all / any** and
  **Not**, and groups can hold groups. For example, *Series equals A* and a group set to
  **Not**, **any** of *pIC50 ≥ 8* and *ID equals C4* finds Series A compounds that are neither
  potent nor C4: `Series equals A and not (pIC50 ≥ 8 or ID equals C4)`. The ✕ at a group's
  right removes it with its conditions. A group with no conditions keeps the search buttons
  disabled.
- To exclude values from several columns, tick **Not** on each condition: *not (Series
  equals A)* and *not (ID contains 5)*.
- A condition with a problem (for example, "between 6 and lots") is marked, and the search
  buttons stay disabled until it is fixed.
- Use the result in one of three ways:
  - **Show** shows only the matching rows. This works together with your filters.
  - **Select** selects the matching rows.
  - **Add to List** adds them to a named list.
- While a search is showing, text and number conditions follow your edits immediately, and
  structure conditions search again when a structure changes.
- The search (with its groups), whether it is showing, and your recent searches are saved
  with the document.

### Structure Search (Ctrl+Shift+M)

Data ▸ Search ▸ Structure Search… searches a structure column.

- **Query.** Draw (paste) or type the query as SMILES, a molfile or InChI. For substructure
  searches you can also type SMARTS, such as `[C,N]`.
  - The pencil opens the sketcher. A plain drawing comes back as SMILES. **Current Row** takes
    the structure in the table's current row.
  - A drawing with query features comes back as SMARTS that keeps them: R-groups, `A`/`Q`/`X`/`M`,
    atom lists, and bond queries. An R-group or `*` matches any atom there.
  - Drawing again reopens your drawing, not the SMARTS.
  - The same applies in the Structure filter card and the structure conditions of Advanced
    Search.
- **Modes:**

  | Mode | Finds |
  |---|---|
  | Substructure | Rows that contain the query. |
  | Exact | Rows that are the query, with or without stereochemistry. |
  | Similarity | Rows at or above a similarity threshold, using a fingerprint and metric you pick. |

- **Matches.** **Filter** shows them alone (as a structure card in the Filters panel, so it
  applies to every row), **Select** selects them, and **Add to List** adds them to the list you
  name. Select and Add to List search **All**, **Shown** or **Selected** rows.
- **Highlighting.** Substructure matches can be highlighted in the table's drawings. Data ▸ Search ▸ Clear
  Highlighting turns that off.
- **Speed.** Searches are fast even on 200k rows: a fingerprint screen rules out most rows
  before any atom-by-atom matching.

---

## 7. Filters

Data ▸ Filters (Ctrl+Shift+F) opens the Filters panel, pinned at the right edge of the window
at a fixed width. Each **card** filters on one column, and
a row is shown only if it passes every enabled card.

| Card | For | Controls |
|---|---|---|
| Number Range | Number columns | A two-handle slider over a histogram of the values, plus min/max boxes. Choose whether blanks are kept. |
| Category | Text or number columns (best with few distinct values) | Tick boxes, one per value, with counts. |
| Text Match | Any column | Contains / equals / starts / ends / regular expression, with **Not**. |
| Structure | Structure columns | Substructure, exact or similarity query; searches in the background after you stop typing. |
| Reaction | Reaction columns | A reaction query (see [Filtering by reaction](#filtering-by-reaction)). |

- **Adding a card.** The **+** button in the panel's header lists the kinds of card. Pick one
  and it starts on a suitable column (for Category, the first column with a short list of
  values); choose its column in the box at the top of the card. Switching the column keeps the
  card's place, on/off and Not. A kind with no suitable column is greyed out. A column
  header's **Filter** adds the best kind of card for that column directly. The
  panel header's other buttons are **Filter Sets** and **Remove All Filters**; point at Remove
  All for how many rows are showing.
- Each card has an on/off tick box, a **¬** (Not) toggle that keeps the rows the card would
  hide, and **×** to remove it. Cards alternate in shade, like table rows. Point at a card's
  column box for its kind and its values (range, or how many categories are listed).
- A category card's list is as tall as its values (up to eight, then it scrolls); the two tick
  boxes beside **Find a value** tick or untick every value. The Structure card's pencil opens
  the sketcher.
- **Filter sets.** The **Filter Sets** menu saves the current cards under a name (such as
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
- **Status bar.** The running job shows at the far left of the status bar. Tools that work
  through molecules or rows one at a time count them as they finish: "Predicting pKa…
  (120/500)", "Generating conformers… (40/500)". Work on the whole set at once, such as a
  plot, a model fit, an embedding, clustering or Matched Molecular Pairs, shows only what is
  being calculated: "Calculating TMAP…", "Calculating UMAP…". Click it to open the Log window
  at its Jobs tab.
- **Log**, at the right end of the menu bar (**Window ▸ Log** on macOS), opens the Log window.
  A red dot beside it shows that a job failed and you have not looked yet. The window opens at **Jobs** while jobs run or a failure is unseen,
  otherwise at **Messages** (see [The log](#16-the-log-and-reporting-problems)).
- The **Jobs** tab lists running, queued and finished jobs with their elapsed time and result.
  The tab's title counts the running and queued jobs, and turns red when one has failed.
  From there you can:
  - cancel a job;
  - remove a queued job, or move it up or down;
  - **Remove** finished, failed or cancelled jobs from the list (select one or several, or
    press Delete); **Clear Finished** removes them all;
  - copy an error;
  - **Show Messages**: the job's log messages, on the Messages tab (or double-click the job).

  Right-click a job for the same commands.
- **Cancelling.** A job stops within about a tenth of a second. While it stops, the status bar
  says "Cancelling …".
- **Partial results.** Tools that work row by row keep what they finished when you cancel:
  - the finished rows are written to the table;
  - the message says "partial results kept";
  - rows not reached are left as they were.

  | Partial results kept on cancel | Discarded on cancel |
  |---|---|
  | Calculate Descriptors | Clustering |
  | Predict pKa | |
  | MedChem plots (descriptor step) | Dimensionality Reduction |
  | Generate Conformers | Self-Organizing Map |
  | Superpose | |
  | Standardize Values | |
  | Similarity Search | Free-Wilson |
  | R-Groups (numbered core) | Hansch |
  | Diverse Subset (the picks made so far) | Matched Pairs |
  | | Activity Cliffs |

  The tools in the right-hand column have no meaningful result until the whole run is complete.
- **Safety.**
  - Closing a document cancels its queued jobs and asks about running ones.
  - Closing the last document window (its ✕) asks first when jobs you started are
    running or queued, in any window, and says how many background tasks (pictures, views
    made ahead, structures being read) stop with them. Background work alone stops without
    asking.
  - When the last document window closes, everything closes with it: the Log, Sketcher and
    Macromolecule Viewer windows, the worker processes and external programs such as AutoDock
    Vina, together with anything they started. Nothing keeps running after chemsuite is
    gone, and a save that was cut short leaves the old file as it was.
  - Results from a job are never written into a different document.
- **Long runs.**
  - The computer is kept awake while jobs run.
  - When a job that took more than 30 seconds finishes or fails while chemsuite is not in
    front, the dock icon (taskbar on Windows) asks for attention and a desktop notification
    names the job and its result. Settings ▸ General ▸ **Notify when a long job finishes** turns
    the notification off.

---

## 9. Plots

Open a plot from the **Plot** menu. It opens to the right of the table. The plot area is
white with dark axes in every theme, so points, lines and labels stay legible in dark mode;
the toolbar and hover cards follow the theme.

| Plot | Shows |
|---|---|
| **Scatter** (Ctrl+Shift+P) | One number column against another, coloured and sized by others. |
| **Line** | Points joined in X order, one line per category. Rows sharing an X value (replicates, a time point) are joined through their mean. |
| **Histogram** | Rows per value range (numbers) or per value (text). |
| **Bar Chart** | The mean, median, sum, minimum or maximum of a number column per group, or the row count, with optional ± SD error bars. |
| **Box Plot** | The spread of a number column per group. |
| **Density** | Where crowded points sit: a 2D histogram of X and Y, shaded by count. Box, lasso and click still select the rows beneath. |
| **Radar** | Three or more number columns as spokes, each scaled to its own range (centre: lowest, rim: highest). Shows the mean of the plotted rows in grey, and each selected row (up to 12) in its own colour. |
| **3D Scatter** | Three number columns in a cube you can turn, coloured and sized by others (see below). |
| **Constellation Plot** | The cores that analog series share, placed in chemical space: each circle is a core, sized by its compounds and coloured by their activity (see below). |
| **Chemical Space Network** | Compounds joined to those they resemble, or form a matched pair with; each cluster laid out by force (see below). For hundreds to a few thousand compounds. |
| **TMAP** | A tree of nearest neighbours: each compound sits on a branch beside its closest relatives (see below). For sets of up to hundreds of thousands. |
| **GTM Landscape** | A trained Generative Topographic Map: a column's landscape over the map, with each compound at its place (see [Generative Topographic Map](#generative-topographic-map)). |
| **3D Landscape** | An activity landscape: a smooth surface of one column (the height, e.g. pIC50) over two others, usually a chemical-space map, with the compounds standing on it (see below). |

**Switching plots.** The cog's **Plot** list holds every kind, 2D and 3D. Changing it redraws
in the same pane: the columns you chose, the title, the point cards and the selection carry
over, and the toolbar keeps the same buttons. Its first button is **Pan and Zoom** on a 2D
plot and **Turn** on a 3D one. A plot fits its pane until you pan or zoom, so a plot drawn
while its pane was still opening never stays small.

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

The table lists the main kinds; Line and Density take X and Y like a scatter plot (Line also
takes Color, Density also takes Bins).

| Setting | Scatter | Histogram | Box plot |
|---|---|---|---|
| Plot | Switch to any 2D kind | | |
| X | Number column | Number or text column | **Group by** (an optional text column) |
| Y | Number column | — | **Values** (a number column) |
| Color | Number (gradient) or text (one colour per category, with a legend) | — | — |
| Size | Number column | — | — |
| Bins | — | Number of bins, or Auto | — |
| Curve | — | Overlay a Normal or Log-normal fit | — |
| Guide | BOILED-Egg, Golden Triangle or LipE lines drawn behind the points | — | — |
| Fit line | Least-squares line; one per category when coloured by text | — | — |

A **bar chart** has **Group by** (a text column), **Values** (a number column, or (none) to
count rows), **Summary** (Mean, Median, Sum, Minimum, Maximum or Count) and **Error bars**.
A **radar plot** has **Axes**: tick the number columns to draw as spokes.

**Colours and sizes.** These rows appear only where they apply:

| Setting | Shown | What it does |
|---|---|---|
| Colors | Color is a number column; always on Density and 3D Landscape | The colour map the values run through: Viridis (the default), Plasma, Inferno, Magma, Cividis, Turbo, Rainbow, Jet, Red–Yellow, Blue–Red or Greys. **Reverse** runs it from high to low. |
| Color range | Color is a number column | The values at the two ends of the colour map. Values beyond an end take that end's colour. A box left empty (**Auto**) follows the data. |
| Size range | Size is a number column | The values given the smallest and the largest point. Values beyond an end take that end's size. A box left empty (**Auto**) follows the data. |
| Point size | Size is a number column | The smallest and largest point, in pixels (1–60; 4 and 18 to start). |

Each setting applies to Scatter, Line, 3D Scatter, Chemical Space Network and TMAP points
(Line takes no Size), and to the cores of a
Constellation Plot (where **Compounds** takes the place of Size range). The colour map also
applies to the Density counts and the 3D Landscape heights. A range end that is not a number goes back to
Auto. These settings are saved with the plot.

Hover over **Fit line** (or **Curve**) to see what was drawn: how many points there are, how
many rows were left out for missing values, and the fit's equation and R² (or μ and σ for a
fitted curve).

**Point cards** show a point's structure and values. Their settings are at the bottom of the
popup:

| Setting | Choices |
|---|---|
| Point cards | **On Hover** (the default): a card follows the mouse. **On Hover and Selected**: as well, each selected point keeps a card beside it. **Selected Points Only**: only the kept cards. **Off**. |
| Show | **Structure** (the drawing), **Label** (the row's first text column, such as its ID) and **Values** (the plotted columns: X, Y, colour and size). |
| More columns | Tick any other columns to show on the cards. |
| Drawing size | The largest drawing: Small, Medium or Large. A smaller molecule gets a smaller card, drawn at the normal bond length so it stays sharp. |

- Kept cards are laid out together: each takes the nearest free spot around its point, never
  on another card or another selected point, with a thin line back to its point when it had
  to move away. They move with the plot when you pan or zoom. When many cards share a small
  plot their drawings get smaller, and in a very crowded plot they show text only. A card disappears when its point is deselected, filtered out or out of view.
- Cards are kept for up to 12 selected points; select more and none are kept (they would
  cover the plot). Hover still works.
- The card settings are saved with the plot.

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
- **Hover** over a point to see its structure and values (see Point cards above).
- The plot shows only the rows the table shows, and redraws when the data, filters or columns
  change.
- Plots are saved with the document and reopen with it.
- Click a bar (bar chart) or a cell (density) to select its rows. On a radar plot, click a
  corner of a selected row's outline to select only that row.

### 3D plots

**Data ▸ Plot ▸ 3D Scatter** draws three number columns (X, Y and Z, set with the cog) in a
cube, on a white background. Each axis is scaled to the cube's side, so columns with very
different ranges share the space; the lowest and highest values are written at the ends of
each axis. **Color** takes a number or text column, **Size** a number column. **Drop lines**
draws a line from each point down to the floor, so heights read at a glance (above 20,000
points, only for the selected ones).

**Data ▸ Plot ▸ 3D Landscape** is an activity landscape (Peltason and Bajorath, 2007). Set **X** and
**Y** to the plane, usually the two columns of a Dimensionality Reduction map, and **Height** to the
activity. The surface at each point of a grid is the weighted mean of the compounds near it,
coloured by height; where no compound is near, the surface is left open rather than made up.
Smooth plateaus are regions where structure changes alter activity little; sharp peaks
beside valleys are activity cliffs. **Grid** sets the cells along each side (Auto picks from
the number of compounds). The compounds stand on the surface and select as on any plot.

The 3D plots use the same toolbar as the 2D ones:

| Button | Name | What it does |
|---|---|---|
| Arrow round a circle | Turn | Drag to turn the cube, right-drag (or Shift-drag) to move it, scroll to zoom, click a point to select it. |
| Dashed box | Box Select | Drag a rectangle to select the points inside it, as seen from where you look. |
| Loop | Lasso Select | Draw round points to select them, as seen from where you look. |
| Four corners | Reset View | Looks at the cube from the start position again. |
| Arrow into a tray | Export | Saves the view as a PNG image. |

- Shift adds to the selection; ⌘/Ctrl toggles. Hover over a point to see its card; point
  cards beside selected points follow them as you turn the cube.
- The plot is saved with the document, with the direction you were looking from.

### Constellation Plot

**Data ▸ Constellation Plot** shows the analog series in a set and how they relate
(Naveja and Medina-Franco, *Front. Chem.* 2019, 7, 510).

- **Cores.** Each compound is cut at its retrosynthetic (RECAP) bonds. Every fragment that
  keeps at least two thirds of the compound's heavy atoms is one of its cores, and so is the
  compound itself. Attachment points become hydrogens. This is the compound–core relationship
  (CCR) method of Naveja et al., *ACS Omega* 2019, 4, 1027.
- **Analog series.** Compounds that share a core, directly or through other compounds, form
  one analog series. Series are numbered from the largest.
- **The plot.** Each circle is a core shared by at least **Min compounds** compounds. Of cores
  that hold exactly the same compounds, only the largest is shown. Circles are placed by the
  cores' Morgan fingerprints (radius 2) on a t-SNE map, so similar cores sit together. A
  circle's size is its number of compounds. Grey lines join cores that share a compound.
- **Color.** Choose a number column and a **Summary** (Mean, Median, Minimum, Maximum or
  Count) to colour each core by its compounds' values. A core whose compounds have no value
  is grey. The Mean of a column of 0 and 1 is the share of actives.

| Setting | What it does |
|---|---|
| Structure | The structure column whose cores are found. |
| Map | t-SNE (the default), UMAP or PCA. Under five cores the map is PCA. |
| Min compounds | Cores shared by fewer compounds are left out (2 to start). |
| Core size | The share of a compound's heavy atoms a core must keep (67 % to start, as in the paper). |
| Series labels | Writes each series' number beside its largest core (up to 400 cores). |

- Cores are found in the background for the rows the table shows. The plot says "Finding
  cores…" until they are ready. Each structure's cores are remembered, so changing a filter
  or Min compounds only regroups and maps them again.
- Click a core, or box or lasso several, to select their compounds' rows. A core holding a
  selected row is ringed.
- Hover over a core to see its structure, its series, its number of compounds and its
  summary value.
- The settings are saved with the plot, and the cores are found again when the document
  reopens.

### Chemical Space Network

**Data ▸ Chemical Space Network** draws each compound as a point and joins the ones that
are alike (Maggiora and Bajorath, *J. Comput.-Aided Mol. Des.* 2014, 28, 795; Zwierzyna et al.,
*J. Comput.-Aided Mol. Des.* 2015, 29, 113).

- **Edges.** **Similarity** joins every pair whose Tanimoto similarity (Morgan fingerprints,
  radius 2) reaches the **Threshold**. **Matched Pairs** joins every matched molecular pair, as
  in Data ▸ SAR ▸ Matched Pairs.
- **Layout.** Each connected cluster is laid out by force, so tightly joined compounds sit
  together, and the clusters are placed in rows, largest first. Compounds with no edge
  (singletons) come after them, in rows; untick **Singletons** to leave them out.
- **Color and Size** work as on a scatter plot. **Edge color** colors each edge by its
  similarity, by the difference in the two compounds' Color values, or by their **SALI**,
  |difference| / (1 − similarity) (see [Activity Cliffs](#activity-cliffs-sali)): a large
  difference between alike compounds is an activity cliff.
- More than 50,000 pairs at the threshold make a network too dense to read, so the plot stops
  and says so: raise the threshold.

| Setting | What it does |
|---|---|
| Structure | The structure column. |
| Edges | Similarity or Matched Pairs. |
| Threshold | The similarity a pair needs (0.65 to start). Similarity edges only. |
| Edge color | Plain, Similarity, Activity Difference or SALI (of the Color column; SALI needs Similarity edges). |
| Singletons | Draws compounds with no edge. |

### TMAP

**Data ▸ TMAP** lays a large set out as a tree (Probst and Reymond, *J. Cheminf.* 2020,
12, 12). Each compound is joined to its most similar compounds (**Neighbours**, by the
**Fingerprint**); the tree keeps the shortest of those links that join every compound once, and
is laid out so that each compound sits on a branch beside its closest relatives. Series show as
branches; distances between far branches mean little.

- **ECFP4** (the default) compares Morgan fingerprints, radius 2, by Tanimoto, as the
  similarity tools do. Up to 10,000 compounds the neighbours are exact. Larger sets find them
  with MinHash and a forest of sorted prefixes, as TMAP does, then check them with the exact
  similarity: about nine in ten of the true ten nearest are found. 50,000 compounds take about
  12 seconds, 250,000 about a minute (and about 1.2 GB of memory).
- **MHFP6** is the fingerprint TMAP was designed with (Probst and Reymond, *J. Cheminf.* 2018,
  10, 66): each compound's substructures up to three bonds from an atom, written as SMILES,
  reduced to 256 MinHash values. Two compounds' similarity is the share of values they agree
  on, an estimate of the overlap of their substructures. Its values go straight into the
  forest at every size. Working out MHFP6 is much slower than ECFP4: about a minute for 50,000
  compounds and four for 250,000.
- If some compounds share no neighbour with the rest, there are several trees, placed side by
  side.

| Setting | What it does |
|---|---|
| Structure | The structure column. |
| Neighbours | How many most similar compounds each is first joined to (10 to start). |
| Fingerprint | **ECFP4** or **MHFP6**: what makes compounds neighbors. |

Both plots:

- are built in the background for the rows the table shows, saying "Finding neighbors…" until
  they are ready, and are built again when a filter changes;
- select rows as a scatter plot does (click, box, lasso), and the table scrolls to them; hover
  shows the row's card and how many edges it has;
- are drawn the same way every time, so a reopened document shows the same picture.

### MedChem plots

The Data ▸ MedChem menu has three plots. Each one first calculates any properties it needs that
are not already columns, then plots them. To calculate them again anyway, tick **Recalculate
existing columns** under **More Options**. Hover over **Plot** for what the plot shows.

| Plot | Axes | What it shows |
|---|---|---|
| BOILED-Egg | WLOGP against TPSA | The white predicts gut absorption; the yolk, brain penetration (Daina & Zoete 2016). |
| Golden Triangle | MW against logP | The region favoured for permeability and clearance (Johnson et al. 2009). |
| Lipophilic Efficiency (LipE) | Potency (pActivity) against logP | Lines of equal LipE = pActivity − logP. |

---

## 10. Chemistry: descriptors, pKa, standardizing, conformers and retrosynthesis

### Calculate Descriptors (Ctrl+Shift+D)

Chemistry ▸ Calculate ▸ Descriptors…

Pick a structure column, the rows (All, Shown or Selected), and the properties you want.
**Find** narrows the list; **Common Set**, **All** and **None** tick or untick many at once.
Each property becomes a number (or text) column. An existing column is never replaced:
calculating MW again adds "MW (2)".

| Group | Descriptors |
|---|---|
| Physicochemical | MW, Exact Mass, cLogP, cLogD 7.4, LogS (ESOL), TPSA, TPSA (N,O,S,P), Molar Refractivity, Labute ASA, Fsp3 |
| Ionisation | pKa (most acidic), pKa (most basic), Charge @ pH 7.4, % Neutral @ pH 7.4 |
| Drug-likeness | QED, SAScore, Ro5 Violations, Ro5 Pass, CNS MPO, BBB Score, AB-MPS |
| Counts | HBD, HBA, NH/OH Count, N+O Count, Rotatable Bonds, Heavy Atoms, Heteroatoms, Valence Electrons, Net Charge |
| Rings | Rings, Aromatic Rings, Aliphatic Rings, Saturated Rings, Spiro Atoms, Bridgehead Atoms |
| Atom Counts | C, N, O, S, P, F, Cl, Br, I |
| Complexity | Bertz Complexity, Balaban J, Hall-Kier Alpha, Kappa 1–3 |
| Partial Charges | Max/Min Partial Charge, Max/Min Abs Partial Charge |
| MOE-type (VSA) | SlogP_VSA1–12, SMR_VSA1–10, PEOE_VSA1–14, EState_VSA1–11, VSA_EState1–10 |
| 3D Shape | PMI 1–3, NPR1, NPR2, PBF, Radius of Gyration, Asphericity, Eccentricity, Spherocity Index, Inertial Shape Factor |
| Identifiers | Canonical SMILES, InChIKey, InChI, Formula, Murcko Scaffold, Generic Scaffold |

- **pKa-based values.** The Ionisation group, cLogD 7.4, CNS MPO, BBB Score and AB-MPS use
  chemsuite pKa (see [Predict pKa](#predict-pka)), which takes seconds per molecule. For CNS MPO and BBB Score, a compound with no basic site is scored as if its most
  basic pKa were 2.
- **SAScore** is the synthetic accessibility score of Ertl and Schuffenhauer (2009): 1 (easy
  to make) to 10 (very hard), from how common the molecule's fragments are in known compounds
  and a penalty for size, stereocentres, macrocycles and unusual ring fusions. Most drug-like
  compounds score 2-4.
- **3D shape descriptors** use a 3D structure when the cell has conformers, with all its
  hydrogens: they carry mass, and without them PMI and PBF can be a tenth to over half too
  low. A structure read without some of them (a pose with only polar hydrogens, say) gets
  the missing ones placed before it is measured; the table is not changed.
- **Every other descriptor** sees a 3D structure as drawn: without its explicit hydrogens or
  its conformers, so it gives what the same compound's 2D cell gives. **InChI** and
  **InChIKey** in particular keep a stereocenter or double bond the structure leaves undefined
  undefined, rather than reading it from the first conformer's geometry; the Calculator, Copy As
  and InChI export give the same. A structure read from a file, taken from a protein or docked
  has its stereo set from its own geometry, so its InChI says which isomer it is.
- **MOE-type (VSA)** are Labute's binned surface-area descriptors (2000), with MOE's names and
  bins, computed by RDKit (MOE itself is not needed; values can differ slightly from MOE's).
  Each sums the surface area of the atoms whose property (Crippen logP, molar refractivity,
  Gasteiger charge or E-state) falls in one range; VSA_EState sums E-states over atoms binned
  by surface area. A column's tooltip gives its range.
- **Fingerprints.** Tick any of the fingerprints the similarity tools use (ECFP4, ECFP6, FCFP4,
  RDKit paths, atom pairs, topological torsions, MACCS keys). Each becomes a text column
  ("ECFP4", "MACCS Keys"…) holding the bits in hexadecimal, bit 0 first (the high bit of the
  first byte), ready to export for modelling. A row with no valid structure is left blank.
- **Cancelling** keeps the rows already calculated; a fingerprint is written only when every
  row has one.

### Predict pKa

Chemistry ▸ Calculate ▸ pKa… predicts pKa values with **chemsuite pKa**, chemsuite's own model:
the Uni-pKa method (Luo et al., JACS Au 2024) on the Uni-Mol network (Zhou et al., ICLR 2023),
trained on pKaHub's aqueous pKa data. It gives a free energy for every protonation state
(microstate) of the whole molecule, so sites that affect each other, as in zwitterions and
polyamines, are handled together. On held-out test sets (SAMPL6–8 and Novartis acids and
bases) its error is RMSE 0.81 pKa units.

- **The model** (88 MB) comes with chemsuite and needs no network. The first prediction
  unpacks it once into chemsuite's cache folder (about 190 MB on disk), and every worker
  process then shares that one copy rather than loading its own.
- **Speed.** About 1 s for a small molecule and 4–20 s for a drug-sized one, which has more
  protonation states. Rows run in parallel across the worker processes, and progress counts
  rows as they finish.
- **Input.** A structure column and the rows. Salts are reduced to the largest fragment and
  charges are neutralized before prediction, so a sodium salt gives its parent acid's values.
- **pH.** Charge, % neutral, cLogD and the major species are given at this pH (default 7.4).
  The protonation states considered start at net charges −1 to +1, and widen, up to −4 and +3,
  while more than 1% of the molecules sit at the end at this pH. So a pKa between charges
  that are not populated at the pH (a dication's, near pH 7) may be listed only at another pH.
- **Columns** (choose any):

  | Column | Meaning |
  |---|---|
  | pKa (most acidic) | The lowest acid pKa. Blank without an acid step. |
  | pKa (most basic) | The highest pKa of a protonated base (the strongest base). Blank without a base step. |
  | pKa Values | Every pKa, lowest first, each marked acid or base, for example glycine's "2.11 (acid), 9.76 (base)". |
  | Charge @ pH, % Neutral @ pH | Mean net charge, and the share of molecules with no site ionized. Zwitterions count as ionized; a permanent charge (a quaternary nitrogen) does not. |
  | cLogD *pH* | Crippen cLogP of the neutral parent plus log10(fraction neutral). |
  | Major Species @ pH | The most populated protonation state, as a structure column. |

  Each run writes new columns and never changes an earlier run's. When a name is already
  taken, all the columns of the run get the same suffix, (2), (3)…, so they read as a set:
  `pKa (most acidic) (2)` and `cLogD 7.4 (2)` come from the same run.

- **Macroscopic pKa values.** Each pKa is a step in net charge, between all the states of one
  charge and all those of the next. It is marked by the atom that loses the most hydrogen
  across the step: a base when that atom is positively charged before it loses it, an acid
  otherwise. So glycine's carboxylic acid is an acid and its ammonium a base.
- **Limits.**
  - A molecule with more than 48 protonation states is left blank; widening stops before it
    would pass 48, so a polyphosphate may lack its last pKa.
  - States that protonate an oxygen or sulfur (a carbonyl, an alcohol) are left out: their
    pKa values lie below about −2, where the model was not trained.
  - Because only the neutral form is assumed to partition, cLogD of zwitterions and of
    strong acids and bases comes out far lower than measured.
  - Treat the values as screening estimates, not measurements.
- **Cancelling** keeps the rows already predicted.

### Sites of Metabolism

Chemistry ▸ Calculate ▸ Sites of Metabolism… ranks the atoms of each structure by how likely a
cytochrome P450 isoform is to metabolize them. The model is chemsuite's own: a gradient boosting
machine over 34 cheap atom descriptors (element, hybridization, ring and valence context,
Gasteiger charge, E-state and Crippen contributions) plus the isoform, trained on the measured
sites of 679 compounds across the nine isoforms. It comes with chemsuite, needs no network, and
scores a molecule in about a millisecond.

- **Isoforms.** CYP3A4, 2D6, 2C9 and 1A2 are on by default; 2A6, 2B6, 2C8, 2C19 and 2E1 are under
  More Options. Each isoform gets its own columns. The model scores them all at once, so asking
  for more of them costs almost nothing.
- **Input.** A structure column and the rows. Hydrogens are not sites: explicit ones are ignored,
  and only heavy atoms are scored. Salts are not reduced for you — run Standardize first if your
  structures have counter-ions, because the model was trained on single, neutralized fragments.
- **Columns** (choose any):

  | Column | Meaning |
  |---|---|
  | Top Site *(isoform)* | The highest-scoring atom, named by its element and its index in the structure, for example `C7`. Atoms that score identically, such as the two methyls of an isopropyl group, are reported by the lower index. |
  | Score *(isoform)* | That atom's score, 0 to 1. |
  | Sites *(isoform)* | The three highest-scoring atoms, best first, each with its score: `C7 0.59, C18 0.04, C19 0.04`. |
  | Applicability *(isoform)* | The share of the top atom's surroundings the model saw while training, 0 to 1. Below 1 means the molecule is unlike the training set at that atom, so trust the ranking less. |

  Each run writes new columns and never changes an earlier run's. When a name is already taken,
  all the columns of the run get the same suffix, (2), (3)…

- **How accurate it is.** Measured over ten random splits that keep every isoform of a molecule
  on one side of the split, so the numbers describe a molecule the model has never seen:

  | | top-1 | top-2 | top-3 | MCC |
  |---|---|---|---|---|
  | All isoforms | 0.59 | 0.73 | 0.81 | 0.48 |
  | CYP3A4 | 0.53 | 0.67 | 0.76 | 0.43 |
  | CYP2D6 | 0.59 | 0.72 | 0.81 | 0.44 |
  | CYP2C9 | 0.61 | 0.73 | 0.83 | 0.47 |
  | CYP1A2 | 0.61 | 0.74 | 0.84 | 0.48 |

  So the measured site is the top-ranked atom about three times in five, and among the top three
  atoms about four times in five. CYP3A4 is the weakest isoform even though it has the most data:
  it is the promiscuous one, with the largest substrates.

- **What counts as a site.** Only *primary* sites of metabolism — the dominant route reported for
  a compound. Published figures for other tools usually count secondary and tertiary annotations
  as hits too, which raises top-1 by about 0.10, and most of them are measured on splits that put
  one isoform of a molecule in training and another in the test set, which raises it further. The
  numbers above do neither, so they are lower than published ones and not comparable with them.
- **Limits.**
  - Scores rank the atoms of one molecule. They are not probabilities you can compare between
    molecules, and they do not say whether the compound is a substrate at all.
  - The isoform matters less than the molecule: the same atom often leads for every isoform.
  - Phase II metabolism (glucuronidation, sulfation), non-P450 enzymes and the structures of the
    metabolites themselves are out of scope.
  - Treat the ranking as a screening aid, not a measurement.
- **On the structures.** When a run finishes, the sites of the most prescribed isoform you
  asked for (CYP3A4 first, then 2D6, 2C9, 1A2, else the first isoform of the run) are marked on
  that structure column's drawings: the ranked Sites column if you asked for it, otherwise the
  Top Site column. That is how you see which atom `C7` means, since the table does not number
  atoms; hover over the cell for a larger drawing with the same atoms marked. Right-click the
  structure column's header for **Highlight Atoms ▸** to mark another isoform's column
  instead, or Data ▸ Search ▸ Clear Highlighting to stop.
- **Cancelling** keeps the rows already found.

### Standardize

Chemistry ▸ Standardize prepares structures for docking, or for any work that needs
one consistent form of each molecule. It runs the usual preparation steps in this order, the
order of the ChEMBL structure pipeline and RDKit's MolStandardize:

| Step | Default | What it does |
|---|---|---|
| Clean up | On | Removes explicit hydrogens, disconnects metals bonded to N, O or F, draws functional groups one way (for example a charge-separated pyridone becomes the pyridone), and puts charges on the strongest acid. |
| Strip salts | On | Keeps the largest organic fragment: a hydrochloride gives the free base, a hydrate the parent. |
| Neutralise | On | Adds or removes hydrogens so that charges cancel. Quaternary nitrogens stay charged. |
| Remove isotopes | Off | Deuterium becomes hydrogen, ¹³C becomes carbon. |
| Canonical tautomer | Off | One consistent tautomer per molecule (RDKit), so duplicates collapse. It is not necessarily the most stable or the bioactive form. If the canonical tautomer would lose any drawn stereo, the drawn tautomer is kept. |
| Protonate at pH | Off | The most populated protonation state at the pH (default 7.4), from chemsuite pKa (see [Predict pKa](#predict-pka) and its limits). It is never a carbanion. Kept counter-ions stay as they are. Without the model, Standardize cannot run with this step on. |
| Form salt | Off | Adds the counter-ion you pick and balances the charge. Runs last, so leaving **Strip salts** on switches one salt for another. |

- **Stereo** is never invented or lost. Each step keeps drawn stereocentres and double-bond
  geometry. Over the 3,024 FDA-approved structures, with every step on, no drawn stereo was
  lost. Undefined stereocentres and double bonds are listed in the notes, so you can
  enumerate or fix them before docking.
- **Output.**
  - **Write to.** **New Column** (the default) writes "*column* (Standardized)", replacing
    it if it is already there, so running again updates it. **Source** overwrites the
    structures in the source column; a row that cannot be standardized keeps its structure.
  - A **Standardization Notes** column is always written (replaced on a rerun). It says what
    changed in each row, for example "salts stripped; neutralised; protonated at pH 7.4; net
    charge +1; 1 undefined stereocentre".
  - A **Salts** column is written whenever **Strip salts** removed something, holding the
    fragments it took (`Cl`, `[Na+]`), so a counter-ion is recorded rather than lost. Rows
    that had nothing stripped are blank, and the column is left out when no row had a salt.
  - Undo restores the table in one step.
- **Form salt** takes the parent to its salt form, the opposite of stripping.
  - **The counter-ion** list holds the common ones: hydrochloride, hydrobromide, sulfate,
    phosphate, mesylate, tosylate, besylate, trifluoroacetate, acetate, citrate, fumarate,
    maleate, succinate, tartrate, oxalate, sodium, potassium, lithium, calcium, magnesium,
    zinc, ammonium, choline, meglumine and tromethamine. **×** sets how many go on.
  - **The charge is balanced.** An acid's ion (chloride, mesylate) protonates the strongest
    base in the molecule; a base's ion (sodium, calcium) deprotonates the strongest acid.
    Charge the molecule already carries counts, so a quaternary ammonium just gains its ion.
  - **The strongest site** is chosen by group: for bases, guanidine, then amidine, aliphatic
    amine, aromatic nitrogen, aniline; for acids, sulfonic acid, then carboxylic acid,
    tetrazole, acyl sulfonamide, phosphonic acid, phenol. It is the usual ordering of these
    groups, not a pKa calculation; [Predict pKa](#predict-pka) does that.
  - **A multiply charged ion repeats the molecule** when one copy cannot balance it, so
    ethylamine and sulfate give the 2:1 salt, and acetic acid and calcium the calcium diacetate.
  - **A molecule with no suitable site is left alone**, and its notes say "no site for a
    hydrochloride salt".
- **Next steps.** For docking, follow with [Generate Conformers](#generate-conformers) for 3D
  coordinates. Standardized structures carry no explicit hydrogens; docking programs add them.
- **Cancelling** keeps the rows already standardized.

### Generate Conformers

**Chemistry ▸ 3D ▸ Generate Conformers…** generates 3D conformers with one of two **Generators**:

- **CONFORGE (CDPKit)**, the default: builds conformers from fragment and torsion libraries
  and scores them with MMFF94s (Seidel et al., *J. Chem. Inf. Model.* 2023, 63, 5549). It
  comes from the `cdpkit` package (LGPL-2.1), which installs with chemsuite wherever CDPKit
  ships one: Windows with Python 3.12, Linux on x86-64, and Macs with Apple silicon. Elsewhere
  (Windows with Python 3.11 or 3.13, Intel Macs) the choice is greyed out and ETKDG is used.
  - It reproduced crystal poses best of the generators tested: on 3,354 drug-like PDB
    ligands (Wang et al., *J. Chem. Inf. Model.* 2023, 63, 6525), the closest of 10
    conformers was 0.84 Å from the crystal pose on average, against 0.99 Å for ETKDG with
    its default settings.
  - CONFORGE gives the same conformers every run, uses its own force field (**Optimize with**
    and **Seed** are greyed out), and always applies the energy window if it is ticked.
  - Macrocycles are sampled with 200 trial conformers each, which takes several seconds per
    molecule (an 18-membered lactam: about 7 s).
- **ETKDG (RDKit)**: distance geometry, then optimisation with MMFF94 (UFF where MMFF has no
  parameters). Where it is the default (no CONFORGE), the **Energy window** starts unticked:
  in the same benchmark the window only dropped conformers close to the crystal pose.

- **Input:** the structure column and **All**, **Shown** or **Selected** rows.
- **Method:** the **Generator**, **Conformers** per molecule, **Distinct above** (the RMSD for
  two conformers to count as distinct), the force field (**Optimize with**), and an **Energy
  window** (tick it to drop conformers more than that many kcal/mol above the lowest). The
  random **Seed** is under **More Options**.
- **Output.** The conformers go into a new 3D structure column, lowest energy first, with
  columns **Conformers** (the count) and **Conformer Energy Range (kcal/mol)**.
  - Each conformer keeps its energy (kcal/mol, from MMFF94, UFF or CONFORGE's MMFF94s), its
    **ΔE** above the molecule's lowest conformer and its Boltzmann **population** at 298.15 K
    over the conformers kept. They are saved with the document, listed in the
    [Browser](#browser), written by [Export Conformers](#export-conformers-sdf) and carried
    onto superposed poses.
  - Energies compare conformers of the same molecule only. The populations come from
    gas-phase force-field energies, which favor intramolecular hydrogen bonds, so read them as
    a rough guide.
  - ETKDG with **Optimize with** set to none stores no energies. Conformers made by an
    earlier version of chemsuite, or read from a file, have none either; generate them again
    to get them. Standardize and Predict pKa drop the energies of a structure they change.
  - Running again writes new columns ("SMILES 3D (2)"); earlier conformers are kept.
- **Cancelling** keeps the molecules already done.
- A structure ETKDG cannot embed, or not within a minute (it tries twice: from its own
  geometry, then from random coordinates), gets no conformers; the message at the end says
  how many structures could not be embedded.

### Export conformers (SDF)

Writes every conformer as its own SD record, carrying the row's other columns. Choose **All**,
**Shown** or **Selected** rows (Selected when rows are selected), then **Export** asks for the
file.

**Hydrogens** sets what is written: **As Stored** (the default) writes each structure as the
table holds it, **Add** writes every hydrogen, placing missing ones on each conformer, and
**Remove** writes heavy atoms only (deuterium and hydrogens that define stereo stay). The
table itself is not changed. Conformers made by Generate Conformers already carry every
hydrogen.

Each record has a **Conformer** field numbering it (when the row has several). Conformers that
carry energies also get **Energy (kcal/mol)**, **Relative Energy (kcal/mol)** (ΔE above the
lowest), **Population (%)** and **Force Field**.

### Superpose

**Chemistry ▸ 3D ▸ Superpose…** moves molecules onto a reference in 3D, or lines up each
row's conformers on its first.

- **Input:**
  - **Structures:** the column to superpose (a 3D column, if there is one, at first), and
    **All**, **Shown** or **Selected** rows.
  - Every 3D conformer a row holds is tried, and the one that fits best is kept. A row with
    only a 2D structure gets conformers first, as Generate Conformers makes them (CONFORGE
    where it is installed, else ETKDG; how many is under **More Options**).
  - **Superpose:** **Onto a reference**, or **Each row's conformers onto its first**.
  - **Reference:** **Current row** (the table's current row: click another row while the
    dialog is open to change it), or a ligand of a protein in the document, in its bound pose
    ("1ATP ▸ E/ATP 355"). For a row, **Reference from** is the column its 3D structure is read
    from; its first conformer is used.
- **Method:**

  | Method | Fits | Scores |
  |---|---|---|
  | Substructure (MCS) | the largest common substructure of each row and the reference, atom on atom | **Core RMSD (Å)**, **Core Atoms** |
  | Substructure (core) | a core you type (SMILES or SMARTS) or draw (**Draw…**) | **Core RMSD (Å)**, **Core Atoms** |
  | Atoms (O3A) | atoms paired by type without a common core (Open3DALIGN; Tosco et al., *J. Comput.-Aided Mol. Des.* 2011, 25, 777) | **O3A Score** (higher is better), **RMSD (Å)** of the paired atoms |
  | Shape | the overlap of the two molecules' volumes (Gaussian atoms, as ROCS; RDKit's shape alignment, from PubChem) | **Shape Tanimoto** (0–1) |
  | Shape + Features | volume plus pharmacophore features: donors, acceptors, rings, charges | **Combo** (0–2: the two Tanimotos added), **Shape Tanimoto**, **Color Tanimoto** |
  | Field (Electrostatics + Shape) | the overlap of the two molecules' electrostatic potentials, with shape | **Field Similarity** (−1 to 1), **Shape Tanimoto** |

  - Symmetric cores are tried every way round, so a phenyl fits whichever way its atoms are
    numbered.
  - **Align on** (Atoms, Shape, Shape + Features and Field; optional): a substructure you type
    (SMILES or SMARTS, at least 3 atoms) or draw. Each conformer is fitted atom on atom on it,
    and the method scores that pose and keeps the best-scoring conformer; O3A pairs the
    other atoms with those pairs held. It adds a **Core RMSD (Å)** column. A row without the
    substructure is fitted as if none were given, with Core RMSD left blank.
  - Field: each atom's partial charge (Gasteiger, or MMFF94 under **More Options**) gives a
    potential, and the overlap of the two potentials, as a Carbó index, is the field
    similarity. It starts from the best shape fits of the three best-fitting conformers,
    and from each turned half over about its axes, so it can tell a pyridine's nitrogen from
    its carbons where shape cannot.
  - **More Options:** the MCS settings (rings match rings only, complete rings only, compare
    elements, a time limit per row), the O3A atom typing (MMFF94 or Crippen), the **Feature
    weight** (Shape + Features) and **Field weight** (Field), each the share of the
    optimization given to features or field, and **Conformers for 2D rows**.
- **Output** (one undo step):
  - a new structure column, "<column> superposed", with each row's best conformer moved onto
    the reference;
  - **Best Conformer** (1 is the first stored) and a column for each score;
  - **Pose ΔE (kcal/mol)**: the kept conformer's energy above the lowest conformer of its
    molecule, a rough measure of strain. It is written when the conformers carry energies
    (from [Generate Conformers](#generate-conformers), or made here for 2D rows) and left
    blank for those that do not;
  - the reference itself, kept with the document so the [Browser](#browser) can show it
    with each superposed row. Rename or copy the column and the link stays.
  - A row with no structure, or that cannot be fitted (the core is not in it), is left blank;
    the message says how many.
- **Each row's conformers onto its first** fits every conformer of a row onto its first, on
  heavy atoms (**Heavy atoms only**, under **More Options**), and writes **Mean Conformer RMSD
  (Å)** and **Max Conformer RMSD (Å)**.
- **Cancelling** keeps the rows already done.

### Browser

The Browser is a panel that shows one row's structure at a time, drawn in **2D** or in
**3D**, starting at the table's current row. Open it with **Data ▸ Browser (Ctrl+Shift+B)**,
at the bottom of the Data menu: a new Browser opens in 2D, and an open one is brought forward
as it is. The **2D | 3D** switch at the left of its header changes only
how the row is drawn: the row, the structure column, the rows browsed and the caption stay as
they are, and a 3D view you turned is still turned when you come back to it.

- **Browsing rows:** ‹ and › under the view (or ← and →, Home and End) step through the rows
  in the table's order; type a number in the row box to jump to it. The table's current row
  moves along (the selection is left alone), and clicking a row in the table shows it in the
  viewer. Which rows are browsed is set in the options menu (gear):
  - **Browse Selected:** only the selected rows.
  - **Browse Visible** (on by default): only the rows the filters and search show.
  - Both on: the selected rows that are shown. Both off: every row in the table.
  - Deselecting the row shown while browsing the selection moves on to the next one.
- **Selected** (the tick box under the view, or Space) adds the row to the selection or takes
  it out.
- **Caption:** column values drawn over the structure, in the corner, as the table shows
  them. Choose the columns in the options menu (an ID or name column is chosen at first, if
  there is one). The caption is the same in 2D and 3D, and is part of an exported picture.
- **Structure column:** chosen in the options menu.
- **Export** saves the view (with its caption) as a PNG.
- The panel, its settings, 2D or 3D, the caption, the row and conformer shown and the 3D view
  are saved with the document. A Structure Viewer or 3D Viewer panel saved by an earlier
  version reopens as a Browser.

**In 2D** the structure is drawn as in the table, larger; ↑ and ↓ step through the rows too.

**In 3D:**

- A row with only a 2D structure gets 3D coordinates in the background (one ETKDG conformer);
  they are not written to the table. They are made ahead for the rows either side of the one
  shown, in 2D as well, so switching to 3D or stepping on usually finds them ready.
- **Style:** Ball and Stick, Stick, Spacefill or Wire. **Hydrogens:** All, Polar (on atoms
  other than carbon) or None. **Color:** by element, or each conformer's carbons in a colour of
  its own (with Overlay Rows, each row's). These options show in the options menu in 3D.
- **Conformers:** for a row with several, ‹ and › in the panel's header (or ↑ and ↓) step
  through them; the counter shows which one, and its ΔE above the lowest when the conformers
  carry energies. **Overlay** draws every conformer at once as thin sticks with polar
  hydrogens. **Align Conformers** (on by default) fits them onto the first on heavy atoms.
- **Conformer Energies** (on by default, beside ›): for conformers made by
  [Generate Conformers](#generate-conformers), a list under the view gives each conformer's
  number, **ΔE (kcal/mol)**, **Population (%)** (Boltzmann, 298.15 K) and **Energy
  (kcal/mol)**. It follows the view: click a line to show that conformer; stepping moves the
  highlight. Select several lines (Shift- or Ctrl-click; ⌘-click on macOS) to overlay just those. Click a
  heading to sort; drag the boundary above the list to resize it. The header tooltips name
  the force field.
- **Overlay Rows** draws every row being browsed at once (all, shown or selected, as
  **Browse Selected** and **Browse Visible** set; up to 200), as thin sticks, each at the
  conformer shown, or all of its conformers with **Overlay** on too. On a column written by
  [Superpose](#superpose) it shows the poses a series takes on the reference; **Color ▸
  Carbons by Conformer or Row** gives each its own carbon color. Rows with only a 2D structure
  are left out, except the one browsed to. Stepping through rows keeps the view.
- **Show Reference** (on by default, in the options menu): a row written by
  [Superpose](#superpose) is drawn with the reference it was fitted to, as thin sticks with
  blue carbons.
- **Mouse:** drag to turn, right-drag, middle-drag or Shift-drag to move, wheel to zoom. Pointing
  at an atom names it (for example "N4"). Ctrl-click (⌘-click on macOS) an atom to add or remove
  its row from the selection (with Overlay Rows, the row that atom belongs to).
- **Measure** (ruler): click two atoms for a distance, three for an angle, four for a dihedral; a
  fifth click starts again.
- **Each row browsed to** is turned square on (its longest direction along the panel's longer
  side), centred and fitted to the panel, and it stays fitted as the panel is resized until you
  turn, move or zoom it. Stepping through conformers keeps the view, so they can be compared.
- **Reset View** turns the molecule square on and fits it again.

3D drawing, here and in the 3D plots, uses each system's own graphics: Metal on macOS,
Direct3D on Windows and OpenGL on Linux. The 3D view is only created when you first switch to
it, so the 2D view works on any computer. On Windows, if the graphics card cannot draw 3D
(common on remote desktops, virtual machines and VDI), Windows draws it in software instead:
slower, but it works. chemsuite checks the graphics when it starts; where 3D cannot be drawn,
the viewer stays in 2D, the 3D switch and the 3D plot types are greyed out, and their
tooltips give the reason. Settings ▸ General ▸ 3D graphics can force software drawing
(Windows and Linux) or turn 3D off.

### Macromolecules

The **Macromolecule** button on the menu bar (Ctrl+Shift+R; on macOS, **Macromolecule ▸
Macromolecule Viewer**) opens a workspace for macromolecules —
proteins, nucleic acids, glycans and their complexes: a window of its own, titled after the
document, with its own menu bar (on macOS the menu bar at the top of the screen switches to the
viewer's menus while it is in front). It is not docked beside the table; it closes with its
document and reopens with it, in the same size and place. The table is only the way in:
everything the viewer does (preparing, minimizing, docking, rescoring) is in its own menus and
the right-click menus:

- **File:** **Open…** (Ctrl+O; PDB or mmCIF files, gzipped or not, and session files,
  `.cspx`), **Fetch from PDB…** (Ctrl+Shift+O; download an entry by its ID, such as 1ATP),
  **Save…** (Ctrl+S), **Save Session…** (Ctrl+Shift+S), **Save Selection…** (PDB or mmCIF),
  **Export Image…** (Ctrl+E, below) and **Clear Session…**. Ctrl+W closes the window.
  - **Save…** writes every shown structure to one PDB or mmCIF file. One structure is written
    whole (every model and alternate location); several are written as one model, and a
    chain whose ID another structure took first is renamed (the status line says how). A
    selection across several structures is written the same way.
  - **Save Session…** writes the whole workspace to a session file: every structure, how each
    looks, saved selections, the focused ligand, labels, measurements, the view (background
    and depth cue), and the preparation reports, poses and records. **Open…** on a session
    file brings it all back. With structures already open it asks: **Replace** them (their
    ligand rows in the table are kept; Undo takes the session's structures out and brings
    back those it replaced) or **Add** the session's beside them (a structure whose ID is taken
    gets a new one, such as 1ATP-2, everywhere the session names it, and a saved selection
    whose name is taken gets `_2`). A session holds no table rows: **Tools ▸ Ligands to
    Table** makes them.
  - **Clear Session…** removes every structure, with the focus, saved selections, labels and
    measurements: an empty workspace. When ligand rows in the table name the structures you
    choose **Delete Rows** or **Keep Rows**, once for all of them. One Undo brings the
    structures (and any rows deleted) back, looking as they did; the focus, saved selections,
    labels and measurements are not brought back.
  - Saving writes a temporary file first, so a save that fails leaves the old file as it was.
- **The structures list:** right-click a structure for **Rename…** (the name shown; the ID the
  table's rows link to stays), **Duplicate** (a copy, `name-copy`, looking the same and keeping
  its preparation record) and **Delete…** (or press Delete). When ligand rows in the table name
  it you choose **Delete Rows** or **Keep Rows**; either way one Undo brings it all back. A
  kept row whose structure is gone shows its **Protein** cell in red in the table (its tooltip
  names the missing structure), and says so on the viewer's status line when you make it
  current.
- **Edit:** **Undo** and **Redo** (Ctrl+Z, Ctrl+Y: the document's last step, named in
  the menu, whether made here or in the table), **Copy** (Ctrl+C: the selected atoms as PDB
  text), **Copy Selection As ▸ PDB / mmCIF**, **Copy Image** (Ctrl+Shift+C: the view at twice
  its size, labels included), and the edits below: **Delete** (Delete or Backspace),
  **Element ▸**, **Charge ▸**, **Rename Residue…**, **Rename Chain…** and **Renumber
  Residues…**.
- **Select:** **All** (Ctrl+A), the usual picks (ligands, binding site, protein, backbone,
  side chains, waters, ions…), **Chain ▸**, **Around ▸**, **Saved Selections ▸**, **Save
  Selection…** and **Clear Selection**.
- **Display:** **Show**, **Hide**, **Color** and **Label** (on the selection), **Hydrogens**
  (on every protein) and **Surface ▸** (below).
- **Tools:** **Prepare…**, **Inspect…**, **Minimise…**, **Dock…** and **Rescore Poses…**
  (below), **Measure ▸**, **Ligands to Table** and **Interactions to Table** (below).
- **View:** **Frame Selection** (F), **Reset View**, **Ligand ▸** (**Focus ▸** a ligand in
  its pocket, **Pocket Surface ▸**, **Clear Focus**), **Interactions**, **Interaction Types ▸**,
  **Background ▸** (the
  theme's, white, grey or black), **Orthographic** (parallel projection: no perspective, so
  sizes compare across depth), **Depth Cue** (fade what lies further back; off draws everything
  fully), whether the **Structures List** shows, **Show Table** (brings the document's table
  forward) and **Full Screen**. Background, projection and depth cue are saved with the
  document.
- **Sequence** (pinned at the far right of the menu bar; on macOS, the last menu): **Sequence
  Viewer** (Ctrl+Shift+L) opens the sequence in a window of its own, and **Rows ▸** chooses
  which of its rows for ligands, glycans, additives, ions and waters show (below).

**Figures.** **File ▸ Export Image…** draws the view again at the size you ask for, rather
than enlarging the screen: give the width in pixels, centimetres or inches with a resolution
(300 dpi to start; it is also written into the file), and the height follows, as the picture
keeps the view's shape so labels and measurements land where they are on screen. **See-through
background** leaves the background out (PNG and TIFF) for a coloured slide. PNG, TIFF and JPEG;
at most 12,000 pixels a side. For a paper, a white background with Depth Cue off and
Orthographic on is a common choice.

Proteins are saved inside the document, with how they look (colours given by hand
included), saved selections, the focused ligand, labels, measurements and the view. These
belong to the document, not the window: closing the viewer keeps them, it opens again as you
left it, and they are saved even while the viewer is closed. Changing how proteins look marks
the document as changed; selecting atoms or turning the camera does not. Closing a document
that holds only proteins still asks to save it. Opening or
fetching one adds only the protein: the table is not touched. Its ligands' bond orders are worked out as it loads (asking the PDB's
chemical component dictionary when needed), and a ligand whose bond orders had to be guessed is
marked ⚠ in the structures list. Ligands become table rows only when you ask:

- **Tools ▸ Ligands to Table** adds a row for each ligand of every protein that has none in
  the table yet. Cofactors (ATP, NAD, FAD, haem, SAM…) get a row only when they are the
  protein's only ligand; **Add Ligand to Table** in the right-click menu adds one.
- **Tools ▸ Interactions to Table**: for each ligand row (the selected rows, or all), how many
  hydrogen bonds, salt bridges, π-stacking, cation–π, hydrophobic, halogen-bond, metal and
  water-bridge contacts its ligand makes, one number column each, and a **Contacts** column
  naming the residues and how (`MET793 HB; LEU788 X`). With no ligand rows yet, it adds them
  first. A row's own pose is used when it has one. It runs in the job queue; cancelling it
  keeps the rows finished, and running it again updates the same columns.

| Column | What it holds |
|---|---|
| Protein | Which protein the ligand belongs to (a second copy of 1ATP becomes 1ATP-2) |
| Ligand ID | Chain, residue name and number, such as `E/ATP 355` |
| Ligand | The ligand as a structure, in its pose from the file (3D kept) |
| Bond Orders | Where the bond orders came from (the first that fits): **SMILES** (one the row already had), **File** (the file's own component tables, as PDB mmCIF files carry), **CCD (downloaded)** (the PDB's chemical component dictionary, fetched once per ligand and kept, when downloads are allowed), **CCD template** (bundled with chemsuite) or **Guessed from geometry**, which is worth checking. A guessed ligand is marked ⚠ in the structures list |

Make a ligand row current in the table and the viewer shows that ligand in its pocket. When
the row's **Ligand** is 3D and no longer sits where the file has it (a docked pose, or a ligand
moved in the sketcher), the viewer shows the row's own pose instead: in sticks with its
pocket and interactions, and the ligand of the file it stands in for faded beside it for
comparison. Editing the current row's ligand updates the pose at once. The pose is not saved
with the viewer; making the row current brings it back.

**The viewer.** The structures list (tick a protein to show or hide it; click a chain or ligand to
select it) beside the 3D view. The sequence has a window of its own (**Sequence ▸ Sequence
Viewer**): a row for each chain, its letters following the full sequence (residues missing
from the model show faint), and under the chains one row each for the ligands, glycans,
additives, ions and waters, every residue in a cell of its own named by its residue name
(`HOH`, `NO3`), numbered above. So one water among hundreds can be found, selected and
deleted. Point at a cell to light the residue up in 3D and read its name and number in the
window's status line. The window closes with the viewer and is saved with the document, open
or not, with the rows you chose.

- **Selecting.** Click selects an atom; on a cartoon, where no atom is drawn, it selects the
  residue. Double-click selects the residue, Alt-click (Option-click) the chain; Shift adds and
  Ctrl (Cmd) toggles; click empty space or press Esc to clear. In the Sequence window click or
  drag across cells, click a row's name to select the whole row (a chain, or every water), or
  use its keys: the arrows move a cursor, Shift+arrows select from where it started, Return
  selects the residue under it, Home and End jump to the ends of a row. In the viewer, Ctrl/⌘+A selects every atom, Ctrl/⌘+C copies the selection as PDB and F frames
  it (a text field keeps its own). The structures list highlights a chain or ligand that is wholly selected and puts
  one partly selected in bold. Type in the selection field (point at it for the syntax): `chain A` (`chain ''` for a
  chain without a name, as docking and MD files often have; the list calls it unnamed),
  `segi PROA` (an MD file's segment: chains told apart only by segment are listed apart, as
  `Chain PROA`, or `Chain PROA/A` when the chain has a name too; a PDB save keeps the
  segments, an mmCIF save does not),
  `resi 10-50`, `resn HIS+ASP`, `name CA`, `ligand`, `water`, `byres (polymer within 5 of
  ligand)` (whole residues near a ligand; `byres` takes what follows it, hence the brackets),
  `and`, `or`, `not`. The **Select** menu has the usual picks (each shows its expression in the
  field), **Around ▸** (the residues or atoms within 3, 5, 8 or 10 Å of the selection, added to
  it) and **Save Selection…**; saved selections are used as `$name`, and `$selection` is what
  is selected now.
- **A protein opens** with its chains as wire, its ligands and ions as sticks and balls, its
  glycans and additives as sticks, and its waters and non-polar hydrogens hidden. Ligand
  carbons are **purple** and ions amber, so a ligand reads at once against the pocket around
  it.
- **Show, Hide, Color** act on the selection: wire, sticks, ball and stick, spacefill,
  cartoon, surface; color by element, chain, secondary structure, B-factor, hydrophobicity or
  a color you choose.
- **Bond orders** show in sticks and ball and stick: a double bond is two thinner sticks side
  by side, a triple bond three, and aromatic rings are drawn with alternating single and double
  bonds (a Kekulé form; a pyrrole-type NH, O or S takes no double bond). Wire stays one line
  per bond.
- **Label ▸ Residues / Atoms** labels what is selected (a residue on its Cα, a ligand in its
  middle); **Clear Labels** removes them. Labels stay legible: nearer ones win where they would
  overlap, faded and sliced-away atoms lose theirs, and at most 300 show at once.
- **Hydrogens** (all, polar or none) applies to every protein. A hydrogen is drawn as the
  atom it is bonded to is shown, so showing
  a residue as sticks shows its hydrogens too.
- **Surface ▸** sets how the surfaces of the proteins with a selection are made and drawn (of
  every visible protein when nothing is selected); its checks show the first of them.
  - **Kind:** **Molecular (SES)**, the surface a water-sized probe rolling over the atoms
    cannot pass (Connolly's; new surfaces are this kind); **Solvent Accessible (SAS)**, where
    the probe's center can go; **Van der Waals**, the atoms' spheres; **Gaussian**, a smooth
    skin over the atoms (quicker, and the kind of every surface in documents saved before
    there were kinds, which open as they were).
  - **Style:** **Solid**, **Mesh** (its wire frame) or **Dots**.
  - **Color:** by any of the schemes above (by hydrophobicity at first), **By
    Electrostatics** (below), or **Single Color…** for one color over the whole surface.
    **Transparency:** opaque to 80%.
  - **By Electrostatics** colors the surface by the Coulomb potential 1.4 Å above it (as
    ChimeraX does), red at −10 through white to blue at +10 kcal/(mol·e), with a dielectric of
    4r and the charges within 15 Å. Amino acids and nucleotides take Amber's charges (Cornell
    et al. 1995, as kept by ff14SB, from PDB2PQR's table) for the protonation states their
    hydrogens show: a prepared His with HD1 and HE2 is charged, Asp with HD2 or Glu with HE2
    neutral, Lys missing a hydrogen neutral. An unprepared protein takes the usual states
    (His neutral, on NE2), its hydrogens' charges on their heavy atoms. Chain ends are charged.
    Ligands, glycans and residues the table lacks take Gasteiger charges: a prepared ligand in
    the protomer [Prepare](#prepare) chose, else in the form of its template, which may be
    neutral (ATP's phosphates are protonated there). Ions take their charge, or their usual
    one when the file gives none (Na⁺, Mg²⁺, Zn²⁺, Cl⁻ …); waters none. The surface is white until the potential arrives (a
    few seconds for a large protein, in the background); carving and Dots keep it.
  - **Carve Around Selection ▸** 3, 5, 8 or 12 Å shows only the part of every protein's surface
    near the selected atoms, whichever protein they belong to (a ligand or a pose); **Show
    Whole Surface** undoes it. Carving cuts the surface already built, so it is immediate.
  - **Probe Radius** (1–2 Å; SES and SAS), **Quality** (Coarse, Normal or Fine: how fine
    the grid it is made on is) and **Smoothing** (Off, Low or High; not for Gaussian)
    rebuild the surface. A molecular surface of a large protein takes a few seconds, in the
    background.
- **Interactions** draws dashed lines from the ligands (or the selection) to what they touch:
  hydrogen bonds (blue), salt bridges (magenta), π-stacking (green), cation–π (orange),
  hydrophobic contacts (grey), halogen bonds (cyan), metal contacts (light gold) and water bridges
  (light blue). Point at an atom to read its interactions in the status line.
  **Interaction Types ▸** ticks which of these kinds are drawn (all, at first), on every
  structure; the choice is saved with the session and kept through Prepare and other edits.
  Ligands to Table, Interactions to Table and the interaction fingerprint still count every kind.
- **Water bridges** need a water with a hydrogen bond to both the query and a partner that is not
  a water or an additive (2.5–3.5 Å each). Each bond must point the right way, judged as for hydrogen bonds below: a
  donor's hydrogen at a D–H···O angle of at least 100° where it has one, otherwise the donor's
  or acceptor's heavy atoms. The angle at the water between its two partners (the donor's
  hydrogen where there is one) must be 71–140°, PLIP's limits (Salentin et al. 2015, after
  Jiang et al. 2005). A water in reach of both but in line with them, or off an acceptor's
  lone pair, bridges nothing.
- **Hydrogen bonds use the hydrogen where the structure has one** — a protein you have
  prepared, an NMR model, a ligand drawn with hydrogens: the H···A distance (up to 2.5 Å, after
  Baker & Hubbard 1984) and the D–H···A angle (at least 100°, as PLIP) decide, and the dashed
  line starts at the hydrogen, which is where the bond is. A hydroxyl, thiol or ammonium
  (Ser, Thr, Tyr, Cys, Lys, a ligand's OH) is the exception: its hydrogen turns freely about
  one bond and was placed rather than seen, so when it points elsewhere the bond is judged on
  the heavy atoms as below, and drawn from them. A residue that carries hydrogens and has none
  on a polar atom does not donate there, so a prepared His donates from the nitrogen its
  tautomer protonated and not from the other. Structures from the PDB rarely carry hydrogens; without them the donor and
  acceptor heavy atoms are judged on their distance (2.5–3.5 Å) and on whether a hydrogen or
  lone pair could point the right way, and the line runs heavy atom to heavy atom. Either way
  the distance reported is the donor–acceptor one.
- **Focus ▸** a ligand shows its pocket (residues within 5 Å as sticks, the waters and ions
  touching it) with its interactions and a see-through skin over the pocket coloured by
  hydrophobicity, fades the rest of the protein and hides other small molecules. A slice of
  depth around the pocket follows the camera, with the back of the pocket fading into the
  background; faded parts of the protein between you and the ligand are cut away, so a ribbon
  never crosses in front of it. **Pocket Surface ▸** in the same menu turns the skin off or on (**Show**) and
  sets its **Style** and **Transparency**; the skin is of its protein's surface kind. **Clear
  Focus** undoes the fading.
- **Measure ▸ Distance / Angle / Dihedral**, then click two, three or four atoms. Esc stops
  measuring.
- **Editing** (the Edit menu, or **Edit Selection ▸** in the right-click menu) changes the
  protein itself; each edit is one step that Undo takes back, here or in the table, and is
  saved with the document. Everything else about the atoms that stay (how they are shown,
  coloured, selected, labelled or measured) stays as it was.
  - **Delete** (Delete or Backspace, here and in the Sequence window) removes the selected
    atoms; an atom goes with its alternate locations. To delete some waters and keep the rest,
    pick them out in the Waters row of the Sequence window (Ctrl/⌘-click adds one at a time)
    and press Delete. Deleting every atom of a protein is refused: use **File ▸ Remove
    Protein**. When table rows name a residue that lost atoms, the status line says how many.
  - **Element ▸** (C, N, O, S, P, F, Cl, Br, I, H or **Other…**) and **Charge ▸** (−2 to +2)
    change the selected atoms; their bonds stay.
  - **Rename Residue…** renames the selected residues (1 to 5 letters or digits); a residue
    outside the protein's chains takes the kind its new name gives it, so a water renamed `NA`
    becomes an ion. **Rename Chain…** renames the chain holding the selection (no other chain
    may have the name). **Renumber Residues…** numbers the selected residues of one chain from
    the number you give, in order; a number another residue of the chain has is refused, and
    the dialog says which. When a ligand that table rows name (**Ligand ID**) is renamed or
    renumbered, those rows follow in the same step.
- **Models.** A file with several models (an NMR ensemble) gets a bar under the view: step
  with the arrows or `[` and `]`, type a model's number, drag the slider, or press Play to loop
  through them. The selection, colours and measurements stay as the atoms move. When several
  shown proteins have models, the list at its left picks which one the bar steps. The model
  shown is saved with the document.
- Right-click an atom to select its residue, chain or surroundings, focus a ligand, switch an
  alternate location (where a file models two residues at one place, such as PRO and SER at
  22 in 1EJG, the letter picks which residue is shown, one at a time), copy a ligand's SMILES, **Add Ligand to Table** (a cofactor, say, or a
  ligand whose row you deleted), **Label ▸** its residue or the atom, start a **Measure From
  Here ▸** distance, angle or dihedral with it as the first atom, **Delete Atom** or
  **Delete Residue** (what was clicked), **Edit Selection ▸** (as in the Edit menu), or **Copy
  Selection As ▸ PDB / mmCIF**.

Surfaces, interactions and the cartoons of very large proteins are built in the background,
so the viewer never stops answering: each appears as soon as it is ready, and until then the
previous one stays on screen.

Everything shown, selected, focused and measured is saved with the document. Fetch from PDB
only sends the ID, over HTTPS; when IT turns downloads off it says so instead.

### Prepare

**Tools ▸ Prepare…** gives a structure the hydrogens, side chains and protonation states it
needs before anything is measured or docked. Crystal structures carry almost no hydrogens, and
X-ray maps cannot tell an Asn's oxygen from its nitrogen, so these are worked out here.

- **Protein** is the one to prepare (the one the viewer is looking at, by default).
- **pH** decides which form each acid, base and ligand takes.
- **Waters:** keep them all, only those near a ligand (within 5 Å, or **More Options ▸ Water
  Distance**), or none. Every water kept gets
  its two hydrogens: with Optimise on, those within 5 Å of a ligand are oriented together with
  the other choices; the rest are turned one at a time towards the N and O around them.
- **Repair** builds side-chain atoms and terminal oxygens the file is missing, from the
  component templates, turning each rebuilt side chain (in 30° steps about its first two
  bonds) to where it crowds its neighbours least. Where the file stops short of the deposited
  sequence, the ends are capped rather than left as charged termini: an **ACE** cap before the
  residue that resumes, an **NME** cap after the one that stops, as a force field needs. This
  holds at a gap inside a chain and at a chain end the sequence goes on beyond (an OXT is built
  only on a chain's true last residue). A residue missing on its own is **built** instead (its
  backbone from its neighbours, its side chain from its template), since two caps would meet
  where it is. Longer loops are not built. A residue deposited without its CA or C (or, where
  a stretch stops, its backbone O) has nothing to build a side chain on or to cap from, so it
  is **removed** (⚠ in the report, under Residues removed) and the residue before or after it
  takes the cap instead: 4AGC ends in a GLN 1169 that is only its N. A file without a deposited sequence (no SEQRES or
  entity sequence) keeps its chain ends as charged termini.
- Ions take their formal charge from the file, else from the component (``FE`` is Fe³⁺, ``FE2``
  Fe²⁺, ``CU1`` Cu⁺) or the element (Na⁺, Mg²⁺, Zn²⁺, Cl⁻). An ion with neither is marked ⚠:
  set its charge by hand (Edit ▸ Charge).
- **Optimise** chooses each Asn/Gln flip, each His tautomer (HID, HIE or HIP) and each
  hydroxyl, thiol and amine orientation by the hydrogen bonds they make, the choices solved
  together rather than one at a time.
- **Ligand forms** (on by default, More Options ▸ Ligand Forms) let each ligand take another
  protonation state or tautomer when its pocket favors it. The forms and their free energies
  in water come from **chemsuite pKa**, chemsuite's own model (the Uni-pKa method, trained on
  pKaHub's aqueous pKa data), which scores every protonation state of the whole molecule, so
  sites that affect each other are handled together; it comes with chemsuite. The report's
  **Ligands** line names the model used. A ligand with more than 48 protonation states (ATP, NADP) also stays as
  drawn, and the report says so on its line.
- **More Options ▸ Titrate** works out each group's pKa in its surroundings and takes the form
  it mostly holds at that pH, for the groups near a ligand (Pocket) or for all of them. It is
  off by default: on the benchmark the textbook forms at pH 7 are right more often.

The rest of **More Options** fine-tunes the run. Each one defaults to what Prepare does
without it, and each has its explanation in its tooltip:

| Setting | Default | What it changes |
|---|---|---|
| **Side Chains** | on | Builds missing side-chain atoms. Off: they stay as the file has them (a residue Prepare builds from nothing still gets its side chain). |
| **Lone Residues** | on | Builds a residue missing on its own. Off: the ends either side are capped instead. |
| **Cap Ends** | on | Puts ACE/NME caps on the cut ends. Off: the ends are left as the file has them. |
| **Alternate Locations** | Most Occupied | Which location an atom modelled more than once keeps: the most occupied, or the first letter. |
| **Water Distance** | 5 Å | How near a ligand a water must be, with Waters ▸ Near Ligands. |
| **Water B-factor Limit** | Off | Drops waters whose B-factor is above it. A water kept by hand in Inspect stays. |
| **Orient Waters** | on | Turns the waters near ligands with the rest of the hydrogen-bond network. |
| **Ligand Forms** | on | Lets a ligand take another protomer or tautomer when its pocket favours it. Off: as drawn. |
| **His** | Auto | Gives every His one tautomer (HID, HIE or HIP) instead of choosing each. A His bound to a metal keeps its own. |
| **Flips** | on | Lets Asn, Gln and His flip. Off: they keep the file's orientation. |
| **⚠ Threshold** | 1 kcal/mol | A choice whose next best is closer than this is marked ⚠, in the report and in Inspect alike. |

The repair settings apply only with **Repair** on. His, Flips, the threshold, Titrate,
Orient Waters and Ligand Forms apply only with **Optimise** on. When a setting is changed from
its default, the report's first line (**Settings**) names it, so the report says how the
structure was made.

The work runs in the background and can be cancelled. What comes back replaces the protein in
**one** step: Undo (Ctrl+Z, in the table or the viewer) puts the structure back as it was.
Everything else about the protein — how it is shown, coloured, selected and labelled — stays.

The **report** opens beside the structure in the viewer's own window — with the viewer closed
it goes beside the table instead — one line per thing found or decided: alternate
locations, chain breaks, gaps and chain ends not modelled, residues built, caps, missing atoms,
covalent links, disulfides and metal sites, each flip, His state and uncommon protonation
state. A line marked ⚠ is worth your eye: a chain break or gap, a residue built from nothing,
a residue with no template, an ion of unknown charge, or a choice whose next best was nearly as
good. Click a
line to select its residue and zoom the viewer onto it, with some of its surroundings; a line
naming two residues (a gap, a chain break, a disulfide) zooms out to show both and selects the
first. A hidden protein's residue is selected but not zoomed to. Hydroxyl and amine hydrogens with nothing to bind to
are free to spin, so they are counted in one line rather than listed. The panel's title gives
the pH and the structure's net charge; right-click (or Ctrl+C) copies the report as text. The
report is saved with the document and comes back with it.

Accuracy, against 262 neutron structures: Asn flips 86%, Gln flips 80%, His orientation 90%,
His tautomer 64% (Reduce: 85/79/88/62).

### Inspect

**Tools ▸ Inspect…** (in the viewer's **Tools** menu) runs the same preparation as Prepare…,
with the same settings, but applies nothing: it lays out everything it found and would decide
in the **Inspect** panel beside the structure, for you to go through first.

- **Issues** are what the file is missing or carries twice, grouped by kind: alternate
  locations, missing atoms, chain breaks, gaps in the sequence, chain ends not modelled,
  residues built, caps, residues with no template, covalent links, disulfides and metal
  sites. A missing side chain has a box: untick it and it
  is left as the file has it, a gap rather than a modelled guess. The others are there to be
  seen.
- **Choices** are every decision, grouped by kind: alternate locations, amide flips, His
  states, protonation states, rotors, ligand forms, crystal waters and water orientations.
  Each line gives the state chosen and how much worse the next best would be. Double-click the
  state to choose another: each state is listed with its cost against the best, so you can see
  how close the call was. A choice made by hand is in bold. Long groups (rotors, waters) come
  folded. A His can always be made HIP; the other protonation states (ASH, GLH, LYN, CYM, TYM,
  a neutral N-terminus) are choices only for the groups that titrate, so set **Titration** to
  **Pocket** or **All** to choose them.
- **Alternate locations** has a line for each residue the file models more than once, with
  the location kept (the highest occupancy, as Prepare… keeps); choose another letter to keep
  that one instead. Each letter is listed with its mean occupancy. Two residues modelled at
  one place make one line, under the first, each letter naming the residue it keeps, so the
  prepared protein has one residue there. With **Repair** off the
  file's locations are left as they are, so there are no such lines.
- **Crystal waters** has a line for each water in the file, **Kept** or **Dropped** as the
  **Waters** setting left it; overrule any one either way.
- **⚠ Only** shows just the lines worth a look, by the same rule as the preparation report:
  near ties, mixed protonation, chain breaks, residues with no template. A water's orientation
  and a hydroxyl with nothing to bind to are never flagged (any answer is as good).
- Click a line to select its residue and zoom onto it, as in the report.
- **Re-survey** solves again with your choices held. Choices are not independent — the
  hydrogen-bond network couples neighbouring sites — so holding one flip can change what the
  residues around it should do, and this shows you before anything is applied.
- **Apply** prepares the structure as reviewed, in **one** undoable step, and the review
  becomes the protein's preparation report, led by what you decided by hand ("Left as in the
  file, by hand", "Chosen by hand: … (proposed: …)") so the document records how its structure
  was made. **Reset** puts every line back as the survey proposed.

What you have decided so far is saved with the document until you apply it. Reopening the
document surveys the structure again and puts your decisions back, so the costs shown are
always worked out from the structure as it is. Decisions name their residues by chain, name
and number, so an edit elsewhere in the structure in the meantime does not move them.

### Minimise

A prepared structure can be relaxed in [OpenMM](https://openmm.org), so that hydrogens, rebuilt
side chains and flipped groups settle against their neighbours while the deposited atoms stay
near where the file put them. OpenMM must be set up first (see
[Setting up OpenMM](#setting-up-openmm)); until it is, Minimise is greyed out and says why.

- **Tools ▸ Prepare…** with **Minimise** ticked leaves the protein as it was and adds two
  new ones in **one** step: **X-prep**, the prepared structure, and **X-min**, the same
  minimised. Minimise needs **Repair** with **Side Chains** and **Cap Ends** on, since a force
  field has no template for half a residue. If the minimisation cannot run, X-prep is still
  added, and its report's first line (⚠) says why.
- **Tools ▸ Minimise…** (in the viewer's **Tools** menu) adds **X-min** for a protein
  already prepared. A protein edited after it was prepared must be prepared again first. A
  minimised protein is never minimised or prepared again: its restraints would hold it to the
  first result instead of the deposit, so start again from the one it came from. Each new
  minimisation of the same protein is another copy (**X-min-2**…). With several models (NMR),
  the model shown is minimised.

Each new protein is shown in place of the one it was made from, which is hidden; Undo removes
the new ones and shows the original again. The minimisation runs in three stages: hydrogens
only; then side chains and waters, with the backbone, ligands and ions held; then the whole
complex, each heavy atom held near its deposited position. Atoms Prepare built are never held.
A transition metal (Zn, Cu, Ni, Co, Fe, Mn, Cd…) keeps its site as Prepare left it: each N, O
or S bound to it stays at its distance, and its other neighbors within 6 Å come no closer.
While the deposited heavy atoms end further than the **RMSD Limit**, the last stage runs again
held twice as hard (up to three times).

| Setting | Default | What it changes |
|---|---|---|
| **Force Field** | ff14SB | ff14SB with TIP3P water, or ff19SB with OPC water. Ligands take OpenFF Sage either way. ff19SB is not offered while crystal waters are kept (OPC's fourth site is not placed yet). |
| **Solvent** | GBn2 | GBn2 implicit solvent, or none (vacuum). Kept crystal waters stay as explicit waters in both. |
| **Restraints** | Normal | How hard the deposited heavy atoms are held: Loose (half), Normal or Firm (twice). |
| **RMSD Limit** | 0.30 Å | The heavy-atom RMSD the run must end within; above it, the report marks ⚠. |
| **Platform** | Auto | Where OpenMM runs: Auto takes the fastest that works (CUDA, OpenCL, CPU). The others are listed once OpenMM has run this session. |

In Prepare these settings are under **Minimisation**; Prepare and Minimise share them, and
they are remembered until chemsuite closes.

The **report** of a minimised protein gives the energy before and after, how far the heavy
atoms moved (RMSD overall, in the pocket — the deposited atoms within 5 Å of a small molecule
— and of the ligand), the closest contact before and after, where the largest force is left,
the five residues that moved most (click one to zoom to it), the force fields and ligand
charges used, the platform, the time and the OpenMM version. A run that ends beyond the RMSD
Limit, or with two atoms closer than 2.2 Å, is marked ⚠. The report, and how the protein was
prepared and minimised, are saved with the document.

### Dock

**Tools ▸ Dock…** docks the document's ligands into a pocket with
[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) 1.2.7, which ships with
chemsuite (Apache-2.0; no download and nothing to install). The poses are kept with the
protein and listed in the **Poses** panel beside the pocket they were docked into, since that
is where the work was done; **Copy to Table** puts the ones you choose into the table, where
each becomes a row linked to its protein, draws in its pocket when its row is current, and
answers plots, filters and SAR as anything else does.

Prepare the protein first: Vina types an atom by the hydrogens on it, so a structure with none
is turned away with a line pointing at **Tools ▸ Prepare…**.

- **Protein** is the receptor; **Ligands** the structure column to dock, over the rows you
  choose (All, Shown or Selected). A row whose structure has no 3D coordinates gets one
  conformer first; a row with no structure at all is skipped and counted.
- **Box ▸ Around** is where the search box comes from: the **Focused Ligand** in the viewer,
  the **Selection** there, the **Whole Protein**, or **Typed** numbers of your own.
  **Padding** is the room left around those atoms. The centre and size are shown either way,
  in ångströms.
- **Scoring** is Vina's own function or **Vinardo**; **Exhaustiveness** how hard the search
  tries (doubling it roughly doubles the time); **Poses** how many to keep per ligand, best
  first. **More Options** holds the **Energy range** (how far above the best pose to keep) and
  the **Seed** (the same seed docks the same way every time).

The run is a background job with progress, and cancelling keeps the ligands that already
finished.

**Save Input Files…** writes the run to a folder instead, to dock with Vina at the command
line. **Command Line ▸ Write** sets how much: **Everything** (`vina.conf` with the search box
and the conditions, `receptor.pdbqt`, one file per ligand in `ligands`, `ligands.csv` saying
which row and conformer each file came from, and `run.sh` / `run.bat`, which run
`vina --config vina.conf --batch ligands --dir out`), **Receptor** (the config and the
receptor) or **Config** (the config alone). Config needs only a box, so it can be saved before
the protein is prepared, or with no protein at all when the box is typed. A row that already
carries 3D conformers keeps them, one file each; a row with none gets **Conformers** of them
made (1 by default). Vina searches rotatable bonds itself, so more than one mainly helps with
the geometry it holds fixed, such as a saturated ring's pucker.

Every run is kept, numbered in the order it was made and labelled with how it was scored
("Run 2 · Rescored (vinardo)"), so a new run never overwrites an earlier one. The **Poses**
panel shows one protein's runs, one at a time: the list at its top chooses which, newest
first, and **Remove Run** forgets the one shown (rows already copied to the table stay). It
groups that run's poses by ligand, each ligand showing its best affinity and its poses under
it: **Pose** (1 is Vina's best), **Affinity (kcal/mol)**, **RMSD l.b.** and
**RMSD u.b.** (Vina's own spread of that pose against the best one). Click a pose to see it in
the pocket. Ctrl+C copies the scores as text, and **Copy to Table** writes the chosen poses as
rows holding **Protein**, **Ligand** (the pose), **Pose**, **Affinity (kcal/mol)**,
**RMSD l.b.**, **RMSD u.b.** and **Docked From** (the row it was docked from) — then run
**Tools ▸ Interactions to Table** on those rows to profile them. The runs are saved with the
document and come back with it, panel and all.

Only the parts of a receptor that AutoDock knows are used: the polymer, waters, ions and
glycans. Ligands, buffers and additives are left out — a ligand is what is being docked — and
the run says which residue names were left out.

**Tools ▸ Rescore Poses…** scores poses that already exist exactly where they stand, with
Vina or Vinardo. **Poses** chooses what to score: **Docked**, the poses of one of the
protein's runs (**Run**, the one shown in the Poses panel to begin with); or a **Table
Column** of structures over the rows you choose. The scores are a new run beside the one they
came from, which is left as it was, so the two functions can be compared pose by pose. Use it
too to score a pose moved by hand in the sketcher.

Docking is not yet covalent, does not move side chains, and docks one receptor at a time.

### Setting up OpenMM

Minimising a prepared complex (and, later, molecular dynamics and binding free energies) runs in
[OpenMM](https://openmm.org) with the Amber ff14SB protein force field and
[OpenFF Sage 2.3](https://openforcefield.org) for ligands, whose partial charges come from
OpenFF's NAGL model. These live in a Python environment of their own, kept warm in the
background once started. **Minimise** (see [Minimise](#minimise)) is the first thing that
uses it.

1. Install [Miniforge](https://conda-forge.org/download/) (or any conda).
2. Create the environment from the file for your system in chemsuite's `envs` folder, for
   example `conda env create --solver libmamba -f envs/openmm-win-64.yml`, which makes
   `chemsuite-openmm` in Miniforge's `envs` folder.
   On Windows and macOS, the file's first lines give one more command, which installs
   `openmmforcefields`. The environment takes about 3.5 GB.
3. In **Settings ▸ External**, set **OpenMM Python** to that environment's `python`
   and press **Test**. It reports the OpenMM, openmmforcefields and OpenFF versions and each
   platform a calculation can run on: **CUDA** on an NVIDIA graphics card, **OpenCL**, **CPU**.

On an NVIDIA card, the environment's `cuda-version` must not be newer than the driver supports
(`nvidia-smi` shows it); a mismatch leaves CUDA out of the Test result. On Windows, an
environment inside WSL works too: set **OpenMM Python (WSL)** to its Linux path (such as
`/home/you/miniforge3/envs/openmm/bin/python`), created from `envs/openmm-linux-64.yml`. It is
used when no **OpenMM Python** is set.

### Retrosynthesis

Chemistry ▸ Reaction ▸ Retrosynthesis plans synthetic routes to each structure with
[AiZynthFinder](https://github.com/MolecularAI/aizynthfinder) (AstraZeneca, MIT). A
neural-network policy proposes reactions, a tree search takes the molecule apart, and a route
counts as **solved** when every precursor it ends on is in a stock file you supply.

Everything runs on your own machine: the policy models, the reaction templates and the stock
file all come from the `config.yml` you point chemsuite at. A config that would send
predictions to a model server (`use_remote_models: true`) or query a stock database on another
host is refused, and nothing about your structures is sent anywhere.

**Setting it up.** AiZynthFinder needs an older RDKit than chemsuite, so it gets its own Python
environment:

1. Create an environment with Python 3.10–3.12, for example
   `conda create -n aizynth python=3.12`.
2. In it, run `python -m pip install aizynthfinder "pyarrow>=19,<26"`.
3. In **Settings ▸ External**, set **AiZynthFinder Python** to that environment's
   `python` and press **Test**: it reports the AiZynthFinder version it found.
4. Press **Download…** on the **AiZynthFinder config** row and choose a folder. chemsuite
   fetches AiZynthFinder's public USPTO policies, ring-breaker policy, filter and the ZINC stock
   (about 1.5 GB, from Zenodo and figshare) in the background (follow it in the Log window's Jobs tab), writes
   a `config.yml` beside them, and sets the config path when it is done. Files already in the
   folder are not fetched again, and a cancelled or interrupted download resumes where it
   stopped. A `config.yml` of your own in that folder is left alone; the new one is then
   `chemsuite_config.yml`. If you already have the data, choose your own `config.yml` with
   **Browse…** instead. Press **Test** on the config row to check it.

Downloading is the only time this feature uses the network, and only to fetch these public
files. When IT turns network use off (policy `network_allowed: false`), Download… is greyed
out and the data has to be put in place by hand.

Until both are set, Retrosynthesis says what is missing and will not run.

- **Input.** A structure column and the rows. The search starts on the selected rows when you
  have a selection, because it costs a minute or more per structure: run a handful of rows, not
  a whole library. The same structure is searched once however many rows hold it.
- **Policies** and **Stocks** list what your `config.yml` offers. Every one is used unless you
  clear some; `ringbreaker` is trained on ring-forming reactions. A route is solved when its
  precursors are in one of the ticked stocks.
- **Time limit**, **Iterations** and **Routes to keep** bound one structure's search: it stops
  at whichever limit comes first, and keeps that many of the best routes.
- **More Options** holds the search **Algorithm** (Monte Carlo by default; Retro\*,
  breadth-first and proof-number search explore the same templates in a different order),
  **Exploration (C)**, the filter policies, the **Longest route** allowed, **Stop at the first
  solved route**, whether a buyable target is still taken apart, and route clustering.
- **Skip rows already searched** leaves rows that have AiZynth values alone. Clear those cells,
  or turn it off, to search them again.
- **Columns.** Each searched row gets the search's statistics:

  | Column | Meaning |
  |---|---|
  | AiZynth Solved | Yes when the best route ends entirely in the stock. |
  | AiZynth Score | The state score of the best route (0–1; 1 is a fully solved, short route). |
  | AiZynth Steps | Reactions in the best route. |
  | AiZynth Routes, AiZynth Solved Routes | How many routes were kept, and how many of those are solved. |
  | AiZynth Time (s), AiZynth First Solution (s) | How long the search took, and when its first solved route appeared. |
  | AiZynth Precursors in Stock / not in Stock | The precursor SMILES of the best route, by availability. |
  | AiZynth Nodes, AiZynth Precursors, AiZynth Policies, AiZynth Stock | The size of the search, and which policies and stocks it used. |

  Anything else your config reports (an extra scorer, a post-processing job) gets its own
  `AiZynth …` column too.
- **Cancelling** keeps the structures already searched.
- **The routes** go to the **Routes** pane (below), which opens with the first result.

#### The Routes pane

Chemistry ▸ Routes shows the routes found in this document. It saves with the document, so a
reopened file still has its routes.

- **Target** picks which searched structure to look at, with its number of routes.
- The **route table** lists each route: solved or not, its score, how many steps it has, and
  how many of its ending molecules are in stock. Extra scores from your config appear as
  further columns. Pick a row to draw that route.
- **The tree** puts the target on the right with its precursors branching to the left, the way
  a retrosynthetic scheme is read. A molecule's card is framed **green** when it is in stock
  and **amber** when it still needs a route of its own; its tooltip names the stock it was
  found in.
- **Step 1** is the earliest reaction and the highest number makes the target. Click a step to
  see its reaction SMILES, the policy that proposed it, that policy's probability, the reaction
  template and its classification underneath; click empty space to close it again. Right-click
  a step to copy the reaction SMILES or draw it in the sketcher.
- **Right-click a molecule** to copy it (SMILES, molfile, InChI, InChIKey or an image) or to
  add it to the table as a new row — handy for turning in-stock precursors into a
  building-block list.
- **Drag** or **scroll** to pan. **⌘-scroll** (Ctrl-scroll on Windows and Linux) or a
  trackpad **pinch** zooms around the pointer, from 20% to 400%. Structures are redrawn for the
  zoom they are shown at, so they stay sharp close up. A route is fitted to the pane when you
  open it and whenever the pane is resized, until you pan or zoom; **double-click empty space**
  to fit it again.

---

## 11. Table and analysis tools

### Statistics

Data ▸ Table ▸ Statistics… opens a panel with five tabs. It works on the shown rows, all rows, or the
selected rows, and recalculates as the data, filters or selection change.

| Tab | What it gives |
|---|---|
| Summary | N, missing, mean, SD, min, 5%, Q1, median, Q3, 95%, max, skew, kurtosis for every number column. |
| Correlation | Pearson r or Spearman ρ matrix. On big tables it runs as a background job when you ask. |
| Outliers | IQR (Tukey fences), Z-score or modified Z (median and MAD), with a threshold; label rows by a column. |
| Curve Fit | Polynomial (choose the degree), exponential, logarithmic or power fit of Y on X. |
| Tests | Welch t-test, paired t-test, Mann-Whitney U, one-sample t-test, Shapiro-Wilk normality. |

### MPO Score

Data ▸ Scoring ▸ MPO Score… combines several properties into one 0–1 score for **All**,
**Shown** or **Selected** rows. Add one row per property (**Add Property**), each with a weight
and one of these desirability shapes:

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

- **Missing values: Score on the other properties.** Rows missing a property are scored on the
  rest (otherwise they get no score).
- **Desirability: A column per property.** Also adds each property's desirability.
- **Load…** and **Save…** keep a profile in a file to reuse.

### Selectivity

Data ▸ Scoring ▸ Selectivity… takes a primary assay and its off-target assays, and works out
log selectivity (1 = 10-fold).

- Tick the **Off-targets**, and say whether the values are **p-Values** (pIC50, pKi…; higher
  is more potent) or **Concentrations** (IC50, Ki… in one unit; lower is more potent).
- Choose **All**, **Shown** or **Selected** rows.
- It adds a selectivity column per off-target and a **Min logSel** column.
- **Narrow Window.** **Flag** adds a **Narrow window** column marking potent compounds (at
  least **Potent at**, or at most for concentrations) whose selectivity is below **Window
  below**.

### Bin Values

Data ▸ Bin Values… groups a number column into bins.

- **Methods:** equal width, equal count, or custom edges.
- **Options:** a fixed width, from/to limits, and open-ended bins below the first edge and above
  the last.
- **Output:** a **range column** ("[5, 6)") and a **class column** with names you choose. Low /
  Medium / High is filled in for you.
- **Rows:** all, shown or selected rows; other rows stay blank.
- A live preview shows each bin's count.

---

### Calculator

Tools ▸ Calculator… (Ctrl+Shift+=) adds a column worked out from a formula. The
formula is checked as you type, and **Preview** shows the result for the first shown rows.
**Insert ▸ Column** and **Insert ▸ Function** put a column name or function at the cursor.

- **Columns:** `[IC50 (nM)]` for any name, or a plain name such as `MW`. Case does not matter
  when the name is unambiguous.
- **Operators:** `+ − * /`, `^` or `**` for powers, `%` for the remainder, comparisons
  `= != <> < <= > >=`, and `and`, `or`, `not`. `+` joins text when either side is text.
- **Conditions:** `if(condition, then, else)`, e.g. `if([pIC50] >= 7, "Active", "Inactive")`.
- **Blanks:** a blank cell gives a blank result, and so do dividing by zero and the log of a
  negative number. Use `ifblank(x, 0)` to fill blanks and `isblank(x)` to test for them.

| Group | Functions |
|---|---|
| Math | `abs sqrt exp ln log log10 log2 floor ceil round(x, digits) sin cos tan`, `min(a, b, …)`, `max(a, b, …)`, `clip(x, low, high)`; constants `pi`, `e` |
| Logic | `if isblank ifblank` |
| Statistics (over the rows calculated) | `mean median std sum count percentile(x, p) rank(x, descending) zscore normalize`; `min(x)` and `max(x)` with one argument |
| Text | `len upper lower trim contains(t, part) concat(a, b, …) number(t) text(x)` |
| Chemistry | `MW ExactMW logP TPSA HBD HBA RotB Rings AromaticRings HeavyAtoms Fsp3 Charge QED SAScore Formula SMILES InChIKey Scaffold`, each taking a structure column, and `substructure(structure, "SMARTS")`, which counts the matches |

Examples: `9 - log10([IC50 (nM)])` gives pIC50 from an IC50 in nM;
`[pIC50] - logP([Structure])` gives lipophilic efficiency (LipE);
`1.37 * [pIC50] / HeavyAtoms([Structure])` gives ligand efficiency in kcal/mol per heavy atom;
and `substructure([Structure], "[NX3;H2]c")` counts primary anilines.

Choose the rows (All, Shown or Selected; other rows stay blank), the new column's name, and
whether to replace a column of that name. Comparisons give True / False text columns. The
column is one undoable step. Large tables are calculated in the background; cancelling keeps
the rows already done, except for formulas with statistics, which need every row at once.
Formulas cannot run code, open files or reach outside the table.

### Standardize Values

**Tools ▸ Standardize Values…** puts values given with a unit, such as assay potencies in one
column and nM, µM or mg/mL in the next, into one new column in a single unit. (Chemistry ▸
Standardize is a different tool: it cleans up structures.)

- **Input:** **Values** (a number column, or a text column with values such as `>10`),
  **Units** (a text column) and **All**, **Shown** or **Selected** rows. A column named like
  "IC50" or "Unit" is chosen at first.
- **Units read**, in any case: M, mM, µM (also uM, μM, mcM), nM, pM, fM, mol/L and the words
  (micromolar…); g/L, mg/mL, mg/L, µg/mL, µg/L, ng/mL, ng/L and pg/mL.
- **Method:** **Convert to** the unit you choose (nM at first). Molar units convert among
  themselves, and mass units among themselves, by powers of ten. Between molar and mass units
  each row uses its molecular weight, from the structure in **Molecular weight from** (as
  drawn, so a salt's counterion counts): g/L = mol/L × g/mol. That field is only active when
  the rows mix the two kinds.
- **Output:** a new number column, "<values> (<unit>)", in one undo step. A qualifier read
  with a value (`>`, `<`, `>=`, `<=`, `~`) goes in a column of its own, "<values> (<unit>)
  Qualifier", next to the number. Rows with no number, a unit that is not read, or no
  structure to weigh are left blank; the message says how many.
- **Cancelling** keeps the rows already done.

### Random Numbers

Tools ▸ Random ▸ Random Numbers… fills a column with random values:

| Distribution | Settings |
|---|---|
| Uniform | Low, high, decimals |
| Whole numbers | Low and high, both included (a dice roll is 1 to 6) |
| Normal | Mean, standard deviation, decimals; optionally kept within low and high |
| Log-normal | Median and spread factor (about two thirds of values fall between median ÷ factor and median × factor); optionally kept within low and high |
| Random order | 1 to the number of rows, each once: sort by it to shuffle the table |

A **seed** makes the draw repeatable: the same seed and settings give the same numbers again.
Leave it empty for a new draw each time.

### Random Rows

Tools ▸ Random ▸ Random Rows… picks a random sample from all, the shown or the
selected rows, either a number of rows or a percentage. It then **selects them**, **adds them
to the selection**, or saves them as a **new list** (for example a test set). A seed picks the
same rows again.

### Notepad

Tools ▸ Notepad opens a new note in a pane of its own, docked beside the table like a plot.
Open it again for another note: they are titled Notepad, Notepad 2, …, and double-clicking
the current tab renames it. The notes are saved with the document, so typing in one marks the
document as needing a save, and they reopen as you left them.

- **Toolbar:** Undo and Redo; font and size; **Bold** (Ctrl+B), *Italic* (Ctrl+I), Underline
  (Ctrl+U) and Strikethrough; Text Color and Highlight; align left, center, right or justify;
  Bullets and Numbering (pressed again, the paragraphs leave the list); Decrease and Increase
  Indent (a list item moves a level out or in); Clear Formatting.
- The buttons show the format where the cursor is. With text selected they change it;
  otherwise they change what you type next.
- Pasting keeps the text's formatting but leaves out pictures.
- Closing a note's pane (its ✕) deletes the note; when it has text, you are asked first.

## 12. Similarity

Diverse Subset and Cluster compare compounds either by **fingerprints** or by **number
columns** (scaled Euclidean distance); Similarity Search uses fingerprints and the other
methods below. The fingerprint options are:

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

### Similarity Search

Chemistry ▸ Similarity ▸ Similarity Search… scores compounds against one or more references,
in several ways at once.

**Search in:**

- **This Table:** **All**, **Shown** or **Selected** rows of a structure column. The scores
  are written as columns.
- **A File:** a library file (SDF, SMILES, CSV and the other formats File ▸ Open reads), such
  as a vendor catalog, screened without opening it. The hits open as a new document, best
  first.

**References** (**Compare with**):

- **Structures:** a list. Type or paste SMILES (or draw it with the pencil) and press **Add**.
  **Current Row** adds the table's current row, and **Ligand** adds a ligand of a protein in
  the document, in its bound pose. The current row is added when the dialog opens. A listed
  row is compared with itself too, so it scores 1.000 on its own row.
- **Selected Rows:** every selected row, followed as the selection changes. A selected row is
  never compared with itself, so each row's score is to its nearest *other* selected compound.
- **Fuse references** (with several): a row's score is its **Best** score against any of
  them, or the **Mean** over all of them. A **Nearest Reference** column names the one each
  row is most like.

**Method:** tick one or more.

| Method | Score column | Score |
|---|---|---|
| Fingerprint (any of those above), with Tanimoto, Dice, Cosine or **Tversky** | e.g. ECFP4 Tanimoto | 0–1 |
| Pharmacophore (Gobbi 2D) | Pharmacophore 2D Tanimoto | 0–1 |
| Substructure: Common Substructure (MCS) | MCS Tanimoto | Atoms of the largest common substructure over the heavy atoms of the two molecules, 0–1 |
| Substructure: Same Scaffold, Same Generic Scaffold | Scaffold Match, Generic Scaffold Match | 1 when the Bemis–Murcko scaffolds match, else 0 |
| 3D: Shape, Shape + Features, Field + Shape | Shape Tanimoto, Shape + Features, Field Similarity | As in [Superpose](#superpose): 0–1, 0–2 and up to 1 |

- **Tversky** is asymmetric: **α** weighs the bits only the reference has, **β** those only the
  row has (More Options). α 1, β 0 asks how much of the reference a row contains.
- **3D** superposes each row onto every reference as [Superpose](#superpose) does and scores
  the best fit. The references need 3D coordinates (a 3D structure column, or a protein's
  ligand). A row with only 2D gets conformers first (**Conformers for 2D rows**, 30), and a row
  with 3D is fitted as stored.
- **Combine** (with two or more methods) adds a **Consensus** column: the **Mean** of the scores
  (each on a 0–1 scale: Shape + Features halved, a negative field score taken as 0), the **Mean
  Rank** of each row among the rows, or the **Normalized Mean** (each method's scores
  stretched to fill 0–1 first). A row needs every method's score to get one.
- **Prefilter** (More Options, 500 by default): the common substructure and 3D methods take
  seconds a row, so they score only the rows the fingerprint ranks best, and leave the rest
  blank. The ranking uses the first fingerprint method chosen, or ECFP4 Tanimoto. **Every
  row** turns it off. **MCS time limit** (2 s a row) keeps the best substructure found by then.

**Output** (a table):

- **Sort: Most similar first** sorts by Consensus, when there is one, else by the first score.
- **Select: At least** selects the rows scoring at least that much, by the same score, on a
  0–1 scale.
- **Keep aligned poses** (with 3D) writes each row's best pose, fitted onto the reference it
  is most like, in "*column* superposed". As with Superpose, the [Browser](#browser) shows
  each pose with its reference, and **Pose ΔE (kcal/mol)** gives the pose conformer's energy
  above its molecule's lowest when the conformers carry energies.

**Output** (a file):

- **Keep:** the best this many compounds (1,000 by default; **Every hit** for all), ranked by
  the fingerprint while the file is read. Every method then scores just those, and the
  prefilter applies to them as it does to a table.
- **At least:** keeps only compounds whose fingerprint similarity reaches this.
- The new document has the file's own columns, the score columns, Consensus and Nearest
  Reference.

Fingerprints of the table's rows are cached, so a second search is faster. Cancelling keeps
the rows already scored. Cancelling a file screen stops reading the file and opens the hits
found so far.

### Diverse Subset

Chemistry ▸ Similarity ▸ Diverse Subset… picks *N* compounds that cover the chemistry, using MaxMin:
each pick is the compound farthest from all the earlier picks.

- **Pick from** all, shown or selected rows, and set **How many**.
- **Start from: Add to the selected rows** starts from the current selection, so the new picks
  complement compounds you have already chosen.
- The picks get a rank column, and **Select: The picked rows** selects them.
- Cancelling keeps the picks made so far. Each of them is as diverse as a full run would make
  it.

### Cluster

Data ▸ Cluster… groups compounds using one of three methods:

| Method | You set |
|---|---|
| Butina | A similarity cutoff. A higher cutoff gives smaller, tighter clusters. |
| k-means | The number of clusters. |
| Hierarchical | The number of clusters and the linkage. |

- Clusters are numbered by size, so cluster 1 is the largest.
- The result is three columns: the cluster number, **… Size**, and **… Centre**, which marks
  each cluster's representative compound.
- Rows that cannot be clustered (no structure, a missing number) get no cluster.
- Butina works on fingerprints only. Ward linkage needs number columns. Hierarchical clustering
  has a row limit. The dialog says so, above its buttons, when one of these applies.
- Clustering again keeps the earlier clustering beside the new one ("Cluster (2)").

---

## 13. SAR and QSAR

### Fragments

Chemistry ▸ Decompose ▸ Fragments… cuts each compound into pieces a chemist could join
again.

- **Cut at.**
  - **BRICS** (Degen et al., *ChemMedChem* 2008) cuts 16 kinds of bond and labels each open
    end with the environment it came from (`[16*]`, `[1*]`), so the pieces can only be put
    back together where the rules allow. [Recombine Fragments](#recombine-fragments) needs
    those labels.
  - **RECAP** (Lewell et al., *J. Chem. Inf. Comput. Sci.* 1998) cuts 11 retrosynthetic bonds
    and writes every open end as a plain `*`. A compound with nothing to cut is its own
    fragment.
- **Smallest piece** leaves out fragments with fewer heavy atoms than the number given
  (3 by default). Attachment points are not counted. To recombine the fragments, set it to 1:
  BRICS cuts an amide into an acyl, an aryl and a one-atom linker (`[5*]N[5*]`), and without
  the linker no amide can be rebuilt.
- **Plain attachment points** writes every open end as `*` for both methods, so the same
  piece counts once however it was cut. Leave it off if you will recombine the fragments.
- **Output.** Either or both:
  - **Fragment columns**: `BRICS 1`, `BRICS 2`… (or `RECAP 1`…), largest fragment first, and
    a `BRICS Fragments` count column. A compound cut into more than 20 pieces writes only
    the largest 20.
  - **Fragment panel**: one row per distinct fragment, with how many compounds hold it and,
    when an **Activity** column is chosen, the mean and standard deviation of those
    compounds. **Select Compounds** selects the rows holding the fragments you pick, so you
    can see what a common piece does to activity. **Copy Fragments** copies the table.
- **Cancelling** keeps the compounds already fragmented.

### Central Ring

Chemistry ▸ Decompose ▸ Central Ring… splits each compound around its own most central ring
system, so there is no core to enter.

- **How the core is chosen.** A ring's centrality is the mean bond-path distance from its
  atoms to every heavy atom; the lowest mean wins. Ties prefer the larger ring, then the ring
  holding the lowest canonical rank. Rings sharing an atom with the winner — fused, bridged
  or spiro — join it, so an indole stays whole; a ring joined only by a bond, such as
  biphenyl's second ring, becomes a substituent instead.
- **Output.** A **Core** column with the numbered attachment points (`c1cc([*:2])ccc1[*:1]`),
  then **R1**, **R2**… labelled as [R-Group Decomposition](#r-group-decomposition) labels
  them (`*Cl`, `*OC`). Positions are numbered in the core's canonical rank order, so the same
  scaffold numbers them the same way every time.
- **Each compound is decomposed on its own**, so two rows share an R number only when they
  share a core. Check the Core column before reading R1 across the whole table; to compare a
  series on one scaffold, use [R-Group Decomposition](#r-group-decomposition) instead.
- A compound with no ring, or a ring system of fewer than three atoms, gets an empty core and
  no R-groups.
- **Cancelling** keeps the compounds already decomposed.

### R-Group Decomposition

Chemistry ▸ Decompose ▸ R-Group Decomposition… splits a series onto a common core.

- **Core.** Enter it as SMILES or SMARTS, or let the tool suggest one.
  - Label attachment points `[*:1]`, `[*:2]`… to fix which position is R1, R2.
  - Without labels, RDKit numbers every open position it finds.
- **Output.** A column per position (R1, R2…) and a match column.
  - Substituents are labelled by SMILES with `*` at the attachment (`*Cl`, `*OC`).
  - **H** means unsubstituted.
- **Suggest** fills in the core the selected compounds share (with fewer than two selected,
  the first rows shown).
- **Matching: Exhaustive** is slower; use it for symmetric cores.
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

- **Min. occurrences** leaves out compounds carrying a substituent seen fewer times (a
  substituent seen once is fitted exactly). **Higher is better** sets which way the Best
  Unmade tab ranks.
- **Contributions** tab: each substituent's effect, both as Fujita-Ban (against a reference
  group) and as Free-Wilson (centred).
- **Best Unmade** tab: the highest-scoring combinations not yet made. Set whether higher activity
  is better, and the minimum number of occurrences a substituent needs.
- Substituents that always occur together cannot be separated by the data. They are reported,
  not given arbitrary values.

### Hansch Analysis

Data ▸ SAR ▸ Hansch… regresses activity on substituent constants:

    activity = a·π + b·π² + ρ·σ + c·(steric) + d

- **Constants.** Published constant tables are not shipped with chemsuite. Point the tool at
  your own CSV with these columns:

      smiles,pi,sigma_m,sigma_p,f,r,es,mr,upsilon,sigma_star

  `smiles` is `H` or an attachment SMILES such as `*Cl`. A blank means the constant is unknown,
  and rows needing it are left out. The header needs `smiles` and at least one of the constant
  names above (in any case and order); other columns, such as a name or a source, are ignored.
  A file whose header names no constant is refused rather than read as empty.
- **Terms.** Choose the hydrophobic (π, π²), electronic (Hammett σ, which positions are meta)
  and steric terms.
- **Coefficients** are either **Summed** over positions (classical Hansch) or fitted **Per
  Position**.

### Matched Molecular Pairs

Data ▸ SAR ▸ Matched Pairs… finds pairs of compounds that differ by one small change
at one place. H → X changes are included.

- Set the **Largest change**, in heavy atoms, and **Min. pairs**, the number of pairs a
  transformation needs to be listed. **Activity** is optional: "(none)" lists the pairs
  without activity shifts.
- Transformations (A → B) are summarized over all their pairs, with the activity shift:
  **Pairs**, **Mean Δ**, **SD**, and how many pairs go **Up** and **Down**. Δ is the "to"
  compound's activity minus the "from" compound's.
- The panel draws each transformation's **From** and **To** groups, the attachment point as a
  wavy line. They sort, and **Copy Transformations** copies them, as SMILES (`*Cl`, `H`).
- Choosing a transformation (or several) lists its pairs below, both compounds drawn the same
  way round on the part they share, with the changed part highlighted (for H, the atom it sits
  on).
- Above the pairs, a strip shows one dot per pair on the Δ axis (green up, red down), with zero
  and the mean marked: a tight cluster is a consistent change, a wide spread is not. Click a dot
  to choose its pair. The strip is hidden when no activity column was chosen.
- Double-click a pair to go to one of its compounds in the table: its "from" columns go to the
  "from" compound, the others to the "to" compound.
- **Current Compound** lists only the transformations of the table's current row, and follows
  it as you move through the table. It is saved with the document.
- **Select Compounds** selects the compounds of the chosen pairs, or of the chosen
  transformations when no pair is chosen.

### Activity Cliffs (SALI)

Data ▸ SAR ▸ Activity Cliffs… lists similar compounds with very different activity.

    SALI = |Δactivity| / (1 − similarity)

- Only pairs at least **Min. similarity** alike are considered.
- **Identical: Skip identical fingerprints** ignores salts and stereoisomers, which would otherwise divide
  by zero.
- Each compound's highest SALI goes in a **Max SALI** column, and the partner it reaches it
  with in **Cliff Partner**.
- **Show Network** opens a [Chemical Space Network](#chemical-space-network) of the rows the
  table shows, as an activity-cliff map (as in DataWarrior): compounds joined at the **Min. similarity**,
  edges colored by SALI, points colored by activity and sized by Max SALI, singletons left
  out. The network compares ECFP4 (Morgan, radius 2) Tanimoto similarity, so with another
  fingerprint or metric its edges can differ from the pairs in the list.

---

### Auto QSAR

Data ▸ QSAR ▸ Auto QSAR… fits the best model it can for one column of measured values, and
chooses the features, the model and its settings itself. Give it the **Target** — the property
or activity to predict — and it does the rest.

- **Rows.** The compounds the model learns from. Rows with no value in the target column, or
  whose structure cannot be read, are left out; the dialog says so before you start if too few
  are left. **20 is the fewest** a model can be judged on.
- **Target.** A number to predict, or a column of labels to assign. **Predict** shows which it
  detected; change it to treat a numeric code (1, 2, 3) as classes rather than as a quantity.
- **Effort.** How many models to try.

  | | What it does |
  |---|---|
  | Quick | A few settings for each of the four cheapest model kinds. Use it to see whether the data supports a model at all. |
  | Standard | Every model kind, two dozen settings each. The default. |
  | Thorough | Eighty settings each, for the last fraction of accuracy once the data is worth it. |

- **Test on.** How the held-out compounds are chosen, which decides what the reported score
  means.

  | | What question it answers |
  |---|---|
  | Scaffolds | Can it predict a chemical series it has never seen? Whole Murcko scaffolds are held out together. The default, and the honest one. |
  | Clusters | The same question with Butina clusters instead of scaffolds, for sets where the scaffolds are all alike. |
  | Random | Can it predict a compound much like ones it has seen? Flattering, and the right choice only when the new compounds really will be close analogs. |
  | A column | Train on the early rows and test on the later ones, ordered by a column you choose — a date, a batch number. What the model will actually be asked to do. |

More Options holds the pool and the arithmetic:

- **Features from.** Descriptors, fingerprints (Morgan counts) and the table's own numeric
  columns. The target, and any column named after it, are never offered to the model.
- **Confidence** is how often the interval, or the set of classes, should contain the truth.
  It is *checked* on held-out compounds, not assumed.
- **Most features** caps what the model may end up using. A model on twenty features is one a
  chemist can argue with.
- **Seed.** The same seed, rows and settings give the same model. Change it to see how much
  the answer depends on where the folds fell.

#### What the columns mean

The columns a fit writes are **out of fold**: every number comes from a model fitted without
that row. The prediction beside a compound is therefore not what the finished model thinks of
a compound it learned from — it is what it would have said had it never seen it. That is the
only honest thing to put in the table, and it is why these numbers differ from
[Apply Model](#apply-model) run on the same rows.

| Column | What it holds |
|---|---|
| `<target> Predicted` | The prediction. For classes, the chosen label. |
| `<target> Confidence` | Classes only: the score for the label it chose. |
| `<target> Low`, `<target> High` | The ends of the interval, at the confidence asked for. |
| `<target> Possible` | Classes only: every label the model cannot rule out. Blank means it cannot place the compound at all at this confidence — a useful warning. |
| `<target> Residual` | Measured minus predicted. For classes, `<target> Correct` is Yes or No. |
| `<target> Distance` | How far the compound sits from the training set. |
| `<target> In Domain` | No when that distance is past the 95th percentile of the training set's own distances: the model is extrapolating, and its interval should not be trusted. |
| `<target> Fold` | Which held-out fold produced the row's numbers. |

#### The QSAR Models panel

A fit opens it, and Data ▸ QSAR ▸ QSAR Models reopens it. Models belong to the **document**:
they are saved with it and are there when you open it again.

| Tab | What it shows |
|---|---|
| Scores | Every metric, measured out of fold. R², Q², RMSE, MAE and Spearman for a number; ROC-AUC, PR-AUC, MCC, balanced accuracy and calibration error for classes. |
| Leaderboard | What else was tried and what it scored. A narrow win means the choice of model hardly matters and the features are doing the work. |
| Features | The features the model ended up using, and where each came from. |
| History | One line per fit or refit: the date, what was done, how many compounds, and the score. |

**Save Model…** writes a `.csqm` file — arrays and text, never code — to use in another
document or send to a colleague. **A model file carries the compounds it was trained on**,
both their keys and their measured values, because that is what makes it refinable. Do not
pass one outside the project if the data is confidential.

**Export for Python…** writes the fitted scikit-learn object as a `.joblib` pickle for a
notebook, beside a manifest naming the versions it needs, the feature list, and a short script
showing how to build those features. It is offered for a model fitted in the current session,
since a reopened model is arrays rather than a scikit-learn object. chemsuite **never reads a
pickle back**: loading one runs whatever code it contains, so Apply Model… accepts `.csqm`
only.

### Apply Model

Data ▸ QSAR ▸ Apply Model… predicts with a model this document already holds, or with a
`.csqm` file (**Open…**). It writes the prediction, the interval and the domain columns for the
rows you choose, skips any whose structure cannot be read, and says how many came out outside
the model's domain.

The features are rebuilt from the model's own recipe, so the document needs a structure column
but nothing else prepared in advance.

### Refine Model

Data ▸ QSAR ▸ Refine Model… fits a model again once more data has been measured. The recipe is
kept — the same target, features and split — so the model's history stays comparable; only the
rows and the **Effort** change.

Compounds are matched by InChIKey. A row the model was already trained on **replaces** its old
value, a new compound is **added**, and a compound the model remembers that is no longer in the
table cannot be refitted, so it is left out and the count is recorded in the History tab.

Refining is a genuine refit, not a nudge: the fill values, the chosen features, the settings
and the calibration are all learned again from the larger set, and the new score is measured
the same way the first one was. **Replace the old model** is on by default; turn it off to keep
both and compare their scores before dropping one.


## 14. Chemical space

### Dimensionality Reduction

Data ▸ Dimensionality Reduction ▸ **PCA…**, **UMAP…**, **t-SNE…** or **PaCMAP…** places
compounds on a 2D or 3D map by fingerprints, number columns, or both. Each dialog opens with
the settings it last ran with; **Defaults** sets every setting back.

**Input**

- **Fingerprint** (or **None**) and any ticked **Columns**: a fingerprint alone maps by its
  bits, columns alone by their scaled values, and a fingerprint with ticked columns uses the
  two together, for rows that have all of them (the fingerprints are reduced to their principal
  components, both blocks are brought to the same spread and weighted, and distances are
  Euclidean).
- **Rows**: the rows to map.
- **Fit on** (not t-SNE): **All rows**, **Shown rows**, **Selected rows** or a **List**. The map
  is made from those rows, and the other rows are then placed onto it, so you can fit on a
  reference set and see where new compounds fall in its space. The status bar says how many
  rows it was fitted on. Fit on at least 5 of the rows mapped.

**Method** and **More Options** (open). The fingerprint and column settings are in More
Options for every method:

| Setting | What it does |
|---|---|
| **Fingerprint size** | 512, 1,024, 2,048 (the default) or 4,096 bits, for Morgan, FCFP, RDKit, atom-pair and torsion fingerprints |
| **Radius** | 1 to 4 bonds, for Morgan and FCFP (each starts from its own: 2 for ECFP4 and FCFP4, 3 for ECFP6) |
| **Counts** | How often each feature occurs, not only whether it does. Counts are compared by Euclidean-type distances, so bit distances (Jaccard and the like) are then not available |
| **Fingerprint weight** | With a fingerprint and columns: the fingerprints' share, 0 to 100% (default 50%); the columns get the rest |
| **Column scaling** | **Z-score** (the default: mean 0, standard deviation 1), **Min-max** (0 to 1), **Robust** (median 0, interquartile range 1, so outliers weigh less) or **None** |

| Method | Method section | More Options |
|---|---|---|
| PCA | **Dimensions**, **Components**: PC1 and PC2 (the default), PC1 and PC3, PC2 and PC3, or PC3 and PC4 | **Write components** (Plotted only, or 3 to 20: writes PC1 to PCn as columns), **Loadings** (number columns without a fingerprint: opens a PCA Loadings panel of each column's weight on each component) |
| UMAP | **Dimensions**, **Neighbors**, **Min. distance**, **Guide by** (None or a column), **Groups** (number guide, default 5), **Guide weight** (default 0.5) | **Distance** (Jaccard, the default, Dice, Rogers-Tanimoto, Hamming, Euclidean, Manhattan, Cosine, Correlation), **Spread**, **Epochs**, **Fit on up to** (default 50,000), **Start from** (Spectral, the default, PCA, Random), **densMAP**, **Density weight** (default 2), **Local connectivity**, **Repulsion**, **Negative samples**, **Learning rate**, **Seed** |
| t-SNE | **Dimensions**, **Perplexity** | **Learning rate** (Automatic or a number), **Iterations**, **Exaggeration**, **PCA dimensions** (default Off), **Distance** (Euclidean, Cosine, Manhattan, Chebyshev; for columns, counts or PCA dimensions), **Algorithm** (Exact, or Barnes-Hut when unticked), **Angle** (Barnes-Hut), **Start from** (Random, the default, or PCA), **Seed** |
| PaCMAP | **Dimensions**, **Neighbors** (default 10), **Variant** (PaCMAP or LocalMAP) | **Mid-near pairs** (default 0.5), **Further pairs** (default 2), **Distance** (Jaccard, the default, Hamming, Euclidean, Manhattan, Angular), **PCA dimensions** (default Off), **Local distance limit** (LocalMAP, default 10), **Learning rate**, **Iterations** (the three phases, 100, 100 and 250), **Start from** (PCA, the default, or Random), **Seed** |

**Output**

- **Color by** (the activity column when there is one, or None) and **Size by** (None or a
  number column): what the plot shows.
- **Plot**: untick **Open a plot** to only write the map's columns.

Notes:

- **Jaccard** (Tanimoto), **Dice**, **Rogers-Tanimoto** and **Hamming** compare fingerprint
  bits, so they need a fingerprint alone, with no columns ticked and **Counts** off (and, for
  PaCMAP, **PCA dimensions** Off). When they don't apply the distance switches to Euclidean,
  and back once they do. PaCMAP cannot place rows by a bit distance, so with **Fit on** choose
  another distance.
- **PCA dimensions** is Off by default, so t-SNE and PaCMAP map the features themselves.
  Setting a number first reduces the features to that many principal components, which is
  faster on large sets (bit distances are then not available).
- **Guide by** makes UMAP supervised: compounds with the same value of the column are pulled
  together. A text column's values are the groups; a number column is cut into **Groups** of
  about equal size. Rows without a value are placed by their features alone. **Guide weight**
  runs from 0 (the features alone shape the map) to 1 (the guide column does).
- **densMAP** keeps how dense each region is: tight groups stay small and loose ones spread
  out. It cannot place rows, so it needs **Fit on** set to all rows and at most **Fit on up
  to** rows.
- **LocalMAP** (Wang et al., 2025) sharpens PaCMAP's local groups by letting go of neighbor
  pairs that end up farther apart on the map than the **Local distance limit**.
- **Exact** t-SNE measures every pair and takes at most 5,000 rows; Barnes-Hut approximates far
  points, more accurately and more slowly at a lower **Angle**.
- **Min. distance** must not be larger than **Spread**.
- **Seed**: the same seed gives the same map every run (UMAP then runs on one processor core).
  Set it to **Random** for a faster map that differs each time.
- PCA's note in the status bar gives the share of the variance on each chosen component. With
  **Write components** the columns are named "PC1", "PC2", … and the plot uses the chosen
  pair. The PCA Loadings panel is saved with the document.
- **PaCMAP** (Wang et al., 2021) balances neighbor, mid-near and further pairs, so it keeps
  local groups like UMAP while distances between groups mean more. With Jaccard or Hamming its
  neighbors are found exactly on the fingerprint bits, which takes longer on very large sets
  than a few PCA dimensions do.
- t-SNE on very large sets (over 50,000 rows) can take many minutes; UMAP and PaCMAP are
  quicker there.
- The map is written to new columns named after the method: "UMAP 1", "UMAP 2" (and
  "UMAP 3" in 3D). A map never replaces a column: the next UMAP map writes "UMAP 1 (1)",
  "UMAP 2 (1)", the one after "(2)", and so on. Rename them in the table if you like.
- A map needs at least 3 rows.
- The job shows as "Calculating UMAP…" (or PCA / t-SNE / PaCMAP) in the status bar until
  the plot opens. UMAP, t-SNE and PaCMAP run in a separate process, so the window stays
  responsive on large sets.

### Self-Organizing Map

Data ▸ Dimensionality Reduction ▸ Self-Organizing Map… trains a grid of cells. Similar
compounds land in the same or neighboring cells. Like the maps, it opens with the settings it
last ran with, and **Defaults** sets them back.

- **Input.** Compare by a fingerprint or by number columns; **Rows**; **Fit on** trains the map
  on those rows only and then places every row in its closest cell.
- **Grid.** Set rows × columns, or tick **Automatic** to let the size follow the number of
  compounds; **Cells** are **Square** or **Hexagonal** (six neighbors at the same distance, so
  groups spread more evenly); **Training passes** sets how long it trains.
- **More Options** (open): the fingerprint settings above, **Neighborhood** (**Gaussian**, whose
  pull fades with distance, or **Bubble**, which pulls every cell within the radius alike),
  **Start radius** (Automatic: half the grid) and **End radius** (default 0.7) in cells,
  **Start from** (**PCA plane**, the default, or **Random compounds**), **Wrap edges** (the grid
  joins opposite edges like a torus, so nothing crowds into the corners) and **Seed**.
- **Output.** Each compound's cell is written to the "SOM Row" and "SOM Column" columns. A
  panel draws the grid as squares or hexagons, colored by the count or by any number column,
  with **Higher is better** for activity. Click a cell to select its compounds; hover over it
  for its count and value.

### Generative Topographic Map

Data ▸ Dimensionality Reduction ▸ Generative Topographic Map… fits a smooth, square sheet of
nodes to the compounds (Bishop, Svensén and Williams, *Neural Comput.* 1998, 10, 215; used for
"chemography" by Gaspar et al., *J. Chem. Inf. Model.* 2015, 55, 84). Unlike the
Self-Organizing Map, each compound is spread over the nodes it fits, not put in one cell. That
gives smooth landscapes of a column over the map, and a measure of how well each compound fits
the map. The dialog opens with the settings it last ran with.

- **Input.** Compare by a fingerprint or by number columns, as for the maps above (reduced to
  at most 50 principal components); **Rows**; **Fit on** trains the map on those rows only
  and then places every row on it.
- **Method.** **Nodes** per side: tick **Automatic** for 10 below 300 compounds, 16 below
  3,000 and 25 above. **Basis functions** per side set how freely the sheet can bend (fewer is
  smoother; Automatic uses a quarter of the nodes, at least 3). **Iterations** is the most
  training passes (100); training stops sooner once the fit stops improving.
- **More Options** (folded): the fingerprint settings, **Basis width** (1.0, in spacings
  between basis functions; wider is smoother) and **Regularization** (0.1; higher keeps the
  sheet from folding into the data).
- **Output.** Three columns: "GTM X" and "GTM Y", each compound's place on the map (−1 to 1),
  and "GTM Log-likelihood", how well the map explains the compound, per dimension. A low
  value means the compound lies outside the space the map was trained on, which makes it a
  check for compounds placed with **Fit on**. The map itself is kept in the document, saved
  with it, and Undo removes it with the columns. A GTM Landscape plot opens on it.

**The GTM Landscape plot** (a plot kind; choose it from **Plot:** in the plot's settings once
a map exists) draws the map as a colored image, with each compound a point at its place.

- **Map.** Which trained map to draw, when the document has more than one (the newest is chosen).
- **Landscape.** **Density** shades where the compounds sit. A number column shows its
  mean over the compounds each node holds, weighted by how much of each compound the node
  holds (the training opens on an activity column when one is found). A text column shows the
  value most of each node's compounds have; mixed nodes are paler. Nodes with few compounds
  fade out, and empty ones are not drawn. **Colors** sets the color map; the landscape always
  spans its own range.
- **Color** and **Size** color and size the points as on a scatter plot. Click, box or lasso
  points to select their rows; hover over one for its values and log-likelihood.

---

## 15. Settings

The **Settings** menu (after Macromolecule) opens each group of settings in a dialog of its own:
**Theme…** (Ctrl+,), **Hotkeys…**, **External…** and **General…**, then **Status Bar**.
Changes preview immediately in every window. **OK** keeps them, **Cancel** puts everything
back, and **Reset** returns that dialog's settings to their defaults. Below them is what used
to be the Help menu: **User Guide**, **Log Folder** and **About**.

| Dialog | Contents |
|---|---|
| Theme | Theme: System (follows your computer's light or dark mode), Light and Dark (fixed, whatever the computer uses), Groovy (a random colorful palette) or a custom theme. Switching repaints every open window and panel in the new colors. Also the application and table fonts, text alignment and structure row height. |
| Hotkeys | Change any command's shortcut. |
| External | Paths to the optional programs chemsuite runs: AiZynthFinder (its Python and config), [OpenMM](#setting-up-openmm) (its Python, or on Windows a Python inside WSL) and, on Windows, the WSL that runs it, each with a **Test** button. Tools that need a missing program are grayed out, and their tooltip says why. |
| General | Whether to show the status bar, processor cores to use (0 = automatic: one per physical core but one, and fewer when the computer is short of memory, so the window and other programs stay responsive; a change applies the next time the workers start), 3D graphics (Automatic, Software on Windows and Linux, or Off; takes effect after a restart), recent files to keep, and whether a long job that finishes in the background posts a desktop notification. Background workers run at below-normal priority with one thread each. They close a minute after the last job ends, giving their memory back (the next job starts them again), and when chemsuite closes, even after a crash. |

Your IT department can fix some settings for everyone. Those are shown grayed out with a
padlock.

---

## 16. The log and reporting problems

- The **Messages** tab of the Log window (**Log**, at the right end of the menu bar; Window ▸
  Log on macOS) shows the session's log. Filter it by level and text. A job's messages
  (Show Messages on the Jobs tab) show with the job's name above the list; click it to show
  every message again. Double-click a message from a job to find the job on the Jobs tab.
- **Settings ▸ Log Folder** opens the folder of log files to attach to a bug report. If
  chemsuite ever closes by itself, `crash.log` there says where it was at that moment.
- **Settings ▸ User Guide** opens this guide; **About** (in the chemsuite menu on macOS) shows
  the version and license.
- If chemsuite hits an unexpected error, it shows a report you can copy. The details are also
  written to the log.

---

## 17. Keyboard shortcuts

You can change any of these in Settings ▸ Hotkeys.

| Action | Shortcut |
|---|---|
| New / Open / Save / Save As | Ctrl+N / Ctrl+O / Ctrl+S / Ctrl+Shift+S |
| Close (clear the document) | Ctrl+F4 |
| Undo / Redo | Ctrl+Z / Ctrl+Y |
| Copy / Paste | Ctrl+C / Ctrl+V |
| Search Rows | Ctrl+F |
| Find and Replace | Ctrl+H |
| Structure Search | Ctrl+Shift+M |
| Select All / Invert / Clear Selection | Ctrl+A / Ctrl+I / Ctrl+Shift+A |
| Delete Selection | Ctrl+Backspace |
| Add Selected to List | Ctrl+L |
| Calculate Descriptors | Ctrl+Shift+D |
| Scatter plot | Ctrl+Shift+P |
| Filters | Ctrl+Shift+F |
| Browser | Ctrl+Shift+B |
| Settings ▸ Theme | Ctrl+, |
| Larger / Smaller Structures / Reset Structure Size | Ctrl+= / Ctrl+- / Ctrl+0 |
