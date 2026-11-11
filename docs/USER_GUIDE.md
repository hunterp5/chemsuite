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
13. [SAR](#13-sar)
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
  | File | New, New from Shown Rows, Open, Save, Export ▸ Table, Selected Rows or Conformers (SDF), Close |
  | Edit | Undo, Copy and Paste, Edit Structure, selecting rows, Lists, Delete Rows, Settings |
  | View | Search Bar, Arrange Panes, Jobs, Log, Status Bar |
  | Tools ▸ Table | Add compounds, rows and columns; Split and Join columns; Calculator; random numbers and rows |
  | Chemistry | Sketcher, Calculate Descriptors, Standardize, Predict pKa, Structure Viewer, Generate Conformers, Split Reactions, Enumerate Products, R-Group Decomposition |
  | Data | Search ▸, Plot ▸, Filters, Browser, then Bin Values, Statistics, Similarity Search, Diverse Subset, Cluster, SAR ▸, Scoring ▸, Chemical Space ▸ |
  | Data ▸ Search | Search Rows, Find and Replace, Advanced Search, Structure and Reaction Search, Clear Highlighting, Show All Rows |
  | Data ▸ Plot | Scatter, Line, Bar Chart, Histogram, Box Plot, Density, Radar, 3D Scatter, 3D Landscape, MedChem ▸ |
  | Protein | Protein Viewer: the workspace for proteins, with its own menus |
  | Help | This guide, the log folder, About |

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
    Every panel (plots, Filters, Advanced Search, Structure Viewer) has the same header
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
    header's empty space) to move the panel. Double-click the current tab to rename a plot.
    When the header is narrow, long tab titles are shortened with "…"; hover over a tab to
    see its full title.
  - Plots leave a little space between the axis titles and the panel's edges.
  - **Arranging panels.** For stacked plots beside the table, open a second plot (it joins
    the first as a tab), then choose Arrange ▸ Below ▸ *the first plot*. You can also drag a
    panel by its header onto the top, bottom or side of another panel to split its space, or
    onto its middle to tab them. **View ▸ Arrange Panes** lays out every docked panel at once:
    **As Tabs**, **Side by Side**, **Stacked** or as a **Grid**.
  - **The Filters panel stays put.** It always sits at the right edge of the window at a fixed
    width, whatever you do with the other panels: new panels, arrangements and moves all go
    to its left. It has no Arrange or Undock button and cannot be dragged elsewhere.

- **Status bar**, on one line:
  - At the left: what the current background job is doing, such as "Generating conformers…
    (120/500)" or "Calculating UMAP…", followed by "· 2 queued" or "· 1 failed" when that
    applies. Click it to open the Jobs window.
  - Next to it: the latest message, such as "Opened compounds.sdf" or "Cancelled:
    Generating conformers; partial results kept".
  - At the far right: the row and selection counts ("200 rows · 3 selected").

  View ▸ Status Bar hides or shows it.

### Tool dialogs

Tool dialogs do not block the window. While one is open you can still select rows, change
filters and scroll the table, and the dialog's row counts follow along. Opening another
document closes any open tool dialogs.

Every tool dialog uses one layout, and so do Split Column, Join Columns, Number Format, Align
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

Use File ▸ Open (Ctrl+O) or File ▸ Open Recent.

| Format | Extensions | Notes |
|---|---|---|
| SD file | `.sdf`, `.sd`, `.mol` | Data fields become columns. A record whose structure does not parse is **kept**: its other fields still load, and the structure cell is marked invalid. Nothing in your file disappears silently. |
| CSV / TSV | `.csv`, `.tsv`, `.tab`, `.txt` | A column headed SMILES, Structure, Mol, Molecule, Canonical_SMILES or Isomeric_SMILES becomes a structure column. A column whose first 20 filled cells are all reaction SMILES (such as USPTO's `rxn_smiles`) becomes a [reaction column](#reactions-in-the-table), even if its name ends in "smiles". |
| SMILES list | `.smi`, `.smiles` | One SMILES per line, with an optional name. |
| Excel | `.xlsx` | The first sheet. Structure and reaction columns are found as for CSV. |
| RD file | `.rdf`, `.rd` | Reactions with their data: a **Reaction** column, and one column per data field (`$DTYPE`). Registry numbers go in an **ID** column. As with SD files, a record whose reaction does not parse is kept and marked invalid. Molecule records (`$MFMT`) go in a **Structure** column. |
| RXN file | `.rxn` | One reaction, as a one-row table. |
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
  - open plots and panels with their settings. Analysis panels (Free-Wilson, Hansch, Matched
    Pairs, Activity Cliffs) rerun their analysis on the same rows when the document opens
    (they show "Recalculating…" meanwhile); Hansch keeps its substituent constants, so the CSV
    file is not needed again;
  - the selected rows, the current cell and where the table was scrolled to;
  - the workspace layout: where each panel is docked (beside, stacked, tabbed, left, right or
    below the table, floating), the size of each pane (where you dragged the dividers
    between the table and the panes, and between panes), which panes are open (Filters,
    Advanced Search, the panels), column widths, and the search bar with its text and
    options.
- A document reopens exactly as you left it: a pane that was closed stays closed, even when
  its filters or search are hiding rows (the status bar still says how many rows are shown).
- Saving is safe. The new file is written next to the old one and swapped in only once it is
  complete, so a crash or a full disk never damages your previous version.
- **File ▸ Export ▸ Table…** writes SDF, CSV/TSV, SMILES, Excel, an RD file or an RXN file.
  **File ▸ Export ▸ Selected Rows…** writes only the selected rows.
  - An **RD file** holds the first reaction column as RXN blocks and every other column as
    data. Reactions read from SMILES have no drawing, so they are laid out as the sketcher lays
    them out (using every processor core). Blank values are left out.
  - An **RXN file** holds one reaction: export a single selected row.
  - CSV, TSV and Excel write reactions as reaction SMILES.

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
  Structure Viewer at that row (in 2D, or as it is when it is already open). Structure cells
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
compounds, so plots and the Structure Viewer are left alone. **Delete** removes the selected
columns in one step you can undo. Sorting is in the header's right-click menu.

**Moving columns.** Drag a column header left or right to move the column (one undoable
step, saved with the document). Hold it near the table's left or right edge to scroll the
table: the closer to the edge, the faster it scrolls, and faster still if you drag past the
edge, so a column can go anywhere in a wide table. After the drop the table stays where you
scrolled to, the moved column stays in view and, if it was selected, stays selected.

Right-click a column header for:

- **Sort ▸ Ascending / Descending.** **Then Ascending / Then Descending** add a second (third…)
  sort key; **Clear** returns rows to their original order.
- **Filter.** Adds a filter card for the column (see [Filters](#7-filters)).
- **Search.** Opens the search bar set to this column, ready to type.
- **Select ▸** picks shown rows by this column's value, top to bottom: **Empty**, **Not Empty**,
  **First Occurrence** (the first row of each value), **Duplicates** (every row whose value
  appears in a row above it, so deleting them keeps one of each) and **All Duplicated Values**
  (every row whose value appears more than once). Structures count as the same when they are
  the same molecule, however they were written.
- **Rename…**, **Duplicate Column** (a copy, with its formatting, to the right),
  **Delete Column**. Delete acts on all selected columns when you right-click one of them.
- **Logarithmic (log₁₀)** (Number columns). Converts the values to log₁₀; choose it again to
  convert them back (10ˣ). The header shows "(log₁₀)". It needs every value above zero. One
  step you can undo.
- **Number Format…** (Number columns). How every value is shown:

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
- **Colour Column…** Colours cell backgrounds: **Low → high** or **Low → middle → high**
  gradients for numbers, or **One colour per value** for any column (the same value always gets
  the same colour). Range ends left empty follow the data. Choose **None** to remove the
  colours.
- **Align Drawings…** (Structure columns). Type, draw or take from the current row a core
  (SMILES or SMARTS); every structure that contains it is drawn with the core in the same place
  and orientation, so a series lines up down the column. Structures without the core keep
  their usual drawing. Only the drawings change, not the structures. **Stop Aligning
  Drawings** brings the usual drawings back. The core is saved with the document's layout.
- **Hover over a structure cell** for a large drawing of the structure with its SMILES below.
- **Freeze Columns up to Here** keeps that column and every column left of it in place while
  the rest of the table scrolls sideways; **Unfreeze Columns** lets them scroll again. Drag the
  columns you want to keep in view (ID, Structure) to the front first. Frozen columns are saved
  with the document.
- Number formats, Logarithmic and colours are saved in the document and can be undone.
- **Hover over a column's header** for a summary of the column: how many cells are filled and
  blank, and for numbers the range, mean ± standard deviation and median, or for text the
  number of distinct values and the commonest one.
- **Column Type ▸ Text / Number / Structure / Reaction.** Converting to Number warns you if
  some values are not numbers; those cells would be cleared. Converting to Structure or
  Reaction reads each cell's text; cells that do not read are kept and marked invalid.

From Tools ▸ Table:

| Command | What it does |
|---|---|
| Add Compound… | Opens the sketcher; draw (or paste) a structure, and **Add** appends it as a new row. |
| Add Rows… | Adds empty rows. |
| Add Column… | Adds a column. |
| Split Column… / Join Columns… | Each with a live preview of the first rows shown. Choose **All**, **Shown** or **Selected** rows; other rows stay blank. Split finds the delimiter itself (**Auto**) or uses one you choose, and keeps every field or only the first, last, largest or smallest; Join can skip blank cells so no delimiters are doubled. |

### Rows

- **Delete Rows (Ctrl+Backspace)** deletes the selected rows. You can undo it.
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
  order, into a new untitled document. Column types are kept.

---

## 4. Structures

- Structures are edited in the **sketcher** only. **Double-click** a structure, press **F2**
  or start typing in the cell, choose **Edit ▸ Edit Structure… (Ctrl+E)**, or right-click and
  choose **Edit in Sketcher…**. OK writes the drawing back as one step you can undo. The
  drawing is stored as a molfile, which keeps your layout and wedges.
- To bring in SMILES, InChI or a molfile, paste it onto the drawing (Ctrl+V), or paste it
  straight into the cell from the table.
- A cell whose text is not a valid structure opens an empty sketcher. Its text stays in the
  cell until you draw a replacement.
- The right-click menu of a structure or reaction cell has **Copy As** (Molfile, SMILES,
  InChI, InChIKey, Image) rather than a plain Copy, so you choose what is copied.
- **Edit ▸ Copy Structure As ▸ Molfile / SMILES / InChI / InChIKey / Image** copies the current
  structure (or the selected structures in its column). A molfile pastes into drawing programs.
- A structure that cannot be read is shown as **invalid**. Its original text is kept, so you
  can fix it.
- Rows are taller when the table has a structure column. Change the height under Settings ▸
  Appearance.

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

Chemistry ▸ Split Reactions… writes the parts of each reaction as structure columns.

- **One per role**: a Reactants, an Agents and a Products column, each holding all of that
  role's components.
- **One per component**: Reactant 1, Reactant 2 … Product 1 …, in the order the reaction
  lists them.
- **Keep atom map numbers** is off by default, so the same molecule compares equal across
  reactions in structure search, duplicates and descriptors.

### Enumerate Products

Chemistry ▸ Enumerate Products… applies a reaction template to building blocks.

- **Template**: an atom-mapped reaction SMARTS, such as
  `[C:1](=[O:2])[OH].[N;H2,H1;!$(NC=O):3]>>[C:1](=[O:2])[N:3]` (amide coupling), or
  **Draw…** and map the atoms with the sketcher's Atom Map tool. Up to four reactants.
- **Reactant 1, 2 …**: a structure column of this table for each, or **From a File…** (SDF,
  CSV or any table chemsuite opens).
- **Combine: Row by row** reacts reactant 1 and reactant 2 from the same row and writes the
  product to a **Product** column, with **Products Found** counting how many products the row
  gave (a building block with two reactive sites gives two).
- **Combine: Every combination** reacts every distinct building block with every other and
  opens the virtual library in a new window, one product per row with its building blocks.
  Up to 2,000,000 combinations.
- Cancelling keeps the products already made.
- Checked against hand-worked products of amide coupling, Suzuki coupling, reductive
  amination and Boc removal, including that stereocentres away from the reaction are kept.

### Filtering by reaction

**Add Filter** on a reaction column gives a reaction card: the same query and **Find** choices
as Reaction Search, applied as a filter together with the other cards.

### The sketcher

Opened by double-clicking a structure cell, Ctrl+E, F2, **Add Compound…**,
**Draw…** in the structure fields of Structure Search and Similarity Search, or
blank with **Chemistry ▸ Sketcher…**, whose **Add to Table** adds the drawing as a new,
selected row.

Drawing tools are along the top, elements down the left and rings along the bottom. OK is
available only when every atom has a valid valence; a line under the drawing says what is
wrong. Undo, Redo, Cut, Copy, Paste, Select All and Clear are on the keyboard (Ctrl/⌘+Z, Y, X,
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
reaction cannot be OK'd into a structure cell; the sketcher asks you to delete the arrow and
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
    (Ctrl+Shift+A)**.
- The selection is shared by the table, the plots and the tools:
  - selecting points in a plot selects the rows;
  - selecting rows highlights their points.
- Many tools can run on all rows, the shown rows or just the selected rows.

### Browser

Data ▸ Browser (Ctrl+Shift+B) opens the [Structure Viewer](#structure-viewer) in 2D: it
steps through rows one at a time with a large drawing, keeping the table's current row in
step. If a Structure Viewer is already open, it is brought forward and switched to 2D.

### Named lists

Lists are saved sets of rows, such as "actives", "to make" or "sent for assay". They are saved
with the document.

- **Edit ▸ Lists ▸ Add Selected to List… (Ctrl+L)** adds the selected rows to a new or existing
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
  | `-Series:B` | Series does not contain B; a `-` excludes any term |
  | `"Mol Wt"<400` | quote a column name that has spaces |

  Column names ignore case. A word before a colon that is not a column name is searched for
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
- **Status bar.** The job's progress shows at the far left of the status bar: "Generating
  conformers… (120/500)", or just "Calculating UMAP…" when there is no count. Click it to open
  the Jobs window.
- **View ▸ Jobs** lists running, queued and finished jobs with their elapsed time and result.
  From there you can:
  - cancel a job;
  - remove a queued job, or move it up or down;
  - **Remove** finished, failed or cancelled jobs from the list (select one or several, or
    press Delete); **Clear Finished** removes them all;
  - copy an error.

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
  | Similarity Search | Free-Wilson |
  | R-Groups (numbered core) | Hansch |
  | Diverse Subset (the picks made so far) | Matched Pairs |
  | | Activity Cliffs |

  The tools in the right-hand column have no meaningful result until the whole run is complete.
- **Safety.**
  - Closing a document cancels its queued jobs and asks about running ones.
  - Quitting with jobs running asks first.
  - Results from a job are never written into a different document.
- **Long runs.**
  - The computer is kept awake while jobs run.
  - When a job that took more than 30 seconds finishes or fails while chemsuite is not in
    front, the dock icon (taskbar on Windows) asks for attention and a desktop notification
    names the job and its result. Settings ▸ General ▸ **Notify when a long job finishes** turns
    the notification off.

---

## 9. Plots

Open a plot from the **Plot** menu. It opens to the right of the table.

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

### MedChem plots

The Data ▸ Plot ▸ MedChem menu has three plots. Each one first calculates any properties it needs that
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
| 3D Shape | PMI 1–3, NPR1, NPR2, PBF, Radius of Gyration, Asphericity, Eccentricity, Spherocity Index, Inertial Shape Factor |
| Identifiers | Canonical SMILES, InChIKey, InChI, Formula, Murcko Scaffold, Generic Scaffold |

- **pKa-based values.** The Ionisation group, cLogD 7.4, CNS MPO, BBB Score and AB-MPS use
  MolGpKa's predicted pKa values (see [Predict pKa](#predict-pka)). For CNS MPO and BBB Score, a
  compound with no basic site is scored as if its most basic pKa were 2.
- **SAScore** is the synthetic accessibility score of Ertl and Schuffenhauer (2009): 1 (easy
  to make) to 10 (very hard), from how common the molecule's fragments are in known compounds
  and a penalty for size, stereocentres, macrocycles and unusual ring fusions. Most drug-like
  compounds score 2-4.
- **3D shape descriptors** use a 3D structure when the cell has conformers.
- **Fingerprints.** Tick any of the fingerprints the similarity tools use (ECFP4, ECFP6, FCFP4,
  RDKit paths, atom pairs, topological torsions, MACCS keys). Each becomes a text column
  ("ECFP4", "MACCS Keys"…) holding the bits in hexadecimal, bit 0 first (the high bit of the
  first byte), ready to export for modelling. A row with no valid structure is left blank.
- **Cancelling** keeps the rows already calculated; a fingerprint is written only when every
  row has one.

### Predict pKa

Chemistry ▸ Predict pKa predicts the pKa of every ionisable site with **MolGpKa**, a
graph-convolutional network (Pan et al., J. Chem. Inf. Model. 2021). Its authors report a mean
absolute error of about 0.6 pKa units on their literature test set. The model ships with
chemsuite and runs on the CPU, at roughly 2 ms per molecule.

- **Input.** A structure column and the rows. Salts are reduced to the largest fragment and
  charges are neutralised before prediction, so a sodium salt gives its parent acid's values.
- **pH.** Charge, % neutral, cLogD and the major species are given at this pH (default 7.4).
  The pKa values do not depend on it.
- **Columns** (choose any):

  | Column | Meaning |
  |---|---|
  | pKa (most acidic) | The lowest acid pKa. Blank without an acidic site. |
  | pKa (most basic) | The highest pKa of a protonated base (the strongest base). Blank without a basic site. |
  | pKa Values | Every site, lowest pKa first, each marked acid or base, for example "4.47 (acid), 9.05 (base)". |
  | Charge @ pH, % Neutral @ pH | Mean net charge, and the share of molecules with no site ionised. Zwitterions count as ionised. |
  | cLogD *pH* | Crippen cLogP of the neutral parent plus log10(fraction neutral). |
  | Major Species @ pH | The dominant protonation state, as a structure column. |

  Each run writes new columns and never changes an earlier run's. When a name is already
  taken, all the columns of the run get the same suffix, (2), (3)…, so they read as a set:
  `pKa (most acidic) (2)` and `cLogD 7.4 (2)` come from the same run.

- **How sites combine.** Each site is predicted, and ionised, independently (Henderson–
  Hasselbalch per site). This is a good approximation when sites are far apart and a rough one
  when they interact.
- **Limits.**
  - Only groups in MolGpKa's site templates are found; methylamine, for example, has no
    matching template.
  - Very small molecules lie outside the training data and can be badly off: triethylamine
    comes out at 1.7 against a measured 10.75, while larger tertiary amines are close.
  - Because only the neutral form is assumed to partition, cLogD of zwitterions and of
    strong acids and bases comes out far lower than measured: glycine gives about −8.5.
  - Treat the values as screening estimates, not measurements.
- **Cancelling** keeps the rows already predicted.

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
| Protonate at pH | Off | Each ionisable site in its majority form at the pH (default 7.4), from MolGpKa's predicted pKa values (see [Predict pKa](#predict-pka) and its limits). Carbon acids are ignored, so no carbanions are made. Kept counter-ions stay as they are. |

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
  - Undo restores the table in one step.
- **Next steps.** For docking, follow with [Generate Conformers](#generate-conformers) for 3D
  coordinates. Standardized structures carry no explicit hydrogens; docking programs add them.
- **Cancelling** keeps the rows already standardized.

### Generate Conformers

Generates 3D conformers with RDKit's ETKDG and optimises them with MMFF94 (UFF where MMFF has no
parameters).

- **Input:** the structure column and **All**, **Shown** or **Selected** rows.
- **Method:** **Conformers** per molecule, **Distinct above** (the RMSD for two conformers to
  count as distinct), the force field (**Optimise with**), and an **Energy window** (tick it to
  drop conformers more than that many kcal/mol above the lowest). The random **Seed** is under
  **More Options**.
- **Output.** The conformers go into a new 3D structure column, lowest energy first, with
  columns **Conformers** (the count) and **Conformer Energy Range (kcal/mol)**.
  - Energies compare conformers of the same molecule only.
  - Running again writes new columns ("SMILES 3D (2)"); earlier conformers are kept.
- **Cancelling** keeps the molecules already done.

### Export conformers (SDF)

Writes every conformer as its own SD record, carrying the row's other columns. Choose **All**,
**Shown** or **Selected** rows (Selected when rows are selected), then **Export** asks for the
file.

### Structure Viewer

The Structure Viewer is a panel that shows one row's structure at a time, drawn in **2D** or
in **3D**, starting at the table's current row. Open it in 2D with **Data ▸ Browser
(Ctrl+Shift+B)** or in 3D with **Chemistry ▸ Structure Viewer (Ctrl+Shift+3)**; either
brings an open viewer forward. The **2D | 3D** switch at the left of its header changes only
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
- **Structure column:** chosen in the options menu. Opened in 3D, a "3D" column is chosen
  first.
- **Export** saves the view (with its caption) as a PNG.
- The panel, its settings, 2D or 3D, the caption, the row and conformer shown and the 3D view
  are saved with the document. A 3D Viewer panel saved by an earlier version reopens as a
  Structure Viewer in 3D.

**In 2D** the structure is drawn as in the table, larger; ↑ and ↓ step through the rows too.

**In 3D:**

- A row with only a 2D structure gets 3D coordinates in the background (one ETKDG conformer);
  they are not written to the table. They are made ahead for the rows either side of the one
  shown, in 2D as well, so switching to 3D or stepping on usually finds them ready.
- **Style:** Ball and Stick, Stick, Spacefill or Wire. **Hydrogens:** All, Polar (on atoms
  other than carbon) or None. **Color:** by element, or each conformer's carbons in a colour of
  its own. These options show in the options menu in 3D.
- **Conformers:** for a row with several, ‹ and › in the panel's header (or ↑ and ↓) step
  through them; the counter shows which one. **Overlay** draws every conformer at once as thin
  sticks with polar hydrogens. **Align Conformers** (on by default) fits them onto the first on
  heavy atoms.
- **Mouse:** drag to turn, right-drag, middle-drag or Shift-drag to move, wheel to zoom. Pointing
  at an atom names it (for example "N4"). Ctrl-click (⌘-click on macOS) an atom to add or remove
  the row from the selection.
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

### Proteins

**Protein ▸ Protein Viewer** (Ctrl+Shift+R) opens a workspace for proteins in its own
window (dock it beside the table with the button in its title bar if you prefer). It has its
own menus:

- **File:** **Open Protein…** (PDB or mmCIF files, gzipped or not), **Fetch from PDB…**
  (download an entry by its ID, such as 1ATP), **Save Selection As…** (PDB or mmCIF) and
  **Close Viewer**. With no protein open, the middle of the window offers Open and Fetch.
- **Edit:** **Copy Selection As ▸ PDB / mmCIF**, **Select All**, **Clear Selection**, **Save
  Selection…**.
- **View:** Show, Hide, Colour and Label, Hydrogens, Surface Transparency, Frame Selection, Reset
  View, and whether the proteins list and the sequence show.
- **Pocket:** Interactions, Focus, Measure and **Extract Ligands to Table** (rows for the
  ligands of proteins that have none yet).

Proteins are saved inside the document. Adding one also adds a row for each of its ligands to
the table, in one step you can undo:

| Column | What it holds |
|---|---|
| Protein | Which protein the ligand belongs to (a second copy of 1ATP becomes 1ATP-2) |
| Ligand ID | Chain, residue name and number, such as `E/ATP 355` |
| Ligand | The ligand as a structure, in its pose from the file (3D kept) |
| Bond Orders | Where the bond orders came from: the PDB's chemical component dictionary (CCD), a SMILES, or **Guessed from geometry**, which is worth checking |

Make a ligand row current in the table and the viewer shows that ligand in its pocket.

**The viewer.** The proteins list (tick a protein to show or hide it; click a chain or ligand to
select it), the 3D view and the sequence underneath.

- **Selecting.** Click selects a residue, Alt-click an atom, double-click its chain; Shift adds
  and Ctrl (Cmd) toggles; click empty space or press Esc to clear. Click or drag across letters
  in the sequence. Type in the selection field (point at it for the syntax): `chain A`,
  `resi 10-50`, `resn HIS+ASP`, `name CA`, `ligand`, `water`, `byres (polymer within 5 of
  ligand)` (whole residues near a ligand; `byres` takes what follows it, hence the brackets),
  `and`, `or`, `not`. **Select ▸** has the usual picks, **Around ▸** (the residues or
  atoms within 3, 5, 8 or 10 Å of the selection, added to it) and **Save Selection…**; saved
  selections are used as `$name`, and `$selection` is what is selected now.
- **Show, Hide, Colour** act on the selection: wire, sticks, ball and stick, spacefill,
  cartoon, surface; colour by element, chain, secondary structure, B-factor, hydrophobicity or
  a colour you choose.
- **Label ▸ Residues / Atoms** labels what is selected (a residue on its Cα, a ligand in its
  middle); **Clear Labels** removes them. Labels stay legible: nearer ones win where they would
  overlap, faded and sliced-away atoms lose theirs, and at most 300 show at once.
- **Hydrogens** (all, polar or none) and **Transparency** (of surfaces, from opaque to 80%)
  apply to every protein. A hydrogen is drawn as the atom it is bonded to is shown, so showing
  a residue as sticks shows its hydrogens too. Both are in the View menu as well.
- **Interactions** draws dashed lines from the ligands (or the selection) to what they touch:
  hydrogen bonds (blue), salt bridges (magenta), π-stacking (green), cation–π (orange),
  hydrophobic contacts (grey), halogen bonds (cyan), metal contacts (purple) and water bridges
  (light blue). Point at an atom to read its interactions in the status line. Structures from
  the PDB rarely carry hydrogens, so these are judged from heavy atoms.
- **Focus ▸** a ligand shows its pocket (residues within 5 Å as sticks, the waters and ions
  touching it) with its interactions and a see-through skin over the pocket coloured by
  hydrophobicity, fades the rest of the protein and hides other small molecules. A slice of
  depth around the pocket follows the camera, with the back of the pocket fading into the
  background. **Pocket Surface** in the same menu turns the skin off or on. **Clear Focus**
  undoes the fading.
- **Measure ▸ Distance / Angle / Dihedral**, then click two, three or four atoms. Esc stops
  measuring.
- **Models.** A file with several models (an NMR ensemble) gets a bar under the view: step
  with the arrows or `[` and `]`, type a model's number, drag the slider, or press Play to loop
  through them. The selection, colours and measurements stay as the atoms move. When several
  shown proteins have models, the list at its left picks which one the bar steps. The model
  shown is saved with the document.
- Right-click an atom to select its residue, chain or surroundings, focus a ligand, switch an
  alternate location, copy a ligand's SMILES, **Add Ligand to Table** (a cofactor, say, or a
  ligand whose row you deleted), **Label ▸** its residue or the atom, start a **Measure From
  Here ▸** distance, angle or dihedral with it as the first atom, or **Copy Selection As ▸
  PDB / mmCIF**.

Surfaces, interactions and the cartoons of very large proteins are built in the background,
so the viewer never stops answering: each appears as soon as it is ready, and until then the
previous one stays on screen.

Everything shown, selected, focused and measured is saved with the document. Fetch from PDB
only sends the ID, over HTTPS; when IT turns downloads off it says so instead.

### Retrosynthesis

Chemistry ▸ Retrosynthesis plans synthetic routes to each structure with
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
3. In **Settings ▸ External Tools**, set **AiZynthFinder Python** to that environment's
   `python` and press **Test**: it reports the AiZynthFinder version it found.
4. Press **Download…** on the **AiZynthFinder config** row and choose a folder. chemsuite
   fetches AiZynthFinder's public USPTO policies, ring-breaker policy, filter and the ZINC stock
   (about 1.5 GB, from Zenodo and figshare) in the background (follow it in View ▸ Jobs), writes
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

Tools ▸ Table ▸ Calculator… (Ctrl+Shift+=) adds a column worked out from a formula. The
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

### Random Numbers

Tools ▸ Table ▸ Random Numbers… fills a column with random values:

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

Tools ▸ Table ▸ Random Rows… picks a random sample from all, the shown or the
selected rows, either a number of rows or a percentage. It then **selects them**, **adds them
to the selection**, or saves them as a **new list** (for example a test set). A seed picks the
same rows again.

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

### Similarity Search

Data ▸ Similarity Search… scores **All**, **Shown** or **Selected** rows against one of two
references (**Compare With**):

- **A Structure** you type or paste, draw (the pencil), or take from the current row
  (**Current Row**); or
- **Selected Rows**. Each row gets its similarity to the nearest selected compound, and
  a column naming that compound. A selected row is never its own nearest neighbour.

Options:

- **Sort: Most similar first** sorts the table by the new column.
- **Select: At least** selects the rows at or above a similarity.

Cancelling keeps the rows already scored.

### Diverse Subset

Data ▸ Diverse Subset… picks *N* compounds that cover the chemistry, using MaxMin:
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

## 13. SAR

### R-Group Decomposition

Chemistry ▸ R-Group Decomposition… splits a series onto a common core.

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
  and rows needing it are left out.
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
- Transformations (A → B) are summarised over all their pairs, with the activity shift.

### Activity Cliffs (SALI)

Data ▸ SAR ▸ Activity Cliffs… lists similar compounds with very different activity.

    SALI = |Δactivity| / (1 − similarity)

- Only pairs at least **Min. similarity** alike are considered.
- **Identical: Skip identical fingerprints** ignores salts and stereoisomers, which would otherwise divide
  by zero.

---

## 14. Chemical space

### Dimensionality Reduction

Data ▸ Chemical Space ▸ Dimensionality Reduction… places compounds on a 2D or 3D map by
fingerprints, number columns, or both. Choose a **Fingerprint** (or **None**) and tick any
**Columns**: a fingerprint alone maps by its bits, columns alone by their z-scores, and a
fingerprint with ticked columns uses the two together, for rows that have all of them: the
fingerprints are reduced to their principal components, the columns are z-scored, and the two
are weighted equally so neither drowns out the other (distances are Euclidean).
**Dimensions: 2D** writes two columns and opens a scatter plot;
**3D** writes three and opens a 3D Scatter (a 3D PCA map is always PC1, PC2 and PC3).

| Method | Options | More Options (folded) |
|---|---|---|
| PCA | **Components**: PC1 and PC2 (the default), PC1 and PC3, PC2 and PC3, or PC3 and PC4 | — |
| UMAP | Neighbours, Min. distance | **Distance** (Jaccard, the default, Dice, Rogers-Tanimoto, Hamming, Euclidean, Manhattan, Cosine, Correlation), **Spread**, **Epochs**, **Fit on up to** (rows; larger sets are fitted on a sample of this size and the rest placed onto that map; default 50,000), **Seed** |
| PaCMAP | Neighbours (default 10) | **Mid-near pairs** (per neighbour pair; default 0.5), **Further pairs** (default 2), **PaCMAP distance** (Jaccard, the default, Hamming, Euclidean, Manhattan, Angular), **PCA dimensions** (default Off), **Seed** |
| t-SNE | Perplexity | **Learning rate** (Automatic or a number), **Iterations**, **Exaggeration** (early exaggeration), **PCA dimensions** (default Off: Jaccard on fingerprint bits; a number runs t-SNE on that many principal components), **Start from** (Random, the default, or PCA), **Seed** |

- **Jaccard** (Tanimoto), **Dice**, **Rogers-Tanimoto** and **Hamming** compare fingerprint
  bits, so they need a fingerprint with no columns ticked (and, for PaCMAP, **PCA dimensions**
  Off). When they don't apply the distance switches to Euclidean, and back once they do.
- **PCA dimensions** is Off by default, so t-SNE and PaCMAP map the features themselves.
  Setting a number first reduces the features to that many principal components, which is
  faster on large sets (bit distances are then not available).
- **Min. distance** must not be larger than **Spread**.
- **Seed**: the same seed gives the same map every run (UMAP then runs on one processor core).
  Set it to **Random** for a faster map that differs each time.
- PCA's note in the status bar gives the share of the variance on each chosen component.
- **PaCMAP** (Wang et al., 2021) balances neighbour, mid-near and further pairs, so it keeps
  local groups like UMAP while distances between groups mean more. With Jaccard or Hamming its
  neighbours are found exactly on the fingerprint bits, which takes longer on very large sets
  than a few PCA dimensions do.
- t-SNE on very large sets (over 50,000 rows) can take many minutes; UMAP and PaCMAP are
  quicker there.
- The map is written to new columns named after the method: "UMAP 1", "UMAP 2" (and
  "UMAP 3" in 3D). A map never replaces a column: the next UMAP map writes "UMAP 1 (1)",
  "UMAP 2 (1)", the one after "(2)", and so on. Rename them in the table if you like.
- The plot opens coloured by the activity column when there is one; change it with the cog.
- A map needs at least 3 rows.
- The job shows as "Calculating UMAP…" (or PCA / t-SNE / PaCMAP) in the status bar.

### Self-Organizing Map

Data ▸ Chemical Space ▸ Self-Organizing Map… trains a grid of cells. Similar compounds land in
the same or neighbouring cells.

- **Grid.** Set rows × columns, or tick **Automatic** to let the size follow the number of
  compounds; **Training passes** sets how long it trains.
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
| External Tools | Paths to optional programs (Open Babel, CONFORGE, Gnina, AmberTools, MAFFT, AiZynthFinder), each with a **Test** button. Tools that need a missing program are greyed out, and their tooltip says why. |
| General | Processor cores to use (0 = automatic: all but one, so the window stays responsive), 3D graphics (Automatic, Software on Windows and Linux, or Off; takes effect after a restart), recent files to keep, and whether a long job that finishes in the background posts a desktop notification. Background workers run at below-normal priority with one thread each, and close when chemsuite closes, even after a crash. |

Your IT department can fix some settings for everyone. Those are shown greyed out with a
padlock.

---

## 16. The log and reporting problems

- **View ▸ Log** shows the session's log. Filter it by level and text.
- **Help ▸ Show Log Folder** opens the folder of log files to attach to a bug report.
- **Help ▸ chemsuite User Guide** opens this guide; **About chemsuite** (in the chemsuite menu
  on macOS) shows the version and licence.
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
| Scatter plot | Ctrl+Shift+P |
| Filters | Ctrl+Shift+F |
| Browser (Structure Viewer in 2D) | Ctrl+Shift+B |
| Structure Viewer in 3D | Ctrl+Shift+3 |
| Settings | Ctrl+, |
