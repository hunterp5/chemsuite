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
10. [Chemistry: descriptors, pKa, standardizing and conformers](#10-chemistry-descriptors-pka-standardizing-and-conformers)
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
  | Chemistry | Sketcher, Calculate Descriptors, Standardize, Predict pKa, Generate Conformers, Split Reactions, Enumerate Products, R-Group Decomposition |
  | Data | Search ▸, Plot ▸, Filters, Browser, then Bin Values, Statistics, Similarity Search, Diverse Subset, Cluster, SAR ▸, Scoring ▸, Chemical Space ▸ |
  | Data ▸ Search | Search Rows, Find and Replace, Advanced Search, Structure and Reaction Search, Clear Highlighting, Show All Rows |
  | Data ▸ Plot | Scatter, Line, Bar Chart, Histogram, Box Plot, Density, Radar, 3D Scatter, 3D Landscape, MedChem ▸ |
  | Help | This guide, the log folder, About |

  - Every command that a tool adds appears in these menus. If a command is greyed out, its
    tooltip says why: for example, there is no structure column, or an external program is
    missing.
- **Search bar.** Hidden until you ask for it with Data ▸ Search ▸ Search Rows (Ctrl+F). See
  [Searching](#6-searching).
- **The table.** One row per compound. Structure columns show drawings.
- **Docked panels.**
  - Filters sits at the right edge at a fixed width, and the Browser docks at the
    side. Advanced Search opens as its own floating window.
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

Tools are moving to one dialog layout. Cluster, Dimensionality Reduction and Predict pKa use it now,
and the other tools will follow. The window title names the tool; every dialog on the layout
reads the same way, top to bottom:

- **Input.** What to work on: structures or number columns, and **Rows**.
  - Rows is **All · Shown · Selected**, with each option's count. The counts update as you
    select rows or change filters.
  - If you choose **Selected** with nothing selected, the dialog says so above its buttons. Select rows in the
    table and the dialog is ready to run.
- **Method.** The tool's settings. Settings that only apply to one method appear when you
  choose that method. Hover over a field or its label for an explanation.
- **More Options.** Optional settings whose defaults suit most runs sit in a folded section.
  Click its arrow or title to open it.
- **Output.** The name of the result column, and **If it exists: Replace it · Add a new
  column**. Adding a new column keeps the earlier result and names the new one "Cluster (2)",
  and so on.
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
    Advanced Search, Browser), column widths, and the search bar with its text and
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
- Right-click a cell for Copy, Paste and Clear Cells. Structure cells offer **Edit in
  Sketcher…** and **Copy As ▸ Molfile / SMILES / InChI / InChIKey / Image** instead of Copy.
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

- **Sort Ascending / Descending.** **Then Sort…** adds a second (third…) sort key.
  **Clear Sort** returns rows to their original order.
- **Filter This Column.** Adds a filter card for the column (see [Filters](#7-filters)).
- **Rename…**, **Duplicate Column** (a copy, with its formatting, to the right),
  **Insert Column Left/Right…**, **Delete Column**. Delete acts on all selected columns when you
  right-click one of them.
- **Split Column…** and **Join Columns…**, each with a live preview.
- **Logarithmic (log₁₀)** (Number columns). Converts the values to log₁₀; choose it again to
  convert them back (10ˣ). The header shows "(log₁₀)". It needs every value above zero. One
  step you can undo.
- **Precision…** (Number columns). The number of decimals shown for every value. The stored
  values are not rounded, so sorting, filters and tools still use the full numbers. "As typed"
  goes back to each cell's own format.
- **Colour Column…** Colours cell backgrounds: **Low → high** or **Low → middle → high**
  gradients for numbers, or **One colour per value** for any column (the same value always gets
  the same colour). Range ends left empty follow the data. Choose **None** to remove the
  colours.
- Precision, Logarithmic and colours are saved in the document and can be undone.
- **Column Type ▸ Text / Number / Structure / Reaction.** Converting to Number warns you if
  some values are not numbers; those cells would be cleared. Converting to Structure or
  Reaction reads each cell's text; cells that do not read are kept and marked invalid.

From Tools ▸ Table:

| Command | What it does |
|---|---|
| Add Compound… | Opens the sketcher; draw (or paste) a structure, and **Add** appends it as a new row. |
| Add Rows… | Adds empty rows. |
| Add Column… | Adds a column. |
| Split Column… / Join Columns… | As above. |

### Rows

- **Delete Rows (Ctrl+Backspace)** deletes the selected rows. You can undo it.
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

- Click, Shift-click and Ctrl-click rows in the table.
  - **Select All (Ctrl+A)**, **Invert Selection (Ctrl+I)**, **Clear Selection
    (Ctrl+Shift+A)**.
- The selection is shared by the table, the plots and the tools:
  - selecting points in a plot selects the rows;
  - selecting rows highlights their points.
- Many tools can run on all rows, the shown rows or just the selected rows.

### Browser

Data ▸ Browser (Ctrl+Shift+B) steps through rows one at a time with a large drawing.

- It browses either the selected rows or every shown row, in the table's order.
- It keeps the table's current row in step with what it shows.
- Its tick box adds or removes the compound from the selection.

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
rows that have a matching cell.

- The **column menu** at the left of the bar chooses where to look: **All Columns** (every text
  and number column) or one column.
- **Substructure search.** Choose a structure column and type a substructure as SMILES or
  SMARTS (for example `c1ccccc1` or `[OH]a`). The table shows the rows that contain it, and
  the matching atoms are highlighted in the drawings. The search runs in the background
  shortly after you stop typing. A query that cannot be read turns the box red; hover over it
  to see why.
- The **Options** button sets **Match Case**, **Whole Cell** or **Regular Expression** (text
  searches only).
- If the chosen column is removed, the bar goes back to All Columns. A renamed column stays
  chosen. The column is saved with the document, with the search text and options.
- Escape clears the search. A second Escape, or the ✕ button, hides the bar.
- Hiding the bar always clears its search, so you never have rows hidden by a search you cannot
  see.

### Find and Replace (Ctrl+H)

A window that steps through matching cells, with **Replace** and **Replace All**. A replacement
is one step you can undo. **Replace All includes hidden rows** controls whether rows hidden by
filters or searches are changed too.

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
  - **Draw…** opens the sketcher. A plain drawing comes back as SMILES.
  - A drawing with query features comes back as SMARTS that keeps them: R-groups, `A`/`Q`/`X`/`M`,
    atom lists, and bond queries. An R-group or `*` matches any atom there.
  - **Draw…** again reopens your drawing, not the SMARTS.
  - The same applies in the Structure filter card and the structure conditions of Advanced
    Search.
- **Modes:**

  | Mode | Finds |
  |---|---|
  | Substructure | Rows that contain the query. |
  | Exact | Rows that are the query, with or without stereochemistry. |
  | Similarity | Rows at or above a similarity threshold, using a fingerprint and metric you pick. |

- **Results.** **Show** them alone (as a structure card in the Filters panel), **Select** them,
  or **Add to List**.
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
| Range | Number columns | A two-handle slider over a histogram of the values, plus min/max boxes. Choose whether blanks are kept. |
| Category | Text columns with few distinct values | Tick boxes, one per value, with counts. |
| Text | Other text | Contains / equals / starts / ends / regular expression, with **Not**. |
| Structure | Structure columns | Substructure, exact or similarity query; searches in the background after you stop typing. |

- Add cards with the **+** button in the panel's header or a column header's **Filter This
  Column**. The header's other buttons are **Filter Sets** and **Remove All Filters**; point at
  Remove All for how many rows are showing.
- Each card has an on/off tick box and a remove button. Point at a card's title for its column
  and its values (range, or how many categories are listed).
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
  - A job that takes more than 30 seconds and finishes while you are elsewhere shows a
    notification.

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
| Arrow round a circle | Turn | Drag to turn the cube, scroll to zoom, click a point to select it. |
| Dashed box | Box Select | Drag a rectangle to select the points inside it, as seen from where you look. |
| Loop | Lasso Select | Draw round points to select them, as seen from where you look. |
| Four corners | Reset View | Looks at the cube from the start position again. |
| Arrow into a tray | Export | Saves the view as a PNG image. |

- Shift adds to the selection; ⌘/Ctrl toggles. Hover over a point to see its card; point
  cards beside selected points follow them as you turn the cube.
- The plot is saved with the document, with the direction you were looking from.

### MedChem plots

The Data ▸ Plot ▸ MedChem menu has three plots. Each one first calculates any properties it needs that
are not already columns, then plots them.

| Plot | Axes | What it shows |
|---|---|---|
| BOILED-Egg | WLOGP against TPSA | The white predicts gut absorption; the yolk, brain penetration (Daina & Zoete 2016). |
| Golden Triangle | MW against logP | The region favoured for permeability and clearance (Johnson et al. 2009). |
| Lipophilic Efficiency (LipE) | Potency (pActivity) against logP | Lines of equal LipE = pActivity − logP. |

---

## 10. Chemistry: descriptors, pKa, standardizing and conformers

### Calculate Descriptors (Ctrl+Shift+D)

Pick a structure column and the properties you want. Each property becomes a number (or text)
column. **Recalculate descriptors that are already columns** controls whether existing results
are overwritten.

| Group | Descriptors |
|---|---|
| Physicochemical | MW, Exact Mass, cLogP, cLogD 7.4, LogS (ESOL), TPSA, TPSA (N,O,S,P), Molar Refractivity, Labute ASA, Fsp3 |
| Ionisation | pKa (most acidic), pKa (most basic), Charge @ pH 7.4, % Neutral @ pH 7.4 |
| Drug-likeness | QED, Ro5 Violations, Ro5 Pass, CNS MPO, BBB Score, AB-MPS |
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
- **3D shape descriptors** use a 3D structure when the cell has conformers.
- **Cancelling** keeps the rows already calculated.

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

- **Settings:** conformers per molecule, a minimum RMSD for two conformers to count as
  distinct, an energy window ("drop conformers more than … kcal/mol above the lowest"), the
  force field, and a random seed.
- **Output.** The conformers go into a new 3D structure column, lowest energy first, with
  columns **Conformers** (the count) and **Conformer Energy Range (kcal/mol)**.
  - Energies compare conformers of the same molecule only.
  - **Replace an existing 3D column of the same name** controls whether a previous result is
    overwritten.
- **Cancelling** keeps the molecules already done.

### Export conformers (SDF)

Writes every conformer as its own SD record, carrying the row's other columns. You can export
all rows or only the selected ones.

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

### Bin Values

Data ▸ Bin Values… groups a number column into bins.

- **Methods:** equal width, equal count, or custom edges.
- **Options:** a fixed width, from/to limits, and open-ended bins below the first edge and above
  the last.
- **Output:** a **range column** ("[5, 6)") and a **class column** with names you choose. Low /
  Medium / High is filled in for you.
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
| Chemistry | `MW ExactMW logP TPSA HBD HBA RotB Rings AromaticRings HeavyAtoms Fsp3 Charge QED Formula SMILES InChIKey Scaffold`, each taking a structure column, and `substructure(structure, "SMARTS")`, which counts the matches |

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

Data ▸ Similarity Search… scores every row (or the shown or selected rows)
against one of two references:

- a **structure** you type or paste; or
- the **selected compounds**. Each row gets its similarity to the nearest selected compound, and
  a column naming that compound. A selected row is never its own nearest neighbour.

Options:

- Sort the table, most similar first.
- Select rows with similarity of at least a threshold.

Cancelling keeps the rows already scored.

### Diverse Subset

Data ▸ Diverse Subset… picks *N* compounds that cover the chemistry, using MaxMin:
each pick is the compound farthest from all the earlier picks.

- **Pick from** all, shown or selected rows.
- You can start from the current selection, so the new picks complement compounds you have
  already chosen.
- The picks get a rank column, and can be selected.
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
- **If it exists: Add a new column** keeps an earlier clustering beside the new one.

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

Data ▸ SAR ▸ Hansch… regresses activity on substituent constants:

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

Data ▸ SAR ▸ Matched Pairs… finds pairs of compounds that differ by one small change
at one place. H → X changes are included.

- Set the largest change, in heavy atoms, and the minimum number of pairs a transformation needs
  to be listed.
- Transformations (A → B) are summarised over all their pairs, with the activity shift.

### Activity Cliffs (SALI)

Data ▸ SAR ▸ Activity Cliffs… lists similar compounds with very different activity.

    SALI = |Δactivity| / (1 − similarity)

- Only pairs at least the **minimum similarity** alike are considered.
- **Skip identical fingerprints** ignores salts and stereoisomers, which would otherwise divide
  by zero.

---

## 14. Chemical space

### Dimensionality Reduction

Data ▸ Chemical Space ▸ Dimensionality Reduction… places compounds on a 2D or 3D map by
fingerprints, number columns, or both. **Compare by: Both** uses the fingerprints and the
ticked number columns together, for rows that have all of them: the fingerprints are reduced to
their principal components, the columns are z-scored, and the two are weighted equally so
neither drowns out the other (distances are Euclidean). **Dimensions: 2D** writes two columns and opens a scatter plot;
**3D** writes three and opens a 3D Scatter (a 3D PCA map is always PC1, PC2 and PC3).

| Method | Options | More Options (folded) |
|---|---|---|
| PCA | **Components**: PC1 and PC2 (the default), PC1 and PC3, PC2 and PC3, or PC3 and PC4 | — |
| UMAP | Neighbours, Min. distance | **Distance** (Automatic, Jaccard, Dice, Rogers-Tanimoto, Euclidean, Manhattan, Cosine, Correlation), **Spread**, **Epochs**, **Fit on up to** (rows; larger sets are fitted on a sample of this size and the rest placed onto that map; default 50,000), **Seed** |
| PaCMAP | Neighbours (default 10) | **Mid-near pairs** (per neighbour pair; default 0.5), **Further pairs** (default 2), **PaCMAP distance** (Euclidean, Manhattan, Angular), **PCA dimensions** (default 50), **Seed** |
| t-SNE | Perplexity | **Learning rate** (Automatic or a number), **Iterations**, **Exaggeration** (early exaggeration), **PCA dimensions** (t-SNE runs on this many principal components; default 50), **Start from** (PCA or Random), **Seed** |

- **Distance: Automatic** is Jaccard (Tanimoto) on fingerprints and Euclidean on number
  columns. Jaccard, Dice and Rogers-Tanimoto need fingerprints.
- **Min. distance** must not be larger than **Spread**.
- **Seed**: the same seed gives the same map every run (UMAP then runs on one processor core).
  Set it to **Random** for a faster map that differs each time.
- PCA's note in the status bar gives the share of the variance on each chosen component.
- **PaCMAP** (Wang et al., 2021) balances neighbour, mid-near and further pairs, so it keeps
  local groups like UMAP while distances between groups mean more. It runs on the PCA
  dimensions of the features and is fast on large sets.
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
| Browser | Ctrl+Shift+B |
| Settings | Ctrl+, |
