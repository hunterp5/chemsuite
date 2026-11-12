# chemsuite architecture: the rules and what enforces them

MolManager, the first attempt, did not fail in one step. It had a window built from about 60
mixins, 271 places where dialogs reached into the window, about 2,000 `getattr` fallbacks and
1,088 blind `except` blocks. Each one arrived as a reasonable-looking change. chemsuite keeps
its shape with **rules that fail a check** and **ratchets that only go down**. A rule written
here and nowhere else would drift too, so every rule below names the test or lint that holds
it.

The full design is in [PLAN.md](PLAN.md); this page is the short list to keep in mind while
changing code.

## Layout

- `chemsuite/core`: the table, documents, jobs and settings, beneath everything else.
- `chemsuite/io` (reading and writing compound files), `analysis`, `sketch`, `ionisation`,
  `molview` and the other domain packages: Qt-free, above `core`.
- `chemsuite/app`: the Qt application.
  - `shell`: what every window uses (start-up, jobs and `JobHost`, errors, notifications,
    settings, themes, commands, glyphs, `ToolWindow`).
  - `widgets`: shared widgets (the dialog kit, `ToolPanel`, flow layout, panes).
  - `render3d`: the `RhiView` base and the shared shaders.
  - `document`: the table window (main window, session, `ToolContext`, panels and docks).
  - `table`: the table view and model, column actions, search and filters.
  - `plots`, `viewer`, `sketcher`: the plot panels, the Structure Viewer, the sketcher.
- `chemsuite/protein`: the Protein Workbench, a subproject. `structure`, `chem` and `look`
  are Qt-free; `ui` is its Qt workbench.
- `chemsuite/tools`: one package per menu feature. `tools/protein` is the workbench's
  bridge to the table.
- `tests` mirrors the code. The GUI test folders (`tests/app`, `tests/tools`,
  `tests/protein/ui`) share the Qt setup in `tests/gui.py`.

## Rules

| Rule | Why | Enforced by |
|---|---|---|
| `core`, `analysis`, `molview`, `sketch`, `protein` (outside `protein/ui`) and the other domain packages never import Qt or the app. | They stay testable, scriptable and fast to import. | `tests/test_import_boundaries.py` |
| `core` imports no domain package (file formats live above it, in `chemsuite/io`, which may use `sketch` to lay out reactions it writes), and the app never imports a tool: it reaches one through its `Tool` declaration (`open_at` for Data ▸ Browser). | Layers only point down; tools plug in. | `test_core_imports_no_domain_package`, `test_the_app_never_imports_a_tool` |
| Tools use the app only through `ToolContext` (`chemsuite.app.document.tools`) and a short list of shared widgets. | Tools never depend on the window, so adding one never means editing it. | `test_tools_use_the_app_only_through_its_tool_api` (`TOOL_APP_IMPORTS`) |
| The Protein Workbench's Qt layer (`protein/ui`) never imports the table window, the table view or a tool; it reaches the document's current row through a `LinkedRow` (`protein/ui/host.py`) its host hands it. Only the bridge, `tools/protein`, imports it. | The workbench grows on its own and will run without a table. | `test_the_workbench_ui_reaches_the_table_only_through_the_bridge` (`panel.py` and `session.py` still read the `Table`: a ratchet until the workspace engine), `test_only_the_bridge_imports_the_workbench_ui` |
| Tools never reach into the window or the context's private members. | MolManager's 271 `app._x` reaches. | `test_tools_never_reach_into_the_window` |
| Tools read the table through `ToolContext` (`structure_inputs`, `table.read`) and write results only through `ctx.apply_results`. | One place handles undo, name clashes, deleted rows, and speed at 200k rows. | `test_tools_use_the_table_through_the_context` (ratchet) |
| A panel that opens from a tool saves and reopens with the document. | Documents reopen exactly as saved. | `tests/app/document/test_session.py` |
| A tool's panel is a `ToolPanel` (`app/widgets/panel.py`); docks, tool windows and the saved session ask it only what that class declares. | No duck-typed protocol: a host never asks whether a panel has something. | `test_every_tool_panel_is_saved_with_the_document`, the `getattr` ratchet |
| A document window's docks and tool windows belong to its `DocumentPanels` (`app/document/panels.py`); nothing else keeps or changes those lists. | MolManager's helpers reached into the window's private state. | review; `window._x` reaches fail `test_tools_never_reach_into_the_window` for tools |
| Every window runs its work through a `JobHost` (`app/shell/job_host.py`): progress, failures, notices, kept partial results, results dropped once the document is gone. | One way to run work, whatever the window. | `tests/app/shell/test_jobs_ui.py` |
| Work runs as jobs (`ctx.run`): the queue lane for heavy work, the quick lane for interactive work; per-row tools keep finished rows on cancel; work the user did not ask for is `background=True`. | The window never blocks; nothing is lost; closing never asks about work nobody started. | `tests/app/shell/test_jobs_ui.py` |
| No blind `except Exception` except at a top-level handler, with a reason. | Errors surface instead of hiding. | ruff `BLE001`, `S110` and the ratchet below |
| No function above complexity 15 without a stated reason. | Branchy code is where bugs hide. | ruff `C90` (exemptions are `# noqa: C901 - reason`; `RUF100` removes stale ones) |
| Dialogs show no helper text; explanations go in tooltips. | A clean, native look. | `CLAUDE.md`, review |

## Ratchets

`tests/test_architecture.py` holds today's counts as budgets. A change may lower a count, and
then lowers the budget in the same commit so it stays lowered. Raising a budget means editing
that file, which a reviewer sees and has to agree to.

- Column internals and direct table writes in tools, per file (`COLUMN_INTERNALS_IN_TOOLS`,
  `DIRECT_WRITES_IN_TOOLS`).
- `getattr`/`hasattr`, `except Exception` and `type: ignore` across the package (`BUDGETS`).
- File length: new files stay under 1,000 lines; the few longer ones have budgets
  (`LONG_FILES`). When a file hits its budget, split it, as `column_actions.py` and
  `list_actions.py` were split out of the main window, rather than raise the number.

## When a rule gets in the way

Change the rule on purpose rather than work around it: edit the allowlist or budget, say why
in the commit message, and update this page if the rule itself changed. A helper on
`ToolContext` (as `structure_inputs` and `workers_for` replaced five copies of the same code)
is usually better than an exemption.

## Decisions

- **Structures are a column kind** (`StructureColumn`), not a separate store: one undo, copy and
  save path for every column. (PLAN.md §2.2)
- **3D draws through QRhi** (Qt's GPU layer): the Structure Viewer and the 3D plots, in the
  system's own API (Direct3D 11, Metal, OpenGL on Linux). Every GPU widget in a window shares
  one API, so a new 3D view must be a `QRhiWidget` too, never a `QOpenGLWidget`. Shaders are
  compiled to `.qsb` with `scripts/build_shaders.py` and checked in. Each package keeps its
  own in a `shaders` folder with its own manifest (`app/render3d`, `protein/ui`); a view names
  its package in `RhiView.shaders`.
- **3D is checked once at start-up** (`chemsuite/app/shell/graphics.py`): anything that draws 3D
  reads `graphics.problem()` to grey itself out with the reason instead of failing when
  opened, and a view that fails when drawing sets it.
- **One container for every file** (`core/container.py`): SQLite with a `meta` table and an
  optional `assets` table, saved atomically, opened only when its version is not newer than
  the reader's. The `.csd` adds its column chunks (`core/document.py`).
- **No QtWebEngine**: PySide6-Essentials only, to keep the install small. (PLAN.md §1)
- **Assets beside the rows.** A document holds immutable assets (`core/assets.py`), such
  as proteins, in an `AssetStore`, changed only by `PutAsset` ops. Inside a table document the
  `Table` holds the store and runs its ops, so they share the table's `History`
  (`core/history.py`) and one undo stack; alone (the protein workspace without a table) the
  store keeps its own history. Readers see either through the `AssetHost` protocol. Rows name an asset in
  a link column (`LINK_COLUMNS`: protein → `Protein`); copying rows to a new document takes the
  assets they name. Adding or removing a protein and its ligand rows is one undo step
  (`add_rows`/`remove_rows` with `extra_ops`).
- **The Protein Workbench is a subproject** (`chemsuite/protein`), so it can grow on its own and
  later open without a table: Qt-free `structure` (the arrays, files, bonds, selections),
  `chem` (ligands, poses, interactions) and `look` (displays, scenes, meshes), and the Qt
  `ui`. Only the bridge, `tools/protein`, knows both the table and the workbench: the
  ligand-row columns and `ligand_rows` (`rows.py`), following the current row (`follow.py`,
  the workbench's `LinkedRow`), and Ligands and Interactions to Table (`to_table.py`).
- **Proteins are Qt-free arrays** (`chemsuite/protein/structure`): an immutable struct-of-arrays
  `Structure` keyed by atom uids (selections, labels, measurements and focus survive edits),
  a `Display` per protein holding how it looks, and a composed `Scene` the view uploads as a
  geometry buffer and a separate look buffer, so a selection change re-uploads colours only.
  The app side (`chemsuite/protein/ui`) is a session object that the tree, sequence, selection
  bar and view read and change, never each other.
- **One `RhiView` base for 3D views** (`app/render3d/rhi_view.py`): camera, picking, slab and
  fog. Meshes fall back to 16-bit index chunks where the GPU lacks 32-bit indices.
- **A workspace of its own is a window, not a dock.** `ToolContext.add_window` hosts a tool's
  widget in a top-level `ToolWindow` (`app/shell/windows.py`, owned by the window's
  `DocumentPanels` in `app/document/panels.py`) with the tool's own menu bar
  and no toolbar: the Protein Viewer. It belongs to the document (closes with it, saved and
  reopened with it, geometry included) and its shortcuts are its own. Docks
  (`add_panel`) stay for panels that sit beside the table.
- **Slow protein parts build in the background, never stale.** Surfaces, interactions, pocket
  skins and large cartoons are `Build`s from `Display.pending()`, run as quick background jobs
  (`protein/ui/builds.py`); a result is kept only if the display is still shown, still holds
  the structure the work was handed, and the part's version has not moved on. The last result
  stays drawn until the new one lands.
- **Every download goes through `core/download.py`**: HTTPS only, audited, refused when IT
  policy turns the network off, with an injectable opener so tests never touch the network.
  Files that can be fetched again (the CCD components) live in `ToolContext.cache_folder`.
- **Heavy simulation runs in its own interpreter** (planned: OpenMM, openfe), driven as a
  JSON-lines worker like `retro/runner.py`, so its pins never meet RDKit's and a native crash
  never takes the app down.
