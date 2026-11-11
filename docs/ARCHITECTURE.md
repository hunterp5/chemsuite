# chemsuite architecture: the rules and what enforces them

MolManager, the first attempt, did not fail in one step. It had a window built from about 60
mixins, 271 places where dialogs reached into the window, about 2,000 `getattr` fallbacks and
1,088 blind `except` blocks. Each one arrived as a reasonable-looking change. chemsuite keeps
its shape with **rules that fail a check** and **ratchets that only go down**. A rule written
here and nowhere else would drift too, so every rule below names the test or lint that holds
it.

The full design is in [PLAN.md](PLAN.md); this page is the short list to keep in mind while
changing code.

## Rules

| Rule | Why | Enforced by |
|---|---|---|
| `core`, `analysis`, `molview`, `sketch` and the other domain packages never import Qt or the app. | They stay testable, scriptable and fast to import. | `tests/test_import_boundaries.py` |
| Tools use the app only through `ToolContext` (`chemsuite.app.tools`) and a short list of shared widgets. | Tools never depend on the window, so adding one never means editing it. | `test_tools_use_the_app_only_through_its_tool_api` (`TOOL_APP_IMPORTS`) |
| Tools never reach into the window or the context's private members. | MolManager's 271 `app._x` reaches. | `test_tools_never_reach_into_the_window` |
| Tools read the table through `ToolContext` (`structure_inputs`, `table.read`) and write results only through `ctx.apply_results`. | One place handles undo, name clashes, deleted rows, and speed at 200k rows. | `test_tools_use_the_table_through_the_context` (ratchet) |
| A panel that opens from a tool saves and reopens with the document. | Documents reopen exactly as saved. | `tests/app/test_session.py` |
| Work runs as jobs (`ctx.run`): the queue lane for heavy work, the quick lane for interactive work; per-row tools keep finished rows on cancel; work the user did not ask for is `background=True`. | The window never blocks; nothing is lost; closing never asks about work nobody started. | `tests/app/test_jobs_ui.py` |
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
- **3D draws through QRhi** (Qt's GPU layer). The Structure Viewer does today. Moving the 3D
  plots off pyqtgraph.opengl onto QRhi is a high priority (ROADMAP.md §7): it frees every view
  to use Direct3D on Windows and Metal on macOS instead of OpenGL everywhere.
- **No QtWebEngine**: PySide6-Essentials only, to keep the install small. (PLAN.md §1)
