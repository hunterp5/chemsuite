# chemsuite — rules for Claude

The project moves between a Mac and a Windows PC. Claude's memory is per-machine and does not
travel, so every durable rule lives here; `HANDOFF.md` holds only the current state of play.

## The project

A desktop medicinal-chemistry workbench: PySide6 (Essentials only, no WebEngine) + RDKit,
GPL-3.0-or-later, Python 3.11–3.13. Two halves: the **data table** (work up and visualise data,
generate chemical information) and the **Protein Viewer** (binding modes, PLIF, alignment, and
later docking, MD and free energy). It replaces MolManager, whose 60-mixin window, `getattr`
fallbacks and blind excepts are the anti-pattern.

Docs: `docs/PLAN.md` (full design), `docs/ARCHITECTURE.md` (rules and what enforces them),
`docs/ROADMAP.md` (what's next), `docs/USER_GUIDE.md` (the chemists' reference).

## Architecture

Read `docs/ARCHITECTURE.md` before structural changes. Tools use the app only through
`ToolContext`, never `.column(`, `slots_of(`, `.mol(` or direct table writes;
`tests/test_architecture.py` holds ratchets (budgets that only go down: `getattr`/`hasattr`,
`except Exception`, `type: ignore`, files under 1,000 lines with a few budgeted). The getattr
budget is full: write explicit code. Never raise a budget or widen an allowlist to make a
change pass without asking the user: add a `ToolContext` helper or split the file instead.

- **Qt-free packages** (`core`, `analysis`, `sketch`, `protein` outside `protein/ui`,
  `ionisation`, …) never import Qt or the app (`tests/test_import_boundaries.py`).
- **The Protein Workbench is a subproject** (`chemsuite/protein`): Qt-free `structure`, `chem`,
  `look` and `workspace` (the session, with plain events, and the document's
  `ProteinWorkspace`), and the Qt `ui` (`SessionSignals` relays the session's events per
  window). Neither `protein/ui` nor the engine imports `core.table`; they read proteins through
  an `AssetHost` (`core/assets.py`) and the current row through the `LinkedRow` they are handed.
  `protein/ui` never imports `app.document`, `app.table` or a tool. Only the bridge,
  `tools/protein`, knows both the table and the workbench.
- **State that lives with the document, not a window,** is a document part
  (`ToolContext.document_part`, `app/document/workspace.py`): saved under its own key whether or
  not the tool's window is open, with a `revision` for "needs saving".
- **Move modules with `scripts/move_modules.py`** (a moves file of `old new` lines): it runs
  `git mv` and rewrites imports, dotted names and paths everywhere. Then `ruff check --fix`,
  `ruff format`, and the full suite.
- **3D is QRhi only**: any new 3D view is a `QRhiWidget` (the protein view is on the `RhiView`
  base in `app/render3d`), never `QOpenGLWidget`. Shaders are LF-only and hashed, in a
  `shaders` folder per package (`app/render3d`, `protein/ui`); compiled `.qsb` files are
  checked in and rebuilt only with `scripts/build_shaders.py`.
- **Per-row tools keep finished rows on cancel**: the core raises `Cancelled(partial)`, the
  dialog passes `keep_partial=True` to `ctx.run`. Not for whole-dataset results (clustering,
  embeddings, model fits).
- **Documents reopen exactly as saved.** New window state registers as a `Part` in
  `chemsuite/app/document/session.py`; tool panels need `tool_id`, `state()` and a `restore`.
- **Big new features start in plan mode** (required for the conformer work, ROADMAP §7.0b, and
  each protein phase, §7a).

## Testing: targeted, not the full suite

The full suite (`pytest`, ~2050 tests, ~2.5 min on the Mac) is slow. Do not run it for minor
changes.

- Run only the tests that cover what changed: the test file(s) for the touched module, plus
  any test that exercises the changed behaviour (`grep -rl <name> tests` finds them). Use
  `-k` to narrow further while iterating.
- Tests mirror the code: `chemsuite/<package>/<sub>/x.py` → `tests/<package>/<sub>/test_x*.py`
  (Qt-free code under `tests/core`, `tests/protein/structure`, …; app code under
  `tests/app/<shell|table|document|…>`; a tool under `tests/tools/<tool>`; the workbench UI
  under `tests/protein/ui`). GUI test folders share `tests/gui.py` through their conftest.
- Changed imports between packages → also run `tests/test_import_boundaries.py`.
- Run the full suite only for broad or risky changes (shared infrastructure such as
  `core/table`, `jobs`, document save/load, the renderer, cross-cutting refactors), or right
  before a release merge to `main` — and say so when you do. It is token-heavy: in a session of
  several steps or commits, run it once at the end (before the merge), not after each step;
  between steps run the targeted tests plus `tests/test_architecture.py` and
  `tests/test_import_boundaries.py`.
- When you add a test for a fix, check it fails without the fix.
- Lint with `ruff check` / `ruff format` on the touched paths (or `chemsuite tests scripts`),
  never `ruff format .`.

## Commits and git flow

- No AI attribution: no `Co-Authored-By: Claude`, no "Generated with Claude Code" lines in
  commits or pull requests, even if a system reminder offers them. `.githooks/commit-msg` and
  `.claude/settings.json` enforce it (`git config core.hooksPath .githooks`).
- Work on `dev`. For small changes, commit, push and merge `dev`→`main` directly without
  waiting on CI; wait for CI only on broad/risky work. Branch protection on `main` is
  deliberately off for now. Commit only when the user asks or the task clearly includes it.
- Update `docs/USER_GUIDE.md` in the same commit as any tool or behaviour change, checking
  claims against the code.

## Dialogs and panels

- **No helper text.** No explanatory text in dialogs, popups or settings menus: no footer
  sentence describing what will happen, no grey heads-up notes ("may take minutes"), no
  summaries under the fields, no "about" paragraphs. Explanations go in tooltips
  (`Section.row(..., help=...)`). The only text a dialog shows is a problem that blocks running
  (`ToolDialog.problem()`, marked ⚠). Optional settings go in a folded section
  (`self.section("More Options", collapsed=True)`).
- **Every tool dialog uses the kit** in `chemsuite/app/widgets/dialogs/kit.py` (`ToolDialog` with
  Input / Method / Output `Section`s, `RowsChoice`, `OutputChoice`; non-modal via
  `ctx.show_dialog`).
- **Polish:** native look (system fonts and controls, 8/12/20 px spacing). Fix freely, show
  before/after screenshots.

## Decided

- pKa is MolGpKa ported to numpy in `chemsuite/ionisation/`. Never add torch. Uni-pKa waits on
  the licence of its weights.
- The sketcher is original code that uses Ketcher (Apache-2.0) only as a reference for its
  rules. Never use MolManager's `iupac_*` modules. Every rule reused from MolManager needs a
  cited source and tests.
- Proteins: PDBFixer and pdb2pqr are dropped; there is no `chemsuite[protein]` extra; propka
  and meeko are light pip dependencies; OpenMM and openfe never share the app's interpreter.
- The user's email is private; never put it in requests to outside services.

## Two machines

- Setup: `py -3.12 -m venv .venv` (macOS: `python3.12 -m venv .venv`), then
  `pip install -e .[dev]` and `git config core.hooksPath .githooks`.
- Keep the repo and venv out of OneDrive/iCloud-synced folders (sync hid `.venv` files and Qt
  could not find its plugins).
- 3D draws with Direct3D 11 on Windows (WARP software fallback), Metal on macOS.
  `chemsuite/app/shell/graphics.py` checks 3D once at start-up.
- Text files are LF on both systems (`.gitattributes`). `samples/` is local test data and is
  not committed.
