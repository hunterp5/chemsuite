# chemsuite — rules for Claude

## Testing: targeted, not the full suite

The full suite (`pytest`, ~1500 tests) is slow. Do not run it for minor changes.

- Run only the tests that cover what changed: the test file(s) for the touched module, plus
  any test that exercises the changed behaviour (`grep -rl <name> tests` finds them). Use
  `-k` to narrow further while iterating.
- Qt-free code (`chemsuite/core`, `analysis`, `sketch`): its own `tests/<package>/…` file.
  App code: the matching `tests/app/test_*_ui.py` / `test_*.py` file.
- Changed imports between packages → also run `tests/test_import_boundaries.py`.
- Run the full suite only for broad or risky changes (shared infrastructure such as
  `core/table`, `jobs`, document save/load, cross-cutting refactors), or right before a
  release merge to `main` — and say so when you do.
- When you add a test for a fix, check it fails without the fix.
- Lint with `ruff check` / `ruff format` on the touched paths (or `chemsuite tests scripts`),
  never `ruff format .`.

## Commits

No AI attribution: no `Co-Authored-By: Claude`, no "Generated with Claude Code" lines in commits
or pull requests. `.githooks/commit-msg` enforces it (`git config core.hooksPath .githooks`).

## Dialogs and panels: no helper text

No explanatory text in dialogs, popups or settings menus: no footer sentence describing what
will happen, no grey heads-up notes ("may take minutes"), no summaries under the fields, no
"about" paragraphs. Explanations go in tooltips (`Section.row(..., help=...)`). The only text
a dialog shows is a problem that blocks running (`ToolDialog.problem()`, marked ⚠). Optional
settings go in a folded section (`self.section("More Options", collapsed=True)`).
