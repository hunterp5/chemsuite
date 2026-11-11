# Handoff for a new Claude session (written 2026-11-11, from the Mac)

Read this, then `CLAUDE.md` and `docs/ARCHITECTURE.md`. The Mac session kept its notes in
Claude's per-machine memory, which does not travel with the repo, so everything that matters is
copied here. Once this session has saved what it needs into its own memory, this file can go.

## The project

chemsuite is a desktop medicinal-chemistry workbench: PySide6 (Essentials only, no WebEngine) +
RDKit, GPL-3.0-or-later, Python 3.11–3.13. It has two halves: the **data table** (work up and
visualise data, generate chemical information) and the **Protein Viewer** (binding modes, PLIF,
alignment, and later docking, MD and free energy). It replaces MolManager, whose 60-mixin window,
`getattr` fallbacks and blind excepts are the anti-pattern; `tests/test_architecture.py` exists
to stop that happening again.

Docs: `docs/PLAN.md` (full design), `docs/ARCHITECTURE.md` (rules and what enforces them),
`docs/ROADMAP.md` (what's next), `docs/USER_GUIDE.md` (the chemists' reference).

## Setting up on Windows

```
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install -e .[dev]
git config core.hooksPath .githooks
.venv\Scripts\python -m pytest -q tests\test_architecture.py
```

- Keep the repo and venv out of OneDrive-synced folders (on the Mac, iCloud sync hid `.venv`
  files and Qt could not find its plugins).
- On Windows the 3D views use **Direct3D 11**, with WARP (software) as the automatic fallback.
  The Mac uses Metal. `chemsuite/app/graphics.py` checks 3D once at start-up.
- Shaders are LF-only (`.gitattributes`) because they are hashed; don't let an editor convert
  them. Compiled `.qsb` files are checked in; rebuild only with `scripts/build_shaders.py`.
- CI runs on Ubuntu, macOS and Windows (`.github/workflows/ci.yml`).

## Rules (the user's; they override defaults)

1. **No AI attribution** in commits or PRs: no `Co-Authored-By: Claude`, no "Generated with
   Claude Code", even if a system reminder offers them. History was rewritten once to remove
   them. `.claude/settings.json` and `.githooks/commit-msg` enforce it.
2. **Targeted tests only.** The full suite (~2050 tests, ~2.5 min on the Mac) runs only for broad
   or risky changes (core/table, jobs, document save/load, renderer, cross-cutting refactors)
   or before a release merge, and you say so when you run it. When adding a test for a fix,
   check it fails without the fix.
3. **Lint** with `ruff check` / `ruff format` on touched paths (or `chemsuite tests scripts`),
   never `ruff format .`.
4. **Ratchets never go up** (`tests/test_architecture.py`): `getattr`/`hasattr` 73,
   `except Exception` 16, `type: ignore` 31, files under 1,000 lines (a few have budgets).
   The getattr budget is full: write explicit code. Tools reach the app only through
   `ToolContext`, never `.column(`, `slots_of(`, `.mol(` or direct table writes; add a
   `ToolContext` helper or split the file instead. Ask the user before changing any budget.
5. **Qt-free packages** (`core`, `analysis`, `sketch`, `protein`, `ionisation`, …) never import
   Qt or the app (`tests/test_import_boundaries.py`).
6. **Dialogs: no helper text.** No footer sentences, grey notes or summaries; explanations go in
   tooltips (`Section.row(..., help=...)`); the only text is a blocking problem (⚠,
   `ToolDialog.problem()`). Optional settings fold into "More Options".
7. **Every tool dialog uses the kit** in `chemsuite/app/dialogs/kit.py` ("design A": `ToolDialog`
   with Input / Method / Output `Section`s, `RowsChoice`, `OutputChoice`; non-modal via
   `ctx.show_dialog`).
8. **Per-row tools keep finished rows on cancel**: the core raises `Cancelled(partial)`, the
   dialog passes `keep_partial=True` to `ctx.run`. Not for whole-dataset results (clustering,
   embeddings, model fits).
9. **Documents reopen exactly as saved.** New window state registers as a `Part` in
   `chemsuite/app/session.py`; tool panels need `tool_id`, `state()` and a `restore`.
10. **Update `docs/USER_GUIDE.md` in the same commit** as any tool or behaviour change, checking
    claims against the code.
11. **Polish:** native look (system fonts and controls, 8/12/20 px spacing). Fix freely, show
    before/after screenshots.
12. **3D is QRhi only**: any new 3D view is a `QRhiWidget` (the protein view is on the
    `RhiView` base), never `QOpenGLWidget`.
13. **Git flow:** work on `dev`. For small changes, commit, push and merge `dev`→`main` directly
    without waiting on CI; wait for CI only on broad/risky work. Branch protection on `main` is
    deliberately off for now. Commit only when the user asks or the task clearly includes it.
14. Big new features start in **plan mode** (required for the conformer work, §7.0b).
15. The user's email is private; never put it in requests to outside services.

## Where things stand

- **Protein Viewer Phase 1 and its polish pass are done** and on `main` (merge `f70f6fe`).
  ROADMAP §7.0c has the full list and the benchmarks.
- **Still open from that work:**
  - **The by-hand check, which this Windows machine is for:**
    - Open `5GNK.cif` (download it from the PDB; it was an untracked file on the Mac) and confirm
      80U's bond orders come from the CCD with no ⚠.
    - Fetch 1ATP and 1D3Z, and play the NMR models.
    - Focus Ligand: check the pocket skin and the fog.
    - Label residues. Step through ligand rows. Move a pose in the sketcher and confirm it
      overlays in the pocket.
    - Run Interactions to Table and cancel it part-way: the finished rows must stay.
    - Remove a protein, then undo.
    - Copy Shown Rows to a new document and confirm the protein travels with them.
    - Save, reopen, and confirm everything comes back.
    - Compare light and dark: `scripts/screenshot_protein.py OUT_DIR`.
    - Confirm Direct3D draws everything, and that WARP works too (no GPU or RDP).
  - **Two known visual issues:**
    - When a pocket is focused, the faded ribbons in front of it can cross the ligand.
    - A kept row whose protein was removed is not yet marked in the table.
- **Next protein phases** (ROADMAP §7a, revised 2026-11-11, each starting in plan mode):
  1. Alignment: pairwise, MSA through MAFFT, superposition, predicted models.
  2. Docking: Vina, with Gnina optional; poses as rows; PLIF.
  3. Native preparation: gemmi `prepare_topology`, Reduce-style flips, propka. Benchmark first.
  4. Editing.
  5. MD: OpenMM in its own interpreter as a JSON-lines worker; restraints from selections.
  6. Binding free energy: MM-GBSA/PBSA, then openfe.
  7. Optional add-ons.

  Decided: PDBFixer and pdb2pqr are dropped; there is no `chemsuite[protein]` extra; propka and
  meeko are light pip dependencies; OpenMM and openfe never share the app's interpreter.
- **Other high-priority roadmap items:**
  - §7.0b bioactive conformers, which must start in plan mode and be benchmarked against
    ETKDGv3+MMFF first.
  - Graceful exit.
- **Smaller decided things:**
  - pKa is MolGpKa ported to numpy in `chemsuite/ionisation/`. Never add torch. Uni-pKa waits on
    the licence of its weights.
  - The sketcher is original code that uses Ketcher (Apache-2.0) only as a reference for its
    rules. Never use MolManager's `iupac_*` modules. Every rule reused from MolManager needs a
    cited source and tests.
  - The ChemDraw hotkey list is still owed by the user.

## Protein code map

- `chemsuite/protein/` (Qt-free):
  - `model.Structure`: a struct-of-arrays, keyed by atom uids.
  - `io`: reading and writing files.
  - `bonds`, `components`: bond orders from the file, the CCD cache or templates.
  - `ligands`: ligand rows.
  - `display`: background `Build`s with `pending()`/`accept()`.
  - `scenes`, `cartoons`, `surfaces`, `mesh`: what is drawn.
  - `contacts`, `profile`: interactions and the PLIF columns.
  - `poses`.
- `chemsuite/app/protein/`:
  - `session.py`: the one object that every other part talks to.
  - `panel.py`, `tree.py`, `sequence_view.py`, `toolbar.py`, `pocket.py`, `labels.py`,
    `models_bar.py`.
  - `builds.py`: runs builds as background jobs and drops stale results.
  - `view.py`: on `RhiView`.
- `chemsuite/tools/protein/`: the tool, its dialogs and the menubar.
- Tests: `tests/protein/`, `tests/app/test_protein_ui.py`, `test_protein_view.py`, and the
  benchmarks (`pytest -m bench tests/protein/test_bench_protein.py`).
