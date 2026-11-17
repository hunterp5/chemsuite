# A unified Similarity Search (plan approved 2026-11-16)

## Context

Similarity Search… today scores rows against one structure (or the nearest selected row) with
one 2D fingerprint. Superposition now gives chemsuite shape, shape + features and field scoring
(`chemsuite/superpose`). ROADMAP §6 names "3D similarity on the same engine" as next. The goal
is one tool that finds similar chemical matter in several ways and in two places:

- **Where it searches:** the open table, or a library file (SDF or SMILES, such as a vendor
  catalog), with the hits opened as a new document.
- **Methods:** 2D fingerprints; substructure (MCS and scaffold); pharmacophore (2D Gobbi now);
  3D shape, features and field, aligned with the superpose engine.
- **References:** one or more; scores fused per row; several methods combined into one
  consensus score.
- **Old dialog:** it is replaced. The new tool keeps the Similarity Search… title, and its
  fingerprint method does everything the old one did, including "nearest selected row".
  Diverse Subset is unchanged.

(Decided with the user: table + files, replace the old dialog, all four method families.)

## Architecture

There is a new **Qt-free engine package, `chemsuite/similarity/`**, added to `QT_FREE` in
`tests/test_import_boundaries.py`. It follows `chemsuite/superpose`'s shape:
model / methods / run. The tool UI stays in `chemsuite/tools/similarity/`. The tool reaches the
app only through `ToolContext`: `structure_inputs`, `run(keep_partial=True)`, `apply_results`,
`open_table`, `workers_for`, `cache_folder`.

```
chemsuite/similarity/
  model.py        Reference(mol, name), Method specs (frozen dataclasses), Search settings,
                  score column names, LOWER_IS_BETTER; validation in __post_init__ (as superpose)
  fingerprint.py  wraps analysis.fingerprints; Tversky; several references
  substructure.py MCS Tanimoto (rdFMCS, timeout) and scaffold match (descriptors/scaffold.py)
  pharmacophore.py Gobbi 2D pharmacophore fingerprint (RDKit Pharm2D), folded to fixed bits
  shape3d.py      delegates to superpose (shape, shape + features, field), best over references
  fusion.py       per-row fusion over references (max | mean); consensus over methods
                  (mean of scores | mean rank | min-max normalized mean)
  run.py          score_cells (table rows) and screen_file (streamed library), both chunked
                  with ordered_map, CancelToken, Cancelled(partial)
```

### Reuse, not rewrite

- **2D:** `analysis/fingerprints.py` already has `KINDS`, `compute_fingerprints`,
  `FingerprintCache`/`cache_for` (table rows), `similarity`, `max_similarity`/`_best`
  (nearest of many references) and `from_counts`/`at_least`.
  - Add **Tversky (α, β)** to `from_counts` and `at_least`.
  - Add a **`pharm2d`** `FingerprintKind` (Gobbi, folded to 4096 bits), so pharmacophore 2D
    rides on the same packed, BLAS-tiled path, with no new similarity code.
- **3D:** `superpose/run.py` has `prepare_reference`, `probe_of` (generates conformers for 2D
  rows), `superpose_one`, the `ALIGNERS` and `Fit.scores`.
  - Add one engine function, `best_over_references(mol, refs, settings) -> (ref_index, Fit)`.
    It prepares each reference once per chunk and generates the probe's conformers once for
    all references. The Superpose dialog is unchanged.
  - Keeping poses: lift the "write superposed column + reference asset" code out of
    `tools/superpose/dialog.py` into a small shared helper both tools call (`blob_with_reference`,
    `PutAsset`), so the Browser shows the reference as it does today.
- **Scaffolds:** `chemsuite/descriptors/scaffold.py` (Bemis–Murcko and generic).
- **Files:** `io/load.read_batches` (streamed), `io/load.TableBuilder` to build the hits
  table, `ctx.open_table(table, message)` to open it as a new document.
- **UI kit:** `ToolDialog`, `Section.row(help=…)`, `RowsChoice`, `Segmented`, `inline`,
  `follow_rows`. `FeatureChooser`'s fingerprint combo (`tools/similarity/common.py`) and the
  reference field (SMILES + draw + Current Row) from today's `search.py`.

### Two-stage search (keeps 3D and MCS affordable)

Fingerprint scores are cheap; MCS and 3D cost seconds per molecule. When an expensive method is
on, stage 1 scores every candidate with the chosen fingerprint, and stage 2 runs MCS/3D only on
the best **N** (More Options ▸ Prefilter, default 500; 0 means everything). Rows outside the
prefilter get blank expensive scores. On a table this is optional. On a file it is always on
when 3D is chosen.

### File screening

`screen_file(path, search, top_k, threshold)` streams the file in batches. It fingerprints
each batch in the process pool (`fingerprint_chunk`) and keeps a bounded min-heap of the top K
(default 1,000) rows that pass the threshold, keeping only those rows' raw records. When the
file is done, stage 2 runs on the survivors. It returns a `Table`: the structure, the file's
own properties, the score columns, Consensus and Nearest Reference. Cancel returns the hits so
far (per-row semantics).

## The dialog (Chemistry ▸ Similarity ▸ Similarity Search…)

Kit layout; no helper text, so explanations go in tooltips:

- **Input:** Search In `Segmented` (This Table | A File). Table: Structure column and
  `RowsChoice`. File: a path field with Browse…. The ⚠ problem shows for a missing or
  unreadable file.
- **References:** a short list. Add from SMILES/drawing, Current Row, Selected Rows, or a
  protein's ligand (as Superpose does, needed for 3D). When several are given, a Fuse row:
  Best (max) | Mean. "Selected Rows" with Best gives today's "Max Similarity to Selection" +
  "Nearest Selected".
- **Methods:** one check row per family, each with its own settings in line:
  - Fingerprint: kind and metric. With Tversky, α and β go in More Options.
  - Substructure: MCS | Scaffold | Generic scaffold.
  - Pharmacophore 2D.
  - 3D: Shape | Shape + Features | Field.
  - Combine: None | Mean | Mean Rank | Normalized Mean, shown when two or more methods are on.
- **Output:**
  - Table: score columns (one per method, plus Consensus, plus Nearest Reference), sort by
    the first score, optionally Select rows at or above a threshold or Add them to a list.
    3D: Keep Aligned Poses writes "<column> superposed" with the reference asset.
  - File: Keep top N and a Threshold; the hits open as a new document.
- **More Options (collapsed):** Prefilter N, MCS timeout, 3D conformers and color/field
  weights (superpose `Settings`), Tversky α/β.

The dialog is split to stay well under 1,000 lines: `search.py` (dialog), `references.py`
(reference list widget), `methods.py` (method rows → engine specs).

## Phases (each a commit on `dev`, targeted tests + architecture/boundaries per step)

1. **Engine, 2D:** `similarity/model.py`, `fingerprint.py`, `fusion.py`, `run.score_cells`;
   Tversky; the `pharm2d` kind. Tests in `tests/similarity/`: identity = 1, Tversky
   asymmetry, several-reference max/mean, consensus ranks, Cancelled partial.
2. **Substructure:** MCS Tanimoto (timeout), scaffold and generic scaffold; the prefilter.
   Tests on known pairs.
3. **3D:** `fits_over_references` in `superpose/run.py` (one embedding for all the
   references) and `similarity/shape3d.py`. **Done 2026-11-16:** shape, shape + features and
   field, fused over the references as the search asks, with a flat reference turned away up
   front by `check_references`. Tests: a molecule against itself scores highest, the best of
   two references is chosen, stored conformers give exactly Superpose's score, and the
   superpose tests still pass. **Keeping the aligned poses moved to phase 5**: it needs the
   worker to carry the winning conformer and the reference's asset id, which belongs with the
   dialog's Keep Aligned Poses option and the helper lifted out of `tools/superpose/dialog.py`.
4. **File screening:** `screen_file` with the heap, two stages and the hits `Table`. Tests on a
   small SDF and SMILES file in `tmp_path`: top K is exact, the threshold works, and cancel
   keeps the hits so far.
5. **Dialog, and the aligned poses:** the shared pose-writing helper out of
   `tools/superpose/dialog.py` (`blob_with_reference`, `PutAsset`), the 3D worker carrying the
   winning pose, and Keep Aligned Poses writing "<column> superposed" so the Browser shows the
   reference. Then replace `tools/similarity/search.py`, update `tool.py` and its tip, and update
   `tests/tools/similarity/test_similarity_ui.py` (table run, several references, the file run
   opens a document, keep poses). Update `docs/USER_GUIDE.md` (Similarity Search section;
   keep the `similarity-search` anchor) and ROADMAP §6 (3D similarity done). Before/after
   screenshots of the dialog.

Before Phase 1: save this plan as `docs/SIMILARITY_PLAN.md` and link it from ROADMAP §6.

**Later (not this plan):** 3D pharmacophore point matching (CDPKit's `PharmacophoreAlignment`
where installed, LGPL, optional); ErG and count-based fingerprints; PubChem/ChEMBL similarity
(network, off under policy).

## Notes on conformers (not in this plan; revisit later)

How the superpose engine, and so 3D similarity, uses conformers today:

- **Reference:** only its first conformer, held rigid (`prepare_reference`). A 2D reference
  is rejected.
- **Probe with 3D:** its stored conformers are used as they are. A row with one 3D pose is
  never sampled further.
- **Probe with 2D:** 30 conformers by default, generated with CONFORGE or ETKDG, and not
  saved except for the winning pose.
- **Picking the best:** every conformer is aligned, and the best primary score wins. Field
  refines only the top 3 by shape and ranks by its optimized objective.

Possible options later:
1. **All reference conformers:** use every conformer of the reference and keep the best pair
   (cost: conformers × conformers).
2. **Always generate conformers:** sample even rows that already hold a 3D pose.
3. **File screening:** 3D embedding of the prefilter's survivors dominates the run time, so
   keep the prefilter on by default.

## Constraints to respect

- Every new code path is RDKit + numpy (BSD); no new dependencies.
- American spelling in all text: normalized, color.
- No `getattr`/`hasattr`; no blind excepts (use superpose's `_FAILURES` tuple pattern).
- Per-row results keep finished rows on cancel; a file screen keeps the hits so far.
- `main_window.py` is untouched (the menu entry comes from `tool.py`).

## Verification

- Per phase: `pytest tests/similarity tests/superpose tests/analysis/test_fingerprints*.py
  tests/tools/similarity tests/test_architecture.py tests/test_import_boundaries.py`, plus
  `ruff check` / `ruff format` on the touched paths.
- At the end: the full suite once, since `analysis/fingerprints.py` and `superpose/run.py` are
  shared, and say so.
- In the app (run skill):
  - Search a sample table by ECFP4 + 3D Shape against two references: Consensus and Nearest
    Reference columns appear, and Keep Aligned Poses shows the reference in the Browser.
  - Screen an SDF from `samples/` with 3D on: the hits open as a new document, sorted.
  - Cancel midway keeps the finished rows.
