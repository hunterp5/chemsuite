# Protein alignment: sequences, MSA, superposition and predicted models

Plan approved 2026-11-26. ROADMAP §7a 1 (Alignment and correspondence). A0, the licences and
the benchmark, is next.

## Context

The Protein Viewer has no way to match residues between chains or proteins. It has no sequence
alignment, no MSA and no protein superposition, and Fetch reaches only RCSB. ROADMAP §7a 1 puts
this phase first because ensembles, PLIF across proteins and FEP mapping all need residue
correspondence.

None of it exists in code yet. What we can build on:
- `structure/sequence.py` (`ChainSequence`, `chain_sequences`, `one_letter`) gives each chain's
  SEQRES sequence, with `residue_at` set to -1 where a residue is not modelled.
- `superpose/kabsch.py` is numpy-only and Qt-free. It has `fit(probe, reference) -> (4x4, rmsd)`.
- gemmi 0.7.5's aligner does global alignment only, with match/mismatch or BLOSUM62 scoring. It
  serves as a test cross-check, not as the engine.

**Decided with the user (2026-11-26):**
- **MSA is chemsuite's own numpy progressive aligner,** structure-aware when coordinates exist.
  This holds on one condition: it must measure close to the industry tools on accuracy and speed
  (gates in A0). MAFFT is not bundled. ROADMAP's "MSA through the registered MAFFT tool" line is
  rewritten.
- **Superposition is sequence-based:** Cα pairs from the alignment, fitted by Kabsch with outlier
  trimming, on a whole chain or on a pocket. A sequence-independent fold alignment comes later.
- **Alignments show in the Sequence window:** aligned columns with gaps, plus a conservation
  track. They are kept with the document, and conservation can also color the 3D view.

**What research found on licences and availability (2026-11-26):**
- **MAFFT:** the core is BSD. Its extensions bundle includes Vienna RNA code, which bars
  commercial use, so the extensions are never used. MAFFT runs as a separate program only; it has
  no wheel.
- **pyFAMSA:** GPL-3.0-or-later, runs in-process from abi3 wheels for win/mac-arm64/linux.
- **MUSCLE 5:** GPL-3.0. A Windows binary is unconfirmed.
- **Clustal Omega:** GPL-2+, with no official Windows build.
- **Kalign 4:** Apache-2.0, with no Windows wheel.
- **T-Coffee:** excluded. Its non-commercial restriction conflicts with GPL.
- **BLOSUM/PAM matrices:** the NCBI files carry no licence text. They are widely shipped (the MIT
  `scoring-matrices` package treats them as "believed public domain").
- **Benchmark sets:** BAliBASE, PREFAB, SABmark, HOMSTRAD and `qscore` state no licence at all.

**Goal:** the user can
- align any chains in the document, pairwise or as an MSA;
- see them column for column in the Sequence window, with conservation;
- superpose proteins on a chain or a pocket, with linked ligand rows and docked runs moving in
  the same undo step;
- fetch AlphaFold DB and PDB-REDO models.

This file is the approved plan, linked from ROADMAP §7a 1 and HANDOFF. Phases are named A0 to A6, which avoids HANDOFF calling this
"Phase 1" (in ROADMAP, Phase 1 is the finished viewer).

## Phases (each one a commit on `dev`, usable on its own)

### A0: licences, the benchmark and the kernel micro-benchmark, first (as P0 was for Prepare)
- **Licence decisions** go in CLAUDE.md "Decided":
  - **Matrices:** the numbers are shipped as JSON with their citations (Henikoff 1992, Dayhoff
    1978) and the NCBI source file named. Each file is pinned by SHA-256, and a test checks it.
  - **Benchmark sets:** fetched by a script into `samples/align/`, never committed or shipped,
    used only for measuring.
  - **Comparators:** MAFFT L-INS-i (the accuracy bar), FAMSA and MUSCLE 5. They live in a
    benchmark env only, never in the app.
- **`scripts/review/align_bench.py`**, plus a `-m bench` test.
  - **Sets:**
    - BAliBASE 3 RV11/RV12, the small equidistant sets, closest to 2–50 chains, plus all six
      categories;
    - SABmark Twilight and HOMSTRAD, for the structure-aware mode;
    - a chemsuite-shaped set: one target from several PDB entries, plus AlphaFold models, plus
      homologs. Its columns should come out essentially exact, with SEQRES-versus-observed
      residues handled.
  - **Measures:** SP and TC on core blocks, computed by our own scorer and cross-checked against
    `qscore` on a few families. Wall time.
  - **Kernel micro-benchmark:** a vectorized affine-gap DP at 400 × 400, and the cost of a
    consistency pass at 50 sequences (1,225 pairs). This decides whether A3 is affordable in
    numpy.
- **Output:** `docs/align/a0-review.md`, with comparator numbers measured locally on the Windows
  PC. Gates are set relative to these runs, not to published figures. The proposed gates, to be
  confirmed in the review:
  - **Accuracy:** on RV11/RV12, mean SP within 2 points and TC within 3 points of L-INS-i.
  - **Baseline:** beats our own ClustalW-style baseline in all six categories.
  - **Structure:** with structures, SABmark Twilight and HOMSTRAD score at least as well as
    sequence-only, and better in the bin under 30% identity.
  - **Speed:** 50 × 400 aa in 10 s or less, and no slower than 3× L-INS-i.
- **If the gates fail:** offer MAFFT and MUSCLE 5 as optional `ExternalTool` entries
  (`core/external_tools.py`), or pyFAMSA as an optional pip dependency. That is the user's call
  at that point.

### A1: the pairwise engine, Qt-free (`chemsuite/macromolecule/align/`, a new package)
- **`matrices.py` and `data/*.json`:** BLOSUM45/62/80 and PAM30/70/250. The 20-letter core is
  fixed, with a documented rule for the ambiguity codes (B, Z, X, U, O, J). Add a `pyproject.toml`
  package-data entry.
- **`pairwise.py`:** Gotoh affine-gap alignment in three modes:
  - global (Needleman–Wunsch);
  - local (Smith–Waterman);
  - end-gap-free, the default for chains with missing termini.

  Rows are vectorized in numpy (prefix max via `np.maximum.accumulate`). It is cancellable through
  `CancelToken`. Output: an `Alignment` with identity, similarity and score.
- **`model.py`:**
  - `AlignedChain(protein, chain, columns)`, where each column holds a residue index or -1;
  - `Alignment` (rows, scores, method, parameters), with `to_json`/`from_json`. Saved keys stay
    stable.
- **Tests** in `tests/macromolecule/align/`:
  - textbook cases;
  - agreement with gemmi's `align_string_sequences` on global BLOSUM62;
  - symmetric and degenerate inputs;
  - a matrix-hash test.

### A2: progressive MSA (`align/tree.py`, `align/profile.py`, `align/msa.py`, `align/conservation.py`)
- **Guide tree:** UPGMA on all-pairs identity. At N ≤ 50 that is 1,225 pairs, so no k-mer
  shortcut is needed; k-mer distances are used only above roughly 200 sequences.
- **Profiles:** sequence weights, position-specific gap penalties, and profile–profile DP. This
  is the step from ClustalW to Clustal Omega, TC 0.37 → 0.55 on BAliBASE.
- **Iterative refinement:** split the tree at each edge, realign the two halves, and keep the
  result if the score improves. About +7 TC in MAFFT's own figures.
- **Conservation:** a per-column score, Jensen–Shannon divergence with a gap penalty.
- **Measured against A0's bench.** Each step is kept only if it helps.

### A3: consistency and structure (only if A0's kernel numbers allow it)
- **Consistency:** extend pair scores through every third sequence (T-Coffee/L-INS-i style) from a
  library of all-pairs alignments. This is the largest accuracy gain, roughly +7–8 TC, and the
  costliest step.
- **Structure-aware:** when the chains have coordinates, each pair is superposed (A4's core). Its
  Cα correspondences, within about 3.5 Å after the fit, enter the library as extra, heavily
  weighted pairs. 3D-Coffee saw +4 to +10 points this way.
- **Gate:** A0's structure gate. If consistency is too slow in numpy, the structure library
  still feeds the A2 profile scores, and the A0 review records that.

### A4: superposition, Qt-free (`align/superpose.py`)
- Cα pairs come from an `Alignment`. The fit reuses `chemsuite/superpose/kabsch.fit`, with
  PyMOL-style iterative trimming (cutoff 2 Å, 5 cycles). It returns
  `Fit(transform, rmsd, pairs used, pairs total)`.
- **Scope:** the whole chain or chains, or a pocket, meaning the reference residues within 6 Å of
  a ligand or the selection. Both values are settable.
- **`transform_structure(s, T)`** calls `s.with_coords(...)` over every model.
- **Tests:**
  - a known rotation is recovered;
  - trimming drops a loop that was moved;
  - the RMSD matches gemmi's `superpose_positions` (× √3).

### A5: in the app
- **Superposition** goes in a new `workspace/superpose.py`, because `workspace/session.py` is 991
  lines. It makes one `session.edit(label, {mobile: moved}, extra_ops)`. The extra ops:
  - **Linked rows' Ligand poses move.** Add a new `LinkedRow.move_poses(protein, T) -> ops` in
    `ui/host.py`, implemented by `TableRow` in `tools/macromolecule/follow.py` with
    `ctx.cell_ops(LIGAND_COLUMN, …)` and `superpose/model.apply_transform`. It is injected like
    `relabel`, the precedent being `session._relabelled` and `tools/macromolecule/tool.py`.
  - **Docked runs and boxes move.** A new `Op` with an inverse transforms `workspace.poses`, so
    undo moves them back.
- **Saved alignments:** `workspace/alignments.py` adds an `alignments` key to
  `MacromoleculeWorkspace` (`_load`, `capture`, `take`, `replace`).
  - Rows whose protein is gone are dropped.
  - The key is optional, so there is no version bump.
  - It bumps `revision`.
- **Sequence window** (`ui/sequence_view.py`, 422 lines; split out `ui/aligned_rows.py`):
  - an aligned mode, with the rows sharing column space and gaps drawn as `-`;
  - a conservation track shaded per column;
  - hover, selection and `show_residue` work through each row's per-column `residue_at`;
  - the Sequence menu gets Alignment ▸ (None or a saved alignment) and Show Conservation;
  - the alignment choice is saved in the window's `state`.
- **3D color:** a Conservation scheme from per-residue values (`look/colours.py`, or
  `Display.colour` per residue).
- **Dialogs**, both built with the ToolDialog kit, with tooltips and no helper text. They are
  added to the viewer's Tools menu through four steps: a `menubar.py` constant, a `panel.py`
  re-export, a `tool.py` command, and a lazily imported opener.
  - `Align Sequences…`: pick chains; Method (Global, Local, End-gap-free or MSA); Matrix; gap
    open and extend under More Options.
  - `Superpose Proteins…`: reference and mobile proteins and chains; Whole Chain or Pocket;
    cutoff and cycles under More Options. RMSD and pairs go to the status line, and the
    alignment it used is saved.

  Both run as jobs through `ctx.run` and can be cancelled.
- **Tests:**
  - `tests/macromolecule/workspace/test_alignments.py` (round trip);
  - `tests/macromolecule/ui/test_sequence_view_aligned.py`;
  - `tests/tools/macromolecule/test_superpose.py`: one Undo restores the protein, the rows'
    poses and the docked runs.

### A6: predicted and re-refined models
- **`structure/fetch.py`:**
  - AlphaFold DB by UniProt accession, CC-BY-4.0, with the attribution kept in `meta`;
  - PDB-REDO by PDB ID, licence checked first.
- **FetchDialog** gets a Source choice. For AlphaFold models, the B-factor field holds pLDDT, so
  the B-factor color scheme shows confidence.
- **Tests** use a fake opener, with no network.

Each phase updates `docs/USER_GUIDE.md` (A5 and A6), and ROADMAP and HANDOFF, in the same commit.

## Files
- **New:**
  - `chemsuite/macromolecule/align/{__init__,matrices,pairwise,model,tree,profile,msa,conservation,superpose}.py`
    and `align/data/`;
  - `macromolecule/workspace/{alignments,superpose}.py`;
  - `macromolecule/ui/aligned_rows.py`;
  - `tools/macromolecule/{align_dialog,superpose_dialog}.py`;
  - `scripts/review/align_bench.py`;
  - `docs/align/a0-review.md`.
- **Changed:**
  - `ui/{sequence_view,menubar,panel,host}.py`;
  - `workspace/workspace.py`;
  - `tools/macromolecule/{tool,follow,dialogs}.py`;
  - `structure/fetch.py`, `look/colours.py`;
  - `pyproject.toml`;
  - `docs/{USER_GUIDE,ROADMAP}.md`, `HANDOFF.md`, `CLAUDE.md` ("Decided");
  - `tests/test_import_boundaries.py`, adding `align` as a Qt-free package.
- **Reused:**
  - `superpose/kabsch.py`, `superpose/model.apply_transform`;
  - `structure/sequence.py`, `Structure.with_coords`;
  - `session.edit(extra_ops=)`, `ctx.cell_ops`;
  - `core/download.fetch_bytes`, `core/jobs` (`CancelToken`, `Cancelled`).

## Verification
- **Per phase:** run the new tests plus `tests/test_architecture.py` and
  `tests/test_import_boundaries.py`, then `ruff check` and `ruff format` on the touched paths.
- **A0, A2 and A3:** run `python scripts/review/align_bench.py` against the gates and record the
  numbers in `docs/align/a0-review.md`.
- **By hand in the app:**
  - fetch 1ATP and the AlphaFold model of PKA (P17612), align them, and superpose on the pocket;
  - the linked rows' poses stay in the pocket, and one Undo restores everything;
  - save and reopen keeps the alignment and the Sequence window's mode.
- **Full suite:** run once, before the merge to `main` (this touches document save/load and undo).
