# Loop modeler for Prepare: closing gaps longer than one residue

Plan approved 2026-11-19. ROADMAP §7a.3 (Repair). Nothing is built yet: L0 is next.

## Context

Prepare builds a residue only when it is missing alone (`repair.fill_single_gaps`,
SIMULATE_PLAN M1). Every longer gap is capped with ACE/NME and marked ⚠. PREPARE_PLAN step 2 and
ROADMAP §7a.3 promised "short loops closed", but no milestone was ever planned for it. A capped
gap leaves the protein with two loose ends: wrong electrostatics near the gap, a hole in the
surface, and docking or MD done on a protein that is not there. Most crystal structures have at
least one short disordered loop, often near a pocket.

**Goal:** close gaps of 2 to N residues (N set by the benchmark, about 8–12) in the Qt-free
Prepare, with:
- geometry that closes exactly onto both anchors;
- no clashes;
- backbone angles in allowed Ramachandran regions;
- side chains from the existing templates.

Each loop is then relaxed by the existing Minimize, which already frees every atom Prepare
built (`Prepared.built`). Longer gaps and missing terminal tails stay capped, as now.

The work is our own code from published algorithms: CCD, Canutescu & Dunbrack 2003; Ramachandran
classes, Lovell et al. 2003. MODELLER is excluded (not free for commercial use). PDBFixer (MIT)
is used only as a benchmark comparator, in the OpenMM environment, never in the app.

## Milestones

### L0 — the benchmark, first (as P0 was for Prepare)

`scripts/review/loop_bench.py`, with data under `samples/loops/` (not committed).
- **Set:** high-resolution (≤ 1.8 Å), nonredundant chains with no gaps and no alternate
  locations on the backbone, picked through RCSB's search API and its 30% sequence clusters
  (CC0). Per length 2, 3, 4, 6, 8, 10 and 12: about 40 segments cut at random (seeded) away
  from chain ends. Cutting deletes the residues and keeps SEQRES, so Prepare sees a real gap.
  Train/test halves by CRC32 of the ID, as `prepare_review.training`.
- **Measures, per loop:**
  - backbone RMSD (N, CA, C, O) to the deposit, with no superposition, since the anchors are
    fixed;
  - closure: junction peptide bond and angles;
  - heavy-atom clashes (< 2.2 Å severe, < 3.0 Å contacts);
  - Ramachandran outliers;
  - time.
- **Comparators:**
  - caps only (today's behavior, as the floor);
  - PDBFixer `addMissingResidues` in `cs-sim-m0` (benchmark-only, like PROPKA and Protoss);
  - each of those, and ours, after our Minimize.
- **Output:** `docs/loops/l0-review.md`, with the gates for L1 set from what PDBFixer and the
  caps score.

### L1 — closure, Qt-free (`chemsuite/macromolecule/prepare/loops.py`, new; `repair.py` is 443 lines)

- **Ramachandran data:** `scripts/build_ramachandran.py` reads φ/ψ from the L0 training set's
  deposits (B-factor filtered). It writes `prepare/data/ramachandran.json`, 10° log-probability
  grids for four classes: general, Gly, Pro, pre-Pro. The file is small and derived from CC0
  PDB data, rebuilt only by the script.
- **Algorithm, per gap** (anchors: residue *a*, open C; residue *b*, open N; *m* missing
  residues named from `Structure.sequence_of(entity)`):
  1. Free torsions: ψ of *a*, φ/ψ of each new residue, φ of *b*. ω is trans; *a*'s O is
     rebuilt in the new peptide plane.
  2. Sample K candidates (about 200, vectorized over candidates in numpy). Each draws φ/ψ
     from the class's grid and builds forward from *a* with `templates.place` and the bond
     constants already in `repair.py`.
  3. Close each one by CCD onto *b*'s N, CA and C, with moves weighted by the Ramachandran
     grid. Keep closures under 0.1 Å RMSD over the three atoms.
  4. Score: clashes against the environment's heavy atoms (backbone and a CB pseudo-atom),
     Ramachandran log-probability, then closure.
  5. For the best few, build side chains with `repair._build` and turn them clear with
     `repair._turned_clear`. Re-score with side chains and keep the best.
  6. Insert the residues with `edit.insert_residue`. Numbering comes from the gap in
     `res_seq` where it fits, else `repair._free_icode`.
- **Deterministic:** the random generator is seeded from the gap's label, so preparing again
  gives the same loop. This matters because documents reopen exactly as saved.
- **Models:** each model of a multi-model file is closed on its own anchors.
- **Lone residues:** a one-residue gap goes through the same path once L0 shows it does at
  least as well as `_fill`. Otherwise `_fill` stays for m = 1.
- **Unit tests** (`tests/macromolecule/prepare/test_loops.py`), on a fixture protein with a
  loop cut out:
  - closes exactly;
  - no severe clash;
  - φ/ψ allowed;
  - same result twice;
  - RMSD under a generous bound;
  - a gap too long for the number of missing residues is capped instead.

### L2 — in Prepare, Inspect and the dialog

- `pipeline.Options`: keep `fill_gaps` (saved key; tooltip "Missing Loops") and add
  `loop_max: int` (the longest gap built; 1 = today's behavior). `settings_of` and
  `options_from` in `tools/macromolecule/prepare_options.py` read and write it, clamped to
  1–12. Older documents take the default.
- `pipeline.prepare`: `repair.fill_single_gaps` becomes `loops.fill_gaps(s, loop_max, cancel)`
  (cancellable; `check(cancel)` per gap). Caps go on what is left.
- **Report:** one `FILLED` finding per gap, "Loop of n residues built", ⚠, selecting the
  loop.
- **Inspect:** each `GAP` row gets a box, as a missing side chain has one. Its label goes into
  `Decisions.repair` (unticked: left capped). `review.decisions` and the saved review carry it.
- **Prepare dialog:** a "Longest Loop" spin box under More Options, beside Lone Residues, with
  its help in a tooltip (no helper text). Minimize's readiness check is unchanged.
- **Default on** for gaps up to the length L0/L3 shows is safe. The decision is recorded in
  the review doc.
- `docs/USER_GUIDE.md`'s Repair paragraph is updated, in the same commit.

### L3 — measured, then default

- L0's benchmark is run on the test half: ours, ours + Minimize, PDBFixer, and PDBFixer +
  Minimize. It goes in `docs/loops/l3-review.md`.
- **Gates for default-on at a length:**
  - every loop closed;
  - no severe clash after Minimize;
  - Ramachandran outliers below 5%;
  - median backbone RMSD no worse than PDBFixer's at that length.

  `loop_max`'s default is the longest length that passes all four.
- One by-hand run in the app: Prepare + Minimize on a structure with an 8-residue gap, looked
  at in the viewer.

## Files

- **New:**
  - `chemsuite/macromolecule/prepare/loops.py`
  - `chemsuite/macromolecule/prepare/data/ramachandran.json`
  - `scripts/build_ramachandran.py`
  - `scripts/review/loop_bench.py`
  - `tests/macromolecule/prepare/test_loops.py`
  - `docs/LOOPS_PLAN.md` (this plan)
  - `docs/loops/`
- **Changed:**
  - `prepare/pipeline.py` (Options, the repair stage, the report)
  - `prepare/repair.py` (shared helpers stay here; the module docstring)
  - `prepare/inspect.py` (FILLED's comment)
  - `prepare/review.py` (gap boxes)
  - `tools/macromolecule/prepare_options.py`
  - `docs/USER_GUIDE.md`, `docs/ROADMAP.md`, `HANDOFF.md`
- **Reused:**
  - `templates.place`
  - `repair._build`, `_turned_clear`, `_crowding`, `_free_icode`, and its bond constants
  - `edit.insert_residue`
  - `inspect.ends` and `_missing_between`
  - `Prepared.built` (the loop is free in Minimize)
  - `simulate.minimise`

## Verification

- `pytest tests/macromolecule/prepare/test_loops.py tests/macromolecule/prepare/test_pipeline.py tests/macromolecule/prepare/test_review.py tests/tools/macromolecule tests/test_architecture.py tests/test_import_boundaries.py`
- `pytest -m openmm tests/simulate`, with `CHEMSUITE_OPENMM_PYTHON` set: Prepare + Minimize
  on a structure with a built loop.
- `python scripts/review/loop_bench.py score --test`: the L3 table.
- `ruff check` and `ruff format` on the touched paths.
- The full suite once, before the merge to `main`, since the change touches Prepare's
  pipeline.
