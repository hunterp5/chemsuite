# Simulation plan: restrained minimisation after Prepare, on a worker built for MD, MM-GBSA and FEP (approved 2026-11-15)


## Context

Protein ▸ Prepare (ROADMAP §7a phase 3, `docs/PREPARE_PLAN.md` P0–P4) repairs the structure,
places hydrogens and picks protein, ligand, water and metal states with the coupled site graph.
It never moves a heavy atom and never assigns force-field parameters, so new hydrogens and
rebuilt side chains sit at ideal geometry, strained against their neighbours, and nothing tells
the chemist whether the complex would survive a simulation.

This adds the last step: a **restrained energy minimisation** of the prepared complex, run in
OpenMM's own interpreter, which is also the first use of the **simulation worker** that phases 5
and 6 (MD, MM-GBSA/PBSA, free energy) will drive. The point of the design is that
parameterisation happens once, is recorded with the document, and is reused unchanged by every
later method — a minimised complex, its MD and its MM-GBSA score all come from the same system.

Decided by the user in planning:
- The worker runs on a **native env where one exists, with an optional WSL interpreter** as a
  second configured path (needed later for openfe, MM-PBSA and AmberTools-only methods).
- **No in-app MMFF tier.** MMFF94 is not validated for proteins and is gas-phase; without the
  worker, Prepare behaves as today and the report says minimisation was skipped.
- The minimised coordinates become a **new model beside the prepared one**, not a replacement.
- Default force field: **ff14SB + OpenFF Sage 2.3 + TIP3P** (openfe's benchmark combination;
  Sage was fitted against TIP3P; GBn2 is validated with ff14SB). ff19SB + OPC is offered too.
- Minimisation is **on by default when the worker is configured**, with a checkbox to turn it off.
- Hardware: NVIDIA GPU on the Windows PC, Apple silicon Mac, a Linux cluster later.
- The environment is a **pinned env file per platform plus a USER_GUIDE section**; Settings ▸
  Test reports the versions it finds. No app-managed install for now.

## What exists, and the gaps this plan must close

- `prepare()` (`macromolecule/prepare/pipeline.py:93`) returns
  `Prepared(structure, states, report, choices, pka, population)`. `states` holds Amber-style
  labels (HID/HIE/HIP, ASH, GLH, LYN, TYM, CYM, "NH2") **but `res_name` is never changed**, and
  `prepare_dialog.py:168-176` throws `states`, `pka`, `population` and `choices` away: only the
  structure and the report rows survive. The worker must therefore be told the states
  explicitly, and must cross-check them against the hydrogens actually present.
- `Structure` (`structure/model.py:79`) holds `coords` as `(k, n, 3)`, so a minimised model fits
  beside the prepared one; `with_coords` (`model.py:226`) replaces coordinates in bulk, and uids
  give a one-to-one atom mapping.
- `io.to_text` writes **no CONECT, SSBOND or LINK records and no bond orders**
  (`structure/io.py:320`), so the topology must travel as data, not as PDB text.
- `ligands.to_mol` (`chem/ligands.py:77`) returns heavy atoms only, with the CCD's charges, not
  the chosen protomer. `ligand_sites._drawn` (`ligand_sites.py:132`) is the code that applies the
  structure's own H counts and charges; a public helper is missing.
- `retro/runner.py` and `retro/worker.py` are the worker pattern to generalise: `Worker`,
  `worker_for`, `child_environment`, `probe`, and the `os.dup2(2, 1)` trick that keeps a
  library's prints out of the reply stream. Its protocol has **one reply per request and no
  in-flight progress**, which MD needs.
- Ratchet budgets (`tests/test_architecture.py:156`) are **exactly full**: `getattr`/`hasattr`
  54/54, `except Exception` 16/16, `type: ignore` 29/29. The regexes scan file text, so tokens
  inside an embedded worker script count. The worker is written with narrow excepts and
  `isinstance` checks; no budget is raised.

### Prepare bugs and gaps found while planning (fixed in M1, each with its own test)

1. **A chain end before a gap has no OXT** (`inspect.py:116-124` only builds OXT for
   `last_in_chain`), so the residue ends at a bare carbonyl carbon that matches no Amber
   template. Gaps get neutral **ACE/NME caps**, which is also what the GAP finding promises.
2. **A truncated chain end is treated as a real terminus.** Compare against SEQRES
   (`seqres`, `entity_seq`) and cap an end that is merely the last observed residue.
3. **Waters outside 5 Å of a ligand have no hydrogens** (`waters.py:146`); every kept water
   needs them for MM.
4. **Ion charges come from the file** (`io.py:110`), which is often 0 while `energy.py` assumes
   +2. Assign from a table by element and oxidation state, and refuse an unknown metal.
5. **`"NH2"` as the label for a neutral N-terminus collides with the NH2 cap residue name** —
   rename the state label.
6. **A public ligand-form helper** (`chem/ligands.py`): the chosen protomer/tautomer as an RDKit
   mol with explicit hydrogens at the Structure's own coordinates, with `uid` props on every
   atom. Built from `to_mol` + the `_drawn` logic + `sites.bond_graph`; never `AddHs`.
7. The docking export (`dock/pdbqt.py:92`) must keep the new caps, or name them in
   `Receptor.dropped`.

## M0 — feasibility gates (no app code; a spike script and a written result)

Each gate can change the design, so none of M1–M5 starts until M0 is written up in
`docs/simulate/m0-review.md`.

1. **Native Windows env.** Build a win-64 env with `openmm`, `openmmforcefields`,
   `openff-toolkit`, `openff-nagl`, `openff-nagl-models` (RDKit backend, no AmberTools), and
   confirm which OpenMM platforms appear (CUDA on the NVIDIA PC). Same for osx-arm64.
2. **Ligand charges without AmberTools.** Sage 2.3 uses AshGC/NAGL charges; openmmforcefields'
   support for that is unconfirmed. The design does not depend on it: the worker computes the
   charges itself (NAGL), then hands the template generator a `Molecule` that already carries
   nonzero charges, which it uses as given. **The charges are stored in the document** so MD and
   MM-GBSA reuse exactly the same ones.
3. **Implicit solvent with a Sage ligand.** Check whether `implicit/gbn2.xml` assigns GB
   parameters to SMIRNOFF/GAFF atom types through `SystemGenerator`. If it does not, the
   fallback is a **restrained explicit TIP3P shell, stripped after minimisation** — arguably the
   more honest route to "simulation-ready" anyway. Decide here, and record the reason.
4. **Template coverage.** Enumerate what ff14SB and ff19SB in OpenMM provide for every state
   Prepare can pick (TYM, CYM, CYX, neutral termini, ACE/NME, HIP…). Any state with no template
   **blocks with a ⚠ problem** and is never silently substituted.
5. **Massless atoms vs restraints.** Confirm whether `LocalEnergyMinimizer` leaves zero-mass
   particles fixed. If not, the hydrogen-only stage uses stiff positional restraints instead.
6. **WSL mechanics.** Path translation (`/mnt/c`), and whether a Job Object kill of `wsl.exe`
   stops the Linux process — the worker also exits when stdin closes, so it cannot be orphaned.

## M1 — the system spec (Qt-free, `chemsuite/simulate/`)

**Done 2026-11-15.** The Prepare fixes 1–7 above, each with its tests, plus what they turned
up:
- **A residue missing alone is built, not capped** (the user's choice, 2026-11-15): two caps
  across a one-residue gap meet where its CA would be (methyls 1.6–1.9 Å apart).
  `repair.fill_single_gaps` builds its backbone from its neighbours' peptide planes, its
  side chain from its template, and reports it with ⚠.
- **Rebuilt side chains are turned clear** of their neighbours (χ1, χ2 in 30° steps): a
  side chain at the template's ideal angles sat 1.6 Å from an Arg, and Meeko read bonds into
  the contact.
- Report lines now select the right residue after caps are put in or waters dropped (they
  were found on the input structure and kept its indices).
- `simulate/spec.py` (`SystemSpec`, `ForceFieldSpec`, `MinimiseSpec` with the M3 stages),
  `build.py` (forms read from the hydrogens; `StateMismatch`; blocking `Problem`s, plus bond
  orders guessed from geometry, nucleic acids, glycans and alternate locations),
  `result.py` (`MinimiseResult`, coordinates by uid, `warnings`);
  `ligands.prepared_form` is fix 6.

A new Qt-free top-level package, registered in `tests/test_import_boundaries.py:37` `QT_FREE`.
It imports `protein.structure`, `protein.chem` and `core.jobs`; never Qt, `core.table` or OpenMM.

- `spec.py`: frozen dataclasses, JSON in and out — `SystemSpec` (atoms: element, name, residue,
  chain, uid; the explicit bond list; per-residue state labels; ligand molblocks with explicit H
  and `uid` atom props; per-ligand charges once computed; ions with assigned charges; disulfides;
  caps), `ForceFieldSpec` (protein FF, small-molecule FF, water model, solvent treatment),
  `MinimiseSpec` (restraint constant, tolerance, max iterations, the stage schedule).
- `build.py`: `system_spec(structure, states, ...) -> SystemSpec` — the Prepared → spec
  translation, including the states cross-check against the hydrogens present (a disagreement is
  a hard error, not a warning).
- `result.py`: `MinimiseResult` (new coordinates by uid, energy before/after, heavy-atom and
  ligand RMSD, worst remaining clash, per-atom max force, the charges used, FF versions, the
  platform and timings).
- Prepare's own `relax.py` is **not** created; `docs/PREPARE_PLAN.md:143` and ROADMAP §7a 3 are
  rewritten to point here.

## M2 — the worker and its runner (`chemsuite/simulate/`)

**Done 2026-11-15.** As below, with:
- `core/jobs/lines.py` (`LineWorker`, `WarmWorkers`, `child_environment`, `is_app_python`);
  retro is a `LineWorker` subclass, its API unchanged. A cancelled request that cannot stop
  is now killed at once (it waited 5 s for a quit the busy worker never read), and a worker
  that dies is reported with its exit code.
- The worker's `minimise` runs the stages already (OpenMM's `MinimizationReporter` gives the
  energy events and the stop between iterations); M3 is the result (RMSD, contacts, what
  moved), the tuning and the benchmark. First run, T4 lysozyme L99A with its ligand, 200
  iterations a stage on the CPU: heavy-atom RMSD 0.24 Å, but surface side-chain tips swing
  up to 1.1 Å in the free side-chain stage (Glu 108 OE1, Lys 124 NZ): M3 should weigh a light
  restraint there.
- `build` gives every residue but a small molecule its template by name (ions too: FE and FE2
  are the same atom to a graph), refuses a molecule with no charges, and checks every
  residue's total charge in the system against the spec's formal charge.
- Env files in `envs/` (8.4.0 pinned, provisionally); the Mac one is untested.
- Tests: `tests/core/test_lines.py` (a plain-Python fake worker, so they run on Windows),
  `tests/simulate/test_runner.py` (a fake OpenMM worker), and `tests/simulate/test_worker.py`
  against a real env (`CHEMSUITE_OPENMM_PYTHON`, marker `openmm`). Retro's own worker tests
  still need a POSIX shell: run `tests/retro` on the Mac before the merge.

- Extract the JSON-lines client from `retro/runner.py` into a shared
  `core/jobs/lines.py` (`LineWorker`, `worker_for`, `child_environment`, `probe` helpers) and
  **move retro onto it**, so there is one implementation, not two. Protocol additions needed by
  MD: **event lines** (progress, energies, frames) interleaved before a request's reply, and a
  **cooperative stop** message so a run can end at a checkpoint instead of being killed.
- `worker.py`: the script run with `python -I -c`, keeping stdout clean via `os.dup2(2, 1)`.
  Ops: `hello` (versions, platforms, GPU), `build` (spec → System, cached), `minimise`,
  later `md`, `mmgbsa`. No `except Exception`, no `hasattr`.
- `runner.py`: the app-side client; one warm worker per interpreter, closed from
  `app/shell/main.py` beside `retro.close_workers()`; `children.start`/`stop` as the rule
  requires (`test_external_programs_start_through_the_registry`).
- `config.py`: two new external tools — `openmm_python` and an optional `openmm_wsl_python` —
  registered in `core/external_tools.py:47`, with a `probe` behind Settings ▸ Test reporting
  OpenMM, OpenFF and force-field versions and the available platforms. `ready_problem`-style
  functions give the dialog its blocking ⚠ text.
- Tests: a **Python fake worker** (the shell-script fake in `tests/retro/fake.py` skips on
  win32, so protocol tests never run on the Windows PC), plus a marker for the few tests that
  need a real OpenMM env.

## M3 — the minimisation protocol

**Done 2026-11-15** (the app side, Qt-free; the UI is M4). As below, with:
- `Prepared.built`: the uids of the heavy atoms Prepare built (filled residues, rebuilt side
  chains, caps, OXT), to be kept with the document beside `states` in M4. They are never
  restrained and are left out of every RMSD, since no deposited position says where they
  belong. An explicit set, not a uid cutoff, so a later **Protein ▸ Minimise…** on a reopened
  structure still knows them.
- `spec.Stage` holds a restraint per class of heavy atom (`backbone`, `side_chains`,
  `ligands`, `waters`, `ions`; `RESTRAINED`) and moves `hydrogens` (heavy atoms massless) or
  `heavy atoms`. Stage 2 now holds ligands and ions at 10 kcal/mol/Å² and waters at 1, so
  nothing unbonded drifts or collapses while the side chains settle. Restraints always pull
  to the deposited start.
- `MinimiseSpec.tighten` (3): while the held heavy atoms end further than `rmsd_limit`, the
  last stage runs again from where it began with its restraints doubled (×2, ×4, ×8), each run
  named in the stages. The worker measures the RMSD and returns it.
- `tolerance` is OpenMM's: an RMS force, not the largest one. The largest left is reported.
- `simulate/minimise.py`: `minimise(worker, structure, system, …, built=)` charges, runs and
  measures (`measure`): RMSD overall / pocket (deposited heavy atoms within 5 Å of a small
  molecule at the start) / ligand, the closest contact before and after (heavy atoms of
  different residues, no ion, not bonded and sharing no bonded neighbour: an amide O and the
  next N are 2.25 Å apart), the five residues that moved most, where the largest force is
  left, the charges and versions. A cancelled run raises `Cancelled` with what it got to,
  measured.
- **Tuning, T4 lysozyme L99A with its ligand and crambin (GBn2, no waters kept, CUDA, under a
  second a stage):** heavy-atom RMSD 0.22 / 0.17 Å, ligand 0.17 Å, pocket 0.18 Å. The surface
  tips that swing up to 1 Å (Glu 108 OE1, Asn 116, Lys 124, Asp 127) move in stage 2. A
  side-chain restraint there of 0 to 2 kcal/mol/Å² changes the largest move by only 0.06 Å:
  exposed carboxylates relaxing in implicit solvent (Glu 108's moves off Gly 77's carbonyl O,
  3.3 Å away in the crystal), not salt bridges forming. Holding them would mean freezing the
  side chains, so stage 2 keeps 1 kcal/mol/Å². Whether kept crystal waters damp it is a
  question for the M5 benchmark.

Staged, so that nothing moves before its neighbours are sane:

1. Hydrogens only (heavy atoms fixed or stiffly restrained — M0 gate 5).
2. Rebuilt and flipped side chains and kept waters, protein backbone restrained.
3. The whole complex with harmonic positional restraints on heavy atoms resolved from the
   deposited coordinates, tightened until heavy-atom RMSD ≤ 0.3 Å; the ligand restrained more
   weakly so it can settle in the pocket.

Reported: energy before and after, RMSD overall / pocket / ligand, the worst clash before and
after, residues that moved most, and every force-field and charge version used. A run that ends
above the RMSD limit, or with a remaining clash, is marked ⚠ rather than quietly accepted.

## M4 — UI, storage and invalidation

**Done 2026-11-15.** As below, except that the result is a **new protein**, not a model of
the old one: Prepare with Minimise ticked leaves the protein as it was and adds `X-prep` and
`X-min` in one undoable step (`ctx.edit`); Protein ▸ Minimise… adds `X-min` (then `X-min-2`…)
for a prepared protein. Each copy names its source (`derived_from` in the asset meta) and is
shown in its place (`WorkbenchSession._swap_sources`, undone both ways). Each keeps its own
`PreparedRecord` (the minimised one with a `Minimisation`: source, force field, system, the
result's numbers) and report. Invalidation is the record's atom digest: a protein edited after
preparation must be prepared again, and a minimised one is never prepared or minimised again.
`MinimiseSettings` (force field, GBn2/vacuum, restraint scale, RMSD limit, platform) is shared
by both dialogs (`tools/macromolecule/minimise_options.py`); ff19SB is refused while crystal waters
are kept (OPC's virtual site is not placed). Prepare's Minimise is off by default until M5.
Code: `tools/macromolecule/minimise.py`, `minimise_dialog.py`; tests
`tests/tools/macromolecule/test_minimise_ui.py` (the fake worker).

- `tools/macromolecule/prepare_dialog.py`: a **Minimise** checkbox in Method (on when the worker is
  configured; its ⚠ problem text when not), and the force-field choice under More Options
  (ff14SB+Sage+TIP3P or ff19SB+OPC). No helper text — explanations go in `help=` tooltips.
- A separate **Protein ▸ Minimise…** command for an already-prepared structure, registered in
  `TOOLS` (`tools/macromolecule/tool.py:358`) and the viewer's command map.
- Run through `ctx.run(..., keep_partial=False)` (a whole-structure result): `job.report` for
  progress, `job.emit_partial` for live energy.
- Applied as **one undoable edit** via `WorkbenchSession.edit` (`workspace/session.py:332`),
  adding the minimised model to `coords`; the models bar names it, and the new model is shown.
- The report extends `PrepareReport`/`Row` (`workspace/reports.py:33`) with the minimisation
  rows, saved under the protein workspace. `Prepared.states`, `pka` and `population` are
  **persisted too** (they are currently dropped), since MD and MM-GBSA need them.
- **Invalidation:** re-running Prepare or editing the structure drops or flags the minimised
  model and its stored system; multi-model (NMR) input minimises the shown model only.
- **Found in M3, to hold in M4:**
  - `Prepared.built` is relative to `prepare()`'s input, and Prepare and Inspect's Apply both
    prepare the structure *shown*, which may be prepared already: a second run would find no
    caps or rebuilt side chains to call built, and restrain them where they were built. Keep
    `built` with the document and take the union with the stored set on every re-prepare
    (uids of atoms since deleted drop out).
  - Minimise from the **prepared** model, never from a minimised one: restraints and RMSD
    are measured from the model passed in, so a second minimisation from the first would
    drift away from the deposit.

## M5 — benchmark and acceptance (`-m bench`, not the regular suite)

On the P0 neutron and named-case sets (`tests/macromolecule/test_bench_prepare.py` alongside):
heavy-atom RMSD to the deposit, clash score, hydrogen-position accuracy before and after
(minimisation must not spoil the states Prepare chose), and ligand RMSD. Plus the real test of
"simulation-ready": a **short explicit-solvent MD (100 ps–1 ns)** on each benchmark complex runs
without instability and keeps the ligand in the pocket. It becomes Prepare's default step only
if the clash score improves, the RMSD limit holds, the neutron marks do not regress, and the MD
runs clean.

**M5a done 2026-11-24** (the minimization metrics; M5b, the MD gate, is planned with §7a 5).
`scripts/review/minimise_review.py run` prepares each of the 262 P0 neutron entries as the
`chemsuite` variant and minimizes it with `MinimiseSettings()` (ff14SB + Sage + TIP3P, GBn2,
CUDA), writing both models from the same run (`out/chemsuite-prep`, `out/chemsuite-min`) and
the numbers (`out/minimise/<ID>.json`); `score` reports them.
`tests/simulate/test_bench_minimise.py` (`-m bench`) reads that run and gates on it. About
1.5 h on the Windows PC. The slow part is Prepare's ligand pKa (up to 10 min on some ligands,
about 13 s per protonation form), not OpenMM (median 7 s).

- **Coverage:** 167 of 262 minimized. 76 fail by design: 23 nonstandard residues (LLP, PCA,
  SUI, CRQ…), 18 hemes (OpenFF's charge model refuses iron), 13 glycans, 11 covalent
  ligands, 7 nucleic acids, 4 ligands with atoms missing. **19 fail for a defect**, the
  test's ratchet:
  - ligand hydrogen missing (6): the carbonic anhydrase sulfonamides bound to Zn (4G0C,
    a named case, 5C8I, 6BBS, 6BC9, 6BCC; "an atom is missing a hydrogen (S)") and 13HO;
  - Prepare picks tyrosinate (TYM), which has no ff14SB template (5): 5TKI, 7T5D, 7T5E,
    8E1W, 8RK1;
  - a chain end left without OXT or cap (3): 1IU6, 3VXF, 5KWF;
  - Mn(III) with no ion template (2): 6KK8, 7KKS;
  - one each: an Arg with 22 atoms for a 24-atom template (2INQ), a ligand topology OpenFF
    refuses (5VG1), NaN coordinates in the worker (6H1M).
  Hemes fail in the worker with OpenFF's raw `ValueError`; build should name them as a
  `Problem` instead.
- **Heavy atoms held** (from the prepared start): RMSD median 0.26, p95 0.30, max 0.41 Å; two
  entries over the 0.3 Å limit (2XQZ, 3INS), nine with a ⚠. Pocket median 0.31, p95 0.62 Å;
  ligand median 0.38, p95 1.04, max 2.79 Å. The ligands that move most are additives in the
  solvent (FMT, MPD, NO2, PEO, D3O in 9KVM, 9KVL, 4AR4, 9BVY); real ligands moving about 1 Å
  are ADP in two kinases (21FE, 7F4X) and 8DG in 9MBE. Whether stage 3 holds ligands tightly
  enough is a question for M5b.
- **Clash score** (our own, Probe-like: pairs more than three bonds apart, hydrogens
  included, overlapping by 0.4 Å or more, H-bonds excepted, per 1000 atoms): median
  28.3 → 2.3, p95 65.7 → 10.2; worse after in no entry.
- **Neutron marks** (the same 167 entries): unchanged. Asn orientation 85.6 → 85.8%, ligand
  polar-H count 95.9 → 95.6%, every other mark equal.
- **Polar hydrogens against the deposit's D**, as the share whose parent→H direction is
  turned more than 60° from the deposit's (drift of the parent atom cannot cause that):
  Ser/Thr/Tyr/Cys 39.0 → 37.6% and waters 47.0 → 41.5% (better); Lys NH₃ 0.9 → 8.9% and
  N-termini 3.7 → 7.2% (worse: surface amines swing in implicit solvent, the CE–NZ axis
  tilting, not only the NH₃ turning); ligands 25.1 → 26.4%. Distance alone misleads here:
  the held atoms drift by about 0.25 Å, carrying their hydrogens with them.
- Named cases: 4PDJ, 7D6G, 2ZYE, 5MO0, 5MOQ, 4JEC and 5T8H all minimize (RMSD 0.24–0.28 Å,
  ligand 0.17–0.57 Å, clashes 69–112 → 3–8); 4G0C fails (above).

Prepare's Minimize stays off by default: M5b is still owed, and the 19 defects come first.

## Later phases this design carries (ROADMAP §7a 5–7, stated here so nothing is painted into a corner)

The `SystemSpec` + worker pair is the single entry point for all of it:
- **MD** (phase 5): solvation and ions, restraints from any selection, minimise/heat/equilibrate/
  produce, checkpoints, GPU; trajectories beside the document (needs a narrow `ToolContext`
  method for the document's folder — `Workspace.path` is not exposed and is `None` when unsaved).
- **MM-GBSA** runs in the same worker with GBn2 and needs no AmberTools; **MM-PBSA**, GIST/3D-RISM
  solvent analysis and **openfe** RBFE are the pieces that need the WSL (or Linux) interpreter.
- Also reachable with this stack, worth considering in that order: **GCMC water placement** and
  WaterMap-style hydration-site thermodynamics from GIST for solvent analysis; **ATM**
  (`ATMForce`, OpenMM ≥ 8.2) for absolute binding free energies without atom mapping;
  **metadynamics** for residence time and cryptic pockets; **hybrid ML/MM** ligand potentials
  (MACE-OFF via `openmm-ml`) where Sage is not good enough; and **constant-pH MD** for the sites
  Prepare marks ambiguous.

## Verification

- `pytest tests/simulate tests/macromolecule/prepare tests/tools/macromolecule tests/test_architecture.py
  tests/test_import_boundaries.py`; the full suite once before the release merge.
- `pytest -m bench tests/simulate/test_bench_minimise.py` against the M5 thresholds.
- `ruff check` / `ruff format` on the touched paths.
- By hand: prepare and minimise 5GNK and a DHFR–MTX complex in the workbench — check the report,
  the model in the models bar, undo as one step, save, reopen; then with the worker unconfigured,
  confirm Prepare still runs and says minimisation was skipped.
- `docs/USER_GUIDE.md` (Prepare, the new command, External Tools at :2240 and the env setup
  section), `docs/ROADMAP.md` §7a 3 and 5, `docs/PREPARE_PLAN.md` P5, and a new
  `docs/SIMULATE_PLAN.md` holding this plan, updated in the same commits.
