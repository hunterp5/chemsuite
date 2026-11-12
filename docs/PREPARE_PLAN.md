# Protein preparation plan: coupled protein–ligand protonation (approved 2026-11-11)

## Context

ROADMAP §7a phase 3 says: CCD templates for missing atoms, gemmi for hydrogens, a Reduce-style
network for His/Asn/Gln flips and rotors, and propka at a chosen pH. That plan treats the
protein and the ligand **separately**. Propka shifts protein pKa values. The ligand's state comes
from Standardize/MolGpKa in water. Neither side sees what the other chose. The user wants to
know whether preparation can make the two aware of each other (microenvironment-aware states
for both), and what the most rigorous route is to structures ready for docking and MD.

**Answer: yes.** This is a solved problem in principle and the reference methods exist:
- **Protoss** (Lippert & Rarey 2009; Bietz et al. 2014, ZBH Hamburg) places hydrogens and picks
  tautomers and protonation states of protein, ligand and waters together. It optimises one
  hydrogen-bond network over all of them.
- **Protonate3D** (Labute, *Proteins* 2009, MOE) picks flip, rotamer, tautomer and ionisation
  states of every group (protein, ligand and solvent) under one energy model. The model covers
  vdW, Coulomb, implicit solvation, and titration and tautomer penalties.
- **Schrödinger PrepWizard** does this in steps rather than jointly: Epik for the ligand,
  PROPKA for the protein, then ProtAssign on the H-bond network that includes the ligand.
- **Constant-pH MD** and **alchemical pKa** calculations are the rigorous tier: states sampled
  or computed with full dynamics. They are too expensive to run on everything. They are worth it
  for one critical site, such as a catalytic Asp dyad or a ligand amine next to an Asp.

None of these is open source in the base-install sense, so we build a Qt-free version of
our own, measure it against them, and keep the expensive tiers for the external worker.

## The model: one coupled "site graph"

Every group whose hydrogens are uncertain becomes a **site** with discrete **states**:
- protein: Asp/Glu (charged, or neutral with H on either O, syn), His (HID/HIE/HIP, each
  flipped or not), Lys/Arg/Tyr/Cys (charged or neutral), Asn/Gln (flipped or not), Ser/Thr/Tyr
  OH, Cys SH and Lys NH₃⁺ rotors, the termini;
- **ligand**: each protomer and tautomer from MolGpKa sites plus RDKit tautomer enumeration,
  pruned by a solution free-energy window, with polar-H rotors on each;
- **waters** kept in the pocket: discrete orientations;
- **metals**: coordinating His/Cys/Asp and ligand atoms take states fixed by coordination
  (thiolate Cys, the His N pointing at the metal, a deprotonated sulfonamide on Zn).

The total free energy of one assignment is
`G = Σ_i G_intrinsic(i, s_i, pH) + Σ_i<j W(i, s_i, j, s_j) + Σ_i E_self(i, s_i)`:
- **Intrinsic (solution) term.** For protein sites: `1.364·(pH − pKa_model)` per proton, from
  model pKa values. For ligand sites: MolGpKa's per-site pKa values
  (`chemsuite/ionisation/molgpka.py`, `predict_sites`) plus a tautomer penalty. This is the
  piece every open tool lacks: ligand states start from a real solution-pKa model.
- **Environment (microenvironment) term.** A PROPKA-style empirical desolvation term from
  burial (atom counts in shells) applies to **every** site, protein and ligand alike, so a
  buried ligand amine pays to stay charged exactly as a buried Lys does.
- **Pairwise term W.** Screened Coulomb (distance-dependent dielectric) between charged
  states; a hydrogen-bond score (Reduce/Protoss geometry: distance, angles,
  donor/acceptor complementarity); an H–H and H–heavy clash penalty. Because W runs across the
  protein–ligand interface, each side's state depends on the other's. Example: methotrexate's
  N1 becomes protonated opposite Asp27 in DHFR.
- **Solver.** Build the interaction graph (sites within ~6 Å). Split it into connected
  clusters (Reduce's cliques). Solve each cluster exactly where it is small: dead-end
  elimination, then branch-and-bound. Use Monte Carlo with annealing where it is large.
  Then estimate **Boltzmann populations** per site by MC at 300 K.
- **Calibration.** Fit the few weights (dielectric, desolvation scale, H-bond weight) on
  experimental protein pKa values (PKAD-2). Check the result against PROPKA3 on the same
  complexes; propka stays a dependency as that cross-check and a fallback.

Output per site: chosen state, population, and the evidence (partners, ΔG to the next state).
A site whose second state is within ~1 kcal/mol is marked **⚠ ambiguous**.

## Pipeline (Qt-free, `chemsuite/protein/prepare/`)

1. **Inspect**: altlocs (take the highest occupancy, or keep both as models); chain breaks;
   missing atoms and residues; non-standard residues; covalent links; metals; disulfides. Each
   finding becomes a report row.
2. **Repair**: missing side-chain atoms and termini/caps from CCD templates
   (`structure/components.py`, `CcdCache`); short loops closed; long gaps capped (ACE/NME)
   and marked ⚠. Ligand chemistry from `chem/ligands.py` (bond orders from the file or the CCD;
   `GUESSED` ligands block until the user confirms the structure).
3. **Hydrogens**: ideal-geometry heavy-atom-anchored placement. Standard residues use our own
   templates; ligands use RDKit on the ligand mol. Benchmark gemmi `prepare_topology` first, as
   already decided. Its monomer-library dependency for ligands may make our own placement the
   better choice.
4. **Waters**: keep those within X Å of the ligand that are well ordered (B-factor, H-bonds
   to ≥2 partners) as candidate structural waters; others removed by choice.
5. **Coupled state assignment** (the model above) at a chosen pH.
6. **Relax**: H-only optimisation in the site-graph score in-process. Restrained heavy-atom
   minimisation (amber14 + OpenFF for the ligand, harmonic restraints, ≤0.3 Å RMSD) runs in the
   OpenMM worker when it is configured (phase 5's worker; skipped with a note otherwise).
7. **Report**: every decision with its evidence; flips; ambiguous sites; buried unsatisfied
   donors/acceptors; clashes; net charge. Shown in the workbench as rows and drawn in the view
   (H-bond partners of the selected site).

## Docking vs MD: what "best" means for each

- **Docking**: the receptor in its most probable state. **Ambiguous sites expand into an
  ensemble** (for example both His tautomers), giving receptor variants for ensemble docking
  (phase 2). Ligands to dock keep their solution protomers/tautomers, each carrying its
  state penalty (`G_intrinsic`), which is added to the docking score. That is the Epik approach.
  Optional: after docking, re-run the coupled model on each pose so the pose and its states
  agree.
- **MD**: fixed states from the coupled model. Run one simulation per state of each ambiguous
  site, or constant-pH MD in the worker for those sites. OpenMM has no mature built-in
  constant-pH class yet (Eastman's `openmm-cph` is an initial implementation, and `protons`
  does not support ligands), so this stays a later add-on. The ligand handed to
  parameterisation is the exact protonated and tautomeric form chosen in the complex, not the
  solution form.
- **FEP (phase 6)**: protonation changes between congeneric ligands are flagged. A ligand
  whose bound state differs from its solution state gets the `1.364·(pH − pKa)` correction
  term recorded with the result.

## Benchmark first (as with the conformer plan); `-m bench`, not the regular suite

- **Hydrogen ground truth: neutron and joint X/N structures** from the PDB, where H/D positions
  are observed. Score His tautomer accuracy, Asn/Gln/His flip accuracy, ligand protonation and
  tautomer accuracy, and polar-H placement. Canonical cases to include as named regression
  tests:
  - DHFR–methotrexate (ligand N1 protonated);
  - carbonic anhydrase II with a sulfonamide (anion on Zn);
  - HIV protease and endothiapepsin/BACE1 Asp dyads with inhibitors (one Asp protonated);
  - a kinase hinge-binder tautomer;
  - a zinc finger (thiolate Cys).
- **pKa**: PKAD-2 experimental protein pKa values (RMSE vs PROPKA3's ≈0.8–1.0).
- **Comparisons** (decided 2026-11-11):
  - local: Reduce (flips), PROPKA3 (pKa), and PDBFixer + PDB2PQR in a separate benchmark venv
    (not a dependency);
  - Protoss outputs from the ProteinsPlus web service, sending **public PDB IDs only**, cached
    under `samples/`;
  - MOE Protonate3D outputs, produced by the user by hand from a script and file list we
    supply, kept as a reference column. Optional (2026-11-12): hard to get, so not awaited.

  Licences of every dataset and tool are checked before use.
- **Acceptance**: it becomes the default only if it matches or beats Reduce on flips and His
  tautomers and PROPKA3 on pKa RMSE, and gets the named ligand cases right. It should also
  reach agreement with Protoss (and Protonate3D, if run) on ligand states in the neutron set;
  where they disagree, the neutron data decide.

## Code layout

- `chemsuite/protein/prepare/` (Qt-free, imports `structure`, `chem`, `ionisation`; never Qt or
  `core.table`):
  - `inspect.py`;
  - `repair.py`;
  - `hydrogens.py`;
  - `sites.py`: site and state enumeration, protein templates in `data/`;
  - `ligand_sites.py`: MolGpKa and tautomers to states;
  - `energy.py`: intrinsic, desolvation, Coulomb, H-bond, clash;
  - `solve.py`: clusters, DEE, branch-and-bound, MC populations;
  - `relax.py`;
  - `report.py`.

  Each file stays under 1,000 lines.
- `structure/edit.py`: add atoms and hydrogens as one undoable edit (`ProteinSession.edit`), so
  Prepare is a single undo step. Phase 4's "protonation states and flips by hand" reuses
  `sites.py` (set a site's state, re-solve its neighbours).
- `protein/ui`: a Prepare dialog on the kit (`ToolDialog`; Input: protein and ligands; Method:
  pH, waters, caps, relax; More Options folded), run as a background job that keeps nothing
  partial (it is a whole-structure result). The report goes in a panel registered as a `Part`.
- Prepared state is saved with the document (workspace part).
- `tools/protein`: once prepared, a ligand row's structure can be updated to the bound state
  chosen.
- Tests: `tests/protein/prepare/` (unit tests per module plus the named cases, using small local
  fixtures), benchmark `tests/protein/test_bench_prepare.py`.

## Phases (each usable alone)

- **P0**: benchmark harness, datasets and licences, and baseline numbers for Reduce, PROPKA3
  and gemmi H placement.
- **P1**: inspect, repair and hydrogens. Protein-only flips and rotors (Reduce parity).
- **P2**: the coupled site graph with protein titration (PROPKA parity on PKAD).
- **P3**: ligand sites, waters and metals in the graph. The named cases pass. Populations and ⚠.
- **P4**: UI, report and docking export (ensemble for ambiguous sites, Meeko PDBQT).
- **P5**: restrained relax (now ROADMAP §7a 3, "Relax after preparation": optional, an in-app
  MMFF pocket tier and an OpenMM tier, planned on its own) and MD hand-off through the OpenMM
  worker. Per-site alternate states for MD.

## Verification

- `pytest tests/protein/prepare tests/test_architecture.py tests/test_import_boundaries.py`.
- `pytest -m bench tests/protein/test_bench_prepare.py`: report vs the acceptance thresholds.
- By hand: prepare 5GNK and a DHFR–MTX complex in the workbench. Check the report, H-bond
  partners drawn, undo as one step, save, reopen.
- `docs/USER_GUIDE.md` updated with each phase, and ROADMAP §7a phase 3 rewritten to this plan.
