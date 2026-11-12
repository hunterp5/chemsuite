# Protein preparation P1: hydrogens, flips, His states and rotors (2026-11-12)

Plan: [`docs/PREPARE_PLAN.md`](../PREPARE_PLAN.md); baselines and pass marks:
[`p0-review.md`](p0-review.md). Code: `chemsuite/protein/prepare/` (`templates`, `hydrogens`,
`sites`, `energy`, `solve`, `inspect`, `repair`, `pipeline`). P1 is done.

## 1. What runs

`pipeline.prepare(structure, Options(ph, waters, repair, optimise))` runs the steps below in
order and returns the structure, the His and Cys states, the site choices and a report
(`inspect.Finding` rows, ⚠ where the user should look). Old atoms keep their uids, so the
workbench can apply the result as one undoable edit (P4).

0. **Inspect and repair** (`inspect.inspect`, `repair`). The report lists:
   - alternate locations;
   - chain breaks, and residues missing from the deposited sequence (⚠; no loops are built,
     the ends stay charged termini);
   - missing atoms, and residues with no template;
   - covalent links, disulfides and metal sites.

   Repair keeps the location shown and drops the others. It builds missing side-chain atoms and
   a chain's last OXT atom by atom from the template's internal coordinates. A rebuilt side chain
   takes the CCD's ideal conformation, which is not checked for clashes. Waters are kept, kept
   within 5 Å of a ligand, or removed.
1. **Hydrogens** (`hydrogens.add_hydrogens`). Polymer residues are built from CCD templates in
   their textbook states (Asp/Glu/Lys/Arg charged, His HIE, NH₃⁺ and carboxylate termini).
   - An atom bonded to another residue keeps one hydrogen fewer per bond. Its hydrogens are
     placed from its real neighbours, so the amide H lies in the peptide plane.
   - A soft metal (Zn, Fe, Cu, …) within 2.9 Å takes a proton's place, as in a Zn thiolate.
   - Ligands get hydrogens from RDKit on the component's own molecule.
   - On ubiquitin's NMR hydrogens (1D3Z), the hydrogens that do not turn land 0.03 Å from the
     deposit on average.
2. **Sites** (`sites.find_sites`):
   - Asn/Gln as built or flipped;
   - His as HID/HIE/HIP, each as built or flipped (a His bonded to another residue is fixed);
   - the polar rotors (OH, SH, NH₃⁺, on protein and ligands) in 15° steps, with a phenol or
     acid OH kept in its plane.
3. **Energy** (`energy.tables`): a weighted sum of terms, so the weights could be fitted on
   precomputed tables:
   - H-bond (distance and D–H···A angle);
   - H–H and H–heavy overlap;
   - lone pair on lone pair;
   - contact of a bare acceptor with water;
   - a metal term (a bare N or O on a metal −1, a carbon +1);
   - for HIP, Coulomb with ε = 4r and burial;
   - priors for a flip, HIP, HID and pH.
4. **Solver** (`solve.optimise`): clusters of touching sites, Goldstein dead-end elimination,
   then a branch-and-bound from a greedy start. It checks against brute force in the tests.
   Each site keeps its gap to its next best state, for the ⚠ of P3.

## 2. Fitting

`prepare_review.py fit` holds the 262 entries in two halves by CRC32 of the ID. It runs a
coordinate search over the weights on the training half and reports both halves. The fitted
weights are the defaults in `energy.Weights`: H-bond 1, clash 160, lone 16, water 1, flip 0,
HIP 3, HID 0.8, pH 1, Coulomb 1, burial 1, metal 2. The search ended at the top of the clash
and lone ranges; going further gained nothing on the test half.

How the terms came in, measured on the test half (Asn / Gln / His orientation / His tautomer):

| Step | Test half |
|---|---|
| H-bond, clash, lone, water, priors | 85.4 / 79.7 / 86.1 / 58.2% |
| own backbone scored (only atoms within 3 bonds left out) | 85.9 / 80.8 / 86.1 / 58.6% |
| Coulomb and burial for HIP | 85.9 / 80.6 / 86.2 / 60.2% |
| His on metals as sites, with the metal term | **85.9 / 80.6 / 87.3 / 63.3%** |
| Reduce on the same half | 84 / 79 / 87 / 61% |

The test half is harder for every method: Reduce scores 86 / 79 / 90 / 64% on the training
half.

## 3. Scores on all 262 entries

`prepare_review.py chemsuite --redo`, then `score`. The step runs `pipeline.prepare` with
repair on, writes PDB files and scores them like the other methods (about 3 min for the
set):

| Question | Pass mark | chemsuite | Reduce | Protoss |
|---|---|---|---|---|
| Asn orientation | 85% | **86.2%** | 85% | – |
| Gln orientation | 79% | **80.4%** | 79% | – |
| His orientation | 88% | **89.6%** | 88% | – |
| His tautomer | 62% | **63.6%** | 62% | 59% |
| His charge | 74% (a P2 mark) | 74% | 74% | 73% |
| Asp / Glu state | textbook | 99.0 / 98.1% | 99 / 98% | 98 / 98% |
| Ligand polar-atom H count | 96% | 96.0% | 91% | 94% |
| Entries answered | 262 | 262 | 262 | 259 |

Half of the full set was used to fit, so the test-half row in §2 is the fair comparison. It
matches or beats Reduce there too. `tests/protein/test_bench_prepare.py` holds these marks.

## 4. Found on the way

- The CCD's `Template.with_hydrogen` dropped the charge list (crash on any N-terminus); fixed.
- P0's harness read a hydrogen as bonded within 1.25 Å, which missed every 1.34 Å S–H. The Cys
  rows of p0-review.md were re-scored.
- Nothing scored a site against its own residue, so a His ND1–H···O=C of its own backbone
  was invisible. Only atoms within three bonds are left out now.
- A His on a metal was skipped as "decided". In the inputs half of them are flipped, which
  left a carbon facing the metal; they are now sites that the metal term orients.

## 5. Open

- His charge stays at the textbook level (74%). In the training half 28% of His are HIP, and
  we call 95 of the 216 HIP neutral. P2's titration model (desolvation, the full Coulomb
  network, PKAD-2) is where that improves.
- HID is still read as HIE in 70 of 185 training cases.
- Seven Tyr on metals are deprotonated by the metal rule but are protonated in the neutron data.
- Waters carry no hydrogens yet, so a water counts as donor and acceptor at once. MD needs
  them (P5).
- Repair puts a rebuilt side chain in its ideal conformation, with no rotamer search, and
  caps no gaps (ACE/NME). Both are left for later.
- Each `Finding.residue` indexes the structure it was found on. Inspection runs on the input,
  so the UI should go by the label.
