# Protein preparation: existing pKa tools, licenses first (2026-11-19)

Follows [`p5-review.md`](p5-review.md) §5: before building Poisson–Boltzmann or ML titration,
run existing implementations through `scripts/review/prepare_states.py` as benchmark-only tools
(the way PROPKA and Protoss are). This page is the survey: what exists, its license, what it
outputs, and the order to try them in. No code yet.

## 1. The bar, stated first

Every candidate below predicts **per-site pKa values** (PypKa and MCCE also give
populations). None decides His tautomer or flip. Those stay with P1's hydrogen-bond optimizer,
which is where P1 already wins (tautomer 63.3%, orientation 87.6%). So a tool can only move the
states benchmark through `pipeline.protonated`'s charge decision: Asp/Glu protonated, Lys
neutral, Tyr/Cys anions, His charged.

- **That is a narrow target.** On PKAD-R at pH 7, 93% of residues are in the usual form. PROPKA
  3.5.1 caught 2 of 21 unusual states and made no wrong switches (p5 §1).
- **PROPKA has never been scored on the neutron set's states.** If PROPKA, at about 1.1 RMSE
  on PKAD, does not beat P1 there, a tool at 0.8–1.0 RMSE needs to be much better on the
  *unusual* states specifically, not on average.
- **What makes a candidate worth its cost:** recall of the unusual states without wrong
  switches. KaML is the only one that reports this directly (§3).

## 2. The candidates

| Tool | Method | License | Use here | Runs on | Outputs |
|---|---|---|---|---|---|
| **KaML-CBtree** (Shen lab, [KaMLs](https://github.com/JanaShenLab/KaMLs)) | CatBoost trees on 37 structure features | MIT, models included | **Benchmark; could ship** (trees in numpy) | Linux (WSL); CPU; DSSP and RIDA binaries | pKa per Asp, Glu, His, Cys, Tyr, Lys |
| **pKAI+** (Bayer, [pKAI](https://github.com/bayer-science-for-a-better-life/pKAI)) | MLP on N/O/S atoms within 15 Å | MIT, TorchScript weights in the repo | **Benchmark; could ship** (MLP in numpy) | Any; torch, CPU | pKa per Asp, Glu, His, Lys, Cys, Tyr, termini |
| **MCCE4** with **NextGenPB** ([MCCE4-Alpha](https://github.com/GunnerLab/MCCE4-Alpha), [NextGenPB](https://github.com/concept-lab/NextGenPB)) | PB + Monte Carlo over conformers | MIT + GPL-3.0 | Benchmark only (separate interpreter, like OpenMM) | Linux/WSL via Apptainer; Windows not supported natively | pKa, occupancy of each conformer (charge and tautomer) at a pH |
| **PypKa** ([PypKa](https://github.com/mms-fcul/PypKa)) | PB + Monte Carlo, with tautomers | LGPL-3.0, **but needs DelPhi** (proprietary; free to registered non-profit labs) | Benchmark only, **if** DelPhi's terms allow it | Linux | pKa, tautomer populations, protonated PDB at a pH |
| **openmm-cph** ([openmm-cph](https://github.com/openmm/openmm-cph)) | Constant-pH MD (explicit solvent, implicit-solvent MC moves) | MIT | Later, in the OpenMM worker (p5 §5.3) | Wherever OpenMM runs | States sampled by simulation |
| KaML-ESM ([kaml-esm](https://github.com/JanaShenLab/kaml-esm)) | ESM embeddings + MLP | **CC BY-NC 4.0**; ESM-C through a paid API | Excluded | 24 GB GPU, 128 GB RAM | pKa |
| pKa-ANI ([pKa-ANI](https://github.com/isayevlab/pKa-ANI)) | ANI representations | **Non-commercial, no sharing; the licensor owns modifications** | Excluded, even as a benchmark | — | — |
| PDB2PKA | PB (APBS) | BSD | Excluded: removed from PDB2PQR 3.x, which is itself dropped | — | — |
| ME-pKa ([ME-pKa](https://github.com/yzjyg215/ME-pKa)) | ESM + structure, multimodal | MIT | Not now: no published weights; needs CHARMM and Chimera to featurize | CUDA | — |

Notes:
- **openmm-cph** is still prerelease ("not yet well validated", no releases, last commit
  2026-08-14) and not in OpenMM itself. It supports `ForceField`-built systems only, which
  matches the worker. It stays the MD-side option and is not a preparation step.
- **H++, Karlsberg+, DeepKa and the DelPhiPKa web server** are web services. Sending structures
  to them is publishing them and they cannot be scripted at benchmark scale, so they are left
  out.

## 3. What each claims

- **KaML-CBtree** (Shen & Shen, 2024–2025, [PMC11601431](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601431/)),
  trained on **PKAD-3** (1,167 experimental values, 247 proteins):
  - RMSE: acids 0.76 against PROPKA3's 1.28, bases 0.79 against 0.95. Per type: Asp 0.75,
    Glu 0.60, His 0.85, Cys 1.50, Tyr 1.24, Lys 0.70.
  - **States:** protonated Asp recall 0.77 (PROPKA3 0.28); protonated Glu precision 1.0
    (0.67); Cys thiolate precision/recall 0.73/0.76 and neutral Lys 1.0/0.8 (PROPKA3 finds
    neither).
  - **Caveat:** the split is random by residue, not by protein, so these numbers are
    optimistic. PKAD-3 contains PKAD-R's residues, so PKAD-R cannot test it. The neutron set
    is independent of its training data and is the fair test.
  - Input must be a single chain without missing atoms. Our prepared models are complete, so
    the bench can feed it one chain at a time (losing inter-chain contacts).
- **pKAI+** (Reis et al., JCTC 2022, [PMC9369009](https://pmc.ncbi.nlm.nih.gov/articles/PMC9369009)):
  - Trained on about 3 million PypKa-computed pKa values, not on experiment; pKAI+ adds a
    penalty on large shifts. RMSE 0.98 on 736 PKAD residues (pKAI 1.15).
  - Heavy atoms only (N, O, S within 15 Å), no hydrogens, 4,008 inputs, three hidden layers.
  - Its shrinkage toward the model pKa is what made it accurate on average, and it is also
    what makes it unlikely to call unusual states. The bench will show which wins.
- **MCCE4 and PypKa** are the "Poisson–Boltzmann titration" family that p5 §5.1 names. They
  replace exactly the Coulomb and desolvation terms P2/P5 found weak. Neither publishes run times on
  the pages read; PB with Monte Carlo is expected to take minutes per protein (MCCE with
  side-chain conformers longer), to be measured. That is fine for a benchmark of about 100
  structures run once, but not interactive.

## 4. Order of work (proposed)

1. **One hook in the harness.** Add a method to `prepare_states.py` that takes per-site pKa
   values from a file (`<ID>.csv`: chain, residue, name, pKa) and scores them through the same
   decision as P1's titration (`pipeline.protonated` at the crystal's pH, tautomer and flip
   from P1). One hook serves every tool below.
2. **PROPKA through it first.** It is already a bench dependency, and it sets the bar.
3. **KaML-CBtree and pKAI+** in one WSL bench venv (catboost, pycaret, torch CPU, DSSP, RIDA;
   nothing reaches the app's environment).
4. **MCCE4 + NextGenPB** in WSL through its Apptainer setup, if step 3 shows that pKa quality
   moves the states score at all.
5. **PypKa** only if you decide DelPhi's registration covers this use (your call: it is free
   to non-profit labs that register with Columbia).

If a tool beats P1 on the test half, shipping it is a separate plan:
- **KaML-CBtree:** CatBoost exports oblivious trees to JSON, which evaluate in numpy in a few
  lines. The work is reimplementing its 37 features Qt-free (DSSP-like secondary structure and
  the disorder flag included). MIT models can ship with a NOTICE, as chemsuite pKa's do.
- **pKAI+:** an MLP, so a numpy port is direct (the chemsuite pKa route; never torch in the
  app). Its featurization is simple.

Data to check before any training of our own (not proposed now): PKAD-3's download terms. Its
site did not load on 2026-11-19, and the paper names no data license.

## 5. Steps 1 and 2 done: the hook, and PROPKA as the bar (2026-11-19)

**The hook.** `prepare_states.py score --tools a,b` adds a column per tool, read from
`samples/prepare/out/pka/<tool>/<ID>.csv` (chain, res_seq, icode, group, pka; group a key of
`sites.MODEL_PKA`, so `NTERM` for an N-terminus).
- Each value gives its group's protonated share at the crystal pH (Henderson–Hasselbalch).
  `pipeline.protonated` turns the share into a form, as titration does.
- A group the file leaves out keeps its form in water at that pH.
- The titrated tables are then solved with those forms fixed, so His tautomer and flip are
  chosen as titration on chooses them.
- Only entries every listed tool covers are scored.
- The cases now carry each site's key (`Case.keys`), so `build --redo` was needed once.
- `prepare_states.py propka` turns `prepare_review.py propka`'s `.pka` files into these CSVs
  (262 entries).

A tool with no values does not reproduce the null column, by design. The null holds every
group in its pH 7 form, while the hook (like the pipeline) uses the form at the crystal pH: at
pH 5.5 every His is charged. Checked on all 15,243 scored sites: each lands in its form in
water at the crystal pH, except the 150 His in entries at pH 6.5 (His's model pKa), where the
tie goes to the charged form, as it does in the pipeline.

**PROPKA 3.5.1 on the test half** (124 entries, default weights, `score --test --tools propka`):

| Question | Null | P1 | Titration | PROPKA |
|---|---|---|---|---|
| Asp state | 99.3% | 99.3% | 98.4% | 95.7% |
| Asp unusual found | 0/13 | 0/13 | 1/13 | 4/13 |
| Asp wrong switches | 0 | 0 | 19 | 73 |
| Glu state | 98.2% | 98.2% | 94.3% | 91.5% |
| Glu unusual found | 0/34 | 0/34 | 17/34 | 14/34 |
| Glu wrong switches | 0 | 0 | 90 | 140 |
| Lys, Tyr, Cys unusual found | 0/40 | 0/40 | 0/40 | 0/40 |
| His charge | 68.7% | **72.6%** | 63.5% | 66.7% |
| His charged found | 0/261 | 72/261 | 169/261 | 62/261 |
| His wrongly charged | 0 | 39 | 212 | 78 |
| His tautomer | 43.2% | **63.3%** | 58.0% | 57.7% |
| His orientation | 50.3% | **87.6%** | 84.8% | 86.4% |
| Objective | 0.817 | **0.893** | 0.866 | 0.864 |

The training half agrees (P1 0.901, PROPKA 0.876).

- **PROPKA loses to P1, and slightly to our own titration.** It finds about a third of the
  unusual Asp and Glu, at about 3 wrong switches for each find.
- **Its His calls are worse than P1's on every His question.** P1 charges a His only where
  its neighbors make that the best state. PROPKA charges fewer correct ones and more wrong
  ones, and each wrong HIP also costs a tautomer.
- **The bar for the next tools is therefore P1's 0.893, not PROPKA's.** For a pKa tool to pass
  it, wrong switches must stay below about one per find. KaML-CBtree's reported precision
  (§3) is the only claim in reach of that. Lys, Tyr and Cys anions were not found by anything
  so far (40 on the test half).
- A margin (a switch only when the pKa clears the pH by some units) is the obvious knob for any
  tool. `SWITCH` already takes a threshold per group, so it can be tried without new code by
  passing thresholds to `decide_external`.

## 6. Step 3: pKAI+, and the PROPKA bar corrected (2026-11-24)

**The runner.** `scripts/review/pka_tools.py` runs in WSL in its own conda env (`pka-bench`:
Python 3.11, pKAI 1.2.0 with torch 1.13.1, numpy < 2, pycaret 3.3.2, catboost 1.2.10,
biopandas). It imports nothing from chemsuite. Every tool gets the input PROPKA got
(`samples/prepare/inputs/<ID>.pdb`), cleaned the same way:
- ATOM records of the 20 standard amino acids only, with no hydrogens and the first
  conformation. Neither tool reads HETATM records by design, so ligands, metals and waters are
  invisible to both.
- Each chain is renumbered 1, 2, … and the values are mapped back. Without this, pKAI drops
  every atom with an insertion code (16 entries have them) and KaML reads numbers as integers.
- Half of the inputs' Asn/Gln/His side chains are flipped at random (as for every method).

**Coverage, checked before scoring** (`prepare_states.py coverage --tools …`). A value whose
key matches no site falls back to the water form without a word, so this is checked first.
- **pKAI+** covers 100% of the Asp, Glu, His, Cys, Tyr and Lys sites. It gives no N-terminus
  values; pKAI 1.2.0 skips termini (a TODO in its code), so they keep their water form.
- The 8 entries without protein (nucleic acids alone, e.g. 1V9G) get empty files.
- The runner loads `pKAI+_model.pt` explicitly; pKAI's own default is the plain pKAI model.

**The PROPKA bar of §5 was off.** Its CSVs dropped about 5% of sites, and those sites were
scored in their water forms:
- PROPKA glues numbers of 1000 and up to the group (`ASP1008 B`). The whitespace split missed
  them.
- PROPKA prints no insertion codes, so `TYR 184` (really 184A in 1NTP) never matched.

`propka_csv` now reads the summary with a regex. A number takes the insertion code of the one
residue of that name it can be. Coverage went from about 95% to 99.6%; what is left is
PROPKA's own merging of residues that share a number (3VXF's light chain). The old CSVs are
kept as `out/pka/propka.before-fix`.

**Scores** (`score --tools propka.before-fix,propka,pkai`; every entry is covered, 138
training, 124 test):

| Test half | P1 | Titration | PROPKA (§5) | PROPKA (fixed) | pKAI+ |
|---|---|---|---|---|---|
| Asp unusual found | 0/13 | 1/13 | 4/13 | 4/13 | 0/13 |
| Asp wrong switches | 0 | 19 | 73 | 76 | 11 |
| Glu unusual found | 0/34 | 17/34 | 14/34 | 17/34 | 1/34 |
| Glu wrong switches | 0 | 90 | 140 | 143 | 9 |
| Lys, Tyr, Cys unusual found | 0/40 | 0/40 | 0/40 | 0/40 | 0/40 |
| His charged found | 72/261 | 169/261 | 62/261 | 61/261 | 115/261 |
| His wrongly charged | 39 | 212 | 78 | 70 | 142 |
| His tautomer | **63.3%** | 58.0% | 57.7% | 58.6% | 57.7% |
| His orientation | **87.6%** | 84.8% | 86.4% | 86.6% | 85.8% |
| Objective | **0.893** | 0.866 | 0.864 | 0.866 | 0.874 |

Training half: P1 0.901, titration 0.855, PROPKA 0.876 → 0.883 after the fix, pKAI+ 0.880.

- **Nothing passes P1.** pKAI+ is the best tool on the test half (0.874) but loses to P1 on
  both halves.
- **The fix moves PROPKA a little** (test 0.864 → 0.866, training 0.876 → 0.883). It finds
  three more unusual Glu and is still at about 3 to 4 wrong switches per find. The conclusion
  of §5 stands.
- **pKAI+ is what its design predicts** (§3). Its shrinkage toward the model pKa leaves Asp
  and Glu almost always in the usual form: 1 find against 20 wrong switches on the test half.
  It does not reach the buried acids.
- **Its His calls cost more than they gain.** pKAI+ charges 115 of 261 truly charged His, but
  also 142 neutral ones; P1 charges 72 and 39. Every wrong HIP also costs a tautomer. A margin
  on His (§5's `SWITCH` knob, chosen on the training half) is the one setting likely to help it.
- Scoring all four built-in methods with three tools took 3.5 h, almost all of it the full
  `titrated` solve. Pass `--fast` when only the external columns are new (fast matches
  titrated to within 0.1 point here).

## 7. Step 3 continued: KaML-CBtree, and margins (2026-11-24)

**Running it.** KaML's repo ships DSSP and RIDA only as prebuilt Linux binaries. The user chose
to run them in WSL; the 32-bit `mkdssp` needs `libc6:i386`. `pka_tools.py kaml` runs each chain
alone, since KaML takes single chains, so inter-chain contacts are lost.
- It loads the two finalized models once per worker and runs DSSP once per chain. KaML runs
  it again for every residue, and the output does not change.
- A residue KaML cannot compute features for (a missing side-chain atom gives an error row with
  one field too many; 42 entries have one) gets no value and keeps its form in water. Blank
  numeric features go to pycaret's imputer as NaN, as they would in training.
- Coverage: 100% of the Asp, Glu, His, Cys, Tyr and Lys sites; no N-termini (KaML has none).

**A bug in KaML's published feature code.** The first run scored 0.774 on the test half, below
the null method (0.817). Features regenerated for 20 of KaML's own training structures (1KSM,
4LZT, …, fetched from the PDB) were compared with its training table, column by column. All
matched closely except the four backbone hydrogen-bond distances (DA1, DA2, DD1, DD2), which
were about 14 Å off: `features.bbh_dis` computes `x1 - y1` where it means `x1 - x2`. The models
were trained on correct distances, so the published prediction script feeds them wrong ones.
- `pka_tools.py` replaces `bbh_dis` with the correct distance; KaML's clone is untouched.
- On those training proteins the fix roughly halves the error (MAE His 0.53 → 0.23,
  Lys 0.86 → 0.48, Asp 0.36 → 0.18).
- A downward bias remains on Lys (−0.47) and Tyr (−0.63). It goes away when each row gets its
  true `group_code`, a categorical feature that the published inference always sets to 0, an
  unseen category. That is how KaML ships, so it is left as is.
- Casting `charge` to a float, as the training table has it, changes nothing.
- The first run's CSVs are kept as `out/pka/kaml.published-distance`.

**Leakage.** 43 of the 262 entries share a UniProt accession with KaML's training rows (mostly
carbonic anhydrase II, 12; lysozyme, 9; RNase A, 4), 3 share a PDB ID, and 24 have no UniProt
DBREF to check. The lists are `out/pka/kaml_overlap.txt` and `kaml_overlap_or_unknown.txt`
for `score --skip`.

**Scores** (`score --fast --test --tools kaml`; objective):

| Test half | Entries | Null | P1 | PROPKA | pKAI+ | KaML |
|---|---|---|---|---|---|---|
| All | 124 | 0.817 | **0.893** | 0.866 | 0.874 | 0.779 |
| Without overlap | 103 | 0.827 | **0.907** | 0.875 | 0.880 | 0.782 |
| Without overlap or unknown | 88 | 0.830 | **0.908** | 0.877 | 0.883 | 0.789 |

Training half: KaML 0.793 (P1 0.901). The leakage does not change the picture.

- **KaML loses on Cys and Tyr.** It switches 84 of 165 usual Cys to thiolate and 177 of 1,141
  usual Tyr to the anion (test half). Its training labels explain the Cys calls: their median
  is 6.3, because the Cys pKas people measure are mostly perturbed active-site ones. KaML's
  predicted medians on these inputs are Cys 6.9, Tyr 8.5, Lys 9.1.
- **On Asp and Glu it is more careful than PROPKA:** Glu 12/34 found with 52 wrong switches
  (PROPKA 17/34, 143); Asp 0/13 with 34. That is still about 4 wrong switches for each find.
- **His:** 105/261 charged found, 178 wrongly charged, tautomer 51.6% (P1 63.3%).

**Margins** (`margin --tools X --jobs 10`). Each group's threshold is searched over the shares
0.5–0.99 on the training half, then the test half is scored once. 90 min for both tools.

| Test half | KaML default | KaML margin | pKAI+ default | pKAI+ margin |
|---|---|---|---|---|
| Thresholds | 0.5 | Cys, Glu, Tyr 0.99; Asp 0.95; His 0.8 | 0.5 | His 0.6 |
| Cys wrong switches | 84 | 6 | 0 | 0 |
| Tyr wrong switches | 177 | 33 | 0 | 0 |
| Glu found / wrong | 12 / 52 | 4 / 27 | 1 / 9 | 1 / 9 |
| His charged found / wrongly | 105 / 178 | 79 / 147 | 115 / 142 | 93 / 133 |
| Objective | 0.779 | 0.851 | **0.874** | 0.869 |

- Margins take most of KaML's damage away, but what is left is still below pKAI+ and P1.
- pKAI+'s His margin gains 0.0003 on the training half and loses 0.005 on the test half: noise.
  Its defaults stand.

## 8. Where this leaves titration (2026-11-24)

**No external pKa tool beats P1**: not PROPKA, pKAI+ or KaML-CBtree, with or without margins,
and not on the entries free of KaML's training proteins. The test-half objectives are P1
0.893, pKAI+ 0.874, PROPKA 0.866, KaML 0.851 with margins.

- **Every tool loses on His.** P1 charges a His only where its hydrogen-bond network makes HIP
  the best state, and so holds the tautomer and orientation scores (63.3%, 87.6%). Every pKa
  tool charges more His wrongly than P1 does, and each wrong HIP also costs a tautomer.
- **The unusual acid states are where the tools could help**, and none finds them cleanly. The
  best ratio is pKAI+'s, which finds almost none. PROPKA and KaML each make about 3–4 wrong
  switches for each find.
- **Lys, Tyr and Cys:** 40 unusual states on the test half. Only KaML finds any (8, at 268
  wrong switches), and its margins take them away again.

Recommendation: keep titration off and ship none of these tools. Step 4 (MCCE4 + NextGenPB) was
conditional on step 3 showing that pKa quality moves the states score. It doesn't for these
three, so it is not worth starting unless the user wants Poisson–Boltzmann checked for its own
sake. Constant-pH MD (openmm-cph) stays the later, MD-side option.
