# Conformer benchmark: crystal poses

How closely each generator reproduces the bioactive pose of 3,354 drug-like PDB ligands
(Wang et al.'s Dataset I: J. Chem. Inf. Model. 2023, 63, 6525–6536,
doi:10.1021/acs.jcim.3c01519). Each generator starts from the ligand's SMILES and makes at most
1, 10, 50 or 250 conformers; a ligand's score is the heavy-atom RMSD, with symmetry, of its
closest conformer. Made with `scripts/review/conformer_bench.py score` (its docstring gives
the steps) on 2026-11-15: Windows PC, 8 worker processes, RDKit 2026.09.2, CDPKit 1.3.0.

| Variant | Generator | Optimised | Energy window | RMS pruning |
|---|---|---|---|---|
| etkdg | ETKDGv3 (the app's default) | MMFF94 | 10 kcal/mol | 0.5 Å |
| etkdg-raw | ETKDGv3 | no | none | none |
| etkdg-prune | ETKDGv3 | no | none | 0.5 Å |
| etkdg-mmff | ETKDGv3 | MMFF94 | none | 0.5 Å |
| conforge | CONFORGE | (its own) | 10 kcal/mol | 0.5 Å |
| conforge-wide | CONFORGE | (its own) | none | 0.5 Å |

The rows marked "(Wang et al.)" are their Table 1.

## Findings

- **CONFORGE is the most accurate generator here** by AMR-R and by coverage within 0.5 and
  1.0 Å, at every count. At 250 conformers conforge-wide (AMR-R 0.582 Å, 49.9% within 0.5 Å)
  is ahead of OMEGA (0.591 Å) and ConfGenX (0.597 Å) in Wang et al.'s table; only torsional
  diffusion (0.574 Å) is lower. Within 2.0 Å every generator covers 98–100% at 50 or more.
- **The app's ETKDG default is the weakest ETKDG variant at 50 and 250 conformers.** Its
  energy window throws most conformers away (31 of 250 kept on average), and MMFF clean-up does
  not bring poses closer to the crystal. ETKDG with neither (etkdg-raw) does best: 0.747 /
  0.623 Å at 50 / 250, against 0.808 / 0.725 Å for the default. This agrees with Wang et al.
- **etkdg-raw reproduces Wang et al.'s RDKit ETKDG figures** (their 1.599 / 1.003 / 0.781 /
  0.648 Å), slightly better at every count.
- **Failures:** ETKDG could not embed 4 ligands (6–8 at 250 conformers); CONFORGE failed on 2.
- **Timing is only roughly comparable.** The ETKDG and CONFORGE runs were made at different
  times on a shared machine. CONFORGE builds its fragment conformers once per ligand and reuses
  them for the later counts, so its 10-, 50- and 250-conformer times are understated: the
  1-conformer time (mean 0.3 s) is closer to what one molecule costs in the app.

## Results

Ligands scored by every variant: 3350 of 3354

AMR-R (Å, mean closest heavy-atom RMSD to the crystal pose; lower is better)

| method | 1 | 10 | 50 | 250 |
|---|---|---|---|---|
| etkdg | 1.544 | 0.985 | 0.808 | 0.725 |
| etkdg-raw | 1.547 | 0.962 | 0.747 | 0.623 |
| etkdg-prune | 1.547 | 0.988 | 0.796 | 0.692 |
| etkdg-mmff | 1.544 | 0.959 | 0.773 | 0.688 |
| conforge | 1.432 | 0.837 | 0.677 | 0.622 |
| conforge-wide | 1.430 | 0.852 | 0.661 | 0.582 |
| ConfGenX (Wang et al.) | 1.510 | 0.952 | 0.714 | 0.597 |
| Conformator (Wang et al.) | 1.549 | 0.969 | 0.783 | 0.669 |
| OMEGA (Wang et al.) | 1.480 | 0.965 | 0.739 | 0.591 |
| RDKit ETKDG (Wang et al.) | 1.599 | 1.003 | 0.781 | 0.648 |
| torsional diffusion (Wang et al.) | 1.742 | 0.999 | 0.731 | 0.574 |

COV-R at 0.5 Å (% of ligands with a conformer within it)

| method | 1 | 10 | 50 | 250 |
|---|---|---|---|---|
| etkdg | 11.0 | 23.5 | 30.6 | 34.0 |
| etkdg-raw | 9.6 | 26.2 | 36.9 | 45.6 |
| etkdg-prune | 9.6 | 22.6 | 30.0 | 35.5 |
| etkdg-mmff | 11.0 | 23.9 | 31.3 | 34.8 |
| conforge | 13.4 | 34.0 | 43.6 | 46.7 |
| conforge-wide | 13.4 | 31.9 | 43.5 | 49.9 |

COV-R at 1.0 Å (% of ligands with a conformer within it)

| method | 1 | 10 | 50 | 250 |
|---|---|---|---|---|
| etkdg | 29.5 | 57.4 | 70.8 | 77.7 |
| etkdg-raw | 27.8 | 57.6 | 73.5 | 82.9 |
| etkdg-prune | 27.8 | 56.8 | 72.2 | 81.5 |
| etkdg-mmff | 29.5 | 58.6 | 73.1 | 81.1 |
| conforge | 35.7 | 67.3 | 80.2 | 85.0 |
| conforge-wide | 35.6 | 66.0 | 81.0 | 87.7 |

COV-R at 2.0 Å (% of ligands with a conformer within it)

| method | 1 | 10 | 50 | 250 |
|---|---|---|---|---|
| etkdg | 72.1 | 93.2 | 97.3 | 99.0 |
| etkdg-raw | 72.1 | 93.6 | 97.8 | 99.3 |
| etkdg-prune | 72.1 | 93.6 | 97.8 | 99.3 |
| etkdg-mmff | 72.1 | 94.5 | 98.2 | 99.6 |
| conforge | 76.2 | 95.6 | 98.0 | 98.7 |
| conforge-wide | 76.3 | 95.5 | 98.5 | 99.4 |

Seconds per ligand (median / mean), conformers made (mean), failures

| method | count | median s | mean s | made | failed |
|---|---|---|---|---|---|
| etkdg | 1 | 0.07 | 0.10 | 1.0 | 4 |
| etkdg | 10 | 0.50 | 0.79 | 5.5 | 4 |
| etkdg | 50 | 1.98 | 3.29 | 14.1 | 4 |
| etkdg | 250 | 7.66 | 14.23 | 30.9 | 6 |
| etkdg-raw | 1 | 0.03 | 0.07 | 1.0 | 4 |
| etkdg-raw | 10 | 0.26 | 0.47 | 10.0 | 4 |
| etkdg-raw | 50 | 1.28 | 1.96 | 50.0 | 4 |
| etkdg-raw | 250 | 6.30 | 9.38 | 250.0 | 8 |
| etkdg-prune | 1 | 0.03 | 0.06 | 1.0 | 4 |
| etkdg-prune | 10 | 0.24 | 0.41 | 7.1 | 4 |
| etkdg-prune | 50 | 1.18 | 1.73 | 24.4 | 4 |
| etkdg-prune | 250 | 5.87 | 8.46 | 74.2 | 6 |
| etkdg-mmff | 1 | 0.07 | 0.10 | 1.0 | 4 |
| etkdg-mmff | 10 | 0.51 | 0.80 | 6.2 | 4 |
| etkdg-mmff | 50 | 1.99 | 3.34 | 18.3 | 4 |
| etkdg-mmff | 250 | 7.75 | 14.85 | 48.0 | 6 |
| conforge | 1 | 0.02 | 0.29 | 1.0 | 2 |
| conforge | 10 | 0.01 | 0.19 | 8.2 | 2 |
| conforge | 50 | 0.02 | 0.20 | 27.6 | 2 |
| conforge | 250 | 0.02 | 0.32 | 72.7 | 2 |
| conforge-wide | 1 | 0.02 | 0.27 | 1.0 | 2 |
| conforge-wide | 10 | 0.02 | 0.17 | 8.9 | 2 |
| conforge-wide | 50 | 0.02 | 0.18 | 35.5 | 2 |
| conforge-wide | 250 | 0.05 | 0.36 | 125.7 | 2 |
