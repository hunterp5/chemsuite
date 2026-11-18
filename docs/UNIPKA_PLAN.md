# Our own Uni-pKa (plan approved 2026-11-12)

chemsuite trains its own model of Uni-pKa's design and runs it in numpy. It replaces MolGpKa:
first for the ligand states in protein preparation (`docs/PREPARE_PLAN.md` P3), then
everywhere once the benchmark shows it is better. This supersedes `docs/PKA_PLAN.md`'s
optional Uni-pKa add-on.

## Why

- **MolGpKa is not good enough for preparation.** On the neutron set's ligands it catches 10 of
  the 57 atoms protonated unlike the CCD (3 of 33 on the test half), and breaks 12 the drawing
  had right. It predicts each site with the rest of the molecule neutral, so it cannot rank
  whole microstates or tautomers. Methotrexate N1, folate N3 and the PLP pyridiniums are missed.
- **Uni-pKa** (Luo et al., JACS Au 2024, doi:10.1021/jacsau.4c00271) gives a free energy for
  each whole microstate.
- **Its published weights cannot be used:**
  - the official ones come only through DP Technology's Bohrium notebook
    (bohrium.com/notebooks/38543442597), which reads (checked 2026-11-12): "Licensing
    Agreement: This work is licensed under the Creative Commons
    Attribution-NonCommercial-ShareAlike 4.0 International License";
  - a container on Zenodo (record 19627026, P. Polishchuk, for EasyDock; the record itself CC
    BY 4.0) says the model file it took from that notebook "is distributed under Apache 2.0".
    The notebook says otherwise, and the Apache-2.0 licence belongs to the code;
  - the Apache-2.0 re-host (Hugging Face `bindresearch/pkasso_unipka`) does not say whether it
    trained them;
  - chemsuite is used inside a company, so non-commercial terms rule it out.

## Name and credit

It is not Uni-pKa, and the app does not call it that, since users would expect the paper's
accuracy. In the app it is **chemsuite pKa** (the Uni-pKa method on the Uni-Mol network). The
user guide and the docstrings cite both:

- **Uni-Mol** (Zhou et al., ICLR 2023; MIT): the network and the pretrained backbone;
- **Uni-pKa** (Luo et al., JACS Au 2024; Apache-2.0): the microstate enumerator and templates,
  the charge embedding, the free-energy training objective and the recipe. The enumerator and
  forward pass are ports of its code, so its notice ships (`data/LICENSE-Uni-pKa.txt`).

The weights are chemsuite's own, trained on pKaHub's compilation (citing pKaHub, and ChEMBL as
CC BY-SA 3.0 requires).

## Licence chain

| Part | Source | Licence | Use |
|---|---|---|---|
| Model and training code | dptech-corp/Uni-pKa | Apache-2.0 | Trained with, in WSL; the forward pass ported to numpy |
| Microstate templates | Uni-pKa's `smarts_pattern.tsv` (full; the simple one kept for tests), from MolGpKa's | Apache-2.0, MIT | Vendored (`chemsuite/ionisation/unipka/data`) |
| Pretrained backbone | Uni-Mol `mol_pre_all_h_220816.pt` (v0.1 release) | MIT | Initial weights |
| Training data (shipped model) | pKaHub (Sipos-Szabó et al., JCIM 2026, doi:10.1021/acs.jcim.6c00107; keserulab/pkahub): Baltruschat's ChEMBL set; the Hunt, Yu, Settimo, Manchester, Caine and Jensen sets; OCHEM; QSARToolbox; with pKaHub's Epik assignments | MIT (the compilation); ChEMBL CC BY-SA 3.0 | The shipped weights cite pKaHub and ChEMBL |
| **Not used** from pKaHub | IUPAC digitised (CC BY-NC), DataWarrior (unknown), the AttenGpKa set (may repeat both) | – | – |
| Test data | SAMPL6–8, Novartis acid/base (in the Uni-pKa repo) | as published | Evaluation only, never trained on |
| **Not used:** Uni-pKa's weights | Bohrium, Hugging Face re-host | unclear / NC | – |
| **Not used:** dwar-iBond for shipping | iBonD (nonprofit academic only), DataWarrior (unknown provenance) | – | A local reference model only, never committed or shipped |
| **Not used:** IUPAC digitised pKa | IUPAC | CC BY-NC 4.0 | – |

## Steps

- **U0, done 2026-11-12:**
  - environment in WSL2: `~/bench/unipka` (Python 3.11, torch 2.5.1 + CUDA 12.4, Uni-Core
    without its CUDA extensions, the Uni-pKa and Uni-Mol clones);
  - Uni-pKa's `dict.txt` equals Uni-Mol's, and the backbone loads (only the charge embedding
    and the head start fresh);
  - `scripts/unipka/train.sh` fine-tunes one fold on one GPU, at about 4 s per update of 32
    molecules.
  - Still to do: the reference fold on dwar-iBond, against the paper's numbers.
- **U1, done 2026-11-12:** `scripts/unipka/pkahub.py` builds the training set from pKaHub's
  `pka_all_datasets.tsv`, which already lists both sides' microstates (Epik-assigned).
  - Dropped: the sources not shipped, temperatures outside 20–30 °C, assignment errors over
    2 units, atoms charged beyond ±1, and repeats spread over more than 1 unit.
  - Held out, by InChIKey skeleton (669): pKaHub's test sets, the Uni-pKa test files and the
    neutron ligands.
  - Result: 22,132 transitions with a mean pKa of 6.36, 2.7 × Uni-pKa's dwar-iBond.
  - A first attempt, a ChEMBL pull with MolGpKa guessing each pKa's transition, was dropped
    for this.
- **U2:** train 5 folds of each model and evaluate on the test sets.
  - Speed: numpy takes about 0.1 s per conformer per fold for a drug-sized microstate, single
    threaded. Upstream's 11 conformers × 5 folds is about 5 s per microstate, too slow for a
    20-microstate ligand.
  - So U2 measures how much accuracy fewer conformers and folds cost, or distils a smaller
    model.
  - **Results, 2026-11-14** (shipped weights, `scripts/unipka/evaluate.py --enumerated`; the
    356 test molecules; RMSE):

    | microstates | templates | 3 conformers + 2D, 1 fold | 3 + 2D, 5 folds | 10 + 2D, 5 folds | failed |
    |---|---|---|---|---|---|
    | given (the transition's own) | – | 0.85 | 0.84 | 0.83 | 0 |
    | enumerated | simple | 0.85 | 0.83 | 0.82 | 24 |
    | enumerated | full | 0.84 | 0.82 | 0.80 | 2 |

    Uni-pKa's paper (its own weights, ChEMBL pretraining): 0.83 on the same sets. The simple
    templates had no site for 24 weak N–H acids (imides, hydantoins, acylated aminothiazoles,
    rhodanines); the full set finds 22 of them, and Novartis acids improve (1.15 → 1.08), so it
    is the bundled default. The full set costs about 21% more microstates.
  - **20 epochs, 2026-11-15** (fold 0 alone, same data and settings; full templates; RMSE):
    the shipped folds all peaked at their 10th epoch, so fold 0 was trained to 20. Validation
    RMSE 0.742 → 0.699 (best at epoch 19).

    | fold 0 | given, 2D only | given, 3 + 2D | enumerated, 3 + 2D | enumerated, 10 + 2D |
    |---|---|---|---|---|
    | 10 epochs (shipped) | 1.02 | 0.85 | 0.84 | 0.83 |
    | 20 epochs | 0.96 | 0.84 | 0.81 | 0.80 |

    One 20-epoch fold matches the five 10-epoch folds (0.82). Novartis acids 1.08 → 1.03,
    SAMPL6 0.74 → 0.64, SAMPL8 0.78 → 0.72; SAMPL7 (20 molecules) 0.82 → 0.94.
  - **All 5 folds at 20 epochs, 2026-11-17** (`samples/unipka/weights/e20.npz`, 440 MB;
    `samples/unipka/e20.eval.txt`, `e20.enum.eval.txt`). Validation RMSE of folds 1–4: 0.760,
    0.711, 0.698, 0.767. External tests, full templates (RMSE):

    | 20 epochs | 2D only | 1 + 2D | 3 + 2D | 10 + 2D |
    |---|---|---|---|---|
    | given, 1 fold | 0.96 | 0.87 | 0.84 | 0.83 |
    | given, 5 folds | 0.93 | 0.86 | 0.83 | 0.82 |
    | enumerated, 1 fold | – | 0.83 | 0.81 | 0.80 |
    | enumerated, 5 folds | – | 0.82 | 0.80 | 0.79 |

    Better than the 10-epoch folds throughout (enumerated, 3 + 2D, 5 folds 0.82 → 0.80) except
    SAMPL7, where 5 folds (0.94–0.99) do worse than 1 (0.90–0.95). MolGpKa on the same sets:
    1.02. Five folds buy only 0.01–0.02 for 5 × the time, and 3 conformers come within 0.01
    of 10.
  - Errors by acid group (`novartis_acid`, SAMPL7): ring N–H acids (RMSE 1.66) and
    sulfonamides (SAMPL7 0.88 too acidic on average) are worst; the training set has about
    590 sulfonamide and 570 ring N–H transitions against 11,900 bases. More data in those
    classes, not label cleaning, is the likely fix.
- **U3:** `chemsuite/ionisation/unipka/`:
  - `microstates.py`: upstream's enumerator, identical on SAMPL6 (with the simple templates),
    with atom maps kept;
  - `conformers.py`: upstream's recipe;
  - `model.py`: matches torch to 1.2e-4 with fp16 weights;
  - `thermo.py`;
  - `Predictor`, which caches by microstate.
  - The weights archive is converted by `scripts/unipka/convert_weights.py` and checked by
    `scripts/unipka/golden.py`.
- **U4:** ligand states in `macromolecule/prepare/ligand_sites.py` from the model's free energies. A
  microstate's energy at pH is `G + ln10 · n · (pH − μ)` kT, with μ the training mean.
  - **Built 2026-11-17** (plan approved that day):
    - The app runs fold 0 of `e20` (`e20_f0.npz`, 1 fold, 3 conformers + 2D), released as
      GPL-3.0-or-later. The release is `pka-weights-1`, with Uni-Mol's MIT license and a NOTICE.
      `ionisation/unipka/weights.py` pins the URL and SHA-256.
    - Prepare and Inspect ask once per session to download it (`tools/macromolecule/pka_model.py`).
      They use MolGpKa until it is downloaded, when it is declined, or when IT policy blocks the
      network.
    - `ligand_pka.LigandPka` is the model's interface; `ligand_sites._model_forms` makes the
      forms. Each form's `intrinsic` is `kT·G − RT ln10·n·μ`; `within`, `attract` and
      `tautomer_cost` are 0, because G already holds them.
    - RDKit's tautomers of the best microstate are scored by the model too.
    - Upstream's enumerator stops one charge short of `lower` and `upper`, so the default
      −2..+2 gives charges −1..+1, and that is what U2 evaluated.
    - Prepare starts there and widens, up to charges −4..+3, while the best form holds the
      lowest or highest charge. The wider range has not been evaluated.
    - A ligand with more than 48 microstates uses MolGpKa and says so in the report. ATP has
      455 microstates, takes 19 s just to enumerate, and the model takes about 0.45 s per
      microstate.
    - The report's **Ligands** line names the model that ran.
  - **Ligand hydrogens on the 262 neutron entries**, at the default weights
    (`prepare_review.py fit_ligands --pka`, against the same run without `--pka`):

    | | MolGpKa | chemsuite pKa |
    |---|---|---|
    | Training half | 97.0% (815/840) | 96.7% (812/840) |
    | Test half | 94.7% (643/679) | 94.7% (643/679) |
    | Test, unlike the CCD | 3 of 33 | 3 of 33 |

    - It is a tie. The prior toward the CCD drawing decides most atoms, and too few atoms differ
      from the drawing to tell the two models apart.
    - Refitting the ligand weights for the model (`ccd` 8 → 6) overfits: training goes to
      96.9%, test drops to 94.3%. So the defaults stay.
    - The case for the model rests on the external pKa tests (0.81 against MolGpKa's 1.02) and on
      its handling of sites that affect each other.
- **U5:** Predict pKa, standardise and the descriptors move over; MolGpKa is removed.

## Notes

- **BLAS threads:** numpy's multithreaded BLAS loses badly on these small matrices when the CPU
  is busy (oversubscription). Workers should run with one BLAS thread.
- **Contiguous embeddings:** a broadcast (zero-stride) embedding array makes numpy skip BLAS;
  the model copies it out.
