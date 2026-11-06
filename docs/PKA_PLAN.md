# Plan: pKa prediction with Uni-pKa (drafted 2026-11-05)

## 1. Goal

Add a **Predict pKa** tool that writes the following for each row:

- macro pKa values;
- the charge and the dominant protomer at a chosen pH;
- a model-based cLogD.

Uni-pKa supplies the free energy of each protonation microstate; chemsuite enumerates the
microstates and applies the thermodynamics. This settles ROADMAP §6's open question about the
"(est.)" descriptors. Those currently guess a basic pKa from the nitrogen count
(`chemsuite/descriptors/medchem.py`).

MolManager already did this, in about 2,100 lines across `molmanager/ionization/`. Its chemistry
is worth porting. Its wiring (mixins, GPU hints, its own result cache) is not.

## 2. What exists upstream

| Piece | Source | Licence | Notes |
|---|---|---|---|
| Uni-pKa model code | dptech-corp/Uni-pKa (Luo et al., JACS Au 2024) | Apache-2.0 | No releases. Its environment is old: Python 3.8, a CUDA Docker image. |
| `unipkainfer` 0.2.0 (PyPI, 2026-07) | Inference-only repackaging of Uni-pKa, Uni-Mol and a pure-Python Uni-Core | Apache-2.0 AND MIT | CPU or single-GPU. Requires `torch>=2.2,<3`, `lmdb`, `pandas`, `pooch`, `platformdirs`. Supports Python 3.10–3.13. This is what MolManager used. |
| Weights (5 folds, `checkpoint_best.pt`) | Hugging Face `bindresearch/pkasso_unipka`, fetched by `unipkainfer` on first use and checked by SHA-256 | Repo says **Apache-2.0** | A third-party re-host (Bind Research, "for use in pKasso"). |
| Original weights | DP Technology's Bohrium notebook | Reported as **CC BY-NC-SA 4.0** | Non-commercial. This conflicts with the re-host's Apache-2.0 tag. |
| Microstate templates | MolGpKa SMARTS (`smarts_pattern.tsv`, `simple_smarts_pattern.tsv`) | MIT | MolManager vendors them, with notices. |

## 3. Decisions needed before building

1. **Weights licence (blocker).** chemsuite is for an internal team at a company, which counts as
   commercial use. If the Bind Research weights really derive from DP Technology's checkpoint, the
   NC licence may still apply whatever the re-host says. This needs a yes from legal, or written
   terms from DP Technology. If the answer is no, the fallback is to train our own weights on the
   public datasets with the Apache-2.0 code. That is a project of its own and is not covered here.
2. **Install size.** The CPU build of torch is about 250–700 MB unpacked, depending on the OS, and
   the weights are about 0.5 GB per fold. The base install budget is 550 MB unpacked (PLAN.md §1),
   so pKa **must be an optional add-on**:
   - `pip install chemsuite[pka]`;
   - a separate installer variant built by IT.

   The tool's `requires=("torch", "unipkainfer")` already greys it out, with the reason, when the
   add-on is missing.
3. **Where the weights come from.** Corporate machines may not reach Hugging Face, and the settings
   have an offline switch. Recommendation:
   - Look in an **IT-managed model folder** first. It is set in settings or the IT policy file, and
     `unipkainfer`'s `model_dir=` points at it.
   - Download only when that folder is empty and the app is online, and say so in the dialog
     first.
   - Ship **one fold** (fold 1, about 0.5 GB). A five-fold ensemble becomes a setting, off by
     default.
4. **GPU.** Start with CPU only, which is simpler for IT and works everywhere. MolManager's CUDA
   detection and hints (`pka_gpu_hint.py`) can come later if CPU speed is not enough.

## 4. How it fits the architecture

- **New Qt-free package `chemsuite/ionisation/`.** Add it to `QT_FREE` in
  `tests/test_import_boundaries.py`. It contains:
  - `enumerate.py`: a port of MolManager's `unipka_enumerator.py`, with the MolGpKa templates
    vendored under `chemsuite/ionisation/data/` along with their licence notices. Plain RDKit.
  - `thermo.py`: pure numpy, ported from `unipka_ensembles.py`. Free energies become macro pKa
    values, populations at a pH, mean charge and pI. It also computes logD as
    `logP + log10(fraction neutral)`.
  - `model.py`: a small `FreeEnergyModel` protocol, `free_energies(mols) -> np.ndarray`.
    - `UnipkaModel` imports `unipkainfer` and `torch` lazily, only when a run starts, so startup
      and the base install never load torch.
    - A `FakeModel` (for example, a sum of per-site contributions) lets the whole pipeline be tested
      without torch.
  - `predict.py`: `predict_pka(mols, model, ph, progress, cancel)`. It works row by row and raises
    `Cancelled(partial)` with the rows it has finished, following the partial-results rule.
- **Model host.** Loading the weights costs seconds and about 1 GB of memory, so the shared process
  pool must not reload them in every worker.
  - One long-lived spawn-context process loads the model once and serves batches of microstates.
    It pins torch's thread count to the core budget.
  - It shuts down with the app, or after it has been idle for a few minutes.
  - MolManager's `unipka_score_host.py` is the template, including its weights-only, memory-mapped
    load.
- **Cache.** Keep free energies in memory, keyed by canonical microstate SMILES. Series share most
  microstates, so a re-run or a change of pH is close to free. A disk cache can wait until someone
  asks for it.
- **Tool `chemsuite/tools/pka/`**, built on the dialog kit:
  - **Input:** structures and Rows.
  - **Method:** pH (default 7.4); "Ensemble of 5 models (slower)".
  - **Output:** choose which columns to write; If it exists.
  - The tool uses `ctx.run(..., keep_partial=True)` and `ctx.apply_results`, and nothing else from
    the app.
- **Output columns**, all optional, with names stated in the dialog footer:
  - `pKa (acidic)`, `pKa (basic)`: the strongest of each kind. A full list such as `pKa values`
    is available as text.
  - `Charge @ pH 7.4`, `% Neutral @ pH 7.4`.
  - `Major Species @ pH 7.4`: a structure column holding the dominant protomer, plus
    `Major Species %`.
  - `cLogD 7.4 (Uni-pKa)`: a separate column, so the fast, model-free "(est.)" descriptor stays
    available for people without the add-on.
- **User guide:** a new "Predict pKa" section explaining the method, what the columns mean, its
  limits (zwitterions, tautomers, metal complexes), and how to install the add-on.

## 5. Phases

| Phase | Work | Done when |
|---|---|---|
| P0 spike (1–2 days) | Install `unipkainfer` and CPU torch in a fresh venv on macOS and Windows. Measure the unpacked size, the model load time and the time per molecule. Run 50 reference compounds through MolManager and record its outputs as golden values. | The numbers are written into this plan; the licence answer is in hand. |
| P1 core | `chemsuite.ionisation`: enumerator port, thermodynamics, `FakeModel`, `predict_pka` with partial results. | Unit tests pass without torch; the golden values match MolManager to 0.01 when the real model is installed (`-m pka` marker). |
| P2 model host | The single model process, cancellation between batches, the free-energy cache, the idle shutdown. | Cancelling stops within one batch; a second run on the same set uses the cache. |
| P3 tool | Dialog, output columns, a status line ("Predicting pKa… (x/y)"), the user guide section. | It works end to end on 1,000 rows; the guide is updated. |
| P4 protomers | `Major Species` structure column, populations, and a "Protomers at pH" tool that expands one row into its states (MolManager's `protomer.py`). | — |
| P5 packaging | The `pka` extra; the IT model-folder setting; the installer variant; CI tests with `FakeModel` on all three OSes plus one Linux job with the real model, cached. | The installer variant stays inside an agreed size budget. |

## 6. Risks

- **Licence:** see §3.1. It can block the whole feature.
- **Speed on CPU:** this is unknown until P0. If it is too slow for 200k rows, offer "Shown" or
  "Selected" rows by default and rely on the cache. GPU support moves forward if needed.
- **torch on corporate Windows:** large wheels and antivirus scanning of `.pyd` files. IT should
  test the installer variant early (P0).
- **`unipkainfer` is a young, single-maintainer package.** Pin it to an exact version. The
  `FreeEnergyModel` protocol keeps it swappable for a vendored copy if needed; it is Apache and MIT
  licensed, so a vendored copy is allowed.
