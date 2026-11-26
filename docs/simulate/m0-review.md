# M0 review: the feasibility gates for minimise-after-Prepare (2026-11-15)

The gates of [`SIMULATE_PLAN.md`](../SIMULATE_PLAN.md) M0, run on the Windows PC (RTX 4070,
driver 560.94) natively and in WSL Ubuntu. Scripts: [`m0_spike.py`](m0_spike.py) (gates 1–5,
run in the OpenMM env) and [`m0_wsl.py`](m0_wsl.py) (gate 6, run with chemsuite's interpreter).
**All six gates pass; the plan's design stands**, with the changes listed at the end.
Not run: osx-arm64 (gate 1 on the Mac is still owed).

## The envs that were tested

| | native Windows (`D:\chemsuite-envs\cs-sim-m0`) | WSL (`~/miniconda3/envs/cs-sim-m0`) |
|---|---|---|
| python | 3.12.16 | 3.12.16 |
| openmm | 8.4.0 | 8.4.0 |
| openmmforcefields | 0.15.1 (pip, no deps, from the GitHub tag) | 0.15.1 (conda) |
| openff-toolkit-base / interchange-base | 0.19.0 / 0.5.5 | 0.19.0 / 0.5.5 |
| openff-nagl / nagl-models | 0.6.1 / 2026.09.0 | 0.6.1 / 2026.09.0 |
| pytorch | 2.14.1 CPU | 2.14.1 CPU |
| rdkit | 2026.09.2 | 2026.03.1 |
| AmberTools | none | 24.8 (a dependency of conda's openmmforcefields) |
| OpenMM platforms | Reference, CPU, CUDA, OpenCL | Reference, CPU, CUDA |

Both from conda-forge with `nodefaults`:
`python=3.12 openmm cuda-version=12.6 pytorch-cpu openff-toolkit-base openff-interchange-base
openff-nagl openff-nagl-models rdkit numpy lxml`, plus `openmmforcefields` (conda on Linux; on
Windows `pip install --no-deps git+https://github.com/openmm/openmmforcefields@0.15.1`).

## Gate 1 — the native env: passes, with four things the env file must do

1. **conda-forge's `openmmforcefields` cannot be installed on Windows**: every build requires
   `ambertools`, which has no win-64 package, and `conda install --no-deps` still solves and
   fails. It is not on PyPI either. It is pure Python, so pip installs it from the GitHub tag
   without dependencies; the one dependency that then has to be named is `lxml`. With that it
   works with no AmberTools at all.
2. **The CUDA version must match the driver.** Left alone, conda installs CUDA 12.9 libraries;
   driver 560.94 supports 12.6, and OpenMM then reports "No compatible CUDA device". Pinning
   `cuda-version=12.6` fixes it, but **OpenMM 8.6 is only built for CUDA ≥ 12.9**, so the pin
   holds OpenMM at 8.4.0, and 8.4.0 needs openmmforcefields 0.15.1 (0.16.0 imports
   `MergedResidue`, which 8.4 lacks). A newer NVIDIA driver would allow 8.6.1 + 0.16.0; gates
   2–5 were also run on that pair (CPU only) and gave the same results.
3. **PyTorch must be the CPU build** (`pytorch-cpu`). NAGL charges take about 1 s per ligand on
   the CPU; the CUDA build makes the env 12 GB instead of 3.5 GB and gains nothing.
4. **The runner must put the env's folders on `PATH` on Windows** (`<env>`, `<env>\Library\bin`,
   `<env>\Scripts`). Started by its `python.exe` alone, as the app does, the CUDA plugin fails
   to load (error 126); with the folders on `PATH` it loads. This belongs in
   `child_environment`.

Also: a conda `solver: classic` setting makes this solve take many minutes; the env docs should
say `--solver libmamba`. Disk: about 3.5 GB per env plus the package cache.

## Gate 2 — ligand charges without AmberTools: passes

- Sage 2.3 (`openff-2.3.0`) carries a `NAGLCharges` section naming
  `openff-gnn-am1bcc-1.0.0.pt`; with `openff-nagl` installed, `SMIRNOFFTemplateGenerator`
  assigns those charges itself.
- A `Molecule` handed in with charges is used **exactly as given** (max difference 0.0 e on
  ibuprofen anion, benzamidinium and imatinib), and the charges NAGL computes directly are
  identical to what the generator computes on its own (0.0 e).

So the plan holds: the worker computes NAGL charges once, they are stored in the document, and
every later system is built from the stored ones.

## Gate 3 — implicit solvent with a Sage ligand: passes, GBn2 stays

`SystemGenerator(forcefields=["amber/ff14SB.xml", "implicit/gbn2.xml"],
small_molecule_forcefield="openff-2.3.0")` builds a peptide + ligand system whose `CustomGBForce`
covers every particle, ligand included, and it minimises. The ligand's GB parameters are set by
element (C, O and H each get one set; OBC2 likewise). The explicit-shell fallback is not needed.

Not checked: GBn2's bonded-atom special cases (a hydrogen on nitrogen, for one) on ligand
atoms; the test ligand has only C, O and C-bound H. Worth a look in M3 before MM-GBSA numbers
are trusted, not a blocker for minimisation.

## Gate 4 — template coverage: the same in ff14SB and ff19SB

| Prepare state | Template |
|---|---|
| HID, HIE, HIP | yes (and N-/C-terminal forms) |
| ASH, GLH, LYN | yes, mid-chain only (no N-/C-terminal form) |
| CYM | yes, mid-chain only |
| CYX | yes (and terminal forms) |
| ACE, NME (and NHE) | yes |
| **TYM** (tyrosinate) | **none** |
| **neutral N-terminus** (the `NH2` state) | **none** |
| **`HIS-`** (imidazolate, a histidine bridging two metals, `hydrogens.py:230`) | **none** |
| **modified residues** in Prepare's templates (ALY, CGU, CME, CSD, CSO, FME, …) | **none** |
| SEP, TPO, PTR (dianion) and S1P, T1P, Y1P (monoanion) | yes, mid-chain only, from openmmforcefields' `amber/phosaa14SB.xml` and `phosaa19SB.xml` (added 2026-11-26; phosaa14SB needs `amber/ff14SB.xml` as its base, whose energies match `amber14/protein.ff14SB.xml` term for term on crambin) |

Per the plan, a state with no template blocks with a ⚠ problem and is never substituted. Two
of these were not in the plan's list: `HIS-`, and the modified residues. ASH/GLH/LYN/CYM on a
chain's first or last residue block too.

## Gate 5 — massless atoms: the minimiser leaves them fixed

With heavy-atom masses set to zero, `LocalEnergyMinimizer` moves heavy atoms 0.0 nm
(Reference, CPU) or 5e-8 nm (CUDA, OpenCL: single-precision round-off) while hydrogens move
up to 0.2 nm. **The system must be built without constraints**: with `HBonds` every platform
raises "A constraint cannot involve a massless particle". The hydrogen-only stage therefore
uses zero masses on an unconstrained system; no stiff restraints are needed.

## Gate 6 — WSL mechanics: passes

- A worker started as `wsl.exe -e python3 -I -c …` through `children.start` answers in 0.3 s
  and round-trips UTF-8 JSON lines.
- `children.stop` (the Job Object kill of `wsl.exe`) ends the Linux process at once, including
  a busy one that never reads stdin.
- Closing stdin ends the worker (exit code 0), so it cannot be orphaned.
- Paths: `wslpath -a` gives `/mnt/c/…`; the worker script itself travels in `-c`, so only data
  paths need translating.

## Changes to the plan

1. **Env files (M2):** pin `cuda-version` to the driver, `pytorch-cpu`, and on Windows install
   openmmforcefields by pip without dependencies. Decide between OpenMM 8.4.0 (works with the
   current driver) and 8.6.1 (needs a newer NVIDIA driver).
2. **`child_environment` (M2):** add the env's `Library\bin` folders to `PATH` on Windows.
3. **Blocking states (M1 `build.py`):** TYM, the neutral N-terminus, `HIS-`, every modified
   residue, and ASH/GLH/LYN/CYM at a chain end.
4. **Stage 1 (M3):** zero masses on heavy atoms, system built with `constraints=None`.
5. **Gate 3's fallback is dropped:** minimisation runs in GBn2.
