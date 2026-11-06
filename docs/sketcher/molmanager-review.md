# Sketcher P0: review of MolManager's sketcher code (2026-11-05)

**Rule (decided by the user).** The new sketcher is an original implementation.
- MolManager's `iupac_*` modules are not used, and neither is anything that imports them.
- Any other MolManager code must be verified against an independent, cited reference before
  any of it, or its data, is reused.

**Reproduce:**

```
python scripts/review/molmanager_review.py \
  <MolManager>/molmanager/ui/sketcher <ketcher>/packages/ketcher-react/src/templates/fg.sdf
```

The MolManager modules are loaded read-only for testing; nothing is copied into chemsuite.

## 1. `alkene_stereo.py`: E/Z labels

**Claim:** E/Z by Cahn–Ingold–Prelog priority.

**How it actually ranks:** with RDKit's `CanonicalRankAtoms`, a canonical numbering that has no
relation to CIP priority.

**Reference:** RDKit's new CIP labeller (`rdCIPLabeler`; Hanson et al., *J. Chem. Inf. Model.*
2018, 58, 1755) on 30 alkenes of known configuration, including tamoxifen, geraniol, nerol,
linoleic acid and sorbic acid. Their 2D coordinates come from the stereo SMILES.

**Result: 28/32 double bonds right, 4 wrong.**

| Molecule | CIP | MolManager |
|---|---|---|
| `C/C(=C\c1ccccc1)C(=O)O` (α-methylcinnamic acid) | E | Z ✗ |
| `CO/N=C/c1ccccc1` (oxime ether) | E | Z ✗ |
| `C/C(I)=C(/Br)Cl` | Z | E ✗ |
| `OC(=O)C(/C)=C(\C)C(C)(C)C` | E | Z ✗ |

Two of these were checked by hand:
- In `C/C(I)=C(/Br)Cl`, I outranks CH₃ and Br outranks Cl, and they are cis, so it is Z.
- In α-methylcinnamic acid, COOH outranks CH₃ and is trans to phenyl, so it is E.

**Verdict: drop.** The sketcher takes E/Z, and R/S, from `rdCIPLabeler`.

## 2. `contracted_labels.py`: abbreviations (105 labels)

**Claim:** contracted labels following IUPAC GR-2.2/2.3, with the first atom of each SMILES as
the attachment point.

**References:**
- Ketcher's `fg.sdf` (attachment atom from its `SAP` record), compared by canonical SMILES of the
  group bonded to `[*]`;
- the isomer each label's name or alias states.

**Results:**
- **34 agree** with Ketcher, 70 are not in Ketcher, and **1 differs**: *Mes*. Here MolManager
  is right (2,4,6-trimethylphenyl) and **Ketcher's record is wrong** (3-methylbenzyl, 8 carbons);
  see `licence-review.md`.
- **Naph is wrong.** Its alias says "1-NAPH", but the SMILES `c1ccc2ccccc2c1` attaches at C2,
  so it is 2-naphthyl.
- **Pr and nPr are defined twice.** Both register the alias "NPR", so whichever is registered
  last wins. The data is not self-consistent.
- The other named isomers (p-tolyl, mesityl, 1-adamantyl, sec-butyl, isobutyl, PMB) are right.

**Verdict: do not reuse as is.** Phase 2 builds the abbreviation table from Ketcher's `fg.sdf`.
- Every entry is checked against an independent definition: formula, and the attachment atom
  given by the name.
- MolManager's 70 extra labels are candidates to add, each one checked the same way. The Naph
  error shows they cannot be taken on trust.

## 3. `valence.py`: implicit hydrogens

**Claim:** valences from RDKit's allowed-valence lists, adjusted for charge.

**Reference:** the 61 species in `valence-review.md`.

**Result: 49/61 right.** It fails on chloride, fluoride, methanethiolate, triiodide,
trimethylsulfonium and diethylzinc (it gives each 2 hydrogens; all have 0). It also gives
*negative* hydrogen counts for hypervalent drawings: nitro −2, chloric acid −4, perchloric
acid −6.

**Verdict: drop.** The sketcher uses RDKit directly (`valence-review.md`).

## 4. The other modules

| Module | Finding | Verdict |
|---|---|---|
| `iupac_labels`, `iupac_orient`, `iupac_rings`, `iupac_style`, `iupac_validate` | Excluded by decision | **Not used** |
| `fused_templates.py` | Imports `iupac_rings._largest_gap_bisector` | **Not used.** Fusing is re-implemented from Ketcher's rules (plan §4.3). |
| `acs_style.py` | Restates RDKit's `SetACS1996Mode` parameters | **Drop.** Call RDKit. |
| `element_colors.py` | Copy of RDKit's default palette | **Drop.** Ask RDKit. |
| `rgroups.py` | R*n* stored as a dummy atom with atom-map *n* | **Drop.** Use RDKit's molfile R-group fields (`RGP`, `_MolFileRLabel`), so R-groups survive a molfile round trip. Map numbers stay an export option for chemsuite's R-group tools. |
| `sketch_grid.py`, `sketch_reactions.py` | Small geometry helpers (30° grid, arrow snapping) | **Rewrite.** New geometry follows Ketcher's rules and is tested on its own. |
| `bonds.py`, `sketch_graph.py`, `constants.py`, `chem.py`, `wildcards.py` | Data structures and helpers tied to the old widget | **Rewrite** in the new model. |
| `widget*.py`, `dialog.py`, `toolbar_glyphs.py`, `customize_elements.py`, `physical_properties.py`, `bioisostere_dialog.py` | Qt user interface; `dialog.py` reaches into the main window 10 times | **Rewrite** |

## Summary

No MolManager sketcher code passes verification as is.
- Its useful content is a list of extra abbreviation names to check and possibly add.
- The rest is replaced by RDKit (valence, CIP, colours, ACS style) or by rules re-implemented
  from Ketcher, each with tests.
