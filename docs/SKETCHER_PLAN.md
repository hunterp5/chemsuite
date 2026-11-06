# Plan: a native structure sketcher (drafted 2026-11-05)

## 1. Goal and ground rules

The goal is a structure sketcher built into chemsuite that chemists can use as they would a
commercial one (ChemDraw, Marvin). It should let them:

- draw and edit a structure cell in place;
- draw queries for structure search and filters;
- add new compounds.

It replaces the text-first editor (`chemsuite/app/structure_editor.py`), and keeps that editor's
paste and type-SMILES paths as one way in.

The rules for this work:

1. **An original implementation.** This is new code written for chemsuite. It is not a
   translation of another program.
2. **Ketcher (EPAM, Apache-2.0) is the reference for chemistry and interaction rules.**
   - We re-implement its rules in Python and cite the source file each rule comes from.
   - Its data files (abbreviations, salts and solvents, templates) may ship as they are,
     with Ketcher's licence and NOTICE (see §3).
3. **Nothing from MolManager's `iupac_*` modules is used.** That covers `iupac_labels`,
   `iupac_orient`, `iupac_rings`, `iupac_style` and `iupac_validate`, and anything that imports
   them.
4. **Every rule is scientifically verified before it ships, including any taken from
   MolManager.** Each rule has a cited source and tests against independent reference data
   (§5). "MolManager did it this way" is not a source.

## 2. Architecture

```
chemsuite/sketch/            Qt-free (added to QT_FREE in tests/test_import_boundaries.py)
  model.py      the drawing: atoms and bonds with 2D coordinates, charges, isotopes,
                radicals, query features, S-groups (abbreviations), R-groups,
                attachment points; later arrows and text
  edits.py      every change as an undoable edit (add bond, merge atoms, attach
                template, …) built from small reversible operations
  valence.py    implicit hydrogens and valence errors (§4.1)
  geometry.py   where new atoms go: bond placement, chains, snapping (§4.2)
  templates.py  attaching rings and fused templates to atoms, bonds or empty canvas (§4.3)
  stereo.py     wedge and hash validity, CIP labels, enhanced stereo (§4.4)
  abbrev.py     functional-group abbreviations: contract and expand (§4.5)
  convert.py    model to and from RDKit and molfile (V2000/V3000), SMILES and SMARTS
  data/         Ketcher's LICENSE and NOTICE (for ctfile.py); later salts-and-solvents.sdf
                and library.sdf after their checks
chemsuite/app/sketcher/      Qt
  canvas.py     the drawing surface, built on QGraphicsView, drawn with chemsuite's depiction
                code so a sketch looks like the table's drawings
  tools.py      the tools (select, bond, chain, ring and template, atom, charge, eraser,
                R-group, …), each turning mouse and keyboard input into edits
  toolbar.py, dialog.py
```

- **Why an own model and not RDKit's editable molecule.** A sketch must hold states RDKit
  rejects or rewrites while you draw: a pentavalent carbon in mid-edit, an abbreviation as one
  node, query atoms, attachment points, arrows. RDKit is the **chemistry oracle**. It sanitises,
  perceives aromaticity, assigns CIP labels and lays out 2D coordinates through `convert.py`.
- **The Qt side holds no chemistry.** Tools only issue edits, so every chemistry rule can be
  tested without a GUI.

## 3. Sources and licences

| Source | Licence | Used for |
|---|---|---|
| Ketcher, github.com/epam/ketcher | Apache-2.0 (with NOTICE) | Rules re-implemented from `ketcher-core` (`domain/entities/atom.ts`, `application/editor/actions/{utils,bond,atomMerge,closelyFusing,aromaticFusing,template,chain}.ts`, `domain/helpers/stereoValidator.ts`, `domain/constants/generics.ts`) and `ketcher-react` (`script/editor/tool/*.ts`). Data: `templates/fg.sdf` (62 abbreviations), `salts-and-solvents.sdf` (135), `library.sdf` (276 templates). Its `__tests__` provide expected behaviour for our tests. |
| BIOVIA CTfile Formats (2020), the valence model | Specification | Implicit hydrogens and valences. Ketcher cites it (`atom.ts`), and molfiles from ChemDraw and others assume it. |
| IUPAC *Graphical representation standards for chemical structure diagrams* (Brecher, Pure Appl. Chem. 2008, 80, 277–410) | Published recommendation | Drawing conventions: bond angles, label placement, contracted labels. Used as a reference text only; no code from MolManager's `iupac_*` modules. |
| IUPAC *Graphical representation of stereochemical configuration* (Brecher, Pure Appl. Chem. 2006, 78, 1897–1970) | Published recommendation | Which wedge and hash drawings are valid or ambiguous. |
| RDKit (BSD-3) | Dependency | Sanitisation, 2D layout, CIP labels (`rdCIPLabeler`, the new labeller of Hanson et al. 2018), molfile, SMILES and CDXML input. |

Ketcher's NOTICE lists further licences for its dependencies. We take only `ketcher-core` and
`ketcher-react` sources and their data files, so only Apache-2.0 applies.
- P0 confirmed this at commit `2b031cc2` (`docs/sketcher/licence-review.md`).
- **Ketcher's data is checked like any other source.** Its "Mes" record is 3-methylbenzyl,
  not mesityl, so every record is checked against an independent definition before it ships.

## 4. Chemistry rules (each needs a verification test before it ships)

### 4.1 Valence and implicit hydrogens

- **Rule (decided in P0, `docs/sketcher/valence-review.md`):** RDKit's valence model, run
  through the same molfile and sanitisation path as the rest of chemsuite.
  - On 61 real reference species RDKit got 59 right and Ketcher's CTfile model 57. Ketcher
    rejects PF₆⁻ and the other hexafluoro counter-ions.
  - A drawing is valid in the sketcher exactly when the table can use it.
- **Shown while drawing:** an atom whose valence cannot be satisfied is marked, with the reason,
  not refused (as Ketcher's `badConn` does).
- **Exporting a molfile** writes an explicit valence on any atom where a CTfile reader would
  count hydrogens differently from RDKit. RDKit's own writer misses some (it leaves PF₆⁻'s
  phosphorus unmarked), so `sketch/ctfile.py` carries a port of Ketcher's `calcValence` to
  find them; a test reads every exported reference species back through it.
- **Tests:** the reference species become a test (`scripts/review/valence_review.py` holds
  them). ClF₃ and IF₇ are recorded as known gaps.

### 4.2 Geometry: where new atoms go

- **New bond from an atom:** Ketcher's `atomForNewBond`:
  - bisect the largest gap between existing bonds;
  - draw chains as 120° zig-zags;
  - make a triple bond or cumulated double bond 180° to the previous bond;
  - reuse an existing atom when the new position lands within 0.1 bond lengths of it.
- **Dragging a bond:** it snaps to 15° steps (Ketcher's `fracAngle`, `FRAC = π/12`), and
  holding Ctrl frees it. The chain tool follows Ketcher's
  `chain.ts`.
- **Verified by:** geometry unit tests (angles to 0.1°, no overlapping atoms) on Ketcher's test
  cases plus our own, and round trips through RDKit that keep the drawn coordinates.

### 4.3 Templates and fusing

- **Rule:** dropping a ring or template on an atom (spiro or substituent), on a bond (fused,
  sharing the edge) or on empty canvas, with atom merging, follows Ketcher's `template.ts`,
  `atomMerge.ts`, `closelyFusing.ts` and `aromaticFusing.ts`. Templates come from `library.sdf`
  and a small ring bar (3- to 8-membered rings and benzene).
- **Verified by:**
  - fusing results compared by canonical SMILES (naphthalene from benzene + benzene,
    indole, steroid core and so on);
  - no overlapping atoms;
  - aromaticity matching RDKit's perception.

### 4.4 Stereochemistry

- **Wedges and hashes:**
  - placement rules and validity checks follow Ketcher's `stereoValidator.ts` and the IUPAC
    2006 recommendations;
  - invalid or ambiguous wedges are flagged.
- **CIP labels** (R/S, E/Z) come from RDKit's new CIP labeller (`rdCIPLabeler`), shown on
  request.
  - MolManager's `alkene_stereo.py` ranks substituents with `CanonicalRankAtoms`. Canonical
    ranks are not CIP priorities, so that module is not reused.
- **Enhanced stereo** (ABS, AND n, OR n groups) follows the V3000 molfile representation.
- **Verified against:**
  - a set of stereocentres and double bonds with literature descriptors (amino acids, drug
    stereocentres, the hard cases from Hanson et al. 2018), checked through a draw, molfile
    and RDKit round trip;
  - Ketcher's stereo tests.

### 4.5 Abbreviations, generic atoms and R-groups

- **Abbreviations** (Ph, Boc, CO₂Me, OTs …) come from chemsuite's own table, **not**
  `fg.sdf`.
  - P2 found 9 unusable records and a wrong Mes in `fg.sdf`, and a wrong iBu in RDKit's list
    (`docs/sketcher/abbreviation-review.md`).
  - They are stored as superatom S-groups, contracted and expanded without changing the
    molecule, and written as `SUP` groups in molfiles.
  - **Verified by:** each entry's formula from its IUPAC name, isomer checks, placement on a
    ring with canonical SMILES, molfile round trips; Ketcher and RDKit as cross-checks
    (`scripts/review/abbreviation_review.py`).
- **Generic and query atoms** (A, Q, X, M and their H variants; atom lists; NOT lists) map to
  SMARTS as in Ketcher's `generics.ts` and the CTfile specification.
  - **Verified by:** substructure searches against reference molecules.
- **R-groups:** `R1`…`Rn` with attachment points, written as V3000 RGP or as `[*:n]` for
  chemsuite's R-group tools.

### 4.6 Clean-up and layout

- **Clean structure:** RDKit's 2D coordinate generation, keeping what the user placed where
  possible (CoordGen, or constrained `Compute2DCoords`). Phase 3 compares both on a benchmark
  set before choosing.
- **Verified by:** on about 500 drug structures, counting overlapping atoms, bond-length spread
  and crossing bonds.

## 5. Verifying anything taken from MolManager

**Done in P0:** `docs/sketcher/molmanager-review.md`, reproducible with
`scripts/review/molmanager_review.py`. **No MolManager sketcher code passes as is.**

| Module | Result | Verdict |
|---|---|---|
| `iupac_*` (5 files), `fused_templates.py` (imports `iupac_rings`) | Excluded by decision | **Not used** |
| `alkene_stereo.py` | 28/32 double bonds agree with CIP; it ranks by canonical order, not CIP | **Drop**; use `rdCIPLabeler` |
| `valence.py` | 49/61 reference species; gives chloride 2 H, and negative H counts | **Drop**; use RDKit |
| `contracted_labels.py` | 34 agree with Ketcher; *Naph* is 2-naphthyl despite its "1-NAPH" alias; duplicate *NPR* alias | **Drop**; its 70 extra names are candidates to add after checking |
| `acs_style.py`, `element_colors.py` | Copies of RDKit settings | **Drop**; call RDKit |
| `rgroups.py` | Atom-map R-groups do not survive a molfile round trip | **Drop**; use molfile R-group fields |
| All other modules | Tied to the old widget, or small geometry helpers | **Rewrite** |

## 6. How it is tested

- **Rule tests** (Qt-free, fast): one file per rule area in §4, each citing its reference data.
- **Round-trip tests:**
  - draw a structure with edit sequences, export a molfile, read it with RDKit and compare the
    canonical SMILES;
  - read a molfile, edit it, export it and check nothing else changed.
- **Ketcher parity tests:** the expectations in Ketcher's `__tests__` for the rules we
  re-implement, rewritten as pytest cases with a note of their origin.
- **Interop corpus:** molfiles and CDXML exported from ChemDraw and Marvin by the team, and
  SMILES from the clipboard. Each must read in, look right and export back unchanged.
- **UI tests** (pytest-qt): tools driven by synthetic mouse and key events, checking the
  resulting model.

## 7. Phases

| Phase | Work | Done when |
|---|---|---|
| P0 ✅ | Licence check, MolManager review, valence review, interaction spec: `docs/sketcher/`. | Done 2026-11-05. Still open: ChemDraw hotkeys checked against the team's ChemDraw (interaction-spec, open item 1). |
| P1 ✅ | Model, edits with undo, valence, geometry, convert; canvas with select, bond, chain, atom, charge and eraser tools; molfile and SMILES copy and paste. | Done 2026-11-05: `chemsuite/sketch/` (plus `render.py` for drawing conventions and `ctfile.py` for export valences) and `chemsuite/app/sketcher/`, opened by **Draw…** in the structure editor. Moved to P2: lasso selection, the rotate handle, and double-click label editing (abbreviations need P2's `fg.sdf`). Ketcher's `0` (any bond) and generic-atom keys wait for P4. |
| P2 | Rings and templates (`library.sdf`), fusing, abbreviations (`fg.sdf`), salts and solvents, clean-up. | Template, abbreviation and layout tests pass. |
| P3 | Stereo: wedges, CIP display, enhanced stereo; layout comparison. | Stereo reference set passes. |
| P4 | Queries: generic atoms, atom lists, bond queries, R-groups, attachment points; hook into structure search and filters. | Query tests pass; search uses drawn queries. |
| P5 | Edit cells in place, add compounds, replace the text editor; user guide section. | Used for all structure entry in the app. |
| Later | Reactions and arrows, text, CDXML export. | — |

## 8. Open questions

- **Interaction model.** Should it follow ChemDraw's hotkeys exactly where the team knows
  them, or Ketcher's? The proposal is ChemDraw's for the keys chemists use daily, and Ketcher's
  for everything else.
- **CDXML.** RDKit can read it, which is enough for paste. Writing CDXML waits until someone
  needs it.
