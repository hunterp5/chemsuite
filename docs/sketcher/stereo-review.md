# Sketcher P3: stereo review (2026-11-06)

This review covers four things:
- the sketcher's R/S and E/Z labels against an independent reference;
- the rules behind the wedge warnings;
- "either" double bonds;
- enhanced stereo.

## 1. CIP labels on the CIP Validation Suite

**Reference:** the CIP Validation Suite, github.com/CIPValidationSuite/ValidationSuite. These
are the 300 structures of Hanson, Mayfield, Vainio, Yerin, Musacchio and Redkin, *J. Chem. Inf.
Model.* 2018, 58, 1755, with the labels the authors agreed.

**Method:** `python scripts/review/cip_review.py compounds_2d.sdf`.
- Each record is opened as a user would paste it (`convert.from_mol`, keeping its drawn
  coordinates and wedges).
- Each is labelled by `stereo.labels`.
- The labels are compared with the suite's tetrahedral (R/S/r/s) and double-bond (E/Z) labels.

**Verdict:**
- **The sketcher's drawing path loses nothing.** On every one of the 266 records with only
  tetrahedral and E/Z units, its labels are exactly those RDKit's labeller (`rdCIPLabeler`)
  gives on the suite's own molfile.
- **233 of 266 agree with the suite.** The 33 that differ are differences of RDKit's labeller
  itself, listed below.
  - Most are centres or bonds it does not perceive as stereogenic at all (no label):
    the rule 6 spiro cases and ring double bonds whose stereo comes from ring substituents.
  - The rest are rule 4 and 5 cases where it misses some centres, or adds labels the suite
    does not give (VS155, VS229).
  - These are upstream limits, recorded here, not hidden.
- **15 suite records are in `tests/sketch/test_stereo_checks.py`.** Their SMILES and labels are
  verbatim from `compounds.smi`. They need sequence rules 3, 4b and 5: Z beats E;
  like/unlike pairs; pseudo-asymmetric r/s.
- **Not labelled by the sketcher:** axial, helical and extended (allene, cumulene) units
  (34 records).

266 records with only tetrahedral and E/Z units: 233 agree, 33 differ. 34 skipped, 0 unreadable.

Of the 33 that differ, 33 give exactly what RDKit's labeller gives on the suite's own molfile.

| Record | Suite (atom + label) | Sketcher | RDKit on the suite molfile |
|---|---|---|---|
| VS132 | 13S 15S | (none) | (none) |
| VS155 | 6Z 7Z 8Z 9Z | 6Z 7Z 8P 9M 13M 15P 17M 19P | 6Z 7Z 8P 9M 13M 15P 17M 19P |
| VS188 | 5E 6E 7E 8E | (none) | (none) |
| VS189 | 7Z 8Z 9Z 10Z | (none) | (none) |
| VS190 | 7E 8E 9E 10E 11E 12E | (none) | (none) |
| VS207 | 16s 17s 18r 19r 20r 21r | 16s 17s 18r 19r | 16s 17s 18r 19r |
| VS214 | 13E 14s 15E 16s | (none) | (none) |
| VS226 | 25r 26r 27s 28s 29s 30s 31r 32r 33r 34r | 25r 26r 27s 28s 29s 30s 31r 32r | 25r 26r 27s 28s 29s 30s 31r 32r |
| VS229 | 8R 9S | 6e 7e 8R 9S | 6e 7e 8R 9S |
| VS246 | 9R | (none) | (none) |
| VS247 | 9R | (none) | (none) |
| VS248 | 9S | (none) | (none) |
| VS279 | 14R 15r 16S 17s 18S | 15r 16r 17s 18s | 15r 16r 17s 18s |
| VS280 | 16R 17R 18R 19R 20R | (none) | (none) |
| VS281 | 17S | (none) | (none) |
| VS282 | 17R | (none) | (none) |
| VS283 | 13S | (none) | (none) |
| VS284 | 13R | (none) | (none) |
| VS285 | 11R | (none) | (none) |
| VS286 | 15S | (none) | (none) |
| VS288 | 11S | (none) | (none) |
| VS289 | 7R | (none) | (none) |
| VS290 | 7S | (none) | (none) |
| VS291 | 9R | (none) | (none) |
| VS292 | 13S | (none) | (none) |
| VS293 | 2Z 5Z 8Z 13Z 14Z 15Z 16S | 2Z 5Z 8Z 13Z 14Z 15Z | 2Z 5Z 8Z 13Z 14Z 15Z |
| VS294 | 9S 10R 11R 12S 13R | 9S 10R 11R 12S | 9S 10R 11R 12S |
| VS295 | 9S 10R 11S 12R 13S | 9S 10R 11S 12R | 9S 10R 11S 12R |
| VS296 | 19S 20S 21S | (none) | (none) |
| VS297 | 17R 18S 19S | (none) | (none) |
| VS298 | 7R 8R 9S | (none) | (none) |
| VS299 | 11S | (none) | (none) |
| VS300 | 28s 29S 30R | (none) | (none) |

Skipped by stereo unit: AT 6, AT TH 1, CT CT4 1, CT TH5 1, CT4 2, CT4 TH 2, HE 2, TH TH3 4, TH3 4, TH5 1, no labels 10

## 2. Wedge warnings

Ketcher's `stereoValidator.ts` (commit `2b031cc2`) checks only neighbour counts and the parity
of implicit hydrogens, which neither the CTfile spec nor IUPAC states. So it is **not** used.

A wedge is instead judged by what RDKit's molfile reader makes of it, the same reader that
decides the molecule (`stereo.wedge_problems`):

| Warning | Rule |
|---|---|
| *Defines nothing* | RDKit reads no stereocentre (or atropisomer axis) at the narrow end, which is also not a possible stereocentre. |
| *Points the wrong way* | As above, but the wide end is a possible stereocentre. IUPAC 2006 (Brecher) puts the narrow end at the stereocentre. |
| *Contradict each other* | The atom is a possible stereocentre but RDKit reads no configuration. An example is a wedge and a hash on opposite bonds, where RDKit warns "opposing bonds have opposite wedging". |
| *Ambiguous* | Turning any one bond at the centre by 10° changes the configuration RDKit reads. This catches three bonds in a half plane, such as a T-shape with the wedge between two collinear plain bonds, and bonds drawn on top of each other. |

- Warnings are marked in amber, explained on hover and in the line under the drawing.
- They never block OK: the molecule is still what RDKit reads.
- Wavy bonds are never flagged.

## 3. "Either" double bonds

- A double bond can be marked *Unknown (crossed bond)*. It is written as CTfile bond stereo 3
  ("cis or trans") and read back the same way.
- RDKit reads it as `STEREOANY`, and SMILES carries no E/Z for it.
- A wavy single bond beside a double bond already means the same thing. RDKit's writer then
  also puts stereo 3 on the double bond, so on reading the sketch keeps only the wavy bond.

## 4. Enhanced stereo

- Stereocentres can be put in **ABS**, **AND n** (racemic: this and its mirror image) or
  **OR n** (one of the two, unknown which) groups. These follow the V3000
  `MDLV30/STEABS`, `STERAC` and `STEREL` collections, as RDKit's `StereoGroup` holds them.
- Only atoms RDKit reads as stereocentres are put in a group.
- A drawing with groups is written as V3000. Explicit valences are kept there as `VAL=`,
  e.g. `VAL=6` on PF₆⁻.
- CXSMILES `|&1:…,o1:…|` reads in. RDKit renumbers groups within each type when it writes.
