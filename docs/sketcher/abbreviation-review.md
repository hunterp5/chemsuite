# Sketcher P2: abbreviations (2026-11-05)

**Question.** Where should the sketcher's abbreviations (Ph, Boc, CO₂Me …) come from?

The plan named Ketcher's `fg.sdf` as the source. P0 found one wrong record there (Mes). P2
checked every record and also RDKit's built-in list. Neither can ship as is.

**Reproduce:**

```
python scripts/review/abbreviation_review.py <ketcher>/packages/ketcher-react/src/templates/fg.sdf
```

## Findings

### Ketcher `fg.sdf` (commit 2b031cc2, 62 records)

**9 records cannot be used:**
- Et, tBu, PO2, PO3, PO4, SO3, SO4 and "ster" are unreadable by RDKit. The causes are an
  old-style CTfile `G` group line, a malformed counts line ("ster"), and invalid valence fields.
- TBDPS fails valence checks.
- SO2 is a sulfonyl drawn with no valence left for its attachment (`*[SH](=O)=O` once read).

**Wrong:**
- **Mes** is 3-methylbenzyl (C₈H₉ through a CH₂), not mesityl (2,4,6-trimethylphenyl).

**Against chemsuite's table:** 40 records agree, 1 differs (Mes), 3 are unusable (Et, tBu,
TBDPS), and 19 of our labels are not in Ketcher.

### RDKit default abbreviations (`rdAbbreviations.GetDefaultAbbreviations()`, 2026.09.1)

**Wrong:**
- **iBu** is `*C(C)CC`, butan-2-yl (sec-butyl). Isobutyl is 2-methylpropyl, `*CC(C)C`.
- **iPent** is `*C(C)CCC`, pentan-2-yl. Isopentyl is 3-methylbutyl. We have no entry for it.

**Against chemsuite's table:** 19 agree, 1 differs (iBu).

## Decision

chemsuite keeps **its own table** (`chemsuite/sketch/abbrev.py`, 63 entries).

Each entry has:
- **its IUPAC substituent name and the formula that name gives, written by hand.** A test
  checks every structure against that formula.
- **one attachment point**, on the atom the name implies. Tests pin the isomers that names
  commonly get wrong: iBu, sBu, Mes, p-Tol, THP at C2, Bpin.
- **a mirrored form** for bonds coming in from the right (MeO₂C, AcO, F₃C).

Ketcher's and RDKit's lists remain cross-checks: the review script reports every disagreement.
`fg.sdf` is **not** shipped, so no correction file is needed.
- Ketcher's `LICENSE` and `NOTICE` stay in `chemsuite/sketch/data/` because of `ctfile.py`.
- Its salts and solvents and its template library still need the same checks when they are
  added.

## How abbreviations behave

**Storage.** An abbreviation is real atoms plus a group record. The molecule, its valence, its
SMILES and its descriptors therefore never depend on whether the label is shown.

**Molfiles.** It is written as a CTfile superatom S-group: `STY SUP`, `SAL` atoms, `SBL`
crossing bond, `SMT` label, `SAP` attachment, and `SDS EXP` when expanded.
- ChemDraw, Marvin and Ketcher read these as abbreviations.
- RDKit reads the full structure.

**Typed labels:**
1. An exact abbreviation comes first, so `Ac`, `Pr` and `Ts` are acetyl, propyl and tosyl as
   in ChemDraw, not actinium, praseodymium and tennessine.
2. Then element symbols, so `Co` is cobalt.
3. Then abbreviations in any case.

Element labels may carry an isotope, the atom's own hydrogens and a charge: `13C`, `NH2`, `O-`.
