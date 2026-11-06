# Sketcher P2: salts and solvents review (2026-11-06)

**Source list:** the 135 records of Ketcher's `salts-and-solvents.sdf` (commit `2b031cc2`).

**What ships:**
- Ketcher's file is **not** shipped. `chemsuite/sketch/salts.py` writes each structure out
  afresh as SMILES, in groups for the menu, with full names and the short names chemists use
  (DMF, THF, TFA).
- Each entry has a `pubchem` lookup name: its name, a synonym or a CAS number.

**Method:** `python scripts/review/salts_review.py KETCHER_SALTS_SDF`.
- Each entry is looked up on PubChem (PUG REST) by its lookup name and compared by InChIKey.
- Each Ketcher record is then matched to an entry by InChIKey.

## Verdict

- **All 135 entries agree with PubChem.**
  - 134 agree by full InChIKey.
  - cis-decalin agrees by connectivity, because PubChem's record for it (CAS 493-01-6) has no
    stereo. `tests/sketch/test_salts.py` checks it is cis: in an embedded, MMFF-optimised
    3D structure, the bridgehead H–C–C–H dihedral is under 90°, never anti.
- **Two Ketcher records are drawn wrongly; both are corrected:**
  - **DBU.** Ketcher's record is `C1=C2CCCCNN2CCC1`, a hydrazine with a C=C bond and an N–N
    bond, under the wrong name "1,5-diazabicyclo[5.4.0]undec-7-ene". DBU is
    1,8-diazabicyclo[5.4.0]undec-7-ene, an amidine: `C1CCC2=NCCCN2CC1`, C₉H₁₆N₂.
  - **Guanidine hydrochloride.** Ketcher draws the cation as H₃N⁺–C(=NH)NH₂, protonated on an
    amino nitrogen. The salt is guanidinium chloride, (H₂N)₂C=NH₂⁺ Cl⁻. InChI normalises
    both to the same key, so the review cannot flag this one; a SMARTS test pins it.
- **Renamed for clarity:**
  - "perfluorocyclic ether" is perfluoro(2-butyltetrahydrofuran) (CAS 335-36-4, matched by
    InChIKey).
  - Names follow current usage: methanesulfonic acid, tert-butanol, propionitrile and so on.
- **Disodium EDTA** is looked up by CAS 139-33-3. Its common names resolve on PubChem to the
  dihydrate.

## Full results

| Entry | chemsuite SMILES | PubChem name | Verdict |
|---|---|---|---|
| acetic acid | `CC(=O)O` | acetic acid | agrees (InChIKey) |
| acetic anhydride | `CC(=O)OC(C)=O` | acetic anhydride | agrees (InChIKey) |
| formic acid | `O=CO` | formic acid | agrees (InChIKey) |
| methanesulfonic acid | `CS(=O)(=O)O` | methanesulfonic acid | agrees (InChIKey) |
| propionic acid | `CCC(=O)O` | propionic acid | agrees (InChIKey) |
| chloroacetic acid | `OC(=O)CCl` | chloroacetic acid | agrees (InChIKey) |
| trichloroacetic acid | `OC(=O)C(Cl)(Cl)Cl` | trichloroacetic acid | agrees (InChIKey) |
| trifluoroacetic acid | `OC(=O)C(F)(F)F` | trifluoroacetic acid | agrees (InChIKey) |
| 1,2-propanediol | `CC(O)CO` | 1,2-propanediol | agrees (InChIKey) |
| 1,3-propanediol | `OCCCO` | 1,3-propanediol | agrees (InChIKey) |
| 1,4-butanediol | `OCCCCO` | 1,4-butanediol | agrees (InChIKey) |
| 1-butanol | `CCCCO` | 1-butanol | agrees (InChIKey) |
| 1-propanol | `CCCO` | 1-propanol | agrees (InChIKey) |
| 2-butanol | `CCC(C)O` | 2-butanol | agrees (InChIKey) |
| 2-ethylhexanol | `CCCCC(CC)CO` | 2-ethyl-1-hexanol | agrees (InChIKey) |
| isopropanol | `CC(C)O` | isopropanol | agrees (InChIKey) |
| 2-methoxyethanol | `COCCO` | 2-methoxyethanol | agrees (InChIKey) |
| 2-pentanol | `CCCC(C)O` | 2-pentanol | agrees (InChIKey) |
| benzyl alcohol | `OCc1ccccc1` | benzyl alcohol | agrees (InChIKey) |
| cyclohexanol | `OC1CCCCC1` | cyclohexanol | agrees (InChIKey) |
| ethanol | `CCO` | ethanol | agrees (InChIKey) |
| ethylene glycol | `OCCO` | ethylene glycol | agrees (InChIKey) |
| glycerol | `OCC(O)CO` | glycerol | agrees (InChIKey) |
| isoamyl alcohol | `CC(C)CCO` | isoamyl alcohol | agrees (InChIKey) |
| isobutanol | `CC(C)CO` | isobutanol | agrees (InChIKey) |
| methanol | `CO` | methanol | agrees (InChIKey) |
| tert-butanol | `CC(C)(C)O` | tert-butanol | agrees (InChIKey) |
| 2,2,2-trifluoroethanol | `OCC(F)(F)F` | 2,2,2-trifluoroethanol | agrees (InChIKey) |
| water | `O` | water | agrees (InChIKey) |
| heavy water | `[2H]O[2H]` | deuterium oxide | agrees (InChIKey) |
| carbon dioxide | `O=C=O` | carbon dioxide | agrees (InChIKey) |
| carbon disulfide | `S=C=S` | carbon disulfide | agrees (InChIKey) |
| benzene | `c1ccccc1` | benzene | agrees (InChIKey) |
| cumene | `CC(C)c1ccccc1` | cumene | agrees (InChIKey) |
| mesitylene | `Cc1cc(C)cc(C)c1` | mesitylene | agrees (InChIKey) |
| toluene | `Cc1ccccc1` | toluene | agrees (InChIKey) |
| o-xylene | `Cc1ccccc1C` | o-xylene | agrees (InChIKey) |
| m-xylene | `Cc1cccc(C)c1` | m-xylene | agrees (InChIKey) |
| p-xylene | `Cc1ccc(C)cc1` | p-xylene | agrees (InChIKey) |
| N,N-dimethylaniline | `CN(C)c1ccccc1` | N,N-dimethylaniline | agrees (InChIKey) |
| pyridine | `c1ccncc1` | pyridine | agrees (InChIKey) |
| triethylamine | `CCN(CC)CC` | triethylamine | agrees (InChIKey) |
| N,N-diisopropylethylamine | `CCN(C(C)C)C(C)C` | N,N-diisopropylethylamine | agrees (InChIKey) |
| 1,8-diazabicyclo[5.4.0]undec-7-ene | `C1CCC2=NCCCN2CC1` | 1,8-diazabicyclo[5.4.0]undec-7-ene | agrees (InChIKey) |
| tris(hydroxymethyl)aminomethane | `NC(CO)(CO)CO` | tromethamine | agrees (InChIKey) |
| acetonitrile | `CC#N` | acetonitrile | agrees (InChIKey) |
| N,N-dimethylacetamide | `CC(=O)N(C)C` | N,N-dimethylacetamide | agrees (InChIKey) |
| N,N-dimethylformamide | `CN(C)C=O` | N,N-dimethylformamide | agrees (InChIKey) |
| dimethyl sulfoxide | `CS(C)=O` | dimethyl sulfoxide | agrees (InChIKey) |
| N,N′-dimethylpropyleneurea | `CN1CCCN(C)C1=O` | 1,3-dimethyl-3,4,5,6-tetrahydro-2(1H)-pyrimidinone | agrees (InChIKey) |
| formamide | `NC=O` | formamide | agrees (InChIKey) |
| N-methyl-2-pyrrolidone | `CN1CCCC1=O` | N-methyl-2-pyrrolidone | agrees (InChIKey) |
| N-methylformamide | `CNC=O` | N-methylformamide | agrees (InChIKey) |
| propionitrile | `CCC#N` | propionitrile | agrees (InChIKey) |
| sulfolane | `O=S1(=O)CCCC1` | sulfolane | agrees (InChIKey) |
| hexamethylphosphoramide | `CN(C)P(=O)(N(C)C)N(C)C` | hexamethylphosphoramide | agrees (InChIKey) |
| hexamethylphosphorous triamide | `CN(C)P(N(C)C)N(C)C` | hexamethylphosphorous triamide | agrees (InChIKey) |
| nitromethane | `C[N+](=O)[O-]` | nitromethane | agrees (InChIKey) |
| amyl acetate | `CCCCCOC(C)=O` | pentyl acetate | agrees (InChIKey) |
| butyl acetate | `CCCCOC(C)=O` | butyl acetate | agrees (InChIKey) |
| dimethyl carbonate | `COC(=O)OC` | dimethyl carbonate | agrees (InChIKey) |
| ethyl acetate | `CCOC(C)=O` | ethyl acetate | agrees (InChIKey) |
| ethyl formate | `CCOC=O` | ethyl formate | agrees (InChIKey) |
| ethyl lactate | `CCOC(=O)C(C)O` | ethyl lactate | agrees (InChIKey) |
| ethyl propionate | `CCOC(=O)CC` | ethyl propionate | agrees (InChIKey) |
| ethylene carbonate | `O=C1OCCO1` | ethylene carbonate | agrees (InChIKey) |
| isobutyl acetate | `CC(=O)OCC(C)C` | isobutyl acetate | agrees (InChIKey) |
| isopropyl acetate | `CC(=O)OC(C)C` | isopropyl acetate | agrees (InChIKey) |
| methyl acetate | `COC(C)=O` | methyl acetate | agrees (InChIKey) |
| methyl formate | `COC=O` | methyl formate | agrees (InChIKey) |
| methyl lactate | `COC(=O)C(C)O` | methyl lactate | agrees (InChIKey) |
| octyl acetate | `CCCCCCCCOC(C)=O` | octyl acetate | agrees (InChIKey) |
| propyl acetate | `CCCOC(C)=O` | propyl acetate | agrees (InChIKey) |
| propylene carbonate | `CC1COC(=O)O1` | propylene carbonate | agrees (InChIKey) |
| tert-butyl acetate | `CC(=O)OC(C)(C)C` | tert-butyl acetate | agrees (InChIKey) |
| 1,2-dimethoxyethane | `COCCOC` | 1,2-dimethoxyethane | agrees (InChIKey) |
| 1,4-dioxane | `C1COCCO1` | 1,4-dioxane | agrees (InChIKey) |
| 2-methyltetrahydrofuran | `CC1CCCO1` | 2-methyltetrahydrofuran | agrees (InChIKey) |
| anisole | `COc1ccccc1` | anisole | agrees (InChIKey) |
| butyl carbitol | `CCCCOCCOCCO` | diethylene glycol monobutyl ether | agrees (InChIKey) |
| cyclopentyl methyl ether | `COC1CCCC1` | cyclopentyl methyl ether | agrees (InChIKey) |
| dibutyl ether | `CCCCOCCCC` | dibutyl ether | agrees (InChIKey) |
| diethyl ether | `CCOCC` | diethyl ether | agrees (InChIKey) |
| diethylene glycol | `OCCOCCO` | diethylene glycol | agrees (InChIKey) |
| diglyme | `COCCOCCOC` | diglyme | agrees (InChIKey) |
| diisopropyl ether | `CC(C)OC(C)C` | diisopropyl ether | agrees (InChIKey) |
| dimethyl ether | `COC` | dimethyl ether | agrees (InChIKey) |
| diphenyl ether | `c1ccc(Oc2ccccc2)cc1` | diphenyl ether | agrees (InChIKey) |
| ethoxybenzene | `CCOc1ccccc1` | ethoxybenzene | agrees (InChIKey) |
| methyl tert-butyl ether | `COC(C)(C)C` | methyl tert-butyl ether | agrees (InChIKey) |
| tert-amyl methyl ether | `CCC(C)(C)OC` | tert-amyl methyl ether | agrees (InChIKey) |
| ethyl tert-butyl ether | `CCOC(C)(C)C` | ethyl tert-butyl ether | agrees (InChIKey) |
| tetrahydrofuran | `C1CCOC1` | tetrahydrofuran | agrees (InChIKey) |
| triethylene glycol | `OCCOCCOCCO` | triethylene glycol | agrees (InChIKey) |
| 1,2,4-trichlorobenzene | `Clc1ccc(Cl)c(Cl)c1` | 1,2,4-trichlorobenzene | agrees (InChIKey) |
| 1,2-dichlorobenzene | `Clc1ccccc1Cl` | 1,2-dichlorobenzene | agrees (InChIKey) |
| 1,2-dichloroethane | `ClCCCl` | 1,2-dichloroethane | agrees (InChIKey) |
| 2,3,4-trifluorotoluene | `Cc1ccc(F)c(F)c1F` | 2,3,4-trifluorotoluene | agrees (InChIKey) |
| carbon tetrachloride | `ClC(Cl)(Cl)Cl` | carbon tetrachloride | agrees (InChIKey) |
| chlorobenzene | `Clc1ccccc1` | chlorobenzene | agrees (InChIKey) |
| chloroform | `ClC(Cl)Cl` | chloroform | agrees (InChIKey) |
| dichloromethane | `ClCCl` | dichloromethane | agrees (InChIKey) |
| fluorobenzene | `Fc1ccccc1` | fluorobenzene | agrees (InChIKey) |
| perfluoro(2-butyltetrahydrofuran) | `FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)OC(F)(F)C(F)(F)C1(F)F` | 335-36-4 | agrees (InChIKey) |
| perfluorocyclohexane | `FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F` | dodecafluorocyclohexane | agrees (InChIKey) |
| perfluorohexane | `FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F` | perfluorohexane | agrees (InChIKey) |
| perfluorotoluene | `Fc1c(F)c(F)c(C(F)(F)F)c(F)c1F` | octafluorotoluene | agrees (InChIKey) |
| trichloroacetonitrile | `N#CC(Cl)(Cl)Cl` | trichloroacetonitrile | agrees (InChIKey) |
| (trifluoromethyl)benzene | `FC(F)(F)c1ccccc1` | benzotrifluoride | agrees (InChIKey) |
| 2-methylpentane | `CCCC(C)C` | 2-methylpentane | agrees (InChIKey) |
| cis-decalin | `C1CC[C@H]2CCCC[C@H]2C1` | 493-01-6 | agrees (connectivity; PubChem record has no stereo) |
| cyclohexane | `C1CCCCC1` | cyclohexane | agrees (InChIKey) |
| isooctane | `CC(C)CC(C)(C)C` | 2,2,4-trimethylpentane | agrees (InChIKey) |
| methylcyclohexane | `CC1CCCCC1` | methylcyclohexane | agrees (InChIKey) |
| methylcyclopentane | `CC1CCCC1` | methylcyclopentane | agrees (InChIKey) |
| heptane | `CCCCCCC` | heptane | agrees (InChIKey) |
| hexane | `CCCCCC` | hexane | agrees (InChIKey) |
| pentane | `CCCCC` | pentane | agrees (InChIKey) |
| 2-pentanone | `CCCC(C)=O` | 2-pentanone | agrees (InChIKey) |
| 3-pentanone | `CCC(=O)CC` | 3-pentanone | agrees (InChIKey) |
| acetone | `CC(C)=O` | acetone | agrees (InChIKey) |
| cyclohexanone | `O=C1CCCCC1` | cyclohexanone | agrees (InChIKey) |
| cyclopentanone | `O=C1CCCC1` | cyclopentanone | agrees (InChIKey) |
| methyl ethyl ketone | `CCC(C)=O` | 2-butanone | agrees (InChIKey) |
| methyl isobutyl ketone | `CC(=O)CC(C)C` | methyl isobutyl ketone | agrees (InChIKey) |
| TRIS hydrochloride | `[NH3+]C(CO)(CO)CO.[Cl-]` | tromethamine hydrochloride | agrees (InChIKey) |
| sodium chloride | `[Na+].[Cl-]` | sodium chloride | agrees (InChIKey) |
| sodium acetate | `CC(=O)[O-].[Na+]` | sodium acetate | agrees (InChIKey) |
| potassium chloride | `[K+].[Cl-]` | potassium chloride | agrees (InChIKey) |
| monopotassium phosphate | `OP(=O)(O)[O-].[K+]` | potassium dihydrogen phosphate | agrees (InChIKey) |
| lithium chloride | `[Li+].[Cl-]` | lithium chloride | agrees (InChIKey) |
| guanidine hydrochloride | `NC(N)=[NH2+].[Cl-]` | guanidine hydrochloride | agrees (InChIKey) |
| disodium phosphate | `OP(=O)([O-])[O-].[Na+].[Na+]` | sodium phosphate dibasic | agrees (InChIKey) |
| disodium EDTA | `OC(=O)CN(CCN(CC(=O)O)CC(=O)[O-])CC(=O)[O-].[Na+].[Na+]` | 139-33-3 | agrees (InChIKey) |
| ammonium sulfate | `[NH4+].[NH4+].[O-]S(=O)(=O)[O-]` | ammonium sulfate | agrees (InChIKey) |

135 entries: 134 agree in full, 1 by connectivity, 0 differ, 0 not found.

| Ketcher record | Ketcher structure | chemsuite entry |
|---|---|---|
| acetic acid | `CC(=O)O` | acetic acid |
| acetic anhydride | `CC(=O)OC(C)=O` | acetic anhydride |
| formic acid | `O=CO` | formic acid |
| methane sulphonic acid | `CS(=O)(=O)O` | methanesulfonic acid |
| propionic acid | `CCC(=O)O` | propionic acid |
| 1,2-propanediol | `CC(O)CO` | 1,2-propanediol |
| 1,3-propanediol | `OCCCO` | 1,3-propanediol |
| 1,4-butanediol | `OCCCCO` | 1,4-butanediol |
| 1-butanol | `CCCCO` | 1-butanol |
| 1-propanol | `CCCO` | 1-propanol |
| 2-butanol | `CCC(C)O` | 2-butanol |
| 2-ethylhexanol | `CCCCC(CC)CO` | 2-ethylhexanol |
| isopropanol | `CC(C)O` | isopropanol |
| 2-methoxyethanol | `COCCO` | 2-methoxyethanol |
| 2-pentanol | `CCCC(C)O` | 2-pentanol |
| benzyl alcohol | `OCc1ccccc1` | benzyl alcohol |
| cyclohexanol | `OC1CCCCC1` | cyclohexanol |
| ethanol | `CCO` | ethanol |
| ethylene glycol | `OCCO` | ethylene glycol |
| glycerol | `OCC(O)CO` | glycerol |
| isoamyl alcohol | `CC(C)CCO` | isoamyl alcohol |
| isobutanol | `CC(C)CO` | isobutanol |
| methanol | `CO` | methanol |
| t-butanol | `CC(C)(C)O` | tert-butanol |
| water | `O` | water |
| carbon dioxide | `O=C=O` | carbon dioxide |
| benzene | `c1ccccc1` | benzene |
| cumene | `CC(C)c1ccccc1` | cumene |
| mesitylene | `Cc1cc(C)cc(C)c1` | mesitylene |
| toluene | `Cc1ccccc1` | toluene |
| n,n-dimethylaniline | `CN(C)c1ccccc1` | N,N-dimethylaniline |
| pyridine | `c1ccncc1` | pyridine |
| triethylamine | `CCN(CC)CC` | triethylamine |
| acetonitrile | `CC#N` | acetonitrile |
| dimethyl acetamide | `CC(=O)N(C)C` | N,N-dimethylacetamide |
| dimethyl formamide | `CN(C)C=O` | N,N-dimethylformamide |
| dimethyl sulfoxide | `CS(C)=O` | dimethyl sulfoxide |
| dimethylpropylene urea | `CN1CCCN(C)C1=O` | N,N′-dimethylpropyleneurea |
| formamide | `NC=O` | formamide |
| n-methyl-2-pyrrolidone | `CN1CCCC1=O` | N-methyl-2-pyrrolidone |
| n-methylformamide | `CNC=O` | N-methylformamide |
| propane nitrile | `CCC#N` | propionitrile |
| sulfolane | `O=S1(=O)CCCC1` | sulfolane |
| amyl acetate | `CCCCCOC(C)=O` | amyl acetate |
| butyl acetate | `CCCCOC(C)=O` | butyl acetate |
| dimethyl carbonate | `COC(=O)OC` | dimethyl carbonate |
| ethyl acetate | `CCOC(C)=O` | ethyl acetate |
| ethyl formate | `CCOC=O` | ethyl formate |
| ethyl lactate | `CCOC(=O)C(C)O` | ethyl lactate |
| ethyl propionate | `CCOC(=O)CC` | ethyl propionate |
| ethylene carbonate | `O=C1OCCO1` | ethylene carbonate |
| isobutyl acetate | `CC(=O)OCC(C)C` | isobutyl acetate |
| isopropyl acetate | `CC(=O)OC(C)C` | isopropyl acetate |
| methyl acetate | `COC(C)=O` | methyl acetate |
| methyl formate | `COC=O` | methyl formate |
| methyl lactate | `COC(=O)C(C)O` | methyl lactate |
| n-octylacetate | `CCCCCCCCOC(C)=O` | octyl acetate |
| propyl acetate | `CCCOC(C)=O` | propyl acetate |
| propylene carbonate | `CC1COC(=O)O1` | propylene carbonate |
| t-butyl acetate | `CC(=O)OC(C)(C)C` | tert-butyl acetate |
| 1,2-dimethoxyethane | `COCCOC` | 1,2-dimethoxyethane |
| 1,4-dioxane | `C1COCCO1` | 1,4-dioxane |
| 2-methyl tetrahydrofuran | `CC1CCCO1` | 2-methyltetrahydrofuran |
| anisole | `COc1ccccc1` | anisole |
| butyl carbitol | `CCCCOCCOCCO` | butyl carbitol |
| cyclopentyl methyl ether | `COC1CCCC1` | cyclopentyl methyl ether |
| dibutyl ether | `CCCCOCCCC` | dibutyl ether |
| diethyl ether | `CCOCC` | diethyl ether |
| diethylene glycol | `OCCOCCO` | diethylene glycol |
| diethylene glycol dimethyl ether | `COCCOCCOC` | diglyme |
| diisopropyl ether | `CC(C)OC(C)C` | diisopropyl ether |
| dimethyl ether | `COC` | dimethyl ether |
| diphenyl ether | `c1ccc(Oc2ccccc2)cc1` | diphenyl ether |
| ethoxybenzene | `CCOc1ccccc1` | ethoxybenzene |
| methyl t-butyl ether | `COC(C)(C)C` | methyl tert-butyl ether |
| t-amyl methyl ether | `CCC(C)(C)OC` | tert-amyl methyl ether |
| t-butylethyl ether | `CCOC(C)(C)C` | ethyl tert-butyl ether |
| tetrahydrofuran | `C1CCOC1` | tetrahydrofuran |
| triethylene glycol | `OCCOCCOCCO` | triethylene glycol |
| 1,2,4-trichlorobenzene | `Clc1ccc(Cl)c(Cl)c1` | 1,2,4-trichlorobenzene |
| 1,2-dichlorobenzene | `Clc1ccccc1Cl` | 1,2-dichlorobenzene |
| 1,2-dichloroethane | `ClCCCl` | 1,2-dichloroethane |
| 2,2,2-trifluoroethanol | `OCC(F)(F)F` | 2,2,2-trifluoroethanol |
| 2,3,4-trifluorotoluene | `Cc1ccc(F)c(F)c1F` | 2,3,4-trifluorotoluene |
| carbon tetrachloride | `ClC(Cl)(Cl)Cl` | carbon tetrachloride |
| chloroacetic acid | `O=C(O)CCl` | chloroacetic acid |
| chlorobenzene | `Clc1ccccc1` | chlorobenzene |
| chloroform | `ClC(Cl)Cl` | chloroform |
| dichloromethane | `ClCCl` | dichloromethane |
| fluorobenzene | `Fc1ccccc1` | fluorobenzene |
| perfluorocyclic ether | `FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)OC(F)(F)C(F)(F)C1(F)F` | perfluoro(2-butyltetrahydrofuran) |
| perfluorocyclohexane | `FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F` | perfluorocyclohexane |
| perfluorohexane | `FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F` | perfluorohexane |
| perfluorotoluene | `Fc1c(F)c(F)c(C(F)(F)F)c(F)c1F` | perfluorotoluene |
| trichloroacetic acid | `O=C(O)C(Cl)(Cl)Cl` | trichloroacetic acid |
| trichloroacetonitrile | `N#CC(Cl)(Cl)Cl` | trichloroacetonitrile |
| trifluoroacetic acid | `O=C(O)C(F)(F)F` | trifluoroacetic acid |
| trifluoromethylbenzene | `FC(F)(F)c1ccccc1` | (trifluoromethyl)benzene |
| 2-methylpentane | `CCCC(C)C` | 2-methylpentane |
| cis-decalin | `C1CC[C@H]2CCCC[C@H]2C1` | cis-decalin |
| cyclohexane | `C1CCCCC1` | cyclohexane |
| isooctane | `CC(C)CC(C)(C)C` | isooctane |
| methylcyclohexane | `CC1CCCCC1` | methylcyclohexane |
| methylcyclopentane | `CC1CCCC1` | methylcyclopentane |
| heptane | `CCCCCCC` | heptane |
| hexane | `CCCCCC` | hexane |
| pentane | `CCCCC` | pentane |
| 2-pentanone | `CCCC(C)=O` | 2-pentanone |
| 3-pentanone | `CCC(=O)CC` | 3-pentanone |
| acetone | `CC(C)=O` | acetone |
| cyclohexanone | `O=C1CCCCC1` | cyclohexanone |
| cyclopentanone | `O=C1CCCC1` | cyclopentanone |
| methyl ethyl ketone | `CCC(C)=O` | methyl ethyl ketone |
| methyl isobutyl ketone | `CC(=O)CC(C)C` | methyl isobutyl ketone |
| carbon disulfide | `S=C=S` | carbon disulfide |
| nitromethane | `C[N+](=O)[O-]` | nitromethane |
| o-xylene | `Cc1ccccc1C` | o-xylene |
| m-xylene | `Cc1cccc(C)c1` | m-xylene |
| p-xylene | `Cc1ccc(C)cc1` | p-xylene |
| n,n-diisopropylethylamine | `CCN(C(C)C)C(C)C` | N,N-diisopropylethylamine |
| 1,5-diazabicyclo[5.4.0]undec-7-ene | `C1=C2CCCCNN2CCC1` | **none** |
| hexamethylphosphoramide | `CN(C)P(=O)(N(C)C)N(C)C` | hexamethylphosphoramide |
| hexamethylphosphorous triamide | `CN(C)P(N(C)C)N(C)C` | hexamethylphosphorous triamide |
| heavy water | `[2H]O[2H]` | heavy water |
| TRIS | `NC(CO)(CO)CO` | tris(hydroxymethyl)aminomethane |
| TRIS HCl | `[Cl-].[NH3+]C(CO)(CO)CO` | TRIS hydrochloride |
| Sodium-chloride | `[Cl-].[Na+]` | sodium chloride |
| Sodium-acetate | `CC(=O)[O-].[Na+]` | sodium acetate |
| Potassium-chloride | `[Cl-].[K+]` | potassium chloride |
| Monopotassium phosphate | `O=P([O-])(O)O.[K+]` | monopotassium phosphate |
| Lithium-chloride | `[Cl-].[Li+]` | lithium chloride |
| Guanidine HCl | `N=C(N)[NH3+].[Cl-]` | guanidine hydrochloride |
| Disodium phosphate | `O=P([O-])([O-])O.[Na+].[Na+]` | disodium phosphate |
| Disodium EDTA | `O=C([O-])CN(CCN(CC(=O)O)CC(=O)O)CC(=O)[O-].[Na+].[Na+]` | disodium EDTA |
| Ammonium-sulfate | `O=S(=O)([O-])[O-].[NH4+].[NH4+]` | ammonium sulfate |
