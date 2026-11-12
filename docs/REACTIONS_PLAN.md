# Plan: reactions in the sketcher and the app (drafted 2026-11-06)

## 1. Goal and ground rules

Chemists should be able to draw, store, view, search and use reactions as they do structures:

- draw a reaction in the sketcher (arrow, plus, reactants, agents, products, atom maps);
- open and save reactions: RXN files, RDF files, and reaction SMILES in any table file;
- keep reactions in a table column that draws them, edits in the sketcher and copies out;
- search reactions by a drawn reaction query, and split or apply them (library enumeration).

The sketcher's ground rules carry over unchanged (SKETCHER_PLAN §1):

- an original implementation;
- Ketcher (EPAM, Apache-2.0, commit `2b031cc2`) as the cited rule reference;
- nothing from MolManager's `iupac_*`;
- every rule verified against independent reference data before it ships.

## 2. Rules taken from references

| Rule | Source | Ours |
|---|---|---|
| Which drawn fragments are reactants and which are products | Ketcher `ketcher-core/src/domain/entities/struct.ts`, `getComponents` and `defineRxnFragmentTypeForAtomset` | **Grouping.** Fragments are grouped into components between the plus signs and the arrow, so a salt drawn without a plus is one component. A component is a reactant if its bounding-box centre is left of the arrow's centre, otherwise a product. **Extension:** a component whose centre lies within the arrow's horizontal span and above or below it is an **agent** (reagent, catalyst, solvent). Ketcher has no agents and saves such a component as a product. |
| One arrow per reaction | Ketcher `getComponents` ("there's just one arrow") | A drawing with more than one arrow is not a single reaction; the sketcher says so and will not save it as one. |
| RXN V2000 / V3000 layout | BIOVIA CTfile Formats (2020), "RXN" and "RDfile" | Written and read by RDKit (`rdChemReactions`); the sketcher writes component molfiles with its own coordinates. |
| Reaction SMILES and SMARTS, agents between the `>` | Daylight Theory Manual §5 (reaction SMILES); RDKit `ReactionFromSmarts(useSmiles=True)` | RDKit reads and writes. Atom maps are `:n` in brackets. |
| Atom-atom mapping numbers | CTfile atom block field 14 (aamap); `M  V30 … AAMAP=` | Each atom carries `mapping` (0 = none). |
| Arrow kinds | Ketcher `rxnArrow.ts` `RxnArrowMode` | **Kept:** plain (open angle), equilibrium, retrosynthetic and failed (crossed). **Not kept:** Ketcher's other drawing variants (half-bows, arcs). |

## 3. Phases

| Phase | Work | Done when |
|---|---|---|
| R1 ✅ | Qt-free reaction model. `Arrow` and `Plus` objects and atom `mapping` in `sketch/model.py`, undoable through `edits.py`. `sketch/reaction.py`: sorts components into reactants, agents and products (§2); converts to and from RDKit's `ChemicalReaction`, RXN (V2000/V3000) and reaction SMILES; lays out a reaction that is read in (components left to right, plus signs between, arrow, agents above the arrow); Clean for reactions. | Done 2026-11-06: of the 50,016 USPTO-50k reactions, 49,986 come back identical through drawing and RXN, 27 are the same molecules by InChIKey, and 3 hold explicit hydrogen counts on phosphorus that the sketcher cannot draw (docs/sketcher/reaction-review.md). |
| R2 ✅ | Sketcher UI. Arrow tool (drag; kinds in its menu), plus tool, atom-map tool (drag from an atom to an atom gives both the next number; Clear Mapping), drawing arrows, plus signs and map numbers; select, move, erase and rotate work on them; copy and paste RXN and reaction SMILES. OK gives a reaction when the drawing has an arrow. | Done 2026-11-06: UI tests drive each tool and the arrow menu; a pasted reaction (RXN or reaction SMILES) copies back out unchanged. Layout leaves room for heteroatom labels, and map numbers keep clear of them. |
| R3 ✅ | Reaction columns. A `reaction` column kind (text plus a parsed RDKit reaction blob): read reaction SMILES and RXN text, draw the reaction in the cell, edit in the sketcher, Copy As (RXN, reaction SMILES, image), save in documents. Open and save RXN and RDF files; detect reaction SMILES in CSV/TSV/XLSX; Convert Column to Reaction. | Done 2026-11-06: all 50,016 USPTO-50k reactions open from CSV (the `rxn_smiles` column is found by content), parse with none invalid, and come back identical through a document, CSV and Excel. Through an RD file (written with fresh layouts, 84 s on all cores) 49,986 come back identical and the same 30 as in R1 differ (27 aromatic forms, 3 P–H counts); every data field comes back identical. Cells draw and edit in the sketcher. |
| R4 ✅ | Using reactions. **Reaction Search** (substructure, with a drawn reaction query; RDKit `HasReactionSubstructMatch`). **Split Reaction** into reactant, agent and product structure columns. **Enumerate Products** applies a reaction (with a drawn or stored template) to building-block columns or files (library enumeration; `RunReactants`), keeping partial results on cancel. Filter cards for reaction columns. User guide section. | Done 2026-11-07: each tool tested against hand-checked reactions (amide coupling, Suzuki, reductive amination, Boc removal; `tests/analysis/test_reactions.py`), UI tests drive search, split, both enumeration modes, a building-block file and the filter card (`tests/tools/reactions/test_reaction_tools_ui.py`); user guide sections Reaction Search, Split Reactions, Enumerate Products and Filtering by reaction. |
| Later | Automatic atom mapping (rxnmapper, an optional add-on: it needs PyTorch); reaction centre highlighting; multi-step schemes (several arrows). | — |

## 4. Open questions

- **Automatic mapping.** RDKit has no atom mapper. rxnmapper (IBM, MIT) is the usual choice
  but needs PyTorch; it would be an optional add-on like Uni-pKa. Manual mapping ships in R2.
- **Reaction columns in plots.** Reactions have no descriptors of their own; R4's Split
  Reaction gives structure columns that do.
