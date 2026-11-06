# Sketcher P0: interaction spec (2026-11-05, draft)

**Decision (user, 2026-11-05):**
- ChemDraw's keys for what chemists do every day: atom labels and bond orders over the
  hovered atom or bond, and Clean Structure.
- Ketcher's keys for everything else.

Ketcher's keys below come from its source at `2b031cc2`
(`ketcher-react/src/script/ui/action/{atoms,tools,index,templates}.ts`). ChemDraw's come only
from public descriptions of its hotkeys, because the official cheat sheet is not publicly
fetchable; **open item 1** checks them against the team's installed ChemDraw.

## 1. Hover keys: act on the atom or bond under the pointer

This is how ChemDraw and Ketcher both work: point at an atom or bond, then press a key.

| Key | Over an atom | Over a bond | Source |
|---|---|---|---|
| `c` `n` `o` `s` `p` `f` `h` `i` | Set the element (C, N, O, S, P, F, H, I) | — | ChemDraw convention; Ketcher `atomCuts` |
| `l` / `b` | Cl / Br | — | Ketcher `atomCuts` (ChemDraw: open item 1) |
| `Shift+s` / `Shift+n` / `Shift+b` | Si / Na / B | — | Ketcher |
| `a` `q` `x` `m` `r` | Generic atoms A, Q, X (halogen), M (metal); R-group | — | Ketcher |
| `d` | Deuterium | — | Ketcher |
| `1` `2` `3` | Add a single, double or triple bond in the best direction (§3) | Set the bond order | Bonds: ChemDraw ("type 3 on a bond: triple") and Ketcher `bondCuts` |
| `4` | — | Aromatic bond | Ketcher (ChemDraw uses 4 for quadruple; we follow Ketcher, since quadruple bonds do not occur in med-chem) |
| `0` | — | Any-bond query | Ketcher |
| `=` / `-` | Charge +1 / −1 | — | Ketcher |
| `Delete`, `Backspace` | Erase | Erase | Ketcher |

Without hovering, the same keys pick the matching tool, as in Ketcher.

## 2. Commands

| Keys (Mod = Ctrl or ⌘) | Action | Source |
|---|---|---|
| `Mod+Shift+K` | Clean structure | ChemDraw. Ketcher uses this key for "copy as KET", a format we do not have. |
| `Mod+Z` / `Mod+Shift+Z`, `Mod+Y` | Undo / redo | Both |
| `Mod+X` `Mod+C` `Mod+V` | Cut, copy, paste (molfile and SMILES; molfile formats from ChemDraw, Marvin and Biovia Draw, as today's editor reads them) | Both |
| `Mod+Shift+M` | Copy as molfile | Ketcher |
| `Mod+A` / `Mod+Shift+A` | Select all / none | Ketcher |
| `Escape`, `Shift+Tab` | Select tool; Escape also cancels the current drag | Ketcher |
| `Alt+H` / `Alt+V` | Flip the selection horizontally / vertically | Ketcher |
| `Mod+G` | S-group (abbreviation, superatom) | Ketcher |
| `t` / `Shift+T` | Cycle the ring templates / open the template library | Ketcher |
| `Alt+E` | Enhanced stereo | Ketcher |

## 3. Mouse behaviour

These are re-implemented from Ketcher (`ketcher-react/src/script/editor/tool/*.ts`), each with a
test of its own.

**Bond tool**

| Action | Result | Ketcher rule |
|---|---|---|
| Click empty canvas | New bond: a single bond at −30°, other types horizontal | `bond.ts` |
| Click an atom | New bond in the largest gap. Zig-zag after one neighbour. 180° when a triple or cumulated double bond needs it. Reuses an atom that lands within 0.1 bond lengths. | `atomForNewBond` |
| Click a bond | Set it to the tool's type. With the single-bond tool, a plain bond cycles single → double → triple. Clicking a wedge or hash with the same tool flips its direction. | `bondChangingAction` (`ketcher-core/…/actions/bond.ts`) |
| Drag from an atom | Bond follows the pointer, snapping to fixed angle steps (Ctrl frees it; the step size comes from Ketcher's `fracAngle` in P1). Releasing on an atom bonds to it. | `bond.ts` |

**Chain tool**
- Dragging draws a zig-zag chain. One bond is added per bond length dragged, and the bond count is shown while dragging. Ctrl frees the angle. (`chain.ts`)

**Ring and template tool**

| Action | Result | Ketcher rule |
|---|---|---|
| Click empty canvas | Place the template | `template.ts` |
| Click an atom | Spiro attachment, or a substituent for templates with an attachment point | `template.ts` |
| Click a bond | Fuse on that bond | `closelyFusing.ts`, `aromaticFusing.ts` |
| Drag | Rotate before placing | `template.ts` |

**Atom tool**
- Click an atom: set its element.
- Click empty canvas: a lone atom.
- Double-click an atom: edit its label as text (element, abbreviation such as `CO2Me`,
  isotope, charge).

**Select tool**
- Rectangle or lasso selection; click a fragment to take all of it.
- Dragging the selection merges atoms that land on other atoms (`atomMerge.ts`).
- Rotate handle.

**Hover**
- The atom or bond under the pointer is highlighted.
- A bad valence (§`valence-review.md`) is marked with the reason in a tooltip.

## 4. Open items

1. **ChemDraw keys.** Check the hover letters (`l`, `b`, and any letters ChemDraw uses for
   abbreviations such as Me, Et and Ph) against the Hotkeys Cheat Sheet of the team's ChemDraw
   version, and settle clashes in favour of ChemDraw.
2. **Angle-snap step.** Read the step from Ketcher's `fracAngle` (P1).
3. **Which ChemDraw abbreviation hotkeys to add** (for example Shift+F for CF₃ in recent
   versions). This depends on item 1.
