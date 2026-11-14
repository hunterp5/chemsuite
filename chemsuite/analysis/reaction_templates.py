# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.

"""Named reaction templates for Enumerate Products.

Each :class:`Template` is a reaction SMARTS with a name for every component, so a dialog can
ask for "Aryl halide" and "Boronic acid" rather than "Reactant 1" and "Reactant 2". The set
follows the robust medicinal-chemistry reactions of Hartenfeller et al., *J. Chem. Inf.
Model.* 2011, 51, 3093, with each SMARTS written here and checked against hand-worked
products (``tests/analysis/test_reaction_templates.py``).

Chemoselectivity lives in the SMARTS: :data:`AMINE` is an amine and not an amide, sulfonamide,
hydrazine or hydroxylamine nitrogen; :data:`OXO` is an aldehyde or ketone and not an acid,
ester or amide. Reagents that are not varied (TosMIC in a Van Leusen, azide in a tetrazole
synthesis) are written into the product rather than asked for.

A component's **role** is what to put in it. A template's components are in the order its
SMARTS lists them, which is the order the enumerator takes its building blocks.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from rdkit import Chem, rdBase

#: An amine nitrogen: not an amide, sulfonamide, hydrazine, hydroxylamine or imine nitrogen.
AMINE = "[N;H2,H1;!$(N[#6]=[O,S,N]);!$(N[S,P]=O);!$(N[#8]);!$(N[#7]);!$(N=*):{n}]"
#: A primary amine, by the same rules.
PRIMARY_AMINE = "[N;H2;!$(N[#6]=[O,S,N]);!$(N[S,P]=O);!$(N[#8]);!$(N[#7]);!$(N=*):{n}]"
#: An aldehyde or ketone carbon: not an acid, ester, amide or acyl halide.
OXO = "[C;!$(C[O,N,S,F,Cl,Br,I]);!$(C=[!O]):{n}]=[O]"
#: A carboxylic acid carbon.
ACID = "[C:{n}](=[O:{o}])[OH]"
#: A leaving group on an aromatic carbon for cross-coupling.
ARYL_LEAVING = "[$(Cl),$(Br),$(I),$(OS(=O)(=O)C(F)(F)F)]"
#: An alkyl leaving group for N- and O-alkylation.
ALKYL_LEAVING = "[$(Cl),$(Br),$(I),$(OS(=O)(=O)c1ccccc1),$(OS(=O)(=O)C)]"
#: A saturated carbon that can be alkylated from (no neighbouring heteroatom or double bond).
ALKYL = "[C;X4;!$(C[O,N,S]);!$(C=*):{n}]"
#: An aromatic carbon activated for SNAr: ortho or para to a ring nitrogen or a nitro group.
ACTIVATED_ARYL = "[c;$(c:n),$(c:c:c:n),$(c:c[N+](=O)[O-]),$(c:c:c:c[N+](=O)[O-]):{n}]"

#: The order groups are offered in.
GROUPS = ("Coupling", "Ring Formation", "Multicomponent", "Deprotection")


@dataclass(frozen=True)
class Template:
    """A named reaction template.

    *roles* names each component of *smarts*, in order. *source* cites where the reaction
    comes from, and *note* says what the template assumes (an implicit reagent, a
    regiochemistry).
    """

    id: str
    name: str
    group: str
    roles: tuple[str, ...]
    smarts: str
    source: str
    note: str = ""

    @property
    def components(self) -> int:
        return len(self.roles)

    def describe(self) -> str:
        """One line for a tooltip: what it makes, what it assumes, and where it comes from."""
        parts = [" + ".join(self.roles)]
        if self.note:
            parts.append(self.note)
        parts.append(self.source)
        return ". ".join(parts)


HARTENFELLER = "Hartenfeller et al., J. Chem. Inf. Model. 2011, 51, 3093"

TEMPLATES: tuple[Template, ...] = (
    # -- couplings --------------------------------------------------------------------------
    Template(
        id="amide",
        name="Amide coupling",
        group="Coupling",
        roles=("Carboxylic acid", "Amine"),
        smarts=f"{ACID.format(n=1, o=2)}.{AMINE.format(n=3)}>>[C:1](=[O:2])[N:3]",
        source=HARTENFELLER,
        note="The acid is coupled as activated (HATU, EDC); the amine must not be an amide N",
    ),
    Template(
        id="amide-acyl-chloride",
        name="Amide from acyl chloride",
        group="Coupling",
        roles=("Acyl chloride", "Amine"),
        smarts=f"[C:1](=[O:2])[Cl].{AMINE.format(n=3)}>>[C:1](=[O:2])[N:3]",
        source=HARTENFELLER,
    ),
    Template(
        id="sulfonamide",
        name="Sulfonamide coupling",
        group="Coupling",
        roles=("Sulfonyl chloride", "Amine"),
        smarts=(f"[S:1](=[O:2])(=[O:3])[Cl].{AMINE.format(n=4)}>>[S:1](=[O:2])(=[O:3])[N:4]"),
        source=HARTENFELLER,
    ),
    Template(
        id="urea",
        name="Urea formation",
        group="Coupling",
        roles=("Amine", "Isocyanate"),
        smarts=f"{AMINE.format(n=1)}.[N:2]=[C:3]=[O:4]>>[N:1][C:3](=[O:4])[N:2]",
        source=HARTENFELLER,
    ),
    Template(
        id="carbamate",
        name="Carbamate formation",
        group="Coupling",
        roles=("Chloroformate", "Amine"),
        smarts=(f"[C:1](=[O:2])([Cl])[O:3].{AMINE.format(n=4)}>>[C:1](=[O:2])([N:4])[O:3]"),
        source=HARTENFELLER,
    ),
    Template(
        id="ester",
        name="Esterification",
        group="Coupling",
        roles=("Carboxylic acid", "Alcohol"),
        smarts=(f"{ACID.format(n=1, o=2)}.[C;!$(C=[O,N,S]):3][OH]>>[C:1](=[O:2])[O][C:3]"),
        source=HARTENFELLER,
    ),
    Template(
        id="suzuki",
        name="Suzuki-Miyaura coupling",
        group="Coupling",
        roles=("Aryl halide or triflate", "Boronic acid or ester"),
        smarts=(f"[c:1]{ARYL_LEAVING}.[c:2][$(B(O)O),$(B1OC(C)(C)C(C)(C)O1)]>>[c:1][c:2]"),
        source="Miyaura & Suzuki, Chem. Rev. 1995, 95, 2457",
        note="Aryl-aryl only; boronic acids and pinacol esters both match",
    ),
    Template(
        id="buchwald",
        name="Buchwald-Hartwig amination",
        group="Coupling",
        roles=("Aryl halide or triflate", "Amine"),
        smarts=f"[c:1]{ARYL_LEAVING}.{AMINE.format(n=2)}>>[c:1][N:2]",
        source="Ruiz-Castillo & Buchwald, Chem. Rev. 2016, 116, 12564",
    ),
    Template(
        id="sonogashira",
        name="Sonogashira coupling",
        group="Coupling",
        roles=("Aryl halide or triflate", "Terminal alkyne"),
        smarts=f"[c:1]{ARYL_LEAVING}.[CH1:2]#[C:3]>>[c:1][C:2]#[C:3]",
        source="Chinchilla & Najera, Chem. Rev. 2007, 107, 874",
    ),
    Template(
        id="chan-lam",
        name="Chan-Lam amination",
        group="Coupling",
        roles=("Boronic acid", "Amine"),
        smarts=f"[c:1][$(B(O)O)].{AMINE.format(n=2)}>>[c:1][N:2]",
        source="Chan et al., Tetrahedron Lett. 1998, 39, 2933",
    ),
    Template(
        id="snar",
        name="SNAr amination",
        group="Coupling",
        roles=("Activated aryl halide", "Amine"),
        smarts=(
            f"{ACTIVATED_ARYL.format(n=1)}[F,Cl,Br,$(OS(=O)(=O)C(F)(F)F)]."
            f"{AMINE.format(n=2)}>>[c:1][N:2]"
        ),
        source=HARTENFELLER,
        note="The ring carbon must be ortho or para to a ring nitrogen or a nitro group",
    ),
    Template(
        id="reductive-amination",
        name="Reductive amination",
        group="Coupling",
        roles=("Aldehyde or ketone", "Amine"),
        smarts=f"{OXO.format(n=1)}.{AMINE.format(n=2)}>>[C:1][N:2]",
        source=HARTENFELLER,
    ),
    Template(
        id="n-alkylation",
        name="N-alkylation",
        group="Coupling",
        roles=("Alkyl halide or sulfonate", "Amine"),
        smarts=(f"{ALKYL.format(n=1)}{ALKYL_LEAVING}.{AMINE.format(n=2)}>>[C:1][N:2]"),
        source=HARTENFELLER,
    ),
    Template(
        id="ether",
        name="Williamson ether synthesis",
        group="Coupling",
        roles=("Phenol", "Alkyl halide or sulfonate"),
        smarts=(f"[c:1][OH].{ALKYL.format(n=2)}{ALKYL_LEAVING}>>[c:1][O][C:2]"),
        source=HARTENFELLER,
    ),
    Template(
        id="click",
        name="Azide-alkyne cycloaddition (CuAAC)",
        group="Coupling",
        roles=("Organic azide", "Terminal alkyne"),
        smarts=("[#6:1][N:2]=[N+:3]=[N-:4].[CH1:5]#[C:6]>>[#6:1][n:2]1[n+0:3][n+0:4][c:6][cH:5]1"),
        source="Rostovtsev et al., Angew. Chem. Int. Ed. 2002, 41, 2596",
        note="Copper-catalysed, so the 1,4-triazole only",
    ),
    # -- ring formations --------------------------------------------------------------------
    Template(
        id="hantzsch-thiazole",
        name="Hantzsch thiazole synthesis",
        group="Ring Formation",
        roles=("alpha-Haloketone", "Thioamide"),
        smarts=(
            "[C:1](=[O:2])[CH2:3][$(Cl),$(Br),$(I)].[NH2:4][C:5]=[S:6]>>[c:1]1[n:4][c:5][s:6][c:3]1"
        ),
        source=HARTENFELLER,
    ),
    Template(
        id="paal-knorr-pyrrole",
        name="Paal-Knorr pyrrole synthesis",
        group="Ring Formation",
        roles=("1,4-Diketone", "Primary amine"),
        smarts=(
            "[C:1](=[O:2])[C:3][C:4][C:5](=[O:6])."
            f"{PRIMARY_AMINE.format(n=7)}>>[c:1]1[c:3][c:4][c:5][n:7]1"
        ),
        source=HARTENFELLER,
    ),
    Template(
        id="pyrazole",
        name="Pyrazole synthesis",
        group="Ring Formation",
        roles=("1,3-Diketone", "Hydrazine"),
        smarts=(
            "[C:1](=[O:2])[CH2:3][C:4](=[O:5]).[NH2:6][NH1:7][*:8]>>"
            "[c:1]1[c:3][c:4][n:7]([*:8])[n:6]1"
        ),
        source=HARTENFELLER,
        note="A substituted hydrazine gives both regioisomers; the template writes one",
    ),
    Template(
        id="benzimidazole",
        name="Benzimidazole synthesis",
        group="Ring Formation",
        roles=("1,2-Diaminobenzene", "Carboxylic acid"),
        smarts=("[c:1]([NH2:2])[c:3]([NH2:4]).[C:5](=[O:6])[OH]>>[c:1]1[n:2][c:5][nH:4][c:3]1"),
        source=HARTENFELLER,
    ),
    Template(
        id="oxadiazole",
        name="1,2,4-Oxadiazole synthesis",
        group="Ring Formation",
        roles=("Amidoxime", "Carboxylic acid"),
        smarts=("[C:1](=[N:2][OH:3])[NH2:4].[C:5](=[O:6])[OH]>>[c:1]1[n:4][c:5][o:3][n:2]1"),
        source=HARTENFELLER,
    ),
    Template(
        id="pictet-spengler",
        name="Pictet-Spengler",
        group="Ring Formation",
        roles=("Tryptamine", "Aldehyde or ketone"),
        smarts=(
            "[cH1:1](:[nH:7]):[c:2][CH2:3][CH2:4][NH2:5]."
            f"{OXO.format(n=6)}>>[nH:7][c:1]1[c:2][CH2:3][CH2:4][NH:5][C:6]1"
        ),
        source=HARTENFELLER,
        note="Closes onto the indole 2-position, giving a tetrahydro-beta-carboline",
    ),
    Template(
        id="tetrazole",
        name="Tetrazole from a nitrile",
        group="Ring Formation",
        roles=("Nitrile",),
        smarts="[C:1]#[N:2]>>[c:1]1[n:2][nH][n][n]1",
        source=HARTENFELLER,
        note="Azide is the implicit reagent (NaN3 or TMSN3)",
    ),
    # -- multicomponent ---------------------------------------------------------------------
    Template(
        id="ugi",
        name="Ugi four-component",
        group="Multicomponent",
        roles=("Carboxylic acid", "Primary amine", "Aldehyde or ketone", "Isocyanide"),
        smarts=(
            f"{ACID.format(n=1, o=2)}.{PRIMARY_AMINE.format(n=3)}.{OXO.format(n=4)}."
            "[C-:5]#[N+:6][*:7]>>[C:1](=[O:2])[N:3][C:4][C+0:5](=O)[N+0;H1:6][*:7]"
        ),
        source="Ugi, Angew. Chem. Int. Ed. 1962, 1, 8",
    ),
    Template(
        id="passerini",
        name="Passerini three-component",
        group="Multicomponent",
        roles=("Carboxylic acid", "Aldehyde or ketone", "Isocyanide"),
        smarts=(
            f"{ACID.format(n=1, o=2)}.{OXO.format(n=3)}.[C-:4]#[N+:5][*:6]>>"
            "[C:1](=[O:2])[O][C:3][C+0:4](=O)[N+0;H1:5][*:6]"
        ),
        source="Passerini, Gazz. Chim. Ital. 1921, 51, 126",
    ),
    Template(
        id="biginelli",
        name="Biginelli dihydropyrimidinone",
        group="Multicomponent",
        roles=("Aldehyde", "beta-Ketoester", "Urea"),
        smarts=(
            f"{OXO.format(n=1)}.[C:5](=[O:6])[CH2:7][C:8](=[O:9])[O:10]."
            "[NH2:11][C:12](=[O:13])[NH2:14]>>"
            "[C:1]1[C:7]([C:8](=[O:9])[O:10])=[C:5][NH:11][C:12](=[O:13])[NH:14]1"
        ),
        source="Kappe, Acc. Chem. Res. 2000, 33, 879",
        note="Thiourea works in the flask but is not matched here",
    ),
    Template(
        id="gbb",
        name="Groebke-Blackburn-Bienayme",
        group="Multicomponent",
        roles=("2-Aminoazine", "Aldehyde", "Isocyanide"),
        smarts=(
            "[NH2:1][c:2][n;R:3].[CH1:4]=[O].[C-:5]#[N+:6][*:7]>>"
            "[n:1]1[c:4][c+0:5]([N+0;H1:6][*:7])[n:3][c:2]1"
        ),
        source="Groebke et al., Synlett 1998, 661",
        note="Gives the 3-amino-imidazo[1,2-a]azine",
    ),
    Template(
        id="van-leusen",
        name="Van Leusen imidazole",
        group="Multicomponent",
        roles=("Aldehyde or ketone", "Primary amine"),
        smarts=(f"{OXO.format(n=1)}.{PRIMARY_AMINE.format(n=2)}>>[c:1]1[n:2][cH][n][cH]1"),
        source="Van Leusen et al., Tetrahedron Lett. 1972, 13, 2369",
        note="TosMIC is the implicit third component; gives the 1,5-disubstituted imidazole",
    ),
    # -- deprotection -----------------------------------------------------------------------
    Template(
        id="boc-removal",
        name="Boc removal",
        group="Deprotection",
        roles=("Boc-protected amine",),
        smarts="[N:1][C](=O)[O][C]([CH3])([CH3])[CH3]>>[N:1]",
        source=HARTENFELLER,
    ),
    Template(
        id="ester-hydrolysis",
        name="Ester hydrolysis",
        group="Deprotection",
        roles=("Ester",),
        smarts="[C:1](=[O:2])[O][C;!$(C=*):3]>>[C:1](=[O:2])[OH]",
        source=HARTENFELLER,
    ),
)

_BY_ID = {t.id: t for t in TEMPLATES}


def by_id(template_id: str) -> Template | None:
    """The template with this id, or None."""
    return _BY_ID.get(template_id)


def in_group(group: str) -> tuple[Template, ...]:
    """The templates of *group*, in their listed order."""
    return tuple(t for t in TEMPLATES if t.group == group)


@cache
def role_queries(template_id: str) -> tuple[Chem.Mol, ...]:
    """Each component's reactant template as a query molecule, for :func:`fits`."""
    template = _BY_ID[template_id]
    with rdBase.BlockLogs():
        from rdkit.Chem import rdChemReactions

        rxn = rdChemReactions.ReactionFromSmarts(template.smarts)
        if rxn is None:
            return ()
        return tuple(rxn.GetReactants())


def fits(template_id: str, mol: Chem.Mol) -> list[int]:
    """Which components of this template *mol* could fill (indices into ``roles``)."""
    queries = role_queries(template_id)
    with rdBase.BlockLogs():
        return [i for i, query in enumerate(queries) if mol.HasSubstructMatch(query)]


def any_fits(template_id: str, component: int, mols: list[Chem.Mol]) -> bool:
    """Whether any of *mols* can fill this component."""
    queries = role_queries(template_id)
    if component >= len(queries):
        return False
    query = queries[component]
    with rdBase.BlockLogs():
        return any(mol is not None and mol.HasSubstructMatch(query) for mol in mols)


__all__ = [
    "GROUPS",
    "TEMPLATES",
    "Template",
    "any_fits",
    "by_id",
    "fits",
    "in_group",
    "role_queries",
]
