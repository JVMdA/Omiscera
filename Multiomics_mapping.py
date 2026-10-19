"""Multi-omics - Cross-omics Pathway Mapping (second tile of the Multi-omics submodule).

QUESTION ANSWERED
-----------------
Joint Pathway Enrichment says WHICH pathways the omics point at together. This
tile goes inside ONE pathway and shows which of its members every omics hit, and
how they connect: enzyme (gene / protein) -> reaction -> metabolite. It is the
established "paint several omics on one pathway" view (PaintOmics 4, Liu et al.
2022, NAR 50:W551; Pathview, Luo & Brouwer 2013, Bioinformatics 29:1830; KEGG
Mapper, Kanehisa & Sato 2020, Protein Sci 29:28). It is DESCRIPTIVE: it computes
no p-value. Which members count as "hit" is decided by the same rule as the
joint test (a plain list is the hit set; for a differential-expression table
Omiscera applies the padj / |log2FC| cutoffs, up and down pooled).

WHERE THE PATHWAY STRUCTURE COMES FROM
--------------------------------------
KEGG: the pathway's KGML file (genes with their reactions, compounds, KEGG's own
layout). Reactome: the pathway's diagram layout plus its entity graph (a
complex or set is opened down to its UniProt / ChEBI leaves). MetaCyc is not
offered: since 2024 it needs a paid BioCyc subscription. Fetching and caching
live in Multiomics_module_retriever; THIS module is pure (no network) so every
rule below is unit-testable offline.

DESIGN RULES (each one a decision, tested in tests/test_multiomics_mapping.py)
------------------------------------------------------------------------------
1. Three states per omics per node, never two. "hit" = significant / in the
   list; "measured" = present in the omics' table but not a hit (only an omics
   with a measured-reference notion - a differential-expression table - can say
   this; a plain list has no "measured, not significant"); "absent" = not in
   the omics' data (or not matched to this node). Conflating measured-not-hit
   with absent would hide that a member was looked at and did not change.
2. A node is hit by an omics if ANY of its members is. A multi-gene KEGG box or
   a Reactome complex lists every member with its own status.
3. Identifier bridges are the ones the rest of the app already uses and each
   match records which bridge it went through. A Reactome chemical is matched on
   its ChEBI ID first and, when Reactome cites a different protonation state of
   the same compound (it lists citrate(3-), CHEBI:16947; a cross-reference
   usually lists citric acid, CHEBI:30769), on the InChIKey SKELETON - the same
   bridge Metabolomics' own Reactome universe uses, labelled as such: UniProt -> KEGG gene
   (conv/uniprot, via Genomics_module.build_kegg_gene_map), Entrez -> KEGG gene,
   gene symbol -> KEGG gene (last resort), KEGG gene -> UniProt for a
   differential-expression table that has no UniProt column, and gene symbol ->
   a Reactome entity's own name only when no UniProt is known. An identifier
   that maps to MORE than one KEGG gene is never guessed (build_kegg_gene_map
   already drops it).
4. A Reactome protein / complex node is matched on its UniProt leaves only and
   a Reactome chemical node on its ChEBI leaves only. A cofactor inside a
   complex (e.g. a 4Fe-4S cluster in ACO2) must not make the complex look like
   a metabolite hit.
5. A "cross-omics reaction link" needs evidence at the reaction: its enzyme hit
   by a gene / protein omics AND a substrate or product hit by Metabolomics, or
   its enzyme hit by two or more gene / protein omics. Direction is reported
   for the differential-expression omics only (a compound list has none); the
   column says "concordant" / "discordant" only when both Transcriptomics and
   Proteomics are hits with a clear direction.
6. Nothing is padded or guessed: an omics left out of the run has no chips; a
   pathway whose structure could not be fetched says so and shows no map.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET

import pandas as pd

from Metabolomics_module import inchikey_skeleton  # read-only import: the first InChIKey block

OMICS_ORDER = ("genomics", "transcriptomics", "proteomics", "metabolomics")
GENE_OMICS = ("genomics", "transcriptomics", "proteomics")
OMICS_LABELS = {"genomics": "Genomics", "transcriptomics": "Transcriptomics",
                "proteomics": "Proteomics", "metabolomics": "Metabolomics"}
MAP_DATABASES = ("KEGG", "Reactome")

_REACTOME_COMPOUND_CLASSES = ("Chemical", "ChemicalDrug")
_REACTOME_GENE_CLASSES = ("Protein", "Complex", "EntitySet", "Gene", "RNA")
_REACTOME_LINK_CLASSES = ("Pathway", "ProcessNode")

MAP_CAVEAT = (
    "Pathway mapping is descriptive: it shows where each omics' hits fall in the pathway and computes no p-value. "
    "Significance comes from each omics' own cutoffs and from Joint Pathway Enrichment."
)


def finite_or_none(value):
    """JSON-safe number: NaN/inf/non-numeric -> None."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


# ---------------------------------------------------------------------------
# Identifier normalisation
# ---------------------------------------------------------------------------

def norm_uniprot(value) -> str:
    """UniProt accession without an isoform suffix ('Q99798-2' -> 'Q99798')."""
    s = str(value or "").strip().upper()
    return s.split("-")[0] if s else ""


def norm_chebi(value) -> str:
    """ChEBI number only ('CHEBI:30031' / 'chebi:30031' / '30031' -> '30031')."""
    m = re.search(r"(\d+)", str(value or ""))
    return m.group(1) if m else ""


def norm_kegg_compound(value) -> str:
    """KEGG compound C-number ('cpd:C00158' -> 'C00158'); '' for anything else (a glycan, a blank)."""
    m = re.search(r"\bC\d{5}\b", str(value or "").upper().replace("CPD:", " "))
    return m.group(0) if m else ""


# ---------------------------------------------------------------------------
# Pathway structure parsers (pure)
#
# Both return the SAME shape so the mapper and the chart never care which
# database a structure came from:
#   {"title", "nodes": [{"id","kind","type","label","x","y","w","h","members": {...}}],
#    "reactions": [{"id","name","reversible","enzymes","substrates","products","pos","segments"}]}
# kind is "gene" (a gene / protein / complex), "compound", "link" (a link to
# another pathway) or "other". y grows DOWNWARD in both databases' own layouts.
# ---------------------------------------------------------------------------

def _f(value, default=0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _polyline(coords) -> list:
    """KGML 'coords' ('x1,y1,x2,y2,...') -> [[x, y], ...] with consecutive duplicates removed (KGML repeats the
    corner points of its rounded lines)."""
    nums = [_f(v, None) for v in str(coords or "").split(",")]
    nums = [v for v in nums if v is not None]
    pts: list = []
    for i in range(0, len(nums) - 1, 2):
        pt = [nums[i], nums[i + 1]]
        if not pts or pts[-1] != pt:
            pts.append(pt)
    return pts


def _flat_ints(points) -> list:
    return [int(round(v)) for pt in points for v in pt]


def parse_kgml_full(xml_text: str) -> dict:
    """Parse a KEGG KGML file keeping what Metabolomics' own parser (parse_kegg_kgml, never edited) drops:
      - GENE entries and the reaction -> enzyme link;
      - the geometry of 'line' graphics. KEGG's global overview maps (Metabolic pathways, ...) draw a gene as a LINE
        along its reaction arrow, with no x / y at all: such a node keeps its polyline in node["path"] (flat ints,
        x1, y1, x2, y2, ...) and its x / y are the middle vertex. A line ORTHOLOG is decoration, returned in
        structure["decor_lines"] and never a node; structure["line_graphics"] says the map is of this kind;
      - GROUP entries (a drawn complex): kind 'group', members = the union of its components' genes;
      - RELATIONS (activation, inhibition, binding, phosphorylation, ...): the lines of a signalling map, which has
        no reactions at all."""
    root = ET.fromstring(xml_text)
    nodes, by_entry, decor, group_components = [], {}, [], {}
    for entry in root.findall("entry"):
        eid = entry.get("id")
        etype = (entry.get("type") or "").lower()
        names = (entry.get("name") or "").split()
        g = entry.find("graphics")
        label = ""
        x = y = w = h = 0.0
        path = None
        line_color = ""
        if g is not None:
            label = (g.get("name") or "").split(",")[0].strip().rstrip(".")
            x, y, w, h = _f(g.get("x")), _f(g.get("y")), _f(g.get("width")), _f(g.get("height"))
            if (g.get("type") or "").lower() == "line":
                pts = _polyline(g.get("coords"))
                if len(pts) >= 2:
                    path = _flat_ints(pts)
                    x, y = pts[len(pts) // 2]
                    line_color = (g.get("fgcolor") or "").upper()
        if etype == "map" and label.upper().startswith("TITLE:"):
            continue  # the map's own caption, not a link to another pathway
        if etype == "ortholog" and path is not None:
            decor.append(path)  # a grey line KEGG draws for an enzyme this organism lacks: shape only
            continue
        members: dict = {}
        if etype == "gene":
            kind = "gene"
            members["kegg_genes"] = [n for n in names if ":" in n and not n.startswith(("ko:", "path:", "cpd:", "gl:", "rn:"))]
        elif etype == "compound":
            kind = "compound"
            members["kegg_compounds"] = [c for c in (norm_kegg_compound(n) for n in names if n.startswith("cpd:")) if c]
        elif etype == "map":
            kind = "link"
        elif etype == "group":
            kind = "group"
            group_components[eid] = [c.get("id") for c in entry.findall("component")]
        else:
            kind = "other"
        node = {"id": f"k{eid}", "kind": kind, "type": etype, "label": label or (names[0] if names else ""),
                "x": x, "y": y, "w": w, "h": h, "members": members}
        if path is not None:
            node["path"] = path
            if line_color:
                node["color"] = line_color  # KEGG colours every line of a global map by its pathway CLASS
        if etype == "map":
            m = re.search(r"(\d{5})$", entry.get("name") or "")
            if m:
                node["pathway"] = m.group(1)  # the map number this label links to
        nodes.append(node)
        by_entry[eid] = node
    for eid, comps in group_components.items():
        parts = [by_entry[c] for c in comps if c in by_entry]
        node = by_entry[eid]
        node["components"] = [p["id"] for p in parts]
        genes: list = []
        for part in parts:
            for gene in part["members"].get("kegg_genes", []):
                if gene not in genes:
                    genes.append(gene)
        node["members"] = {"kegg_genes": genes}
        node["label"] = "+".join(p["label"] for p in parts[:4] if p["label"]) + ("+..." if len(parts) > 4 else "")
    reactions = []
    for rx in root.findall("reaction"):
        eid = rx.get("id")
        names = [n.replace("rn:", "") for n in (rx.get("name") or "").split()]
        subs = [f"k{s.get('id')}" for s in rx.findall("substrate") if s.get("id") in by_entry]
        prods = [f"k{p.get('id')}" for p in rx.findall("product") if p.get("id") in by_entry]
        reactions.append({"id": f"rx{len(reactions)}", "name": " ".join(names) or f"reaction {eid}",
                          "reversible": (rx.get("type") or "").lower() == "reversible",
                          "enzymes": [f"k{eid}"] if eid in by_entry else [],
                          "substrates": subs, "products": prods, "pos": None, "segments": []})
    relations = []
    for rel in root.findall("relation"):
        rtype = (rel.get("type") or "").lower()
        a, b = rel.get("entry1"), rel.get("entry2")
        if rtype == "maplink" or a not in by_entry or b not in by_entry:
            continue
        subs = [(sv.get("name") or "").lower() for sv in rel.findall("subtype")]
        subs = [t for t in subs if t and t not in ("compound", "hidden compound")]
        relations.append({"a": f"k{a}", "b": f"k{b}", "rtype": rtype, "types": subs or [rtype or "relation"]})
    return {"title": root.get("title") or "", "nodes": nodes, "reactions": reactions, "relations": relations,
            "decor_lines": decor, "line_graphics": any("path" in n for n in nodes)}


def _reactome_leaves(db_id, graph_by_id: dict, seen: "set | None" = None) -> list:
    """Every leaf entity under a diagram node's entity (a complex or set is opened
    through its children, recursively; cycles are guarded)."""
    seen = seen if seen is not None else set()
    if db_id in seen:
        return []
    seen.add(db_id)
    node = graph_by_id.get(db_id)
    if node is None:
        return []
    children = node.get("children") or []
    if not children:
        return [node]
    out: list = []
    for child in children:
        out.extend(_reactome_leaves(child, graph_by_id, seen))
    return out


def _strip_compartment(name: str) -> str:
    return re.sub(r"\s*\[[^\]]*\]\s*$", "", str(name or "")).strip()


def parse_reactome_diagram(diagram: dict, graph: dict, keep_reaction_db_ids: "set | None" = None) -> dict:
    """Combine a Reactome diagram layout (nodes + reaction edges) with its entity
    graph (each node's leaves, with UniProt / ChEBI identifiers).

    keep_reaction_db_ids: when the wanted pathway has no diagram of its own (it is
    drawn inside an ancestor's), only the reactions in this set - and the nodes
    they touch - are kept."""
    graph_by_id = {n.get("dbId"): n for n in (graph.get("nodes") or [])}
    keep = {int(i) for i in keep_reaction_db_ids} if keep_reaction_db_ids is not None else None

    reactions, touched = [], set()
    for edge in diagram.get("edges") or []:
        rdb = edge.get("reactomeId")
        if keep is not None and rdb not in keep:
            continue
        ins = [f"r{i.get('id')}" for i in edge.get("inputs") or []]
        outs = [f"r{i.get('id')}" for i in edge.get("outputs") or []]
        cats = [f"r{i.get('id')}" for i in edge.get("catalysts") or []]
        segs = [[_f(s.get("from", {}).get("x")), _f(s.get("from", {}).get("y")),
                 _f(s.get("to", {}).get("x")), _f(s.get("to", {}).get("y"))] for s in edge.get("segments") or []]
        pos = edge.get("position") or {}
        reactions.append({"id": f"rx{edge.get('id')}", "name": str(edge.get("displayName") or ""),
                          "reversible": str(edge.get("reactionType") or "").lower() == "reversible",
                          "enzymes": cats, "substrates": ins, "products": outs,
                          "pos": [_f(pos.get("x")), _f(pos.get("y"))] if pos else None, "segments": segs})
        touched.update(ins + outs + cats)

    nodes = []
    for d in diagram.get("nodes") or []:
        nid = f"r{d.get('id')}"
        if keep is not None and nid not in touched:
            continue
        rclass = str(d.get("renderableClass") or "")
        prop = d.get("prop") or {}
        pos = d.get("position") or {}
        leaves = _reactome_leaves(d.get("reactomeId"), graph_by_id)
        uniprot, chebi, symbols = [], [], []
        for leaf in leaves:
            ref, ident = str(leaf.get("referenceType") or ""), leaf.get("identifier")
            if ref == "ReferenceGeneProduct" and ident:
                uniprot.append(norm_uniprot(ident))
                symbols.append(_strip_compartment(leaf.get("displayName")).upper())
            elif ref == "ReferenceMolecule" and ident:
                chebi.append(norm_chebi(ident))
        uniprot, chebi = sorted({u for u in uniprot if u}), sorted({c for c in chebi if c})
        if rclass in _REACTOME_LINK_CLASSES:
            kind = "link"
        elif rclass in _REACTOME_COMPOUND_CLASSES or (rclass in _REACTOME_GENE_CLASSES and chebi and not uniprot):
            # A set / complex whose every leaf is a small molecule is a metabolite node.
            kind = "compound"
        elif rclass in _REACTOME_GENE_CLASSES:
            kind = "gene"
        else:
            kind = "other"
        members = {"uniprot": uniprot, "symbols": sorted({s for s in symbols if s})} if kind == "gene" else \
                  {"chebi": chebi} if kind == "compound" else {}
        node = {"id": nid, "kind": kind, "type": rclass,
                "label": _strip_compartment(d.get("displayName")),
                "x": _f(pos.get("x")), "y": _f(pos.get("y")),
                "w": _f(prop.get("width")), "h": _f(prop.get("height")), "members": members}
        if kind == "link":
            node["reactome_id"] = d.get("reactomeId")  # the sub-pathway this box stands for
        nodes.append(node)
    return {"title": str(diagram.get("displayName") or ""), "nodes": nodes, "reactions": reactions,
            "relations": [], "decor_lines": [], "line_graphics": False}


def is_overview(structure: dict) -> bool:
    """A pathway whose diagram holds only links to sub-pathways (Reactome's 'Disease', 'Metabolism', 'Immune System',
    ...): there is no gene or compound to put a hit on, so nothing can be mapped on it."""
    kinds = {n["kind"] for n in structure.get("nodes") or []}
    return "link" in kinds and not (kinds & {"gene", "compound", "group", "other"})


# ---------------------------------------------------------------------------
# Per-omics feature tables. A "feature" is one gene / protein / compound as the
# user measured it:
#   {"omics","label","hit","log2fc","padj","symbol","uniprots": set, "kegg_genes": set,
#    "kegg_via","uniprot_via","kegg_compounds": set, "chebi": set}
# ---------------------------------------------------------------------------

def invert_acc_to_gene(acc_to_gene: dict) -> dict:
    """{kegg_gene: {accession}} from KEGG's own, already de-ambiguated acc_to_gene."""
    out: dict = {}
    for acc, gene in (acc_to_gene or {}).items():
        out.setdefault(gene, set()).add(norm_uniprot(acc))
    return out


def _bridge_gene(feature: dict, gene_map: "dict | None", gene_to_acc: "dict | None") -> dict:
    """Attach KEGG gene IDs (and, when no UniProt is known, UniProt accessions) to
    one gene / protein feature, recording which bridge produced each."""
    gm = gene_map or {}
    acc_to_gene, entrez_to_gene, symbol_to_gene = (gm.get("acc_to_gene") or {}, gm.get("entrez_to_gene") or {},
                                                   gm.get("symbol_to_gene") or {})
    if not feature["kegg_genes"]:
        for acc in sorted(feature["uniprots"]):
            g = acc_to_gene.get(acc)
            if g:
                feature["kegg_genes"].add(g)
                feature["kegg_via"] = "UniProt -> KEGG gene"
        if not feature["kegg_genes"] and feature.get("entrez"):
            g = entrez_to_gene.get(str(feature["entrez"]))
            if g:
                feature["kegg_genes"].add(g)
                feature["kegg_via"] = "Entrez Gene ID -> KEGG gene"
        if not feature["kegg_genes"] and feature.get("symbol"):
            g = symbol_to_gene.get(str(feature["symbol"]).upper())
            if g:
                feature["kegg_genes"].add(g)
                feature["kegg_via"] = "gene symbol -> KEGG gene"
    else:
        feature["kegg_via"] = feature.get("kegg_via") or "KEGG gene ID"
    if not feature["uniprots"] and feature["kegg_genes"] and gene_to_acc:
        for g in feature["kegg_genes"]:
            feature["uniprots"] |= gene_to_acc.get(g, set())
        if feature["uniprots"]:
            feature["uniprot_via"] = "KEGG gene -> UniProt"
    elif feature["uniprots"]:
        feature["uniprot_via"] = "UniProt"
    return feature


def _new_gene_feature(omics: str, label: str, hit: bool, log2fc=None, padj=None, symbol: str = "",
                      uniprots=(), kegg_genes=(), entrez: str = "") -> dict:
    return {"omics": omics, "label": label, "hit": bool(hit), "log2fc": finite_or_none(log2fc),
            "padj": finite_or_none(padj), "symbol": str(symbol or "").upper(),
            "uniprots": {u for u in (norm_uniprot(x) for x in uniprots) if u},
            "kegg_genes": {str(k).strip() for k in kegg_genes if str(k).strip()},
            "entrez": str(entrez or "").strip(), "kegg_via": "", "uniprot_via": "",
            "kegg_compounds": set(), "chebi": set()}


def gene_features_genomics(final_df: pd.DataFrame, gene_map: "dict | None" = None) -> list:
    """Genomics is a plain list: every usable row is a hit. final_df: the resolved
    table with ambiguous / unresolved rows already removed."""
    id_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else "Ensembl Gene ID"
    gene_to_acc = invert_acc_to_gene((gene_map or {}).get("acc_to_gene"))
    out = []
    for _, row in final_df.iterrows():
        label = str(row.get(id_col, "")).strip()
        feat = _new_gene_feature(
            "genomics", label, True, symbol=label,
            uniprots=[row.get("UniProt ID", "")] if str(row.get("UniProt ID", "")).strip() else [],
            kegg_genes=[row.get("KEGG ID", "")] if str(row.get("KEGG ID", "")).strip() else [],
            entrez=row.get("NCBI Gene ID", ""))
        out.append(_bridge_gene(feat, gene_map, gene_to_acc))
    return out


def gene_features_de(omics: str, df: pd.DataFrame, significant: set, gene_map: "dict | None" = None) -> list:
    """Transcriptomics / Proteomics: every row of the table is MEASURED; a hit is a
    member of `significant` (symbols for Transcriptomics, accessions for
    Proteomics - the same sets the joint test uses)."""
    gene_to_acc = invert_acc_to_gene((gene_map or {}).get("acc_to_gene"))
    out = []
    for _, row in df.iterrows():
        symbol = str(row.get("gene_symbol", "")).strip()
        if omics == "proteomics":
            raw = str(row.get("accession", "")).strip()
            acc = norm_uniprot(raw)
            if not acc:
                continue
            hit = raw in significant or acc in significant
            feat = _new_gene_feature(omics, symbol or acc, hit, row.get("log2FC"), row.get("padj"),
                                     symbol=symbol, uniprots=[acc])
        else:
            if not symbol:
                continue
            hit = symbol.upper() in significant
            feat = _new_gene_feature(omics, symbol, hit, row.get("log2FC"), row.get("padj"),
                                     symbol=symbol, entrez=row.get("entrez_id", ""))
        out.append(_bridge_gene(feat, gene_map, gene_to_acc))
    return out


def compound_features(matched_df: pd.DataFrame) -> list:
    """Metabolomics is a plain list of matched compounds; class-level-only rows
    (lipid shorthand such as 'PC 34:1') are not a specific structure and are
    left out - the same rule the pathway tests apply."""
    df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
    if "class_level_only" in df.columns:
        df = df[~df["class_level_only"].astype(bool)]
    out = []
    for _, row in df.iterrows():
        kegg, chebi = norm_kegg_compound(row.get("kegg_id")), norm_chebi(row.get("chebi_id"))
        skeleton = inchikey_skeleton(row.get("inchikey", "")) if "inchikey" in df.columns else ""
        if not kegg and not chebi and not skeleton:
            continue  # no pathway-database identifier: it can never sit on a pathway map
        out.append({"omics": "metabolomics", "label": str(row.get("Name", "") or row.get("compound_key")).strip(),
                    "hit": True, "log2fc": None, "padj": None, "symbol": "", "uniprots": set(), "kegg_genes": set(),
                    "entrez": "", "kegg_via": "", "uniprot_via": "",
                    "kegg_compounds": {kegg} if kegg else set(), "chebi": {chebi} if chebi else set(),
                    "skeleton": skeleton or ""})
    return out


# ---------------------------------------------------------------------------
# Mapping features onto a structure
# ---------------------------------------------------------------------------

def _index(features: list, getter) -> dict:
    idx: dict = {}
    for feat in features:
        for key in getter(feat):
            idx.setdefault(key, []).append(feat)
    return idx


def _direction(members: list):
    signs = {1 if m["log2fc"] > 0 else -1 for m in members if m["hit"] and m["log2fc"] not in (None, 0)}
    if not signs:
        return None
    return "up" if signs == {1} else "down" if signs == {-1} else "mixed"


def structure_chebi_ids(structure: dict) -> set:
    """Every ChEBI number on the structure's metabolite nodes (for the skeleton lookup)."""
    return {c for n in structure["nodes"] if n["kind"] == "compound" for c in n["members"].get("chebi", [])}


def map_omics_to_structure(structure: dict, features_by_omics: dict, database: str,
                           chebi_skeletons: "dict | None" = None) -> list:
    """Return the structure's nodes with a per-omics status attached:
    node["status"][omics] = {"state": "hit"|"measured"|"absent", "direction": ..., "members": [...]}.
    Only the omics in features_by_omics are mapped; the others stay absent.

    chebi_skeletons: {ChEBI number: InChIKey skeleton} for the structure's chemicals (Reactome only);
    without it a chemical can only be matched on its exact ChEBI ID."""
    chebi_skeletons = chebi_skeletons or {}
    indexes: dict = {}
    for omics, feats in features_by_omics.items():
        if omics == "metabolomics":
            indexes[omics] = {"cpd": _index(feats, lambda f: f["kegg_compounds"]) if database == "KEGG" else {},
                              "chebi": _index(feats, lambda f: f["chebi"]) if database == "Reactome" else {},
                              "skeleton": (_index(feats, lambda f: [f["skeleton"]] if f.get("skeleton") else [])
                                           if database == "Reactome" else {})}
        elif database == "KEGG":
            indexes[omics] = {"gene": _index(feats, lambda f: f["kegg_genes"])}
        else:
            indexes[omics] = {"uniprot": _index(feats, lambda f: f["uniprots"]),
                              # symbol is a last resort, only for a feature with NO UniProt of its own
                              "symbol": _index([f for f in feats if not f["uniprots"] and f["symbol"]],
                                               lambda f: [f["symbol"]])}
    out = []
    for node in structure["nodes"]:
        status = {o: {"state": "absent", "direction": None, "members": []} for o in OMICS_ORDER}
        for omics, idx in indexes.items():
            matched: dict = {}
            if node["kind"] in ("gene", "group") and omics in GENE_OMICS:
                m = node["members"]
                if database == "KEGG":
                    for g in m.get("kegg_genes", []):
                        for f in idx["gene"].get(g, []):
                            matched.setdefault(id(f), (f, f["kegg_via"] or "KEGG gene ID"))
                else:
                    for u in m.get("uniprot", []):
                        for f in idx["uniprot"].get(u, []):
                            matched.setdefault(id(f), (f, f["uniprot_via"] or "UniProt"))
                    for s in m.get("symbols", []):
                        for f in idx["symbol"].get(s, []):
                            matched.setdefault(id(f), (f, "gene symbol (no UniProt known)"))
            elif node["kind"] == "compound" and omics == "metabolomics":
                m = node["members"]
                if database == "KEGG":
                    for c in m.get("kegg_compounds", []):
                        for f in idx["cpd"].get(c, []):
                            matched.setdefault(id(f), (f, "KEGG compound ID"))
                else:
                    for c in m.get("chebi", []):
                        for f in idx["chebi"].get(c, []):
                            matched.setdefault(id(f), (f, "ChEBI ID"))
                    # A different protonation state of the same compound: same InChIKey skeleton. Exact ChEBI wins.
                    for c in m.get("chebi", []):
                        skel = chebi_skeletons.get(c)
                        for f in idx["skeleton"].get(skel, []) if skel else []:
                            matched.setdefault(id(f), (f, "InChIKey skeleton (another protonation state)"))
            if not matched:
                continue
            members = [{"label": f["label"], "hit": f["hit"], "log2fc": f["log2fc"], "padj": f["padj"], "via": via}
                       for f, via in matched.values()]
            members.sort(key=lambda r: (not r["hit"], str(r["label"])))
            hit = any(r["hit"] for r in members)
            status[omics] = {"state": "hit" if hit else "measured",
                             "direction": _direction(members) if omics in ("transcriptomics", "proteomics") else None,
                             "members": members}
        out.append({**node, "status": status})
    return out


def _omics_phrase(status: dict, omics: str) -> str:
    d = status[omics]["direction"]
    return OMICS_LABELS[omics] + (f" ({d})" if d else "")


def reaction_links(structure: dict, mapped_nodes: list, selected: "list | tuple") -> list:
    """Reactions with cross-omics evidence (design rule 5)."""
    by_id = {n["id"]: n for n in mapped_nodes}
    has_met = "metabolomics" in selected
    rows = []
    for rx in structure["reactions"]:
        enz_nodes = [by_id[i] for i in rx["enzymes"] if i in by_id and by_id[i]["kind"] == "gene"]
        enz_hit_omics = [o for o in GENE_OMICS if o in selected and any(n["status"][o]["state"] == "hit" for n in enz_nodes)]
        compounds = []
        if has_met:
            for role, ids in (("substrate", rx["substrates"]), ("product", rx["products"])):
                for i in ids:
                    n = by_id.get(i)
                    if n and n["kind"] == "compound" and n["status"]["metabolomics"]["state"] == "hit":
                        for m in n["status"]["metabolomics"]["members"]:
                            if m["hit"]:
                                compounds.append(f"{m['label']} ({role})")
        if not ((enz_hit_omics and compounds) or len(enz_hit_omics) >= 2):
            continue
        parts, dirs = [], {}
        for o in enz_hit_omics:
            sset = {n["status"][o]["direction"] for n in enz_nodes if n["status"][o]["state"] == "hit"} - {None}
            dirs[o] = sset.pop() if len(sset) == 1 else ("mixed" if sset else None)
            parts.append(OMICS_LABELS[o] + (f" ({dirs[o]})" if dirs[o] else ""))
        t, p = dirs.get("transcriptomics"), dirs.get("proteomics")
        if t in ("up", "down") and p in ("up", "down"):
            concord = "concordant" if t == p else "discordant"
        else:
            concord = "n/a"
        enz_labels = []
        for n in enz_nodes:
            for o in enz_hit_omics:
                for m in n["status"][o]["members"]:
                    if m["hit"] and m["label"] not in enz_labels:
                        enz_labels.append(m["label"])
        rows.append({
            "Reaction": rx["name"], "Enzyme (hit)": ", ".join(enz_labels),
            "Enzyme omics": "; ".join(parts), "Metabolites hit": "; ".join(dict.fromkeys(compounds)),
            "Evidence": ("enzyme + metabolite" if (enz_hit_omics and compounds) else "enzyme in several omics")
                        + (" + enzyme in several omics" if (enz_hit_omics and compounds and len(enz_hit_omics) >= 2) else ""),
            "Transcriptomics vs Proteomics direction": concord,
            "Reversible": "yes" if rx["reversible"] else "no",
        })
    return rows


MEMBER_COLUMNS = ["Pathway node", "Box ID", "Node type", "Member", "Omics", "Status", "Direction", "log2FC", "padj", "Matched by"]
LINK_COLUMNS = ["Reaction", "Enzyme (hit)", "Enzyme omics", "Metabolites hit", "Evidence",
                "Transcriptomics vs Proteomics direction", "Reversible"]


def member_rows(mapped_nodes: list, selected: "list | tuple") -> list:
    rows = []
    for n in mapped_nodes:
        if n["kind"] == "group":
            continue  # a drawn complex repeats its components' members, which have rows of their own
        for o in OMICS_ORDER:
            if o not in selected:
                continue
            st = n["status"][o]
            for m in st["members"]:
                rows.append({"Pathway node": n["label"], "Box ID": n["id"], "Node type": "compound" if n["kind"] == "compound" else "gene / protein",
                             "Member": m["label"], "Omics": OMICS_LABELS[o],
                             "Status": "hit" if m["hit"] else "measured, not significant",
                             "Direction": (("up" if m["log2fc"] > 0 else "down") if (m["hit"] and m["log2fc"]) else ""),
                             "log2FC": m["log2fc"], "padj": m["padj"], "Matched by": m["via"]})
    return rows


def summarize_mapping(mapped_nodes: list, links: list, selected: "list | tuple", interactions: "list | None" = None) -> dict:
    """Counts for the overview. Boxes are counted per BOX (a gene KEGG draws five times is five boxes), members per
    distinct member (that same gene is ONE member hit, not five). A group (drawn complex) is not a box of its own."""
    boxes = [n for n in mapped_nodes if n["kind"] != "group"]
    per_omics = {}
    for o in selected:
        hit_labels = {m["label"] for n in boxes for m in n["status"][o]["members"] if m["hit"]}
        per_omics[o] = {"label": OMICS_LABELS[o],
                        "nodes_hit": sum(1 for n in boxes if n["status"][o]["state"] == "hit"),
                        "nodes_measured_only": sum(1 for n in boxes if n["status"][o]["state"] == "measured"),
                        "members_hit": len(hit_labels)}
    return {"omics": list(selected), "per_omics": per_omics, "n_nodes": len(boxes),
            "n_nodes_with_hit": sum(1 for n in boxes if any(n["status"][o]["state"] == "hit" for o in selected)),
            "n_links": len(links), "n_interactions": len(interactions or [])}


INTERACTION_COLUMNS = ["From", "To", "Relation", "From hit by", "To hit by"]


def _hit_summary(node: dict, selected) -> "tuple[list, str]":
    """(omics that hit this box, 'Genomics; Transcriptomics (up)') - '' omics when none did."""
    hit = [o for o in selected if node["status"][o]["state"] == "hit"]
    return hit, "; ".join(_omics_phrase(node["status"], o) for o in hit)


def _box_name(node: dict, selected) -> str:
    """The box's label, or - for a group, whose label is built from its components - the hit members' own names."""
    names: list = []
    for o in selected:
        for m in node["status"][o]["members"]:
            if m["hit"] and m["label"] not in names:
                names.append(m["label"])
    if node["kind"] == "group" and names:
        return node["label"] + " (" + ", ".join(names[:4]) + (", ..." if len(names) > 4 else "") + ")"
    return node["label"]


def interaction_links(structure: dict, mapped_nodes: list, selected: "list | tuple") -> list:
    """Relations of a signalling map (activation, inhibition, binding, phosphorylation, ...) whose BOTH partners were
    hit by some omics. KEGG's relation is the database's own statement; that both ends were hit is what is shown
    here, descriptive and not a test."""
    by_id = {n["id"]: n for n in mapped_nodes}
    rows = []
    for rel in structure.get("relations") or []:
        a, b = by_id.get(rel["a"]), by_id.get(rel["b"])
        if a is None or b is None:
            continue
        hit_a, text_a = _hit_summary(a, selected)
        hit_b, text_b = _hit_summary(b, selected)
        if not hit_a or not hit_b:
            continue
        rows.append({"From": _box_name(a, selected), "To": _box_name(b, selected),
                     "Relation": ", ".join(rel["types"]), "From hit by": text_a, "To hit by": text_b})
    return rows


# ---------------------------------------------------------------------------
# Pathway picker (coverage, no statistics)
# ---------------------------------------------------------------------------

PICKER_COLUMNS_FIXED = ["Term ID", "Term Name", "Omics with a hit", "Total hits", "Joint combined FDR"]


def build_pathway_picker(inputs: dict, selected: "list | tuple", joint_fdr: "dict | None" = None) -> list:
    """One row per pathway that at least one omics hit, ranked by how many omics
    hit it and then by the number of hits. `inputs` is the per-omics dict the
    joint test already builds (term_map / list_set / universe_set)."""
    joint_fdr = joint_fdr or {}
    ids: list = []
    names: dict = {}
    for o in selected:
        for tid, info in (inputs.get(o, {}).get("term_map") or {}).items():
            if tid not in names:
                ids.append(tid)
                names[tid] = str(info.get("name") or tid)
    rows = []
    for tid in ids:
        row = {"Term ID": tid, "Term Name": names[tid]}
        n_omics = total = 0
        for o in selected:
            inp = inputs.get(o) or {}
            info = (inp.get("term_map") or {}).get(tid)
            genes = set(info["genes"]) if info else set()
            k = len(genes & set(inp.get("list_set") or ()))
            K = len(genes & set(inp.get("universe_set") or ()))
            row[f"{OMICS_LABELS[o]}: hit / annotated"] = f"{k} / {K}"
            if k:
                n_omics += 1
                total += k
        if not total:
            continue
        row["Omics with a hit"] = n_omics
        row["Total hits"] = total
        row["Joint combined FDR"] = finite_or_none(joint_fdr.get(tid))
        rows.append(row)
    rows.sort(key=lambda r: (-r["Omics with a hit"], -r["Total hits"], r["Term Name"]))
    return rows


def picker_columns(selected: "list | tuple") -> list:
    return (["Term ID", "Term Name", "Omics with a hit", "Total hits"]
            + [f"{OMICS_LABELS[o]}: hit / annotated" for o in selected] + ["Joint combined FDR"])


def relabel_kegg_compounds(structure: dict, names: dict) -> None:
    """KGML labels a compound by its C-number; show KEGG's own primary name instead
    (names: {C-number: name}). A compound with no known name keeps the C-number."""
    for node in structure["nodes"]:
        if node["kind"] == "compound":
            for c in node["members"].get("kegg_compounds", []):
                if names.get(c):
                    node["label"] = str(names[c]).split(";")[0].strip()
                    break


def overview_rows(subpathways: list, picker_rows: list, selected: "list | tuple") -> list:
    """The sub-pathways of an overview pathway with the picker's hit counts. The picker lists only pathways some omics
    hit, so a sub-pathway it does not list has no hit (shown as 0 / -: the annotated count is not known here)."""
    by_id = {r["Term ID"]: r for r in picker_rows}
    rows = []
    for sp in subpathways:
        r = by_id.get(sp["stId"])
        row = {"Term ID": sp["stId"], "Term Name": sp["name"],
               "Omics with a hit": r["Omics with a hit"] if r else 0, "Total hits": r["Total hits"] if r else 0}
        for o in selected:
            col = f"{OMICS_LABELS[o]}: hit / annotated"
            row[col] = r[col] if r else "0 / -"
        row["Joint combined FDR"] = r.get("Joint combined FDR") if r else None
        rows.append(row)
    rows.sort(key=lambda r: (-r["Omics with a hit"], -r["Total hits"], r["Term Name"]))
    return rows


# ---------------------------------------------------------------------------
# Pathway-class colours (KEGG's global overview maps)
#
# KEGG paints every line of a global map in the colour of the pathway CLASS it belongs to (Carbohydrate metabolism,
# Lipid metabolism, ...), and stores that colour in the line's own graphics. The colour carries no class NAME, so the
# name comes from KEGG's own classification (br08901, map number -> class) by a vote: a colour is named by the class most of
# its genes' pathways belong to. The global and overview maps themselves (01100 ...) are not a class a gene belongs to
# and take no part in the vote.
# ---------------------------------------------------------------------------

OVERVIEW_CLASS = "Global and overview maps"
OTHER_CLASS = "Other"


def parse_brite_classes(text: str) -> dict:
    """KEGG's pathway classification (br08901) -> {map number: class}, for the METABOLISM branch only. Its levels: 'A'
    (Metabolism, Genetic Information Processing, Human Diseases, ...), 'B' (Carbohydrate metabolism, ...), 'C' (one pathway
    each, '00010  Glycolysis / ...'). The class a pathway is coloured by is the B level. Pathways of the other A branches are
    left out on purpose: a metabolic map's genes also sit on disease pathways ('Neurodegenerative disease' shares many
    metabolic enzymes) and those must not outvote the metabolic class the line is coloured by."""
    out: dict = {}
    a = b = ""
    for line in str(text or "").splitlines():
        if line.startswith("A"):
            a = line[1:].strip()
            b = ""
        elif line.startswith("B"):
            b = line[1:].strip()
        elif line.startswith("C") and b and a == "Metabolism":
            parts = line[1:].split()
            if parts and parts[0].isdigit():
                out[parts[0]] = b
    return out


def assign_color_classes(nodes: list, gene_to_maps: dict, map_to_class: dict) -> dict:
    """{line colour: class name}: the class most of the colour's genes' pathways belong to (ties: alphabetical, so the
    answer never depends on dict order); a colour none of whose genes is on any classified pathway is 'Other'."""
    votes: dict = {}
    for n in nodes:
        color = n.get("color")
        if not color or n.get("kind") != "gene":
            continue
        tally = votes.setdefault(color, {})
        for gene in n["members"].get("kegg_genes", []):
            for num in gene_to_maps.get(gene, ()):
                cls = map_to_class.get(num)
                if cls and cls != OVERVIEW_CLASS:
                    tally[cls] = tally.get(cls, 0) + 1
    out = {}
    for color, tally in votes.items():
        out[color] = sorted(tally.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] if tally else OTHER_CLASS
    return out


def class_legend(nodes: list, color_class: dict) -> list:
    """[{name, color, n_lines, n_hit}] sorted by name. A class drawn in several KEGG colours takes its most frequent one as
    its default colour (the user can recolour the class as a whole)."""
    by_cls: dict = {}
    for n in nodes:
        if n.get("kind") != "gene" or not n.get("color"):
            continue
        cls = color_class.get(n["color"], OTHER_CLASS)
        e = by_cls.setdefault(cls, {"name": cls, "colors": {}, "n_lines": 0, "n_hit": 0})
        e["colors"][n["color"]] = e["colors"].get(n["color"], 0) + 1
        e["n_lines"] += 1
        if any(v["state"] == "hit" for v in n["status"].values()):
            e["n_hit"] += 1
    out = []
    for cls, e in by_cls.items():
        color = sorted(e["colors"].items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        out.append({"name": cls, "color": color, "n_lines": e["n_lines"], "n_hit": e["n_hit"]})
    out.sort(key=lambda r: r["name"])
    return out
