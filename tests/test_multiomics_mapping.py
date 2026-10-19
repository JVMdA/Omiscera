# -*- coding: utf-8 -*-
"""Unit tests for Multiomics_mapping (Cross-omics Pathway Mapping engine) and the structure
fetchers in Multiomics_module_retriever, offline - small embedded KGML / Reactome fixtures and
mocked HTTP.

What is under test (each rule is mutation-tested; a mutant breaks the named check):
  - the KGML parser keeps GENE entries and the reaction -> enzyme link (the existing Metabolomics
    parser drops both, which is why this is a separate parser), reversible flags, glycans excluded;
  - the Reactome parser opens a complex / set down to its UniProt / ChEBI leaves, treats a set of
    only small molecules as a metabolite, guards cycles, strips an isoform suffix, and restricts to
    a sub-pathway's own reactions;
  - three states per omics per node - hit / measured / absent - never two;
  - a node is hit if ANY member is (multi-gene KEGG box, Reactome complex);
  - every match records the identifier bridge it used;
  - a cofactor inside a Reactome complex never makes the complex a metabolite hit;
  - a reaction link needs enzyme-hit + metabolite-hit, or enzyme-hit in two omics; nothing else;
  - the fetchers never cache a transport error, do cache a definite "nodata", and draw a pathway
    that has no diagram of its own from its NEAREST ancestor, restricted to its own reactions.
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Multiomics_mapping as mm  # noqa: E402
import Multiomics_module_retriever as mr  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def node(structure, nid):
    return next(n for n in structure["nodes"] if n["id"] == nid)


def mapped_node(mapped, nid):
    return next(n for n in mapped if n["id"] == nid)


# ---------------------------------------------------------------------------
KGML = """<?xml version="1.0"?>
<pathway name="path:hsa99999" org="hsa" number="99999" title="Test cycle">
  <entry id="1" name="hsa:100 hsa:101" type="gene" reaction="rn:R00001">
    <graphics name="AAA, BBB" x="100" y="100" width="46" height="17"/></entry>
  <entry id="2" name="hsa:200" type="gene" reaction="rn:R00002">
    <graphics name="CCC, X" x="200" y="100" width="46" height="17"/></entry>
  <entry id="3" name="ko:K00001" type="ortholog" reaction="rn:R00004">
    <graphics name="K00001" x="300" y="100" width="46" height="17"/></entry>
  <entry id="4" name="cpd:C00001" type="compound"><graphics name="C00001" x="100" y="200" width="8" height="8"/></entry>
  <entry id="5" name="cpd:C00002" type="compound"><graphics name="C00002" x="150" y="200" width="8" height="8"/></entry>
  <entry id="6" name="cpd:C00003" type="compound"><graphics name="C00003" x="200" y="200" width="8" height="8"/></entry>
  <entry id="7" name="gl:G00001" type="compound"><graphics name="G00001" x="250" y="200" width="8" height="8"/></entry>
  <entry id="8" name="path:hsa00010" type="map"><graphics name="Glycolysis" x="300" y="200" width="40" height="20"/></entry>
  <entry id="9" name="hsa:300" type="gene" reaction="rn:R00005">
    <graphics name="DDD" x="400" y="100" width="46" height="17"/></entry>
  <entry id="11" name="path:hsa99999" type="map"><graphics name="TITLE:Test cycle" x="50" y="20" width="40" height="20"/></entry>
  <entry id="10" name="hsa:400" type="gene" reaction="rn:R00006">
    <graphics name="EEE" x="500" y="100" width="46" height="17"/></entry>
  <reaction id="1" name="rn:R00001" type="reversible"><substrate id="4" name="cpd:C00001"/><product id="5" name="cpd:C00002"/></reaction>
  <reaction id="2" name="rn:R00002" type="irreversible"><substrate id="5" name="cpd:C00002"/><product id="6" name="cpd:C00003"/></reaction>
  <reaction id="3" name="rn:R00004" type="irreversible"><substrate id="6" name="cpd:C00003"/><product id="4" name="cpd:C00001"/></reaction>
  <reaction id="9" name="rn:R00005" type="irreversible"><substrate id="5" name="cpd:C00002"/><product id="6" name="cpd:C00003"/></reaction>
  <reaction id="10" name="rn:R00006" type="irreversible"><substrate id="4" name="cpd:C00001"/><product id="5" name="cpd:C00002"/></reaction>
</pathway>"""

print("=== KGML parser keeps genes and the reaction -> enzyme link ===")
S = mm.parse_kgml_full(KGML)
kinds = {}
for n in S["nodes"]:
    kinds[n["kind"]] = kinds.get(n["kind"], 0) + 1
check(kinds == {"gene": 4, "other": 1, "compound": 4, "link": 1}, f"node kinds: 4 genes, 1 ortholog (other), 4 compound entries, 1 map link ({kinds})")
check(node(S, "k1")["members"]["kegg_genes"] == ["hsa:100", "hsa:101"],
      "SENSITIVITY: a multi-gene KGML box keeps EVERY gene (the Metabolomics parser drops gene entries altogether)")
check(node(S, "k4")["members"]["kegg_compounds"] == ["C00001"] and node(S, "k7")["members"]["kegg_compounds"] == [],
      "a compound entry yields its C-number; a glycan (gl:G...) yields none, so it can never match a compound")
check(all(not n["label"].startswith("TITLE:") for n in S["nodes"]) and len(S["nodes"]) == 10,
      "the map's own TITLE caption is not a box (it is not a link to another pathway)")
check(node(S, "k3")["kind"] == "other" and "kegg_genes" not in node(S, "k3")["members"],
      "an ortholog (ko:) entry is drawn but carries no gene members")
check(node(S, "k1")["label"] == "AAA" and (node(S, "k1")["x"], node(S, "k1")["y"], node(S, "k1")["w"]) == (100.0, 100.0, 46.0),
      "label and KEGG's own x / y / width are kept")
rx = {r["name"]: r for r in S["reactions"]}
check(rx["R00001"]["enzymes"] == ["k1"] and rx["R00001"]["substrates"] == ["k4"] and rx["R00001"]["products"] == ["k5"],
      "SENSITIVITY: a reaction links its ENZYME node to its substrate and product nodes (entry ids, not compound ids)")
check(rx["R00001"]["reversible"] is True and rx["R00002"]["reversible"] is False, "reversible flag read from the reaction type")
check(len(S["reactions"]) == 5, "every reaction element is kept")

print("=== feature tables and identifier bridges ===")
GENE_MAP = {"acc_to_gene": {"P2": "hsa:200", "P9": "hsa:300"},
            "entrez_to_gene": {"100": "hsa:100", "101": "hsa:101", "200": "hsa:200", "300": "hsa:300"},
            "symbol_to_gene": {"EEE": "hsa:400"}}
genomics_df = pd.DataFrame({"Gene Symbol": ["AAA"], "UniProt ID": ["P1"], "NCBI Gene ID": [""], "KEGG ID": ["hsa:100"]})
t_df = pd.DataFrame({"gene_symbol": ["BBB", "CCC", "DDD", "EEE", "ZZZ"], "entrez_id": ["101", "200", "300", "", ""],
                     "log2FC": [2.5, 3.0, 0.1, 1.5, -1.0], "padj": [0.001, 0.002, 0.9, 0.01, 0.5]})
p_df = pd.DataFrame({"accession": ["P2", "P9"], "gene_symbol": ["CCC", "DDD"], "log2FC": [-2.0, 0.2], "padj": [0.003, 0.8]})
matched_df = pd.DataFrame({"compound_key": ["k1", "k2", "k3", "k4"], "Name": ["Cmp1", "Cmp3", "ShorthandLipid", "NoIds"],
                           "kegg_id": ["C00001", "cpd:C00003", "C00002", ""], "chebi_id": ["CHEBI:30031", "", "", ""],
                           "class_level_only": [False, False, True, False]})
F = {"genomics": mm.gene_features_genomics(genomics_df, GENE_MAP),
     "transcriptomics": mm.gene_features_de("transcriptomics", t_df, {"BBB", "CCC", "EEE"}, GENE_MAP),
     "proteomics": mm.gene_features_de("proteomics", p_df, {"P2"}, GENE_MAP),
     "metabolomics": mm.compound_features(matched_df)}
tf = {f["label"]: f for f in F["transcriptomics"]}
check(tf["BBB"]["kegg_genes"] == {"hsa:101"} and tf["BBB"]["kegg_via"] == "Entrez Gene ID -> KEGG gene",
      "Transcriptomics: Entrez -> KEGG gene, and the bridge is recorded")
check(tf["EEE"]["kegg_genes"] == {"hsa:400"} and tf["EEE"]["kegg_via"] == "gene symbol -> KEGG gene",
      "a gene with no Entrez ID falls back to the gene symbol, labelled as the weakest bridge")
check(tf["ZZZ"]["kegg_genes"] == set(), "a gene KEGG cannot place is NOT guessed")
pf = {f["label"]: f for f in F["proteomics"]}
check(pf["CCC"]["kegg_genes"] == {"hsa:200"} and pf["CCC"]["kegg_via"] == "UniProt -> KEGG gene" and pf["CCC"]["hit"] is True
      and pf["DDD"]["hit"] is False, "Proteomics: UniProt -> KEGG gene; hit only for accessions in the significant set")
check(F["genomics"][0]["kegg_via"] == "KEGG gene ID" and F["genomics"][0]["hit"] is True,
      "Genomics: the assigned KEGG ID is used directly and every listed gene is a hit")
check(len(F["metabolomics"]) == 2 and {f["label"] for f in F["metabolomics"]} == {"Cmp1", "Cmp3"},
      "SENSITIVITY: class-level-only (lipid shorthand) and unmatched compounds are not features")
check(tf["BBB"]["hit"] and not tf["DDD"]["hit"] and tf["DDD"]["log2fc"] == 0.1, "Transcriptomics hit = member of the significant set; others stay measured")

print("=== KEGG mapping: three states, any-member rule, bridges ===")
M = mm.map_omics_to_structure(S, F, "KEGG")
n1, n2, n9 = mapped_node(M, "k1"), mapped_node(M, "k2"), mapped_node(M, "k9")
check(n1["status"]["genomics"]["state"] == "hit" and n1["status"]["transcriptomics"]["state"] == "hit"
      and n1["status"]["proteomics"]["state"] == "absent",
      "SENSITIVITY: a two-gene box is hit by Genomics (hsa:100 only) AND Transcriptomics (hsa:101 only) - ANY member counts")
check(n1["status"]["transcriptomics"]["direction"] == "up" and n1["status"]["genomics"]["direction"] is None,
      "direction is reported for a differential-expression omics only (a list has none)")
check(n9["status"]["transcriptomics"]["state"] == "measured" and n9["status"]["proteomics"]["state"] == "measured"
      and n9["status"]["genomics"]["state"] == "absent",
      "SENSITIVITY: measured-but-not-significant is its own state, distinct from absent")
check(n2["status"]["transcriptomics"]["direction"] == "up" and n2["status"]["proteomics"]["direction"] == "down",
      "the same node carries a different direction per omics")
check([m["via"] for m in n1["status"]["transcriptomics"]["members"]] == ["Entrez Gene ID -> KEGG gene"],
      "every matched member records its bridge")
c1, c2, c3 = mapped_node(M, "k4"), mapped_node(M, "k5"), mapped_node(M, "k6")
check(c1["status"]["metabolomics"]["state"] == "hit" and c3["status"]["metabolomics"]["state"] == "hit"
      and c2["status"]["metabolomics"]["state"] == "absent",
      "compounds: matched by KEGG C-number ('cpd:C00003' normalised); the class-level one (C00002) is absent")
check(mapped_node(M, "k3")["status"]["genomics"]["state"] == "absent" and mapped_node(M, "k8")["status"]["transcriptomics"]["state"] == "absent",
      "an ortholog entry and a map link can never be hit")
check(mapped_node(M, "k7")["status"]["metabolomics"]["state"] == "absent", "a glycan can never be a Metabolomics hit")
M_only_g = mm.map_omics_to_structure(S, {"genomics": F["genomics"]}, "KEGG")
check(mapped_node(M_only_g, "k1")["status"]["transcriptomics"]["state"] == "absent",
      "an omics that is not part of the run stays absent everywhere")

# a box with TWO genes, one a hit and one only measured: the box is a hit (any member), and both members are listed
two = pd.DataFrame({"gene_symbol": ["BBB", "AAB"], "entrez_id": ["101", "100"], "log2FC": [2.5, 0.1], "padj": [0.001, 0.9]})
M_two = mm.map_omics_to_structure(S, {"transcriptomics": mm.gene_features_de("transcriptomics", two, {"BBB"}, GENE_MAP)}, "KEGG")
st_two = mapped_node(M_two, "k1")["status"]["transcriptomics"]
check(st_two["state"] == "hit" and [m["hit"] for m in st_two["members"]] == [True, False],
      "SENSITIVITY: a two-gene box with one hit gene and one merely-measured gene is HIT (any member), and lists both, hit first")
only_measured = mm.map_omics_to_structure(S, {"transcriptomics": mm.gene_features_de("transcriptomics", two, set(), GENE_MAP)}, "KEGG")
check(mapped_node(only_measured, "k1")["status"]["transcriptomics"]["state"] == "measured",
      "...and with no hit gene the same box is only 'measured'")

print("=== cross-omics reaction links ===")
links = mm.reaction_links(S, M, ["genomics", "transcriptomics", "proteomics", "metabolomics"])
by_rx = {r["Reaction"]: r for r in links}
check(set(by_rx) == {"R00001", "R00002", "R00006"},
      f"links: R00001 (enzyme in 2 omics + metabolite), R00002 (enzyme + metabolite) and R00006 (symbol-bridged enzyme + metabolite) - got {sorted(by_rx)}")
check("enzyme in several omics" in by_rx["R00001"]["Evidence"] and "enzyme + metabolite" in by_rx["R00001"]["Evidence"],
      "R00001 carries both kinds of evidence")
check(by_rx["R00002"]["Transcriptomics vs Proteomics direction"] == "discordant",
      "SENSITIVITY: Transcriptomics up and Proteomics down on the same enzyme is reported 'discordant'")
check(by_rx["R00001"]["Transcriptomics vs Proteomics direction"] == "n/a", "no concordance claim without both DE omics")
check("R00005" not in by_rx, "SENSITIVITY: an enzyme only MEASURED (never a hit) with no metabolite is not a link")
check("R00004" not in by_rx, "SENSITIVITY: an ortholog-only reaction is not a link")
only_enzyme = mm.reaction_links(S, M, ["genomics", "metabolomics"])
check({r["Reaction"] for r in only_enzyme} == {"R00001"},
      "with only Genomics + Metabolomics selected, only a reaction whose enzyme AND metabolite are hit is a link")
no_met = mm.reaction_links(S, M, ["genomics", "transcriptomics", "proteomics"])
check({r["Reaction"] for r in no_met} == {"R00001", "R00002"},
      "without Metabolomics, a link needs the enzyme hit in two or more omics")
check(by_rx["R00001"]["Reversible"] == "yes" and by_rx["R00002"]["Reversible"] == "no", "reversibility is carried into the table")
# compound hit but enzyme not hit anywhere: no link
M_no_enz = mm.map_omics_to_structure(S, {"metabolomics": F["metabolomics"]}, "KEGG")
check(mm.reaction_links(S, M_no_enz, ["genomics", "metabolomics"]) == [],
      "SENSITIVITY: a hit metabolite whose enzyme nobody hit is not a cross-omics link")

print("=== member table and summary ===")
rows = mm.member_rows(M, ["genomics", "transcriptomics", "proteomics", "metabolomics"])
status_set = {(r["Omics"], r["Member"], r["Status"]) for r in rows}
check(("Transcriptomics", "DDD", "measured, not significant") in status_set and ("Transcriptomics", "BBB", "hit") in status_set,
      "member rows keep measured-not-significant members next to hits")
check(all(r["Direction"] in ("up", "down", "") for r in rows) and
      next(r for r in rows if r["Member"] == "BBB")["Direction"] == "up", "direction column only for hits with a log2FC")
sm = mm.summarize_mapping(M, links, ["genomics", "transcriptomics", "proteomics", "metabolomics"])
check(sm["per_omics"]["transcriptomics"]["nodes_hit"] == 3 and sm["n_links"] == 3 and sm["n_nodes"] == 10,
      f"summary counts nodes hit per omics and links ({sm['per_omics']['transcriptomics']}, links {sm['n_links']})")

print("=== KEGG compound relabel ===")
S2 = mm.parse_kgml_full(KGML)
mm.relabel_kegg_compounds(S2, {"C00001": "H2O; Water", "C00002": ""})
check(node(S2, "k4")["label"] == "H2O" and node(S2, "k5")["label"] == "C00002",
      "a compound shows KEGG's primary name; one without a name keeps its C-number")

# ---------------------------------------------------------------------------
print("=== Reactome parser ===")
GRAPH = {"nodes": [
    {"dbId": 10, "stId": "R-HSA-10", "displayName": "AAA [cytosol]", "schemaClass": "EntityWithAccessionedSequence",
     "identifier": "P11111", "referenceType": "ReferenceGeneProduct"},
    {"dbId": 20, "displayName": "AAA:BBB [cytosol]", "schemaClass": "Complex", "children": [21, 22, 23]},
    {"dbId": 21, "displayName": "BBB [cytosol]", "identifier": "P33333-2", "referenceType": "ReferenceGeneProduct"},
    {"dbId": 22, "displayName": "FES [cytosol]", "identifier": "33722", "referenceType": "ReferenceMolecule"},
    {"dbId": 23, "displayName": "loop", "children": [20]},
    {"dbId": 30, "displayName": "Cit [mitochondrial matrix]", "identifier": "30031", "referenceType": "ReferenceMolecule"},
    {"dbId": 40, "displayName": "Sugars [cytosol]", "children": [41, 42]},
    {"dbId": 41, "displayName": "Glc", "identifier": "4167", "referenceType": "ReferenceMolecule"},
    {"dbId": 42, "displayName": "Fru", "identifier": "28757", "referenceType": "ReferenceMolecule"},
    {"dbId": 50, "displayName": "CCC [cytosol]", "identifier": "P44444", "referenceType": "ReferenceGeneProduct"},
]}
DIAGRAM = {"displayName": "Test pathway", "nodes": [
    {"id": 1, "reactomeId": 10, "displayName": "AAA [cytosol]", "renderableClass": "Protein", "position": {"x": 100, "y": 100}, "prop": {"width": 40, "height": 20}},
    {"id": 2, "reactomeId": 20, "displayName": "AAA:BBB [cytosol]", "renderableClass": "Complex", "position": {"x": 200, "y": 100}, "prop": {"width": 60, "height": 20}},
    {"id": 3, "reactomeId": 30, "displayName": "Cit [mitochondrial matrix]", "renderableClass": "Chemical", "position": {"x": 100, "y": 200}, "prop": {"width": 30, "height": 20}},
    {"id": 4, "reactomeId": 40, "displayName": "Sugars [cytosol]", "renderableClass": "EntitySet", "position": {"x": 200, "y": 200}, "prop": {"width": 30, "height": 20}},
    {"id": 5, "reactomeId": 999, "displayName": "Other pathway", "renderableClass": "ProcessNode", "position": {"x": 300, "y": 200}, "prop": {"width": 30, "height": 20}},
    {"id": 6, "reactomeId": 50, "displayName": "CCC [cytosol]", "renderableClass": "Protein", "position": {"x": 400, "y": 100}, "prop": {"width": 40, "height": 20}},
], "edges": [
    {"id": 100, "reactomeId": 1001, "displayName": "AAA makes Cit", "reactionType": "transition", "position": {"x": 150, "y": 150},
     "segments": [{"from": {"x": 140, "y": 150}, "to": {"x": 150, "y": 150}}, {"from": {"x": 150, "y": 150}, "to": {"x": 160, "y": 150}}],
     "inputs": [{"id": 4}], "outputs": [{"id": 3}], "catalysts": [{"id": 1}]},
    {"id": 101, "reactomeId": 1002, "displayName": "complex does X", "reactionType": "reversible", "position": {"x": 250, "y": 150},
     "segments": [], "inputs": [{"id": 3}], "outputs": [{"id": 4}], "catalysts": [{"id": 2}]},
    {"id": 102, "reactomeId": 1003, "displayName": "CCC does Y", "position": {"x": 400, "y": 150}, "segments": [],
     "inputs": [{"id": 3}], "outputs": [{"id": 4}], "catalysts": [{"id": 6}]},
]}
R = mm.parse_reactome_diagram(DIAGRAM, GRAPH)
rk = {n["id"]: n for n in R["nodes"]}
check(rk["r1"]["kind"] == "gene" and rk["r1"]["members"]["uniprot"] == ["P11111"], "a protein node resolves to its UniProt")
check(rk["r2"]["kind"] == "gene" and rk["r2"]["members"]["uniprot"] == ["P33333"],
      "SENSITIVITY: a complex is opened to its leaves (a cycle in the entity graph is guarded) and an isoform suffix is stripped (P33333-2 -> P33333)")
check("chebi" not in rk["r2"]["members"], "SENSITIVITY: the cofactor ChEBI leaf inside a complex is NOT kept as a metabolite member")
check(rk["r3"]["kind"] == "compound" and rk["r3"]["members"]["chebi"] == ["30031"], "a Chemical node resolves to its ChEBI number")
check(rk["r4"]["kind"] == "compound" and rk["r4"]["members"]["chebi"] == ["28757", "4167"],
      "SENSITIVITY: a set whose every leaf is a small molecule is a metabolite node with all its ChEBI members")
check(rk["r5"]["kind"] == "link", "a ProcessNode is a link, never mapped")
check((rk["r1"]["x"], rk["r1"]["y"], rk["r1"]["w"], rk["r1"]["h"]) == (100.0, 100.0, 40.0, 20.0), "position is the node centre; width / height from prop")
rxr = {r["name"]: r for r in R["reactions"]}
check(rxr["AAA makes Cit"]["enzymes"] == ["r1"] and rxr["AAA makes Cit"]["substrates"] == ["r4"] and rxr["AAA makes Cit"]["products"] == ["r3"],
      "a reaction's inputs / outputs / catalysts become substrates / products / enzymes")
check(rxr["AAA makes Cit"]["pos"] == [150.0, 150.0] and len(rxr["AAA makes Cit"]["segments"]) == 2 and rxr["complex does X"]["reversible"] is True,
      "reaction position, drawn segments and reversibility are kept")
Rk = mm.parse_reactome_diagram(DIAGRAM, GRAPH, keep_reaction_db_ids={1001})
check([r["name"] for r in Rk["reactions"]] == ["AAA makes Cit"] and {n["id"] for n in Rk["nodes"]} == {"r1", "r3", "r4"},
      "SENSITIVITY: a sub-pathway keeps ONLY its own reactions and the nodes they touch (no link nodes, no other pathway's entities)")

print("=== Reactome mapping: bridges and the cofactor rule ===")
RGENE = {"acc_to_gene": {"P33333": "hsa:777"}, "entrez_to_gene": {"555": "hsa:777"}, "symbol_to_gene": {}}
rg = pd.DataFrame({"Gene Symbol": ["AAA"], "UniProt ID": ["P11111"], "NCBI Gene ID": [""], "KEGG ID": [""]})
rt = pd.DataFrame({"gene_symbol": ["BBB-GENE", "CCC"], "entrez_id": ["555", ""], "log2FC": [2.0, 1.2], "padj": [0.01, 0.02]})
rp = pd.DataFrame({"accession": ["P44444"], "gene_symbol": ["CCC"], "log2FC": [-1.5], "padj": [0.04]})
rm = pd.DataFrame({"compound_key": ["a", "b", "c"], "Name": ["Citrate", "Iron-sulfur", "Glucose"],
                   "kegg_id": ["", "", ""], "chebi_id": ["CHEBI:30031", "CHEBI:33722", "4167"], "class_level_only": [False] * 3})
RF = {"genomics": mm.gene_features_genomics(rg, RGENE),
      "transcriptomics": mm.gene_features_de("transcriptomics", rt, {"BBB-GENE", "CCC"}, RGENE),
      "proteomics": mm.gene_features_de("proteomics", rp, {"P44444"}, RGENE),
      "metabolomics": mm.compound_features(rm)}
RM = mm.map_omics_to_structure(R, RF, "Reactome")
check(mapped_node(RM, "r1")["status"]["genomics"]["state"] == "hit"
      and mapped_node(RM, "r1")["status"]["genomics"]["members"][0]["via"] == "UniProt",
      "Genomics matches a Reactome protein on its own UniProt")
tm_members = mapped_node(RM, "r2")["status"]["transcriptomics"]["members"]
check(len(tm_members) == 1 and tm_members[0]["label"] == "BBB-GENE" and tm_members[0]["via"] == "KEGG gene -> UniProt",
      "Transcriptomics (no UniProt column) reaches a Reactome complex through Entrez -> KEGG gene -> UniProt, recorded as such")
r6 = mapped_node(RM, "r6")["status"]
check(r6["transcriptomics"]["members"][0]["via"] == "gene symbol (no UniProt known)" and r6["proteomics"]["members"][0]["via"] == "UniProt",
      "a gene with no UniProt AND no KEGG gene is matched by symbol and labelled so; Proteomics matches on UniProt")
# a protein WITH a UniProt must never be matched by its symbol to a different protein that happens to share the name
rs = mm.map_omics_to_structure(R, {"proteomics": mm.gene_features_de(
    "proteomics", pd.DataFrame({"accession": ["P99999"], "gene_symbol": ["CCC"], "log2FC": [1.5], "padj": [0.01]}), {"P99999"}, RGENE)},
    "Reactome")
check(mapped_node(rs, "r6")["status"]["proteomics"]["state"] == "absent",
      "SENSITIVITY: a protein that HAS a UniProt (P99999) is never matched to another protein by a shared gene symbol (CCC)")
check(mapped_node(RM, "r3")["status"]["metabolomics"]["state"] == "hit"
      and mapped_node(RM, "r3")["status"]["metabolomics"]["members"][0]["via"] == "ChEBI ID",
      "Metabolomics matches a Reactome chemical on its ChEBI ID ('CHEBI:30031' normalised)")
check(mapped_node(RM, "r2")["status"]["metabolomics"]["state"] == "absent",
      "SENSITIVITY: the iron-sulfur cofactor in the complex does NOT make the complex a Metabolomics hit")
check(mapped_node(RM, "r4")["status"]["metabolomics"]["state"] == "hit",
      "a set of sugars is a hit when ANY of its ChEBI members is in the list (glucose)")
rlinks = mm.reaction_links(R, RM, ["genomics", "transcriptomics", "proteomics", "metabolomics"])
rl = {r["Reaction"]: r for r in rlinks}
check(set(rl) == {"AAA makes Cit", "complex does X", "CCC does Y"}, f"Reactome links found: {sorted(rl)}")
check(rl["CCC does Y"]["Transcriptomics vs Proteomics direction"] == "discordant",
      "the symbol-matched Transcriptomics hit (up) and the Proteomics hit (down) on one enzyme are 'discordant'")

# ---------------------------------------------------------------------------
print("=== Reactome: another protonation state is bridged by InChIKey skeleton ===")
sk_m = pd.DataFrame({"compound_key": ["c1", "c2", "c3"], "Name": ["Citrate", "Glucose", "Mystery"],
                     "kegg_id": ["", "", ""], "chebi_id": ["CHEBI:30769", "", ""],
                     "inchikey": ["KRKNYBCHXYNGOX-UHFFFAOYSA-N", "WQZGKKKJIJFFOK-GASJEMHNSA-N", ""],
                     "class_level_only": [False, False, False]})
SKF = {"metabolomics": mm.compound_features(sk_m)}
check({f["label"] for f in SKF["metabolomics"]} == {"Citrate", "Glucose"},
      "a compound with an InChIKey but no ChEBI / KEGG ID is still a feature; one with no identifier at all is not")
# Reactome cites citrate(3-) = CHEBI:30031 here, the cross-reference cites citric acid = CHEBI:30769: different ChEBI IDs
SKELS = {"30031": "KRKNYBCHXYNGOX", "4167": "WQZGKKKJIJFFOK", "28757": "BJHIKXHVCXFQLS"}
bridged = mm.map_omics_to_structure(R, SKF, "Reactome", SKELS)
cit = mapped_node(bridged, "r3")["status"]["metabolomics"]
check(cit["state"] == "hit" and cit["members"][0]["via"] == "InChIKey skeleton (another protonation state)",
      "SENSITIVITY: citrate (ChEBI 30769) reaches Reactome's citrate(3-) (ChEBI 30031) through the InChIKey skeleton, and says so")
check(mapped_node(bridged, "r4")["status"]["metabolomics"]["state"] == "hit",
      "...a set is hit through ANY member's skeleton (glucose inside the sugar set)")
unbridged = mm.map_omics_to_structure(R, SKF, "Reactome")
check(mapped_node(unbridged, "r3")["status"]["metabolomics"]["state"] == "absent",
      "SENSITIVITY: without the skeleton lookup the same compounds are NOT matched (exact ChEBI differs) - the bridge is what finds them")
wrong = mm.map_omics_to_structure(R, SKF, "Reactome", {"30031": "AAAAAAAAAAAAAA", "4167": "BBBBBBBBBBBBBB", "28757": "CCCCCCCCCCCCCC"})
check(mapped_node(wrong, "r3")["status"]["metabolomics"]["state"] == "absent",
      "SENSITIVITY: a different skeleton is never a match (the bridge does not conflate different compounds)")
both = mm.map_omics_to_structure(R, {"metabolomics": mm.compound_features(pd.DataFrame({
    "compound_key": ["c1"], "Name": ["Cit"], "kegg_id": [""], "chebi_id": ["30031"], "inchikey": ["KRKNYBCHXYNGOX-UHFFFAOYSA-N"],
    "class_level_only": [False]}))}, "Reactome", SKELS)
check(mapped_node(both, "r3")["status"]["metabolomics"]["members"][0]["via"] == "ChEBI ID",
      "an exact ChEBI match is reported as exact, not as the weaker skeleton bridge")
check(mm.structure_chebi_ids(R) == {"30031", "4167", "28757"}, "the structure's chemicals are listed for the skeleton lookup (not the cofactor inside a complex)")
kegg_bridge = mm.map_omics_to_structure(S, SKF, "KEGG", SKELS)
check(all(n["status"]["metabolomics"]["state"] == "absent" for n in kegg_bridge),
      "KEGG never uses the skeleton bridge (it matches on the C-number)")

print("=== pathway picker ===")
inputs = {
    "genomics": {"term_map": {"P1": {"name": "Path one", "genes": {"a", "b", "c"}}, "P2": {"name": "Path two", "genes": {"d", "e"}},
                              "P3": {"name": "Untouched", "genes": {"x", "y"}},
                              "P4": {"name": "Big one-omics", "genes": {"q1", "q2", "q3", "q4", "q5"}}},
                 "list_set": {"a", "b", "d", "q1", "q2", "q3", "q4", "q5"},
                 "universe_set": {"a", "b", "c", "d", "e", "x", "y", "q1", "q2", "q3", "q4", "q5"}},
    "metabolomics": {"term_map": {"P1": {"name": "Path one", "genes": {"m1", "m2"}}, "P3": {"name": "Untouched", "genes": {"m5"}}},
                     "list_set": {"m1"}, "universe_set": {"m1", "m2", "m5"}},
}
picker = mm.build_pathway_picker(inputs, ["genomics", "metabolomics"], {"P1": 0.01})
check([r["Term ID"] for r in picker] == ["P1", "P4", "P2"],
      "SENSITIVITY: ranked by number of omics hitting FIRST (P1: 2 omics, 3 hits) then total hits (P4: 1 omics, 5 hits, before P2: 1 omics, 1 hit); a pathway nobody hit is not listed")
check(picker[0]["Omics with a hit"] == 2 and picker[0]["Total hits"] == 3 and picker[0]["Genomics: hit / annotated"] == "2 / 3"
      and picker[0]["Metabolomics: hit / annotated"] == "1 / 2", "per-omics 'hit / annotated' counts")
check(picker[0]["Joint combined FDR"] == 0.01 and picker[1]["Joint combined FDR"] is None and picker[2]["Joint combined FDR"] is None,
      "the joint FDR is shown only where a joint run produced one")
check(picker[2]["Metabolomics: hit / annotated"] == "0 / 0", "an omics with no annotation in a pathway shows 0 / 0")
check(mm.picker_columns(["genomics", "proteomics"])[-1] == "Joint combined FDR" and "Proteomics: hit / annotated" in mm.picker_columns(["genomics", "proteomics"]),
      "picker columns follow the selected omics")

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
print("=== KEGG overview maps: genes drawn as LINES ===")
KGML_LINES = """<?xml version="1.0"?>
<pathway name="path:hsa88888" org="hsa" number="88888" title="Line map">
  <entry id="1" name="hsa:10" type="gene" reaction="rn:R1"><graphics name="LNA, X" type="line" coords="0,0,0,0,100,0,100,0,100,50,100,50"/></entry>
  <entry id="2" name="hsa:20 hsa:21" type="gene" reaction="rn:R2"><graphics name="LNB" type="line" coords="100,50,100,100,200,100"/></entry>
  <entry id="3" name="ko:K1" type="ortholog" reaction="rn:R3"><graphics name="K1" type="line" coords="200,100,300,100"/></entry>
  <entry id="4" name="cpd:C00001" type="compound"><graphics name="C00001" type="circle" x="0" y="0" width="8" height="8"/></entry>
  <entry id="5" name="cpd:C00002" type="compound"><graphics name="C00002" type="circle" x="100" y="50" width="8" height="8"/></entry>
  <entry id="6" name="cpd:C00003" type="compound"><graphics name="C00003" type="circle" x="200" y="100" width="8" height="8"/></entry>
  <entry id="7" name="cpd:C00004" type="compound"><graphics name="C00004" type="circle" x="300" y="100" width="8" height="8"/></entry>
  <reaction id="1" name="rn:R1" type="irreversible"><substrate id="4" name="cpd:C00001"/><product id="5" name="cpd:C00002"/></reaction>
  <reaction id="2" name="rn:R2" type="irreversible"><substrate id="5" name="cpd:C00002"/><product id="6" name="cpd:C00003"/></reaction>
  <reaction id="3" name="rn:R3" type="irreversible"><substrate id="6" name="cpd:C00003"/><product id="7" name="cpd:C00004"/></reaction>
</pathway>"""
SL = mm.parse_kgml_full(KGML_LINES)
check(SL["line_graphics"] is True and not mm.parse_kgml_full(KGML)["line_graphics"],
      "a map whose genes are lines says so; a box map does not")
nl = {n["id"]: n for n in SL["nodes"]}
check(nl["k1"]["path"] == [0, 0, 100, 0, 100, 50] and (nl["k1"]["x"], nl["k1"]["y"]) == (100.0, 0.0),
      "SENSITIVITY: a line gene keeps its polyline (KGML's repeated corner points removed) and sits at the middle vertex - "
      "not at (0, 0), where every gene of KEGG's global maps used to pile up")
check(nl["k2"]["path"] == [100, 50, 100, 100, 200, 100] and (nl["k2"]["x"], nl["k2"]["y"]) == (100.0, 100.0)
      and nl["k2"]["members"]["kegg_genes"] == ["hsa:20", "hsa:21"] and nl["k2"]["label"] == "LNB",
      "...with its genes and label")
check("k3" not in nl and SL["decor_lines"] == [[200, 100, 300, 100]],
      "a line ORTHOLOG is decoration (kept as a grey shape), never a node a hit could land on")
check(all(n["kind"] != "gene" or "path" in n for n in SL["nodes"]) and len(SL["nodes"]) == 6,
      "every gene of a line map has a path; 2 genes + 4 compounds are the nodes")
check(next(r for r in SL["reactions"] if r["name"] == "R3")["enzymes"] == [],
      "a reaction whose enzyme is a (dropped) ortholog line has no enzyme node, not a dangling id")
LF = {"genomics": mm.gene_features_genomics(pd.DataFrame({"Gene Symbol": ["LNB"], "UniProt ID": [""], "NCBI Gene ID": [""],
                                                         "KEGG ID": ["hsa:21"]}))}
ML = mm.map_omics_to_structure(SL, LF, "KEGG")
check(mapped_node(ML, "k2")["status"]["genomics"]["state"] == "hit" and mapped_node(ML, "k1")["status"]["genomics"]["state"] == "absent"
      and "path" in mapped_node(ML, "k2"), "a line gene is matched like a box gene, and keeps its path through the mapping")

print("=== KEGG signalling maps: groups and relations ===")
KGML_SIG = """<?xml version="1.0"?>
<pathway name="path:hsa77777" org="hsa" number="77777" title="Signalling map">
  <entry id="1" name="hsa:1" type="gene"><graphics name="AAA, A1" type="rectangle" x="100" y="100" width="46" height="17"/></entry>
  <entry id="2" name="hsa:2" type="gene"><graphics name="BBB" type="rectangle" x="200" y="100" width="46" height="17"/></entry>
  <entry id="3" name="hsa:3" type="gene"><graphics name="CCC" type="rectangle" x="300" y="100" width="46" height="17"/></entry>
  <entry id="4" name="hsa:4" type="gene"><graphics name="DDD" type="rectangle" x="360" y="100" width="46" height="17"/></entry>
  <entry id="5" name="undefined" type="group"><graphics type="rectangle" x="335" y="100" width="130" height="30"/><component id="3"/><component id="4"/></entry>
  <entry id="6" name="cpd:C00001" type="compound"><graphics name="C00001" type="circle" x="100" y="200" width="8" height="8"/></entry>
  <entry id="7" name="hsa:1" type="gene"><graphics name="AAA" type="rectangle" x="500" y="300" width="46" height="17"/></entry>
  <relation entry1="1" entry2="2" type="PPrel"><subtype name="activation" value="--&gt;"/></relation>
  <relation entry1="2" entry2="5" type="PPrel"><subtype name="inhibition" value="--|"/><subtype name="phosphorylation" value="+p"/></relation>
  <relation entry1="1" entry2="6" type="PCrel"><subtype name="compound" value="6"/></relation>
  <relation entry1="1" entry2="99" type="PPrel"><subtype name="activation" value="--&gt;"/></relation>
  <relation entry1="5" entry2="5" type="maplink"><subtype name="compound" value="x"/></relation>
</pathway>"""
SS = mm.parse_kgml_full(KGML_SIG)
ns = {n["id"]: n for n in SS["nodes"]}
check(ns["k5"]["kind"] == "group" and ns["k5"]["components"] == ["k3", "k4"] and ns["k5"]["members"]["kegg_genes"] == ["hsa:3", "hsa:4"]
      and ns["k5"]["label"] == "CCC+DDD",
      "SENSITIVITY: a group (a drawn complex) holds its components' genes and is labelled from them (KEGG gives it no name)")
check(len(SS["relations"]) == 3 and {r["rtype"] for r in SS["relations"]} == {"pprel", "pcrel"},
      f"relations are read; one to a missing entry and a map link are skipped ({len(SS['relations'])} kept)")
check(SS["relations"][1]["types"] == ["inhibition", "phosphorylation"] and SS["relations"][0]["types"] == ["activation"]
      and SS["relations"][2]["types"] == ["pcrel"],
      f"every subtype is kept in order; the 'compound' subtype is not a type; none falls back to the relation type ({[r['types'] for r in SS['relations']]})")
gdf = pd.DataFrame({"Gene Symbol": ["AAA", "BBB", "DDD"], "UniProt ID": ["", "", ""], "NCBI Gene ID": ["", "", ""],
                    "KEGG ID": ["hsa:1", "hsa:2", "hsa:4"]})
MS = mm.map_omics_to_structure(SS, {"genomics": mm.gene_features_genomics(gdf)}, "KEGG")
check(mapped_node(MS, "k5")["status"]["genomics"]["state"] == "hit"
      and [m["label"] for m in mapped_node(MS, "k5")["status"]["genomics"]["members"]] == ["DDD"],
      "SENSITIVITY: a group is hit if ANY component is (DDD hit, CCC not)")
ilinks = mm.interaction_links(SS, MS, ["genomics"])
check(len(ilinks) == 2 and [(r["From"], r["Relation"]) for r in ilinks] == [("AAA", "activation"), ("BBB", "inhibition, phosphorylation")],
      f"interaction links need BOTH partners hit: AAA->BBB and BBB->group; AAA->compound (compound not hit) is not one ({[(r['From'], r['To']) for r in ilinks]})")
check(ilinks[1]["To"] == "CCC+DDD (DDD)" and ilinks[1]["To hit by"] == "Genomics" and ilinks[0]["From hit by"] == "Genomics",
      "a group partner is named with the member that was hit; each side says which omics hit it")
check(mm.interaction_links(SS, mm.map_omics_to_structure(SS, {"genomics": mm.gene_features_genomics(gdf.iloc[[0]])}, "KEGG"), ["genomics"]) == [],
      "SENSITIVITY: with only one partner hit there is no interaction link")
check(mm.interaction_links({"relations": []}, MS, ["genomics"]) == [] and mm.interaction_links({}, MS, ["genomics"]) == [],
      "a map with no relations has none")
check(all(r["Pathway node"] != "CCC+DDD" for r in mm.member_rows(MS, ["genomics"])),
      "a group is not a row of the member table (its components have their own rows)")
check(mm.member_rows(MS, ["genomics"])[0]["Box ID"].startswith("k"), "member rows carry the box id, so repeated boxes can be told apart")
sm_s = mm.summarize_mapping(MS, [], ["genomics"], ilinks)
check(sm_s["n_nodes"] == 6 and sm_s["n_interactions"] == 2 and sm_s["per_omics"]["genomics"]["nodes_hit"] == 4,
      f"a group is not counted as a box of its own: 6 boxes, and 4 of them hit (AAA's two boxes, BBB, DDD) - got "
      f"{sm_s['n_nodes']} / {sm_s['per_omics']['genomics']['nodes_hit']}")
check(sm_s["per_omics"]["genomics"]["members_hit"] == 3,
      f"SENSITIVITY: AAA drawn in two boxes is ONE member hit (AAA, BBB, DDD = 3), not two ({sm_s['per_omics']['genomics']['members_hit']})")

print("=== pathway-class colours of KEGG's global maps ===")
KGML_COLOR = """<?xml version="1.0"?>
<pathway name="path:hsa01100" org="hsa" number="01100" title="Colour map">
  <entry id="1" name="hsa:10" type="gene"><graphics name="A1" type="line" fgcolor="#b3b3e6" coords="0,0,100,0"/></entry>
  <entry id="2" name="hsa:20" type="gene"><graphics name="B1" type="line" fgcolor="#B3B3E6" coords="0,10,100,10"/></entry>
  <entry id="3" name="hsa:30" type="gene"><graphics name="C1" type="line" fgcolor="#80CCCC" coords="0,20,100,20"/></entry>
  <entry id="4" name="hsa:40" type="gene"><graphics name="D1" type="line" fgcolor="#FFCC66" coords="0,30,100,30"/></entry>
  <entry id="5" name="hsa:50" type="gene"><graphics name="E1" type="line" fgcolor="#6C63F6" coords="0,40,100,40"/></entry>
  <entry id="6" name="hsa:60" type="gene"><graphics name="F1" type="rectangle" fgcolor="#000000" x="10" y="10" width="46" height="17"/></entry>
  <entry id="7" name="path:hsa00010" type="map"><graphics name="Glycolysis" fgcolor="none" type="roundrectangle" x="50" y="50" width="40" height="20"/></entry>
  <entry id="8" name="path:hsa01100" type="map"><graphics name="Metabolic pathways" fgcolor="none" type="roundrectangle" x="60" y="60" width="40" height="20"/></entry>
</pathway>"""
SC = mm.parse_kgml_full(KGML_COLOR)
cn = {n["id"]: n for n in SC["nodes"]}
check(cn["k1"]["color"] == "#B3B3E6" and cn["k2"]["color"] == "#B3B3E6" and cn["k3"]["color"] == "#80CCCC",
      "a line gene keeps KEGG's colour for it (upper-cased, so '#b3b3e6' and '#B3B3E6' are one colour)")
check("color" not in cn["k6"], "a box gene has no line colour")
check(cn["k7"]["pathway"] == "00010" and cn["k8"]["pathway"] == "01100", "a pathway label keeps the map number it links to")
BRITE = """+C	Map number
!
AMetabolism
B  Global and overview maps
C    01100  Metabolic pathways
B  Carbohydrate metabolism
C    00010  Glycolysis / Gluconeogenesis
C    00020  Citrate cycle (TCA cycle)
B  Lipid metabolism
C    00061  Fatty acid biosynthesis
B  Amino acid metabolism
C    00250  Alanine, aspartate and glutamate metabolism
AGenetic Information Processing
B  Translation
C    03010  Ribosome
AHuman Diseases
B  Neurodegenerative disease
C    05010  Alzheimer disease
"""
classes_ = mm.parse_brite_classes(BRITE)
check(classes_ == {"01100": "Global and overview maps", "00010": "Carbohydrate metabolism", "00020": "Carbohydrate metabolism",
                   "00061": "Lipid metabolism", "00250": "Amino acid metabolism"},
      f"KEGG's classification is read at its B level, for the METABOLISM branch only ({classes_})")
check("03010" not in classes_ and "05010" not in classes_,
      "SENSITIVITY: pathways of other branches (Ribosome, Alzheimer disease) are left out - a metabolic map's genes also sit on disease "
      "pathways, which once named a colour 'Neurodegenerative disease'")
g2m = {"hsa:10": {"00010", "00020", "01100", "05010"}, "hsa:20": {"00020"}, "hsa:30": {"00061"}, "hsa:40": {"01100"},
       "hsa:50": {"00010", "00061"}}
cc = mm.assign_color_classes(SC["nodes"], g2m, classes_)
check(cc["#B3B3E6"] == "Carbohydrate metabolism" and cc["#80CCCC"] == "Lipid metabolism", "a colour is named by the class its genes' pathways belong to")
check(cc["#FFCC66"] == "Other",
      "SENSITIVITY: a colour whose only pathway is the overview map itself is 'Other' (the global map is not a class a gene belongs to)")
check(cc["#6C63F6"] == "Carbohydrate metabolism",
      "SENSITIVITY: a tie (one pathway of each of two classes) goes to the alphabetically first class, never to dict order")
check("#000000" not in cc and len(cc) == 4, "a box gene's colour is not a line colour")
for n in SC["nodes"]:
    n["status"] = {}
SC["nodes"][0]["status"] = {"genomics": {"state": "hit", "direction": None, "members": []}}
leg = mm.class_legend(SC["nodes"], cc)
by_name = {r["name"]: r for r in leg}
check([r["name"] for r in leg] == ["Carbohydrate metabolism", "Lipid metabolism", "Other"],
      f"the legend lists each class once, sorted ({[r['name'] for r in leg]})")
check(by_name["Carbohydrate metabolism"]["n_lines"] == 3 and by_name["Carbohydrate metabolism"]["n_hit"] == 1
      and by_name["Lipid metabolism"]["n_hit"] == 0,
      "...with how many of its lines there are and how many were hit")
check(by_name["Carbohydrate metabolism"]["color"] == "#B3B3E6",
      "a class drawn in several KEGG colours takes its most frequent one as the default (#B3B3E6 on 2 lines, #6C63F6 on 1)")

print("=== Reactome overview pathways ===")
OV_DIAGRAM = {"displayName": "Disease", "nodes": [
    {"id": 1, "reactomeId": 111, "displayName": "Diseases of metabolism", "renderableClass": "ProcessNode", "position": {"x": 100, "y": 100}, "prop": {"width": 80, "height": 30}},
    {"id": 2, "reactomeId": 222, "displayName": "Diseases of immune system", "renderableClass": "ProcessNode", "position": {"x": 300, "y": 100}, "prop": {"width": 80, "height": 30}},
], "edges": []}
OV = mm.parse_reactome_diagram(OV_DIAGRAM, {"nodes": []})
check(mm.is_overview(OV) and [n["reactome_id"] for n in OV["nodes"]] == [111, 222],
      "a diagram of only ProcessNodes is an overview, and each link keeps the pathway it stands for")
check(not mm.is_overview(R) and not mm.is_overview({"nodes": []}) and not mm.is_overview(mm.parse_kgml_full(KGML)),
      "SENSITIVITY: a real diagram, an empty one and a KEGG map are not overviews")
check(not mm.is_overview(mm.parse_reactome_diagram(dict(OV_DIAGRAM, nodes=OV_DIAGRAM["nodes"] + [DIAGRAM["nodes"][0]]), GRAPH)),
      "one real entity among the links makes it a pathway that can be mapped")
orows = mm.overview_rows([{"stId": "R-1", "name": "Zeta"}, {"stId": "R-2", "name": "Alpha"}, {"stId": "R-3", "name": "Beta"}],
                         [{"Term ID": "R-3", "Term Name": "Beta", "Omics with a hit": 2, "Total hits": 3, "Genomics: hit / annotated": "2 / 9",
                           "Joint combined FDR": 0.01},
                          {"Term ID": "R-1", "Term Name": "Zeta", "Omics with a hit": 1, "Total hits": 8, "Genomics: hit / annotated": "8 / 20",
                           "Joint combined FDR": None}], ["genomics"])
check([r["Term ID"] for r in orows] == ["R-3", "R-1", "R-2"] and orows[2]["Total hits"] == 0 and orows[2]["Genomics: hit / annotated"] == "0 / -"
      and orows[0]["Genomics: hit / annotated"] == "2 / 9",
      "sub-pathways ranked by omics count then hits; one the picker does not list has no hit (0 / -)")

print("=== fetchers: caching and the no-diagram fallback (HTTP mocked) ===")
store = {}
calls = {"kegg": 0, "json": []}
mr._gm._load_universe_cache = lambda db, key: store.get((db, key))
mr._gm._save_universe_cache = lambda db, key, data: store.__setitem__((db, key), data)


def fake_kegg(status, text):
    def _get(url, verbose=False):
        calls["kegg"] += 1
        return status, text
    return _get


mr._gm._kegg_safe_get = fake_kegg("error", "boom")
out = mr.fetch_pathway_structure("KEGG", "hsa99999")
check(out["status"] == "error" and out["structure"] is None and not store, "SENSITIVITY: a KEGG transport error is reported and NEVER cached")
mr._gm._kegg_safe_get = fake_kegg("nodata", "")
out = mr.fetch_pathway_structure("KEGG", "hsa99998")
check(out["status"] == "nodata" and ("MultiomicsMap_KEGG_kgml_v4", "hsa99998") in store, "a definite 'no map' is cached")
mr._gm._kegg_safe_get = fake_kegg("ok", KGML)
out = mr.fetch_pathway_structure("KEGG", "hsa99999")
n_before = calls["kegg"]
again = mr.fetch_pathway_structure("KEGG", "hsa99999")
check(out["status"] == "ok" and len(out["structure"]["nodes"]) == 10 and again["status"] == "ok" and calls["kegg"] == n_before,
      "a good map is parsed, cached, and the second request makes no network call")
mr._gm._kegg_safe_get = fake_kegg("ok", "<not xml")
check(mr.fetch_pathway_structure("KEGG", "hsa99997", force_refresh=True)["status"] == "error", "an unreadable KGML is an error, not a crash")
check(mr._reactome_reaction_db_ids([{"dbId": 5}, 7, {"dbId": 9, "x": 1}, "junk"]) == {5, 7, 9},
      "contained-event ids: Reactome's bare-integer repeats count as well as full objects")

chain_pathway = {"hasDiagram": False, "displayName": "Sub", "stId": "R-HSA-5"}
responses = {
    "query/R-HSA-5": ("ok", chain_pathway),
    "event/R-HSA-5/ancestors": ("ok", [[{"stId": "R-HSA-5", "hasDiagram": False},
                                        {"stId": "R-HSA-6", "displayName": "Parent", "hasDiagram": True},
                                        {"stId": "R-HSA-7", "displayName": "Grand", "hasDiagram": True}]]),
    "pathway/R-HSA-5/containedEvents": ("ok", [{"dbId": 1001}, 1002]),
    "R-HSA-6.json": ("ok", DIAGRAM), "R-HSA-6.graph.json": ("ok", GRAPH),
}


def fake_json(url):
    calls["json"].append(url)
    for tail, resp in responses.items():
        if url.endswith(tail):
            return resp
    return "nodata", None


mr._http_get_json = fake_json
out = mr.fetch_pathway_structure("Reactome", "R-HSA-5")
check(out["status"] == "ok" and out["restricted_to"] == "Parent (R-HSA-6)",
      "SENSITIVITY: a pathway with no diagram is drawn from its NEAREST ancestor that has one (not the grand-parent)")
check({r["name"] for r in out["structure"]["reactions"]} == {"AAA makes Cit", "complex does X"} and out["structure"]["title"] == "Sub",
      "...restricted to the pathway's own reactions, titled with the pathway itself")
check(any(u.endswith("R-HSA-6.json") for u in calls["json"]) and not any(u.endswith("R-HSA-7.json") for u in calls["json"]),
      "only the nearest ancestor's diagram is downloaded")
n_calls = len(calls["json"])
mr.fetch_pathway_structure("Reactome", "R-HSA-5")
check(len(calls["json"]) == n_calls, "a fetched Reactome structure is cached")
responses["query/R-HSA-8"] = ("error", "down")
check(mr.fetch_pathway_structure("Reactome", "R-HSA-8")["status"] == "error" and ("MultiomicsMap_Reactome_diagram_v1", "R-HSA-8") not in store,
      "a Reactome transport error is never cached")
responses["query/R-HSA-9"] = ("ok", {"hasDiagram": False})
responses["event/R-HSA-9/ancestors"] = ("ok", [[{"stId": "R-HSA-9", "hasDiagram": False}]])
o9 = mr.fetch_pathway_structure("Reactome", "R-HSA-9")
check(o9["status"] == "nodata" and "no diagram" in o9["note"], "no diagram anywhere above: an explicit 'nodata' note, no invented map")
try:
    mr.fetch_pathway_structure("WikiPathways", "x")
    check(False, "an unsupported database raises")
except ValueError:
    check(True, "an unsupported database raises (MetaCyc is not offered)")

print("=== fetchers: overview pathways ===")
responses.update({
    "query/R-HSA-60": ("ok", {"hasDiagram": True, "displayName": "Disease"}),
    "R-HSA-60.json": ("ok", OV_DIAGRAM), "R-HSA-60.graph.json": ("ok", {"nodes": []}),
    "query/111": ("ok", {"stId": "R-HSA-111", "displayName": "Diseases of metabolism"}),
    "query/222": ("ok", {"stId": "R-HSA-222", "displayName": "Diseases of immune system"}),
})
ov = mr.fetch_pathway_structure("Reactome", "R-HSA-60")
check(ov["status"] == "overview" and ov["subpathways"] == [{"stId": "R-HSA-111", "name": "Diseases of metabolism"},
                                                         {"stId": "R-HSA-222", "name": "Diseases of immune system"}]
      and "only as links to its 2 sub-pathways" in ov["note"],
      "SENSITIVITY: an overview pathway is reported as one, with its sub-pathways resolved to stable ids and names")
n_calls = len(calls["json"])
mr.fetch_pathway_structure("Reactome", "R-HSA-60")
check(len(calls["json"]) == n_calls and ("MultiomicsMap_Reactome_diagram_v2", "R-HSA-60") in store, "an overview is cached")
responses["query/R-HSA-61"] = ("ok", {"hasDiagram": True, "displayName": "Disease 2"})
responses["R-HSA-61.json"] = ("ok", OV_DIAGRAM)
responses["R-HSA-61.graph.json"] = ("ok", {"nodes": []})
responses["query/111"] = ("error", "down")
ov2 = mr.fetch_pathway_structure("Reactome", "R-HSA-61")
check(ov2["status"] == "error" and ("MultiomicsMap_Reactome_diagram_v2", "R-HSA-61") not in store,
      "SENSITIVITY: if a sub-pathway cannot be looked up, the overview is an error and is NOT cached (it would be a half list for good)")
responses["query/111"] = ("ok", {"stId": "R-HSA-111", "displayName": "Diseases of metabolism"})
responses["query/222"] = ("ok", {"stId": "R-HSA-111", "displayName": "Duplicate"})
ov3 = mr.fetch_pathway_structure("Reactome", "R-HSA-61", force_refresh=True)
check([x["stId"] for x in ov3["subpathways"]] == ["R-HSA-111"], "a sub-pathway reached through two link boxes is listed once")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
