# -*- coding: utf-8 -*-
"""Unit tests for Pathway Topology Analysis (Metabolomics_module.
run_pathway_topology and its building blocks), offline. Only KEGG's network
edges are faked - the enrichment path itself (build_kegg_compound_universe,
_eligible_keys, _drop_zero_overlap, _run_fisher_for_namespace /
_run_wilcoxon_for_namespace) is the REAL code test_metabolomics_ora.py
already exercises, reused here unchanged.

THE HEADLINE GUARANTEE under test: run_pathway_topology calls the exact same
building blocks run_ora / run_quantitative_enrichment use for their own KEGG
database entry, so a topology run and the corresponding ORA/Quantitative KEGG
run must report IDENTICAL p-values/FDRs for the same input and config - only
the added Impact/Matched-in-graph/Graph-nodes/Topology-note columns are new.
Checked directly against a live run_ora call on the same fixture.

THE SECOND THING under test: Brandes' betweenness (this module's own
implementation, no networkx) is checked against a brute-force shortest-path
enumeration on random small directed graphs - the same discipline used for
every other from-scratch statistic in this codebase (see e.g.
test_regulon_engine.py's closed-form-vs-brute-force check).

THE THIRD THING under test: KEGG's reaction graph is organism-specific in its
EDGES even though compound membership is not (verified live 2026-09-29:
hsa00020/sce00020/eco00020 share 23 compounds but have 22/21/21 reactions) -
so fetch_pathway_graph must be fetched and cached PER organism PER map, never
reused across organisms the way build_kegg_organism_pathway_index's map INDEX
correctly is.
"""
import itertools
import os
import random
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Metabolomics_module as mm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# ---------------------------------------------------------------------------
# parse_kegg_kgml
# ---------------------------------------------------------------------------
print("=== parse_kegg_kgml: nodes, reversible/irreversible edges, glycans, zero-reaction maps ===")

_KGML_BASIC = """<?xml version="1.0"?>
<pathway name="path:hsa00020" org="hsa" number="00020">
  <entry id="1" name="cpd:C00022" type="compound">
    <graphics name="C00022" x="100" y="200" width="8" height="8"/>
  </entry>
  <entry id="2" name="cpd:C00024" type="compound">
    <graphics name="C00024" x="150" y="220" width="8" height="8"/>
  </entry>
  <entry id="3" name="cpd:C00036" type="compound">
    <graphics name="C00036" x="180" y="250" width="8" height="8"/>
  </entry>
  <entry id="4" name="gl:G00001" type="compound">
    <graphics name="G00001" x="300" y="100" width="8" height="8"/>
  </entry>
  <reaction id="10" name="rn:R00001" type="reversible">
    <substrate id="1" name="cpd:C00022"/>
    <product id="2" name="cpd:C00024"/>
  </reaction>
  <reaction id="11" name="rn:R00002" type="irreversible">
    <substrate id="2" name="cpd:C00024"/>
    <substrate id="4" name="gl:G00001"/>
    <product id="3" name="cpd:C00036"/>
  </reaction>
</pathway>
"""
graph = mm.parse_kegg_kgml(_KGML_BASIC)
check(set(graph["nodes"]) == {"C00022", "C00024", "C00036", "G00001"},
      f"every compound entry becomes a node, glycan (gl:) kept (got {sorted(graph['nodes'])})")
check(graph["nodes"]["C00022"]["x"] == 100.0 and graph["nodes"]["C00022"]["y"] == 200.0,
      "KGML x/y coordinates are read through")
check(graph["n_reactions"] == 2, f"n_reactions counts <reaction> elements (got {graph['n_reactions']})")
edge_set = {tuple(e) for e in graph["edges"]}
check(("C00022", "C00024") in edge_set and ("C00024", "C00022") in edge_set,
      "a REVERSIBLE reaction adds both directions")
check(("C00024", "C00036") in edge_set and ("G00001", "C00036") in edge_set,
      "a multi-substrate reaction adds one edge per substrate->product pair")
check(("C00036", "C00024") not in edge_set and ("C00036", "G00001") not in edge_set,
      "an IRREVERSIBLE reaction does NOT add the reverse direction")

_KGML_DUP_ENTRY = """<?xml version="1.0"?>
<pathway name="path:test" org="hsa" number="00000">
  <entry id="1" name="cpd:C00001 cpd:C00002" type="compound">
    <graphics name="C00001, C00002" x="50" y="60" width="8" height="8"/>
  </entry>
</pathway>
"""
graph_dup = mm.parse_kegg_kgml(_KGML_DUP_ENTRY)
check(set(graph_dup["nodes"]) == {"C00001", "C00002"} and graph_dup["nodes"]["C00001"]["x"] == 50.0,
      "an entry naming several ids adds each as its own node, sharing that entry's coordinates")

_KGML_ZERO_REACTIONS = """<?xml version="1.0"?>
<pathway name="path:hsa00190" org="hsa" number="00190">
  <entry id="1" name="cpd:C00002" type="compound"><graphics name="C00002" x="1" y="1"/></entry>
  <entry id="2" name="cpd:C00008" type="compound"><graphics name="C00008" x="2" y="2"/></entry>
</pathway>
"""
graph_zero = mm.parse_kegg_kgml(_KGML_ZERO_REACTIONS)
check(len(graph_zero["nodes"]) == 2 and graph_zero["edges"] == [] and graph_zero["n_reactions"] == 0,
      "a real pathway with entries but no <reaction> elements returns nodes with an empty edge list, not an error")


# ---------------------------------------------------------------------------
# _brandes_betweenness vs brute-force shortest-path enumeration
# ---------------------------------------------------------------------------
print("=== _brandes_betweenness: matches brute-force enumeration on random digraphs ===")


def _brute_force_betweenness(nodes, adjacency):
    from collections import deque
    centrality = {v: 0.0 for v in nodes}
    n = len(nodes)
    if n < 3:
        return centrality
    for s in nodes:
        dist = {v: None for v in nodes}
        dist[s] = 0
        preds = {v: [] for v in nodes}
        queue = deque([s])
        while queue:
            v = queue.popleft()
            for w in adjacency.get(v, []):
                if dist[w] is None:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:
                    preds[w].append(v)
        for t in nodes:
            if t == s or dist[t] is None:
                continue
            paths = []

            def backtrack(node, tail):
                if node == s:
                    paths.append([s] + tail[::-1])
                    return
                for p in preds[node]:
                    backtrack(p, tail + [node])

            backtrack(t, [])
            for path in paths:
                for mid in path[1:-1]:
                    centrality[mid] += 1.0 / len(paths)
    norm = (n - 1) * (n - 2)
    if norm > 0:
        for v in centrality:
            centrality[v] /= norm
    return centrality


random.seed(20260929)
_all_match = True
for _trial in range(25):
    _n = random.randint(4, 8)
    _nodes = [f"n{i}" for i in range(_n)]
    _adj = {v: [] for v in _nodes}
    for a, b in itertools.permutations(_nodes, 2):
        if random.random() < 0.35:
            _adj[a].append(b)
    got = mm._brandes_betweenness(_nodes, _adj)
    expected = _brute_force_betweenness(_nodes, _adj)
    if any(abs(got[v] - expected[v]) > 1e-9 for v in _nodes):
        _all_match = False
        check(False, f"trial {_trial} mismatch: got={got} expected={expected}")
        break
check(_all_match, "Brandes betweenness matches brute-force shortest-path enumeration on 25 random digraphs")

# A graph with < 3 nodes has no meaningful betweenness (no node can lie
# BETWEEN two others) - every centrality is 0, not undefined/NaN.
check(mm._brandes_betweenness(["a", "b"], {"a": ["b"], "b": []}) == {"a": 0.0, "b": 0.0},
      "a 2-node graph has zero betweenness for every node")


# ---------------------------------------------------------------------------
# node_importance / pathway_impact
# ---------------------------------------------------------------------------
print("=== node_importance / pathway_impact ===")
_chain_graph = {"nodes": {"A": {}, "B": {}, "C": {}}, "edges": [["A", "B"], ["B", "C"]]}
imp_out = mm.node_importance(_chain_graph, "outdegree")
check(imp_out == {"A": 1.0, "B": 1.0, "C": 0.0}, f"out-degree counts outgoing edges only (got {imp_out})")
imp_bw = mm.node_importance(_chain_graph, "betweenness")
# Normalised by (n-1)(n-2) = 2: the middle node lies on the ONE shortest
# path (A,C), so its raw count 1 becomes 1/2 - see _brandes_betweenness's
# own docstring for this exact worked example.
check(imp_bw["B"] == 0.5 and imp_bw["A"] == 0.0 and imp_bw["C"] == 0.0,
      f"the middle node of a 3-chain has normalised betweenness 0.5, the ends 0 (got {imp_bw})")

impact_full, wd_full = mm.pathway_impact(imp_out, {"A", "B", "C"})
check(impact_full == 1.0 and wd_full, "impact is 1.0 when every node is matched")
impact_none, wd_none = mm.pathway_impact(imp_out, set())
check(impact_none == 0.0 and wd_none, "impact is 0.0 (well-defined) when nothing is matched but the graph has weight")
impact_a, wd_a = mm.pathway_impact(imp_out, {"A"})
check(abs(impact_a - 0.5) < 1e-9 and wd_a, f"impact is the matched share of total importance (got {impact_a})")

_zero_graph = {"nodes": {"A": {}, "B": {}}, "edges": []}
imp_zero = mm.node_importance(_zero_graph, "betweenness")
impact_zero, wd_zero = mm.pathway_impact(imp_zero, {"A"})
check(impact_zero == 0.0 and wd_zero is False,
      "SENSITIVITY: a zero-reaction map's impact is reported as 0.0 but flagged NOT well-defined, "
      "distinct from a genuinely-computed zero")
check(mm.node_importance({"nodes": {}}, "betweenness") == {}, "an empty graph has no importances, not an error")
check(mm.node_importance({"nodes": {"A": {}}}, "outdegree") == {"A": 0.0}, "a single node has out-degree 0")


# ---------------------------------------------------------------------------
# fetch_pathway_graph: caching, nodata, error, and organism-specificity
# ---------------------------------------------------------------------------
print("=== fetch_pathway_graph: caching, nodata vs error, per-organism (not per-map) ===")

_FAKE_CACHE_STORE = {}
mm._load_universe_cache = lambda db, key: _FAKE_CACHE_STORE.get((db, key))
mm._save_universe_cache = lambda db, key, data: _FAKE_CACHE_STORE.__setitem__((db, key), data)

_KGML_HSA = _KGML_BASIC  # 2 reactions, defined above
_KGML_SCE = """<?xml version="1.0"?>
<pathway name="path:sce00020" org="sce" number="00020">
  <entry id="1" name="cpd:C00022" type="compound"><graphics name="C00022" x="10" y="10"/></entry>
  <entry id="2" name="cpd:C00024" type="compound"><graphics name="C00024" x="20" y="20"/></entry>
</pathway>
"""  # 0 reactions - deliberately different from the human graph above

_fetch_calls = []  # KGML fetches only - list/compound (the name index) is a separate, independently-cached artifact.
_FAKE_COMPOUND_NAMES = "C00022\tPyruvate; 2-Oxopropanoic acid\nC00024\tAcetyl-CoA\nC00036\tOxaloacetate\n"


def _fake_kegg_get_kgml(url, verbose=False, **kwargs):
    if url.endswith("/list/compound"):
        return "ok", _FAKE_COMPOUND_NAMES
    _fetch_calls.append(url)
    if url.endswith("/get/hsa00020/kgml"):
        return "ok", _KGML_HSA
    if url.endswith("/get/sce00020/kgml"):
        return "ok", _KGML_SCE
    if url.endswith("/get/hsa99999/kgml"):
        return "nodata", ""
    if url.endswith("/get/hsa88888/kgml"):
        return "error", ""
    return "error", ""


mm._kegg_safe_get = _fake_kegg_get_kgml

g_hsa = mm.fetch_pathway_graph("hsa00020")
check(g_hsa["status"] == "ok" and g_hsa["n_reactions"] == 2, f"human TCA graph fetched OK (got {g_hsa['status']}, {g_hsa['n_reactions']} reactions)")
check(g_hsa["nodes"]["C00022"]["label"] == "Pyruvate",
      f"a node known to KEGG's list/compound gets its REAL name, not the bare id (got {g_hsa['nodes']['C00022']['label']!r})")
check(g_hsa["nodes"]["G00001"]["label"] == "G00001",
      "a node absent from the (test) name index falls back to the bare id, not an error")
g_sce = mm.fetch_pathway_graph("sce00020")
check(g_sce["status"] == "ok" and g_sce["n_reactions"] == 0,
      "SENSITIVITY: yeast's SAME map number has a DIFFERENT reaction graph (0 reactions) - "
      "proving the fetch is keyed per organism, not just per map number")

check(len(_fetch_calls) == 2, f"two real fetches so far (got {len(_fetch_calls)})")
g_hsa_again = mm.fetch_pathway_graph("hsa00020")
check(len(_fetch_calls) == 2, "a second fetch of the SAME organism+map is served from cache, not re-fetched")
check(g_hsa_again["n_reactions"] == 2, "the cached graph round-trips correctly")

g_nodata = mm.fetch_pathway_graph("hsa99999")
check(g_nodata["status"] == "nodata" and g_nodata["nodes"] == {},
      "a map with no KGML diagram (KEGG's 404) is reported as 'nodata', not an error")
_calls_before = len(_fetch_calls)
mm.fetch_pathway_graph("hsa99999")
check(len(_fetch_calls) == _calls_before, "a 'nodata' result is cached too, so it is not re-fetched every run")

g_error = mm.fetch_pathway_graph("hsa88888")
check(g_error["status"] == "error", "a transient failure is reported as 'error'")
_calls_before2 = len(_fetch_calls)
mm.fetch_pathway_graph("hsa88888")
check(len(_fetch_calls) > _calls_before2, "SENSITIVITY: an 'error' result is NEVER cached - it is retried next time")


# ---------------------------------------------------------------------------
# run_pathway_topology end-to-end: p-value identity with run_ora, both input
# modes, and a genuine Impact computed through the real pipeline.
# ---------------------------------------------------------------------------
print("=== run_pathway_topology: identical p-values to run_ora's own KEGG block ===")


def _row(name, key, kegg_id):
    return {"Name": name, "compound_key": key, "match_tier": "kegg", "class_level_only": False,
            "chebi_id": "", "kegg_id": kegg_id, "inchikey": key, "pubchem_cid": "", "hmdb_id": "", "lm_id": "",
            "super_class": "", "main_class": "", "sub_class": "",
            "lm_category": "", "lm_main_class": "", "lm_sub_class": ""}


TCA_NAMES = ["Citrate", "Isocitrate", "2-Oxoglutarate", "Succinate", "Fumarate", "Malate"]
TCA_KEGG = ["C00022", "C00024", "C00036", "C00042", "C00122", "C00149"]
TCA_ROWS = [_row(TCA_NAMES[i], f"KEY_TCA{i}", TCA_KEGG[i]) for i in range(6)]
NOISE_ROWS = [_row(f"NOISE{i}", f"KEY_NOISE{i}", f"C9{i:04d}") for i in range(4)]
MATCHED = pd.DataFrame(TCA_ROWS + NOISE_ROWS, columns=mm.MATCHED_COLUMNS)

BY_KEGG = {row["kegg_id"]: row["compound_key"] for row in TCA_ROWS + NOISE_ROWS}
# 20 extra database-only compounds so the "database" universe scope differs
# from the "experiment" one, exactly as test_metabolomics_ora.py's fixture does.
for i in range(20):
    BY_KEGG[f"C8{i:04d}"] = f"KEY_EXTRA{i}"
CROSSWALK = {"by_chebi": {}, "by_kegg": BY_KEGG, "by_skeleton": {}}

FAKE_LINKS = "\n".join([f"cpd:{c}\tpath:map00020" for c in TCA_KEGG]
                       + [f"cpd:{c}\tpath:map00650" for c in [f"C9{i:04d}" for i in range(4)] + [f"C8{i:04d}" for i in range(20)]])
FAKE_NAMES = "path:map00020\tCitrate cycle (TCA cycle)\npath:map00650\tButanoate metabolism\n"
FAKE_ORG_LIST = "hsa00020\tCitrate cycle (TCA cycle) - Homo sapiens (human)\nhsa00650\tButanoate metabolism - Homo sapiens (human)\n"


# A dedicated graph carrying all 6 TCA compounds (unlike _KGML_BASIC above,
# built for the parser-level tests with only 3 of them plus a glycan) - a
# simple cycle, so every compound is reachable from every other.
_KGML_TCA_FULL = "\n".join(
    ['<?xml version="1.0"?>', '<pathway name="path:hsa00020" org="hsa" number="00020">']
    + [f'  <entry id="{i + 1}" name="cpd:{c}" type="compound"><graphics name="{c}" x="{i * 40}" y="{i * 20}"/></entry>'
       for i, c in enumerate(TCA_KEGG)]
    + [f'  <reaction id="{50 + i}" name="rn:R{i}" type="reversible">'
       f'<substrate id="{i + 1}" name="cpd:{TCA_KEGG[i]}"/>'
       f'<product id="{(i + 1) % 6 + 1}" name="cpd:{TCA_KEGG[(i + 1) % 6]}"/></reaction>'
       for i in range(6)]
    + ["</pathway>"]
)


def _fake_kegg_get_full(url, verbose=False, **kwargs):
    if "link/pathway/compound" in url:
        return "ok", FAKE_LINKS
    if url.endswith("/list/pathway/hsa"):
        return "ok", FAKE_ORG_LIST
    if url.endswith("/list/pathway"):
        return "ok", FAKE_NAMES
    if url.endswith("/get/hsa00020/kgml"):
        return "ok", _KGML_TCA_FULL
    if url.endswith("/get/hsa00650/kgml"):
        return "ok", _KGML_ZERO_REACTIONS
    return "error", ""


mm.resolve_kegg_organism_code = lambda *a, **k: {
    "organism_code": "hsa", "organism_name": "Homo sapiens", "verified": True, "note": None, "caveat": None, "candidates": [],
}
mm._kegg_safe_get = _fake_kegg_get_full
_FAKE_CACHE_STORE.clear()

HUMAN = {"taxon_id": "9606", "taxon_name": "Homo sapiens"}

topo_list = mm.run_pathway_topology(
    MATCHED, CROSSWALK, input_mode="list", measure="betweenness",
    min_term_size=3, max_term_size=500, **HUMAN,
)
block_topo = topo_list["tables"]["KEGG"]
check(block_topo["n_terms_tested"] == 2, f"both maps tested (got {block_topo['n_terms_tested']})")
check(block_topo.get("measure") == "betweenness", "block carries the measure used")
check("graphs" in block_topo and set(block_topo["graphs"]) == {"hsa00020", "hsa00650"},
      f"a graph is fetched for every TESTED pathway (got {list(block_topo.get('graphs', {}))})")

ora_result = mm.run_ora(
    MATCHED, CROSSWALK, databases=["KEGG"], min_term_size=3, max_term_size=500, **HUMAN,
)
ora_table = ora_result["tables"]["KEGG"]["table"]
topo_table = block_topo["table"]
merged = topo_table.merge(ora_table, on="Term ID", suffixes=("_topo", "_ora"))
check(len(merged) == 2, "both tested terms matched up between the two tables")
check((merged["p-value_topo"] == merged["p-value_ora"]).all() and (merged["FDR_topo"] == merged["FDR_ora"]).all(),
      "SENSITIVITY: p-values AND FDRs are byte-identical between run_pathway_topology and run_ora's own KEGG block")

tca_row = topo_table[topo_table["Term ID"] == "hsa00020"].iloc[0]
check(tca_row["Impact"] > 0, f"the TCA map (real reactions, all 6 compounds matched) has positive Impact (got {tca_row['Impact']})")
check(tca_row["Matched in graph"] == 6, f"all 6 TCA compounds are matched in the graph (got {tca_row['Matched in graph']})")
butanoate_row = topo_table[topo_table["Term ID"] == "hsa00650"].iloc[0]
check(butanoate_row["Impact"] == 0.0 and "no reactions" in str(butanoate_row["Topology note"]).lower(),
      f"the zero-reaction map's Impact is 0 AND flagged, not silently indistinguishable from a real zero "
      f"(got Impact={butanoate_row['Impact']}, note={butanoate_row['Topology note']!r})")
check("Overlapping Compounds" in topo_table.columns and "Citrate" in str(tca_row["Overlapping Compounds"]),
      "the member column is relabelled to the user's own compound Names, same as every other tile")


print("=== run_pathway_topology: concentration-table (Wilcoxon) mode ===")
rng = np.random.default_rng(7)
conc_names = TCA_NAMES + [f"NOISE{i}" for i in range(4)]
values = pd.DataFrame(rng.normal(10, 1, size=(12, len(conc_names))), columns=conc_names,
                      index=[f"S{i}" for i in range(12)])
values.iloc[6:, 0] += 3.0  # planted shift in Citrate
groups = pd.Series(["A"] * 6 + ["B"] * 6, index=values.index)
stats_df, _stats_meta = mm.compute_compound_statistics(values, groups, "A", "B", log_transform=False)

topo_conc = mm.run_pathway_topology(
    MATCHED, CROSSWALK, input_mode="conc", stats_df=stats_df, metric="log2fc", measure="outdegree",
    min_term_size=3, max_term_size=500, **HUMAN,
)
block_conc = topo_conc["tables"]["KEGG"]
check(block_conc.get("measure") == "outdegree", "concentration-mode block carries its own chosen measure")
check("Effect (mean difference)" in block_conc["table"].columns,
      "concentration mode uses the Wilcoxon/burden column shape (Effect, not Fold Enrichment)")

quant_result = mm.run_quantitative_enrichment(
    MATCHED, CROSSWALK, stats_df, databases=["KEGG"], metric="log2fc",
    min_term_size=3, max_term_size=500, **HUMAN,
)
quant_table = quant_result["tables"]["KEGG"]["table"]
merged_conc = block_conc["table"].merge(quant_table, on="Term ID", suffixes=("_topo", "_quant"))
check(len(merged_conc) > 0 and (merged_conc["p-value_topo"] == merged_conc["p-value_quant"]).all(),
      "SENSITIVITY: concentration-mode p-values are identical to run_quantitative_enrichment's own KEGG block too")


print("=== run_pathway_topology: no-genome organism -> empty with a note, never a fallback ===")
mm.resolve_kegg_organism_code = lambda *a, **k: {
    "organism_code": "", "organism_name": "", "verified": False,
    "note": 'KEGG has no genome for "Nonexistus fakeus".', "caveat": None, "candidates": [],
}
topo_none = mm.run_pathway_topology(
    MATCHED, CROSSWALK, input_mode="list", taxon_id="0", taxon_name="Nonexistus fakeus", min_term_size=3,
)
block_none = topo_none["tables"]["KEGG"]
check(block_none["n_terms_tested"] == 0 and block_none.get("note"),
      f"no KEGG genome -> zero terms tested with a note explaining why (got note={block_none.get('note')!r})")
check(block_none["table"].empty and "Impact" in block_none["table"].columns,
      "the empty table still carries the topology-specific columns (so the frontend's column set never varies)")


print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
