# -*- coding: utf-8 -*-
"""Direct unit tests for the Cross-omics Pathway Mapping routes and jobs in API.py: no live
server, no network - imports API.py and drives the route functions / jobs in-process with
hand-built state and stubs (same technique as test_multiomics_session.py).

THE THINGS UNDER TEST
  1. State: a new session carries session["multiomics"]["mapping"]; configuring / re-uploading /
     re-resolving drops it (a mapping built on yesterday's uploads must never survive them).
  2. Validation: out-of-order requests, MetaCyc / unknown databases, malformed pathway ids (the id
     goes into a URL - '../x' must be refused) and bad cutoffs are refused with a message.
  3. The picker and the run use the SAME per-omics hit rule as Joint Pathway Enrichment (the
     shared _multiomics_prepare_omics / _multiomics_database_inputs), resolve the KEGG organism
     ONCE, skip Genomics' KEGG-id pass for Reactome, and never read the four single-omics tiles.
  4. A pathway whose structure cannot be fetched returns a result that says why and no map.
  5. Payloads are strict-JSON safe (no NaN / inf), and the CSV downloads carry the last result.
"""
import json
import os
import sys

import pandas as pd
from fastapi import HTTPException

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import API  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def refused(fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except HTTPException as exc:
        return exc.status_code, str(exc.detail)
    return None, ""


def snapshot(session):
    def norm(x):
        if isinstance(x, pd.DataFrame):
            return ("df", x.shape, tuple(x.columns), x.astype(str).values.tolist())
        if isinstance(x, dict):
            return {k: norm(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [norm(v) for v in x]
        if isinstance(x, set):
            return ("set", sorted(map(str, x)))
        return x
    return {k: norm(session[k]) for k in ("ora", "transcriptomics", "proteomics", "metabolomics")}


# ---------------------------------------------------------------------------
print("=== state ===")
state = API.new_session_state()
check(state["multiomics"]["mapping"] == {"config": None, "pickers": {}, "results": None},
      "a new session carries an empty mapping slot")
check(API._fresh_multiomics_mapping() is not API._fresh_multiomics_mapping()
      and API._fresh_multiomics_mapping()["pickers"] is not API._fresh_multiomics_mapping()["pickers"],
      "each call returns independent dicts")
mx0 = API._fresh_multiomics_state()
mx0["joint_config"], mx0["joint_results"] = {"a": 1}, {"b": 2}
mx0["mapping"] = {"config": {"x": 1}, "pickers": {"KEGG": {}}, "results": {"y": 1}}
API._multiomics_clear_results(mx0)
check(mx0["joint_config"] is None and mx0["joint_results"] is None and mx0["mapping"] == API._fresh_multiomics_mapping(),
      "SENSITIVITY: _multiomics_clear_results drops the joint AND the mapping state")

# ---------------------------------------------------------------------------
print("=== routes refuse out-of-order or invalid requests ===")
SID = "mo-map-session-test"
API.SESSIONS.pop(SID, None)
check(refused(API.multiomics_mapping_configure, {}, x_session_id=SID)[0] == 400, "configure before choosing omics is refused")
API.multiomics_select({"omics_types": ["genomics", "transcriptomics", "metabolomics"]}, x_session_id=SID)
check(refused(API.multiomics_mapping_configure, {}, x_session_id=SID)[0] == 400, "configure before resolving is refused")
S = API.get_session(SID)
MX = S["multiomics"]
MX["taxon_id"], MX["taxon_name"] = "9606", "Homo sapiens"
MX["inputs"]["genomics"].update({"kind": "list", "resolved_df": pd.DataFrame({
    "Gene Symbol": ["CS", "ACO2"], "UniProt ID": ["O75390", "Q99798"], "NCBI Gene ID": ["", ""], "KEGG ID": ["", ""],
    "Resolution_status": ["resolved", "resolved"]})})
MX["inputs"]["transcriptomics"].update({"kind": "de", "resolved_df": pd.DataFrame({
    "gene_symbol": ["CS", "ACO2", "IDH3A"], "entrez_id": ["1431", "50", "3419"], "log2FC": [2.0, -1.5, 0.1], "padj": [0.01, 0.02, 0.9]})})
MX["inputs"]["metabolomics"].update({"kind": "list", "resolved_df": pd.DataFrame({"Name": ["Citrate"]}), "matched_df": pd.DataFrame({
    "Name": ["Citrate"], "compound_key": ["K1"], "kegg_id": ["C00158"], "chebi_id": ["CHEBI:30769"], "class_level_only": [False]})})
check(refused(API.multiomics_mapping_pathways, {"database": "KEGG"}, x_session_id=SID)[0] == 400,
      "listing pathways before configuring is refused")
check(refused(API.multiomics_mapping_run, {"database": "KEGG", "pathway_id": "hsa00020"}, x_session_id=SID)[0] == 400,
      "mapping before configuring is refused")
check(refused(API.multiomics_mapping_configure, {"de_cutoffs": {"transcriptomics": {"padj": 0}}}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_mapping_configure, {"de_cutoffs": {"transcriptomics": {"padj": "x"}}}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_mapping_configure, {"de_cutoffs": {"transcriptomics": {"padj": 2}}}, x_session_id=SID)[0] == 400,
      "a padj cutoff of 0, above 1 or a non-number is refused")
cfg = API.multiomics_mapping_configure({"de_cutoffs": {"transcriptomics": {"padj": 0.1, "lfc": -2}}}, x_session_id=SID)["config"]
check(cfg == {"de_cutoffs": {"transcriptomics": {"padj": 0.1, "lfc": 2.0}}}
      and MX["mapping"]["config"] == cfg,
      f"cutoffs are stored for the DE omics that take part only, |log2FC| as a magnitude (got {cfg})")
MX["mapping"]["pickers"]["KEGG"] = {"stale": 1}
MX["mapping"]["results"] = {"stale": 1}
API.multiomics_mapping_configure({}, x_session_id=SID)
check(MX["mapping"]["pickers"] == {} and MX["mapping"]["results"] is None
      and MX["mapping"]["config"]["de_cutoffs"] == {"transcriptomics": {"padj": 0.05, "lfc": 1.0}},
      "re-configuring starts clean (no stale picker / result) and falls back to the default cutoffs")
check(refused(API.multiomics_joint_configure, {"de_cutoffs": {"transcriptomics": {"padj": 0}}}, x_session_id=SID)[0] == 400,
      "REGRESSION: Joint's configure still validates cutoffs (it now shares the helper)")

submitted = []
_real_submit = API.submit_job
API.submit_job = lambda fn, *a, **k: (submitted.append((fn.__name__, a)) or "job-1")
check(refused(API.multiomics_mapping_pathways, {"database": "MetaCyc"}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_mapping_pathways, {}, x_session_id=SID)[0] == 400,
      "MetaCyc / no database is refused (MetaCyc needs a paid BioCyc subscription)")
for bad in ("../../etc/passwd", "hsa00020/../x", "R-HSA-71403", "", "HSA00020", "hsa0002", "hsa00020;x"):
    check(refused(API.multiomics_mapping_run, {"database": "KEGG", "pathway_id": bad}, x_session_id=SID)[0] == 400,
          f"SENSITIVITY: KEGG pathway id {bad!r} is refused before it can reach a URL")
for bad in ("hsa00020", "R-HSA-", "R-HSA-71403/../..", "r-hsa-71403", "R-HS-71403"):
    check(refused(API.multiomics_mapping_run, {"database": "Reactome", "pathway_id": bad}, x_session_id=SID)[0] == 400,
          f"Reactome pathway id {bad!r} is refused")
check(submitted == [], "no job was started for any refused request")
out = API.multiomics_mapping_run({"database": "KEGG", "pathway_id": "hsa00020", "pathway_name": "TCA"}, x_session_id=SID)
out2 = API.multiomics_mapping_run({"database": "Reactome", "pathway_id": "R-MMU-71403"}, x_session_id=SID)
out3 = API.multiomics_mapping_pathways({"database": "Reactome"}, x_session_id=SID)
check(out == out2 == out3 == {"job_id": "job-1"} and [s[0] for s in submitted] == [
    "_job_run_multiomics_mapping", "_job_run_multiomics_mapping", "_job_multiomics_mapping_pathways"],
    "valid ids (human and mouse Reactome, KEGG organism code) start the right jobs")
API.submit_job = _real_submit

# ---------------------------------------------------------------------------
print("=== picker + run jobs: shared hit rule, one KEGG organism, strict JSON ===")
mo, mm, mr = API.multiomics_engine, API.multiomics_mapping, API.multiomics_retriever
seen = {"kegg_calls": [], "resolver_calls": 0, "ids": [], "fetch": []}


def _fake_resolver(taxon_id, taxon_name, verbose=False, force_refresh=False, confirmed_organism_code=""):
    seen["resolver_calls"] += 1
    return seen["resolver_result"]


def _term_input(prefix):
    return {"term_map": {"p1": {"name": "Pathway one", "genes": {prefix + "1", prefix + "2", prefix + "3"}},
                         "p2": {"name": "Pathway two", "genes": {prefix + "7"}}},
            "list_set": {prefix + "1", prefix + "2"}, "universe_set": {prefix + "1", prefix + "2", prefix + "3", prefix + "7"},
            "labels": None, "note": None, "caveat": None}


KGML_STRUCT = {"title": "Citrate cycle", "nodes": [
    {"id": "k1", "kind": "gene", "type": "gene", "label": "CS", "x": 10.0, "y": 10.0, "w": 40.0, "h": 17.0,
     "members": {"kegg_genes": ["hsa:1431"]}},
    {"id": "k2", "kind": "gene", "type": "gene", "label": "ACO2", "x": 100.0, "y": 10.0, "w": 40.0, "h": 17.0,
     "members": {"kegg_genes": ["hsa:50"]}},
    {"id": "k3", "kind": "compound", "type": "compound", "label": "C00158", "x": 10.0, "y": 80.0, "w": 8.0, "h": 8.0,
     "members": {"kegg_compounds": ["C00158"]}}],
    "reactions": [{"id": "rx0", "name": "R00351", "reversible": False, "enzymes": ["k1"], "substrates": [],
                   "products": ["k3"], "pos": None, "segments": []}]}
KGML_STRUCT["nodes"] += [
    # the same gene drawn a second time, as a LINE (a global-map style box); and a group (a drawn complex)
    {"id": "k4", "kind": "gene", "type": "gene", "label": "CS", "x": 50.0, "y": 50.0, "w": 0.0, "h": 0.0,
     "members": {"kegg_genes": ["hsa:1431"]}, "path": [0, 0, 50, 50, 100, 0]},
    {"id": "k5", "kind": "group", "type": "group", "label": "CS+ACO2", "x": 55.0, "y": 10.0, "w": 130.0, "h": 30.0,
     "members": {"kegg_genes": ["hsa:1431", "hsa:50"]}, "components": ["k1", "k2"]},
    {"id": "k6", "kind": "gene", "type": "gene", "label": "NOHIT", "x": 300.0, "y": 300.0, "w": 40.0, "h": 17.0,
     "members": {"kegg_genes": ["hsa:999"]}}]
KGML_STRUCT["relations"] = [{"a": "k1", "b": "k2", "rtype": "pprel", "types": ["activation"]},
                            {"a": "k2", "b": "k3", "rtype": "pcrel", "types": ["pcrel"]},
                            {"a": "k2", "b": "k6", "rtype": "pprel", "types": ["inhibition"]}]
KGML_STRUCT["decor_lines"] = [[0, 0, 5, 5]]
KGML_STRUCT["line_graphics"] = False
REACT_STRUCT = {"title": "Reactome path", "nodes": [
    {"id": "r1", "kind": "gene", "type": "Protein", "label": "CS", "x": 10.0, "y": 10.0, "w": 40.0, "h": 20.0,
     "members": {"uniprot": ["O75390"], "symbols": ["CS"]}}], "reactions": []}


def fake_fetch(database, pathway_id, force_refresh=False):
    seen["fetch"].append((database, pathway_id))
    if pathway_id == "R-HSA-60":
        return {"status": "overview", "note": "Reactome draws \"Disease\" only as links to its 2 sub-pathways.",
                "structure": {"title": "Disease", "nodes": [], "reactions": []},
                "subpathways": [{"stId": "R-HSA-222", "name": "Quiet sub-pathway"}, {"stId": "p1", "name": "Pathway one"}]}
    if pathway_id == "hsa99999":
        return {"status": "nodata", "structure": None, "note": "KEGG has no pathway map (KGML) for hsa99999."}
    if pathway_id == "R-HSA-8":
        return {"status": "error", "structure": None, "note": "Could not reach Reactome."}
    if database == "KEGG":
        return {"status": "ok", "structure": json.loads(json.dumps(KGML_STRUCT)), "note": None}
    return {"status": "ok", "structure": json.loads(json.dumps(REACT_STRUCT)), "note": None, "restricted_to": "Parent (R-HSA-6)"}


_saved = (mo.build_genomics_input, mo.build_transcriptomics_input, mo.build_proteomics_input, mo.build_metabolomics_input,
          API.transcriptomics_engine.resolve_kegg_organism_code, mr.assign_genomics_kegg_ids,
          API.metabolomics_engine.build_compound_crosswalk, mr.select_transcriptomics, mr.fetch_pathway_structure,
          API.genomics_ora.build_kegg_gene_map, API.metabolomics_engine.build_kegg_compound_name_index)
mo.build_genomics_input = lambda df, db, tid, tname, kegg, **k: (seen["kegg_calls"].append(("genomics", db, kegg)) or _term_input("g"))
mo.build_transcriptomics_input = lambda udf, sig, db, tname, kegg, **k: (seen["kegg_calls"].append(("transcriptomics", db, kegg)) or _term_input("t"))
mo.build_metabolomics_input = lambda mdf, cw, db, tid, tname, kegg, kname="", **k: (seen["kegg_calls"].append(("metabolomics", db, kegg)) or _term_input("m"))
API.transcriptomics_engine.resolve_kegg_organism_code = _fake_resolver
mr.assign_genomics_kegg_ids = lambda df, code: (seen["ids"].append(code) or df.assign(**{"KEGG ID": ["hsa:1431", "hsa:50"]}))
API.metabolomics_engine.build_compound_crosswalk = lambda verbose=False, **k: {}
mr.select_transcriptomics = lambda df, p, l: {"universe_df": df, "significant": {"CS", "ACO2"}, "n_tested": 3, "n_significant": 2}
mr.fetch_pathway_structure = fake_fetch
API.genomics_ora.build_kegg_gene_map = lambda org, **k: {"acc_to_gene": {"O75390": "hsa:1431", "Q99798": "hsa:50"},
                                                         "entrez_to_gene": {"1431": "hsa:1431", "50": "hsa:50"}, "symbol_to_gene": {}}
API.metabolomics_engine.build_kegg_compound_name_index = lambda **k: {"C00158": "Citrate"}

seen["resolver_result"] = {"organism_code": "hsa", "organism_name": "Homo sapiens (human)", "note": None, "caveat": None}
API.multiomics_mapping_configure({"de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0}}}, x_session_id=SID)
MX["joint_results"] = {"blocks": {"KEGG": {"table": pd.DataFrame({"Term ID": ["p1"], "Combined FDR": [0.003]})}}}
before = snapshot(S)
pick = API._job_multiomics_mapping_pathways("job-p", S, "KEGG")
rows = pick["table"]["rows"]
check([r["Term ID"] for r in rows] == ["p1"] and rows[0]["Omics with a hit"] == 3 and rows[0]["Total hits"] == 6,
      f"the picker lists only pathways some omics hit, ranked (got {[r['Term ID'] for r in rows]})")
check(rows[0]["Joint combined FDR"] == 0.003 and pick["has_joint_fdr"] is True,
      "a joint run's FDR for that database is shown beside the pathway")
check(seen["resolver_calls"] == 1 and seen["ids"] == ["hsa"] and {c for _, d, c in seen["kegg_calls"] if d == "KEGG"} == {"hsa"},
      "SENSITIVITY: one KEGG organism resolution; the one code reaches every KEGG builder and Genomics' id pass")
check(pick["table"]["columns"] == mm.picker_columns(["genomics", "transcriptomics", "metabolomics"]),
      "picker columns follow the selected omics")
seen["resolver_calls"] = 0; seen["ids"].clear(); seen["kegg_calls"].clear()
pick_r = API._job_multiomics_mapping_pathways("job-pr", S, "Reactome")
check(seen["ids"] == [], "SENSITIVITY: Genomics' KEGG-gene-ID pass is NOT run for Reactome (it would only cost network calls)")
check(pick_r["has_joint_fdr"] is False and rows and pick_r["database"] == "Reactome", "no joint FDR for a database that was not run")
json.dumps(pick, allow_nan=False)
json.dumps(pick_r, allow_nan=False)
check(True, "picker payloads are strict-JSON safe")

seen["resolver_calls"] = 0
res = API._job_run_multiomics_mapping("job-r", S, "KEGG", "hsa00020", "TCA given by the picker")
check(res["status"] == "ok" and res["pathway_name"] == "Citrate cycle" and res["omics"] == ["genomics", "transcriptomics", "metabolomics"],
      "the run maps the fetched structure for the selected omics, titled from the structure")
cs = next(n for n in res["nodes"] if n["id"] == "k1")
check(cs["status"]["genomics"]["state"] == "hit" and cs["status"]["transcriptomics"]["state"] == "hit"
      and cs["status"]["transcriptomics"]["direction"] == "up" and "proteomics" not in cs["status"],
      "the node carries a status per SELECTED omics only (Proteomics is not part of this run)")
cit = next(n for n in res["nodes"] if n["id"] == "k3")
check(cit["label"] == "Citrate" and cit["status"]["metabolomics"]["state"] == "hit",
      "a KEGG compound box is relabelled with KEGG's name and matched on its C-number")
check(res["summary"]["n_links"] == 1 and res["links"]["rows"][0]["Reaction"] == "R00351"
      and "CS" in res["links"]["rows"][0]["Enzyme (hit)"], "the CS -> citrate reaction link is found")
check(seen["resolver_calls"] == 1 and seen["fetch"][-1] == ("KEGG", "hsa00020"), "one KEGG organism resolution; the right pathway fetched")
json.dumps(res, allow_nan=False)
check(True, "the mapping payload is strict-JSON safe")
check(len(res["relations"]) == 3 and res["decor_lines"] == [[0, 0, 5, 5]] and res["line_graphics"] is False,
      "the result carries the map's relations and decoration lines")
byid = {n["id"]: n for n in res["nodes"]}
check(byid["k4"]["path"] == [0, 0, 50, 50, 100, 0] and byid["k5"]["components"] == ["k1", "k2"] and "path" not in byid["k1"],
      "SENSITIVITY: a line gene keeps its polyline and a group its components through the payload; a box carries neither")
check(byid["k6"]["status"] == {} and set(byid["k1"]["status"]) == {"genomics", "transcriptomics"} and "metabolomics" in byid["k3"]["status"],
      "SENSITIVITY: only the omics that touched a box carry a status ('absent' is left out - a global map has ~5,000 boxes)")
check(byid["k5"]["status"]["genomics"]["state"] == "hit", "a group is hit when a component is")
inter = res["interactions"]["rows"]
check([(r["From"], r["To"], r["Relation"]) for r in inter] == [("CS", "ACO2", "activation"), ("ACO2", "Citrate", "pcrel")]
      and res["summary"]["n_interactions"] == 2 and res["interactions"]["columns"] == mm.INTERACTION_COLUMNS,
      f"interaction links need both partners hit (ACO2 -> NOHIT is not one): {[(r['From'], r['To']) for r in inter]}")
check(res["summary"]["n_nodes"] == 5 and res["summary"]["per_omics"]["genomics"]["nodes_hit"] == 3,
      f"a group is not a box: 5 boxes, CS (twice) and ACO2 hit by Genomics (got {res['summary']['n_nodes']} / {res['summary']['per_omics']['genomics']['nodes_hit']})")
dli = API.multiomics_mapping_download("interactions", x_session_id=SID)
check(dli.media_type == "text/csv" and b"From hit by" in dli.body and b"pcrel" in dli.body, "the interaction links download as CSV too")
check(API.multiomics_mapping_download("anything-else", x_session_id=SID).body.startswith(b"Pathway node,Box ID"),
      "an unknown table name falls back to the member table, never to a path")
dl = API.multiomics_mapping_download("members", x_session_id=SID)
dl2 = API.multiomics_mapping_download("links", x_session_id=SID)
check(dl.media_type == "text/csv" and b"Matched by" in dl.body and b"R00351" in dl2.body,
      "the members and links CSVs carry the last mapped pathway")
check(snapshot(S) == before, "ISOLATION: picker and run left the four single-omics session sub-dicts untouched")

seen["ids"].clear()
res_r = API._job_run_multiomics_mapping("job-rr", S, "Reactome", "R-HSA-71403", "")
check(seen["ids"] == [], "SENSITIVITY: mapping a REACTOME pathway does not run Genomics' KEGG-gene-ID pass (it matches on UniProt) even though a KEGG organism is resolved for the bridge")
check(res_r["status"] == "ok" and res_r["restricted_to"] == "Parent (R-HSA-6)" and res_r["pathway_name"] == "Reactome path",
      "a Reactome pathway drawn from an ancestor reports which one")
seen["resolver_result"] = {"organism_code": "", "organism_name": "", "caveat": None, "note": "Several KEGG genomes match."}
res_n = API._job_run_multiomics_mapping("job-rn", S, "Reactome", "R-HSA-71403", "")
check(res_n["status"] == "ok" and "gene symbol" in (res_n["note"] or ""),
      "SENSITIVITY: with no KEGG organism, Reactome still maps and SAYS Transcriptomics fell back to gene-symbol matching")
try:
    API._job_run_multiomics_mapping("job-rk", S, "KEGG", "hsa00020", "")
    check(False, "KEGG with no organism raises")
except ValueError as exc:
    check("Several KEGG genomes" in str(exc), "KEGG with no resolvable organism fails with the organism note (no guessing)")
seen["resolver_result"] = {"organism_code": "hsa", "organism_name": "x", "note": None, "caveat": None}

res_x = API._job_run_multiomics_mapping("job-x", S, "KEGG", "hsa99999", "Nope")
check(res_x["status"] == "nodata" and res_x["nodes"] == [] and "hsa99999" in res_x["note"] and res_x["pathway_name"] == "Nope",
      "an unknown pathway returns a result that says why, with no map")
check(MX["mapping"]["results"] is None and refused(API.multiomics_mapping_download, "members", x_session_id=SID)[0] == 400,
      "...and the previous result no longer downloads as if it were this pathway")
res_e = API._job_run_multiomics_mapping("job-e", S, "Reactome", "R-HSA-8", "")
check(res_e["status"] == "error" and res_e["nodes"] == [] and res_e["note"], "a transport error is reported, not hidden")

# ---------------------------------------------------------------------------
print("=== a global map: line colours are named from KEGG's classification ===")
LINE_STRUCT = {"title": "Metabolic pathways", "nodes": [
    {"id": "k1", "kind": "gene", "type": "gene", "label": "CS", "x": 50.0, "y": 50.0, "w": 0.0, "h": 0.0,
     "members": {"kegg_genes": ["hsa:1431"]}, "path": [0, 50, 100, 50], "color": "#B3B3E6"},
    {"id": "k2", "kind": "gene", "type": "gene", "label": "ACO2", "x": 50.0, "y": 80.0, "w": 0.0, "h": 0.0,
     "members": {"kegg_genes": ["hsa:50"]}, "path": [0, 80, 100, 80], "color": "#80CCCC"},
    {"id": "k3", "kind": "gene", "type": "gene", "label": "NOHIT", "x": 50.0, "y": 110.0, "w": 0.0, "h": 0.0,
     "members": {"kegg_genes": ["hsa:999"]}, "path": [0, 110, 100, 110], "color": "#B3B3E6"},
    {"id": "k4", "kind": "link", "type": "map", "label": "Glycolysis", "x": 200.0, "y": 50.0, "w": 40.0, "h": 20.0,
     "members": {}, "pathway": "00010"},
    {"id": "k5", "kind": "link", "type": "map", "label": "Ribosome", "x": 200.0, "y": 90.0, "w": 40.0, "h": 20.0,
     "members": {}, "pathway": "03010"}],
    "reactions": [], "relations": [], "decor_lines": [], "line_graphics": True}
_fetch_prev = mr.fetch_pathway_structure
_classes_prev, _univ_prev = mr.fetch_pathway_classes, API.genomics_ora.build_kegg_universe
mr.fetch_pathway_structure = lambda database, pid, force_refresh=False: (
    {"status": "ok", "structure": json.loads(json.dumps(LINE_STRUCT)), "note": None} if pid == "hsa01100" else _fetch_prev(database, pid))
mr.fetch_pathway_classes = lambda **k: {"00010": "Carbohydrate metabolism", "00061": "Lipid metabolism", "01100": "Global and overview maps"}
API.genomics_ora.build_kegg_universe = lambda org, **k: {"pathways": {
    "hsa00010": {"genes": {"hsa:1431", "hsa:999"}}, "hsa00061": {"genes": {"hsa:50"}},
    "hsa01100": {"genes": {"hsa:1431", "hsa:50", "hsa:999"}}}}
gl = API._job_run_multiomics_mapping("job-glob", S, "KEGG", "hsa01100", "")
byn = {n["id"]: n for n in gl["nodes"]}
check(gl["line_graphics"] is True and [c["name"] for c in gl["classes"]] == ["Carbohydrate metabolism", "Lipid metabolism"],
      f"a line map's result lists its pathway classes ({[c['name'] for c in gl['classes']]})")
cb = {c["name"]: c for c in gl["classes"]}
check(cb["Carbohydrate metabolism"] == {"name": "Carbohydrate metabolism", "color": "#B3B3E6", "n_lines": 2, "n_hit": 1}
      and cb["Lipid metabolism"]["n_hit"] == 1,
      f"...with KEGG's colour, the number of lines and the number hit ({cb['Carbohydrate metabolism']})")
check(byn["k1"]["cls"] == "Carbohydrate metabolism" and byn["k2"]["cls"] == "Lipid metabolism" and byn["k3"]["cls"] == "Carbohydrate metabolism",
      "each line gene carries its class")
check(byn["k4"]["cls"] == "Carbohydrate metabolism" and "cls" not in byn["k5"],
      "a pathway label carries its class where KEGG classifies it as metabolism (Ribosome: none)")
json.dumps(gl, allow_nan=False)
check(True, "the global-map payload is strict-JSON safe")
box = API._job_run_multiomics_mapping("job-box", S, "KEGG", "hsa00020", "")
check(box["classes"] == [] and all("cls" not in n for n in box["nodes"]), "a map made of boxes has no class colouring")
rr = API._job_run_multiomics_mapping("job-box2", S, "Reactome", "R-HSA-71403", "")
check(rr["classes"] == [], "neither has a Reactome diagram")
mr.fetch_pathway_structure = _fetch_prev
mr.fetch_pathway_classes, API.genomics_ora.build_kegg_universe = _classes_prev, _univ_prev

print("=== an overview pathway lists its sub-pathways with the picker's hit counts ===")
res_ov = API._job_run_multiomics_mapping("job-ov", S, "Reactome", "R-HSA-60", "Disease")
sp = res_ov["subpathways"]
check(res_ov["status"] == "overview" and res_ov["nodes"] == [] and res_ov["pathway_name"] == "Disease" and "2 sub-pathways" in res_ov["note"],
      "an overview result has no map, says why, and carries the pathway's name")
check([r["Term ID"] for r in sp["rows"]] == ["p1", "R-HSA-222"] and sp["rows"][0]["Total hits"] == 6 and sp["rows"][1]["Total hits"] == 0
      and sp["columns"] == mm.picker_columns(["genomics", "transcriptomics", "metabolomics"]),
      f"SENSITIVITY: the sub-pathway the picker lists carries its hit counts and ranks first; the quiet one has 0 ({[r['Term ID'] for r in sp['rows']]})")
check(MX["mapping"]["results"] is None, "an overview is not a result to download")
json.dumps(res_ov, allow_nan=False)
check(True, "the overview payload is strict-JSON safe")
MX["mapping"]["pickers"].clear()
seen["resolver_calls"] = 0
res_ov2 = API._job_run_multiomics_mapping("job-ov2", S, "Reactome", "R-HSA-60", "Disease")
check([r["Term ID"] for r in res_ov2["subpathways"]["rows"]] == ["p1", "R-HSA-222"] and "rows" in MX["mapping"]["pickers"]["Reactome"],
      "SENSITIVITY: reaching an overview from a Joint result (no picker built yet) builds the counts itself and keeps them")

print("=== re-resolving drops the mapping ===")
API.multiomics_mapping_configure({}, x_session_id=SID)
MX["mapping"]["pickers"]["KEGG"] = {"x": 1}
_rr, _rm = mr.resolve_upload, API._resolve_and_match_compound_names
mr.resolve_upload = lambda omics_type, df, taxon_id, identifier_type=None, progress_callback=None: {
    "resolved_df": MX["inputs"][omics_type]["resolved_df"], "meta": {"n_input_rows": len(df)}}
API._resolve_and_match_compound_names = lambda job_id, df, progress_start=0.0, progress_span=1.0, dedupe_names=True: (
    df.copy(), MX["inputs"]["metabolomics"]["matched_df"], {"n_input_rows": len(df)})
for o in ("genomics", "transcriptomics", "metabolomics"):
    MX["inputs"][o]["df"] = pd.DataFrame({"x": [1]})
API._job_resolve_multiomics("job-rs", S, "9606", "Homo sapiens")
check(MX["mapping"] == API._fresh_multiomics_mapping(), "SENSITIVITY: a successful re-resolve drops the mapping config, pickers and result")
mr.resolve_upload, API._resolve_and_match_compound_names = _rr, _rm

(mo.build_genomics_input, mo.build_transcriptomics_input, mo.build_proteomics_input, mo.build_metabolomics_input,
 API.transcriptomics_engine.resolve_kegg_organism_code, mr.assign_genomics_kegg_ids,
 API.metabolomics_engine.build_compound_crosswalk, mr.select_transcriptomics, mr.fetch_pathway_structure,
 API.genomics_ora.build_kegg_gene_map, API.metabolomics_engine.build_kegg_compound_name_index) = _saved
API.SESSIONS.pop(SID, None)

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
