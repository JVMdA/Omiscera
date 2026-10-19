# -*- coding: utf-8 -*-
"""Direct unit tests for the Network-based Integration routes and job in API.py: no live server, no network - the
route functions and the job are driven in-process with hand-built state, a synthetic STRING graph and stubs.

THE THINGS UNDER TEST
  1. State: a new session carries session["multiomics"]["network"]; re-uploading / re-resolving / re-configuring drops it.
  2. Validation: out-of-order requests, one omics, an organism STRING does not have, and every parameter's range are refused
     with a message; nothing half-valid is stored.
  3. The job: seeds from every omics (equal weight), both networks in one result, a strict-JSON payload, the per-omics
     tables, and the four single-omics sub-dicts untouched.
  4. Downloads carry the last result; an unknown table or a missing result is a 400.
"""
import gzip
import json
import os
import sys
import tempfile

import numpy as np
import pandas as pd
from fastapi import HTTPException

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import API  # noqa: E402
import Genomics_module as gm  # noqa: E402

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
check(state["multiomics"]["network"] == {"config": None, "results": None}, "a new session carries an empty network slot")
check(API._fresh_multiomics_network() is not API._fresh_multiomics_network(), "each call returns an independent dict")
mx0 = API._fresh_multiomics_state()
mx0["network"] = {"config": {"x": 1}, "results": {"y": 1}}
API._multiomics_clear_results(mx0)
check(mx0["network"] == API._fresh_multiomics_network(), "SENSITIVITY: _multiomics_clear_results drops the network state too")

# ---------------------------------------------------------------------------
print("=== routes refuse out-of-order or invalid requests ===")
SID = "mo-net-session-test"
API.SESSIONS.pop(SID, None)
check(refused(API.multiomics_network_configure, {}, x_session_id=SID)[0] == 400, "configure before choosing omics is refused")
check(refused(API.multiomics_network_run, x_session_id=SID)[0] == 400, "run before configuring is refused")
API.multiomics_select({"omics_types": ["genomics", "transcriptomics", "proteomics", "metabolomics"]}, x_session_id=SID)
check(refused(API.multiomics_network_configure, {}, x_session_id=SID)[0] == 400, "configure before resolving is refused")

# ----- a synthetic STRING-like graph of 60 proteins in 3 communities -----
TAX = "9995"
N = 60
names = [f"Q{i}" for i in range(N)]
PID = {n: f"{TAX}.{n}" for n in names}
rng = np.random.default_rng(5)
edges = []
for c in range(3):
    b0 = c * 20
    for i in range(20):
        for j in range(i + 1, 20):
            if rng.random() < 0.45:
                edges.append((names[b0 + i], names[b0 + j], int(rng.integers(750, 999))))
    for i in range(19):
        edges.append((names[b0 + i], names[b0 + i + 1], 900))
links = []
for a, b, w in edges:
    links += [(PID[a], PID[b], w), (PID[b], PID[a], w)]


def _gz(path, header, rows, sep):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(str(x) for x in r) + "\n")


TMP = tempfile.TemporaryDirectory()
L, I, A = (os.path.join(TMP.name, f) for f in ("l.gz", "i.gz", "a.gz"))
_gz(L, ["protein1", "protein2", "combined_score"], links, " ")
_gz(I, ["#string_protein_id", "preferred_name", "protein_size", "annotation"], [(PID[n], f"GENE{i}", "100", "x") for i, n in enumerate(names)], "\t")
_gz(A, ["#string_protein_id", "alias", "source"], [(PID[n], f"ACC{i}", "UniProt_AC") for i, n in enumerate(names)], "\t")
gnet = API.genomics_net
_real = (gnet._string_file_paths, gnet.is_taxon_in_string, gm.build_go_universe, API._multiomics_prepare_omics, API._multiomics_database_inputs,
         API.multiomics_retriever.fetch_reactome_reaction_maps, API.multiomics_retriever.fetch_chebi_skeletons_cached)
gnet._string_file_paths = lambda taxon_id: {"links": L, "info": I, "aliases": A}
gm.build_go_universe = lambda *a, **k: {"biological_process": {"GO:1": {"name": "community A", "genes": {f"ACC{i}" for i in range(20)}}},
                                        "molecular_function": {}, "cellular_component": {}}


def _rows_p(sig):
    return pd.DataFrame([{"input_id": f"r{i}", "accession": f"ACC{i}", "gene_symbol": f"GENE{i}", "protein_name": "",
                          "log2FC": 3.0 if i in sig else 0.05, "pvalue": 1e-8 if i in sig else 0.8,
                          "padj": 1e-6 if i in sig else 0.9, "status": "resolved"} for i in range(N)])


def _rows_t(sig):
    return pd.DataFrame([{"input_id": f"t{i}", "gene_symbol": f"GENE{i}", "ensembl_id": "", "entrez_id": "",
                          "log2FC": -2.5 if i in sig else 0.1, "pvalue": 1e-7 if i in sig else 0.7,
                          "padj": 1e-5 if i in sig else 0.8, "status": "resolved"} for i in range(N)])


S = API.get_session(SID)
MX = S["multiomics"]
MX["taxon_id"], MX["taxon_name"] = TAX, "Test organism"
MX["inputs"]["genomics"].update({"kind": "list", "resolved_df": pd.DataFrame({
    "Gene Symbol": [f"GENE{i}" for i in (45, 46, 47, 55)], "UniProt ID": [f"ACC{i}" for i in (45, 46, 47, 55)],
    "Resolution_status": ["resolved", "resolved", "resolved", "ambiguous"]})})
MX["inputs"]["transcriptomics"].update({"kind": "de", "resolved_df": _rows_t(set(range(3, 9)))})
MX["inputs"]["proteomics"].update({"kind": "de", "resolved_df": _rows_p(set(range(0, 6)))})
MX["inputs"]["metabolomics"].update({"kind": "list", "resolved_df": pd.DataFrame({"Name": ["metabA", "metabB"]}), "matched_df": pd.DataFrame({
    "Name": ["metabA", "metabB"], "compound_key": ["K1", "K2"], "kegg_id": ["", ""], "chebi_id": ["10", "11"], "class_level_only": [False, False]})})

gnet.is_taxon_in_string = lambda taxon_id: False
code, detail = refused(API.multiomics_network_configure, {}, x_session_id=SID)
check(code == 400 and "STRING has no protein-interaction network" in detail and "Test organism" in detail, f"an organism STRING does not have is refused by name ({detail[:70]})")
gnet.is_taxon_in_string = lambda taxon_id: True
MX["taxon_id"] = ""
check(refused(API.multiomics_network_configure, {}, x_session_id=SID)[0] == 400, "no organism is refused")
MX["taxon_id"] = TAX
for bad in ({"confidence": 123}, {"confidence": "x"}, {"n_permutations": 200}, {"n_permutations": 10 ** 9}, {"fdr_alpha": 0}, {"fdr_alpha": 2},
            {"module_database": "MSigDB"}, {"min_module_size": 1}, {"max_hits_for_modules": 10}, {"max_reactions": 2}, {"max_reactions": 5000}, {"max_reaction_proteins": 1}, {"max_reaction_proteins": 9999},
            {"pathway_top_n": 1}, {"pathway_top_n": 500}, {"pathway_min_jaccard": 0}, {"pathway_min_jaccard": 1.5}, {"pathway_merge_jaccard": 0.1}, {"pathway_merge_jaccard": 1.5},
            {"de_cutoffs": {"transcriptomics": {"padj": 0}}}):
    check(refused(API.multiomics_network_configure, bad, x_session_id=SID)[0] == 400, f"refused: {bad}")
check(MX["network"]["config"] is None, "nothing half-valid was stored")
cfg = API.multiomics_network_configure({"n_permutations": 500, "confidence": 700, "min_module_size": 3, "fdr_alpha": 0.5, "module_database": "GO",
                                        "pathway_top_n": 20, "pathway_min_jaccard": 0.2, "max_reactions": 40, "max_reaction_proteins": 25}, x_session_id=SID)["config"]
check(cfg["n_permutations"] == 500 and cfg["max_reactions"] == 40 and cfg["max_reaction_proteins"] == 25 and cfg["restart_prob"] == gnet.DEFAULT_RESTART_PROB and cfg["seed"] == 42
      and set(cfg["de_cutoffs"]) == {"transcriptomics", "proteomics"} and MX["network"]["config"] == cfg, "a valid config is stored (DE cutoffs for the two DE omics only)")
MX["network"]["results"] = {"stale": 1}
API.multiomics_network_configure({"fdr_alpha": 0.5, "n_permutations": 500}, x_session_id=SID)
check(MX["network"]["results"] is None, "re-configuring drops the previous result")
API.multiomics_network_configure({"n_permutations": 500, "fdr_alpha": 0.5, "min_module_size": 3, "pathway_top_n": 20, "pathway_min_jaccard": 0.2, "max_reactions": 40}, x_session_id=SID)

submitted = []
_real_submit = API.submit_job
API.submit_job = lambda fn, *a, **k: (submitted.append(fn.__name__) or "job-n")
check(API.multiomics_network_run(x_session_id=SID) == {"job_id": "job-n"} and submitted == ["_job_multiomics_network"], "run starts the network job")
API.submit_job = _real_submit

# ---------------------------------------------------------------------------
print("=== the job ===")
API.multiomics_retriever.fetch_reactome_reaction_maps = lambda taxon_id, taxon_name, force_refresh=False: {
    "status": "ok", "note": None, "chebi_reactions": {"10": ["RA"], "11": ["RB"]}, "reaction_proteins": {"RA": ["ACC40", "ACC41"], "RB": ["ACC42"]},
    "reaction_names": {}}
API.multiomics_retriever.fetch_chebi_skeletons_cached = lambda ids: {}
TERMS = {"R-T1": {"name": "Community A pathway", "genes": {f"GENE{i}" for i in range(0, 12)}},
         "R-T2": {"name": "Community A overlap", "genes": {f"GENE{i}" for i in range(4, 14)}},
         "R-T3": {"name": "Community C pathway", "genes": {f"GENE{i}" for i in range(44, 50)}}}
PTERMS = {k: {"name": v["name"], "genes": {f"ACC{g[4:]}" for g in v["genes"]}} for k, v in TERMS.items()}
UNIV = {f"GENE{i}" for i in range(N)}


def _fake_prepare(mx, selected, cutoffs, kegg_code, job_id=None, assign_kegg_ids=True):
    return {o: {} for o in selected}, {}


def _fake_inputs(mx, selected, prep, database, kegg_code, kegg_name, kegg_note, kegg_caveat):
    out = {}
    if "genomics" in selected:
        out["genomics"] = {"term_map": PTERMS, "list_set": {f"ACC{i}" for i in (45, 46, 47)}, "universe_set": {f"ACC{i}" for i in range(N)},
                           "labels": {f"ACC{i}": f"GENE{i}" for i in range(N)}}
    if "transcriptomics" in selected:
        out["transcriptomics"] = {"term_map": TERMS, "list_set": {f"GENE{i}" for i in range(3, 9)}, "universe_set": UNIV, "labels": None}
    if "proteomics" in selected:
        out["proteomics"] = {"term_map": PTERMS, "list_set": {f"ACC{i}" for i in range(0, 6)}, "universe_set": {f"ACC{i}" for i in range(N)},
                             "labels": {f"ACC{i}": f"GENE{i}" for i in range(N)}}
    if "metabolomics" in selected:
        out["metabolomics"] = {"term_map": {"R-T1": {"name": "Community A pathway", "genes": {"K1"}}}, "list_set": {"K1"}, "universe_set": {"K1", "K2", "K3"},
                               "labels": {"K1": "metabA"}}
    return out


API._multiomics_prepare_omics, API._multiomics_database_inputs = _fake_prepare, _fake_inputs
before = snapshot(S)
progress = []
_real_update = API.update_job
API.update_job = lambda jid, **kw: progress.append(kw.get("progress"))
try:
    out = API._job_multiomics_network("job-run", S)
finally:
    API.update_job = _real_update
json.dumps(out, allow_nan=False)
check(True, "the whole payload is strict-JSON safe")
check(out["omics_used"] == ["genomics", "transcriptomics", "proteomics", "metabolomics"], f"all four omics took part ({out['omics_used']})")
_ps = [p for p in progress if p is not None]
check(len(_ps) >= 5 and _ps[-1] == 1.0 and all(0.0 <= p <= 1.0 for p in _ps) and 0.85 in _ps, f"progress is reported, stays within 0-1 and ends at 1.0 ({len(_ps)} updates)")
prop = out["propagation"]
tbl = {r["Gene Symbol"]: r for r in prop["table"]["rows"]}
check(tbl["GENE4"]["Seed omics"] == "Transcriptomics + Proteomics" and tbl["GENE46"]["Seed omics"] == "Genomics" and tbl["GENE41"]["Seed omics"] == "Metabolomics",
      "each seed names the omics that seeded it (a metabolite's reaction seeds ACC41)")
check(abs(sum(r["Seed score"] for r in prop["table"]["rows"] if r["Seed score"] is not None) - out["seeds"]["n_seeds"]) < 1e-6, "the seed weights sum to the number of seeds")
check(tbl["GENE55"]["Seed"] is False and tbl["GENE47"]["Seed"] is True and out["omics"][0]["In the seed list"] == 3,
      "SENSITIVITY: an AMBIGUOUS Genomics row is never a seed (only the 3 resolved genes are)")
check(sum(1 for r in out["omics"] if r["Used"]) == 4 and abs(sum(r["Share of the combined weight"] for r in out["omics"]) - 1) < 1e-12, "four omics, equal shares")
check(out["metabolites"]["columns"][:2] == ["Compound", "Matched ChEBI"] and len(out["metabolites"]["rows"]) == 2, "the metabolite -> enzyme table is in the payload")
check(prop["edges"] and all({"source", "target", "weight"} <= set(e) for e in prop["edges"]), "the sub-network edges are shipped for the existing figures")
check(out["modules"] is not None and "Seeds from Genomics" in out["modules"]["table"]["columns"], "modules carry per-omics seed counts")
pw = out["pathway_network"]
check([n["id"] for n in pw["nodes"]] and pw["omics"] and all("x" in n and "y" in n for n in pw["nodes"]), f"the pathway network has nodes with positions ({[n['id'] for n in pw['nodes']]})")
check(pw["canvas"]["width"] >= 1100 and all(0 < n["x"] < pw["canvas"]["width"] and 0 < n["y"] < pw["canvas"]["height"] for n in pw["nodes"]),
      "the engine sizes the canvas and every node is inside it (positions are pixels)")
check(any({e["a"], e["b"]} == {"R-T1", "R-T2"} for e in pw["edges"]), "the two overlapping Community-A pathways are joined")
check(snapshot(S) == before, "ISOLATION: the four single-omics sub-dicts are untouched")
check(S["multiomics"]["joint_results"] is None and S["multiomics"]["mapping"] == API._fresh_multiomics_mapping() and S["multiomics"]["correlation"]["results"] is None,
      "...and Joint, Mapping and Correlation")
res = MX["network"]["results"]
check(res is not None and set(res["tables"]) == {"hits", "modules", "metabolites", "omics", "pathway_nodes", "pathway_edges"}, "the DataFrames are kept in the session")

print("=== downloads ===")
for which, header in [("hits", b"STRING ID,Gene Symbol"), ("modules", b"Module ID,Size"), ("metabolites", b"Compound,Matched ChEBI"), ("omics", b"Omics,Entities uploaded"),
                      ("pathway_nodes", b"Term ID,Pathway"), ("pathway_edges", b"Pathway A,Pathway B")]:
    dl = API.multiomics_network_download(which, x_session_id=SID)
    check(dl.media_type == "text/csv" and dl.body.startswith(header) and "multiomics_network_" in dl.headers.get("content-disposition", ""), f"'{which}' downloads as CSV with its columns and a multiomics_network_ file name")
check("multiomics_network_hits.csv" in API.multiomics_network_download("hits", x_session_id=SID).headers["content-disposition"], "the hit table's file name is multiomics_network_hits.csv")
mid = out["modules"]["table"]["rows"][0]["Module ID"]
e = out["modules"]["enrichment"].get(mid)
if e:
    dle = API.multiomics_network_download("enrichment", key=mid, x_session_id=SID)
    check(dle.media_type == "text/csv" and len(dle.body) > 20, "a module's enrichment downloads")
check(refused(API.multiomics_network_download, "nope", x_session_id=SID)[0] == 400, "an unknown table is refused")
check(refused(API.multiomics_network_download, "enrichment", "M99", x_session_id=SID)[0] == 400, "an unknown module is refused")
check(refused(API.multiomics_network_download, "hits", x_session_id="mo-net-none")[0] == 400, "no result is a 400")

print("=== a run with too few omics, and clearing ===")
MX["selected_omics"] = ["transcriptomics", "proteomics"]
API.multiomics_network_configure({"n_permutations": 500, "fdr_alpha": 0.5}, x_session_id=SID)
API._multiomics_clear_results(MX)
check(MX["network"] == API._fresh_multiomics_network(), "a re-upload / re-resolve drops the network result")
MX["selected_omics"] = ["transcriptomics"]
check(refused(API.multiomics_network_configure, {}, x_session_id=SID)[0] == 400, "ONE selected omics cannot be configured")

(gnet._string_file_paths, gnet.is_taxon_in_string, gm.build_go_universe, API._multiomics_prepare_omics, API._multiomics_database_inputs,
 API.multiomics_retriever.fetch_reactome_reaction_maps, API.multiomics_retriever.fetch_chebi_skeletons_cached) = _real
API.SESSIONS.pop(SID, None)
print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
