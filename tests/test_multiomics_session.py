# -*- coding: utf-8 -*-
"""Direct unit tests for the Multi-omics session state and jobs in API.py: no
live server, no network - imports API.py and drives the route functions and
jobs in-process with hand-built state and stubs (same technique as
test_metabolomics_session.py).

THE THINGS UNDER TEST
  1. ISOLATION (the user's explicit concern): Multi-omics has its own
     session["multiomics"] and never reads or writes session["ora"] /
     ["transcriptomics"] / ["proteomics"] / ["metabolomics"] - a full
     select -> resolve -> configure -> run pass leaves all four byte-identical.
  2. ONE KEGG organism per run: resolve_kegg_organism_code is called exactly
     ONCE however many omics take part, and the SAME code reaches every omics'
     KEGG builder (commit 27ce033 was two paths resolving it separately and
     disagreeing). A resolution that cannot produce a code (no genome /
     ambiguous) makes every omics' KEGG input an explanatory note, while
     Reactome is untouched.
  3. A failed resolve commits nothing (taxon, tables and results all stay as
     they were), the same commit-together discipline as _job_resolve_de.
  4. Route validation: <2 omics, an unselected omics, configure before resolve,
     a bad method, bad cutoffs are all refused with a message.
"""
import copy
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
    """Deep-comparable view of the four single-omics sub-dicts."""
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
print("=== state factory ===")
state = API.new_session_state()
check("multiomics" in state and set(state["multiomics"]["inputs"]) == {"genomics", "transcriptomics", "proteomics", "metabolomics"},
      "a new session carries a multiomics sub-dict with one input slot per omics")
check(state["multiomics"]["selected_omics"] == [] and state["multiomics"]["joint_results"] is None
      and state["multiomics"]["kegg_organism_confirmation"] is None,
      "...empty: nothing selected, no results, no KEGG confirmation")
check(API._fresh_multiomics_state() is not API._fresh_multiomics_state()
      and API._fresh_multiomics_state()["inputs"]["genomics"] is not API._fresh_multiomics_state()["inputs"]["genomics"],
      "each call returns independent dicts (no shared mutable default)")

# ---------------------------------------------------------------------------
print("=== select: 2-4 omics, keeps what is still selected ===")
SID = "mo-session-test"
API.SESSIONS.pop(SID, None)
code, detail = refused(API.multiomics_select, {"omics_types": ["genomics"]}, x_session_id=SID)
check(code == 400 and "two" in detail, f"one omics is refused ({code}: {detail})")
check(refused(API.multiomics_select, {"omics_types": []}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_select, {"omics_types": ["genomics", "bogus"]}, x_session_id=SID)[0] == 400,
      "none / an unknown name that leaves fewer than two are refused")
out = API.multiomics_select({"omics_types": ["metabolomics", "genomics", "genomics"]}, x_session_id=SID)
check(out["selected"] == ["genomics", "metabolomics"], f"selection is de-duplicated and kept in canonical order (got {out['selected']})")
S = API.get_session(SID)
S["multiomics"]["inputs"]["genomics"]["df"] = pd.DataFrame({"Gene Symbol": ["A"]})
S["multiomics"]["joint_config"] = {"databases": ["KEGG"]}
API.multiomics_select({"omics_types": ["genomics", "proteomics"]}, x_session_id=SID)
S = API.get_session(SID)
check(S["multiomics"]["inputs"]["genomics"]["df"] is not None and S["multiomics"]["joint_config"] is None
      and S["multiomics"]["selected_omics"] == ["genomics", "proteomics"],
      "re-selecting keeps an omics' upload but drops the config that described the old combination")

# ---------------------------------------------------------------------------
print("=== routes refuse out-of-order or invalid requests ===")
API.SESSIONS.pop(SID, None)
check(refused(API.multiomics_resolve, {"taxon_id": "9606"}, x_session_id=SID)[0] == 400,
      "resolve before choosing omics is refused")
API.multiomics_select({"omics_types": ["genomics", "metabolomics"]}, x_session_id=SID)
code, detail = refused(API.multiomics_resolve, {"taxon_id": "9606"}, x_session_id=SID)
check(code == 400 and "Genomics" in detail and "Metabolomics" in detail, f"resolve before uploading names what is missing ({detail})")
S = API.get_session(SID)
for o in ("genomics", "metabolomics"):
    S["multiomics"]["inputs"][o]["df"] = pd.DataFrame({"x": [1]})
check(refused(API.multiomics_resolve, {"taxon_id": ""}, x_session_id=SID)[0] == 400,
      "resolve without an organism is refused")
check(refused(API.multiomics_joint_configure, {}, x_session_id=SID)[0] == 400,
      "configure before every omics is resolved is refused")
check(refused(API.multiomics_joint_run, x_session_id=SID)[0] == 400, "run before configure is refused")
for o in ("genomics", "metabolomics"):
    S["multiomics"]["inputs"][o]["resolved_df"] = pd.DataFrame({"x": [1]})
check(refused(API.multiomics_joint_configure, {"databases": ["GO"]}, x_session_id=SID)[0] == 400,
      "a database outside KEGG / Reactome is refused")
check(refused(API.multiomics_joint_configure, {"method": "stouffer"}, x_session_id=SID)[0] == 400,
      "Stouffer is refused at the API too, not only in the engine")
check(refused(API.multiomics_joint_configure, {"min_term_size": "x"}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_joint_configure, {"min_term_size": 10, "max_term_size": 5}, x_session_id=SID)[0] == 400,
      "non-numeric / inverted term-size limits are refused")
ok = API.multiomics_joint_configure({"databases": ["KEGG", "Reactome"], "method": "tippett", "min_term_size": 5},
                                    x_session_id=SID)["config"]
check(ok["databases"] == ["KEGG", "Reactome"] and ok["method"] == "tippett" and ok["min_term_size"] == 5
      and ok["max_term_size"] == 500 and ok["de_cutoffs"] == {},
      f"a valid config is stored (500 default max, no DE cutoffs when no DE omics are selected): {ok}")
API.multiomics_select({"omics_types": ["genomics", "transcriptomics"]}, x_session_id=SID)
S = API.get_session(SID)
S["multiomics"]["inputs"]["transcriptomics"]["resolved_df"] = pd.DataFrame({"x": [1]})
check(refused(API.multiomics_joint_configure, {"de_cutoffs": {"transcriptomics": {"padj": 0, "lfc": 1}}}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_joint_configure, {"de_cutoffs": {"transcriptomics": {"padj": "x"}}}, x_session_id=SID)[0] == 400,
      "a padj cutoff of 0 or a non-number is refused")
cfg = API.multiomics_joint_configure({"de_cutoffs": {"transcriptomics": {"padj": 0.1, "lfc": -2}}}, x_session_id=SID)["config"]
check(cfg["de_cutoffs"] == {"transcriptomics": {"padj": 0.1, "lfc": 2.0}},
      f"cutoffs are stored per DE omics; |log2FC| is taken as a magnitude (got {cfg['de_cutoffs']})")

# ---------------------------------------------------------------------------
print("=== resolve job: stubbed resolvers, commit-together ===")
API.SESSIONS.pop(SID, None)
API.multiomics_select({"omics_types": ["genomics", "metabolomics"]}, x_session_id=SID)
S = API.get_session(SID)
S["multiomics"]["inputs"]["genomics"].update({"df": pd.DataFrame({"Gene Symbol": ["CS", "ACO2", "ZZZ"]}),
                                               "kind": "list", "identifier_type": "gene_symbol"})
S["multiomics"]["inputs"]["metabolomics"].update({"df": pd.DataFrame({"Name": ["Citrate", "Pyruvate"]}), "kind": "list"})
S["multiomics"]["taxon_id"], S["multiomics"]["taxon_name"] = "OLD", "Old organism"
S["multiomics"]["joint_results"] = {"stale": True}

_real_resolve_upload = API.multiomics_retriever.resolve_upload
_real_match = API._resolve_and_match_compound_names


def _fake_resolve_upload(omics_type, df, taxon_id, identifier_type=None, progress_callback=None):
    return {"resolved_df": pd.DataFrame({"Gene Symbol": df["Gene Symbol"], "UniProt ID": ["U1", "U2", ""],
                                         "KEGG ID": ["", "", ""], "Resolution_status": ["resolved", "resolved", "not_found"]}),
            "meta": {"identifier_type": identifier_type, "n_input_rows": len(df)}}


matched_ok = pd.DataFrame({"Name": ["Citrate", "Pyruvate"], "compound_key": ["K1", "K2"]})
API.multiomics_retriever.resolve_upload = _fake_resolve_upload
API._resolve_and_match_compound_names = lambda job_id, df, progress_start=0.0, progress_span=1.0, dedupe_names=True: (
    df.copy(), matched_ok, {"n_input_rows": len(df)})
before = snapshot(S)
res = API._job_resolve_multiomics("job-x", S, "9606", "Homo sapiens")
check(S["multiomics"]["taxon_id"] == "9606" and S["multiomics"]["taxon_name"] == "Homo sapiens"
      and S["multiomics"]["inputs"]["genomics"]["resolved_df"] is not None
      and S["multiomics"]["inputs"]["metabolomics"]["matched_df"] is not None,
      "a successful resolve commits the organism and both tables together")
check(S["multiomics"]["joint_results"] is None and S["multiomics"]["joint_config"] is None
      and S["multiomics"]["kegg_organism_confirmation"] is None,
      "...and drops results/config/KEGG confirmation that described the previous resolve")
check(res["omics"]["genomics"]["n_usable"] == 2 and res["omics"]["metabolomics"]["n_usable"] == 2,
      f"the summary counts usable rows per omics (genomics drops the unresolved gene): {res['omics']['genomics']['n_usable']}")
check(snapshot(S) == before, "ISOLATION: resolving left session['ora'/'transcriptomics'/'proteomics'/'metabolomics'] untouched")

# a failure in the SECOND omics must commit nothing
S["multiomics"]["taxon_id"], S["multiomics"]["taxon_name"] = "KEEP", "Keep organism"
old_gen = S["multiomics"]["inputs"]["genomics"]["resolved_df"]
S["multiomics"]["joint_results"] = {"keep": True}


def _boom(job_id, df, **kw):
    raise RuntimeError("PubChem unreachable")


API._resolve_and_match_compound_names = _boom
try:
    API._job_resolve_multiomics("job-x", S, "10090", "Mus musculus")
    check(False, "a resolver failure propagates")
except RuntimeError:
    check(True, "a resolver failure propagates (the job reports an error)")
check(S["multiomics"]["taxon_id"] == "KEEP" and S["multiomics"]["inputs"]["genomics"]["resolved_df"] is old_gen
      and S["multiomics"]["joint_results"] == {"keep": True},
      "SENSITIVITY: after a failure part-way, the organism, the earlier omics' table and the results are exactly as before")
API.multiomics_retriever.resolve_upload = _real_resolve_upload
API._resolve_and_match_compound_names = _real_match

# ---------------------------------------------------------------------------
print("=== run job: ONE KEGG organism per run, shared by every omics ===")
API.SESSIONS.pop(SID, None)
API.multiomics_select({"omics_types": ["genomics", "transcriptomics", "proteomics", "metabolomics"]}, x_session_id=SID)
S = API.get_session(SID)
MX = S["multiomics"]
MX["taxon_id"], MX["taxon_name"] = "9606", "Homo sapiens"
MX["inputs"]["genomics"].update({"kind": "list", "resolved_df": pd.DataFrame({
    "Gene Symbol": ["CS"], "UniProt ID": ["U1"], "KEGG ID": [""], "Resolution_status": ["resolved"]})})
MX["inputs"]["transcriptomics"].update({"kind": "de", "resolved_df": pd.DataFrame({
    "gene_symbol": ["CS"], "entrez_id": ["1"], "log2FC": [2.0], "padj": [0.01]})})
MX["inputs"]["proteomics"].update({"kind": "de", "resolved_df": pd.DataFrame({
    "accession": ["P1"], "gene_symbol": ["CS"], "log2FC": [2.0], "padj": [0.01]})})
MX["inputs"]["metabolomics"].update({"kind": "list", "matched_df": pd.DataFrame({
    "Name": ["Citrate"], "compound_key": ["K1"], "class_level_only": [False]})})
MX["joint_config"] = {"databases": ["KEGG", "Reactome"], "method": "fisher", "min_term_size": 1, "max_term_size": 500,
                      "de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0}, "proteomics": {"padj": 0.05, "lfc": 1.0}}}

seen = {"kegg_calls": [], "resolver_calls": 0, "ids": []}


def _fake_kegg_resolver(taxon_id, taxon_name, verbose=False, force_refresh=False, confirmed_organism_code=""):
    seen["resolver_calls"] += 1
    seen["confirmed"] = confirmed_organism_code
    return seen["resolver_result"]


def _term_input(omics):
    return {"term_map": {"p1": {"name": "Pathway one", "genes": {"m1", "m2", "m3"}}},
            "list_set": {"m1", "m2"}, "universe_set": {"m1", "m2", "m3", "m4"}, "labels": None,
            "note": None, "caveat": None}


mo = API.multiomics_engine
_saved = (mo.build_genomics_input, mo.build_transcriptomics_input, mo.build_proteomics_input,
          mo.build_metabolomics_input, API.transcriptomics_engine.resolve_kegg_organism_code,
          API.multiomics_retriever.assign_genomics_kegg_ids, API.metabolomics_engine.build_compound_crosswalk,
          API.multiomics_retriever.select_transcriptomics, API.multiomics_retriever.select_proteomics)
mo.build_genomics_input = lambda df, db, tid, tname, kegg, **k: (seen["kegg_calls"].append(("genomics", db, kegg)) or _term_input("g"))
mo.build_transcriptomics_input = lambda udf, sig, db, tname, kegg, **k: (seen["kegg_calls"].append(("transcriptomics", db, kegg)) or _term_input("t"))
mo.build_proteomics_input = lambda rdf, sig, db, tid, tname, kegg, **k: (seen["kegg_calls"].append(("proteomics", db, kegg)) or _term_input("p"))
mo.build_metabolomics_input = lambda mdf, cw, db, tid, tname, kegg, kname="", **k: (seen["kegg_calls"].append(("metabolomics", db, kegg)) or _term_input("m"))
API.transcriptomics_engine.resolve_kegg_organism_code = _fake_kegg_resolver
API.multiomics_retriever.assign_genomics_kegg_ids = lambda df, code: (seen["ids"].append(code) or df)
API.metabolomics_engine.build_compound_crosswalk = lambda verbose=False, **k: {}
API.multiomics_retriever.select_transcriptomics = lambda df, p, l: {"universe_df": df, "significant": {"CS"}, "n_tested": 1, "n_significant": 1}
API.multiomics_retriever.select_proteomics = lambda df, p, l: {"significant": {"P1"}, "n_tested": 1, "n_significant": 1}

seen["resolver_result"] = {"organism_code": "hsa", "organism_name": "Homo sapiens (human)", "note": None, "caveat": None}
before = snapshot(S)
out = API._job_run_multiomics_joint("job-y", S)
check(seen["resolver_calls"] == 1,
      f"SENSITIVITY: the KEGG organism is resolved exactly ONCE for a four-omics, two-database run (got {seen['resolver_calls']})")
kegg_codes = {(o, c) for o, d, c in seen["kegg_calls"] if d == "KEGG"}
check(kegg_codes == {("genomics", "hsa"), ("transcriptomics", "hsa"), ("proteomics", "hsa"), ("metabolomics", "hsa")},
      f"...and that one code reaches every omics' KEGG builder (got {sorted(kegg_codes)})")
check(seen["ids"] == ["hsa"], "...and Genomics' KEGG-gene-ID pass (a separate step) used it too")
check(set(out["blocks"]) == {"KEGG", "Reactome"} and out["summary"]["kegg_organism_code"] == "hsa"
      and out["summary"]["omics"] == ["genomics", "transcriptomics", "proteomics", "metabolomics"],
      "the result names the organism and the omics that took part")
check(out["blocks"]["KEGG"]["omics_combined"] == ["genomics", "transcriptomics", "proteomics", "metabolomics"]
      and out["blocks"]["KEGG"]["n_terms_combined"] == 1,
      "all four omics combine (stub pathway hit by every omics)")
check(all("Per-omics p-values are combined as if independent" in (b["caveat"] or "") for b in out["blocks"].values()),
      "every block carries the independence caveat")
check(snapshot(S) == before, "ISOLATION: a full run left the four single-omics session sub-dicts untouched")
check(MX["joint_results"]["blocks"]["KEGG"]["table"].shape[0] == 1, "the DataFrames stay in session state (only the payload is JSON)")

print("=== run job: Brown's method is the default, and its correlation reaches the payload ===")
check(API.multiomics_engine.COMBINE_METHOD_DEFAULT == "brown", "the engine's default method is Brown's")
MX["joint_config"]["method"] = "brown"
seen["kegg_calls"].clear(); seen["resolver_calls"] = 0
seen["resolver_result"] = {"organism_code": "hsa", "organism_name": "Homo sapiens (human)", "note": None, "caveat": None}
out_b = API._job_run_multiomics_joint("job-brown", S)
bk = out_b["blocks"]["KEGG"]
check(bk["method"] == "brown" and "Brown" in out_b["summary"]["method_label"],
      f"the run records Brown's method and labels it ({out_b['summary']['method_label']})")
check(bk["correlation"] is not None and bk["correlation"]["applied"] is True and len(bk["correlation"]["pairs"]) == 6,
      "the serialized block carries the correlation: six pairs for four omics, marked as applied")
check(bk["correlation"]["fallback"] is True and "fewer than 20" in (bk["note"] or ""),
      "SENSITIVITY: with a single jointly tested pathway Brown falls back to Fisher and the BLOCK NOTE says so (never silent)")
check("Brown's method" in bk["caveat"] and "as if independent" not in bk["caveat"],
      "SENSITIVITY: a Brown run carries Brown's caveat - it must NOT claim the omics were combined as if independent")
import json as _json  # noqa: E402
check(len(_json.dumps(out_b["blocks"]["KEGG"]["correlation"], allow_nan=False)) > 10,
      "the correlation block is strict-JSON safe (no NaN / inf)")
MX["joint_config"]["method"] = "fisher"
out_f = API._job_run_multiomics_joint("job-fisher-again", S)
check(out_f["blocks"]["KEGG"]["correlation"]["applied"] is False and "as if independent" in out_f["blocks"]["KEGG"]["caveat"],
      "a Fisher run still REPORTS the correlation but marks it not applied, and keeps the independence caveat")
MX["joint_config"]["method"] = "fisher"

print("=== run job: no usable KEGG organism -> KEGG explains itself, Reactome still runs ===")
seen["kegg_calls"].clear(); seen["resolver_calls"] = 0; seen["ids"].clear()
seen["resolver_result"] = {"organism_code": "", "organism_name": "", "caveat": None,
                           "note": "Several KEGG genomes match this organism; confirm one."}
out2 = API._job_run_multiomics_joint("job-z", S)
check(seen["resolver_calls"] == 1 and seen["ids"] == [], "still one resolution; no KEGG-gene-ID pass without a code")
check(not any(d == "KEGG" for _, d, _ in seen["kegg_calls"]),
      "SENSITIVITY: no omics' KEGG builder is called when there is no organism (no silent fallback to generic maps)")
kb = out2["blocks"]["KEGG"]
check(kb["n_terms_combined"] == 0 and "Several KEGG genomes" in (kb["note"] or ""),
      f"the KEGG block is empty and says why (note: {(kb['note'] or '')[:90]})")
check(out2["blocks"]["Reactome"]["n_terms_combined"] == 1, "Reactome is unaffected and still produces a joint result")
MX["kegg_organism_confirmation"] = {"taxon_id": "9606", "taxon_name": "Homo sapiens", "organism_code": "bsu",
                                    "organism_name": "x"}
API._job_run_multiomics_joint("job-c", S)
check(seen["confirmed"] == "bsu", "a confirmation the user made for THIS organism is forwarded to the one resolution")
MX["kegg_organism_confirmation"] = {"taxon_id": "OTHER", "taxon_name": "Homo sapiens", "organism_code": "bsu", "organism_name": "x"}
API._job_run_multiomics_joint("job-d", S)
check(seen["confirmed"] == "", "a confirmation made for a DIFFERENT organism is ignored (stale)")

(mo.build_genomics_input, mo.build_transcriptomics_input, mo.build_proteomics_input, mo.build_metabolomics_input,
 API.transcriptomics_engine.resolve_kegg_organism_code, API.multiomics_retriever.assign_genomics_kegg_ids,
 API.metabolomics_engine.build_compound_crosswalk, API.multiomics_retriever.select_transcriptomics,
 API.multiomics_retriever.select_proteomics) = _saved

print("=== download ===")
resp = API.multiomics_joint_download("KEGG", "combined", x_session_id=SID)
check(resp.media_type == "text/csv" and b"Combined FDR" in resp.body, "the combined table downloads as CSV")
check(refused(API.multiomics_joint_download, "GO", "combined", x_session_id=SID)[0] == 400,
      "a database that was not run is refused")
API.SESSIONS.pop(SID, None)

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
