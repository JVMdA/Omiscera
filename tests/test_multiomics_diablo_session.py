# -*- coding: utf-8 -*-
"""Direct unit tests for the DIABLO / Multiblock Integration routes and jobs in API.py: no server, no network - imports API.py
and drives the route functions / jobs in-process, with the identifier resolvers and the pathway builders stubbed (the real ones
are exercised by live_e2e_multiomics_diablo.py).

THE THINGS UNDER TEST
  1. State: session["multiomics"]["diablo"] is its own slot; choosing omics for Joint keeps it; a change of the Correlation
     data drops the config, result, enrichment and cached pathway sets (the fingerprint), re-resolving identical data does not.
  2. Validation: nothing before the tables are matched, a group column REQUIRED (unlike MOFA), group sources (a file column,
     a sample sheet), every setting range, a cross-validation fold count larger than the smallest group refused.
  3. The run job on the REAL example files (planted factors recovered): component search, keepX search, the cross-validated
     error and its permutation p-value, strict JSON (no NaN - the known failure mode found live: AVE is NaN for a block with
     missing values and must become null, not break the response), selected features, downloads, isolation from the four
     single-omics sub-dicts and from Joint / Correlation.
  4. The enrichment job: validation, the real builders called with the right inputs, the 'unavailable' path, downloads.
"""
import asyncio
import json
import os
import sys

import numpy as np
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


class FakeUpload:
    def __init__(self, data, filename="x.csv"):
        self._d, self.filename = data, filename

    async def read(self):
        return self._d


def upload(omics, data, sid):
    return asyncio.run(API.multiomics_correlation_upload(omics, FakeUpload(data), "auto", x_session_id=sid))


def ex(name):
    return open(os.path.join(_ROOT, "Data3", name), "rb").read()


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


T_EX, P_EX, M_EX = ex("example_mofa_transcriptomics.csv"), ex("example_mofa_proteomics.csv"), ex("example_mofa_metabolomics.csv")
_t = pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_transcriptomics.csv"))
_p = pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_proteomics.csv"))
GENE_OF = {"transcriptomics": {g: g for g in _t["Gene"]}, "proteomics": dict(zip(_p["Protein IDs"], _t["Gene"]))}
ENTREZ = {g: str(1000 + i) for i, g in enumerate(_t["Gene"])}


def fake_resolve(omics_type, ids, taxon_id, progress_callback=None):
    mp = {i: GENE_OF[omics_type][i] for i in ids if i in GENE_OF[omics_type]}
    out = {"gene_of": mp, "n_features": len(ids), "n_mapped": len(mp), "n_unmapped": len(ids) - len(mp), "identifier_type": "stub"}
    if omics_type == "transcriptomics":
        out["entrez_of"] = {i: ENTREZ[mp[i]] for i in mp}
    else:
        out["acc_of"] = {i: i for i in mp}
    return out


API.multiomics_retriever.resolve_matrix_features = fake_resolve

print("=== state ===")
state = API.new_session_state()
check(state["multiomics"]["diablo"] == {"config": None, "sheet": None, "results": None, "enrich": {}, "metab_match": None, "inputs": {}},
      "a new session carries an empty DIABLO slot")
check(API._fresh_multiomics_diablo()["enrich"] is not API._fresh_multiomics_diablo()["enrich"], "each call returns independent dicts")

SID = "mo-diablo-test"
API.SESSIONS.pop(SID, None)
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "proteomics", "metabolomics"]}, x_session_id=SID)
S = API.get_session(SID)
CS, DS = S["multiomics"]["correlation"], S["multiomics"]["diablo"]

print("=== nothing before the tables are matched ===")
opt = API.multiomics_diablo_options(x_session_id=SID)
check(opt["ready"] is False and opt["samples"] == [] and opt["group_options"] == [], "options say 'not ready' and offer no groups")
check(refused(API.multiomics_diablo_configure, {}, x_session_id=SID)[0] == 400, "configure is refused before the tables are matched")
check(refused(API.multiomics_diablo_run, x_session_id=SID)[0] == 400, "run is refused before configuring")
check(refused(lambda: asyncio.run(API.multiomics_diablo_groups(FakeUpload(b"a,b\n1,2\n"), x_session_id=SID)))[0] == 400, "a sample sheet is refused before the tables are matched")
check(refused(API.multiomics_diablo_enrich, {"omics": "transcriptomics"}, x_session_id=SID)[0] == 400, "enrichment is refused before a result exists")
check(refused(API.multiomics_diablo_download, "cv", x_session_id=SID)[0] == 400, "download is refused before a result exists")

upload("transcriptomics", T_EX, SID)
upload("proteomics", P_EX, SID)
upload("metabolomics", M_EX, SID)
before = snapshot(S)
API._job_correlation_resolve("job-d1", S, "9606", "Homo sapiens", False)
check(len(CS["aligned"]["samples"]) == 36, "36 samples matched")
opt = API.multiomics_diablo_options(x_session_id=SID)
check(opt["ready"] and len(opt["samples"]) == 36 and [o["key"] for o in opt["group_options"]] == ["metabolomics:Group"]
      and opt["group_options"][0]["groups"] == ["Control", "Treated"] and opt["group_options"][0]["usable"],
      "the Group column of the metabolite file is offered, and it is usable (2+ groups of 3+ samples)")
check(opt["defaults"]["design_value"] == 0.1 and opt["min_samples_per_group"] == 3 and opt["max_components"] == 5,
      "the options carry the defaults and the sample / component limits")

print("=== configure: a group column is REQUIRED (unlike MOFA) ===")
code, detail = refused(API.multiomics_diablo_configure, {}, x_session_id=SID)
check(code == 400 and "group" in detail.lower(), f"SENSITIVITY: no group column at all is refused, with the reason ({detail[:80]})")
code, detail = refused(API.multiomics_diablo_configure, {"group_column": "nope:col"}, x_session_id=SID)
check(code == 400, "an unknown group column is refused")

print("=== configure: validation ===")
for bad, label in (({"group_column": "metabolomics:Group", "top_variance": 10}, "top variance 10"),
                   ({"group_column": "metabolomics:Group", "design_value": 1.5}, "design value 1.5"),
                   ({"group_column": "metabolomics:Group", "max_components": 0}, "max components 0"),
                   ({"group_column": "metabolomics:Group", "n_folds": 1}, "n_folds 1"),
                   ({"group_column": "metabolomics:Group", "n_folds": 20}, "n_folds too many"),
                   ({"group_column": "metabolomics:Group", "n_repeats": 0}, "n_repeats 0"),
                   ({"group_column": "metabolomics:Group", "distance": "nope"}, "unknown distance"),
                   ({"group_column": "metabolomics:Group", "keepX_grid": [1]}, "keepX grid value too small"),
                   ({"group_column": "metabolomics:Group", "keepX_grid": []}, "empty keepX grid"),
                   ({"group_column": "metabolomics:Group", "n_folds": 19}, "more folds than the smallest group (18) can support")):
    code, detail = refused(API.multiomics_diablo_configure, bad, x_session_id=SID)
    check(code == 400 and detail, f"SENSITIVITY: {label} is refused ({detail[:70]})")
check(DS["config"] is None, "...and a refused setting commits nothing")
out = API.multiomics_diablo_configure({"group_column": "metabolomics:Group", "max_components": 2, "keepX_grid": [5, 10, 20],
                                       "n_folds": 5, "n_repeats": 3, "top_variance": 150}, x_session_id=SID)
cfg = out["config"]
check(cfg["group_column"] == "metabolomics:Group" and cfg["max_components"] == 2 and cfg["keepX_grid"] == [5, 10, 20]
      and cfg["n_folds"] == 5 and cfg["n_repeats"] == 3 and cfg["design_value"] == 0.1 and cfg["distance"] == "centroids.dist"
      and cfg["seed"] == 42 and out["n_samples"] == 36 and out["groups"] == ["Control", "Treated"],
      "a valid configuration is stored with the defaults filled in")

print("=== run job on the real example files ===")
pay = API._job_diablo_run("job-d2", S)
json.dumps(pay, allow_nan=False)
check(True, "SENSITIVITY: the payload is strict-JSON safe (AVE's NaN for a block with missing values must serialise as null)")
check(any(v is None for v in pay["ave"]["proteomics"]), "...specifically, the proteomics AVE (missing values) came through as null, not dropped or crashing")
check(pay["n_components"] >= 1 and pay["n_components"] <= 2 and len(pay["component_search"]) == 2, f"{pay['n_components']} component(s) chosen out of 2 offered")
names = pay["component_names"]
check(len(names) == pay["n_components"], "one name per kept component")
check(pay["cv"]["BER"] < 0.1, f"SENSITIVITY: the cross-validated error on the real planted treatment signal is small ({pay['cv']['BER']:.3f})")
check(pay["permutation"]["p_value"] < 0.2, f"the permutation test flags it as unlikely by chance (p={pay['permutation']['p_value']:.3f})")
check(all(pay["n_selected"][o] > 0 for o in pay["omics"]), f"every omics selected at least one feature ({pay['n_selected']})")
check(set(pay["selected_features"]["transcriptomics"]) <= set(w["feature"] for w in pay["top_weights"]["transcriptomics"][names[0]]),
      "the selected features are among the reported top weights")
check(pay["groups"] and set(pay["groups"].values()) == {"Control", "Treated"}, "every sample's group is reported")
check(snapshot(S) == before, "ISOLATION: the four single-omics sub-dicts are untouched")
check(CS["results"] is None, "ISOLATION: Correlation's own result slot is untouched")

print("=== downloads ===")
for tb, head in (("cv", "Metric"), ("component_search", "Components"), (f"weights:{pay['omics'][0]}", "Feature"), (f"components:{pay['omics'][0]}", "Sample")):
    r = API.multiomics_diablo_download(tb, x_session_id=SID)
    first = r.body.decode().splitlines()[0]
    check(r.media_type == "text/csv" and first.split(",")[0] == head, f"download '{tb}' is a CSV starting with {head}")
check(refused(API.multiomics_diablo_download, "nonsense", x_session_id=SID)[0] == 400, "an unknown table is refused")
check(refused(API.multiomics_diablo_download, "enrichment:transcriptomics", x_session_id=SID)[0] == 400, "an enrichment that was not run is refused")

print("=== sample sheet overrides / adds a group source, and drops the result ===")
samples = pay["samples"]
sheet = "id,Batch\n" + "\n".join(f"{s},{'b1' if i % 3 == 0 else 'b2'}" for i, s in enumerate(samples)) + "\n"
check(DS["results"] is not None, "(a result exists before the sheet)")
res = asyncio.run(API.multiomics_diablo_groups(FakeUpload(sheet.encode(), "sheet.csv"), x_session_id=SID))
check(res["n_matched"] == 36 and res["columns"] == ["Batch"] and [o["key"] for o in res["group_options"]] == ["metabolomics:Group", "sheet:Batch"],
      f"the sheet's Batch column is offered alongside the file column ({res['columns']})")
check(DS["results"] is None and DS["enrich"] == {}, "a new sample sheet drops the old result")

print("=== the result is tied to the data it was fitted on ===")
API.multiomics_diablo_configure({"group_column": "metabolomics:Group", "max_components": 1, "n_folds": 3, "n_repeats": 2}, x_session_id=SID)
API._job_diablo_run("job-d3", S)
check(DS["results"] is not None and DS["config"] is not None, "(a result and a config exist)")
low = pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_proteomics.csv"))
low.iloc[0, 1] = low.iloc[0, 1] + 1.0
upload("proteomics", low.to_csv(index=False).encode(), SID)
API.multiomics_diablo_options(x_session_id=SID)
check(DS["results"] is None and DS["config"] is None and DS["sheet"] is None and DS["enrich"] == {} and DS["inputs"] == {},
      "SENSITIVITY: a re-uploaded table clears the DIABLO config, result, enrichment, cached pathway sets and sheet")
API._job_correlation_resolve("job-d4", S, "9606", "Homo sapiens", False)
check(refused(API.multiomics_diablo_run, x_session_id=SID)[0] == 400, "...and after re-matching, a run needs a new configuration")
upload("proteomics", P_EX, SID)
API._job_correlation_resolve("job-d5", S, "9606", "Homo sapiens", False)

print("=== Joint selection keeps the DIABLO state ===")
API.multiomics_diablo_configure({"group_column": "metabolomics:Group", "max_components": 1, "n_folds": 3, "n_repeats": 1}, x_session_id=SID)
API.multiomics_select({"omics_types": ["transcriptomics", "proteomics"]}, x_session_id=SID)
check(S["multiomics"]["diablo"]["config"] is not None and S["multiomics"]["diablo"] is DS, "choosing omics for Joint does not wipe DIABLO's state")

print("=== enrichment ===")
API.multiomics_diablo_configure({"group_column": "metabolomics:Group", "max_components": 1, "n_folds": 3, "n_repeats": 1, "keepX_grid": [5, 10]}, x_session_id=SID)
payE = API._job_diablo_run("job-d6", S)
for bad in ({"omics": "genomics"}, {"omics": "transcriptomics", "min_size": 1}, {"omics": "transcriptomics", "min_size": 50, "max_size": 20}):
    check(refused(API.multiomics_diablo_enrich, bad, x_session_id=SID)[0] == 400, f"SENSITIVITY: enrichment with {list(bad.values())[1:]} is refused")
check("job_id" in API.multiomics_diablo_enrich({"omics": "transcriptomics"}, x_session_id=SID), "a valid request starts a job")

wdf = DS["results"]["tables"]["weights:transcriptomics"]
selected = set(payE["selected_features"]["transcriptomics"])
label_of = dict(zip(wdf["Feature"].astype(str), wdf["Label"].astype(str)))
fake_terms = {"T_SEL": {"name": "sel", "genes": {label_of[f] for f in selected}},
             "T_OTHER": {"name": "other", "genes": set(label_of[f] for f in wdf["Feature"] if f not in selected) & set(list(label_of.values())[:40])}}
seen = {"calls": 0}


def fake_inputs(session, omics, job_id=None):
    seen["calls"] += 1
    feats = [str(c) for c in DS["results"]["blocks"][omics].columns]
    return [label_of.get(f) for f in feats], fake_terms, feats


_real_inputs = API._diablo_enrichment_inputs
API._diablo_enrichment_inputs = fake_inputs
outE = API._job_diablo_enrich("job-e1", S, "transcriptomics", 1, 500)
json.dumps(outE, allow_nan=False)
rowsE = outE["table"]["rows"] if "rows" in outE["table"] else outE["table"]["data"]
check(outE["n_tested"] == 2 and rowsE[0]["Term ID"] == "T_SEL" and rowsE[0]["FDR"] < 0.2, f"the selected-feature pathway is found first ({[r['Term ID'] for r in rowsE]})")
check(outE["n_selected"] == len(selected) and outE["n_selected_matched"] <= outE["n_selected"], "the payload counts the selected features and how many matched a key")
outE2 = API._job_diablo_enrich("job-e2", S, "transcriptomics", 1, 500)
check(seen["calls"] == 1, "SENSITIVITY: the pathway sets of an omics are built once per result (a second enrichment call reuses them)")
dl = API.multiomics_diablo_download("enrichment:transcriptomics", x_session_id=SID)
check(dl.body.decode().splitlines()[0].startswith("Term ID,Term Name"), "the enrichment is downloadable once it has been run")


def unavailable(session, omics, job_id=None):
    raise API._DiabloEnrichmentUnavailable("No pathway for these.")


API._diablo_enrichment_inputs = unavailable
DS["inputs"] = {}
outU = API._job_diablo_enrich("job-e3", S, "metabolomics", 10, 500)
check(outU["n_tested"] == 0 and outU["note"] == "No pathway for these.", "when no pathway can be built the job returns a note, not an error")
API._diablo_enrichment_inputs = _real_inputs

print("=== the real input builders are called with the right inputs ===")
cap = {}


def fake_build_t(universe_df, selected_arg, db, taxon_name, extra):
    cap["t"] = (universe_df.copy(), db, taxon_name)
    return {"term_map": {"X": {"name": "x", "genes": {"CS"}}}, "note": None}


_bt = API.multiomics_engine.build_transcriptomics_input
API.multiomics_engine.build_transcriptomics_input = fake_build_t
DS["inputs"] = {}
keysT, tmT, featsT = API._diablo_enrichment_inputs(S, "transcriptomics")
feats_t = [str(c) for c in DS["results"]["blocks"]["transcriptomics"].columns]
check(cap["t"][1] == "Reactome" and cap["t"][2] == "Homo sapiens" and set(cap["t"][0].columns) == {"gene_symbol", "entrez_id"} and len(cap["t"][0]) == len(feats_t),
      "transcripts: the universe is the model's features with their gene symbol and Entrez id, Reactome, the chosen species")
check(keysT == feats_t and featsT == feats_t, "transcript keys are the gene symbols, background is every feature that entered the model")
API.multiomics_engine.build_transcriptomics_input = _bt

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
