# -*- coding: utf-8 -*-
"""Direct unit tests for the MOFA / Factor-based Integration routes and jobs in API.py: no server, no network - imports API.py and
drives the route functions / jobs in-process, with the identifier resolvers and the pathway builders stubbed (the real ones are exercised
by live_e2e_multiomics_mofa.py).

THE THINGS UNDER TEST
  1. State: session["multiomics"]["mofa"] is its own slot; choosing omics for Joint keeps it; a change of the Correlation data drops
     the config, result, enrichment and sample sheet (the fingerprint), re-resolving identical data does not.
  2. Validation: nothing before the tables are matched, every setting range, an unknown group column, a bad sample sheet.
  3. The run job on the REAL MOFA example files: planted factors recovered, payload strict-JSON safe, group test, 25 top weights,
     downloads, sample-sheet groups, nothing leaking into the single-omics sub-dicts / Joint / Correlation.
  4. The enrichment job: validation, collapse of features sharing a gene, the planted pathway found, the "unavailable" path,
     downloads - and the real builders are called with the right inputs.
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
TRUTH = pd.read_csv(os.path.join(_ROOT, "tests", "fixtures", "mofa_example_truth.csv")).set_index("Sample")
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
check(state["multiomics"]["mofa"] == {"config": None, "results": None, "sheet": None, "enrich": {}, "metab_match": None, "inputs": {}},
      "a new session carries an empty MOFA slot")
check(API._fresh_multiomics_mofa()["enrich"] is not API._fresh_multiomics_mofa()["enrich"], "each call returns independent dicts")

SID = "mo-mofa-test"
API.SESSIONS.pop(SID, None)
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "proteomics", "metabolomics"]}, x_session_id=SID)
S = API.get_session(SID)
CS, MS = S["multiomics"]["correlation"], S["multiomics"]["mofa"]

print("=== nothing before the tables are matched ===")
opt = API.multiomics_mofa_options(x_session_id=SID)
check(opt["ready"] is False and opt["samples"] == [] and opt["group_options"] == [], "options say 'not ready' and offer no groups")
check(refused(API.multiomics_mofa_configure, {}, x_session_id=SID)[0] == 400, "configure is refused before the tables are matched")
check(refused(API.multiomics_mofa_run, x_session_id=SID)[0] == 400, "run is refused before configuring")
check(refused(lambda: asyncio.run(API.multiomics_mofa_groups(FakeUpload(b"a,b\n1,2\n"), x_session_id=SID)))[0] == 400, "a sample sheet is refused before the tables are matched")
check(refused(API.multiomics_mofa_enrich, {"omics": "transcriptomics", "factor": "Factor1"}, x_session_id=SID)[0] == 400, "enrichment is refused before a result exists")
check(refused(API.multiomics_mofa_download, "factors", x_session_id=SID)[0] == 400, "download is refused before a result exists")

upload("transcriptomics", T_EX, SID)
upload("proteomics", P_EX, SID)
upload("metabolomics", M_EX, SID)
before = snapshot(S)
API._job_correlation_resolve("job-m1", S, "9606", "Homo sapiens", False)
check(len(CS["aligned"]["samples"]) == 36 and CS["resolved"]["transcriptomics"]["entrez_of"]["CS"] == ENTREZ["CS"]
      if "CS" in ENTREZ else len(CS["aligned"]["samples"]) == 36, "36 samples matched; the resolved Entrez ids reach the Correlation state")
check(CS["resolved"]["proteomics"]["acc_of"], "protein accessions are kept next to the gene symbols (the enrichment needs them)")
opt = API.multiomics_mofa_options(x_session_id=SID)
check(opt["ready"] and len(opt["samples"]) == 36 and [o["key"] for o in opt["group_options"]] == ["metabolomics:Group"]
      and opt["group_options"][0]["groups"] == ["Control", "Treated"] and opt["group_options"][0]["n_samples"] == 36,
      "the Group column of the metabolite file is offered as a sample-group source")
check(opt["defaults"]["factors"] == 10 and opt["min_samples"] == 6 and opt["warn_samples"] == 15, "the options carry the defaults and the sample limits")

print("=== configure: validation ===")
for bad, label in (({"factors": 1}, "factors 1"), ({"factors": 26}, "factors 26"), ({"factors": "x"}, "factors text"),
                   ({"convergence": "slowest"}, "convergence"), ({"top_variance": 10}, "top variance 10"), ({"top_variance": 9000}, "top variance 9000"),
                   ({"n_refits": 9}, "refits 9"), ({"n_refits": -1}, "refits -1"), ({"min_r2": 50}, "min R2 50"),
                   ({"stability_cut": 2}, "stability cut 2"), ({"group_column": "nope:col"}, "unknown group column")):
    code, detail = refused(API.multiomics_mofa_configure, bad, x_session_id=SID)
    check(code == 400 and detail, f"SENSITIVITY: {label} is refused ({detail[:60]})")
check(MS["config"] is None, "...and a refused setting commits nothing")
out = API.multiomics_mofa_configure({"factors": 10, "n_refits": 2, "group_column": "metabolomics:Group", "log2": False}, x_session_id=SID)
cfg = out["config"]
check(cfg["factors"] == 10 and cfg["n_refits"] == 2 and cfg["group_column"] == "metabolomics:Group" and cfg["convergence"] == "fast"
      and cfg["scale_views"] is True and cfg["top_variance"] == 2000 and cfg["seed"] == 42 and out["n_samples"] == 36,
      "a valid configuration is stored with the defaults filled in")

print("=== run job on the real example files ===")
pay = API._job_mofa_run("job-m2", S)
json.dumps(pay, allow_nan=False)
check(True, "the payload is strict-JSON safe")
names = pay["factor_values"]["names"]
Z = pd.DataFrame(pay["factor_values"]["values"], index=pay["factor_values"]["samples"], columns=names).loc[TRUTH.index]
rec = {t: max(abs(np.corrcoef(Z[f], TRUTH[t])[0, 1]) for f in names) for t in ("treatment", "cholesterol", "lipids")}
check(all(v > 0.95 for v in rec.values()), f"SENSITIVITY: the three planted factors are recovered through the API ({ {k: round(v, 3) for k, v in rec.items()} })")
check(pay["converged"] and pay["n_factors"] == len(names) <= 10 and pay["n_factors_requested"] == 10, f"{pay['n_factors']} of 10 factors kept, converged")
r2 = np.array(pay["r2"]["values"])
check(r2.shape == (3, len(names)) and pay["r2"]["omics"] == ["Transcriptomics", "Proteomics", "Metabolomics"] and (r2 >= 0).all() and r2.sum(axis=1).max() < 100,
      "variance explained is omics x factors, in percent")
check(all(len(pay["top_weights"][o][f]) == 25 for o in pay["omics"] for f in names), "25 top-weighted features per factor and omics")
tw = pay["top_weights"]["transcriptomics"][names[0]]
check(all(abs(tw[i]["weight"]) >= abs(tw[i + 1]["weight"]) for i in range(24)) and tw[0]["label"] == tw[0]["feature"], "they are sorted by |weight|, labelled with the gene symbol")
_wt = MS["results"]["tables"]["weights:transcriptomics"]
_exp = _wt.reindex(_wt[names[0]].abs().sort_values(ascending=False).index)["Feature"].tolist()[:25]
check([r["feature"] for r in tw] == _exp and any(r["weight"] < 0 for r in tw) and any(r["weight"] > 0 for r in tw),
      "SENSITIVITY: the 25 shown are the largest by |weight| (both signs), exactly as the full table ranks them")
check(pay["group_column"] == "metabolomics:Group" and pay["groups"]["Control_1"] == "Control" and pay["groups"]["Treated_18"] == "Treated", "the groups of the chosen column are returned per sample")
gt = pay["group_tests"]
check(gt["test"] == "Mann-Whitney U" and min(r["FDR"] for r in gt["rows"]) < 1e-4 and sorted(r["FDR"] for r in gt["rows"])[1] > 0.05,
      "exactly one factor separates the groups (the treatment factor)")
check(pay["stability"]["n_refits"] == 2 and not any(p["unstable"] for p in pay["stability"]["per_factor"][:3]), "stability was checked with 2 refits")
check(snapshot(S) == before, "ISOLATION: the four single-omics sub-dicts are untouched")
check(CS["results"] is None, "ISOLATION: Correlation's own result slot is untouched")

print("=== downloads ===")
for tb, head in (("factors", "Sample"), ("r2", "Omics"), ("stability", "Factor"), ("groups", "Factor"), ("weights:transcriptomics", "Feature")):
    r = API.multiomics_mofa_download(tb, x_session_id=SID)
    first = r.body.decode().splitlines()[0]
    check(r.media_type == "text/csv" and first.split(",")[0] == head, f"download '{tb}' is a CSV starting with {head}")
rw = pd.read_csv(__import__("io").StringIO(API.multiomics_mofa_download("weights:proteomics", x_session_id=SID).body.decode()))
check(len(rw) == pay["n_features"]["proteomics"] and "Label" in rw.columns, "the weights table covers every feature that entered the model")
check(refused(API.multiomics_mofa_download, "nonsense", x_session_id=SID)[0] == 400, "an unknown table is refused")
check(refused(API.multiomics_mofa_download, "enrichment:transcriptomics:Factor1:signed", x_session_id=SID)[0] == 400, "an enrichment that was not run is refused")

print("=== sample sheet ===")
bad_ids = b"id,Batch\nzz1,a\nzz2,b\nzz3,a\n"
code, detail = refused(lambda: asyncio.run(API.multiomics_mofa_groups(FakeUpload(bad_ids), x_session_id=SID)))
check(code == 400 and "Only 0" in detail, f"SENSITIVITY: a sheet whose ids match no sample is refused ({detail[:60]})")
samples = pay["samples"]
unique_only = "id,Name\n" + "\n".join(f"{s},n{i}" for i, s in enumerate(samples)) + "\n"
code, detail = refused(lambda: asyncio.run(API.multiomics_mofa_groups(FakeUpload(unique_only.encode()), x_session_id=SID)))
check(code == 400 and "group label" in detail, "a sheet with only unique-valued columns offers no group")
sheet = "id,Batch,Name\n" + "\n".join(f"{s},{'b1' if i % 3 == 0 else 'b2'},n{i}" for i, s in enumerate(samples)) + "\n"
check(MS["results"] is not None, "(a result exists before the sheet)")
res = asyncio.run(API.multiomics_mofa_groups(FakeUpload(sheet.encode(), "sheet.csv"), x_session_id=SID))
check(res["n_matched"] == 36 and res["columns"] == ["Batch"] and [o["key"] for o in res["group_options"]] == ["metabolomics:Group", "sheet:Batch"],
      f"the sheet's Batch column is offered; 'Name' (all different) is not ({res['columns']})")
check(MS["results"] is None and MS["enrich"] == {}, "a new sample sheet drops the old result (it was coloured by other groups)")
API.multiomics_mofa_configure({"factors": 10, "n_refits": 0, "group_column": "sheet:Batch"}, x_session_id=SID)
pay2 = API._job_mofa_run("job-m3", S)
check(pay2["groups"]["Control_1"] == "b1" and pay2["groups"]["Control_2"] == "b2" and pay2["group_tests"]["rows"], "the run is coloured and tested by the sheet's column")
check(pay2["stability"]["n_refits"] == 0, "n_refits = 0 skips the stability refits")
API.multiomics_mofa_configure({"factors": 10, "n_refits": 0, "group_column": ""}, x_session_id=SID)
pay3 = API._job_mofa_run("job-m4", S)
check(pay3["groups"] is None and pay3["group_tests"]["rows"] == [], "with no group column nothing is tested")
check(np.allclose(pay2["factor_values"]["values"], pay3["factor_values"]["values"]), "SENSITIVITY: groups are for interpretation only - the fitted factors are identical with and without them")

print("=== enrichment ===")
API.multiomics_mofa_configure({"factors": 10, "n_refits": 0, "group_column": "metabolomics:Group"}, x_session_id=SID)
payE = API._job_mofa_run("job-m5", S)
fnames = payE["factor_values"]["names"]
for bad in ({"omics": "genomics", "factor": fnames[0]}, {"omics": "transcriptomics", "factor": "Factor99"}, {"omics": "transcriptomics", "factor": fnames[0], "mode": "x"},
            {"omics": "transcriptomics", "factor": fnames[0], "min_size": 1}, {"omics": "transcriptomics", "factor": fnames[0], "min_size": 50, "max_size": 20}):
    check(refused(API.multiomics_mofa_enrich, bad, x_session_id=SID)[0] == 400, f"SENSITIVITY: enrichment with {list(bad.values())[1:]} is refused")
check("job_id" in API.multiomics_mofa_enrich({"omics": "transcriptomics", "factor": fnames[0]}, x_session_id=SID), "a valid request starts a job")

# which factor is the treatment factor, and which genes belong to it (the file's treatment genes are the ones with a large weight on it)
f_treat = max(fnames, key=lambda f: abs(np.corrcoef(pd.DataFrame(payE["factor_values"]["values"], index=payE["factor_values"]["samples"], columns=fnames).loc[TRUTH.index, f],
                                                    TRUTH["treatment"])[0, 1]))
wdf = MS["results"]["tables"]["weights:transcriptomics"]
top_genes = wdf.reindex(wdf[f_treat].abs().sort_values(ascending=False).index)["Label"].tolist()
planted = set(top_genes[:40])
other = [g for g in wdf["Label"] if g not in planted]
# the database keys are NOT the labels (a real run has InChIKeys / accessions): the job must show the labels in "Top members"
KP = "K_"
fake_terms = {"T_planted": {"name": "planted", "genes": {KP + g for g in planted}}, "T_other": {"name": "other", "genes": {KP + g for g in other[:40]}},
              "T_other2": {"name": "other2", "genes": {KP + g for g in other[40:80]}}}
seen = {"calls": 0}


def fake_inputs(session, omics, job_id=None):
    seen["omics"] = omics
    seen["calls"] += 1
    return [KP + str(l) for l in wdf["Label"]], fake_terms


_real_inputs = API._mofa_enrichment_inputs
API._mofa_enrichment_inputs = fake_inputs
outE = API._job_mofa_enrich("job-e1", S, "transcriptomics", f_treat, "absolute", 10, 500)
json.dumps(outE, allow_nan=False)
rowsE = outE["table"]["rows"] if "rows" in outE["table"] else outE["table"]["data"]
check(outE["n_tested"] == 3 and rowsE[0]["Term ID"] == "T_planted" and rowsE[0]["FDR"] < 1e-6, f"the planted pathway is found first on the treatment factor ({f_treat}), absolute mode")
check(outE["factor"] == f_treat and outE["mode"] == "absolute" and outE["n_features_matched"] == len(wdf), "the payload names the factor / mode and counts the matched features")
tm = rowsE[0]["Top members"]
check(KP not in tm and all(g in planted for g in tm.split(", ")), f"SENSITIVITY: 'Top members' shows the features' own labels, not the database keys ({tm})")
outE2 = API._job_mofa_enrich("job-e1b", S, "transcriptomics", fnames[0], "signed", 10, 500)
check(seen["calls"] == 1 and "transcriptomics" in MS["inputs"], "SENSITIVITY: the pathway sets of an omics are built once per result (a second factor / test reuses them)")
dl = API.multiomics_mofa_download(f"enrichment:transcriptomics:{f_treat}:absolute", x_session_id=SID)
check(dl.body.decode().splitlines()[0].startswith("Term ID,Term Name"), "the enrichment is downloadable once it has been run")
# features sharing a gene count once: duplicate every key, the pathway test must not change
MS["inputs"].clear()
dup_keys = [KP + str(l) for l in wdf["Label"]]
dup_keys[1] = dup_keys[0]
API._mofa_enrichment_inputs = lambda session, omics, job_id=None: (dup_keys, fake_terms)
outD = API._job_mofa_enrich("job-e2", S, "transcriptomics", f_treat, "absolute", 10, 500)
rowsD = outD["table"]["rows"] if "rows" in outD["table"] else outD["table"]["data"]
check(outD["n_features_matched"] == len(wdf) - 1, "SENSITIVITY: two features of one gene count once in the pathway test")


def unavailable(session, omics, job_id=None):
    raise API._MofaEnrichmentUnavailable("No pathway for these.")


API._mofa_enrichment_inputs = unavailable
outU = API._job_mofa_enrich("job-e3", S, "metabolomics", fnames[0], "signed", 10, 500)
check(outU["n_tested"] == 0 and outU["note"] == "No pathway for these.", "when no pathway can be built the job returns a note, not an error")
API._mofa_enrichment_inputs = _real_inputs

print("=== the real input builders are called with the right inputs ===")
cap = {}


def fake_build_t(universe_df, selected, db, taxon_name, extra):
    cap["t"] = (universe_df.copy(), selected, db, taxon_name)
    return {"term_map": {"X": {"name": "x", "genes": {"CS"}}}, "note": None}


def fake_build_p(df, selected, db, taxon_id, taxon_name, extra):
    cap["p"] = (df.copy(), selected, db, taxon_id, taxon_name)
    return {"term_map": {"X": {"name": "x", "genes": {"P1"}}}, "note": None}


_bt, _bp = API.multiomics_engine.build_transcriptomics_input, API.multiomics_engine.build_proteomics_input
API.multiomics_engine.build_transcriptomics_input, API.multiomics_engine.build_proteomics_input = fake_build_t, fake_build_p
keysT, tmT = API._mofa_enrichment_inputs(S, "transcriptomics")
feats_t = [str(c) for c in MS["results"]["views"]["transcriptomics"].columns]
udf = cap["t"][0]
check(cap["t"][2] == "Reactome" and cap["t"][3] == "Homo sapiens" and set(udf.columns) == {"gene_symbol", "entrez_id"} and len(udf) == len(feats_t)
      and dict(zip(udf["gene_symbol"], udf["entrez_id"])) == {g: ENTREZ[g] for g in feats_t}, "transcripts: the universe is the model's features with their gene symbol and Entrez id, Reactome, the chosen species")
check(keysT == feats_t, "transcript keys are the gene symbols")
keysP, tmP = API._mofa_enrichment_inputs(S, "proteomics")
feats_p = [str(c) for c in MS["results"]["views"]["proteomics"].columns]
check(cap["p"][2] == "Reactome" and cap["p"][3] == "9606" and set(cap["p"][0].columns) == {"accession", "gene_symbol"} and len(cap["p"][0]) == len(feats_p)
      and keysP == feats_p, "proteins: accessions + gene symbols, Reactome, taxon 9606; keys are the accessions")
CS["resolved"]["transcriptomics"]["entrez_of"] = {}
code = None
try:
    API._mofa_enrichment_inputs(S, "transcriptomics")
except API._MofaEnrichmentUnavailable as exc:
    code = str(exc)
check(code and "Entrez" in code, "transcripts without any Entrez id say why no pathway can be built")
API.multiomics_engine.build_transcriptomics_input, API.multiomics_engine.build_proteomics_input = _bt, _bp

print("=== the result is tied to the data it was fitted on ===")
payF = API._job_mofa_run("job-m6", S)
check(MS["results"] is not None and MS["config"] is not None, "(a result and a config exist)")
API._job_correlation_resolve("job-m7", S, "9606", "Homo sapiens", False)
API.multiomics_mofa_options(x_session_id=SID)
check(MS["results"] is not None and MS["config"] is not None, "re-resolving IDENTICAL data keeps the MOFA result")
MS["inputs"]["transcriptomics"] = ([], {})
low = pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_proteomics.csv"))
low.iloc[0, 1] = low.iloc[0, 1] + 1.0
upload("proteomics", low.to_csv(index=False).encode(), SID)
API.multiomics_mofa_options(x_session_id=SID)
check(MS["results"] is None and MS["config"] is None and MS["sheet"] is None and MS["enrich"] == {} and MS["inputs"] == {},
      "SENSITIVITY: a re-uploaded table (resolution dropped) clears the MOFA config, result, enrichment, cached pathway sets and sheet")
API._job_correlation_resolve("job-m8", S, "9606", "Homo sapiens", False)
check(refused(API.multiomics_mofa_run, x_session_id=SID)[0] == 400, "...and after re-matching, a run needs a new configuration")
fp_a = API._mofa_fingerprint(CS)
upload("proteomics", P_EX, SID)
API._job_correlation_resolve("job-m9", S, "9606", "Homo sapiens", False)
check(API._mofa_fingerprint(CS) != fp_a, "SENSITIVITY: the fingerprint differs for different data")
ren = dict(CS["aligned"])
ren["samples"] = ["renamed"] + list(CS["aligned"]["samples"][1:])
check(API._mofa_fingerprint({"aligned": ren}) != API._mofa_fingerprint(CS), "SENSITIVITY: the fingerprint also differs when only a sample's name differs")

print("=== Joint selection keeps the MOFA state ===")
API.multiomics_mofa_configure({"factors": 4, "n_refits": 0}, x_session_id=SID)
API.multiomics_select({"omics_types": ["transcriptomics", "proteomics"]}, x_session_id=SID)
check(S["multiomics"]["mofa"]["config"] is not None and S["multiomics"]["mofa"] is MS, "choosing omics for Joint does not wipe MOFA's state")

print("=== too few samples ===")
S2 = "mo-mofa-test-2"
API.SESSIONS.pop(S2, None)
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "metabolomics"]}, x_session_id=S2)
SS = API.get_session(S2)
SS["multiomics"]["correlation"]["aligned"] = {"samples": ["a", "b", "c"], "mats": {}, "report": {}, "warnings": []}
code, detail = refused(API.multiomics_mofa_configure, {}, x_session_id=S2)
check(code == 400 and "at least 6 samples" in detail and "3 were matched" in detail, f"SENSITIVITY: fewer than 6 matched samples is refused with the count ({detail[:70]})")
API.SESSIONS.pop(S2, None)
API.SESSIONS.pop(SID, None)

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
