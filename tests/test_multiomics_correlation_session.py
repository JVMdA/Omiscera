# -*- coding: utf-8 -*-
"""Direct unit tests for the Correlation-based Integration routes and jobs in API.py: no server, no network - imports
API.py and drives the route functions / jobs in-process, with the identifier resolvers stubbed (the real ones are
exercised by live_e2e_multiomics_correlation.py).

THE THINGS UNDER TEST
  1. State: session["multiomics"]["correlation"] is its own slot; choosing omics for Joint does not wipe it; a re-upload
     drops the resolution, the alignment and the result.
  2. Validation: Genomics refused (a list cannot be correlated across samples), <2 omics, an unselected omics, a file that
     is not a matrix, no organism, bad numbers.
  3. The resolve job: transcript / protein features are resolved with their own mapping, metabolites are NOT (nothing needs
     their database ids), samples are matched across the files, a failure commits nothing, a too-small overlap says which
     file is probably the wrong way round.
  4. The run job on the REAL example files: the payload is strict-JSON safe, the planted structure is in it, the scatter
     endpoint returns the pair's values with the correlation recomputed, downloads carry the ids, and nothing leaks into
     the four single-omics sub-dicts or into Joint / Mapping.
"""
import asyncio
import json
import os
import sys

import pandas as pd
from fastapi import HTTPException

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import API  # noqa: E402

mc = API.multiomics_correlation
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


def upload(omics, data, sid, orientation="auto"):
    return asyncio.run(API.multiomics_correlation_upload(omics, FakeUpload(data), orientation, x_session_id=sid))


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


# the example files: row i of the RNA table is the same gene as row i of the protein table (the generator writes them in step)
T_EX, P_EX = ex("example_corr_transcriptomics.csv"), ex("example_corr_proteomics.csv")
M_EX = ex("example_corr_metabolomics.csv")
_t = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_transcriptomics.csv"))
_p = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_proteomics.csv"))
GENE_OF = {"transcriptomics": {g: g for g in _t["Gene"]}, "proteomics": dict(zip(_p["Protein IDs"], _t["Gene"]))}
calls = []


def fake_resolve(omics_type, ids, taxon_id, progress_callback=None):
    calls.append((omics_type, list(ids), taxon_id))
    mp = {i: GENE_OF[omics_type][i] for i in ids if i in GENE_OF[omics_type]}
    return {"gene_of": mp, "n_features": len(ids), "n_mapped": len(mp), "n_unmapped": len(ids) - len(mp), "identifier_type": "stub"}


_real_resolve = API.multiomics_retriever.resolve_matrix_features
API.multiomics_retriever.resolve_matrix_features = fake_resolve

print("=== state ===")
state = API.new_session_state()
cs0 = state["multiomics"]["correlation"]
check(cs0 == {"selected": [], "taxon_id": "", "taxon_name": "", "uploads": {}, "resolved": {}, "aligned": None, "config": None, "results": None},
      "a new session carries an empty correlation slot, separate from Joint and Mapping")
check(API._fresh_multiomics_correlation() is not API._fresh_multiomics_correlation()
      and API._fresh_multiomics_correlation()["uploads"] is not API._fresh_multiomics_correlation()["uploads"], "each call returns independent dicts")

print("=== select ===")
SID = "mo-corr-test"
API.SESSIONS.pop(SID, None)
code, detail = refused(API.multiomics_correlation_select, {"omics_types": ["genomics", "transcriptomics"]}, x_session_id=SID)
check(code == 400 and "gene list" in detail, f"SENSITIVITY: Genomics is refused, with the reason ({detail[:70]})")
check(refused(API.multiomics_correlation_select, {"omics_types": ["transcriptomics"]}, x_session_id=SID)[0] == 400
      and refused(API.multiomics_correlation_select, {"omics_types": []}, x_session_id=SID)[0] == 400, "fewer than two omics is refused")
out = API.multiomics_correlation_select({"omics_types": ["metabolomics", "transcriptomics", "proteomics", "proteomics"]}, x_session_id=SID)
check(out["selected"] == ["transcriptomics", "proteomics", "metabolomics"], f"the selection is de-duplicated and in canonical order ({out['selected']})")
S = API.get_session(SID)
CS = S["multiomics"]["correlation"]

print("=== upload ===")
before = snapshot(S)
check(refused(lambda: asyncio.run(API.multiomics_correlation_upload("genomics", FakeUpload(T_EX), "auto", x_session_id=SID)))[0] == 400,
      "an omics that is not selected is refused")
res = upload("transcriptomics", T_EX, SID)
check(res["orientation"] == "columns" and res["n_samples"] == 12 and res["n_features"] == 299 and res["first_samples"][0] == "Control_1",
      f"RNA table: samples in the columns, detected and reported ({res['orientation_reason']})")
res = upload("metabolomics", M_EX, SID)
check(res["orientation"] == "rows" and res["set_aside_columns"] == ["Group"] and res["n_features"] == 42,
      "metabolite table: samples in the rows, the Group column set aside")
res = upload("proteomics", P_EX, SID)
check(res["n_samples"] == 12 and res["n_features"] == 299, "protein table read")
code, detail = (None, "")
try:
    asyncio.run(API.multiomics_correlation_upload("proteomics", FakeUpload(b"just,text\nhello,world\n"), "auto", x_session_id=SID))
except HTTPException as exc:
    code, detail = exc.status_code, str(exc.detail)
check(code == 400 and detail, f"a file that is not a matrix is refused with a message ({detail[:60]})")
check(CS["uploads"]["proteomics"]["df"].shape == (12, 299), "...and the earlier valid upload is untouched by the failed one")
res = upload("transcriptomics", T_EX, SID, orientation="rows")
check(res["orientation"] == "rows" and res["n_samples"] == 299, "the orientation you set wins over the guess (here deliberately wrong)")
upload("transcriptomics", T_EX, SID)
check(snapshot(S) == before, "ISOLATION: uploads left the four single-omics sub-dicts untouched")

print("=== resolve job ===")
check(refused(API.multiomics_correlation_resolve, {"taxon_id": ""}, x_session_id=SID)[0] == 400, "resolve without an organism is refused")
S2 = "mo-corr-test-2"
API.SESSIONS.pop(S2, None)
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "proteomics"]}, x_session_id=S2)
check(refused(API.multiomics_correlation_resolve, {"taxon_id": "9606"}, x_session_id=S2)[0] == 400, "resolve before uploading names what is missing")
API.SESSIONS.pop(S2, None)
CS["resolved"] = {}
CS["aligned"] = {"stale": True}
CS["results"] = {"stale": True}
calls.clear()
API.multiomics_correlation_resolve  # noqa: B018  (route exists)
pay = API._job_correlation_resolve("job-r", S, "9606", "Homo sapiens", False)
check([c[0] for c in calls] == ["transcriptomics", "proteomics"] and calls[0][2] == "9606",
      f"SENSITIVITY: transcript and protein features are resolved, once each, with the chosen organism - the metabolites are NOT ({[c[0] for c in calls]})")
check(calls[0][1][:2] == ["CS", "ACO2"] and len(calls[0][1]) == 299, "the resolver gets the matrix's feature ids as written")
check(pay["n_shared"] == 12 and pay["omics"]["transcriptomics"]["n_mapped"] == 299 and pay["omics"]["metabolomics"]["n_mapped"] is None
      and pay["omics"]["proteomics"]["n_used"] == 12 and pay["omics"]["proteomics"]["dropped"] == [],
      "the alignment card: 12 shared samples, nothing dropped, resolution counts for T and P, none for metabolites")
check(CS["aligned"]["samples"][0] == "Control_1" and CS["results"] is None and set(CS["resolved"]) == {"transcriptomics", "proteomics", "metabolomics"}
      and CS["taxon_id"] == "9606", "the job commits the alignment and drops the earlier result")
json.dumps(pay, allow_nan=False)
check(True, "the alignment payload is strict-JSON safe")
before_fail = CS["aligned"]
short = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_proteomics.csv"))
short = short.rename(columns={"Control_1": "X1", "Control_2": "X2", "Control_3": "X3", "Control_4": "X4", "Control_5": "X5", "Control_6": "X6",
                              "Treated_1": "Y1", "Treated_2": "Y2"})
upload("proteomics", short.to_csv(index=False).encode(), SID)
try:
    API._job_correlation_resolve("job-r2", S, "9606", "Homo sapiens", False)
    msg = ""
except ValueError as exc:
    msg = str(exc)
check("Only 4 sample(s) are in every file" in msg, f"SENSITIVITY: too few shared samples fails with the count ({msg[:80]})")
check(CS["aligned"] is None, "a re-upload had already dropped the old alignment (nothing stale is kept)")
upload("proteomics", P_EX, SID)
CS["aligned"] = before_fail
transposed = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_proteomics.csv")).set_index("Protein IDs").T.reset_index()
upload("proteomics", transposed.to_csv(index=False).encode(), SID, orientation="columns")
try:
    API._job_correlation_resolve("job-r3", S, "9606", "Homo sapiens", False)
    msg = ""
except ValueError as exc:
    msg = str(exc)
check("other orientation" in msg and "Proteomics" in msg, f"SENSITIVITY: a protein file the wrong way round is diagnosed ({msg[:110]})")
upload("proteomics", P_EX, SID)
low = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_proteomics.csv"))
low = low.rename(columns={c: c.lower() for c in low.columns if c != "Protein IDs"})
upload("proteomics", low.to_csv(index=False).encode(), SID)
try:
    API._job_correlation_resolve("job-r4", S, "9606", "Homo sapiens", False)
    msg = ""
except ValueError as exc:
    msg = str(exc)
check("sample(s)" in msg, "lower-case sample ids do not match by default...")
pay_ci = API._job_correlation_resolve("job-r5", S, "9606", "Homo sapiens", True)
check(pay_ci["n_shared"] == 12, "...and do with the case-insensitive option")
upload("proteomics", P_EX, SID)
pay = API._job_correlation_resolve("job-r6", S, "9606", "Homo sapiens", False)
check(CS["aligned"] is not None, "back to a good state")
dup_t = pd.concat([_t, _t[_t["Gene"] == "CS"]], ignore_index=True)        # the CS row twice (two probes of one gene)
upload("transcriptomics", dup_t.to_csv(index=False).encode(), SID)
calls.clear()
API._job_correlation_resolve("job-dup", S, "9606", "Homo sapiens", False)
t_ids = calls[0][1]
gof = CS["resolved"]["transcriptomics"]["gene_of"]
check(t_ids.count("CS") == 1 and "CS (2)" not in t_ids and gof.get("CS (2)") == "CS" and gof.get("CS") == "CS"
      and CS["resolved"]["transcriptomics"]["n_mapped"] == 300,
      "SENSITIVITY: a repeated id renamed 'CS (2)' is resolved by its original id, so BOTH copies reach the gene (300 of 300 mapped)")
upload("transcriptomics", T_EX, SID)
pay = API._job_correlation_resolve("job-r7", S, "9606", "Homo sapiens", False)

print("=== configure ===")
S3 = "mo-corr-test-3"
API.SESSIONS.pop(S3, None)
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "proteomics"]}, x_session_id=S3)
check(refused(API.multiomics_correlation_configure, {}, x_session_id=S3)[0] == 400, "configure before the samples are matched is refused")
API.SESSIONS.pop(S3, None)
cfg = API.multiomics_correlation_configure({}, x_session_id=SID)["config"]
check(cfg == mc.CONFIG_DEFAULTS, "no payload gives the documented defaults")
for bad in ({"method": "kendall"}, {"fdr": 0}, {"fdr": 2}, {"rho_cut": 1.5}, {"top_variance": 10}, {"top_variance": "x"}, {"min_samples": 3}, {"heat_k": 500}, {"net_k": 1},
            {"min_samples": 20}):
    check(refused(API.multiomics_correlation_configure, bad, x_session_id=SID)[0] == 400, f"refused: {bad}")
cfg = API.multiomics_correlation_configure({"method": "pearson", "fdr": 0.1, "rho_cut": 0.5, "log2": True, "top_variance": 100}, x_session_id=SID)["config"]
check(cfg["method"] == "pearson" and cfg["fdr"] == 0.1 and cfg["log2"] is True and cfg["top_variance"] == 100 and CS["config"] == cfg, "a valid config is stored")
CS["results"] = {"stale": True}
API.multiomics_correlation_configure({}, x_session_id=SID)
check(CS["results"] is None, "re-configuring drops the previous result")

print("=== run job on the example files ===")
check(refused(API.multiomics_correlation_run, x_session_id="mo-corr-none")[0] == 400, "run before configure is refused")
before = snapshot(S)
out = API._job_correlation_run("job-run", S)
json.dumps(out, allow_nan=False)
check(True, "the whole result payload is strict-JSON safe")
check(set(out["blocks"]) == {"transcriptomics|proteomics", "transcriptomics|metabolomics", "proteomics|metabolomics"} and out["samples"][0] == "Control_1" and len(out["samples"]) == 12,
      "three blocks for three omics, over the 12 shared samples")
conc = {r["Gene"]: r for r in out["concordance"]["table"]["rows"]}
planted = [g for g in _t["Gene"]][:15]
n_sig_pl = sum(1 for g in planted if conc[g]["FDR"] < 0.05)
check(all(conc[g]["rho"] > 0.5 for g in planted) and n_sig_pl >= 8,
      f"BIOLOGY (planted): the 15 concordant genes are concordant - rho > 0.5 for all, FDR < 0.05 for {n_sig_pl} (at 12 samples with a missing cell or two, not every one clears BH over 299 genes)")
neg = [g for g in _t["Gene"]][15:19]
check(all(conc[g]["rho"] < -0.3 for g in neg), "the planted DISCORDANT genes come out negative")
bg = [g for g in _t["Gene"]][29:]
check(sum(1 for g in bg if conc[g]["FDR"] < 0.05) <= 3 and abs(out["concordance"]["median_rho"]) < 0.3,
      f"the 270 background genes are not concordant (few significant, median rho {out['concordance']['median_rho']:.2f})")
tm = out["blocks"]["transcriptomics|metabolomics"]
check(any(r["Feature A"] == "CS" and r["Feature B"] == "Citrate" for r in tm["table"]["rows"]), "BIOLOGY (planted): CS is correlated with citrate")
check(out["blocks"]["transcriptomics|proteomics"]["n_tests"] > 40000 and tm["n_tests"] > 5000, "the block sizes are the tested pairs")
check(out["blocks"]["transcriptomics|proteomics"]["table"]["columns"] == ["Feature A", "Feature B", "Omics A", "Omics B", "rho", "p-value", "FDR", "n"]
      and len(out["blocks"]["transcriptomics|proteomics"]["ids"]) == len(out["blocks"]["transcriptomics|proteomics"]["table"]["rows"]),
      "the pair table has its columns and one id pair per row (for the scatter viewer)")
check(out["alignment"]["n_shared"] == 12 and out["feature_report"]["proteomics"]["n_too_many_missing"] == 0, "the alignment and the feature report travel with the result")
check(out["network"]["nodes"] and out["network"]["edges"] and out["network"]["crossings"] >= 0, "a network layout is returned")
check(snapshot(S) == before, "ISOLATION: a full run left the four single-omics sub-dicts untouched")
check(S["multiomics"]["joint_results"] is None and S["multiomics"]["mapping"] == API._fresh_multiomics_mapping(), "...and Joint and Mapping")

print("=== scatter, download ===")
sc = API.multiomics_correlation_scatter(kind="gene", gene="CS", x_session_id=SID)
check(sc["samples"][0] == "Control_1" and len(sc["x"]) == 12 and abs(sc["rho"] - conc["CS"]["rho"]) < 1e-9 and sc["n"] == conc["CS"]["n"],
      f"the gene scatter returns the 12 values per sample and the SAME rho and n as the table ({sc['rho']:.3f}, n = {sc['n']})")
check(sc["x_label"] == "Transcript CS" and sc["y_label"] == "Protein CS" and all(isinstance(v, float) for v in sc["x"]),
      "SENSITIVITY: the axis labels are their own keys - they once overwrote the x / y value arrays (the response had 13 'values')")
row = next(r for r in tm["table"]["rows"] if r["Feature A"] == "CS" and r["Feature B"] == "Citrate")
ids = tm["ids"][tm["table"]["rows"].index(row)]
sc2 = API.multiomics_correlation_scatter(kind="pair", block="transcriptomics|metabolomics", a=ids[0], b=ids[1], x_session_id=SID)
check(abs(sc2["rho"] - row["rho"]) < 1e-9 and abs(sc2["p"] - row["p-value"]) < 1e-9 * max(1, row["p-value"]) + 1e-15,
      "the pair scatter recomputes the table's rho and p exactly")
check(sc2["x_label"] == "Transcriptomics: CS" and sc2["y_label"] == "Metabolomics: Citrate", f"pair scatter axes: {sc2['x_label']} / {sc2['y_label']}")
acc = _p["Protein IDs"][_t["Gene"] == "CS"].iloc[0]
sc3 = API.multiomics_correlation_scatter(kind="pair", block="proteomics|metabolomics", a=acc, b="Citrate", x_session_id=SID)
check(sc3["x_label"] == f"Proteomics: CS ({acc})", f"a protein axis names its gene AND its accession ({sc3['x_label']})")
# a protein with a missing value: the scatter's rho is the exact pairwise Spearman the table uses
miss_acc = next(a for a in _p["Protein IDs"] if _p.set_index("Protein IDs").loc[a].isna().any() and a in CS["results"]["mats"]["proteomics"].columns)
sc4 = API.multiomics_correlation_scatter(kind="pair", block="proteomics|metabolomics", a=miss_acc, b="Citrate", x_session_id=SID)
Rm, Nm = mc.correlate_all(CS["results"]["mats"]["proteomics"][[miss_acc]].to_numpy(), CS["results"]["mats"]["metabolomics"][["Citrate"]].to_numpy(), "spearman", 5)
check(sc4["n"] == 11 and abs(sc4["rho"] - Rm[0, 0]) < 1e-12, f"a pair with a missing value: the scatter (n = {sc4['n']}) and the block agree exactly")
check(refused(API.multiomics_correlation_scatter, "pair", "transcriptomics|metabolomics", "nope", "Citrate", "", x_session_id=SID)[0] == 404
      and refused(API.multiomics_correlation_scatter, "gene", "", "", "", "NOPE", x_session_id=SID)[0] == 404, "an unknown feature is a 404")
check(refused(API.multiomics_correlation_scatter, "gene", "", "", "", "CS", x_session_id="mo-corr-none")[0] == 400, "no result is a 400")
dl = API.multiomics_correlation_download("concordance", x_session_id=SID)
check(dl.media_type == "text/csv" and dl.body.startswith(b"Gene,rho,p-value,FDR,n"), "the concordance table downloads as CSV")
dl2 = API.multiomics_correlation_download("transcriptomics|metabolomics", x_session_id=SID)
check(b"Feature A id" in dl2.body and b"CS" in dl2.body, "a block downloads with the feature ids")
check(refused(API.multiomics_correlation_download, "nope", x_session_id=SID)[0] == 400, "an unknown table is refused")

print("=== options reach the engine ===")
API.multiomics_correlation_configure({"log2": True}, x_session_id=SID)
out_l = API._job_correlation_run("job-l", S)
check(out_l["config"]["log2"] is True and out_l["concordance"]["n_common"] == out["concordance"]["n_common"], "log2 runs (the example values are positive)")
neg_t = pd.read_csv(os.path.join(_ROOT, "Data3", "example_corr_transcriptomics.csv"))
neg_t.iloc[0, 1] = -5
upload("transcriptomics", neg_t.to_csv(index=False).encode(), SID)
API._job_correlation_resolve("job-n", S, "9606", "Homo sapiens", False)
API.multiomics_correlation_configure({"log2": True}, x_session_id=SID)
try:
    API._job_correlation_run("job-n2", S)
    msg = ""
except ValueError as exc:
    msg = str(exc)
check("negative" in msg, f"log2 of negative values fails with a message ({msg[:60]})")
check(S["multiomics"]["correlation"]["results"] is None, "a failed run leaves no half-built result")
upload("transcriptomics", T_EX, SID)
check(CS["results"] is None and CS["aligned"] is None and "transcriptomics" not in CS["resolved"], "SENSITIVITY: a re-upload drops the resolution, the alignment and the result")
API._job_correlation_resolve("job-r", S, "9606", "Homo sapiens", False)
API.multiomics_correlation_configure({}, x_session_id=SID)
API._job_correlation_run("job-r2", S)
check(CS["results"] is not None and CS["aligned"] is not None, "(a fresh alignment and result exist before the omics are re-chosen)")
API.multiomics_correlation_select({"omics_types": ["transcriptomics", "proteomics"]}, x_session_id=SID)
check(CS["results"] is None and CS["aligned"] is None and set(CS["uploads"]) == {"transcriptomics", "proteomics"} and "metabolomics" not in CS["resolved"],
      "SENSITIVITY: choosing different omics drops the old alignment and result, and the upload of an omics no longer chosen")
select = API.multiomics_select({"omics_types": ["genomics", "proteomics"]}, x_session_id=SID)
check(S["multiomics"]["correlation"] is CS and set(CS["uploads"]) == {"transcriptomics", "proteomics"} and CS["selected"] == ["transcriptomics", "proteomics"],
      "choosing omics for Joint leaves the correlation slot alone (uploads and choice untouched)")
check(S["multiomics"]["correlation"] is CS, "SENSITIVITY: the very same correlation dict survives Joint's select (it used to rebuild the whole multiomics state)")

API.multiomics_retriever.resolve_matrix_features = _real_resolve
API.SESSIONS.pop(SID, None)
print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
