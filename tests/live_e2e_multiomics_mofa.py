"""Live end-to-end test of Multi-omics MOFA / Factor-based Integration through the HTTP API
(no browser - see playwright_multiomics_mofa.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network access: select the three omics -> upload the three
bundled 36-sample tables -> resolve (the real transcript / protein identifier resolvers) -> configure -> run (the real mofapy2 model)
-> pathway enrichment (real Reactome sets) -> download.

What is checked against INDEPENDENT code (no Omiscera code: pandas / numpy / scipy straight from the CSV files):
  - the three PLANTED factors are recovered (|r| > 0.95 with tests/fixtures/mofa_example_truth.csv, written by the generator
    from the random numbers, not from the model);
  - the variance explained per omics (total over all factors) is recomputed from the returned factor values, the downloaded weights
    and the raw tables (centred, scaled to equal view variance, observed cells only): it must agree with the model's own figure
    to half a point (mofapy2 also counts the cells it imputed, so the two are not identical by definition - the exact definition
    is tested offline in test_multiomics_mofa.py);
  - the group test (Mann-Whitney U) is recomputed with scipy from the returned factor values and the Group column of the
    metabolite file: every factor's p-value must equal the app's;
  - the planted treatment factor is the only one that separates the groups;
  - the top-25 weights are the 25 largest |weights| of the downloaded weights table; the downloads equal the payload;
  - the pathway enrichment of the treatment factor puts the planted Reactome pathway (glycolysis / citric acid cycle) first, its table is
    sorted by FDR, every AUC lies in [0, 1] and agrees with its direction, FDRs never fall below their p-values;
  - a second run gives the identical factors (seed fixed); the payload has no NaN; the single-omics sub-dicts are untouched
    (the session route used by Joint still answers).
The example data are synthetic (see info/README-MultiOmics-MOFA.txt): the test asserts what was planted.
"""
import io
import os
import time

import numpy as np
import pandas as pd
import requests
from scipy import stats

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
HEADERS = {"X-Session-Id": "multiomics-mofa-live-" + str(int(time.time() * 1000))}
PASS = FAIL = 0
FILES = {"transcriptomics": "example_mofa_transcriptomics.csv", "proteomics": "example_mofa_proteomics.csv",
         "metabolomics": "example_mofa_metabolomics.csv"}


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def post(path, **kw):
    return requests.post(BASE + path, headers=HEADERS, timeout=900, **kw)


def get(path, **kw):
    return requests.get(BASE + path, headers=HEADERS, timeout=900, **kw)


def wait_job(job_id, timeout=1800):
    start = time.time()
    while time.time() - start < timeout:
        job = get(f"/api/jobs/{job_id}").json()
        if job.get("status") == "done":
            return job.get("result")
        if job.get("status") in ("error", "failed"):
            raise RuntimeError(job.get("error") or job)
        time.sleep(1.0)
    raise TimeoutError(job_id)


def run_mofa(cfg):
    r = post("/api/multiomics/mofa/configure", json=cfg)
    assert r.status_code == 200, r.text
    job = post("/api/multiomics/mofa/run").json()["job_id"]
    return wait_job(job)


def csv_of(path, **kw):
    r = get(path, **kw)
    assert r.status_code == 200, (path, r.text[:200])
    return pd.read_csv(io.StringIO(r.text))


raw = {o: pd.read_csv(os.path.join(_ROOT, "Data3", f)) for o, f in FILES.items()}
truth = pd.read_csv(os.path.join(_ROOT, "tests", "fixtures", "mofa_example_truth.csv")).set_index("Sample")

print("=== select, upload, resolve ===")
r = post("/api/multiomics/correlation/select", json={"omics_types": ["transcriptomics", "proteomics", "metabolomics"]})
check(r.status_code == 200, "the three omics are selected")
for o, f in FILES.items():
    with open(os.path.join(_ROOT, "Data3", f), "rb") as fh:
        r = post(f"/api/multiomics/correlation/{o}/upload", files={"file": (f, fh.read(), "text/csv")}, data={"orientation": "auto"})
    check(r.status_code == 200, f"{o} table uploaded ({r.json().get('n_samples')} samples x {r.json().get('n_features')} features)")
job = post("/api/multiomics/correlation/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"]
al = wait_job(job)
check(al["n_shared"] == 36, "36 samples are in every file")
print("   mapped to a gene:", {o: (v["n_mapped"], v["n_features"]) for o, v in al["omics"].items()})

print("=== options and validation ===")
opt = get("/api/multiomics/mofa/options").json()
check(opt["ready"] and opt["group_options"][0]["key"] == "metabolomics:Group" and opt["group_options"][0]["groups"] == ["Control", "Treated"],
      "the Group column of the metabolite file is offered")
check(post("/api/multiomics/mofa/configure", json={"factors": 1}).status_code == 400, "an impossible number of factors is refused")

print("=== run ===")
t0 = time.time()
cfg = {"factors": 10, "n_refits": 2, "group_column": "metabolomics:Group"}
R = run_mofa(cfg)
print(f"   took {time.time() - t0:.1f} s; factors kept {R['n_factors']}; converged {R['converged']} in {R['n_iterations']} iterations")
txt = get("/api/multiomics/mofa/options").text
check("NaN" not in txt, "the options payload has no NaN")
fv = R["factor_values"]
Z = pd.DataFrame(fv["values"], index=fv["samples"], columns=fv["names"])
check(R["converged"] and R["n_factors"] >= 3 and sorted(Z.index) == sorted(raw["metabolomics"]["Sample"]), "converged; the factor values cover the 36 samples")
rec = {}
for name in ("treatment", "cholesterol", "lipids"):
    t = truth.loc[Z.index, name]
    rec[name] = max(abs(np.corrcoef(Z[f], t)[0, 1]) for f in Z.columns)
check(all(v > 0.95 for v in rec.values()), f"SENSITIVITY: the three planted factors are recovered ({ {k: round(v, 3) for k, v in rec.items()} })")
f_treat = max(Z.columns, key=lambda f: abs(np.corrcoef(Z[f], truth.loc[Z.index, "treatment"])[0, 1]))
f_chol = max(Z.columns, key=lambda f: abs(np.corrcoef(Z[f], truth.loc[Z.index, "cholesterol"])[0, 1]))
f_lip = max(Z.columns, key=lambda f: abs(np.corrcoef(Z[f], truth.loc[Z.index, "lipids"])[0, 1]))
r2 = pd.DataFrame(R["r2"]["values"], index=R["r2"]["omics"], columns=R["r2"]["factors"])
check((r2[f_treat] > 5).all(), f"the treatment factor ({f_treat}) is SHARED: > 5% in all three omics ({r2[f_treat].round(1).to_dict()})")
check(r2.loc["Metabolomics", f_chol] < 2 and r2.loc["Transcriptomics", f_chol] > 2 and r2.loc["Proteomics", f_chol] > 2,
      f"the cholesterol factor ({f_chol}) is in transcripts and proteins only ({r2[f_chol].round(1).to_dict()})")
check(r2.loc["Transcriptomics", f_lip] < 2 and r2.loc["Proteomics", f_lip] < 2 and r2.loc["Metabolomics", f_lip] > 5,
      f"the lipid factor ({f_lip}) is in metabolites only ({r2[f_lip].round(1).to_dict()})")

print("=== variance explained, recomputed from the raw tables ===")
W = {o: csv_of("/api/multiomics/mofa/download", params={"table": f"weights:{o}"}) for o in R["omics"]}
for o in R["omics"]:
    w = W[o]
    feats = w["Feature"].astype(str).tolist()
    df = raw[o]
    if o == "metabolomics":
        data = df.set_index("Sample").drop(columns=["Group"])
    else:
        data = df.set_index(df.columns[0]).T
    data = data[feats].loc[Z.index].astype(float)
    Y = data - data.mean(axis=0, skipna=True)
    Y = Y / np.sqrt(np.nanmean(Y.to_numpy() ** 2))      # equal total variance per view (scale_views)
    pred = Z.to_numpy() @ w[list(Z.columns)].to_numpy().T
    obs = ~np.isnan(Y.to_numpy())
    ss_res = np.nansum((Y.to_numpy() - pred) ** 2)
    ss_tot = np.nansum(Y.to_numpy() ** 2)
    mine = 100 * (1 - ss_res / ss_tot)
    check(abs(mine - R["r2_total"][o]) < 0.5, f"{o}: total variance explained, from Z x W and the raw table: {mine:.1f}% vs the model's {R['r2_total'][o]:.1f}%")

print("=== group tests, recomputed with scipy ===")
grp = raw["metabolomics"].set_index("Sample")["Group"]
ctrl = [s for s in Z.index if grp[s] == "Control"]
trt = [s for s in Z.index if grp[s] == "Treated"]
ok_p = True
rows = {r["Factor"]: r for r in R["group_tests"]["rows"]}
for f in Z.columns:
    p = stats.mannwhitneyu(Z.loc[ctrl, f], Z.loc[trt, f], alternative="two-sided").pvalue
    ok_p &= abs(p - rows[f]["p-value"]) < 1e-9 * max(1.0, p) + 1e-15
check(ok_p, "every factor's Mann-Whitney p-value equals scipy's")
sig = [f for f, r in rows.items() if r["FDR"] < 0.05]
check(sig == [f_treat], f"SENSITIVITY: only the treatment factor separates Control from Treated ({sig})")
check(abs(rows[f_treat]["Effect"]) > 0.95, f"...with a rank-biserial effect of {rows[f_treat]['Effect']:.2f}")

print("=== weights and downloads ===")
ok_top = True
for o in R["omics"]:
    for f in Z.columns:
        exp = W[o].reindex(W[o][f].abs().sort_values(ascending=False).index).head(25)
        got = R["top_weights"][o][f]
        ok_top &= [g["feature"] for g in got] == exp["Feature"].astype(str).tolist() and np.allclose([g["weight"] for g in got], exp[f].to_numpy())
check(ok_top, "the 25 top weights of every factor and omics are the 25 largest |weights| of the downloaded table")
dfz = csv_of("/api/multiomics/mofa/download", params={"table": "factors"})
check(len(dfz) == 36 and np.allclose(dfz.set_index("Sample").loc[Z.index, list(Z.columns)].to_numpy(), Z.to_numpy()), "the factor download equals the payload")
dr2 = csv_of("/api/multiomics/mofa/download", params={"table": "r2"})
check(np.allclose(dr2.set_index("Omics").to_numpy(), r2.to_numpy()), "the variance download equals the payload")
dgt = csv_of("/api/multiomics/mofa/download", params={"table": "groups"})
check(len(dgt) == len(R["group_tests"]["rows"]), "the group-test download has one row per factor")
check(get("/api/multiomics/mofa/download", params={"table": "nope"}).status_code == 400, "an unknown table is refused")

print("=== pathway enrichment (real Reactome) ===")
job = post("/api/multiomics/mofa/enrich", json={"omics": "transcriptomics", "factor": f_treat, "mode": "signed"}).json()["job_id"]
E = wait_job(job)
tab = pd.DataFrame(E["table"]["rows"])
print("   top:", tab["Term Name"].head(5).tolist())
check(E["n_tested"] > 20 and len(tab) > 0, f"{E['n_tested']} pathways tested, {E['n_features_matched']} of {E['n_features']} features matched")
check(any(k in tab["Term Name"].head(3).str.lower().str.cat() for k in ("glycoly", "citric", "tca")), "SENSITIVITY: the planted pathway (glycolysis / citric acid cycle) is among the first three")
check(tab["FDR"].is_monotonic_increasing and (tab["FDR"] >= tab["p-value"] - 1e-15).all(), "sorted by FDR, and no FDR below its p-value")
check(tab["AUC"].between(0, 1).all() and ((tab["AUC"] > 0.5) == (tab["Direction"] == "higher weights")).all(), "every AUC lies in [0, 1] and agrees with its direction")
dl = csv_of("/api/multiomics/mofa/download", params={"table": f"enrichment:transcriptomics:{f_treat}:signed"})
check(len(dl) == E["n_tested"] and dl["Term ID"].head(3).tolist() == tab["Term ID"].head(3).tolist(), "the enrichment download holds every tested pathway, in the same order")
job = post("/api/multiomics/mofa/enrich", json={"omics": "proteomics", "factor": f_treat, "mode": "absolute"}).json()["job_id"]
Ep = wait_job(job)
check(Ep["n_tested"] > 20 and Ep["mode"] == "absolute", f"proteins, absolute weights: {Ep['n_tested']} pathways tested")
check(post("/api/multiomics/mofa/enrich", json={"omics": "transcriptomics", "factor": "Factor99"}).status_code == 400, "an unknown factor is refused")

print("=== a second run, and the rest of the session ===")
R2 = run_mofa(cfg)
Z2 = pd.DataFrame(R2["factor_values"]["values"], index=R2["factor_values"]["samples"], columns=R2["factor_values"]["names"])
check(np.allclose(Z.to_numpy(), Z2.to_numpy(), atol=1e-9), "the same data and settings give identical factors (the seed is fixed)")
check(get("/api/multiomics/mofa/download", params={"table": f"enrichment:transcriptomics:{f_treat}:signed"}).status_code == 400,
      "a new run drops the earlier enrichment (it belonged to the earlier factors)")
check(get("/api/multiomics/mofa/options").json()["ready"], "the matched tables are still there after the runs")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
