"""Live end-to-end test of Multi-omics Correlation-based Integration through the HTTP API
(no browser - see playwright_multiomics_correlation.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network access: select the three omics -> upload the
three bundled matched-sample tables -> resolve (the real transcript / protein identifier resolvers) -> run -> download.

Every number is checked against an INDEPENDENT recomputation straight from the uploaded CSV files, written here without any
Omiscera code (pandas average ranks, elementwise numpy Pearson over the samples both features have, the t distribution for the
p-value from scipy.stats.t, Benjamini-Hochberg by its definition):
  - EVERY gene of the same-gene concordance table (rho, n) - the pairing of transcript and protein is taken from the file's
    row order (row i of the RNA table and of the protein table are the same gene by construction), not from the app's own
    identifier mapping, so a wrong mapping would show up as a wrong rho;
  - the COMPLETE RNA x metabolite block: all genes x all compounds are recomputed, BH-corrected, and the set of pairs with
    FDR <= 0.05 and |rho| >= 0.6 must equal the app's reported set EXACTLY (nothing missed, nothing invented);
  - the planted structure must be recovered: CS - citrate, IDH3A - isocitrate, MDH2 - malate, and the planted discordant genes
    come out NEGATIVE;
  - every download equals the payload table it came from; the scatter endpoint returns the same rho as the table.
The example data are synthetic (see info/README-MultiOmics-Correlation.txt): the test asserts what was planted.
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
HEADERS = {"X-Session-Id": "multiomics-corr-live-" + str(int(time.time() * 1000))}
PASS = FAIL = 0
FILES = {"transcriptomics": "example_corr_transcriptomics.csv", "proteomics": "example_corr_proteomics.csv",
         "metabolomics": "example_corr_metabolomics.csv"}


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
        if job.get("status") == "error":
            raise RuntimeError(job.get("error") or job)
        time.sleep(2)
    raise TimeoutError(job_id)


# ---------------------------------------------------------------- independent statistics (no Omiscera code)
def spearman_matrix(A, B):
    """A: samples x p, B: samples x q (DataFrames, NaN = missing). The textbook pairwise-complete Spearman, one pair at a
    time (as R's cor(use = "pairwise.complete.obs", method = "spearman")): keep the samples BOTH have, rank each over
    those samples (pandas average ranks), Pearson of the ranks written out by hand (np.corrcoef crashes here).
    Returns rho, n arrays (p x q)."""
    av, bv = A.to_numpy(float), B.to_numpy(float)
    rho = np.full((av.shape[1], bv.shape[1]), np.nan)
    nn = np.zeros(rho.shape, int)
    for i in range(av.shape[1]):
        for j in range(bv.shape[1]):
            ok = ~np.isnan(av[:, i]) & ~np.isnan(bv[:, j])
            nn[i, j] = int(ok.sum())
            if ok.sum() < 3:
                continue
            x = pd.Series(av[ok, i]).rank(method="average").to_numpy()
            y = pd.Series(bv[ok, j]).rank(method="average").to_numpy()
            x, y = x - x.mean(), y - y.mean()
            den = np.sqrt((x * x).sum() * (y * y).sum())
            rho[i, j] = (x * y).sum() / den if den > 0 else np.nan
    return rho, nn


def pvalues(r, n):
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    p = np.where(np.abs(r) >= 1.0, 0.0, p)
    return np.where(n >= 5, p, np.nan)


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    adj = pv[order] * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


# ---------------------------------------------------------------- raw tables, read independently
raw_t = pd.read_csv(os.path.join(_ROOT, "Data3", FILES["transcriptomics"]), index_col=0)
raw_p = pd.read_csv(os.path.join(_ROOT, "Data3", FILES["proteomics"]), index_col=0)
raw_m = pd.read_csv(os.path.join(_ROOT, "Data3", FILES["metabolomics"]), index_col=0).drop(columns=["Group"])
T, P, M = raw_t.T, raw_p.T, raw_m          # samples x features
samples = list(M.index)
assert list(T.index) == samples == list(P.index)

print("=== select, upload, resolve ===")
r = post("/api/multiomics/correlation/select", json={"omics_types": ["transcriptomics", "proteomics", "metabolomics"]})
check(r.status_code == 200, "three omics selected")
for o, f in FILES.items():
    r = post(f"/api/multiomics/correlation/{o}/upload", files={"file": (f, open(os.path.join(_ROOT, "Data3", f), "rb"))}, data={"orientation": "auto"})
    check(r.status_code == 200, f"{o}: uploaded ({r.json().get('n_samples')} samples x {r.json().get('n_features')} features, {r.json().get('orientation')})")
r = post("/api/multiomics/correlation/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens", "case_insensitive": False})
al = wait_job(r.json()["job_id"])
check(al["n_shared"] == 12 and al["samples"] == samples, "12 shared samples, in the file's order")
check(al["omics"]["transcriptomics"]["n_mapped"] >= 280 and al["omics"]["proteomics"]["n_mapped"] >= 280,
      f"identifiers mapped by the real resolvers (RNA {al['omics']['transcriptomics']['n_mapped']}, protein {al['omics']['proteomics']['n_mapped']} of 299)")

print("=== run ===")
check(post("/api/multiomics/correlation/configure", json={"method": "spearman"}).status_code == 200, "configured (Spearman, defaults)")
t0 = time.time()
res = wait_job(post("/api/multiomics/correlation/run", json={}).json()["job_id"])
print(f"   run took {time.time() - t0:.1f}s")
cfg = res["config"]
check(cfg["fdr"] == 0.05 and cfg["rho_cut"] == 0.6, "default cutoffs FDR 0.05 and |rho| 0.6")

print("=== same-gene concordance vs independent recomputation (every gene) ===")
conc = res["concordance"]["table"]["rows"]
check(len(conc) >= 270, f"{len(conc)} genes tested")
sym = list(raw_t.index)
paired_rho = {}
for i in range(len(sym)):
    r1, n1 = spearman_matrix(T.iloc[:, [i]], P.iloc[:, [i]])
    paired_rho[sym[i]] = (r1[0, 0], n1[0, 0])
worst = 0.0
bad = []
for row in conc:
    g = row["Gene"]
    if g not in paired_rho:
        bad.append((g, "not in the RNA file"))
        continue
    rr, nn = paired_rho[g]
    if int(row["n"]) != int(nn) or abs(float(row["rho"]) - rr) > 5e-4:
        bad.append((g, row["rho"], rr, row["n"], nn))
    worst = max(worst, abs(float(row["rho"]) - rr))
check(not bad, f"every gene's rho and n equal the independent recomputation (largest difference {worst:.2e}; {len(bad)} disagree: {bad[:3]})")
p_ind = pvalues(np.array([paired_rho[x["Gene"]][0] for x in conc]), np.array([paired_rho[x["Gene"]][1] for x in conc]))
check(np.allclose(p_ind, [float(x["p-value"]) for x in conc], rtol=2e-3, atol=1e-12), "every gene's p-value equals the t-distribution p-value")
q_ind = bh(p_ind)
check(np.allclose(q_ind, [float(x["FDR"]) for x in conc], rtol=3e-3, atol=1e-12), "every gene's FDR equals Benjamini-Hochberg over the genes tested")
byg = {x["Gene"]: x for x in conc}
check(float(byg["CS"]["rho"]) > 0.9 and float(byg["IDH3A"]["rho"]) > 0.9, "PLANTED: CS and IDH3A transcripts track their proteins (rho > 0.9)")
check(res["concordance"]["n_sig_pos"] >= 8 and res["concordance"]["n_sig_neg"] >= 2,
      f"PLANTED: concordant ({res['concordance']['n_sig_pos']}) and discordant ({res['concordance']['n_sig_neg']}) genes are both recovered")
neg_genes = [x["Gene"] for x in conc if float(x["rho"]) < -0.8 and float(x["FDR"]) < 0.05]
check(len(neg_genes) >= 2 and all(paired_rho[g][0] < -0.8 for g in neg_genes), f"the discordant genes are negative in the raw files too ({neg_genes})")

print("=== the complete RNA x metabolite block, rebuilt independently ===")
blk = res["blocks"]["transcriptomics|metabolomics"]
Tk = T.loc[:, T.isna().mean() <= 0.2]
Tk = Tk.loc[:, Tk.std(ddof=0) > 0]
Mk = M.loc[:, M.isna().mean() <= 0.2]
Mk = Mk.loc[:, Mk.std(ddof=0) > 0]
rho_b, n_b = spearman_matrix(Tk, Mk)
p_b = pvalues(rho_b, n_b)
q_b = bh(p_b.ravel()).reshape(p_b.shape)
want = {(Tk.columns[i], Mk.columns[j]) for i, j in zip(*np.where((q_b <= 0.05) & (np.abs(rho_b) >= 0.6)))}
got = {(a, b) for a, b in blk["ids"]}
check(blk["n_tests"] == rho_b.size, f"the number of tests is the full block ({blk['n_tests']} = {Tk.shape[1]} genes x {Mk.shape[1]} compounds)")
check(got == want, f"the reported set EQUALS the independent one - {len(got)} pairs, none missed, none invented (missing {sorted(want - got)[:3]}, extra {sorted(got - want)[:3]})")
check(blk["n_sig_fdr"] == int((q_b <= 0.05).sum()), f"the number significant at the FDR agrees ({blk['n_sig_fdr']})")
for gene, met in [("CS", "Citrate"), ("IDH3A", "Isocitrate"), ("MDH2", "Malate")]:
    check((gene, met) in got, f"PLANTED: {gene} - {met} is reported")
check(abs(blk["expected_false"] - float(np.sum(q_b[(q_b <= 0.05) & (np.abs(rho_b) >= 0.6)]))) < 1e-6, "the expected-false count is the sum of the reported q-values")

print("=== the complete protein x metabolite block (proteins WITH missing values), rebuilt independently ===")
blk_pm = res["blocks"]["proteomics|metabolomics"]
Pk = P.loc[:, P.isna().mean() <= 0.2]
Pk = Pk.loc[:, Pk.std(ddof=0) > 0]
rho_pm, n_pm = spearman_matrix(Pk, Mk)
p_pm = pvalues(rho_pm, n_pm)
q_pm = bh(p_pm.ravel()).reshape(p_pm.shape)
want_pm = {(Pk.columns[i], Mk.columns[j]) for i, j in zip(*np.where((q_pm <= 0.05) & (np.abs(rho_pm) >= 0.6)))}
got_pm = {(a, b) for a, b in blk_pm["ids"]}
check(int((n_pm < 12).sum()) > 0, f"(the block really has pairs with fewer than 12 samples: {int((n_pm < 12).sum())})")
check(got_pm == want_pm and blk_pm["n_sig_fdr"] == int((q_pm <= 0.05).sum()),
      f"the reported set EQUALS the per-pair-ranked independent one ({len(got_pm)} pairs; missing {sorted(want_pm - got_pm)[:3]}, extra {sorted(got_pm - want_pm)[:3]})")

print("=== the other two blocks: every reported pair recomputed ===")
for key, (X, Y) in {"transcriptomics|proteomics": (T, P), "proteomics|metabolomics": (P, M)}.items():
    b = res["blocks"][key]
    bad = 0
    for (ia, ib), row in zip(b["ids"], b["table"]["rows"]):
        x, y = X[ia], Y[ib]
        rr, nn = spearman_matrix(pd.DataFrame({"x": x}), pd.DataFrame({"y": y}))
        bad += int(abs(rr[0, 0] - float(row["rho"])) > 5e-4 or int(nn[0, 0]) != int(row["n"]) or abs(float(row["rho"])) < 0.6 or float(row["FDR"]) > 0.05)
    check(bad == 0, f"{key}: all {len(b['ids'])} reported pairs match the raw files and pass both cutoffs")

print("=== downloads and the scatter endpoint ===")
for table, rows in [("concordance", conc), ("transcriptomics|metabolomics", blk["table"]["rows"])]:
    d = get("/api/multiomics/correlation/download", params={"table": table})
    df = pd.read_csv(io.StringIO(d.content.decode("utf-8-sig")))
    check(d.status_code == 200 and len(df) == len(rows) and "rho" in df.columns, f"download '{table}': {len(df)} rows = the table shown")
check(get("/api/multiomics/correlation/download", params={"table": "nope"}).status_code == 400, "an unknown table is refused")
sc = get("/api/multiomics/correlation/scatter", params={"kind": "gene", "gene": "CS"}).json()
check(abs(sc["rho"] - float(byg["CS"]["rho"])) < 5e-4 and sc["n"] == int(byg["CS"]["n"]) and len(sc["x"]) == 12 and sc["x_label"] == "Transcript CS",
      "the gene scatter returns the table's rho and n, 12 values, its own axis labels")
sp = get("/api/multiomics/correlation/scatter", params={"kind": "pair", "block": "transcriptomics|metabolomics", "a": "CS", "b": "Citrate"}).json()
rr, _ = spearman_matrix(pd.DataFrame({"x": T["CS"]}), pd.DataFrame({"y": M["Citrate"]}))
check(abs(sp["rho"] - rr[0, 0]) < 1e-6, f"the pair scatter's rho equals the independent one ({sp['rho']:.4f})")
check(get("/api/multiomics/correlation/scatter", params={"kind": "pair", "block": "transcriptomics|metabolomics", "a": "NOPE", "b": "Citrate"}).status_code == 404,
      "an unknown feature is a 404")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
