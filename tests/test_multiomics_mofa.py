# -*- coding: utf-8 -*-
"""Offline tests of Multiomics_mofa.py (Multi-omics MOFA / Factor-based Integration). Needs mofapy2 and a numpy with a working BLAS
(the project environment, see requirements.txt).

Each number is checked against an INDEPENDENT reference or the planted truth:
  - recovery of PLANTED factors (a factor shared by three omics, one by two, one specific to a single omics, and noise factors that
    must be dropped), the R2 pattern (shared factor explains variance in every omics, specific ones only in their own);
  - the variance explained equals an exact re-implementation of mofapy2's own definition computed from the model's Z, W and the data it saw;
  - factor order (by total variance explained), sign (largest weight positive), determinism for one seed;
  - stability: a real factor reappears across seeds (~1), a made-up one does not and is flagged;
  - group tests against scipy and a by-hand Kruskal-Wallis, BH against Genomics_module's own;
  - the weight-set enrichment against scipy's Mann-Whitney per pathway (ties included), BH, the size gates, the keyless features left out,
    a planted pathway found in the right direction, and a null where p-values are uniform.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm  # noqa: E402
import Multiomics_mofa as mf  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def raises(fn, *a, **k):
    try:
        fn(*a, **k)
    except Exception as exc:
        return type(exc).__name__, str(exc)
    return None, ""


# ------------------------------------------------------------------ planted data
rng = np.random.RandomState(1)
N = 48
samples = [f"S{i}" for i in range(N)]
z_all = rng.randn(N)          # shared by the three omics (the "treatment")
z_tp = rng.randn(N)           # transcripts and proteins only
z_m = rng.randn(N)            # metabolites only


def sparse(p, frac):
    return rng.randn(p) * (rng.rand(p) < frac)


pT, pP, pM = 120, 90, 40
wT_all, wT_tp = sparse(pT, 0.35), sparse(pT, 0.3)
wP_all, wP_tp = sparse(pP, 0.35), sparse(pP, 0.3)
wM_all, wM_m = sparse(pM, 0.5), sparse(pM, 0.5)
T = np.outer(z_all, wT_all) + np.outer(z_tp, wT_tp) + 0.6 * rng.randn(N, pT)
P = np.outer(z_all, wP_all) + np.outer(z_tp, wP_tp) + 0.6 * rng.randn(N, pP)
M = np.outer(z_all, wM_all) + np.outer(z_m, wM_m) + 0.6 * rng.randn(N, pM)
P[rng.rand(N, pP) < 0.05] = np.nan
MATS = {"transcriptomics": pd.DataFrame(T, index=samples, columns=[f"t{i}" for i in range(pT)]),
        "proteomics": pd.DataFrame(P, index=samples, columns=[f"p{i}" for i in range(pP)]),
        "metabolomics": pd.DataFrame(M, index=samples, columns=[f"m{i}" for i in range(pM)])}
CFG = {"factors": 8, "convergence": "fast", "top_variance": 2000, "n_refits": 2, "seed": 7}

print("=== preparing the views ===")
raw = MATS["proteomics"].copy()
raw["const"] = 3.0
raw["gappy"] = np.where(np.arange(N) < 30, np.nan, 1.0 + np.arange(N))        # 62% missing
views, reps = mf.prepare_views({"proteomics": raw, "metabolomics": MATS["metabolomics"]}, mf.CONFIG_DEFAULTS)
check("const" not in views["proteomics"].columns and "gappy" not in views["proteomics"].columns and reps["proteomics"]["n_constant"] == 1
      and reps["proteomics"]["n_too_many_missing"] == 1, "a constant feature and one missing in 62% of the samples are dropped and counted")
v2, _ = mf.prepare_views(MATS, dict(mf.CONFIG_DEFAULTS, top_variance=50))
check(all(v2[o].shape[1] == min(50, MATS[o].shape[1]) for o in v2) and v2["transcriptomics"].shape[1] == 50 and v2["metabolomics"].shape[1] == 40,
      "the most variable 50 features are kept per omics (an omics with fewer keeps them all)")
top_var = MATS["transcriptomics"].var(axis=0, ddof=0).sort_values(ascending=False).index[:50]
check(set(v2["transcriptomics"].columns) == set(top_var), "SENSITIVITY: they are the 50 with the largest variance")
name_e, msg_e = raises(mf.prepare_views, {"proteomics": pd.DataFrame({"a": [1.0, 1.0, 1.0, 1.0], "b": [2.0, 2.0, 2.0, 2.0]})}, mf.CONFIG_DEFAULTS)
check(name_e == "ValueError" and "fewer than two usable features" in msg_e, "an omics with nothing usable left is refused with the reason")
name_e, msg_e = raises(mf.fit_mofa, {"metabolomics": MATS["metabolomics"].iloc[:4]}, mf.CONFIG_DEFAULTS)
check(name_e == "ValueError" and "at least 6 samples" in msg_e, f"fewer than {mf.MIN_SAMPLES} samples is refused ({msg_e[:50]})")
check(raises(mf.fit_mofa, MATS, dict(CFG, factors=1))[0] == "ValueError" and raises(mf.fit_mofa, MATS, dict(CFG, convergence="x"))[0] == "ValueError", "bad settings are refused")

print("=== fitting: the planted factors come back ===")
res = mf.run_mofa(MATS, CFG, groups=None)
fit, Z = res["fit"], res["factors"]
K = Z.shape[1]
check(2 <= K <= 5, f"MOFA keeps {K} factors from 8 offered: the unused ones are dropped, not kept as noise")
rec = {}
for name, truth in (("shared", z_all), ("transcript+protein", z_tp), ("metabolite", z_m)):
    r = [abs(np.corrcoef(Z.iloc[:, k], truth)[0, 1]) for k in range(K)]
    rec[name] = (max(r), int(np.argmax(r)))
check(all(v[0] > 0.9 for v in rec.values()), f"SENSITIVITY: every planted factor is recovered, |r| with the truth {[round(v[0], 3) for v in rec.values()]}")
check(len({v[1] for v in rec.values()}) == 3, "...each by a different factor")
r2 = res["r2"]
fs, ftp, fm = (f"Factor{rec[n][1] + 1}" for n in ("shared", "transcript+protein", "metabolite"))
check(all(r2.loc[OL, fs] > 5 for OL in r2.index), f"the SHARED factor explains > 5% in every omics ({r2[fs].round(1).to_dict()})")
check(r2.loc["Metabolomics", ftp] < 2 and r2.loc["Transcriptomics", ftp] > 5 and r2.loc["Proteomics", ftp] > 5,
      f"the transcript+protein factor explains variance in those two omics only ({r2[ftp].round(1).to_dict()})")
check(r2.loc["Transcriptomics", fm] < 2 and r2.loc["Proteomics", fm] < 2 and r2.loc["Metabolomics", fm] > 5, f"the metabolite factor is specific to metabolites ({r2[fm].round(1).to_dict()})")

print("=== order, sign, determinism ===")
tot = r2.sum(axis=0).to_numpy()
check(all(tot[i] >= tot[i + 1] - 1e-12 for i in range(len(tot) - 1)), f"factors are ordered by total variance explained ({tot.round(1).tolist()})")
ok_sign = True
for j in range(K):
    best = max(((abs(fit["W"][o][:, j]).max(), fit["W"][o][int(np.argmax(abs(fit["W"][o][:, j]))), j]) for o in fit["views"]), key=lambda t: t[0])
    ok_sign &= best[1] > 0
check(ok_sign, "SENSITIVITY: every factor's largest weight (over all omics) is positive - the arbitrary sign is fixed")
res2 = mf.run_mofa(MATS, CFG, groups=None)
check(np.allclose(res["factors"].to_numpy(), res2["factors"].to_numpy(), atol=1e-9) and all(np.allclose(res["fit"]["W"][o], res2["fit"]["W"][o], atol=1e-9) for o in res["omics"]),
      "the same data and seed give the same factors and weights")
res3 = mf.run_mofa(MATS, dict(CFG, seed=8), groups=None)
rec3 = [max(abs(np.corrcoef(res3["factors"].iloc[:, k], t)[0, 1]) for k in range(res3["factors"].shape[1])) for t in (z_all, z_tp, z_m)]
check(min(rec3) > 0.9, f"another seed finds the same planted factors (|r| {[round(float(x), 3) for x in rec3]}) - the result does not hinge on the seed")

print("=== variance explained = mofapy2's own definition, recomputed ===")
Zf = res["fit"]["Z"]
worst = 0.0
for i, o in enumerate(fit["views"]):
    Y = fit["Y"][o]
    miss = MATS[o][list(res["views"][o].columns)].isna().to_numpy()
    SS = float((Y ** 2).sum())
    for k in range(K):
        pred = np.outer(Zf[:, k], fit["W"][o][:, k])
        pred[miss] = 0.0
        worst = max(worst, abs(100 * (1 - float(((Y - pred) ** 2).sum()) / SS) - r2.iloc[i, k]))
    allp = Zf @ fit["W"][o].T
    allp[miss] = 0.0
    worst = max(worst, abs(100 * (1 - float(((Y - allp) ** 2).sum()) / SS) - res["r2_total"][o]))
check(worst < 1e-8, f"R2 per factor and in total equal the exact formula from Z, W and the data the model saw (largest difference {worst:.1e} points)")
check(all(0 < res["r2_total"][o] < 100 for o in res["omics"]), f"total variance explained per omics: { {o: round(v, 1) for o, v in res['r2_total'].items()} }")
obs_only = {}
for o in res["omics"]:
    y = res["views"][o].to_numpy(float)
    y = (y - np.nanmean(y, axis=0))
    y = y / np.sqrt(np.nanvar(y))                       # unit total variance, observed cells only
    ok = ~np.isnan(y)
    pred = np.where(ok, Zf @ fit["W"][o].T, 0.0)
    obs_only[o] = 100 * (1 - float(((np.where(ok, y, 0) - pred) ** 2).sum()) / float((y[ok] ** 2).sum()))
check(all(abs(obs_only[o] - res["r2_total"][o]) < 6 for o in obs_only), f"(sanity) from the observed cells alone the totals agree within a few points ({ {o: round(v, 1) for o, v in obs_only.items()} })")

print("=== stability ===")
st = res["stability"]
check(st["n_refits"] == 2 and all(p["stability"] is not None for p in st["per_factor"]), "the model was refitted with 2 other seeds")
sdict = {p["factor"]: p["stability"] for p in st["per_factor"]}
check(all(sdict[f] > 0.9 for f in (fs, ftp, fm)), f"the real factors reappear in every refit ({ {k: round(sdict[k], 3) for k in (fs, ftp, fm)} })")
fake = {"Z": np.random.RandomState(99).randn(N, 3)}
st_fake = mf.stability(res["views"], dict(CFG, n_refits=2), fake)
check(all(p["unstable"] for p in st_fake["per_factor"]) and max(p["stability"] for p in st_fake["per_factor"]) < 0.7,
      f"SENSITIVITY: factors that are not in the data are flagged unstable (stability {[round(p['stability'], 2) for p in st_fake['per_factor']]})")
_real_fit = mf.fit_mofa


def _fit_that_fails_on_refits(views, cfg=None, seed=None):
    if seed is not None:
        raise ValueError("no factor survived")
    return _real_fit(views, cfg, seed)


mf.fit_mofa = _fit_that_fails_on_refits
try:
    st_fail = mf.stability(res["views"], dict(CFG, n_refits=2), res["fit"])
finally:
    mf.fit_mofa = _real_fit
check(st_fail["n_refits"] == 2 and all(p["unstable"] and p["stability"] == 0.0 for p in st_fail["per_factor"]) and [f["n_factors"] for f in st_fail["fits"]] == [0, 0],
      "SENSITIVITY: a refit that finds no factor at all does not abort the run - the reference factors are then flagged unstable")
no_refit = mf.run_mofa(MATS, dict(CFG, n_refits=0))
check(no_refit["stability"]["n_refits"] == 0 and all(p["stability"] is None and not p["unstable"] for p in no_refit["stability"]["per_factor"]), "n_refits = 0 skips the check and says nothing about it")

print("=== warnings ===")
small = {o: m.iloc[:12] for o, m in MATS.items()}
res_small = mf.run_mofa(small, dict(CFG, n_refits=0))
check(any("Only 12 samples" in w for w in res_small["warnings"]), "fewer than 15 samples warns that the weights are unreliable")
check(not any("Only" in w for w in res["warnings"]), "48 samples gives no sample-size warning")
res_one = mf.run_mofa({"transcriptomics": MATS["transcriptomics"]}, dict(CFG, n_refits=0))
check(any("Only one omics" in w for w in res_one["warnings"]), "a single omics is allowed but says there is nothing shared")

print("=== group tests ===")
grp = np.array(["A"] * 24 + ["B"] * 24)
Zt = np.column_stack([np.where(grp == "A", 1.0, -1.0) + 0.3 * rng.randn(N), rng.randn(N), rng.randn(N)])
gt = mf.group_tests(Zt, grp)
ref_p = [stats.mannwhitneyu(Zt[grp == "A", k], Zt[grp == "B", k], alternative="two-sided").pvalue for k in range(3)]
check(gt["test"] == "Mann-Whitney U" and np.allclose([r["p-value"] for r in gt["rows"]], ref_p, rtol=1e-9), "two groups: Mann-Whitney p equals scipy's")
check(np.allclose([r["FDR"] for r in gt["rows"]], gm.benjamini_hochberg(ref_p), rtol=1e-9), "FDR equals Genomics_module's Benjamini-Hochberg")
u = float(sum((a > b) + 0.5 * (a == b) for a in Zt[grp == "A", 0] for b in Zt[grp == "B", 0]))
check(abs(gt["rows"][0]["Effect"] - (2 * u / (24 * 24) - 1)) < 1e-12 and gt["rows"][0]["Effect"] > 0.9, f"effect = rank-biserial correlation counted by hand ({gt['rows'][0]['Effect']:.3f}); A is higher on factor 1")
check(gt["rows"][0]["FDR"] < 1e-6 and gt["rows"][1]["FDR"] > 0.05, "the planted group factor is significant, the noise factors are not")
grp3 = np.array(["A"] * 16 + ["B"] * 16 + ["C"] * 16)
Z3 = np.column_stack([np.repeat([0.0, 1.0, 3.0], 16) + 0.5 * rng.randn(N), rng.randn(N)])
g3 = mf.group_tests(Z3, grp3)
allz = Z3[:, 0]
rk = stats.rankdata(allz)
H = 12.0 / (N * (N + 1)) * sum(rk[grp3 == x].sum() ** 2 / 16 for x in "ABC") - 3 * (N + 1)
tie = 1 - float((np.unique(allz, return_counts=True)[1] ** 3 - np.unique(allz, return_counts=True)[1]).sum()) / (N ** 3 - N)
check(g3["test"] == "Kruskal-Wallis" and abs(stats.chi2.sf(H / tie, 2) - g3["rows"][0]["p-value"]) < 1e-12, "three groups: Kruskal-Wallis p equals the by-hand formula")
miss = grp.astype(object).copy()
miss[:4] = ""
gm_ = mf.group_tests(Zt, miss)
check(gm_["counts"] == {"A": 20, "B": 24}, "samples without a group are left out")
check(mf.group_tests(Zt, np.array(["A"] * N))["rows"] == [] and "Fewer than two groups" in mf.group_tests(Zt, np.array(["A"] * N))["note"], "one group cannot be tested, and the note says so")
check(mf.group_tests(Zt, np.array(["A"] * 47 + ["B"]))["rows"] == [], "a group of one sample is dropped")
res_g = mf.run_mofa(MATS, dict(CFG, n_refits=0), groups={s: ("treated" if z_all[i] > 0 else "control") for i, s in enumerate(samples)})
best_row = min(res_g["group_tests"]["rows"], key=lambda r: r["p-value"])
check(best_row["Factor"] == fs and best_row["FDR"] < 1e-6, f"through the whole pipeline: the shared factor ({fs}) is the one separating the groups (FDR {best_row['FDR']:.1e})")

print("=== weight-set enrichment ===")
rg = np.random.RandomState(5)
G = 600
keys = [f"G{i}" for i in range(G)]
w = rg.randn(G)
w[:40] += 2.0                                                    # planted: the first 40 genes have high weights
w = np.round(w, 1)                                               # lots of ties
terms = {"T_planted": {"name": "planted", "genes": {f"G{i}" for i in range(0, 40)}},
         "T_neg": {"name": "planted low", "genes": {f"G{i}" for i in range(500, 540)}},
         "T_small": {"name": "too small", "genes": {"G1", "G2", "G3"}}}
w[500:540] -= 2.0
for t in range(60):
    terms[f"T_rand{t}"] = {"name": f"random {t}", "genes": {keys[i] for i in rg.choice(np.arange(40, 500), size=rg.randint(12, 80), replace=False)}}
en = mf.weight_set_enrichment(w, keys, terms, min_size=10, max_size=500)
check("T_small" not in set(en["Term ID"]), "a pathway below the minimum size is not tested")
ref = {}
for tid, info in terms.items():
    mem = np.array([k in info["genes"] for k in keys])
    if mem.sum() < 10:
        continue
    ref[tid] = stats.mannwhitneyu(w[mem], w[~mem], alternative="two-sided", method="asymptotic", use_continuity=True)
got = en.set_index("Term ID")
check(set(got.index) == set(ref) and all(abs(got.loc[t, "p-value"] - ref[t].pvalue) < 1e-9 * max(1, ref[t].pvalue) + 1e-15 for t in ref),
      f"SENSITIVITY: all {len(ref)} p-values equal scipy's asymptotic Mann-Whitney, with ties and continuity correction")
check(all(abs(got.loc[t, "AUC"] - ref[t].statistic / (got.loc[t, "Members"] * (G - got.loc[t, "Members"]))) < 1e-12 for t in ref), "AUC = U / (n1 n2)")
check(np.allclose(got["FDR"].to_numpy(), gm.benjamini_hochberg(got["p-value"].to_numpy()), rtol=1e-9), "FDR is Benjamini-Hochberg over the tested pathways")
check(got.loc["T_planted", "Direction"] == "higher weights" and got.loc["T_planted", "FDR"] < 1e-10 and got.loc["T_neg", "Direction"] == "lower weights" and got.loc["T_neg", "FDR"] < 1e-10,
      "the planted high-weight and low-weight pathways are found, in the right direction")
check(en.iloc[0]["FDR"] <= en.iloc[-1]["FDR"] and en["FDR"].is_monotonic_increasing, "results are sorted by FDR")
en_abs = mf.weight_set_enrichment(w, keys, terms, min_size=10, mode="absolute").set_index("Term ID")
refa = stats.mannwhitneyu(np.abs(w[:40]), np.abs(w[40:]), alternative="greater", method="asymptotic", use_continuity=True)
check(abs(en_abs.loc["T_planted", "p-value"] - refa.pvalue) < 1e-9 * max(1, refa.pvalue) + 1e-15 and en_abs.loc["T_planted", "Direction"] == "larger |weight|", "'absolute' mode: one-sided test on |weight| equals scipy's")
keys_gap = [k if i % 3 else None for i, k in enumerate(keys)]
en_gap = mf.weight_set_enrichment(w, keys_gap, terms, min_size=5)
members_expected = sum(1 for i in range(40) if i % 3)          # the planted pathway's genes that still have an identifier
check(en_gap.set_index("Term ID").loc["T_planted", "Members"] == members_expected == 26,
      f"SENSITIVITY: features without an identifier are left out of the test (the planted pathway has {members_expected} testable members, not 40)")
check(mf.weight_set_enrichment(w[:15], keys[:15], terms, min_size=10).empty, "too few features to test gives an empty table, not an error")
check(raises(mf.weight_set_enrichment, w, keys, terms, 10, 500, "x")[0] == "ValueError", "an unknown mode is refused")
pnull = []
for rep in range(8):
    wn = np.random.RandomState(100 + rep).randn(G)
    tm = {f"N{t}": {"name": "n", "genes": {keys[i] for i in np.random.RandomState(rep * 1000 + t).choice(G, size=40, replace=False)}} for t in range(120)}
    pnull += mf.weight_set_enrichment(wn, keys, tm, 10)["p-value"].tolist()
ks = stats.kstest(pnull, "uniform").statistic
check(ks < 0.04 and 0.03 < np.mean(np.array(pnull) < 0.05) < 0.07, f"NULL CALIBRATION: with random weights and random pathways the p-values are uniform (KS {ks:.3f}; {np.mean(np.array(pnull) < 0.05):.1%} below 0.05 over {len(pnull)})")

print("=== several features for one gene ===")
cw, ck = mf.collapse_by_key(np.array([0.5, -2.0, 1.0, 0.3, 4.0, -0.1]), ["A", "A", "B", None, "B", ""])
check(ck == ["A", None, "B", ""], f"each gene counts once; features without a key stay, in their original order ({ck})")
check(np.allclose(cw, [-2.0, 0.3, 4.0, -0.1]), f"the feature with the largest |weight| represents its gene, sign kept ({cw.tolist()})")
cw2, ck2 = mf.collapse_by_key(np.array([1.0, -1.0]), ["A", "A"])
check(ck2 == ["A"] and len(cw2) == 1 and abs(cw2[0]) == 1.0, "a tie keeps the first feature, one entry")
cw3, ck3 = mf.collapse_by_key(np.array([3.0, 1.0, 2.0]), ["X", "Y", "Z"])
check(np.allclose(cw3, [3.0, 1.0, 2.0]) and ck3 == ["X", "Y", "Z"], "all-distinct keys are unchanged, order kept")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
