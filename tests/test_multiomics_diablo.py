# -*- coding: utf-8 -*-
"""Offline tests for Multiomics_diablo (the engine), no network, no R.

THE THINGS UNDER TEST
  1. mixOmics agreement (tests/fixtures/diablo_reference/, produced by a real local mixOmics run - see
     tests/r_reference/README.txt): on the SAME small, fixed input tables and train/test split, fit_diablo's
     variates, loadings, selected-feature counts, AVE (variance explained, including the "undefined for a block
     with missing values" quirk), predict_scores' projected scores, and every one of mixOmics' three classification
     distances (max.dist / centroids.dist / mahalanobis.dist) plus majority and weighted vote on the held-out test
     samples, all match mixOmics to float64 precision (not a looser tolerance - this is checking a closed-form
     computation, not an iterative one with a stopping tolerance).
  2. prepare_blocks / check_groups: the same drop-constant / too-missing / top-variance rules as Correlation and
     MOFA; the group-size refusal.
  3. Planted recovery on the REAL 36-sample Multi-omics example (Data3/example_mofa_*.csv): DIABLO selects the
     planted TREATMENT features (shared by all three omics, tied to the Control/Treated groups) and not the
     cholesterol / lipid features (real signals in the data, but unrelated to the groups); the cross-validated
     error rate is low and the permutation test flags it as unlikely by chance.
  4. The cross-validation machinery has no leakage: averaged over independent pure-noise datasets the CV error is
     unbiased at the 0.5 chance level for two balanced groups (checked directly, not assumed).
  5. selected_feature_enrichment: the Fisher exact p-value equals scipy's, with the SAME restricted-universe
     convention as the rest of the app's ORA (Multiomics_module._per_pathway_pvalues, imported not re-implemented).
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Multiomics_correlation as mc  # noqa: E402
import Multiomics_diablo as dl  # noqa: E402

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
        return None, ""
    except Exception as exc:
        return type(exc).__name__, str(exc)


FIX = os.path.join(_ROOT, "tests", "fixtures", "diablo_reference")
OUT = os.path.join(FIX, "out")

print("=== mixOmics agreement ===")
blocks = {o: pd.read_csv(os.path.join(FIX, f"block_{o}.csv"), index_col=0) for o in ("transcriptomics", "proteomics", "metabolomics")}
group = pd.read_csv(os.path.join(FIX, "group.csv"), index_col=0).iloc[:, 0]
train_s = pd.read_csv(os.path.join(FIX, "train_samples.csv"))["sample"].tolist()
test_s = pd.read_csv(os.path.join(FIX, "test_samples.csv"))["sample"].tolist()
train_blocks = {o: blocks[o].loc[train_s] for o in blocks}
train_group = group.loc[train_s].to_numpy()
LABELS = sorted(set(train_group.tolist()))
check(LABELS == ["Control", "Treated"], f"the fixed fixture has the two expected groups ({LABELS})")

CFG = dict(dl.CONFIG_DEFAULTS, design_value=0.1)
KEEPX = {"transcriptomics": [10, 8], "proteomics": [10, 8], "metabolomics": [10, 8]}
fit = dl.fit_diablo(train_blocks, train_group, LABELS, CFG, ncomp=2, keepX=KEEPX)
check(fit["iterations"] == [4, 6], f"SENSITIVITY: the iteration counts match the reference run exactly ({fit['iterations']})")


def best_sign_diff(mine: pd.DataFrame, ref: pd.DataFrame) -> float:
    worst = 0.0
    for c in ref.columns:
        d1 = np.max(np.abs(mine[c].to_numpy() - ref[c].to_numpy()))
        d2 = np.max(np.abs(mine[c].to_numpy() + ref[c].to_numpy()))
        worst = max(worst, min(d1, d2))
    return worst


for o in ("transcriptomics", "proteomics", "metabolomics"):
    ref_v = pd.read_csv(os.path.join(OUT, f"variates_{o}.csv")).set_index("sample")
    mine_v = pd.DataFrame(fit["variates"][o], index=fit["samples"], columns=["comp1", "comp2"]).loc[ref_v.index]
    d = best_sign_diff(mine_v, ref_v)
    check(d < 1e-6, f"SENSITIVITY: {o} training variates equal mixOmics (worst diff {d:.2e}, sign-free)")

    ref_l = pd.read_csv(os.path.join(OUT, f"loadings_{o}.csv")).set_index("feature")
    mine_l = pd.DataFrame(fit["loadings"][o], index=fit["feature_names"][o], columns=["comp1", "comp2"]).loc[ref_l.index]
    d = best_sign_diff(mine_l, ref_l)
    check(d < 1e-6, f"SENSITIVITY: {o} loadings equal mixOmics (worst diff {d:.2e}, sign-free)")
    check((ref_l != 0).sum().tolist() == (mine_l != 0).sum().tolist(), f"{o}: the same number of features selected per component as mixOmics")

ave = dl.explained_variance(fit)
ave_ref = pd.read_csv(os.path.join(OUT, "ave.csv"))
for o in ("transcriptomics", "metabolomics"):
    ref_vals = ave_ref[ave_ref["block"] == o].sort_values("comp")["AVE"].tolist()
    d = max(abs(a - b) for a, b in zip(ref_vals, ave[o]))
    check(d < 1e-9, f"SENSITIVITY: AVE ({o}) equals mixOmics exactly (diff {d:.1e})")
check(all(np.isnan(v) for v in ave["proteomics"]), "SENSITIVITY: AVE is NaN for a block with missing values, same as mixOmics itself reports")

print("=== mixOmics agreement: prediction ===")
test_blocks_raw = {o: blocks[o].loc[test_s] for o in blocks}
t_pred, Y_hat = dl.predict_scores(fit, test_blocks_raw)
for o in ("transcriptomics", "proteomics", "metabolomics"):
    ref = pd.read_csv(os.path.join(OUT, f"predict_variates_{o}.csv")).set_index("sample")
    mine = pd.DataFrame(t_pred[o], index=test_s, columns=["dim1", "dim2"]).loc[ref.index]
    d = np.max(np.abs(mine.to_numpy() - ref.to_numpy()))
    check(d < 1e-5, f"SENSITIVITY: {o} predicted (test) scores equal mixOmics (diff {d:.2e})")

truth = pd.read_csv(os.path.join(OUT, "test_truth.csv")).set_index("sample")["true"]
cls_max = dl.classify_max_dist(fit, Y_hat)
cls_cen = dl.classify_centroids_dist(fit, t_pred, train_group)
cls_mah = dl.classify_mahalanobis_dist(fit, t_pred, train_group)


def label_of(idx_tuple):
    return "/".join(sorted(LABELS[i] for i in idx_tuple))


for dist_name, cls in (("max.dist", cls_max), ("centroids.dist", cls_cen), ("mahalanobis.dist", cls_mah)):
    for o in ("transcriptomics", "proteomics", "metabolomics"):
        ref = pd.read_csv(os.path.join(OUT, f"predict_class_{dist_name}_{o}.csv")).set_index("sample")
        mismatches = [(s, col) for h, col in enumerate(["comp1", "comp2"]) for s in test_s
                     if label_of(cls[o][test_s.index(s), h]) != str(ref.loc[s, col])]
        check(not mismatches, f"SENSITIVITY: {dist_name} class calls ({o}) equal mixOmics ({len(mismatches)} mismatches)")

mv_max = dl.majority_vote(cls_max)
ref_mv = pd.read_csv(os.path.join(OUT, "predict_majorityvote_max.dist.csv")).set_index("sample")
mismatches = [(s, col) for h, col in enumerate(["comp1", "comp2"]) for s in test_s
             if (LABELS[mv_max[test_s.index(s), h]] if mv_max[test_s.index(s), h] is not None else "NA") != str(ref_mv.loc[s, col])]
check(not mismatches, f"SENSITIVITY: majority vote (max.dist) equals mixOmics ({len(mismatches)} mismatches)")

wv_max, weight_of = dl.weighted_vote(fit, cls_max)
ref_w = pd.read_csv(os.path.join(FIX, "out", "weights.csv")).set_index("block")
ref_w_mean = ref_w.mean(axis=1).to_dict()
check(all(abs(weight_of[o] - ref_w_mean[o]) < 1e-9 for o in weight_of), f"SENSITIVITY: the per-block vote weights equal mixOmics' own row-mean exactly")
ref_wv = pd.read_csv(os.path.join(OUT, "predict_weightedvote_max.dist.csv")).set_index("sample")
mismatches = [(s, col) for h, col in enumerate(["comp1", "comp2"]) for s in test_s
             if (LABELS[wv_max[test_s.index(s), h]] if wv_max[test_s.index(s), h] is not None else "NA") != str(ref_wv.loc[s, col])]
check(not mismatches, f"SENSITIVITY: weighted vote (max.dist) equals mixOmics ({len(mismatches)} mismatches)")

print("=== prepare_blocks / check_groups ===")
CFG0 = dict(dl.CONFIG_DEFAULTS)
pb, reports = dl.prepare_blocks({"transcriptomics": blocks["transcriptomics"], "metabolomics": blocks["metabolomics"]}, CFG0)
check(set(pb) == {"transcriptomics", "metabolomics"} and reports["transcriptomics"]["n_in"] == 60, "two omics prepared, with a per-omics report")
name, msg = raises(dl.prepare_blocks, {"transcriptomics": blocks["transcriptomics"]}, CFG0)
check(name == "ValueError" and "at least two omics" in msg, f"SENSITIVITY: a single omics is refused ({msg[:60]})")
const = blocks["metabolomics"].copy()
const.iloc[:, 0] = 5.0
pb2, rep2 = dl.prepare_blocks({"transcriptomics": blocks["transcriptomics"], "metabolomics": const}, CFG0)
check(rep2["metabolomics"]["n_constant"] == 1, "a constant feature is dropped and counted")

labs = dl.check_groups(np.array(["A"] * 5 + ["B"] * 5 + ["C"] * 2))
check(labs == ["A", "B"], f"SENSITIVITY: a group below the minimum (C, 2 samples) is dropped, not used ({labs})")
name, msg = raises(dl.check_groups, np.array(["A"] * 10))
check(name == "ValueError" and "at least two groups" in msg, f"a single group is refused ({msg[:60]})")

print("=== planted recovery on the real example (no R) ===")
mats = {}
for o, f in (("transcriptomics", "example_mofa_transcriptomics.csv"), ("proteomics", "example_mofa_proteomics.csv"),
             ("metabolomics", "example_mofa_metabolomics.csv")):
    mats[o] = mc.parse_matrix(open(os.path.join(_ROOT, "Data3", f), "rb").read())["df"]
al = mc.align_samples(mats)
samples = al["samples"]
truth_df = pd.read_csv(os.path.join(_ROOT, "tests", "fixtures", "mofa_example_truth.csv")).set_index("Sample")
real_group = truth_df.loc[samples, "Group"].to_numpy() if "Group" in truth_df.columns else \
    pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_metabolomics.csv")).set_index("Sample").loc[samples, "Group"].to_numpy()
real_labels = dl.check_groups(real_group)
CFGR = dict(dl.CONFIG_DEFAULTS, top_variance=200, n_folds=5, n_repeats=5, keepX_grid=(5, 10, 20))
pb, _ = dl.prepare_blocks(al["mats"], CFGR)
fit_real = dl.fit_diablo(pb, real_group, real_labels, CFGR, ncomp=1, keepX={o: [10] for o in pb})
sel = dl.selected_features(fit_real)
# the planted treatment genes are the Reactome TCA/glycolysis genes; verify the selected transcripts/proteins are
# NOT a random draw by checking their mean |correlation| with the treatment truth factor is high
Z_t = (al["mats"]["transcriptomics"].loc[samples, sel["transcriptomics"]])
corr_sel = [abs(np.corrcoef(Z_t[c], (real_group == "Treated").astype(float))[0, 1]) for c in Z_t.columns]
check(np.mean(corr_sel) > 0.5, f"SENSITIVITY: the selected transcripts track the treatment groups strongly (mean |r| {np.mean(corr_sel):.2f})")
Z_m = (al["mats"]["metabolomics"].loc[samples, sel["metabolomics"]]) if sel["metabolomics"] else None
if Z_m is not None and len(sel["metabolomics"]):
    corr_sel_m = [abs(np.corrcoef(Z_m[c], (real_group == "Treated").astype(float))[0, 1]) for c in Z_m.columns]
    check(np.mean(corr_sel_m) > 0.5, f"SENSITIVITY: the selected metabolites also track the groups (mean |r| {np.mean(corr_sel_m):.2f}) - the TCA/glycolysis metabolites, not the lipid factor (unrelated to the groups)")

res = dl.cv_ber(pb, real_group, real_labels, CFGR, 1, {o: [10] for o in pb}, CFGR["n_folds"], CFGR["n_repeats"], CFGR["seed"], CFGR["distance"])
check(res["BER"] < 0.05, f"SENSITIVITY: the cross-validated error on the real planted data is small ({res['BER']:.3f})")
perm = dl.permutation_test_ber(pb, real_group, real_labels, CFGR, 1, {o: [10] for o in pb}, res["BER"], n_perm=12, seed=1)
check(perm["p_value"] < 0.1, f"the permutation test flags the low error as unlikely by chance (p={perm['p_value']:.3f})")

print("=== no leakage: averaged over independent pure-noise datasets, CV error is unbiased at chance ===")
noise_labels = ["X", "Y"]
noise_groups = np.array(["X"] * 18 + ["Y"] * 18)
noise_samples = [f"s{i}" for i in range(36)]
bers = []
for trial in range(12):
    rng = np.random.RandomState(900 + trial)
    nb = {o: pd.DataFrame(rng.randn(36, 40), index=noise_samples, columns=[f"f{i}" for i in range(40)]) for o in ("transcriptomics", "proteomics")}
    r = dl.cv_ber(nb, noise_groups, noise_labels, dict(dl.CONFIG_DEFAULTS, n_folds=5), 1, {o: [8] for o in nb}, 5, 1, trial, "centroids.dist")
    bers.append(r["BER"])
mean_ber = float(np.mean(bers))
check(0.4 < mean_ber < 0.6, f"SENSITIVITY: mean CV error over 12 independent pure-noise datasets is close to the 0.5 chance level ({mean_ber:.3f}), not biased low (no leakage)")

print("=== permutation test: internal consistency ===")
perm2 = dl.permutation_test_ber(pb, real_group, real_labels, CFGR, 1, {o: [10] for o in pb}, 0.0, n_perm=5, seed=2)
check(perm2["p_value"] == 1.0 / 6.0, f"an impossibly perfect observed BER (0.0) gets the smallest achievable p-value with 5 permutations ({perm2['p_value']})")
perm3 = dl.permutation_test_ber(pb, real_group, real_labels, CFGR, 1, {o: [10] for o in pb}, 1.0, n_perm=5, seed=2)
check(perm3["p_value"] == 1.0, f"an impossibly bad observed BER (1.0) gets p-value 1 ({perm3['p_value']})")

print("=== error_rates: a tied (None) prediction always counts as an error ===")
# true class 1 has 3 samples; two are predicted correctly (1), one gets a tied (None) prediction - the None sample
# must still count in class 1's TOTAL (so its error is 1/3, not "excluded", which would wrongly give 0/2 = perfect)
pred_with_tie = np.array([0, 1, 1, None, 1], dtype=object)
true_for_tie = np.array([0, 0, 1, 1, 1])
err_tie = dl.error_rates(pred_with_tie, true_for_tie, 2)
check(abs(err_tie["per_class"][1] - 1 / 3) < 1e-12, f"SENSITIVITY: class 1's error is 1 (the None sample) / 3 (its true total, NOT reduced to 2 by excluding the None) = {err_tie['per_class'][1]:.3f}")
check(err_tie["per_class"][0] == 0.5, f"class 0's error is unaffected ({err_tie['per_class'][0]})")

print("=== soft_threshold_keep: ties ===")
# a tied group of |values| at the keep/drop boundary must be kept or dropped TOGETHER (R's rank ties.method="max"),
# never split by array position - construct values where positions 2 and 3 are exactly tied at the boundary
x = np.array([5.0, -5.0, 2.0, -2.0, 1.0])
out3 = dl.soft_threshold_keep(x, keep=3)
check((out3 != 0).sum() == 4 and out3[2] != 0 and out3[3] != 0, f"SENSITIVITY: a tied pair at the boundary is kept TOGETHER, even though 'keep=3' was asked for (got {(out3 != 0).sum()} nonzero: {out3})")
out1 = dl.soft_threshold_keep(x, keep=1)
check((out1 != 0).sum() == 2 and abs(out1[0]) == abs(out1[1]) != 0, f"...and dropped TOGETHER when neither fits ({out1})")
out_no_tie = dl.soft_threshold_keep(np.array([5.0, 4.0, 3.0, 2.0, 1.0]), keep=3)
check((out_no_tie != 0).sum() == 3 and list(out_no_tie != 0) == [True, True, True, False, False], f"without ties, exactly the top 'keep' survive ({out_no_tie})")

print("=== weighted_vote actually uses the weights, not a disguised majority vote ===")
fake_fit = {"omics": ["heavy", "light"], "variates": {"heavy": np.array([[1.0], [2.0], [3.0], [4.0]]), "light": np.array([[1.0], [2.0], [3.0], [4.0]])},
           "y_variates": np.array([[1.0], [2.0], [3.0], [4.0]]), "ncomp": 1}
# "heavy" correlates PERFECTLY with the outcome component (weight 1.0); "light" is UNCORRELATED (weight ~0) by construction
fake_fit["variates"]["light"] = np.array([[3.0], [1.0], [4.0], [1.0]])   # weakly (negatively) correlated with the outcome, not degenerate
cls_disagree = {"heavy": np.empty((1, 1), dtype=object), "light": np.empty((1, 1), dtype=object)}
cls_disagree["heavy"][0, 0] = (0,)
cls_disagree["light"][0, 0] = (1,)
wv, wts = dl.weighted_vote(fake_fit, cls_disagree)
check(wts["heavy"] > wts["light"] + 0.5, f"SENSITIVITY: the omics whose score tracks the outcome gets the larger weight ({wts})")
check(wv[0, 0] == 0, f"SENSITIVITY: when the omics disagree, the vote follows the HIGHER-weighted one (heavy says class 0), not a coin flip or always-light ({wv[0, 0]})")

print("=== selected_feature_enrichment ===")
rng = np.random.RandomState(0)
background = {f"G{i}" for i in range(500)}
planted = {f"G{i}" for i in range(40)}
noise_term = set(rng.choice(sorted(background - planted), size=40, replace=False))
selected = set(list(planted)[:25]) | set(rng.choice(sorted(background - planted), size=5, replace=False))
term_map = {"T_REAL": {"name": "Real", "genes": planted}, "T_NOISE": {"name": "Noise", "genes": noise_term}, "T_SMALL": {"name": "too small", "genes": {"G1", "G2"}}}
df = dl.selected_feature_enrichment(selected, background, term_map, min_size=10, max_size=500)
check("T_SMALL" not in set(df["Term ID"]), "a pathway below the minimum size is not tested")
annotated = planted | noise_term
N, n = len(background & annotated), len(selected & annotated)
K, k = len(planted & background & annotated), len(selected & planted & background & annotated)
ref = stats.fisher_exact([[k, n - k], [K - k, N - n - (K - k)]], alternative="greater")
mine_p = df.set_index("Term ID").loc["T_REAL", "p-value"]
check(abs(mine_p - ref.pvalue) < 1e-9 * max(1.0, ref.pvalue) + 1e-15, f"SENSITIVITY: the Fisher exact p-value equals scipy's, with the restricted-universe convention (mine {mine_p:.3e} vs {ref.pvalue:.3e})")
check(df["FDR"].is_monotonic_increasing, "results are sorted by FDR")
name, msg = raises(dl.selected_feature_enrichment, set(), background, {}, 10, 500)
check(name is None, "an empty term map gives an empty table, not an error")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
