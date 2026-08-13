"""MEASUREMENT HARNESS, not a test suite. Decides whether empirical-null
calibration is worth building for TF activity inference, before any engine
code is written (see plan: 'measurement-first', with a kill criterion).

Reuses tr._ulm_score directly (already tested/verified) rather than
reimplementing the statistic - only the CALIBRATION logic here is throwaway/
standalone, deliberately, so this script proves the idea before it becomes
production code.

Three scenarios + a real-data check:
  A. Independent, no signal   - sanity check the harness itself is unbiased.
  B. Co-regulated, no signal  - THE measurement that decides this plan.
  C. Real signal present      - power check: calibration must not erase it.
  D. lambda-hat on the user's actual example dataset.
"""
import os
import sys

import numpy as np
from scipy.special import ndtr, stdtr

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx
import Transcriptomics_module_retriever as retr
import Transcriptomics_regulon as tr

ALPHA = 0.05


def robust_lambda(scores):
    """Standalone (not-yet-production) empirical-null estimate: median centre,
    MAD-based scale, clamped to never go below 1.0 (theory)."""
    scores = np.asarray(scores, dtype=float)
    center = float(np.median(scores))
    mad = float(np.median(np.abs(scores - center)))
    sigma_hat = mad / 0.6745 if mad > 0 else 1.0
    lam = max(1.0, sigma_hat)
    return center, lam


def calibrated_pvalues(scores, center, lam):
    z = (np.asarray(scores, dtype=float) - center) / lam
    return 2.0 * ndtr(-np.abs(z))


def make_regulons(rng, n_genes, n_tfs, size_lo=5, size_hi=100, pos_frac_range=(0.3, 0.7), pool=None):
    """pool: optional array of gene indices to draw targets from (for scenario
    B, restricts a TF's targets to ONE correlated block, matching how a real
    TF's targets can cluster inside one co-regulated module)."""
    regulons = []
    universe = np.arange(n_genes) if pool is None else pool
    for _ in range(n_tfs):
        size = rng.integers(size_lo, size_hi + 1)
        size = min(size, len(universe))
        targets = rng.choice(universe, size=size, replace=False)
        pos_frac = rng.uniform(*pos_frac_range)
        n_pos = max(1, int(round(size * pos_frac)))
        n_pos = min(n_pos, size - 1) if size > 1 else size
        pos = targets[:n_pos].tolist()
        neg = targets[n_pos:].tolist()
        regulons.append((pos, neg))
    return regulons


def run_scenario(name, y, regulons, known_active_idx=None):
    n = len(y)
    yc = y - y.mean()
    syy = float(np.dot(yc, yc))
    scores, pvals = [], []
    for pos, neg in regulons:
        score, p, n_used, beta1 = tr._ulm_score(pos, neg, yc, syy, n)
        if score is None:
            continue
        scores.append(score)
        pvals.append(p)
    scores = np.array(scores)
    pvals = np.array(pvals)

    fpr_raw = float(np.mean(pvals < ALPHA))
    center, lam = robust_lambda(scores)
    cal_pvals = calibrated_pvalues(scores, center, lam)
    fpr_cal = float(np.mean(cal_pvals < ALPHA))

    print(f"\n=== Scenario {name} ===")
    print(f"  n_tfs_tested = {len(scores)}")
    print(f"  raw parametric FPR at alpha=0.05   = {fpr_raw:.4f}  (nominal target: {ALPHA})")
    print(f"  median score = {center:.4f}   MAD-based sigma_hat = {lam:.4f}   lambda_used = {lam:.4f}")
    print(f"  calibrated FPR at alpha=0.05        = {fpr_cal:.4f}")

    if known_active_idx is not None:
        active_scores = scores[known_active_idx]
        active_p_raw = pvals[known_active_idx]
        active_p_cal = cal_pvals[known_active_idx]
        power_raw = float(np.mean(active_p_raw < ALPHA))
        power_cal = float(np.mean(active_p_cal < ALPHA))
        print(f"  KNOWN-ACTIVE TFs: n={len(known_active_idx)}  power_raw={power_raw:.3f}  power_calibrated={power_cal:.3f}")

    return {"fpr_raw": fpr_raw, "fpr_cal": fpr_cal, "lambda": lam, "n": len(scores)}


rng = np.random.default_rng(2026)
N_GENES = 3000
N_TFS = 400

# --- Scenario A: independent, no signal --------------------------------
y_a = rng.normal(size=N_GENES)
regulons_a = make_regulons(rng, N_GENES, N_TFS)
res_a = run_scenario("A (independent, no signal)", y_a, regulons_a)

# --- Scenario B: co-regulated, no signal --------------------------------
# 15 blocks of 60 genes, each sharing a per-block random common factor
# (simulates real co-regulation: genes driven by the same upstream event end
# up with correlated values in a single-summary-statistic dataset). NO TF
# has a true effect - "no signal" - but each TF's targets are drawn entirely
# from ONE block, so genuinely co-move for reasons outside that TF.
n_blocks, block_size = 15, 60
y_b = rng.normal(size=N_GENES)
block_bounds = []
for b in range(n_blocks):
    lo = b * block_size
    hi = lo + block_size
    shift = rng.normal(scale=2.0)
    y_b[lo:hi] += shift + rng.normal(scale=0.4, size=block_size)
    block_bounds.append((lo, hi))

regulons_b = []
for _ in range(N_TFS):
    lo, hi = block_bounds[rng.integers(0, n_blocks)]
    pool = np.arange(lo, hi)
    regulons_b.extend(make_regulons(rng, N_GENES, 1, size_lo=5, size_hi=min(40, block_size), pool=pool))
res_b = run_scenario("B (co-regulated targets, no true signal)", y_b, regulons_b)

# --- Scenario C: real signal present, still co-regulated background -----
y_c = y_b.copy()
n_active = 30
active_regulons = []
active_shift = 4.0
for i in range(n_active):
    size = rng.integers(15, 40)
    targets = rng.choice(N_GENES, size=size, replace=False)
    n_pos = size // 2
    pos, neg = targets[:n_pos].tolist(), targets[n_pos:].tolist()
    # Genuinely move these targets in the predicted direction - a REAL effect
    # layered on top of the same co-regulated-block background as Scenario B.
    y_c[pos] += active_shift + rng.normal(scale=0.5, size=len(pos))
    y_c[neg] -= active_shift + rng.normal(scale=0.5, size=len(neg))
    active_regulons.append((pos, neg))
regulons_c = regulons_b + active_regulons
known_active_idx = np.arange(len(regulons_b), len(regulons_c))
res_c = run_scenario("C (real signal + co-regulated background)", y_c, regulons_c, known_active_idx)

# --- Scenario D: lambda-hat on the user's REAL data ---------------------
print("\n=== Scenario D (real data: TRRUST + DoRothEA + CollecTRI) ===")
with open(os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"), "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
resolved_df, _ = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)
real_result = tr.run_regulon_activity(
    resolved_df, ["TRRUST", "DoRothEA", "CollecTRI"], "9606",
    metric=tx.GSEA_RANKING_METRIC_DEFAULT, min_regulon_size=5, n_permutations=200, seed=42,
)
for db, block in real_result["tables"].items():
    rows = block["table"]
    if rows.empty:
        print(f"  {db}: no TFs tested")
        continue
    scores = rows["Activity Score"].to_numpy(dtype=float)
    center, lam = robust_lambda(scores)
    raw_p = rows["p-value (parametric)"].to_numpy(dtype=float)
    fpr_raw = float(np.mean(raw_p < ALPHA))
    print(f"  {db}: n_tfs={len(scores)}  median={center:.4f}  lambda_hat={lam:.4f}  "
          f"frac_p<0.05={fpr_raw:.4f}")

# --- Summary / kill-criterion evaluation --------------------------------
print("\n" + "=" * 70)
print("KILL-CRITERION EVALUATION")
print("=" * 70)
print(f"Scenario A (sanity):     raw FPR = {res_a['fpr_raw']:.4f}, lambda = {res_a['lambda']:.4f}")
print(f"Scenario B (the test):   raw FPR = {res_b['fpr_raw']:.4f}, calibrated FPR = {res_b['fpr_cal']:.4f}, lambda = {res_b['lambda']:.4f}")
material = res_b["fpr_raw"] >= 0.08
print(f"\nMaterial inflation in Scenario B (raw FPR >= 0.08)? {material}")
if material:
    improved = res_b["fpr_cal"] < res_b["fpr_raw"] - 0.02
    print(f"Does calibration measurably improve it (cal FPR < raw FPR - 0.02)? {improved}")
    print(f"\n=> {'PROCEED to Step 2' if improved else 'DO NOT PROCEED - inflation real but calibration does not fix it, need a different approach'}")
else:
    print("\n=> KILL CRITERION HIT: inflation is not material. DO NOT PROCEED to Step 2.")
    print("   README's existing caveat stands; this script is the evidence.")
