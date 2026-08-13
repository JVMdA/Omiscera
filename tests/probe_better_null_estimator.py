"""FEASIBILITY PROBE, before planning any fix: is the lambda ambiguity
actually reducible, or is it inherent?

MAD fails under signal because contamination in the tails pushes the median
absolute deviation UP the null's own CDF (with 30% contamination the 50th
percentile of |dev| lands at the ~71st percentile of the null).

Candidate: a TRUNCATED estimator that uses only the central window |z| < c,
where the null dominates, and corrects analytically for the truncation.
For a centred normal truncated to |x| < c:
    E[x^2 | |x|<c] = sigma^2 * (1 - 2*c*phi(c/sigma) / (sigma*(2*Phi(c/sigma)-1)))
Solved by fixed-point iteration on sigma.

Compared against MAD on the SAME contamination grid, at the realistic TF
counts this module actually produces (TRRUST 134, DoRothEA 129, CollecTRI
494). Truth is sigma = 1.0 in every case, so the correct answer is lambda
= 1.0 - anything above that is the estimator being fooled by signal.
"""
import os
import sys

import numpy as np
from scipy.special import ndtr

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_regulon as tr


def _phi(x):
    return np.exp(-0.5 * x * x) / np.sqrt(2.0 * np.pi)


def truncated_sigma(scores, c_mult=1.5, iters=60):
    """Truncated-normal MLE of the null scale, using only |x - median| < c.
    c is set as c_mult * an initial MAD estimate, so the window adapts to the
    data rather than being a fixed absolute cutoff."""
    x = np.asarray(scores, dtype=float)
    x = x[np.isfinite(x)]
    center = float(np.median(x))
    d = x - center
    mad = float(np.median(np.abs(d)))
    if mad <= 0:
        return center, 1.0
    sigma = mad / 0.6745
    for _ in range(iters):
        c = c_mult * sigma
        inner = d[np.abs(d) < c]
        if inner.size < 10:
            break
        m2 = float(np.mean(inner ** 2))
        a = c / sigma
        # E[x^2 | |x|<c] / sigma^2 for a centred normal
        denom = (2.0 * ndtr(a) - 1.0)
        if denom <= 1e-12:
            break
        factor = 1.0 - (2.0 * a * _phi(a)) / denom
        if factor <= 1e-12:
            break
        new_sigma = np.sqrt(m2 / factor)
        if not np.isfinite(new_sigma) or new_sigma <= 0:
            break
        if abs(new_sigma - sigma) < 1e-9:
            sigma = new_sigma
            break
        sigma = new_sigma
    return center, float(sigma)


print("TRUTH in every row below is sigma = 1.0 (no inflation). Any lambda > 1")
print("is the estimator being fooled by REAL SIGNAL.\n")
print(f"{'n_TFs':>6} {'%active':>8} {'MAD lambda':>12} {'Truncated lambda':>18}   verdict")
print("-" * 70)

for n_tfs in (134, 494):
    for frac in (0.0, 0.10, 0.25, 0.40, 0.55):
        mad_runs, trunc_runs = [], []
        for rep in range(40):
            rng = np.random.default_rng(1000 + rep)
            s = rng.normal(size=n_tfs)                    # TRUE sigma = 1
            n_act = int(frac * n_tfs)
            if n_act:
                s[:n_act] += rng.choice([-1, 1], size=n_act) * rng.uniform(4, 9, size=n_act)
            _, lam_mad, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
            _, sig_tr = truncated_sigma(s)
            mad_runs.append(lam_mad)
            trunc_runs.append(max(1.0, sig_tr))
        mad_m = float(np.mean(mad_runs))
        tr_m = float(np.mean(trunc_runs))
        better = tr_m < mad_m - 0.03
        verdict = "truncated BETTER" if better else ("similar" if abs(tr_m - mad_m) <= 0.03 else "truncated WORSE")
        print(f"{n_tfs:>6} {int(frac*100):>7}% {mad_m:>12.3f} {tr_m:>18.3f}   {verdict}")
    print()

print("Reference: how far each estimator is from the TRUTH (1.0) matters more")
print("than which is smaller - an estimator that under-shoots below 1.0 would")
print("be clamped to 1.0 anyway, which is safe.")
