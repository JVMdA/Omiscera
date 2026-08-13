"""Two conditions the truncated estimator MUST satisfy before it can replace
MAD. Fixing the false-alarm case is worthless if it breaks the case
calibration exists for.

  1. It must still RECOVER GENUINE INFLATION. Under true co-regulation the
     whole distribution is scaled; the estimator has to report that scale,
     not flatten it to 1.0. If it can't, calibration stops working entirely.
  2. It must hold up under CONTINUOUS signal, not just the cleanly-separated
     |4-9| clumps used in the first probe. Real TF activity is graded, so
     some genuinely-active TFs sit INSIDE the truncation window.
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


print("=== CONDITION 1: does it still recover GENUINE inflation? ===")
print("(no signal at all; the whole null is truly scaled by the stated factor)")
print(f"{'true sigma':>11} {'MAD lambda':>12} {'Truncated lambda':>18}   recovers?")
print("-" * 62)
ok1 = True
for true_sigma in (1.0, 1.3, 1.6, 2.0, 3.0):
    mads, trs = [], []
    for rep in range(40):
        rng = np.random.default_rng(500 + rep)
        s = rng.normal(scale=true_sigma, size=300)
        _, lam_mad, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
        _, sig_tr = truncated_sigma(s)
        mads.append(lam_mad); trs.append(max(1.0, sig_tr))
    m, t = float(np.mean(mads)), float(np.mean(trs))
    rel_err = abs(t - true_sigma) / true_sigma
    good = rel_err < 0.08
    ok1 = ok1 and good
    print(f"{true_sigma:>11.1f} {m:>12.3f} {t:>18.3f}   {'YES' if good else 'NO  (err %.0f%%)' % (rel_err*100)}")

print("\n=== CONDITION 2: CONTINUOUS (graded) signal, not separated clumps ===")
print("(active TFs drawn from a wide normal, so many land inside the window)")
print(f"{'%active':>8} {'MAD lambda':>12} {'Truncated lambda':>18}   (truth = 1.0)")
print("-" * 62)
for frac in (0.0, 0.10, 0.25, 0.40):
    mads, trs = [], []
    for rep in range(40):
        rng = np.random.default_rng(900 + rep)
        s = rng.normal(size=300)
        n_act = int(frac * 300)
        if n_act:
            # graded effect sizes: many small, some large - the realistic shape
            s[:n_act] += rng.normal(scale=3.0, size=n_act)
        _, lam_mad, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
        _, sig_tr = truncated_sigma(s)
        mads.append(lam_mad); trs.append(max(1.0, sig_tr))
    m, t = float(np.mean(mads)), float(np.mean(trs))
    print(f"{int(frac*100):>7}% {m:>12.3f} {t:>18.3f}   truncated is {'better' if t < m - 0.03 else 'similar'}")

print("\n=== VERDICT ===")
print(f"Condition 1 (recovers genuine inflation): {'PASS' if ok1 else 'FAIL - estimator unusable'}")
print("Condition 2: see table - continuous signal is harder than separated clumps for BOTH,")
print("so the honest claim is 'materially better', not 'solved'.")
