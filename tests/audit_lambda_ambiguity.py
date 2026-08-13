"""THE CRITICAL AMBIGUITY: lambda cannot, by itself, tell "the null is
inflated by co-regulation" apart from "a lot of TFs are genuinely active".

RE-RUN after Part 1 replaced plain MAD/0.6745 with the truncated-normal MLE
(_truncated_normal_sigma) as the null-scale estimator. This file originally
(pre-fix) found lambda 1.74 (pure inflation) vs 1.57 (30% active, CLUMPED
signal) - i.e. INDISTINGUISHABLE. The whole point of Part 1 was to reduce
that; re-running with updated numbers rather than just re-passing the old
assertions, per this session's own rule that a claimed improvement must be
re-measured, not assumed.

RESULT, measured below: for the ORIGINAL clumped-signal fixture, the fix
works better than the probes even predicted - lambda now SEPARATES the two
scenarios cleanly (~1.6-1.7 vs ~1.0-1.1). That is real progress, not
smoothed over.

But the ambiguity is not eliminated, only narrowed. A second, harder
scenario is added below: WEAK genuine inflation (sigma ~1.2-1.3) vs
MODERATE graded/continuous signal (~30-50% active) - both land lambda in
the same ~1.0-1.4 band. That is where the surviving ambiguity now lives,
and it is where the QQ-plot hypothesis (does the plot's shape distinguish
them?) is re-tested, since the original clump-vs-inflation question this
file was built around is no longer the hard case.
"""
import os
import sys

import numpy as np
from scipy.special import ndtri

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_regulon as tr

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def central_slope(scores, frac=0.5):
    """Least-squares slope of the CENTRAL `frac` of a QQ plot (observed vs
    theoretical normal quantiles) - i.e. what the bulk of the points do,
    ignoring the tails."""
    s = np.sort(np.asarray(scores, dtype=float))
    n = len(s)
    theo = ndtri((np.arange(n) + 0.5) / n)
    lo = int(n * (1 - frac) / 2)
    hi = n - lo
    x, y = theo[lo:hi], s[lo:hi]
    xc, yc = x - x.mean(), y - y.mean()
    return float(np.dot(xc, yc) / np.dot(xc, xc))


def mean_lambda(gen, n_reps=30, seed0=0):
    lams = []
    for rep in range(n_reps):
        rng = np.random.default_rng(seed0 + rep)
        s = gen(rng)
        _, lam, _, _ = tr._empirical_null_calibration(s)
        lams.append(lam)
    return float(np.mean(lams))


N = 600

print("=== Scenario 1: UNIFORM INFLATION (co-regulation) - no TF truly active ===")
true_infl = 1.6
rng1 = np.random.default_rng(2027)
s_infl = rng1.normal(scale=true_infl, size=N)
_, lam_infl, _, _ = tr._empirical_null_calibration(s_infl)
slope_infl = central_slope(s_infl)
print(f"    lambda (truncated)    = {lam_infl:.3f}")
print(f"    central-50% QQ slope  = {slope_infl:.3f}")
print(f"    ratio central/lambda  = {slope_infl / lam_infl:.3f}")
check(abs(slope_infl / lam_infl - 1.0) < 0.15,
      f"under uniform inflation the BULK follows the fitted lambda line (ratio {slope_infl/lam_infl:.3f} ~ 1)")

print("\n=== Scenario 2: REAL SIGNAL - true null, 30% of TFs genuinely active (CLUMPED effect sizes) ===")
rng2 = np.random.default_rng(2027)
s_sig = rng2.normal(size=N)                       # true null sd = 1.0
n_act = int(0.30 * N)
s_sig[:n_act] += rng2.choice([-1, 1], size=n_act) * rng2.uniform(4, 9, size=n_act)
_, lam_sig, _, _ = tr._empirical_null_calibration(s_sig)
slope_sig = central_slope(s_sig)
print(f"    lambda (truncated)    = {lam_sig:.3f}")
print(f"    central-50% QQ slope  = {slope_sig:.3f}   <- the naive central-slope statistic is STILL fooled")
print(f"    ratio central/lambda  = {slope_sig / lam_sig:.3f}")
# The naive central-slope statistic is a diagnostic-only quantity this probe
# computes for comparison - it is NOT what the app's QQ plot draws (the app
# draws the fitted line at slope=lambda_used, confirmed exactly equal in
# playwright_regulon.py). It still suffers the same rank-axis-compression
# artifact as before (tail contamination pushes central RANKS outward
# regardless of which sigma estimator lambda now uses), so it still reads
# ~1.5 here even though lambda itself is now accurate (~1.0).
check(slope_sig > 1.3,
      f"the NAIVE central-slope statistic is still fooled by clumped signal (got {slope_sig:.3f}) - "
      f"this is a property of the raw statistic, not of lambda, and is why the app does not use it")

print("\n=== Can lambda ALONE tell Scenario 1 and 2 apart now? ===")
print(f"    Scenario 1 (pure inflation, no signal):        lambda = {lam_infl:.3f}")
print(f"    Scenario 2 (no inflation, 30% active, clumped): lambda = {lam_sig:.3f}")
separated_clump = abs(lam_infl - lam_sig) > 0.5
check(separated_clump,
      f"UPDATED RESULT (was 'indistinguishable' with plain MAD: 1.74 vs 1.57): with the truncated estimator "
      f"the two lambdas are now clearly separated ({lam_infl:.2f} vs {lam_sig:.2f}) for this CLEANLY-SEPARATED "
      f"clump-signal fixture - Part 1 works, and works better than the feasibility probes alone predicted")

print("\n=== Scenario 3 (the surviving hard case): WEAK genuine inflation vs MODERATE graded signal ===")
print("    Both land in a similar lambda band - this is where the ambiguity now actually lives.")


def graded_gen(frac):
    def gen(rng):
        s = rng.normal(size=N)
        n = int(frac * N)
        if n:
            s[:n] += rng.normal(scale=3.0, size=n)     # graded, not clumped: many small effects, some large
        return s
    return gen


lam_weak_infl_mean = mean_lambda(lambda rng: rng.normal(scale=1.3, size=N), n_reps=40, seed0=6000)
lam_graded45_mean = mean_lambda(graded_gen(0.45), n_reps=40, seed0=7000)
print(f"    weak genuine inflation (true sigma=1.3), mean lambda over 40 reps: {lam_weak_infl_mean:.3f}")
print(f"    graded signal, 45% active, mean lambda over 40 reps:               {lam_graded45_mean:.3f}")
overlap_graded = abs(lam_weak_infl_mean - lam_graded45_mean) < 0.35
check(overlap_graded,
      f"CONFIRMED: unlike the clump case above, weak inflation and moderate GRADED signal still land close "
      f"({lam_weak_infl_mean:.2f} vs {lam_graded45_mean:.2f}) - lambda alone remains genuinely ambiguous here, "
      f"exactly the graded-signal caveat _truncated_normal_sigma's docstring states")

print("\n=== Re-testing the QQ-plot hypothesis, specifically for Scenario 3 (the case that still needs it) ===")
# An OBSERVATION, not a validated diagnostic: with the new estimator, the
# ratio of the naive central-slope statistic to lambda tends to stay near
# 1.0 for genuine inflation (any strength) but rises above it under graded
# contamination. Reported here because it appeared unprompted while
# re-measuring this scenario, and because a probe that stayed silent about
# a lead it stumbled onto would be the same failure mode as smoothing over
# a bad number - but it is NOT proposed as a fix: single-run ranges still
# overlap at the boundary (measured separately, not shown here), it has
# only been tested on this one synthetic fixture shape, and turning it into
# an actual app feature would need its own dedicated validation this
# session was not scoped to do.
ratios_infl, ratios_graded = [], []
for rep in range(40):
    rng_a = np.random.default_rng(9000 + rep)
    s_a = rng_a.normal(scale=1.3, size=N)
    _, lam_a, _, _ = tr._empirical_null_calibration(s_a)
    ratios_infl.append(central_slope(s_a) / lam_a)

    rng_b = np.random.default_rng(9000 + rep)
    s_b = graded_gen(0.45)(rng_b)
    _, lam_b, _, _ = tr._empirical_null_calibration(s_b)
    ratios_graded.append(central_slope(s_b) / lam_b)
ratio_infl_mean = float(np.mean(ratios_infl))
ratio_graded_mean = float(np.mean(ratios_graded))
print(f"    mean central-slope/lambda ratio, weak genuine inflation: {ratio_infl_mean:.3f}")
print(f"    mean central-slope/lambda ratio, graded 45% active:      {ratio_graded_mean:.3f}")
print("    (an unvalidated lead, not a claim - see comment above)")

print("\n=== Conclusion for the README (updated after Part 1) ===")
print("    Lambda alone now SEPARATES pure inflation from cleanly-separated-clump signal - the case this")
print("    file originally showed was ambiguous under plain MAD is no longer ambiguous under the truncated")
print("    estimator. The ambiguity is narrowed, not gone: it survives specifically at the boundary between")
print("    WEAK genuine inflation and MODERATE graded (continuous) signal, where both still produce a")
print("    similarly elevated lambda. The QQ plot's naive central slope does not resolve this (it inherits")
print("    the same rank-axis-compression artifact as before). Therefore: calibrated FDR remains the")
print("    CONSERVATIVE default, the parametric column stays visible beside it, and the README should state")
print("    the narrowed claim precisely - 'materially reduces, does not eliminate' - rather than either the")
print("    old blanket ambiguity or an overclaim that Part 1 solved it.")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
