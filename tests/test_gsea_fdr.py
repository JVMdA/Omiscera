"""Offline correctness tests for the two FDR/q-value columns GSEA reports
side by side in every result table (Transcriptomics_module.py, also shared
verbatim by the Regulon activity module - see run_gsea's docstring on
multiple-testing families):

  - "FDR q-value" (_pooled_fdr) - the GSEA-paper pooled-permutation FDR
    (Subramanian et al. 2005 / fgsea-style): compares each term's observed
    NES against a null pool drawn from EVERY tested term's own permutation
    null, not a plain correction on nominal p. Positive- and negative-NES
    terms are independent multiple-testing families, matching desktop GSEA.
  - "FDR (BH)" (benjamini_hochberg) - the classic Benjamini-Hochberg
    procedure on nominal p-values, offered as a more familiar alternative
    (see the module's own comment on why it can be anti-conservative for
    correlated gene sets - this file does not re-litigate that tradeoff,
    only that the arithmetic itself is correct).

This file exists because _pooled_fdr's own docstring claims its monotonicity
"is verified explicitly in this module's unit tests" - an audit found no
such test anywhere in tests/. The vectorized implementation was cross-checked
by hand against a literal, unoptimized reimplementation of the GSEA paper's
formula and found to match exactly (0.0 max diff) before this file was
written, so these tests encode that same reference rather than merely
re-describing the vectorized code's own logic back at itself - a test that
just re-derives what the code already computes would pass even if the code
were wrong in a way the "reference" repeated by construction.
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import numpy as np

import Transcriptomics_module as tm

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# ===========================================================================
# _pooled_fdr: literal brute-force reference, independent of the vectorized
# implementation under test (no shared helper functions, no numpy tricks -
# just the paper's definition written out in plain loops).
# ===========================================================================

def _brute_force_pooled_fdr(observed_nes, size_of_term, null_nes_by_size, size_counts):
    S = len(observed_nes)
    q = np.zeros(S)
    for i in range(S):
        nes_i = observed_nes[i]
        sign = 1 if nes_i >= 0 else -1

        same_sign_obs = [x for x in observed_nes if (x >= 0) == (nes_i >= 0)]
        if sign > 0:
            obs_ge = sum(1 for x in same_sign_obs if x >= nes_i)
        else:
            obs_ge = sum(1 for x in same_sign_obs if x <= nes_i)
        obs_frac = obs_ge / len(same_sign_obs) if same_sign_obs else 0.0

        total_null = 0
        null_ge = 0
        for size, null_arr in null_nes_by_size.items():
            cnt = size_counts.get(size, 0)
            if cnt == 0:
                continue
            side_null = [x for x in null_arr if (x >= 0) == (sign > 0)]
            if not side_null:
                continue
            total_null += cnt * len(side_null)
            if sign > 0:
                null_ge += cnt * sum(1 for x in side_null if x >= nes_i)
            else:
                null_ge += cnt * sum(1 for x in side_null if x <= nes_i)
        null_frac = (null_ge / total_null) if total_null > 0 else 1.0

        q[i] = min(1.0, null_frac / obs_frac) if obs_frac > 0 else 1.0

    pos_idx = [i for i in range(S) if observed_nes[i] >= 0]
    neg_idx = [i for i in range(S) if observed_nes[i] < 0]
    for idx_group, keyf in (
        (pos_idx, lambda i: observed_nes[i]),
        (neg_idx, lambda i: -observed_nes[i]),
    ):
        idx_sorted = sorted(idx_group, key=keyf)  # ascending |NES|
        running_min = 1.0
        for i in idx_sorted:
            running_min = min(running_min, q[i])
            q[i] = running_min
    return q


def _random_scenario(rng, n_terms=40, sizes=(10, 20, 50, 100), nes_scale=1.3, null_scale=1.0, n_null=500):
    """A GSEA-shaped scenario: n_terms observed NES spread across a handful
    of size buckets, each bucket carrying its own permutation-null NES pool
    (exactly the shape _run_gsea_for_sets builds - one null per unique size,
    reused/reweighted across every term that size, not one null per term)."""
    size_of_term = rng.choice(sizes, size=n_terms)
    size_counts = {s: int(np.sum(size_of_term == s)) for s in sizes}
    null_nes_by_size = {s: rng.normal(0, null_scale, size=n_null) for s in sizes}
    observed_nes = rng.normal(0, nes_scale, size=n_terms)
    return observed_nes, size_of_term, null_nes_by_size, size_counts


print("=== _pooled_fdr matches the literal GSEA-paper formula exactly ===")
for seed in (1, 2, 3, 42, 12345):
    rng = np.random.default_rng(seed)
    observed_nes, size_of_term, null_nes_by_size, size_counts = _random_scenario(rng)
    q_real = tm._pooled_fdr(observed_nes, size_of_term, null_nes_by_size, size_counts)
    q_ref = _brute_force_pooled_fdr(observed_nes, size_of_term, null_nes_by_size, size_counts)
    max_diff = float(np.max(np.abs(q_real - q_ref)))
    check(max_diff < 1e-9,
          f"seed={seed}: vectorized _pooled_fdr matches the independent brute-force "
          f"reference to floating-point precision (max abs diff {max_diff:.2e})")

print("\n=== monotonicity: q must not increase as |NES| increases, within each sign ===")
for seed in (1, 2, 3):
    rng = np.random.default_rng(seed)
    observed_nes, size_of_term, null_nes_by_size, size_counts = _random_scenario(rng, n_terms=60)
    q = tm._pooled_fdr(observed_nes, size_of_term, null_nes_by_size, size_counts)

    pos = observed_nes >= 0
    order_pos = np.argsort(observed_nes[pos])  # ascending NES = ascending |NES| on the positive side
    q_pos_ordered = q[pos][order_pos]
    check(bool(np.all(np.diff(q_pos_ordered) <= 1e-12)),
          f"seed={seed}: positive-NES family - q is non-increasing as NES increases "
          f"(a more strongly enriched term never gets a WORSE q than a weaker one)")

    neg = ~pos
    order_neg = np.argsort(-observed_nes[neg])  # ascending |NES| for the negative side
    q_neg_ordered = q[neg][order_neg]
    check(bool(np.all(np.diff(q_neg_ordered) <= 1e-12)),
          f"seed={seed}: negative-NES family - q is non-increasing as |NES| increases")

    # SENSITIVITY: a backward running-min (the wrong direction, easy to get
    # backwards by accident - see _pooled_fdr's own docstring) would violate
    # this on a scenario with real variation in q. Confirm the check isn't
    # vacuously true by verifying q actually varies across ranks.
    check(float(np.std(q_pos_ordered)) > 1e-6,
          "SENSITIVITY: q genuinely varies across the positive-side ranks in this scenario "
          "(otherwise the monotonicity check above would pass trivially even with a broken "
          "direction, since a constant sequence is 'monotone' either way)")

print("\n=== bounds: every q in [0, 1] ===")
rng = np.random.default_rng(99)
observed_nes, size_of_term, null_nes_by_size, size_counts = _random_scenario(rng, n_terms=80, nes_scale=2.5)
q = tm._pooled_fdr(observed_nes, size_of_term, null_nes_by_size, size_counts)
check(bool(np.all((q >= 0.0) & (q <= 1.0))),
      f"all q-values land in [0, 1] even for very extreme NES (got range "
      f"[{q.min():.4f}, {q.max():.4f}])")

print("\n=== edge cases ===")
# No observed terms of one sign at all - the other sign's family must still
# compute cleanly and not error out or leak into the missing sign.
rng = np.random.default_rng(5)
all_pos_nes = np.abs(rng.normal(1.0, 0.5, size=10))
sizes_all = np.full(10, 20)
null_one_size = {20: rng.normal(0, 1, size=300)}
counts_one_size = {20: 10}
q_all_pos = tm._pooled_fdr(all_pos_nes, sizes_all, null_one_size, counts_one_size)
check(len(q_all_pos) == 10 and np.all((q_all_pos >= 0) & (q_all_pos <= 1)),
      f"an entirely positive-NES scenario (no negative-side terms at all) computes cleanly "
      f"(got {len(q_all_pos)} q-values, range [{q_all_pos.min():.3f}, {q_all_pos.max():.3f}])")

# A null pool that has no same-sign values at all for one term's sign -
# denom_null_total stays 0 for that side, which must fall back to q=1
# rather than raising or producing NaN/inf.
rng = np.random.default_rng(6)
one_nes = np.array([1.5, -1.5])
one_size = np.array([20, 20])
strictly_positive_null = {20: np.abs(rng.normal(1, 1, size=200))}  # no negative null values at all
one_count = {20: 2}
q_no_neg_null = tm._pooled_fdr(one_nes, one_size, strictly_positive_null, one_count)
check(q_no_neg_null[1] == 1.0,
      f"SENSITIVITY: a negative observed NES with NO negative null values to compare against "
      f"(denom_null_total=0 for that side) falls back to q=1.0 (non-significant), not NaN or "
      f"an unhandled division error (got {q_no_neg_null[1]!r})")
check(not np.isnan(q_no_neg_null).any(), "no NaN anywhere in the output")

print("\n=== benjamini_hochberg: the sibling 'FDR (BH)' column ===")


def _reference_bh(pvalues):
    """Independent, unoptimized reference implementation of the textbook
    Benjamini-Hochberg procedure, written without looking at
    Transcriptomics_module.benjamini_hochberg's own code."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    ranked = [pvalues[i] for i in order]
    q = [ranked[i] * m / (i + 1) for i in range(m)]
    for i in range(m - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = [min(x, 1.0) for x in q]
    out = [0.0] * m
    for rank_pos, orig in enumerate(order):
        out[orig] = q[rank_pos]
    return out


for seed in (1, 2, 3):
    rng = np.random.default_rng(seed)
    pvals = rng.uniform(0, 1, size=30).tolist() + [0.0005, 0.001, 0.5, 1.0]
    got = tm.benjamini_hochberg(pvals)
    expected = _reference_bh(pvals)
    max_diff = max(abs(a - b) for a, b in zip(got, expected))
    check(max_diff < 1e-10,
          f"seed={seed}: benjamini_hochberg matches an independent reference implementation "
          f"exactly (max abs diff {max_diff:.2e})")

pvals = [0.5, 0.001, 0.3, 0.01, 0.9, 0.02]
q_bh = tm.benjamini_hochberg(pvals)
order = sorted(range(len(pvals)), key=lambda i: pvals[i])
q_sorted = [q_bh[i] for i in order]
check(all(q_sorted[i] <= q_sorted[i + 1] + 1e-12 for i in range(len(q_sorted) - 1)),
      "q is non-decreasing when sorted by ascending p-value")
check(all(q_bh[i] >= pvals[i] - 1e-12 for i in range(len(pvals))),
      "the two BH invariants hold: q is never below its own nominal p-value")
check(tm.benjamini_hochberg([]) == [], "an empty p-value list returns an empty list, not an error")


print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
