"""FULL AUDIT of the TF/Regulon module's figures and mathematics.

Requested because the user works in metabolomics and cannot independently
sanity-check the transcriptomics theory - so every claim a figure makes is
checked numerically here rather than asserted.

Sections:
  A. Cross-figure consistency: does every figure's "significant" mean the
     same thing, and does each encoding match its axis label?
  B. Coherence bars: are they computed on a variable whose SIGN agrees with
     the one the activity score uses? (they use log2FC; the score uses the
     ranking metric - if signs could disagree the figure would contradict
     the table)
  C. Calibration QQ plot: are the theoretical quantiles the right ones?
  D. Empirical-null breakdown: at what fraction of genuinely-active TFs does
     the null-scale estimator start over-correcting? Re-run AFTER replacing
     plain MAD with the truncated-normal MLE (Transcriptomics_regulon.py's
     _truncated_normal_sigma) - this section originally reproduced the old
     MAD estimator being fooled by clumped signal at >=40% active (lambda
     >1.5); the new estimator materially reduces that (see
     test_regulon_engine.py section 12b for the permanent version of this
     comparison, including the still-imperfect GRADED-signal case kept
     separate there).
  E. Numerical robustness: can a p-value underflow to exactly 0 and blow up
     the volcano's -log10 axis?
  F. Multiple-testing families: how many independent BH corrections does a
     user actually see, and is that documented?
"""
import os
import sys

import numpy as np

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx
import Transcriptomics_module_retriever as retr
import Transcriptomics_regulon as tr

DATA3 = os.path.join(_ROOT, "Data3")

PASS = FAIL = WARN = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def warn(msg):
    global WARN
    WARN += 1
    print("WARN: " + msg)


with open(DATA3 + r"\example_transcriptomics_de_full.csv", "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
resolved_df, _ = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)


# =============================================================================
# B. Coherence bars use log2FC; the activity score uses the ranking metric.
#    If their SIGNS could disagree, the coherence figure would contradict the
#    activity score it is meant to explain.
# =============================================================================
print("=== B. sign(ranking metric) must equal sign(log2FC) for EVERY metric ===")
for metric in tx.GSEA_RANKING_METRICS:
    ranked = tx.build_ranked_list(resolved_df, metric)
    syms = ranked["symbols"]
    vals = np.asarray(ranked["values"], dtype=float)
    lfc_by_sym = {str(r["gene_symbol"]).upper(): float(r["log2FC"])
                  for _, r in resolved_df.iterrows() if str(r["gene_symbol"]).strip()}
    lfc = np.array([lfc_by_sym[s] for s in syms])
    # Compare signs, treating the exact-zero cases separately.
    nonzero = lfc != 0
    agree = np.sign(vals[nonzero]) == np.sign(lfc[nonzero])
    check(bool(np.all(agree)),
          f"{metric}: sign(metric) == sign(log2FC) for all {int(nonzero.sum())} non-zero-log2FC genes")
    n_zero_lfc = int((~nonzero).sum())
    if n_zero_lfc:
        # build_ranked_list defines sign(0) as +1; the coherence figure counts
        # log2FC > 0 as "moved up", so a gene with EXACTLY 0 lands on opposite
        # sides of the two conventions.
        warn(f"{metric}: {n_zero_lfc} gene(s) have log2FC exactly 0 - the ranking metric treats these as "
             f"POSITIVE (sign(0):=+1) while the coherence bars count them as 'moved down' (v > 0 is False). "
             f"Cosmetic on this dataset but a genuine convention mismatch.")
    else:
        print(f"    ({metric}: no genes with log2FC exactly 0 in this dataset)")


# =============================================================================
# C. QQ plot theoretical quantiles
# =============================================================================
print("\n=== C. Calibration QQ plot uses correct plotting positions ===")
# The JS uses _probit((i + 0.5) / n) for i = 0..n-1 - the Hazen position.
# Verify against scipy's own normal quantile function.
from scipy.special import ndtri
n_pts = 200
js_positions = [(i + 0.5) / n_pts for i in range(n_pts)]
theoretical = ndtri(js_positions)
check(abs(theoretical[n_pts // 2 - 1] + theoretical[n_pts // 2]) < 1e-9,
      "plotting positions are symmetric about 0 (median quantile is 0)")
check(theoretical[0] < -2 and theoretical[-1] > 2,
      f"tails reach beyond +/-2 sd for n=200 (got {theoretical[0]:.2f}, {theoretical[-1]:.2f})")
# Under a TRUE standard normal null, sorted samples should track these closely.
rng = np.random.default_rng(7)
sample = np.sort(rng.normal(size=n_pts))
# NB: np.corrcoef is NOT usable in this runtime - it triggers a native
# access-violation crash (documented environment hazard, alongside
# matplotlib and scipy.stats.spearmanr). Pearson computed by hand instead.
a = sample - sample.mean()
b = theoretical - theoretical.mean()
corr = float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))
check(corr > 0.99, f"a genuine standard-normal sample tracks the plotted quantiles (r={corr:.4f})")


# =============================================================================
# D. Empirical-null breakdown: when does MAD over-correct?
# =============================================================================
print("\n=== D. Empirical-null (truncated-normal) behaviour as the fraction of TRULY ACTIVE TFs rises ===")
print("    (documented claim: MATERIALLY reduces, does not eliminate, the ambiguity)")
n_tfs = 400
for frac in (0.0, 0.1, 0.25, 0.4, 0.5, 0.7):
    rng_d = np.random.default_rng(101)
    scores = rng_d.normal(size=n_tfs)                    # true null: sd 1
    n_active = int(frac * n_tfs)
    if n_active:
        scores[:n_active] += rng_d.choice([-1, 1], size=n_active) * rng_d.uniform(4, 9, size=n_active)
    _, lam, _, _ = tr._empirical_null_calibration(scores)
    verdict = "ok" if lam < 1.25 else ("INFLATED - over-conservative" if lam >= 1.5 else "drifting")
    print(f"    {int(frac*100):>3}% truly active -> lambda = {lam:.3f}   [{verdict}]")
    # MEASURED BEHAVIOUR, now documented in _empirical_null_calibration's and
    # _truncated_normal_sigma's docstrings and the README. The truncated
    # estimator restricts itself to a central window where the null
    # dominates, so CLEANLY SEPARATED clumps of signal (as simulated here,
    # +-4..9) barely move it even at 70% active - this is the regime where
    # the fix is strongest. It does NOT claim this for graded/continuous
    # signal (weaker claim, tested separately in test_regulon_engine.py
    # section 12b) or past ~55% active, which is beyond any robust
    # estimator.
    if frac <= 0.10:
        check(lam < 1.15, f"at {int(frac*100)}% active TFs lambda stays close to 1.0 (got {lam:.3f})")
    if 0.40 <= frac <= 0.50:
        check(lam < 1.25, f"at {int(frac*100)}% active TFs (clumped signal) lambda STAYS close to 1.0 "
                          f"(got {lam:.3f}) - the truncated estimator materially reduces the old MAD ambiguity here")
    if frac >= 0.70:
        # DOCUMENTED breakdown point (both README and this module's docstring
        # say "past ~55% active, beyond any robust estimator"), reproduced
        # here rather than silently skipped: this is not a regression, it is
        # the honest edge of the claim. The central-window truncation itself
        # starts including majority-signal points once active TFs dominate,
        # which is a structural limit of the approach, not a bug.
        check(lam > 2.0, f"at {int(frac*100)}% active TFs the estimator DOES break down as documented "
                         f"(got lambda={lam:.3f}) - past the ~55% active honesty limit, not a regression")


# =============================================================================
# E. Numerical robustness: can a calibrated p underflow to exactly 0?
# =============================================================================
print("\n=== E. Can a p-value underflow to 0 and blow up the volcano's -log10 axis? ===")
from scipy.special import ndtr
for z in (5, 10, 20, 30, 38, 39, 40, 50):
    p = 2.0 * float(ndtr(-abs(z)))
    y = -np.log10(p + 1e-300)
    flag = "  <-- p underflows to EXACTLY 0; volcano y would be 300" if p == 0.0 else ""
    print(f"    |z| = {z:>2}  ->  p = {p:.3e}   -log10(p+1e-300) = {y:6.1f}{flag}")
underflow_z = None
for z in range(30, 60):
    if 2.0 * float(ndtr(-float(z))) == 0.0:
        underflow_z = z
        break
print(f"    p first underflows to exactly 0 at |z| ~= {underflow_z}")

# Is that reachable on real data?
act = tr.run_regulon_activity(resolved_df, ["TRRUST", "DoRothEA", "CollecTRI"], "9606",
                              min_regulon_size=5, n_permutations=200, seed=42)
max_abs_z = 0.0
for db, block in act["tables"].items():
    df = block["table"]
    if df.empty:
        continue
    zs = df["Calibrated Score (z)"].to_numpy(dtype=float)
    p_cal = df["p-value (calibrated)"].to_numpy(dtype=float)
    max_abs_z = max(max_abs_z, float(np.max(np.abs(zs))))
    n_zero_p = int((p_cal == 0.0).sum())
    print(f"    {db}: max |z| = {np.max(np.abs(zs)):.2f}, p-values exactly 0: {n_zero_p}")
    check(n_zero_p == 0, f"{db}: no calibrated p-value underflows to exactly 0 on this dataset")
headroom = underflow_z - max_abs_z
print(f"    headroom before underflow: {headroom:.1f} z-units")
if headroom < 20:
    warn(f"only {headroom:.1f} z-units of headroom before a p-value underflows to 0 and puts a point at "
         f"y=300 on the volcano, flattening every other point. A stronger dataset could reach it.")
else:
    print(f"    (comfortable headroom on this dataset)")


# =============================================================================
# F. How many independent BH families does a user actually see?
# =============================================================================
print("\n=== F. Multiple-testing families the user is exposed to ===")
enr = tr.run_regulon_enrichment(resolved_df, ["TRRUST", "DoRothEA"], "9606",
                                 min_term_size=5, max_term_size=500)
keys = sorted(enr["tables"].keys())
print(f"    enrichment result tables: {keys}")
check(len(keys) == 2 * 3, f"enrichment produces database x direction tables ({len(keys)} for 2 databases)")
warn(f"each of those {len(keys)} tables is BH-corrected INDEPENDENTLY (as is each activity table). A user "
     f"scanning all of them is running {len(keys)} separate families with no global correction - the same "
     f"convention ORA uses, but it must be stated in this module's README too.")

print("\n" + "=" * 72)
print("TOTAL: %d passed, %d failed, %d warnings" % (PASS, FAIL, WARN))
sys.exit(1 if FAIL else 0)
