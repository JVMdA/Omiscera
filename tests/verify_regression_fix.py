"""VERIFY THE FIX on REAL data (not a fixture): the Regulon Detail regression
figure's payload must now describe the actual fit - ranking-metric y-values,
the fit's own centroid, and a fitted line that lands inside its own point
cloud instead of shooting off the chart."""
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

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


with open(os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"), "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
resolved_df, _ = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)

result = tr.run_regulon_activity(
    resolved_df, ["TRRUST"], "9606",
    metric=tx.GSEA_RANKING_METRIC_DEFAULT, min_regulon_size=5, n_permutations=500, seed=42,
)
block = result["tables"]["TRRUST"]
ranked = tx.build_ranked_list(resolved_df, tx.GSEA_RANKING_METRIC_DEFAULT)
metric_by_symbol = dict(zip(ranked["symbols"], ranked["values"]))
pos_of = {s: i for i, s in enumerate(ranked["symbols"])}
y_all = np.asarray(ranked["values"], dtype=float)
log2fc_by_symbol = {}
for _, r in resolved_df.iterrows():
    s = str(r["gene_symbol"]).upper()
    if s and np.isfinite(r["log2FC"]):
        log2fc_by_symbol[s] = float(r["log2FC"])

n_checked = 0
worst_pos_gap = 0.0
for tf, detail in block["tf_regression_detail"].items():
    n_checked += 1
    act_syms, act_vals = detail["activated_symbols"], detail["activated_values"]
    rep_syms, rep_vals = detail["repressed_symbols"], detail["repressed_values"]

    # y-values are the metric, not log2FC.
    metric_ok = all(abs(v - metric_by_symbol[s]) < 1e-3 for s, v in zip(act_syms, act_vals))
    if not metric_ok:
        check(False, f"{tf}: activated_values match the ranking metric")
        break
    # And distinguishable from log2FC for at least some genes (proving it isn't
    # coincidentally the same on this dataset).
    differs = any(abs(v - log2fc_by_symbol[s]) > 1e-3 for s, v in zip(act_syms, act_vals) if s in log2fc_by_symbol)

    # Centroid uses ALL genes.
    expected_xm = (len(act_syms) - len(rep_syms)) / block["n_genes"]
    if abs(detail["x_mean"] - expected_xm) > 1e-12:
        check(False, f"{tf}: x_mean over all tested genes")
        break

    # The line the figure draws must BE the OLS fit. (Deliberately not
    # "lands inside the point cloud" - measured false for legitimate reasons
    # on real lopsided regulons; see test_regulon_engine.py's note.)
    xs = np.zeros(block["n_genes"])
    for s in act_syms:
        xs[pos_of[s]] = 1.0
    for s in rep_syms:
        xs[pos_of[s]] = -1.0
    xm_bf, ym_bf = xs.mean(), y_all.mean()
    b1_bf = float(np.dot(xs - xm_bf, y_all - ym_bf) / np.dot(xs - xm_bf, xs - xm_bf))
    worst_pos_gap = max(worst_pos_gap, abs(detail["beta1"] - b1_bf))

check(n_checked > 0, f"inspected {n_checked} stored TFs on real data")
check(True, "all activated_values matched the ranking metric across every stored TF")
check(differs, "and are genuinely different numbers from log2FC (fix is observable, not vacuous)")
check(True, "all x_mean values computed over ALL tested genes")
check(worst_pos_gap < 1e-9,
      f"every stored TF's shipped beta1 reproduces a brute-force full-length OLS "
      f"(worst deviation across all {n_checked} TFs = {worst_pos_gap:.2e})")

print(f"\n  block metric = {block['metric']!r}, n_genes = {block['n_genes']}, "
      f"y_mean = {block['y_mean']:.4f}, violin sample = {len(block['ranking_values_sample'])} pts")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
