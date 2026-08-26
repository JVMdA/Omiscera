# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_module.run_quantitative_enrichment, entirely
offline (only "Custom database", same reasoning as test_proteomics_ora.py).

The single most important thing this test proves is the reason this module
exists as a separate change from Genomics' Gene-set Burden Test: log2FC is
SIGNED, so a term shifted DOWN must be exactly as detectable as one shifted
UP. That requires alternative="two-sided" in
Genomics_module._run_wilcoxon_for_namespace (added specifically for this
milestone, default still "greater" so the Burden Test itself is untouched -
see test_genomics_burden_alternative_unchanged below). A down-shifted term
is therefore the LOAD-BEARING fixture here: it must come out significant
under two-sided, and this test explicitly proves it would NOT under the old
one-sided default, by calling the shared engine both ways on the identical
data.
"""
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as genomics_ora  # noqa: E402
import Proteomics_module as pm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# =============================================================================
# Fixture: a resolved protein table with REAL continuous log2FC signal (not
# ORA's significant/non-significant binary - Quantitative Enrichment has no
# cutoff, every quantified protein contributes its score).
#   - SET_UP:    10 proteins, log2FC ~ N(+2.5, 0.3)  -> should enrich UP
#   - SET_DOWN:  10 proteins, log2FC ~ N(-2.5, 0.3)  -> should enrich DOWN
#   - SET_NOISE: 10 proteins, log2FC ~ N(0, 0.3), mixed in with background -> no real signal
#   - background: 170 proteins, log2FC ~ N(0, 0.3)
#   - a handful of NaN-score rows (unquantified in one condition), to prove
#     they are dropped and counted, not silently zero-filled.
#   - one deliberately DUPLICATED accession (two rows, same accession,
#     different log2FC) to prove _build_protein_scores mean-aggregates
#     rather than letting the last row silently win.
# =============================================================================
rng = np.random.default_rng(42)
N_BACKGROUND = 170
accessions = [f"P{i:05d}" for i in range(200)]

set_up = accessions[0:10]
set_down = accessions[10:20]
set_noise = accessions[20:30]
background = accessions[30:200]

rows = []
for acc in set_up:
    rows.append({"accession": acc, "log2FC": float(rng.normal(2.5, 0.3)), "pvalue": 1e-4, "padj": 1e-3})
for acc in set_down:
    rows.append({"accession": acc, "log2FC": float(rng.normal(-2.5, 0.3)), "pvalue": 1e-4, "padj": 1e-3})
for acc in set_noise + background:
    rows.append({"accession": acc, "log2FC": float(rng.normal(0.0, 0.3)), "pvalue": 0.7, "padj": 0.9})

# A few NaN rows: unquantified proteins that still resolved an accession.
for acc in ["P90001", "P90002", "P90003"]:
    rows.append({"accession": acc, "log2FC": float("nan"), "pvalue": float("nan"), "padj": float("nan")})

# Duplicate accession: P00000 (already in set_up) appears a SECOND time with
# a different log2FC - the two should mean-aggregate, not silently pick one.
rows.append({"accession": "P00000", "log2FC": 1.5, "pvalue": 1e-4, "padj": 1e-3})

for r in rows:
    r["input_id"] = r["accession"]
    r["gene_symbol"] = f"GENE_{r['accession']}"
    r["protein_name"] = ""
    r["status"] = "mapped"
resolved_df = pd.DataFrame(rows)

# Custom GMT: SET_UP, SET_DOWN, SET_NOISE, plus 50 single-member filler
# terms drawn from the flat background (no real shift either direction).
lines = [
    "SET_UP\tReal UP shift\t" + "\t".join(set_up),
    "SET_DOWN\tReal DOWN shift\t" + "\t".join(set_down),
    "SET_NOISE\tNo real shift\t" + "\t".join(set_noise),
]
filler_pool = background[:50]
for i, acc in enumerate(filler_pool):
    lines.append(f"FILLER_{i}\tflat background\t{acc}\t{filler_pool[(i + 1) % len(filler_pool)]}")
gmt_text = "\n".join(lines)
custom_db = genomics_ora.load_custom_geneset_db(gmt_text.encode("utf-8"), id_type="accession")
check(custom_db["n_terms"] == 3 + len(filler_pool),
      f"fixture custom DB has {3 + len(filler_pool)} terms (got {custom_db['n_terms']})")


# =============================================================================
# _build_protein_scores: metric formulas, NaN dropping, duplicate aggregation
# =============================================================================
print("=== _build_protein_scores ===")
scores, n_dropped = pm._build_protein_scores(resolved_df, metric="log2fc")
check(n_dropped == 3, f"the 3 NaN-score rows are dropped and counted (got n_dropped={n_dropped})")
check("P90001" not in scores and "P90002" not in scores and "P90003" not in scores,
      "NaN-score accessions do not appear in the score dict at all")
# P00000 appears twice: once inside set_up's N(2.5, 0.3) draw, once as the
# explicit 1.5 duplicate row - the stored score must be their mean, not
# either raw value alone.
p00000_rows = resolved_df[resolved_df["accession"] == "P00000"]["log2FC"].to_numpy()
check(len(p00000_rows) == 2, "fixture sanity: P00000 really does appear on 2 rows")
expected_mean = float(np.mean(p00000_rows))
check(abs(scores["P00000"] - expected_mean) < 1e-9,
      f"duplicate accession P00000 is mean-aggregated ({expected_mean:.4f}), not last-row-wins "
      f"(got {scores['P00000']:.4f})")
check(len(scores) == len(set(resolved_df["accession"])) - 3,  # -3 for the NaN rows' accessions
      f"one score per DISTINCT accession, duplicates collapsed (got {len(scores)} scores)")

scores_p, _ = pm._build_protein_scores(resolved_df, metric="pi_value")
check(scores_p["P00000"] != scores["P00000"], "pi_value differs from log2fc for a non-trivial p-value row")


# =============================================================================
# n_dropped counts DISTINCT PROTEINS with no score, not dropped ROWS - a real
# bug caught in review: resolve_protein_table does NOT deduplicate by
# accession, so an accession can have one valid row and one NaN row (partial
# duplicate). A row-level count would report that as 1 dropped, even though
# the protein IS scored (via its surviving row) and 0 proteins were actually
# lost - or, worse, it could report 2 dropped for a fully-NaN duplicate that
# is really only 1 missing protein. Both must be counted per DISTINCT
# accession, matching "Proteins scored" (also distinct-accession-based) so
# the two summary numbers the UI shows side by side are on the same basis.
# =============================================================================
print("=== n_dropped counts distinct un-scored PROTEINS, not dropped rows ===")
dup_rows = pd.DataFrame([
    {"accession": "DUPA", "log2FC": 2.0, "pvalue": 0.01, "padj": 0.05},   # DUPA: 1 valid + 1 NaN row
    {"accession": "DUPA", "log2FC": float("nan"), "pvalue": float("nan"), "padj": float("nan")},
    {"accession": "DUPB", "log2FC": float("nan"), "pvalue": float("nan"), "padj": float("nan")},  # DUPB: both rows NaN
    {"accession": "DUPB", "log2FC": float("nan"), "pvalue": float("nan"), "padj": float("nan")},
    {"accession": "SOLO", "log2FC": 1.0, "pvalue": 0.02, "padj": 0.05},   # SOLO: single valid row
])
dup_scores, dup_dropped = pm._build_protein_scores(dup_rows, metric="log2fc")
check(set(dup_scores.keys()) == {"DUPA", "SOLO"},
      f"DUPA is scored via its surviving row, DUPB has no score at all (got {sorted(dup_scores.keys())})")
check(dup_dropped == 1,
      f"exactly 1 DISTINCT protein (DUPB) has no score - NOT 3, which is what a row-level count "
      f"(DUPA's 1 NaN row + DUPB's 2 NaN rows) would wrongly report (got n_dropped={dup_dropped})")
check(len(dup_scores) + dup_dropped == dup_rows["accession"].nunique(),
      "scored + dropped exactly partitions the distinct accessions - the invariant the UI's "
      "'Proteins scored' / 'Dropped' metric cards rely on")

print("=== SENSITIVITY: the old row-level formula genuinely overcounts on this fixture ===")
# Reproduce the OLD (buggy) formula directly: n_input (valid-accession ROWS)
# minus the count of ROWS surviving the finite-score filter - exactly what
# _build_protein_scores computed before this fix.
n_input_rows = len(dup_rows)
n_finite_rows = dup_rows["log2FC"].apply(lambda v: v == v).sum()  # v == v is False for NaN, no numpy import needed here
old_style_n_dropped = n_input_rows - n_finite_rows
check(old_style_n_dropped == 3,
      f"confirmed: the OLD row-based formula gives 3 (1 from DUPA's NaN row + 2 from DUPB's) on this "
      f"fixture, not the correct 1 - proving the fix changed a real, reachable result "
      f"(got {old_style_n_dropped})")


# =============================================================================
# run_quantitative_enrichment: the two-sided fix is load-bearing
# =============================================================================
print("=== run_quantitative_enrichment: UP and DOWN both detected ===")
result = pm.run_quantitative_enrichment(
    resolved_df, databases=["Custom database"], taxon_id="9606", taxon_name="Homo sapiens",
    metric="log2fc", min_term_size=1, custom_db=custom_db, verbose=False,
)
block = result["tables"]["Custom database"]
table = block["table"]
check(not table.empty, "Custom database block has rows")
check(block["n_terms_tested"] == 3 + len(filler_pool),
      f"every size-eligible term is tested (no significance cutoff to exclude any) "
      f"(got {block['n_terms_tested']})")

row_up = table[table["Term ID"] == "SET_UP"].iloc[0]
row_down = table[table["Term ID"] == "SET_DOWN"].iloc[0]
row_noise = table[table["Term ID"] == "SET_NOISE"].iloc[0]

check(row_up["FDR"] < 0.01, f"SET_UP (positive shift) is significant (FDR={row_up['FDR']:.2e})")
check(row_up["Effect (mean difference)"] > 1.0, "SET_UP's effect is positive and large")

check(row_down["FDR"] < 0.01,
      f"SET_DOWN (negative shift) is ALSO significant under two-sided (FDR={row_down['FDR']:.2e}) "
      f"- this is the fix this milestone adds")
check(row_down["Effect (mean difference)"] < -1.0, "SET_DOWN's effect is negative and large")

check(row_noise["FDR"] > 0.05, f"SET_NOISE (no real shift) is NOT significant (FDR={row_noise['FDR']:.2e})")

# The proteomics member-column rename/format, same convention as ORA.
check(pm.OVERLAP_COLUMN in table.columns, f"{pm.OVERLAP_COLUMN!r} column present (Burden engine's 'Genes' renamed)")
check("Genes" not in table.columns, "the raw 'Genes' column name from the shared engine does not leak through")
sample_member = table[table["Term ID"] == "SET_UP"][pm.OVERLAP_COLUMN].iloc[0]
check("(" in sample_member and ")" in sample_member,
      f"members are in 'ACCESSION (SYMBOL)' format, matching ORA (got e.g. {sample_member!r})")


# =============================================================================
# THE load-bearing sensitivity check: prove two-sided is not vacuous by
# reverting to one-sided "greater" on the identical data and showing
# SET_DOWN would have been missed.
# =============================================================================
print("=== SENSITIVITY: one-sided 'greater' would have missed SET_DOWN ===")
universe_genes = set(resolved_df["accession"]) - {"P90001", "P90002", "P90003"}
term_map = {tid: {"name": info["name"], "genes": info["genes"]} for tid, info in custom_db["terms"].items()}

two_sided_block = genomics_ora._run_wilcoxon_for_namespace(
    term_map, scores, universe_genes, True, 1, None, alternative="two-sided",
)
one_sided_block = genomics_ora._run_wilcoxon_for_namespace(
    term_map, scores, universe_genes, True, 1, None,  # alternative defaults to "greater"
)
t_two = two_sided_block["table"]
t_one = one_sided_block["table"]
p_down_two = t_two[t_two["Term ID"] == "SET_DOWN"]["p-value"].iloc[0]
p_down_one = t_one[t_one["Term ID"] == "SET_DOWN"]["p-value"].iloc[0]
check(p_down_two < 1e-3,
      f"SET_DOWN's two-sided p-value is small (got {p_down_two:.2e})")
check(p_down_one > 0.99,
      f"the SAME data under the OLD one-sided 'greater' default gives SET_DOWN a p-value near 1 "
      f"(got {p_down_one:.4f}) - proving the fix genuinely changes the result, not just adds a parameter")
p_up_two = t_two[t_two["Term ID"] == "SET_UP"]["p-value"].iloc[0]
p_up_one = t_one[t_one["Term ID"] == "SET_UP"]["p-value"].iloc[0]
check(p_up_one < 1e-3, "SET_UP (a genuinely positive shift) is still caught under the OLD one-sided default")
check(abs(p_up_two - 2 * p_up_one) < 1e-6 or p_up_two < 1e-6,
      "for a one-directional signal, the two-sided p-value is ~2x the one-sided p-value (or both floored) - "
      "sane relationship between the two, not an unrelated computation")


# =============================================================================
# The Genomics Burden Test itself must be UNCHANGED (default alternative
# still "greater") - this milestone must not silently alter Genomics' own
# gene-set burden test just because it added a parameter to the shared
# engine.
# =============================================================================
print("=== Genomics' own Burden Test default is unchanged ===")
default_block = genomics_ora._run_wilcoxon_for_namespace(
    term_map, scores, universe_genes, True, 1, None,  # no `alternative` kwarg at all
)
check(default_block["table"].equals(one_sided_block["table"]),
      "calling _run_wilcoxon_for_namespace with NO alternative kwarg reproduces the explicit "
      "alternative='greater' call exactly - the new parameter's default is genuinely inert")


# =============================================================================
# No significance cutoff: every quantified protein with a finite score
# contributes, regardless of its own p-value/padj.
# =============================================================================
print("=== no cutoff: n_scored reflects ALL finite-score proteins, not a significant subset ===")
check(result["summary"]["n_scored"] == len(scores),
      f"summary n_scored matches _build_protein_scores' own count (got {result['summary']['n_scored']})")
check(result["summary"]["n_scored_dropped"] == 3,
      f"summary n_scored_dropped reports the 3 un-scored proteins, counted per DISTINCT accession "
      f"not per row (see the dedicated check above) (got {result['summary']['n_scored_dropped']})")
check(result["summary"]["metric"] == "log2fc", "summary echoes back the requested metric")


# =============================================================================
# Default metric is log2FC (Perseus' own convention), NOT GSEA's default
# metric (signed_log_pvalue) - a real, deliberate divergence from the GSEA
# constant this module otherwise reuses.
# =============================================================================
print("=== default metric ===")
check(pm.QUANT_METRIC_DEFAULT == "log2fc", "Quantitative Enrichment defaults to log2FC")
result_default = pm.run_quantitative_enrichment(
    resolved_df, databases=["Custom database"], taxon_id="9606", taxon_name="Homo sapiens",
    min_term_size=1, custom_db=custom_db, verbose=False,  # metric omitted entirely
)
check(result_default["summary"]["metric"] == "log2fc", "omitting `metric` uses log2FC, not GSEA's own default")


print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    sys.exit(1)
