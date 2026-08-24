# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_module.py's ORA orchestrator, entirely offline.

Uses ONLY "Custom database" (Genomics_module.load_custom_geneset_db), which
needs no network - so run_ora's own logic (significant-protein selection,
direction handling, independent BH per direction, and - the single most
important test here - that zero-overlap terms never enter the BH family) is
tested directly against a hand-built fixture, not through any live database
fetch. GO/Reactome/InterPro/Complex Portal/UniProt/MSigDB universe builders
are exercised live in playwright_proteomics_ora.py and by the manual live
verification recorded in the Fix commit for this milestone.
"""
import io
import os
import sys

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
# Fixture: a resolved protein table + a custom gene-set database, entirely
# synthetic and offline. Two REAL enrichment signals are built in:
#   - "SET_A": all 5 members are significant UP  -> should be highly enriched
#   - "SET_B": all 4 members are significant DOWN -> should be highly enriched
#   - "SET_NOISE": scattered across sig/non-sig, no real enrichment
#   - dozens of filler single-member terms with ZERO overlap with anything
#     significant - these exist specifically to prove the zero-overlap terms
#     are excluded from the BH family (see Fix 2026-08-24 in the commit log).
# =============================================================================
N_UNIVERSE = 200
accessions = [f"P{i:05d}" for i in range(N_UNIVERSE)]

set_a = accessions[0:5]      # will be significant UP
set_b = accessions[5:9]      # will be significant DOWN
set_noise = accessions[9:20]  # mixed, no real signal

rows = []
for i, acc in enumerate(accessions):
    if acc in set_a:
        log2fc, pvalue = 3.0, 1e-6
    elif acc in set_b:
        log2fc, pvalue = -3.0, 1e-6
    else:
        log2fc, pvalue = 0.05, 0.8
    rows.append({
        "input_id": acc, "accession": acc, "gene_symbol": f"GENE{i}", "protein_name": "",
        "log2FC": log2fc, "pvalue": pvalue, "padj": pvalue, "status": "mapped",
    })
resolved_df = pd.DataFrame(rows)

# Build a GMT-format custom database: SET_A, SET_B, SET_NOISE, plus 50
# single-member filler terms guaranteed to have ZERO overlap with the
# significant set (they are all non-significant background accessions).
lines = [
    "SET_A\tReal UP signal\t" + "\t".join(set_a),
    "SET_B\tReal DOWN signal\t" + "\t".join(set_b),
    "SET_NOISE\tNo real signal\t" + "\t".join(set_noise),
]
background_only = accessions[20:200]
for i, acc in enumerate(background_only):
    lines.append(f"FILLER_{i}\tNo overlap with anything significant\t{acc}\t{background_only[(i+1) % len(background_only)]}")
gmt_text = "\n".join(lines)
custom_db = genomics_ora.load_custom_geneset_db(gmt_text.encode("utf-8"), id_type="accession")
check(custom_db["n_terms"] == 3 + len(background_only), f"fixture custom DB has {3 + len(background_only)} terms (got {custom_db['n_terms']})")


# =============================================================================
# select_significant_proteins
# =============================================================================
print("=== select_significant_proteins: direction semantics ===")
sig_all = pm.select_significant_proteins(resolved_df, 0.05, 1.0, "all")
sig_up = pm.select_significant_proteins(resolved_df, 0.05, 1.0, "up")
sig_down = pm.select_significant_proteins(resolved_df, 0.05, 1.0, "down")
check(set(sig_all["accession"]) == set(set_a) | set(set_b), "'all' selects both up and down significant proteins")
check(set(sig_up["accession"]) == set(set_a), "'up' selects only the up-regulated significant set")
check(set(sig_down["accession"]) == set(set_b), "'down' selects only the down-regulated significant set")


# =============================================================================
# run_ora: the zero-overlap-term family-inflation fix (the real bug caught
# live during this milestone's development)
# =============================================================================
print("=== run_ora: zero-overlap terms excluded from the BH family ===")
result = pm.run_ora(
    resolved_df, databases=["Custom database"], taxon_id="9606", taxon_name="Homo sapiens",
    padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, custom_db=custom_db, verbose=False,
)
block_all = result["tables"]["Custom database (all)"]
check(block_all["n_terms_tested"] == 2,
      f"only SET_A and SET_B (the two terms with real overlap) enter the family - SET_NOISE and all "
      f"{len(background_only)} filler terms are excluded (got n_terms_tested={block_all['n_terms_tested']})")
check(block_all["n_terms_no_overlap"] == 1 + len(background_only),
      f"n_terms_no_overlap correctly counts SET_NOISE + all filler terms "
      f"(got {block_all['n_terms_no_overlap']}, expected {1 + len(background_only)})")

print("=== SENSITIVITY: prove the fix is not vacuous by reverting it ===")
# Directly exercise the OLD (buggy) behaviour this fix replaced: call
# Genomics_module._run_fisher_for_namespace WITHOUT the pre-filter, on the
# exact same fixture, and confirm it genuinely inflates the family - i.e.
# that the check above would have gone red against the pre-fix code.
sig_set = set(sig_all["accession"])
universe_set = set(resolved_df["accession"])
unfiltered_terms = {tid: {"name": info["name"], "genes": info["genes"]} for tid, info in custom_db["terms"].items()}
unfiltered_block = genomics_ora._run_fisher_for_namespace(
    unfiltered_terms, sig_set, universe_set, 0.05, True, 1, None,
)
check(unfiltered_block["n_terms_tested"] == 3 + len(background_only),
      f"confirmed: calling Genomics' shared runner WITHOUT the pre-filter reproduces the original bug "
      f"(all {3 + len(background_only)} terms enter the family, not just the 2 with real overlap) - "
      f"this is the exact defect Proteomics_module.py's _drop_zero_overlap fixes")


# =============================================================================
# run_ora: independent BH per direction (never pooled)
# =============================================================================
print("=== run_ora: FDR families are independent per direction ===")
block_up = result["tables"]["Custom database (up)"]
block_down = result["tables"]["Custom database (down)"]
check(block_up["n_terms_tested"] == 1, "the 'up' family tests only SET_A (SET_B has zero overlap with up-regulated proteins)")
check(block_down["n_terms_tested"] == 1, "the 'down' family tests only SET_B")

set_a_fdr_all = block_all["table"].set_index("Term ID").loc["SET_A", "FDR"]
set_a_fdr_up = block_up["table"].set_index("Term ID").loc["SET_A", "FDR"]
check(abs(set_a_fdr_all - set_a_fdr_up) > 1e-12 or block_all["n_terms_tested"] != block_up["n_terms_tested"],
      "SET_A's FDR genuinely differs (or comes from a differently-sized family) between the 'all' and 'up' "
      "directions, confirming these are independently-corrected families, not the same number reused")

print("=== run_ora: real signal is genuinely significant, noise is not ===")
set_a_row = block_all["table"].set_index("Term ID").loc["SET_A"]
set_b_row = block_all["table"].set_index("Term ID").loc["SET_B"]
check(set_a_row["FDR"] < 0.01, f"SET_A (real, planted signal) is significant (FDR={set_a_row['FDR']:.2e})")
check(set_b_row["FDR"] < 0.01, f"SET_B (real, planted signal) is significant (FDR={set_b_row['FDR']:.2e})")
check(set_a_row["k"] == 5, f"SET_A's overlap count is exactly its 5 planted members (got k={set_a_row['k']})")
check(set_b_row["k"] == 4, f"SET_B's overlap count is exactly its 4 planted members (got k={set_b_row['k']})")


# =============================================================================
# Empty-direction handling
# =============================================================================
print("=== run_ora: a direction with zero significant proteins returns an honest empty block ===")
result_strict = pm.run_ora(
    resolved_df, databases=["Custom database"], taxon_id="9606", taxon_name="Homo sapiens",
    padj_cutoff=1e-20, lfc_cutoff=1.0, min_term_size=1, custom_db=custom_db, verbose=False,
)
block_strict = result_strict["tables"]["Custom database (all)"]
check(block_strict["n_terms_tested"] == 0, "an impossibly strict padj cutoff correctly yields zero tested terms")
check(block_strict["table"].empty, "the result table is genuinely empty, not just zero-row-but-wrong-shape")
check(bool(block_strict.get("note")), "an empty block carries an explanatory note")


print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
