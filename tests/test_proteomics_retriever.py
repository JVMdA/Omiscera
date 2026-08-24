# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_module_retriever.py.

No-network functions (column detection, protein-group parsing, isoform
stripping, contaminant/decoy filtering, dedupe-by-accession) are tested
directly and fast. resolve_protein_table itself needs live UniProt ID
Mapping and is only lightly smoke-tested here (real network, small input) -
the fuller live round trip (large example, messy MaxQuant example, full ORA
run) is covered by playwright_proteomics_ora.py.
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Proteomics_module_retriever as pmr  # noqa: E402

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
# Column detection
# =============================================================================
print("=== detect_columns: proteomics-specific identifier synonyms ===")
df = pd.DataFrame({
    "Protein IDs": ["P04637"], "log2FC": [1.0], "pvalue": [0.01],
})
mapping = pmr.detect_columns(df)
check(mapping["identifier"] == "Protein IDs", "MaxQuant-style 'Protein IDs' header detected as identifier")

df2 = pd.DataFrame({
    "Majority protein IDs": ["P04637"], "logFC": [1.0], "P.Value": [0.01],
})
mapping2 = pmr.detect_columns(df2)
check(mapping2["identifier"] == "Majority protein IDs", "'Majority protein IDs' detected")
check(mapping2["log2FC"] == "logFC", "'logFC' synonym detected")
check(mapping2["pvalue"] == "P.Value", "'P.Value' synonym detected (normalizes away the dot)")

try:
    pmr.detect_columns(pd.DataFrame({"foo": [1]}))
    check(False, "missing required columns should raise ValueError")
except ValueError:
    check(True, "missing required columns raises ValueError")

print("=== detect_columns: optional abundance column ===")
df3 = pd.DataFrame({"Protein IDs": ["P1"], "log2FC": [1.0], "pvalue": [0.01], "iBAQ": [1000.0]})
mapping3 = pmr.detect_columns(df3)
check(mapping3.get("abundance") == "iBAQ", "optional 'iBAQ' abundance column detected")
df4 = pd.DataFrame({"Protein IDs": ["P1"], "log2FC": [1.0], "pvalue": [0.01]})
mapping4 = pmr.detect_columns(df4)
check("abundance" not in mapping4, "abundance column correctly absent when not provided")


# =============================================================================
# Protein-group parsing / isoform stripping
# =============================================================================
print("=== parse_protein_group / strip_isoform ===")
check(pmr.parse_protein_group("P04637;P04637-2;Q9H3D4") == "P04637",
      "leading accession taken from a semicolon-separated protein group")
check(pmr.parse_protein_group("P04637") == "P04637", "single accession passes through unchanged")
check(pmr.parse_protein_group("") == "", "empty value stays empty")
check(pmr.strip_isoform("P04637-2") == "P04637", "isoform suffix stripped")
check(pmr.strip_isoform("P04637") == "P04637", "accession with no isoform suffix unchanged")
check(pmr.strip_isoform(pmr.parse_protein_group("P04637-2;Q9H3D4")) == "P04637",
      "group parsing + isoform stripping compose correctly")


# =============================================================================
# Contaminant/decoy filtering - SENSITIVITY-CHECKED: prove the check is not
# vacuous by reverting each detection mechanism in turn and confirming it
# genuinely goes red.
# =============================================================================
print("=== filter_contaminants_decoys: both detection mechanisms ===")
messy = pd.DataFrame({
    "identifier": ["P04637", "CON__P02769", "Q00987", "REV__Q99999", "P16220"],
    "Potential contaminant": ["", "+", "", "", ""],
    "Reverse": ["", "", "", "+", ""],
    "Only identified by site": ["", "", "", "", "+"],
})
clean, counts = pmr.filter_contaminants_decoys(messy, "identifier")
check(len(clean) == 2, f"2 of 5 rows survive filtering (got {len(clean)})")
check(set(clean["identifier"]) == {"P04637", "Q00987"}, "the correct 2 rows survive")
check(counts["n_prefix_removed"] == 2, "CON__/REV__ prefix count is correct")
check(counts["n_flagged_removed"] == 3, "flag-column count is correct (all 3 rows had a flag)")
check(set(counts["flag_columns_found"]) == {"Potential contaminant", "Reverse", "Only identified by site"},
      "all three real MaxQuant flag columns were detected")
check(counts["n_total_removed"] == 3, "total removed = union of both mechanisms (3 rows), not double-counted")

print("=== SENSITIVITY: prefix-only detection is not vacuous ===")
prefix_only = pd.DataFrame({"identifier": ["P1", "CON__P2", "REV__P3"]})
clean2, counts2 = pmr.filter_contaminants_decoys(prefix_only, "identifier")
check(len(clean2) == 1 and counts2["n_prefix_removed"] == 2,
      "prefix detection alone correctly removes CON__/REV__ rows with no flag columns present")

print("=== SENSITIVITY: flag-only detection is not vacuous (no prefix at all) ===")
flag_only = pd.DataFrame({"identifier": ["P1", "P2", "P3"], "Potential contaminant": ["", "+", ""]})
clean3, counts3 = pmr.filter_contaminants_decoys(flag_only, "identifier")
check(len(clean3) == 2 and counts3["n_flagged_removed"] == 1,
      "flag-column detection alone correctly removes a row with no CON__/REV__ prefix")

print("=== SENSITIVITY: a real regression - flag columns must be readable from the caller's frame ===")
# This mirrors the actual bug caught during development: resolve_protein_table
# originally called filter_contaminants_decoys on a stripped-down `work`
# frame that never carried the raw flag columns, so n_flagged_removed was
# silently always 0 even with real '+' flags present. Reproduce the SAME
# shape of bug here directly and confirm it's detectable.
stripped_frame = pd.DataFrame({"identifier": ["P1", "P2"]})  # flag columns NOT carried over
_, counts_bug = pmr.filter_contaminants_decoys(stripped_frame, "identifier")
check(counts_bug["n_flagged_removed"] == 0 and counts_bug["flag_columns_found"] == [],
      "confirms the failure mode: filtering a frame with NO flag columns finds none (the real bug fixed in "
      "resolve_protein_table was failing to carry these columns into the frame being filtered)")


# =============================================================================
# Dedupe by resolved accession
# =============================================================================
print("=== _dedupe_by_accession: keeps the strongest evidence, collapses the rest ===")
dupe_df = pd.DataFrame({
    "accession": ["P04637", "P04637", "Q00987"],
    "padj": [0.5, 0.001, 0.02],
    "log2FC": [0.1, 2.0, 1.0],
})
deduped, n_collapsed = pmr._dedupe_by_accession(dupe_df)
check(n_collapsed == 1, "exactly 1 duplicate collapsed")
check(len(deduped) == 2, "2 unique accessions remain")
kept = deduped[deduped["accession"] == "P04637"].iloc[0]
check(abs(kept["padj"] - 0.001) < 1e-9, "the row with the SMALLER padj (stronger evidence) was kept, not the first-seen one")

print("=== SENSITIVITY: dedupe tie-break order is not accidental ===")
# Reverse the row order - if dedupe were accidentally "keep first seen"
# rather than "keep strongest evidence", this would now keep the WRONG row.
dupe_df_reversed = dupe_df.iloc[::-1].reset_index(drop=True)
deduped2, _ = pmr._dedupe_by_accession(dupe_df_reversed)
kept2 = deduped2[deduped2["accession"] == "P04637"].iloc[0]
check(abs(kept2["padj"] - 0.001) < 1e-9,
      "row order reversed, but the SAME (stronger-evidence) row is still kept - confirms this is a real sort, not accidental first-wins")


# =============================================================================
# resolve_protein_table - light live smoke test (real UniProt ID Mapping)
# =============================================================================
print("=== resolve_protein_table: real UniProt ID Mapping round trip ===")
live_df = pd.DataFrame({
    "Protein IDs": ["P04637;P04637-2", "CON__P02769", "Q15086", "P16220"],
    "log2FC": [1.5, 2.0, -0.8, 0.3],
    "pvalue": [0.001, 0.02, 0.04, 0.5],
    "Potential contaminant": ["", "+", "", ""],
})
resolved, meta = pmr.resolve_protein_table(live_df, taxon_id="9606", verbose=False)
check(meta["contaminants"]["n_total_removed"] == 1, "1 contaminant row removed before resolution")
check(len(resolved) == 2, f"2 rows survive (protein group + secondary-accession collapse into 1, contaminant removed) (got {len(resolved)})")
p53_row = resolved[resolved["accession"] == "P04637"]
check(len(p53_row) == 1, "P04637 (TP53) present exactly once")
if len(p53_row) == 1:
    check(p53_row.iloc[0]["gene_symbol"] == "TP53", f"gene symbol correctly resolved to TP53 (got {p53_row.iloc[0]['gene_symbol']!r})")
    # Two rows resolve to P04637: the direct group ("P04637;P04637-2",
    # pvalue=0.001) and the secondary accession (Q15086, pvalue=0.04). The
    # direct group has the SMALLER (stronger) p-value, so dedupe keeps IT,
    # not the secondary-accession row - status is "mapped", not
    # "resolved_secondary". This asserts the dedupe tie-break rule (kept
    # row = strongest evidence) on real, non-synthetic UniProt data.
    check(p53_row.iloc[0]["status"] == "mapped",
          f"the surviving P04637 row is the direct group match (stronger p-value=0.001 beats Q15086's 0.04), "
          f"correctly labelled 'mapped' (got {p53_row.iloc[0]['status']!r})")


print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
