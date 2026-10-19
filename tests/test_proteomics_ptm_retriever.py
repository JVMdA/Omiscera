# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_ptm_retriever.py, offline (UniProt ID Mapping
monkeypatched). Covers what is UNIQUE to this module - the site-level input
contract that Proteomics_module_retriever.py cannot express at all (no
residue/position anywhere in that module's RESOLVED_COLUMNS).

THE single most important test in this file: isoform preservation. Every
other Proteomics retriever collapses "P04637-2" to "P04637" via
strip_isoform, because GO/Reactome/InterPro annotate the canonical entry.
This module must NOT do that for site IDENTITY - see the module's own
docstring for why silently merging isoform-specific sites is the one
mistake this whole module exists to avoid.
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Proteomics_ptm_retriever as pr  # noqa: E402

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
print("=== normalize_site: every accepted input shape ===")
# =============================================================================
check(pr.normalize_site(raw_site="S473") == ("S", 473), "bare 1-letter+digits ('S473')")
check(pr.normalize_site(raw_site="Ser473") == ("S", 473), "3-letter code ('Ser473')")
check(pr.normalize_site(raw_site="pS473") == ("S", 473), "phospho-notation prefix ('pS473')")
check(pr.normalize_site(raw_site="s473") == ("S", 473), "lowercase residue letter")
check(pr.normalize_site(raw_site="K652") == ("K", 652), "non-STY residue (lysine, e.g. acetylation/ubiquitination)")
check(pr.normalize_site(residue="S", position="473") == ("S", 473), "separate residue+position columns")
check(pr.normalize_site(residue="S", position="473;450") == ("S", 473),
      "MaxQuant's semicolon-joined multi-protein position - only the FIRST is used, "
      "matching the leading-accession convention for protein groups")
check(pr.normalize_site(raw_site="bogus") == (None, None), "unparseable token returns (None, None)")
check(pr.normalize_site(raw_site="Xyz99") == (None, None), "unknown 3-letter code returns (None, None)")
check(pr.normalize_site(residue=None, position=None) == (None, None), "no input at all returns (None, None)")

# =============================================================================
print("=== detect_columns: site XOR residue+position, both accepted shapes ===")
# =============================================================================
combined_df = pd.DataFrame({"Protein IDs": ["P1"], "Site": ["S1"], "log2FC": [1.0], "pvalue": [0.01]})
mapping = pr.detect_columns(combined_df)
check(mapping.get("site") == "Site", f"combined site column detected (got {mapping})")

split_df = pd.DataFrame({"Proteins": ["P1"], "Amino acid": ["S"], "Positions within proteins": ["1"],
                         "log2FC": [1.0], "pvalue": [0.01]})
mapping2 = pr.detect_columns(split_df)
check(mapping2.get("residue") == "Amino acid" and mapping2.get("position") == "Positions within proteins",
      f"MaxQuant-style split residue+position columns detected (got {mapping2})")

missing_df = pd.DataFrame({"Protein IDs": ["P1"], "log2FC": [1.0], "pvalue": [0.01]})
try:
    pr.detect_columns(missing_df)
    check(False, "a table with NEITHER site NOR residue+position raises")
except ValueError as exc:
    check("site" in str(exc), f"the error names the missing site/residue+position requirement (got: {exc})")

# =============================================================================
print("=== resolve_site_table: ISOFORM PRESERVATION - the headline test ===")
# =============================================================================
_real_idmapping = pr._idmapping_resolve
pr._idmapping_resolve = lambda accs, verbose=False: (
    {a.upper(): {"accession": a.upper(), "gene_symbol": "TP53" if a.upper() == "P04637" else "X", "protein_name": "t"}
     for a in accs}, 0,
)

iso_df = pd.DataFrame({
    "Protein IDs": ["P04637", "P04637-2", "Q00987"],
    "Site": ["S315", "S315", "S166"],
    "log2FC": [2.41, 1.9, 1.12],
    "pvalue": [0.001, 0.002, 0.03],
})
resolved, meta = pr.resolve_site_table(iso_df, taxon_id="9606", localization_threshold=None)
site_ids = set(resolved["site_id"])
check("P04637_S315" in site_ids and "P04637-2_S315" in site_ids,
      f"the canonical AND the isoform-specific site at the SAME residue number both survive as "
      f"DISTINCT rows, not collapsed into one (got {sorted(site_ids)})")
iso_row = resolved[resolved["accession"] == "P04637-2"].iloc[0]
check(iso_row["canonical_accession"] == "P04637" and iso_row["isoform"] == "2",
      f"the isoform row's canonical_accession strips the suffix (for external-database matching only) "
      f"while accession itself keeps it (got accession={iso_row['accession']!r}, "
      f"canonical={iso_row['canonical_accession']!r}, isoform={iso_row['isoform']!r})")
check(iso_row["canonical_site_id"] == "P04637_S315",
      f"canonical_site_id collapses to the SAME key the canonical row also uses - this is the "
      f"deliberate, SEPARATE key for matching against canonical-only external databases (iPTMnet), "
      f"never used as this row's own identity (got {iso_row['canonical_site_id']!r})")

# =============================================================================
print("=== resolve_site_table: protein-group ambiguity flag ===")
# =============================================================================
group_df = pd.DataFrame({
    "Protein IDs": ["P04637;P04637-2", "Q00987"],
    "Site": ["S315", "S166"],
    "log2FC": [2.0, 1.0], "pvalue": [0.01, 0.02],
})
resolved_g, _ = pr.resolve_site_table(group_df, taxon_id="9606", localization_threshold=None)
check(resolved_g.loc[resolved_g["protein_group"] == "P04637;P04637-2", "ambiguous_group"].iloc[0] == True,
      "a multi-member protein group is flagged ambiguous_group=True")
check(resolved_g.loc[resolved_g["protein_group"] == "Q00987", "ambiguous_group"].iloc[0] == False,
      "a single-member group is NOT flagged ambiguous")
check(resolved_g.loc[resolved_g["protein_group"] == "P04637;P04637-2", "accession"].iloc[0] == "P04637",
      "the LEADING member of the group is used as the row's accession")

drop_resolved, drop_meta = pr.resolve_site_table(group_df, taxon_id="9606", localization_threshold=None,
                                                 drop_ambiguous_groups=True)
check(len(drop_resolved) == 1 and drop_resolved.iloc[0]["protein_group"] == "Q00987",
      f"drop_ambiguous_groups=True removes the ambiguous row entirely (got {len(drop_resolved)} rows)")
check(drop_meta["dropped_ambiguous_groups"] is True, "meta reports the drop setting honestly")

# =============================================================================
print("=== resolve_site_table: multiplicity collapse, NEVER averaged ===")
# =============================================================================
multi_df = pd.DataFrame({
    "Protein IDs": ["P04637", "P04637", "P04637"],
    "Site": ["S315", "S315", "S315"],
    "log2FC": [5.0, 1.0, 3.0],  # a mean of these (3.0) must NEVER appear
    "pvalue": [0.5, 0.001, 0.3],
    "Multiplicity": [1, 2, 3],
})
multi_resolved, multi_meta = pr.resolve_site_table(multi_df, taxon_id="9606", localization_threshold=None)
check(len(multi_resolved) == 1, f"three multiplicity-class rows for the SAME site collapse to exactly one (got {len(multi_resolved)})")
check(abs(multi_resolved.iloc[0]["log2FC"] - 1.0) < 1e-9,
      f"the collapsed row keeps the SMALLEST-pvalue row's log2FC (1.0), not an average of 5.0/1.0/3.0 "
      f"(got {multi_resolved.iloc[0]['log2FC']})")
check(multi_meta["n_duplicates_collapsed"] == 2, f"2 collapsed rows reported (got {multi_meta['n_duplicates_collapsed']})")

# =============================================================================
print("=== resolve_site_table: localization-probability filter ===")
# =============================================================================
loc_df = pd.DataFrame({
    "Protein IDs": ["P04637", "Q00987"],
    "Site": ["S315", "S166"],
    "log2FC": [2.0, 1.0], "pvalue": [0.01, 0.02],
    "Localization prob": [0.9, 0.5],
})
loc_resolved, loc_meta = pr.resolve_site_table(loc_df, taxon_id="9606", localization_threshold=0.75)
check(len(loc_resolved) == 1 and loc_resolved.iloc[0]["accession"] == "P04637",
      f"the sub-threshold row (0.5 < 0.75) is dropped, the above-threshold one (0.9) survives "
      f"(got {len(loc_resolved)} rows: {list(loc_resolved['accession'])})")
check(loc_meta["n_dropped_localization"] == 1 and loc_meta["has_localization_prob"] is True,
      f"meta reports the drop count and that the column was present (got {loc_meta['n_dropped_localization']}, "
      f"{loc_meta['has_localization_prob']})")

no_loc_resolved, no_loc_meta = pr.resolve_site_table(iso_df, taxon_id="9606", localization_threshold=0.75)
check(len(no_loc_resolved) == 3 and no_loc_meta["has_localization_prob"] is False,
      f"with NO localization_prob column at all, the filter is a no-op (nothing dropped) and "
      f"has_localization_prob is honestly False, not silently True (got {len(no_loc_resolved)} rows, "
      f"has_localization_prob={no_loc_meta['has_localization_prob']})")

# =============================================================================
print("=== resolve_site_table: unparseable sites are dropped and counted, not crashed on ===")
# =============================================================================
bad_site_df = pd.DataFrame({
    "Protein IDs": ["P04637", "Q00987"],
    "Site": ["S315", "not_a_site"],
    "log2FC": [2.0, 1.0], "pvalue": [0.01, 0.02],
})
bad_resolved, bad_meta = pr.resolve_site_table(bad_site_df, taxon_id="9606", localization_threshold=None)
check(len(bad_resolved) == 1 and bad_meta["n_unparsed_site"] == 1,
      f"the unparseable site is dropped and counted (got {len(bad_resolved)} rows, "
      f"n_unparsed_site={bad_meta['n_unparsed_site']})")

# =============================================================================
print("=== resolve_site_table: contaminant/decoy rows still filtered (reused from protein-level) ===")
# =============================================================================
contam_df = pd.DataFrame({
    "Protein IDs": ["CON__P123", "P04637"],
    "Site": ["S1", "S315"],
    "log2FC": [5.0, 2.0], "pvalue": [0.01, 0.01],
})
contam_resolved, contam_meta = pr.resolve_site_table(contam_df, taxon_id="9606", localization_threshold=None)
check(len(contam_resolved) == 1 and contam_meta["contaminants"]["n_prefix_removed"] == 1,
      f"the CON__-prefixed row is removed before resolution (got {len(contam_resolved)} rows)")

pr._idmapping_resolve = _real_idmapping

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
