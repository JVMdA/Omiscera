# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_ptm.py, offline (iPTMnet fetch monkeypatched to
a small fixture file - see FAKE_PTM_ROWS below). Covers what is UNIQUE to
this module - the shared regulon numerics (_ulm_score/_null_slopes_shared/
_permutation_p_value/_empirical_null_calibration) are already covered by
test_regulon_engine.py and are imported here unchanged, never retested.

THE HEADLINE CLAIM this file demonstrates end to end: a kinase whose known
substrate sites are systematically UP-regulated gets detected as ACTIVE
(positive score, small calibrated FDR), a kinase whose sites are noise does
not - on a synthetic fixture with a KNOWN, planted ground truth, not just
"the code ran without crashing".
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Proteomics_ptm as ptm  # noqa: E402

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
# Fixture: a tiny iPTMnet ptm.txt, no header row, real column order.
# KIN1 phosphorylates 6 sites on proteins P00000-P00005 (will be made
# strongly UP in the synthetic resolved table below).
# KIN2 phosphorylates 6 DIFFERENT sites on P00006-P00011 (will be noise).
# KIN3's substrate is on a DIFFERENT organism (must be excluded).
# KIN4's relationship is ACETYLATION, not phosphorylation (must be excluded).
# =============================================================================
# DELIBERATELY NOT a real taxon id (e.g. "9606"). build_kinase_substrate_universe's
# parsed-JSON cache (cache/ora_universe/iPTMnet_v{version}_{taxon_id}.json) is
# SHARED with real production runs and keyed only on taxon_id - using the real
# human id here with force_refresh=True would silently overwrite a real user's
# cached human kinase universe with this test's 2-kinase fixture. Bug found and
# fixed 2026-09-21 after exactly this happened during interactive testing.
TAX = "0002"
FIXTURE_ROWS = []
for i in range(6):
    FIXTURE_ROWS.append(("PHOSPHORYLATION", "test", f"P{i:05d}", f"GENE{i}", "Homo sapiens (Human)",
                        f"S{100+i}", "KIN1AC", "KIN1", "note", "1"))
for i in range(6, 12):
    FIXTURE_ROWS.append(("PHOSPHORYLATION", "test", f"P{i:05d}", f"GENE{i}", "Homo sapiens (Human)",
                        f"S{100+i}", "KIN2AC", "KIN2", "note", "1"))
FIXTURE_ROWS.append(("PHOSPHORYLATION", "test", "P99999", "GENEX", "Mus musculus (Mouse)",
                    "S1", "KIN3AC", "KIN3", "note", "1"))
FIXTURE_ROWS.append(("ACETYLATION", "test", "P00000", "GENE0", "Homo sapiens (Human)",
                    "K1", "KIN4AC", "KIN4", "note", "1"))

_TMP_DIR = tempfile.mkdtemp(prefix="omiscera_ptm_test_")
FIXTURE_PATH = os.path.join(_TMP_DIR, "fake_iptmnet_ptm.txt")
with open(FIXTURE_PATH, "w", encoding="utf-8") as fh:
    for r in FIXTURE_ROWS:
        fh.write("\t".join(r) + "\n")

_real_ensure = ptm._ensure_iptmnet_file
ptm._ensure_iptmnet_file = lambda url, local_path, refresh_days=90, verbose=False: (
    FIXTURE_PATH if local_path.endswith("iptmnet_ptm.txt") else None
)


def _make_resolved_sites(n_background=40, seed=0):
    """6 KIN1 sites strongly UP, 6 KIN2 sites at noise level, plus a
    background of unrelated sites so the ranked list has enough entries for
    a meaningful null."""
    rows = []
    for i in range(6):
        rows.append({
            "input_id": f"r{i}", "accession": f"P{i:05d}", "canonical_accession": f"P{i:05d}",
            "isoform": "", "protein_group": f"P{i:05d}", "ambiguous_group": False,
            "gene_symbol": f"GENE{i}", "protein_name": "t", "residue": "S", "position": 100 + i,
            "site_id": f"P{i:05d}_S{100+i}", "canonical_site_id": f"P{i:05d}_S{100+i}",
            "multiplicity": np.nan, "localization_prob": np.nan,
            "log2FC": 3.0, "pvalue": 1e-6, "padj": 1e-5, "status": "mapped",
        })
    for i in range(6, 12):
        rows.append({
            "input_id": f"r{i}", "accession": f"P{i:05d}", "canonical_accession": f"P{i:05d}",
            "isoform": "", "protein_group": f"P{i:05d}", "ambiguous_group": False,
            "gene_symbol": f"GENE{i}", "protein_name": "t", "residue": "S", "position": 100 + i,
            "site_id": f"P{i:05d}_S{100+i}", "canonical_site_id": f"P{i:05d}_S{100+i}",
            "multiplicity": np.nan, "localization_prob": np.nan,
            "log2FC": 0.01, "pvalue": 0.9, "padj": 0.95, "status": "mapped",
        })
    rng = np.random.default_rng(seed)
    for i in range(n_background):
        rows.append({
            "input_id": f"bg{i}", "accession": f"B{i:05d}", "canonical_accession": f"B{i:05d}",
            "isoform": "", "protein_group": f"B{i:05d}", "ambiguous_group": False,
            "gene_symbol": f"BGGENE{i}", "protein_name": "t", "residue": "S", "position": 200 + i,
            "site_id": f"B{i:05d}_S{200+i}", "canonical_site_id": f"B{i:05d}_S{200+i}",
            "multiplicity": np.nan, "localization_prob": np.nan,
            "log2FC": float(rng.normal(0, 0.3)), "pvalue": float(rng.uniform(0.1, 0.9)),
            "padj": float(rng.uniform(0.2, 0.95)), "status": "mapped",
        })
    return pd.DataFrame(rows)


RESOLVED_SITES = _make_resolved_sites()

# =============================================================================
print("=== build_kinase_substrate_universe: organism + ptm_type filters ===")
# =============================================================================
universe = ptm.build_kinase_substrate_universe(TAX, taxon_name="Homo sapiens", force_refresh=True)
check(set(universe["terms"].keys()) == {"KIN1AC", "KIN2AC"},
      f"only the human, phosphorylation-type kinases survive - KIN3 (mouse) and KIN4 (acetylation) "
      f"are excluded (got {sorted(universe['terms'].keys())})")
check(universe["terms"]["KIN1AC"]["genes"] == {f"P{i:05d}_S{100+i}" for i in range(6)},
      f"KIN1's substrate set is exactly its 6 real sites (got {universe['terms']['KIN1AC']['genes']})")
check(universe["terms"]["KIN1AC"]["name"] == "KIN1", "the display name is the enzyme_genename, not the accession")

no_name = ptm.build_kinase_substrate_universe("0000", taxon_name="")
check(no_name["terms"] == {} and no_name["caveat"] is not None,
      f"missing taxon_name is refused with a caveat, not silently empty-and-unexplained "
      f"(got caveat={no_name['caveat']!r})")

wrong_organism = ptm.build_kinase_substrate_universe("0001", taxon_name="Danio rerio", force_refresh=True)
check(wrong_organism["terms"] == {} and wrong_organism["caveat"] is not None,
      f"an organism with NO matching rows in the fixture gets an empty universe + caveat, not an "
      f"exception (got caveat={wrong_organism['caveat']!r})")

# =============================================================================
print("=== build_ranked_site_list: site_id vs canonical_site_id, descending order ===")
# =============================================================================
ranked = ptm.build_ranked_site_list(RESOLVED_SITES, metric="log2fc")
check(ranked["n_sites"] == len(RESOLVED_SITES), f"every resolved site is ranked (got {ranked['n_sites']})")
check(ranked["site_ids"][0] in {f"P{i:05d}_S{100+i}" for i in range(6)},
      f"the top-ranked (highest log2FC) site is one of KIN1's strongly-up sites (got {ranked['site_ids'][0]})")
check(ranked["canonical_site_ids"] == ranked["site_ids"],
      "with no isoforms in this fixture, site_ids and canonical_site_ids coincide (sanity check)")

# =============================================================================
print("=== run_kinase_activity: THE HEADLINE CLAIM - active kinase detected, noise kinase is not ===")
# =============================================================================
result = ptm.run_kinase_activity(RESOLVED_SITES, TAX, taxon_name="Homo sapiens", metric="log2fc",
                                 min_kinase_sites=3, n_permutations=200, seed=42, force_refresh=True)
check(result["table"] is not None, f"a real result table was produced (note={result.get('note')})")
if result["table"] is not None:
    tbl = result["table"]
    check(len(tbl) == 2, f"exactly KIN1 and KIN2 were tested (got {len(tbl)})")
    kin1 = tbl[tbl["Kinase"] == "KIN1"].iloc[0]
    kin2 = tbl[tbl["Kinase"] == "KIN2"].iloc[0]
    check(kin1["Activity Score"] > kin2["Activity Score"] and kin1["Activity Score"] > 5,
          f"KIN1 (substrates strongly UP) scores meaningfully higher than KIN2 (noise) "
          f"(KIN1={kin1['Activity Score']:.3f}, KIN2={kin2['Activity Score']:.3f})")
    check(kin1["FDR (calibrated)"] < 0.05,
          f"KIN1 clears a real significance threshold on its calibrated FDR (got {kin1['FDR (calibrated)']:.4g})")
    check(kin1["Substrate Sites Tested"] == 6, f"KIN1 tested on all 6 of its matched sites (got {kin1['Substrate Sites Tested']})")
    check("KIN1" in result["kinase_regression_detail"], "per-kinase regression detail stored for the top kinase")
    check(set(result["kinase_regression_detail"]["KIN1"]["substrate_site_ids"]) == {f"P{i:05d}_S{100+i}" for i in range(6)},
          "the stored detail's substrate_site_ids are this kinase's OWN sites (site_id, not canonical_site_id)")

# min_kinase_sites gate
too_strict = ptm.run_kinase_activity(RESOLVED_SITES, TAX, taxon_name="Homo sapiens", metric="log2fc",
                                     min_kinase_sites=100, n_permutations=50, seed=1)
check(too_strict["table"].empty and too_strict["note"] is not None,
      f"a min_kinase_sites gate no kinase can clear returns an explained empty result (empty DataFrame, "
      f"not None - the house convention), not a crash (got note={too_strict['note']!r})")

# too few sites overall
tiny = _make_resolved_sites(n_background=0).head(2)
degenerate = ptm.run_kinase_activity(tiny, TAX, taxon_name="Homo sapiens", min_kinase_sites=1, n_permutations=50)
check(degenerate["table"].empty and "Too few tested sites" in (degenerate["note"] or ""),
      f"fewer than 3 tested sites is refused with an explanation (got note={degenerate['note']!r})")

# =============================================================================
print("=== run_kinase_substrate_ora: KIN1's sites are enriched among up-regulated sites ===")
# =============================================================================
ora = ptm.run_kinase_substrate_ora(RESOLVED_SITES, TAX, taxon_name="Homo sapiens",
                                   padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, force_refresh=True)
check("Kinase substrates (up)" in ora["tables"], "an 'up' direction table exists")
up_tbl = ora["tables"]["Kinase substrates (up)"]["table"]
check(not up_tbl.empty and up_tbl.iloc[0]["Term ID"] == "KIN1AC",
      f"KIN1 is the (only) enriched kinase among up-regulated sites (got "
      f"{up_tbl.to_dict('records') if not up_tbl.empty else None})")
down_tbl = ora["tables"]["Kinase substrates (down)"]["table"]
check(down_tbl.empty, f"no kinase is enriched among down-regulated sites (none exist in this fixture) (got {len(down_tbl)} rows)")

# =============================================================================
print("=== select_significant_sites ===")
# =============================================================================
sig = ptm.select_significant_sites(RESOLVED_SITES, padj_cutoff=0.05, lfc_cutoff=1.0, direction="up")
check(len(sig) == 6 and set(sig["site_id"]) == {f"P{i:05d}_S{100+i}" for i in range(6)},
      f"exactly the 6 planted up-regulated sites are selected (got {len(sig)})")

ptm._ensure_iptmnet_file = _real_ensure

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
