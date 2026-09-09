"""LIVE network test for KEGG PATHWAY in Genomics ORA + Gene-set Burden Test:
build_kegg_universe and build_kegg_gene_map (Genomics_module.py), and
Genomics_module_retriever.resolve_kegg_ids (the retriever-side gene mapping).

Does NOT need the app's own uvicorn server - only real internet access to
KEGG's REST API (rest.kegg.jp) and UniProt/NCBI (via the real retriever),
same pattern as smoke_regulon_fetchers.py and live_msigdb_urls.py.

WHY THE UNIVERSE NUMBERS HERE DIFFER FROM AN EARLIER VERSION OF THIS FILE:
build_kegg_universe used to be keyed on UniProt accessions. It is now keyed on
KEGG GENE IDs (see that function's docstring) - accession space inflated
counts via UniProt fan-out (one KEGG gene can have several accessions) and
could never reach the ~9% of pathway genes with no UniProt cross-reference at
all. This is a deliberate, verified improvement, not a regression: human
N drops 10,873 -> 9,424 and hsa00010's K drops 107 -> 67, both because they
now count real genes once each instead of accessions.

WHAT THIS FILE PROVES, beyond the offline tests in
test_kegg_universe_parsing.py / test_kegg_gene_map.py / test_kegg_retriever_mapping.py:
  - The whole chain (KEGG's real REST responses -> parsing -> universe/gene-map
    assembly -> the retriever's real UniProt/NCBI resolution -> the app's own
    Fisher/Wilcoxon engines) produces the correct biological answer on a KNOWN
    input, not just a well-shaped one.
  - KEGG's non-model-organism coverage is real, not just claimed (chickpea),
    AND that the retriever's three-bridge gene mapping genuinely recovers
    genes UniProt alone could not place (the ARF2/ARF6 case).
  - The real failure paths (an organism KEGG has no data for at all; a real
    organism with pathways but no NCBI Entrez cross-reference) degrade to an
    honest caveat/flag instead of raising or silently caching garbage.
"""
import os
import random
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm
import Genomics_module_retriever as gr
import database_catalog as dc

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def _clear_cache(cache_db, org):
    path = gm._universe_cache_file_path(cache_db, org)
    try:
        os.remove(path)
    except OSError:
        pass
    return path


print("=== build_kegg_universe('hsa'): real human KEGG PATHWAY build, KEGG-gene-ID space ===")
cache_path = _clear_cache(f"KEGG_v{gm.KEGG_UNIVERSE_CACHE_VERSION}", "hsa")
t0 = time.time()
hsa = gm.build_kegg_universe("hsa", verbose=True, force_refresh=True)
elapsed = time.time() - t0
check(330 <= len(hsa["pathways"]) <= 400,
      f"human pathway count in the expected real range (got {len(hsa['pathways'])}, measured 372)")
check(8500 <= hsa["universe_size"] <= 10500,
      "universe_size counts real KEGG GENES once each (measured 9,424), not UniProt accessions "
      f"(the old, superseded number was 10,873) (got {hsa['universe_size']})")
check(hsa.get("caveat") is None, f"a real, successful human build carries no caveat (got {hsa.get('caveat')!r})")
check(os.path.exists(cache_path), "a successful build writes its cache file to disk")
print(f"  (cold build took {elapsed:.1f}s)")

t1 = time.time()
hsa_cached = gm.build_kegg_universe("hsa", force_refresh=False)
cached_elapsed = time.time() - t1
check(cached_elapsed < 0.5, f"a cached re-read is fast, not a re-fetch (got {cached_elapsed:.3f}s)")
sample_pid = next(iter(hsa_cached["pathways"]))
check(isinstance(hsa_cached["pathways"][sample_pid]["genes"], set),
      "a cache-loaded universe's 'genes' values are real sets, not JSON lists")
sample_gene = next(iter(hsa_cached["pathways"][sample_pid]["genes"]))
check(sample_gene.startswith("hsa:"),
      f"pathway members are KEGG gene IDs (hsa:NNNN), not UniProt accessions (got {sample_gene!r})")

print("=== build_kegg_universe('cam'): non-model organism (chickpea) has real data ===")
cam = gm.build_kegg_universe("cam", verbose=False, force_refresh=True)
check(len(cam["pathways"]) > 100, f"chickpea has real KEGG pathway coverage (got {len(cam['pathways'])}, measured 163)")
check(cam["universe_size"] > 5000, f"chickpea universe is non-trivial (got {cam['universe_size']})")
check(cam.get("caveat") is None, "chickpea build carries no caveat")

print("=== build_kegg_universe('vg'): real total-failure path degrades honestly ===")
# Verified live: KEGG's "vg" pseudo-organism genuinely has list/pathway -> 200
# with 0 rows, AND link/pathway -> HTTP 400 - a real organism code with no
# usable pathway data, not a fabricated edge case.
vg_cache_path = _clear_cache(f"KEGG_v{gm.KEGG_UNIVERSE_CACHE_VERSION}", "vg")
vg = gm.build_kegg_universe("vg", verbose=True, force_refresh=True)
check(vg["pathways"] == {}, f"an organism with no KEGG pathway data returns an empty universe, not an exception (got {len(vg['pathways'])} pathways)")
check(bool(vg.get("caveat")), "the empty result carries a non-empty caveat explaining why")
check(not os.path.exists(vg_cache_path), "a failed/empty build is NOT cached (no TTL on this layer)")

print("=== build_kegg_gene_map: human vs chickpea - the measured symbol-coverage asymmetry, locked in ===")
hsa_map = gm.build_kegg_gene_map("hsa", force_refresh=True)
cam_map = gm.build_kegg_gene_map("cam", force_refresh=True)
hsa_symbol_ratio = len(hsa_map["symbol_to_gene"]) / max(len(hsa_map["entrez_to_gene"]), 1)
cam_symbol_ratio = len(cam_map["symbol_to_gene"]) / max(len(cam_map["entrez_to_gene"]), 1)
check(hsa_symbol_ratio > 0.80,
      f"human: symbol coverage is high (measured ~94%), the symbol bridge is genuinely useful there "
      f"(got ratio={hsa_symbol_ratio:.3f})")
check(cam_symbol_ratio < 0.10,
      f"chickpea: symbol coverage is near-nothing (measured ~1.6-1.9%) - locks in why the symbol bridge "
      f"is a low-priority LAST resort, not the primary route, for non-model organisms (got ratio={cam_symbol_ratio:.3f})")
check(hsa_map["has_entrez"] and cam_map["has_entrez"], "both organisms have a real Entrez cross-reference")

print("=== build_kegg_gene_map('syn'): a REAL organism with pathways but NO Entrez cross-reference ===")
# Verified live (independently of this suite): Synechocystis (syn) has 128
# real pathways and 3,564 real UniProt conversions, but conv/ncbi-geneid
# returns HTTP 400 for it - a genuine, not fabricated, has_entrez=False case.
syn_cache_path = _clear_cache(f"KEGG_GENEMAP_v{gm.KEGG_GENEMAP_CACHE_VERSION}", "syn")
syn_map = gm.build_kegg_gene_map("syn", verbose=True, force_refresh=True)
check(syn_map["has_entrez"] is False, f"syn genuinely has no Entrez cross-reference (got has_entrez={syn_map['has_entrez']})")
check(len(syn_map["acc_to_gene"]) > 500, f"syn's accession bridge still works fine independently (got {len(syn_map['acc_to_gene'])})")
check(os.path.exists(syn_cache_path), "a build with a missing OPTIONAL map is still cached (has_entrez=False is a real, permanent property)")

print("=== Entrez identity guard: KEGG gene IDs really are Entrez IDs for a RefSeq-sourced organism ===")
sample_pairs = list(cam_map["entrez_to_gene"].items())[:200]
matches = sum(1 for entrez, gene in sample_pairs if gene.split(":", 1)[1] == entrez)
check(len(sample_pairs) > 0 and matches / len(sample_pairs) >= 0.95,
      f"at least 95% of sampled chickpea genes have identical numeric KEGG-gene-ID/Entrez-ID parts "
      f"(got {matches}/{len(sample_pairs)}) - the assumption this whole feature leans on, checked live")

print("=== database_catalog.search_contexts('KEGG', ...) still resolves 'hsa' (the picker this feature reuses) ===")
frame = dc.search_contexts("KEGG", "Homo sapiens", limit=5, allow_online=True)
ids = list(frame["context_id"]) if not frame.empty else []
check("hsa" in ids, f"searching 'Homo sapiens' finds the hsa organism code (got {ids})")

print("=== Non-model known-answer (headline): real chickpea gene symbols through the REAL retriever ===")
# The exact case that motivated this whole feature: ARF2/ARF6 have NO
# reference proteome for chickpea (a non-model organism), so UniProt's own
# tie-break cannot disambiguate them and they come back "ambiguous" with an
# empty accession under resolve_gene_table alone.
import pandas as pd
CHICKPEA_SYMBOLS = ["ARF1", "NAC01", "UGT85J2", "HSFA4C", "ARF2", "ARF6", "ARF3"]
genes_df = pd.DataFrame({"Gene Symbol": CHICKPEA_SYMBOLS})
resolved = gr.resolve_gene_table(genes_df, "gene_symbol", taxon_id="3827", verbose=False)

still_ambiguous_before = resolved[resolved["Resolution_status"] == "ambiguous"]
check(len(still_ambiguous_before) >= 1,
      f"at least one real chickpea gene is genuinely ambiguous under UniProt resolution alone "
      f"(got {len(still_ambiguous_before)}: {still_ambiguous_before['Gene Symbol'].tolist()}) - "
      "confirms this test exercises the real problem, not an already-solved input")

mapped = gr.resolve_kegg_ids(resolved, "Gene Symbol", "cam", verbose=False)
still_ambiguous_after = mapped[mapped["Resolution_status"] == "ambiguous"]
check(len(still_ambiguous_after) < len(still_ambiguous_before),
      f"at least one previously-ambiguous gene was recovered by KEGG mapping "
      f"(before: {len(still_ambiguous_before)}, after: {len(still_ambiguous_after)})")
check((mapped["KEGG ID"] != "").sum() == len(mapped),
      f"EVERY one of the {len(mapped)} real chickpea genes now has a KEGG ID - the headline result "
      f"(got {(mapped['KEGG ID'] != '').sum()}/{len(mapped)})")
for _, row in mapped.iterrows():
    bridge = ("accession" if row["Resolution_status"] != "resolved_via_kegg"
              else ("disambiguation/entrez" if row["KEGG ID"] else "unmatched"))
    print(f"  {row['Gene Symbol']}: KEGG ID={row['KEGG ID']!r}, status={row['Resolution_status']!r}, bridge~={bridge}")

print("=== Known-answer check: 10 canonical human glycolytic enzymes, through the REAL retriever ===")
GLYCOLYTIC_SYMBOLS = ["HK1", "GPI", "PFKM", "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM"]
glyco_df = pd.DataFrame({"Gene Symbol": GLYCOLYTIC_SYMBOLS})
glyco_resolved = gr.resolve_gene_table(glyco_df, "gene_symbol", taxon_id="9606", verbose=False)
glyco_mapped = gr.resolve_kegg_ids(glyco_resolved, "Gene Symbol", "hsa", verbose=False)
check((glyco_mapped["KEGG ID"] != "").all(),
      f"all 10 glycolytic genes got a real KEGG ID via the retriever (got "
      f"{(glyco_mapped['KEGG ID'] != '').sum()}/10)")

kegg_id_list = glyco_mapped["KEGG ID"].tolist()
kegg_id_labels = dict(zip(glyco_mapped["KEGG ID"], glyco_mapped["Gene Symbol"]))
ora = gm.run_ora_for_database(
    "KEGG", "hsa", [], list_kegg_ids=kegg_id_list, kegg_id_labels=kegg_id_labels,
    min_term_size=3, max_term_size=500,
)
kegg_block = ora["KEGG"]
check(kegg_block["n_input_dropped"] == 0,
      f"ALL 10 real, retriever-resolved KEGG gene IDs are present in the KEGG universe "
      f"(got n_input_dropped={kegg_block['n_input_dropped']})")

table = kegg_block["table"]
top = table.iloc[0]
check(top["Term ID"] == "hsa00010", f"hsa00010 (Glycolysis) ranks #1 by FDR (got top hit {top['Term ID']!r})")
check(top["Term Name"] == "Glycolysis / Gluconeogenesis",
      f"suffix stripper removed ' - Homo sapiens (human)' without eating the internal ' / ' (got {top['Term Name']!r})")
check(float(top["FDR"]) < 1e-10, f"top hit is overwhelmingly significant (got FDR={top['FDR']!r})")
check(int(top["k"]) >= 9, f"at least 9 of the 10 input genes fall in hsa00010 (got k={top['k']})")
check(50 <= int(top["K"]) <= 90,
      f"hsa00010's REAL gene-count term size (measured 67, down from the superseded accession-space "
      f"107) (got K={top['K']})")
check(set(str(top["Overlapping Genes"]).split(", ")) == set(GLYCOLYTIC_SYMBOLS[:int(top['k'])]) or
      all(g in GLYCOLYTIC_SYMBOLS for g in str(top["Overlapping Genes"]).split(", ")),
      f"Overlapping Genes shows the user's OWN gene symbols (via kegg_id_labels), never a raw "
      f"KEGG gene ID like 'hsa:3098' (got {top['Overlapping Genes']!r})")

term_ids = set(table["Term ID"])
check("hsa01100" not in term_ids,
      f"hsa01100 (Metabolic pathways, a KEGG overview map) is EXCLUDED from the results table - "
      f"max_term_size handles it rather than a build-time exclusion (n_terms_excluded_min_size={kegg_block['n_terms_excluded_min_size']})")
check(kegg_block["n_terms_excluded_min_size"] > 0, "at least one term (the overview map) was excluded for exceeding max_term_size")

print("=== Cross-database sanity: the SAME 10 (accession-resolved) genes also top Reactome's own Glycolysis pathway ===")
glyco_accessions = [a for a in glyco_resolved["UniProt ID"].tolist() if a]
reactome = gm.run_ora_for_database(
    "Reactome", "9606", glyco_accessions, taxon_name="Homo sapiens",
    min_term_size=3, max_term_size=500,
)
reactome_top = reactome["Reactome"]["table"].iloc[0]
check("70171" in str(reactome_top["Term ID"]) or "glycolysis" in str(reactome_top["Term Name"]).lower(),
      "Reactome, an INDEPENDENT database in a DIFFERENT ID space, also ranks Glycolysis on the same "
      f"genes - agreement across independent spaces is what proves this is right, not accidentally empty "
      f"(got {reactome_top['Term ID']!r}, {reactome_top['Term Name']!r})")

print("=== Burden Test known-answer: same 10 glycolytic KEGG gene IDs scored high vs. sampled low ===")
random.seed(0)
universe_genes = set()
for info in gm.build_kegg_universe("hsa")["pathways"].values():
    universe_genes |= info["genes"]
background_sample = random.sample(sorted(universe_genes - set(kegg_id_list)), 500)
kegg_scores = {g: 1.0 for g in kegg_id_list}
kegg_scores.update({g: 0.0 for g in background_sample})

burden = gm.run_burden_for_database(
    "KEGG", "hsa", {}, score_by_kegg_id=kegg_scores, kegg_id_labels=kegg_id_labels,
    min_term_size=3, max_term_size=500,
)
burden_table = burden["KEGG"]["table"]
glyco_row = burden_table[burden_table["Term ID"] == "hsa00010"]
check(not glyco_row.empty, "hsa00010 appears in the Burden Test result table")
if not glyco_row.empty:
    row = glyco_row.iloc[0]
    check(float(row["FDR"]) < 0.01, f"hsa00010 is significant in the burden test (got FDR={row['FDR']!r})")
    check(float(row["Effect (mean difference)"]) > 0,
          f"hsa00010's effect direction is positive (higher-scored genes concentrate there) (got {row['Effect (mean difference)']!r})")
    # hsa00010 (K=67 real genes) has far more members than the 10 labeled
    # glycolytic genes plus the 500 randomly sampled (deliberately unlabeled)
    # background genes cover - most term members legitimately fall back to a
    # raw "hsa:NNNN" ID here, which is correct, not a bug. The real assertion
    # is narrower: every one of the 10 KNOWN, labeled genes must appear by
    # its real symbol, not by its raw KEGG ID.
    burden_gene_tokens = set(str(row["Genes"]).split(", "))
    labeled_present = {g for g in GLYCOLYTIC_SYMBOLS if g in burden_gene_tokens}
    check(len(labeled_present) >= 9,
          f"at least 9 of the 10 labeled glycolytic genes appear in Genes BY THEIR REAL SYMBOL, not a "
          f"raw KEGG gene ID (got {sorted(labeled_present)})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
