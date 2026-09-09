"""Offline unit tests for Genomics_module_retriever.resolve_kegg_ids and the
retriever changes that feed it (Entrez ID retention, harvest-for-ambiguous).

Written after a real bug was caught by hand during development: a gene
recovered ONLY via the NCBI Entrez bridge (no UniProt accession at all) got a
"KEGG ID" filled in but its Resolution_status stayed "ambiguous" - which
means it would NEVER reach final_df (API.py's table1_df = resolved_df[~
ambiguous_mask] excludes every ambiguous row unconditionally), so the KEGG ID
would be silently unreachable by any analysis. Fixed by promoting any
previously-ambiguous row that gets a KEGG ID (via ANY bridge) to a new status,
"resolved_via_kegg". These tests exist specifically to guard that promotion,
since it is easy to reintroduce silently (e.g. by adding a new bridge and
forgetting the promotion step).

No server needed - _safe_uniprot_search / _ncbi_gene_search / build_kegg_gene_map
are all monkeypatched (no network).
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import pandas as pd

import API
import Genomics_module as gm
import Genomics_module_retriever as gr

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


print("=== _resolve_via_geneid: NCBIGeneID kept even when the bridge stays ambiguous ===")
_real_search = gr._safe_uniprot_search
gr._safe_uniprot_search = lambda query, verbose=False: (
    [{"primaryAccession": "P1", "uniProtkbId": "P1_HUMAN", "proteinDescription": {},
      "organism": {"taxonId": 9606}, "extraAttributes": {}},
     {"primaryAccession": "P2", "uniProtkbId": "P2_HUMAN", "proteinDescription": {},
      "organism": {"taxonId": 9606}, "extraAttributes": {}}]
    if "reviewed:true" in query else []
)
try:
    result = gr._resolve_via_geneid("999", taxon_for_ref="", verbose=False)
    check(result["Resolution_status"] == "ambiguous",
          f"two hits stays ambiguous, not silently collapsed (got {result['Resolution_status']})")
    check(result["NCBIGeneID"] == "999",
          f"SENSITIVITY: the Entrez ID is kept even though the UniProt bridge itself stayed ambiguous "
          f"(got {result['NCBIGeneID']!r}) - this is the exact fix for the bug where it used to be thrown away")
finally:
    gr._safe_uniprot_search = _real_search

print("=== resolve_gene_symbol: harvest-only Entrez lookup for an ambiguous symbol ===")
gr._safe_uniprot_search = lambda query, verbose=False: (
    [{"primaryAccession": "P1", "uniProtkbId": "P1_X", "proteinDescription": {},
      "organism": {"taxonId": 1}, "extraAttributes": {}},
     {"primaryAccession": "P2", "uniProtkbId": "P2_X", "proteinDescription": {},
      "organism": {"taxonId": 1}, "extraAttributes": {}}]
    if "reviewed:false" in query else []
)
gr._ncbi_gene_search = lambda term, verbose=False: ["12345"]
gr._NCBI_GENEID_CACHE.clear()
try:
    result = gr.resolve_gene_symbol("AMBIGSYM", "1", verbose=False)
    check(result["Resolution_status"] == "ambiguous", "status is untouched by the harvest step")
    check(result["Candidates"], "Candidates is untouched (still populated) by the harvest step")
    check(result["NCBIGeneID"] == "12345", f"a single esearch hit is harvested into NCBIGeneID (got {result['NCBIGeneID']!r})")
finally:
    gr._safe_uniprot_search = _real_search

print("=== resolve_gene_symbol: two esearch hits -> NCBIGeneID stays empty (never guesses) ===")
gr._safe_uniprot_search = lambda query, verbose=False: (
    [{"primaryAccession": "P1", "uniProtkbId": "P1_X", "proteinDescription": {},
      "organism": {"taxonId": 1}, "extraAttributes": {}},
     {"primaryAccession": "P2", "uniProtkbId": "P2_X", "proteinDescription": {},
      "organism": {"taxonId": 1}, "extraAttributes": {}}]
    if "reviewed:false" in query else []
)
gr._ncbi_gene_search = lambda term, verbose=False: ["111", "222"]
gr._NCBI_GENEID_CACHE.clear()
try:
    result = gr.resolve_gene_symbol("AMBIGSYM2", "1", verbose=False)
    check(result["Resolution_status"] == "ambiguous", "status is untouched")
    check(result["NCBIGeneID"] == "",
          f"SENSITIVITY: two esearch hits means real ambiguity (likely paralogs) - the bridge stays "
          f"silent rather than guessing gene_ids[0] (got {result['NCBIGeneID']!r})")
finally:
    gr._safe_uniprot_search = _real_search
    gr._ncbi_gene_search = gr._ncbi_gene_search  # restore happens naturally via module reload boundary in other tests

print("=== resolve_kegg_ids: direct bridges (accession / entrez / symbol), first hit wins ===")
_real_build_gene_map = gm.build_kegg_gene_map
gm.build_kegg_gene_map = lambda org, verbose=False, force_refresh=False: {
    "acc_to_gene": {"ACCA": "org:1"},
    "entrez_to_gene": {"555": "org:2"},
    "symbol_to_gene": {"SYMB": "org:3"},
    "has_entrez": True, "has_symbols": True, "fetched_at": "", "caveat": None,
}
try:
    df = pd.DataFrame([
        {"Gene Symbol": "GENEA", "UniProt ID": "ACCA", "Entry Name": "", "Protein Name": "",
         "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "resolved_unreviewed", "Candidates": ""},
        {"Gene Symbol": "GENEB", "UniProt ID": "", "Entry Name": "", "Protein Name": "",
         "NCBI Gene ID": "555", "KEGG ID": "", "Resolution_status": "ambiguous", "Candidates": "[]"},
        {"Gene Symbol": "SYMB", "UniProt ID": "", "Entry Name": "", "Protein Name": "",
         "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "not_found", "Candidates": ""},
    ])
    mapped = gr.resolve_kegg_ids(df, "Gene Symbol", "org", verbose=False)
    row_a = mapped[mapped["Gene Symbol"] == "GENEA"].iloc[0]
    row_b = mapped[mapped["Gene Symbol"] == "GENEB"].iloc[0]
    row_c = mapped[mapped["Gene Symbol"] == "SYMB"].iloc[0]
    check(row_a["KEGG ID"] == "org:1", f"accession bridge fires for a resolved row (got {row_a['KEGG ID']!r})")
    check(row_a["Resolution_status"] == "resolved_unreviewed",
          "a row that was NOT ambiguous is never touched by the promotion step")
    check(row_b["KEGG ID"] == "org:2", f"entrez bridge fires for an ambiguous row with a harvested Entrez ID (got {row_b['KEGG ID']!r})")
    check(row_b["Resolution_status"] == "resolved_via_kegg",
          f"CRITICAL: the entrez-only match promotes the row out of 'ambiguous' so it can reach final_df "
          f"(got {row_b['Resolution_status']!r}) - this is the exact bug this test file exists to guard")
    check(row_c["KEGG ID"] == "org:3", f"symbol bridge fires as a last resort (got {row_c['KEGG ID']!r})")

    print("=== resolve_kegg_ids: accession-candidate disambiguation (no direct bridge fires) ===")
    df2 = pd.DataFrame([{
        "Gene Symbol": "GENED", "UniProt ID": "", "Entry Name": "", "Protein Name": "",
        "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "ambiguous",
        "Candidates": '[{"accession": "NOMATCH", "entry_name": "N1", "protein_name": "P1"}, '
                      '{"accession": "ACCA", "entry_name": "E_A", "protein_name": "Protein A"}]',
    }])
    mapped2 = gr.resolve_kegg_ids(df2, "Gene Symbol", "org", verbose=False)
    row_d = mapped2.iloc[0]
    check(row_d["KEGG ID"] == "org:1",
          f"exactly one candidate (ACCA) is in KEGG's accession map -> disambiguated (got {row_d['KEGG ID']!r})")
    check(row_d["UniProt ID"] == "ACCA",
          f"disambiguation ALSO fills the real UniProt ID, not just KEGG ID (got {row_d['UniProt ID']!r})")
    check(row_d["Resolution_status"] == "resolved_via_kegg", "promoted out of ambiguous")
    check(row_d["Candidates"] == "", "Candidates is cleared once a decision has been made")

    print("=== resolve_kegg_ids: SENSITIVITY - two matching candidates -> left ambiguous, never guessed ===")
    df3 = pd.DataFrame([{
        "Gene Symbol": "GENEE", "UniProt ID": "", "Entry Name": "", "Protein Name": "",
        "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "ambiguous",
        "Candidates": '[{"accession": "ACCA", "entry_name": "E1", "protein_name": "P1"}, '
                      '{"accession": "ACCB", "entry_name": "E2", "protein_name": "P2"}]',
    }])
    gm.build_kegg_gene_map = lambda org, verbose=False, force_refresh=False: {
        "acc_to_gene": {"ACCA": "org:1", "ACCB": "org:9"},  # BOTH candidates map to a (different) real gene
        "entrez_to_gene": {}, "symbol_to_gene": {},
        "has_entrez": False, "has_symbols": False, "fetched_at": "", "caveat": None,
    }
    mapped3 = gr.resolve_kegg_ids(df3, "Gene Symbol", "org", verbose=False)
    row_e = mapped3.iloc[0]
    check(row_e["KEGG ID"] == "", f"two DIFFERENT genes among the candidates -> no disambiguation (got {row_e['KEGG ID']!r})")
    check(row_e["Resolution_status"] == "ambiguous", "row correctly stays ambiguous rather than guessing")
finally:
    gm.build_kegg_gene_map = _real_build_gene_map

print("=== resolve_kegg_ids: empty organism code is a no-op (KEGG ID stays empty for every row) ===")
df4 = pd.DataFrame([{"Gene Symbol": "X", "UniProt ID": "ACC1", "Entry Name": "", "Protein Name": "",
                      "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "resolved", "Candidates": ""}])
mapped4 = gr.resolve_kegg_ids(df4, "Gene Symbol", "", verbose=False)
check((mapped4["KEGG ID"] == "").all(), "an empty KEGG organism code leaves every KEGG ID blank, no error")

print("=== resolve_kegg_ids: does not mutate the input DataFrame ===")
df5 = pd.DataFrame([{"Gene Symbol": "X", "UniProt ID": "ACC1", "Entry Name": "", "Protein Name": "",
                      "NCBI Gene ID": "", "KEGG ID": "", "Resolution_status": "resolved", "Candidates": ""}])
df5_copy_before = df5.copy()
gr.resolve_kegg_ids(df5, "Gene Symbol", "", verbose=False)
check(df5.equals(df5_copy_before), "the input DataFrame is left unchanged (a new frame is returned)")

print("=== _job_run_ora: KEGG organism guard - session value wins over a mismatched config taxon_id ===")
# Regression test for a real bug found by audit: the config screen used to
# carry its OWN independent KEGG organism picker (dbCfg.taxonId), separate
# from the one the retrieval step actually resolved KEGG IDs against
# (session["ora"]["kegg_organism_code"]). The two could silently diverge -
# genes resolved as "hsa" tested against an "mmu" universe, an empty result
# with no warning. Fixed by (a) removing the second picker in static/app.js
# and (b) this server-side guard, which must never trust the config screen's
# taxon_id for KEGG over the session's own record of what the gene IDs were
# actually resolved for.
_real_build_kegg_universe = gm.build_kegg_universe
captured_taxon_ids = []


def _fake_build_kegg_universe(taxon_id, verbose=False, force_refresh=False):
    captured_taxon_ids.append(taxon_id)
    return {"pathways": {"p1": {"name": "Test", "genes": {"hsa:1"}}},
            "universe_size": 1, "fetched_at": "", "caveat": None}


gm.build_kegg_universe = _fake_build_kegg_universe
try:
    final_df = pd.DataFrame([{
        "Gene Symbol": "GENEA", "UniProt ID": "ACC1", "Entry Name": "", "Protein Name": "",
        "NCBI Gene ID": "", "KEGG ID": "hsa:1", "Resolution_status": "resolved", "Candidates": "",
    }])
    session = {"ora": API._fresh_ora_state()}
    session["ora"]["final_df"] = final_df
    # The organism the gene IDs were ACTUALLY resolved for (set by
    # _job_resolve_kegg_ids in the real flow).
    session["ora"]["kegg_organism_code"] = "hsa"
    # A DELIBERATELY mismatched config-screen organism, simulating the
    # divergence the removed second picker used to allow.
    session["ora"]["ora_config"] = {
        "databases": ["KEGG"],
        "per_database": {"KEGG": {"taxon_id": "mmu", "taxon_name": "", "include_iea": True}},
    }

    result = API._job_run_ora("fake-job-kegg-guard", session)
    check(captured_taxon_ids == ["hsa"],
          f"build_kegg_universe was called with the SESSION's organism ('hsa'), never the mismatched "
          f"config-screen one ('mmu') (got {captured_taxon_ids})")
    kegg_block = result.get("KEGG", {}).get("KEGG", {})
    check(bool(kegg_block.get("caveat")) and "did not match" in kegg_block.get("caveat", ""),
          f"a mismatch is surfaced as a caveat, not silently corrected (got {kegg_block.get('caveat')!r})")

    print("=== _job_run_ora: KEGG organism guard - no caveat when the two agree ===")
    captured_taxon_ids.clear()
    session2 = {"ora": API._fresh_ora_state()}
    session2["ora"]["final_df"] = final_df
    session2["ora"]["kegg_organism_code"] = "hsa"
    session2["ora"]["ora_config"] = {
        "databases": ["KEGG"],
        "per_database": {"KEGG": {"taxon_id": "hsa", "taxon_name": "", "include_iea": True}},
    }
    result2 = API._job_run_ora("fake-job-kegg-agree", session2)
    check(captured_taxon_ids == ["hsa"], f"agreeing organisms still build for 'hsa' (got {captured_taxon_ids})")
    kegg_block2 = result2.get("KEGG", {}).get("KEGG", {})
    check(not kegg_block2.get("caveat"), f"no spurious caveat when config and session already agree (got {kegg_block2.get('caveat')!r})")

    print("=== _job_run_burden: the identical guard ===")
    captured_taxon_ids.clear()
    session3 = {"ora": API._fresh_ora_state()}
    session3["ora"]["final_df"] = final_df
    session3["ora"]["kegg_organism_code"] = "hsa"
    session3["ora"]["variant_scores"] = {"genea": 1.0}
    session3["ora"]["burden_config"] = {
        "databases": ["KEGG"],
        "per_database": {"KEGG": {"taxon_id": "mmu", "taxon_name": "", "include_iea": True}},
    }
    result3 = API._job_run_burden("fake-job-kegg-burden-guard", session3)
    check(captured_taxon_ids == ["hsa"],
          f"burden's build_kegg_universe call also used the session organism, not the mismatched config "
          f"one (got {captured_taxon_ids})")
    kegg_block3 = result3.get("KEGG", {}).get("KEGG", {})
    check(bool(kegg_block3.get("caveat")) and "did not match" in kegg_block3.get("caveat", ""),
          f"burden also surfaces the mismatch as a caveat (got {kegg_block3.get('caveat')!r})")
finally:
    gm.build_kegg_universe = _real_build_kegg_universe

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
