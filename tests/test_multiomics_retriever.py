# -*- coding: utf-8 -*-
"""Unit tests for Multiomics_module_retriever (the thin glue file), offline.

What is under test is that this file adds NO resolution logic of its own:
  - every upload / resolve call is routed to that omics' OWN existing retriever
    function with the arguments that omics' own tile would pass (checked with
    spies on the real modules), so an ID resolves exactly as it does in the
    single-omics tile;
  - the per-omics list rule: Transcriptomics / Proteomics significant sets are
    exactly what that module's own selection function returns for direction
    "all" (up AND down pooled), against a reference set of the genes /
    proteins TESTED; Genomics rows that are ambiguous or unresolved are left
    out, never guessed.
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module_retriever as gr  # noqa: E402
import Multiomics_module_retriever as mr  # noqa: E402
import Proteomics_module as pm  # noqa: E402
import Proteomics_module_retriever as pr  # noqa: E402
import Transcriptomics_module as tm  # noqa: E402
import Transcriptomics_module_retriever as tr  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


calls = []

# ---------------------------------------------------------------------------
print("=== parse_upload routes to each omics' own loader ===")
_orig = (gr.load_gene_table_from_bytes, gr.deduplicate_gene_table, tr.load_de_table_from_bytes, tr.detect_columns,
         pr.load_de_table_from_bytes, pr.detect_columns)
gr.load_gene_table_from_bytes = lambda b: (calls.append(("g.load", b)) or (pd.DataFrame({"Gene Symbol": ["A", "A", "B"]}), "gene_symbol"))
gr.deduplicate_gene_table = lambda df, t: (calls.append(("g.dedup", t)) or (df.drop_duplicates(), {"duplicates_removed": 1}))
tr.load_de_table_from_bytes = lambda b: (calls.append(("t.load", b)) or pd.DataFrame({"gene": ["A"], "lfc": [1.0], "p": [0.1]}))
tr.detect_columns = lambda df: (calls.append(("t.detect",)) or {"identifier": "gene", "log2FC": "lfc", "pvalue": "p"})
pr.load_de_table_from_bytes = lambda b: (calls.append(("p.load", b)) or pd.DataFrame({"prot": ["P1"]}))
pr.detect_columns = lambda df: (calls.append(("p.detect",)) or {"identifier": "prot"})

g = mr.parse_upload("genomics", b"gg")
check(("g.load", b"gg") in calls and ("g.dedup", "gene_symbol") in calls and g["identifier_type"] == "gene_symbol"
      and g["kind"] == "list" and len(g["df"]) == 2 and g["n_duplicates_removed"] == 1,
      "Genomics: Genomics' own loader + de-duplicator, the detected identifier type kept, duplicates counted")
t = mr.parse_upload("transcriptomics", b"tt")
check(("t.load", b"tt") in calls and ("t.detect",) in calls and t["kind"] == "de" and t["column_mapping"]["pvalue"] == "p",
      "Transcriptomics: its own DE-table loader + column detection (a table with log2FC and p-value is REQUIRED)")
p = mr.parse_upload("proteomics", b"pp")
check(("p.load", b"pp") in calls and ("p.detect",) in calls and p["kind"] == "de",
      "Proteomics: its own loader + column detection")
try:
    mr.parse_upload("metabolomics", b"x")
    check(False, "metabolomics is not parsed here")
except ValueError:
    check(True, "Metabolomics is refused here - it is parsed by API's existing shared loader, not duplicated")
try:
    mr.parse_upload("bogus", b"x")
    check(False, "an unknown omics is refused")
except ValueError:
    check(True, "an unknown omics type is refused")
(gr.load_gene_table_from_bytes, gr.deduplicate_gene_table, tr.load_de_table_from_bytes, tr.detect_columns,
 pr.load_de_table_from_bytes, pr.detect_columns) = _orig

# ---------------------------------------------------------------------------
print("=== resolve_upload routes with the arguments that omics' own tile passes ===")
calls.clear()
_orig2 = (gr.resolve_gene_table, tr.resolve_de_table, pr.resolve_protein_table)
gr.resolve_gene_table = lambda df, id_type, taxon_id=None, verbose=False: (
    calls.append(("g.resolve", id_type, taxon_id)) or pd.DataFrame({
        "Gene Symbol": ["A", "B", "C", "D"], "UniProt ID": ["U1", "U2", "", ""],
        "Resolution_status": ["resolved", "ambiguous", "not_found", "resolved_via_ncbi"]}))
tr.resolve_de_table = lambda df, taxon, identifier_type=None, progress_callback=None, verbose=False: (
    calls.append(("t.resolve", taxon, identifier_type)) or (pd.DataFrame({"gene_symbol": ["A"]}), {"n_mapped": 1}))
pr.resolve_protein_table = lambda df, taxon, progress_callback=None, verbose=False: (
    calls.append(("p.resolve", taxon)) or (pd.DataFrame({"accession": ["P1"]}), {"n_mapped": 1}))

out = mr.resolve_upload("genomics", pd.DataFrame(), "9606", identifier_type="gene_symbol")
check(("g.resolve", "gene_symbol", "9606") in calls, "Genomics gene symbols: resolved with the organism, like /api/ora/genes/resolve")
mr.resolve_upload("genomics", pd.DataFrame(), "9606", identifier_type="ensembl_gene_id")
check(("g.resolve", "ensembl_gene_id", None) in calls,
      "Genomics Ensembl IDs: NO taxon passed (Genomics' own rule - Ensembl IDs are organism-specific already)")
check(out["meta"]["n_ambiguous"] == 1 and out["meta"]["n_not_found"] == 1 and out["meta"]["n_resolved"] == 2,
      f"the summary counts ambiguous / unresolved genes that will be left out (got {out['meta']})")
mr.resolve_upload("transcriptomics", pd.DataFrame(), "9606")
check(("t.resolve", "9606", None) in calls, "Transcriptomics: identifier type auto-detected (None), exactly as its own resolver allows")
mr.resolve_upload("proteomics", pd.DataFrame(), "9606")
check(("p.resolve", "9606") in calls, "Proteomics: its own UniProt ID-mapping resolver")
try:
    mr.resolve_upload("metabolomics", pd.DataFrame(), "9606")
    check(False, "metabolomics is not resolved here")
except ValueError:
    check(True, "Metabolomics is resolved through API's existing _resolve_and_match_compound_names, not here")
gr.resolve_gene_table, tr.resolve_de_table, pr.resolve_protein_table = _orig2

print("=== Genomics: ambiguous and unresolved rows are never guessed ===")
resolved = pd.DataFrame({"Gene Symbol": list("ABCDE"),
                         "Resolution_status": ["resolved", "ambiguous", "not_found", "user_provided", "resolved_via_kegg"]})
usable = mr.usable_genomics_rows(resolved)
check(list(usable["Gene Symbol"]) == ["A", "D", "E"] and len(resolved) == 5,
      "ambiguous and not_found rows are dropped, every other status kept (and the input is not mutated)")
_orig3 = gr.resolve_kegg_ids
gr.resolve_kegg_ids = lambda df, id_col, org, verbose=False: (calls.append(("g.kegg", id_col, org)) or df)
mr.assign_genomics_kegg_ids(pd.DataFrame({"Gene Symbol": ["A"]}), "hsa")
mr.assign_genomics_kegg_ids(pd.DataFrame({"Ensembl Gene ID": ["E"]}), "mmu")
check(("g.kegg", "Gene Symbol", "hsa") in calls and ("g.kegg", "Ensembl Gene ID", "mmu") in calls,
      "the KEGG-gene-ID pass is Genomics' own resolve_kegg_ids, run with the single organism code it is given")
gr.resolve_kegg_ids = _orig3

# ---------------------------------------------------------------------------
print("=== DE selection: pooled up + down, reference = genes/proteins TESTED ===")
de = pd.DataFrame({
    "gene_symbol": ["UPG", "DOWNG", "FLAT", "LOWP", "BORDER", "", "NOPADJ"],
    "entrez_id": ["1", "2", "3", "4", "5", "6", "7"],
    "log2FC": [2.5, -2.5, 0.2, 3.0, 1.0, 5.0, 4.0],
    "padj": [0.001, 0.001, 0.001, 0.4, 0.05, 0.001, float("nan")],
})
sel = mr.select_transcriptomics(de, 0.05, 1.0)
check(sel["significant"] == {"UPG", "DOWNG", "BORDER"},
      f"up AND down both selected, boundary values inclusive, weak/flat/NaN-padj/unmapped excluded (got {sorted(sel['significant'])})")
own = set(tm.select_degs(de, 0.05, 1.0, "all")["gene_symbol"].astype(str).str.upper())
check(sel["significant"] == own, "SENSITIVITY: identical to Transcriptomics' own select_degs(direction='all')")
up_only = set(tm.select_degs(de, 0.05, 1.0, "up")["gene_symbol"].astype(str).str.upper())
check("DOWNG" not in up_only and "DOWNG" in sel["significant"],
      "...and a direction-split selection would have lost the down-regulated gene (why pooling matters)")
check(sel["n_tested"] == 6 and sel["n_significant"] == 3 and len(sel["universe_df"]) == 6,
      f"the reference set is every TESTED gene with a symbol (6), not a database (got n_tested={sel['n_tested']})")
sel2 = mr.select_transcriptomics(de, 0.0005, 2.0)
check(sel2["significant"] == set(), "tighter cutoffs narrow it (and an empty selection is returned, not an error)")

dp = pd.DataFrame({
    "accession": ["P1", "P2", "P3", "P4", ""],
    "gene_symbol": ["a", "b", "c", "d", "e"],
    "log2FC": [2.0, -2.0, 0.1, 3.0, 5.0], "padj": [0.01, 0.01, 0.01, 0.9, 0.01],
})
sp = mr.select_proteomics(dp, 0.05, 1.0)
check(sp["significant"] == {"P1", "P2"} and sp["n_tested"] == 4,
      f"Proteomics: pooled up + down, reference = measured accessions (got {sorted(sp['significant'])}, tested {sp['n_tested']})")
check(sp["significant"] == set(pm.select_significant_proteins(dp, 0.05, 1.0, "all")["accession"].astype(str)),
      "SENSITIVITY: identical to Proteomics' own select_significant_proteins(direction='all')")

print("=== matrix feature resolution keeps the ids the pathway step needs (MOFA enrichment) ===")
_orig_map, _orig_idm = tr.map_identifiers, pr._idmapping_resolve


def _fake_map(df, id_type=None, taxon_id=None, progress_callback=None, verbose=False):
    out = df.copy()
    out["gene_symbol"] = ["tp53", "brca1", ""]
    out["entrez_id"] = ["7157", "", "nan"]
    out["status"] = ["mapped", "mapped", "unmapped"]
    return out, {}


tr.map_identifiers = _fake_map
rt = mr.resolve_matrix_features("transcriptomics", ["A", "B", "C"], "9606")
tr.map_identifiers = _orig_map
check(rt["gene_of"] == {"A": "TP53", "B": "BRCA1"} and rt["entrez_of"] == {"A": "7157"},
      f"transcripts: the Entrez id comes with the gene symbol; a missing / 'nan' id is simply absent ({rt['entrez_of']})")
pr._idmapping_resolve = lambda accs, verbose=False: ({"P04637": {"gene_symbol": "tp53"}}, 0)
rp = mr.resolve_matrix_features("proteomics", ["P04637-2", "sp|P04637|X;Q9ZZZ9", "junk"], "9606")
pr._idmapping_resolve = _orig_idm
check(rp["acc_of"].get("P04637-2") == "P04637" and rp["gene_of"].get("P04637-2") == "TP53" and "junk" in rp["acc_of"],
      f"proteins: the cleaned accession (isoform / protein-group stripped) is kept next to the symbol ({rp['acc_of']})")

print("=== constants ===")
check(mr.OMICS_TYPES == ("genomics", "transcriptomics", "proteomics", "metabolomics") and mr.DE_OMICS == ("transcriptomics", "proteomics"),
      "the four omics and the two DE-table omics")
check(mr.PADJ_CUTOFF_DEFAULT == tm.PADJ_CUTOFF_DEFAULT == 0.05 and mr.LFC_CUTOFF_DEFAULT == tm.LFC_CUTOFF_DEFAULT == 1.0,
      "default cutoffs equal the Transcriptomics module's own defaults")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
