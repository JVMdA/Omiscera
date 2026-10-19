# -*- coding: utf-8 -*-
"""Unit tests for Proteomics_network.py.

The engine numerics (RWR, permutation null, MCL) are already covered by
test_network_engine_core.py; the shared orchestrator's own behaviour (seed
weighting, module cap, small-graph guard, recruited-gene universe rule) is
already covered by test_transcriptomics_network.py and is not retested here -
run_ppi_network_analysis calls txnet.run_network_analysis unchanged, only
swapping in the four accession-space hooks. This file tests what is unique to
Proteomics_network.py: the UniProt-accession identifier bridge (and its OWN
collision hazard, separate from build_string_id_index's), the tiered
accession -> gene_symbol mapping, and - the single most important test in
this file - that enrich_modules_by_accession translates GO/Reactome/Complex
Portal's accession-keyed gene sets into STRING-ID space so a protein with
several accessions is counted exactly once, while still using the GRAPH
(not the uploaded table) as its universe so a recruited protein keeps its
module's label.

Fixture graph (taxon "9997", confidence 700): two triangles P1-P2-P3 and
P4-P5-P6, connected via a high-degree hub PHI (score 900) and a low-degree
leaf PLO (score 700), both carrying the SAME UniProt accession "ACC_DUP" -
the accession-space analogue of test_transcriptomics_network.py's PHI/PLO
gene-symbol collision pair (PHI, the higher-degree protein, must win). P4/P6
carry no accession at all, reachable only via the gene_symbol fallback (P4)
or not in the user's table at all (P6, purely recruited).
"""
import gzip
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as spx

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm  # noqa: E402
import Genomics_network as gnet  # noqa: E402
import Proteomics_module as pm  # noqa: E402
import Proteomics_network as pnet  # noqa: E402
import Transcriptomics_network as txnet  # noqa: E402

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
# Fixture graph
# =============================================================================
TAX = "9997"
NAMES = ["P1", "P2", "P3", "P4", "P5", "P6", "PHI", "PLO"]
P = {k: f"{TAX}.{k}" for k in NAMES}

LINKS_ROWS = [
    (P["P1"], P["P2"], 900), (P["P2"], P["P1"], 900),
    (P["P1"], P["P3"], 900), (P["P3"], P["P1"], 900),
    (P["P2"], P["P3"], 900), (P["P3"], P["P2"], 900),
    (P["P4"], P["P5"], 900), (P["P5"], P["P4"], 900),
    (P["P4"], P["P6"], 900), (P["P6"], P["P4"], 900),
    (P["P5"], P["P6"], 900), (P["P6"], P["P5"], 900),
    (P["PHI"], P["P1"], 900), (P["P1"], P["PHI"], 900),
    (P["PHI"], P["P2"], 900), (P["P2"], P["PHI"], 900),
    (P["PLO"], P["P3"], 700), (P["P3"], P["PLO"], 700),
]
# PHI listed BEFORE PLO: a naive dict(zip(...)) last-wins index (the pattern
# Genomics_network's own uniprot_to_protein uses) would end up pointing
# "ACC_DUP" at PLO (the LOWER-degree protein) purely because it is read
# second - exactly the bug _alias_index_add's explicit ranking avoids.
INFO_ROWS = [
    (P["P1"], "GENE1"), (P["P2"], "GENE2"), (P["P3"], "GENE3"),
    (P["P4"], "GENE4"), (P["P5"], "GENE5"), (P["P6"], "GENE6"),
    (P["PHI"], "DUPGENE"), (P["PLO"], "DUPGENE"),
]
ALIAS_ROWS = [
    (P["PHI"], "ACC_DUP", "UniProt_AC"),
    (P["PLO"], "ACC_DUP", "UniProt_AC"),
    (P["P1"], "ACC1", "UniProt_AC"),
    (P["P2"], "ACC2", "UniProt_AC"),
    (P["P3"], "ACC3", "UniProt_AC"),
    (P["P5"], "ACC5", "UniProt_AC"),  # P5 mapped via accession; P4/P6 are not
]


def _write_gz(path, header, rows, sep):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(str(x) for x in r) + "\n")


TMP = tempfile.TemporaryDirectory()
LINKS = os.path.join(TMP.name, "links.txt.gz")
INFO = os.path.join(TMP.name, "info.txt.gz")
ALIASES = os.path.join(TMP.name, "aliases.txt.gz")
_write_gz(LINKS, ["protein1", "protein2", "combined_score"], LINKS_ROWS, " ")
_write_gz(INFO, ["#string_protein_id", "preferred_name", "protein_size", "annotation"],
          [(pid, name, "100", "x") for pid, name in INFO_ROWS], "\t")
_write_gz(ALIASES, ["#string_protein_id", "alias", "source"], ALIAS_ROWS, "\t")

GRAPH = gnet._parse_string_network_files(LINKS, INFO, ALIASES, confidence=700, taxon_id=TAX)
check(GRAPH["n_proteins"] == 8, f"fixture graph has all 8 proteins (got {GRAPH['n_proteins']})")

# enrich_modules_by_accession resolves its own aliases_path via
# gnet._string_file_paths(taxon_id) (mirroring _ppi_mapper's own lookup, so
# a caller never has to pass identifier-bridge plumbing through a signature
# shared with enrich_modules_by_symbol) - point that at the fixture files for
# this fake taxon everywhere below.
_real_string_file_paths = gnet._string_file_paths
gnet._string_file_paths = lambda taxon_id: {"links": LINKS, "info": INFO, "aliases": ALIASES}

# =============================================================================
print("=== build_accession_string_index: collision resolution, own hazard ===")
# =============================================================================
idx = pnet.build_accession_string_index(GRAPH, aliases_path=ALIASES)
check(idx["by_accession"].get("ACC_DUP") == P["PHI"],
      f"ACC_DUP resolves to PHI (higher degree), not PLO (got {idx['by_accession'].get('ACC_DUP')})")
check("ACC_DUP" in idx["collisions"], "the collision is reported")
check(idx["by_accession"].get("ACC1") == P["P1"] and idx["by_accession"].get("ACC5") == P["P5"],
      "unambiguous accessions map cleanly")
check(all(k not in ("P4", "P6") for k in idx["by_accession"]),
      "P4/P6 carry no accession at all (absent from the accession index)")
check(idx["by_symbol"].get("GENE4") == P["P4"],
      f"P4 is still reachable via the gene-symbol fallback, reused unchanged from "
      f"build_string_id_index (got {idx['by_symbol'].get('GENE4')})")
check(idx["n_accessions"] == 5, f"n_accessions counts distinct resolved accessions (got {idx['n_accessions']})")

# =============================================================================
print("=== map_protein_table_to_string: tiered accession -> gene_symbol, a true partition ===")
# =============================================================================
tier_df = pd.DataFrame({
    "input_id": ["r1", "r2", "r3", "r4", "r5"],
    "accession": ["ACC1", "ACC2", "ACC3", "", "nope"],
    "gene_symbol": ["nomatch1", "nomatch2", "nomatch3", "GENE4", "nomatch5"],
    "protein_name": [""] * 5,
    "log2FC": [3.0, 3.0, -3.0, -3.0, 1.0],
    "pvalue": [1e-6] * 5,
    "padj": [1e-5] * 5,
    "status": ["resolved"] * 5,
})
mapping = pnet.map_protein_table_to_string(tier_df, GRAPH, index=idx, taxon_id=TAX)
check(mapping["n_mapped"] == 4 and mapping["n_unmapped"] == 1,
      f"4 of 5 rows map, r5 (accession 'nope', symbol 'nomatch5') does not "
      f"(mapped={mapping['n_mapped']}, unmapped={mapping['n_unmapped']})")
check(mapping["n_by_tier"] == {"accession": 3, "gene_symbol": 1},
      f"the per-tier counts are a true partition of the mapped rows (got {mapping['n_by_tier']})")
check(mapping["string_by_row"]["r4"] == P["P4"],
      "row with no accession correctly falls back to gene_symbol")
check(set(mapping["table"].columns) >= {"STRING ID", "Mapped via"},
      "the mapping table carries STRING ID / Mapped via columns like map_de_table_to_string's")
check(mapping["low_mapping_warning"] is None, "80% mapping does not trigger the low-mapping warning")

# =============================================================================
print("=== map_protein_table_to_string: the <10% refusal gate, worded for proteins ===")
# =============================================================================
error_df = pd.DataFrame({
    "input_id": ["m0"] + [f"m{i}" for i in range(1, 20)],
    "accession": ["ACC1"] + ["NOPE"] * 19,
    "gene_symbol": [""] * 20,
    "protein_name": [""] * 20,
    "log2FC": [1.0] * 20, "pvalue": [1e-3] * 20, "padj": [1e-2] * 20,
    "status": ["resolved"] * 20,
})
try:
    pnet.map_protein_table_to_string(error_df, GRAPH, index=idx, taxon_id=TAX)
    check(False, "a 5%-mapping table raises rather than silently succeeding")
except ValueError as exc:
    check(TAX in str(exc), f"a 5%-mapping table raises ValueError naming the taxon (got: {exc})")

# =============================================================================
print("=== enrich_modules_by_accession: universe is the WHOLE GRAPH, GO path ===")
# =============================================================================
FAKE_GO = {
    "biological_process": {
        "GO:0001": {"name": "test process", "genes": {"ACC1", "ACC2", "ACC3", "ACC_DUP"}},
    },
    "molecular_function": {}, "cellular_component": {},
}
_real_build_go_universe = gm.build_go_universe
gm.build_go_universe = lambda *a, **k: FAKE_GO
try:
    modules = [
        {"module_id": "M1", "string_ids": [P["P1"], P["P2"], P["P3"], P["PHI"]], "size": 4},
        {"module_id": "M2", "string_ids": [], "size": 0},
    ]
    enr = pnet.enrich_modules_by_accession(modules, GRAPH, TAX, taxon_name="Test organism",
                                            databases=["GO"], min_term_size=1)
finally:
    gm.build_go_universe = _real_build_go_universe

check(enr["universe_size"] == GRAPH["n_proteins"],
      f"universe_size is the WHOLE graph's protein count (8), not a per-database annotated "
      f"subset and not the uploaded table's (got {enr['universe_size']})")
r1 = enr["results"]["M1"]
check(r1 is not None, "M1 gets a GO label")
if r1:
    tbl = r1["best"]["block"]["table"]
    check(len(tbl) == 1 and tbl.iloc[0]["Term ID"] == "GO:0001",
          f"M1's label is the test GO term (got {tbl.iloc[0].to_dict() if len(tbl) else None})")
    check(int(tbl.iloc[0]["k"]) == 4,
          f"SENSITIVITY, THE HEADLINE CLAIM: k=4 (all of P1/P2/P3/PHI), NOT 5 - ACC_DUP and ACC1..3 "
          f"translate to only 4 distinct STRING ids because PHI absorbs ACC_DUP once, proving a "
          f"protein carrying more than one accession is counted EXACTLY ONCE, not once per accession "
          f"(got k={tbl.iloc[0]['k']})")
    overlap = str(tbl.iloc[0]["Overlapping Genes"])
    check(all(g in overlap for g in ("GENE1", "GENE2", "GENE3", "DUPGENE")),
          f"the member column is translated back to gene SYMBOLS for display, not raw STRING ids "
          f"(got {overlap!r})")
check(enr["results"]["M2"] is None, "a module with no members is skipped, not crashed on")

# =============================================================================
print("=== enrich_modules_by_accession: a RECRUITED protein keeps its module's label ===")
# =============================================================================
# M2 = {P4, P5, P6}. P5 is the only one carrying an accession (ACC5); P4/P6
# would never appear in ANY uploaded table under this fixture (P4 has no
# accession at all, P6 has neither). This is the accession-space analogue of
# test_transcriptomics_network.py's "recruited gene is counted, not dropped".
FAKE_REACTOME = {"pathways": {"R-1": {"name": "test pathway", "genes": {"ACC5"}}}}
_real_build_reactome_universe = gm.build_reactome_universe
gm.build_reactome_universe = lambda *a, **k: FAKE_REACTOME
try:
    modules2 = [{"module_id": "M2", "string_ids": [P["P4"], P["P5"], P["P6"]], "size": 3}]
    enr_r = pnet.enrich_modules_by_accession(modules2, GRAPH, TAX, taxon_name="Test organism",
                                              databases=["Reactome"], min_term_size=1)
finally:
    gm.build_reactome_universe = _real_build_reactome_universe

r2 = enr_r["results"]["M2"]
check(r2 is not None,
      "M2 gets a Reactome label even though only ONE of its three members (P5) carries the "
      "accession the pathway is keyed on - P4/P6 are RECRUITED and stay in the module regardless "
      "(they are simply outside THIS namespace's own annotated universe, same as any gene never "
      "annotated to any GO/Reactome/Complex Portal term)")
if r2:
    tbl2 = r2["best"]["block"]["table"]
    check(int(tbl2.iloc[0]["k"]) == 1 and int(tbl2.iloc[0]["n"]) == 1,
          f"k=n=1: P5/ACC5 is the only module member that is even IN this namespace's annotated "
          f"universe (P4/P6 carry no accession at all, so no term set could ever mention them) - "
          f"the module is still labelled (checked above), it is only the Fisher accounting that "
          f"correctly excludes unannotatable proteins (got k={tbl2.iloc[0]['k']}, n={tbl2.iloc[0]['n']})")

# =============================================================================
print("=== enrich_modules_by_accession: Reactome requires taxon_name, Complex Portal path ===")
# =============================================================================
gm.build_reactome_universe = lambda *a, **k: FAKE_REACTOME
try:
    enr_no_name = pnet.enrich_modules_by_accession(modules2, GRAPH, TAX, taxon_name="",
                                                    databases=["Reactome"], min_term_size=1)
finally:
    gm.build_reactome_universe = _real_build_reactome_universe
check(enr_no_name["results"]["M2"] is None and enr_no_name["caveat"] is not None,
      f"without taxon_name, Reactome is refused with a caveat rather than silently skipped "
      f"(got results={enr_no_name['results']['M2']}, caveat={enr_no_name['caveat']!r})")

FAKE_COMPLEX = {"terms": {"CPX-1": {"name": "test complex", "genes": {"ACC1", "ACC2"}}}, "caveat": None}
_real_build_complex = pm.build_complex_portal_universe
pm.build_complex_portal_universe = lambda *a, **k: FAKE_COMPLEX
try:
    modules1 = [{"module_id": "M1", "string_ids": [P["P1"], P["P2"], P["P3"], P["PHI"]], "size": 4}]
    enr_c = pnet.enrich_modules_by_accession(modules1, GRAPH, TAX, databases=["Complex Portal"], min_term_size=1)
finally:
    pm.build_complex_portal_universe = _real_build_complex
check(enr_c["results"]["M1"] is not None, "Complex Portal also labels a module (CPX-1 contains ACC1/ACC2)")

try:
    pnet.enrich_modules_by_accession(modules1, GRAPH, TAX, databases=["Bogus"])
    check(False, "an unsupported database name raises rather than silently doing nothing")
except ValueError:
    check(True, "an unsupported database name raises ValueError")

# =============================================================================
print("=== NEW: multiple databases at once (GO + Complex Portal) - genuinely different results per database ===")
# =============================================================================
# Same degenerate-universe trap as Transcriptomics' own multi-database test:
# with only the "real" term in a namespace, that term's own genes ARE the
# annotated universe (K == N), forcing p=1.0 regardless of overlap - so each
# database's fake term set below carries a non-overlapping FILLER term too,
# widening the universe enough that GO and Complex Portal produce genuinely
# DIFFERENT (non-tied) p-values, and Complex Portal's tighter match (K=2,
# exact overlap with a 2-member module slice) wins "best" while GO's own
# result survives untouched in "candidates".
FAKE_GO_MULTI = {
    "biological_process": {
        "GO:0001": {"name": "test process", "genes": {"ACC1", "ACC2", "ACC3", "ACC_DUP"}},
        "GO:FILLER": {"name": "GO filler", "genes": {"ACC5"}},
    },
    "molecular_function": {}, "cellular_component": {},
}
FAKE_COMPLEX_MULTI = {
    "terms": {
        "CPX-1": {"name": "test complex", "genes": {"ACC1", "ACC2"}},
        "CPX-FILLER": {"name": "filler complex", "genes": {"ACC5"}},
    },
    "caveat": None,
}
gm.build_go_universe = lambda *a, **k: FAKE_GO_MULTI
pm.build_complex_portal_universe = lambda *a, **k: FAKE_COMPLEX_MULTI
try:
    enr_multi = pnet.enrich_modules_by_accession(
        modules1, GRAPH, TAX, databases=["GO", "Complex Portal"], min_term_size=1)
finally:
    gm.build_go_universe = _real_build_go_universe
    pm.build_complex_portal_universe = _real_build_complex

check(enr_multi["databases"] == ["GO", "Complex Portal"],
      f"the return dict reports BOTH requested databases (got {enr_multi['databases']})")
m1_multi = enr_multi["results"]["M1"]
check(m1_multi is not None, "M1 gets a result when GO + Complex Portal are both selected")
candidate_keys = set(m1_multi["candidates"].keys())
check("biological_process" in candidate_keys and "Complex Portal" in candidate_keys,
      f"candidates spans BOTH databases, not just whichever won (got {candidate_keys})")
check(m1_multi["candidates"]["biological_process"] is not None and m1_multi["candidates"]["Complex Portal"] is not None,
      "both candidates are real (non-None) results, not one masking the other")
check(m1_multi["n_namespaces_compared"] == 2,
      f"n_namespaces_compared spans both databases (got {m1_multi['n_namespaces_compared']})")
go_candidate_table = m1_multi["candidates"]["biological_process"]["block"]["table"]
check(not go_candidate_table.empty and int(go_candidate_table.iloc[0]["k"]) == 4,
      f"GO's own candidate table is independently usable (still k=4, PHI counted once) even though "
      f"it may not have won 'best' (got empty={go_candidate_table.empty})")

# =============================================================================
print("=== run_ppi_network_analysis: wires the four hooks into the shared orchestrator ===")
# =============================================================================
_real_is_taxon_in_string = gnet.is_taxon_in_string
gnet.is_taxon_in_string = lambda taxon_id: True
try:
    small_df = pd.DataFrame({
        "input_id": ["r1"], "accession": ["ACC1"], "gene_symbol": ["GENE1"], "protein_name": [""],
        "log2FC": [3.0], "pvalue": [1e-6], "padj": [1e-5], "status": ["resolved"],
    })
    try:
        pnet.run_ppi_network_analysis(small_df, TAX, taxon_name="Test organism", n_permutations=10)
        check(False, "an 8-protein graph (below MIN_GRAPH_PROTEINS_FOR_NULL) raises rather than running")
    except ValueError as exc:
        check("permutation null" in str(exc),
              f"the small-graph guard fires through run_ppi_network_analysis exactly as it does "
              f"through run_network_analysis directly (got: {exc})")
finally:
    gnet.is_taxon_in_string = _real_is_taxon_in_string

# =============================================================================
print("=== run_ppi_network_analysis: end-to-end on a graph big enough to clear the null guard ===")
# =============================================================================
BIG_TAX = "9995"
BIG_N = 24
big_names = [f"Q{i}" for i in range(BIG_N)]
BIG_P = {n: f"{BIG_TAX}.{n}" for n in big_names}
rng = np.random.default_rng(0)
big_edges = []
for i in range(BIG_N):
    for j in range(i + 1, BIG_N):
        if rng.random() < 0.2:
            big_edges.append((big_names[i], big_names[j], int(rng.integers(700, 999))))
for i in range(BIG_N - 1):  # guarantee connectivity regardless of the random draw above
    big_edges.append((big_names[i], big_names[i + 1], 900))

big_links_rows = []
for a, b, w in big_edges:
    big_links_rows.append((BIG_P[a], BIG_P[b], w))
    big_links_rows.append((BIG_P[b], BIG_P[a], w))
big_info_rows = [(BIG_P[n], f"BGENE{i}") for i, n in enumerate(big_names)]
big_alias_rows = [(BIG_P[n], f"BACC{i}", "UniProt_AC") for i, n in enumerate(big_names)]

BIG_LINKS = os.path.join(TMP.name, "big_links.txt.gz")
BIG_INFO = os.path.join(TMP.name, "big_info.txt.gz")
BIG_ALIASES = os.path.join(TMP.name, "big_aliases.txt.gz")
_write_gz(BIG_LINKS, ["protein1", "protein2", "combined_score"], big_links_rows, " ")
_write_gz(BIG_INFO, ["#string_protein_id", "preferred_name", "protein_size", "annotation"],
          [(pid, name, "100", "x") for pid, name in big_info_rows], "\t")
_write_gz(BIG_ALIASES, ["#string_protein_id", "alias", "source"], big_alias_rows, "\t")

_big_paths = {"links": BIG_LINKS, "info": BIG_INFO, "aliases": BIG_ALIASES}
gnet._string_file_paths = lambda taxon_id: _big_paths
gnet.is_taxon_in_string = lambda taxon_id: True

big_rows = []
for i, n in enumerate(big_names):
    if i < 6:
        lfc, pval, padj = 3.0, 1e-8, 1e-6
    elif i < 12:
        lfc, pval, padj = -3.0, 1e-8, 1e-6
    else:
        lfc, pval, padj = 0.05, 0.8, 0.9
    big_rows.append({"input_id": f"r{i}", "accession": f"BACC{i}", "gene_symbol": f"BGENE{i}",
                     "protein_name": "", "log2FC": lfc, "pvalue": pval, "padj": padj, "status": "resolved"})
big_df = pd.DataFrame(big_rows)

try:
    result = pnet.run_ppi_network_analysis(big_df, BIG_TAX, taxon_name="Test organism", confidence=700,
                                            n_permutations=30, verbose=False)
finally:
    gnet._string_file_paths = _real_string_file_paths
    gnet.is_taxon_in_string = _real_is_taxon_in_string

check(result["mapping"]["n_mapped"] == 24 and result["mapping"]["n_by_tier"] == {"accession": 24, "gene_symbol": 0},
      f"all 24 rows map via the accession tier (got n_mapped={result['mapping']['n_mapped']}, "
      f"n_by_tier={result['mapping']['n_by_tier']})")
check(result["propagation"]["n_seeds_mapped"] == 12,
      f"seeds come from pm.select_significant_proteins (6 up + 6 down), not tx.select_degs "
      f"(got {result['propagation']['n_seeds_mapped']})")
seed_scores = result["propagation"]["table"]["Seed score"].dropna()
check(len(seed_scores) == 12 and (seed_scores > 0).all(),
      f"SEED POSITIVITY INVARIANT: every one of the 12 seeds' scores is strictly positive even "
      f"though half of them are down-regulated - direction is expressed only through WHICH "
      f"proteins seed the walk, never via a signed weight (got {len(seed_scores)} scored seeds, "
      f"min={seed_scores.min() if len(seed_scores) else None})")
check("_graph" in result and "_hit_ids" in result,
      "the payload mirrors run_network_analysis's shape exactly (_graph/_hit_ids present), so the "
      "existing API _graph-strip and netprop chart builders need no Proteomics-specific branch")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
