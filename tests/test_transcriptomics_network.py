# -*- coding: utf-8 -*-
"""Unit tests for Transcriptomics_network.py, which had ZERO unit coverage.

The engine numerics (RWR, permutation null, MCL) are already covered by
test_network_engine_core.py and are NOT retested here. This file tests the
code that is unique to this module: the identifier bridge, seed weighting,
the module-detection scale guard, and - the single most important test in
this file - that per-module enrichment uses the GRAPH's genes as its universe,
not the DE table's, so a RECRUITED gene is not silently dropped from its own
module's label.

Fixture graph (taxon "9998", confidence 700): two triangles P1-P2-P3 and
P4-P5-P6; a collision pair PHI/PLO both named "DUPGENE" (PHI the higher-degree,
correct winner); P9/P10 named so their SYMBOL never matches the DE table,
forcing the Ensembl/Entrez fallback tiers to be exercised for real.
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

import Genomics_network as gnet  # noqa: E402
import Transcriptomics_module as tx  # noqa: E402
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
TAX = "9998"
NAMES = ["P1", "P2", "P3", "P4", "P5", "P6", "PHI", "PLO", "P9", "P10"]
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
    (P["P9"], P["P4"], 900), (P["P4"], P["P9"], 900),
    (P["P10"], P["P5"], 900), (P["P5"], P["P10"], 900),
]
# PHI listed BEFORE PLO: a naive dict(zip(...)) last-wins index would end up
# pointing "DUPGENE" at PLO (the LOWER-degree protein) purely because it is
# read second - exactly the bug build_string_id_index's explicit ranking
# avoids.
INFO_ROWS = [
    (P["P1"], "GENE1"), (P["P2"], "GENE2"), (P["P3"], "GENE3"),
    (P["P4"], "GENE4"), (P["P5"], "GENE5"), (P["P6"], "GENE6"),
    (P["PHI"], "DUPGENE"), (P["PLO"], "DUPGENE"),
    (P["P9"], "EGENE9"), (P["P10"], "EGENE10"),
]
ALIAS_ROWS = [
    (P["P9"], "ENSG00000009", "Ensembl_gene"),
    (P["P10"], "999009", "UniProt_DR_GeneID"),
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
check(GRAPH["n_proteins"] == 10, f"fixture graph has all 10 proteins (got {GRAPH['n_proteins']})")


def _graph_from_dense(A, prefix="N"):
    """A fully-connected (or otherwise dense) graph, bypassing STRING files -
    used only for the module-detection cap/giant-warning tests, which are
    pure graph structure and need no identifiers."""
    n_nodes = A.shape[0]
    raw = spx.csr_matrix(A)
    ids = [f"{TAX}.{prefix}{i}" for i in range(n_nodes)]
    return {
        "raw_adjacency": raw,
        "protein_index": {pid: i for i, pid in enumerate(ids)},
        "index_protein": ids,
        "protein_gene": {pid: f"G{i}" for i, pid in enumerate(ids)},
        "n_proteins": n_nodes,
    }


# =============================================================================
print("=== build_string_id_index: symbol/ensembl/entrez tiers, collision determinism ===")
# =============================================================================
idx = txnet.build_string_id_index(GRAPH, aliases_path=ALIASES)

check(idx["by_symbol"].get("GENE3") == P["P3"] and idx["by_symbol"].get("GENE4") == P["P4"],
      "non-colliding symbols map to the right protein")
check(idx["by_symbol"].get("DUPGENE") == P["PHI"],
      f"the collision winner is PHI, the HIGHER-degree protein (got {idx['by_symbol'].get('DUPGENE')})")
check("DUPGENE" in idx["collisions"] and idx["collisions"]["DUPGENE"][0]["kept"] == P["PHI"]
      and idx["collisions"]["DUPGENE"][0]["dropped"] == P["PLO"],
      f"the collision is reported, not silently resolved (got {idx['collisions'].get('DUPGENE')})")
check(idx["n_symbols"] == 9, f"9 distinct symbols after the collision collapses 2 proteins into 1 (got {idx['n_symbols']})")
check(idx["by_ensembl"].get("ENSG00000009") == P["P9"], "the Ensembl alias resolves to P9")
check(idx["by_entrez"].get("999009") == P["P10"], "the Entrez alias resolves to P10")

# The broken stand-in, run for real against the SAME graph, to prove the
# check above is not vacuous.
naive = {}
for pid, name in GRAPH["protein_gene"].items():
    naive[str(name).strip().upper()] = pid
check(naive["DUPGENE"] == P["PLO"],
      f"SENSITIVITY: dict(zip(...))-style last-wins picks PLO - the WRONG, lower-degree protein "
      f"(got {naive['DUPGENE']}) - confirming build_string_id_index's explicit ranking is the thing "
      f"doing the work, not an accident of iteration order")

# =============================================================================
print("=== map_de_table_to_string: tiering is a true partition, symbol beats ensembl ===")
# =============================================================================
def _row(input_id, gene_symbol="", ensembl_id="", entrez_id="", log2fc=2.0, padj=0.001, pvalue=0.0001):
    return {"input_id": input_id, "gene_symbol": gene_symbol, "ensembl_id": ensembl_id,
            "entrez_id": entrez_id, "log2FC": log2fc, "padj": padj, "pvalue": pvalue}


tier_df = pd.DataFrame([
    _row("r1", gene_symbol="GENE1"),
    _row("r2", gene_symbol="nomatch_sym", ensembl_id="ENSG00000009"),
    _row("r3", gene_symbol="nomatch_sym2", entrez_id="999009"),
    _row("r4", gene_symbol="nomatch_sym3", ensembl_id="nope", entrez_id="nope"),
    # gene_symbol matches AND ensembl_id also (coincidentally) matches a
    # different protein's alias - symbol must win, ensembl must never be consulted.
    _row("r5", gene_symbol="GENE4", ensembl_id="ENSG00000009"),
])
tier_map = txnet.map_de_table_to_string(tier_df, GRAPH, aliases_path=ALIASES, taxon_id=TAX)
check(tier_map["n_mapped"] == 4 and tier_map["n_unmapped"] == 1,
      f"4 of 5 rows map, r4 does not (mapped={tier_map['n_mapped']}, unmapped={tier_map['n_unmapped']})")
check(tier_map["n_by_tier"] == {"gene_symbol": 2, "ensembl_id": 1, "entrez_id": 1},
      f"the per-tier counts are a true partition of the mapped rows (got {tier_map['n_by_tier']})")
r5 = tier_map["table"].set_index("input_id").loc["r5"]
check(r5["STRING ID"] == P["P4"] and r5["Mapped via"] == "gene_symbol",
      f"SENSITIVITY: r5's gene_symbol ('GENE4') wins over its ALSO-matching ensembl_id - "
      f"got STRING ID {r5['STRING ID']!r} via {r5['Mapped via']!r} (would be P9/ensembl_id if the "
      f"tier order were broken)")
check(tier_map["low_mapping_warning"] is None, "80% mapping does not trigger the low-mapping warning")

# =============================================================================
print("=== map_de_table_to_string: the <10% refusal and the 10-50% warning gates ===")
# =============================================================================
error_df = pd.DataFrame(
    [_row("m0", gene_symbol="GENE1")] + [_row(f"m{i}", gene_symbol="NOPE") for i in range(1, 20)])
try:
    txnet.map_de_table_to_string(error_df, GRAPH, aliases_path=ALIASES, taxon_id=TAX)
    check(False, "a 5%-mapping table raises rather than silently succeeding")
except ValueError as exc:
    check(TAX in str(exc), f"a 5%-mapping table raises ValueError naming the taxon (got: {exc})")

warn_df = pd.DataFrame(
    [_row(f"w{i}", gene_symbol=g) for i, g in enumerate(["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6"])]
    + [_row(f"w{i}", gene_symbol="NOPE") for i in range(6, 20)])
warn_map = txnet.map_de_table_to_string(warn_df, GRAPH, aliases_path=ALIASES, taxon_id=TAX)
check(abs(warn_map["fraction_mapped"] - 0.30) < 1e-9, f"30% mapping computed correctly (got {warn_map['fraction_mapped']:.3f})")
check(warn_map["low_mapping_warning"] is not None,
      "a 30%-mapping table (between the 10% refusal and the 50% warning line) runs but sets the warning")

# =============================================================================
print("=== build_seed_weights: DEG set matches tx.select_degs; duplicate-protein rows sum ===")
# =============================================================================
seed_df = pd.DataFrame([
    _row("r1", gene_symbol="G1", log2fc=2.0, padj=0.001),
    _row("r2", gene_symbol="G1b", log2fc=2.5, padj=0.002),   # maps to the SAME protein as r1
    _row("r3", gene_symbol="G2", log2fc=-3.0, padj=1e-10),
])
seed_map = {"string_by_row": {"r1": "9998.W1", "r2": "9998.W1", "r3": "9998.W2"}}

for direction in ("all", "up", "down"):
    real_degs = tx.select_degs(seed_df, tx.PADJ_CUTOFF_DEFAULT, tx.LFC_CUTOFF_DEFAULT, direction)
    w = txnet.build_seed_weights(seed_df, seed_map, direction=direction, weighting="equal")
    check(w["n_degs"] == len(real_degs),
          f"direction={direction!r}: n_degs ({w['n_degs']}) matches tx.select_degs directly ({len(real_degs)}), "
          f"so the two modules can never disagree about what a DEG is")

w_all = txnet.build_seed_weights(seed_df, seed_map, direction="all", weighting="equal")
check(w_all["n_degs"] == 3 and w_all["n_seeds"] == 2,
      f"3 DEG rows collapse onto 2 STRING proteins (n_degs={w_all['n_degs']}, n_seeds={w_all['n_seeds']})")
check(abs(w_all["weights"]["9998.W1"] - 2.0) < 1e-9,
      f"SENSITIVITY (documents actual behaviour): two DEG rows mapping to ONE protein are SUMMED "
      f"(1.0+1.0), not deduplicated to 1.0 - 'equal' weighting is per-DEG-row, not per-protein "
      f"(got {w_all['weights']['9998.W1']})")
check(abs(w_all["weights"]["9998.W2"] - 1.0) < 1e-9, "a protein with a single DEG row weighs exactly 1.0")

# =============================================================================
print("=== build_seed_weights: abs_metric is strictly positive, and SYMMETRIC in sign ===")
# =============================================================================
sym_df = pd.DataFrame([
    _row("r1", gene_symbol="G3", log2fc=3.0, padj=1e-10, pvalue=1e-10),
    _row("r2", gene_symbol="G4", log2fc=-3.0, padj=1e-10, pvalue=1e-10),
])
sym_map = {"string_by_row": {"r1": "9998.W1", "r2": "9998.W2"}}
sym_w = txnet.build_seed_weights(sym_df, sym_map, weighting="abs_metric", metric="signed_log_pvalue")
check(all(v > 0 for v in sym_w["weights"].values()), "every abs_metric weight is strictly positive")
check(abs(sym_w["weights"]["9998.W1"] - sym_w["weights"]["9998.W2"]) < 1e-9,
      f"log2FC=+3 and log2FC=-3 at the SAME |p-value| get the SAME weight "
      f"(W1={sym_w['weights']['9998.W1']:.4f}, W2={sym_w['weights']['9998.W2']:.4f}) - direction is expressed "
      f"only through seed SELECTION, never as a sign in the weight")

# SENSITIVITY: this is exactly what the engine's shift branch would break, if
# a caller ever passed a SIGNED weight instead. Demonstrated with the real
# shared function, not a re-implementation.
signed = {"9998.W1": 3.0, "9998.W2": -3.0}
shifted, _ = gnet._map_seeds_to_indices(
    {"9998.W1": "9998.W1", "9998.W2": "9998.W2"}, {"9998.W1": 0, "9998.W2": 1}, signed)
check(abs(shifted[0] - shifted[1]) > 1e-6,
      f"SENSITIVITY: had a signed metric been used directly (skipping build_seed_weights' abs()), "
      f"_map_seeds_to_indices' shift-positive branch would fire and BREAK the symmetry just proven above "
      f"(shifted weights {shifted[0]:.2f} vs {shifted[1]:.2f}, not equal) - this is why every weight this "
      f"module produces is asserted strictly positive")

# =============================================================================
print("=== build_seed_weights: winsorization caps the ratio at exactly 20 ===")
# =============================================================================
ratio_df = pd.DataFrame([
    _row(f"r{i}", gene_symbol=f"G{5+i}", log2fc=lfc, padj=0.001)
    for i, lfc in enumerate([1.0, 1.0, 1.0, 1.0, 100.0])
])
ratio_map = {"string_by_row": {f"r{i}": f"9998.W{i}" for i in range(5)}}
rw = txnet.build_seed_weights(ratio_df, ratio_map, weighting="abs_metric", metric="log2fc")
check(rw["n_winsorized"] == 1, f"exactly the one outlier (log2FC=100) is winsorized (got {rw['n_winsorized']})")
check(abs(rw["weight_ratio"] - 20.0) < 1e-9,
      f"the floor (0.05) and the cap (1.0 after normalizing to the 99th percentile) bound the max/min "
      f"weight ratio at EXACTLY 20 (got {rw['weight_ratio']:.4f})")

# =============================================================================
print("=== detect_modules_capped: the cap keeps the top-N by Z-score, and reports it ===")
# =============================================================================
# Row order deliberately does NOT match Z-score rank - the true top-4
# (P1..P4, Z=9..6) are scattered through the table, not in the first 4 rows.
# Otherwise "top-4 by Z-score" and "first 4 rows" would be indistinguishable
# and the check below would be vacuous.
_cap_row_order = ["P10", "P1", "P9", "P2", "PHI", "P3", "PLO", "P4", "P6", "P5"]
_cap_z = {"P1": 9, "P2": 8, "P3": 7, "P4": 6, "P5": 5, "P6": 4, "PHI": 3, "PLO": 2, "P9": 1, "P10": 0}
cap_hits = pd.DataFrame({
    "STRING ID": [P[k] for k in _cap_row_order],
    "Z-score": [_cap_z[k] for k in _cap_row_order],
})
captured_ids = []
_real_detect = gnet.detect_modules


def _spy_detect(hit_string_ids, graph, **kw):
    captured_ids.append(list(hit_string_ids))
    return _real_detect(hit_string_ids, graph, **kw)


gnet.detect_modules = _spy_detect
try:
    modules, cap_info = txnet.detect_modules_capped(cap_hits, GRAPH, min_module_size=2, max_hits=4)
finally:
    gnet.detect_modules = _real_detect

check(cap_info["capped"] is True and cap_info["n_used_for_modules"] == 4,
      f"a 10-hit table capped at max_hits=4 reports capped=True, n_used=4 (got {cap_info})")
check(set(captured_ids[0]) == {P["P1"], P["P2"], P["P3"], P["P4"]},
      f"the 4 hits actually passed to clustering are the top-4 by Z-score, not the first 4 rows "
      f"(got {captured_ids[0]})")

# =============================================================================
print("=== detect_modules_capped: a single dense cluster sets the giant-module warning ===")
# =============================================================================
n_clique = 6
A = np.ones((n_clique, n_clique)) * 0.9
np.fill_diagonal(A, 0.0)
clique_graph = _graph_from_dense(A, prefix="C")
clique_hits = pd.DataFrame({
    "STRING ID": [f"{TAX}.C{i}" for i in range(n_clique)],
    "Z-score": list(range(n_clique, 0, -1)),
})
clique_modules, clique_cap_info = txnet.detect_modules_capped(clique_hits, clique_graph, min_module_size=2)
check(len(clique_modules) == 1 and clique_modules[0]["size"] == n_clique,
      f"a single fully-connected clique clusters into ONE module of all {n_clique} genes "
      f"(got {[m['size'] for m in clique_modules]})")
check(clique_cap_info["giant_module_warning"] is not None,
      "and that triggers the giant-module warning, since one module then holds 100% of the clustered genes")

# =============================================================================
print("=== run_network_analysis: the small-graph guard refuses BEFORE any mapping work ===")
# =============================================================================
_real_build_network = gnet.build_string_network
gnet.build_string_network = lambda *a, **kw: {"n_proteins": 5}
try:
    try:
        txnet.run_network_analysis(pd.DataFrame([_row("x", gene_symbol="G1")]), TAX)
        check(False, "a 5-protein graph (below MIN_GRAPH_PROTEINS_FOR_NULL) raises rather than running")
    except ValueError as exc:
        check("permutation null" in str(exc),
              f"and the message explains WHY (degree bins degenerate below this size): {exc}")
finally:
    gnet.build_string_network = _real_build_network

# =============================================================================
print("=== enrich_modules_by_symbol: THE UNIVERSE - a recruited gene is counted, not dropped ===")
# =============================================================================
FAKE_GO_SETS = {
    "biological_process": {
        "GO:TEST": {"name": "Test term", "genes": ["GENE1", "GENE2", "GENE3", "GENE4"]},
    },
    "molecular_function": {},
    "cellular_component": {},
    "fetched_at": "2026-01-01T00:00:00",
}
_real_fetch_go_sets = tx.fetch_go_sets
tx.fetch_go_sets = lambda *a, **kw: FAKE_GO_SETS
try:
    # Module M1 = {GENE1, GENE2, GENE3}. GENE1/GENE2 are "seeds" (present in a
    # DE table below); GENE3 is RECRUITED - it would never appear in any DE
    # table, since the whole point of this module is finding genes the user
    # never uploaded.
    modules = [{"module_id": "M1", "string_ids": [P["P1"], P["P2"], P["P3"]], "size": 3}]
    enr = txnet.enrich_modules_by_symbol(modules, GRAPH, TAX, database="GO", min_term_size=1)
finally:
    tx.fetch_go_sets = _real_fetch_go_sets

check(enr["universe_size"] == 9,
      f"universe_size is the GRAPH's distinct symbol count (9, matching build_string_id_index's n_symbols "
      f"above), NOT the DE table's (got {enr['universe_size']})")
block = enr["results"]["M1"]["block"]
check(not block["table"].empty, "module M1 got a real enrichment result, not an empty/skipped block")
top = block["table"].iloc[0]
check(int(top["k"]) == 3 and "GENE3" in top["Overlapping Genes"],
      f"THE HEADLINE CLAIM: the recruited gene GENE3 IS counted in k and listed in 'Overlapping Genes' "
      f"(k={top['k']}, overlap={top['Overlapping Genes']!r})")
check(list(block["table"].columns) == tx.RESULT_COLUMNS,
      f"the emitted table's columns are exactly tx.RESULT_COLUMNS, so every existing ORA chart/downstream "
      f"contract needs no changes (got {list(block['table'].columns)})")

# SENSITIVITY: re-run the SAME real engine call (_run_ora_for_sets) with the
# universe swapped to what the DE-table-based "obvious default" would have
# used - the tested genes only, which by construction excludes the recruited
# gene. This is the exact silent failure this design decision prevents.
wrong_universe = {"GENE1", "GENE2"}  # "genes tested in the DE table" - GENE3 excluded by definition
right_block = tx._run_ora_for_sets(FAKE_GO_SETS["biological_process"], {"GENE1", "GENE2", "GENE3", "GENE4"},
                                   {"GENE1", "GENE2", "GENE3"}, 1, 2000, None)
wrong_block = tx._run_ora_for_sets(FAKE_GO_SETS["biological_process"], wrong_universe,
                                   {"GENE1", "GENE2", "GENE3"}, 1, 2000, None)
right_k = int(right_block["table"].iloc[0]["k"])
wrong_k = int(wrong_block["table"].iloc[0]["k"])
check(right_k == 3 and wrong_k == 2,
      f"SENSITIVITY: with the graph universe k=3 (GENE3 counted); with the DE-table universe k={wrong_k} "
      f"and GENE3 silently VANISHES from its own module's enrichment - proving enrich_modules_by_symbol's "
      f"choice of universe is load-bearing, not cosmetic")
check("GENE3" not in wrong_block["table"].iloc[0]["Overlapping Genes"],
      "confirmed: GENE3 is absent from 'Overlapping Genes' under the wrong universe")

TMP.cleanup()
print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
