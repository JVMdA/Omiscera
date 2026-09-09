"""Offline unit tests for the KEGG PATHWAY universe builder added to
Genomics_module.py (Genomics ORA + Gene-set Burden Test).

These guard the parsing logic in isolation from network access, and exist
because KEGG's own REST responses contain two traps that a naive
implementation would get wrong silently:

  1. Pathway names carry a " - <Organism name>" suffix, but several REAL
     pathway names also contain an internal " - " (e.g. "Glycosaminoglycan
     biosynthesis - keratan sulfate", "Drug metabolism - cytochrome P450") -
     verified against the real human pathway list. A naive
     rsplit(" - ", 1) mangles those into the wrong name. The fix strips only
     the single most-common " - <tail>" across all rows for that organism.
  2. One KEGG gene can map to several UniProt accessions (verified live:
     hsa:1 -> P04217 AND V9HWD8). Every accession must be kept, not just one -
     see build_kegg_universe's own docstring for why collapsing would be
     wrong, not just less thorough.

No server needed - pure functions, hand-built fixture strings and a
monkeypatched build_kegg_universe for the dispatcher checks (no network).
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


print("=== _parse_kegg_pathway_list: suffix stripping ===")
pathway_text = (
    "hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)\n"
    "hsa00533\tGlycosaminoglycan biosynthesis - keratan sulfate - Homo sapiens (human)\n"
    "hsa01100\tMetabolic pathways - Homo sapiens (human)\n"
)
names = gm._parse_kegg_pathway_list(pathway_text)
check(names.get("hsa00010") == "Glycolysis / Gluconeogenesis",
      f"internal ' / ' survives suffix stripping (got {names.get('hsa00010')!r})")
check(names.get("hsa00533") == "Glycosaminoglycan biosynthesis - keratan sulfate",
      "SENSITIVITY: internal ' - ' survives - the exact case a naive rsplit(' - ', 1) "
      f"would mangle into 'Glycosaminoglycan biosynthesis' (got {names.get('hsa00533')!r})")
check(names.get("hsa01100") == "Metabolic pathways",
      f"a name with only the organism suffix strips cleanly (got {names.get('hsa01100')!r})")

no_suffix_text = "p1\tFoo\np2\tBar\n"
check(gm._parse_kegg_pathway_list(no_suffix_text) == {"p1": "Foo", "p2": "Bar"},
      "no common suffix across rows -> nothing stripped")

check(gm._parse_kegg_pathway_list("") == {}, "empty input -> empty dict")

print("=== _parse_kegg_gene_pathway_links ===")
link_text = "hsa:10327\tpath:hsa00010\nhsa:124\tpath:hsa00010\nmalformed_line_no_tab\n"
links = gm._parse_kegg_gene_pathway_links(link_text)
check(links == [("hsa:10327", "hsa00010"), ("hsa:124", "hsa00010")],
      f"'path:' prefix stripped and malformed rows skipped (got {links})")

print("=== _parse_kegg_uniprot_conv ===")
conv_text = "hsa:1\tup:P04217\nhsa:1\tup:V9HWD8\nhsa:2\tup:P00001\n"
conv = gm._parse_kegg_uniprot_conv(conv_text)
check(conv.get("hsa:1") == {"P04217", "V9HWD8"},
      f"'up:' prefix stripped and gene->accession fan-out preserved as a set (got {conv.get('hsa:1')})")
check(conv.get("hsa:2") == {"P00001"}, "single-accession gene parses to a one-element set")

print("=== universe assembly: a pathway whose every gene lacks a UniProt mapping is dropped ===")
links2 = [("g1", "p1"), ("g2", "p1"), ("g3", "p2")]  # g3 has no UniProt mapping below
conv2 = {"g1": {"A1"}, "g2": {"A2"}}
pathways: dict = {}
annotated: set = set()
for gene, pid in links2:
    accs = conv2.get(gene)
    if not accs:
        continue
    bucket = pathways.setdefault(pid, {"name": pid, "genes": set()})
    bucket["genes"] |= accs
    annotated |= accs
check("p1" in pathways and pathways["p1"]["genes"] == {"A1", "A2"},
      "pathway with mapped genes keeps their union of accessions")
check("p2" not in pathways,
      "pathway whose only gene has no UniProt mapping is dropped entirely, not kept empty")
check(len(annotated) == 2, f"universe_size counts the accession union, not the link-row count (got {len(annotated)})")

print("=== cache round-trip: 'genes' keys survive as sets (not lists) ===")
synthetic_universe = {
    "pathways": {"p1": {"name": "Test Pathway", "genes": {"A1", "A2"}}},
    "universe_size": 2, "fetched_at": "2026-01-01T00:00:00+00:00", "caveat": None,
}
safe = gm._universe_to_json_safe(synthetic_universe)
check(isinstance(safe["pathways"]["p1"]["genes"], list),
      "to_json_safe converts the set to a plain list for JSON serialization")
restored = gm._universe_from_json_safe(safe)
check(isinstance(restored["pathways"]["p1"]["genes"], set),
      "from_json_safe rehydrates the 'genes' key back into a set")
check(restored["pathways"]["p1"]["genes"] == {"A1", "A2"}, "round-trip preserves the exact accession set")

print("=== run_ora_for_database('KEGG', ...): correct shape and k/K/n/N in KEGG-gene-ID space (no network) ===")
# Note: the universe is keyed on KEGG gene IDs (org:12345-style), NOT UniProt
# accessions - this is the re-keying this file's own build_kegg_universe
# tests above already exercise via _parse_kegg_gene_pathway_links. list_kegg_ids
# (not list_accessions, which KEGG no longer reads at all) is what the
# dispatcher now matches against.
_real_build_kegg_universe = gm.build_kegg_universe
_synthetic = {
    "pathways": {
        "p1": {"name": "Test Pathway 1", "genes": {"org:1", "org:2", "org:3", "org:4", "org:5"}},
        "p2": {"name": "Test Pathway 2", "genes": {"org:6", "org:7"}},
    },
    "universe_size": 7, "fetched_at": "x", "caveat": None,
}
gm.build_kegg_universe = lambda *a, **k: _synthetic
try:
    ora_result = gm.run_ora_for_database(
        "KEGG", "hsa", [], list_kegg_ids=["org:1", "org:2", "org:3"],
        kegg_id_labels={"org:1": "GENE_A", "org:2": "GENE_B", "org:3": "GENE_C"},
    )
    check(set(ora_result.keys()) == {"KEGG"}, f"ORA result keyed exactly {{'KEGG'}} (got {set(ora_result.keys())})")
    table = ora_result["KEGG"]["table"]
    row = table[table["Term ID"] == "p1"].iloc[0]
    check(int(row["k"]) == 3 and int(row["K"]) == 5 and int(row["n"]) == 3 and int(row["N"]) == 7,
          f"k/K/n/N computed correctly against the synthetic universe (got k={row['k']}, K={row['K']}, "
          f"n={row['n']}, N={row['N']})")
    check(row["Overlapping Genes"] == "GENE_A, GENE_B, GENE_C",
          f"SENSITIVITY: Overlapping Genes shows the user's OWN identifiers via kegg_id_labels, not raw "
          f"KEGG gene IDs like 'org:1' (got {row['Overlapping Genes']!r})")

    print("=== run_ora_for_database('KEGG', ...): a gene with no label falls back to its raw ID, never dropped ===")
    ora_result_partial = gm.run_ora_for_database(
        "KEGG", "hsa", [], list_kegg_ids=["org:1", "org:2", "org:3"],
        kegg_id_labels={"org:1": "GENE_A"},  # org:2 and org:3 deliberately unlabeled
    )
    row_partial = ora_result_partial["KEGG"]["table"][ora_result_partial["KEGG"]["table"]["Term ID"] == "p1"].iloc[0]
    check(row_partial["Overlapping Genes"] == "GENE_A, org:2, org:3",
          f"an unlabeled gene falls back to its raw KEGG ID rather than vanishing from the list "
          f"(got {row_partial['Overlapping Genes']!r})")

    print("=== run_burden_for_database('KEGG', ...): correct shape and relabeling (no network) ===")
    kegg_scores = {"org:1": 1.0, "org:2": 1.0, "org:3": 1.0, "org:4": 0.0, "org:5": 0.0, "org:6": 0.0, "org:7": 0.0}
    burden_result = gm.run_burden_for_database(
        "KEGG", "hsa", {}, score_by_kegg_id=kegg_scores,
        kegg_id_labels={"org:1": "GENE_A", "org:2": "GENE_B", "org:3": "GENE_C"},
    )
    check(set(burden_result.keys()) == {"KEGG"}, f"Burden result keyed exactly {{'KEGG'}} (got {set(burden_result.keys())})")
    check(list(burden_result["KEGG"]["table"].columns) == gm.BURDEN_RESULT_COLUMNS,
          "Burden result table uses the standard BURDEN_RESULT_COLUMNS shape")
    burden_row = burden_result["KEGG"]["table"][burden_result["KEGG"]["table"]["Term ID"] == "p1"].iloc[0]
    # Burden's "Genes" lists EVERY scored gene in the term (not just an
    # "input" subset like ORA's Overlapping Genes) - p1 has 5 members, all
    # scored, so all 5 appear: the 3 labeled ones plus org:4/org:5 falling
    # back to their raw ID (no label supplied for those two).
    check(burden_row["Genes"] == "GENE_A, GENE_B, GENE_C, org:4, org:5",
          f"Burden's 'Genes' column is also relabeled via kegg_id_labels, with unlabeled genes falling "
          f"back to their raw KEGG ID rather than vanishing (got {burden_row['Genes']!r})")
finally:
    gm.build_kegg_universe = _real_build_kegg_universe

print("=== unknown database still raises (KEGG branch didn't swallow the fallthrough) ===")
try:
    gm.run_ora_for_database("NOT_A_REAL_DATABASE", "hsa", ["A1"])
    check(False, "run_ora_for_database should have raised ValueError")
except ValueError:
    check(True, "run_ora_for_database still raises ValueError for an unknown database")

try:
    gm.run_burden_for_database("NOT_A_REAL_DATABASE", "hsa", {"A1": 1.0})
    check(False, "run_burden_for_database should have raised ValueError")
except ValueError:
    check(True, "run_burden_for_database still raises ValueError for an unknown database")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
