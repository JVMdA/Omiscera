"""Offline unit tests for build_kegg_gene_map (Genomics_module.py) - the
lookup tables (UniProt accession / NCBI Entrez ID / gene symbol -> KEGG gene
ID) the retriever uses to assign each resolved gene a KEGG gene ID.

Written after re-keying build_kegg_universe's pathway sets from UniProt
accessions to KEGG gene IDs (measured live: accession space inflated human
counts via UniProt fan-out - N=10,873 accessions vs N=9,424 real genes,
hsa00010 K=107 vs K=67 - and could never reach the ~9% of pathway genes with
no UniProt cross-reference at all). These tests guard the ONE new invariant a
lookup table needs that a pathway gene-set does not: a key claimed by more
than one KEGG gene must never be guessed at - it is dropped from the map
entirely, in both directions (accession bridge and Entrez bridge) and for the
symbol bridge.

No server needed - pure functions, hand-built fixture strings and a
monkeypatched _kegg_safe_get for the builder checks (no network).
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


print("=== _parse_kegg_ncbi_geneid_conv ===")
geneid_text = "cam:101488213\tncbi-geneid:101488213\ncam:101488214\tncbi-geneid:101488214\nmalformed\n"
geneid_conv = gm._parse_kegg_ncbi_geneid_conv(geneid_text)
check(geneid_conv == {"cam:101488213": "101488213", "cam:101488214": "101488214"},
      f"'ncbi-geneid:' prefix stripped and malformed rows skipped (got {geneid_conv})")

print("=== _parse_kegg_gene_symbols ===")
gene_list_text = (
    "hsa:1\tCDS\t19:...\tA1BG, A1B, ABG; alpha-1-B glycoprotein\n"
    "cam:101488545\tCDS\t1:...\tuncharacterized protein\n"
    "hsa:3098\tCDS\t10:...\tHK1; hexokinase 1\n"
)
symbols = gm._parse_kegg_gene_symbols(gene_list_text)
check(symbols.get("hsa:1") == ["A1BG", "A1B", "ABG"],
      f"multi-symbol field parses all tokens before the ';' (got {symbols.get('hsa:1')})")
check("cam:101488545" not in symbols,
      "a bare description with no symbols (no comma-separated tokens before ';') yields nothing")
check(symbols.get("hsa:3098") == ["HK1"], f"single-symbol field parses correctly (got {symbols.get('hsa:3098')})")

print("=== _parse_kegg_gene_symbols: SENSITIVITY - a description containing ';' must not fake a symbol ===")
# A real KEGG-shaped hazard: a long, space-containing description with an
# embedded ';' but no genuine symbol before it - must NOT be read as symbols.
hazard_text = "cam:999\tCDS\t1:...\tputative uncharacterized protein; contains a colon-like note\n"
hazard = gm._parse_kegg_gene_symbols(hazard_text)
check("cam:999" not in hazard,
      f"a space-containing description before ';' is rejected, not misread as symbols (got {hazard.get('cam:999')})")

too_many_text = "cam:998\tCDS\t1:...\tA, B, C, D, E, F; too many aliases\n"
too_many = gm._parse_kegg_gene_symbols(too_many_text)
check("cam:998" not in too_many, "more than 5 comma-separated tokens is rejected as implausible for real gene aliases")

print("=== _kegg_invert_unique ===")
forward_acc = {"hsa:1": {"P1", "P2"}, "hsa:2": {"P3"}, "hsa:3": {"P2"}}  # P2 claimed by hsa:1 AND hsa:3
inverted = gm._kegg_invert_unique(forward_acc)
check(inverted.get("P1") == "hsa:1", "an accession claimed by exactly one gene is kept")
check(inverted.get("P3") == "hsa:2", "a second, independent single-claim accession is kept")
check("P2" not in inverted,
      "SENSITIVITY: an accession claimed by TWO different genes is dropped from the map entirely, "
      f"never guessed at (got {inverted.get('P2')!r})")

forward_scalar = {"hsa:1": "E1", "hsa:2": "E2"}
inverted_scalar = gm._kegg_invert_unique(forward_scalar)
check(inverted_scalar == {"E1": "hsa:1", "E2": "hsa:2"}, "scalar (non-set) forward values are also handled correctly")

print("=== build_kegg_gene_map: monkeypatched _kegg_safe_get (no network) ===")
_real_kegg_safe_get = gm._kegg_safe_get
# Fake, never-real KEGG organism codes - a genuine organism code here would
# WRITE this fake fixture data to the real on-disk universe cache (force_
# refresh only skips READING a stale cache, it does not skip the save at the
# end of a successful build), silently poisoning that organism for every
# later real use on this machine. Real short KEGG codes are always 3-4 plain
# letters, so these fictitious 6+ char codes can never collide with one.
FAKE_ORG = "xtestorg1"
FAKE_ORG_EMPTY = "xtestorg2"


def _fake_kegg_safe_get(url, verbose=False):
    if "conv/uniprot" in url:
        return "ok", "cam:1\tup:ACC1\ncam:2\tup:ACC2\ncam:3\tup:ACC2\n"  # ACC2 claimed by two genes
    if "conv/ncbi-geneid" in url:
        return "ok", "cam:1\tncbi-geneid:111\ncam:2\tncbi-geneid:222\n"
    if "list/" in url:
        return "ok", "cam:1\tCDS\t1:...\tSYM1; description one\ncam:2\tCDS\t1:...\tSYM2, SYM2B; description two\n"
    return "nodata", ""


gm._kegg_safe_get = _fake_kegg_safe_get
try:
    result = gm.build_kegg_gene_map(FAKE_ORG, force_refresh=True)
    check(result["acc_to_gene"] == {"ACC1": "cam:1"},
          f"ACC2 (claimed by cam:2 AND cam:3) is dropped from acc_to_gene entirely (got {result['acc_to_gene']})")
    check(result["entrez_to_gene"] == {"111": "cam:1", "222": "cam:2"},
          f"unique Entrez IDs both map correctly (got {result['entrez_to_gene']})")
    check(result["symbol_to_gene"].get("SYM1") == "cam:1" and result["symbol_to_gene"].get("SYM2B") == "cam:2",
          f"symbol aliases map to their owning gene (got {result['symbol_to_gene']})")
    check(result["has_entrez"] is True and result["has_symbols"] is True, "has_entrez/has_symbols reflect real data")
    check(result["caveat"] is None, "a successful build carries no caveat")

    print("=== build_kegg_gene_map: no key is ever a plain dict->set that would corrupt on JSON round-trip ===")
    safe = gm._universe_to_json_safe(result)
    restored = gm._universe_from_json_safe(safe)
    check(restored["acc_to_gene"] == result["acc_to_gene"],
          "acc_to_gene survives the universe cache's JSON round-trip as a plain dict (not turned into a set)")
    check(isinstance(restored["acc_to_gene"]["ACC1"], str), "values stay plain strings, not sets, after round-trip")

    print("=== build_kegg_gene_map: an organism with no UniProt data gets a caveat and is not cached ===")
    gm._kegg_safe_get = lambda url, verbose=False: ("nodata", "")
    empty_result = gm.build_kegg_gene_map(FAKE_ORG_EMPTY, force_refresh=True)
    check(empty_result["acc_to_gene"] == {}, "no UniProt data -> empty acc_to_gene")
    check(bool(empty_result["caveat"]), "a non-empty caveat explains why")
finally:
    gm._kegg_safe_get = _real_kegg_safe_get
    # This fixture data must never survive on disk: even under a fictitious
    # org code, a leftover file is exactly the kind of stale test artifact
    # that caused a real bug during this feature's own development (a prior,
    # unguarded version of this test wrote fake data under the REAL "cam"
    # code and silently broke live chickpea results until caught by hand).
    for org in (FAKE_ORG, FAKE_ORG_EMPTY):
        cache_db = f"KEGG_GENEMAP_v{gm.KEGG_GENEMAP_CACHE_VERSION}"
        path = gm._universe_cache_file_path(cache_db, org)
        if os.path.exists(path):
            os.remove(path)

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
