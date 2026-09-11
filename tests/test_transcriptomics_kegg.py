"""Offline unit tests for KEGG in Transcriptomics ORA + GSEA.

Transcriptomics is gene-SYMBOL-native, while KEGG's pathway membership is
keyed by KEGG gene ID. That makes KEGG's arrival here a two-bridge problem
that nothing else in this module has, and each bridge has a way to fail
silently that these tests exist to prevent:

  1. ORGANISM. KEGG identifies organisms by its own code ("hsa"), and
     publishes no taxon -> code endpoint. resolve_kegg_organism_code derives
     the code from the run's organism NAME and verifies it against KEGG's own
     "TAXONOMY  TAX:<id>" line. A name-only match would be a real hazard:
     verified live, species-level taxon 1423 (Bacillus subtilis) matches 14
     strain-level KEGG genomes and NONE carries taxon 1423 - picking the
     first would silently test the wrong strain. This is the same class of
     bug as commit 27ce033 (KEGG organism divergence in Genomics ORA), so it
     gets the same treatment: make it impossible, then test that it is.

  2. IDENTIFIER. KEGG gene -> Entrez -> this run's gene symbols. The Entrez
     leg comes from KEGG's own conv/ncbi-geneid table; the symbol leg comes
     from the resolved DE table itself (the universe here IS the user's
     measured genes, so no external mapping call is needed - the same bridge
     Reactome already uses, now shared via _entrez_sets_to_symbol_sets).

Every "KEGG returned nothing" path must produce an explanatory note rather
than an empty table with no reason, and must never be mistaken for
_run_ora_for_sets' generic "your DEGs aren't annotated" message - that would
blame the user's gene list for a service or coverage problem. Those notes are
asserted here too.

No network: every KEGG call is monkeypatched, and the module's on-disk cache
is redirected to a temp dir so a test run can neither read nor poison the
real cache (the exact mistake fixed in tests/test_kegg_gene_map.py earlier).
"""
import os
import shutil
import sys
import tempfile

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import pandas as pd

import Transcriptomics_module as tm

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# --- Redirect the on-disk cache for the whole run. resolve_kegg_organism_code
# and fetch_kegg_sets both cache, and a shared cache between tests would let
# one test's fixture leak into another's (and, worse, into real user data).
_TMP_CACHE = tempfile.mkdtemp(prefix="omiscera_kegg_test_")
_REAL_CACHE_DIR = tm.ORA_UNIVERSE_CACHE_DIR
tm.ORA_UNIVERSE_CACHE_DIR = _TMP_CACHE


def _reset_cache():
    """Drop every cached blob between scenarios."""
    for name in os.listdir(_TMP_CACHE):
        os.remove(os.path.join(_TMP_CACHE, name))


# Real list/genome row shapes, including the two that matter: KEGG's yeast
# genome is named without the strain suffix the app's organism carries, and
# B. subtilis has many strain-level genomes.
GENOME_TEXT = (
    "T01001\thsa; Homo sapiens (human)\n"
    "T01002\tmmu; Mus musculus (house mouse)\n"
    "T02819\tcam; Cicer arietinum (chickpea)\n"
    "T00005\tsce; Saccharomyces cerevisiae (budding yeast)\n"
    "T00030\tbsu; Bacillus subtilis subsp. subtilis 168\n"
    "T00624\tbss; Bacillus subtilis subsp. spizizenii 43F\n"
)
GENOME_TAXA = {
    "T01001": "9606", "T01002": "10090", "T02819": "3827",
    "T00005": "559292", "T00030": "224308", "T00624": "1052585",
}


def _fake_kegg_get(genome_text=GENOME_TEXT, taxa=None, conv_text="", conv_status="ok"):
    """Build a _genomics_kegg_get stand-in over hand-built fixtures."""
    taxa = GENOME_TAXA if taxa is None else taxa

    def _get(url, verbose=False):
        if url == tm.KEGG_LIST_GENOME_URL:
            return ("ok", genome_text) if genome_text is not None else ("error", "boom")
        if url.startswith("https://rest.kegg.jp/get/genome:"):
            tnum = url.rsplit(":", 1)[1]
            tax = taxa.get(tnum)
            if tax is None:
                return "nodata", ""
            return "ok", f"ENTRY       {tnum}            Genome\nTAXONOMY    TAX:{tax}\n  LINEAGE   x\n"
        if "conv/ncbi-geneid" in url:
            return conv_status, conv_text
        raise AssertionError("unexpected KEGG URL in an offline test: " + url)

    return _get


# ===========================================================================
print("=== _parse_kegg_genome_list ===")
rows = tm._parse_kegg_genome_list(GENOME_TEXT)
check(len(rows) == 6, f"all 6 rows parsed (got {len(rows)})")
check(rows[0] == ("T01001", "hsa", "Homo sapiens (human)"),
      f"code and name split on the '; ' inside one tab field (got {rows[0]!r})")
check(tm._parse_kegg_genome_list("garbage with no tab\n\n") == [],
      "rows without a tab are skipped, not half-parsed")

check(tm._kegg_strip_common_name("Homo sapiens (human)") == "Homo sapiens",
      "trailing common name stripped")
check(tm._kegg_strip_common_name("Bacillus subtilis subsp. subtilis 168")
      == "Bacillus subtilis subsp. subtilis 168",
      "a name with no parenthesised common name is left alone")


# ===========================================================================
print("\n=== resolve_kegg_organism_code: the exact-match happy path ===")
_reset_cache()
tm._genomics_kegg_get = _fake_kegg_get()
res = tm.resolve_kegg_organism_code("9606", "Homo sapiens")
check(res["organism_code"] == "hsa", f"Homo sapiens/9606 -> 'hsa' (got {res['organism_code']!r})")
check(res["verified"] is True, "taxon-verified against KEGG's own TAXONOMY TAX: line")
check(res["note"] is None and res["caveat"] is None, "a clean resolution carries no note and no caveat")

print("\n=== resolve_kegg_organism_code: genus+species fallback (real yeast case) ===")
_reset_cache()
res = tm.resolve_kegg_organism_code("559292", "Saccharomyces cerevisiae S288C")
check(res["organism_code"] == "sce",
      "the app's strain-suffixed 'Saccharomyces cerevisiae S288C' still matches KEGG's "
      f"unsuffixed genome name via genus+species (got {res['organism_code']!r})")
check(res["verified"] is True, "and it is still taxon-verified, not merely name-matched")

print("\n=== resolve_kegg_organism_code: several strains, none carrying the taxon ===")
_reset_cache()
res = tm.resolve_kegg_organism_code("1423", "Bacillus subtilis")
check(res["organism_code"] == "",
      "SENSITIVITY: species-level taxon 1423 matches 2 strain-level genomes and none "
      "carries taxon 1423 -> resolves NOTHING rather than silently picking a strain "
      f"(got {res['organism_code']!r})")
check(res["note"] and "strain-level" in res["note"],
      f"and says why, naming the strain-level cause (got {res['note']!r})")
check(res["note"] and "bsu" in res["note"] and "bss" in res["note"],
      "listing the actual candidates, so the user can see what it could not choose between")
check(len(res["candidates"]) == 2,
      f"the ambiguous case exposes a STRUCTURED candidate list too, not just the note text - "
      f"this is what the frontend picker/the confirm endpoint read (got {res['candidates']})")
check({c["organism_code"] for c in res["candidates"]} == {"bsu", "bss"},
      f"candidates carry the real organism codes (got {res['candidates']})")

print("\n=== resolve_kegg_organism_code: candidates field is empty outside the ambiguous case ===")
_reset_cache()
clean = tm.resolve_kegg_organism_code("9606", "Homo sapiens")
check(clean["candidates"] == [], "a cleanly-verified resolution carries no candidates")
_reset_cache()
single = tm.resolve_kegg_organism_code("99999", "Nonexistus fakeus")
check(single["candidates"] == [], "'no genome at all' carries no candidates either (nothing to pick)")

print("\n=== resolve_kegg_organism_code: confirmed_organism_code, valid candidate ===")
_reset_cache()
confirmed = tm.resolve_kegg_organism_code("1423", "Bacillus subtilis", confirmed_organism_code="bsu")
check(confirmed["organism_code"] == "bsu",
      f"a confirmed code that IS one of this run's own candidates is accepted (got {confirmed['organism_code']!r})")
check(confirmed["verified"] is False,
      "still reported as unverified - confirmation is a manual override of an ambiguity KEGG's "
      "own taxon record never resolved, not a TAX-verified match")
check(confirmed["caveat"] and "manually confirmed" in confirmed["caveat"],
      f"the caveat says it was a manual confirmation, not automatic (got {confirmed['caveat']!r})")
check(confirmed["candidates"] == [], "once confirmed, no more candidates to pick from")
check(confirmed["note"] is None, "a successful confirmation carries no refusal note")

print("\n=== resolve_kegg_organism_code: confirmed_organism_code, NOT a real candidate ===")
_reset_cache()
rejected = tm.resolve_kegg_organism_code("1423", "Bacillus subtilis", confirmed_organism_code="hsa")
check(rejected["organism_code"] == "",
      "SENSITIVITY: a confirmed code that is NOT among this run's own name-matched candidates "
      "is silently ignored, falling back to the plain ambiguous-refusal result, rather than "
      f"being trusted (got {rejected['organism_code']!r}) - this is what stops confirmation "
      "from becoming an independent organism source unrelated to the actual taxon_name")
check(len(rejected["candidates"]) == 2, "the real candidate list is still reported so the caller can retry")

print("\n=== resolve_kegg_organism_code: confirmed_organism_code is a no-op outside the ambiguous case ===")
_reset_cache()
ignored = tm.resolve_kegg_organism_code("9606", "Homo sapiens", confirmed_organism_code="bsu")
check(ignored["organism_code"] == "hsa" and ignored["verified"] is True,
      "a confirmed code is only ever consulted in the genuinely-ambiguous branch - when the "
      f"organism already resolves cleanly, it is ignored entirely (got {ignored})")

print("\n=== resolve_kegg_organism_code: confirmation is never cached ===")
_reset_cache()
tm.resolve_kegg_organism_code("1423", "Bacillus subtilis", confirmed_organism_code="bsu")
check(not os.listdir(_TMP_CACHE),
      "a confirmed result is never written to the on-disk cache (it's tied to one specific "
      "confirmed_organism_code, and sharing the plain cache key with the unconfirmed case "
      f"would risk serving it back for either one) (cache dir: {os.listdir(_TMP_CACHE)})")
# And a plain (unconfirmed) call afterward still sees the real ambiguity, not a leaked
# confirmed result from the call above.
_reset_cache()
tm.resolve_kegg_organism_code("1423", "Bacillus subtilis", confirmed_organism_code="bsu")
still_ambiguous = tm.resolve_kegg_organism_code("1423", "Bacillus subtilis")
check(still_ambiguous["organism_code"] == "" and len(still_ambiguous["candidates"]) == 2,
      f"a later unconfirmed call for the same organism is unaffected by an earlier confirmation "
      f"(got {still_ambiguous['organism_code']!r}, candidates={still_ambiguous['candidates']})")

print("\n=== run_ora/run_gsea: kegg_confirmed_organism_code reaches a real result ===")
_reset_cache()
DE_BSU = pd.DataFrame({
    "input_id": ["G1", "G2", "G3"], "gene_symbol": ["G1", "G2", "G3"],
    "entrez_id": ["111", "222", "333"], "ensembl_id": [""] * 3,
    "log2FC": [3.0, 3.0, -3.0], "pvalue": [1e-6] * 3, "padj": [1e-5] * 3,
    "status": ["resolved"] * 3,
})
BSU_SETS = {
    "pathways": {"bsu00010": {"name": "Glycolysis", "genes": ["111", "222"]}},
    "n_pathways_before_bridge": 1, "fetch_failed": False, "caveat": None,
}
_fetch_before_confirm_tests = tm.fetch_kegg_sets  # saved/restored locally: _real_fetch
                                                   # (used elsewhere below) isn't defined yet
tm.fetch_kegg_sets = lambda *a, **k: BSU_SETS
try:
    confirmed_ora = tm.run_ora(
        DE_BSU, databases=["KEGG"], taxon_id="1423", taxon_name="Bacillus subtilis",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, max_term_size=500,
        use_length_correction=False, kegg_confirmed_organism_code="bsu",
    )
    unconfirmed_ora = tm.run_ora(
        DE_BSU, databases=["KEGG"], taxon_id="1423", taxon_name="Bacillus subtilis",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, max_term_size=500,
        use_length_correction=False,
    )
    confirmed_gsea = tm.run_gsea(
        DE_BSU, databases=["KEGG"], taxon_id="1423", taxon_name="Bacillus subtilis",
        min_term_size=1, max_term_size=500, n_permutations=100, seed=1,
        kegg_confirmed_organism_code="bsu",
    )
finally:
    tm.fetch_kegg_sets = _fetch_before_confirm_tests
check(len(confirmed_ora["tables"]["KEGG (up)"]["table"]) > 0,
      f"run_ora with a valid confirmed code produces real KEGG results "
      f"(got {len(confirmed_ora['tables']['KEGG (up)']['table'])} rows)")
check(len(unconfirmed_ora["tables"]["KEGG (up)"]["table"]) == 0,
      "SENSITIVITY: the SAME run without a confirmed code stays refused (still ambiguous) - "
      f"proving the confirmed path genuinely changed the outcome, not something else "
      f"(got {len(unconfirmed_ora['tables']['KEGG (up)']['table'])} rows)")
unconfirmed_note = unconfirmed_ora["tables"]["KEGG (up)"]["note"] or ""
check("matches" in unconfirmed_note and "KEGG organisms" in unconfirmed_note,
      f"and the unconfirmed run's note is the real ambiguity refusal, not the confirmed "
      f"path's caveat text (got {unconfirmed_note!r})")
check(len(confirmed_gsea["tables"]["KEGG"]["table"]) > 0,
      f"run_gsea with a valid confirmed code also produces real KEGG results "
      f"(got {len(confirmed_gsea['tables']['KEGG']['table'])} rows)")

print("\n=== resolve_kegg_organism_code: no KEGG genome at all ===")
_reset_cache()
res = tm.resolve_kegg_organism_code("99999", "Nonexistus fakeus")
check(res["organism_code"] == "" and res["note"],
      "an organism KEGG has never sequenced resolves to nothing, with a note")
check("no genome" in res["note"], f"note names the real cause (got {res['note']!r})")

print("\n=== resolve_kegg_organism_code: single name match, taxon disagrees ===")
_reset_cache()
tm._genomics_kegg_get = _fake_kegg_get(
    genome_text="T01001\thsa; Homo sapiens (human)\n", taxa={"T01001": "63221"},
)
res = tm.resolve_kegg_organism_code("9606", "Homo sapiens")
check(res["organism_code"] == "hsa",
      "one unambiguous name match is still used when the taxon disagrees (a strain/"
      "assembly difference is the common cause, and refusing would be unhelpful)")
check(res["verified"] is False and res["caveat"],
      "but it is flagged unverified AND carries a caveat - never accepted silently")
check("9606" in res["caveat"], f"the caveat names this run's taxon (got {res['caveat']!r})")

print("\n=== resolve_kegg_organism_code: transient failure is not cached ===")
_reset_cache()
tm._genomics_kegg_get = _fake_kegg_get(genome_text=None)
res = tm.resolve_kegg_organism_code("9606", "Homo sapiens")
check(res["organism_code"] == "" and "Re-run" in (res["note"] or ""),
      "an unreachable KEGG produces a retryable note, not a 'no data' verdict")
check(not os.listdir(_TMP_CACHE),
      "SENSITIVITY: nothing was cached, so the next run retries instead of being told "
      f"'KEGG is unavailable' for {tm.KEGG_REFRESH_DAYS} days (cache dir: "
      f"{os.listdir(_TMP_CACHE)})")

print("\n=== resolve_kegg_organism_code: no organism name ===")
_reset_cache()
res = tm.resolve_kegg_organism_code("9606", "   ")
check(res["organism_code"] == "" and "resolved organism name" in (res["note"] or ""),
      "a blank organism name is refused with an actionable note, not guessed at")


# ===========================================================================
print("\n=== fetch_kegg_sets: KEGG gene IDs -> Entrez ===")
_reset_cache()

FAKE_UNIVERSE = {
    "pathways": {
        "hsa00010": {"name": "Glycolysis / Gluconeogenesis",
                     "genes": {"hsa:3098", "hsa:2821", "hsa:5230", "hsa:NOCONV"}},
        "hsa04110": {"name": "Cell cycle", "genes": {"hsa:983", "hsa:1017"}},
        "hsa99999": {"name": "Untranslatable", "genes": {"hsa:NOCONV"}},
    },
    "universe_size": 6,
    "caveat": None,
}
CONV_TEXT = (
    "hsa:3098\tncbi-geneid:3098\n"
    "hsa:2821\tncbi-geneid:2821\n"
    "hsa:5230\tncbi-geneid:5230\n"
    "hsa:983\tncbi-geneid:983\n"
    "hsa:1017\tncbi-geneid:1017\n"
)

tm._genomics_build_kegg_universe = lambda org, verbose=False, force_refresh=False: FAKE_UNIVERSE
tm._genomics_kegg_get = _fake_kegg_get(conv_text=CONV_TEXT)

sets = tm.fetch_kegg_sets("hsa")
check(sets["pathways"]["hsa00010"]["genes"] == ["2821", "3098", "5230"],
      f"KEGG gene IDs converted to bare Entrez IDs, sorted (got {sets['pathways']['hsa00010']['genes']!r})")
check(sets["pathways"]["hsa00010"]["name"] == "Glycolysis / Gluconeogenesis",
      "pathway name carried through unchanged")
check("hsa99999" not in sets["pathways"],
      "SENSITIVITY: a pathway whose every gene lacks an Entrez cross-reference is DROPPED, "
      "not kept as an empty set that would inflate 'terms tested'")
check(sets["n_pathways_before_bridge"] == 3,
      f"the pre-bridge pathway count is kept, so 'KEGG has no pathways' and 'none survived "
      f"the bridge' stay distinguishable (got {sets['n_pathways_before_bridge']})")

print("\n=== fetch_kegg_sets: organism with no Entrez cross-references ===")
_reset_cache()
tm._genomics_kegg_get = _fake_kegg_get(conv_text="")
sets = tm.fetch_kegg_sets("xyz")
check(sets["pathways"] == {} and sets["n_pathways_before_bridge"] == 3,
      "pathways exist upstream but none survives the bridge")
check("no NCBI Entrez cross-references" in (sets["caveat"] or ""),
      f"and the caveat says exactly that, not 'no pathways' (got {sets['caveat']!r})")
check(not os.listdir(_TMP_CACHE), "an empty build is never cached")

print("\n=== fetch_kegg_sets: conv endpoint unreachable ===")
_reset_cache()
tm._genomics_kegg_get = _fake_kegg_get(conv_status="error", conv_text="boom")
sets = tm.fetch_kegg_sets("hsa")
check(sets["fetch_failed"] is True, "a transport failure is reported as fetch_failed")
check("Re-run" in (sets["caveat"] or ""), "with a retryable caveat")
check(not os.listdir(_TMP_CACHE), "and is never cached")

print("\n=== fetch_kegg_sets: blank organism code ===")
sets = tm.fetch_kegg_sets("")
check(sets["pathways"] == {} and "No KEGG organism code" in (sets["caveat"] or ""),
      "a blank organism code short-circuits with a caveat and no network call")


# ===========================================================================
print("\n=== _entrez_sets_to_symbol_sets: the shared Reactome/KEGG bridge ===")
UNIVERSE_DF = pd.DataFrame({
    "gene_symbol": ["gapdh", "Eno1", "PKM", "CDK1", "ORPHAN"],
    "entrez_id":   ["2597",  "2023", "5315", "983",  ""],
})
bridged = tm._entrez_sets_to_symbol_sets(
    {"p1": {"name": "P One", "genes": ["2597", "2023", "999999"]},
     "p2": {"name": "P Two", "genes": ["983"]},
     "p3": {"name": "P Three", "genes": ["999999"]},
     "p4": {"name": "P Four", "genes": [""]}},
    UNIVERSE_DF,
)
check(bridged["p1"]["genes"] == ["ENO1", "GAPDH"],
      f"symbols are upper-cased and sorted, unmatched Entrez IDs dropped (got {bridged['p1']['genes']!r})")
check("p3" not in bridged, "a set that translates to nothing is dropped, not kept empty")
check("p4" not in bridged,
      "SENSITIVITY: a blank Entrez ID does NOT match the blank entrez_id of gene ORPHAN - "
      "without the `and e` guard every empty-Entrez gene would join every set carrying a "
      "blank ID")
check(tm._entrez_sets_to_symbol_sets(None, UNIVERSE_DF) == {},
      "None pathways (a failed fetch) bridges to {} rather than raising")


# ===========================================================================
# End-to-end through run_ora / run_gsea, still fully offline.
# ===========================================================================
DE_DF = pd.DataFrame({
    "input_id":    ["GAPDH", "ENO1", "PKM", "CDK1", "CCNB1", "TP53", "ORPHAN"],
    "gene_symbol": ["GAPDH", "ENO1", "PKM", "CDK1", "CCNB1", "TP53", "ORPHAN"],
    "entrez_id":   ["2597", "2023", "5315", "983", "891", "7157", ""],
    "ensembl_id":  [""] * 7,
    "log2FC":      [3.0, 2.5, 2.2, -2.4, -2.6, 0.1, 1.0],
    "pvalue":      [1e-6, 1e-5, 1e-5, 1e-6, 1e-5, 0.9, 0.5],
    "padj":        [1e-5, 1e-4, 1e-4, 1e-5, 1e-4, 0.95, 0.6],
    "status":      ["resolved"] * 7,
})
KEGG_SETS_HSA = {
    "pathways": {
        "hsa00010": {"name": "Glycolysis / Gluconeogenesis", "genes": ["2597", "2023", "5315"]},
        "hsa04110": {"name": "Cell cycle", "genes": ["983", "891", "7157"]},
    },
    "n_pathways_before_bridge": 2, "fetch_failed": False, "caveat": None,
}

_real_resolve = tm.resolve_kegg_organism_code
_real_fetch = tm.fetch_kegg_sets


def _run_ora_kegg(resolve_result, fetch_result):
    tm.resolve_kegg_organism_code = lambda *a, **k: resolve_result
    tm.fetch_kegg_sets = lambda *a, **k: fetch_result
    try:
        return tm.run_ora(
            DE_DF, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
            padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=2, max_term_size=500,
            use_length_correction=False,
        )
    finally:
        tm.resolve_kegg_organism_code = _real_resolve
        tm.fetch_kegg_sets = _real_fetch


OK_ORG = {"organism_code": "hsa", "organism_name": "Homo sapiens",
          "verified": True, "note": None, "caveat": None}

print("\n=== run_ora with KEGG: real terms come back ===")
result = _run_ora_kegg(OK_ORG, KEGG_SETS_HSA)
keys = sorted(result["tables"].keys())
check(keys == ["KEGG (all)", "KEGG (down)", "KEGG (up)"],
      f"one independently-corrected family per direction, same as Reactome (got {keys})")
up = result["tables"]["KEGG (up)"]["table"]
check(len(up) >= 1, f"the up direction found terms (got {len(up)} rows)")
check("hsa00010" in set(up["Term ID"].astype(str)),
      f"glycolysis is enriched among the up-regulated genes (got {sorted(set(up['Term ID'].astype(str)))})")
down = result["tables"]["KEGG (down)"]["table"]
check("hsa04110" in set(down["Term ID"].astype(str)),
      "and cell cycle among the down-regulated ones - the two directions are genuinely "
      f"different results, not one repeated (got {sorted(set(down['Term ID'].astype(str)))})")

print("\n=== run_ora with KEGG: every 'nothing to report' path explains itself ===")
NO_ORG = {"organism_code": "", "organism_name": "", "verified": False,
          "note": "KEGG has no genome for \"Nonexistus fakeus\".", "caveat": None}
res = _run_ora_kegg(NO_ORG, KEGG_SETS_HSA)
block = res["tables"]["KEGG (all)"]
check(len(block["table"]) == 0 and "no genome" in (block["note"] or ""),
      f"unresolvable organism -> empty block carrying the resolver's own note (got {block['note']!r})")

NO_PATHWAYS = {"pathways": {}, "n_pathways_before_bridge": 0, "fetch_failed": False, "caveat": None}
res = _run_ora_kegg(OK_ORG, NO_PATHWAYS)
note = res["tables"]["KEGG (all)"]["note"] or ""
check("no pathway data" in note and "Homo sapiens" in note,
      f"organism resolved but KEGG has no pathways -> its own distinct note (got {note!r})")

NO_BRIDGE = {"pathways": {}, "n_pathways_before_bridge": 12, "fetch_failed": False,
             "caveat": "KEGG publishes no NCBI Entrez cross-references for organism 'xyz'."}
res = _run_ora_kegg(OK_ORG, NO_BRIDGE)
note = res["tables"]["KEGG (all)"]["note"] or ""
check("Entrez cross-references" in note,
      f"pathways exist but none converts to Entrez -> the fetcher's caveat, NOT 'no pathway "
      f"data' (got {note!r})")

FOREIGN = {"pathways": {"mmu00010": {"name": "Glycolysis", "genes": ["111111", "222222"]}},
           "n_pathways_before_bridge": 1, "fetch_failed": False, "caveat": None}
res = _run_ora_kegg(OK_ORG, FOREIGN)
note = res["tables"]["KEGG (all)"]["note"] or ""
check("identifier-mapping gap" in note,
      "SENSITIVITY: pathways exist and convert to Entrez, but none of THIS experiment's "
      "genes carries a matching Entrez ID -> named as an identifier gap, explicitly not as "
      f"'your genes aren't annotated' (got {note!r})")

FAILED = {"pathways": {}, "n_pathways_before_bridge": 0, "fetch_failed": True,
          "caveat": "Could not reach KEGG REST for organism 'hsa'. Re-run to retry."}
res = _run_ora_kegg(OK_ORG, FAILED)
note = res["tables"]["KEGG (all)"]["note"] or ""
check("Re-run" in note,
      f"a transient failure is reported as retryable, never as 'no data' (got {note!r})")

print("\n=== run_ora with KEGG: an unverified organism caveat reaches the results ===")
UNVERIFIED = dict(OK_ORG, verified=False,
                  caveat="KEGG organism \"hsa\" was matched by name, but KEGG records a "
                         "different NCBI taxon for it than this run's (9606).")
res = _run_ora_kegg(UNVERIFIED, KEGG_SETS_HSA)
block = res["tables"]["KEGG (up)"]
check(len(block["table"]) >= 1, "the run still produces results")
check("different NCBI taxon" in (block["note"] or ""),
      f"but the organism caveat is attached to every direction's block (got {block['note']!r})")

print("\n=== run_ora: adding KEGG does not disturb the other databases ===")
tm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
tm.fetch_kegg_sets = lambda *a, **k: KEGG_SETS_HSA
_real_reactome = tm.fetch_reactome_sets
tm.fetch_reactome_sets = lambda *a, **k: {
    "pathways": {"R-HSA-70171": {"name": "Glycolysis", "genes": ["2597", "2023", "5315"]}},
    "matched_via_genus_species_fallback": False,
}
try:
    both = tm.run_ora(
        DE_DF, databases=["Reactome", "KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=2, max_term_size=500,
        use_length_correction=False,
    )
finally:
    tm.resolve_kegg_organism_code = _real_resolve
    tm.fetch_kegg_sets = _real_fetch
    tm.fetch_reactome_sets = _real_reactome

check(len(both["tables"]) == 6,
      f"Reactome x 3 directions + KEGG x 3 directions = 6 families (got {len(both['tables'])})")
r_up = both["tables"]["Reactome (up)"]["table"]
check("R-HSA-70171" in set(r_up["Term ID"].astype(str)),
      "Reactome still works through the now-SHARED _entrez_sets_to_symbol_sets bridge - "
      f"the refactor did not change its output (got {sorted(set(r_up['Term ID'].astype(str)))})")

progress: list = []
tm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
tm.fetch_kegg_sets = lambda *a, **k: KEGG_SETS_HSA
try:
    tm.run_ora(
        DE_DF, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=2, max_term_size=500,
        use_length_correction=False,
        progress_callback=lambda step, total, label: progress.append((step, total)),
    )
finally:
    tm.resolve_kegg_organism_code = _real_resolve
    tm.fetch_kegg_sets = _real_fetch
check(progress and progress[-1][0] == progress[-1][1] == 3,
      f"steps_total accounts for KEGG exactly, so the progress bar reaches 1.0 and no "
      f"further (got {progress[-1] if progress else None})")


# ===========================================================================
print("\n=== run_gsea with KEGG ===")


def _run_gsea_kegg(resolve_result, fetch_result, databases=("KEGG",)):
    tm.resolve_kegg_organism_code = lambda *a, **k: resolve_result
    tm.fetch_kegg_sets = lambda *a, **k: fetch_result
    try:
        return tm.run_gsea(
            DE_DF, databases=list(databases), taxon_id="9606", taxon_name="Homo sapiens",
            min_term_size=2, max_term_size=500, n_permutations=200, seed=7,
        )
    finally:
        tm.resolve_kegg_organism_code = _real_resolve
        tm.fetch_kegg_sets = _real_fetch


res = _run_gsea_kegg(OK_ORG, KEGG_SETS_HSA)
check(list(res["tables"].keys()) == ["KEGG"],
      f"one family, no direction suffix - NES sign carries direction (got {list(res['tables'])})")
table = res["tables"]["KEGG"]["table"]
check(len(table) == 2, f"both KEGG terms were tested (got {len(table)})")
nes_by_term = dict(zip(table["Term ID"].astype(str), table["NES"]))
check(nes_by_term.get("hsa00010", 0) > 0 and nes_by_term.get("hsa04110", 0) < 0,
      "glycolysis (up-regulated genes) gets a positive NES and cell cycle (down-regulated) "
      f"a negative one - the sign really tracks the ranking (got {nes_by_term})")

res = _run_gsea_kegg(NO_ORG, KEGG_SETS_HSA)
check(len(res["tables"]["KEGG"]["table"]) == 0 and "no genome" in (res["tables"]["KEGG"]["note"] or ""),
      "GSEA's unresolvable-organism path carries the same explanatory note as ORA's")

res = _run_gsea_kegg(UNVERIFIED, KEGG_SETS_HSA)
check("different NCBI taxon" in (res["tables"]["KEGG"]["note"] or ""),
      "and the unverified-organism caveat also reaches GSEA's block")

print("\n=== run_gsea: KEGG is covered by the 'fewer than 2 ranked genes' early-out ===")
EMPTY_DE = DE_DF.iloc[:1].copy()
tm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
tm.fetch_kegg_sets = lambda *a, **k: KEGG_SETS_HSA
try:
    res = tm.run_gsea(EMPTY_DE, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
                      min_term_size=2, n_permutations=50, seed=7)
finally:
    tm.resolve_kegg_organism_code = _real_resolve
    tm.fetch_kegg_sets = _real_fetch
check("KEGG" in res["tables"] and "Fewer than 2 genes" in (res["tables"]["KEGG"]["note"] or ""),
      "SENSITIVITY: KEGG gets its own labelled empty block in the early-out, rather than "
      f"vanishing from the results entirely (got {list(res['tables'])})")


# ===========================================================================
print("\n=== API allow-lists actually offer KEGG ===")
import API  # noqa: E402  (imported late: it is slow and only needed for this check)

check("KEGG" in API.TRANSCRIPTOMICS_ORA_DATABASES,
      "Transcriptomics ORA accepts KEGG (without this the configure endpoint silently "
      "drops it from the requested databases)")
check("KEGG" in API.TRANSCRIPTOMICS_GSEA_DATABASES, "Transcriptomics GSEA accepts KEGG")

app_js = open(os.path.join(_ROOT, "static", "app.js"), encoding="utf-8").read()
check('const TRANSCRIPTOMICS_ORA_DATABASES = ["GO", "Reactome", "KEGG", "MSigDB", "Custom database"];' in app_js,
      "the frontend offers KEGG on both Transcriptomics config screens")
check('"Reactome", "KEGG",' in app_js,
      "and both result-family orderings place KEGG explicitly (an unlisted family sorts to "
      "999 and lands after 'Custom database')")


shutil.rmtree(_TMP_CACHE, ignore_errors=True)
tm.ORA_UNIVERSE_CACHE_DIR = _REAL_CACHE_DIR

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
