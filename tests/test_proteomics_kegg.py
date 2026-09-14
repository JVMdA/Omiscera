# -*- coding: utf-8 -*-
"""Offline unit tests for KEGG in Proteomics ORA + Quantitative Enrichment.

Proteomics is UniProt-accession-native (see Proteomics_module_retriever.py),
so unlike Transcriptomics (gene-symbol-native, needed an Entrez bridge) KEGG
here needs no identifier bridge at all: the retriever's own "accession"
column is already the exact key KEGG's own conv/uniprot table uses. This
file guards the two things that ARE new:

  1. build_kegg_accession_sets' accession -> gene -> pathway join. The
     critical property: going gene -> accession must reuse the SAME
     (already-deduplicated) acc_to_gene map build_kegg_gene_map produces,
     never KEGG's raw, un-deduplicated conv/uniprot data directly - that
     naive direction is exactly the accession fan-out bug that made Genomics
     move its OWN KEGG universe off accession space in the first place
     (verified there: human N=10,873 accessions vs N=9,424 real genes).
  2. That the organism-resolution/confirmation machinery - built for and
     already tested against Transcriptomics - was genuinely GENERALIZED
     (not duplicated) to serve Proteomics too: resolve_kegg_organism_code
     is reused unchanged (no Proteomics-specific copy), and
     _resolve_kegg_for_proteomics produces the same four distinct "nothing
     to report" notes run_ora/run_quantitative_enrichment need.

No network: every KEGG call is monkeypatched, and the on-disk cache is
redirected to a temp dir so a run can neither read nor poison the real
cache (the mistake fixed in tests/test_kegg_gene_map.py in an earlier
session).
"""
import os
import shutil
import sys
import tempfile

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import pandas as pd

import Genomics_module as gm  # noqa: E402
import Proteomics_module as pm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# --- Redirect the on-disk cache for the whole run (build_kegg_gene_map and
# build_kegg_universe both live in Genomics_module and cache via its
# ORA_UNIVERSE_CACHE_DIR).
_TMP_CACHE = tempfile.mkdtemp(prefix="omiscera_proteomics_kegg_test_")
_REAL_CACHE_DIR = gm.ORA_UNIVERSE_CACHE_DIR
gm.ORA_UNIVERSE_CACHE_DIR = _TMP_CACHE

_real_resolve = pm.resolve_kegg_organism_code
_real_gene_map = pm.build_kegg_gene_map
_real_universe = pm.build_kegg_universe


def _reset():
    for name in os.listdir(_TMP_CACHE):
        os.remove(os.path.join(_TMP_CACHE, name))
    pm.resolve_kegg_organism_code = _real_resolve
    pm.build_kegg_gene_map = _real_gene_map
    pm.build_kegg_universe = _real_universe


# ===========================================================================
print("=== build_kegg_accession_sets: the accession -> gene -> pathway join ===")
_reset()

GENE_MAP = {
    "acc_to_gene": {
        "P1": "hsa:1", "P2": "hsa:2",  # both survive -> both land under hsa00010
        "P3": "hsa:3",                  # lands under hsa04110
        "P9": "hsa:9",                  # gene never on any pathway
        # P4 deliberately absent: build_kegg_gene_map already dropped it as
        # ambiguous (mapped to >1 KEGG gene) - never guessed at here either.
    },
    "entrez_to_gene": {}, "symbol_to_gene": {},
    "has_entrez": False, "has_symbols": False,
    "fetched_at": "", "caveat": None,
}
UNIVERSE = {
    "pathways": {
        "hsa00010": {"name": "Glycolysis / Gluconeogenesis", "genes": {"hsa:1", "hsa:2"}},
        "hsa04110": {"name": "Cell cycle", "genes": {"hsa:3"}},
        "hsa99999": {"name": "Untranslatable", "genes": {"hsa:404"}},  # gene with no accession at all
    },
    "universe_size": 3, "fetched_at": "", "caveat": None,
}
pm.build_kegg_gene_map = lambda org, verbose=False, force_refresh=False: GENE_MAP
pm.build_kegg_universe = lambda org, verbose=False, force_refresh=False: UNIVERSE

result = pm.build_kegg_accession_sets("hsa", {"P1", "P2", "P3", "P4", "P5", "P9"})
check(result["pathways"]["hsa00010"]["genes"] == {"P1", "P2"},
      f"two measured accessions sharing a pathway's gene BOTH appear (isoforms are not fan-out - "
      f"see module docstring) (got {result['pathways']['hsa00010']['genes']!r})")
check(result["pathways"]["hsa04110"]["genes"] == {"P3"},
      f"a single-accession pathway resolves correctly (got {result['pathways']['hsa04110']['genes']!r})")
check("hsa99999" not in result["pathways"],
      "SENSITIVITY: a pathway whose only gene has no accession behind it in acc_to_gene is DROPPED, "
      "not kept as an empty set that would inflate 'terms tested'")
check(result["n_measured"] == 6, f"n_measured counts the caller's full measured set (got {result['n_measured']})")
check(result["n_mapped"] == 4,
      f"n_mapped counts measured accessions with a KEGG gene (P1,P2,P3,P9; P4/P5 excluded) "
      f"(got {result['n_mapped']})")
check(result["n_on_pathway"] == 3,
      f"n_on_pathway counts measured accessions whose gene sits on a pathway (P1,P2,P3 - P9's gene "
      f"hsa:9 is on no pathway) (got {result['n_on_pathway']})")
check(result["fetch_failed"] is False and result["caveat"] is None, "a clean fetch carries no caveat")

print("\n=== build_kegg_accession_sets: an accession dropped as ambiguous never appears ===")
_reset()
# P_AMBIG deliberately absent from acc_to_gene (as build_kegg_gene_map would
# do for a real ambiguous accession) - confirm no code path here re-derives
# it from raw conv/uniprot data behind the scenes.
pm.build_kegg_gene_map = lambda org, verbose=False, force_refresh=False: {
    "acc_to_gene": {"P1": "hsa:1"}, "entrez_to_gene": {}, "symbol_to_gene": {},
    "has_entrez": False, "has_symbols": False, "fetched_at": "", "caveat": None,
}
pm.build_kegg_universe = lambda org, verbose=False, force_refresh=False: {
    "pathways": {"hsa00010": {"name": "Glycolysis", "genes": {"hsa:1"}}},
    "universe_size": 1, "fetched_at": "", "caveat": None,
}
result2 = pm.build_kegg_accession_sets("hsa", {"P1", "P_AMBIG"})
check(result2["pathways"]["hsa00010"]["genes"] == {"P1"},
      f"the ambiguous accession never appears in any pathway's gene set (got "
      f"{result2['pathways']['hsa00010']['genes']!r})")
check(result2["n_mapped"] == 1, f"n_mapped correctly excludes the ambiguous accession (got {result2['n_mapped']})")

print("\n=== build_kegg_accession_sets: fetch failure propagates, not cached ===")
_reset()
pm.build_kegg_gene_map = lambda org, verbose=False, force_refresh=False: {
    "acc_to_gene": {}, "entrez_to_gene": {}, "symbol_to_gene": {},
    "has_entrez": False, "has_symbols": False, "fetched_at": "",
    "caveat": "Could not reach KEGG REST for organism 'hsa'. Re-run to retry.",
}
pm.build_kegg_universe = lambda org, verbose=False, force_refresh=False: {
    "pathways": {}, "universe_size": 0, "fetched_at": "", "caveat": None,
}
result3 = pm.build_kegg_accession_sets("hsa", {"P1"})
check(result3["fetch_failed"] is True, "a 'Could not reach' caveat from either half sets fetch_failed")
check("Re-run" in (result3["caveat"] or ""), "and the caveat is retryable, not a 'no data' verdict")

print("\n=== build_kegg_accession_sets: blank organism code short-circuits ===")
_reset()
result4 = pm.build_kegg_accession_sets("", {"P1", "P2"})
check(result4["pathways"] == {} and "No KEGG organism code" in (result4["caveat"] or ""),
      "a blank organism code short-circuits with a caveat and no network call")
check(result4["n_measured"] == 2, "n_measured is still reported even in the short-circuit case")


# ===========================================================================
print("\n=== _resolve_kegg_for_proteomics: the four distinct 'nothing to report' notes ===")

OK_ORG = {"organism_code": "hsa", "organism_name": "Homo sapiens",
          "verified": True, "note": None, "caveat": None, "candidates": []}
HSA_SETS = {
    "pathways": {"hsa00010": {"name": "Glycolysis", "genes": {"P1", "P2"}}},
    "n_measured": 5, "n_mapped": 2, "n_on_pathway": 2,
    "fetch_failed": False, "caveat": None,
}

_reset()
pm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
pm.build_kegg_accession_sets = lambda *a, **k: HSA_SETS
pathways, caveat, note = pm._resolve_kegg_for_proteomics("9606", "Homo sapiens", {"P1", "P2", "P3", "P4", "P5"})
check(note is None and pathways == HSA_SETS["pathways"],
      "a clean resolution + real pathways returns them with no note")
check(caveat and "2 of your 5" in caveat,
      f"and a coverage caveat reporting real mapped/measured counts (got {caveat!r})")

NO_ORG = {"organism_code": "", "organism_name": "", "verified": False,
          "note": "KEGG has no genome for \"Nonexistus\".", "caveat": None, "candidates": []}
pm.resolve_kegg_organism_code = lambda *a, **k: NO_ORG
pathways, caveat, note = pm._resolve_kegg_for_proteomics("0", "Nonexistus", {"P1"})
check(pathways == {} and note and "no genome" in note,
      f"organism unresolvable -> the resolver's own note, unchanged (got {note!r})")

pm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
pm.build_kegg_accession_sets = lambda *a, **k: {
    "pathways": {}, "n_measured": 5, "n_mapped": 0, "n_on_pathway": 0,
    "fetch_failed": False, "caveat": None,
}
pathways, caveat, note = pm._resolve_kegg_for_proteomics("9606", "Homo sapiens", {"P1", "P2", "P3", "P4", "P5"})
check(note and "None of your quantified proteins" in note,
      f"organism resolved but ZERO measured proteins map to any KEGG gene -> its own distinct note, "
      f"not confused with 'no pathway data' (got {note!r})")

pm.build_kegg_accession_sets = lambda *a, **k: {
    "pathways": {}, "n_measured": 5, "n_mapped": 3, "n_on_pathway": 0,
    "fetch_failed": False, "caveat": None,
}
pathways, caveat, note = pm._resolve_kegg_for_proteomics("9606", "Homo sapiens", {"P1", "P2", "P3", "P4", "P5"})
check(note and "3 of your 5" in note and "none of those genes sit on any KEGG pathway" in note,
      f"SENSITIVITY: proteins DID map to KEGG genes, but none of those genes sit on a pathway - "
      f"named as a coverage gap distinct from 'nothing mapped at all' (got {note!r})")

pm.build_kegg_accession_sets = lambda *a, **k: {
    "pathways": {}, "n_measured": 5, "n_mapped": 5, "n_on_pathway": 5,
    "fetch_failed": True, "caveat": "Could not reach KEGG REST for organism 'hsa'. Re-run to retry.",
}
pathways, caveat, note = pm._resolve_kegg_for_proteomics("9606", "Homo sapiens", {"P1", "P2", "P3", "P4", "P5"})
check(note and "Re-run" in note, f"a transient KEGG failure is reported as retryable (got {note!r})")

print("\n=== _resolve_kegg_for_proteomics: unverified/manual-confirmation caveat reaches the caller ===")
UNVERIFIED = dict(OK_ORG, verified=False,
                  caveat="KEGG organism \"hsa\" was matched by name, but KEGG records a different "
                         "NCBI taxon for it than this run's (9606).")
pm.resolve_kegg_organism_code = lambda *a, **k: UNVERIFIED
pm.build_kegg_accession_sets = lambda *a, **k: HSA_SETS
pathways, caveat, note = pm._resolve_kegg_for_proteomics("9606", "Homo sapiens", {"P1", "P2", "P3", "P4", "P5"})
check(note is None and caveat and "different NCBI taxon" in caveat and "2 of your 5" in caveat,
      f"organism-resolution caveat and the coverage caveat are BOTH combined onto the success path "
      f"(got {caveat!r})")

_reset()


# ===========================================================================
print("\n=== run_ora / run_quantitative_enrichment: end-to-end through the real functions ===")

DF = pd.DataFrame({
    "input_id": ["P1", "P2", "P3", "P4", "P5", "P6", "P7"],
    "accession": ["P1", "P2", "P3", "P4", "P5", "P6", "P7"],
    "gene_symbol": ["G1", "G2", "G3", "G4", "G5", "G6", "G7"],
    "protein_name": [""] * 7,
    "log2FC": [3.0, 3.0, -3.0, -3.0, 0.05, 0.05, 0.05],
    "pvalue": [1e-6, 1e-6, 1e-6, 1e-6, 0.9, 0.9, 0.9],
    "padj": [1e-5, 1e-5, 1e-5, 1e-5, 0.95, 0.95, 0.95],
    "status": ["resolved"] * 7,
})
KEGG_SETS = {
    "pathways": {
        "hsa00010": {"name": "Glycolysis", "genes": {"P1", "P2"}},
        "hsa04110": {"name": "Cell cycle", "genes": {"P3", "P4"}},
    },
    "n_measured": 7, "n_mapped": 4, "n_on_pathway": 4,
    "fetch_failed": False, "caveat": None,
}


def _run_ora_kegg(resolve_result, sets_result):
    pm.resolve_kegg_organism_code = lambda *a, **k: resolve_result
    pm.build_kegg_accession_sets = lambda *a, **k: sets_result
    try:
        return pm.run_ora(
            DF, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
            padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, max_term_size=500,
        )
    finally:
        _reset()


result = _run_ora_kegg(OK_ORG, KEGG_SETS)
keys = sorted(result["tables"].keys())
check(keys == ["KEGG (all)", "KEGG (down)", "KEGG (up)"],
      f"one independently-corrected family per direction, same as every other Proteomics ORA "
      f"database (got {keys})")
up = result["tables"]["KEGG (up)"]["table"]
check(len(up) >= 1, f"the up direction found terms (got {len(up)} rows)")
check("hsa00010" in set(up["Term ID"].astype(str)),
      f"glycolysis is enriched among the up-regulated proteins (got {sorted(set(up['Term ID'].astype(str)))})")
down = result["tables"]["KEGG (down)"]["table"]
check("hsa04110" in set(down["Term ID"].astype(str)),
      "and cell cycle among the down-regulated ones - genuinely different results per direction "
      f"(got {sorted(set(down['Term ID'].astype(str)))})")
overlap_col = up.loc[up["Term ID"] == "hsa00010", "Overlapping Proteins"]
check(len(overlap_col) and "P1" in str(overlap_col.iloc[0]) and "G1" in str(overlap_col.iloc[0]),
      f"members show as 'ACCESSION (SYMBOL)', the same normalization every other Proteomics ORA "
      f"database gets (_finalize_member_column) (got {str(overlap_col.iloc[0]) if len(overlap_col) else None!r})")
check("4 of your 7" in (result["tables"]["KEGG (up)"].get("caveat") or ""),
      f"the coverage caveat reaches the stored block (got {result['tables']['KEGG (up)'].get('caveat')!r})")

print("\n=== run_ora: KEGG's empty path carries its own note, not the generic 'no overlap' one ===")
result_no_org = _run_ora_kegg(NO_ORG, KEGG_SETS)
block = result_no_org["tables"]["KEGG (all)"]
check(len(block["table"]) == 0 and "no genome" in (block["note"] or ""),
      f"unresolvable organism -> empty block carrying the resolver's own note (got {block['note']!r})")

print("\n=== run_ora: adding KEGG does not disturb Reactome's own output ===")
pm.resolve_kegg_organism_code = lambda *a, **k: OK_ORG
pm.build_kegg_accession_sets = lambda *a, **k: KEGG_SETS
_real_build_reactome = pm.build_reactome_universe
pm.build_reactome_universe = lambda *a, **k: {
    "pathways": {"R-HSA-1": {"name": "Test Reactome pathway", "genes": {"P1", "P2", "P3"}}},
    "universe_size": 3, "fetched_at": "", "caveat": None,
}
try:
    both = pm.run_ora(
        DF, databases=["Reactome", "KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=1, max_term_size=500,
    )
finally:
    pm.build_reactome_universe = _real_build_reactome
    _reset()
check(len(both["tables"]) == 6,
      f"Reactome x 3 directions + KEGG x 3 directions = 6 families (got {len(both['tables'])})")
r_up = both["tables"]["Reactome (up)"]["table"]
check("R-HSA-1" in set(r_up["Term ID"].astype(str)),
      f"Reactome still finds its own term, unaffected by the new KEGG block "
      f"(got {sorted(set(r_up['Term ID'].astype(str)))})")

print("\n=== run_quantitative_enrichment: single KEGG family, no direction suffix ===")


def _run_quant_kegg(resolve_result, sets_result):
    pm.resolve_kegg_organism_code = lambda *a, **k: resolve_result
    pm.build_kegg_accession_sets = lambda *a, **k: sets_result
    try:
        return pm.run_quantitative_enrichment(
            DF, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
            min_term_size=1, max_term_size=500,
        )
    finally:
        _reset()


quant = _run_quant_kegg(OK_ORG, KEGG_SETS)
check(list(quant["tables"].keys()) == ["KEGG"],
      f"one family, no direction suffix - the signed statistic already captures direction "
      f"(got {list(quant['tables'])})")
qtbl = quant["tables"]["KEGG"]["table"]
check(len(qtbl) == 2, f"both KEGG terms were tested (got {len(qtbl)})")
check("4 of your 7" in (quant["tables"]["KEGG"].get("caveat") or ""),
      "the same coverage caveat reaches Quantitative Enrichment's block too")

quant_no_org = _run_quant_kegg(NO_ORG, KEGG_SETS)
check(len(quant_no_org["tables"]["KEGG"]["table"]) == 0
      and "no genome" in (quant_no_org["tables"]["KEGG"]["note"] or ""),
      "Quantitative Enrichment's unresolvable-organism path carries the same note ORA's does")


# ===========================================================================
print("\n=== API allow-lists actually offer KEGG ===")
import API  # noqa: E402  (imported late: slow, only needed for this check)

check("KEGG" in API.proteomics_engine.PROTEOMICS_ORA_DATABASES,
      "Proteomics ORA/Quant both read from PROTEOMICS_ORA_DATABASES (same list, see "
      "run_quantitative_enrichment's own docstring) - KEGG reaches both from one edit")

app_js = open(os.path.join(_ROOT, "static", "app.js"), encoding="utf-8").read()
check('const PROTEOMICS_ORA_DATABASES = [\n    "GO", "Reactome", "KEGG", "MSigDB",' in app_js,
      "the frontend offers KEGG on both Proteomics config screens (shared PROTEOMICS_ORA_DATABASES)")
check('"Reactome", "KEGG", "UniProt disease"' in app_js,
      "PROTEOMICS_ORA_FAMILY_ORDER places KEGG explicitly (an unlisted family sorts to 999)")

check(hasattr(API, "_kegg_confirmed_code"), "API._kegg_confirmed_code exists (generalized helper)")
check(not hasattr(API, "_transcriptomics_kegg_confirmed_code"),
      "the OLD Transcriptomics-only name is gone, not left as dead code alongside the new one")
routes = [r.path for r in API.app.routes if "/kegg/" in getattr(r, "path", "")]
check("/api/proteomics/kegg/resolve-organism" in routes and "/api/proteomics/kegg/confirm-organism" in routes,
      f"both Proteomics KEGG endpoints are registered (got {sorted(routes)})")
check("/api/transcriptomics/kegg/resolve-organism" in routes and "/api/transcriptomics/kegg/confirm-organism" in routes,
      "and the original Transcriptomics endpoints are still registered too")


shutil.rmtree(_TMP_CACHE, ignore_errors=True)
gm.ORA_UNIVERSE_CACHE_DIR = _REAL_CACHE_DIR

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
