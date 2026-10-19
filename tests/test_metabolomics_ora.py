# -*- coding: utf-8 -*-
"""Unit tests for Metabolomics_module.run_ora (compound-LIST input, database
universe) and run_quantitative_enrichment (concentration-table input,
measured-compound universe), offline. Only the network edges are faked - the
KEGG downloads and the universe cache - so the REAL universe builders, the
REAL overview-map exclusion, and the shared Fisher/Wilcoxon engines all run.

THE HEADLINE FIX under test: KEGG reserves map01000-map01099 as global
"chemical structure" overview maps that aggregate compounds across huge
swaths of metabolism rather than naming a real pathway - verified live
2026-09-22 that these outranked a genuinely planted TCA-cycle signal (k==K
on the overview maps: every annotated compound in the universe was also
significant, because their real content is "most of central metabolism").
The fixture plants the identical failure mode, small enough to assert on
exactly, and drives it through the real builder (an earlier version of this
file replaced the builder with a fake that re-implemented the exclusion,
which tested the fake, not the code).

THE SECOND THING under test: universe_scope. A plain compound list has no
"measured" set - if the universe were the user's own list, every term would
contain all of it and Fisher's test would be degenerate - so list tiles test
against the DATABASE's reference set ("database"), while the concentration
tile keeps the narrower measured set ("experiment"). Both must coexist and
each result must say which one it used.

THE THIRD THING under test: KEGG is organism-specific. KEGG serves no
organism-specific compound links (verified live: link/compound/hsa is empty)
and a map's compounds do not vary by organism (hsa00020 == sce00020 ==
map00020), so the organism selects WHICH MAPS are tested - the reference-map
links are one shared download and list/pathway/<org> says which maps an
organism has. Terms are reported under the organism's own ids (hsa00020).
With no KEGG genome, or an ambiguous organism the user has not resolved, KEGG is
EMPTY with a note - never a silent fallback to the generic reference maps. The
user can also pick a KEGG organism EXPLICITLY (a search box like Metabolic
Connection's); an explicit code is used as given and the organism derivation is
not consulted at all.
"""
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Metabolomics_module as mm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def _row(name, key, kegg_id, main_class):
    return {"Name": name, "compound_key": key, "match_tier": "kegg", "class_level_only": False,
            "chebi_id": "", "kegg_id": kegg_id,
            "super_class": "Organic acids", "main_class": main_class, "sub_class": main_class,
            "lm_category": "", "lm_main_class": "", "lm_sub_class": ""}


# Six TCA compounds (one deliberately named with a comma), six unrelated
# "noise" compounds - the compounds the user has.
TCA_NAMES = ["TCA0", "TCA1", "TCA2", "TCA3", "TCA4", "2,3-Comma name"]
TCA_ROWS = [_row(TCA_NAMES[i], f"KEY_TCA{i}", f"C{10000+i}", "TCA acids") for i in range(6)]
NOISE_ROWS = [_row(f"NOISE{i}", f"KEY_NOISE{i}", f"C{20000+i}", "Other acids") for i in range(6)]
MEASURED = pd.DataFrame(TCA_ROWS + NOISE_ROWS, columns=mm.MATCHED_COLUMNS)
LIST_ONLY = pd.DataFrame(TCA_ROWS, columns=mm.MATCHED_COLUMNS)

# Compounds the DATABASE knows about that the user does not have - what makes
# the database-scope universe larger than the measured one.
N_EXTRA = 30
EXTRA_KEYS = [f"KEY_EXTRA{i}" for i in range(N_EXTRA)]
EXTRA_KEGG = [f"C{30000+i}" for i in range(N_EXTRA)]

BY_KEGG = {f"C{10000+i}": f"KEY_TCA{i}" for i in range(6)}
BY_KEGG.update({f"C{20000+i}": f"KEY_NOISE{i}" for i in range(6)})
BY_KEGG.update(dict(zip(EXTRA_KEGG, EXTRA_KEYS)))
BY_KEY = {}
for i in range(6):
    BY_KEY[f"KEY_TCA{i}"] = {"main_class": "TCA acids", "sub_class": "TCA acids"}
    BY_KEY[f"KEY_NOISE{i}"] = {"main_class": "Other acids", "sub_class": "Other acids"}
for k in EXTRA_KEYS:
    BY_KEY[k] = {"main_class": "Other acids", "sub_class": "Other acids"}
CROSSWALK = {"by_kegg": BY_KEGG, "by_chebi": {}, "by_skeleton": {}, "by_key": BY_KEY, "sources": []}

TCA_KEGG = [f"C{10000+i}" for i in range(6)]
NOISE_KEGG = [f"C{20000+i}" for i in range(6)]
# Reference maps: 00020 (TCA), 00640 and 00650 (noise + the database-only
# compounds), plus two OVERVIEW maps holding everything. Human "has" 00020 and
# 00640 but NOT 00650; yeast has 00650 but NOT 00020 - so the same shared link
# file must yield different term sets per organism.
FAKE_LINKS = "\n".join(
    [f"cpd:{c}\tpath:map00020" for c in TCA_KEGG]
    + [f"cpd:{c}\tpath:map00640" for c in NOISE_KEGG + EXTRA_KEGG]
    + [f"cpd:{c}\tpath:map00650" for c in NOISE_KEGG + EXTRA_KEGG]
    + [f"cpd:{c}\tpath:map01062" for c in TCA_KEGG + NOISE_KEGG]
    + [f"cpd:{c}\tpath:map01100" for c in TCA_KEGG + NOISE_KEGG + EXTRA_KEGG]
)
FAKE_NAMES = ("path:map00020\tCitrate cycle (TCA cycle)\npath:map00640\tPropanoate metabolism\n"
              "path:map00650\tButanoate metabolism\n"
              "path:map01062\tBiosynthesis of terpenoids and steroids\npath:map01100\tMetabolic pathways\n")
FAKE_ORG_LISTS = {
    "hsa": ("hsa00020\tCitrate cycle (TCA cycle) - Homo sapiens (human)\n"
            "hsa00640\tPropanoate metabolism - Homo sapiens (human)\n"
            "hsa01062\tBiosynthesis of terpenoids and steroids - Homo sapiens (human)\n"
            "hsa01100\tMetabolic pathways - Homo sapiens (human)\n"),
    "sce": ("sce00650\tButanoate metabolism - Saccharomyces cerevisiae (budding yeast)\n"
            "sce01100\tMetabolic pathways - Saccharomyces cerevisiae (budding yeast)\n"),
}


def _fake_kegg_get(url, verbose=False, **kwargs):
    if "link/pathway/compound" in url:
        return "ok", FAKE_LINKS
    for org, text in FAKE_ORG_LISTS.items():
        if url.endswith("/list/pathway/" + org):
            return "ok", text
    if url.endswith("/list/pathway/nomaps"):
        return "ok", ""
    return "ok", FAKE_NAMES


def _resolver(code, name, note=None, caveat=None, candidates=()):
    return lambda *a, **k: {"organism_code": code, "organism_name": name, "verified": bool(code),
                            "note": note, "caveat": caveat, "candidates": list(candidates)}


# resolve_kegg_organism_code reaches KEGG through Transcriptomics_module's own
# alias of _kegg_safe_get - a SEPARATE binding that faking mm._kegg_safe_get
# does not intercept - so the resolver itself is stubbed (the pattern
# tests/test_proteomics_kegg.py uses). What is under test here is what this
# module does with its answer.
mm.resolve_kegg_organism_code = _resolver("hsa", "Homo sapiens")
HUMAN = {"taxon_id": "9606", "taxon_name": "Homo sapiens"}


mm._kegg_safe_get = _fake_kegg_get
mm._load_universe_cache = lambda *a, **k: None
mm._save_universe_cache = lambda *a, **k: None
_NOT_TESTED = {"terms": {}, "universe_size": 0, "fetched_at": "", "caveat": "not tested here"}
mm.build_reactome_compound_universe = lambda *a, **k: dict(_NOT_TESTED)
mm.build_smpdb_universe = lambda *a, **k: dict(_NOT_TESTED)
mm.build_pathbank_universe = lambda *a, **k: dict(_NOT_TESTED)

print("=== KEGG_COMPOUND_EXCLUDED_MAPS: the overview-map numbers are excluded ===")
check("01062" in mm.KEGG_COMPOUND_EXCLUDED_MAPS, "the chemical-structure-map range (01000-01099) is excluded")
check("01100" in mm.KEGG_COMPOUND_EXCLUDED_MAPS,
      "the shared gene-side overview range (imported from database_catalog) is also excluded")
check("00020" not in mm.KEGG_COMPOUND_EXCLUDED_MAPS, "a real, specific pathway map is NOT excluded")

print("=== names: commas are stripped from NAME: keys (they would tear a member list apart) ===")
check("," not in mm._norm_name("Glucose, 6-phosphate"),
      f"a comma inside a compound name never reaches a compound key (got {mm._norm_name('Glucose, 6-phosphate')!r})")
check(mm._norm_name("1,3-Diaminopropane") != mm._norm_name("13-Diaminopropane"),
      "...but the digits either side of it are not merged together")

print("=== run_ora: THE HEADLINE FIX - overview maps excluded, real pathway wins ===")
result = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, **HUMAN)
table = result["tables"]["KEGG"]["table"]
check(list(result["tables"]) == ["KEGG"], f"one result family per database, no direction suffix (got {list(result['tables'])})")
check(not table.empty, "the KEGG family produced results")
if not table.empty:
    term_ids = set(table["Term ID"])
    check("hsa01062" not in term_ids and "hsa01100" not in term_ids,
          f"the two overview maps never entered the result table at all - excluded before "
          f"testing, not merely outranked (got: {sorted(term_ids)})")
    check(table.iloc[0]["Term ID"] == "hsa00020",
          f"the real TCA cycle pathway is the TOP hit once overview maps are excluded (got {table.iloc[0]['Term ID']!r})")
    check(table.iloc[0]["k"] == 6 and table.iloc[0]["K"] == 6,
          f"the TCA hit's k/K are exactly the 6 list compounds (got k={table.iloc[0]['k']}, K={table.iloc[0]['K']})")
    check(table.iloc[0]["FDR"] < 0.05, f"the planted signal is significant (FDR={table.iloc[0]['FDR']})")

print("=== run_ora: the overlap column shows the user's own names, semicolon-joined ===")
if not table.empty:
    members = table.iloc[0]["Overlapping Compounds"]
    check("Overlapping Compounds" in table.columns and "Overlapping Genes" not in table.columns,
          "the member column is named for compounds, not genes")
    check("KEY_TCA" not in members, f"raw internal compound keys never reach the table (got {members!r})")
    check("2,3-Comma name" in members.split("; "),
          f"a name containing a comma survives intact as ONE member (got {members!r})")
    check(len(members.split("; ")) == 6, f"all six members are listed (got {members.split('; ')})")

print("=== run_ora: summary describes a list run, with no significance split ===")
summ = result["summary"]
check(summ["universe_scope"] == "database", f"summary states its universe scope (got {summ.get('universe_scope')!r})")
check(summ["n_input_compounds"] == 6 and summ["n_matched"] == 6, f"counts match the fixture (got {summ})")
check("n_sig_up" not in summ and "n_sig_all" not in summ and "n_sig_down" not in summ,
      "no up/down/all significance counts - a plain list has no direction")

print("=== universe_scope: database universe is the DATABASE's set, not the user's list ===")
uni_db = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database", **HUMAN)
uni_exp = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="experiment", **HUMAN)
check(uni_db["universe_size"] == 6 + 6 + N_EXTRA,
      f"database scope covers every compound annotated to a real pathway (got {uni_db['universe_size']})")
check(uni_exp["universe_size"] == 6,
      f"experiment scope is confined to the user's own compounds (got {uni_exp['universe_size']})")
check(len(uni_db["terms"]["hsa00640"]["genes"]) == 6 + N_EXTRA and "hsa00640" not in uni_exp["terms"],
      "K for a pathway counts the whole database under database scope; under experiment scope a term "
      "with none of the user's compounds is simply absent")
check(uni_exp["terms"]["hsa00020"]["genes"] == set(LIST_ONLY["compound_key"]),
      "SENSITIVITY: under experiment scope a list-only input makes the TCA term equal the whole "
      "universe (K==N, so Fisher is degenerate) - the reason list tiles must not use it")

print("=== run_ora: a background list narrows the database universe ===")
bg = {f"KEY_TCA{i}" for i in range(6)} | {f"KEY_NOISE{i}" for i in range(6)}
result_bg = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, background_compounds=bg, **HUMAN)
check(result_bg["summary"]["background_applied"] is True, "summary records that a background was applied")
check(result_bg["tables"]["KEGG"]["universe_size"] == 12,
      f"universe narrowed to the 12 background compounds (got {result_bg['tables']['KEGG']['universe_size']})")
check(result["tables"]["KEGG"]["universe_size"] == 6 + 6 + N_EXTRA,
      f"without a background the full database universe is used (got {result['tables']['KEGG']['universe_size']})")

print("=== run_ora: Chemical Class recovers the planted signal from the crosswalk taxonomy ===")
result2 = mm.run_ora(LIST_ONLY, CROSSWALK, ["Chemical Class"], min_term_size=3)
chem = result2["tables"]["Chemical Class"]["table"]
check(not chem.empty and chem.iloc[0]["Term Name"] == "TCA acids",
      f"Chemical Class independently recovers 'TCA acids' (got {chem.iloc[0]['Term Name'] if not chem.empty else None})")

print("=== run_ora: an empty database list produces empty tables, not a crash ===")
check(mm.run_ora(LIST_ONLY, CROSSWALK, [])["tables"] == {}, "an empty database list produces an empty tables dict")

print("=== KEGG is organism-specific: the organism selects which MAPS are tested ===")
check(str(table.iloc[0]["Term Name"]) == "Citrate cycle (TCA cycle)",
      f"terms carry the organism's own name with KEGG's ' - Homo sapiens (human)' suffix stripped "
      f"(got {table.iloc[0]['Term Name']!r})")
check("hsa00650" not in set(uni_db["terms"]),
      "a reference map the organism LACKS (00650 - not in human's list) is not a term")
check(all(tid.startswith("hsa") for tid in uni_db["terms"]),
      f"every term id names the organism (got {sorted(uni_db['terms'])})")
check("hsa01062" not in uni_db["terms"] and "hsa01100" not in uni_db["terms"],
      "the overview maps stay excluded AFTER organism re-keying (the exclusion runs on the reference id; "
      "on 'hsa01100' the 'map' prefix it strips would be gone and they would silently return)")
check("Homo sapiens (hsa)" in (uni_db["caveat"] or "") and "2 of 3 reference maps" in (uni_db["caveat"] or ""),
      f"the caveat names the organism and how many reference maps it has (got {uni_db['caveat']!r})")
check("does not vary by organism" in (uni_db["caveat"] or ""),
      "the caveat does NOT claim species-specific compound content - only which maps are tested")

mm.resolve_kegg_organism_code = _resolver("sce", "Saccharomyces cerevisiae")
uni_sce = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database",
                                          taxon_id="4932", taxon_name="Saccharomyces cerevisiae")
check(set(uni_sce["terms"]) == {"sce00650"},
      f"a second organism yields a DIFFERENT term set from the SAME shared link data "
      f"(yeast: {sorted(uni_sce['terms'])}; human: {sorted(uni_db['terms'])}) - this is the check that "
      f"fails if the builder is ever made organism-blind again")
check(uni_sce["universe_size"] == 6 + N_EXTRA and uni_db["universe_size"] == 6 + 6 + N_EXTRA,
      f"the universe follows the organism's maps (yeast {uni_sce['universe_size']}, human {uni_db['universe_size']})")
mm.resolve_kegg_organism_code = _resolver("hsa", "Homo sapiens")

print("=== KEGG: no genome / ambiguous organism -> EMPTY with a note, never a fallback ===")
mm.resolve_kegg_organism_code = _resolver(
    "", "", note='KEGG has no genome for "Nonexistus".')
res_none = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, taxon_id="0", taxon_name="Nonexistus")
blk = res_none["tables"]["KEGG"]
check(blk["table"].empty and "no genome" in (blk["note"] or "") and "KEGG organism box" in (blk["note"] or ""),
      f"no KEGG genome -> an empty table whose note points at the KEGG organism box (got {blk['note']!r})")
check("MSigDB" not in (blk["note"] or ""), "...and carries no gene-module advice ('try GO or MSigDB')")
check(blk["n_terms_tested"] == 0 and blk["universe_size"] == 0,
      "SENSITIVITY: nothing was tested - no silent fallback to the generic reference maps")

mm.resolve_kegg_organism_code = _resolver(
    "", "", note='"Bacillus subtilis" matches 14 KEGG organisms; pick one explicitly.',
    candidates=[{"organism_code": "bsu", "organism_name": "B. subtilis 168"}])
res_amb = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, taxon_id="1423", taxon_name="Bacillus subtilis")
_amb_note = res_amb["tables"]["KEGG"]["note"] or ""
check(res_amb["tables"]["KEGG"]["table"].empty and "KEGG organism box" in _amb_note and "bsu" in _amb_note,
      f"an ambiguous organism with nothing chosen is empty, with a note that names the candidates and says "
      f"how to fix it (got {_amb_note!r})")
check("MSigDB" not in _amb_note and "try GO" not in _amb_note,
      "the note carries no gene-module advice ('GO and MSigDB...') - it means nothing in Metabolomics")

print("=== KEGG: an EXPLICIT organism (the KEGG search box) is used as given ===")
def _boom(*a, **k):
    raise AssertionError("the organism derivation must not be consulted for an explicit KEGG organism")
mm.resolve_kegg_organism_code = _boom
uni_pick = mm.build_kegg_compound_universe(
    LIST_ONLY, CROSSWALK, universe_scope="database",
    kegg_organism_code="sce", kegg_organism_name="Saccharomyces cerevisiae (budding yeast)")
check(set(uni_pick["terms"]) == {"sce00650"},
      f"an explicit code selects that organism's maps with NO taxon given at all (got {sorted(uni_pick['terms'])})")
check("chosen explicitly" in (uni_pick["caveat"] or "") and "(budding yeast) (sce)" in (uni_pick["caveat"] or ""),
      f"the caveat says it was chosen explicitly and names the organism (got {uni_pick['caveat']!r})")
check(mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database",
                                      kegg_organism_code="SCE")["terms"].keys() == uni_pick["terms"].keys(),
      "the code is case-insensitive, and needs no organism name")

mm.resolve_kegg_organism_code = _resolver(
    "", "", note='"Bacillus subtilis" matches 14 KEGG organisms; pick one explicitly.',
    candidates=[{"organism_code": "bsu", "organism_name": "B. subtilis 168"}])
res_win = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, taxon_id="1423", taxon_name="Bacillus subtilis",
                     kegg_organism_code="hsa", kegg_organism_name="Homo sapiens (human)")
check(not res_win["tables"]["KEGG"]["table"].empty
      and all(tid.startswith("hsa") for tid in res_win["tables"]["KEGG"]["table"]["Term ID"]),
      "an explicit KEGG organism unblocks an ambiguous run: run_ora forwards it through _build_universes to the builder")
res_q = mm.run_quantitative_enrichment(
    MEASURED, CROSSWALK, pd.DataFrame({"Name": TCA_NAMES + [f"NOISE{i}" for i in range(6)],
                                       "log2fc": [2.5] * 6 + [0.0] * 6}),
    ["KEGG"], min_term_size=3, taxon_id="1423", taxon_name="Bacillus subtilis",
    kegg_organism_code="hsa", kegg_organism_name="Homo sapiens (human)")
check(not res_q["tables"]["KEGG"]["table"].empty,
      "run_quantitative_enrichment forwards the explicit organism too")
_saved_get3 = mm._kegg_safe_get
mm._kegg_safe_get = lambda url, verbose=False, **k: (("ok", "") if url.endswith("/list/pathway/zzz") else _saved_get3(url))
res_zzz = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database", kegg_organism_code="zzz")
mm._kegg_safe_get = _saved_get3
check(res_zzz["terms"] == {} and "no pathway maps for organism" in (res_zzz.get("note") or ""),
      f"an explicit code KEGG has no maps for is empty with its own note (got {res_zzz.get('note')!r})")
mm.resolve_kegg_organism_code = _resolver("hsa", "Homo sapiens")

mm.resolve_kegg_organism_code = _resolver(
    "hsa", "Homo sapiens", caveat="KEGG organism hsa was matched by name but records a different taxon.")
res_cav = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, **HUMAN)
check("different taxon" in (res_cav["tables"]["KEGG"]["caveat"] or ""),
      "an organism-resolution caveat survives into a successful result")
mm.resolve_kegg_organism_code = _resolver(
    "", "", note="KEGG requires a resolved organism name.", caveat="an organism-resolution warning")
res_cn = mm.run_ora(LIST_ONLY, CROSSWALK, ["KEGG"], min_term_size=3, taxon_id="", taxon_name="")
check(res_cn["tables"]["KEGG"]["caveat"] == "an organism-resolution warning"
      and "resolved organism name" in res_cn["tables"]["KEGG"]["note"],
      "the caveat travels BESIDE the note on an empty block instead of being dropped (_empty_block(note, caveat))")
mm.resolve_kegg_organism_code = _resolver("hsa", "Homo sapiens")

_saved_get = mm._kegg_safe_get
mm._kegg_safe_get = lambda url, verbose=False, **k: (("ok", "") if url.endswith("/list/pathway/hsa") else _saved_get(url))
res_nomaps = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database", **HUMAN)
check(res_nomaps["terms"] == {} and "no pathway maps" in (res_nomaps.get("note") or ""),
      f"a genome with no listed maps is its own diagnosis (got {res_nomaps.get('note')!r})")
mm._kegg_safe_get = lambda url, verbose=False, **k: (("error", "") if url.endswith("/list/pathway/hsa") else _saved_get(url))
res_down = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database", **HUMAN)
check(res_down["terms"] == {} and "Could not reach KEGG" in (res_down.get("note") or ""),
      f"KEGG being unreachable is reported as that, not as 'no overlap' (got {res_down.get('note')!r})")
mm._kegg_safe_get = _saved_get

# Compounds map to KEGG ids, but every map holding them is absent for the organism.
mm._kegg_safe_get = lambda url, verbose=False, **k: (
    ("ok", "hsa09999\tSomething else - Homo sapiens (human)\n") if url.endswith("/list/pathway/hsa") else _saved_get(url))
res_lack = mm.build_kegg_compound_universe(LIST_ONLY, CROSSWALK, universe_scope="database", **HUMAN)
check(res_lack["terms"] == {} and "none of the KEGG maps" in (res_lack.get("note") or ""),
      f"'your compounds have KEGG ids but none of their maps exist for this organism' is distinct from "
      f"'no compound had a KEGG id' (got {res_lack.get('note')!r})")
mm._kegg_safe_get = _saved_get

print("=== build_kegg_organism_pathway_index: ids are read from the 5-digit TAIL ===")
# KEGG organism codes are not always 3 characters (a fixed-width slice would be
# wrong), and this repo's two KEGG list parsers disagree about the "path:"
# prefix - so the number must come from the tail whatever the prefix or code.
_ids_text = ("path:ecoj00020\tCitrate cycle (TCA cycle) - Escherichia coli J53 (E. coli)\n"
             "ecoj00010\tGlycolysis / Gluconeogenesis - Escherichia coli J53 (E. coli)\n"
             "ecoj00000\tNot a five digit tail? actually five zeros - Escherichia coli J53 (E. coli)\n"
             "ecoj-x\tNo digits at all - Escherichia coli J53 (E. coli)\n")
_saved_get2 = mm._kegg_safe_get
mm._kegg_safe_get = lambda url, verbose=False, **k: ("ok", _ids_text)
_idx = mm.build_kegg_organism_pathway_index("ecoj", force_refresh=True)
mm._kegg_safe_get = _saved_get2
check(set(_idx["names"]) == {"00020", "00010", "00000"},
      f"4-character organism codes and a 'path:' prefix both yield the bare 5-digit number; a row with no "
      f"5-digit tail is dropped (got {sorted(_idx['names'])})")
check(_idx["names"]["00020"] == "Citrate cycle (TCA cycle)",
      f"the organism suffix is stripped from the name (got {_idx['names']['00020']!r})")
check(mm.build_kegg_organism_pathway_index("")["names"] == {}, "an empty organism code yields an empty index, not a request")

print("=== run_quantitative_enrichment: Wilcoxon over the MEASURED universe ===")
rng = np.random.default_rng(0)
stats = pd.DataFrame({
    "Name": TCA_NAMES + [f"NOISE{i}" for i in range(6)],
    "log2fc": [float(rng.normal(2.5, 0.2)) for _ in range(6)] + [float(rng.normal(0, 0.2)) for _ in range(6)],
})
qres = mm.run_quantitative_enrichment(MEASURED, CROSSWALK, stats, ["KEGG"], metric="log2fc", min_term_size=3, **HUMAN)
qtable = qres["tables"]["KEGG"]["table"]
check(qres["summary"]["universe_scope"] == "experiment", "the concentration tile keeps the measured-compound universe")
check(not qtable.empty, "Quant produced results without any significance cutoff")
if not qtable.empty:
    tca_row = qtable[qtable["Term ID"] == "hsa00020"]
    check(not tca_row.empty and tca_row.iloc[0]["Effect (mean difference)"] > 2.0,
          f"the TCA term's effect reflects the planted +2.5 log2FC shift "
          f"(got {tca_row.iloc[0]['Effect (mean difference)'] if not tca_row.empty else None})")
    check("hsa01062" not in set(qtable["Term ID"]), "the overview-map exclusion applies identically to Quant")
    check("Compounds in set" in qtable.columns and "Genes in set" not in qtable.columns,
          "burden columns are renamed for compounds")
    check("2,3-Comma name" in tca_row.iloc[0]["Overlapping Compounds"].split("; "),
          "a comma-containing name survives intact in Quant's member column too")
check(qres["summary"]["n_scored"] == 12, f"every compound with a finite score is counted (got {qres['summary']['n_scored']})")
check(qres["tables"]["KEGG"]["universe_size"] == 12,
      f"Quant's universe is the 12 measured compounds, NOT the database's {6 + 6 + N_EXTRA} "
      f"(got {qres['tables']['KEGG']['universe_size']})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
