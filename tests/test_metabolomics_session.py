# -*- coding: utf-8 -*-
"""Direct unit tests for two Metabolomics session-state fixes in API.py: no
live server, no network - imports API.py directly and drives it with
hand-built session state and direct in-process job calls, the same way
tests/test_organism_consistency.py already does for the analogous
Transcriptomics bug.

FIX 1 - repeated compound names. The LIST tiles (ORA/Chemical Class/Lipid
Class) never deduplicated an exact-repeated Name: both rows were resolved,
matched, and counted in "Matched", wasting a lookup and inflating the
metric card for no reason - a plain name list has no concept of "two
different measurements" of one compound the way the concentration table's
columns do (see _job_resolve_metabolites_conc's own DELIBERATE choice to
keep repeats, tested here too, so a future change cannot silently start
deduping adduct columns).

FIX 2 - cross-tile staleness. The list and concentration upload paths share
ONE organism (taxon_id/taxon_name), but each resolve job only cleared its
OWN downstream results. Resolving the concentration table for a new
organism left ORA's already-computed, KEGG/Reactome/SMPDB/PathBank-derived
results on record for the organism it was previously resolved under -
silently wrong if ever re-downloaded. The guard added mirrors this one on
each path, and is asymmetric on purpose: Chemical Class / Lipid Class never
use the organism (Metabolomics_module._build_universes), so a taxon change
must NOT touch them.
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import API  # noqa: E402

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
# Fixtures: fast, deterministic stand-ins for the network-touching pieces
# (resolve_compound_table, build_compound_crosswalk, match_compounds) - the
# same technique test_organism_consistency.py uses for resolve_de_table.
# =============================================================================
def _fake_resolve_compound_table(df, verbose=False):
    return pd.DataFrame({
        "Name": df["Name"].tolist(),
        "InChIKey": ["FAKEKEY-" + str(n).upper() for n in df["Name"]],
        "KEGG ID": [""] * len(df), "PubChem ID": [""] * len(df), "SMILES": [""] * len(df),
        "HMDB ID": [""] * len(df), "Formula": [""] * len(df),
        "Mapping_status": ["mapped"] * len(df),
    })


def _fake_match_compounds(resolved_df, crosswalk, user_id_mask=None):
    matched = resolved_df.copy()
    matched["compound_key"] = matched["Name"].astype(str).str.upper()
    matched["match_tier"] = "name"
    matched["class_level_only"] = False
    for col in API.metabolomics_engine.MATCHED_COLUMNS:
        if col not in matched.columns:
            matched[col] = ""
    meta = {
        "n_input": len(matched), "n_matched": len(matched), "n_unmatched": 0,
        "n_matched_by_tier": {"name": len(matched)}, "n_skeleton_only": 0,
        "n_class_level_only": 0,
    }
    return matched, meta


_old_resolve = API.metabolomics_retriever.resolve_compound_table
_old_crosswalk = API.metabolomics_engine.build_compound_crosswalk
_old_match = API.metabolomics_engine.match_compounds
API.metabolomics_retriever.resolve_compound_table = _fake_resolve_compound_table
API.metabolomics_engine.build_compound_crosswalk = lambda verbose=False, progress_callback=None: {}
API.metabolomics_engine.match_compounds = _fake_match_compounds


def _conc_df(names, n_samples=12):
    # Built by POSITION, not a {name: ...} dict comprehension: a dict would
    # collapse two columns that share a literal name down to one, which is
    # exactly the repeated-adduct case this fixture needs to keep as two.
    idx = [f"S{i}" for i in range(n_samples)]
    df = pd.DataFrame([[1.0] * len(names) for _ in range(n_samples)], index=idx, columns=range(len(names)))
    df.columns = list(names)
    return df


try:
    # =========================================================================
    # 1a. LIST path: an exact-duplicate Name is removed, once, and counted.
    # =========================================================================
    print("=== _dedupe_metabolomics_names: direct cases ===")
    deduped, n = API._dedupe_metabolomics_names(pd.DataFrame({"Name": ["Citrate", "Citrate", "Caffeine"]}))
    check(list(deduped["Name"]) == ["Citrate", "Caffeine"] and n == 1,
          f"one exact repeat is removed and counted (got names={list(deduped['Name'])}, n={n})")
    deduped2, n2 = API._dedupe_metabolomics_names(pd.DataFrame({"Name": ["Citrate", "citrate ", "Caffeine"]}))
    check(n2 == 1, f"case/whitespace-insensitive (reuses metabolomics_retriever's own key) (got n={n2})")
    deduped3, n3 = API._dedupe_metabolomics_names(pd.DataFrame({"Name": ["A", "B", "C"]}))
    check(n3 == 0 and len(deduped3) == 3, "no duplicates: nothing removed, nothing miscounted")

    print("=== _job_resolve_metabolites_list: the LIST path dedupes, and reports both counts ===")
    session = {"metabolomics": API._fresh_metabolomics_state()}
    session["metabolomics"]["list_df"] = pd.DataFrame({"Name": ["Citrate", "Citrate", "Caffeine"]})
    API._job_resolve_metabolites_list("job-1", session, "9606", "Homo sapiens")
    matched = session["metabolomics"]["list_matched_df"]
    meta = session["metabolomics"]["list_meta"]
    check(len(matched) == 2, f"the duplicate row never reaches resolution/matching (got {len(matched)} rows, not 3)")
    check(meta["n_duplicate_names_removed"] == 1, f"...and the removal is counted (got {meta.get('n_duplicate_names_removed')})")
    check(meta["n_input_rows"] == 3, f"n_input_rows still reports the ORIGINAL upload count, 3 (got {meta['n_input_rows']})")
    check(meta["n_input"] == 2, f"n_input (from match_compounds) reflects the DEDUPED count, 2 (got {meta['n_input']})")

    # =========================================================================
    # 1b. CONCENTRATION path: the SAME repeated Name is deliberately kept.
    # =========================================================================
    print("=== _job_resolve_metabolites_conc: the CONCENTRATION path never dedupes (a repeat is a real feature) ===")
    session_c = {"metabolomics": API._fresh_metabolomics_state()}
    values = _conc_df(["Malic acid", "Malic acid", "Citric acid"])  # two Malic acid FEATURES
    session_c["metabolomics"]["conc_values_df"] = values
    API._job_resolve_metabolites_conc("job-2", session_c, "9606", "Homo sapiens")
    matched_c = session_c["metabolomics"]["conc_matched_df"]
    check(len(matched_c) == 3,
          f"SENSITIVITY: all three columns are resolved as three rows - the list-path dedup is NOT applied here "
          f"(got {len(matched_c)} rows; a wrongly-shared dedup would give 2)")
    check("n_duplicate_names_removed" not in session_c["metabolomics"]["conc_meta"]
          or session_c["metabolomics"]["conc_meta"].get("n_duplicate_names_removed", 0) == 0,
          "the concentration path's meta reports no removals (dedup truly never ran)")

    # =========================================================================
    # 2. Cross-tile staleness: a taxon change on ONE path clears the OTHER
    #    path's organism-dependent results, and only those.
    # =========================================================================
    print("=== cross-tile staleness: resolving CONC for a new organism clears ORA (not Chemical/Lipid Class) ===")
    session2 = {"metabolomics": API._fresh_metabolomics_state()}
    session2["metabolomics"]["list_df"] = pd.DataFrame({"Name": ["Citrate"]})
    session2["metabolomics"]["conc_values_df"] = _conc_df(["Citrate"])
    API._job_resolve_metabolites_list("job-3", session2, "9606", "Homo sapiens")
    # Pretend every tile has already produced results, under human.
    for key in ("ora_results", "chemclass_results", "lipidclass_results"):
        session2["metabolomics"][key] = {"tables": {}, "summary": {}}
    API._job_resolve_metabolites_conc("job-4", session2, "9606", "Homo sapiens")
    for key in ("ora_results", "chemclass_results", "lipidclass_results"):
        session2["metabolomics"][key] = {"tables": {}, "summary": {}}
    session2["metabolomics"]["quant_results"] = {"tables": {}, "summary": {}}

    # Re-resolving CONC for the SAME organism must not touch anything cross-tile.
    API._job_resolve_metabolites_conc("job-5", session2, "9606", "Homo sapiens")
    check(session2["metabolomics"]["ora_results"] is not None,
          "re-resolving conc for the SAME organism leaves ORA's results alone")

    # Resolving CONC for mouse: ORA is organism-dependent and must be cleared.
    API._job_resolve_metabolites_conc("job-6", session2, "10090", "Mus musculus")
    check(session2["metabolomics"]["ora_results"] is None,
          "SENSITIVITY: resolving conc for a NEW organism clears ORA's now-stale results")
    check(session2["metabolomics"]["chemclass_results"] is not None,
          "...but Chemical Class is untouched - it never used the organism")
    check(session2["metabolomics"]["lipidclass_results"] is not None,
          "...and neither did Lipid Class")

    print("=== cross-tile staleness: resolving the LIST path for a new organism clears Quant, symmetrically ===")
    session2["metabolomics"]["quant_results"] = {"tables": {}, "summary": {}}  # pretend Quant ran again, for mouse
    API._job_resolve_metabolites_list("job-7", session2, "10090", "Mus musculus")
    check(session2["metabolomics"]["quant_results"] is not None,
          "re-resolving the list path for the SAME organism (mouse again) leaves Quant's results alone")
    session2["metabolomics"]["quant_results"] = {"tables": {}, "summary": {}}
    API._job_resolve_metabolites_list("job-8", session2, "10116", "Rattus norvegicus")
    check(session2["metabolomics"]["quant_results"] is None,
          "SENSITIVITY: resolving the list path for a NEW organism clears Quant's now-stale results")

    print("=== the flip side: quant_stats_df is NEVER cleared by an organism change (it does not depend on one) ===")
    session3 = {"metabolomics": API._fresh_metabolomics_state()}
    session3["metabolomics"]["list_df"] = pd.DataFrame({"Name": ["Citrate"]})
    session3["metabolomics"]["conc_values_df"] = _conc_df(["Citrate"])
    API._job_resolve_metabolites_conc("job-9", session3, "9606", "Homo sapiens")
    session3["metabolomics"]["quant_stats_df"] = pd.DataFrame({"Name": ["Citrate"], "log2fc": [1.0]})
    session3["metabolomics"]["quant_stats_meta"] = {"group_a": "A", "group_b": "B"}
    API._job_resolve_metabolites_list("job-10", session3, "10090", "Mus musculus")
    check(session3["metabolomics"]["quant_stats_df"] is not None,
          "a list-path organism change does not touch quant_stats_df/meta (pure per-compound statistics)")
finally:
    API.metabolomics_retriever.resolve_compound_table = _old_resolve
    API.metabolomics_engine.build_compound_crosswalk = _old_crosswalk
    API.metabolomics_engine.match_compounds = _old_match

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
