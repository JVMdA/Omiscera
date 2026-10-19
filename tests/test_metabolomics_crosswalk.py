# -*- coding: utf-8 -*-
"""Unit tests for Metabolomics_module's identifier-bridging layer, offline
(every bulk fetch monkeypatched to a small fixture written to a temp dir).

THE HEADLINE TESTS of this milestone. The crosswalk is what every downstream
result depends on: if a compound reaches the wrong entry, or silently fails to
reach one, every table below it is wrong in a way no statistical test can
reveal. So these verify the join ORDER, the fallback TIERS, and - most
importantly - the two places where a convenient match would be a dishonest
one: stereochemistry-insensitive skeleton matching, and lipid sum-composition
names that cannot support a pathway claim.
"""
import json
import os
import sys
import tempfile
import zipfile

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


# =============================================================================
# Fixtures. Deliberately small and hand-checkable.
#
#   PALMITIC   full record, present in BOTH sources (merge path)
#   CHEBIONLY  reachable ONLY via a ChEBI ID - the Reactome case, which has no
#              other route in
#   LIPIDONLY  reachable ONLY via an LM_ID
#   STEREO_A / STEREO_B  two genuine stereoisomers sharing a skeleton, so a
#              skeleton match must never be preferred over an exact one
# =============================================================================
_TMP = tempfile.mkdtemp(prefix="omiscera_metab_test_")

_SDF_RECORDS = [
    # (LM_ID, NAME, CATEGORY, MAIN_CLASS, SUB_CLASS, INCHIKEY, PUBCHEM, KEGG, HMDB, CHEBI, SYNONYMS)
    ("LMFA01010001", "Palmitic acid", "Fatty Acyls [FA]", "Fatty Acids and Conjugates [FA01]",
     "Straight chain fatty acids [FA0101]", "IPCSVZSSVZVIGE-UHFFFAOYSA-N", "985", "C00249",
     "HMDB0000220", "15756", "Cetylic acid; Palmitate"),
    ("LMGP01010999", "Lipid Only Compound", "Glycerophospholipids [GP]",
     "Glycerophosphocholines [GP01]", "", "LIPIDONLYKEY-UHFFFAOYSA-N", "", "", "", "", ""),
    ("LMST00000001", "Stereo A", "Sterol Lipids [ST]", "Sterols [ST01]", "",
     "SKELSHAREDXX-AAAAAAAAAA-N", "1111", "", "", "", ""),
    ("LMST00000002", "Stereo B", "Sterol Lipids [ST]", "Sterols [ST01]", "",
     "SKELSHAREDXX-BBBBBBBBBB-N", "2222", "", "", "", ""),
]

_sdf_path = os.path.join(_TMP, "lipidmaps.sdf.zip")
_sdf_text = []
for lm_id, name, cat, main, sub, ikey, cid, kegg, hmdb, chebi, syn in _SDF_RECORDS:
    _sdf_text.append(lm_id + "\n  TEST\n\nM  END\n")
    for tag, value in (("LM_ID", lm_id), ("NAME", name), ("CATEGORY", cat),
                       ("MAIN_CLASS", main), ("SUB_CLASS", sub), ("INCHI_KEY", ikey),
                       ("PUBCHEM_CID", cid), ("KEGG_ID", kegg), ("HMDB_ID", hmdb),
                       ("CHEBI_ID", chebi), ("SYNONYMS", syn)):
        if value:
            _sdf_text.append("> <%s>\n%s\n\n" % (tag, value))
    _sdf_text.append("$$$$\n")
with zipfile.ZipFile(_sdf_path, "w") as zf:
    zf.writestr("structures.sdf", "".join(_sdf_text))

_refmet_path = os.path.join(_TMP, "refmet.json")
with open(_refmet_path, "w", encoding="utf-8") as fh:
    json.dump({
        "Row1": {"name": "Palmitic acid", "pubchem_cid": "985", "chebi_id": "15756",
                 "kegg_id": "C00249", "hmdb_id": "HMDB0000220",
                 "inchi_key": "IPCSVZSSVZVIGE-UHFFFAOYSA-N", "refmet_id": "RM1",
                 "super_class": "Fatty Acyls", "main_class": "Fatty acids",
                 "sub_class": "Long-chain fatty acids"},
        "Row2": {"name": "ChEBI Only Compound", "chebi_id": "99999", "refmet_id": "RM2",
                 "super_class": "Organic acids", "main_class": "TCA acids",
                 "sub_class": "TCA acids"},
        "Row3": {"name": "PC 34:1", "pubchem_cid": "5497103", "chebi_id": "64517",
                 "inchi_key": "SHORTHANDKEY-UHFFFAOYSA-N", "refmet_id": "RM3",
                 "super_class": "Glycerophospholipids", "main_class": "PC", "sub_class": "PC"},
    }, fh)


def _fake_ensure(url, local_name, refresh_days=30, verbose=False, force_refresh=False, timeout=300):
    if "lipidmaps" in url.lower():
        return _sdf_path
    if "refmet" in url.lower():
        return _refmet_path
    return None            # ChEBI: simulate the unreachable-host case


mm._ensure_bulk_file = _fake_ensure
mm._ensure_chebi_files = lambda verbose=False, force_refresh=False: (None, None)
# Never touch the real on-disk universe cache from a test - a previous run of
# the real app must not be able to satisfy (or be corrupted by) these fixtures.
mm._load_universe_cache = lambda database, taxon_id: None
mm._save_universe_cache = lambda database, taxon_id, data: None

CROSSWALK = mm.build_compound_crosswalk(force_refresh=True)

# =============================================================================
print("=== crosswalk build: sources merge, missing source degrades honestly ===")
# =============================================================================
check(set(CROSSWALK["sources"]) == {"LIPID MAPS", "RefMet"},
      f"only the reachable sources are reported as contributing (got {CROSSWALK['sources']})")
check(CROSSWALK["caveat"] and "ChEBI" in CROSSWALK["caveat"],
      f"the unreachable source is NAMED in the caveat rather than silently dropped "
      f"(got {CROSSWALK['caveat']!r})")

palmitic_key = CROSSWALK["by_inchikey"].get("IPCSVZSSVZVIGE-UHFFFAOYSA-N")
palmitic = CROSSWALK["by_key"].get(palmitic_key, {})
check(palmitic.get("lm_main_class") == "Fatty Acids and Conjugates"
      and palmitic.get("main_class") == "Fatty acids",
      "a compound present in BOTH sources merges into one entry carrying BOTH taxonomies "
      f"(got lipid={palmitic.get('lm_main_class')!r}, chemical={palmitic.get('main_class')!r})")
check(palmitic.get("lm_main_class") and "[" not in palmitic.get("lm_main_class"),
      "LIPID MAPS' bracketed class codes are stripped from the display label")
check(CROSSWALK["by_name"].get("cetylic acid") == palmitic_key,
      "a SYNONYM resolves to the same compound as its primary name")

# =============================================================================
print("=== match tiers: order, and the two routes with no alternative ===")
# =============================================================================
df = pd.DataFrame([
    {"Name": "Palmitic acid", "InChIKey": "IPCSVZSSVZVIGE-UHFFFAOYSA-N",
     "PubChem ID": "985", "KEGG ID": "C00249"},
    {"Name": "ChEBI Only Compound", "ChEBI ID": "99999"},
    {"Name": "Lipid Only Compound", "LIPID MAPS ID": "LMGP01010999"},
    {"Name": "Unresolvable Nonsense"},
])
matched, meta = mm.match_compounds(df, CROSSWALK)
tiers = dict(zip(matched["Name"], matched["match_tier"]))

check(tiers.get("Palmitic acid") == "inchikey",
      f"the MOST specific available identifier wins - InChIKey, not the PubChem/KEGG also "
      f"present on the same row (got {tiers.get('Palmitic acid')!r})")
check(tiers.get("ChEBI Only Compound") == "chebi",
      f"a compound reachable ONLY by ChEBI ID still resolves - this is the Reactome case, "
      f"which has no other route in (got {tiers.get('ChEBI Only Compound')!r})")
check(tiers.get("Lipid Only Compound") == "lm_id",
      f"a lipid reachable ONLY by LM_ID still resolves (got {tiers.get('Lipid Only Compound')!r})")
check(meta["n_unmatched"] == 1 and meta["n_matched"] == 3,
      f"an unresolvable compound is counted, not silently dropped or crashed on "
      f"(matched={meta['n_matched']}, unmatched={meta['n_unmatched']})")
check(sum(meta["n_matched_by_tier"].values()) == meta["n_matched"],
      f"tier counts PARTITION the matched rows - they must sum to n_matched exactly, or "
      f"the summary double-counts (sum={sum(meta['n_matched_by_tier'].values())}, "
      f"n_matched={meta['n_matched']})")

# =============================================================================
print("=== InChIKey skeleton fallback: preferred LAST, and reported honestly ===")
# =============================================================================
stereo = pd.DataFrame([
    # Exact match available -> must NOT collapse onto the other stereoisomer
    {"Name": "Stereo A", "InChIKey": "SKELSHAREDXX-AAAAAAAAAA-N"},
    {"Name": "Stereo B", "InChIKey": "SKELSHAREDXX-BBBBBBBBBB-N"},
    # No exact match -> skeleton fallback, and it must be COUNTED
    {"Name": "Stereo C", "InChIKey": "SKELSHAREDXX-CCCCCCCCCC-N"},
])
smatched, smeta = mm.match_compounds(stereo, CROSSWALK)
keys = dict(zip(smatched["Name"], smatched["compound_key"]))
stiers = dict(zip(smatched["Name"], smatched["match_tier"]))
check(keys["Stereo A"] != keys["Stereo B"],
      "two genuine stereoisomers with exact keys do NOT collapse onto one compound - the "
      "skeleton index must never override an available exact match")
check(stiers["Stereo A"] == "inchikey" and stiers["Stereo B"] == "inchikey",
      f"both exact rows report the exact tier (got {stiers})")
check(stiers["Stereo C"] == "skeleton" and smeta["n_skeleton_only"] == 1,
      f"the row with no exact match falls back to the skeleton AND is counted separately, "
      f"so a user can see how much rests on stereochemistry-insensitive matching "
      f"(tier={stiers['Stereo C']!r}, n_skeleton_only={smeta['n_skeleton_only']})")

# =============================================================================
print("=== lipid shorthand: class yes, pathway no ===")
# =============================================================================
shorthand = pd.DataFrame([
    {"Name": "PC 34:1"},
    {"Name": "TG 52:2"},
    {"Name": "Palmitic acid", "InChIKey": "IPCSVZSSVZVIGE-UHFFFAOYSA-N"},
])
hmatched, hmeta = mm.match_compounds(shorthand, CROSSWALK)
by_name = {r["Name"]: r for _, r in hmatched.iterrows()}

check(bool(by_name["PC 34:1"]["class_level_only"]),
      "a shorthand name that DOES exist in RefMet is still class-level-only: RefMet gives "
      "bulk species a representative ChEBI ID, and letting that through would assert that "
      "one sum composition takes part in a specific reaction")
check(bool(by_name["TG 52:2"]["class_level_only"])
      and by_name["TG 52:2"]["lm_main_class"] == "Triradylglycerols",
      f"a shorthand name absent from every source still gets its class from the prefix "
      f"(got {by_name['TG 52:2']['lm_main_class']!r})")
check(by_name["PC 34:1"]["lm_main_class"] == "Glycerophosphocholines",
      f"a RefMet-matched shorthand row still gets LIPID MAPS fields filled from the prefix, "
      f"without which the Lipid Class analysis would be empty for it "
      f"(got {by_name['PC 34:1']['lm_main_class']!r})")
check(not bool(by_name["Palmitic acid"]["class_level_only"]),
      "a discrete compound is NOT flagged class-level-only")

pathway_keys = mm._pathway_capable_keys(hmatched)
check(len(pathway_keys) == 1,
      f"only the discrete compound is pathway-capable; both shorthand rows are excluded "
      f"(got {len(pathway_keys)} of 3)")
lipid_universe = mm.build_lipid_class_universe(hmatched, level="lm_main_class")
lipid_members = set()
for info in lipid_universe["terms"].values():
    lipid_members |= info["genes"]
check(len(lipid_members) == 3,
      f"...but ALL THREE reach the Lipid Class universe, which is the analysis their names "
      f"can actually support (got {len(lipid_members)})")

# A user who pins a structure themselves keeps pathway capability.
mask = pd.Series([True, False, False])
pmatched, _ = mm.match_compounds(shorthand, CROSSWALK, user_id_mask=mask)
check(not bool(pmatched.iloc[0]["class_level_only"]),
      "a shorthand row where the USER supplied a structure identifier is exempt - they have "
      "pinned one structure, unlike a retriever-guessed representative")

# =============================================================================
print("=== class universes, experiment scope: narrowed to measured compounds ===")
# =============================================================================
chem = mm.build_chemical_class_universe(matched, level="main_class")
chem_members = set()
for info in chem["terms"].values():
    chem_members |= info["genes"]
check(chem_members <= set(matched["compound_key"].astype(str)),
      "under experiment scope the class universe contains ONLY compounds from this experiment - "
      "the concentration tile's universe must stay the measured set")

# =============================================================================
print("=== class universes, database scope: the WHOLE crosswalk, not the user's rows ===")
# =============================================================================
one = pd.DataFrame([{"Name": "Palmitic acid", "InChIKey": "IPCSVZSSVZVIGE-UHFFFAOYSA-N"}])
one_matched, _ = mm.match_compounds(one, CROSSWALK)
chem_db = mm.build_chemical_class_universe(one_matched, CROSSWALK, level="main_class", universe_scope="database")
chem_exp = mm.build_chemical_class_universe(one_matched, CROSSWALK, level="main_class", universe_scope="experiment")
check(chem_db["universe_size"] == len(CROSSWALK["by_key"]) - sum(
          1 for e in CROSSWALK["by_key"].values() if not e.get("main_class")),
      f"database scope covers every crosswalk compound that has a RefMet class "
      f"(got {chem_db['universe_size']})")
check(chem_exp["universe_size"] == 1,
      f"the same one-compound list under experiment scope collapses to a 1-compound universe "
      f"(got {chem_exp['universe_size']}) - which is why list tiles must use database scope")
check(chem_db["universe_size"] > chem_exp["universe_size"], "database universe is strictly larger than the list itself")

lip_db = mm.build_lipid_class_universe(hmatched, CROSSWALK, level="lm_main_class", universe_scope="database")
glyco = lip_db["terms"].get("LIPIDCLASS:Glycerophosphocholines", {"genes": set()})["genes"]
pc_key = str(by_name["PC 34:1"]["compound_key"])
check(pc_key in glyco,
      "a lipid SHORTHAND compound (no LIPID MAPS structure of its own) is unioned into its own class "
      "under database scope - otherwise a shorthand-heavy list would score its best analysis as empty")
lmsd_members = {k for k, e in CROSSWALK["by_key"].items() if e.get("lm_main_class") == "Glycerophosphocholines"}
check(lmsd_members <= glyco, "...and the class still contains every real LIPID MAPS member of it")

# =============================================================================
print("=== _crosswalk_id_maps / _id_to_key_maps(scope='database') read the whole crosswalk ===")
# =============================================================================
by_chebi_all, by_kegg_all, by_skel_all = mm._crosswalk_id_maps(CROSSWALK)
check(by_kegg_all.get("C00249") == palmitic_key, "KEGG id -> key covers compounds the user never uploaded")
check(mm._id_to_key_maps(CROSSWALK, scope="database")[1] == by_kegg_all,
      "_id_to_key_maps(scope='database') is the crosswalk's own KEGG map")
_, kegg_exp = mm._id_to_key_maps(one_matched, scope="experiment")
check(set(kegg_exp) <= {"C00249"}, f"experiment scope only knows the user's own compounds (got {sorted(kegg_exp)})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
