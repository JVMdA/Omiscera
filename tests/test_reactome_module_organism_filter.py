"""Regression test for a real bug found while auditing organism-generality
across the app: Genomics_network.enrich_modules (the per-module functional
enrichment step behind Network Propagation / Module Detection's Reactome
option) called Genomics_module.build_reactome_universe WITHOUT taxon_name.

build_reactome_universe filters by organism NAME, not taxon_id (unlike GO,
which filters by taxon_id and is unaffected). With taxon_name omitted, the
filter never runs, so the "universe" silently pools EVERY organism's
Reactome pathways - measured live against the real bulk file: 91,531
accessions across all species vs 12,155 for human alone, a >7x inflation
that deflates fold-enrichment and corrupts p-values. It hit non-model
organisms hardest (their true pathway set is small against the pooled
all-species one) and was cached under the SAME (database, taxon_id) key
every other caller (Genomics ORA, Burden Test) uses for that taxon - so a
single module-detection run for taxon X could poison Reactome results for
every other analysis on taxon X afterward.

THE FIX: enrich_modules now accepts taxon_name and threads it through to
build_reactome_universe (Genomics_network.py); the API.py call site
(_job_netprop_run) now passes it. The universe cache key was also bumped
(REACTOME_UNIVERSE_CACHE_VERSION) so any pre-fix poisoned cache entry is
never re-served, without needing a manual purge.

No server needed - direct function calls against the real, local Reactome
bulk file (already cached on disk; no network call is made unless that
cache has expired past REACTOME_REFRESH_DAYS).
"""
import os
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm
import Genomics_network as gn

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


TAXON = "9606"

print("=== organism filter genuinely restricts the universe (real Reactome data) ===")
# NOTE: build_reactome_universe's cache key is (database, taxon_id) only -
# NOT taxon_name - so calling it with taxon_name="" and force_refresh=True
# would overwrite the real, shared cache file for taxon 9606 with a
# genuinely wrong (all-species) universe. Guard with try/finally: whatever
# happens, the LAST write for this taxon must be the correct, human-only
# one, so any real analysis reading this cache afterward gets real data.
try:
    unnamed = gm.build_reactome_universe(TAXON, taxon_name="", force_refresh=True)
    named = gm.build_reactome_universe(TAXON, taxon_name="Homo sapiens", force_refresh=True)
finally:
    gm.build_reactome_universe(TAXON, taxon_name="Homo sapiens", force_refresh=True)

check(named["universe_size"] > 0, f"human-filtered universe is non-empty (got {named['universe_size']})")
check(named["universe_size"] < 20000,
      f"human-filtered universe stays in the real, measured human-only range "
      f"(got {named['universe_size']}, real measurement was 12,155)")
check(unnamed["universe_size"] > 3 * named["universe_size"],
      f"an EMPTY taxon_name pools every organism's pathways together, at least 3x larger than the "
      f"human-only universe for the exact same taxon_id (human={named['universe_size']}, "
      f"unfiltered={unnamed['universe_size']}, real measurement was a >7x inflation)")

print("=== enrich_modules threads taxon_name through to build_reactome_universe ===")
# Spy on the real function (no gene-set data needed - modules=[] means the
# per-module loop never runs, but the universe is still built beforehand,
# exercising the exact line this bug lived in).
captured = {}
real_fn = gm.build_reactome_universe


def _spy(taxon_id, taxon_name="", **kwargs):
    captured["taxon_id"] = taxon_id
    captured["taxon_name"] = taxon_name
    return {"pathways": {}, "universe_size": 0, "fetched_at": ""}


gn.genomics_ora.build_reactome_universe = _spy
try:
    gn.enrich_modules([], {"uniprot_to_protein": {}}, TAXON, taxon_name="Homo sapiens", database="Reactome")
    check(captured.get("taxon_name") == "Homo sapiens",
          f"enrich_modules passes ITS taxon_name argument through unchanged "
          f"(got {captured.get('taxon_name')!r})")
finally:
    gn.genomics_ora.build_reactome_universe = real_fn

print("=== SENSITIVITY: the wiring check above is not vacuous ===")
# Re-run with the EXACT pre-fix call signature at the API.py call site: no
# taxon_name argument passed at all (falls back to the parameter's own
# default of "", which is what silently caused the bug before this fix -
# the API.py caller simply never supplied the argument).
captured.clear()
gn.genomics_ora.build_reactome_universe = _spy
try:
    gn.enrich_modules([], {"uniprot_to_protein": {}}, TAXON, database="Reactome")
    check(captured.get("taxon_name") == "",
          f"omitting taxon_name (the pre-fix API.py call shape) genuinely reproduces the empty-name bug "
          f"(got {captured.get('taxon_name')!r}) - proving the check above is not vacuous")
finally:
    gn.genomics_ora.build_reactome_universe = real_fn

print("=== cache-version bump: a pre-existing poisoned entry under the OLD key is never re-served ===")
old_style_db = "Reactome"  # the pre-fix cache_db string, with no version suffix
old_style_path = gm._universe_cache_file_path(old_style_db, TAXON)
# This is the SAME file Genomics ORA/Burden legitimately use for GO's own
# "Reactome" cache_db elsewhere in the app under the pre-fix naming, so back
# up and restore it rather than just deleting - a crash mid-test must not
# leave a real cache entry gone or poisoned for the next real user.
had_backup = os.path.exists(old_style_path)
backup_path = old_style_path + ".test_backup"
if had_backup:
    os.replace(old_style_path, backup_path)

try:
    poisoned = {
        "pathways": {"FAKE_POISONED_PATHWAY": {"name": "poisoned", "genes": {"FAKE1", "FAKE2"}}},
        "universe_size": 999999,
        "fetched_at": "2000-01-01T00:00:00+00:00",
    }
    gm._save_universe_cache(old_style_db, TAXON, poisoned)
    reread_old = gm._load_universe_cache(old_style_db, TAXON)
    check(reread_old is not None and reread_old["universe_size"] == 999999,
          "sanity: the old-style cache file was actually written and is independently readable")

    fresh = gm.build_reactome_universe(TAXON, taxon_name="Homo sapiens", force_refresh=False)
    check(fresh["universe_size"] != 999999,
          f"build_reactome_universe ignores the old-key poisoned entry entirely - the version-bumped "
          f"cache key means it was never even looked at (got universe_size={fresh['universe_size']})")
    check("FAKE_POISONED_PATHWAY" not in fresh["pathways"],
          "the poisoned pathway does not leak into a real (post-fix) universe")
finally:
    try:
        os.remove(old_style_path)
    except OSError:
        pass
    if had_backup:
        os.replace(backup_path, old_style_path)

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
