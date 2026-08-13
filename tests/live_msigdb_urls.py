"""LIVE network test: does every pinned MSigDB URL in MSIGDB_COLLECTIONS
still resolve, for every organism it claims to support, and does the
resulting gene-set count look like MSigDB's own published figures (or, for
the couple of collections with no single published number, this session's
own directly-measured baseline)? This is the check that turns a rotted
version pin into a loud, specific failure instead of the kind of silent
2020-vintage drift the OLD Hallmark-via-Enrichr fetcher had for years.

Needs internet. Makes ~14 real downloads (one per available collection x
organism pair) - the two biggest (Perturbations/CGP, Immunologic) are a few
MB each but still fast, since these are static files on Broad's own server,
not a live query API.
"""
import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# (expected_count, tolerance_fraction) - tight where MSigDB's own collections
# pages give an exact published number (verified live before this session's
# plan was written), looser where only this suite's OWN direct fetch is the
# baseline (mouse Phenotype/MPT and mouse WikiPathways are not broken out on
# MSigDB's summary pages the way the human ones are - 92 and 202 are what
# were actually measured here, not numbers copied from a webpage).
#
# Every pair now has a baseline: an earlier version of this file left mouse
# WikiPathways at None ("no baseline on file"), which meant a silent collapse
# to, say, 3 sets would still have passed. The README's own coverage table is
# generated from these same numbers, so a drift here is also a docs drift.
EXPECTED = {
    ("hallmark", "9606"): (50, 0.05),
    ("hallmark", "10090"): (50, 0.05),
    ("perturbations", "9606"): (3555, 0.15),
    ("perturbations", "10090"): (984, 0.15),
    ("wikipathways", "9606"): (925, 0.15),
    ("wikipathways", "10090"): (202, 0.25),
    ("pid", "9606"): (196, 0.15),
    ("phenotype", "9606"): (5793, 0.15),
    ("phenotype", "10090"): (92, 0.25),
    ("oncogenic", "9606"): (189, 0.05),
    ("immunologic", "9606"): (5219, 0.15),
    ("immunologic", "10090"): (787, 0.15),
    ("cell_type", "9606"): (866, 0.15),
    ("cell_type", "10090"): (233, 0.15),
}

print(f"=== live fetch: every (collection, organism) pair MSIGDB_COLLECTIONS claims is available (v{tx.MSIGDB_VERSION}) ===")
n_checked = 0
for collection_key, info in tx.MSIGDB_COLLECTIONS.items():
    for taxon_id in ("9606", "10090"):
        if not tx.msigdb_collection_available(collection_key, taxon_id):
            continue
        n_checked += 1
        result = tx.fetch_msigdb_sets(collection_key, taxon_id, verbose=False, force_refresh=True)
        label = f"{collection_key}/{taxon_id}"
        check(result["fetch_failed"] is False, f"{label}: pinned URL resolves (fetch_failed=False)")
        n_sets = len(result["sets"])
        check(n_sets > 0, f"{label}: at least one gene set was parsed (got {n_sets})")

        expected, tol = EXPECTED.get((collection_key, taxon_id), (None, None))
        if expected is not None:
            lo, hi = expected * (1 - tol), expected * (1 + tol)
            check(lo <= n_sets <= hi,
                  f"{label}: set count {n_sets} is within {tol:.0%} of the expected ~{expected} "
                  f"(range [{lo:.0f}, {hi:.0f}])")
        else:
            print(f"    ({label}: no baseline count on file, got {n_sets} sets - informational only)")

check(n_checked == 14, f"exactly 14 available (collection, organism) pairs were exercised (got {n_checked}) - "
                        f"if this changes, MSIGDB_COLLECTIONS' availability map changed and EXPECTED above needs "
                        f"a matching update, not just a passing re-run")


# =============================================================================
# Upstream-coverage watch: this is the mechanism that makes "if MSigDB adds
# another organism, will the app catch it?" actually true, rather than just
# hoped for. MSIGDB_TAXON_SUFFIX (Transcriptomics_module.py) hardcodes
# {"9606": "Hs", "10090": "Mm"} because that is genuinely all MSigDB has ever
# published (verified live: every release directory 2022.1 through the
# current pin lists only .Hs/.Mm). If Broad ever adds a third organism
# suffix, this check goes red - not silently invisible - which is the
# trigger to widen MSIGDB_TAXON_SUFFIX rather than something nobody notices.
# =============================================================================
print(f"=== upstream-coverage watch: MSigDB v{tx.MSIGDB_VERSION}'s own release listing still covers only Hs/Mm ===")
try:
    resp = tx.HTTP_SESSION.get(tx.MSIGDB_BASE_URL + "/", headers=tx.OMISCERA_HEADERS, timeout=60)
    check(resp.status_code == 200, f"MSigDB release index page is reachable (got {resp.status_code})")
    if resp.status_code == 200:
        dirs = set(re.findall(r'href="' + re.escape(tx.MSIGDB_VERSION) + r'\.([A-Za-z]+)/"', resp.text))
        check(len(dirs) > 0,
              f"the current pinned version ({tx.MSIGDB_VERSION}) appears in the release index at all "
              f"(found organism suffixes: {sorted(dirs)}) - if this is empty, the index page's HTML "
              f"structure changed and this check needs updating, not the organism gate")
        check(dirs == {"Hs", "Mm"},
              f"MSigDB v{tx.MSIGDB_VERSION} still publishes exactly {{Hs, Mm}} and nothing else "
              f"(got {sorted(dirs)}) - if this ever grows, MSIGDB_TAXON_SUFFIX in Transcriptomics_module.py "
              f"is now artificially narrower than upstream and should be widened")
except Exception as exc:
    check(False, f"upstream-coverage watch request failed: {exc}")

# SENSITIVITY for the watch above: it depends on a live server's real
# content, which this suite doesn't control, so prove the PARSING/COMPARISON
# logic itself would catch a third organism by feeding it synthetic HTML
# containing one, using the exact same regex and equality check.
synthetic_html = (
    f'<a href="{tx.MSIGDB_VERSION}.Hs/">{tx.MSIGDB_VERSION}.Hs/</a>\n'
    f'<a href="{tx.MSIGDB_VERSION}.Mm/">{tx.MSIGDB_VERSION}.Mm/</a>\n'
    f'<a href="{tx.MSIGDB_VERSION}.Zz/">{tx.MSIGDB_VERSION}.Zz/</a>\n'  # hypothetical 3rd organism
)
synthetic_dirs = set(re.findall(r'href="' + re.escape(tx.MSIGDB_VERSION) + r'\.([A-Za-z]+)/"', synthetic_html))
check(synthetic_dirs == {"Hs", "Mm", "Zz"},
      f"sanity: the regex correctly extracts all three organism suffixes from synthetic HTML (got {synthetic_dirs})")
check(synthetic_dirs != {"Hs", "Mm"},
      f"a synthetic release listing WITH a third organism genuinely fails the {{Hs, Mm}} equality check "
      f"(got {synthetic_dirs}) - proving the real upstream-coverage watch above is not vacuous")


# =============================================================================
# Sensitivity: deliberately break the version pin, confirm the live check
# above would actually go red - a test that can't fail is not evidence.
# =============================================================================
print("=== sensitivity: deliberately break MSIGDB_VERSION, confirm the live fetch fails loudly ===")
_source_before = None
with open(tx.__file__, "rb") as fh:
    _source_before = fh.read()

old_version = tx.MSIGDB_VERSION
tx.MSIGDB_VERSION = "1900.1"  # a version that has certainly never existed on the release server
try:
    broken_result = tx.fetch_msigdb_sets("hallmark", "9606", verbose=False, force_refresh=True)
finally:
    tx.MSIGDB_VERSION = old_version

with open(tx.__file__, "rb") as fh:
    _source_after = fh.read()
check(_source_after == _source_before,
      "the module's source file on disk was never touched by this sensitivity test (monkeypatch only, in-memory)")
check(broken_result["fetch_failed"] is True,
      "a bogus/rotted version pin makes the fetch fail (404) - proving the live check above would "
      "actually catch a real rotted pin, not just pass unconditionally")

# Confirm the real pin still works immediately after, i.e. this test didn't
# leave any cached bad state behind for MSIGDB_VERSION's real value.
restored_result = tx.fetch_msigdb_sets("hallmark", "9606", verbose=False, force_refresh=True)
check(restored_result["fetch_failed"] is False,
      "the real, current MSIGDB_VERSION still resolves immediately after the sensitivity check")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
