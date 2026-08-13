"""Guards against a real class of drift risk found while auditing organism
generality across the app: several per-database organism allowlists are
DUPLICATED between Python (the source of truth - TRRUST_URL_BY_TAXON,
OMNIPATH_TAXA, CHEA3_TAXON in Transcriptomics_regulon.py; MSIGDB_TAXON_SUFFIX,
MSIGDB_COLLECTIONS in Transcriptomics_module.py) and JS (the UI mirror -
TRANSCRIPTOMICS_REGULON_TAXA/_regulonDbAvailableForTaxon,
_msigdbAvailableForTaxon/_msigdbCollectionAvailableForTaxon in app.js). If
the Python side ever gains a taxon (e.g. TRRUST publishes a rat file), the
JS mirror will NOT pick it up automatically - it has to be edited by hand -
and nothing today would say so. This test compares the two live against
each other instead of against a hardcoded literal, so a change on the
Python side that isn't mirrored in JS fails here immediately, not months
later when a user's checkbox is mysteriously greyed out for an organism
that actually has data.

(playwright_msigdb.py already has one such check - MSigDB's collection
KEYS/ORDER match exactly - but against a hardcoded expected list, not the
live Python constant, and it doesn't cover per-taxon availability at all.
This file is complementary, not a replacement.)

Needs a live server (to load app.js) - no external network calls.
"""
import os
import sys

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx
import Transcriptomics_regulon as tr

BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


PROBE_TAXA = ["9606", "10090", "10116", "3702"]  # human, mouse, rat, Arabidopsis


def regulon_py_available(db: str, taxon: str) -> bool:
    if db == "ChEA3":
        return taxon == tr.CHEA3_TAXON
    if db == "TRRUST":
        return taxon in tr.TRRUST_URL_BY_TAXON
    return taxon in tr.OMNIPATH_TAXA  # DoRothEA, CollecTRI


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.goto(BASE, wait_until="networkidle")

    print("=== Regulon: JS taxon/database mirror matches Python source of truth ===")
    trrust_taxa = set(tr.TRRUST_URL_BY_TAXON.keys())
    omnipath_taxa = set(tr.OMNIPATH_TAXA)
    check(trrust_taxa == omnipath_taxa,
          f"TRRUST_URL_BY_TAXON and OMNIPATH_TAXA cover the exact same taxa in Python - the JS mirror "
          f"assumes this (ONE shared taxon map covers TRRUST/DoRothEA/CollecTRI) and would silently be "
          f"wrong for whichever source diverges first (TRRUST={trrust_taxa}, OmniPath={omnipath_taxa})")

    js_regulon_taxa = set(page.evaluate("Object.keys(TRANSCRIPTOMICS_REGULON_TAXA)"))
    check(js_regulon_taxa == omnipath_taxa,
          f"TRANSCRIPTOMICS_REGULON_TAXA (app.js) mirrors OMNIPATH_TAXA / TRRUST_URL_BY_TAXON exactly "
          f"(js={sorted(js_regulon_taxa)}, python={sorted(omnipath_taxa)})")

    js_databases = page.evaluate("TRANSCRIPTOMICS_REGULON_DATABASES")
    check(js_databases == tr.REGULON_DATABASES,
          f"TRANSCRIPTOMICS_REGULON_DATABASES (app.js) matches REGULON_DATABASES exactly, same order "
          f"(js={js_databases}, python={tr.REGULON_DATABASES})")

    mismatches = []
    for db in tr.REGULON_DATABASES:
        for taxon in PROBE_TAXA:
            py_avail = regulon_py_available(db, taxon)
            js_avail = page.evaluate(f"_regulonDbAvailableForTaxon({db!r}, {taxon!r})")
            if js_avail != py_avail:
                mismatches.append((db, taxon, "python=" + str(py_avail), "js=" + str(js_avail)))
    check(len(mismatches) == 0,
          f"_regulonDbAvailableForTaxon (app.js) agrees with Python's ACTUAL per-source gating for every "
          f"database x probed taxon (incl. rat/10116 and Arabidopsis/3702, which Python correctly rejects "
          f"for all four regulon sources) - got mismatches: {mismatches}")

    print("=== MSigDB: JS taxon/collection mirror matches Python source of truth ===")
    js_msigdb_taxa = {t: page.evaluate(f"_msigdbAvailableForTaxon({t!r})") for t in PROBE_TAXA}
    py_msigdb_taxa = {t: (t in tx.MSIGDB_TAXON_SUFFIX) for t in PROBE_TAXA}
    check(js_msigdb_taxa == py_msigdb_taxa,
          f"_msigdbAvailableForTaxon (app.js) matches MSIGDB_TAXON_SUFFIX for every taxon probed "
          f"(js={js_msigdb_taxa}, python={py_msigdb_taxa})")

    collection_mismatches = []
    for key in tx.MSIGDB_COLLECTIONS:
        for taxon in ("9606", "10090"):
            py_avail = tx.msigdb_collection_available(key, taxon)
            js_avail = page.evaluate(f"""() => {{
                const c = TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.find(x => x.key === {key!r});
                return c ? _msigdbCollectionAvailableForTaxon(c, {taxon!r}) : null;
            }}""")
            if js_avail != py_avail:
                collection_mismatches.append((key, taxon, "python=" + str(py_avail), "js=" + str(js_avail)))
    check(len(collection_mismatches) == 0,
          f"_msigdbCollectionAvailableForTaxon (app.js) agrees with msigdb_collection_available for every "
          f"collection x (human, mouse) pair - got mismatches: {collection_mismatches}")

    print("=== SENSITIVITY: the exhaustive regulon cross-check is not vacuous ===")
    # The exact bug class this exists to catch: a gate that always says
    # "available", silently offering a database for an organism it has no
    # real data for.
    fake_mismatches = [
        (db, taxon) for db in tr.REGULON_DATABASES for taxon in PROBE_TAXA
        if True != regulon_py_available(db, taxon)  # noqa: E712 - "always available" stand-in
    ]
    check(len(fake_mismatches) > 0,
          f"a deliberately-broken 'always available' gate DOES produce mismatches against the real "
          f"Python source of truth (got {len(fake_mismatches)} of {len(tr.REGULON_DATABASES) * len(PROBE_TAXA)}) "
          f"- proving the cross-check above is not vacuous")

    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
