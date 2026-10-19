# -*- coding: utf-8 -*-
"""Unit tests for Metabolomics_retriever.py, offline.

This file is a DELIBERATE fork of test_compound_resolver.py (the resolver
actually wired into the separate, existing "Metabolic Connection" feature) -
kept permanently independent so this new module can adapt it freely without
risking that feature. See its own module docstring for the full rationale.

Most of the resolution logic already carries its own offline, mocked-network
test functions inline (test_user_provided_identifier_preservation, etc.) -
this file calls those directly (they raise AssertionError on failure) rather
than re-deriving equivalent coverage, and adds tests for what is UNIQUE to
this file as a fork: cache-path isolation from test_compound_resolver.py, and
a real crash this fork fixed (an empty cache CSV raising EmptyDataError and
killing an entire resolve run - found live 2026-09-21, the empty-cache state
is legitimate, not an error).
"""
import os
import sys

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Metabolomics_retriever as mr  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def run_inline(fn, label):
    try:
        fn(verbose=False)
        check(True, label)
    except AssertionError as exc:
        check(False, f"{label} - assertion failed: {exc}")
    except Exception as exc:
        check(False, f"{label} - raised {type(exc).__name__}: {exc}")


# =============================================================================
print("=== cache isolation from test_compound_resolver.py - the whole point of the fork ===")
# =============================================================================
try:
    import test_compound_resolver as tcr
    check(mr.CACHE_DIR != tcr.CACHE_DIR,
          f"CACHE_DIR differs from the Metabolic Connection resolver's own "
          f"(mr={mr.CACHE_DIR!r}, tcr={tcr.CACHE_DIR!r})")
    for name in ("PUBCHEM_CACHE_PATH", "KEGG_CACHE_PATH", "PUBCHEM_SYNONYMS_CACHE_PATH",
                "HMDB_CACHE_PATH", "FINAL_RESULT_CACHE_PATH"):
        mr_path, tcr_path = getattr(mr, name), getattr(tcr, name)
        check(mr_path != tcr_path,
              f"{name} points at a different file than the live Metabolic Connection "
              f"resolver's own cache, so the two can never silently share (or corrupt) "
              f"on-disk state despite having no shared imports "
              f"(mr={mr_path!r}, tcr={tcr_path!r})")
except ImportError:
    check(False, "test_compound_resolver.py could not be imported for comparison")

# =============================================================================
print("=== _read_cache_csv: an empty cache is a legitimate state, not a crash ===")
# =============================================================================
_probe_dir = mr.CACHE_DIR
os.makedirs(_probe_dir, exist_ok=True)
_empty_path = os.path.join(_probe_dir, "_test_empty_cache_probe.csv")
with open(_empty_path, "w", encoding="utf-8") as fh:
    pass  # 0 bytes - the exact shape that previously raised pandas.errors.EmptyDataError
try:
    result = mr._read_cache_csv(_empty_path)
    check(result is None,
          f"a 0-byte cache file returns None (a legitimate 'nothing cached yet' state), "
          f"not a raised exception that would kill the whole resolve run (got {result!r})")
finally:
    os.remove(_empty_path)

_missing_path = os.path.join(_probe_dir, "_test_nonexistent_cache_probe.csv")
check(mr._read_cache_csv(_missing_path) is None,
      "a cache file that does not exist at all also returns None, not an exception")

_corrupt_path = os.path.join(_probe_dir, "_test_corrupt_cache_probe.csv")
with open(_corrupt_path, "w", encoding="utf-8") as fh:
    fh.write("this is not,valid\ncsv\"content\x00broken")
try:
    result = mr._read_cache_csv(_corrupt_path)
    check(result is None or hasattr(result, "columns"),
          "a corrupt cache file never crashes the caller, even if pandas can partially "
          f"parse it (got {type(result).__name__})")
finally:
    os.remove(_corrupt_path)

check(mr.ensure_cache_loaded() is None,
      "ensure_cache_loaded() runs cleanly against whatever real cache state exists on "
      "disk right now (the fix's real-world entry point)")

# =============================================================================
print("=== resolve_compound_table contract (offline, mocked network) ===")
# =============================================================================
run_inline(mr.test_user_provided_identifier_preservation,
          "user-supplied PubChem ID/SMILES/InChIKey/HMDB ID/Formula survive resolution "
          "unmodified, and OUTPUT_COLUMNS shape holds")
run_inline(mr.test_final_cache_preserves_user_identifiers,
          "the final-result cache never overwrites a user-supplied identifier, including "
          "the not_found/exact-cache-hit edge case")
run_inline(mr.test_hmdb_pubchem_enrichment,
          "HMDB IDs can be enriched from PubChem synonym data")
run_inline(mr.test_hmdb_normalization,
          "HMDB ID normalization/deduplication (HMDB0000162 / HMDB00162 / hmdb162 all "
          "collapse to one canonical form)")

# =============================================================================
print("=== OUTPUT_COLUMNS is the exact contract the new module's join logic depends on ===")
# =============================================================================
check(list(mr.OUTPUT_COLUMNS) == ["Name", "KEGG ID", "PubChem ID", "SMILES", "InChIKey",
                                  "HMDB ID", "Formula", "Mapping_status"],
      f"OUTPUT_COLUMNS carries exactly the identifiers Metabolomics_module.match_compounds "
      f"joins on (PubChem/InChIKey/KEGG/HMDB) - a change here would silently break that "
      f"join (got {list(mr.OUTPUT_COLUMNS)})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
