"""Direct unit tests for the MSigDB fetcher added to Transcriptomics_module.py
(fetch_msigdb_sets, msigdb_collection_available, _msigdb_url, _parse_msigdb_gmt,
_msigdb_selected_keys). No live network for most of this file - the live-URL
check (do the pinned URLs still resolve, do the counts look right) lives
separately in tests/live_msigdb_urls.py, since that one genuinely needs the
internet and takes longer.

Mirrors test_regulon_engine.py's fake-HTTP-session idiom, and its cache-
isolation discipline: an earlier session in THIS project once had a test
write fake data straight into the real on-disk ORA universe cache, and the
live app silently served it for the rest of that session (30-day refresh
window). Every fetch below is redirected to a throwaway temp directory for
the whole file, restored at the end - never assume a "test" fetch is safe
just because it looks like a unit test.
"""
import os
import shutil
import sys
import tempfile

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx

_REAL_CACHE_DIR = tx.ORA_UNIVERSE_CACHE_DIR
_TEST_CACHE_DIR = tempfile.mkdtemp(prefix="omiscera_msigdb_test_cache_")
tx.ORA_UNIVERSE_CACHE_DIR = _TEST_CACHE_DIR

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


class _FakeResp:
    def __init__(self, text=None, status_code=200):
        self.status_code = status_code
        self.content = (text or "").encode("utf-8")


class _FakeSession:
    """Returns a fixed response (or raises) regardless of the URL - fine
    here since these tests only ever check ONE fetch_msigdb_sets call at a
    time, unlike the multi-collection integration tests in
    regression_ora_gsea.py which need real per-URL routing."""
    def __init__(self, text=None, status_code=200, raise_exc=None):
        self._text = text
        self._status_code = status_code
        self._raise_exc = raise_exc

    def get(self, *a, **k):
        if self._raise_exc is not None:
            raise self._raise_exc
        return _FakeResp(self._text, self._status_code)


# =============================================================================
# 1. MSIGDB_COLLECTIONS / MSIGDB_TAXON_SUFFIX - the availability map itself
# =============================================================================
print("=== MSIGDB_COLLECTIONS: shape and the explicit (never-derived) availability map ===")
check(len(tx.MSIGDB_COLLECTIONS) == 8, f"8 collections offered (got {len(tx.MSIGDB_COLLECTIONS)})")
for key, info in tx.MSIGDB_COLLECTIONS.items():
    check("label" in info and "human_file" in info and "mouse_file" in info,
          f"{key!r} has label/human_file/mouse_file keys")
    check(bool(info["human_file"]), f"{key!r} always has a human file (got {info['human_file']!r})")

print("=== msigdb_collection_available: known true/false cases ===")
check(tx.msigdb_collection_available("hallmark", "9606") is True, "hallmark available for human")
check(tx.msigdb_collection_available("hallmark", "10090") is True, "hallmark available for mouse")
check(tx.msigdb_collection_available("pid", "10090") is False, "PID has NO mouse file - correctly unavailable")
check(tx.msigdb_collection_available("oncogenic", "10090") is False, "Oncogenic has NO mouse file - correctly unavailable")
check(tx.msigdb_collection_available("cell_type", "10090") is True, "Cell type DOES have a mouse file (m8.all)")
check(tx.msigdb_collection_available("hallmark", "7955") is False, "no organism outside human/mouse is available (zebrafish)")
check(tx.msigdb_collection_available("bogus_key", "9606") is False, "an unknown collection key is never available")

print("=== _msigdb_url: exact construction, matches the live-verified pattern ===")
check(tx._msigdb_url("hallmark", "9606") ==
      f"{tx.MSIGDB_BASE_URL}/{tx.MSIGDB_VERSION}.Hs/h.all.v{tx.MSIGDB_VERSION}.Hs.symbols.gmt",
      f"human hallmark URL matches the pinned-version pattern (got {tx._msigdb_url('hallmark', '9606')})")
check(tx._msigdb_url("hallmark", "10090") ==
      f"{tx.MSIGDB_BASE_URL}/{tx.MSIGDB_VERSION}.Mm/mh.all.v{tx.MSIGDB_VERSION}.Mm.symbols.gmt",
      f"mouse hallmark URL uses the 'mh.all' filename, not 'h.all' (got {tx._msigdb_url('hallmark', '10090')})")
check(tx._msigdb_url("pid", "10090") is None, "no URL is ever constructed for an unavailable collection/organism pair")
check(tx._msigdb_url("bogus_key", "9606") is None, "no URL for an unknown collection key")
check(tx._msigdb_url("hallmark", "7955") is None, "no URL for an organism MSigDB doesn't cover")

print("=== _msigdb_selected_keys: canonical order, silently drops unknown keys ===")
check(tx._msigdb_selected_keys(["cell_type", "hallmark", "bogus"]) == ["hallmark", "cell_type"],
      f"returns in MSIGDB_COLLECTIONS' own canonical order, regardless of input order, dropping unknowns "
      f"(got {tx._msigdb_selected_keys(['cell_type', 'hallmark', 'bogus'])})")
check(tx._msigdb_selected_keys(None) == [], "None is treated as no selection, not an error")
check(tx._msigdb_selected_keys([]) == [], "empty selection stays empty")


# =============================================================================
# 2. _parse_msigdb_gmt - GMT parsing correctness on a hand-built fixture
# =============================================================================
print("=== _parse_msigdb_gmt: parsing ===")
fixture_gmt = (
    "SET_A\thttps://msigdb.example/SET_A\tgene1\tGENE2\t gene3 \n"
    "\n"
    "SET_B\thttps://msigdb.example/SET_B\tgene4\n"
    "TOO_SHORT\tonly_one_field\n"
    "SET_C\thttps://msigdb.example/SET_C\t\t\n"
)
parsed = tx._parse_msigdb_gmt(fixture_gmt)
check(set(parsed.keys()) == {"SET_A", "SET_B"},
      f"only well-formed, non-empty-gene lines become sets (got {sorted(parsed.keys())})")
check(parsed["SET_A"]["genes"] == ["GENE1", "GENE2", "GENE3"],
      f"genes are upper-cased, whitespace-trimmed, and de-duplicated/sorted (got {parsed['SET_A']['genes']})")
check(parsed["SET_A"]["name"] == "SET_A", "set 'name' matches its own id")
check("TOO_SHORT" not in parsed, "a line with fewer than 3 tab-separated fields is rejected, not crashed on")
check("SET_C" not in parsed, "a line whose gene fields are all blank produces no set (empty genes list)")

dup_gmt = "SET_D\turl\tgeneX\tGENEX\tgenex\n"
check(tx._parse_msigdb_gmt(dup_gmt)["SET_D"]["genes"] == ["GENEX"],
      "case-insensitive duplicate gene symbols collapse to one entry")


# =============================================================================
# 3. fetch_msigdb_sets - the failure contract (never raises, never poisons
#    the cache on a bad fetch) + the success/caching path
# =============================================================================
print("=== fetch_msigdb_sets: HTTP failure -> fetch_failed=True, cache untouched ===")
old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(status_code=503)
try:
    r_http_fail = tx.fetch_msigdb_sets("oncogenic", "9606", force_refresh=True)
finally:
    tx.HTTP_SESSION = old_session
check(r_http_fail["fetch_failed"] is True, "a non-200 response sets fetch_failed=True")
check(r_http_fail["sets"] == {}, "and returns an empty sets dict, not a partial/garbage one")
cache_path_after_fail = tx._universe_cache_path(f"msigdb_oncogenic_9606_{tx.MSIGDB_VERSION}")
check(not os.path.exists(cache_path_after_fail), "a failed fetch never writes to the on-disk cache")

print("=== fetch_msigdb_sets: exception during fetch -> fetch_failed=True, no crash ===")
old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(raise_exc=ConnectionError("simulated network failure"))
try:
    r_exc = tx.fetch_msigdb_sets("oncogenic", "9606", force_refresh=True)
finally:
    tx.HTTP_SESSION = old_session
check(r_exc["fetch_failed"] is True, "a raised exception is caught and reported as fetch_failed, not propagated")

print("=== fetch_msigdb_sets: unavailable collection/organism -> degrades safely, does not fetch ===")
old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(text="SHOULD_NOT_BE_PARSED\turl\tG1\n", status_code=200)
try:
    r_unavailable = tx.fetch_msigdb_sets("pid", "10090", force_refresh=True)
finally:
    tx.HTTP_SESSION = old_session
check(r_unavailable["fetch_failed"] is True, "calling fetch_msigdb_sets for an unavailable pair still fails safely")
check(r_unavailable["sets"] == {}, "and never returns sets even if a (misconfigured) session would have responded")

print("=== fetch_msigdb_sets: success path, caching, and force_refresh ===")
fake_gmt = "AKT_SIGNALING\turl\tAKT1\tAKT2\tAKT3\n" "MTOR_SIGNALING\turl\tMTOR\tRPTOR\n"
old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(text=fake_gmt, status_code=200)
try:
    r_ok = tx.fetch_msigdb_sets("oncogenic", "9606", force_refresh=True)
finally:
    tx.HTTP_SESSION = old_session
check(r_ok["fetch_failed"] is False, "a genuine 200 with parseable content is NOT marked failed")
check(set(r_ok["sets"].keys()) == {"AKT_SIGNALING", "MTOR_SIGNALING"}, "both fixture sets are parsed")

cache_path_ok = tx._universe_cache_path(f"msigdb_oncogenic_9606_{tx.MSIGDB_VERSION}")
check(os.path.exists(cache_path_ok), "a genuinely successful, non-empty fetch DOES write to the cache")

# Now simulate a session outage - if caching worked, this must NOT reach the
# (broken) session at all.
old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(status_code=500)
try:
    r_cached = tx.fetch_msigdb_sets("oncogenic", "9606", force_refresh=False)
finally:
    tx.HTTP_SESSION = old_session
check(r_cached["fetch_failed"] is False and set(r_cached["sets"].keys()) == {"AKT_SIGNALING", "MTOR_SIGNALING"},
      "a subsequent call with force_refresh=False serves the cached result without touching HTTP_SESSION again")

old_session = tx.HTTP_SESSION
tx.HTTP_SESSION = _FakeSession(status_code=500)
try:
    r_forced = tx.fetch_msigdb_sets("oncogenic", "9606", force_refresh=True)
finally:
    tx.HTTP_SESSION = old_session
check(r_forced["fetch_failed"] is True, "force_refresh=True bypasses the cache and genuinely re-fetches")

print("=== cache key includes MSIGDB_VERSION - a version bump self-invalidates ===")
old_version = tx.MSIGDB_VERSION
tx.MSIGDB_VERSION = "1999.9"
try:
    cache_path_new_version = tx._universe_cache_path(f"msigdb_oncogenic_9606_{tx.MSIGDB_VERSION}")
finally:
    tx.MSIGDB_VERSION = old_version
check(cache_path_new_version != cache_path_ok,
      "a different MSIGDB_VERSION produces a DIFFERENT cache path, so a version bump can never accidentally "
      "serve a stale-version cache entry")


# =============================================================================
# 4. Sensitivity: deliberately break msigdb_collection_available's mapping,
#    confirm the availability tests above would have caught it
# =============================================================================
print("=== sensitivity: deliberately break the mouse-availability map, confirm it's caught ===")
_source_before = None
with open(tx.__file__, "rb") as fh:
    _source_before = fh.read()

_original_available = tx.msigdb_collection_available


def _broken_available(collection_key, taxon_id):
    # DELIBERATELY WRONG: claims every collection is available for mouse,
    # the exact bug this module's design note warns against (a c->m
    # substitution assumption instead of the explicit per-collection map).
    info = tx.MSIGDB_COLLECTIONS.get(collection_key)
    if info is None:
        return False
    return str(taxon_id).strip() in ("9606", "10090")


tx.msigdb_collection_available = _broken_available
broken_pid_mouse = tx.msigdb_collection_available("pid", "10090")
broken_oncogenic_mouse = tx.msigdb_collection_available("oncogenic", "10090")
tx.msigdb_collection_available = _original_available

with open(tx.__file__, "rb") as fh:
    _source_after = fh.read()
check(_source_after == _source_before,
      "the module's source file on disk was never touched by this sensitivity test (monkeypatch only, in-memory)")
check(broken_pid_mouse is True and broken_oncogenic_mouse is True,
      "the broken version WOULD wrongly claim PID/Oncogenic are available for mouse - proving the real "
      "availability tests above are not vacuous (they assert False for exactly these two)")


tx.ORA_UNIVERSE_CACHE_DIR = _REAL_CACHE_DIR
shutil.rmtree(_TEST_CACHE_DIR, ignore_errors=True)

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
