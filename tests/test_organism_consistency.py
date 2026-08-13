"""Direct unit tests for the organism-consistency fix in API.py: no live
server, no network - imports API.py directly and drives it with hand-built
session state and direct in-process calls, the same way test_example_files.py
and audit_regulon_untested.py already do.

WHY THIS FIX EXISTS. Two pre-existing bugs, found by a post-implementation
audit of the MSigDB work, that predate MSigDB and are more dangerous than
anything it introduced:

  (a) taxon_name used to be committed to the session SYNCHRONOUSLY in the
      resolve endpoint, while taxon_id was only committed inside the resolve
      JOB, on success. run_job does no session rollback on failure. So: a
      successful human resolve, then a re-resolve as mouse whose JOB fails,
      left the session holding the NEW taxon_name beside the OLD taxon_id and
      the OLD resolved table.
  (b) Changing the organism selector cleared nothing downstream - the
      resolved table, ORA/GSEA/regulon configs and results all stayed on
      screen describing the PREVIOUS organism.

Both are dangerous specifically because they fail SILENTLY: taxon_id drives
GO/QuickGO, MSigDB and the BioMart gene-length correction, while taxon_name
drives Reactome - so a drifted session doesn't error, it returns one result
set quietly mixing two organisms, for a user who cannot independently referee
the biology.

The fix has three layers, each tested separately below:
  1. _job_resolve_de now commits taxon_id AND taxon_name together, only on
     success (Part 1) - tested by forcing a failure and checking NOTHING
     partially committed.
  2. _require_consistent_taxon (Part 2) is the server-side invariant that
     makes the bug class unreachable regardless of cause - tested directly,
     and tested as WIRED INTO each of the three configure endpoints
     separately, since "the function exists" and "the function is actually
     called by all three routes" are different claims.
  3. The frontend clears downstream state on both an organism change and a
     new upload (Parts 3-4) - covered in tests/playwright_msigdb.py, which
     can drive the real UI; not duplicated here.
"""
import os
import sys

import pandas as pd

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import API

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
# 1. _require_consistent_taxon, unit-tested directly
# =============================================================================
print("=== _require_consistent_taxon: direct cases ===")


def _session_with(meta_taxon, current_taxon):
    tx = API._fresh_transcriptomics_state()
    if meta_taxon is not None:
        tx["meta"] = {"taxon": meta_taxon}
    tx["taxon_id"] = current_taxon
    return {"transcriptomics": tx}


# Agreeing taxa -> no exception.
try:
    API._require_consistent_taxon(_session_with("9606", "9606"))
    check(True, "agreeing taxa (9606 == 9606) do not raise")
except Exception as exc:
    check(False, f"agreeing taxa unexpectedly raised: {exc}")

# Disagreeing taxa -> HTTPException(400), message names BOTH values.
try:
    API._require_consistent_taxon(_session_with("9606", "10090"))
    check(False, "disagreeing taxa (9606 vs 10090) should have raised")
except API.HTTPException as exc:
    check(exc.status_code == 400, f"disagreeing taxa raise HTTP 400 (got {exc.status_code})")
    check("9606" in exc.detail and "10090" in exc.detail,
          f"error message names BOTH the resolved-for and current taxon (got: {exc.detail!r})")

# Never resolved (no meta yet) -> no exception - a DIFFERENT, pre-existing
# check ("resolved table required") handles this case with a better message;
# duplicating it here would produce a worse, confusing one.
try:
    API._require_consistent_taxon(_session_with(None, "9606"))
    check(True, "no meta yet (never resolved) does not raise - handled elsewhere")
except Exception as exc:
    check(False, f"no-meta case unexpectedly raised: {exc}")

# Resolved, but current taxon_id somehow blank -> no exception (defensive;
# should not be reachable via the UI, but must degrade safely, not crash).
try:
    API._require_consistent_taxon(_session_with("9606", ""))
    check(True, "blank current taxon_id does not raise (defensive - not a genuine mismatch to report)")
except Exception as exc:
    check(False, f"blank-taxon_id case unexpectedly raised: {exc}")

# Empty session / missing "transcriptomics" key entirely -> no crash.
try:
    API._require_consistent_taxon({})
    check(True, "a session with no 'transcriptomics' key at all does not crash")
except Exception as exc:
    check(False, f"empty session unexpectedly raised: {exc}")


# =============================================================================
# 2. Atomic commit: _job_resolve_de leaves NOTHING partially updated on failure
# =============================================================================
print("=== _job_resolve_de: a failed resolve commits nothing (Part 1) ===")

_FIXTURE_DF_HUMAN = pd.DataFrame({"gene_symbol": ["TP53"], "log2FC": [1.0], "padj": [0.01]})
_FIXTURE_META_HUMAN = {"taxon": "9606", "n_mapped": 1, "n_unmapped": 0}


def _fake_resolve_de_table(de_df, taxon_id, identifier_type=None, progress_callback=None, verbose=False):
    if taxon_id == "9606":
        return _FIXTURE_DF_HUMAN, _FIXTURE_META_HUMAN
    raise RuntimeError("simulated resolve failure (e.g. a real mygene.info outage)")


old_resolve_fn = API.transcriptomics_retriever.resolve_de_table
API.transcriptomics_retriever.resolve_de_table = _fake_resolve_de_table
try:
    session = {"transcriptomics": API._fresh_transcriptomics_state()}
    session["transcriptomics"]["de_df"] = pd.DataFrame({"x": [1]})

    # First call: succeeds, as human.
    API._job_resolve_de("fake-job-1", session, "gene_symbol", "9606", "Homo sapiens")
    check(session["transcriptomics"]["taxon_id"] == "9606", "first (successful) resolve commits taxon_id")
    check(session["transcriptomics"]["taxon_name"] == "Homo sapiens", "and taxon_name")
    check(session["transcriptomics"]["meta"]["taxon"] == "9606", "and meta")

    # Second call: fails, as mouse. Wrapped the way run_job wraps every job
    # function (API.py:619-625) - _job_resolve_de itself does not catch.
    try:
        API._job_resolve_de("fake-job-2", session, "gene_symbol", "10090", "Mus musculus")
        check(False, "the second call was supposed to raise (fixture simulates a resolve failure)")
    except RuntimeError:
        pass

    check(session["transcriptomics"]["taxon_id"] == "9606",
          f"AFTER the failed mouse resolve, taxon_id is still the original human value "
          f"(got {session['transcriptomics']['taxon_id']!r})")
    check(session["transcriptomics"]["taxon_name"] == "Homo sapiens",
          f"AFTER the failed mouse resolve, taxon_name is still 'Homo sapiens' - the bug this "
          f"guards was exactly a NEW taxon_name surviving beside an OLD taxon_id "
          f"(got {session['transcriptomics']['taxon_name']!r})")
    check(session["transcriptomics"]["resolved_df"] is _FIXTURE_DF_HUMAN,
          "AFTER the failed mouse resolve, resolved_df is still the original human DataFrame")
    check(session["transcriptomics"]["meta"] is _FIXTURE_META_HUMAN,
          "AFTER the failed mouse resolve, meta is still the original human meta")
finally:
    API.transcriptomics_retriever.resolve_de_table = old_resolve_fn


print("=== sensitivity: the atomic-commit test is not vacuous ===")
# Re-run the same scenario against a reimplementation of the PRE-FIX ordering
# (taxon_name committed unconditionally, before the resolve call that can
# fail) and confirm it DOES exhibit the drift the fixed version does not.
_source_before = None
with open(API.__file__, "rb") as fh:
    _source_before = fh.read()


def _broken_job_resolve_de(job_id, session, identifier_type, taxon_id, taxon_name):
    # DELIBERATELY BROKEN: mirrors the old bug - taxon_name committed first,
    # unconditionally, before the call that can raise.
    session["transcriptomics"]["taxon_name"] = taxon_name
    resolved_df, meta = API.transcriptomics_retriever.resolve_de_table(
        session["transcriptomics"]["de_df"], taxon_id, identifier_type=identifier_type,
    )
    session["transcriptomics"]["taxon_id"] = taxon_id
    session["transcriptomics"]["resolved_df"] = resolved_df
    session["transcriptomics"]["meta"] = meta


API.transcriptomics_retriever.resolve_de_table = _fake_resolve_de_table
try:
    broken_session = {"transcriptomics": API._fresh_transcriptomics_state()}
    broken_session["transcriptomics"]["de_df"] = pd.DataFrame({"x": [1]})
    _broken_job_resolve_de("fake-job-3", broken_session, "gene_symbol", "9606", "Homo sapiens")
    try:
        _broken_job_resolve_de("fake-job-4", broken_session, "gene_symbol", "10090", "Mus musculus")
    except RuntimeError:
        pass
    broken_drifted = (
        broken_session["transcriptomics"]["taxon_name"] == "Mus musculus"
        and broken_session["transcriptomics"]["taxon_id"] == "9606"
    )
finally:
    API.transcriptomics_retriever.resolve_de_table = old_resolve_fn

with open(API.__file__, "rb") as fh:
    _source_after = fh.read()
check(_source_after == _source_before,
      "the module's source file on disk was never touched by this sensitivity test (monkeypatch only, in-memory)")
check(broken_drifted,
      f"the PRE-FIX ordering really does leave taxon_name='Mus musculus' beside taxon_id='9606' after a failed "
      f"resolve - proving the atomic-commit test above is not vacuous "
      f"(got taxon_name={broken_session['transcriptomics']['taxon_name']!r}, "
      f"taxon_id={broken_session['transcriptomics']['taxon_id']!r})")


# =============================================================================
# 3. Every configure endpoint refuses a drifted session (wiring, not just
#    the guard function in isolation)
# =============================================================================
print("=== _require_consistent_taxon is actually wired into all three configure endpoints ===")


def _drifted_session(session_id):
    session = API.get_session(session_id)
    session["transcriptomics"] = API._fresh_transcriptomics_state()
    session["transcriptomics"]["resolved_df"] = pd.DataFrame({"gene_symbol": ["A1BG"]})
    session["transcriptomics"]["meta"] = {"taxon": "9606"}
    session["transcriptomics"]["taxon_id"] = "10090"  # deliberately drifted
    return session_id


def _consistent_session(session_id):
    session = API.get_session(session_id)
    session["transcriptomics"] = API._fresh_transcriptomics_state()
    session["transcriptomics"]["resolved_df"] = pd.DataFrame({"gene_symbol": ["A1BG"]})
    session["transcriptomics"]["meta"] = {"taxon": "9606"}
    session["transcriptomics"]["taxon_id"] = "9606"
    return session_id


ora_sid = _drifted_session("test-organism-consistency-ora")
try:
    API.transcriptomics_ora_configure({"databases": ["GO"]}, ora_sid)
    check(False, "ORA configure should have refused a drifted session")
except API.HTTPException as exc:
    check(exc.status_code == 400 and "Organism mismatch" in exc.detail,
          f"ORA /configure refuses a drifted session (got {exc.status_code}: {exc.detail!r})")

gsea_sid = _drifted_session("test-organism-consistency-gsea")
try:
    API.transcriptomics_gsea_configure({"databases": ["GO"]}, gsea_sid)
    check(False, "GSEA configure should have refused a drifted session")
except API.HTTPException as exc:
    check(exc.status_code == 400 and "Organism mismatch" in exc.detail,
          f"GSEA /configure refuses a drifted session (got {exc.status_code}: {exc.detail!r})")

regulon_sid = _drifted_session("test-organism-consistency-regulon")
try:
    # "analyses" omitted deliberately - the endpoint defaults it to
    # TRANSCRIPTOMICS_REGULON_ANALYSES, and the organism guard fires before
    # any of that validation runs, so a minimal payload is enough here.
    API.transcriptomics_regulon_configure({"databases": ["TRRUST"]}, regulon_sid)
    check(False, "Regulon configure should have refused a drifted session")
except API.HTTPException as exc:
    check(exc.status_code == 400 and "Organism mismatch" in exc.detail,
          f"Regulon /configure refuses a drifted session (got {exc.status_code}: {exc.detail!r})")

# And a CONSISTENT session must NOT be refused by the guard - confirms the
# three checks above are catching the drift specifically, not e.g. missing
# required payload fields that would 400 either way.
ok_ora_sid = _consistent_session("test-organism-consistency-ora-ok")
try:
    API.transcriptomics_ora_configure({"databases": ["GO"]}, ok_ora_sid)
    check(True, "ORA /configure accepts a CONSISTENT session (got past the guard)")
except API.HTTPException as exc:
    check("Organism mismatch" not in exc.detail,
          f"a consistent session must not fail on the ORGANISM guard specifically "
          f"(got a different, acceptable 400: {exc.detail!r})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
