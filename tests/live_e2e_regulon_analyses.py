"""Regression test for the ANALYSIS-SELECTION envelope bug.

Found by inspection and confirmed against the live API: the three regulon job
workers returned two different shapes. _job_run_transcriptomics_regulon_both
returned {"enrichment": {...}, "activity": {...}}, but the two
single-analysis workers returned a FLAT {"tables", "summary"}. The results
screen reads results.enrichment / results.activity, so unticking either
analysis on the config screen produced a COMPLETELY BLANK results screen -
no panes, empty section rail, no error anywhere.

Every prior test configured BOTH analyses, so the path was never exercised.
This asserts all three selections return a usable envelope.
"""
import os
import io
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA3 = os.path.join(_ROOT, "Data3")

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def poll(job_id, headers, timeout=300):
    t0 = time.time()
    while time.time() - t0 < timeout:
        j = requests.get(f"{BASE}/api/jobs/{job_id}", headers=headers, timeout=30).json()
        if j.get("status") in ("done", "error"):
            return j
        time.sleep(1.0)
    raise TimeoutError


def run_with(analyses, session_id):
    headers = {"X-Session-Id": session_id}
    with open(DATA3 + r"\example_transcriptomics_de_full.csv", "rb") as fh:
        content = fh.read()
    requests.post(f"{BASE}/api/transcriptomics/upload", headers=headers,
                  files={"file": ("ex.csv", io.BytesIO(content), "text/csv")},
                  data={"identifier_type": "gene_symbol"}, timeout=120).raise_for_status()
    j = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=headers,
                      json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60).json()
    poll(j["job_id"], headers)
    requests.post(f"{BASE}/api/transcriptomics/regulon/configure", headers=headers,
                  json={"databases": ["TRRUST"], "analyses": analyses,
                        "min_regulon_size": 5, "min_term_size": 5}, timeout=30).raise_for_status()
    j = requests.post(f"{BASE}/api/transcriptomics/regulon/run", headers=headers, timeout=30).json()
    done = poll(j["job_id"], timeout=300, headers=headers)
    assert done["status"] == "done", f"job failed: {done.get('error')}"
    return done.get("result") or {}


print("=== enrichment ONLY ===")
res = run_with(["enrichment"], "analyses-enr-only")
check("enrichment" in res, "result carries an 'enrichment' key (the results screen renders that pane)")
check("activity" not in res, "and no 'activity' key, since it was not requested")
enr_tables = (res.get("enrichment") or {}).get("tables") or {}
check(len(enr_tables) > 0, f"enrichment tables present ({len(enr_tables)} keys)")
check(any(k for k in res.keys()), "at least ONE pane will render - not a blank screen")

print("=== activity ONLY ===")
res = run_with(["activity"], "analyses-act-only")
check("activity" in res, "result carries an 'activity' key")
check("enrichment" not in res, "and no 'enrichment' key, since it was not requested")
act_tables = (res.get("activity") or {}).get("tables") or {}
check(len(act_tables) > 0, f"activity tables present ({len(act_tables)} keys)")
trrust = act_tables.get("TRRUST") or {}
check(trrust.get("n_tfs_tested", 0) > 0, f"activity-only run produced real results ({trrust.get('n_tfs_tested')} TFs)")
check("lambda_used" in trrust, "calibration diagnostics survive the activity-only path too")

print("=== BOTH ===")
res = run_with(["enrichment", "activity"], "analyses-both")
check("enrichment" in res and "activity" in res, "both keys present")
check(len(((res.get("activity") or {}).get("tables") or {})) > 0, "activity tables present in the both-path")
check(len(((res.get("enrichment") or {}).get("tables") or {})) > 0, "enrichment tables present in the both-path")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
