"""End-to-end HTTP test of the new TF/Regulon API against a running uvicorn
(127.0.0.1:8009), exercising the REAL upload -> resolve -> configure -> run
-> download pipeline (no mocking), using the app's own example
transcriptomics DE file."""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA3 = os.path.join(_ROOT, "Data3")
SESSION_ID = "regulon-e2e-test-session"
HEADERS = {"X-Session-Id": SESSION_ID}

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def poll_job(job_id, timeout=180):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = requests.get(f"{BASE}/api/jobs/{job_id}", headers=HEADERS, timeout=30)
        r.raise_for_status()
        j = r.json()
        if j.get("status") in ("done", "error"):
            return j
        time.sleep(1.0)
    raise TimeoutError(f"job {job_id} did not finish within {timeout}s")


print("=== upload + resolve (real mygene.info call, human) ===")
with open(f"{DATA3}\\example_transcriptomics_de_full.csv", "rb") as fh:
    r = requests.post(
        f"{BASE}/api/transcriptomics/upload",
        headers=HEADERS,
        files={"file": ("example_transcriptomics_de_full.csv", fh, "text/csv")},
        data={"identifier_type": "gene_symbol"},
        timeout=60,
    )
check(r.status_code == 200, f"upload succeeded (got {r.status_code}: {r.text[:300]})")

r = requests.post(
    f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
    json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=30,
)
check(r.status_code == 200, f"resolve accepted (got {r.status_code}: {r.text[:300]})")
job = poll_job(r.json()["job_id"])
check(job["status"] == "done", f"resolve job finished (status={job['status']}, error={job.get('error')})")

print("=== regulon configure (both analyses, all four databases) ===")
r = requests.post(
    f"{BASE}/api/transcriptomics/regulon/configure", headers=HEADERS,
    json={
        "databases": ["TRRUST", "DoRothEA", "CollecTRI", "ChEA3"],
        "analyses": ["enrichment", "activity"],
        "dorothea_levels": "AB",
        "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
        "min_term_size": 3, "min_regulon_size": 3,
    },
    timeout=30,
)
check(r.status_code == 200, f"configure accepted (got {r.status_code}: {r.text[:300]})")

print("=== regulon run (both) ===")
r = requests.post(f"{BASE}/api/transcriptomics/regulon/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"run accepted (got {r.status_code}: {r.text[:300]})")
job = poll_job(r.json()["job_id"], timeout=240)
check(job["status"] == "done", f"regulon job finished (status={job['status']}, error={job.get('error')})")

result = job.get("result") or {}
enrichment_tables = (result.get("enrichment") or {}).get("tables") or {}
activity_tables = (result.get("activity") or {}).get("tables") or {}
check(len(enrichment_tables) > 0, f"enrichment produced result tables (got {len(enrichment_tables)} keys)")
check(len(activity_tables) == 4, f"activity produced one block per requested database (got {len(activity_tables)})")

chea3_activity_note = (activity_tables.get("ChEA3") or {}).get("note")
check(chea3_activity_note is not None and "mode-of-regulation" in chea3_activity_note,
      f"ChEA3 activity block explains the unsigned skip over HTTP too (got: {chea3_activity_note})")

any_enrichment_rows = any(
    len((block.get("table") or {}).get("rows", [])) > 0 for block in enrichment_tables.values()
)
print(f"    (enrichment tables have any rows: {any_enrichment_rows} - informational, real data can legitimately be empty)")

any_activity_rows = any(
    len((block.get("table") or {}).get("rows", [])) > 0 for block in activity_tables.values()
)
print(f"    (activity tables have any rows: {any_activity_rows} - informational)")

print("=== enrichment CSV download ===")
some_enrich_key = next(iter(enrichment_tables.keys()))
r = requests.get(
    f"{BASE}/api/transcriptomics/regulon/download", headers=HEADERS,
    params={"key": some_enrich_key}, timeout=30,
)
check(r.status_code == 200, f"enrichment download succeeded for key {some_enrich_key!r} (got {r.status_code})")
check(r.headers.get("content-type", "").startswith("text/csv"), "enrichment download is a CSV response")

print("=== activity CSV download ===")
some_activity_key = next(iter(activity_tables.keys()))
r = requests.get(
    f"{BASE}/api/transcriptomics/regulon/download-activity", headers=HEADERS,
    params={"key": some_activity_key}, timeout=30,
)
check(r.status_code == 200, f"activity download succeeded for key {some_activity_key!r} (got {r.status_code})")

print("=== unknown key -> 400, not a crash ===")
r = requests.get(
    f"{BASE}/api/transcriptomics/regulon/download", headers=HEADERS,
    params={"key": "Not A Real Key"}, timeout=30,
)
check(r.status_code == 400, f"unknown table key returns 400 cleanly (got {r.status_code})")

print("=== configure without databases -> 400 ===")
r = requests.post(
    f"{BASE}/api/transcriptomics/regulon/configure", headers=HEADERS,
    json={"databases": [], "analyses": ["enrichment"]}, timeout=30,
)
check(r.status_code == 400, f"empty database list rejected with 400 (got {r.status_code})")

print("=== configure without a resolved table (fresh session) -> 400 ===")
r = requests.post(
    f"{BASE}/api/transcriptomics/regulon/configure", headers={"X-Session-Id": "regulon-e2e-empty-session"},
    json={"databases": ["TRRUST"], "analyses": ["enrichment"]}, timeout=30,
)
check(r.status_code == 400, f"no resolved table rejected with 400 (got {r.status_code})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
