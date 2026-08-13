"""End-to-end через the REAL API on the NEW Ensembl example: download the
example exactly as the "Example file" button does, upload it, resolve, and run
TF/Regulon - asserting the activity table comes back NON-EMPTY.

This is the precise failure the old 18-gene ensembl example produced (min
regulon size 5 -> no TF clears the bar -> empty activity table), so it is
asserted against directly rather than assumed fixed.
"""
import io
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
SESSION_ID = "example-ensembl-e2e"
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


def poll_job(job_id, timeout=300):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = requests.get(f"{BASE}/api/jobs/{job_id}", headers=HEADERS, timeout=30)
        r.raise_for_status()
        j = r.json()
        if j.get("status") in ("done", "error"):
            return j
        time.sleep(1.0)
    raise TimeoutError("job did not finish")


print("=== download the Ensembl example exactly as the UI button does ===")
r = requests.get(f"{BASE}/api/transcriptomics/example", params={"identifier_type": "ensembl_gene"}, timeout=60)
check(r.status_code == 200, f"example endpoint returns 200 (got {r.status_code})")
content = r.content
n_rows = len(content.decode("utf-8").strip().splitlines()) - 1
print(f"    downloaded {n_rows} data rows")
check(n_rows > 1000, f"the served Ensembl example is the FULL one ({n_rows} rows, was 18 before this fix)")
first_col = content.decode("utf-8").splitlines()[0].split(",")[0].lstrip("﻿")
check(first_col == "Ensembl Gene ID", f"header is the Ensembl ID column (got {first_col!r})")

print("=== upload + resolve it ===")
r = requests.post(
    f"{BASE}/api/transcriptomics/upload", headers=HEADERS,
    files={"file": ("example_transcriptomics_ensembl_full.csv", io.BytesIO(content), "text/csv")},
    data={"identifier_type": "ensembl_gene"}, timeout=120,
)
check(r.status_code == 200, f"upload accepted (got {r.status_code}: {r.text[:250]})")

r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
                   json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
check(r.status_code == 200, f"resolve accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished (status={job['status']}, error={job.get('error')})")
n_resolved = len((job.get("result") or {}).get("resolved_table", {}).get("rows", []))
print(f"    resolved rows: {n_resolved}")
check(n_resolved > 1000, f"resolved table is large ({n_resolved} genes)")

print("=== configure + run TF/Regulon ===")
r = requests.post(f"{BASE}/api/transcriptomics/regulon/configure", headers=HEADERS,
                   json={"databases": ["TRRUST"], "analyses": ["enrichment", "activity"],
                         "min_regulon_size": 5, "min_term_size": 5}, timeout=30)
check(r.status_code == 200, f"configure accepted (got {r.status_code}: {r.text[:250]})")
r = requests.post(f"{BASE}/api/transcriptomics/regulon/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"run accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"], timeout=300)
check(job["status"] == "done", f"regulon job finished (status={job['status']}, error={job.get('error')})")

result = job.get("result") or {}
activity_tables = (result.get("activity") or {}).get("tables") or {}
trrust = activity_tables.get("TRRUST") or {}
n_tfs = trrust.get("n_tfs_tested", 0)
print(f"    TRRUST activity: {n_tfs} TFs tested, lambda={trrust.get('lambda_used')}")
check(n_tfs > 0,
      f"THE FIX: the Ensembl example now produces a NON-EMPTY activity table ({n_tfs} TFs) - "
      f"the old 18-gene example produced 0")

enrichment_tables = (result.get("enrichment") or {}).get("tables") or {}
enr_rows = len(((enrichment_tables.get("TRRUST (all)") or {}).get("table") or {}).get("rows", []))
print(f"    TRRUST (all) enrichment: {enr_rows} TFs tested")
check(enr_rows > 0, f"and a non-empty target-enrichment table ({enr_rows} TFs)")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
