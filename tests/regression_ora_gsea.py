"""Targeted regression check for the Transcriptomics ORA and GSEA endpoints,
extended over two sessions:
  - Originally: does the GO path still work after adding the TF/Regulon
    module (which shares the session-state dict and API.py module).
  - Then: does the GO path still work after adding MSigDB, and does the new
    MSigDB path itself work end to end through the REAL HTTP API - configure
    validation ("MSigDB" requires msigdb_collections), the live job worker,
    and a real MSigDB download - a different layer than the direct
    Transcriptomics_module.run_ora() call in test_msigdb_fetcher.py or the
    browser-level check in playwright_regulon.py."""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
SESSION_ID = "regression-ora-gsea-session"
HEADERS = {"X-Session-Id": SESSION_ID}
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


print("=== upload + resolve ===")
with open(f"{DATA3}\\example_transcriptomics_de_full.csv", "rb") as fh:
    r = requests.post(
        f"{BASE}/api/transcriptomics/upload", headers=HEADERS,
        files={"file": ("example_transcriptomics_de_full.csv", fh, "text/csv")},
        data={"identifier_type": "gene_symbol"}, timeout=60,
    )
check(r.status_code == 200, f"upload succeeded (got {r.status_code})")
r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
                   json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=30)
check(r.status_code == 200, f"resolve accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished (status={job['status']}, error={job.get('error')})")

print("=== ORA (pre-existing, untouched code path) ===")
r = requests.post(f"{BASE}/api/transcriptomics/ora/configure", headers=HEADERS,
                   json={"databases": ["GO"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0, "min_term_size": 5}, timeout=30)
check(r.status_code == 200, f"ORA configure succeeded (got {r.status_code}: {r.text[:200]})")
r = requests.post(f"{BASE}/api/transcriptomics/ora/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"ORA run accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"], timeout=180)
check(job["status"] == "done", f"ORA job finished (status={job['status']}, error={job.get('error')})")
check(len((job.get("result") or {}).get("tables") or {}) > 0, "ORA produced result tables")

print("=== GSEA (pre-existing, untouched code path) ===")
r = requests.post(f"{BASE}/api/transcriptomics/gsea/configure", headers=HEADERS,
                   json={"databases": ["GO"], "min_term_size": 15, "n_permutations": 500}, timeout=30)
check(r.status_code == 200, f"GSEA configure succeeded (got {r.status_code}: {r.text[:200]})")
r = requests.post(f"{BASE}/api/transcriptomics/gsea/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"GSEA run accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"], timeout=180)
check(job["status"] == "done", f"GSEA job finished (status={job['status']}, error={job.get('error')})")
check(len((job.get("result") or {}).get("tables") or {}) > 0, "GSEA produced result tables")

print("=== MSigDB configure validation: 'MSigDB' requires msigdb_collections ===")
r = requests.post(f"{BASE}/api/transcriptomics/ora/configure", headers=HEADERS,
                   json={"databases": ["MSigDB"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0, "min_term_size": 5}, timeout=30)
check(r.status_code == 400, f"selecting MSigDB with no collections is rejected (got {r.status_code})")

print("=== ORA with MSigDB (real fetch: oncogenic collection, 189 sets - small and fast) ===")
r = requests.post(f"{BASE}/api/transcriptomics/ora/configure", headers=HEADERS,
                   json={"databases": ["MSigDB"], "msigdb_collections": ["oncogenic"],
                         "padj_cutoff": 0.05, "lfc_cutoff": 1.0, "min_term_size": 5}, timeout=30)
check(r.status_code == 200, f"ORA+MSigDB configure succeeded (got {r.status_code}: {r.text[:200]})")
r = requests.post(f"{BASE}/api/transcriptomics/ora/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"ORA+MSigDB run accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"], timeout=180)
check(job["status"] == "done", f"ORA+MSigDB job finished (status={job['status']}, error={job.get('error')})")
ora_msigdb_tables = (job.get("result") or {}).get("tables") or {}
check("MSigDB: Oncogenic (all)" in ora_msigdb_tables,
      f"result key uses the new 'MSigDB: <collection label> (direction)' convention (got {sorted(ora_msigdb_tables.keys())})")
ora_summary = (job.get("result") or {}).get("summary") or {}
check(isinstance(ora_summary.get("n_independent_families"), int),
      f"ORA summary reports n_independent_families for the multiple-testing caption (got {ora_summary.get('n_independent_families')!r})")

print("=== GSEA with MSigDB (real fetch: oncogenic collection) ===")
r = requests.post(f"{BASE}/api/transcriptomics/gsea/configure", headers=HEADERS,
                   json={"databases": ["MSigDB"], "msigdb_collections": ["oncogenic"],
                         "min_term_size": 5, "n_permutations": 200}, timeout=30)
check(r.status_code == 200, f"GSEA+MSigDB configure succeeded (got {r.status_code}: {r.text[:200]})")
r = requests.post(f"{BASE}/api/transcriptomics/gsea/run", headers=HEADERS, timeout=30)
check(r.status_code == 200, f"GSEA+MSigDB run accepted (got {r.status_code})")
job = poll_job(r.json()["job_id"], timeout=180)
check(job["status"] == "done", f"GSEA+MSigDB job finished (status={job['status']}, error={job.get('error')})")
gsea_msigdb_tables = (job.get("result") or {}).get("tables") or {}
check("MSigDB: Oncogenic" in gsea_msigdb_tables,
      f"GSEA result key has no direction suffix, unlike ORA's (got {sorted(gsea_msigdb_tables.keys())})")
gsea_summary = (job.get("result") or {}).get("summary") or {}
check(isinstance(gsea_summary.get("n_independent_families"), int),
      f"GSEA summary reports n_independent_families (got {gsea_summary.get('n_independent_families')!r})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
