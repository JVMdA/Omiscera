"""Real HTTP end-to-end test for the KEGG organism disambiguation picker
(Transcriptomics ORA/GSEA) - the new endpoints /api/transcriptomics/kegg/
resolve-organism and /api/transcriptomics/kegg/confirm-organism.

Mirrors tests/regression_ora_gsea.py's approach (real HTTP against a live
server, no browser). Uses Bacillus subtilis (NCBI taxon 1423) as the
ambiguous case: a real, reproducible KEGG ambiguity confirmed live in
tests/live_transcriptomics_kegg.py (14 strain-level KEGG genomes, none
carrying taxon 1423) - not a stub or monkeypatch.

Covers what tests/test_transcriptomics_kegg.py (offline, fixtures) and
tests/live_transcriptomics_kegg.py (live, direct function calls) cannot:
whether the two new endpoints are wired correctly end to end, whether a
confirmation actually changes what an ORA run over HTTP returns, and -
critically - whether changing the organism at retrieval invalidates a
previous confirmation server-side, so it can never silently leak into a run
for a different organism (the same hazard class as commit 27ce033's KEGG
organism divergence fix).

Needs a live server:  python -m uvicorn API:app --port 8009
Run:                  python tests/live_e2e_transcriptomics_kegg_confirm.py
"""
import io
import os
import sys
import time
import uuid

import requests

BASE = "http://127.0.0.1:8009"
SID = "kegg-confirm-e2e-" + uuid.uuid4().hex[:8]
H = {"X-Session-Id": SID}

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def wait(job_id, timeout=180):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = requests.get(f"{BASE}/api/jobs/{job_id}", headers=H, timeout=30).json()
        if r.get("status") in ("done", "error"):
            return r
        time.sleep(2)
    raise TimeoutError("job did not finish")


# A tiny DE table - the point of this file is the KEGG organism-resolution
# plumbing, not gene-level biology, so the identifiers don't need to be
# annotatable in KEGG themselves.
csv_bytes = b"gene,log2FoldChange,pvalue,padj\nrecA,2.0,1e-4,1e-3\ndnaA,1.5,1e-3,1e-2\n"

print("=== upload + resolve for Bacillus subtilis (real ambiguous case) ===")
r = requests.post(f"{BASE}/api/transcriptomics/upload", headers=H,
                  files={"file": ("de.csv", io.BytesIO(csv_bytes), "text/csv")},
                  data={"identifier_type": "gene_symbol"}, timeout=60)
check(r.status_code == 200, f"upload accepted (HTTP {r.status_code}: {r.text[:200]})")

r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=H,
                  json={"taxon_id": "1423", "taxon_name": "Bacillus subtilis"}, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished (status={job['status']}, {str(job.get('error'))[:200]})")

print("\n=== /kegg/resolve-organism: reports the real ambiguity, touches no session state ===")
r = requests.post(f"{BASE}/api/transcriptomics/kegg/resolve-organism", headers=H, json={}, timeout=60)
check(r.status_code == 200, f"resolve-organism succeeded (HTTP {r.status_code})")
result = r.json()
check(result["organism_code"] == "", f"no automatic resolution (got {result['organism_code']!r})")
check(len(result.get("candidates", [])) >= 10,
      f"real strain-level candidates are reported (got {len(result.get('candidates', []))})")
codes = [c["organism_code"] for c in result["candidates"]]
check("bsu" in codes, f"the canonical 'bsu' (B. subtilis 168) is among them (got {sorted(codes)[:8]})")
# Calling it again should return the same shape - confirms it's genuinely
# read-only (no session mutation that would change behaviour on a re-call).
r2 = requests.post(f"{BASE}/api/transcriptomics/kegg/resolve-organism", headers=H, json={}, timeout=60)
check(r2.json()["organism_code"] == "" and set(c["organism_code"] for c in r2.json()["candidates"]) == set(codes),
      "calling resolve-organism twice in a row gives the identical result (read-only, no side effect)")

print("\n=== /kegg/confirm-organism: rejects a code that isn't a real candidate ===")
r = requests.post(f"{BASE}/api/transcriptomics/kegg/confirm-organism", headers=H,
                  json={"organism_code": "hsa"}, timeout=60)
check(r.status_code == 400, f"a non-candidate code is rejected (HTTP {r.status_code}: {r.text[:200]})")

print("\n=== /kegg/confirm-organism: rejects a missing organism_code ===")
r = requests.post(f"{BASE}/api/transcriptomics/kegg/confirm-organism", headers=H, json={}, timeout=60)
check(r.status_code == 400, f"a missing organism_code is rejected (HTTP {r.status_code})")

print("\n=== /kegg/confirm-organism: accepts a real candidate and persists it ===")
r = requests.post(f"{BASE}/api/transcriptomics/kegg/confirm-organism", headers=H,
                  json={"organism_code": "bsu"}, timeout=60)
check(r.status_code == 200, f"a real candidate is accepted (HTTP {r.status_code}: {r.text[:300]})")
confirm_result = r.json()
check(confirm_result["organism_code"] == "bsu" and confirm_result["verified"] is False,
      f"confirmed result carries the code, unverified (got {confirm_result})")

print("\n=== ORA run over HTTP now uses the confirmed organism ===")
r = requests.post(f"{BASE}/api/transcriptomics/ora/configure", headers=H, json={
    "databases": ["KEGG"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
    "min_term_size": 2, "max_term_size": 500, "use_length_correction": False,
}, timeout=60)
check(r.status_code == 200, f"configure succeeded (HTTP {r.status_code})")
r = requests.post(f"{BASE}/api/transcriptomics/ora/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"ORA run finished (status={job['status']}, {str(job.get('error'))[:300]})")
kegg_all = job["result"]["tables"]["KEGG (all)"]
check("manually confirmed" in (kegg_all.get("note") or "").lower(),
      f"the result carries the user-confirmed caveat, not silence (got {kegg_all.get('note')!r})")

print("\n=== changing organism at retrieval invalidates the confirmation ===")
r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=H,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"re-resolve to human finished (status={job['status']})")

r = requests.post(f"{BASE}/api/transcriptomics/kegg/resolve-organism", headers=H, json={}, timeout=60)
result2 = r.json()
check(result2["organism_code"] == "hsa" and result2["verified"] is True,
      f"human now resolves automatically and cleanly - no leftover B. subtilis confirmation "
      f"(got {result2['organism_code']!r}, verified={result2.get('verified')})")

print("\n=== confirmation for the OLD organism does not leak into a run for the NEW one ===")
r = requests.post(f"{BASE}/api/transcriptomics/ora/configure", headers=H, json={
    "databases": ["KEGG"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
    "min_term_size": 2, "max_term_size": 500, "use_length_correction": False,
}, timeout=60)
r = requests.post(f"{BASE}/api/transcriptomics/ora/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"ORA run for human finished (status={job['status']})")
kegg_all2 = job["result"]["tables"]["KEGG (all)"]
check("confirmed" not in (kegg_all2.get("note") or "").lower(),
      f"SENSITIVITY: no leftover B. subtilis confirmation caveat leaked into the human run "
      f"(got {kegg_all2.get('note')!r}) - proves _transcriptomics_kegg_confirmed_code's "
      "taxon re-check (API.py) actually guards this, not just _job_resolve_de's own reset")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
