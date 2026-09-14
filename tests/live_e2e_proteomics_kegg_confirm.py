"""Real HTTP end-to-end test for KEGG in Proteomics ORA/Quantitative
Enrichment, including the strain-disambiguation picker - the SAME
generalized endpoints (/api/proteomics/kegg/*) and confirmation machinery
already proven for Transcriptomics
(tests/live_e2e_transcriptomics_kegg_confirm.py), now exercised through
Proteomics' own session key and job runners.

Mirrors that file's structure closely on purpose: the point of this test is
to prove the GENERALIZATION (Transcriptomics-only code made module-generic)
didn't silently narrow or break anything for the second module that now
depends on it. Uses the same real, reproducible KEGG ambiguity (Bacillus
subtilis, NCBI taxon 1423 - 14 strain-level KEGG genomes, none carrying that
taxon - confirmed live in tests/live_transcriptomics_kegg.py and again here).

Needs a live server:  python -m uvicorn API:app --port 8009
Run:                  python tests/live_e2e_proteomics_kegg_confirm.py
"""
import io
import sys
import time
import uuid

import requests

BASE = "http://127.0.0.1:8009"
SID = "proteomics-kegg-confirm-e2e-" + uuid.uuid4().hex[:8]
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


# Real human accessions with real KEGG pathway membership (same ones
# verified in tests/live_proteomics_kegg.py) - real biology, small universe,
# fast test.
GLYCOLYSIS = ["P04406", "P06733", "P14618", "P04075"]  # GAPDH, ENO1, PKM, ALDOA
CELL_CYCLE = ["P06493", "P14635"]  # CDK1, CCNB1
OTHER = ["P60709", "P68371", "P62826", "P07437", "P11021", "P08238", "P30101"]

lines = ["accession,log2FC,pvalue,padj"]
for a in GLYCOLYSIS:
    lines.append(f"{a},3.0,1e-8,1e-7")
for a in CELL_CYCLE:
    lines.append(f"{a},-3.0,1e-8,1e-7")
for a in OTHER:
    lines.append(f"{a},0.05,0.9,0.95")
csv_bytes = ("\n".join(lines) + "\n").encode("utf-8")

print("=== upload + resolve (real human accessions) ===")
r = requests.post(f"{BASE}/api/proteomics/upload", headers=H,
                  files={"file": ("de.csv", io.BytesIO(csv_bytes), "text/csv")}, timeout=60)
check(r.status_code == 200, f"upload accepted (HTTP {r.status_code}: {r.text[:200]})")

r = requests.post(f"{BASE}/api/proteomics/resolve", headers=H,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished (status={job['status']}, {str(job.get('error'))[:200]})")

print("\n=== ORA: configure + run with KEGG (clean, automatic organism) ===")
r = requests.post(f"{BASE}/api/proteomics/ora/configure", headers=H, json={
    "databases": ["KEGG"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
    "min_term_size": 2, "max_term_size": 500,
}, timeout=60)
check(r.status_code == 200 and r.json()["config"]["databases"] == ["KEGG"],
      f"KEGG survives the allow-list filter over HTTP (HTTP {r.status_code}: {r.text[:300]})")
r = requests.post(f"{BASE}/api/proteomics/ora/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"ORA run finished (status={job['status']}, {str(job.get('error'))[:300]})")
tables = job["result"]["tables"]
check(sorted(tables.keys()) == ["KEGG (all)", "KEGG (down)", "KEGG (up)"],
      f"three direction families came back over HTTP (got {sorted(tables.keys())})")
up_ids = [row["Term ID"] for row in tables["KEGG (up)"]["table"]["rows"]]
check("hsa00010" in up_ids, f"glycolysis is in the up results (got {up_ids[:6]})")
down_ids = [row["Term ID"] for row in tables["KEGG (down)"]["table"]["rows"]]
check("hsa04110" in down_ids, f"cell cycle is in the down results (got {down_ids[:6]})")
overlap_key = "Overlapping Proteins"
up_row = next((row for row in tables["KEGG (up)"]["table"]["rows"] if row["Term ID"] == "hsa00010"), None)
check(up_row is not None and "P04406" in up_row.get(overlap_key, "") and "GAPDH" in up_row.get(overlap_key, ""),
      f"members show as ACCESSION (SYMBOL) over HTTP too (got {up_row.get(overlap_key) if up_row else None!r})")

print("\n=== Quantitative Enrichment: configure + run with KEGG ===")
r = requests.post(f"{BASE}/api/proteomics/quant/configure", headers=H, json={
    "databases": ["KEGG"], "metric": "log2fc", "min_term_size": 2, "max_term_size": 500,
}, timeout=60)
check(r.status_code == 200, f"quant configure succeeded (HTTP {r.status_code}: {r.text[:300]})")
r = requests.post(f"{BASE}/api/proteomics/quant/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"quant run finished (status={job['status']}, {str(job.get('error'))[:300]})")
qtables = job["result"]["tables"]
check(list(qtables.keys()) == ["KEGG"], f"one KEGG family, no direction suffix (got {list(qtables.keys())})")
check(len(qtables["KEGG"]["table"]["rows"]) > 0, "quant found real KEGG terms")

print("\n=== switch to Bacillus subtilis (real ambiguous case) ===")
r = requests.post(f"{BASE}/api/proteomics/resolve", headers=H,
                  json={"taxon_id": "1423", "taxon_name": "Bacillus subtilis"}, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"resolve to B. subtilis finished (status={job['status']})")

print("\n=== /api/proteomics/kegg/resolve-organism: reports the real ambiguity ===")
r = requests.post(f"{BASE}/api/proteomics/kegg/resolve-organism", headers=H, json={}, timeout=60)
check(r.status_code == 200, f"resolve-organism succeeded (HTTP {r.status_code})")
result = r.json()
check(result["organism_code"] == "", f"no automatic resolution (got {result['organism_code']!r})")
codes = [c["organism_code"] for c in result.get("candidates", [])]
check(len(codes) >= 10 and "bsu" in codes,
      f"real strain-level candidates reported, including canonical 'bsu' (got {len(codes)}: {codes[:5]})")

print("\n=== /api/proteomics/kegg/confirm-organism: rejects a non-candidate, accepts a real one ===")
r = requests.post(f"{BASE}/api/proteomics/kegg/confirm-organism", headers=H,
                  json={"organism_code": "hsa"}, timeout=60)
check(r.status_code == 400, f"a non-candidate code is rejected (HTTP {r.status_code})")

r = requests.post(f"{BASE}/api/proteomics/kegg/confirm-organism", headers=H,
                  json={"organism_code": "bsu"}, timeout=60)
check(r.status_code == 200, f"a real candidate is accepted (HTTP {r.status_code}: {r.text[:300]})")
confirm_result = r.json()
check(confirm_result["organism_code"] == "bsu" and confirm_result["verified"] is False,
      f"confirmed result carries the code, unverified (got {confirm_result})")

print("\n=== ORA run now uses the confirmed organism ===")
r = requests.post(f"{BASE}/api/proteomics/ora/configure", headers=H, json={
    "databases": ["KEGG"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
    "min_term_size": 2, "max_term_size": 500,
}, timeout=60)
r = requests.post(f"{BASE}/api/proteomics/ora/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"ORA run for B. subtilis finished (status={job['status']})")
kegg_all = job["result"]["tables"]["KEGG (all)"]
check("manually confirmed" in (kegg_all.get("note") or "").lower(),
      f"the result carries the manual-confirmation caveat (got {kegg_all.get('note')!r})")

print("\n=== changing organism back to human invalidates the confirmation ===")
r = requests.post(f"{BASE}/api/proteomics/resolve", headers=H,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"re-resolve to human finished (status={job['status']})")

r = requests.post(f"{BASE}/api/proteomics/kegg/resolve-organism", headers=H, json={}, timeout=60)
result2 = r.json()
check(result2["organism_code"] == "hsa" and result2["verified"] is True,
      f"human resolves automatically and cleanly again, no leftover B. subtilis confirmation "
      f"(got {result2['organism_code']!r}, verified={result2.get('verified')})")

r = requests.post(f"{BASE}/api/proteomics/ora/configure", headers=H, json={
    "databases": ["KEGG"], "padj_cutoff": 0.05, "lfc_cutoff": 1.0,
    "min_term_size": 2, "max_term_size": 500,
}, timeout=60)
r = requests.post(f"{BASE}/api/proteomics/ora/run", headers=H, timeout=60)
job = wait(r.json()["job_id"])
check(job["status"] == "done", f"ORA run for human finished (status={job['status']})")
kegg_all2 = job["result"]["tables"]["KEGG (all)"]
check("confirmed" not in (kegg_all2.get("note") or "").lower(),
      f"SENSITIVITY: no leftover B. subtilis confirmation caveat leaked into the human run "
      f"(got {kegg_all2.get('note')!r})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
