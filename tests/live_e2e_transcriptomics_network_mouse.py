"""Evaluation Part 2: a full real run on a SECOND organism (mouse, taxon
10090), kept in its own file so a STRING mouse-network download failure is
diagnosable on its own rather than aborting the human-only suite.

Directly answers "does this work for every organism, not just the one it was
demoed on" - and is also the first time the Ensembl/Entrez fallback tiers and
the mapping-coverage gates get exercised against a NON-human STRING
nomenclature, where they matter most in practice.

Data3/example_transcriptomics_mouse_test.csv is a small hand-built table of
well-known mouse genes (Trp53, Myc, Cdkn1a, ...) with synthetic but
plausible log2FC/p-values - real identifiers, manufactured statistics, same
convention as this session's other synthetic-but-real fixtures.
"""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
    if os.path.basename(os.path.dirname(os.path.abspath(__file__))) == "tests" \
    else r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
DATA3 = os.path.join(_ROOT, "Data3")
SESSION_ID = "tx-network-mouse"
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


def poll(job_id, timeout=1200):
    t0 = time.time()
    last = ""
    while time.time() - t0 < timeout:
        r = requests.get(f"{BASE}/api/jobs/{job_id}", headers=HEADERS, timeout=30)
        r.raise_for_status()
        j = r.json()
        if j.get("message") and j["message"] != last:
            last = j["message"]
            print(f"    .. {j.get('progress', 0)*100:5.1f}% {last}")
        if j.get("status") in ("done", "error"):
            return j
        time.sleep(1.5)
    raise TimeoutError("job did not finish")


print("=== upload + resolve a MOUSE table (real mygene.info, taxon 10090) ===")
with open(os.path.join(DATA3, "example_transcriptomics_mouse_test.csv"), "rb") as fh:
    r = requests.post(f"{BASE}/api/transcriptomics/upload", headers=HEADERS,
                      files={"file": ("mouse_example.csv", fh, "text/csv")},
                      data={"identifier_type": "gene_symbol"}, timeout=120)
check(r.status_code == 200, f"upload succeeded ({r.status_code})")

r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
                  json={"taxon_id": "10090", "taxon_name": "Mus musculus"}, timeout=60)
check(r.status_code == 200, f"resolve accepted ({r.status_code}: {r.text[:300]})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished ({job['status']}, {job.get('error')})")

print("=== configure + run (this triggers the mouse STRING download if not already cached) ===")
r = requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
                  json={"confidence": 700, "n_permutations": 500, "module_detection": True,
                        "min_module_size": 2, "enrichment_database": "GO",
                        "direction": "all", "weighting": "equal"},
                  timeout=120)
check(r.status_code == 200, f"configure accepted for taxon 10090 ({r.status_code}: {r.text[:300]})")

r = requests.post(f"{BASE}/api/transcriptomics/network/run", headers=HEADERS, timeout=60)
check(r.status_code == 200, f"run accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"the mouse network job finished ({job['status']}, {job.get('error')})")

res = job.get("result") or {}
prop = res.get("propagation") or {}
mp = res.get("mapping") or {}
sd = res.get("seeds") or {}
print(f"    taxon {res.get('taxon_id')} ({res.get('taxon_name')})")
print(f"    mapping {mp.get('n_mapped')}/{mp.get('n_input')} = {mp.get('fraction_mapped')}")
print(f"    n_by_tier {mp.get('n_by_tier')}")
print(f"    seeds {sd.get('n_seeds')} | hits {prop.get('n_hits')} "
      f"(seed {prop.get('n_hits_seed')}, recruited {prop.get('n_recruited')}) | "
      f"n_proteins {prop.get('n_proteins')}")

check(res.get("taxon_id") == "10090", f"the run actually used taxon 10090, not a leftover human session ({res.get('taxon_id')})")
check((prop.get("n_proteins") or 0) > 1000,
      f"a real mouse STRING network was built, not a stub ({prop.get('n_proteins')} proteins)")
check((mp.get("n_mapped") or 0) >= 20,
      f"the great majority of this small hand-built table mapped to mouse STRING proteins "
      f"({mp.get('n_mapped')}/{mp.get('n_input')})")
check((sd.get("n_seeds") or 0) > 0, f"at least one mouse gene became a seed ({sd.get('n_seeds')})")
check((prop.get("n_hits") or 0) > 0, f"the mouse run produced significant proteins ({prop.get('n_hits')})")

mods = res.get("modules")
if mods:
    mtbl = mods.get("table") or {}
    check(len(mtbl.get("rows") or []) > 0, f"modules were detected on the mouse network ({len(mtbl.get('rows') or [])})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
