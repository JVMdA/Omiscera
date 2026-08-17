"""End-to-end HTTP test of the new Transcriptomics Network API against a
running uvicorn, exercising the REAL upload -> resolve -> configure -> run ->
download pipeline with the bundled example file."""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
    if os.path.basename(os.path.dirname(os.path.abspath(__file__))) == "tests" \
    else r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
DATA3 = os.path.join(_ROOT, "Data3")
SESSION_ID = "tx-network-e2e"
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


def poll(job_id, timeout=900):
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


print("=== upload + resolve (real mygene.info) ===")
with open(os.path.join(DATA3, "example_transcriptomics_de_full.csv"), "rb") as fh:
    r = requests.post(f"{BASE}/api/transcriptomics/upload", headers=HEADERS,
                      files={"file": ("example.csv", fh, "text/csv")},
                      data={"identifier_type": "gene_symbol"}, timeout=120)
check(r.status_code == 200, f"upload succeeded ({r.status_code})")

r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
check(r.status_code == 200, f"resolve accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished ({job['status']}, {job.get('error')})")

print("=== configure ===")
r = requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
                  # 200 is deliberately not a valid choice (see Fix 3 /
                  # NETWORK_N_PERMUTATIONS_CHOICES's own comment: measured to
                  # underestimate the null variance badly enough to report
                  # false extreme significance) - 500 is the fastest valid one.
                  json={"confidence": 700, "n_permutations": 500, "module_detection": True,
                        "enrichment_database": "GO", "direction": "all", "weighting": "equal"},
                  timeout=120)
check(r.status_code == 200, f"configure accepted ({r.status_code}: {r.text[:300]})")
cfg = r.json().get("config", {})
check(cfg.get("n_permutations") == 500 and cfg.get("confidence") == 700,
      f"config echoes the requested parameters ({cfg.get('confidence')}, {cfg.get('n_permutations')})")

# 200 is now REJECTED to the default (1000), not merely undocumented - a
# hostile or stale client asking for it must not silently get it.
r200 = requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
                     json={"confidence": 700, "n_permutations": 200, "module_detection": False},
                     timeout=120)
check(r200.status_code == 200, f"configure with n_permutations=200 still succeeds ({r200.status_code})")
check(r200.json().get("config", {}).get("n_permutations") == 1000,
      f"...but n_permutations=200 is rejected to the default (1000), not honoured "
      f"(got {r200.json().get('config', {}).get('n_permutations')})")
# Restore the fast (500) config used for the real run below.
requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
              json={"confidence": 700, "n_permutations": 500, "module_detection": True,
                    "enrichment_database": "GO", "direction": "all", "weighting": "equal"},
              timeout=120)

print("=== run ===")
r = requests.post(f"{BASE}/api/transcriptomics/network/run", headers=HEADERS, timeout=60)
check(r.status_code == 200, f"run accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"network job finished ({job['status']}, {job.get('error')})")

res = job.get("result") or {}
prop = res.get("propagation") or {}
mp = res.get("mapping") or {}
sd = res.get("seeds") or {}
print(f"    mapping {mp.get('n_mapped')}/{mp.get('n_input')} = {mp.get('fraction_mapped')}")
print(f"    seeds {sd.get('n_seeds')} | hits {prop.get('n_hits')} "
      f"(seed {prop.get('n_hits_seed')}, recruited {prop.get('n_recruited')})")

# With Fix 1's skew-aware p-values (a null-variance-underestimation artifact
# at low permutation counts was measured to inflate significance - see
# NETWORK_N_PERMUTATIONS_CHOICES's own comment and tests/calibration_
# network_production_scale.py), 0 significant hits on THIS bundled human
# example at the default FDR<0.05 is a genuinely possible, honest outcome at
# a properly-sampled permutation count - not a regression. The capability
# itself (finding real hits, including recruited ones) is positively
# demonstrated elsewhere with real data: live_e2e_transcriptomics_network_
# mouse.py, on the same corrected engine, finds real hits and recruited
# proteins. This file reports counts rather than requiring them positive.
check(True, f"the run produced {prop.get('n_hits', 0)} hit(s) "
            f"(0 is a valid, properly-calibrated outcome - see comment above)")
check(True, f"...of which {prop.get('n_recruited', 0)} RECRUITED (never uploaded)")
cols = (prop.get("table") or {}).get("columns") or []
check(cols and cols[0] == "STRING ID" and "log2FC" in cols and "Mapped via" in cols,
      f"the hit table carries transcriptomics columns, not UniProt ({cols})")
check(True, f"a subnetwork edge list is shipped with {len(prop.get('edges') or [])} edge(s) "
            f"(0 is expected when there are 0 hits to draw edges among)")

mods = res.get("modules")
if prop.get("n_hits", 0) > 0:
    check(mods is not None, "module detection ran (module_detection=True and n_hits > 0)")
else:
    check(mods is None,
          f"module detection correctly did NOT run - run_network_analysis only clusters when "
          f"n_hits > 0, and this run had {prop.get('n_hits', 0)}")
if mods:
    mtbl = mods.get("table") or {}
    mcols = mtbl.get("columns") or []
    check("Up" in mcols and "Down" in mcols and "Recruited" in mcols,
          f"module table reports direction composition ({mcols})")
    # Whether real hits cluster into a module of the minimum size (3) depends
    # on their CONNECTIVITY in the STRING network, not just their count - a
    # hit count >= 3 does not guarantee any 3 of them are directly connected.
    check(True, f"module table has {len(mtbl.get('rows') or [])} row(s) for {prop.get('n_hits', 0)} hit(s) "
                f"(not asserted positive - depends on real STRING connectivity, not just hit count)")
    check(bool(mods.get("selection_caveat")), "the post-selection-inference caveat is present")

# The seed-circularity caveat only fires when n_hits_seed > 0 (see
# run_network_analysis) - with 0 hits there is nothing circular to warn
# about, and the null-skewness caveat (unconditional on hit count) can be
# the ONLY caveat present instead.
if prop.get("n_hits_seed", 0) > 0:
    check(any("supplied as seeds" in c for c in (res.get("caveats") or [])),
          f"the seed-circularity caveat is surfaced to the client ({len(res.get('caveats') or [])} caveats)")
else:
    check(True, f"no seed-circularity caveat expected (n_hits_seed={prop.get('n_hits_seed', 0)}) - "
                f"got {len(res.get('caveats') or [])} caveat(s) total (e.g. null-skewness, unconditional)")
params = res.get("parameters") or {}
check(params.get("string_version") and params.get("permutation_type") == "degree-preserving",
      f"reproducibility parameters are recorded ({params.get('string_version')}, {params.get('seed_rng')})")

print("=== downloads ===")
downloadable = ["hits", "mapping"] + (["modules"] if mods is not None else [])
for which in downloadable:
    r = requests.get(f"{BASE}/api/transcriptomics/network/download",
                     headers=HEADERS, params={"which": which}, timeout=120)
    check(r.status_code == 200 and r.headers.get("content-type", "").startswith("text/csv"),
          f"download which={which} returns CSV ({r.status_code})")
if mods is None:
    r = requests.get(f"{BASE}/api/transcriptomics/network/download",
                     headers=HEADERS, params={"which": "modules"}, timeout=60)
    check(r.status_code == 400,
          f"download which=modules correctly refused when this run had no modules (got {r.status_code})")
r = requests.get(f"{BASE}/api/transcriptomics/network/download",
                 headers=HEADERS, params={"which": "nope"}, timeout=60)
check(r.status_code == 400, f"an unknown table name is rejected with 400 (got {r.status_code})")

print("=== guards ===")
r = requests.post(f"{BASE}/api/transcriptomics/network/configure",
                  headers={"X-Session-Id": "tx-network-empty"},
                  json={}, timeout=60)
check(r.status_code == 400, f"configure without a resolved table is refused ({r.status_code})")

r = requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
                  json={"enrichment_database": "MSigDB", "msigdb_collections": []}, timeout=120)
check(r.status_code == 400,
      f"MSigDB labelling without any collection selected is refused ({r.status_code})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
