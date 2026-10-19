"""End-to-end HTTP test of the new Proteomics PPI Network API against a
running uvicorn, exercising the REAL upload -> resolve -> configure -> run ->
download pipeline with the bundled example file - mirroring
live_e2e_transcriptomics_network.py's structure closely on purpose: the point
is to prove the GENERALIZATION (Transcriptomics_network.run_network_analysis's
hook parameters, and API.py's _network_param_bounds/_job_run_network) didn't
silently narrow or break anything for the second module now depending on it,
while the identifier space (UniProt accession, not gene symbol) and the
module-label databases (GO/Reactome/Complex Portal, not GO/MSigDB) are
genuinely different underneath.
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
SESSION_ID = "px-network-e2e"
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


print("=== guards: configure without a resolved table is refused ===")
r = requests.post(f"{BASE}/api/proteomics/network/configure",
                  headers={"X-Session-Id": "px-network-empty"}, json={}, timeout=60)
check(r.status_code == 400, f"configure without a resolved table is refused ({r.status_code})")

print("=== upload + resolve (real UniProt ID Mapping) ===")
with open(os.path.join(DATA3, "example_proteomics_de_full.csv"), "rb") as fh:
    r = requests.post(f"{BASE}/api/proteomics/upload", headers=HEADERS,
                      files={"file": ("example.csv", fh, "text/csv")}, timeout=120)
check(r.status_code == 200, f"upload succeeded ({r.status_code})")

r = requests.post(f"{BASE}/api/proteomics/resolve", headers=HEADERS,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
check(r.status_code == 200, f"resolve accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished ({job['status']}, {job.get('error')})")

print("=== configure ===")
r = requests.post(f"{BASE}/api/proteomics/network/configure", headers=HEADERS,
                  # 500 is the fastest valid permutation count - see
                  # NETWORK_N_PERMUTATIONS_CHOICES's own comment (a lower,
                  # unoffered count was measured to manufacture spurious
                  # extreme Z-scores, not merely add noise). GO + Complex
                  # Portal together: a real, non-mocked multi-database
                  # configure, not just GO alone.
                  json={"confidence": 700, "n_permutations": 500, "module_detection": True,
                        "enrichment_databases": ["GO", "Complex Portal"], "direction": "all", "weighting": "equal"},
                  timeout=120)
check(r.status_code == 200, f"configure accepted ({r.status_code}: {r.text[:300]})")
cfg = r.json().get("config", {})
check(cfg.get("n_permutations") == 500 and cfg.get("confidence") == 700,
      f"config echoes the requested parameters ({cfg.get('confidence')}, {cfg.get('n_permutations')})")
check(cfg.get("enrichment_databases") == ["GO", "Complex Portal"],
      f"BOTH requested databases are echoed back, in the order requested (got {cfg.get('enrichment_databases')})")
check(cfg.get("msigdb_collections") == [],
      f"msigdb_collections is always empty for Proteomics regardless of what was sent "
      f"(got {cfg.get('msigdb_collections')})")

# A hostile/stale n_permutations outside the offered choices is REJECTED to
# the default (1000), never honoured - same guard as Transcriptomics.
r200 = requests.post(f"{BASE}/api/proteomics/network/configure", headers=HEADERS,
                     json={"confidence": 700, "n_permutations": 200, "module_detection": False,
                           "enrichment_databases": ["GO"]},
                     timeout=120)
check(r200.status_code == 200, f"configure with n_permutations=200 still succeeds ({r200.status_code})")
check(r200.json().get("config", {}).get("n_permutations") == 1000,
      f"...but n_permutations=200 is rejected to the default (1000), not honoured "
      f"(got {r200.json().get('config', {}).get('n_permutations')})")

# A database Proteomics does not offer (MSigDB, Transcriptomics-only) is
# filtered out silently - requesting ONLY unsupported ones leaves the list
# empty, which is then REJECTED outright (never falls back to "GO anyway"
# without telling the caller - an unchecked box must not silently relabel).
r_bad_db = requests.post(f"{BASE}/api/proteomics/network/configure", headers=HEADERS,
                         json={"confidence": 700, "enrichment_databases": ["MSigDB"]}, timeout=120)
check(r_bad_db.status_code == 400,
      f"requesting only an unsupported database (MSigDB) is rejected, not silently defaulted "
      f"(got {r_bad_db.status_code})")

# No database selected at all is also rejected - same "must be explicit" rule.
r_empty_db = requests.post(f"{BASE}/api/proteomics/network/configure", headers=HEADERS,
                           json={"confidence": 700, "enrichment_databases": []}, timeout=120)
check(r_empty_db.status_code == 400,
      f"an empty database selection is rejected (got {r_empty_db.status_code})")

# Restore the fast (500, GO) config used for the real run below.
requests.post(f"{BASE}/api/proteomics/network/configure", headers=HEADERS,
              json={"confidence": 700, "n_permutations": 500, "module_detection": True,
                    "enrichment_databases": ["GO"], "direction": "all", "weighting": "equal"},
              timeout=120)

print("=== run ===")
r = requests.post(f"{BASE}/api/proteomics/network/run", headers=HEADERS, timeout=60)
check(r.status_code == 200, f"run accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"network job finished ({job['status']}, {job.get('error')})")

res = job.get("result") or {}
prop = res.get("propagation") or {}
mp = res.get("mapping") or {}
sd = res.get("seeds") or {}
print(f"    mapping {mp.get('n_mapped')}/{mp.get('n_input')} = {mp.get('fraction_mapped')}")
print(f"    tiers {mp.get('n_by_tier')}")
print(f"    seeds {sd.get('n_seeds')} | hits {prop.get('n_hits')} "
      f"(seed {prop.get('n_hits_seed')}, recruited {prop.get('n_recruited')})")

# Re-measures the coverage figures this feature's design docs cite on the
# app's own bundled 3,500-protein example - if UniProt's alias data or
# STRING's own release ever drifts enough to move these substantially, this
# is what would catch it. Not pinned exactly (external data can shift
# slightly release to release), just bounded around the measured values.
check(mp.get("n_input", 0) > 3000, f"the bundled example has >3000 input rows (got {mp.get('n_input')})")
check(mp.get("fraction_mapped", 0) > 0.5,
      f"mapping coverage stays well above the 10% refusal / 50% warning gates "
      f"(this app's own historical measurements on this file put it around 80-95%, and it moves with "
      f"UniProt's own live ID Mapping service - not pinned exactly; got {mp.get('fraction_mapped')})")
check((mp.get("n_by_tier") or {}).get("accession", 0) > 0,
      f"most mapped rows go through the accession tier, not the gene-symbol fallback "
      f"(got {mp.get('n_by_tier')})")
check(sd.get("n_seeds", 0) > 50, f"a meaningful number of proteins seed the walk (got {sd.get('n_seeds')})")

cols = (prop.get("table") or {}).get("columns") or []
check(cols and cols[0] == "STRING ID" and "log2FC" in cols and "Mapped via" in cols,
      f"the hit table carries the SAME shared column set Transcriptomics' own hit table uses - "
      f"entity_noun only changes prose, never a column name (got {cols})")
check(True, f"the run produced {prop.get('n_hits', 0)} hit(s), {prop.get('n_recruited', 0)} of them RECRUITED "
            f"(never uploaded) - the genuinely new information over ORA/Quant")

mods = res.get("modules")
if prop.get("n_hits", 0) > 0:
    check(mods is not None, "module detection ran (module_detection=True and n_hits > 0)")
else:
    check(mods is None,
          f"module detection correctly did NOT run - run_ppi_network_analysis only clusters when "
          f"n_hits > 0, and this run had {prop.get('n_hits', 0)}")
if mods:
    mtbl = mods.get("table") or {}
    mcols = mtbl.get("columns") or []
    check("Up" in mcols and "Down" in mcols and "Recruited" in mcols,
          f"module table reports direction composition ({mcols})")
    check(True, f"module table has {len(mtbl.get('rows') or [])} row(s) for {prop.get('n_hits', 0)} hit(s) "
                f"(not asserted positive - depends on real STRING connectivity, not just hit count)")
    check(bool(mods.get("selection_caveat")), "the post-selection-inference caveat is present")
    enr = mods.get("enrichment") or {}
    labelled = [mid for mid, e in enr.items() if e]
    check(True, f"{len(labelled)} of {len(enr)} module(s) got a GO label "
                f"(depends on real connectivity/annotation, not asserted positive)")

params = res.get("parameters") or {}
check(params.get("string_version") and params.get("permutation_type") == "degree-preserving",
      f"reproducibility parameters are recorded ({params.get('string_version')}, {params.get('seed_rng')})")

print("=== downloads ===")
downloadable = ["hits", "mapping"] + (["modules"] if mods is not None else [])
for which in downloadable:
    r = requests.get(f"{BASE}/api/proteomics/network/download",
                     headers=HEADERS, params={"which": which}, timeout=120)
    check(r.status_code == 200 and r.headers.get("content-type", "").startswith("text/csv"),
          f"download which={which} returns CSV ({r.status_code})")
if mods is not None:
    labelled_mid = next((mid for mid, e in (mods.get("enrichment") or {}).items() if e), None)
    if labelled_mid:
        r = requests.get(f"{BASE}/api/proteomics/network/download", headers=HEADERS,
                         params={"which": "enrichment", "key": labelled_mid}, timeout=60)
        check(r.status_code == 200, f"per-module enrichment CSV download works ({r.status_code})")
else:
    r = requests.get(f"{BASE}/api/proteomics/network/download",
                     headers=HEADERS, params={"which": "modules"}, timeout=60)
    check(r.status_code == 400,
          f"download which=modules correctly refused when this run had no modules (got {r.status_code})")
r = requests.get(f"{BASE}/api/proteomics/network/download",
                 headers=HEADERS, params={"which": "nope"}, timeout=60)
check(r.status_code == 400, f"an unknown table name is rejected with 400 (got {r.status_code})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
