"""End-to-end HTTP test of the new Proteomics PTM API against a running
uvicorn: real upload -> resolve (real UniProt ID Mapping) -> configure ->
run (real iPTMnet fetch, kinase activity + substrate ORA) -> download, on
the bundled example file (real human accessions/sites, built from iPTMnet's
own data - see Data3/example_proteomics_phosphosites.csv's own generation:
25 real MAPK1 substrate sites simulated strongly UP, 15 real AKT1 substrate
sites simulated DOWN, 260 background sites as noise)."""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
    if os.path.basename(os.path.dirname(os.path.abspath(__file__))) == "tests" \
    else r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
DATA3 = os.path.join(_ROOT, "Data3")
SESSION_ID = "px-ptm-e2e"
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


def poll(job_id, timeout=300):
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
        time.sleep(1.0)
    raise TimeoutError("job did not finish")


print("=== guards: configure without a resolved table is refused ===")
r = requests.post(f"{BASE}/api/proteomics/ptm/configure",
                  headers={"X-Session-Id": "px-ptm-empty"}, json={}, timeout=60)
check(r.status_code == 400, f"configure without a resolved table is refused ({r.status_code})")

print("=== upload + resolve (real UniProt ID Mapping) ===")
with open(os.path.join(DATA3, "example_proteomics_phosphosites.csv"), "rb") as fh:
    r = requests.post(f"{BASE}/api/proteomics/ptm/upload", headers=HEADERS,
                      files={"file": ("example.csv", fh, "text/csv")}, timeout=120)
check(r.status_code == 200, f"upload succeeded ({r.status_code})")
check(r.json().get("n_rows") == 300, f"300 rows parsed (got {r.json().get('n_rows')})")

r = requests.post(f"{BASE}/api/proteomics/ptm/resolve", headers=HEADERS,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
check(r.status_code == 200, f"resolve accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished ({job['status']}, {job.get('error')})")
meta = job["result"]["metrics"]
print(f"    mapped {meta.get('n_mapped')}/{meta.get('n_input_rows')}, "
      f"localization-dropped {meta.get('n_dropped_localization')}")
check(meta.get("n_mapped", 0) > 100, f"a meaningful number of sites resolved (got {meta.get('n_mapped')})")

print("=== isoform preservation check (real data) ===")
resolved_table = job["result"]["resolved_table"]
cols = resolved_table["columns"]
check("accession" in cols and "canonical_accession" in cols and "site_id" in cols and "canonical_site_id" in cols,
      f"resolved table carries the isoform-aware columns (got {cols})")
rows = resolved_table["rows"]
isoform_rows = [r for r in rows if r.get("isoform")]
check(len(isoform_rows) > 0, f"at least one isoform-suffixed site survived resolution (got {len(isoform_rows)})")
if isoform_rows:
    r0 = isoform_rows[0]
    check(r0["accession"] != r0["canonical_accession"] and r0["accession"].startswith(r0["canonical_accession"]),
          f"an isoform row's accession differs from its canonical_accession ({r0['accession']} vs {r0['canonical_accession']})")

print("=== configure ===")
r = requests.post(f"{BASE}/api/proteomics/ptm/configure", headers=HEADERS,
                  json={"padj_cutoff": 0.05, "lfc_cutoff": 1.0, "min_term_size": 3,
                        "min_kinase_sites": 5, "metric": "log2fc", "n_permutations": 500},
                  timeout=60)
check(r.status_code == 200, f"configure accepted ({r.status_code}: {r.text[:300]})")
cfg = r.json().get("config", {})
check(cfg.get("n_permutations") == 500, f"config echoes n_permutations (got {cfg.get('n_permutations')})")

print("=== run (real iPTMnet fetch - cached from generating the example, so should be fast) ===")
r = requests.post(f"{BASE}/api/proteomics/ptm/run", headers=HEADERS, timeout=60)
check(r.status_code == 200, f"run accepted ({r.status_code})")
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"PTM job finished ({job['status']}, {job.get('error')})")

result = job["result"]
check("activity" in result and "ora" in result, f"result carries both activity and ora keys (got {list(result.keys())})")

activity_tables = result["activity"]["tables"]
check("iPTMnet" in activity_tables, f"activity result wrapped under 'iPTMnet' key (got {list(activity_tables.keys())})")
act_block = activity_tables["iPTMnet"]
act_rows = act_block["table"]["rows"]
print(f"    kinases tested: {act_block.get('n_kinases_tested')}, top rows:")
for row in act_rows[:5]:
    print(f"      {row.get('Kinase')}: score={row.get('Activity Score'):.3f} "
          f"FDR(calibrated)={row.get('FDR (calibrated)'):.3g} n={row.get('Substrate Sites Tested')}")

check(act_block.get("n_kinases_tested", 0) > 0, f"at least one kinase tested (got {act_block.get('n_kinases_tested')})")
check(len(act_rows) > 0 and act_rows[0]["Kinase"] == "MAPK1",
      f"MAPK1 (the kinase whose substrates were simulated strongly UP) is the TOP-ranked kinase "
      f"(got top={act_rows[0]['Kinase'] if act_rows else None})")
if act_rows:
    check(act_rows[0]["Activity Score"] > 0, f"MAPK1's activity score is positive (got {act_rows[0]['Activity Score']})")
    check(act_rows[0]["FDR (calibrated)"] < 0.01,
          f"MAPK1's calibrated FDR is genuinely significant (got {act_rows[0]['FDR (calibrated)']})")
akt1_row = next((r for r in act_rows if r["Kinase"] == "AKT1"), None)
check(akt1_row is not None and akt1_row["Activity Score"] < 0,
      f"AKT1 (substrates simulated DOWN) has a negative activity score (got {akt1_row})")

ora_tables = result["ora"]["tables"]
up_key = next((k for k in ora_tables if "(up)" in k), None)
check(up_key is not None, f"an 'up' direction ORA table exists (got {list(ora_tables.keys())})")
if up_key:
    up_rows = ora_tables[up_key]["table"]["rows"]
    check(len(up_rows) > 0 and up_rows[0]["Term Name"] == "MAPK1",
          f"MAPK1 is the top enriched kinase among up-regulated sites (got {up_rows[0] if up_rows else None})")

print("=== downloads ===")
r = requests.get(f"{BASE}/api/proteomics/ptm/download", headers=HEADERS, params={"which": "activity"}, timeout=60)
check(r.status_code == 200 and r.headers.get("content-type", "").startswith("text/csv"),
      f"activity CSV download works ({r.status_code})")
r = requests.get(f"{BASE}/api/proteomics/ptm/download", headers=HEADERS,
                 params={"which": "ora", "key": up_key}, timeout=60)
check(r.status_code == 200 and r.headers.get("content-type", "").startswith("text/csv"),
      f"ora CSV download works ({r.status_code})")
r = requests.get(f"{BASE}/api/proteomics/ptm/download/resolved.csv", headers=HEADERS, timeout=60)
check(r.status_code == 200, f"resolved-sites CSV download works ({r.status_code})")
r = requests.get(f"{BASE}/api/proteomics/ptm/download", headers=HEADERS, params={"which": "bogus"}, timeout=60)
check(r.status_code == 400, f"an unknown table name is rejected with 400 (got {r.status_code})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
