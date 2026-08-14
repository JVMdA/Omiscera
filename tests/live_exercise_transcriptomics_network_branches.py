"""Evaluation Part 2: exercises code paths in the Transcriptomics Network
module that the original build (tests/live_e2e_transcriptomics_network.py)
never actually ran even once - abs_metric seed weighting, direction=up/down,
MSigDB module labelling, the enrichment CSV download, and the Ensembl/Entrez
identifier fallback tiers. Each of these was previously "tested" only in the
sense that the code existed and didn't obviously error; this file runs them
for real, against a live server and the REAL cached human STRING network, and
asserts on the actual output.

The Ensembl/Entrez fallback needs identifiers that are real (so STRING's own
alias file resolves them) but whose SYMBOL does not match STRING's preferred
name (so the symbol tier is skippped and the fallback genuinely fires). Real
STRING protein ids/aliases for TP53 and MYC, read directly from the already-
cached human aliases file, are used for this - see the comment at that
section for how they were obtained.
"""
import os
import sys
import time

import requests

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
    if os.path.basename(os.path.dirname(os.path.abspath(__file__))) == "tests" \
    else r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
sys.path.insert(0, _ROOT)
DATA3 = os.path.join(_ROOT, "Data3")
SESSION_ID = "tx-network-branches"
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


def configure_and_run(config):
    r = requests.post(f"{BASE}/api/transcriptomics/network/configure", headers=HEADERS,
                      json=config, timeout=120)
    if r.status_code != 200:
        return None, r
    r2 = requests.post(f"{BASE}/api/transcriptomics/network/run", headers=HEADERS, timeout=60)
    if r2.status_code != 200:
        return None, r2
    job = poll(r2.json()["job_id"])
    return job, r2


print("=== upload + resolve (real mygene.info) ===")
with open(os.path.join(DATA3, "example_transcriptomics_de_full.csv"), "rb") as fh:
    r = requests.post(f"{BASE}/api/transcriptomics/upload", headers=HEADERS,
                      files={"file": ("example.csv", fh, "text/csv")},
                      data={"identifier_type": "gene_symbol"}, timeout=120)
check(r.status_code == 200, f"upload succeeded ({r.status_code})")
r = requests.post(f"{BASE}/api/transcriptomics/resolve", headers=HEADERS,
                  json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}, timeout=60)
job = poll(r.json()["job_id"])
check(job["status"] == "done", f"resolve finished ({job['status']}, {job.get('error')})")

# =============================================================================
print("=== abs_metric seed weighting: never run before, now for real ===")
# =============================================================================
job, resp = configure_and_run({"confidence": 700, "n_permutations": 100, "module_detection": False,
                               "direction": "all", "weighting": "abs_metric", "metric": "pi_value"})
check(job is not None and job["status"] == "done",
      f"abs_metric run completed ({resp.status_code if job is None else job['status']})")
if job:
    seeds = (job.get("result") or {}).get("seeds") or {}
    check(seeds.get("weighting") == "abs_metric" and seeds.get("winsor_q") is not None,
          f"the run actually used abs_metric weighting and recorded a real winsorization quantile "
          f"(weighting={seeds.get('weighting')}, winsor_q={seeds.get('winsor_q')})")
    check((seeds.get("weight_ratio") or 0) <= 20.0 + 1e-6,
          f"the max/min weight ratio is bounded at <=20 on real data (got {seeds.get('weight_ratio')})")

# =============================================================================
print("=== direction=up / direction=down: never run before, now for real ===")
# =============================================================================
up_job, _ = configure_and_run({"confidence": 700, "n_permutations": 100, "module_detection": False,
                               "direction": "up", "weighting": "equal"})
down_job, _ = configure_and_run({"confidence": 700, "n_permutations": 100, "module_detection": False,
                                 "direction": "down", "weighting": "equal"})
check(up_job is not None and up_job["status"] == "done", "direction=up run completed")
check(down_job is not None and down_job["status"] == "done", "direction=down run completed")
if up_job and down_job:
    n_up = ((up_job.get("result") or {}).get("seeds") or {}).get("n_seeds", 0)
    n_down = ((down_job.get("result") or {}).get("seeds") or {}).get("n_seeds", 0)
    n_all = ((job.get("result") or {}).get("seeds") or {}).get("n_seeds", 0) if job else 0
    check(0 < n_up < n_all and 0 < n_down < n_all,
          f"up ({n_up}) and down ({n_down}) each seed a real, STRICT SUBSET of the all-direction seed "
          f"count ({n_all}) - proving direction actually restricts seed selection, not just relabels it")
    up_params = ((up_job.get("result") or {}).get("parameters") or {})
    check(up_params.get("direction") == "up", "the recorded parameters honestly reflect direction=up")

# =============================================================================
print("=== MSigDB module labelling: never run before, now for real ===")
# =============================================================================
msig_job, msig_resp = configure_and_run({
    "confidence": 700, "n_permutations": 100, "module_detection": True,
    "min_module_size": 3, "direction": "all", "weighting": "equal",
    "enrichment_database": "MSigDB", "msigdb_collections": ["hallmark"],
})
check(msig_job is not None and msig_job["status"] == "done",
      f"an MSigDB-labelled run completed ({msig_resp.status_code if msig_job is None else msig_job['status']}, "
      f"{(msig_job or {}).get('error')})")
mod_key = None
if msig_job:
    mods = (msig_job.get("result") or {}).get("modules") or {}
    rows = (mods.get("table") or {}).get("rows") or []
    check(len(rows) > 0, f"modules were detected under MSigDB labelling ({len(rows)})")
    labelled = [row for row in rows if row.get("Enrichment database", "").startswith("MSigDB")]
    check(len(labelled) > 0,
          f"at least one module's top term actually came from MSigDB, not an empty/GO fallback "
          f"({len(labelled)}/{len(rows)} modules labelled)")
    if rows:
        mod_key = rows[0]["Module ID"]

# =============================================================================
print("=== download?which=enrichment&key=<module>: never called before, now for real ===")
# =============================================================================
if mod_key:
    r = requests.get(f"{BASE}/api/transcriptomics/network/download", headers=HEADERS,
                     params={"which": "enrichment", "key": mod_key}, timeout=120)
    check(r.status_code == 200 and r.headers.get("content-type", "").startswith("text/csv"),
          f"enrichment CSV for module {mod_key!r} downloads (status {r.status_code})")
    check(len(r.text.splitlines()) >= 1, "the downloaded CSV has at least a header row")
else:
    check(False, "no module id was available to exercise the enrichment download (see MSigDB run above)")

r = requests.get(f"{BASE}/api/transcriptomics/network/download", headers=HEADERS,
                 params={"which": "enrichment", "key": "NOT_A_REAL_MODULE"}, timeout=60)
check(r.status_code == 400, f"an unknown module id in the enrichment download is rejected (got {r.status_code})")

# =============================================================================
print("=== Ensembl / Entrez fallback tiers: real STRING aliases, symbols that CANNOT match ===")
# =============================================================================
# TP53 (STRING protein 9606.ENSP00000269305, Ensembl gene ENSG00000141510) and
# MYC (9606.ENSP00000478887, Entrez 4609) - both read directly from the real
# cached 9606 aliases file. Each row's gene_symbol is a value STRING will
# never recognise, so the symbol tier cannot fire; only the alias fallback can
# resolve them.
import pandas as pd  # noqa: E402
import Genomics_network as gnet  # noqa: E402
import Transcriptomics_network as txnet  # noqa: E402

real_graph = gnet.build_string_network("9606", confidence=700)
real_paths = gnet._string_file_paths("9606")
fallback_df = pd.DataFrame([
    {"input_id": "f1", "gene_symbol": "NOTAREALSYMBOL_TP53", "ensembl_id": "ENSG00000141510",
     "entrez_id": "", "log2FC": 2.0, "padj": 0.001, "pvalue": 0.0001},
    {"input_id": "f2", "gene_symbol": "NOTAREALSYMBOL_MYC", "ensembl_id": "", "entrez_id": "4609",
     "log2FC": -2.0, "padj": 0.001, "pvalue": 0.0001},
])
fallback_map = txnet.map_de_table_to_string(fallback_df, real_graph, aliases_path=real_paths.get("aliases"),
                                            taxon_id="9606")
check(fallback_map["n_by_tier"]["ensembl_id"] == 1 and fallback_map["n_by_tier"]["entrez_id"] == 1
      and fallback_map["n_by_tier"]["gene_symbol"] == 0,
      f"both rows are unreachable by symbol and mapped ONLY via their real alias tier "
      f"(got n_by_tier={fallback_map['n_by_tier']})")
tp53_row = fallback_map["table"].set_index("input_id").loc["f1"]
myc_row = fallback_map["table"].set_index("input_id").loc["f2"]
check(tp53_row["STRING ID"] == "9606.ENSP00000269305" and tp53_row["Mapped via"] == "ensembl_id",
      f"TP53 resolved via its real Ensembl gene id to the correct STRING protein (got {tp53_row['STRING ID']} "
      f"via {tp53_row['Mapped via']})")
check(myc_row["STRING ID"] == "9606.ENSP00000478887" and myc_row["Mapped via"] == "entrez_id",
      f"MYC resolved via its real Entrez id to the correct STRING protein (got {myc_row['STRING ID']} "
      f"via {myc_row['Mapped via']})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
