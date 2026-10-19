"""Live end-to-end test of the five Metabolomics tiles through the HTTP API
(no browser - see playwright_metabolomics.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network
access: upload -> resolve (PubChem/KEGG name resolution + the local
LIPID MAPS/RefMet/ChEBI crosswalk) -> configure -> run -> download, for

    ORA (KEGG + SMPDB), Chemical Class, Lipid Class      [Data3/example_metabolomics_list.csv]
    Quantitative Enrichment                              [Data3/example_metabolomics_concentrations.csv]
    Pathway Topology Analysis                            EITHER of the two, user's choice

and asserts on the BIOLOGY the bundled examples plant, not just on shapes:
the list is mostly TCA-cycle/glycolysis intermediates plus a few lipids, and
the concentration table has TCA/glycolysis compounds UP and lipids DOWN in
the Treated group. Each database must independently recover that, which is
what makes agreement between them meaningful (nothing here is hardcoded to
a specific database's output).

KEGG is organism-specific: KEGG has no organism-specific compound links, and a
map's compounds do not vary by organism, so the organism selects WHICH MAPS are
tested (KEGG's list/pathway/<org>). Covered here against real KEGG: human terms
are hsa-prefixed, yeast gets a different (smaller) universe from the same shared
link file, an organism KEGG has no genome for is EMPTY with a note (never a
fallback to generic maps), an ambiguous organism (B. subtilis, 14 KEGG
strains) is empty until the user picks a KEGG organism explicitly (the KEGG
search box: /api/metabolomics/kegg/check-organism + kegg_organism_code in the
run config), and an explicit pick wins over the derived organism.

Also covers the input-contract errors that motivated the rework: a list with
no Name column is rejected with a message that says so (not a KeyError), a
concentration table with one group is rejected, and the layout/group
overrides on the Quantitative upload work.

Pathway Topology Analysis (KEGG only, either input mode): checked directly
against a fresh ORA/Quantitative KEGG run on the SAME input - the p-values
and FDRs must be byte-identical (run_pathway_topology calls the exact same
building blocks internally), so any drift between the two paths fails loudly
here rather than silently producing two different "correct" answers.
"""
import os
import re
import time

import requests

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
HEADERS = {"X-Session-Id": "metab-live-e2e"}
PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def post(path, **kwargs):
    return requests.post(BASE + path, headers=HEADERS, timeout=600, **kwargs)


def get(path, **kwargs):
    return requests.get(BASE + path, headers=HEADERS, timeout=600, **kwargs)


def wait_job(job_id, timeout=1800):
    start = time.time()
    while time.time() - start < timeout:
        job = get(f"/api/jobs/{job_id}").json()
        if job.get("status") == "done":
            return job.get("result")
        if job.get("status") == "error":
            raise RuntimeError(job.get("error") or job)
        time.sleep(2)
    raise TimeoutError(job_id)


def example(name):
    return open(os.path.join(_ROOT, "Data3", name), "rb").read()


def run_tile(kind, config):
    r = post(f"/api/metabolomics/{kind}/configure", json=config)
    check(r.status_code == 200, f"[{kind}] configure accepted ({r.status_code} {r.text[:120] if r.status_code != 200 else ''})")
    r = post(f"/api/metabolomics/{kind}/run")
    check(r.status_code == 200, f"[{kind}] run started")
    return wait_job(r.json()["job_id"])


print("=== input-contract errors are explained, not crashes ===")
r = post("/api/metabolomics/list/upload", files={"file": ("bad.csv", b"Compound_ID,Value\nx,1\n")})
check(r.status_code == 400 and "Name" in r.text, f"a list with no Name column is rejected naming the column (got {r.status_code}: {r.text[:100]})")
one_group = b"Sample,Group,A,B\ns1,X,1,2\ns2,X,2,3\ns3,X,3,4\n"
r = post("/api/metabolomics/quant/upload", files={"file": ("one.csv", one_group)})
check(r.status_code == 400 and ("two" in r.text.lower() or "group" in r.text.lower()),
      f"a concentration table with a single group is rejected (got {r.status_code}: {r.text[:120]})")
r = post("/api/metabolomics/quant/upload?orientation=sideways", files={"file": ("x.csv", one_group)})
check(r.status_code == 400, "an unknown orientation override is rejected")

print("=== list path: upload + resolve (shared by three tiles) ===")
r = post("/api/metabolomics/list/upload", files={"file": ("list.csv", example("example_metabolomics_list.csv"))})
check(r.status_code == 200 and r.json()["n_rows"] >= 50, f"list uploaded ({r.json().get('n_rows')} compounds)")
check(not any(c.lower() in ("log2fc", "pvalue", "padj") for c in r.json()["columns"]),
      "the example list carries no statistics columns")

print("=== a repeated compound name in the uploaded list is removed, once, and reported ===")
example_bytes = example("example_metabolomics_list.csv")
# Duplicate the first data row using its OWN original terminator bytes
# (splitlines(keepends=True) keeps them) rather than hard-coding one, so this
# does not depend on whether the example file is LF- or CRLF-terminated.
data_lines = example_bytes.splitlines(keepends=True)
dup_bytes = example_bytes + data_lines[1]
r_dup = post("/api/metabolomics/list/upload", files={"file": ("list_dup.csv", dup_bytes)})
check(r_dup.status_code == 200 and r_dup.json()["n_rows"] >= 61,
      f"the upload itself keeps the duplicate row (it is not de-duplicated until resolve) (got {r_dup.json().get('n_rows')})")
resolve_dup = wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
meta_dup = resolve_dup["metrics"]
check(meta_dup.get("n_duplicate_names_removed") == 1,
      f"the resolve step reports exactly one duplicate removed (got {meta_dup.get('n_duplicate_names_removed')})")
check(meta_dup["n_input_rows"] == r_dup.json()["n_rows"],
      "n_input_rows still reports the ORIGINAL upload count, including the duplicate")
check(meta_dup["n_input"] == meta_dup["n_input_rows"] - 1,
      "n_input (compounds actually resolved) is exactly one less than the upload")

print("=== list path: upload + resolve (shared by three tiles) ===")
resolve = wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
meta = resolve["metrics"]
check(meta["n_matched"] >= 50, f"most compounds matched ({meta['n_matched']} of {meta['n_input']})")
check(sum(meta["n_matched_by_tier"].values()) == meta["n_matched"], "tier counts partition the matched rows")

print("=== Over-Representation Analysis ===")
ora = run_tile("ora", {"databases": ["KEGG", "SMPDB"], "min_term_size": 3})
check(ora["summary"]["universe_scope"] == "database", "list tiles test against the database universe")
kegg = ora["tables"]["KEGG"]
top = kegg["table"]["rows"][0]
check(float(top["FDR"]) < 0.01, f"KEGG's top term is strongly significant ({top['Term Name']!r}, FDR={top['FDR']:.2g})")
# Organism ids are <code><5 digits> and the code is not always 3 characters, so
# the map number is taken from the tail. (The old check stripped "map", which
# would have become a silent no-op on "hsa..." ids and could never fail again.)
num = re.search(r"(\d{5})$", top["Term ID"]).group(1)
check("overview" not in top["Term Name"].lower() and not num.startswith("010"),
      f"the top term is a real pathway, not a KEGG overview map (got {top['Term ID']})")
check(top["Term ID"].startswith("hsa"), f"the term id names the organism (got {top['Term ID']})")
check(all(r["Term ID"].startswith("hsa") for r in kegg["table"]["rows"]), "EVERY KEGG term is a human map")
check(kegg["universe_size"] > 1000,
      f"the universe is the organism's KEGG reference set, not the 60-compound list (N={kegg['universe_size']})")
check("Homo sapiens (hsa)" in (kegg.get("caveat") or ""), "the KEGG caveat names the organism used")
human_universe = kegg["universe_size"]
smpdb = ora["tables"]["SMPDB"]
check(smpdb["n_terms_tested"] > 0, f"SMPDB tested real terms ({smpdb['n_terms_tested']})")
check("Overlapping Compounds" in top and "KEY" not in top["Overlapping Compounds"] and top["Overlapping Compounds"],
      f"members are the user's own names, not internal keys (got {top['Overlapping Compounds'][:60]!r})")
dl = get("/api/metabolomics/ora/download", params={"key": "KEGG"})
check(dl.ok and dl.text.splitlines()[0].count("Overlapping Compounds") == 1, "result CSV downloads")

print("=== KEGG organism resolution routes ===")
r = post("/api/metabolomics/kegg/resolve-organism")
check(r.status_code == 200 and r.json()["organism_code"] == "hsa",
      f"the derived KEGG organism for human is hsa (got {r.status_code} {r.text[:80]})")

print("=== KEGG is organism-specific: a second organism gets a different universe from the same link file ===")
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "4932", "taxon_name": "Saccharomyces cerevisiae"}).json()["job_id"])
yeast = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3})["tables"]["KEGG"]
check(yeast["n_terms_tested"] > 0 and all(r["Term ID"].startswith("sce") for r in yeast["table"]["rows"]),
      f"yeast terms are sce-prefixed ({yeast['n_terms_tested']} tested)")
check(0 < yeast["universe_size"] < human_universe,
      f"yeast has a strictly smaller KEGG universe than human (yeast {yeast['universe_size']} < human {human_universe}) - "
      f"fails if KEGG is ever made organism-blind again")

print("=== KEGG: an organism KEGG has no genome for is EMPTY with a note, never a fallback ===")
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "0", "taxon_name": "Nonexistus fakeus"}).json()["job_id"])
none_block = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3})["tables"]["KEGG"]
check(none_block["n_terms_tested"] == 0 and not none_block["table"]["rows"] and none_block["universe_size"] == 0,
      "no KEGG genome -> nothing tested (no silent fallback to the generic reference maps)")
check("no genome" in (none_block.get("note") or "").lower(), f"...and it says why (got {none_block.get('note')!r:.120})")

print("=== KEGG: an ambiguous organism is empty until a KEGG organism is picked ===")
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "1423", "taxon_name": "Bacillus subtilis"}).json()["job_id"])
amb = post("/api/metabolomics/kegg/resolve-organism").json()
check(amb["organism_code"] == "" and len(amb["candidates"]) > 1,
      f"B. subtilis is genuinely ambiguous ({len(amb['candidates'])} candidates)")
unpicked = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3})["tables"]["KEGG"]
check(unpicked["n_terms_tested"] == 0 and (unpicked.get("note") or ""), "nothing picked -> KEGG empty with a note")

print("=== KEGG organism box: check-organism validates what the user picked ===")
good = post("/api/metabolomics/kegg/check-organism", json={"organism_code": "bsu"})
check(good.status_code == 200 and good.json()["organism_code"] == "bsu" and good.json()["n_maps"] > 50
      and "subtilis" in good.json()["organism_name"],
      f"a real KEGG organism returns its name and map count (got {good.status_code} {good.text[:120]})")
check(post("/api/metabolomics/kegg/check-organism", json={"organism_code": "zzzz"}).status_code == 400,
      "a code KEGG has no maps for is rejected, not accepted and left to return an empty table")
check(post("/api/metabolomics/kegg/check-organism", json={"organism_code": "../etc"}).status_code == 400,
      "a malformed code is rejected before it can reach a URL or cache path")
check(post("/api/metabolomics/ora/configure", json={"databases": ["KEGG"], "kegg_organism_code": "zzzz"}).status_code == 400,
      "configure re-validates the explicit organism (the run cannot be started with a bad one)")

print("=== KEGG: an explicit pick unblocks the ambiguous organism and wins over the derived one ===")
picked = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3, "kegg_organism_code": "bsu"})["tables"]["KEGG"]
check(picked["n_terms_tested"] > 0 and all(r["Term ID"].startswith("bsu") for r in picked["table"]["rows"]),
      f"picking bsu makes KEGG test bsu maps ({picked['n_terms_tested']} terms)")
check("chosen explicitly" in (picked.get("caveat") or ""), "the caveat says the organism was chosen explicitly")
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
cross = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3, "kegg_organism_code": "sce"})["tables"]["KEGG"]
check(cross["n_terms_tested"] > 0 and all(r["Term ID"].startswith("sce") for r in cross["table"]["rows"]),
      "with human chosen above, an explicit yeast pick still tests yeast maps (the explicit choice wins)")
check(0 < cross["universe_size"] < human_universe, "...with yeast's own, smaller universe")

# Back to human for the remaining tiles (a new resolve also clears the confirmation).
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
r = post("/api/metabolomics/kegg/resolve-organism")
check(r.json()["organism_code"] == "hsa", "back on human")

print("=== Chemical Class Enrichment ===")
cc = run_tile("chemclass", {"min_term_size": 3, "chemical_class_level": "sub_class"})
rows = cc["tables"]["Chemical Class"]["table"]["rows"]
check(rows and "TCA" in rows[0]["Term Name"], f"RefMet recovers the planted TCA signal (top={rows[0]['Term Name'] if rows else None})")
check(float(rows[0]["FDR"]) < 0.01, f"...significantly (FDR={rows[0]['FDR']:.2g})")

print("=== Lipid Class / Structural Enrichment ===")
lc = run_tile("lipidclass", {"min_term_size": 3, "lipid_class_level": "lm_main_class"})
lrows = lc["tables"]["Lipid Class"]["table"]["rows"]
check(lc["tables"]["Lipid Class"]["n_terms_tested"] > 0, "LIPID MAPS classes tested")
check(any("Glycerophosphocholines" in r["Term Name"] for r in lrows),
      "the example's phosphatidylcholine shorthand names (PC 34:1, ...) reach their class even though "
      "no single LIPID MAPS structure carries those names")

print("=== cross-tile staleness: resolving Quant for a DIFFERENT organism clears ORA's cached results ===")
# The Chemical Class / Lipid Class runs just above each did their own
# list/resolve for human (list resolve unconditionally clears ORA's cached
# results, not only on an organism change), so ORA has to be re-run here to
# get a genuine "downloadable before anything else changes" baseline.
run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3})
dl_before = get("/api/metabolomics/ora/download", params={"key": "KEGG"})
check(dl_before.ok, "ORA's results are still downloadable before anything else changes the organism")
r = post("/api/metabolomics/quant/upload", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
check(r.status_code == 200, "concentration table re-uploaded")
wait_job(post("/api/metabolomics/quant/resolve", json={"taxon_id": "10090", "taxon_name": "Mus musculus"}).json()["job_id"])
dl_after = get("/api/metabolomics/ora/download", params={"key": "KEGG"})
check(dl_after.status_code == 400,
      f"SENSITIVITY: after resolving Quant for a NEW organism (mouse), ORA's now-stale human results are no "
      f"longer downloadable - it must be re-run for the current organism, not silently served from before "
      f"(got {dl_after.status_code})")
# The Quantitative Enrichment section right below re-uploads and re-resolves
# the concentration table for human on its own, so no restore is needed here.

print("=== Quantitative Enrichment (concentration table) ===")
r = post("/api/metabolomics/quant/upload", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
info = r.json()
check(r.status_code == 200 and info["orientation"] == "standard" and info["group_column"] == "Group",
      f"layout and group column auto-detected (got {info.get('orientation')!r}, {info.get('group_column')!r})")
check(sorted(info["group_levels"]) == ["Control", "Treated"], "two groups found")
r = post("/api/metabolomics/quant/upload?group_column=Nope", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
check(r.status_code == 400 and "Nope" in r.text,
      f"an unknown group-column override is rejected, not silently replaced by the auto-detected one ({r.status_code}: {r.text[:80]})")
r = post("/api/metabolomics/quant/upload?group_column=Group", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
check(r.status_code == 200 and r.json()["group_column"] == "Group", "a valid group-column override is accepted")
r = post("/api/metabolomics/quant/upload", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
qres = wait_job(post("/api/metabolomics/quant/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
check(qres["metrics"]["n_matched"] >= 35, f"compound headers resolved ({qres['metrics']['n_matched']} matched)")

bad = post("/api/metabolomics/quant/configure", json={"group_a": "Treated", "group_b": "Treated", "databases": ["KEGG"]})
check(bad.status_code == 400, "identical groups are rejected")
cfg = post("/api/metabolomics/quant/configure", json={
    "group_a": "Treated", "group_b": "Control", "log_transform": True, "min_per_group": 3,
    "metric": "log2fc", "databases": ["KEGG", "Chemical Class", "Lipid Class"], "min_term_size": 3})
check(cfg.status_code == 200, f"configure computes the statistics ({cfg.status_code})")
sm = cfg.json()["stats_meta"]
check(sm["n_compounds_tested"] > 30 and sm["group_a"] == "Treated", f"statistics computed by Omiscera (tested {sm['n_compounds_tested']})")
stats = get("/api/metabolomics/quant/download/statistics.csv")
check(stats.ok and stats.text.splitlines()[0].startswith("Name,n_a,n_b,mean_a,mean_b,log2fc"), "statistics CSV is downloadable")
result = wait_job(post("/api/metabolomics/quant/run").json()["job_id"])
check(result["summary"]["universe_scope"] == "experiment", "the concentration tile keeps the measured universe")
check(len(result["compound_log2fc"]) > 30, "per-compound log2FC is returned for network colouring")
qk = result["tables"]["KEGG"]["table"]["rows"]
check(qk and float(qk[0]["Effect (mean difference)"]) > 0, "KEGG's top term is shifted UP in Treated (planted TCA/glycolysis increase)")
ql = result["tables"]["Lipid Class"]["table"]["rows"]
check(ql and min(float(r_["Effect (mean difference)"]) for r_ in ql) < 0,
      "a Lipid Class term is shifted DOWN in Treated (planted lipid decrease) - the test is two-sided")

print("=== Pathway Topology Analysis: list mode (Fisher, identical p-values to ORA's own KEGG) ===")
r = post("/api/metabolomics/list/upload", files={"file": ("list.csv", example("example_metabolomics_list.csv"))})
check(r.status_code == 200, "list re-uploaded")
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
topo_list = run_tile("topology", {"input_mode": "list", "measure": "betweenness", "min_term_size": 3, "max_term_size": 500})
tblock = topo_list["tables"]["KEGG"]
check(tblock["n_terms_tested"] > 0 and tblock.get("measure") == "betweenness",
      f"KEGG tested ({tblock['n_terms_tested']} terms), measure recorded")
trows = tblock["table"]["rows"]
check(all("Impact" in row and "Matched in graph" in row and "Graph nodes" in row for row in trows[:5]),
      "every row carries the topology columns")
check(all(0.0 <= float(row["Impact"]) <= 1.0 for row in trows), "Impact is always in [0, 1]")
check(set(tblock.get("graphs", {})) == {row["Term ID"] for row in trows},
      "a graph is present for every TESTED pathway, and only those")
top_graph = tblock["graphs"][trows[0]["Term ID"]]
check(top_graph["status"] in ("ok", "nodata", "error") and "nodes" in top_graph,
      f"the top pathway's graph carries a real KEGG fetch status (got {top_graph['status']})")

ora_again = run_tile("ora", {"databases": ["KEGG"], "min_term_size": 3, "max_term_size": 500})
ora_by_id = {row["Term ID"]: row for row in ora_again["tables"]["KEGG"]["table"]["rows"]}
mismatches = [row["Term ID"] for row in trows if row["Term ID"] in ora_by_id
              and (row["p-value"] != ora_by_id[row["Term ID"]]["p-value"] or row["FDR"] != ora_by_id[row["Term ID"]]["FDR"])]
check(not mismatches,
      f"SENSITIVITY: p-values/FDRs are byte-identical to ORA's own KEGG block for the same input/config "
      f"(mismatches: {mismatches})")

dl = get("/api/metabolomics/topology/download", params={"key": "KEGG"})
check(dl.ok and dl.text.splitlines()[0].startswith("Term ID,Term Name"), "result table is downloadable")

print("=== Pathway Topology Analysis: concentration mode (Wilcoxon, identical p-values to Quant's own KEGG) ===")
r = post("/api/metabolomics/quant/upload", files={"file": ("conc.csv", example("example_metabolomics_concentrations.csv"))})
check(r.status_code == 200, "concentration table re-uploaded")
wait_job(post("/api/metabolomics/quant/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
topo_conc = run_tile("topology", {
    "input_mode": "conc", "measure": "outdegree", "min_term_size": 3, "max_term_size": 500,
    "group_a": "Treated", "group_b": "Control", "log_transform": True, "min_per_group": 3, "metric": "log2fc",
})
cblock = topo_conc["tables"]["KEGG"]
check(cblock["n_terms_tested"] > 0 and cblock.get("measure") == "outdegree",
      f"KEGG tested in concentration mode ({cblock['n_terms_tested']} terms), own measure recorded")
check("Effect (mean difference)" in cblock["table"]["rows"][0], "concentration mode uses the Wilcoxon/burden columns")

quant_kegg = result["tables"]["KEGG"]["table"]["rows"]
quant_by_id = {row["Term ID"]: row for row in quant_kegg}
crows = cblock["table"]["rows"]
mismatches_c = [row["Term ID"] for row in crows if row["Term ID"] in quant_by_id
                and (row["p-value"] != quant_by_id[row["Term ID"]]["p-value"] or row["FDR"] != quant_by_id[row["Term ID"]]["FDR"])]
check(not mismatches_c,
      f"SENSITIVITY: concentration-mode p-values/FDRs are byte-identical to Quantitative Enrichment's own KEGG "
      f"block for the same input/config (mismatches: {mismatches_c})")

stats_dl = get("/api/metabolomics/topology/download/statistics.csv")
check(stats_dl.ok, "topology's own statistics CSV (concentration mode) is downloadable")

print("=== Pathway Topology Analysis: organism with no KEGG genome -> empty with a note ===")
r = post("/api/metabolomics/list/upload", files={"file": ("list.csv", example("example_metabolomics_list.csv"))})
wait_job(post("/api/metabolomics/list/resolve",
              json={"taxon_id": "0", "taxon_name": "Nonexistus fakeus"}).json()["job_id"])
topo_none = run_tile("topology", {"input_mode": "list", "min_term_size": 3})
none_block = topo_none["tables"]["KEGG"]
check(none_block["n_terms_tested"] == 0 and none_block.get("note"),
      f"no KEGG genome -> nothing tested, with a note (got {none_block.get('note')!r})")
# Restore human for anything that runs after this in a future edit.
wait_job(post("/api/metabolomics/list/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
