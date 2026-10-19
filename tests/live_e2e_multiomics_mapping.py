"""Live end-to-end test of Multi-omics Cross-omics Pathway Mapping through the HTTP API
(no browser - see playwright_multiomics_mapping.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network access: select all
four omics -> upload the four bundled examples -> resolve each omics with its OWN resolver ->
configure -> list the pathways -> map KEGG hsa00020 and Reactome R-HSA-71403 / R-HSA-70171.

Every claim about the biology is checked against an INDEPENDENT recomputation from the raw
files: KEGG's own KGML and gene list are fetched here and the expected Genomics hits per map are
rebuilt without any Omiscera code; the differential-expression numbers shown on a node are
compared with the uploaded CSV row for that gene. The bundled Genomics / Metabolomics examples
plant TCA-cycle + glycolysis genes and intermediates, so the TCA cycle must come out as one of
the best-covered pathways and citrate synthase (CS) must be linked to citrate.
"""
import json
import os
import re
import time

import pandas as pd
import requests

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
_RUN = str(int(time.time() * 1000))
HEADERS = {"X-Session-Id": "multiomics-map-live-" + _RUN}
PASS = FAIL = 0
OMICS = ["genomics", "transcriptomics", "proteomics", "metabolomics"]


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def post(path, **kwargs):
    return requests.post(BASE + path, headers=HEADERS, timeout=900, **kwargs)


def get(path, **kwargs):
    return requests.get(BASE + path, headers=HEADERS, timeout=900, **kwargs)


def wait_job(job_id, timeout=2400):
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


ABSENT = {"state": "absent", "direction": None, "members": []}


def S(node, omics):
    """A box's status for one omics. The payload carries only the omics that touched the box (a global map has ~5,000
    boxes), so a missing key IS 'absent'."""
    return node["status"].get(omics, ABSENT)


def node_by_label(result, label, kind=None):
    return [n for n in result["nodes"] if (n["label"] == label) and (kind is None or n["kind"] == kind)]


def hit_labels(result, omics, kind="gene"):
    out = set()
    for n in result["nodes"]:
        if n["kind"] != kind:
            continue
        st = n["status"].get(omics)
        if st and st["state"] == "hit":
            out |= {m["label"] for m in st["members"] if m["hit"]}
    return out


# ---------------------------------------------------------------------------
print("=== guard rails before anything is set up ===")
r = post("/api/multiomics/mapping/configure", json={})
check(r.status_code == 400, f"configuring before the omics are chosen/resolved is refused ({r.status_code})")
r = post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa00020"})
check(r.status_code == 400, "mapping before configuring is refused")

print("=== set up four omics ===")
r = post("/api/multiomics/select", json={"omics_types": OMICS})
check(r.status_code == 200 and r.json()["selected"] == OMICS, "all four omics selected")
files = {"genomics": "example_multiomics_genes.csv", "transcriptomics": "example_transcriptomics_de_full.csv",
         "proteomics": "example_proteomics_de_full.csv", "metabolomics": "example_metabolomics_list.csv"}
for o, f in files.items():
    r = post(f"/api/multiomics/{o}/upload", files={"file": (f, example(f))})
    check(r.status_code == 200, f"{o} example uploaded")
res = wait_job(post("/api/multiomics/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
check(all(res["omics"][o]["n_usable"] > 0 for o in OMICS), "every omics resolved with usable rows "
      + ", ".join(f"{o}={res['omics'][o]['n_usable']}" for o in OMICS))

r = post("/api/multiomics/mapping/configure", json={"de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0},
                                                                    "proteomics": {"padj": 0.05, "lfc": 1.0}}})
check(r.status_code == 200 and r.json()["config"]["de_cutoffs"]["transcriptomics"]["padj"] == 0.05, "mapping configured")
r = post("/api/multiomics/mapping/configure", json={"de_cutoffs": {"transcriptomics": {"padj": 5, "lfc": 1}}})
check(r.status_code == 400, "a padj cutoff above 1 is refused")
r = post("/api/multiomics/mapping/configure", json={"de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0},
                                                                    "proteomics": {"padj": 0.05, "lfc": 1.0}}})

r = post("/api/multiomics/mapping/pathways", json={"database": "MetaCyc"})
check(r.status_code == 400, "MetaCyc is not offered (needs a paid BioCyc subscription)")
r = post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "../../etc/passwd"})
check(r.status_code == 400, "a malformed pathway id is refused before it can reach a URL")
r = post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "R-HSA-71403"})
check(r.status_code == 400, "a Reactome id is not accepted for KEGG")

# ---------------------------------------------------------------------------
print("=== pathway pickers ===")
kegg_pick = wait_job(post("/api/multiomics/mapping/pathways", json={"database": "KEGG"}).json()["job_id"])
reactome_pick = wait_job(post("/api/multiomics/mapping/pathways", json={"database": "Reactome"}).json()["job_id"])
for name, pick in (("KEGG", kegg_pick), ("Reactome", reactome_pick)):
    rows = pick["table"]["rows"]
    check(pick["n_pathways"] == len(rows) and rows, f"[{name}] picker lists {len(rows)} pathways hit by at least one omics")
    check(all(r["Total hits"] > 0 for r in rows), f"[{name}] every listed pathway has at least one hit")
    key = [(-r["Omics with a hit"], -r["Total hits"], r["Term Name"]) for r in rows]
    check(key == sorted(key), f"[{name}] ranked by omics-with-a-hit, then total hits")
    check(pick["table"]["columns"][:4] == ["Term ID", "Term Name", "Omics with a hit", "Total hits"]
          and all(f"{o.capitalize()}: hit / annotated" in pick["table"]["columns"] for o in OMICS),
          f"[{name}] columns include every omics' hit / annotated count")
def k_of(r, omics):
    return int(r[f"{omics.capitalize()}: hit / annotated"].split(" / ")[0])


tca = next(r for r in kegg_pick["table"]["rows"] if r["Term ID"] == "hsa00020")
rank = kegg_pick["table"]["rows"].index(tca) + 1
check(k_of(tca, "genomics") >= 8 and k_of(tca, "metabolomics") >= 5 and tca["Omics with a hit"] >= 2,
      f"BIOLOGY: the planted TCA genes and intermediates hit KEGG's TCA cycle (Genomics {k_of(tca, 'genomics')}, "
      f"Metabolomics {k_of(tca, 'metabolomics')}; hit by {tca['Omics with a hit']} omics, rank {rank} of {len(kegg_pick['table']['rows'])} - "
      f"ranked by omics count first, and the synthetic DE examples do not plant TCA signal)")
check(all(r["Omics with a hit"] >= tca["Omics with a hit"] for r in kegg_pick["table"]["rows"][:rank - 1]),
      "every pathway ranked above it is hit by at least as many omics")
rtca = next(r for r in reactome_pick["table"]["rows"] if r["Term ID"] == "R-HSA-71403")
check(k_of(rtca, "genomics") >= 5 and k_of(rtca, "metabolomics") >= 1 and rtca["Omics with a hit"] >= 2,
      f"BIOLOGY: Reactome's citric acid cycle is hit by Genomics ({k_of(rtca, 'genomics')}) and Metabolomics "
      f"({k_of(rtca, 'metabolomics')}) - {rtca['Omics with a hit']} omics")
check(all(r["Joint combined FDR"] is None for r in kegg_pick["table"]["rows"]) and not kegg_pick["has_joint_fdr"],
      "no joint FDR is shown before a joint run exists")

# ---------------------------------------------------------------------------
print("=== KEGG hsa00020 (citrate cycle) ===")
kegg = wait_job(post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa00020"}).json()["job_id"])
check(kegg["status"] == "ok" and kegg["omics"] == OMICS and kegg["pathway_name"].startswith("Citrate cycle"),
      f"mapped: {kegg['pathway_name']}")
check(kegg["summary"]["n_nodes"] == len(kegg["nodes"]) and len(kegg["reactions"]) == 22,
      f"{len(kegg['nodes'])} boxes, {len(kegg['reactions'])} reactions (KEGG's own map has 22)")

raw_kgml = requests.get("https://rest.kegg.jp/get/hsa00020/kgml", timeout=60).text
gene_ids_in_map = set(re.findall(r'name="(hsa:\d+(?: hsa:\d+)*)" type="gene"', raw_kgml))
gene_ids_in_map = {g for entry in gene_ids_in_map for g in entry.split()}
listing = requests.get("https://rest.kegg.jp/list/hsa", timeout=120).text
symbols_of = {}
for line in listing.splitlines():
    parts = line.split("\t")
    if len(parts) >= 4 and parts[0] in gene_ids_in_map:
        symbols_of[parts[0]] = {s.strip().upper() for s in parts[3].split(";")[0].split(",")}
genomics_symbols = {s.strip().upper() for s in pd.read_csv(os.path.join(_ROOT, "Data3", files["genomics"]))["Gene Symbol"]}
expected = {s for syms in symbols_of.values() for s in syms} & genomics_symbols
got = {s.upper() for s in hit_labels(kegg, "genomics")}
check(len(expected) >= 8, f"INDEPENDENT: {len(expected)} of the 29 Genomics genes sit on KEGG's own TCA map")
check(got == expected, f"INDEPENDENT: the Genomics hits on the map equal the genes found by re-reading KEGG's KGML + gene list "
      f"(got {len(got)}, expected {len(expected)}; missing {sorted(expected - got)[:5]}, extra {sorted(got - expected)[:5]})")

mets = hit_labels(kegg, "metabolomics", "compound")
check({"Citrate", "Succinate", "Fumarate"} <= {m.split(";")[0].strip().capitalize() for m in mets} or len(mets) >= 6,
      f"BIOLOGY: planted TCA intermediates are hit compounds on the map ({sorted(mets)[:8]})")
check(all(any(n["label"] for n in node_by_label(kegg, n["label"])) for n in kegg["nodes"] if n["kind"] == "compound"),
      "compound boxes are labelled (KEGG's name, not a bare C-number, where KEGG has one)")
labels = [n["label"] for n in kegg["nodes"] if n["kind"] == "compound"]
check(sum(1 for l in labels if re.fullmatch(r"C\d{5}", l)) <= len(labels) // 4,
      f"most compound boxes carry a real name ({sum(1 for l in labels if re.fullmatch(r'C' + chr(92) + 'd{5}', l))} bare C-numbers of {len(labels)})")

# three states + direction, checked against the uploaded DE table
t_csv = pd.read_csv(os.path.join(_ROOT, "Data3", files["transcriptomics"])).set_index("Gene Symbol")
checked_dir = bad_dir = bad_num = n_measured = 0
for n in kegg["nodes"]:
    st = S(n, "transcriptomics")
    if st["state"] == "measured":
        n_measured += 1
    for m in st["members"]:
        row = t_csv.loc[m["label"]] if m["label"] in t_csv.index else None
        if row is None:
            continue
        row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
        if m["log2fc"] is None or abs(m["log2fc"] - float(row["log2FoldChange"])) > 1e-6:
            bad_num += 1
        is_hit = float(row["padj"]) <= 0.05 and abs(float(row["log2FoldChange"])) >= 1.0
        if is_hit != m["hit"]:
            bad_num += 1
        if m["hit"]:
            checked_dir += 1
    if st["state"] == "hit" and st["direction"] in ("up", "down"):
        hits = [m for m in st["members"] if m["hit"]]
        if not all((m["log2fc"] > 0) == (st["direction"] == "up") for m in hits) and len({m["log2fc"] > 0 for m in hits}) == 1:
            bad_dir += 1
check(bad_num == 0, f"INDEPENDENT: every Transcriptomics member's log2FC equals the uploaded row and 'hit' equals padj <= 0.05 & |log2FC| >= 1 "
      f"({bad_num} mismatches over the members shown)")
check(bad_dir == 0, "a node's direction agrees with the sign of its hit members")
check(n_measured >= 1, f"SENSITIVITY: measured-but-not-significant boxes exist ({n_measured}) - distinct from absent")

links = kegg["links"]["rows"]
cs = [r for r in links if "CS" in [g.strip() for g in r["Enzyme (hit)"].split(",")]]
check(cs and any("itrate" in r["Metabolites hit"] for r in cs),
      f"BIOLOGY: citrate synthase (CS) is linked to citrate ({[(r['Reaction'], r['Metabolites hit'][:40]) for r in cs][:2]})")
check(all(r["Evidence"] for r in links) and len(links) == kegg["summary"]["n_links"], f"{len(links)} cross-omics reaction links, each with its evidence")
check(all(("enzyme" in r["Evidence"]) for r in links), "every link names its enzyme evidence")
mem = kegg["members"]["rows"]
check(all(m["Matched by"] for m in mem) and {"hit", "measured, not significant"} >= {m["Status"] for m in mem},
      "every member row records its identifier bridge and a status of hit / measured, not significant")
check(any(m["Matched by"] == "UniProt -> KEGG gene" for m in mem if m["Omics"] == "Proteomics"), "Proteomics reaches KEGG genes through UniProt")
check(any(m["Matched by"].startswith("Entrez") for m in mem if m["Omics"] == "Transcriptomics"), "Transcriptomics reaches KEGG genes through Entrez")

dl = get("/api/multiomics/mapping/download?table=members")
dfm = pd.read_csv(__import__("io").StringIO(dl.text))
check(dl.status_code == 200 and len(dfm) == len(mem) and list(dfm.columns) == kegg["members"]["columns"], "members CSV = the payload, same columns")
dl2 = get("/api/multiomics/mapping/download?table=links")
check(dl2.status_code == 200 and len(pd.read_csv(__import__("io").StringIO(dl2.text))) == len(links), "links CSV = the payload")

# ---------------------------------------------------------------------------
print("=== cutoffs change the hits ===")
post("/api/multiomics/mapping/configure", json={"de_cutoffs": {"transcriptomics": {"padj": 1e-12, "lfc": 3.0},
                                                               "proteomics": {"padj": 1e-12, "lfc": 3.0}}})
strict = wait_job(post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa00020"}).json()["job_id"])
n_loose = sum(1 for n in kegg["nodes"] if S(n, "transcriptomics")["state"] == "hit")
n_strict = sum(1 for n in strict["nodes"] if S(n, "transcriptomics")["state"] == "hit")
check(n_strict <= n_loose and strict["summary"]["per_omics"]["genomics"]["nodes_hit"] == kegg["summary"]["per_omics"]["genomics"]["nodes_hit"],
      f"a stricter cutoff never adds Transcriptomics hits ({n_loose} -> {n_strict}) and does not touch the list omics")
n_meas_strict = sum(1 for n in strict["nodes"] if S(n, "transcriptomics")["state"] == "measured")
check(n_meas_strict >= n_measured, f"...the genes it drops become 'measured', not 'absent' ({n_measured} -> {n_meas_strict})")
post("/api/multiomics/mapping/configure", json={"de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0},
                                                               "proteomics": {"padj": 0.05, "lfc": 1.0}}})

# ---------------------------------------------------------------------------
print("=== Reactome R-HSA-71403 (citric acid cycle) ===")
rx = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": "R-HSA-71403"}).json()["job_id"])
check(rx["status"] == "ok" and rx["restricted_to"] is None and "citric acid cycle" in rx["pathway_name"].lower(),
      f"mapped from its own diagram: {rx['pathway_name']}")
kinds = {k: sum(1 for n in rx["nodes"] if n["kind"] == k) for k in ("gene", "compound")}
check(kinds["gene"] >= 10 and kinds["compound"] >= 20, f"diagram boxes: {kinds}")
check({"CS", "ACO2"} & hit_labels(rx, "genomics") or {"CS", "ACO2"} & {m["label"] for n in rx["nodes"] for m in S(n, "genomics")["members"]},
      "BIOLOGY: Genomics TCA genes (CS / ACO2) are found on the Reactome diagram")
check(any(S(n, "metabolomics")["state"] == "hit" for n in rx["nodes"] if n["kind"] == "compound"),
      "BIOLOGY: Metabolomics compounds are matched to Reactome chemicals by ChEBI")
cof = [n for n in rx["nodes"] if n["kind"] == "gene" and S(n, "metabolomics")["state"] != "absent"]
check(not cof, "SENSITIVITY: no protein / complex box is ever a Metabolomics hit (cofactors inside complexes are ignored)")
check(any(r["Metabolites hit"] and r["Enzyme (hit)"] for r in rx["links"]["rows"]), f"{len(rx['links']['rows'])} Reactome reaction links with an enzyme and a metabolite")
check(any(m["Matched by"] == "UniProt" for m in rx["members"]["rows"]), "Reactome matches Genomics / Proteomics on UniProt")
check(rx["caveat"] and "computes no p-value" in rx["caveat"], "the descriptive caveat is part of the result")

print("=== Reactome sub-pathway without its own diagram (Glycolysis, R-HSA-70171) ===")
gly = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": "R-HSA-70171"}).json()["job_id"])
check(gly["status"] == "ok" and gly["restricted_to"] == "Glucose metabolism (R-HSA-70326)" and gly["pathway_name"] == "Glycolysis",
      f"drawn inside its nearest ancestor with a diagram: {gly['restricted_to']}")
check(len(gly["reactions"]) < 40 and len(gly["nodes"]) < 120,
      f"restricted to glycolysis' own reactions ({len(gly['reactions'])} reactions, {len(gly['nodes'])} boxes), not the whole parent diagram")
check(any(S(n, "genomics")["state"] == "hit" for n in gly["nodes"]), "BIOLOGY: planted glycolysis genes are hit on the glycolysis map")

print("=== KEGG's global overview map hsa01100: genes are LINES ===")
import xml.etree.ElementTree as _ET

raw01100 = requests.get("https://rest.kegg.jp/get/hsa01100/kgml", timeout=120).text
root01100 = _ET.fromstring(raw01100)
entries = root01100.findall("entry")
gene_lines = [e for e in entries if e.get("type") == "gene"]
expected_nodes = sum(1 for e in entries if e.get("type") in ("gene", "compound", "map")
                     and not ((e.find("graphics").get("name") or "").upper().startswith("TITLE:")))
g01 = wait_job(post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa01100"}).json()["job_id"])
check(g01["status"] == "ok" and g01["line_graphics"] is True and len(g01["nodes"]) == expected_nodes,
      f"hsa01100 is mapped as a line map: {len(g01['nodes'])} boxes = the {expected_nodes} gene / compound / link entries of KEGG's own file")
gnodes = [n for n in g01["nodes"] if n["kind"] == "gene"]
check(len(gnodes) == len(gene_lines) == 1342 and all("path" in n and len(n["path"]) >= 4 for n in gnodes),
      f"SENSITIVITY: all {len(gnodes)} genes carry a polyline of at least two points (before: none, so every gene sat at (0, 0))")
check(sum(1 for n in gnodes if n["x"] == 0 and n["y"] == 0) == 0 and len({(n["x"], n["y"]) for n in gnodes}) > 1000,
      f"SENSITIVITY: no gene at (0, 0); {len({(n['x'], n['y']) for n in gnodes})} different positions - the pile in the top-left corner is gone")
raw_first = {e.get("id"): [float(v) for v in e.find("graphics").get("coords").split(",")[:2]] for e in gene_lines}
n_first_ok = sum(1 for n in gnodes if [float(n["path"][0]), float(n["path"][1])] == raw_first[n["id"][1:]])
check(n_first_ok == len(gnodes), f"INDEPENDENT: every polyline starts where KEGG's own coords start ({n_first_ok} of {len(gnodes)})")
check(len(g01["decor_lines"]) == sum(1 for e in entries if e.get("type") == "ortholog") == 2536 and len(g01["reactions"]) == len(root01100.findall("reaction")),
      f"the {len(g01['decor_lines'])} ortholog lines are decoration and all {len(g01['reactions'])} reactions are kept")
check(g01["relations"] == [] and g01["interactions"]["rows"] == [] and g01["summary"]["n_interactions"] == 0, "a metabolic map has no relations")
gsyms = {}
for line in requests.get("https://rest.kegg.jp/list/hsa", timeout=120).text.splitlines():
    parts = line.split("\t")
    if len(parts) >= 4:
        gsyms[parts[0]] = {x.strip().upper() for x in parts[3].split(";")[0].split(",")}
on_map = set()
for e in gene_lines:
    for gid in e.get("name").split():
        on_map |= gsyms.get(gid, set())
exp01100 = on_map & genomics_symbols
got01100 = {x.upper() for x in hit_labels(g01, "genomics")}
check(len(exp01100) >= 20 and got01100 == exp01100,
      f"INDEPENDENT: the Genomics hits on the global map equal the genes found by re-reading KEGG's file ({len(got01100)} of {len(genomics_symbols)}; "
      f"missing {sorted(exp01100 - got01100)[:4]}, extra {sorted(got01100 - exp01100)[:4]})")
check(len(hit_labels(g01, "metabolomics", "compound")) >= 20 and g01["summary"]["per_omics"]["metabolomics"]["members_hit"] >= 20,
      f"{g01['summary']['per_omics']['metabolomics']['members_hit']} Metabolomics compounds are hit on the global map")
check(len(json.dumps(g01)) < 6_000_000, f"the payload of the biggest map stays small enough to send ({len(json.dumps(g01)) / 1e6:.1f} MB)")

print("=== hsa01100: the lines are coloured by pathway CLASS, named from KEGG's own classification ===")
br_text = requests.get("https://rest.kegg.jp/get/br:br08901", timeout=60).text
br_class, cur_a, cur_b = {}, "", ""
for ln in br_text.splitlines():            # an independent reading of the classification: A / B / C lines
    if ln.startswith("A"):
        cur_a, cur_b = ln[1:].strip(), ""
    elif ln.startswith("B"):
        cur_b = ln[1:].strip()
    elif ln.startswith("C") and cur_a == "Metabolism" and cur_b:
        br_class[ln[1:].split()[0]] = cur_b
cls_names = {c["name"] for c in g01["classes"]}
check(8 <= len(g01["classes"]) <= 14 and cls_names <= set(br_class.values()) | {"Other"},
      f"SENSITIVITY: {len(g01['classes'])} pathway classes, every one a class of KEGG's own metabolism classification ({sorted(cls_names)})")
raw_color = {e.get("id"): e.find("graphics").get("fgcolor").upper() for e in gene_lines}
check(all(n.get("cls") in cls_names for n in gnodes) and sum(c["n_lines"] for c in g01["classes"]) == len(gnodes),
      f"every one of the {len(gnodes)} line genes carries a class and is counted in exactly one")
by_color_cls = {}
for n in gnodes:
    by_color_cls.setdefault(raw_color[n["id"][1:]], set()).add(n["cls"])
check(all(len(v) == 1 for v in by_color_cls.values()), f"INDEPENDENT: one colour is one class - every one of KEGG's {len(by_color_cls)} line colours maps to a single class")
check({c["color"] for c in g01["classes"]} <= set(raw_color.values()), "each class's default colour is one KEGG itself uses on the map")
tca = [n for n in gnodes if n["label"] == "CS"]
check(tca and {n["cls"] for n in tca} == {br_class["00020"]} == {"Carbohydrate metabolism"},
      f"INDEPENDENT: citrate synthase (a TCA-cycle enzyme, hsa00020) is coloured as Carbohydrate metabolism, as br08901 says ({sorted({n['cls'] for n in tca})})")
raw_links = {e.get("id"): e.get("name") for e in entries if e.get("type") == "map" and not (e.find("graphics").get("name") or "").upper().startswith("TITLE:")}
lnodes = [n for n in g01["nodes"] if n["kind"] == "link"]
exp_cls = {nid[1:]: br_class.get(raw_links[nid[1:]][-5:]) for nid in (n["id"] for n in lnodes)}
check(all(n.get("cls") == exp_cls[n["id"][1:]] for n in lnodes) and sum(1 for n in lnodes if n.get("cls")) > 100,
      f"INDEPENDENT: all {len(lnodes)} pathway labels carry the class br08901 gives their map ({sum(1 for n in lnodes if n.get('cls'))} metabolic ones)")
check(sum(c["n_hit"] for c in g01["classes"]) == sum(1 for n in gnodes if any(v["state"] == "hit" for v in n["status"].values())) > 20,
      "the per-class hit counts add up to the genes that were hit")

print("=== KEGG signalling map hsa04210 (apoptosis): relations and groups ===")
raw04210 = requests.get("https://rest.kegg.jp/get/hsa04210/kgml", timeout=60).text
root04210 = _ET.fromstring(raw04210)
ids04210 = {e.get("id") for e in root04210.findall("entry")}
exp_rel = [r for r in root04210.findall("relation") if r.get("type") != "maplink" and r.get("entry1") in ids04210 and r.get("entry2") in ids04210]
g04 = wait_job(post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa04210"}).json()["job_id"])
check(g04["status"] == "ok" and not g04["line_graphics"] and g04["reactions"] == [] and len(g04["relations"]) == len(exp_rel) == 135,
      f"INDEPENDENT: all {len(exp_rel)} relations of KEGG's own file are in the result (a signalling map has no reactions)")
exp_pairs = sorted((r.get("entry1"), r.get("entry2")) for r in exp_rel)
check(sorted((r["a"][1:], r["b"][1:]) for r in g04["relations"]) == exp_pairs, "...between the same pairs of boxes")
sub_names = {st.get("name") for r in exp_rel for st in r.findall("subtype")}
check({"activation", "inhibition", "phosphorylation"} <= {t for r in g04["relations"] for t in r["types"]}
      and {t for r in g04["relations"] for t in r["types"]} <= (sub_names | {"pprel"}), "the subtypes are KEGG's own")
groups = [n for n in g04["nodes"] if n["kind"] == "group"]
raw_groups = {e.get("id"): [c.get("id") for c in e.findall("component")] for e in root04210.findall("entry") if e.get("type") == "group"}
check(len(groups) == len(raw_groups) == 4 and all(n["components"] == ["k" + c for c in raw_groups[n["id"][1:]]] for n in groups),
      f"INDEPENDENT: the {len(groups)} groups hold exactly the components KEGG lists")
by_id4 = {n["id"]: n for n in g04["nodes"]}
def hit4(n):
    return any(v["state"] == "hit" for v in n["status"].values())
check(all(hit4(n) == any(hit4(by_id4[c]) for c in n["components"]) for n in groups), "SENSITIVITY: a group is hit exactly when one of its components is")
exp_int = sum(1 for r in g04["relations"] if hit4(by_id4[r["a"]]) and hit4(by_id4[r["b"]]))
check(exp_int == len(g04["interactions"]["rows"]) == g04["summary"]["n_interactions"] and exp_int > 0,
      f"INDEPENDENT: the {exp_int} interaction links are exactly the relations with BOTH ends hit (recounted from the payload)")
check(all(r["From hit by"] and r["To hit by"] and r["Relation"] for r in g04["interactions"]["rows"]), "each interaction row says who hit each end and the relation")
dl_i = get("/api/multiomics/mapping/download?table=interactions")
check(dl_i.status_code == 200 and len(pd.read_csv(__import__("io").StringIO(dl_i.text))) == exp_int, "the interaction links download as CSV")
mem4 = g04["members"]["rows"]
check(all(m["Box ID"].startswith("k") for m in mem4) and not any(m["Pathway node"].startswith("undefined") for m in mem4),
      "member rows carry the box id; no row for a group")
box_counts = {}
for n in g04["nodes"]:
    if n["kind"] == "gene":
        box_counts.setdefault(n["label"], 0)
        box_counts[n["label"]] += 1
multi = [l for l, c in box_counts.items() if c > 1 and any(hit4(n) for n in g04["nodes"] if n["label"] == l)]
if multi:
    lab = multi[0]
    n_hit_members = len({m["Member"] for m in mem4 if m["Member"] == lab and m["Status"] == "hit"})
    check(n_hit_members == 1 and g04["summary"]["per_omics"]["transcriptomics"]["members_hit"] <= g04["summary"]["per_omics"]["transcriptomics"]["nodes_hit"],
          f"SENSITIVITY: {lab}, drawn in {box_counts[lab]} boxes, is ONE member hit (members hit <= boxes hit)")

print("=== Reactome overview pathways (Disease) list their sub-pathways ===")
dis = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": "R-HSA-1643685"}).json()["job_id"])
check(dis["status"] == "overview" and dis["nodes"] == [] and dis["pathway_name"] == "Disease" and len(dis["subpathways"]["rows"]) == 12,
      f"SENSITIVITY: 'Disease' is reported as an overview with its 12 sub-pathways (it used to come back as a blank chart)")
subs = dis["subpathways"]["rows"]
check(all(r["Term ID"].startswith("R-HSA-") for r in subs) and "Diseases of metabolism" in {r["Term Name"] for r in subs}
      and {r["Term Name"] for r in subs} >= {"Diseases of programmed cell death", "Diseases of Immune System"},
      "the sub-pathways are named and carry Reactome stable ids")
keyed = [(-r["Omics with a hit"], -r["Total hits"], r["Term Name"]) for r in subs]
check(keyed == sorted(keyed), "ranked by omics hitting, then by hits")
pick_rows = {r["Term ID"]: r for r in reactome_pick["table"]["rows"]}
check(all((r["Total hits"] == pick_rows[r["Term ID"]]["Total hits"]) if r["Term ID"] in pick_rows else (r["Total hits"] == 0) for r in subs),
      "INDEPENDENT: each sub-pathway's hit count equals the picker's row for the same pathway (0 where the picker lists none)")
def find_diagram(rows, depth=0, path=()):
    """Reactome nests overviews (Disease > Diseases of metabolism > ...): follow the sub-pathways down until one has a
    diagram of its own. Returns (path of names, result) for the first pathway that maps."""
    for r in rows:
        out = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": r["Term ID"], "pathway_name": r["Term Name"]}).json()["job_id"])
        if out["status"] == "ok":
            return path + (r["Term Name"],), out
        if out["status"] == "overview" and depth < 3:
            found = find_diagram(out["subpathways"]["rows"][:4], depth + 1, path + (r["Term Name"],))
            if found:
                return found
    return None


mapped_sub = find_diagram(subs[:3])
check(mapped_sub is not None and len(mapped_sub[1]["nodes"]) > 5,
      f"SENSITIVITY: following the sub-pathways down reaches a real diagram (Reactome nests its overviews): "
      f"{' > '.join(mapped_sub[0]) if mapped_sub else None} ({len(mapped_sub[1]['nodes']) if mapped_sub else 0} boxes)")
for sid in ("R-HSA-1430728", "R-HSA-168256"):
    o = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": sid}).json()["job_id"])
    check(o["status"] == "overview" and o["subpathways"]["rows"], f"{o['pathway_name']} ({sid}) is an overview too, with {len(o['subpathways']['rows'])} sub-pathways")

print("=== missing / unknown pathways are explained ===")
bad = wait_job(post("/api/multiomics/mapping/run", json={"database": "Reactome", "pathway_id": "R-HSA-999999999"}).json()["job_id"])
check(bad["status"] == "nodata" and bad["nodes"] == [] and bad["note"], f"an unknown Reactome pathway: no map, with a reason ({bad['note']})")
bad2 = wait_job(post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa99999"}).json()["job_id"])
check(bad2["status"] == "nodata" and bad2["nodes"] == [], f"an unknown KEGG pathway: no map, with a reason ({bad2['note']})")

# ---------------------------------------------------------------------------
print("=== re-uploading clears the mapping ===")
post("/api/multiomics/genomics/upload", files={"file": ("genes.csv", example("example_multiomics_genes.csv"))})
r = post("/api/multiomics/mapping/run", json={"database": "KEGG", "pathway_id": "hsa00020"})
check(r.status_code == 400, "after a re-upload the mapping must be configured again")
r = get("/api/multiomics/mapping/download?table=members")
check(r.status_code == 400, "...and there is nothing to download")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
