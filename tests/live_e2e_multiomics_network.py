"""Live end-to-end test of Multi-omics Network-based Integration through the HTTP API
(no browser - see playwright_multiomics_network.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network access: four omics uploaded and resolved with
their own resolvers -> configure -> run on the real human STRING network (500 permutations) and Reactome's real files.

Every claim is checked against an INDEPENDENT recomputation from the raw files, written here without Omiscera code:
  - the SEEDS of each omics: gene symbols / UniProt accessions looked up straight in STRING's own info / alias files, the
    DEGs picked from the uploaded CSVs at the cutoffs (proteomics' padj by Benjamini-Hochberg written out here);
  - the METABOLITE -> ENZYME join: Reactome's two bulk files are streamed again here and joined for every compound the app
    used, with the same two limits, to the same reaction and protein counts;
  - the COMBINED WEIGHTS: equal weight per omics, recomputed from the seed sets;
  - the RANDOM WALK: the STRING graph is rebuilt from the raw links file (score >= 700, column-normalised) and the walk is
    iterated to convergence; every protein's propagated score must equal the app's;
  - the pathway network's edges whose members the payload lists in full.
The bundled examples are synthetic (they plant TCA-cycle genes and metabolites); the test asserts what they plant.
"""
import gzip
import io
import os
import time

import numpy as np
import pandas as pd
import requests
import scipy.sparse as sp

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
HEADERS = {"X-Session-Id": "multiomics-net-live-" + str(int(time.time() * 1000))}
OMICS = ["genomics", "transcriptomics", "proteomics", "metabolomics"]
LABEL = {"genomics": "Genomics", "transcriptomics": "Transcriptomics", "proteomics": "Proteomics", "metabolomics": "Metabolomics"}
STRING_DIR = os.path.join(_ROOT, "cache", "string")
PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def post(path, **kw):
    return requests.post(BASE + path, headers=HEADERS, timeout=900, **kw)


def get(path, **kw):
    return requests.get(BASE + path, headers=HEADERS, timeout=900, **kw)


def wait_job(job_id, timeout=3000):
    start, last = time.time(), ""
    while time.time() - start < timeout:
        job = get(f"/api/jobs/{job_id}").json()
        if job.get("message") != last:
            last = job.get("message")
            print(f"   [{time.time() - start:5.0f}s] {job.get('progress', 0):.2f} {last}", flush=True)
        if job.get("status") == "done":
            return job.get("result")
        if job.get("status") == "error":
            raise RuntimeError(job.get("error") or job)
        time.sleep(3)
    raise TimeoutError(job_id)


def example(name):
    return os.path.join(_ROOT, "Data3", name)


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    order = np.argsort(pv)
    m = len(pv)
    adj = np.minimum.accumulate((pv[order] * m / (np.arange(m) + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


print("=== guard rails before anything is set up ===")
check(post("/api/multiomics/network/configure", json={}).status_code == 400, "configure before choosing / resolving omics is refused")
check(post("/api/multiomics/network/run", json={}).status_code == 400, "run before configuring is refused")

print("=== four omics, real resolvers ===")
check(post("/api/multiomics/select", json={"omics_types": OMICS}).status_code == 200, "four omics selected")
files = {"genomics": "example_multiomics_genes.csv", "transcriptomics": "example_transcriptomics_de_full.csv",
         "proteomics": "example_proteomics_de_full.csv", "metabolomics": "example_metabolomics_list.csv"}
for o, f in files.items():
    check(post(f"/api/multiomics/{o}/upload", files={"file": (f, open(example(f), "rb").read())}).status_code == 200, f"{o} uploaded")
res = wait_job(post("/api/multiomics/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"}).json()["job_id"])
check(all(res["omics"][o]["n_usable"] > 0 for o in OMICS), "every omics resolved")

print("=== configuration is validated ===")
for bad in ({"confidence": 123}, {"n_permutations": 77}, {"module_database": "MSigDB"}, {"fdr_alpha": 5}, {"max_reactions": 1},
            {"max_reaction_proteins": 1}, {"pathway_min_jaccard": 0}):
    check(post("/api/multiomics/network/configure", json=bad).status_code == 400, f"refused: {bad}")
r = post("/api/multiomics/network/configure", json={"n_permutations": 500, "pathway_top_n": 100})
CONFIG = r.json()["config"]
check(r.status_code == 200 and r.json()["config"]["n_permutations"] == 500 and r.json()["config"]["max_reactions"] == 60
      and r.json()["config"]["max_reaction_proteins"] == 30, "a valid configuration is stored with the documented defaults")

print("=== run (real STRING network, real Reactome files) ===")
t0 = time.time()
out = wait_job(post("/api/multiomics/network/run", json={}).json()["job_id"])
print(f"   run took {time.time() - t0:.0f}s")
prop = out["propagation"]
check(out["omics_used"] == OMICS, f"all four omics took part ({out['omics_used']})")
check(prop["n_proteins"] > 15000, f"the human STRING network at confidence 700 ({prop['n_proteins']} proteins)")

print("=== independent rebuild of STRING from the raw files ===")
links = pd.read_csv(os.path.join(STRING_DIR, "9606.protein.links.v12.0.txt.gz"), sep=" ", dtype={"protein1": str, "protein2": str})
links = links[links["combined_score"] >= 700]
info = pd.read_csv(os.path.join(STRING_DIR, "9606.protein.info.v12.0.txt.gz"), sep="\t", dtype=str)
info.columns = [c.lstrip("#").strip() for c in info.columns]
alias = pd.read_csv(os.path.join(STRING_DIR, "9606.protein.aliases.v12.0.txt.gz"), sep="\t", dtype=str)
alias.columns = [c.lstrip("#").strip() for c in alias.columns]
prot_ids = pd.unique(pd.concat([links["protein1"], links["protein2"]], ignore_index=True))
pidx = {p: i for i, p in enumerate(prot_ids)}
n = len(prot_ids)
W = sp.csr_matrix((links["combined_score"].to_numpy(float) / 1000.0, (links["protein1"].map(pidx).to_numpy(), links["protein2"].map(pidx).to_numpy())), shape=(n, n))
W = W.maximum(W.T)
colsum = np.asarray(W.sum(axis=0)).ravel()
Wn = sp.csr_matrix(W @ sp.diags(1.0 / np.where(colsum == 0, 1.0, colsum)))
name_of = dict(zip(info["string_protein_id"], info["preferred_name"]))
by_name = {}
for pid, nm in name_of.items():
    if pid in pidx:
        by_name.setdefault(str(nm).upper(), []).append(pid)
by_acc = {}
for a, pid in zip(alias.loc[alias["source"] == "UniProt_AC", "alias"], alias.loc[alias["source"] == "UniProt_AC", "string_protein_id"]):
    if pid in pidx:
        by_acc.setdefault(a.upper(), set()).add(pid)
check(n == prop["n_proteins"], f"the independent graph has the same {n} proteins")

hits_csv = pd.read_csv(io.StringIO(get("/api/multiomics/network/download", params={"which": "hits"}).content.decode("utf-8-sig")))
check(len(hits_csv) == n, f"the hit-table download carries every protein ({len(hits_csv)})")
seed_rows = hits_csv[hits_csv["Seed"] == True]  # noqa: E712
app_sets = {o: set(seed_rows.loc[seed_rows["Seed omics"].str.contains(LABEL[o]), "STRING ID"]) for o in OMICS}


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else 1.0


print("=== seeds, per omics, from the raw uploads ===")
g_syms = [s.strip().upper() for s in pd.read_csv(example(files["genomics"]))["Gene Symbol"].astype(str)]
exp_g = {pid for s in g_syms for pid in by_name.get(s, [])[:1]}
check(jaccard(app_sets["genomics"], exp_g) >= 0.9, f"Genomics seeds agree with STRING's own names (Jaccard {jaccard(app_sets['genomics'], exp_g):.3f}; {len(app_sets['genomics'])} vs {len(exp_g)})")
tde = pd.read_csv(example(files["transcriptomics"]))
tsel = tde[(tde["padj"] <= 0.05) & (tde["log2FoldChange"].abs() >= 1.0)]
exp_t = {pid for s in tsel["Gene Symbol"].astype(str).str.upper() for pid in by_name.get(s, [])[:1]}
check(jaccard(app_sets["transcriptomics"], exp_t) >= 0.9, f"Transcriptomics seeds agree with the DEGs of the uploaded table (Jaccard {jaccard(app_sets['transcriptomics'], exp_t):.3f}; {len(app_sets['transcriptomics'])} vs {len(exp_t)})")
pde = pd.read_csv(example(files["proteomics"]))
pde["padj"] = bh(pde["pvalue"].to_numpy())
psel = pde[(pde["padj"] <= 0.05) & (pde["log2FC"].abs() >= 1.0)]
exp_p = {sorted(by_acc.get(a.split(";")[0].upper(), set()))[0] for a in psel["Protein IDs"].astype(str) if by_acc.get(a.split(";")[0].upper())}
check(jaccard(app_sets["proteomics"], exp_p) >= 0.9, f"Proteomics seeds agree with the DEPs of the uploaded table, padj by Benjamini-Hochberg (Jaccard {jaccard(app_sets['proteomics'], exp_p):.3f}; {len(app_sets['proteomics'])} vs {len(exp_p)})")

print("=== metabolites -> enzymes, re-joined from Reactome's raw files ===")
met = out["metabolites"]["rows"]


def stream(url):
    # Reactome's server answers 5xx now and then and sometimes closes a large download half-way (a 522 and a 'Response ended
    # prematurely' each ended a whole run): the WHOLE read is repeated, a few times, then it fails loudly
    for attempt in range(5):
        try:
            resp = requests.get(url, stream=True, timeout=300)
            resp.raise_for_status()
            rows = []
            for raw in resp.iter_lines():
                if raw:
                    parts = raw.decode("utf-8", "replace").split("\t")
                    if len(parts) >= 6 and parts[5].strip() == "Homo sapiens":
                        rows.append((parts[0].strip(), parts[1].strip()))
            return rows
        except requests.exceptions.RequestException as exc:
            print(f"   (Reactome download failed: {exc}; retry {attempt + 1}/5)")
            time.sleep(20)
    raise RuntimeError("Reactome's file could not be downloaded after 5 attempts")


chebi_rx, rx_prot = {}, {}
for c, rx in stream("https://reactome.org/download/current/ChEBI2ReactomeReactions.txt"):
    chebi_rx.setdefault(c.replace("CHEBI:", ""), set()).add(rx)
for acc, rx in stream("https://reactome.org/download/current/UniProt2ReactomeReactions.txt"):
    rx_prot.setdefault(rx, set()).add(acc.upper())
bad, exp_m, n_checked = [], set(), 0
for row in met:
    if not row["Matched ChEBI"]:
        continue
    ids = [x.strip() for x in row["Matched ChEBI"].split(",")]
    rxs = set().union(*[chebi_rx.get(i, set()) for i in ids])
    fam = {r for r in rxs if len(rx_prot.get(r, ())) > 30}
    accs = set().union(*[rx_prot.get(r, set()) for r in rxs if r not in fam])
    n_checked += 1
    if row["Reactions"] != len(rxs) or row["Family reactions skipped"] != len(fam) or row["Proteins"] != len(accs):
        bad.append((row["Compound"], row["Reactions"], len(rxs), row["Family reactions skipped"], len(fam), row["Proteins"], len(accs)))
    if row["Status"] == "used":
        for a in accs:
            for pid in sorted(by_acc.get(a, set()))[:1]:
                exp_m.add(pid)
check(n_checked >= 30 and not bad, f"every matched compound's reaction count, skipped family reactions and protein count equal the independent join ({n_checked} compounds; {bad[:2]})")
check(all((r["Status"] == "used") == (0 < r["Reactions"] <= 60 and r["Proteins"] > 0) for r in met), "a compound is used exactly when it passes both limits and has proteins")
hubs = [r for r in met if r["Status"].startswith("left out: in")]
check(len(hubs) >= 3 and {"2-Oxoglutarate", "Succinate"} <= {r["Compound"] for r in hubs}, f"BIOLOGY: busy intermediates are left out as hubs ({[r['Compound'] for r in hubs]})")
check(jaccard(app_sets["metabolomics"], exp_m) >= 0.9, f"Metabolomics seeds equal the independent join's proteins (Jaccard {jaccard(app_sets['metabolomics'], exp_m):.3f}; {len(app_sets['metabolomics'])} vs {len(exp_m)})")
check(max(r["Proteins"] for r in met if r["Status"] == "used") <= 300 and sum(r["Proteins"] for r in met if r["Status"] == "used") < 2500,
      "no single metabolite floods the seed list (the family-reaction guard works on real data)")

print("=== combined weights: equal weight per omics ===")
sets = {o: app_sets[o] for o in OMICS}
union = set().union(*sets.values())
exp_w = {p: len(union) / 4.0 * sum(1.0 / len(sets[o]) for o in OMICS if p in sets[o]) for p in union}
got_w = dict(zip(seed_rows["STRING ID"], seed_rows["Seed score"]))
check(set(exp_w) == set(got_w) and max(abs(exp_w[p] - got_w[p]) for p in exp_w) < 1e-9,
      f"every seed's weight equals n x (1/K) x sum of 1/|S_omics| (largest difference {max(abs(exp_w[p] - got_w[p]) for p in exp_w):.1e}; {len(union)} seeds)")
only = {o: [p for p in sets[o] if all(p not in sets[q] for q in OMICS if q != o)] for o in OMICS}
ok_ratio = True
for o1 in OMICS:
    for o2 in OMICS:
        if o1 < o2 and only[o1] and only[o2]:
            ratio = got_w[only[o1][0]] / got_w[only[o2][0]]
            ok_ratio &= abs(ratio - len(sets[o2]) / len(sets[o1])) < 1e-9
            ok_ratio &= len({round(got_w[p], 12) for p in only[o1]}) == 1
check(ok_ratio and all(only[o] for o in OMICS), f"SENSITIVITY: a protein seeded by ONE omics weighs exactly 1/(that omics' seed count) as much as another omics' - the 29 genes outweigh the 552 transcripts {len(sets['transcriptomics']) / len(sets['genomics']):.1f}-fold each")
check(abs(sum(got_w.values()) - len(union)) < 1e-6 and out["seeds"]["n_seeds"] == len(union), "the weights sum to the number of distinct seeds")

print("=== the random walk, iterated independently on the rebuilt graph ===")
ALPHA = 0.5          # the engine's documented restart probability (Genomics_network.DEFAULT_RESTART_PROB)
check(out["parameters"]["restart_prob"] == ALPHA and out["parameters"]["confidence"] == 700, "the run reports restart probability 0.5 and STRING confidence 700")
s = np.zeros(n)
for p, w in got_w.items():
    s[pidx[p]] = w
s /= s.sum()
x = s.copy()
for _ in range(2000):
    x2 = (1 - ALPHA) * (Wn @ x) + ALPHA * s
    done = np.abs(x2 - x).sum() < 1e-13
    x = x2
    if done:
        break
app_score = dict(zip(hits_csv["STRING ID"], hits_csv["Propagated score"]))
ind = np.array([x[pidx[p]] for p in hits_csv["STRING ID"]])
app = np.array([app_score[p] for p in hits_csv["STRING ID"]])
rel = np.abs(ind - app) / np.maximum(np.abs(app), 1e-300)
check(np.abs(ind - app).max() < 1e-8 and np.median(rel) < 1e-6, f"every protein's propagated score equals the independent random walk (largest absolute difference {np.abs(ind - app).max():.1e}, median relative {np.median(rel):.1e})")
top = hits_csv.sort_values("Propagated score", ascending=False).head(10)
check(set(top["Gene Symbol"]) & {"CS", "ACO2", "IDH3A", "MDH2", "SDHB", "FH", "PDHA1", "DLD", "DLST", "SUCLG1"}, f"BIOLOGY: planted TCA-cycle proteins are among the ten highest propagated scores ({list(top['Gene Symbol'])})")

print("=== the significance machinery ===")
check(hits_csv["FDR"].between(0, 1).all() and (hits_csv["Z-score"].notna()).all(), "every protein has a Z-score and an FDR in [0, 1]")
sig = hits_csv[hits_csv["FDR"] < 0.05]
check(prop["n_hits"] == len(sig) and prop["n_hits_seed"] == int(sig["Seed"].sum()), f"the reported hit counts equal the table's ({prop['n_hits']} significant, {prop['n_hits_seed']} of them seeds)")
check(prop["n_recruited"] <= prop["n_hits"] - prop["n_hits_seed"] and any("circular" in c.lower() or "expected rather than a discovery" in c for c in out["caveats"]),
      "the recruited count never exceeds the non-seed hits, and the seed-circularity caveat is shown")
check(any("Equal weight per omics" in c for c in out["caveats"]) and any("association" in c for c in out["caveats"]), "the equal-weight rule and the metabolite caveat are stated")
mods = out["modules"]
if mods and mods["table"]["rows"]:
    mt = pd.DataFrame(mods["table"]["rows"])
    cols = [f"Seeds from {LABEL[o]}" for o in OMICS]
    check(set(cols) <= set(mt.columns) and (mt["Seeds"] <= mt[cols].sum(axis=1)).all(), "per module, the seeds never exceed the sum of the per-omics seed counts")
    members = {m: set(v) for m, v in mods["membership"].items()}
    ok = True
    for _, r in mt.iterrows():
        mem = members[r["Module ID"]]
        for o in OMICS:
            ok &= int(r[f"Seeds from {LABEL[o]}"]) == len(mem & sets[o])
        ok &= int(r["Seeds"]) == len(mem & union)
    check(ok, f"each module's per-omics seed counts equal the intersection of its members with the independent seed sets ({len(mt)} module(s))")
else:
    print("   (no module at this size - nothing to check)")

print("=== the Reactome pathway network ===")
pw = out["pathway_network"]
check(pw["nodes"] and pw["omics"], f"the pathway network has {len(pw['nodes'])} pathways joined by {len(pw['edges'])} links in {len(pw['groups'])} group(s)")
ids = {nd["id"] for nd in pw["nodes"]}
rank = next((i + 1 for i, x in enumerate(pw["nodes"]) if x["id"] == "R-HSA-71403"), None)
print(f"   Reactome's citric acid cycle (R-HSA-71403): rank {rank} of {len(pw['nodes'])} drawn pathways" if rank else "   Reactome's citric acid cycle is not among the drawn pathways")
shared = [x for x in pw["nodes"] if x["hits"].get("genomics", 0) > 0 and x["hits"].get("metabolomics", 0) > 0]
check(len(shared) >= 3, f"BIOLOGY: pathways hit by BOTH the Genomics list and the Metabolomics list (the planted shared signal) are in the network ({len(shared)}: {[x['name'][:28] for x in shared[:3]]})")
link_thr = CONFIG["pathway_min_jaccard"]
merge_thr = CONFIG["pathway_merge_jaccard"]
check(link_thr == 0.4 and merge_thr == 0.7, "the documented defaults: link at Jaccard 0.4, merge near-duplicates at 0.7")
check(all(e["a"] in ids and e["b"] in ids and link_thr - 1e-12 <= e["jaccard"] < merge_thr and e["shared"] >= 1 for e in pw["edges"]),
      "every link joins two drawn pathways and lies between the link threshold and the merge threshold (nothing more similar than that is drawn twice)")
check(pw["n_merged"] > 0 and all(isinstance(x["merged"], list) for x in pw["nodes"]), f"near-identical pathways were merged ({pw['n_merged']}) and each drawn pathway lists what it stands for")
nd = {x["id"]: x for x in pw["nodes"]}
complete = [e for e in pw["edges"] if all(len(nd[t]["members"].get(o, [])) == nd[t]["hits"].get(o, 0) for t in (e["a"], e["b"]) for o in pw["omics"])]
okj = all(abs(e["jaccard"] - (lambda a, b: len(a & b) / len(a | b))(
    {("M:" + m if o == "metabolomics" else m) for o in pw["omics"] for m in nd[e["a"]]["members"][o]},
    {("M:" + m if o == "metabolomics" else m) for o in pw["omics"] for m in nd[e["b"]]["members"][o]})) < 1e-9 for e in complete)
check(complete and okj, f"for the {len(complete)} links whose members the payload lists in full, the Jaccard is recomputed exactly")
pos = np.array([[x["x"], x["y"]] for x in pw["nodes"]])
bw, bh = pw["canvas"]["box_width"], pw["canvas"]["box_height"]
clear = all(abs(pos[i, 0] - pos[j, 0]) >= bw - 1e-6 or abs(pos[i, 1] - pos[j, 1]) >= bh - 1e-6 for i in range(len(pos)) for j in range(i + 1, len(pos)))
check(clear and (pos[:, 0] > 0).all() and (pos[:, 0] < pw["canvas"]["width"]).all(), f"no two pathway label boxes overlap on the {pw['canvas']['width']}x{pw['canvas']['height']} canvas")
check(pw["nodes"][0]["combined_p"] <= pw["nodes"][-1]["combined_p"], "pathways are the ones with the best combined p-values, best first")

print("=== downloads equal the payload ===")
for which, rows in [("omics", len(out["omics"])), ("pathway_nodes", len(pw["nodes"])), ("pathway_edges", len(pw["edges"]))]:
    d = get("/api/multiomics/network/download", params={"which": which})
    df = pd.read_csv(io.StringIO(d.content.decode("utf-8-sig"))) if d.content.strip() else pd.DataFrame()
    check(d.status_code == 200 and len(df) == rows, f"download '{which}' has {len(df)} rows = the table shown")
d = get("/api/multiomics/network/download", params={"which": "metabolites"})
check(len(pd.read_csv(io.StringIO(d.content.decode("utf-8-sig")))) == len(met), "download 'metabolites' has every compound")
if mods and mods["table"]["rows"]:
    check(len(pd.read_csv(io.StringIO(get("/api/multiomics/network/download", params={"which": "modules"}).content.decode("utf-8-sig")))) == len(mods["table"]["rows"]), "download 'modules' has every module")
check(get("/api/multiomics/network/download", params={"which": "nope"}).status_code == 400, "an unknown table is refused")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
