"""Live end-to-end test of Multi-omics Joint Pathway Enrichment through the HTTP
API (no browser - see playwright_multiomics.py for the UI).

Drives the REAL flow against a running server on port 8009 with real network
access: choose omics -> upload one file per omics -> resolve every omics with
its OWN resolver (UniProt / mygene.info / UniProt ID mapping / PubChem + the
crosswalk) -> KEGG organism -> configure -> run -> download.

It asserts on the BIOLOGY the bundled examples plant, not just on shapes: the
Genomics example is TCA-cycle + glycolysis enzymes and the Metabolomics example
is TCA-cycle + glycolysis intermediates (plus a few lipids), so the joint
result must put KEGG's TCA cycle / glycolysis at the top, with a combined
p-value STRONGER than either omics alone - the whole point of combining.

THE HEADLINE GUARANTEE, checked against the real single-omics tiles in a second
session: the per-omics p-value inside a joint result equals the p-value that
omics' OWN ORA tile reports for the same input (Metabolomics ORA's KEGG and
Genomics ORA's KEGG + Reactome). Nothing here is allowed to drift from them.

Also covered: a four-omics run with DE tables (Transcriptomics + Proteomics,
where Omiscera applies the cutoffs itself, up and down pooled), both combination
methods, the partial-coverage table, the CSV downloads, the one-organism rule
(KEGG organism reported once, "hsa"), and the input-contract errors.
"""
import io
import os
import re
import time

import pandas as pd
import requests

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
# A NEW session per run: the server keeps sessions in memory, and re-choosing omics keeps their uploads (by design),
# so a fixed id would make "before uploading" false on the second run against the same server.
_RUN = str(int(time.time() * 1000))
HEADERS = {"X-Session-Id": "multiomics-live-e2e-" + _RUN}
CMP = {"X-Session-Id": "multiomics-live-compare-" + _RUN}
PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def post(path, headers=HEADERS, **kwargs):
    return requests.post(BASE + path, headers=headers, timeout=900, **kwargs)


def get(path, headers=HEADERS, **kwargs):
    return requests.get(BASE + path, headers=headers, timeout=900, **kwargs)


def wait_job(job_id, headers=HEADERS, timeout=2400):
    start = time.time()
    while time.time() - start < timeout:
        job = get(f"/api/jobs/{job_id}", headers=headers).json()
        if job.get("status") == "done":
            return job.get("result")
        if job.get("status") == "error":
            raise RuntimeError(job.get("error") or job)
        time.sleep(2)
    raise TimeoutError(job_id)


def example(name):
    return open(os.path.join(_ROOT, "Data3", name), "rb").read()


def rows_of(payload):
    return payload["rows"]


def row(rows, term_id):
    return next((r for r in rows if r["Term ID"] == term_id), None)


def close(a, b, tol=1e-9):
    """RELATIVE tolerance. An absolute one (tol * max(1, |b|)) cannot tell 1e-40 from 1e-30 - every comparison of two
    tiny p-values would pass whatever they are. Both exactly 0.0 counts as equal."""
    if a is None or b is None:
        return False
    a, b = float(a), float(b)
    if a == b:
        return True
    return abs(a - b) <= tol * max(abs(a), abs(b))


# ---------------------------------------------------------------------------
print("=== input-contract errors are explained, not crashes ===")
r = post("/api/multiomics/select", json={"omics_types": ["genomics"]})
check(r.status_code == 400 and "two" in r.text, f"one omics is refused ({r.status_code}: {r.text[:90]})")
r = post("/api/multiomics/select", json={"omics_types": ["genomics", "metabolomics"]})
check(r.status_code == 200 and r.json()["selected"] == ["genomics", "metabolomics"], "genomics + metabolomics selected")
r = post("/api/multiomics/proteomics/upload", files={"file": ("p.csv", b"a,b\n1,2\n")})
check(r.status_code == 400 and "not one of the selected" in r.text, "uploading for an omics that was not selected is refused")
r = post("/api/multiomics/genomics/upload", files={"file": ("bad.csv", b"Compound,Value\nx,1\n")})
check(r.status_code == 400, f"a gene file with no gene column is rejected with a message ({r.text[:100]})")
r = post("/api/multiomics/metabolomics/upload", files={"file": ("bad.csv", b"Compound_ID,Value\nx,1\n")})
check(r.status_code == 400 and "Name" in r.text, "a compound file with no Name column is rejected naming the column")
r = post("/api/multiomics/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"})
check(r.status_code == 400, "resolving before uploading is refused")
r = post("/api/multiomics/joint/configure", json={})
check(r.status_code == 400, "configuring before resolving is refused")
r = get("/api/multiomics/example?omics_type=bogus")
check(r.status_code == 404, "an unknown example type is a 404")
for o in ("genomics", "transcriptomics", "proteomics", "metabolomics"):
    check(get(f"/api/multiomics/example?omics_type={o}").status_code == 200, f"example file served for {o}")

# ---------------------------------------------------------------------------
print("=== two omics: Genomics + Metabolomics ===")
r = post("/api/multiomics/genomics/upload", files={"file": ("genes.csv", example("example_multiomics_genes.csv"))})
check(r.status_code == 200 and r.json()["n_rows"] == 29 and r.json()["identifier_type"] == "gene_symbol",
      f"gene list uploaded, identifier type detected ({r.json().get('n_rows')} genes, {r.json().get('identifier_type')})")
r = post("/api/multiomics/metabolomics/upload", files={"file": ("m.csv", example("example_metabolomics_list.csv"))})
check(r.status_code == 200 and r.json()["n_rows"] >= 50 and r.json()["kind"] == "list",
      f"compound list uploaded ({r.json().get('n_rows')} compounds)")

r = post("/api/multiomics/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"})
check(r.status_code == 200, "resolve started")
res = wait_job(r.json()["job_id"])
gsum, msum = res["omics"]["genomics"], res["omics"]["metabolomics"]
check(gsum["n_usable"] >= 25, f"Genomics: most genes resolved by Genomics' own resolver ({gsum['n_usable']} of 29 usable)")
check(msum["n_usable"] >= 40, f"Metabolomics: compounds matched through the existing crosswalk ({msum['n_usable']} usable)")
check(res["taxon_id"] == "9606" and res["taxon_name"] == "Homo sapiens", "the ONE organism is recorded")

r = post("/api/multiomics/kegg/resolve-organism")
kegg = r.json()
check(r.status_code == 200 and kegg["organism_code"] == "hsa" and not kegg.get("note"),
      f"the shared resolver derives KEGG's organism for the run (got {kegg.get('organism_code')})")

r = post("/api/multiomics/joint/configure", json={"databases": ["KEGG", "Reactome"], "method": "fisher", "min_term_size": 3})
check(r.status_code == 200 and r.json()["config"]["max_term_size"] == 500, "configure accepted")
r = post("/api/multiomics/joint/run")
check(r.status_code == 200, "run started")
out = wait_job(r.json()["job_id"])
check(set(out["blocks"]) == {"KEGG", "Reactome"} and out["summary"]["kegg_organism_code"] == "hsa"
      and out["summary"]["omics"] == ["genomics", "metabolomics"],
      f"both databases ran; the ONE KEGG organism is reported ({out['summary']['kegg_organism_code']})")
check("Fisher" in out["summary"]["method_label"], "the method is named in the result")

kb, rb = out["blocks"]["KEGG"], out["blocks"]["Reactome"]
krows = rows_of(kb["table"])
check(kb["n_terms_combined"] > 5 and kb["omics_combined"] == ["genomics", "metabolomics"] and not kb["omics_excluded"],
      f"KEGG: pathways combined across both omics ({kb['n_terms_combined']})")
check(all("Per-omics p-values are combined as if independent" in (b["caveat"] or "") for b in out["blocks"].values()),
      "every block carries the independence caveat")
top = {r_["Term ID"]: r_ for r_ in krows if r_["Combined FDR"] is not None and r_["Combined FDR"] < 0.05}
for tid, name in (("hsa00020", "TCA cycle"), ("hsa00010", "glycolysis")):
    rw = top.get(tid)
    check(rw is not None, f"BIOLOGY: KEGG {name} ({tid}) is jointly significant (Combined FDR < 0.05)")
    if rw:
        best_single = min(rw["Genomics: p-value"], rw["Metabolomics: p-value"])
        check(rw["Combined p-value"] < best_single,
              f"BIOLOGY: {name}: the joint p-value ({rw['Combined p-value']:.2e}) is stronger than either omics alone "
              f"(best single {best_single:.2e})")
check(rb["n_terms_combined"] > 0 and any(re.search(r"citric|tca|glycoly|pyruvate", r_["Term Name"], re.I)
                                         and r_["Combined FDR"] < 0.05 for r_ in rows_of(rb["table"])[:25]),
      "BIOLOGY: Reactome puts a TCA / glycolysis / pyruvate pathway in its top results too")
fdrs = [r_["Combined FDR"] for r_ in krows]
check(fdrs == sorted(fdrs), "the combined table is sorted by Combined FDR")
check(len(set(r_["Term ID"] for r_ in krows)) == len(krows), "no pathway appears twice")
check(all(r_["Genomics: k"] is not None and r_["Metabolomics: k"] is not None for r_ in krows),
      "every combined row has both omics' counts (partial pathways live in their own table)")
check(any(r_["Metabolomics: k"] == 0 and r_["Genomics: k"] > 0 for r_ in krows),
      "a pathway hit by only one omics is still combined, the silent omics entering as p = 1.0 (k = 0 kept)")
ptable = rows_of(kb["partial_table"])
check(all(r_["Coverage"].startswith("partial") for r_ in ptable),
      f"the partial-coverage table is labelled ({len(ptable)} pathways tested in only one omics)")
check(all(r_["Combined p-value"] is None for r_ in ptable), "...and never carries a combined p-value")

r = get("/api/multiomics/joint/download?database=KEGG&table=combined")
df = pd.read_csv(io.StringIO(r.text))
check(r.status_code == 200 and "Combined FDR" in df.columns and len(df) == kb["n_terms_combined"],
      "the combined table downloads as CSV with every row")
check(get("/api/multiomics/joint/download?database=KEGG&table=partial").status_code == 200,
      "the partial table downloads too")
check(get("/api/multiomics/joint/download?database=GO").status_code == 400, "a database that was not run is refused")

# ---------------------------------------------------------------------------
print("=== SENSITIVITY: per-omics p-values equal that omics' OWN ORA tile (separate session) ===")
r = post("/api/metabolomics/list/upload", headers=CMP, files={"file": ("m.csv", example("example_metabolomics_list.csv"))})
r = post("/api/metabolomics/list/resolve", headers=CMP, json={"taxon_id": "9606", "taxon_name": "Homo sapiens"})
wait_job(r.json()["job_id"], CMP)
post("/api/metabolomics/ora/configure", headers=CMP, json={"databases": ["KEGG", "Reactome"], "min_term_size": 3})
mora = wait_job(post("/api/metabolomics/ora/run", headers=CMP).json()["job_id"], CMP)
m_kegg = rows_of(mora["tables"]["KEGG"]["table"])
checked = 0
for tid in [r_["Term ID"] for r_ in krows if (r_["Metabolomics: k"] or 0) > 0][:12]:
    own = row(m_kegg, tid)
    if own is not None:
        checked += 1
        ours = row(krows, tid)["Metabolomics: p-value"]
        check(close(ours, own["p-value"]), f"Metabolomics ORA tile == joint, KEGG {tid}: {own['p-value']:.4g} vs {ours:.4g}")
check(checked >= 5, f"enough Metabolomics KEGG pathways compared ({checked})")
m_rea = rows_of(mora["tables"]["Reactome"]["table"])
rchecked = 0
for r_ in rows_of(rb["table"])[:15]:
    if (r_["Metabolomics: k"] or 0) > 0:
        own = row(m_rea, r_["Term ID"])
        if own is not None:
            rchecked += 1
            check(close(r_["Metabolomics: p-value"], own["p-value"]),
                  f"Metabolomics ORA tile == joint, Reactome {r_['Term ID']}")
check(rchecked >= 3, f"enough Metabolomics Reactome pathways compared ({rchecked})")

r = post("/api/ora/genes/upload", headers=CMP, files={"file": ("g.csv", example("example_multiomics_genes.csv"))},
         data={"identifier_type": "gene_symbol"})
r = post("/api/ora/genes/resolve", headers=CMP, json={"taxon_id": "9606"})
gres = wait_job(r.json()["job_id"], CMP)
r = post("/api/ora/genes/resolve-kegg", headers=CMP, json={"kegg_organism_code": "hsa"})
kres = wait_job(r.json()["job_id"], CMP)
if kres.get("auto_final") is None:
    post("/api/ora/genes/discard-all", headers=CMP)
r = post("/api/ora/analysis/configure", headers=CMP, json={
    "databases": ["KEGG", "Reactome"],
    "per_database": {"KEGG": {"taxon_id": "hsa"}, "Reactome": {"taxon_id": "9606", "taxon_name": "Homo sapiens"}}})
check(r.status_code == 200, "Genomics ORA configured for comparison")
gora = wait_job(post("/api/ora/analysis/run", headers=CMP).json()["job_id"], CMP)
g_kegg = rows_of(gora["KEGG"]["KEGG"]["table"])
gchecked = 0
for r_ in krows:
    if (r_["Genomics: k"] or 0) > 0:
        own = row(g_kegg, r_["Term ID"])
        if own is not None:
            gchecked += 1
            check(close(r_["Genomics: p-value"], own["p-value"]),
                  f"Genomics ORA tile == joint, KEGG {r_['Term ID']}: {own['p-value']:.4g} vs {r_['Genomics: p-value']:.4g}")
        if gchecked >= 12:
            break
check(gchecked >= 5, f"enough Genomics KEGG pathways compared ({gchecked})")
g_rea = rows_of(gora["Reactome"]["Reactome"]["table"])
grchecked = 0
for r_ in rows_of(rb["table"]):
    if (r_["Genomics: k"] or 0) > 0:
        own = row(g_rea, r_["Term ID"])
        if own is not None:
            grchecked += 1
            check(close(r_["Genomics: p-value"], own["p-value"]), f"Genomics ORA tile == joint, Reactome {r_['Term ID']}")
        if grchecked >= 8:
            break
check(grchecked >= 3, f"enough Genomics Reactome pathways compared ({grchecked})")

# ---------------------------------------------------------------------------
print("=== Tippett (minimum-p) is a real second method ===")
r = post("/api/multiomics/joint/configure", json={"databases": ["KEGG"], "method": "tippett", "min_term_size": 3})
tout = wait_job(post("/api/multiomics/joint/run").json()["job_id"])
trows = rows_of(tout["blocks"]["KEGG"]["table"])
t_tca = row(trows, "hsa00020")
check("Tippett" in trows[0]["Method"] and t_tca is not None,
      "Tippett run: the method is recorded on every row and TCA is still present")
if t_tca:
    # -expm1(k * log1p(-p)) == 1 - (1 - p)^k without the float cancellation (for p ~ 1e-40 the naive form is 0.0)
    import math
    expect = -math.expm1(2 * math.log1p(-min(t_tca["Genomics: p-value"], t_tca["Metabolomics: p-value"])))
    check(close(t_tca["Combined p-value"], expect, tol=1e-6), f"Tippett's combined p = 1-(1-min p)^2 ({t_tca['Combined p-value']:.3g})")
check(post("/api/multiomics/joint/configure", json={"method": "stouffer"}).status_code == 400,
      "Stouffer is refused by the API (it is invalid for one-sided enrichment p-values)")

# ---------------------------------------------------------------------------
print("=== four omics: DE tables for Transcriptomics + Proteomics (cutoffs applied by Omiscera) ===")
r = post("/api/multiomics/select", json={"omics_types": ["genomics", "transcriptomics", "proteomics", "metabolomics"]})
check(r.status_code == 200, "all four omics selected (Genomics + Metabolomics uploads are kept)")
r = post("/api/multiomics/transcriptomics/upload", files={"file": ("t.csv", example("example_transcriptomics_de_full.csv"))})
check(r.status_code == 200 and r.json()["kind"] == "de" and r.json()["column_mapping"].get("pvalue"),
      f"a DE table is accepted; its columns were detected ({r.json().get('column_mapping')})")
r = post("/api/multiomics/transcriptomics/upload", files={"file": ("t.csv", example("example_multiomics_genes.csv"))})
check(r.status_code == 400, "a plain gene list is REJECTED for Transcriptomics (a DE table with log2FC + p-value is required)")
post("/api/multiomics/transcriptomics/upload", files={"file": ("t.csv", example("example_transcriptomics_de_full.csv"))})
r = post("/api/multiomics/proteomics/upload", files={"file": ("p.csv", example("example_proteomics_de_full.csv"))})
check(r.status_code == 200 and r.json()["kind"] == "de", "the protein-level DE table is accepted")
r = post("/api/multiomics/resolve", json={"taxon_id": "9606", "taxon_name": "Homo sapiens"})
res4 = wait_job(r.json()["job_id"])
check(set(res4["omics"]) == {"genomics", "transcriptomics", "proteomics", "metabolomics"}
      and all(v["n_usable"] > 0 for v in res4["omics"].values()),
      "all four omics resolved with their own resolvers: " + ", ".join(f"{k}={v['n_usable']}" for k, v in res4["omics"].items()))
r = post("/api/multiomics/joint/configure", json={
    "databases": ["KEGG", "Reactome"], "method": "fisher", "min_term_size": 5,
    "de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0}, "proteomics": {"padj": 0.05, "lfc": 1.0}}})
check(r.status_code == 200 and set(r.json()["config"]["de_cutoffs"]) == {"transcriptomics", "proteomics"}, "four-omics config accepted")
out4 = wait_job(post("/api/multiomics/joint/run").json()["job_id"])
osum = out4["summary"]["omics_summary"]
check(osum["transcriptomics"]["n_in_list"] > 0 and osum["transcriptomics"]["n_tested"] > osum["transcriptomics"]["n_in_list"]
      and osum["proteomics"]["n_in_list"] > 0,
      f"significant sets come from the DE tables (tx {osum['transcriptomics']['n_in_list']} of {osum['transcriptomics']['n_tested']}, "
      f"prot {osum['proteomics']['n_in_list']} of {osum['proteomics']['n_tested']}), the reference is what was TESTED")
b4 = out4["blocks"]["KEGG"]
check(b4["n_terms_combined"] > 0 and set(b4["omics_combined"]) == {"genomics", "transcriptomics", "proteomics", "metabolomics"},
      f"KEGG combines all four omics ({b4['n_terms_combined']} pathways)")
r4 = rows_of(b4["table"])
check(set(x["N omics combined"] for x in r4) == {4}, "every combined row used four p-values (2k = 8 degrees of freedom)")
check(all(f"{o}: p-value" in b4["columns"] for o in ("Genomics", "Transcriptomics", "Proteomics", "Metabolomics")),
      "the table carries every omics' own columns beside the combined ones")
check(out4["summary"]["kegg_organism_code"] == "hsa", "still ONE KEGG organism for four omics")

# ---------------------------------------------------------------------------
print("=== Brown's method (the default): correlated omics are not double-counted ===")
import math as _math  # noqa: E402
from scipy.stats import chi2 as _chi2  # noqa: E402

r = post("/api/multiomics/joint/configure", json={
    "databases": ["KEGG", "Reactome"], "min_term_size": 5,
    "de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0}, "proteomics": {"padj": 0.05, "lfc": 1.0}}})
check(r.status_code == 200 and r.json()["config"]["method"] == "brown", "SENSITIVITY: with no method given the API picks Brown's")
brown = wait_job(post("/api/multiomics/joint/run").json()["job_id"])
post("/api/multiomics/joint/configure", json={
    "databases": ["KEGG", "Reactome"], "method": "fisher", "min_term_size": 5,
    "de_cutoffs": {"transcriptomics": {"padj": 0.05, "lfc": 1.0}, "proteomics": {"padj": 0.05, "lfc": 1.0}}})
fisher4 = wait_job(post("/api/multiomics/joint/run").json()["job_id"])
OM = ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]
check("Brown" in brown["summary"]["method_label"], f"the result names Brown's method ({brown['summary']['method_label']})")

for db in ("KEGG", "Reactome"):
    bb, ff = brown["blocks"][db], fisher4["blocks"][db]
    corr = bb["correlation"]
    pair = {(p_["a"], p_["b"]): p_ for p_ in corr["pairs"]}
    check(len(corr["pairs"]) == 6 and corr["applied"] and not corr["fallback"],
          f"[{db}] four omics -> six estimated pairwise correlations, applied, estimated over {corr['n_pathways']} pathways")
    gm = pair[("genomics", "metabolomics")]
    check(gm["rho"] > 0.1,
          f"[{db}] BIOLOGY: Genomics and Metabolomics, which carry the SAME planted TCA/glycolysis biology, correlate across "
          f"pathways (rho = {gm['rho']:.2f}) - the dependence plain Fisher ignores. (The example Transcriptomics and "
          f"Proteomics tables are independent synthetic files: their rho = "
          f"{pair[('transcriptomics', 'proteomics')]['rho']:.2f} is reported, not asserted.)")
    check(all(0.0 <= p_["rho"] <= 1.0 for p_ in corr["pairs"]), f"[{db}] every correlation used lies in [0, 1] (rhos: "
          + ", ".join(f"{p_['a'][:3]}-{p_['b'][:3]} {p_['rho']:.2f}" for p_ in corr["pairs"]) + ")")
    check(1.0 <= corr["effective_omics"] < 4.0, f"[{db}] the four omics count as {corr['effective_omics']:.2f} independent ones")
    check(ff["correlation"]["applied"] is False and abs(ff["correlation"]["pairs"][0]["rho"] - corr["pairs"][0]["rho"]) < 1e-12,
          f"[{db}] a Fisher run reports the SAME correlation but does not apply it")
    k = len(bb["omics_combined"])
    rho_sum = sum(p_["rho"] for p_ in corr["pairs"])
    var = 4 * k + 8 * rho_sum
    c_, f_ = var / (4 * k), 8 * k * k / var
    check(close(corr["c"], c_) and close(corr["f"], f_), f"[{db}] reported c = {corr['c']:.3f}, f = {corr['f']:.3f} follow from the correlations")
    bad_p = bad_dir = 0
    frows = {r_["Term ID"]: r_ for r_ in ff["table"]["rows"]}
    for r_ in bb["table"]["rows"]:
        ps = [min(max(r_[f"{o}: p-value"], 1e-300), 1.0) for o in OM]
        X = -2 * sum(_math.log(p_) for p_ in ps)
        if not close(r_["Combined p-value"], float(_chi2.sf(X / c_, f_)), 1e-6):
            bad_p += 1
        fr = frows.get(r_["Term ID"])
        if fr and fr["Combined p-value"] <= 0.28 and r_["Combined p-value"] < fr["Combined p-value"] * (1 - 1e-9):
            bad_dir += 1
    check(bad_p == 0, f"[{db}] SENSITIVITY: all {len(bb['table']['rows'])} combined p-values equal an independent chi-squared "
          f"recomputation from the four per-omics p-values and the reported c, f ({bad_p} mismatches)")
    check(bad_dir == 0, f"[{db}] SENSITIVITY: wherever Fisher's p <= 0.28, Brown's p is never smaller than Fisher's ({bad_dir} violations)")
    check(set(frows) == set(r_["Term ID"] for r_ in bb["table"]["rows"]), f"[{db}] Brown and Fisher combine the same pathways")

# the case that motivated Brown: expression layers dominate the Fisher ranking
rb_top = fisher4["blocks"]["Reactome"]["table"]["rows"][:5]
brown_rows = {r_["Term ID"]: r_ for r_ in brown["blocks"]["Reactome"]["table"]["rows"]}
weaker = [r_["Combined p-value"] < brown_rows[r_["Term ID"]]["Combined p-value"] for r_ in rb_top]
check(all(weaker),
      "BIOLOGY: Fisher's five strongest Reactome hits (RNA + protein agreeing, e.g. "
      + rb_top[0]["Term Name"][:40] + ") all become LESS extreme under Brown, by "
      + ", ".join(f"{brown_rows[r_['Term ID']]['Combined p-value'] / max(r_['Combined p-value'], 1e-300):.0e}x" for r_ in rb_top[:3]))
kb = {r_["Term ID"]: i for i, r_ in enumerate(brown["blocks"]["KEGG"]["table"]["rows"])}
check(kb.get("hsa00020", 99) < 5 and brown["blocks"]["KEGG"]["table"]["rows"][kb["hsa00020"]]["Combined FDR"] < 0.05,
      f"BIOLOGY: the TCA cycle (planted in the genes AND the compounds) is still a top-5, significant KEGG pathway under Brown "
      f"(rank {kb.get('hsa00020', 99) + 1})")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
