# -*- coding: utf-8 -*-
"""Unit tests for Multi-omics Joint Pathway Enrichment (Multiomics_module),
offline - fixtures only, no network.

THE HEADLINE GUARANTEE under test: the per-omics p-value inside a joint result
is the SAME number that omics' own ORA engine produces for the same input
(checked against the exact function each single-omics module uses internally:
Genomics_module._run_fisher_for_namespace, Transcriptomics_module.
_run_ora_for_sets, Metabolomics_module._drop_zero_overlap +
_run_fisher_for_namespace - the same pair Proteomics uses).

THE DESIGN RULES (each is a decision with a failure mode, mutation-tested):
  1. k == 0 pathways keep p = 1.0 for the combination (single-omics tiles drop
     them from their own FDR family; the combination must not).
  2. A pathway enters the combined family only if some omics has k > 0.
  3. Only pathways present in EVERY testable omics are combined; the rest go
     to a separate partial-coverage table, never padded with p = 1.0.
  4. An omics with no annotated hit is left out; < 2 testable omics -> empty
     result with a note.
  5. ONE Benjamini-Hochberg pass, over the combined p-values only.
  6. Stouffer is not offered (its z = -inf at p = 1 wipes out a real hit).
"""
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm  # noqa: E402
import Metabolomics_module as mm  # noqa: E402
import Multiomics_module as mo  # noqa: E402
import Proteomics_module as pm  # noqa: E402
import Transcriptomics_module as tm  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


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
print("=== combine_omics_pvalues: the combination math ===")
ps = [0.01, 0.04]
x2 = -2 * (math.log(0.01) + math.log(0.04))
f = mo.combine_omics_pvalues(ps, "fisher")
check(close(f["combined_p"], chi2.sf(x2, 4)),
      f"Fisher equals a hand-computed -2*sum(ln p) against chi-squared with 2k=4 df (got {f['combined_p']:.6g})")
check(close(f["statistic"], x2), "Fisher's statistic is -2*sum(ln p)")
check(close(mo.combine_omics_pvalues(ps, "tippett")["combined_p"], 1 - (1 - 0.01) ** 2),
      "Tippett equals 1-(1-min p)^k")
check(close(mo.combine_omics_pvalues([1e-10, 1.0])["combined_p"],
            mo.combine_omics_pvalues([1e-10])["combined_p"] * 1.0, tol=1.0)
      and mo.combine_omics_pvalues([1e-10, 1.0])["combined_p"] < 1e-8,
      "SENSITIVITY: a p = 1.0 omics does not wipe out a real hit under Fisher "
      f"(got {mo.combine_omics_pvalues([1e-10, 1.0])['combined_p']:.3g})")
check(mo.combine_omics_pvalues([1e-10, 1.0], "tippett")["combined_p"] < 1e-9,
      "...nor under Tippett")
from scipy.stats import combine_pvalues  # noqa: E402
stouffer_p = combine_pvalues([1e-10, 1.0 - 1e-16], method="stouffer")[1]
check(stouffer_p > 0.5,
      f"WHY Stouffer is not offered: the same two p-values give p = {stouffer_p:.2f} under Stouffer")
try:
    mo.combine_omics_pvalues([0.01, 0.02], "stouffer")
    check(False, "Stouffer is rejected")
except ValueError:
    check(True, "Stouffer is rejected (ValueError) - only Fisher and Tippett are offered")
check(math.isfinite(mo.combine_omics_pvalues([0.0, 0.5])["statistic"]),
      "a p-value that underflowed to 0.0 is floored so the statistic stays finite")
check(mo.combine_omics_pvalues([0.3])["combined_p"] > 0 and mo.combine_omics_pvalues([0.3])["n_combined"] == 1,
      "a single p-value combines to itself (n_combined = 1)")

# ---------------------------------------------------------------------------
print("=== _per_pathway_pvalues: reuses Genomics' engine, but keeps k == 0 ===")
TERMS = {
    "T1": {"name": "alpha pathway", "genes": {"a", "b", "c", "d"}},
    "T2": {"name": "beta pathway", "genes": {"e", "f", "g", "h"}},
    "T3": {"name": "gamma pathway", "genes": {"i", "j", "k", "l"}},
    "T4": {"name": "tiny", "genes": {"a", "e"}},
    "T5": {"name": "huge", "genes": set("abcdefghijklm")},
}
UNIVERSE = set("abcdefghijklmnop")
LIST = {"a", "b", "c", "e"}
res = mo._per_pathway_pvalues(TERMS, LIST, UNIVERSE, 3, 10)
check(set(res["terms"]) == {"T1", "T2", "T3"},
      f"the size gate drops T4 (K=2 < 3) and T5 (K=13 > 10) (got {sorted(res['terms'])})")
check(res["terms"].get("T3", {}).get("k") == 0 and res["terms"].get("T3", {}).get("p") == 1.0,
      "SENSITIVITY (rule 1): a zero-hit pathway is KEPT with p = 1.0")
# Genomics' run_ora_for_database builds the reference set as the union of EVERY
# pathway's genes (including ones the size gate then drops) and hands that in.
ALL_ANNOTATED = set().union(*[v["genes"] for v in TERMS.values()])
old = gm._run_fisher_for_namespace(TERMS, LIST, ALL_ANNOTATED, 0.05, False, 3, 10)
check("T3" not in set(old["table"]["Term ID"]),
      "...whereas Genomics' own single-omics engine drops that same pathway (the contrast this rule exists for)")
for tid in ("T1", "T2"):
    row = old["table"][old["table"]["Term ID"] == tid].iloc[0]
    check(close(res["terms"][tid]["p"], row["p-value"]),
          f"SENSITIVITY: {tid}'s p-value equals Genomics' own _run_fisher_for_namespace (got {res['terms'][tid]['p']:.6g})")
check(res["N"] == 13 and res["n"] == 4,
      f"N = reference set restricted to database-annotated genes, n likewise (got N={res['N']}, n={res['n']})")
check(res["terms"]["T1"]["overlap"] == {"a", "b", "c"}, "overlap is the term's members that are in the list")
as_list = {k: {"name": v["name"], "genes": sorted(v["genes"])} for k, v in TERMS.items()}
res_l = mo._per_pathway_pvalues(as_list, LIST, UNIVERSE, 3, 10)
check({k: v["p"] for k, v in res_l["terms"].items()} == {k: v["p"] for k, v in res["terms"].items()},
      "REGRESSION: Transcriptomics' symbol sets carry 'genes' as a sorted LIST - handled identically to a set")
check(mo._per_pathway_pvalues({}, LIST, UNIVERSE, 3, 10)["terms"] == {}
      and mo._per_pathway_pvalues(TERMS, set(), UNIVERSE, 3, 10)["n"] == 0,
      "an empty database / an empty list yields no terms / n = 0 without raising")

# Transcriptomics' own engine, same input
tx = tm._run_ora_for_sets({k: v for k, v in TERMS.items()}, UNIVERSE, LIST, 3, 10, None)
for tid in ("T1", "T2"):
    row = tx["table"][tx["table"]["Term ID"] == tid].iloc[0]
    check(close(res["terms"][tid]["p"], row["p-value"], tol=1e-7),
          f"SENSITIVITY: {tid}'s p-value equals Transcriptomics' own _run_ora_for_sets "
          f"(hypergeometric == one-sided Fisher; got {res['terms'][tid]['p']:.6g} vs {row['p-value']:.6g})")
# Metabolomics' / Proteomics' pair: _drop_zero_overlap then _run_fisher_for_namespace
universe_keys = ALL_ANNOTATED  # Metabolomics' run_ora likewise unions every term before filtering
filtered, _, _ = mm._drop_zero_overlap(TERMS, LIST, universe_keys, 3, 10)
mb = gm._run_fisher_for_namespace(filtered, LIST, universe_keys, 0.05, True, 0, None)
for tid in ("T1", "T2"):
    row = mb["table"][mb["table"]["Term ID"] == tid].iloc[0]
    check(close(res["terms"][tid]["p"], row["p-value"]),
          f"SENSITIVITY: {tid}'s p-value equals Metabolomics'/Proteomics' own path (got {row['p-value']:.6g})")

# ---------------------------------------------------------------------------
print("=== run_joint_pathway_enrichment: partition, combination, ONE FDR pass ===")
A_TERMS = {
    "P1": {"name": "both hit", "genes": {"a1", "a2", "a3", "a4", "a5"}},
    "P2": {"name": "only A hit", "genes": {"a6", "a7", "a8", "a9"}},
    "P3": {"name": "nobody hit", "genes": {"a10", "a11", "a12", "a13"}},
    "P4": {"name": "only in A's database", "genes": {"a1", "a2", "a14", "a15"}},
}
B_TERMS = {
    "P1": {"name": "both hit", "genes": {"b1", "b2", "b3", "b4", "b5"}},
    "P2": {"name": "only A hit", "genes": {"b6", "b7", "b8", "b9"}},
    "P3": {"name": "nobody hit", "genes": {"b10", "b11", "b12", "b13"}},
}
A_UNI = {f"a{i}" for i in range(1, 20)}
B_UNI = {f"b{i}" for i in range(1, 16)}
inputs = {
    "genomics": {"term_map": A_TERMS, "list_set": {"a1", "a2", "a3", "a6", "a7"}, "universe_set": A_UNI,
                 "labels": {"a1": "GENE1"}, "note": None, "caveat": None},
    "metabolomics": {"term_map": B_TERMS, "list_set": {"b1", "b2", "b3"}, "universe_set": B_UNI,
                     "labels": {"b1": "Pyruvate, total", "b2": "Citrate"}, "note": None, "caveat": None},
}
jr = mo.run_joint_pathway_enrichment(inputs, "KEGG", "fisher", 3, 500)
tbl, ptbl = jr["table"], jr["partial_table"]
a_p = mo._per_pathway_pvalues(A_TERMS, inputs["genomics"]["list_set"], A_UNI, 3, 500)["terms"]
b_p = mo._per_pathway_pvalues(B_TERMS, inputs["metabolomics"]["list_set"], B_UNI, 3, 500)["terms"]
check(set(tbl["Term ID"]) == {"P1", "P2"},
      f"combined table = pathways present in BOTH omics with a hit somewhere (got {sorted(tbl['Term ID'])})")
check(jr["n_terms_no_overlap"] == 1,
      "SENSITIVITY (rule 2): the nobody-hit pathway is counted as no-overlap, never enters the family")
check(list(ptbl["Term ID"]) == ["P4"] and "Metabolomics" in ptbl.iloc[0]["Coverage"],
      "SENSITIVITY (rule 3): a pathway tested in only ONE omics goes to the partial table, naming who lacks it")
check(pd.isna(ptbl.iloc[0]["Metabolomics: p-value"]) and not pd.isna(ptbl.iloc[0]["Genomics: p-value"]),
      "...and the missing omics' columns stay empty - never a synthetic p = 1.0")
p1 = tbl[tbl["Term ID"] == "P1"].iloc[0]
check(close(p1["Combined p-value"], mo.combine_omics_pvalues([a_p["P1"]["p"], b_p["P1"]["p"]])["combined_p"]),
      "P1's combined p-value is Fisher over its two per-omics p-values")
p2 = tbl[tbl["Term ID"] == "P2"].iloc[0]
check(p2["Metabolomics: k"] == 0 and close(p2["Metabolomics: p-value"], 1.0)
      and close(p2["Combined p-value"], mo.combine_omics_pvalues([a_p["P2"]["p"], 1.0])["combined_p"]),
      "SENSITIVITY (rule 1): P2, hit only by Genomics, combines Genomics' p with Metabolomics' p = 1.0")
pv = list(tbl["Combined p-value"])
expected_fdr = gm.benjamini_hochberg(pv)
check(all(close(a, b) for a, b in zip(tbl["Combined FDR"], expected_fdr)),
      "SENSITIVITY (rule 5): Combined FDR is ONE BH pass over the combined p-values "
      f"(m = {len(pv)}, got {[round(x, 6) for x in tbl['Combined FDR']]})")
check(not tbl.columns.duplicated().any() and "Genomics: p-value" in tbl.columns
      and "Metabolomics: Overlapping members" in tbl.columns and "Combined FDR" in tbl.columns,
      "the table carries every per-omics column beside the combined ones")
check(tbl.iloc[0]["Combined FDR"] <= tbl.iloc[-1]["Combined FDR"], "sorted by Combined FDR ascending")
check("Pyruvate, total" in str(tbl[tbl["Term ID"] == "P1"].iloc[0]["Metabolomics: Overlapping members"]),
      "a compound name containing a comma survives intact in the members column ('; '-joined)")
check("GENE1" in str(p1["Genomics: Overlapping members"]),
      "members show the user's own label, not the internal key")
check(jr["omics_combined"] == ["genomics", "metabolomics"] and jr["omics_excluded"] == {},
      "both omics were testable and combined")
check(mo.JOINT_CAVEAT.startswith("Per-omics p-values are combined as if independent"),
      "the independence caveat is stated")
jt = mo.run_joint_pathway_enrichment(inputs, "KEGG", "tippett", 3, 500)
check(all("Tippett" in m for m in jt["table"]["Method"]), "the chosen method is recorded on every row")
check(close(jt["table"].set_index("Term ID").loc["P1", "Combined p-value"],
            1 - (1 - min(a_p["P1"]["p"], b_p["P1"]["p"])) ** 2),
      "Tippett's result equals 1-(1-min p)^2")

print("=== run_joint_pathway_enrichment: untestable omics and too few omics ===")
no_hit = dict(inputs)
no_hit["metabolomics"] = dict(inputs["metabolomics"], list_set={"zz1", "zz2"})
nj = mo.run_joint_pathway_enrichment(no_hit, "KEGG", "fisher", 3, 500)
check(nj["table"].empty and "metabolomics" in nj["omics_excluded"] and nj["note"] and "Fewer than two" in nj["note"],
      "SENSITIVITY (rule 4): an omics with no annotated hit is left out; one omics left -> empty table + a note")
three = dict(inputs)
three["proteomics"] = {"term_map": A_TERMS, "list_set": {"a1", "a2", "a6"}, "universe_set": A_UNI,
                       "labels": None, "note": None, "caveat": None}
t3 = mo.run_joint_pathway_enrichment(three, "KEGG", "fisher", 3, 500)
check(t3.get("omics_combined") == ["genomics", "proteomics", "metabolomics"]
      and set(t3["table"]["N omics combined"]) == {3},
      "three omics combine with 2k = 6 degrees of freedom (N omics combined = 3)")
dropped = dict(three)
dropped["proteomics"] = dict(three["proteomics"], note="Could not reach KEGG; re-run to retry.")
d3 = mo.run_joint_pathway_enrichment(dropped, "KEGG", "fisher", 3, 500)
check("proteomics" in d3["omics_excluded"] and d3["omics_combined"] == ["genomics", "metabolomics"]
      and "left out" in (d3["note"] or "") and not d3["table"].empty,
      "an omics that came back with a note is excluded with that reason shown; the rest still combine")
try:
    mo.run_joint_pathway_enrichment({"genomics": inputs["genomics"]}, "KEGG")
    check(False, "a single omics is rejected")
except ValueError:
    check(True, "fewer than two omics selected is rejected outright")

# ---------------------------------------------------------------------------
print("=== per-omics input builders: each calls that omics' OWN builders ===")
calls = {}


def _fake_reactome(taxon_id, taxon_name="", include_iea=True, verbose=False):
    calls["reactome"] = (taxon_id, taxon_name)
    return {"pathways": {"R-1": {"name": "Reactome one", "genes": {"U1", "U2", "U3"}}}, "caveat": None}


def _fake_kegg(org, verbose=False, force_refresh=False):
    calls.setdefault("kegg", []).append(org)
    return {"pathways": {"hsa1": {"name": "Kegg one", "genes": {"hsa:1", "hsa:2", "hsa:3"}}}, "caveat": None}


mo.build_reactome_universe, mo.build_kegg_universe = _fake_reactome, _fake_kegg
gdf = pd.DataFrame({"Gene Symbol": ["G1", "G2", "G3"], "UniProt ID": ["U1", "U2", ""],
                    "KEGG ID": ["hsa:1", "", "hsa:3"]})
gr = mo.build_genomics_input(gdf, "Reactome", "9606", "Homo sapiens", "hsa")
check(gr["list_set"] == {"U1", "U2"} and gr["labels"]["U1"] == "G1" and calls["reactome"] == ("9606", "Homo sapiens"),
      "Genomics + Reactome: matches on UniProt ID, labels back to the user's symbol, organism passed through")
gk = mo.build_genomics_input(gdf, "KEGG", "9606", "Homo sapiens", "hsa")
check(gk["list_set"] == {"hsa:1", "hsa:3"} and calls["kegg"] == ["hsa"] and gk["universe_set"] == {"hsa:1", "hsa:2", "hsa:3"},
      "Genomics + KEGG: matches on KEGG gene ID (not UniProt), universe = the pathway genes, the ONE organism code used")
check("organism name" in mo.build_genomics_input(gdf, "Reactome", "9606", "", "hsa")["note"]
      and "KEGG organism" in mo.build_genomics_input(gdf, "KEGG", "9606", "x", "")["note"],
      "no taxon name (Reactome) / no KEGG organism (KEGG) -> an explanatory note, never a silent empty result")

tm_fetch = {"pathways": {"hsa1": {"name": "Kegg one", "genes": {"1", "2", "3"}}}, "fetch_failed": False}
tm.fetch_kegg_sets = lambda org, verbose=False, force_refresh=False: (calls.setdefault("tm_kegg", []).append(org) or tm_fetch)
tm.fetch_reactome_sets = lambda name, include_iea=True, verbose=False, force_refresh=False: (
    calls.__setitem__("tm_reactome", name) or {"pathways": {"R-1": {"name": "R one", "genes": {"1", "2"}}}})
udf = pd.DataFrame({"gene_symbol": ["Aa", "Bb", "Cc", "Dd"], "entrez_id": ["1", "2", "3", "9"]})
tr = mo.build_transcriptomics_input(udf, {"AA", "BB"}, "KEGG", "Homo sapiens", "hsa")
check(tr["list_set"] == {"AA", "BB"} and tr["universe_set"] == {"AA", "BB", "CC", "DD"}
      and calls["tm_kegg"] == ["hsa"] and set(tr["term_map"]["hsa1"]["genes"]) == {"AA", "BB", "CC"},
      "Transcriptomics: Entrez-keyed pathways bridged to symbols; universe = the genes TESTED; the shared organism code used")
trr = mo.build_transcriptomics_input(udf, {"AA"}, "Reactome", "Homo sapiens", "hsa")
check(calls["tm_reactome"] == "Homo sapiens" and trr["list_set"] == {"AA"},
      "Transcriptomics + Reactome uses the shared taxon name")
nm = mo.build_transcriptomics_input(pd.DataFrame({"gene_symbol": ["Aa"], "entrez_id": [""]}), {"AA"}, "KEGG", "x", "hsa")
check(nm["note"] and not nm["term_map"], "no Entrez ID in common -> an explanatory note")

pm.build_kegg_accession_sets = lambda org, measured, verbose=False, force_refresh=False: (
    calls.setdefault("pm_kegg", []).append(org) or
    {"pathways": {"hsa1": {"name": "K", "genes": {"P1", "P2", "P3"}}}, "fetch_failed": False, "caveat": None})
pdf = pd.DataFrame({"accession": ["P1", "P2", "P3", "P9"], "gene_symbol": ["a", "b", "c", "d"]})
pr = mo.build_proteomics_input(pdf, {"P1", "P2"}, "KEGG", "9606", "Homo sapiens", "hsa")
check(pr["list_set"] == {"P1", "P2"} and pr["universe_set"] == {"P1", "P2", "P3", "P9"}
      and calls["pm_kegg"] == ["hsa"] and pr["labels"]["P1"] == "A",
      "Proteomics: universe = the proteins MEASURED, KEGG via the module's own accession sets with the shared code")

mdf = pd.DataFrame({"compound_key": ["K1", "K2", "K3", ""], "Name": ["Pyruvate", "Citrate", "Lipid", "x"],
                    "class_level_only": [False, False, True, False]})
mm._build_universes = lambda matched, cw, dbs, **kw: (
    calls.__setitem__("mm_kw", kw) or
    {dbs[0]: {"terms": {"hsa1": {"name": "K", "genes": {"K1", "K2", "K4"}}}, "caveat": "scope note"}})
mr = mo.build_metabolomics_input(mdf, {}, "KEGG", "9606", "Homo sapiens", "hsa", "Homo sapiens (human)")
check(mr["list_set"] == {"K1", "K2"} and mr["universe_set"] == {"K1", "K2", "K4"}
      and calls["mm_kw"]["universe_scope"] == "database" and calls["mm_kw"]["kegg_organism_code"] == "hsa",
      "Metabolomics: list = pathway-eligible compounds (class-level-only lipids excluded), universe = the database's "
      "own annotated set, the shared KEGG organism code used")
check(mr["labels"]["K1"] == "Pyruvate" and mr["caveat"] == "scope note",
      "...labels are the user's own compound names; the universe caveat travels along")

print("=== joint_result_columns / helpers ===")
cols = mo.joint_result_columns(["genomics", "metabolomics"])
check(cols[:3] == ["Term ID", "Term Name", "Database"] and cols[-1] == "Coverage" and "Genomics: k" in cols,
      "column order: ids, per-omics blocks in OMICS_TYPES order, combined columns last")
check(mo.finite_or_none(float("nan")) is None and mo.finite_or_none(2) == 2.0 and mo.finite_or_none("x") is None,
      "finite_or_none maps NaN / junk to None")

# ===========================================================================
# BROWN'S METHOD - Fisher corrected for CORRELATED omics (the default)
#
# Fisher's method assumes independent p-values. RNA and protein from the same samples are not independent, so
# Fisher counts one signal twice and the joint p-value comes out too small. Brown (1975) matches the mean and
# variance of Fisher's statistic to a scaled chi-squared; the "empirical" version (Poole et al. 2016) estimates the
# correlation from the data itself. Properties asserted below: rho = 0 IS Fisher; an omics combined with a copy of
# itself returns its own p; a positive rho never makes the result MORE significant; negative rho is clamped.
# ===========================================================================
print("=== brown_combine: the arithmetic, against an independent chi-squared formula ===")
ps3 = [1e-6, 1e-3, 0.2]
rhos3 = [0.5, 0.2, 0.8]
X3 = -2 * sum(math.log(p) for p in ps3)
var3 = 12 + 8 * sum(rhos3)
c3, f3 = var3 / 12, 2 * 36 / var3
b3 = mo.brown_combine(ps3, rhos3)
check(close(b3["combined_p"], chi2.sf(X3 / c3, f3)) and close(b3["c"], c3) and close(b3["f"], f3) and close(b3["statistic"], X3),
      f"three omics, hand-derived: Var = 4k + 8*sum(rho) = {var3}, c = {c3:.3f}, f = {f3:.3f}, p = {b3['combined_p']:.4g}")
check(close(mo.brown_combine(ps3, [0.8, 0.5, 0.2])["combined_p"], b3["combined_p"]),
      "only the SUM of the pairwise correlations matters (order is irrelevant)")
for k in (2, 3, 4, 5):
    pk = [10 ** -(i + 2) for i in range(k)]
    zero = mo.brown_combine(pk, [0.0] * (k * (k - 1) // 2))
    fish = mo.combine_omics_pvalues(pk, "fisher")
    check(close(zero["combined_p"], fish["combined_p"]) and close(zero["c"], 1.0) and close(zero["f"], 2 * k),
          f"SENSITIVITY: with every rho = 0, Brown IS Fisher (k={k}: c = 1, f = 2k = {2 * k}, p identical)")
cp = mo.brown_combine([1e-8, 1e-8], [1.0])
check(close(cp["combined_p"], 1e-8) and close(cp["c"], 2.0) and close(cp["f"], 2.0),
      f"SENSITIVITY: an omics combined with an exact copy of itself (rho = 1) returns its OWN p, 1e-8 (c = 2, f = 2) - "
      f"got {cp['combined_p']:.3g}")
check(mo.combine_omics_pvalues([1e-8, 1e-8], "fisher")["combined_p"] < 1e-14,
      "...whereas Fisher claims ~4e-15 for that same duplicated evidence - the double counting Brown removes")
import random as _random  # noqa: E402
_rng = _random.Random(7)
worst = 1.0
for _ in range(400):
    k = _rng.choice([2, 3, 4])
    pk = [10 ** -_rng.uniform(0, 12) for _ in range(k)]
    rk = [_rng.uniform(0, 1) for _ in range(k * (k - 1) // 2)]
    fis = mo.combine_omics_pvalues(pk, "fisher")["combined_p"]
    bro = mo.brown_combine(pk, rk)["combined_p"]
    if fis <= 0.28 and fis > 1e-250 and bro < fis * (1 - 1e-9):
        check(False, f"Brown more significant than Fisher where it matters: {pk} {rk} {bro} < {fis}")
        break
else:
    check(True, "SENSITIVITY: across 400 random cases with rho in [0, 1], wherever Fisher's p <= 0.28 (every result that could "
                "be significant) Brown's p is NEVER smaller than Fisher's")
# The claim is NOT "everywhere": past ~0.28 the tails cross. State that honestly and pin where.
_cross = [(p1, p2) for p1 in (1.0, 0.9, 0.5, 0.3) for p2 in (1.0, 0.9, 0.5, 0.3)
          if mo.brown_combine([p1, p2], [0.9])["combined_p"] < mo.combine_omics_pvalues([p1, p2], "fisher")["combined_p"] * (1 - 1e-9)]
check(len(_cross) > 0 and all(mo.combine_omics_pvalues(list(c), "fisher")["combined_p"] > 0.28 for c in _cross),
      f"HONESTY: the crossover is real ({len(_cross)} grid cases where Brown < Fisher at rho = 0.9) but ONLY where Fisher's p > 0.28 - "
      "between two non-significant p-values")
_lo = [(a, b, r) for a in (1e-2, 1e-3, 1e-5) for b in (0.5, 0.05, 1e-4) for r in (0.1, 0.5, 1.0)
       if mo.combine_omics_pvalues([a, b], "fisher")["combined_p"] <= 0.28]
check(all(mo.brown_combine([a, b], [r])["combined_p"] >= mo.combine_omics_pvalues([a, b], "fisher")["combined_p"] * (1 - 1e-9) for a, b, r in _lo),
      f"...and below 0.28 Brown >= Fisher holds on a dense grid ({len(_lo)} cases)")
check(all(0.0 <= mo.brown_combine([_rng.random() for _ in range(3)], [0.3, 0.3, 0.3])["combined_p"] <= 1.0 for _ in range(50)),
      "the combined p-value is always a probability")
try:
    mo.brown_combine([0.1, 0.2, 0.3], [0.5])
    check(False, "a wrong number of correlations is refused")
except ValueError:
    check(True, "k omics need exactly k(k-1)/2 correlations - anything else is refused")
try:
    mo.combine_omics_pvalues([0.1, 0.2], "brown")
    check(False, "Brown without correlations is refused")
except ValueError:
    check(True, "method 'brown' without correlations is refused (never silently Fisher)")
check(close(mo.combine_omics_pvalues([0.01, 0.04], "brown", [0.4])["combined_p"], mo.brown_combine([0.01, 0.04], [0.4])["combined_p"]),
      "combine_omics_pvalues dispatches 'brown'")
check(math.isfinite(mo.brown_combine([0.0, 0.5], [0.3])["statistic"]), "a p-value that underflowed to 0.0 is floored")

print("=== _pearson / estimate_omics_correlation ===")
import numpy as _np  # noqa: E402
_x, _y = _np.random.RandomState(1).rand(50), _np.random.RandomState(2).rand(50)
def _pure_pearson(x, y):
    """Reference Pearson in plain Python. numpy.corrcoef / numpy.cov are NOT usable as the reference here: they hard-crash
    the interpreter in this environment (a broken BLAS call - the same bug as scipy.stats.spearmanr, Windows fatal
    exception 0xc06d007f), which is exactly why the engine's _pearson uses only elementwise numpy."""
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))


check(close(mo._pearson(_x, _y), _pure_pearson(list(_x), list(_y))), "hand-written Pearson equals an independent pure-Python Pearson")
check(close(mo._pearson(_x, 3 * _x + 1), 1.0) and close(mo._pearson(_x, -2 * _x), -1.0), "perfect positive / negative linear relations give +1 / -1")
check(math.isnan(mo._pearson([1, 1, 1], [1, 2, 3])), "a constant vector has no correlation (NaN, handled by the caller)")


def _per_from_X(**xs):
    """per[o]['terms'][t]['p'] from X = -2 ln p vectors."""
    n = len(next(iter(xs.values())))
    return {o: {"terms": {f"P{i:03d}": {"p": math.exp(-v / 2)} for i, v in enumerate(vec)}} for o, vec in xs.items()}, \
        [f"P{i:03d}" for i in range(n)]


_g = _np.random.RandomState(11)
_n = 300
_s, _a, _b = _g.exponential(size=_n), _g.exponential(size=_n), _g.exponential(size=_n)
per_c, ids_c = _per_from_X(genomics=list(3 * _s + _a), transcriptomics=list(3 * _s + _b))
est = mo.estimate_omics_correlation(per_c, ["genomics", "transcriptomics"], ids_c)
check(0.85 <= est["pairs"][0]["rho"] <= 0.95 and not est["fallback"],
      f"two omics sharing a signal (true correlation 0.90) are estimated at {est['pairs'][0]['rho']:.2f}")
per_i, ids_i = _per_from_X(genomics=list(_a), transcriptomics=list(_b))
est_i = mo.estimate_omics_correlation(per_i, ["genomics", "transcriptomics"], ids_i)
check(abs(est_i["pairs"][0]["rho"]) < 0.15, f"two independent omics are estimated near 0 ({est_i['pairs'][0]['rho']:.2f})")
per_n, ids_n = _per_from_X(genomics=list(_a), transcriptomics=list(_np.maximum(_a.max() - _a, 0)))
est_n = mo.estimate_omics_correlation(per_n, ["genomics", "transcriptomics"], ids_n)
check(est_n["pairs"][0]["rho_raw"] < -0.5 and est_n["pairs"][0]["rho"] == 0.0 and any("negative" in x for x in est_n["notes"]),
      f"SENSITIVITY: a NEGATIVE estimate ({est_n['pairs'][0]['rho_raw']:.2f}) is clamped to 0 and says so - otherwise Brown "
      "could come out MORE significant than Fisher")
per_s, ids_s = _per_from_X(genomics=list(_a[:19]), transcriptomics=list(_a[:19] * 2))
est_s = mo.estimate_omics_correlation(per_s, ["genomics", "transcriptomics"], ids_s)
check(est_s["fallback"] and est_s["rhos"] == [0.0] and any("Only 19" in x for x in est_s["notes"]),
      f"SENSITIVITY: fewer than {mo.BROWN_MIN_PATHWAYS} pathways (19) -> rho = 0 (plain Fisher) WITH a note, even for "
      "perfectly correlated data")
per_20, ids_20 = _per_from_X(genomics=list(_a[:20]), transcriptomics=list(_a[:20] * 2))
check(not mo.estimate_omics_correlation(per_20, ["genomics", "transcriptomics"], ids_20)["fallback"],
      "...and exactly 20 pathways is enough")
per_k, ids_k = _per_from_X(genomics=[0.0] * 40, transcriptomics=list(_a[:40]))
est_k = mo.estimate_omics_correlation(per_k, ["genomics", "transcriptomics"], ids_k)
check(est_k["rhos"] == [0.0] and any("no variation" in x for x in est_k["notes"]),
      "an omics with no variation (p = 1 everywhere) cannot correlate: rho = 0 with a note, not a crash")
per_3, ids_3 = _per_from_X(genomics=list(3 * _s + _a), transcriptomics=list(3 * _s + _b), proteomics=list(_g.exponential(size=_n)))
est_3 = mo.estimate_omics_correlation(per_3, ["genomics", "transcriptomics", "proteomics"], ids_3)
check([(p["a"], p["b"]) for p in est_3["pairs"]] == [("genomics", "transcriptomics"), ("genomics", "proteomics"), ("transcriptomics", "proteomics")]
      and len(est_3["rhos"]) == 3, "three omics give three pairs in a fixed order")

print("=== run_joint_pathway_enrichment under Brown ===")


def _mk_inputs(n_path, shared, seed=5):
    """Two omics over n_path pathways; omics B's hits track omics A's when shared > 0."""
    g = _np.random.RandomState(seed)
    # A LARGE gene universe on purpose: with few genes almost every pathway overlaps SOME hit by chance, so "pathways with
    # a hit" and "all pathways" coincide and a test cannot tell a correlation estimated on all of them from one estimated
    # only on the hit ones (a mutant doing the latter survived the first version of this fixture).
    universe_a = {f"a{i}" for i in range(4000)}
    universe_b = {f"b{i}" for i in range(4000)}
    tm_a, tm_b, hit_a, hit_b = {}, {}, set(), set()
    # sorted(), NEVER list(set): a set of strings iterates in an order that depends on PYTHONHASHSEED, which made which
    # genes were 'hit' differ between runs and this fixture fail about one run in twenty-five.
    for t in range(n_path):
        ma = {f"a{i}" for i in g.choice(4000, 12, replace=False)}
        mb = {f"b{i}" for i in g.choice(4000, 12, replace=False)}
        tm_a[f"T{t}"] = {"name": f"Pathway {t}", "genes": ma}
        tm_b[f"T{t}"] = {"name": f"Pathway {t}", "genes": mb}
        if t % 4 == 0:                       # a genuinely hit pathway
            hit_a |= set(sorted(ma)[:6])
            hit_b |= set(sorted(mb)[:6 if shared else 0])
        elif t % 7 == 0:
            hit_b |= set(sorted(mb)[:4])
    return {"genomics": {"term_map": tm_a, "list_set": hit_a, "universe_set": universe_a, "labels": None, "note": None, "caveat": None},
            "metabolomics": {"term_map": tm_b, "list_set": hit_b, "universe_set": universe_b, "labels": None, "note": None, "caveat": None}}


inp = _mk_inputs(60, shared=True)
rb = mo.run_joint_pathway_enrichment(inp, "KEGG", "brown", 5, 500)
rf = mo.run_joint_pathway_enrichment(inp, "KEGG", "fisher", 5, 500)
cb = rb["correlation"]
check(cb is not None and cb["applied"] and cb["c"] is not None and cb["f"] is not None and len(cb["pairs"]) == 1,
      f"Brown's block reports the correlation it used (rho = {cb['pairs'][0]['rho']:.2f}, c = {cb['c']:.2f}, f = {cb['f']:.2f})")
_n_hit = len(rb["table"]) + rb["n_terms_no_overlap"] * 0
check(cb["n_pathways"] == 60 and len(rb["table"]) < 40,
      f"the correlation is estimated over EVERY jointly tested pathway ({cb['n_pathways']} of 60) although only {len(rb['table'])} "
      "have a hit anywhere (selecting hit pathways would distort it)")
_per_chk = {o: mo._per_pathway_pvalues(inp[o]["term_map"], inp[o]["list_set"], inp[o]["universe_set"], 5, 500) for o in inp}
_ids_chk = sorted(set(_per_chk["genomics"]["terms"]) & set(_per_chk["metabolomics"]["terms"]))
_xa = [-2 * math.log(max(_per_chk["genomics"]["terms"][t]["p"], 1e-300)) for t in _ids_chk]
_xb = [-2 * math.log(max(_per_chk["metabolomics"]["terms"][t]["p"], 1e-300)) for t in _ids_chk]
check(len(_ids_chk) == 60 and close(cb["pairs"][0]["rho_raw"], _pure_pearson(_xa, _xb)),
      f"SENSITIVITY: rho is the Pearson correlation of -2 ln p over ALL 60 jointly tested pathways, recomputed independently "
      f"(engine {cb['pairs'][0]['rho_raw']:.4f}, independent {_pure_pearson(_xa, _xb):.4f})")
check(rf["correlation"] is not None and not rf["correlation"]["applied"] and rf["correlation"]["c"] is None
      and abs(rf["correlation"]["pairs"][0]["rho"] - cb["pairs"][0]["rho"]) < 1e-12,
      "Fisher's block still REPORTS the correlation (information - the dependence it ignores) but does not apply it")
fb = rf["table"].set_index("Term ID")
bb = rb["table"].set_index("Term ID")
check(set(fb.index) == set(bb.index) and len(bb) > 5, "Brown and Fisher combine the same pathways")
check(all(bb.loc[t, "Combined p-value"] >= fb.loc[t, "Combined p-value"] * (1 - 1e-9) for t in bb.index if fb.loc[t, "Combined p-value"] <= 0.28),
      "SENSITIVITY: on every pathway where Fisher's p <= 0.28 (the only ones that could be significant), Brown's p is at least Fisher's")
rho_used = cb["pairs"][0]["rho"]
check(all(close(bb.loc[t, "Combined p-value"], mo.brown_combine([bb.loc[t, "Genomics: p-value"], bb.loc[t, "Metabolomics: p-value"]], [rho_used])["combined_p"])
          for t in bb.index),
      "every row's combined p-value is Brown's formula applied to ITS per-omics p-values and the one estimated rho")
check(all(close(bb.loc[t, "Statistic"], -2 * (math.log(max(bb.loc[t, "Genomics: p-value"], 1e-300)) + math.log(max(bb.loc[t, "Metabolomics: p-value"], 1e-300))))
          for t in bb.index), "the Statistic column stays Fisher's X for every method")
if cb["pairs"][0]["rho"] > 0.05:
    check(sum(bb["Combined p-value"] > fb["Combined p-value"] * 1.0001) > 0,
          "with a real positive correlation Brown is strictly more conservative somewhere")
check(all(close(a, b) for a, b in zip(bb["Combined FDR"], gm.benjamini_hochberg(list(bb["Combined p-value"])))),
      "SENSITIVITY: Brown's FDR is still ONE Benjamini-Hochberg pass over the combined p-values")
check("Brown" in bb["Method"].iloc[0] and set(bb["N omics combined"]) == {2}, "the method is named on every row")

small = _mk_inputs(10, shared=True)
sb = mo.run_joint_pathway_enrichment(small, "KEGG", "brown", 5, 500)
sf = mo.run_joint_pathway_enrichment(small, "KEGG", "fisher", 5, 500)
check(sb["correlation"]["fallback"] and "fewer than 20" in (sb["note"] or ""),
      "SENSITIVITY: with only 10 jointly tested pathways Brown falls back to plain Fisher AND the block's note says so")
check(len(sb["table"]) == len(sf["table"]) and all(close(a, b) for a, b in zip(sb["table"]["Combined p-value"], sf["table"]["Combined p-value"])),
      "...and the fallback result is exactly Fisher's")
check(rb["correlation"]["effective_omics"] is not None and 1.0 <= rb["correlation"]["effective_omics"] <= 2.0 + 1e-9,
      "the effective number of independent omics lies between 1 and the number combined")
check(mo.run_joint_pathway_enrichment({"genomics": inp["genomics"], "metabolomics": dict(inp["metabolomics"], list_set={"zz"})},
                                      "KEGG", "brown", 5, 500)["correlation"] is None,
      "fewer than two usable omics -> no correlation block (nothing to estimate)")
check(mo.COMBINE_METHOD_DEFAULT == "brown" and mo.COMBINE_METHODS == ("brown", "fisher", "tippett"),
      "Brown is the default method; Fisher and Tippett remain; Stouffer is still absent")
check("Brown" in mo.COMBINE_METHOD_LABELS["brown"] and "independent" in mo.COMBINE_METHOD_LABELS["fisher"],
      "labels say which method corrects for correlation and which assumes independence")
check(mo.joint_caveat("brown") == mo.JOINT_CAVEAT_BROWN and mo.joint_caveat("fisher") == mo.JOINT_CAVEAT
      and mo.joint_caveat("tippett") == mo.JOINT_CAVEAT
      and "Per-omics p-values are combined as if independent" not in mo.JOINT_CAVEAT_BROWN
      and "somewhat conservative" in mo.JOINT_CAVEAT_BROWN,
      "each method carries the caveat that is TRUE for it (Brown must not claim independence)")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
