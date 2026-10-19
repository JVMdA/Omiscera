"""Multi-omics - Joint Pathway Enrichment (first tile of the Multi-omics submodule).

QUESTION ANSWERED
-----------------
Two or more omics (Genomics, Transcriptomics, Proteomics, Metabolomics) each
give a list of significant features. For every KEGG / Reactome pathway: do the
omics TOGETHER point at it? Each omics is tested on its own against that
pathway's own members (a one-sided Fisher exact test, the same math every
single-omics ORA tile here uses), and the per-omics p-values are then combined
into ONE joint p-value per pathway.

METHOD, AND WHERE IT COMES FROM
-------------------------------
Default: Brown's method (Brown 1975, Biometrics 31:987; empirical version of
Poole et al. 2016, Bioinformatics 32:i430). It starts from Fisher's (1932)
statistic X = -2 * sum(ln p_i) but corrects for correlated omics (e.g. RNA and
protein from the same samples), which plain Fisher counts twice: with k omics
and pairwise correlations rho_ij, Var[X] = 4k + 8*sum(rho_ij), c = Var/4k,
f = 8k^2/Var, p = chi2.sf(X/c, f). rho_ij is the Pearson correlation of
-2*ln(p) between two omics across all jointly tested pathways (negative values
clamped to 0; under 20 pathways it falls back to rho = 0, i.e. Fisher, with a
note). Plain Fisher (chi-squared, 2k df, assumes independence) stays available.
Merging per-dataset p-values this way is also how ActivePathways combines omics
(Paczkowska et al. 2020, Nat Commun 11:735). It is deliberately NOT the
"pool every hit into one bigger contingency table" approach: combining each
omics' own p-value keeps every omics' evidence separate and visible next to
the joint result, which a pooled count would hide. The code here is Omiscera's
own implementation of those general statistical methods (the combination
itself uses scipy's chi-squared survival function (Fisher/Tippett via
scipy.stats.combine_pvalues or closed form; Brown via chi2.sf), so no
distribution is re-written by hand; contrast with
Metabolomics_module._brandes_betweenness, which is hand-written only because
no graph library exists in this environment).

SECOND METHOD: Tippett (1931), the minimum-p test, answers a different
question ("is AT LEAST ONE omics strongly enriched here?") where Fisher answers
"does evidence accumulate across omics?". Stouffer's Z method is deliberately
NOT offered: it treats a p-value near 1 as strong evidence AGAINST the
pathway (z = -inf at p = 1), but a one-sided over-representation p-value of
1.0 only means "no enrichment seen" - and an omics with no hits in a pathway
legitimately has p = 1.0. Measured: one omics at p = 1e-10 plus one at p = 1.0
combines to p = 0.90 under Stouffer (a real hit wiped out) versus 2.4e-9 under
Fisher. Fisher and Tippett are both unaffected by p = 1.

DESIGN RULES (each one a decision, tested in tests/test_multiomics_joint.py)
----------------------------------------------------------------------------
1. A pathway where an omics has ZERO hits keeps that omics' p-value as 1.0
   in the combination (_per_pathway_pvalues keeps k == 0 terms). Single-omics
   tiles drop k == 0 terms from their own FDR family (clusterProfiler
   convention) and that is right for a standalone table - but "this omics saw
   nothing here" is real evidence when another omics did, so dropping it would
   silently delete it.
2. A pathway enters the COMBINED family only if at least one omics has k > 0.
   A pathway nobody hit says nothing, and counting it would inflate the BH
   family exactly like the k == 0 inflation the single-omics modules already
   documented fixing.
3. A pathway is combined only when it passes the term-size gate in EVERY
   testable omics. A pathway present in only some omics is reported in a
   separate "partial coverage" table and never padded with a synthetic p = 1.0.
4. An omics with NONE of its hits annotated in a database is "not testable"
   there (every p-value would be undefined, not 1.0) and is left out of that
   database's combination, with a note. Fewer than two testable omics -> an
   empty combined table with a note, never a one-omics "combination".
5. ONE Benjamini-Hochberg pass, over the combined p-values only. The per-omics
   p-values shown beside it are informational and not corrected again: BH is a
   property of one family of tests screened together, and that family here is
   "the pathways jointly tested".
6. The KEGG organism is resolved ONCE by the caller (the shared
   Transcriptomics_module.resolve_kegg_organism_code) and the SAME code is
   handed to every omics' KEGG builder - never resolved twice (commit 27ce033
   was exactly two code paths resolving the organism independently and
   disagreeing). Reactome is aligned through one shared taxon name.

ISOLATION
---------
This module only IMPORTS from the four omics engines; nothing in them is
edited. Their universe builders and statistics primitives are called exactly as
their own ORA entry points call them.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from Genomics_module import (
    benjamini_hochberg,
    bonferroni,
    build_kegg_universe,
    build_reactome_universe,
    fisher_exact_pvalue,
    fold_enrichment,
    format_term_name,
)

try:  # pragma: no cover - scipy is a hard dependency of the single-omics engines too
    from scipy.stats import chi2 as _chi2
    from scipy.stats import combine_pvalues as _scipy_combine_pvalues
except Exception:  # pragma: no cover
    _chi2 = None
    _scipy_combine_pvalues = None


OMICS_TYPES = ("genomics", "transcriptomics", "proteomics", "metabolomics")
OMICS_LABELS = {
    "genomics": "Genomics",
    "transcriptomics": "Transcriptomics",
    "proteomics": "Proteomics",
    "metabolomics": "Metabolomics",
}
JOINT_DATABASES = ("KEGG", "Reactome")
COMBINE_METHODS = ("brown", "fisher", "tippett")
COMBINE_METHOD_DEFAULT = "brown"
COMBINE_METHOD_LABELS = {
    "brown": "Brown's method (Fisher, corrected for correlated omics)",
    "fisher": "Fisher's combined probability test (assumes independent omics)",
    "tippett": "Tippett's minimum-p test",
}
# Below this many jointly tested pathways the correlation between two omics' p-values is too noisy to estimate, and
# Brown's method falls back to rho = 0 (plain Fisher) WITH a note - never silently.
BROWN_MIN_PATHWAYS = 20
MIN_TERM_SIZE_DEFAULT = 3
MAX_TERM_SIZE_DEFAULT = 500
FDR_ALPHA_DEFAULT = 0.05
_P_FLOOR = 1e-300  # keeps ln(p) finite when a Fisher test underflows to 0.0

JOINT_CAVEAT = (
    "Per-omics p-values are combined as if independent. The omics measure related biology "
    "and share pathways, so treat the combined p-value and FDR as a ranking of pathways, "
    "not as an exact error rate."
)
JOINT_CAVEAT_BROWN = (
    "Per-omics p-values are combined with Brown's method: how strongly the omics agree across pathways is "
    "estimated from this run's own p-values, and correlated omics count as less than one independent piece of "
    "evidence each. That estimate includes real shared biology, so the combined p-value and FDR are somewhat "
    "conservative - treat them as a ranking of pathways, not as an exact error rate."
)


def joint_caveat(method: str) -> str:
    """The statement every result carries about what the combination assumes."""
    return JOINT_CAVEAT_BROWN if method == "brown" else JOINT_CAVEAT


# ---------------------------------------------------------------------------
# Statistics (pure functions - no I/O, unit-testable offline)
# ---------------------------------------------------------------------------

def _per_pathway_pvalues(term_map: dict, list_set: set, universe_set: set,
                         min_term_size: int, max_term_size) -> dict:
    """One-sided Fisher exact p-value for every pathway passing the size gate,
    INCLUDING pathways with zero hits (p = 1.0) - see design rule 1.

    Same contingency convention as Genomics_module._run_fisher_for_namespace:
    N = reference set restricted to genes annotated somewhere in this database,
    n = the list restricted to it, K = the term restricted to it, k = overlap.
    fisher_exact_pvalue / fold_enrichment are imported, not re-implemented.
    """
    # set(...) throughout: Transcriptomics' symbol sets carry "genes" as a
    # sorted LIST (Transcriptomics_module._entrez_sets_to_symbol_sets), the
    # other omics carry sets.
    annotated: set = set()
    for info in term_map.values():
        annotated |= set(info["genes"])
    annotated_universe = set(universe_set) & annotated
    in_universe = set(list_set) & annotated_universe
    N, n = len(annotated_universe), len(in_universe)

    terms: dict = {}
    n_excluded_size = 0
    if N > 0:
        for term_id, info in term_map.items():
            members = set(info["genes"]) & annotated_universe
            K = len(members)
            if K < min_term_size or (max_term_size and K > max_term_size):
                n_excluded_size += 1
                continue
            overlap = members & in_universe
            k = len(overlap)
            terms[term_id] = {
                "name": info.get("name", term_id), "k": k, "K": K, "n": n, "N": N,
                "p": fisher_exact_pvalue(k, n, K, N), "fold": fold_enrichment(k, n, K, N),
                "overlap": overlap,
            }
    return {
        "terms": terms, "N": N, "n": n, "n_list": len(list_set),
        "n_excluded_size": n_excluded_size, "n_terms_in_database": len(term_map),
    }


def brown_combine(pvalues: "list[float]", rhos: "list[float]") -> dict:
    """Brown's (1975) extension of Fisher's method for CORRELATED p-values.

    Fisher's statistic is X = -2 * sum(ln p_i). Under H0 and independence E[X] = 2k and Var[X] = 4k. When the tests
    are positively correlated the variance is larger: Var[X] = 4k + 8 * sum_{i<j} rho_ij (each -2 ln p has variance
    4, so a covariance of 4*rho_ij). Brown matches X's mean and variance to a scaled chi-squared:
        c = Var / (2 E),   f = 2 E^2 / Var,   p = P(chi2_f >= X / c).
    rhos: the k*(k-1)/2 pairwise correlations (any order - only their sum matters). With every rho = 0 this is
    EXACTLY Fisher (c = 1, f = 2k); an omics combined with an exact copy of itself (rho = 1) returns its own p
    (c = 2, f = 2) instead of the doubled evidence Fisher would claim.
    For rho >= 0, Brown's p is at least Fisher's wherever Fisher's p is below ~0.28 (checked numerically for 2-5
    omics and every rho; the smallest crossover was 0.285). Above that the two tails cross and Brown can be slightly
    SMALLER - but that is a difference between two clearly non-significant p-values, so it never changes a conclusion.
    Returns {"statistic" (= Fisher's X), "combined_p", "n_combined", "method", "c", "f"}.
    """
    if _chi2 is None:
        raise RuntimeError("scipy is required for brown_combine")
    ps = [min(max(float(p), _P_FLOOR), 1.0) for p in pvalues]
    if not ps:
        raise ValueError("No p-values to combine.")
    k = len(ps)
    if len(rhos) != k * (k - 1) // 2:
        raise ValueError(f"Brown's method over {k} p-values needs {k * (k - 1) // 2} pairwise correlations, got {len(rhos)}.")
    stat = -2.0 * sum(math.log(p) for p in ps)
    mean = 2.0 * k
    var = 4.0 * k + 8.0 * sum(float(r) for r in rhos)
    c = var / (2.0 * mean)
    f = 2.0 * mean * mean / var
    combined = float(_chi2.sf(stat / c, f))
    return {"statistic": float(stat), "combined_p": min(max(combined, 0.0), 1.0), "n_combined": k,
            "method": "brown", "c": float(c), "f": float(f)}


def combine_omics_pvalues(pvalues: "list[float]", method: str = "fisher", rhos: "list[float] | None" = None) -> dict:
    """Combine per-omics p-values for ONE pathway.

    method: "brown" (Brown's method - needs `rhos`, the pairwise correlations), "fisher" (X2 = -2*sum(ln p),
    df = 2k, assumes independence) or "tippett" (minimum p). The default here is "fisher" because Brown needs the
    correlations; the user-facing default (COMBINE_METHOD_DEFAULT) is "brown".
    Returns {"statistic", "combined_p", "n_combined", "method"} (+ "c", "f" for Brown).
    """
    if method not in COMBINE_METHODS:
        raise ValueError(f"Unknown combination method {method!r}; use one of {COMBINE_METHODS}.")
    if method == "brown":
        if rhos is None:
            raise ValueError("Brown's method needs the pairwise correlations (rhos).")
        return brown_combine(pvalues, rhos)
    if _scipy_combine_pvalues is None:
        raise RuntimeError("scipy is required for combine_omics_pvalues")
    ps = [min(max(float(p), _P_FLOOR), 1.0) for p in pvalues]
    if not ps:
        raise ValueError("No p-values to combine.")
    stat, combined = _scipy_combine_pvalues(ps, method=method)
    return {
        "statistic": float(stat), "combined_p": float(min(max(combined, 0.0), 1.0)),
        "n_combined": len(ps), "method": method,
    }


def _pearson(x, y) -> float:
    """Plain Pearson correlation (numpy only - the codebase avoids scipy.stats.spearmanr, which crashes in this
    environment). NaN when either vector is constant."""
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    xs, ys = xs - xs.mean(), ys - ys.mean()
    denom = math.sqrt(float((xs * xs).sum()) * float((ys * ys).sum()))
    return float((xs * ys).sum() / denom) if denom > 0 else float("nan")


def estimate_omics_correlation(per: dict, testable: "list[str]", family_ids: "list[str]") -> dict:
    """How strongly each PAIR of omics agree across pathways, estimated from this run's own p-values (the
    "empirical" Brown's method, Poole et al. 2016) - the only way to get it for omics that share no features
    (genes vs metabolites).

    For each omics, the vector over `family_ids` of X = -2 ln p (p floored at 1e-300); rho_ij is the Pearson
    correlation of two such vectors. `family_ids` is every pathway tested in EVERY usable omics - the whole jointly
    tested population, not just the pathways with a hit, because choosing only hit pathways would distort the
    estimate. rho is clamped to [0, 1]: a negative estimate (noise, or two omics pulling different ways) is set to
    0. A negative rho would SHRINK the variance below Fisher's and make Brown claim MORE evidence than independent
    omics - the opposite of what a correction for dependence is for - so wherever Fisher's p is below ~0.28 (every
    result that could be significant) Brown's p is never smaller than Fisher's.
    Fewer than BROWN_MIN_PATHWAYS pathways -> every rho is 0 (plain Fisher) with an explicit note.
    Returns {"pairs": [{"a", "b", "rho", "rho_raw", "n_pathways"}], "rhos": [...], "n_pathways", "fallback", "notes"}.
    """
    ids = sorted(family_ids)
    n = len(ids)
    pairs, notes = [], []
    fallback = n < BROWN_MIN_PATHWAYS
    if fallback:
        notes.append(f"Only {n} pathway(s) were tested in every omics - fewer than {BROWN_MIN_PATHWAYS}, too few to "
                     "estimate how the omics correlate - so the omics are treated as independent (plain Fisher).")
    x = {o: [-2.0 * math.log(min(max(per[o]["terms"][t]["p"], _P_FLOOR), 1.0)) for t in ids] for o in testable}
    for i, a in enumerate(testable):
        for b in testable[i + 1:]:
            raw = float("nan") if fallback else _pearson(x[a], x[b])
            if fallback:
                rho = 0.0
            elif not math.isfinite(raw):
                rho = 0.0
                notes.append(f"{OMICS_LABELS[a]} and {OMICS_LABELS[b]}: one of them has no variation across pathways, so "
                             "their correlation cannot be estimated and is taken as 0.")
            elif raw < 0:
                rho = 0.0
                notes.append(f"{OMICS_LABELS[a]} and {OMICS_LABELS[b]}: estimated correlation {raw:.2f} is negative and "
                             "was set to 0 (a negative value would make Brown's method claim MORE evidence than independent omics).")
            else:
                rho = min(raw, 1.0)
            pairs.append({"a": a, "b": b, "rho": rho, "rho_raw": raw if math.isfinite(raw) else None, "n_pathways": n})
    return {"pairs": pairs, "rhos": [p["rho"] for p in pairs], "n_pathways": n, "fallback": fallback, "notes": notes}


def _members_text(overlap: set, labels: "dict | None") -> str:
    """'; '-joined display labels - compound names contain commas, so a comma
    join (the gene modules' convention) would be ambiguous here."""
    if not overlap:
        return ""
    shown = {(labels or {}).get(m, m) for m in overlap}
    return "; ".join(sorted(str(s) for s in shown if str(s).strip()))


def joint_result_columns(omics: "list[str]") -> "list[str]":
    """Column order for a joint table over the given omics, in OMICS_TYPES order."""
    cols = ["Term ID", "Term Name", "Database"]
    for o in omics:
        label = OMICS_LABELS[o]
        cols += [f"{label}: k", f"{label}: K", f"{label}: p-value", f"{label}: Overlapping members"]
    cols += ["Combined p-value", "Combined FDR", "Combined Bonferroni", "Statistic", "Method",
             "N omics combined", "Coverage"]
    return cols


def run_joint_pathway_enrichment(
    inputs: dict,
    database: str,
    method: str = COMBINE_METHOD_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE_DEFAULT,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
) -> dict:
    """Combine per-omics pathway enrichment for ONE database.

    inputs: {omics_type: {"term_map": {term_id: {"name", "genes": set}},
                           "list_set": set, "universe_set": set,
                           "labels": {member: display} | None,
                           "note": str | None, "caveat": str | None}}
    Returns {"table", "partial_table", "per_omics": {omics: summary},
             "n_terms_combined", "n_terms_partial", "n_terms_no_overlap",
             "omics_combined", "omics_excluded", "note", "caveat", "method",
             "correlation": {"pairs", "n_pathways", "fallback", "applied", "c", "f", "effective_omics", "notes"}
                            | None (None when fewer than two omics were usable)}.
    """
    omics = [o for o in OMICS_TYPES if o in inputs]
    if len(omics) < 2:
        raise ValueError("Joint pathway enrichment needs at least two omics.")

    per: dict = {}
    for o in omics:
        spec = inputs[o]
        per[o] = _per_pathway_pvalues(
            spec.get("term_map") or {}, spec.get("list_set") or set(),
            spec.get("universe_set") or set(), min_term_size, max_term_size)

    notes: list = []
    excluded: dict = {}
    for o in omics:
        spec_note = inputs[o].get("note")
        if spec_note:
            excluded[o] = spec_note
        elif not per[o]["terms"]:
            excluded[o] = (f"{OMICS_LABELS[o]} has no {database} pathways passing the term-size gate "
                           "for this run.")
        elif per[o]["n"] == 0:
            excluded[o] = (f"None of the {OMICS_LABELS[o]} hits are annotated in {database}, so "
                           f"{OMICS_LABELS[o]} cannot be tested here.")
    testable = [o for o in omics if o not in excluded]
    for o, why in excluded.items():
        notes.append(f"{OMICS_LABELS[o]} was left out of the {database} combination: {why}")

    columns = joint_result_columns(omics)
    per_omics_summary = {
        o: {"n_list": per[o]["n_list"], "n_in_universe": per[o]["n"], "universe_size": per[o]["N"],
            "n_terms_tested": len(per[o]["terms"]), "n_excluded_size": per[o]["n_excluded_size"],
            "testable": o in testable}
        for o in omics
    }
    caveats = [c for c in (inputs[o].get("caveat") for o in omics) if c]
    base = {
        "table": pd.DataFrame(columns=columns), "partial_table": pd.DataFrame(columns=columns),
        "per_omics": per_omics_summary, "n_terms_combined": 0, "n_terms_partial": 0,
        "n_terms_no_overlap": 0, "omics_combined": list(testable), "omics_excluded": excluded,
        "note": None, "caveat": " ".join(dict.fromkeys(caveats)) or None,
        "method": method, "columns": columns, "correlation": None,
    }
    if len(testable) < 2:
        base["note"] = " ".join(notes + [
            f"Fewer than two omics could be tested against {database}, so no joint result is produced."
        ])
        return base

    all_ids = set()
    for o in testable:
        all_ids |= set(per[o]["terms"])
    full_ids = [t for t in all_ids if all(t in per[o]["terms"] for o in testable)]
    partial_ids = [t for t in all_ids if t not in set(full_ids)]

    def _row(term_id: str, present: "list[str]", coverage: str) -> dict:
        name = next((per[o]["terms"][term_id]["name"] for o in present), term_id)
        row = {"Term ID": term_id, "Term Name": format_term_name(name), "Database": database,
               "Coverage": coverage}
        for o in omics:
            label = OMICS_LABELS[o]
            hit = per[o]["terms"].get(term_id) if o in testable else None
            row[f"{label}: k"] = hit["k"] if hit else None
            row[f"{label}: K"] = hit["K"] if hit else None
            row[f"{label}: p-value"] = hit["p"] if hit else None
            row[f"{label}: Overlapping members"] = (
                _members_text(hit["overlap"], inputs[o].get("labels")) if hit else None)
        return row

    # How the omics correlate across pathways: ALWAYS estimated and reported (so a user running plain Fisher can see the
    # dependence it ignores), and used for the combination only under Brown's method.
    correlation = estimate_omics_correlation(per, testable, full_ids)
    rhos = correlation["rhos"] if method == "brown" else None
    if rhos is not None:
        k_all = len(testable)
        var_all = 4.0 * k_all + 8.0 * sum(rhos)
        correlation["c"] = var_all / (4.0 * k_all)
        correlation["f"] = 8.0 * k_all * k_all / var_all
        correlation["effective_omics"] = correlation["f"] / 2.0
    else:
        correlation["c"] = correlation["f"] = correlation["effective_omics"] = None
    correlation["applied"] = method == "brown"

    rows, n_no_overlap = [], 0
    for term_id in full_ids:
        if not any(per[o]["terms"][term_id]["k"] > 0 for o in testable):
            n_no_overlap += 1
            continue
        combo = combine_omics_pvalues([per[o]["terms"][term_id]["p"] for o in testable], method, rhos)
        row = _row(term_id, testable, "all omics")
        row["Combined p-value"] = combo["combined_p"]
        row["Statistic"] = combo["statistic"]
        row["Method"] = COMBINE_METHOD_LABELS[method]
        row["N omics combined"] = combo["n_combined"]
        rows.append(row)

    if rows:
        pvals = [r["Combined p-value"] for r in rows]
        for r, q, b in zip(rows, benjamini_hochberg(pvals), bonferroni(pvals)):
            r["Combined FDR"], r["Combined Bonferroni"] = q, b
    table = pd.DataFrame(rows, columns=columns)
    if not table.empty:
        table = table.sort_values(by=["Combined FDR", "Combined p-value"], kind="stable").reset_index(drop=True)

    partial_rows = []
    for term_id in partial_ids:
        present = [o for o in testable if term_id in per[o]["terms"]]
        if not any(per[o]["terms"][term_id]["k"] > 0 for o in present):
            continue
        missing = [OMICS_LABELS[o] for o in testable if o not in present]
        row = _row(term_id, present, "partial: not tested in " + ", ".join(missing))
        row["N omics combined"] = len(present)
        partial_rows.append(row)
    partial_table = pd.DataFrame(partial_rows, columns=columns)
    if not partial_table.empty:
        best_p = partial_table[[f"{OMICS_LABELS[o]}: p-value" for o in testable]].min(axis=1)
        partial_table = partial_table.loc[best_p.sort_values(kind="stable").index].reset_index(drop=True)

    base.update({
        "table": table, "partial_table": partial_table, "n_terms_combined": len(rows),
        "n_terms_partial": len(partial_rows), "n_terms_no_overlap": n_no_overlap,
        "note": " ".join(notes + (correlation["notes"] if method == "brown" else [])) or None,
        "correlation": correlation,
    })
    return base


# ---------------------------------------------------------------------------
# Per-omics input assembly - calls each omics' OWN universe builders, exactly
# the way that omics' own ORA entry point calls them (read-only imports).
# ---------------------------------------------------------------------------

def _empty_input(note: str, caveat: "str | None" = None) -> dict:
    return {"term_map": {}, "list_set": set(), "universe_set": set(), "labels": None,
            "note": note, "caveat": caveat}


def build_genomics_input(final_df: pd.DataFrame, database: str, taxon_id: str, taxon_name: str,
                         kegg_organism_code: str, include_iea: bool = True,
                         force_refresh: bool = False, verbose: bool = False) -> dict:
    """final_df: Genomics resolved table (UniProt ID, optional KEGG ID, Gene Symbol /
    Ensembl Gene ID), ambiguous / unresolved rows already removed by the caller."""
    id_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else "Ensembl Gene ID"
    if database == "Reactome":
        if not str(taxon_name).strip():
            return _empty_input("Reactome needs an organism name; none was available for this run.")
        universe = build_reactome_universe(taxon_id, taxon_name=taxon_name, include_iea=include_iea,
                                           verbose=verbose)
        key_col = "UniProt ID"
    elif database == "KEGG":
        if not kegg_organism_code:
            return _empty_input("No KEGG organism could be resolved for this run.")
        universe = build_kegg_universe(kegg_organism_code, verbose=verbose, force_refresh=force_refresh)
        key_col = "KEGG ID"
    else:
        raise ValueError(f"Unsupported database {database!r}")
    term_map = universe.get("pathways") or {}
    if not term_map:
        return _empty_input(universe.get("caveat") or f"No {database} pathway data available for this run.")
    if key_col not in final_df.columns:
        return _empty_input(f"The Genomics table has no '{key_col}' column.")
    keys = final_df[key_col].astype(str).str.strip()
    labels = {k: str(v).strip() for k, v in zip(keys, final_df[id_col].astype(str)) if k and str(v).strip()}
    universe_set: set = set()
    for info in term_map.values():
        universe_set |= set(info["genes"])
    return {"term_map": term_map, "list_set": {k for k in keys if k}, "universe_set": universe_set,
            "labels": labels, "note": None, "caveat": universe.get("caveat")}


def build_transcriptomics_input(universe_df: pd.DataFrame, deg_symbols: set, database: str,
                                taxon_name: str, kegg_organism_code: str, include_iea: bool = True,
                                force_refresh: bool = False, verbose: bool = False) -> dict:
    """Symbol-space sets, built exactly as Transcriptomics_module.run_ora builds
    them: the database's own Entrez-keyed pathways bridged to symbols through
    the resolved table's entrez_id column; universe = the genes tested."""
    import Transcriptomics_module as tm  # local: keeps this module importable without it

    universe_symbols = set(universe_df["gene_symbol"].astype(str).str.upper())
    if database == "Reactome":
        if not str(taxon_name).strip():
            return _empty_input("Reactome needs an organism name; none was available for this run.")
        fetched = tm.fetch_reactome_sets(taxon_name, include_iea=include_iea, verbose=verbose,
                                         force_refresh=force_refresh)
    elif database == "KEGG":
        if not kegg_organism_code:
            return _empty_input("No KEGG organism could be resolved for this run.")
        fetched = tm.fetch_kegg_sets(kegg_organism_code, verbose=verbose, force_refresh=force_refresh)
        if fetched.get("fetch_failed"):
            return _empty_input(fetched.get("caveat") or "Could not reach KEGG for this run; re-run to retry.")
    else:
        raise ValueError(f"Unsupported database {database!r}")
    term_map = tm._entrez_sets_to_symbol_sets(fetched.get("pathways"), universe_df)
    if not term_map:
        return _empty_input(
            fetched.get("caveat")
            or f"{database} has no pathway that could be matched to your resolved genes "
               "(no Entrez Gene IDs in common).")
    return {"term_map": term_map, "list_set": set(deg_symbols) & universe_symbols,
            "universe_set": universe_symbols, "labels": None, "note": None,
            "caveat": fetched.get("caveat")}


def build_proteomics_input(resolved_df: pd.DataFrame, sig_accessions: set, database: str,
                           taxon_id: str, taxon_name: str, kegg_organism_code: str,
                           include_iea: bool = True, force_refresh: bool = False,
                           verbose: bool = False) -> dict:
    """Accession-space sets, built as Proteomics_module.run_ora builds them:
    universe = the proteins measured. KEGG goes through the module's own
    build_kegg_accession_sets with the ONE shared organism code."""
    import Proteomics_module as pm

    universe_df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""]
    universe_set = set(universe_df["accession"].astype(str))
    labels = {a: s for a, s in zip(universe_df["accession"].astype(str),
                                   universe_df["gene_symbol"].astype(str).str.upper()) if s.strip()}
    if database == "Reactome":
        if not str(taxon_name).strip():
            return _empty_input("Reactome needs an organism name; none was available for this run.")
        universe = build_reactome_universe(taxon_id, taxon_name=taxon_name, include_iea=include_iea,
                                           verbose=verbose)
        term_map, caveat = universe.get("pathways") or {}, universe.get("caveat")
    elif database == "KEGG":
        if not kegg_organism_code:
            return _empty_input("No KEGG organism could be resolved for this run.")
        kegg = pm.build_kegg_accession_sets(kegg_organism_code, universe_set, verbose=verbose,
                                            force_refresh=force_refresh)
        if kegg.get("fetch_failed"):
            return _empty_input(kegg.get("caveat") or "Could not reach KEGG for this run; re-run to retry.")
        term_map, caveat = kegg.get("pathways") or {}, kegg.get("caveat")
    else:
        raise ValueError(f"Unsupported database {database!r}")
    if not term_map:
        return _empty_input(caveat or f"No {database} pathway data available for this run.")
    return {"term_map": term_map, "list_set": set(sig_accessions) & universe_set,
            "universe_set": universe_set, "labels": labels, "note": None, "caveat": caveat}


def build_metabolomics_input(matched_df: pd.DataFrame, crosswalk: dict, database: str,
                             taxon_id: str, taxon_name: str, kegg_organism_code: str,
                             kegg_organism_name: str = "", force_refresh: bool = False,
                             verbose: bool = False) -> dict:
    """Compound-key-space sets, built as Metabolomics_module.run_ora builds
    them for a plain list: universe = the database's own annotated compounds."""
    import Metabolomics_module as mm

    universes = mm._build_universes(
        matched_df, crosswalk, [database], universe_scope="database", taxon_id=taxon_id,
        taxon_name=taxon_name, kegg_organism_code=kegg_organism_code,
        kegg_organism_name=kegg_organism_name, force_refresh=force_refresh, verbose=verbose)
    universe = universes.get(database) or {}
    terms = universe.get("terms") or {}
    if not terms:
        return _empty_input(universe.get("note") or universe.get("caveat")
                            or f"No {database} annotations available for this run.",
                            caveat=universe.get("caveat") if universe.get("note") else None)
    universe_set: set = set()
    for info in terms.values():
        universe_set |= set(info["genes"])
    return {"term_map": terms, "list_set": mm._eligible_keys(matched_df, database) & universe_set,
            "universe_set": universe_set, "labels": mm._key_to_name_map(matched_df),
            "note": None, "caveat": universe.get("caveat")}


def finite_or_none(value):
    """JSON-safe number: NaN/inf -> None."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None
