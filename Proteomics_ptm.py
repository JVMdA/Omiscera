# -*- coding: utf-8 -*-
"""
Proteomics_ptm.py — PTM / Phosphosite Enrichment engine for Omiscera.

Answers a question none of the other three Proteomics tiles can: given a
SITE-level (not protein-level) phosphoproteomics result, which KINASES look
active, and are my regulated sites enriched for any kinase's known
substrates? Consumes Proteomics_ptm_retriever.resolve_site_table's output.

DATA SOURCE: iPTMnet release 6.2 (research.bioinformatics.udel.edu/iptmnet),
bulk files verified live 2026-09-15:
  ptm.txt   https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ptm.txt
            44 MB, tab-separated, NO header row. Columns in order:
            ptm_type, source, substrate_UniProtAC, substrate_genename,
            organism, site, enzyme_UniProtAC, enzyme_genename, note, pmid.
            ~23,000 enzyme-substrate-site triples across 11 PTM types.
  score.txt https://research.bioinformatics.udel.edu/iptmnet_data/files/current/score.txt
            26 MB, columns: substrate_UniProtAC, site, enzyme_UniProtAC,
            ptm_type, score - a per-relationship confidence, used only as an
            OPTIONAL minimum-confidence filter (default off).
This is currently the ONLY source in this app of kinase->substrate-SITE
relationships; dbPTM (2025, biomics.lab.nycu.edu.tw/dbPTM) has far deeper
site coverage but no clean enzyme table - see the module docstring in
Stage 3 (not yet built) for its narrower, supporting role.

LICENCE: iPTMnet is CC BY-NC-SA 4.0 (NonCommercial, ShareAlike) - stricter
than every other source this app uses (STRING/GO/Reactome/DISEASES/Complex
Portal are all CC BY 4.0 or CC0). This is acceptable for THIS repository
specifically because cache/ is gitignored (~1.4GB, never committed - see
.gitignore's own comment) and every fetcher re-downloads at each user's own
install: Omiscera's GitHub repo ships CODE, not DATA, so it never
redistributes iPTMnet's dataset itself. The NC/ShareAlike terms bind each
user's OWN use of the data they download, not this codebase. Any user
running Omiscera for commercial purposes must independently confirm this
tile's use complies with iPTMnet's licence - stated plainly in
info/README-PTM-proteomics.txt and the in-app caption, not just here.
dbPTM (Stage 3) states no licence at all, only a citation request.

WHY THE REGULON ENGINE'S NUMERICS ARE IMPORTED, NOT GENERALIZED: "kinase ->
substrate sites, scored by those sites' log2FC" is structurally the TF/
Regulon tile's "TF -> target genes" problem (Transcriptomics_regulon.py).
Its pure-numeric core (_ulm_score/_null_slopes_shared/_permutation_p_value/
_empirical_null_calibration/_signature_key) takes only integer positions
into an array plus a centred value vector - nothing gene-shaped - so it is
imported directly rather than hooked or copied. run_regulon_activity ITSELF
is NOT reused/generalized: every kinase->site edge here is +1 (iPTMnet has
no "dephosphorylating kinase" - a phosphatase is a different enzyme entirely,
not a sign on the same edge), so n_neg is always 0 - the signed/unsigned
database split and "Activated/Repressed Targets" columns that function's
prose is built around do not apply. Generalizing it would need ~6 hooks over
a function whose genuinely shared body is a fraction of its length, for a
caller whose payload, figures and prose all differ anyway - the precedent
this avoids over-applying is Transcriptomics_network.run_network_analysis,
which WAS worth generalizing because both callers emit the identical payload
for the identical client figures. Consequence, deliberately: this module
makes ZERO edits to Transcriptomics_regulon.py - its diff stays empty, so
there is no regression surface to reason about, only to gate on
(`git diff --stat Transcriptomics_regulon.py` empty; test_regulon_engine.py
green).
"""

import os
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.special import ndtr

import Genomics_module as genomics_ora
import Proteomics_module as pm
import Transcriptomics_module as tx
from Transcriptomics_regulon import (
    _ulm_score, _null_slopes_shared, _permutation_p_value,
    _empirical_null_calibration, _signature_key,
)
from Proteomics_ptm_retriever import normalize_site

try:
    import requests as _requests
    HTTP_SESSION = _requests.Session()
except ImportError:
    _requests = None
    HTTP_SESSION = None

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PTM_CACHE_DIR = os.path.join(pm.CACHE_DIR, "proteomics_ptm")
PTM_REFRESH_DAYS = 90  # 44MB/26MB global files, released infrequently - matches STRING's 90, not the 30 most Proteomics bulk files use

IPTMNET_PTM_URL = "https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ptm.txt"
IPTMNET_SCORE_URL = "https://research.bioinformatics.udel.edu/iptmnet_data/files/current/score.txt"
IPTMNET_PTM_LOCAL = os.path.join(PTM_CACHE_DIR, "iptmnet_ptm.txt")
IPTMNET_SCORE_LOCAL = os.path.join(PTM_CACHE_DIR, "iptmnet_score.txt")
IPTMNET_UNIVERSE_CACHE_VERSION = 1

DIRECTIONS = pm.DIRECTIONS

KINASE_ACTIVITY_METRIC_DEFAULT = tx.GSEA_RANKING_METRIC_DEFAULT
MIN_KINASE_SITES_DEFAULT = 5
PTM_N_PERMUTATIONS_DEFAULT = 1000
PTM_N_PERMUTATIONS_CHOICES = (500, 1000, 2000, 5000)
PTM_SEED_DEFAULT = 42
PTM_MAX_CURVES_STORED = 50
PTM_MAX_NULL_SAMPLES = 500
PTM_MIN_KINASES_FOR_CALIBRATION = 30

ACTIVITY_RESULT_COLUMNS = [
    "Kinase", "Substrate Sites Tested", "Activity Score",
    "Calibrated Score (z)", "p-value (calibrated)", "FDR (calibrated)",
    "p-value (permutation)", "FDR (permutation)", "p-value (parametric)", "FDR (parametric)",
]

MIN_TERM_SIZE_DEFAULT = 3
MAX_TERM_SIZE_DEFAULT = None  # kinases have far fewer substrates than GO terms have genes; no catch-all risk

# ---------------------------------------------------------------------------
# Significant-site selection (mirrors Proteomics_module.select_significant_proteins)
# ---------------------------------------------------------------------------

def select_significant_sites(resolved_sites_df: pd.DataFrame, padj_cutoff: float, lfc_cutoff: float,
                             direction: str = "all") -> pd.DataFrame:
    """Select significant PHOSPHOSITES from the resolved site table. Same
    semantics as Proteomics_module.select_significant_proteins, keyed on
    site_id instead of accession."""
    df = resolved_sites_df[resolved_sites_df["site_id"].astype(str).str.strip() != ""].copy()
    mask = (df["padj"] <= padj_cutoff) & (df["log2FC"].abs() >= lfc_cutoff)
    if direction == "up":
        mask &= df["log2FC"] >= lfc_cutoff
    elif direction == "down":
        mask &= df["log2FC"] <= -lfc_cutoff
    return df[mask]


# ---------------------------------------------------------------------------
# iPTMnet bulk-file cache (mirrors Genomics_network._ensure_string_file's
# streaming download/cache pattern - these files are STRING-sized, not the
# small single-request downloads most Proteomics fetchers use)
# ---------------------------------------------------------------------------

def _ensure_iptmnet_file(url: str, local_path: str, refresh_days: float = PTM_REFRESH_DAYS,
                         verbose: bool = False) -> "str | None":
    if os.path.exists(local_path):
        age_days = (time.time() - os.path.getmtime(local_path)) / 86400.0
        if age_days < refresh_days:
            if verbose:
                print(f"[iPTMnet] Using cached {os.path.basename(local_path)} ({age_days:.1f} days old)")
            return local_path
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    if HTTP_SESSION is None:
        return local_path if os.path.exists(local_path) else None
    try:
        if verbose:
            print(f"[iPTMnet] Downloading {url} ...")
        resp = HTTP_SESSION.get(url, headers=genomics_ora.OMISCERA_HEADERS, timeout=180, stream=True)
        resp.raise_for_status()
        tmp_path = local_path + ".tmp"
        downloaded = 0
        with open(tmp_path, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1024 * 256):
                fh.write(chunk)
                downloaded += len(chunk)
        os.replace(tmp_path, local_path)
        if verbose:
            print(f"[iPTMnet] Downloaded {downloaded / 1e6:.1f} MB -> {local_path}")
        return local_path
    except Exception as exc:
        if verbose:
            print(f"[iPTMnet] fetch failed: {exc}")
    # A stale cached copy beats nothing, matching every other fetcher here.
    return local_path if os.path.exists(local_path) else None


def _organism_matches(organism_field: str, taxon_name: str) -> bool:
    """iPTMnet keys ptm.txt by SPECIES NAME, e.g. 'Homo sapiens (Human)', not
    NCBI taxon id - so this module needs taxon_name, not just taxon_id (the
    same reason build_reactome_universe needs it). Matched as a prefix,
    case-insensitive: 'homo sapiens (human)'.startswith('homo sapiens')."""
    name = str(taxon_name).strip().lower()
    if not name:
        return False
    return str(organism_field).strip().lower().startswith(name)


def build_kinase_substrate_universe(taxon_id: str, taxon_name: str = "", verbose: bool = False,
                                    force_refresh: bool = False) -> dict:
    """Kinase -> phosphosite universe from iPTMnet's ptm.txt, filtered to
    ptm_type=='Phosphorylation' and this organism.

    Returns {"terms": {enzyme_UniProtAC: {"name": enzyme_genename or the
    accession itself, "genes": set[canonical_site_id]}}, "universe_size",
    "fetched_at", "caveat": str|None}. The "genes" key name (not "sites") is
    DELIBERATE, not a naming slip - Genomics_module._run_fisher_for_namespace
    and its cache (de)serializers hard-key on "genes" (see
    Genomics_module.py:1565 - only that literal key round-trips through the
    set<->list JSON cache), the same reason every accession-keyed Proteomics
    universe already reuses "genes" for UniProt accessions rather than
    literal gene symbols. canonical_site_id (not site_id) is used because
    iPTMnet substrate accessions are canonical-only - see
    Proteomics_ptm_retriever's module docstring.
    """
    cache_db = f"iPTMnet_v{IPTMNET_UNIVERSE_CACHE_VERSION}"
    if not force_refresh:
        cached = genomics_ora._load_universe_cache(cache_db, taxon_id)
        if cached is not None:
            return cached

    if not str(taxon_name).strip():
        return {
            "terms": {}, "universe_size": 0,
            "fetched_at": datetime.now(timezone.utc).isoformat(),
            "caveat": (
                "iPTMnet indexes kinase-substrate relationships by organism NAME, and none was "
                "available for this run, so kinase activity/enrichment is unavailable."
            ),
        }

    path = _ensure_iptmnet_file(IPTMNET_PTM_URL, IPTMNET_PTM_LOCAL, verbose=verbose, refresh_days=(0 if force_refresh else PTM_REFRESH_DAYS))
    terms: dict = {}
    universe: set = set()
    caveat = None
    if path is None:
        caveat = "Could not fetch iPTMnet for this run (transient service error). Re-run to retry."
    else:
        matched_organism = False
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                ptm_type, _source, substrate_ac, _substrate_gene, organism, site, enzyme_ac, enzyme_gene = fields[:8]
                if ptm_type.strip().upper() != "PHOSPHORYLATION":
                    continue
                if not _organism_matches(organism, taxon_name):
                    continue
                matched_organism = True
                enzyme_ac = enzyme_ac.strip()
                if not enzyme_ac or not substrate_ac.strip():
                    continue  # no kinase attribution -> useless for activity/ORA (Stage 3's dbPTM flag will use these)
                res1, pos = normalize_site(raw_site=site)
                if res1 is None:
                    continue
                canonical_site_id = f"{substrate_ac.strip()}_{res1}{pos}"
                bucket = terms.setdefault(enzyme_ac, {"name": enzyme_gene.strip() or enzyme_ac, "genes": set()})
                bucket["genes"].add(canonical_site_id)
                universe.add(canonical_site_id)
        if not matched_organism:
            caveat = (
                f"No iPTMnet phosphorylation data matched organism \"{taxon_name}\" - iPTMnet covers "
                "human, mouse, rat, and a limited set of other model organisms, so kinase activity/"
                "enrichment is unavailable for this organism."
            )

    result = {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": datetime.now(timezone.utc).isoformat(), "caveat": caveat,
    }
    if terms:
        genomics_ora._save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[iPTMnet universe] not caching: empty (taxon={taxon_id})")
    return result


def _load_iptmnet_scores(verbose: bool = False, force_refresh: bool = False) -> "dict | None":
    """{(substrate_AC, canonical_site_id, enzyme_AC): float(score)} from
    score.txt, for the OPTIONAL minimum-confidence filter. Not organism-
    filtered (score.txt carries no organism column) - callers intersect
    against an already organism-filtered term map, so an unmatched key is
    simply never looked up. None on fetch failure (caller must treat that as
    'filter unavailable this run', never as 'every relationship scored 0')."""
    path = _ensure_iptmnet_file(IPTMNET_SCORE_URL, IPTMNET_SCORE_LOCAL, verbose=verbose,
                                refresh_days=(0 if force_refresh else PTM_REFRESH_DAYS))
    if path is None:
        return None
    scores: dict = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            substrate_ac, site, enzyme_ac, _ptm_type, score = fields[:5]
            res1, pos = normalize_site(raw_site=site)
            if res1 is None:
                continue
            try:
                score_val = float(score)
            except ValueError:
                continue
            scores[(substrate_ac.strip(), f"{substrate_ac.strip()}_{res1}{pos}", enzyme_ac.strip())] = score_val
    return scores


def _apply_score_filter(terms: dict, min_score: "float | None", scores: "dict | None") -> dict:
    """Drop (kinase, site) edges below min_score. A site with no score.txt
    entry is KEPT (score.txt is a confidence add-on over ptm.txt's own
    relationships, not a strict subset - dropping unscored edges would be a
    silent, undocumented narrowing of the universe)."""
    if min_score is None or scores is None:
        return terms
    filtered: dict = {}
    for enzyme_ac, info in terms.items():
        kept = {
            site for site in info["genes"]
            if scores.get((site.rsplit("_", 1)[0], site, enzyme_ac), min_score) >= min_score
        }
        if kept:
            filtered[enzyme_ac] = {"name": info["name"], "genes": kept}
    return filtered


# ---------------------------------------------------------------------------
# Ranked site list (mirrors Transcriptomics_module.build_ranked_list, keyed
# on site identity instead of gene symbol)
# ---------------------------------------------------------------------------

def build_ranked_site_list(resolved_sites_df: pd.DataFrame, metric: str = KINASE_ACTIVITY_METRIC_DEFAULT) -> dict:
    """Descending-sorted, deterministically-tie-broken ranked SITE list, for
    kinase activity inference. metric options: see tx.GSEA_RANKING_METRICS
    (identical set - the ranking-metric math is entirely identifier-agnostic).

    Returns {"site_ids": [str,...] (this row's OWN identity, descending
    metric order), "canonical_site_ids": [str,...] (aligned; the key used
    to match against iPTMnet's kinase sets - see Proteomics_ptm_retriever's
    module docstring for why the two differ), "values": np.ndarray,
    "metric", "n_sites", "n_input", "n_dropped_padj", "n_ties",
    "zero_crossing_rank"}.
    """
    if metric not in tx.GSEA_RANKING_METRICS:
        metric = KINASE_ACTIVITY_METRIC_DEFAULT
    df = resolved_sites_df[resolved_sites_df["site_id"].astype(str).str.strip() != ""].copy()
    n_input = len(df)

    n_dropped_padj = 0
    if metric == "signed_log_padj":
        before = len(df)
        df = df[df["padj"].notna()]
        n_dropped_padj = before - len(df)

    log2fc = df["log2FC"].to_numpy(dtype=float)
    sign = np.where(log2fc >= 0, 1.0, -1.0)
    pvalue = np.clip(df["pvalue"].to_numpy(dtype=float), tx.GSEA_MIN_P, 1.0)

    if metric == "signed_log_pvalue":
        values = sign * (-np.log10(pvalue))
    elif metric == "log2fc":
        values = log2fc
    elif metric == "pi_value":
        values = log2fc * (-np.log10(pvalue))
    else:  # signed_log_padj
        padj = np.clip(df["padj"].to_numpy(dtype=float), tx.GSEA_MIN_P, 1.0)
        values = sign * (-np.log10(padj))

    df = df.assign(_metric=values)
    df = df.sort_values(by=["_metric", "site_id"], ascending=[False, True], kind="mergesort")
    site_ids = df["site_id"].tolist()
    canonical_site_ids = df["canonical_site_id"].tolist()
    metric_values = df["_metric"].to_numpy(dtype=float)
    n_ties = int(len(metric_values) - len(np.unique(metric_values)))
    zero_crossing_rank = int(np.sum(metric_values > 0))

    return {
        "site_ids": site_ids, "canonical_site_ids": canonical_site_ids, "values": metric_values,
        "metric": metric, "n_sites": len(site_ids), "n_input": n_input,
        "n_dropped_padj": n_dropped_padj, "n_ties": n_ties,
        "zero_crossing_rank": zero_crossing_rank,
    }


# ---------------------------------------------------------------------------
# A1 — Kinase activity inference (ULM / KSEA-style)
# ---------------------------------------------------------------------------

def _empty_activity_block(note: str) -> dict:
    return {
        "table": pd.DataFrame(columns=ACTIVITY_RESULT_COLUMNS), "n_kinases_tested": 0, "note": note,
        "null_slopes_by_signature": {}, "kinase_regression_detail": {},
        "ranking_values_sample": [], "y_mean": 0.0, "n_sites": 0, "metric": "",
        "lambda_used": 1.0, "null_center": 0.0, "calibration_note": None,
    }


def run_kinase_activity(
    resolved_sites_df: pd.DataFrame,
    taxon_id: str,
    taxon_name: str = "",
    metric: str = KINASE_ACTIVITY_METRIC_DEFAULT,
    min_kinase_sites: int = MIN_KINASE_SITES_DEFAULT,
    min_score: "float | None" = None,
    n_permutations: int = PTM_N_PERMUTATIONS_DEFAULT,
    seed: int = PTM_SEED_DEFAULT,
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Kinase activity via ULM (KSEA-style): each kinase's known substrate
    sites' ranking-metric values are regressed against the rest of the
    experiment. A positive score means that kinase's substrates moved up
    (more phosphorylated) relative to everything else tested - i.e. the
    kinase appears MORE active; negative means less active.

    Unlike TF activity (Transcriptomics_regulon.run_regulon_activity), there
    is no signed/unsigned database split: every kinase->site edge from
    iPTMnet is +1 by construction (a dephosphorylating enzyme is a
    phosphatase, a DIFFERENT enzyme, not a sign on the same edge) - so
    "Activated/Repressed Targets" has no meaning here and is not reported;
    "Substrate Sites Tested" replaces it.

    THREE p-values/scores per kinase (see ACTIVITY_RESULT_COLUMNS), same
    three-column design and same caveats as run_regulon_activity - read that
    function's own docstring for the full statistical justification of each
    (parametric assumes independent, normally-distributed ranking-metric
    residuals; calibrated rescales by an empirical null estimated from the
    BULK of this run's own kinase scores; permutation is a deliberately
    UNcalibrated gene-label-permutation diagnostic). ADDITIONAL caveat
    specific to kinase substrate sets, stated in the UI and README: known
    substrate counts are wildly uneven across kinases (CK2/CDK1/PKA/mTOR
    carry far more curated sites than most other kinases), so statistical
    power is NOT uniform across the table and a top hit can reflect
    annotation depth rather than biology - read the table alongside how many
    sites each kinase had tested, not the rank alone.

    min_score: optional minimum iPTMnet score.txt confidence (0-1-ish scale,
    see that file); None (default) applies no filter. A site absent from
    score.txt is KEPT, not dropped - see _apply_score_filter's docstring.
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    _progress(0, 3, "Loading iPTMnet kinase-substrate data...")
    universe = build_kinase_substrate_universe(taxon_id, taxon_name=taxon_name, verbose=verbose, force_refresh=force_refresh)
    if universe["caveat"] and not universe["terms"]:
        return {"table": None, "n_kinases_tested": 0, "note": universe["caveat"],
               "null_slopes_by_signature": {}, "kinase_regression_detail": {},
               "ranking_values_sample": [], "y_mean": 0.0, "n_sites": 0, "metric": metric,
               "lambda_used": 1.0, "null_center": 0.0, "calibration_note": None,
               "summary": {"n_sites": 0, "metric": metric, "min_kinase_sites": min_kinase_sites,
                          "n_permutations": n_permutations, "seed": seed}}

    terms = universe["terms"]
    if min_score is not None:
        _progress(1, 3, "Applying confidence filter...")
        scores = _load_iptmnet_scores(verbose=verbose, force_refresh=force_refresh)
        terms = _apply_score_filter(terms, min_score, scores)

    ranked = build_ranked_site_list(resolved_sites_df, metric)
    n = ranked["n_sites"]
    y = np.asarray(ranked["values"], dtype=float)
    y_mean = float(y.mean()) if n else 0.0
    yc = y - y_mean
    syy = float(np.dot(yc, yc))
    rng = np.random.default_rng(seed)

    if n < 3:
        block = _empty_activity_block("Too few tested sites to fit kinase activity for this run.")
        block["summary"] = {"n_sites": n, "metric": metric, "min_kinase_sites": min_kinase_sites,
                            "n_permutations": n_permutations, "seed": seed}
        return block

    position_by_canonical: dict = {}
    for i, csid in enumerate(ranked["canonical_site_ids"]):
        position_by_canonical.setdefault(csid, []).append(i)

    _progress(2, 3, "Fitting kinase activity...")
    rows = []
    kinase_positions: dict = {}
    for enzyme_ac, info in terms.items():
        positions = sorted({p for csid in info["genes"] if csid in position_by_canonical for p in position_by_canonical[csid]})
        if len(positions) < min_kinase_sites:
            continue
        score, p_parametric, n_used, beta1 = _ulm_score(positions, [], yc, syy, n)
        if score is None:
            continue
        rows.append({
            "Kinase": info["name"], "Substrate Sites Tested": n_used,
            "Activity Score": score, "p-value (parametric)": p_parametric,
            "_beta1": beta1, "_signature": (len(positions), 0), "_enzyme_ac": enzyme_ac,
        })
        kinase_positions[enzyme_ac] = positions

    if not rows:
        block = _empty_activity_block(
            f"No kinase met the minimum substrate-site count ({min_kinase_sites} sites present in "
            "this experiment) with a usable fit."
        )
        block["summary"] = {"n_sites": n, "metric": metric, "min_kinase_sites": min_kinase_sites,
                            "n_permutations": n_permutations, "seed": seed}
        return block

    signatures_needed = {r["_signature"] for r in rows}
    null_by_signature = _null_slopes_shared(signatures_needed, yc, n, n_permutations, rng)
    for row in rows:
        null_arr = null_by_signature[row["_signature"]]
        row["p-value (permutation)"] = _permutation_p_value(abs(row["_beta1"]), null_arr)

    parametric_fdrs = tx.benjamini_hochberg([r["p-value (parametric)"] for r in rows])
    permutation_fdrs = tx.benjamini_hochberg([r["p-value (permutation)"] for r in rows])
    for row, p_fdr, perm_fdr in zip(rows, parametric_fdrs, permutation_fdrs):
        row["FDR (parametric)"] = p_fdr
        row["FDR (permutation)"] = perm_fdr

    null_center, lambda_used, n_for_calibration, calibration_note = _empirical_null_calibration(
        [r["Activity Score"] for r in rows], min_tfs=PTM_MIN_KINASES_FOR_CALIBRATION,
    )
    for row in rows:
        z = (row["Activity Score"] - null_center) / lambda_used
        row["Calibrated Score (z)"] = z
        row["p-value (calibrated)"] = max(2.0 * float(ndtr(-abs(z))), tx.GSEA_MIN_P)
    calibrated_fdrs = tx.benjamini_hochberg([r["p-value (calibrated)"] for r in rows])
    for row, cal_fdr in zip(rows, calibrated_fdrs):
        row["FDR (calibrated)"] = cal_fdr

    rows.sort(key=lambda r: r["FDR (calibrated)"])

    null_slopes_by_signature: dict = {}
    kinase_regression_detail: dict = {}
    for row in rows[:PTM_MAX_CURVES_STORED]:
        sig = row["_signature"]
        sig_key = _signature_key(*sig)
        if sig_key not in null_slopes_by_signature:
            null_slopes_by_signature[sig_key] = tx._thin_null(null_by_signature[sig], target=PTM_MAX_NULL_SAMPLES)
        positions = kinase_positions[row["_enzyme_ac"]]
        n_pos_i = len(positions)
        kinase_regression_detail[row["Kinase"]] = {
            "signature": sig_key, "beta1": row["_beta1"], "x_mean": n_pos_i / n,
            "substrate_site_ids": [ranked["site_ids"][i] for i in positions],
            "substrate_values": [round(float(y[i]), 4) for i in positions],
        }

    for row in rows:
        del row["_beta1"]
        del row["_signature"]
        del row["_enzyme_ac"]

    df = pd.DataFrame(rows, columns=ACTIVITY_RESULT_COLUMNS)
    return {
        "table": df, "n_kinases_tested": len(rows), "note": None,
        "null_slopes_by_signature": null_slopes_by_signature,
        "kinase_regression_detail": kinase_regression_detail,
        "ranking_values_sample": tx._thin_null(y, target=400),
        "y_mean": float(y_mean), "n_sites": n, "metric": metric,
        "lambda_used": lambda_used, "null_center": null_center, "calibration_note": calibration_note,
        "summary": {"n_sites": n, "metric": metric, "min_kinase_sites": min_kinase_sites,
                   "n_permutations": n_permutations, "seed": seed},
    }


# ---------------------------------------------------------------------------
# A2 — Kinase substrate-set ORA (near-free once A1 exists - same term map,
# fed through the shared Transcriptomics ORA engine, which needs only sets
# of strings and has no gene-symbol-specific logic)
# ---------------------------------------------------------------------------

def run_kinase_substrate_ora(
    resolved_sites_df: pd.DataFrame,
    taxon_id: str,
    taxon_name: str = "",
    padj_cutoff: float = pm.PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = pm.LFC_CUTOFF_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE_DEFAULT,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    min_score: "float | None" = None,
    force_refresh: bool = False,
    verbose: bool = False,
) -> dict:
    """Over-representation of significant sites among each kinase's known
    substrates, one BH-corrected family per direction (all/up/down),
    mirroring Proteomics_module.run_ora's per-direction shape but for a
    single database. Tested and universe sets are canonical_site_ids (the
    space iPTMnet itself publishes) - see build_ranked_site_list's docstring
    for why site_id and canonical_site_id are kept distinct everywhere else
    in this module; ORA is the one place collapsing to canonical is exactly
    right, since a kinase's substrate LIST cannot distinguish isoforms
    iPTMnet never annotated in the first place.
    """
    universe = build_kinase_substrate_universe(taxon_id, taxon_name=taxon_name, verbose=verbose, force_refresh=force_refresh)
    terms = universe["terms"]
    if min_score is not None and terms:
        scores = _load_iptmnet_scores(verbose=verbose, force_refresh=force_refresh)
        terms = _apply_score_filter(terms, min_score, scores)

    universe_sites = set(resolved_sites_df["canonical_site_id"].astype(str))
    tables: dict = {}
    for direction in DIRECTIONS:
        sig_df = select_significant_sites(resolved_sites_df, padj_cutoff, lfc_cutoff, direction)
        sig_sites = set(sig_df["canonical_site_id"].astype(str)) & universe_sites
        block = tx._run_ora_for_sets(terms, universe_sites, sig_sites, min_term_size, max_term_size, None)
        tables[f"Kinase substrates ({direction})"] = block

    return {"tables": tables, "universe_size": len(universe_sites), "caveat": universe["caveat"]}
