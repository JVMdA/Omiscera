# -*- coding: utf-8 -*-
"""
Transcriptomics statistical analysis module for Omiscera — Transcription
Factor / Regulon.

Standalone in the sense that it defines no new statistics engine of its own
where one already exists: a regulon (a TF's set of target genes) satisfies the
same {term_id: {"name": str, "genes": set}} contract every other Transcriptomics
gene-set database uses, so target ENRICHMENT reuses Transcriptomics_module's
existing ORA engine unchanged (see run_regulon_enrichment below). The one truly
new statistic this file adds is TF ACTIVITY inference (run_regulon_activity),
which uses the signed activation/repression edges TRRUST/DoRothEA/CollecTRI
record and ChEA3 does not.

Two analyses:
  A. TARGET ENRICHMENT (run_regulon_enrichment) — "is this TF's target set
     over-represented in my DE genes?" Signs are ignored; every database
     (including unsigned ChEA3) participates.
  B. ACTIVITY INFERENCE (run_regulon_activity) — "does the DIRECTION of this
     TF's targets' movement match what activating/repressing this TF would
     predict?" Univariate linear model (ULM, the same method decoupleR uses).
     Only signed databases (TRRUST, DoRothEA, CollecTRI) participate; ChEA3's
     GMT libraries carry no mode-of-regulation information at all.

Four regulon sources, one canonical shape per source:
    {tf_symbol: {"name": tf_symbol, "genes": set[str], "weights": {target: +1|-1}}}
  - "genes": EVERY target, regardless of whether its mode is known — this is
    what target enrichment uses (Analysis A doesn't care about sign).
  - "weights": ONLY targets with a known mode (Activation/Repression) — this
    is what activity inference uses (Analysis B). An edge whose direction the
    source itself could not determine (TRRUST "Unknown", or any OmniPath edge
    where is_stimulation/is_inhibition are both true or both false) is
    correctly a member of "genes" but absent from "weights": including it as
    +1 by default would be asserting a direction the source never claimed.

  - TRRUST — bulk TSV (grnpedia.org), human + mouse. Columns: TF, target,
    mode (Activation/Repression/Unknown), PMIDs. No header row.
  - DoRothEA / CollecTRI — both fetched from the OmniPath web service
    (https://omnipathdb.org/interactions), human + mouse, confirmed live:
    `?types=transcriptional&datasets=dorothea&dorothea_levels=A&genesymbols=1
    &organisms=9606` and `?types=transcriptional&datasets=collectri&...` both
    resolve. Direction comes from the `is_stimulation`/`is_inhibition` fields
    (present in the live response despite not being documented on the public
    docs page — verified directly, not assumed). CollecTRI's `source` field
    can carry a protein COMPLEX rather than a single TF (e.g.
    "COMPLEX:P05412_P17535" / source_genesymbol "JUN_JUND"); these rows are
    skipped, since a complex is not a single gene this app's gene-centric
    pipeline can otherwise represent.
  - ChEA3 — six primary GMT libraries downloaded from
    https://maayanlab.cloud/chea3/assets/tflibs/<library>.gmt (not ChEA3's
    live scoring API, which returns its own MeanRank/TopRank ranking rather
    than a plain gene set — mixing that ranking into this app's own
    Fisher/BH-FDR table would not be statistically coherent). Human only.
    Confirmed live format: `TF_CELLLINE_GENOME\\tgene1\\tgene2\\t...` (no GMT
    description column). The same TF appears in several rows, one per cell
    line (e.g. STAT1_HELAS3_HG19, STAT1_K562_HG19, ...) — this parser takes
    the first "_"-delimited token as the TF symbol and unions targets across
    every row for that TF into one regulon. Unsigned: "weights" is always {}.

"""

import math
import os
import re
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from scipy.special import ndtr, stdtr

import Transcriptomics_module as tx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REGULON_CACHE_DIR = os.path.join(tx.CACHE_DIR, "transcriptomics_regulons")
REGULON_REFRESH_DAYS = 30

TRRUST_URL_BY_TAXON = {
    "9606": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "10090": "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv",
}

OMNIPATH_INTERACTIONS_URL = "https://omnipathdb.org/interactions"
OMNIPATH_TAXA = ("9606", "10090")

CHEA3_GMT_BASE_URL = "https://maayanlab.cloud/chea3/assets/tflibs/"
CHEA3_LIBRARIES = (
    "ARCHS4_Coexpression", "ENCODE_ChIP-seq", "Enrichr_Queries",
    "GTEx_Coexpression", "Literature_ChIP-seq", "ReMap_ChIP-seq",
)
CHEA3_TAXON = "9606"  # human only

REGULON_DATABASES = ["TRRUST", "DoRothEA", "CollecTRI", "ChEA3"]
SIGNED_DATABASES = {"TRRUST", "DoRothEA", "CollecTRI"}  # ChEA3 carries no mode-of-regulation data

DOROTHEA_CONFIDENCE_DEFAULT = "ABC"
DOROTHEA_CONFIDENCE_CHOICES = ("A", "AB", "ABC", "ABCD", "ABCDE")

MIN_REGULON_SIZE_DEFAULT = 5
MIN_TERM_SIZE_DEFAULT = 5
MAX_TERM_SIZE_DEFAULT = 500

# Permutation null for activity inference - mirrors GSEA's own naming and
# defaults (Transcriptomics_module.py's GSEA_N_PERMUTATIONS_*/GSEA_SEED_DEFAULT/
# GSEA_MAX_CURVES_STORED) so a user who has already tuned GSEA settings finds
# the same vocabulary here.
REGULON_N_PERMUTATIONS_DEFAULT = 1000
REGULON_N_PERMUTATIONS_CHOICES = (500, 1000, 2000, 5000)
REGULON_SEED_DEFAULT = 42
REGULON_MAX_CURVES_STORED = 50   # per database - top TFs by calibrated FDR get a stored null + target list
REGULON_MAX_NULL_SAMPLES = 500   # per distinct (n_pos, n_neg) signature, for the null-fan figure

# Empirical-null calibration (genomic-control style) for activity inference.
# Measured before this was built (tests/calibration_regulon_activity.py):
# on the user's real data, the RAW parametric p-value is inflated by a factor
# lambda ~ 1.45-1.67 across TRRUST/DoRothEA/CollecTRI (22-36% of TFs clear
# p<0.05 under a pure-noise simulation with the SAME co-regulated-block
# structure, against a 5% nominal target) - a real, measured effect, not a
# theoretical worry. Calibration reduces but does NOT eliminate this
# inflation (measured: ~19% -> ~11% false-positive rate in the simulation
# that motivated this) - see ACTIVITY_RESULT_COLUMNS' calibrated columns and
# the module README for how to read it.
REGULON_MIN_TFS_FOR_CALIBRATION = 30   # below this, a median/MAD estimate has too little support to trust

# "Signed Targets Tested" (not "Regulon Size"): this counts the TF's targets
# that are BOTH present in the user's data AND carry a direction, i.e. the
# genes the activity score was actually computed from. It always equals
# Activated Targets + Repressed Targets. Deliberately NOT called "Regulon
# Size", because the enrichment table's K means the TF's whole target set in
# the universe - two different senses of "size" sitting in adjacent tables.
ACTIVITY_RESULT_COLUMNS = [
    "TF", "Signed Targets Tested", "Activated Targets", "Repressed Targets", "Activity Score",
    "Calibrated Score (z)", "p-value (calibrated)", "FDR (calibrated)",
    "p-value (permutation)", "FDR (permutation)", "p-value (parametric)", "FDR (parametric)",
]

DIRECTIONS = tx.DIRECTIONS


# ---------------------------------------------------------------------------
# On-disk cache — sets/weights round-tripped through JSON (sets are not
# natively JSON-serializable; genes are stored sorted for a deterministic file)
# ---------------------------------------------------------------------------

def _cache_path(key: str) -> str:
    os.makedirs(REGULON_CACHE_DIR, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", key)
    return os.path.join(REGULON_CACHE_DIR, safe + ".json")


def _serialize_regulons(regulons: dict) -> dict:
    return {
        tf: {"name": info["name"], "genes": sorted(info["genes"]), "weights": info["weights"]}
        for tf, info in regulons.items()
    }


def _deserialize_regulons(data: dict) -> dict:
    return {
        tf: {
            "name": info["name"],
            "genes": set(info["genes"]),
            "weights": {k: int(v) for k, v in info["weights"].items()},
        }
        for tf, info in data.items()
    }


def _empty_fetch_result(signed: bool, note: "str | None" = None, fetch_failed: bool = False) -> dict:
    return {
        "regulons": {}, "signed": signed,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "fetch_failed": fetch_failed, "note": note,
    }


def _load_cached_fetch(cache_path: str):
    cached = tx._load_json_cache(cache_path, REGULON_REFRESH_DAYS)
    if cached is None:
        return None
    cached = dict(cached)
    cached["regulons"] = _deserialize_regulons(cached["regulons"])
    return cached


def _save_fetch_result(cache_path: str, result: dict) -> None:
    to_cache = dict(result)
    to_cache["regulons"] = _serialize_regulons(result["regulons"])
    tx._save_json_cache(cache_path, to_cache)


# ---------------------------------------------------------------------------
# TRRUST — bulk TSV, human + mouse
# ---------------------------------------------------------------------------

def fetch_trrust_regulons(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """{tf: {"name","genes","weights"}} from TRRUST's bulk TSV. signed=True:
    "Unknown"-mode edges are kept in "genes" but excluded from "weights"."""
    taxon_id = str(taxon_id)
    cache_path = _cache_path(f"trrust_{taxon_id}")
    if not force_refresh:
        cached = _load_cached_fetch(cache_path)
        if cached is not None:
            return cached

    url = TRRUST_URL_BY_TAXON.get(taxon_id)
    if url is None:
        return _empty_fetch_result(True, note=f"TRRUST only covers human and mouse; no data for taxon {taxon_id}.")
    if tx.HTTP_SESSION is None:
        return _empty_fetch_result(True, note="requests is not installed.", fetch_failed=True)

    try:
        resp = tx.HTTP_SESSION.get(url, headers=tx.OMISCERA_HEADERS, timeout=tx.DEFAULT_TIMEOUT * 2)
        if resp.status_code != 200:
            return _empty_fetch_result(True, note=f"TRRUST fetch failed (HTTP {resp.status_code}).", fetch_failed=True)
        text = resp.content.decode("utf-8", errors="replace")
    except Exception as exc:
        if verbose:
            print(f"[TRRUST] fetch failed: {exc}")
        return _empty_fetch_result(True, note=f"TRRUST fetch failed: {exc}", fetch_failed=True)

    regulons: dict = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        tf = parts[0].strip().upper()
        target = parts[1].strip().upper()
        mode = parts[2].strip()
        if not tf or not target:
            continue
        bucket = regulons.setdefault(tf, {"name": tf, "genes": set(), "weights": {}})
        bucket["genes"].add(target)
        if mode == "Activation":
            bucket["weights"][target] = 1
        elif mode == "Repression":
            bucket["weights"][target] = -1
        # "Unknown" -> target stays in "genes" only, no weight assigned.

    result = _empty_fetch_result(True)
    result["regulons"] = regulons
    if regulons:
        _save_fetch_result(cache_path, result)
    elif verbose:
        print("[TRRUST] not caching: no regulons parsed")
    return result


# ---------------------------------------------------------------------------
# DoRothEA / CollecTRI — both via the OmniPath web service
# ---------------------------------------------------------------------------

def _build_regulons_from_omnipath_rows(rows: "list[dict]") -> dict:
    regulons: dict = {}
    for row in rows:
        source = str(row.get("source") or "")
        if source.startswith("COMPLEX:"):
            continue  # a TF complex, not a single gene - see module docstring
        tf = str(row.get("source_genesymbol") or "").strip().upper()
        target = str(row.get("target_genesymbol") or "").strip().upper()
        if not tf or not target:
            continue
        bucket = regulons.setdefault(tf, {"name": tf, "genes": set(), "weights": {}})
        bucket["genes"].add(target)
        is_stim = bool(row.get("is_stimulation"))
        is_inhib = bool(row.get("is_inhibition"))
        if is_stim and not is_inhib:
            bucket["weights"][target] = 1
        elif is_inhib and not is_stim:
            bucket["weights"][target] = -1
        # both true or both false -> direction ambiguous, no weight assigned.
    return regulons


def _fetch_omnipath_dataset(
    dataset: str, taxon_id: str, dorothea_levels: "str | None", verbose: bool = False,
) -> dict:
    if tx.HTTP_SESSION is None:
        return _empty_fetch_result(True, note="requests is not installed.", fetch_failed=True)
    params = {
        "types": "transcriptional",
        "organisms": taxon_id,
        "genesymbols": "1",
        "datasets": dataset,
        "fields": "sources",
        "format": "json",
    }
    if dataset == "dorothea" and dorothea_levels:
        params["dorothea_levels"] = ",".join(list(dorothea_levels))
    try:
        resp = tx.HTTP_SESSION.get(
            OMNIPATH_INTERACTIONS_URL, params=params, headers=tx.OMISCERA_HEADERS,
            timeout=tx.DEFAULT_TIMEOUT * 4,
        )
        if resp.status_code != 200:
            return _empty_fetch_result(True, note=f"OmniPath fetch failed (HTTP {resp.status_code}).", fetch_failed=True)
        rows = resp.json()
    except Exception as exc:
        if verbose:
            print(f"[OmniPath:{dataset}] fetch failed: {exc}")
        return _empty_fetch_result(True, note=f"OmniPath fetch failed: {exc}", fetch_failed=True)

    if not isinstance(rows, list):
        return _empty_fetch_result(True, note="OmniPath returned an unexpected response shape.", fetch_failed=True)

    result = _empty_fetch_result(True)
    result["regulons"] = _build_regulons_from_omnipath_rows(rows)
    return result


def fetch_dorothea_regulons(
    taxon_id: str, confidence_levels: str = DOROTHEA_CONFIDENCE_DEFAULT,
    verbose: bool = False, force_refresh: bool = False,
) -> dict:
    taxon_id = str(taxon_id)
    if taxon_id not in OMNIPATH_TAXA:
        return _empty_fetch_result(True, note=f"DoRothEA only covers human and mouse; no data for taxon {taxon_id}.")
    cache_key = f"dorothea_{taxon_id}_{confidence_levels}"
    cache_path = _cache_path(cache_key)
    if not force_refresh:
        cached = _load_cached_fetch(cache_path)
        if cached is not None:
            return cached
    result = _fetch_omnipath_dataset("dorothea", taxon_id, confidence_levels, verbose=verbose)
    if result["regulons"]:
        _save_fetch_result(cache_path, result)
    return result


def fetch_collectri_regulons(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    taxon_id = str(taxon_id)
    if taxon_id not in OMNIPATH_TAXA:
        return _empty_fetch_result(True, note=f"CollecTRI only covers human and mouse; no data for taxon {taxon_id}.")
    cache_path = _cache_path(f"collectri_{taxon_id}")
    if not force_refresh:
        cached = _load_cached_fetch(cache_path)
        if cached is not None:
            return cached
    result = _fetch_omnipath_dataset("collectri", taxon_id, None, verbose=verbose)
    if result["regulons"]:
        _save_fetch_result(cache_path, result)
    return result


# ---------------------------------------------------------------------------
# ChEA3 — six primary GMT libraries, human only, unsigned
# ---------------------------------------------------------------------------

def fetch_chea3_regulons(
    libraries: "list[str] | None" = None, verbose: bool = False, force_refresh: bool = False,
) -> dict:
    libraries = list(libraries) if libraries else list(CHEA3_LIBRARIES)
    libraries = [lib for lib in libraries if lib in CHEA3_LIBRARIES] or list(CHEA3_LIBRARIES)
    cache_path = _cache_path("chea3_" + "_".join(sorted(libraries)))
    if not force_refresh:
        cached = _load_cached_fetch(cache_path)
        if cached is not None:
            return cached

    if tx.HTTP_SESSION is None:
        return _empty_fetch_result(False, note="requests is not installed.", fetch_failed=True)

    regulons: dict = {}
    any_ok = False
    errors = []
    for lib in libraries:
        try:
            resp = tx.HTTP_SESSION.get(
                CHEA3_GMT_BASE_URL + lib + ".gmt", headers=tx.OMISCERA_HEADERS, timeout=tx.DEFAULT_TIMEOUT * 4,
            )
            if resp.status_code != 200:
                errors.append(f"{lib}: HTTP {resp.status_code}")
                continue
            text = resp.content.decode("utf-8", errors="replace")
        except Exception as exc:
            errors.append(f"{lib}: {exc}")
            if verbose:
                print(f"[ChEA3] {lib} fetch failed: {exc}")
            continue
        for line in text.splitlines():
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            # "TF_CELLLINE_GENOME" - the same TF appears once per cell line;
            # targets are unioned across every such row (see module docstring).
            tf = parts[0].split("_")[0].strip().upper()
            genes = {g.strip().upper() for g in parts[1:] if g.strip()}
            if not tf or not genes:
                continue
            bucket = regulons.setdefault(tf, {"name": tf, "genes": set(), "weights": {}})
            bucket["genes"].update(genes)
        any_ok = True

    if not any_ok:
        note = "ChEA3 fetch failed for every library (" + "; ".join(errors) + ")." if errors else "ChEA3 fetch failed."
        return _empty_fetch_result(False, note=note, fetch_failed=True)

    result = _empty_fetch_result(False)
    result["regulons"] = regulons
    if errors and verbose:
        print(f"[ChEA3] partial fetch, {len(errors)} librar(y/ies) failed: {errors}")
    if regulons:
        _save_fetch_result(cache_path, result)
    return result


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def fetch_regulons(
    source: str, taxon_id: str,
    dorothea_levels: str = DOROTHEA_CONFIDENCE_DEFAULT,
    chea3_libraries: "list[str] | None" = None,
    verbose: bool = False, force_refresh: bool = False,
) -> dict:
    if source == "TRRUST":
        return fetch_trrust_regulons(taxon_id, verbose=verbose, force_refresh=force_refresh)
    if source == "DoRothEA":
        return fetch_dorothea_regulons(taxon_id, dorothea_levels, verbose=verbose, force_refresh=force_refresh)
    if source == "CollecTRI":
        return fetch_collectri_regulons(taxon_id, verbose=verbose, force_refresh=force_refresh)
    if source == "ChEA3":
        # ChEA3's libraries use human gene symbols only (like MSigDB Hallmarks
        # in Transcriptomics_module.run_ora) - for any other organism the
        # symbols would silently fail to match rather than erroring, so this
        # is gated explicitly instead of returning a near-empty, misleading
        # enrichment/activity result.
        if str(taxon_id).strip() != CHEA3_TAXON:
            return _empty_fetch_result(
                False,
                note=f"ChEA3 uses human gene symbols and is only available for human (taxon {CHEA3_TAXON}); "
                     f"this run's organism is taxon {taxon_id}.",
            )
        return fetch_chea3_regulons(chea3_libraries, verbose=verbose, force_refresh=force_refresh)
    raise ValueError(f"Unknown regulon source: {source}")


# ---------------------------------------------------------------------------
# Analysis A — TF target enrichment (no new statistics: reuses
# Transcriptomics_module._run_ora_for_sets exactly as GO/Reactome/Hallmark do)
# ---------------------------------------------------------------------------

def _empty_ora_block(note: str, n_universe_total: int) -> dict:
    return {
        "table": pd.DataFrame(columns=tx.RESULT_COLUMNS),
        "n_terms_tested": 0, "n_terms_excluded_min_size": 0, "n_terms_no_overlap": 0,
        "universe_size": 0, "n_deg_in_universe": 0,
        "n_universe_total": n_universe_total, "n_deg_total": 0,
        "note": note,
    }


def run_regulon_enrichment(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    padj_cutoff: float = tx.PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = tx.LFC_CUTOFF_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE_DEFAULT,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    dorothea_levels: str = DOROTHEA_CONFIDENCE_DEFAULT,
    chea3_libraries: "list[str] | None" = None,
    background_symbols: "set | None" = None,
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """TF target enrichment: each selected regulon database's TFs are treated
    as plain gene sets (signs ignored) and tested exactly like GO/Reactome/
    Hallmark via Transcriptomics_module._run_ora_for_sets - same universe
    convention (genes tested in this experiment), same per-database BH-FDR,
    same term-size gate. No length-bias correction: regulon membership has no
    documented relationship to transcript length the way RNA-seq detection
    does, so none is applied here (unlike run_ora)."""
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universe_df = resolved_df[resolved_df["gene_symbol"].astype(str).str.strip() != ""]
    universe_symbols_all = set(universe_df["gene_symbol"].astype(str).str.upper())
    n_universe_before_background = len(universe_symbols_all)
    if background_symbols:
        universe_symbols_all = universe_symbols_all & {str(s).strip().upper() for s in background_symbols}

    directions = {}
    for direction_key in DIRECTIONS:
        deg_df = tx.select_degs(resolved_df, padj_cutoff, lfc_cutoff, direction_key)
        raw_deg_set = set(deg_df["gene_symbol"].astype(str).str.upper())
        directions[direction_key] = raw_deg_set & universe_symbols_all

    tables: dict = {}
    steps_total = max(len(databases) * len(DIRECTIONS), 1)
    step = 0

    for db in databases:
        fetched = fetch_regulons(
            db, taxon_id, dorothea_levels=dorothea_levels, chea3_libraries=chea3_libraries,
            verbose=verbose, force_refresh=force_refresh,
        )
        term_map = {
            tf: {"name": info["name"], "genes": info["genes"]}
            for tf, info in fetched["regulons"].items()
        }
        for direction_key in DIRECTIONS:
            step += 1
            _progress(step, steps_total, f"{db} ({direction_key})")
            result_key = f"{db} ({direction_key})"
            deg_symbols = directions[direction_key]
            if fetched.get("fetch_failed"):
                tables[result_key] = _empty_ora_block(
                    f"Could not fetch {db} for this run (transient service error). Re-run to retry.",
                    len(universe_symbols_all),
                )
                continue
            if not term_map:
                tables[result_key] = _empty_ora_block(
                    fetched.get("note") or f"No {db} regulons are available for this organism.",
                    len(universe_symbols_all),
                )
                continue
            if not deg_symbols:
                tables[result_key] = _empty_ora_block(
                    "No DEGs in this direction at the current thresholds.", len(universe_symbols_all),
                )
                continue
            tables[result_key] = tx._run_ora_for_sets(
                term_map, universe_symbols_all, deg_symbols, min_term_size, max_term_size, None,
            )

    summary = {
        "n_universe": len(universe_symbols_all),
        "n_universe_before_background": n_universe_before_background,
        "background_applied": bool(background_symbols),
        "n_deg_all": len(directions["all"]),
        "n_deg_up": len(directions["up"]),
        "n_deg_down": len(directions["down"]),
    }
    return {"tables": tables, "summary": summary}


# ---------------------------------------------------------------------------
# Analysis B — TF activity inference (univariate linear model, decoupleR's
# "ulm" method). The one genuinely new statistic in this module.
# ---------------------------------------------------------------------------

def _ulm_score(
    pos_positions: "list[int]", neg_positions: "list[int]",
    yc: np.ndarray, syy: float, n: int,
) -> "tuple[float, float, int, float] | tuple[None, None, int, None]":
    """Closed-form ULM t-statistic and two-sided PARAMETRIC p-value for one
    TF's regulon, derived analytically from ordinary least squares rather
    than looping over every gene per TF (verified equivalent to a naive
    per-gene OLS fit - see test_regulon_engine.py):

    x_i in {+1, -1, 0} (activated target / repressed target / not a target).
    Because x is 0 outside the regulon:
        sum(x)      = n_pos - n_neg
        sum(x^2)    = n_pos + n_neg = n_total
        sum(x*y_c)  = sum(y_c over activated targets) - sum(y_c over repressed targets)
    which gives SXX and SXY without ever touching the n-n_total non-target
    genes individually. SSE then follows from the OLS identity
    SSE = SYY - beta1^2 * SXX (since SXY = beta1 * SXX), so the residual sum
    of squares - which in principle depends on every gene - is recovered
    without a second per-gene pass either.

    This parametric p-value assumes the ranking-metric RESIDUALS are
    approximately normal with constant variance, an assumption real ranking-
    metric distributions (e.g. signed -log10(p), often a spike near zero plus
    a long tail) can violate. run_regulon_activity ALSO reports a permutation
    p-value (see _null_slopes_shared below) built empirically from the actual
    observed distribution, with no normality assumption.

    IMPORTANT, checked directly rather than assumed: permutation is NOT a
    strict fix for co-regulated (correlated) targets specifically. Gene-label
    permutation draws its null from RANDOM same-sized subsets of the ranked
    list, which are typically uncorrelated - so it under-represents how
    extreme a genuinely co-regulated regulon's statistic can be, in the same
    direction as the parametric test's own weakness (this mirrors the
    documented limitation of gene-label permutation for GSEA itself - see
    run_gsea's module docstring). The two p-values are reported side by side
    as two different estimates with two different, known limitations, not as
    one superseding the other - see the module README.

    Returns (score, p_value, n_targets_used, beta1); (None, None,
    n_targets_used, None) if the fit is degenerate (fewer than 2 usable
    targets, or zero variance in x). beta1 is returned because the caller
    needs the raw slope (not just the standardized score) to look it up
    against a permutation null of slopes - see _null_slopes_shared, which
    relies on beta1 being a strictly monotonic function of |score| at a
    fixed (n_pos, n_neg) (since SSE = SYY - beta1^2*SXX depends only on
    beta1 for a fixed signature) to make storing null SLOPES alone sufficient
    - verified in test_regulon_engine.py.
    """
    n_pos, n_neg = len(pos_positions), len(neg_positions)
    n_total = n_pos + n_neg
    if n_total < 2:
        return None, None, n_total, None
    sum_x = n_pos - n_neg
    x_mean = sum_x / n
    sxx = n_total - n * (x_mean ** 2)
    if sxx <= 1e-12:
        return None, None, n_total, None
    sum_xy_c = (float(yc[pos_positions].sum()) if pos_positions else 0.0) \
        - (float(yc[neg_positions].sum()) if neg_positions else 0.0)
    beta1 = sum_xy_c / sxx
    sse = syy - (beta1 ** 2) * sxx
    df_resid = n - 2
    if df_resid <= 0 or sse <= 0:
        return None, None, n_total, None
    se = math.sqrt(sse / df_resid / sxx)
    if se <= 0:
        return None, None, n_total, None
    score = beta1 / se
    # Numerically stable two-sided p: stdtr(df, -x) = 1 - stdtr(df, x) exactly
    # by symmetry, avoiding the cancellation error of computing 1 - stdtr(df, x)
    # directly for large |score| (same rationale as this codebase's existing
    # ndtr-based p-values elsewhere).
    p_value = 2.0 * float(stdtr(df_resid, -abs(score)))
    return score, p_value, n_total, beta1


def _null_slopes_shared(
    signatures, yc: np.ndarray, n: int, n_permutations: int,
    rng: "np.random.Generator", chunk: int = 200,
) -> dict:
    """Permutation null of the ULM slope (beta1), one array per (n_pos, n_neg)
    signature, computed together and sharing one random draw per permutation
    trial across every signature - the same technique
    Transcriptomics_module._null_es_shared uses for GSEA (see that function's
    docstring for the full justification).

    Under gene-label permutation, the null distribution of beta1 for a
    signature depends ONLY on (n_pos, n_neg) - not on which genes are chosen
    - because sum_xy_c is just "sum of a random n_pos-subset of yc minus sum
    of a disjoint random n_neg-subset", and sxx is fixed by (n_pos, n_neg, n)
    alone. So one shared random draw per trial serves every signature: a
    smaller signature just takes a further contiguous prefix of a larger
    one's random subset (identical idea to GSEA's per-size sub-prefixing).

    Only slopes are stored (not the t-statistic) - see _ulm_score's docstring
    for why that is sufficient: an empirical p-value on |beta1| is identical
    to one on |t| at a fixed signature.

    Returns {(n_pos, n_neg): (n_permutations,) float array of null beta1}.
    """
    sig_list = sorted({(int(p_), int(g_)) for p_, g_ in signatures})
    if not sig_list:
        return {}
    max_total = max(p_ + g_ for p_, g_ in sig_list)
    results = {sig: np.empty(n_permutations, dtype=float) for sig in sig_list}
    done = 0
    while done < n_permutations:
        p = min(chunk, n_permutations - done)
        rand_matrix = rng.random((p, n))
        # Same off-by-one rationale as _null_es_shared: argpartition's kth is
        # a 0-indexed position INTO the array, so the boundary that must be
        # correct for the smallest max_total elements is max_total - 1.
        idx_max = np.argpartition(rand_matrix, max_total - 1, axis=1)[:, :max_total]
        # argpartition's within-partition order is unspecified, not random -
        # an explicit reshuffle is required before taking a smaller prefix
        # (same as _null_es_shared).
        shuffle_key = rng.random((p, max_total))
        shuffle_order = np.argsort(shuffle_key, axis=1)
        idx_shuffled = np.take_along_axis(idx_max, shuffle_order, axis=1)
        yc_shuffled = yc[idx_shuffled]              # (p, max_total)
        cs = np.cumsum(yc_shuffled, axis=1)         # (p, max_total)
        zeros_col = np.zeros(p)
        for (n_pos, n_neg) in sig_list:
            total = n_pos + n_neg
            sum_pos = cs[:, n_pos - 1] if n_pos > 0 else zeros_col
            sum_total = cs[:, total - 1]
            sum_neg = sum_total - sum_pos
            sum_xy_c = sum_pos - sum_neg
            x_mean = (n_pos - n_neg) / n
            sxx = total - n * (x_mean ** 2)
            beta1_null = (sum_xy_c / sxx) if sxx > 1e-12 else zeros_col
            results[(n_pos, n_neg)][done:done + p] = beta1_null
        done += p
    return results


def _permutation_p_value(observed_abs_beta1: float, null_beta1: np.ndarray) -> float:
    """Empirical two-sided p-value with the standard +1/+1 floor (matches
    Transcriptomics_module._normalize_es's convention), so a reported p can
    never be exactly 0 regardless of how extreme the observed slope is."""
    null_abs = np.abs(null_beta1)
    count_ge = int(np.sum(null_abs >= observed_abs_beta1))
    return (count_ge + 1) / (null_abs.size + 1)


def _phi(x: float) -> float:
    """Standard normal density. Three-term expression, no new dependency."""
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _truncated_normal_sigma(deviations: np.ndarray, sigma0: float, c_mult: float = 1.5, iters: int = 60) -> float:
    """Truncated-normal MLE of a centred normal's scale, using only the
    CENTRAL window |x| < c = c_mult * sigma (fixed-point iteration on sigma).

    WHY: plain MAD is fooled by real signal, not just by true inflation - see
    _empirical_null_calibration's docstring for the measured numbers. Signal
    contamination in the TAILS pushes the median absolute deviation up the
    null's own CDF. Restricting the estimate to a central window where the
    null dominates, then correcting analytically for having thrown the tails
    away, is dramatically less fooled by that contamination while still
    recovering GENUINE inflation (both measured in
    tests/probe_better_null_estimator.py and
    tests/probe_truncated_validity.py before this was written - the
    second probe is the one that matters: an estimator that fixed the false-
    alarm case by failing to recover true inflation would break the feature
    calibration exists for, and this one does not).

    For a centred normal truncated to |x| < c:
        E[x^2 | |x|<c] = sigma^2 * (1 - 2a*phi(a) / (2*Phi(a) - 1)),  a = c/sigma
    Solved for sigma by fixed-point iteration, seeded from the MAD estimate.
    Falls back to sigma0 (the seed) if the window ever holds too few points
    or the iteration produces a non-finite/non-positive value - a small
    database must degrade to the plain-MAD estimate, not to garbage.
    """
    sigma = sigma0
    for _ in range(iters):
        c = c_mult * sigma
        inner = deviations[np.abs(deviations) < c]
        if inner.size < 10:
            return sigma0
        m2 = float(np.mean(inner ** 2))
        a = c / sigma
        denom = 2.0 * float(ndtr(a)) - 1.0
        if denom <= 1e-12:
            return sigma0
        factor = 1.0 - (2.0 * a * _phi(a)) / denom
        if factor <= 1e-12:
            return sigma0
        new_sigma = math.sqrt(m2 / factor)
        if not math.isfinite(new_sigma) or new_sigma <= 0:
            return sigma0
        if abs(new_sigma - sigma) < 1e-9:
            return new_sigma
        sigma = new_sigma
    return sigma


def _empirical_null_calibration(scores: "list[float] | np.ndarray", min_tfs: int = REGULON_MIN_TFS_FOR_CALIBRATION):
    """Genomic-control-style empirical-null calibration for the parametric
    ULM t-statistics of every TF tested in one database run.

    WHY THIS IS VALID FOR THIS STATISTIC SPECIFICALLY: under the global null
    (no TF has a real effect), _ulm_score's t-statistic is Student's-t(n-2)
    for EVERY TF regardless of regulon size or activated/repressed split -
    the degrees of freedom depend only on n (the number of tested genes),
    which is the same for every TF in a run. That makes every TF's score
    exchangeable under one shared null, which is what makes pooling them
    into a single empirical null legitimate - this is not true of every
    statistic, it is checked to be true of this one (see
    test_regulon_engine.py's exchangeability test).

    MEASURED, not assumed (tests/calibration_regulon_activity.py): most
    TFs in a real run are not truly active, so the BULK of the observed
    score distribution approximates the true null - and on real data that
    bulk is measurably wider than Student's-t(n-2) predicts, consistent with
    target co-regulation inflating the naive variance estimate. Calibration
    divides that inflation back out.

    Center: median of the observed scores. Scale: a TRUNCATED-NORMAL MLE
    (_truncated_normal_sigma), seeded from and falling back to plain
    MAD/0.6745 - see that function's docstring for why plain MAD alone was
    replaced and what "materially better, not perfect" means in practice.

    MEASURED (tests/probe_better_null_estimator.py and
    probe_truncated_validity.py), at the TF counts this module actually
    produces (~130-500):
                              MAD lambda   truncated lambda   (truth)
        0% truly active        ~1.02          ~1.02            1.0
        10%                    ~1.11-1.13     ~1.03-1.06       1.0
        25%                    ~1.41-1.44     ~1.04-1.08       1.0
        40%                    ~2.03-2.04     ~1.04-1.10       1.0
        genuine inflation, sigma=1.6            ~1.56          1.6  (still recovered)
        genuine inflation, sigma=3.0            ~2.92          3.0  (still recovered)
        40% active, GRADED (continuous, not clumped) effect sizes:
                                ~1.43          ~1.16            1.0

    CONSEQUENCE - lambda is STILL AMBIGUOUS, just less so
    (audit_lambda_ambiguity.py). A large lambda can still mean (a) genuine
    null inflation from co-regulated targets, or (b) many TFs genuinely
    active - the truncated estimator narrows the gap between those two
    readings (see the graded-signal row above: ~1.16 instead of ~1.43) but
    does not close it to 1.0. Both estimators degrade past ~55% truly-active
    TFs, beyond any robust estimator and beyond any plausible experiment.
    Calibration remains the CONSERVATIVE choice, not the provably correct
    one. The parametric columns are kept beside it precisely so a user with
    independent reason to expect widespread real TF activity can read those
    instead. The module README states this, with the numbers above.

    lambda is CLAMPED to >= 1.0 (t(n-2)'s own sd is ~1.0003 at realistic n,
    so 1.0 is the correct "no inflation detected" floor): calibration can
    only ever make a result MORE conservative, never less. If the empirical
    spread came out narrower than theory, dividing by it would manufacture
    significance out of a noise estimate rather than correct for real
    inflation - this is the genomic-control convention, not a bespoke choice.

    Returns (center, lambda_used, n_used, note). note is None when
    calibration was actually applied; otherwise it explains why lambda was
    left at 1.0 (too few TFs, or a degenerate/zero-MAD score distribution).
    """
    arr = np.asarray(scores, dtype=float)
    arr = arr[np.isfinite(arr)]
    n_used = int(arr.size)
    if n_used < min_tfs:
        return 0.0, 1.0, n_used, (
            f"Only {n_used} TF(s) were tested in this database - fewer than the {min_tfs} needed for a "
            "reliable empirical-null estimate, so calibration was skipped (lambda = 1.0, uncalibrated)."
        )
    center = float(np.median(arr))
    deviations = arr - center
    mad = float(np.median(np.abs(deviations)))
    if mad <= 0:
        return center, 1.0, n_used, (
            "The observed TF scores had zero spread (MAD = 0) in this run, which cannot support an "
            "empirical-null estimate - calibration was skipped (lambda = 1.0, uncalibrated)."
        )
    mad_sigma = mad / 0.6745
    sigma_hat = _truncated_normal_sigma(deviations, mad_sigma)
    lambda_used = max(1.0, sigma_hat)
    return center, lambda_used, n_used, None


def _empty_activity_block(note: str) -> dict:
    return {
        "table": pd.DataFrame(columns=ACTIVITY_RESULT_COLUMNS), "n_tfs_tested": 0, "note": note,
        "null_slopes_by_signature": {}, "tf_regression_detail": {},
        "ranking_values_sample": [], "y_mean": 0.0, "n_genes": 0, "metric": "",
        "lambda_used": 1.0, "null_center": 0.0, "calibration_note": None,
    }


def _signature_key(n_pos: int, n_neg: int) -> str:
    """JSON object keys are always strings - this is the single format used
    both when writing null_slopes_by_signature and when a TF's
    tf_regression_detail entry points into it, so the two always agree."""
    return f"{n_pos}_{n_neg}"


def run_regulon_activity(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    metric: str = tx.GSEA_RANKING_METRIC_DEFAULT,
    min_regulon_size: int = MIN_REGULON_SIZE_DEFAULT,
    dorothea_levels: str = DOROTHEA_CONFIDENCE_DEFAULT,
    n_permutations: int = REGULON_N_PERMUTATIONS_DEFAULT,
    seed: int = REGULON_SEED_DEFAULT,
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """TF activity inference via ULM. For each TF in a SIGNED database, the
    gene-level ranking metric (build_ranked_list - the same statistic and the
    same user-facing metric choice as GSEA) is regressed on the TF's signed
    target-weight vector across every tested gene. A positive score means the
    TF's activated targets moved up and repressed targets moved down relative
    to the rest of the experiment - i.e. the TF appears ACTIVE; negative means
    it appears repressed. ChEA3 (unsigned) is not run here - see
    _empty_activity_block's note for that database.

    THREE p-values/scores are reported per TF (see ACTIVITY_RESULT_COLUMNS):
      - "p-value (parametric)" from _ulm_score's t-test - decoupleR's own
        "ulm" method. Assumes every gene's ranking value is an independent
        draw, which co-regulated regulon targets can violate.
      - "Calibrated Score (z)" / "p-value (calibrated)" / "FDR (calibrated)"
        from _empirical_null_calibration: the parametric scores are
        rescaled by an empirical-null lambda estimated from the BULK of this
        run's own TF scores (most TFs are not truly active, so they
        approximate the true null). THIS is the primary column - what
        sorts/highlights the table and drives chart filtering. Measured
        (tests/calibration_regulon_activity.py) to meaningfully reduce,
        but NOT eliminate, the excess false-positive rate the parametric
        p-value carries under target co-regulation - report it as an
        improvement, not a complete fix.
      - "p-value (permutation)" from a gene-label permutation null
        (_null_slopes_shared), kept as a separate, deliberately UNcalibrated
        diagnostic - it makes no normality assumption, but is not immune to
        correlated targets either (checked directly, not assumed: an early
        version of this docstring claimed permutation alone "fixes" the
        parametric p-value's co-regulation weakness; a test built to
        demonstrate that showed it does not reliably hold). A TF where the
        calibrated and permutation columns disagree sharply is worth a
        closer look in the Regulon Detail figures, not a reason to trust one
        column over the other.

    For the REGULON_MAX_CURVES_STORED top TFs per database (by calibrated
    FDR), the block also stores a thinned null (for the regression-plot fan)
    and each such TF's activated/repressed target SYMBOLS (for the point
    cloud) - not for every TF, mirroring GSEA's GSEA_MAX_CURVES_STORED
    (storing a full curve for every one of potentially hundreds of TFs would
    be far more payload than any results screen needs).
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    ranked = tx.build_ranked_list(resolved_df, metric)
    symbols = ranked["symbols"]
    n = len(symbols)
    position_by_symbol = {s: i for i, s in enumerate(symbols)}
    y = np.asarray(ranked["values"], dtype=float)
    y_mean = float(y.mean()) if n else 0.0
    yc = y - y_mean
    syy = float(np.dot(yc, yc))
    rng = np.random.default_rng(seed)

    # Keys are pre-inserted in the CALLER's order (not signed-then-unsigned)
    # so the client's Object.keys(tables) - which drives which database its
    # dropdown shows first - reflects the order the user selected databases
    # in, rather than always defaulting to an unsigned (always-empty) one
    # first regardless of what was picked. Caught via a Playwright test that
    # asserted the activity chart toggles actually appear: ChEA3 (unsigned,
    # always 0 TFs) was silently winning the "first key" slot every time.
    tables: dict = {db: None for db in databases}
    signed_dbs = [db for db in databases if db in SIGNED_DATABASES]
    unsigned_dbs = [db for db in databases if db not in SIGNED_DATABASES]
    steps_total = max(len(signed_dbs), 1)
    step = 0

    for db in unsigned_dbs:
        tables[db] = _empty_activity_block(
            f"{db} carries no mode-of-regulation (activation/repression) information, so TF "
            "activity cannot be inferred from it - see Target Enrichment instead."
        )

    for db in signed_dbs:
        step += 1
        _progress(step, steps_total, db)
        fetched = fetch_regulons(db, taxon_id, dorothea_levels=dorothea_levels, verbose=verbose, force_refresh=force_refresh)
        if fetched.get("fetch_failed"):
            tables[db] = _empty_activity_block(f"Could not fetch {db} for this run (transient service error). Re-run to retry.")
            continue
        if not fetched["regulons"]:
            tables[db] = _empty_activity_block(fetched.get("note") or f"No {db} regulons are available for this organism.")
            continue
        if n < 3:
            tables[db] = _empty_activity_block("Too few tested genes to fit TF activity for this run.")
            continue

        # Pass 1: parametric fit for every candidate TF (cheap, closed-form).
        # tf_positions keeps the position lists around for pass 2 - the
        # permutation null needs to know every row's (n_pos, n_neg)
        # signature before it can be computed, and every row's beta1 to look
        # up its own permutation p-value afterwards.
        rows = []
        tf_positions: dict = {}
        for tf, info in fetched["regulons"].items():
            weights = info.get("weights") or {}
            if not weights:
                continue
            pos_positions = [position_by_symbol[t] for t, w in weights.items() if w > 0 and t in position_by_symbol]
            neg_positions = [position_by_symbol[t] for t, w in weights.items() if w < 0 and t in position_by_symbol]
            if len(pos_positions) + len(neg_positions) < min_regulon_size:
                continue
            score, p_parametric, n_used, beta1 = _ulm_score(pos_positions, neg_positions, yc, syy, n)
            if score is None:
                continue
            rows.append({
                "TF": tf, "Signed Targets Tested": n_used,
                "Activated Targets": len(pos_positions), "Repressed Targets": len(neg_positions),
                "Activity Score": score, "p-value (parametric)": p_parametric,
                "_beta1": beta1, "_signature": (len(pos_positions), len(neg_positions)),
            })
            tf_positions[tf] = (pos_positions, neg_positions)

        if not rows:
            tables[db] = _empty_activity_block(
                f"No {db} TF met the minimum regulon size ({min_regulon_size} targets present in "
                "this experiment) with a usable signed fit."
            )
            continue

        # Pass 2: permutation null, shared across every distinct signature
        # among the surviving rows (see _null_slopes_shared's docstring).
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

        # Empirical-null calibration - see _empirical_null_calibration's
        # docstring for why this statistic specifically supports pooling
        # every TF's score into one null, and calibration_regulon_activity.py
        # for the measurement that justified building this at all.
        null_center, lambda_used, n_for_calibration, calibration_note = _empirical_null_calibration(
            [r["Activity Score"] for r in rows]
        )
        for row in rows:
            z = (row["Activity Score"] - null_center) / lambda_used
            row["Calibrated Score (z)"] = z
            # Floored at GSEA_MIN_P (Transcriptomics_module.py:1596, the same
            # constant GSEA's own p-values are clamped with) so an extreme z
            # cannot underflow to a literal 0 - which would put a point at
            # y=300 on the -log10(p) volcano and flatten every other point
            # (measured reachable at |z|~=38; not on the bundled example
            # dataset, but a stronger real dataset could get there).
            row["p-value (calibrated)"] = max(2.0 * float(ndtr(-abs(z))), tx.GSEA_MIN_P)
        calibrated_fdrs = tx.benjamini_hochberg([r["p-value (calibrated)"] for r in rows])
        for row, cal_fdr in zip(rows, calibrated_fdrs):
            row["FDR (calibrated)"] = cal_fdr

        # Sorted/highlighted by the CALIBRATED FDR - the primary column (see
        # this function's docstring). Parametric and permutation remain for
        # comparison; the permutation p-value is deliberately left
        # uncalibrated (see the module README).
        rows.sort(key=lambda r: r["FDR (calibrated)"])

        # Stored detail (thinned null + target symbols) only for the top TFs
        # by calibrated FDR - see this function's docstring.
        null_slopes_by_signature: dict = {}
        tf_regression_detail: dict = {}
        for row in rows[:REGULON_MAX_CURVES_STORED]:
            sig = row["_signature"]
            sig_key = _signature_key(*sig)
            if sig_key not in null_slopes_by_signature:
                null_slopes_by_signature[sig_key] = tx._thin_null(null_by_signature[sig], target=REGULON_MAX_NULL_SAMPLES)
            pos_positions, neg_positions = tf_positions[row["TF"]]
            n_pos_i, n_neg_i = len(pos_positions), len(neg_positions)
            # Ship the RANKING-METRIC value per target - NOT the client's
            # log2FC lookup. These are the exact y-values the fit consumed;
            # an earlier version of the regression figure plotted log2FC
            # instead, which is a different variable entirely whenever the
            # metric isn't "log2fc" (the default is signed -log10(p)), so
            # the drawn slope did not match the drawn points at all.
            # x_mean is likewise the fit's own (n_pos-n_neg)/n over ALL
            # tested genes, not over the targets - the OLS centroid the
            # fitted and null lines must pivot through.
            tf_regression_detail[row["TF"]] = {
                "signature": sig_key,
                "beta1": row["_beta1"],
                "x_mean": (n_pos_i - n_neg_i) / n,
                "activated_symbols": [symbols[i] for i in pos_positions],
                "activated_values": [round(float(y[i]), 4) for i in pos_positions],
                "repressed_symbols": [symbols[i] for i in neg_positions],
                "repressed_values": [round(float(y[i]), 4) for i in neg_positions],
            }

        for row in rows:
            del row["_beta1"]
            del row["_signature"]

        df = pd.DataFrame(rows, columns=ACTIVITY_RESULT_COLUMNS)
        tables[db] = {
            "table": df, "n_tfs_tested": len(rows), "note": None,
            "null_slopes_by_signature": null_slopes_by_signature,
            "tf_regression_detail": tf_regression_detail,
            # Quantile-thinned sample of EVERY tested gene's ranking value, so
            # the regression figure can draw the x=0 (non-target) mass the fit
            # is actually anchored by. Reuses _thin_null purely as a
            # distribution-preserving thinner (it is just evenly-spaced order
            # statistics - see its docstring). Targets are <1% of genes, so
            # this all-genes distribution is visually indistinguishable from
            # any single TF's true non-target set; the figure labels it "all
            # tested genes" rather than claiming otherwise.
            "ranking_values_sample": tx._thin_null(y, target=400),
            "y_mean": float(y_mean),
            "n_genes": n,
            "metric": metric,
            "lambda_used": lambda_used,
            "null_center": null_center,
            "calibration_note": calibration_note,
        }

    summary = {
        "n_genes": n, "metric": metric, "min_regulon_size": min_regulon_size,
        "n_permutations": n_permutations, "seed": seed,
    }
    return {"tables": tables, "summary": summary}

