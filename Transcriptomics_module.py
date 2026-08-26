# -*- coding: utf-8 -*-
"""
Transcriptomics statistical analysis module for Omiscera — Over-Representation
Analysis (ORA).

Standalone: does not import any Genomics_* module. Tests whether the genes in
a differential-expression (DE) result (produced upstream by
Transcriptomics_module_retriever.py) are over-represented in GO terms, Reactome
pathways, or MSigDB gene sets (8 curated collections, human + mouse - see
MSIGDB_COLLECTIONS).

Guiding principles
-------------------
  - Universe = the genes actually tested in this experiment (every mapped row
    of the resolved DE table), not a generic genome-wide background. A gene
    that was never expressed/detected in this experiment could not have been
    called DE regardless of its biology, so it does not belong in the
    background either.
  - Identifier spine = gene SYMBOL. Each database is fetched in whatever ID
    space it is natively keyed by (GO -> symbol, via QuickGO's own SYMBOL
    column; Reactome -> Entrez; MSigDB -> symbol) and reported back in symbol
    space, since the resolver already gives every gene its symbol, Entrez ID,
    and Ensembl ID together. No cross-database ID translation file is needed.
  - Statistical test: Wallenius non-central hypergeometric (GOseq-style gene-
    length bias correction — longer genes accumulate more reads and are more
    likely to be called DE independent of biology), falling back to the plain
    hypergeometric test when per-gene length data is unavailable for the
    organism or scikit-learn is not installed.
  - Multiple-testing families are never pooled across databases/collections -
    each direction of each database (and each MSigDB collection) gets its
    own independent Benjamini-Hochberg correction, same as Genomics ORA and
    the Regulon module. Selecting many MSigDB collections in one run can
    produce up to 24 independent families; both run_ora's and run_gsea's
    summary report n_independent_families so the UI can say so plainly.

Module layout
-------------
  1. Constants
  2. Stats engine              (pure functions, no I/O)
  3. On-disk JSON cache helpers (shared by every fetcher)
  4. Gene-set fetchers          (GO, Reactome, MSigDB)
  5. Gene-length fetcher + PWF fitting (length-bias correction)
  6. DEG selection + universe
  7. ORA orchestrator (run_ora)
"""

import csv
import io
import json
import math
import os
import re
import time
from datetime import datetime, timezone
from typing import Sequence

import numpy as np
import pandas as pd

import go_hierarchy
from Genomics_module import format_term_name

try:
    import requests as _requests
    HTTP_SESSION = _requests.Session()
except ImportError:
    _requests = None
    HTTP_SESSION = None

try:
    from scipy.stats import hypergeom as _scipy_hypergeom
    from scipy.stats import nchypergeom_wallenius as _scipy_wallenius
except ImportError:
    _scipy_hypergeom = None
    _scipy_wallenius = None

try:
    from sklearn.isotonic import IsotonicRegression as _IsotonicRegression
except ImportError:
    _IsotonicRegression = None

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OMISCERA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositorio do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "*/*",
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
ORA_UNIVERSE_CACHE_DIR = os.path.join(CACHE_DIR, "transcriptomics_ora_universe")

DEFAULT_TIMEOUT = 30
DEFAULT_RETRIES = 3

# GO — QuickGO/GOA downloadSearch (native GO annotation source). The TSV's own
# SYMBOL column is used directly as gene identity — no UniProt bridge needed.
QUICKGO_DOWNLOAD_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
QUICKGO_DOWNLOAD_LIMIT = 2_000_000
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GO_OBO_PATH = os.path.join(CACHE_DIR, "go-basic.obo")  # species-independent; shared filename is fine
GO_OBO_REFRESH_DAYS = 30
GO_REFRESH_DAYS = 30
GO_NAMESPACES = ("biological_process", "molecular_function", "cellular_component")

# Reactome — NCBI2Reactome_All_Levels.txt bulk file (native source, Entrez-keyed).
REACTOME_URL = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"
REACTOME_PATH = os.path.join(CACHE_DIR, "transcriptomics_ncbi2reactome.txt")
REACTOME_REFRESH_DAYS = 30

# MSigDB — native download, gene-symbol-keyed, human + mouse. VERIFIED LIVE
# before this was built: https://data.broadinstitute.org/gsea-msigdb/msigdb/
# release/ is openly listable, NOT registration-gated as an earlier version
# of this comment assumed (that led to fetching Hallmark only, via Enrichr's
# redistribution - stale at 2020 data and human-only; replaced by this).
# Main collection is CC BY 4.0 (c) Broad Institute/MIT/UC Regents; three
# subsets carry incompatible license terms and are deliberately never
# fetched here (KEGG_LEGACY/KEGG_MEDICUS/BioCarta - see
# info/README-ORA-transcriptomics.txt "MSigDB collections and licensing").
#
# Version is PINNED, not auto-discovered: a silent upstream bump would change
# every user's results with no record of why. Bump MSIGDB_VERSION deliberately;
# tests/live_msigdb_urls.py fails loudly if the pinned URLs stop resolving.
MSIGDB_VERSION = "2026.1"
MSIGDB_BASE_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
MSIGDB_REFRESH_DAYS = 30

# Mouse collection numbering does NOT mirror human (verified against
# MSigDB's own human/mouse collections pages) - e.g. mouse M8 IS cell type
# like human C8, but mouse M3 is regulatory targets, not a C3 analog of
# anything offered here, and PID/Oncogenic have no mouse file at all. So
# availability is this EXPLICIT per-collection map, never a c->m string
# substitution. mouse_file=None means genuinely absent upstream, not merely
# unfetched. Excluded from this list entirely (see the README section above
# for why): C5:GO and C2:CP:REACTOME (duplicate Omiscera's own GO/Reactome
# with a different source/version), C3:TFT (overlaps the Regulon module),
# C1/C4 (niche), and the three license carve-outs.
MSIGDB_COLLECTIONS = {
    "hallmark":      {"label": "Hallmark",              "human_file": "h.all",              "mouse_file": "mh.all"},
    "perturbations": {"label": "Perturbations (CGP)",   "human_file": "c2.cgp",              "mouse_file": "m2.cgp"},
    "wikipathways":  {"label": "WikiPathways",          "human_file": "c2.cp.wikipathways",  "mouse_file": "m2.cp.wikipathways"},
    "pid":           {"label": "PID",                   "human_file": "c2.cp.pid",           "mouse_file": None},
    "phenotype":     {"label": "Phenotype (HPO / MPT)", "human_file": "c5.hpo",              "mouse_file": "m5.mpt"},
    "oncogenic":     {"label": "Oncogenic",              "human_file": "c6.all",              "mouse_file": None},
    "immunologic":   {"label": "Immunologic",            "human_file": "c7.all",              "mouse_file": "m7.all"},
    "cell_type":     {"label": "Cell type",              "human_file": "c8.all",              "mouse_file": "m8.all"},
}
# taxon -> (URL organism suffix, which MSIGDB_COLLECTIONS field holds the
# filename for that organism).
MSIGDB_TAXON_SUFFIX = {"9606": ("Hs", "human_file"), "10090": ("Mm", "mouse_file")}

# Ensembl BioMart — bulk per-gene transcript length, used only for GOseq-style
# length-bias correction. The longest transcript per gene is used as that
# gene's representative length: a standard, practical proxy for the length-vs-
# detectability trend (goseq's own default, total exonic length from a full
# transcript database, needs exon-level structure this app does not otherwise
# fetch — the longest-transcript length is a documented simplification, not a
# silent one).
BIOMART_URL = "https://www.ensembl.org/biomart/martservice"
BIOMART_DATASET_BY_TAXON = {
    "9606": "hsapiens_gene_ensembl",
    "10090": "mmusculus_gene_ensembl",
    "10116": "rnorvegicus_gene_ensembl",
    "7955": "drerio_gene_ensembl",
    "7227": "dmelanogaster_gene_ensembl",
    "6239": "celegans_gene_ensembl",
    "559292": "scerevisiae_gene_ensembl",
}
GENE_LENGTH_CACHE_NAME_TMPL = "transcriptomics_gene_lengths_{taxon}.csv"
GENE_LENGTH_REFRESH_DAYS = 90  # genome builds/gene models change rarely
BIOMART_CHUNK_SIZE = 300

# clusterProfiler's defaults (the standard transcriptomics ORA reference):
# terms with too few genes are statistically unstable (one shared gene can
# swing k/K wildly); terms with hundreds of genes tend to be uninformatively
# generic ("protein binding"). Both remain user-editable per run.
MIN_TERM_SIZE_DEFAULT = 10
MAX_TERM_SIZE_DEFAULT = 500
PADJ_CUTOFF_DEFAULT = 0.05
LFC_CUTOFF_DEFAULT = 1.0

RESULT_COLUMNS = [
    "Term ID", "Term Name", "k", "K", "n", "N",
    "Fold Enrichment", "p-value", "FDR", "Bonferroni", "Overlapping Genes",
]

DIRECTIONS = ("all", "up", "down")


# ---------------------------------------------------------------------------
# Stats engine (pure functions, no I/O — unit-testable without network)
# ---------------------------------------------------------------------------

def hypergeometric_pvalue(k: int, n: int, K: int, N: int) -> float:
    """One-sided (over-representation) hypergeometric p-value: P(X >= k).

    N = universe size, K = genes in this set (within the universe),
    n = number of DEGs (within the universe), k = DEGs in this set.
    Equivalent to a one-sided ('greater') Fisher's exact test on the same
    2x2 table.
    """
    if _scipy_hypergeom is None:
        raise RuntimeError("scipy is required for hypergeometric_pvalue")
    if k <= 0:
        return 1.0
    return float(_scipy_hypergeom.sf(k - 1, N, K, n))


def wallenius_pvalue(k: int, n: int, K: int, N: int, odds: float) -> float:
    """Length-bias-corrected over-representation p-value: P(X >= k) under the
    Wallenius non-central hypergeometric distribution with the given odds
    (bias) ratio. odds=1.0 is mathematically equivalent to the plain
    hypergeometric test (verified in the engine's unit tests).

    scipy's nchypergeom_wallenius(M, n, N, odds) uses M=population size,
    n=Type-I objects in the population, N=number of draws — which collide
    with this module's own N/K/n/k naming, so the mapping here is:
    M<-N (universe), scipy's n<-K (set size), scipy's N<-n (number of DEGs).
    """
    if _scipy_wallenius is None:
        raise RuntimeError("scipy is required for wallenius_pvalue")
    if k <= 0:
        return 1.0
    # Both ends clamped: scipy silently saturates sf() to 0.0/1.0 at extreme
    # odds rather than raising, so an unclamped odds ratio would fail silently
    # too - _category_odds's averaged weights realistically never approach
    # these bounds, but the clamp costs nothing and removes the edge case.
    #
    # The two bounds are NOT symmetric in risk. Verified empirically: raising
    # odds monotonically drives sf() UP toward 1.0 (conservative - a high
    # odds ratio means the category's genes were "expected" to be drawn, so a
    # given overlap becomes less surprising). scipy saturates to exactly 1.0
    # by odds~100 in practice and only starts returning nan above ~1e15. So a
    # clamp that binds on the high side fails safe (p -> 1, never a false
    # positive). The 1e3 clamp previously used here bound well before that
    # safe saturation point, which would make a genuinely large odds ratio
    # (e.g. from a length-skewed category with few PWF-covered genes) read as
    # LESS conservative than it should - i.e. clamping down made the p-value
    # smaller than the correct, larger one. Raised to 1e9: still 1e6 x more
    # headroom than any realistic _category_odds output, and comfortably
    # below scipy's nan threshold, so it is a safety net rather than a bound
    # that actually shapes results.
    odds = max(1e-6, min(float(odds), 1e9))
    return float(_scipy_wallenius.sf(k - 1, N, K, n, odds))


def fold_enrichment(k: int, n: int, K: int, N: int) -> float:
    """Observed/expected ratio: (k/n) / (K/N)."""
    if n <= 0 or N <= 0 or K <= 0:
        return float("nan")
    expected = n * (K / N)
    if expected <= 0:
        return float("nan")
    return k / expected


def benjamini_hochberg(pvalues: list) -> list:
    """Benjamini-Hochberg FDR correction. Returns adjusted q-values in the
    same order as the input list."""
    m = len(pvalues)
    if m == 0:
        return []
    indexed = sorted(range(m), key=lambda i: pvalues[i])
    ranked_p = [pvalues[i] for i in indexed]
    q_at_rank = [0.0] * m
    for rank, p in enumerate(ranked_p, start=1):
        q_at_rank[rank - 1] = p * m / rank
    running_min = 1.0
    for i in range(m - 1, -1, -1):
        running_min = min(running_min, q_at_rank[i])
        q_at_rank[i] = min(running_min, 1.0)
    result = [0.0] * m
    for rank_pos, original_idx in enumerate(indexed):
        result[original_idx] = q_at_rank[rank_pos]
    return result


def bonferroni(pvalues: list) -> list:
    """Bonferroni correction: min(p * m, 1.0) for each p-value. Stricter and
    more conservative than BH-FDR; reported alongside it for users who prefer
    family-wise error control (matches Genomics_module.bonferroni)."""
    m = len(pvalues)
    return [min(p * m, 1.0) for p in pvalues]


PWF_SIGNIFICANCE_ALPHA = 0.05


def _rankdata(values: np.ndarray) -> np.ndarray:
    """Average ranks (1-indexed, ties get the mean rank of their block) -
    the same convention scipy.stats.rankdata's default 'average' method uses.
    Pure numpy; needed because this file avoids scipy.stats.spearmanr (see
    _spearman_corr_and_pvalue)."""
    n = len(values)
    order = np.argsort(values, kind="mergesort")
    sorted_vals = values[order]
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    return ranks


def _spearman_corr_and_pvalue(x: "Sequence[float]", y: "Sequence[float]") -> "tuple[float, float]":
    """Spearman rank correlation (rho) and its two-sided p-value under
    H0: rho=0, computed with plain numpy + math.erfc - NOT scipy.stats.spearmanr,
    which is verified (in this project's runtime) to trigger a native
    access-violation crash.

    rho itself is just the Pearson correlation of the two rank arrays.
    Significance uses the standard Fisher z-transform large-sample
    approximation (z = arctanh(rho) * sqrt(n-3), asymptotically standard
    normal under H0), which is accurate for the sample sizes real gene lists
    have (hundreds to tens of thousands of genes) and needs nothing beyond
    the standard library's math.erfc for the normal tail.
    """
    n = len(x)
    if n < 4:
        return 0.0, 1.0
    rx = _rankdata(np.asarray(x, dtype=float))
    ry = _rankdata(np.asarray(y, dtype=float))
    rx_c = rx - rx.mean()
    ry_c = ry - ry.mean()
    denom = math.sqrt(float(np.sum(rx_c ** 2)) * float(np.sum(ry_c ** 2)))
    if denom <= 0:
        return 0.0, 1.0  # one array is constant (e.g. every gene DE, or none) - no correlation is definable
    rho = float(np.sum(rx_c * ry_c)) / denom
    rho = max(-0.999999, min(0.999999, rho))  # guard arctanh(+/-1) = +/-inf
    z = math.atanh(rho) * math.sqrt(n - 3)
    p_value = math.erfc(abs(z) / math.sqrt(2.0))  # two-sided: erfc already covers both tails at |z|
    return rho, p_value


def fit_pwf(is_deg: "list[bool]", lengths: "list[float]") -> "list[float] | None":
    """Fit a Probability Weighting Function: P(DE) as a monotone function of
    gene length, via isotonic regression (goseq-style). Returns one weight per
    gene (aligned to the input order), or None if no statistically
    significant length/DE relationship is detected for this gene list -
    callers should then skip length correction and use the plain
    (uncorrected) test instead, rather than fitting a spurious direction to
    noise.

    The monotone DIRECTION is not fit-and-compared (an earlier version tried
    both directions and kept whichever had lower residual error - on truly
    unbiased data that comparison is decided by sampling noise, and was
    observed to pick the wrong direction ~half the time). Instead, the sign
    of the gene-length/DE-status Spearman correlation decides the direction
    directly, and its significance decides whether to apply a correction at
    all: standard full-length RNA-seq protocols give longer genes more reads
    and thus more power to be called DE (positive correlation -> increasing
    PWF), but 3'-tagged protocols (QuantSeq, most droplet single-cell) give
    ~1 read per transcript regardless of length, so the correlation can be
    flat or negative for those - both are handled correctly by this test
    rather than assumed away.

    Note: this deliberately does NOT use scikit-learn's own
    `increasing="auto"` option (which internally calls scipy.stats.spearmanr,
    the function verified to crash here) or scipy.stats.spearmanr directly -
    see _spearman_corr_and_pvalue for the crash-free replacement.

    Raises RuntimeError if scikit-learn is unavailable — callers must catch
    this and fall back to the plain hypergeometric test rather than silently
    skip length correction.
    """
    if _IsotonicRegression is None:
        raise RuntimeError("scikit-learn is required for PWF fitting")
    x = np.asarray(lengths, dtype=float)
    y = np.asarray(is_deg, dtype=float)
    rho, p_value = _spearman_corr_and_pvalue(lengths, is_deg)
    if p_value >= PWF_SIGNIFICANCE_ALPHA or rho == 0.0:
        return None
    iso = _IsotonicRegression(y_min=1e-6, y_max=1.0 - 1e-6, out_of_bounds="clip", increasing=(rho > 0))
    weights = iso.fit_transform(x, y)
    return weights.tolist()


def _category_odds(weights: dict, set_genes: set, universe_genes: set) -> float:
    """Convert PWF weights into a per-category Wallenius odds ratio: the mean
    weight of genes IN the category divided by the mean weight of genes NOT in
    it (goseq-style bias-to-odds conversion). Returns 1.0 (no bias) when there
    isn't enough length data to compute a meaningful ratio.
    """
    in_w = [weights[g] for g in set_genes if g in weights]
    out_genes = universe_genes - set_genes
    out_w = [weights[g] for g in out_genes if g in weights]
    if not in_w or not out_w:
        return 1.0
    mean_in = sum(in_w) / len(in_w)
    mean_out = sum(out_w) / len(out_w)
    if mean_out <= 0:
        return 1.0
    return mean_in / mean_out


# ---------------------------------------------------------------------------
# On-disk JSON cache helpers (shared by all three gene-set fetchers)
# ---------------------------------------------------------------------------

def _universe_cache_path(cache_key: str) -> str:
    os.makedirs(ORA_UNIVERSE_CACHE_DIR, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", cache_key)
    return os.path.join(ORA_UNIVERSE_CACHE_DIR, safe + ".json")


def _load_json_cache(path: str, refresh_days: float):
    if not os.path.exists(path):
        return None
    age_days = (time.time() - os.path.getmtime(path)) / 86400
    if age_days >= refresh_days:
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return None


def _save_json_cache(path: str, data: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh)


# ---------------------------------------------------------------------------
# GO gene sets (QuickGO, symbol-keyed)
# ---------------------------------------------------------------------------

def _quickgo_safe_get(params: dict, verbose: bool = False):
    if HTTP_SESSION is None:
        return None
    headers = dict(OMISCERA_HEADERS)
    headers["Accept"] = "text/tsv"
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            response = HTTP_SESSION.get(
                QUICKGO_DOWNLOAD_URL, params=params, headers=headers, timeout=DEFAULT_TIMEOUT * 4,
            )
            if response.status_code in (200, 201):
                return response.content
            if response.status_code in (429, 500, 502, 503, 504):
                last_error = f"HTTP {response.status_code}"
                if attempt < DEFAULT_RETRIES:
                    time.sleep(1.0)
                    continue
            last_error = f"HTTP {response.status_code}"
            break
        except Exception as exc:
            last_error = str(exc)
            if attempt < DEFAULT_RETRIES:
                time.sleep(0.5)
                continue
            break
    if verbose:
        print(f"[QuickGO] request failed: {last_error}")
    return None


def _parse_quickgo_tsv_symbol(content: bytes, include_iea: bool) -> dict:
    """Parse a QuickGO downloadSearch TSV for one aspect, keyed by gene SYMBOL
    (not UniProt accession). Returns {go_id: {"name": str, "genes": set[SYMBOL]}}.
    """
    term_map: dict = {}
    lines = content.decode("utf-8-sig", errors="replace").splitlines()
    if not lines:
        return term_map
    header = [h.strip().upper() for h in lines[0].split("\t")]
    col = {h: i for i, h in enumerate(header)}

    def _find(*candidates):
        for c in candidates:
            if c in col:
                return col[c]
        return None

    db_col = _find("GENE PRODUCT DB", "DB")
    sym_col = _find("SYMBOL", "GENE PRODUCT SYMBOL")
    go_id_col = _find("GO TERM", "GO ID", "GO TERM ID")
    go_name_col = _find("GO NAME", "GO TERM NAME")
    ev_col = _find("GO EVIDENCE CODE", "EVIDENCE CODE", "GO EVIDENCE")

    if sym_col is None or go_id_col is None:
        return term_map

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if db_col is not None and db_col < len(parts):
            if parts[db_col].strip() not in ("UniProtKB", ""):
                continue
        if sym_col >= len(parts):
            continue
        symbol = parts[sym_col].strip().upper()
        if not symbol:
            continue
        if not include_iea and ev_col is not None and ev_col < len(parts):
            if parts[ev_col].strip() == "IEA":
                continue
        if go_id_col >= len(parts):
            continue
        go_id = parts[go_id_col].strip()
        if not go_id:
            continue
        go_name = parts[go_name_col].strip() if go_name_col is not None and go_name_col < len(parts) else go_id
        bucket = term_map.setdefault(go_id, {"name": go_name, "genes": set()})
        bucket["genes"].add(symbol)

    return term_map


def _ensure_go_obo_file(verbose: bool = False) -> str:
    if os.path.exists(GO_OBO_PATH):
        age_days = (time.time() - os.path.getmtime(GO_OBO_PATH)) / 86400
        if age_days < GO_OBO_REFRESH_DAYS:
            return GO_OBO_PATH
    os.makedirs(CACHE_DIR, exist_ok=True)
    assert HTTP_SESSION is not None, "requests not installed"
    if verbose:
        print(f"[GO OBO] Downloading {GO_OBO_URL} ...")
    resp = HTTP_SESSION.get(GO_OBO_URL, headers=OMISCERA_HEADERS, timeout=180, stream=True)
    resp.raise_for_status()
    # Stream to a temp file and only replace the real path once fully
    # downloaded - an interrupted write straight to GO_OBO_PATH would leave a
    # truncated file with a FRESH mtime, reused as "fresh" for GO_OBO_REFRESH_DAYS.
    tmp_path = GO_OBO_PATH + ".tmp"
    with open(tmp_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
    os.replace(tmp_path, GO_OBO_PATH)
    return GO_OBO_PATH


def fetch_go_sets(taxon_id: str, include_iea: bool = True, verbose: bool = False, force_refresh: bool = False) -> dict:
    """Fetch GO gene sets for a taxon, keyed by gene SYMBOL, one dict per
    namespace.

    Each namespace's direct QuickGO annotations are propagated up the GO
    hierarchy (is_a/part_of ancestors, via go_hierarchy.propagate_term_map)
    before being returned, so a gene annotated only to a specific child term
    also counts toward its parent terms - matching clusterProfiler/goseq/
    topGO/DAVID. Without this, every non-leaf GO term would be systematically
    under-populated, since QuickGO's downloadSearch API returns direct
    annotations only.

    Returns:
        {
          "biological_process": {go_id: {"name": str, "genes": [symbol, ...]}},
          "molecular_function": {...},
          "cellular_component": {...},
          "fetched_at": iso str,
        }
    """
    cache_suffix = "iea" if include_iea else "noiea"
    # "_v2" marks the cache format that includes ancestor propagation - bumped
    # so an on-disk cache written before this feature existed is never reused.
    cache_path = _universe_cache_path(f"go_{taxon_id}_{cache_suffix}_v2")
    if not force_refresh:
        cached = _load_json_cache(cache_path, GO_REFRESH_DAYS)
        # Defense in depth, mirroring the write-side guard below: only ever
        # written when propagated=True and no namespace failed, but reject
        # a blob missing that stamp rather than trust it blindly.
        if cached is not None and cached.get("propagated") is True:
            return cached

    # A network hiccup on the 32 MB OBO must degrade to direct-annotation-only
    # results, not abort the whole ORA run (which would also discard Reactome/
    # Hallmark/Custom-database results for no reason connected to them).
    hierarchy = None
    try:
        obo_path = _ensure_go_obo_file(verbose=verbose)
        hierarchy = go_hierarchy.parse_go_obo(obo_path, verbose=verbose)
    except Exception as exc:
        if verbose:
            print(f"[GO] WARNING: could not load/parse GO hierarchy ({exc}); using direct annotations only")
    go_names = hierarchy["names"] if hierarchy is not None else {}

    result: dict = {ns: {} for ns in GO_NAMESPACES}
    failed_namespaces = []
    for ns in GO_NAMESPACES:
        params = {
            "taxonId": taxon_id,
            "aspect": ns,
            "geneProductType": "protein",
            "downloadLimit": QUICKGO_DOWNLOAD_LIMIT,
        }
        if verbose:
            print(f"[GO] Fetching {ns} (taxon={taxon_id}) from QuickGO...")
        raw = _quickgo_safe_get(params, verbose=verbose)
        if raw is None:
            failed_namespaces.append(ns)
            continue
        term_map = _parse_quickgo_tsv_symbol(raw, include_iea)
        for go_id, info in term_map.items():
            info["name"] = go_names.get(go_id, info["name"])
        if hierarchy is not None:
            term_map = go_hierarchy.propagate_term_map(term_map, hierarchy)
        for info in term_map.values():
            info["genes"] = sorted(info["genes"])
        result[ns] = term_map
        if verbose:
            print(f"[GO] {ns}: {len(term_map)} terms" + (" (after ancestor propagation)" if hierarchy is not None else " (direct annotations only)"))

    result["fetched_at"] = datetime.now(timezone.utc).isoformat()
    result["fetch_failed_namespaces"] = failed_namespaces
    result["propagated"] = hierarchy is not None
    if failed_namespaces or hierarchy is None:
        # A transient QuickGO failure, or a failure to load/parse the GO
        # hierarchy, must not poison the cache with a partial or unpropagated
        # universe for GO_REFRESH_DAYS - serve this run's (possibly degraded)
        # result but let the NEXT run retry from scratch.
        if verbose:
            print(f"[GO] not caching: failed_namespaces={failed_namespaces}, propagated={hierarchy is not None}")
    else:
        _save_json_cache(cache_path, result)
    return result


# ---------------------------------------------------------------------------
# Reactome gene sets (native bulk file, Entrez-keyed)
# ---------------------------------------------------------------------------

def _ensure_reactome_file(verbose: bool = False) -> str:
    if os.path.exists(REACTOME_PATH):
        age_days = (time.time() - os.path.getmtime(REACTOME_PATH)) / 86400
        if age_days < REACTOME_REFRESH_DAYS:
            if verbose:
                print(f"[Reactome] Using cached file ({age_days:.1f} days old)")
            return REACTOME_PATH
    os.makedirs(CACHE_DIR, exist_ok=True)
    assert HTTP_SESSION is not None, "requests not installed"
    if verbose:
        print(f"[Reactome] Downloading {REACTOME_URL} ...")
    resp = HTTP_SESSION.get(REACTOME_URL, headers=OMISCERA_HEADERS, timeout=120, stream=True)
    resp.raise_for_status()
    # Stream to a temp file, then replace atomically - an interrupted 98 MB
    # download straight to REACTOME_PATH would leave a truncated file with a
    # fresh mtime, reused as "fresh" for REACTOME_REFRESH_DAYS.
    tmp_path = REACTOME_PATH + ".tmp"
    with open(tmp_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
    os.replace(tmp_path, REACTOME_PATH)
    return REACTOME_PATH


def _genus_species(name: str) -> str:
    """First two whitespace-separated tokens, lowercased - e.g. 'Saccharomyces
    cerevisiae S288C' -> 'saccharomyces cerevisiae'. Used as a fallback match
    when a strain-level organism name (as returned by UniProt) doesn't equal
    Reactome's own (genus+species-only) organism string exactly."""
    return " ".join(name.strip().lower().split()[:2])


def _parse_reactome_file(tsv_path: str, target_org: str, target_genus_species: str, include_iea: bool) -> tuple:
    """Parse a Reactome NCBI2Reactome_All_Levels-style TSV (columns:
    entrez_id, pathway_id, url, pathway_name, evidence, organism - tab-
    separated, no header), matching rows to target_org exactly first, with a
    genus+species fallback used only if the exact match yields nothing.
    Extracted from fetch_reactome_sets (identical logic, just given a name)
    so it's independently unit-testable without a real 98 MB download.

    Returns (pathways, used_fallback) where pathways is
    {pathway_id: {"name": str, "genes": [entrez_id, ...]}} (genes SORTED,
    matching fetch_reactome_sets' return shape).
    """
    pathways: dict = {}
    pathways_fallback: dict = {}
    with open(tsv_path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            entrez_id, pathway_id, _url, pathway_name, evidence, organism = parts[:6]
            entrez_id = entrez_id.strip()
            pathway_id = pathway_id.strip()
            organism_lower = organism.strip().lower()
            if not entrez_id or not pathway_id:
                continue
            if not include_iea and evidence.strip().upper() == "IEA":
                continue
            if organism_lower == target_org:
                bucket = pathways.setdefault(pathway_id, {"name": pathway_name.strip(), "genes": set()})
                bucket["genes"].add(entrez_id)
            elif target_genus_species and _genus_species(organism_lower) == target_genus_species:
                bucket = pathways_fallback.setdefault(pathway_id, {"name": pathway_name.strip(), "genes": set()})
                bucket["genes"].add(entrez_id)

    used_fallback = False
    if not pathways and pathways_fallback:
        pathways = pathways_fallback
        used_fallback = True

    for info in pathways.values():
        info["genes"] = sorted(info["genes"])
    return pathways, used_fallback


def fetch_reactome_sets(
    taxon_name: str, include_iea: bool = True, verbose: bool = False, force_refresh: bool = False,
) -> dict:
    """Fetch Reactome pathway gene sets for an organism (matched by scientific
    name), keyed by NCBI Entrez Gene ID.

    Matching is exact-string first; if that yields zero pathways (e.g. a
    strain-level name like "Saccharomyces cerevisiae S288C" that Reactome
    itself only records as "Saccharomyces cerevisiae"), falls back to
    matching on genus+species alone. The fallback is reported via
    "matched_via_genus_species_fallback" so callers can attach a caveat.

    include_iea: the bulk file's evidence column carries "TAS" (curated,
    Traceable Author Statement) and "IEA" (Electronic, e.g. orthology-
    projected annotations) exactly like GO's evidence codes - confirmed by
    inspecting the actual downloaded file. When False, IEA rows are excluded,
    mirroring Genomics' build_reactome_universe. Cache key encodes include_iea
    so IEA/no-IEA universes are stored separately, same convention as GO.

    Returns {"pathways": {pathway_id: {"name": str, "genes": [entrez_id, ...]}},
    "fetched_at": iso str, "matched_via_genus_species_fallback": bool}.

    Raises ValueError for a blank taxon_name: an earlier version silently
    pooled every organism's pathways together in that case (matching nothing
    was treated as matching everything) - a real bug once live, fixed by
    having run_ora refuse to call this without a resolved organism name. That
    guard lives in the caller, so this is defense in depth for any other
    caller, not the normal way this case is reached.
    """
    if not taxon_name.strip():
        raise ValueError("fetch_reactome_sets requires a non-blank taxon_name.")

    cache_suffix = "iea" if include_iea else "noiea"
    cache_key = f"reactome_{taxon_name.strip().lower().replace(' ', '_')}_{cache_suffix}"
    cache_path = _universe_cache_path(cache_key)
    if not force_refresh:
        cached = _load_json_cache(cache_path, REACTOME_REFRESH_DAYS)
        if cached is not None:
            return cached

    tsv_path = _ensure_reactome_file(verbose=verbose)
    target_org = taxon_name.strip().lower()
    target_genus_species = _genus_species(taxon_name)
    pathways, used_fallback = _parse_reactome_file(tsv_path, target_org, target_genus_species, include_iea)

    result = {
        "pathways": pathways,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "matched_via_genus_species_fallback": used_fallback,
    }
    # An empty pathway set is ambiguous - genuinely no Reactome data for this
    # organism, or a fetch/parse problem - so don't cache it and let the next
    # run try again, rather than serving "zero pathways" for REACTOME_REFRESH_DAYS.
    if pathways:
        _save_json_cache(cache_path, result)
    elif verbose:
        print(f"[Reactome] not caching: no pathways found for org={taxon_name!r}")
    if verbose:
        print(f"[Reactome] org={taxon_name!r}: {len(pathways)} pathways (fallback={used_fallback})")
    return result


# ---------------------------------------------------------------------------
# MSigDB gene sets (native download, symbol-keyed, human + mouse)
# ---------------------------------------------------------------------------

def msigdb_collection_available(collection_key: str, taxon_id: str) -> bool:
    """True iff this collection has a real MSigDB file for this organism.
    Consults MSIGDB_COLLECTIONS's explicit per-collection map only - never
    derives a filename by substituting the organism code (see that
    constant's comment for why that would be wrong)."""
    info = MSIGDB_COLLECTIONS.get(collection_key)
    if info is None:
        return False
    suffix = MSIGDB_TAXON_SUFFIX.get(str(taxon_id).strip())
    if suffix is None:
        return False
    _, file_key = suffix
    return bool(info.get(file_key))


def _msigdb_url(collection_key: str, taxon_id: str) -> "str | None":
    info = MSIGDB_COLLECTIONS.get(collection_key)
    if info is None:
        return None
    suffix = MSIGDB_TAXON_SUFFIX.get(str(taxon_id).strip())
    if suffix is None:
        return None
    org_suffix, file_key = suffix
    base_name = info.get(file_key)
    if not base_name:
        return None
    return (
        f"{MSIGDB_BASE_URL}/{MSIGDB_VERSION}.{org_suffix}/"
        f"{base_name}.v{MSIGDB_VERSION}.{org_suffix}.symbols.gmt"
    )


def _parse_msigdb_gmt(text: str) -> dict:
    """Parse a standard GMT file (name \\t description/card-url \\t gene1 \\t
    gene2 ...), one gene set per line. The description field is MSigDB's own
    gene-set-card URL - discarded here since no result table has a link
    column yet; a cheap future enhancement, not a silent drop.

    Returns {set_name: {"name": str, "genes": [symbol, ...]}}.
    """
    sets: dict = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        set_name = parts[0].strip()
        genes = sorted({g.strip().upper() for g in parts[2:] if g.strip()})
        if set_name and genes:
            sets[set_name] = {"name": set_name, "genes": genes}
    return sets


def fetch_msigdb_sets(collection_key: str, taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """Fetch one MSigDB collection for one organism, gene-symbol-keyed.

    Caller is responsible for checking msigdb_collection_available(...) first
    and supplying a proper "not available for this organism" note instead of
    calling this when it's False - this function still degrades safely
    (fetch_failed=True) if called anyway, but its failure note is generic.

    Files are downloaded in one request and parsed straight to the small
    {name: {...}} dict this returns - the same in-memory approach the former
    Hallmark fetcher used - because what gets cached to disk is that small
    parsed JSON (via _save_json_cache), never the raw GMT text at a fixed
    path, so the tmp-file+os.replace atomic-write pattern GO/Reactome use for
    their large persistent raw-file caches does not apply here.

    Returns {"sets": {set_name: {"name": str, "genes": [symbol, ...]}},
    "fetched_at": iso str, "fetch_failed": bool}.
    """
    cache_key = f"msigdb_{collection_key}_{taxon_id}_{MSIGDB_VERSION}"
    cache_path = _universe_cache_path(cache_key)
    if not force_refresh:
        cached = _load_json_cache(cache_path, MSIGDB_REFRESH_DAYS)
        if cached is not None:
            return cached

    url = _msigdb_url(collection_key, taxon_id)
    sets: dict = {}
    fetch_ok = False
    if url is not None and HTTP_SESSION is not None:
        try:
            resp = HTTP_SESSION.get(url, headers=OMISCERA_HEADERS, timeout=120)
            if resp.status_code == 200:
                sets = _parse_msigdb_gmt(resp.content.decode("utf-8", errors="replace"))
                fetch_ok = True
            elif verbose:
                print(f"[MSigDB] {collection_key}/{taxon_id}: HTTP {resp.status_code}")
        except Exception as exc:
            if verbose:
                print(f"[MSigDB] {collection_key}/{taxon_id}: fetch failed: {exc}")

    result = {"sets": sets, "fetched_at": datetime.now(timezone.utc).isoformat(), "fetch_failed": not fetch_ok}
    # Same cache-poisoning guard every fetcher in this module uses: only a
    # genuinely successful, non-empty fetch is cached.
    if fetch_ok and sets:
        _save_json_cache(cache_path, result)
    elif verbose:
        print(f"[MSigDB] {collection_key}/{taxon_id}: not caching (fetch_ok={fetch_ok}, n_sets={len(sets)})")
    if verbose:
        print(f"[MSigDB] {collection_key}/{taxon_id}: {len(sets)} sets")
    return result


def _msigdb_selected_keys(msigdb_collections: "list[str] | None") -> "list[str]":
    """Validated, canonically-ordered (MSIGDB_COLLECTIONS' own order) list of
    collection keys actually selected - silently drops any unrecognized key
    rather than raising, the same way an unknown `databases` entry is
    already dropped by the API layer's whitelist (defense in depth here,
    not the primary validation)."""
    requested = set(msigdb_collections or [])
    return [k for k in MSIGDB_COLLECTIONS if k in requested]


# ---------------------------------------------------------------------------
# User-uploaded background (reference) gene list — always gene-SYMBOL-keyed.
# Optional: this module's default universe is already "genes tested in this
# experiment" (see module docstring), which is the RNA-seq-appropriate
# background without any upload needed. This is an escape hatch for a user
# who wants to restrict further, e.g. to a targeted gene panel.
# ---------------------------------------------------------------------------

def load_user_background_list(content: bytes) -> set:
    """Parse an uploaded background gene list (gene symbols).

    Accepts a plain one-symbol-per-line file, or a CSV/TSV with a recognizable
    gene-symbol column header.
    """
    text = content.decode("utf-8-sig", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return set()

    id_like_headers = ("gene symbol", "gene_symbol", "symbol", "gene")
    first_line_lower = lines[0].lower()
    if any(h in first_line_lower for h in id_like_headers) and ("," in lines[0] or "\t" in lines[0]):
        delimiter = "\t" if "\t" in lines[0] else ","
        reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
        id_col = None
        for col in reader.fieldnames or []:
            if any(h in col.strip().lower() for h in id_like_headers):
                id_col = col
                break
        if id_col:
            return {
                str(row.get(id_col, "")).strip().upper()
                for row in reader
                if str(row.get(id_col, "")).strip()
            }

    symbols = set()
    for line in lines:
        token = line.split(",")[0].split("\t")[0].strip()
        if token and token.lower() not in id_like_headers:
            symbols.add(token.upper())
    return symbols


# ---------------------------------------------------------------------------
# Custom (user-uploaded) gene-set database — always gene-SYMBOL-keyed, since
# this module is symbol-native throughout (unlike Genomics, there is no
# accession/Ensembl/BLAST id-type choice to make here).
# ---------------------------------------------------------------------------

_CUSTOM_TERM_ALIASES = (
    "term_id", "term id", "termid", "term", "set", "set_id", "setid",
    "pathway", "pathway_id", "category", "gene_set", "geneset", "id",
    "gene_set_id", "gene set id",
)
_CUSTOM_NAME_ALIASES = (
    "term_name", "term name", "termname", "name", "description",
    "label", "pathway_name", "set_name", "gene_set_name", "gene set name",
)
_CUSTOM_GENE_ALIASES = (
    "gene", "gene_id", "geneid", "symbol", "gene symbol", "gene_symbol",
)


def load_custom_geneset_db(content: bytes) -> dict:
    """Parse a user-uploaded custom gene-set database into the ORA universe
    shape, always in gene-symbol space (uppercased, for case-insensitive
    matching against the resolved DE table's own gene_symbol column).

    Accepts two formats:
      * GMT: one gene set per line, tab-separated ->
            Term_ID <TAB> Description <TAB> Gene1 <TAB> Gene2 ...
      * Tidy TSV/CSV: a header naming a term column and a gene column (plus an
        optional term-name column); one (term, gene) membership per row.

    Returns {"terms": {term_id: {"name": str, "genes": set}},
             "n_terms": int, "universe_size": int}.
    Raises ValueError on an empty/unparseable file.
    """
    text = content.decode("utf-8-sig", errors="replace")
    lines = [ln for ln in (l.rstrip("\r\n") for l in text.splitlines()) if ln.strip()]
    if not lines:
        raise ValueError("The uploaded database file is empty.")

    terms: dict = {}

    header = lines[0]
    delim = "\t" if "\t" in header else ("," if "," in header else None)
    is_tidy = False
    if delim:
        cols = [c.strip().lower() for c in header.split(delim)]
        is_tidy = any(c in _CUSTOM_TERM_ALIASES for c in cols) and any(c in _CUSTOM_GENE_ALIASES for c in cols)

    if is_tidy:
        assert delim is not None  # is_tidy is only True when the `if delim:` guard above passed
        reader = csv.DictReader(io.StringIO(text), delimiter=delim)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]

        def _pick(aliases):
            for f in fieldnames:
                if f.strip().lower() in aliases:
                    return f
            return None

        term_col = _pick(_CUSTOM_TERM_ALIASES)
        name_col = _pick(_CUSTOM_NAME_ALIASES)
        gene_col = _pick(_CUSTOM_GENE_ALIASES)
        if not term_col or not gene_col:
            raise ValueError("Tidy file must have both a term column and a gene column.")
        for row in reader:
            tid = str(row.get(term_col, "") or "").strip()
            gene = str(row.get(gene_col, "") or "").strip().upper()
            if not tid or not gene:
                continue
            name = str(row.get(name_col, "") or "").strip() if name_col else ""
            bucket = terms.setdefault(tid, {"name": name or tid, "genes": set()})
            if name and bucket["name"] == tid:
                bucket["name"] = name
            bucket["genes"].add(gene)
    else:
        # GMT form: Term <TAB> Description <TAB> gene1 <TAB> gene2 ...
        for ln in lines:
            parts = ln.split("\t")
            if len(parts) < 3 and "," in ln and len(ln.split(",")) >= 3:
                parts = ln.split(",")  # tolerate comma-separated GMT-ish lines
            if len(parts) < 3:
                continue
            tid = parts[0].strip()
            if not tid:
                continue
            desc = parts[1].strip()
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            if not genes:
                continue
            bucket = terms.setdefault(tid, {"name": desc or tid, "genes": set()})
            bucket["genes"].update(genes)

    if not terms:
        raise ValueError(
            "No gene sets could be parsed. Expected GMT (Term<TAB>Description<TAB>Gene1<TAB>...) "
            "or a tidy table with term and gene columns."
        )

    universe: set = set()
    for info in terms.values():
        universe |= info["genes"]

    return {"terms": terms, "n_terms": len(terms), "universe_size": len(universe)}


# ---------------------------------------------------------------------------
# Gene lengths (Ensembl BioMart) — for length-bias correction only
# ---------------------------------------------------------------------------

def _biomart_fetch_lengths(ensembl_ids: "list[str]", dataset: str, verbose: bool = False) -> dict:
    """Query Ensembl BioMart for transcript lengths of the given gene IDs, in
    chunks. Returns {ensembl_gene_id: max_transcript_length}."""
    lengths: dict = {}
    assert HTTP_SESSION is not None, "requests not installed"
    ids = list(ensembl_ids)
    for i in range(0, len(ids), BIOMART_CHUNK_SIZE):
        chunk = ids[i:i + BIOMART_CHUNK_SIZE]
        id_filter = ",".join(chunk)
        xml_query = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            "<!DOCTYPE Query>"
            '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">'
            f'<Dataset name="{dataset}" interface="default">'
            f'<Filter name="ensembl_gene_id" value="{id_filter}"/>'
            '<Attribute name="ensembl_gene_id"/>'
            '<Attribute name="transcript_length"/>'
            "</Dataset></Query>"
        )
        try:
            resp = HTTP_SESSION.get(BIOMART_URL, params={"query": xml_query}, headers=OMISCERA_HEADERS, timeout=60)
            if resp.status_code != 200:
                if verbose:
                    print(f"[BioMart] HTTP {resp.status_code}")
                continue
            for line in resp.text.splitlines():
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                gid, length_str = parts[0].strip(), parts[1].strip()
                try:
                    length = float(length_str)
                except ValueError:
                    continue
                if gid not in lengths or length > lengths[gid]:
                    lengths[gid] = length
        except Exception as exc:
            if verbose:
                print(f"[BioMart] chunk fetch failed: {exc}")
            continue
    return lengths


def fetch_gene_lengths(ensembl_ids: "list[str]", taxon_id: str, verbose: bool = False) -> dict:
    """Bulk-fetch per-gene lengths (bp) for length-bias correction. Returns
    {ensembl_gene_id: length}. Returns {} (empty) if the organism has no known
    BioMart dataset or the request fails — callers must treat that as "length
    correction unavailable" and fall back to the plain hypergeometric test,
    never raise.
    """
    dataset = BIOMART_DATASET_BY_TAXON.get(str(taxon_id).strip())
    if not dataset or HTTP_SESSION is None or not ensembl_ids:
        return {}

    cache_path = os.path.join(CACHE_DIR, GENE_LENGTH_CACHE_NAME_TMPL.format(taxon=taxon_id))
    cached_lengths: dict = {}
    if os.path.exists(cache_path):
        age_days = (time.time() - os.path.getmtime(cache_path)) / 86400
        if age_days < GENE_LENGTH_REFRESH_DAYS:
            try:
                df = pd.read_csv(cache_path, dtype=str)
                cached_lengths = {row["ensembl_gene_id"]: float(row["length"]) for _, row in df.iterrows()}
            except Exception:
                cached_lengths = {}

    wanted = set(g for g in ensembl_ids if g)
    missing = [g for g in wanted if g not in cached_lengths]
    if missing:
        if verbose:
            print(f"[BioMart] fetching lengths for {len(missing)} gene(s), dataset={dataset}")
        fetched = _biomart_fetch_lengths(missing, dataset, verbose=verbose)
        cached_lengths.update(fetched)
        try:
            os.makedirs(CACHE_DIR, exist_ok=True)
            pd.DataFrame(
                [{"ensembl_gene_id": g, "length": v} for g, v in cached_lengths.items()]
            ).to_csv(cache_path, index=False)
        except Exception:
            pass

    return {g: cached_lengths[g] for g in ensembl_ids if g in cached_lengths}


# ---------------------------------------------------------------------------
# DEG selection
# ---------------------------------------------------------------------------

def select_degs(resolved_df: pd.DataFrame, padj_cutoff: float, lfc_cutoff: float, direction: str = "all") -> pd.DataFrame:
    """Select DEGs from the resolved DE table.

    direction: "all" (either direction), "up" (log2FC >= +lfc_cutoff), or
    "down" (log2FC <= -lfc_cutoff). Rows with no gene symbol (unmapped) or a
    missing padj are never selected (NaN comparisons evaluate False).
    """
    df = resolved_df[resolved_df["gene_symbol"].astype(str).str.strip() != ""].copy()
    mask = (df["padj"] <= padj_cutoff) & (df["log2FC"].abs() >= lfc_cutoff)
    if direction == "up":
        mask &= df["log2FC"] >= lfc_cutoff
    elif direction == "down":
        mask &= df["log2FC"] <= -lfc_cutoff
    return df[mask]


# ---------------------------------------------------------------------------
# Per-database-family ORA runner (shared by GO/Reactome/Hallmark)
# ---------------------------------------------------------------------------

def _run_ora_for_sets(
    term_map: dict,
    universe_symbols: set,
    deg_symbols: set,
    min_term_size: int,
    max_term_size,
    pwf_weights: "dict | None",
) -> dict:
    """Fisher's-exact-equivalent hypergeometric (or Wallenius, if pwf_weights
    is given) + BH-FDR, one row per gene set passing the term-size gate.

    Terms with ZERO overlapping DEGs (k=0) never enter the p-value/FDR family —
    they are not "tested", they're skipped, matching clusterProfiler (the
    standard transcriptomics reference), which only runs the hypergeometric
    test on terms containing at least one query gene. Including them (an
    earlier version of this function did) inflates the family size m for every
    real test in the same table, which measurably worsens every other term's
    FDR without adding information — a term with 0 overlap is never going to
    be reported as enriched regardless of its (fixed) p=1.0.

    The hypergeometric universe (N) and DEG count (n) are restricted to genes
    ANNOTATED in this database (present in at least one term). A gene with no
    annotation here can never fall inside any term, so counting it in N/n would
    bias every p-value; this matches the clusterProfiler/goseq convention. The
    raw (unrestricted) tested-gene and DEG totals are reported separately for
    transparency.
    """
    n_universe_total = len(universe_symbols)
    n_deg_total = len(deg_symbols)

    # Annotated universe: tested genes that appear in at least one term here.
    annotated_genes: set = set()
    for info in term_map.values():
        annotated_genes.update(info["genes"])
    annotated_universe = universe_symbols & annotated_genes

    N = len(annotated_universe)
    n = len(deg_symbols & annotated_universe)

    if n <= 0:
        # No DEGs in this direction are annotated in this database at all -
        # every term would be degenerate (k is necessarily 0 too, since
        # overlap ⊆ deg_symbols ∩ annotated_universe = ∅). An earlier version
        # emitted one p=NaN/fold=NaN row per surviving term here instead of
        # returning early - for a large family like GO:BP that meant
        # thousands of NaN rows, and running benjamini_hochberg over them
        # silently turned every FDR into 1.0 (Python's min(1.0, nan) == 1.0),
        # producing a huge, meaningless table. Return a proper empty block
        # with a note instead; n/N are still reported for transparency.
        return {
            "table": pd.DataFrame(columns=RESULT_COLUMNS),
            "n_terms_tested": 0,
            "n_terms_excluded_min_size": 0,
            "n_terms_no_overlap": 0,
            "universe_size": N,
            "n_deg_in_universe": 0,
            "n_universe_total": n_universe_total,
            "n_deg_total": n_deg_total,
            "note": "None of the DEGs in this direction are annotated in this database, so no terms could be tested.",
        }

    rows = []
    n_excluded_min_size = 0
    n_terms_no_overlap = 0
    for term_id, info in term_map.items():
        # term genes ⊆ annotated_genes, so intersecting with the annotated
        # universe is the same set as intersecting with the full universe.
        set_genes = set(info["genes"]) & annotated_universe
        K = len(set_genes)
        if K < min_term_size or (max_term_size and K > max_term_size):
            n_excluded_min_size += 1
            continue
        overlap = set_genes & deg_symbols
        k = len(overlap)
        if k == 0:
            n_terms_no_overlap += 1
            continue
        elif pwf_weights:
            odds = _category_odds(pwf_weights, set_genes, annotated_universe)
            p_val, fold = wallenius_pvalue(k, n, K, N, odds), fold_enrichment(k, n, K, N)
        else:
            p_val, fold = hypergeometric_pvalue(k, n, K, N), fold_enrichment(k, n, K, N)
        rows.append({
            "Term ID": term_id, "Term Name": format_term_name(info.get("name", term_id)),
            "k": k, "K": K, "n": n, "N": N,
            "Fold Enrichment": fold, "p-value": p_val,
            # Comma-joined (not "; ") to match the frontend's shared _oraGeneSet
            # parser (static/app.js), which the reused Genomics chart builders
            # (enrichment map, gene-concept network) depend on to split this
            # column back into a gene set.
            "Overlapping Genes": ", ".join(sorted(overlap)),
        })

    if rows:
        pvals = [r["p-value"] for r in rows]
        fdrs = benjamini_hochberg(pvals)
        bonfs = bonferroni(pvals)
        for row, fdr, bonf in zip(rows, fdrs, bonfs):
            row["FDR"] = fdr
            row["Bonferroni"] = bonf

    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    if not df.empty:
        df = df.sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)

    return {
        "table": df,
        "n_terms_tested": len(rows),
        "n_terms_excluded_min_size": n_excluded_min_size,
        "n_terms_no_overlap": n_terms_no_overlap,
        "universe_size": N,
        "n_deg_in_universe": n,
        "n_universe_total": n_universe_total,
        "n_deg_total": n_deg_total,
    }


# ---------------------------------------------------------------------------
# ORA orchestrator
# ---------------------------------------------------------------------------

def run_ora(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    taxon_name: str,
    padj_cutoff: float = PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = LFC_CUTOFF_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE_DEFAULT,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    use_length_correction: bool = True,
    include_iea: bool = True,
    force_refresh: bool = False,
    custom_db: "dict | None" = None,
    background_symbols: "set | None" = None,
    msigdb_collections: "list[str] | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run ORA for the requested databases against the resolved DE table.

    databases: subset of ["GO", "Reactome", "MSigDB", "Custom database"].
    msigdb_collections: required (non-empty) when "MSigDB" is selected - a
        list of MSIGDB_COLLECTIONS keys. Each becomes its OWN independently
        BH-corrected result family (like "Reactome (up)"), e.g.
        "MSigDB: Hallmark (up)", "MSigDB: Cell type (up)". A collection
        unavailable for this run's organism (msigdb_collection_available)
        gets an explanatory empty block instead of a silent zero-result table.
    custom_db: required when "Custom database" is selected - the dict returned
        by load_custom_geneset_db (already parsed by the caller at upload
        time, same pattern as Genomics' custom database). Reuses the exact
        same _run_ora_for_sets engine as GO/Reactome/MSigDB, including this
        module's "universe = tested genes" convention (the custom database's
        own genes are intersected with the resolved DE table's tested genes,
        not used as an independent universe).
    background_symbols: optional gene-symbol set (from load_user_background_list).
        This module's default universe is already "genes tested in this
        experiment" (see module docstring) - the RNA-seq-appropriate
        background with no upload needed. When supplied, the universe is
        additionally restricted to the intersection with this set, for a
        user who wants a narrower reference (e.g. a targeted panel) than
        "everything tested." Applied once, globally, before any per-database
        block - all databases (including a custom one) share this narrowed
        universe, unlike Genomics where background applies per-database.
    include_iea: GO only. IEA (electronic, uncurated) annotations are included
        by default; set False to restrict GO to curated (non-IEA) evidence.
        Reactome's own bulk file also carries an evidence-code column, so
        include_iea is applied there too - it does NOT carry an IEA concept
        for MSigDB, which has no evidence-code data at all.
    force_refresh: bypass and overwrite the on-disk gene-set cache for GO/
        Reactome/MSigDB. Use this to bust a cache poisoned by a past
        transient fetch failure (normally that no longer happens - see
        fetch_go_sets/fetch_msigdb_sets, which now refuse to cache a
        partial/failed result in the first place - but this remains the
        manual escape hatch for an already-poisoned on-disk cache).

    Returns:
        {
          "tables": {result_key: {"table": DataFrame, "n_terms_tested": int,
                                   "n_terms_excluded_min_size": int,
                                   "n_terms_no_overlap": int,
                                   "universe_size": int, "n_deg_in_universe": int,
                                   "length_correction_used": bool,
                                   "length_correction_reason": str,
                                   "note": str | None}, ...},
          "summary": {"n_universe", "n_deg_all", "n_deg_up", "n_deg_down",
                      "length_correction_used", "length_correction_reason"},
        }
    result_key examples: "GO: Biological Process (all)", "Reactome (up)". The
    top-level summary's length_correction_used/reason reflect the "all"
    direction specifically (the primary/default view); each table block also
    carries its own direction-specific values, since "up"/"down"/"all" are
    different gene lists that are now fit and corrected independently (see
    fit_pwf).
    """
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
        # Also restrict the DataFrame used below for gene-length lookups, so
        # the length-bias PWF is fit on the SAME (background-restricted)
        # gene population the statistical test actually uses - not the wider
        # unrestricted set of tested genes.
        universe_df = universe_df[universe_df["gene_symbol"].astype(str).str.upper().isin(universe_symbols_all)]

    directions = {}
    n_deg_before_background = {}
    for direction_key in DIRECTIONS:
        deg_df = select_degs(resolved_df, padj_cutoff, lfc_cutoff, direction_key)
        raw_deg_set = set(deg_df["gene_symbol"].astype(str).str.upper())
        n_deg_before_background[direction_key] = len(raw_deg_set)
        # Restrict to the (possibly background-restricted) universe, matching
        # universe_symbols_all above - without this, a DEG excluded by the
        # user's background list still counted toward n_deg_total/the summary
        # cards, while n_deg_in_universe correctly excluded it, making
        # "DEGs not annotated" (= n_deg_total - n_deg_in_universe) silently
        # count "excluded by background" as "not annotated in this database".
        # When no background is supplied this is a no-op: select_degs already
        # filters to non-blank gene_symbol rows, the same filter
        # universe_symbols_all uses, so raw_deg_set is already a subset.
        directions[direction_key] = raw_deg_set & universe_symbols_all

    # --- Length-bias correction setup: fetch gene lengths ONCE (the expensive
    # network call), then fit a SEPARATE PWF per direction - "all"/"up"/"down"
    # are different gene lists with different detection power, so a bias model
    # fit on one and reused for the others would be a mismatch (goseq refits
    # per list). Each direction that yields an essentially flat PWF (no
    # meaningful bias for that particular list) falls back to the plain
    # hypergeometric test for that direction specifically, rather than forcing
    # a degenerate correction.
    pwf_weights_by_direction: dict = {d: None for d in DIRECTIONS}
    reason_by_direction: dict = {d: "" for d in DIRECTIONS}

    if not use_length_correction:
        for d in DIRECTIONS:
            reason_by_direction[d] = "Length correction disabled."
    elif str(taxon_id) not in BIOMART_DATASET_BY_TAXON:
        for d in DIRECTIONS:
            reason_by_direction[d] = f"No gene-length source available for taxon {taxon_id}."
    elif _IsotonicRegression is None:
        for d in DIRECTIONS:
            reason_by_direction[d] = "scikit-learn is not installed."
    else:
        ensembl_by_symbol = dict(zip(
            universe_df["gene_symbol"].astype(str).str.upper(),
            universe_df["ensembl_id"].astype(str),
        ))
        ensembl_ids = [v for v in ensembl_by_symbol.values() if v]
        lengths_by_ensembl = fetch_gene_lengths(ensembl_ids, taxon_id, verbose=verbose)
        symbols_with_length = [sym for sym, ens in ensembl_by_symbol.items() if ens and ens in lengths_by_ensembl]
        if len(symbols_with_length) < len(universe_symbols_all) * 0.5:
            sparse_reason = (
                f"Gene lengths found for only {len(symbols_with_length)}/"
                f"{len(universe_symbols_all)} universe genes; too sparse to fit reliably."
            )
            for d in DIRECTIONS:
                reason_by_direction[d] = sparse_reason
        else:
            lengths_ordered = [lengths_by_ensembl[ensembl_by_symbol[s]] for s in symbols_with_length]
            for direction_key in DIRECTIONS:
                deg_set = directions[direction_key]
                if not deg_set:
                    reason_by_direction[direction_key] = "No DEGs in this direction at the current thresholds."
                    continue
                is_deg_ordered = [s in deg_set for s in symbols_with_length]
                try:
                    weights = fit_pwf(is_deg_ordered, lengths_ordered)
                except Exception as exc:
                    reason_by_direction[direction_key] = f"PWF fitting failed: {exc}"
                    continue
                if weights is None:
                    reason_by_direction[direction_key] = (
                        "No statistically significant gene-length bias detected for this "
                        "gene list; using the standard (uncorrected) test."
                    )
                    continue
                pwf_weights_by_direction[direction_key] = dict(zip(symbols_with_length, weights))

    length_used_by_direction = {d: pwf_weights_by_direction[d] is not None for d in DIRECTIONS}

    tables: dict = {}
    # GO increments `step` once per NAMESPACE x DIRECTION (3x3=9), and MSigDB
    # once per SELECTED COLLECTION x DIRECTION, not once per direction like
    # every other database - a plain len(databases)*len(DIRECTIONS) would
    # undercount both, so progress (step/steps_total) overshoots 1.0 and the
    # bar freezes at 100% before the run actually finishes.
    n_msigdb_selected = len(_msigdb_selected_keys(msigdb_collections))
    steps_total = max(
        sum(
            len(GO_NAMESPACES) if db == "GO" else
            n_msigdb_selected if db == "MSigDB" else
            1
            for db in databases
        ) * len(DIRECTIONS),
        1,
    )
    step = 0

    def _store(result_key: str, direction_key: str, block: dict) -> None:
        block["length_correction_used"] = length_used_by_direction[direction_key]
        block["length_correction_reason"] = reason_by_direction[direction_key]
        tables[result_key] = block

    def _empty_block(note: str) -> dict:
        return {
            "table": pd.DataFrame(columns=RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0, "n_terms_no_overlap": 0,
            "universe_size": 0, "n_deg_in_universe": 0,
            "n_universe_total": len(universe_symbols_all), "n_deg_total": 0,
            "note": note,
        }

    if "GO" in databases:
        go_sets = fetch_go_sets(taxon_id, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
        failed_ns = set(go_sets.get("fetch_failed_namespaces") or [])
        for ns in GO_NAMESPACES:
            ns_label = ns.replace("_", " ").title()
            ns_note = (
                f"Could not fetch GO {ns_label} from QuickGO for this run (transient "
                "service error); this namespace's terms are unavailable. Re-run the "
                "analysis to retry." if ns in failed_ns else None
            )
            for direction_key in DIRECTIONS:
                deg_symbols = directions[direction_key]
                step += 1
                _progress(step, steps_total, f"GO {ns_label} ({direction_key})")
                result_key = f"GO: {ns_label} ({direction_key})"
                if not deg_symbols:
                    _store(result_key, direction_key, _empty_block("No DEGs in this direction at the current thresholds."))
                    continue
                if ns in failed_ns:
                    assert ns_note is not None  # set above whenever `ns in failed_ns`
                    _store(result_key, direction_key, _empty_block(ns_note))
                    continue
                _store(result_key, direction_key, _run_ora_for_sets(
                    go_sets.get(ns, {}), universe_symbols_all, deg_symbols,
                    min_term_size, max_term_size,
                    pwf_weights_by_direction[direction_key],
                ))

    if "Reactome" in databases:
        if not taxon_name.strip():
            # Without a resolved organism name, Reactome's own file has no
            # species column to filter on - matching that used to mean
            # silently pooling every organism's pathways together (a real
            # bug: taxon_name was never sent by the frontend). Refuse instead.
            for direction_key in DIRECTIONS:
                step += 1
                _progress(step, steps_total, f"Reactome ({direction_key})")
                _store(f"Reactome ({direction_key})", direction_key, _empty_block(
                    "Reactome requires a resolved organism name to filter pathways by "
                    "species; none was available for this run. Re-run retrieval and "
                    "confirm an organism was selected, then re-run this analysis."
                ))
        else:
            reactome = fetch_reactome_sets(taxon_name, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
            fallback_note = None
            if reactome.get("matched_via_genus_species_fallback"):
                fallback_note = (
                    f"No exact Reactome organism match for \"{taxon_name}\"; matched via "
                    f"genus + species (\"{_genus_species(taxon_name)}\") instead."
                )
            # Reactome pathways are Entrez-keyed natively; translate to the
            # universe's symbol space via the resolved table's own entrez_id ->
            # gene_symbol mapping, so every database is tested in the same space.
            entrez_to_symbol = dict(zip(
                universe_df["entrez_id"].astype(str), universe_df["gene_symbol"].astype(str).str.upper()
            ))
            reactome_symbol_sets: dict = {}
            for pw_id, info in reactome.get("pathways", {}).items():
                syms = {entrez_to_symbol[e] for e in info["genes"] if e in entrez_to_symbol and e}
                if syms:
                    reactome_symbol_sets[pw_id] = {"name": info["name"], "genes": sorted(syms)}
            # A genuinely empty pathway set means this organism/setting has NO
            # Reactome data at all; an empty reactome_symbol_sets DESPITE
            # non-empty pathways means Reactome has data but none of it could
            # be translated to symbols (no experiment gene carries a matching
            # Entrez ID) - two different, real causes, checked in that order
            # so neither one masquerades as _run_ora_for_sets' generic
            # "DEGs aren't annotated" note, which would misattribute a
            # database/coverage problem to the gene list. Reactome's curated
            # (non-IEA) evidence is essentially human-only (verified:
            # mouse/rat/zebrafish/cow bulk file rows are 100% electronic), so
            # the IEA-excluded case is the single most likely way a
            # non-human run hits the "no pathways at all" branch.
            no_data_note = None
            if not reactome.get("pathways"):
                if not include_iea:
                    no_data_note = (
                        f"No Reactome pathways found for \"{taxon_name}\" with electronic (IEA) "
                        "annotations excluded. Reactome's curated (non-IEA) coverage is "
                        "essentially human-only - for other organisms, try including electronic "
                        "annotations."
                    )
                else:
                    no_data_note = f"No Reactome pathway data is available for \"{taxon_name}\"."
            elif not reactome_symbol_sets:
                no_data_note = (
                    f"Reactome has pathway data for \"{taxon_name}\", but none of your resolved "
                    "genes carry a matching Entrez Gene ID, so no pathway could be translated to "
                    "your gene symbols. This is an identifier-mapping gap, not a sign that your "
                    "genes aren't biologically annotated - try GO or MSigDB, which match on gene "
                    "symbol instead."
                )
            for direction_key in DIRECTIONS:
                deg_symbols = directions[direction_key]
                step += 1
                _progress(step, steps_total, f"Reactome ({direction_key})")
                result_key = f"Reactome ({direction_key})"
                if no_data_note:
                    combined = (fallback_note + " " + no_data_note) if fallback_note else no_data_note
                    _store(result_key, direction_key, _empty_block(combined))
                    continue
                if not deg_symbols:
                    empty_note = "No DEGs in this direction at the current thresholds."
                    if fallback_note:
                        empty_note = fallback_note + " " + empty_note
                    _store(result_key, direction_key, _empty_block(empty_note))
                    continue
                block = _run_ora_for_sets(
                    reactome_symbol_sets, universe_symbols_all, deg_symbols,
                    min_term_size, max_term_size,
                    pwf_weights_by_direction[direction_key],
                )
                # Append rather than overwrite: _run_ora_for_sets may already
                # have set its own note (e.g. "no DEGs annotated in this
                # database"), and silently replacing it would hide that real
                # explanation behind the organism-match caveat.
                if fallback_note:
                    block["note"] = (block["note"] + " " + fallback_note) if block.get("note") else fallback_note
                _store(result_key, direction_key, block)

    if "MSigDB" in databases:
        # Each SELECTED COLLECTION is its own independently BH-corrected
        # result family, exactly like "Reactome (up)" - see
        # msigdb_collection_available's docstring for why a collection's
        # per-organism availability is an explicit map, never derived.
        for collection_key in _msigdb_selected_keys(msigdb_collections):
            info = MSIGDB_COLLECTIONS[collection_key]
            label = f"MSigDB: {info['label']}"
            available = msigdb_collection_available(collection_key, taxon_id)
            msigdb_data = fetch_msigdb_sets(collection_key, taxon_id, verbose=verbose, force_refresh=force_refresh) if available else {"sets": {}}
            msigdb_failed = available and msigdb_data.get("fetch_failed")
            for direction_key in DIRECTIONS:
                deg_symbols = directions[direction_key]
                step += 1
                _progress(step, steps_total, f"{label} ({direction_key})")
                result_key = f"{label} ({direction_key})"
                if not available:
                    _store(result_key, direction_key, _empty_block(
                        f"MSigDB's {info['label']} collection has no gene sets published for this "
                        f"run's organism (taxon {taxon_id})."
                    ))
                    continue
                if msigdb_failed:
                    _store(result_key, direction_key, _empty_block(
                        f"Could not fetch MSigDB {info['label']} gene sets for this run (transient "
                        "service error); this collection's terms are unavailable. Re-run the "
                        "analysis to retry."
                    ))
                    continue
                if not deg_symbols:
                    _store(result_key, direction_key, _empty_block("No DEGs in this direction at the current thresholds."))
                    continue
                _store(result_key, direction_key, _run_ora_for_sets(
                    msigdb_data.get("sets", {}), universe_symbols_all, deg_symbols,
                    min_term_size, max_term_size,
                    pwf_weights_by_direction[direction_key],
                ))

    if "Custom database" in databases:
        custom_terms = (custom_db or {}).get("terms", {})
        for direction_key in DIRECTIONS:
            deg_symbols = directions[direction_key]
            step += 1
            _progress(step, steps_total, f"Custom database ({direction_key})")
            result_key = f"Custom database ({direction_key})"
            if not custom_terms:
                _store(result_key, direction_key, _empty_block(
                    "No custom gene-set database was uploaded for this run."
                ))
                continue
            if not deg_symbols:
                _store(result_key, direction_key, _empty_block("No DEGs in this direction at the current thresholds."))
                continue
            _store(result_key, direction_key, _run_ora_for_sets(
                custom_terms, universe_symbols_all, deg_symbols,
                min_term_size, max_term_size,
                pwf_weights_by_direction[direction_key],
            ))

    # How many result families in `tables` actually TESTED terms and carry
    # their own independent Benjamini-Hochberg correction - excludes empty
    # blocks (no DEGs in that direction, a database unavailable/failed for
    # this organism), since those carry no p-values and so contribute
    # nothing to the multiple-comparisons picture. This is what the results
    # screen reports so a user who scans every table for the best hit can
    # see how many independent corrections that scan actually crosses - the
    # same caveat already documented for the Regulon module, sharpened here
    # because selecting several MSigDB collections can push this well past
    # what GO+Reactome+one-database ever produced.
    n_independent_families = sum(1 for block in tables.values() if block.get("n_terms_tested", 0) > 0)

    summary = {
        "n_universe": len(universe_symbols_all),
        "n_universe_before_background": n_universe_before_background,
        "background_applied": bool(background_symbols),
        "n_deg_all": len(directions["all"]),
        "n_deg_up": len(directions["up"]),
        "n_deg_down": len(directions["down"]),
        "n_deg_all_before_background": n_deg_before_background["all"],
        "n_deg_up_before_background": n_deg_before_background["up"],
        "n_deg_down_before_background": n_deg_before_background["down"],
        "length_correction_used": length_used_by_direction["all"],
        "length_correction_reason": reason_by_direction["all"],
        "n_independent_families": n_independent_families,
    }
    return {"tables": tables, "summary": summary}


# ---------------------------------------------------------------------------
# GSEA (pre-ranked) — gene-set (gene-label) permutation.
#
# This app's Transcriptomics input is a differential-expression SUMMARY table
# (log2FC/pvalue/padj per gene) — there is no expression matrix and no sample
# class labels anywhere in the app (see Transcriptomics_module_retriever.py's
# module docstring and RESOLVED_COLUMNS). Classic GSEA permutes sample
# PHENOTYPE labels; that is impossible here. This is therefore GSEAPreranked-
# equivalent (the same approach fgsea and clusterProfiler::GSEA take from a
# ranked list): the null is built by permuting which genes fall in a set
# (gene-label permutation), not by reshuffling samples.
#
# That constraint is also what makes this tractable at GO scale. Under gene-
# label permutation, the null distribution of the enrichment score for a
# gene set depends ONLY on that set's size, never on which genes are in it -
# so one null is computed per distinct size and shared across every term of
# that size, in every database, in one run. Combined with an O(m) closed-form
# for the enrichment score (never materializing an N-length running-sum
# vector) and sharing one random draw per permutation trial across every
# size (see _null_es_shared), a full ~9,700-term GO Biological Process run at
# 1000 permutations takes single-digit seconds, not the ~10^11-operation,
# hours-long cost a naive per-term/per-permutation loop would need.
#
# LIMITATION to surface to users (not just here): gene-label permutation
# ignores inter-gene correlation within a set, which is known to be
# anti-conservative for co-regulated gene sets - reported FDR q-values can be
# optimistic for such sets. This is the standard, accepted limitation of
# every pre-ranked GSEA tool (fgsea, clusterProfiler), not a defect specific
# to this implementation, but it is real and is documented in the GSEA
# README and the results screen.
# ---------------------------------------------------------------------------

GSEA_RANKING_METRICS = {
    "signed_log_pvalue": "Signed -log10(p-value)",
    "log2fc": "log2 Fold-Change",
    "pi_value": "π-value (log2FC × -log10 p)",
    "signed_log_padj": "Signed -log10(padj)",
}
GSEA_RANKING_METRIC_DEFAULT = "signed_log_pvalue"
GSEA_MIN_TERM_SIZE_DEFAULT = 15
GSEA_MAX_TERM_SIZE_DEFAULT = 500
GSEA_N_PERMUTATIONS_DEFAULT = 1000
GSEA_N_PERMUTATIONS_CHOICES = (500, 1000, 2000, 5000)
GSEA_WEIGHT_DEFAULT = 1.0
GSEA_SEED_DEFAULT = 42
GSEA_MIN_P = 1e-300  # clamp before log10 so a reported p=0 doesn't produce inf
GSEA_MAX_CURVES_STORED = 50  # per family - full-resolution curves for every tested term would be far too much payload at GO scale
GSEA_MAX_NULL_SAMPLES = 500  # per distinct gene-set size, for the null-distribution figure

GSEA_RESULT_COLUMNS = [
    "Term ID", "Term Name", "Size", "ES", "NES", "p-value",
    "FDR q-value", "FDR (BH)", "Leading Edge Size",
    "Leading Edge Fraction", "Rank at Max", "Leading Edge Genes",
]
# Bonferroni is deliberately NOT included here (unlike RESULT_COLUMNS): on a
# permutation p-value floored at 1/(n_permutations+1) across thousands of
# gene sets it is degenerate (a column of 1.00s), so it is omitted rather
# than shipped as a misleadingly precise-looking but meaningless number.


def build_ranked_list(resolved_df: pd.DataFrame, metric: str = GSEA_RANKING_METRIC_DEFAULT) -> dict:
    """Build a descending-sorted, deterministically-tie-broken ranked gene
    list from the resolved DE table, for pre-ranked GSEA.

    metric options (GSEA_RANKING_METRICS): signed_log_pvalue (default,
    sign(log2FC) * -log10(p)), log2fc (raw effect size), pi_value
    (log2FC * -log10(p), Xiao et al.), signed_log_padj (like the default but
    using the FDR-adjusted p-value; rows with a NaN padj are dropped and
    counted, since padj is nullable per-row in the resolved table unlike
    pvalue/log2FC).

    Ties (equal metric values) break on ascending gene symbol, so the same
    input always produces the same ranking regardless of row order - GSEA's
    permutation test only has meaning if re-running it is reproducible.

    Returns {"symbols": [str,...] (descending metric order),
    "values": np.ndarray (aligned to symbols), "metric": str,
    "n_genes": int, "n_input": int, "n_dropped_padj": int, "n_ties": int,
    "zero_crossing_rank": int (count of genes with metric > 0)}.
    """
    if metric not in GSEA_RANKING_METRICS:
        metric = GSEA_RANKING_METRIC_DEFAULT
    df = resolved_df[resolved_df["gene_symbol"].astype(str).str.strip() != ""].copy()
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.upper()
    n_input = len(df)

    n_dropped_padj = 0
    if metric == "signed_log_padj":
        before = len(df)
        df = df[df["padj"].notna()]
        n_dropped_padj = before - len(df)

    log2fc = df["log2FC"].to_numpy(dtype=float)
    # sign(0) is 0 in numpy, which would silently strand a gene with exactly
    # log2FC==0 in the middle of the ranking regardless of its p-value -
    # define the boundary case as +1 explicitly instead.
    sign = np.where(log2fc >= 0, 1.0, -1.0)
    pvalue = np.clip(df["pvalue"].to_numpy(dtype=float), GSEA_MIN_P, 1.0)

    if metric == "signed_log_pvalue":
        values = sign * (-np.log10(pvalue))
    elif metric == "log2fc":
        values = log2fc
    elif metric == "pi_value":
        values = log2fc * (-np.log10(pvalue))
    else:  # signed_log_padj
        padj = np.clip(df["padj"].to_numpy(dtype=float), GSEA_MIN_P, 1.0)
        values = sign * (-np.log10(padj))

    df = df.assign(_metric=values)
    df = df.sort_values(by=["_metric", "gene_symbol"], ascending=[False, True], kind="mergesort")
    symbols = df["gene_symbol"].tolist()
    metric_values = df["_metric"].to_numpy(dtype=float)
    n_ties = int(len(metric_values) - len(np.unique(metric_values)))
    zero_crossing_rank = int(np.sum(metric_values > 0))

    return {
        "symbols": symbols, "values": metric_values, "metric": metric,
        "n_genes": len(symbols), "n_input": n_input,
        "n_dropped_padj": n_dropped_padj, "n_ties": n_ties,
        "zero_crossing_rank": zero_crossing_rank,
    }


def _es_batch(idx_sorted: np.ndarray, abs_weights_full: np.ndarray, n_genes: int):
    """Vectorized weighted running-sum enrichment score for a batch of gene
    sets, each given as a row of ascending 0-indexed hit positions.

    idx_sorted: (P, m) int array, hit positions sorted ascending per row (P
        independent gene sets/permutation draws, all of the same size m).
    abs_weights_full: (N,) array of |ranking metric|^weight for the WHOLE
        ranked list (computed once per run, shared across every call).
    Returns (es: (P,) float, rank_at_max: (P,) int 0-indexed, sign: (P,) +-1).

    Closed-form O(m) evaluation - never materializes an N-length running-sum
    vector. The running sum only changes value exactly at a hit (a jump up,
    P_hit increases by weight/N_R) or drifts linearly down between hits
    (P_miss), so its extremes occur only exactly at a hit (a candidate
    maximum, "just after" the jump) or exactly one rank before one (a
    candidate minimum, "just before" the jump) - hand-derived and cross-
    checked against a brute-force O(N) reference in this module's unit
    tests, including the off-by-one between the max branch's rank (AT the
    hit) and the min branch's rank (one BEFORE the hit).
    """
    P, m = idx_sorted.shape
    N = n_genes
    w = abs_weights_full[idx_sorted]
    C = np.cumsum(w, axis=1)
    N_R = C[:, -1].copy()
    N_R[N_R <= 0] = 1e-12
    t_prime = np.arange(m)[None, :]
    misses_before = idx_sorted - t_prime
    denom_miss = max(N - m, 1)
    es_max_component = C / N_R[:, None] - misses_before / denom_miss
    C_shifted = np.concatenate([np.zeros((P, 1)), C[:, :-1]], axis=1)
    es_min_component = C_shifted / N_R[:, None] - misses_before / denom_miss

    max_idx = np.argmax(es_max_component, axis=1)
    es_max = es_max_component[np.arange(P), max_idx]
    min_idx = np.argmin(es_min_component, axis=1)
    es_min = es_min_component[np.arange(P), min_idx]

    take_max = es_max >= -es_min
    es = np.where(take_max, es_max, es_min)
    rank_at_max_branch = idx_sorted[np.arange(P), max_idx]
    rank_at_min_branch = idx_sorted[np.arange(P), min_idx] - 1
    rank_at_max = np.where(take_max, rank_at_max_branch, np.maximum(rank_at_min_branch, 0))
    sign = np.where(es >= 0, 1, -1)
    return es, rank_at_max.astype(int), sign


def _null_es_shared(sizes_needed, abs_weights_full: np.ndarray, n_genes: int,
                     n_permutations: int, rng: "np.random.Generator", chunk: int = 200) -> dict:
    """Permutation null of the enrichment score, one array per gene-set size
    in sizes_needed, computed together and sharing one random draw per
    permutation trial across every size.

    Under gene-label permutation the null distribution for a given size
    depends ONLY on that size (see this section's module-level docstring) -
    the single fact that makes this whole approach tractable. Sharing one
    base random subset per permutation trial across sizes (a smaller size
    just takes a further random sub-prefix of a larger one) mirrors how
    classic (phenotype-permutation) GSEA reuses the SAME permutation index
    across every gene set, and was verified in this module's unit tests to
    produce a null statistically indistinguishable (matching mean/std well
    within sampling noise, and an unbiased nominal-p-value calibration under
    a pure-noise input) from independently redrawing per size, while being
    roughly an order of magnitude faster - the dominant cost of a GSEA run.

    Returns {size: (n_permutations,) float array of null ES values}.
    """
    N = n_genes
    sizes_sorted = sorted({int(s) for s in sizes_needed})
    if not sizes_sorted:
        return {}
    max_size = sizes_sorted[-1]
    results = {s: np.empty(n_permutations, dtype=float) for s in sizes_sorted}
    done = 0
    while done < n_permutations:
        p = min(chunk, n_permutations - done)
        rand_matrix = rng.random((p, N))
        # kth = max_size - 1 (not max_size): np.argpartition's kth is a valid
        # 0-indexed position INTO the array (max N-1), and the smallest
        # `max_size` elements are exactly positions [0, max_size), so the
        # partition boundary the call must guarantee correct is at index
        # max_size - 1. Passing max_size itself is off by one and raises
        # "kth out of bounds" whenever a gene set's size equals the entire
        # ranked universe (max_size == N) - a real, if rare, case with a
        # small ranked list, caught by a live run rather than any unit test
        # (which always used N large enough that max_size < N).
        idx_max = np.argpartition(rand_matrix, max_size - 1, axis=1)[:, :max_size]
        # argpartition's within-partition order is unspecified (not random) -
        # an explicit reshuffle is required before taking a further, smaller
        # prefix, or a smaller size's "random" subset would be systematically
        # biased rather than genuinely uniform (verified in this module's
        # unit tests by comparing against independently-drawn nulls).
        shuffle_key = rng.random((p, max_size))
        shuffle_order = np.argsort(shuffle_key, axis=1)
        idx_max_shuffled = np.take_along_axis(idx_max, shuffle_order, axis=1)
        for s in sizes_sorted:
            sub_sorted = np.sort(idx_max_shuffled[:, :s], axis=1)
            es, _rank, _sign = _es_batch(sub_sorted, abs_weights_full, N)
            results[s][done:done + p] = es
        done += p
    return results


def _normalize_es(observed_es: float, null_es: np.ndarray):
    """NES = ES / mean(same-signed null) - makes different-sized gene sets'
    enrichment scores comparable. Returns (nes, null_nes, nominal_p).
    nominal_p uses the +1/+1 correction so it is never exactly 0 (resolution
    floor is 1/(n_permutations+1) - a strongly enriched set cannot score
    below that regardless of how extreme its ES is, at a given permutation
    count)."""
    pos = null_es[null_es > 0]
    neg = null_es[null_es < 0]
    mean_pos = float(pos.mean()) if pos.size else 1e-12
    mean_neg = float(abs(neg.mean())) if neg.size else 1e-12
    null_nes = np.where(null_es >= 0, null_es / mean_pos, null_es / mean_neg)
    if observed_es >= 0:
        nes = observed_es / mean_pos
        same_sign_null = null_es[null_es >= 0]
        ge = int(np.sum(same_sign_null >= observed_es))
        p = (ge + 1) / (same_sign_null.size + 1)
    else:
        nes = observed_es / mean_neg
        same_sign_null = null_es[null_es <= 0]
        le = int(np.sum(same_sign_null <= observed_es))
        p = (le + 1) / (same_sign_null.size + 1)
    return float(nes), null_nes, float(p)


def _pooled_fdr(observed_nes: np.ndarray, size_of_term: np.ndarray,
                 null_nes_by_size: dict, size_counts: dict) -> np.ndarray:
    """GSEA-style pooled-permutation FDR q-value: the null pool for the
    denominator is drawn from EVERY tested term's own (size-shared) null -
    not a plain Benjamini-Hochberg correction on nominal p - which is what
    makes q-values comparable across gene sets of different sizes and to
    published GSEA/fgsea output. Positive-NES and negative-NES gene sets are
    two separate multiple-testing families (matching desktop GSEA); q is
    only guaranteed monotone WITHIN each sign, never enforced jointly across
    both signs combined (verified explicitly in this module's unit tests,
    which check each sign family separately).
    """
    S = len(observed_nes)
    q = np.zeros(S, dtype=float)
    pos_mask = observed_nes >= 0
    neg_mask = ~pos_mask

    def _one_side(mask, sign):
        if not np.any(mask):
            return
        obs = observed_nes[mask]
        obs_sorted_idx = np.argsort(obs if sign > 0 else -obs)
        obs_sorted = obs[obs_sorted_idx]
        if sign > 0:
            obs_ge_counts = len(obs_sorted) - np.searchsorted(obs_sorted, obs_sorted, side="left")
        else:
            neg_obs_sorted = np.sort(-obs)
            obs_ge_counts = len(neg_obs_sorted) - np.searchsorted(neg_obs_sorted, -obs_sorted, side="left")
        n_obs_side = len(obs)

        numerator_counts = np.zeros(len(obs_sorted), dtype=float)
        denom_null_total = 0.0
        for size, null_arr in null_nes_by_size.items():
            cnt = size_counts.get(size, 0)
            if cnt == 0:
                continue
            side_null = null_arr[null_arr >= 0] if sign > 0 else null_arr[null_arr <= 0]
            if side_null.size == 0:
                continue
            denom_null_total += cnt * side_null.size
            if sign > 0:
                sorted_null = np.sort(side_null)
                counts_ge = len(sorted_null) - np.searchsorted(sorted_null, obs_sorted, side="left")
            else:
                sorted_null_mag = np.sort(-side_null)
                counts_ge = len(sorted_null_mag) - np.searchsorted(sorted_null_mag, -obs_sorted, side="left")
            numerator_counts += cnt * counts_ge

        if denom_null_total <= 0:
            q_side = np.ones(len(obs_sorted))
        else:
            null_frac = numerator_counts / denom_null_total
            obs_frac = obs_ge_counts / max(n_obs_side, 1)
            with np.errstate(divide="ignore", invalid="ignore"):
                q_side = null_frac / np.where(obs_frac > 0, obs_frac, np.nan)
            q_side = np.nan_to_num(q_side, nan=1.0)
            q_side = np.clip(q_side, 0.0, 1.0)
            # Forward running minimum over ascending-|NES| order: q must not
            # increase as |NES| increases. A backward running min gives the
            # opposite, wrong-direction monotonicity - a bug caught only by
            # an explicit per-sign-family test in this module's unit tests,
            # not by a naive combined-sign check (which a single interleaved
            # monotone ordering across both signs would incorrectly demand).
            q_side = np.minimum.accumulate(q_side)
        q[np.where(mask)[0][obs_sorted_idx]] = q_side

    _one_side(pos_mask, 1)
    _one_side(neg_mask, -1)
    return q


def _leading_edge(hit_positions_sorted, rank_at_max: int, es_sign: int, symbols_by_position: dict) -> "list[str]":
    """Leading-edge subset: set members at rank <= rank_at_max (ES>=0, top-
    enriched) or rank >= rank_at_max (ES<0, bottom-enriched)."""
    if es_sign >= 0:
        edge_positions = [p for p in hit_positions_sorted if p <= rank_at_max]
    else:
        edge_positions = [p for p in hit_positions_sorted if p >= rank_at_max]
    return sorted(symbols_by_position[p] for p in edge_positions)


def _full_running_sum(hit_positions_sorted, abs_weights_full: np.ndarray, n_genes: int) -> np.ndarray:
    """O(N) full running-sum curve for ONE term - only ever called for the
    small number of top-ranked terms actually stored for the enrichment plot
    (GSEA_MAX_CURVES_STORED), never for every tested term (thousands of full
    N-length curves would be far too much payload/compute at GO scale)."""
    N = n_genes
    hit_set = set(hit_positions_sorted)
    m = len(hit_positions_sorted)
    N_R = float(sum(abs_weights_full[i] for i in hit_positions_sorted))
    if N_R <= 0:
        N_R = 1e-12
    denom_miss = max(N - m, 1)
    rs = np.empty(N, dtype=float)
    running = 0.0
    for j in range(N):
        if j in hit_set:
            running += abs_weights_full[j] / N_R
        else:
            running -= 1.0 / denom_miss
        rs[j] = running
    return rs


def _downsample_curve(values: np.ndarray, keep_rank: int, target_points: int = 2000):
    """Downsample a length-N curve to ~target_points for a compact API
    payload, always keeping the exact endpoints and `keep_rank` (the ES peak,
    or the ranking-metric zero-crossing) so a plotted feature is never lost
    to downsampling. Returns (ranks: list[int], values: list[float])."""
    N = len(values)
    if N <= target_points:
        idx = np.arange(N)
    else:
        idx = np.linspace(0, N - 1, target_points).astype(int)
        keep_rank = max(0, min(int(keep_rank), N - 1))
        idx = np.unique(np.concatenate([idx, [0, N - 1, keep_rank]]))
    return idx.tolist(), np.asarray(values)[idx].tolist()


def _thin_null(null_nes: np.ndarray, target: int = GSEA_MAX_NULL_SAMPLES) -> "list[float]":
    """Evenly-spaced ORDER STATISTICS of a permutation null - a deterministic
    quantile thinning whose empirical CDF matches the full null to within
    1/target, so a histogram of it is faithful to the full-resolution null
    used for the actual p-value/FDR computation.

    Deliberately NOT a random subsample: drawing from the run's `rng` here
    would advance its state and change every downstream family's results,
    silently breaking same-seed reproducibility. Sorting first and taking an
    evenly-spaced index set needs no randomness at all.
    """
    arr = np.sort(np.asarray(null_nes, dtype=float))
    if arr.size > target:
        arr = arr[np.linspace(0, arr.size - 1, target).astype(int)]
    return np.round(arr, 4).tolist()


def _run_gsea_for_sets(
    term_map: dict,
    ranked: dict,
    min_term_size: int,
    max_term_size,
    n_permutations: int,
    weight: float,
    rng: "np.random.Generator",
) -> dict:
    """Run pre-ranked GSEA for one gene-set family against the ranked list.

    term_map: {term_id: {"name": str, "genes": iterable of SYMBOLS}} - genes
        may be a list or a set (both fetchers and the custom-DB loader are
        supported; see this module's other gene-set functions).
    ranked: the dict returned by build_ranked_list.

    Terms whose overlap with the ranked list falls outside
    [min_term_size, max_term_size] are excluded before any permutation work -
    matching run_ora's term-size gate.

    Returns {"table": DataFrame (GSEA_RESULT_COLUMNS), "n_terms_tested": int,
    "n_terms_excluded_min_size": int, "universe_size": int,
    "curves": {term_id: {"ranks": [...], "values": [...], "rank_at_max": int,
    "sign": int, "hit_ranks": [...]}} for the top GSEA_MAX_CURVES_STORED terms
    by FDR q-value,
    "null_nes_by_size": {str(size): [thinned null NES samples, ascending]}
    for each distinct size among those stored curves - shared by size, so two
    equal-sized terms have the identical entry}.
    """
    symbols = ranked["symbols"]
    values = ranked["values"]
    n_genes = len(symbols)
    position_by_symbol = {s: i for i, s in enumerate(symbols)}
    symbols_by_position = {i: s for i, s in enumerate(symbols)}
    abs_weights_full = np.abs(values) ** weight

    per_term_hit_positions = {}
    n_excluded = 0
    sizes_needed = []
    for term_id, info in term_map.items():
        genes = info.get("genes")
        gene_set = genes if isinstance(genes, set) else set(genes)
        hit_positions = sorted(position_by_symbol[g] for g in gene_set if g in position_by_symbol)
        size = len(hit_positions)
        if size < min_term_size or (max_term_size and size > max_term_size):
            n_excluded += 1
            continue
        per_term_hit_positions[term_id] = (info.get("name", term_id), hit_positions)
        sizes_needed.append(size)

    if not per_term_hit_positions:
        return {
            "table": pd.DataFrame(columns=GSEA_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": n_excluded,
            "universe_size": n_genes, "curves": {}, "null_nes_by_size": {},
        }

    null_cache = _null_es_shared(sizes_needed, abs_weights_full, n_genes, n_permutations, rng)

    observed = []
    for term_id, (name, hit_positions) in per_term_hit_positions.items():
        size = len(hit_positions)
        idx = np.array(hit_positions, dtype=int)[None, :]
        es_arr, rank_arr, sign_arr = _es_batch(idx, abs_weights_full, n_genes)
        es_val, rank_val, sign_val = float(es_arr[0]), int(rank_arr[0]), int(sign_arr[0])
        nes_val, _null_nes, p_val = _normalize_es(es_val, null_cache[size])
        observed.append({
            "term_id": term_id, "name": name, "size": size,
            "es": es_val, "nes": nes_val, "nominal_p": p_val,
            "rank_at_max": rank_val, "sign": sign_val, "hit_positions": hit_positions,
        })

    observed_nes_arr = np.array([o["nes"] for o in observed])
    size_arr = np.array([o["size"] for o in observed])
    size_counts: dict = {}
    for s in size_arr:
        size_counts[int(s)] = size_counts.get(int(s), 0) + 1

    null_nes_by_size = {}
    for size, null_arr in null_cache.items():
        pos = null_arr[null_arr > 0]
        neg = null_arr[null_arr < 0]
        mean_pos = float(pos.mean()) if pos.size else 1e-12
        mean_neg = float(abs(neg.mean())) if neg.size else 1e-12
        null_nes_by_size[size] = np.where(null_arr >= 0, null_arr / mean_pos, null_arr / mean_neg)

    q_values = _pooled_fdr(observed_nes_arr, size_arr, null_nes_by_size, size_counts)
    nominal_ps = [o["nominal_p"] for o in observed]
    fdr_bh = benjamini_hochberg(nominal_ps)

    rows = []
    for o, q, bh in zip(observed, q_values, fdr_bh):
        edge_genes = _leading_edge(o["hit_positions"], o["rank_at_max"], o["sign"], symbols_by_position)
        rows.append({
            "Term ID": o["term_id"], "Term Name": format_term_name(o["name"]), "Size": o["size"],
            "ES": o["es"], "NES": o["nes"], "p-value": o["nominal_p"],
            "FDR q-value": float(q), "FDR (BH)": bh,
            "Leading Edge Size": len(edge_genes),
            "Leading Edge Fraction": (len(edge_genes) / o["size"]) if o["size"] else 0.0,
            "Rank at Max": o["rank_at_max"],
            # Comma-joined to match the frontend's shared _oraGeneSet parser
            # (static/app.js), so the enrichment-map/gene-concept-network
            # builders work on leading-edge genes with no changes needed.
            "Leading Edge Genes": ", ".join(edge_genes),
        })

    df = pd.DataFrame(rows, columns=GSEA_RESULT_COLUMNS)
    curves: dict = {}
    # Keyed by str(size) - JSON object keys are always strings, and the
    # client must look this up via String(row["Size"]). Only populated for
    # sizes that actually have a stored curve (<= GSEA_MAX_CURVES_STORED
    # terms), not every distinct size in this family, since the null
    # histogram is only ever shown alongside a term that has a curve.
    null_nes_samples: dict = {}
    if not df.empty:
        df = df.assign(_abs_nes=df["NES"].abs())
        df = df.sort_values(by=["FDR q-value", "_abs_nes"], ascending=[True, False], kind="stable")
        df = df.drop(columns=["_abs_nes"]).reset_index(drop=True)

        hit_positions_by_term = {o["term_id"]: o["hit_positions"] for o in observed}
        rank_sign_by_term = {o["term_id"]: (o["rank_at_max"], o["sign"]) for o in observed}
        size_by_term = {o["term_id"]: o["size"] for o in observed}
        for term_id in df["Term ID"].head(GSEA_MAX_CURVES_STORED):
            hp = hit_positions_by_term[term_id]
            rank_at_max, sign_val = rank_sign_by_term[term_id]
            rs = _full_running_sum(hp, abs_weights_full, n_genes)
            ranks_ds, values_ds = _downsample_curve(rs, rank_at_max)
            curves[term_id] = {
                "ranks": ranks_ds, "values": values_ds,
                "rank_at_max": rank_at_max, "sign": sign_val,
                "hit_ranks": hp,
            }
            # The null is shared by SIZE (that is the whole point of
            # _null_es_shared - see its docstring) - two equal-sized terms
            # legitimately reuse the identical histogram with a different
            # observed-NES marker, so this is only computed once per
            # distinct size encountered among the stored curves.
            size_key = str(size_by_term[term_id])
            if size_key not in null_nes_samples:
                null_nes_samples[size_key] = _thin_null(null_nes_by_size[size_by_term[term_id]])

    return {
        "table": df,
        "n_terms_tested": len(rows),
        "n_terms_excluded_min_size": n_excluded,
        "universe_size": n_genes,
        "curves": curves,
        "null_nes_by_size": null_nes_samples,
    }


def run_gsea(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    taxon_name: str,
    metric: str = GSEA_RANKING_METRIC_DEFAULT,
    min_term_size: int = GSEA_MIN_TERM_SIZE_DEFAULT,
    max_term_size=GSEA_MAX_TERM_SIZE_DEFAULT,
    n_permutations: int = GSEA_N_PERMUTATIONS_DEFAULT,
    weight: float = GSEA_WEIGHT_DEFAULT,
    seed: int = GSEA_SEED_DEFAULT,
    include_iea: bool = True,
    force_refresh: bool = False,
    custom_db: "dict | None" = None,
    background_symbols: "set | None" = None,
    msigdb_collections: "list[str] | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run pre-ranked GSEA for the requested databases against the resolved
    DE table's ranked gene list (see this section's module-level docstring
    for why gene-set/gene-label permutation is the only option here, and
    what that implies).

    Bidirectional by construction: NES sign IS the direction (positive =
    enriched toward the top/up-regulated end, negative = toward the bottom/
    down-regulated end) - unlike run_ora there is no separate all/up/down
    loop, and result keys carry no direction suffix.

    databases: subset of ["GO", "Reactome", "MSigDB", "Custom database"].
    msigdb_collections: required (non-empty) when "MSigDB" is selected - see
        run_ora's docstring; identical contract here, just without the
        all/up/down split (each selected collection is still its own family,
        e.g. "MSigDB: Cell type", but with no direction suffix).
    metric: one of GSEA_RANKING_METRICS.
    seed: reseeds an independent RNG for this run only - reruns with the same
        inputs and seed reproduce identical results.

    Returns:
        {"tables": {family_label: <_run_gsea_for_sets block, or an
                    _empty_block(note) when the family couldn't run>, ...},
         "ranking_profile": {"ranks": [...], "values": [...], "n_genes": int,
                              "zero_crossing_rank": int} - the whole ranked
                              list's metric curve, downsampled once per run,
                              for the frontend's ranked-metric strip,
         "ranking_full": {"symbols": [...], "values": [...], "n_genes": int} -
                          the FULL-resolution ranked list (every gene, not
                          downsampled). Every stored curve's "hit_ranks"
                          indexes directly into this array pair - this is
                          what lets the ridgeline/leading-edge-heatmap
                          figures show each gene set's member metric values
                          without the client re-deriving the ranking (it
                          cannot do so exactly: background-list filtering
                          happens above, before build_ranked_list, and the
                          client only ever sees a filename/count summary of
                          that list, never the symbols themselves),
         "summary": {"n_genes_ranked", "n_input", "n_dropped_padj", "n_ties",
                     "zero_crossing_rank", "metric", "n_permutations", "seed",
                     "n_independent_families"}}
    family_label examples: "GO: Biological Process", "Reactome",
    "MSigDB: Hallmark", "MSigDB: Cell type", "Custom database".
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universe_df = resolved_df[resolved_df["gene_symbol"].astype(str).str.strip() != ""]
    if background_symbols:
        bg = {str(s).strip().upper() for s in background_symbols}
        universe_df = universe_df[universe_df["gene_symbol"].astype(str).str.upper().isin(bg)]

    ranked = build_ranked_list(universe_df, metric=metric)
    rng = np.random.default_rng(seed)

    def _empty_block(note: str) -> dict:
        return {
            "table": pd.DataFrame(columns=GSEA_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
            "universe_size": ranked["n_genes"], "curves": {}, "null_nes_by_size": {}, "note": note,
        }

    tables: dict = {}
    msigdb_keys_selected = _msigdb_selected_keys(msigdb_collections)
    steps_total = max(
        sum(
            len(GO_NAMESPACES) if db == "GO" else
            len(msigdb_keys_selected) if db == "MSigDB" else
            1
            for db in databases
        ),
        1,
    )
    step = 0

    if ranked["n_genes"] < 2:
        for db in databases:
            if db == "GO":
                labels = [f"GO: {ns.replace('_', ' ').title()}" for ns in GO_NAMESPACES]
            elif db == "MSigDB":
                labels = [f"MSigDB: {MSIGDB_COLLECTIONS[k]['label']}" for k in msigdb_keys_selected]
            else:
                labels = [db]
            for label in labels:
                step += 1
                _progress(step, steps_total, f"GSEA: {label}")
                tables[label] = _empty_block(
                    "Fewer than 2 genes have a computable ranking metric for this run; GSEA cannot proceed."
                )
        return {
            "tables": tables,
            "ranking_profile": {"ranks": [], "values": [], "n_genes": 0, "zero_crossing_rank": 0},
            "ranking_full": {"symbols": [], "values": [], "n_genes": 0},
            "summary": {
                "n_genes_ranked": 0, "n_input": ranked["n_input"], "n_dropped_padj": ranked["n_dropped_padj"],
                "n_ties": 0, "zero_crossing_rank": 0, "metric": metric,
                "n_permutations": n_permutations, "seed": seed, "n_independent_families": 0,
            },
        }

    if "GO" in databases:
        go_sets = fetch_go_sets(taxon_id, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
        failed_ns = set(go_sets.get("fetch_failed_namespaces") or [])
        for ns in GO_NAMESPACES:
            ns_label = ns.replace("_", " ").title()
            step += 1
            _progress(step, steps_total, f"GSEA: GO {ns_label}")
            result_key = f"GO: {ns_label}"
            if ns in failed_ns:
                tables[result_key] = _empty_block(
                    f"Could not fetch GO {ns_label} from QuickGO for this run (transient "
                    "service error); this namespace's terms are unavailable. Re-run the "
                    "analysis to retry."
                )
                continue
            tables[result_key] = _run_gsea_for_sets(
                go_sets.get(ns, {}), ranked, min_term_size, max_term_size, n_permutations, weight, rng,
            )

    if "Reactome" in databases:
        step += 1
        _progress(step, steps_total, "GSEA: Reactome")
        if not taxon_name.strip():
            tables["Reactome"] = _empty_block(
                "Reactome requires a resolved organism name to filter pathways by species; "
                "none was available for this run. Re-run retrieval and confirm an organism "
                "was selected, then re-run this analysis."
            )
        else:
            reactome = fetch_reactome_sets(taxon_name, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
            fallback_note = None
            if reactome.get("matched_via_genus_species_fallback"):
                fallback_note = (
                    f"No exact Reactome organism match for \"{taxon_name}\"; matched via "
                    f"genus + species (\"{_genus_species(taxon_name)}\") instead."
                )
            entrez_to_symbol = dict(zip(
                universe_df["entrez_id"].astype(str), universe_df["gene_symbol"].astype(str).str.upper()
            ))
            reactome_symbol_sets: dict = {}
            for pw_id, info in reactome.get("pathways", {}).items():
                syms = {entrez_to_symbol[e] for e in info["genes"] if e in entrez_to_symbol and e}
                if syms:
                    reactome_symbol_sets[pw_id] = {"name": info["name"], "genes": sorted(syms)}
            no_data_note = None
            if not reactome.get("pathways"):
                no_data_note = (
                    f"No Reactome pathways found for \"{taxon_name}\" with electronic (IEA) "
                    "annotations excluded. Reactome's curated (non-IEA) coverage is "
                    "essentially human-only - for other organisms, try including electronic "
                    "annotations." if not include_iea else
                    f"No Reactome pathway data is available for \"{taxon_name}\"."
                )
            elif not reactome_symbol_sets:
                no_data_note = (
                    f"Reactome has pathway data for \"{taxon_name}\", but none of your resolved "
                    "genes carry a matching Entrez Gene ID, so no pathway could be translated to "
                    "your gene symbols. This is an identifier-mapping gap, not a sign that your "
                    "genes aren't biologically annotated - try GO or MSigDB, which match on gene "
                    "symbol instead."
                )
            if no_data_note:
                combined = (fallback_note + " " + no_data_note) if fallback_note else no_data_note
                tables["Reactome"] = _empty_block(combined)
            else:
                block = _run_gsea_for_sets(
                    reactome_symbol_sets, ranked, min_term_size, max_term_size, n_permutations, weight, rng,
                )
                if fallback_note:
                    block["note"] = fallback_note
                tables["Reactome"] = block

    if "MSigDB" in databases:
        for collection_key in msigdb_keys_selected:
            info = MSIGDB_COLLECTIONS[collection_key]
            label = f"MSigDB: {info['label']}"
            step += 1
            _progress(step, steps_total, f"GSEA: {label}")
            available = msigdb_collection_available(collection_key, taxon_id)
            if not available:
                tables[label] = _empty_block(
                    f"MSigDB's {info['label']} collection has no gene sets published for this "
                    f"run's organism (taxon {taxon_id})."
                )
                continue
            msigdb_data = fetch_msigdb_sets(collection_key, taxon_id, verbose=verbose, force_refresh=force_refresh)
            if msigdb_data.get("fetch_failed"):
                tables[label] = _empty_block(
                    f"Could not fetch MSigDB {info['label']} gene sets for this run (transient "
                    "service error); this collection's terms are unavailable. Re-run the "
                    "analysis to retry."
                )
            else:
                tables[label] = _run_gsea_for_sets(
                    msigdb_data.get("sets", {}), ranked, min_term_size, max_term_size, n_permutations, weight, rng,
                )

    if "Custom database" in databases:
        step += 1
        _progress(step, steps_total, "GSEA: Custom database")
        custom_terms = (custom_db or {}).get("terms", {})
        if not custom_terms:
            tables["Custom database"] = _empty_block("No custom gene-set file was uploaded for this run.")
        else:
            tables["Custom database"] = _run_gsea_for_sets(
                custom_terms, ranked, min_term_size, max_term_size, n_permutations, weight, rng,
            )

    ranking_ranks, ranking_values = _downsample_curve(ranked["values"], ranked["zero_crossing_rank"])

    # See run_ora's identical computation for why empty blocks are excluded.
    n_independent_families = sum(1 for block in tables.values() if block.get("n_terms_tested", 0) > 0)

    summary = {
        "n_genes_ranked": ranked["n_genes"],
        "n_input": ranked["n_input"],
        "n_dropped_padj": ranked["n_dropped_padj"],
        "n_ties": ranked["n_ties"],
        "zero_crossing_rank": ranked["zero_crossing_rank"],
        "metric": metric,
        "n_permutations": n_permutations,
        "seed": seed,
        "n_independent_families": n_independent_families,
    }
    return {
        "tables": tables,
        "ranking_profile": {
            "ranks": ranking_ranks, "values": ranking_values,
            "n_genes": ranked["n_genes"], "zero_crossing_rank": ranked["zero_crossing_rank"],
        },
        # Full-resolution ranked list (every gene) - see this function's
        # docstring for why this must be shipped from the backend rather
        # than reconstructed client-side. np.round(...).tolist() (not a bare
        # .tolist()) yields native Python floats at a fixed 4 dp - np.float64
        # is not JSON-serializable, and 4 dp is far finer than the metric's
        # own meaningful precision without materially changing the sort
        # order that build_ranked_list already fixed via its tie-break.
        "ranking_full": {
            "symbols": ranked["symbols"],
            "values": np.round(ranked["values"], 4).tolist(),
            "n_genes": ranked["n_genes"],
        },
        "summary": summary,
    }
