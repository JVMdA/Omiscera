# -*- coding: utf-8 -*-
"""
Genomics statistical analysis module for Omiscera.

Hosts Over-Representation Analysis (ORA) and, going forward, other genomics
statistical analyses. All HTTP calls to external databases (GO, Reactome,
GWAS Catalog) and all statistics for these analyses live in this file.

ORA tests whether a confirmed UniProt ID list (produced upstream by
Genomics_module_retriever.py — via Gene Symbol, Ensembl Gene ID, or
Homology/BLAST resolution; this module does not care which) is
over-represented in GO terms, Reactome pathways, or GWAS Catalog trait
associations versus chance, using a one-sided Fisher's exact test with
Benjamini-Hochberg FDR correction.

Module layout
-------------
  1. Constants
  2. Stats engine            (pure functions, no I/O)
  3. Shared UniProt bulk-fetch + symbol/accession mapping (reused across
     analyses and across GO/Reactome/GWAS universe builders)
  4. ORA universe builders   (GO, Reactome, GWAS Catalog)
  5. ORA universe cache
  6. ORA orchestrator        (run_ora_for_database)


Guiding principle (Part 8 revision)
------------------------------------
Each database's annotations come from THAT database's own authoritative source:
  - GO       : QuickGO/GOA downloadSearch API  (EBI)
  - Reactome : UniProt2Reactome_All_Levels.txt  (Reactome.org)
  - GWAS     : GWAS Catalog bulk associations   (EBI FTP)

UniProt is used ONLY for (a) species selection -> NCBI taxon ID/name (upstream,
in API.py) and (b) gene resolution (Genomics_module_retriever.py). It is NOT
the source of any GO term, pathway, or trait annotation in this module.

GWAS universe is gene-symbol-native: the catalog uses gene symbols (MAPPED_GENE)
and so does the ORA test for GWAS — no UniProt accession bridge. GO/Reactome
universes are UniProt-accession-keyed (their files already use accessions), so
they consume the resolved input accession list.
"""

import csv
import io
import json
import os
import re
import time
import zipfile
from datetime import datetime, timezone

import pandas as pd

import go_hierarchy

try:
    import requests as _requests
    HTTP_SESSION = _requests.Session()
except ImportError:
    _requests = None
    HTTP_SESSION = None

try:
    from scipy.stats import fisher_exact as _scipy_fisher_exact
except ImportError:
    _scipy_fisher_exact = None

try:
    from scipy.stats import mannwhitneyu as _scipy_mannwhitneyu
except ImportError:
    _scipy_mannwhitneyu = None

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# UniProt (still used for gene resolution upstream; kept for SYMBOL_FIELDS reference)
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PAGE_SIZE = 500
SYMBOL_FIELDS = "accession,gene_names"

# GO — QuickGO/GOA downloadSearch (native GO annotation source)
# NOTE (verified live): the downloadSearch `aspect` param requires the FULL
# namespace name (biological_process/molecular_function/cellular_component);
# single letters (P/F/C) are rejected with HTTP 400. Accept header must be
# `text/tsv`. The TSV carries NO GO-term-name column (includeFields=goName is
# ignored), so term names are sourced separately from go-basic.obo.
QUICKGO_DOWNLOAD_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch"
QUICKGO_DOWNLOAD_LIMIT = 2_000_000

GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GO_OBO_LOCAL_NAME = "go-basic.obo"
GO_OBO_REFRESH_DAYS = 30

GO_NAMESPACE_MAP = {"P": "biological_process", "F": "molecular_function", "C": "cellular_component"}
GO_NAMESPACES = ("biological_process", "molecular_function", "cellular_component")

# Reactome — UniProt2Reactome_All_Levels.txt bulk file (native Reactome source)
REACTOME_UNIPROT_URL = "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
REACTOME_LOCAL_NAME = "reactome_uniprot2reactome.txt"
REACTOME_REFRESH_DAYS = 30
# Bumped once: a caller that omitted taxon_name (Genomics_network.enrich_modules,
# fixed to always pass it now) used to build+cache an ALL-SPECIES universe
# under the SAME (database, taxon_id) cache key every other caller uses -
# measured live: 91,531 accessions pooled across every organism vs 12,155 for
# human alone, a >7x inflation silently shared with ORA/Burden's own Reactome
# runs for that taxon. Bumping this invalidates any such pre-fix entry rather
# than requiring a manual cache purge; bump again if this class of bug recurs.
REACTOME_UNIVERSE_CACHE_VERSION = 2

# GWAS Catalog — bulk associations file (native GWAS source; gene-symbol-native)
GWAS_RELEASE_DIR = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
GWAS_ASSOC_FILENAME = "gwas-catalog-associations_ontology-annotated-full.zip"
GWAS_LOCAL_TSV_NAME = "gwas_catalog_associations.tsv"
GWAS_REFRESH_DAYS = 30
GWAS_HUMAN_TAXON_ID = "9606"

MIN_TERM_SIZE = 3
# clusterProfiler's own maxGSSize default. With GO ancestor propagation now
# in effect, near-root terms can carry tens of thousands of genes (e.g.
# "biological_process" itself covers essentially the whole annotated
# proteome) - untested, those would report FDR-significant, fold~1.0
# "enrichment" for a completely uninformative category. Still user-editable
# per run via the config screen's term-size fields.
MAX_TERM_SIZE_DEFAULT = 500
FDR_ALPHA_DEFAULT = 0.05

OMISCERA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositorio do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "application/json",
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
ORA_UNIVERSE_CACHE_DIR = os.path.join(CACHE_DIR, "ora_universe")
ORA_UNIVERSE_INDEX_PATH = os.path.join(ORA_UNIVERSE_CACHE_DIR, "_index.csv")
GWAS_BULK_TSV_PATH = os.path.join(CACHE_DIR, GWAS_LOCAL_TSV_NAME)
REACTOME_BULK_PATH = os.path.join(CACHE_DIR, REACTOME_LOCAL_NAME)
GO_OBO_PATH = os.path.join(CACHE_DIR, GO_OBO_LOCAL_NAME)

DEFAULT_TIMEOUT = 30
DEFAULT_RETRIES = 3

ORA_RESULT_COLUMNS = [
    "Term ID", "Term Name", "k", "K", "n", "N",
    "Fold Enrichment", "p-value", "FDR", "Bonferroni", "Overlapping Genes",
]

BURDEN_RESULT_COLUMNS = [
    "Term ID", "Term Name", "Genes in set", "Scored genes (universe)",
    "Mean score (set)", "Mean score (rest)", "Effect (mean difference)",
    "p-value", "FDR", "Bonferroni", "Genes",
]


# ---------------------------------------------------------------------------
# Stats engine (pure functions, no I/O — unit-testable without network)
# ---------------------------------------------------------------------------

def _ora_contingency(k: int, n: int, K: int, N: int) -> "list[list[int]]":
    """Standard ORA 2x2 contingency table: [[k, n-k], [K-k, N-n-(K-k)]]."""
    return [[k, n - k], [K - k, N - n - (K - k)]]


def fisher_exact_pvalue(k: int, n: int, K: int, N: int) -> float:
    """One-sided ('greater' = over-representation) Fisher's exact test p-value."""
    if _scipy_fisher_exact is None:
        raise RuntimeError("scipy is required for fisher_exact_pvalue")
    table = _ora_contingency(k, n, K, N)
    _stat, p_val = _scipy_fisher_exact(table, alternative="greater")  # type: ignore[misc]
    return float(p_val)  # type: ignore[arg-type]


def fold_enrichment(k: int, n: int, K: int, N: int) -> float:
    """Observed/expected ratio: (k/n) / (K/N)."""
    if n <= 0 or N <= 0 or K <= 0:
        return float("nan")
    expected = n * (K / N)
    if expected <= 0:
        return float("nan")
    return k / expected


def benjamini_hochberg(pvalues: list) -> list:
    """Benjamini-Hochberg FDR correction (manual implementation, no statsmodels).

    Returns adjusted q-values in the same order as the input list.
    """
    m = len(pvalues)
    if m == 0:
        return []

    indexed = sorted(range(m), key=lambda i: pvalues[i])
    ranked_p = [pvalues[i] for i in indexed]

    q_at_rank = [0.0] * m
    for rank, p in enumerate(ranked_p, start=1):
        q_at_rank[rank - 1] = p * m / rank

    # Enforce monotonicity: q[i] = min(q[i], q[i+1], ..., q[m-1]) via cumulative min from the end
    running_min = 1.0
    for i in range(m - 1, -1, -1):
        running_min = min(running_min, q_at_rank[i])
        q_at_rank[i] = min(running_min, 1.0)

    result = [0.0] * m
    for rank_pos, original_idx in enumerate(indexed):
        result[original_idx] = q_at_rank[rank_pos]
    return result


def bonferroni(pvalues: list) -> list:
    """Bonferroni correction: min(p * m, 1.0) for each p-value."""
    m = len(pvalues)
    return [min(p * m, 1.0) for p in pvalues]


# ---------------------------------------------------------------------------
# Shared UniProt bulk-fetch (cursor-paginated) — now used ONLY by
# _taxon_symbol_to_accession_map, which is available for future analyses that
# may need a symbol->accession bridge. GO/Reactome/GWAS no longer use it.
# ---------------------------------------------------------------------------

def _fetch_uniprot_annotations_bulk(taxon_id: str, fields: str, verbose: bool = False) -> list:
    """Cursor-paginate through reviewed UniProt entries for a taxon.

    Returns a list of raw UniProt entry dicts (each with "primaryAccession"
    and whatever cross-reference/field data was requested). Walks the
    Link: rel="next" response header until exhausted.
    """
    if HTTP_SESSION is None:
        if verbose:
            print("requests library not available; cannot fetch UniProt bulk data")
        return []

    entries = []
    url = UNIPROT_SEARCH_URL
    params = {
        "query": f"organism_id:{taxon_id} AND reviewed:true",
        "fields": fields,
        "format": "json",
        "size": UNIPROT_PAGE_SIZE,
    }
    page_num = 0
    while url:
        last_error = ""
        response = None
        for attempt in range(DEFAULT_RETRIES + 1):
            try:
                response = HTTP_SESSION.get(
                    url,
                    params=params if page_num == 0 else None,
                    headers=OMISCERA_HEADERS,
                    timeout=DEFAULT_TIMEOUT,
                )
                if response.status_code in (200, 201):
                    break
                if response.status_code in (429, 500, 502, 503, 504):
                    last_error = f"HTTP {response.status_code}"
                    if attempt < DEFAULT_RETRIES:
                        time.sleep(1.0)
                        continue
                last_error = f"HTTP {response.status_code}"
                response = None
                break
            except Exception as exc:
                last_error = str(exc)
                response = None
                if attempt < DEFAULT_RETRIES:
                    time.sleep(0.5)
                    continue
                break

        if response is None:
            if verbose:
                print(f"[UniProt bulk fetch] page {page_num} failed: {last_error}")
            break

        try:
            page_results = response.json().get("results", [])
        except Exception as exc:
            if verbose:
                print(f"[UniProt bulk fetch] JSON parse failed on page {page_num}: {exc}")
            break

        entries.extend(page_results)
        page_num += 1
        if verbose and page_num % 10 == 0:
            print(f"[UniProt bulk fetch] page {page_num}, {len(entries)} entries so far")

        next_link = response.links.get("next", {})
        url = next_link.get("url", "") or None
        # After the first request the cursor URL carries its own query string.
        params = None

    if verbose:
        print(f"[UniProt bulk fetch] taxon={taxon_id} fields={fields} -> {len(entries)} entries, {page_num} page(s)")
    return entries


def _taxon_symbol_to_accession_map(taxon_id: str, verbose: bool = False) -> dict:
    """Build a gene-symbol -> UniProt accession map for a taxon from ONE bulk fetch.

    Used to bridge GWAS Catalog's MAPPED_GENE symbols to accessions without any
    per-gene API calls. Reusable by any future analysis needing the same bridge.
    A symbol that maps to multiple accessions keeps the first seen (consistent,
    deterministic; ambiguity here doesn't need interactive resolution the way
    Genomics_module_retriever's gene-symbol resolver does, since this is only
    used to build a background universe, not to assign a single confident ID).
    """
    entries = _fetch_uniprot_annotations_bulk(taxon_id, SYMBOL_FIELDS, verbose=verbose)
    symbol_map = {}
    for entry in entries:
        accession = str(entry.get("primaryAccession", "")).strip()
        if not accession:
            continue
        genes = entry.get("genes", []) or []
        for gene in genes:
            primary = (gene.get("geneName") or {}).get("value", "")
            if primary:
                symbol_map.setdefault(primary.upper(), accession)
            for synonym in gene.get("synonyms", []) or []:
                syn_value = (synonym or {}).get("value", "")
                if syn_value:
                    symbol_map.setdefault(syn_value.upper(), accession)
    return symbol_map


# ---------------------------------------------------------------------------
# ORA universe builders
# ---------------------------------------------------------------------------

def _quickgo_safe_get(params: dict, verbose: bool) -> "bytes | None":
    """Single GET against QUICKGO_DOWNLOAD_URL with retry. Returns raw bytes or None."""
    if HTTP_SESSION is None:
        return None
    headers = dict(OMISCERA_HEADERS)
    headers["Accept"] = "text/tsv"
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            resp = HTTP_SESSION.get(
                QUICKGO_DOWNLOAD_URL, params=params, headers=headers, timeout=180, stream=True,
            )
            if resp.status_code == 200:
                return resp.content
            last_error = f"HTTP {resp.status_code}"
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < DEFAULT_RETRIES:
                time.sleep(2.0)
                continue
            break
        except Exception as exc:
            last_error = str(exc)
            if attempt < DEFAULT_RETRIES:
                time.sleep(1.0)
    if verbose:
        print(f"[QuickGO] request failed: {last_error}")
    return None


def _parse_quickgo_tsv(content: bytes, ns: str, include_iea: bool) -> tuple:
    """Parse a QuickGO downloadSearch TSV response for one aspect.

    Returns (term_map, symbol_map):
      term_map   = {go_id: {"name": str, "genes": set[accession]}}
      symbol_map = {SYMBOL_UPPER: set[accession]}  (from the SYMBOL column, used for
                   the gene-symbol fallback when an input gene has no UniProt ID)
    Strips UniProtKB isoform suffixes (e.g. P12345-2 -> P12345).
    """
    term_map: dict = {}
    symbol_map: dict = {}
    lines = content.decode("utf-8-sig", errors="replace").splitlines()
    if not lines:
        return term_map, symbol_map

    header = [h.strip().upper() for h in lines[0].split("\t")]
    col = {h: i for i, h in enumerate(header)}

    # Locate columns robustly — QuickGO TSV header names can vary slightly
    def _find(*candidates):
        for c in candidates:
            if c in col:
                return col[c]
        return None

    db_col      = _find("GENE PRODUCT DB", "DB")
    id_col      = _find("GENE PRODUCT ID", "DB OBJECT ID", "GENE PRODUCT ID (UNIPROT ACCESSION)")
    sym_col     = _find("SYMBOL", "GENE PRODUCT SYMBOL")
    go_id_col   = _find("GO TERM", "GO ID", "GO TERM ID")
    go_name_col = _find("GO NAME", "GO TERM NAME")
    ev_col      = _find("GO EVIDENCE CODE", "EVIDENCE CODE", "GO EVIDENCE")

    if id_col is None or go_id_col is None:
        return term_map, symbol_map

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")

        # Filter to UniProtKB entries only
        if db_col is not None and db_col < len(parts):
            if parts[db_col].strip() not in ("UniProtKB", ""):
                continue

        if id_col >= len(parts):
            continue
        accession = parts[id_col].strip()
        if not accession:
            continue
        # Strip isoform suffix
        if "-" in accession:
            accession = accession.split("-")[0]

        # Filter IEA if requested
        if not include_iea and ev_col is not None and ev_col < len(parts):
            if parts[ev_col].strip() == "IEA":
                continue

        if go_id_col >= len(parts):
            continue
        go_id = parts[go_id_col].strip()
        if not go_id:
            continue

        go_name = (
            parts[go_name_col].strip()
            if go_name_col is not None and go_name_col < len(parts)
            else go_id
        )

        bucket = term_map.setdefault(go_id, {"name": go_name, "genes": set()})
        bucket["genes"].add(accession)

        # Record symbol -> accession for the fallback path
        if sym_col is not None and sym_col < len(parts):
            sym = parts[sym_col].strip().upper()
            if sym:
                symbol_map.setdefault(sym, set()).add(accession)

    return term_map, symbol_map


def _ensure_go_obo_file(verbose: bool = False) -> str:
    """Download (if needed) and cache go-basic.obo. GO terms are species-independent,
    so this file is fetched once total (not per taxon). Refreshed after GO_OBO_REFRESH_DAYS."""
    if os.path.exists(GO_OBO_PATH):
        age_days = (time.time() - os.path.getmtime(GO_OBO_PATH)) / 86400
        if age_days < GO_OBO_REFRESH_DAYS:
            if verbose:
                print(f"[GO OBO] Using cached go-basic.obo ({age_days:.1f} days old)")
            return GO_OBO_PATH
        if verbose:
            print(f"[GO OBO] Cache expired ({age_days:.1f} days), re-downloading...")

    os.makedirs(CACHE_DIR, exist_ok=True)
    assert HTTP_SESSION is not None, "requests not installed"
    if verbose:
        print(f"[GO OBO] Downloading {GO_OBO_URL} ...")
    resp = HTTP_SESSION.get(GO_OBO_URL, headers=OMISCERA_HEADERS, timeout=180, stream=True)
    resp.raise_for_status()
    # Stream to a temp file and only replace the real path once fully
    # downloaded - an interrupted write straight to GO_OBO_PATH would leave a
    # truncated file with a FRESH mtime, which the age check above would then
    # happily reuse for GO_OBO_REFRESH_DAYS.
    tmp_path = GO_OBO_PATH + ".tmp"
    downloaded = 0
    with open(tmp_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
            downloaded += len(chunk)
    os.replace(tmp_path, GO_OBO_PATH)
    if verbose:
        print(f"[GO OBO] Downloaded {downloaded / 1e6:.1f} MB -> {GO_OBO_PATH}")
    return GO_OBO_PATH


def build_go_universe(
    taxon_id: str,
    include_iea: bool = True,
    verbose: bool = False,
    force_refresh: bool = False,
) -> dict:
    """Build the GO-annotated universe from QuickGO/GOA (native GO annotation source).

    Queries EBI QuickGO downloadSearch API — one call per namespace. The TSV has no
    term-name column, so names are overlaid from go-basic.obo. Each namespace's
    direct annotations are then propagated up the GO hierarchy (is_a/part_of
    ancestors, via go_hierarchy.propagate_term_map), so a gene annotated only to
    a specific child term also counts toward its parent terms — matching
    clusterProfiler/goseq/topGO/DAVID; QuickGO's API otherwise returns direct
    annotations only, which would leave every non-leaf term under-populated.

    Returns {namespace: {go_id: {"name", "genes": set}}, "universe_size", "fetched_at",
    "propagated": bool, "caveat": str|None}.
    Cache key encodes include_iea so IEA/no-IEA universes are stored separately.
    """
    cache_db = "GO_v2" if include_iea else "GO_noiea_v2"  # "_v2" = includes ancestor propagation
    if not force_refresh:
        cached = _load_universe_cache(cache_db, taxon_id)
        # A cached _v2 blob is only trustworthy if it was actually propagated -
        # one saved after a propagation failure (see below) is stamped
        # propagated=False and must be treated as a miss, or a degraded
        # (unpropagated) universe would be served indefinitely (no TTL, and
        # force_refresh is not reachable from the running app).
        if cached is not None and cached.get("propagated") is True:
            return cached

    namespaces = {ns: {} for ns in GO_NAMESPACES}
    all_annotated_genes: set = set()
    symbol_accs: dict = {}   # SYMBOL_UPPER -> set[accession]  (all candidates)
    acc_count: dict = {}     # accession -> annotation count  (canonical picker)
    failed_namespaces = []

    for ns in GO_NAMESPACES:
        params = {
            "taxonId": taxon_id,
            "aspect": ns,   # FULL namespace name required (letters P/F/C are rejected)
            "geneProductType": "protein",
            "downloadLimit": QUICKGO_DOWNLOAD_LIMIT,
        }
        if verbose:
            print(f"[GO universe] Fetching {ns} (taxon={taxon_id}, iea={include_iea}) from QuickGO...")
        raw = _quickgo_safe_get(params, verbose=verbose)
        if raw is None:
            failed_namespaces.append(ns)
            if verbose:
                print(f"[GO universe] WARNING: failed to fetch {ns}; universe will be empty for this namespace")
            continue
        term_map, sym_map = _parse_quickgo_tsv(raw, ns, include_iea)
        namespaces[ns] = term_map
        # Symbol/accession bookkeeping is done on DIRECT annotations only,
        # before ancestor propagation below - propagation would otherwise
        # inflate an accession's apparent annotation count across many
        # ancestor terms without it reflecting any new real evidence, which
        # would corrupt the "most-annotated = canonical accession" picker.
        for info in term_map.values():
            all_annotated_genes |= info["genes"]
            for acc in info["genes"]:
                acc_count[acc] = acc_count.get(acc, 0) + 1
        for sym, accs in sym_map.items():
            symbol_accs.setdefault(sym, set()).update(accs)
        if verbose:
            print(f"[GO universe] {ns}: {len(term_map)} terms")

    # Reduce each symbol to ONE representative accession — the most-annotated one,
    # which is almost always the canonical (reviewed) protein rather than a TrEMBL
    # fragment/isoform. Mapping a symbol to ALL its accessions would inflate the
    # input gene set and corrupt the statistics.
    symbol_index: dict = {}
    for sym, accs in symbol_accs.items():
        symbol_index[sym] = max(accs, key=lambda a: acc_count.get(a, 0))

    # Overlay real GO term names from go-basic.obo (TSV has no name column),
    # then propagate each namespace's direct annotations to their ancestors.
    # Ancestor propagation never introduces a gene the namespace didn't
    # already have somewhere (an ancestor's genes are always a subset of
    # genes already counted in all_annotated_genes above), so universe_size
    # is unaffected and is correctly computed from the pre-propagation loop.
    propagated_ok = False
    try:
        obo_path = _ensure_go_obo_file(verbose=verbose)
        hierarchy = go_hierarchy.parse_go_obo(obo_path, verbose=verbose)
        go_names = hierarchy["names"]
        for ns in GO_NAMESPACES:
            for go_id, info in namespaces[ns].items():
                info["name"] = go_names.get(go_id, go_id)
            namespaces[ns] = go_hierarchy.propagate_term_map(namespaces[ns], hierarchy)
        propagated_ok = True
    except Exception as exc:
        if verbose:
            print(f"[GO universe] WARNING: could not load/propagate GO hierarchy ({exc}); using direct annotations only")

    result = {
        "biological_process": namespaces["biological_process"],
        "molecular_function": namespaces["molecular_function"],
        "cellular_component": namespaces["cellular_component"],
        # symbol -> single representative accession; used for symbol fallback
        "symbol_index": symbol_index,
        "universe_size": len(all_annotated_genes),
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "propagated": propagated_ok,
    }
    if failed_namespaces:
        result["caveat"] = (
            f"Could not fetch GO data for: {', '.join(failed_namespaces)} (transient QuickGO "
            "error). Those namespace(s) are empty for this run - re-run to retry."
        )
    elif not propagated_ok:
        result["caveat"] = (
            "Could not load GO's hierarchy for this run; results reflect only direct "
            "annotations (non-leaf terms may be under-populated). Re-run to retry."
        )

    # Only cache a genuinely complete, propagated universe - a namespace fetch
    # failure or a propagation failure must not poison the cache (no TTL, and
    # force_refresh is not reachable from the running app), matching the
    # equivalent guard in Transcriptomics_module.fetch_go_sets.
    if not failed_namespaces and propagated_ok:
        _save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[GO universe] not caching: failed_namespaces={failed_namespaces}, propagated={propagated_ok}")
    if verbose:
        print(f"[GO universe] symbol_index: {len(symbol_index)} symbols")
    return result


# ---------------------------------------------------------------------------
# Reactome bulk-file helper
# ---------------------------------------------------------------------------

def _ensure_reactome_file(verbose: bool = False) -> str:
    """Download (if needed) and cache the Reactome UniProt2Reactome_All_Levels.txt.

    Re-downloaded when older than REACTOME_REFRESH_DAYS.
    Returns the path to the local file.
    """
    if os.path.exists(REACTOME_BULK_PATH):
        age_days = (time.time() - os.path.getmtime(REACTOME_BULK_PATH)) / 86400
        if age_days < REACTOME_REFRESH_DAYS:
            if verbose:
                print(f"[Reactome] Using cached file ({age_days:.1f} days old)")
            return REACTOME_BULK_PATH
        if verbose:
            print(f"[Reactome] Cache expired ({age_days:.1f} days), re-downloading...")

    os.makedirs(CACHE_DIR, exist_ok=True)
    assert HTTP_SESSION is not None, "requests not installed"
    if verbose:
        print(f"[Reactome] Downloading {REACTOME_UNIPROT_URL} ...")
    resp = HTTP_SESSION.get(REACTOME_UNIPROT_URL, headers=OMISCERA_HEADERS, timeout=120, stream=True)
    resp.raise_for_status()
    # Stream to a temp file, then replace atomically - an interrupted 98 MB
    # download straight to REACTOME_BULK_PATH would leave a truncated file
    # with a fresh mtime, reused as "fresh" for REACTOME_REFRESH_DAYS.
    tmp_path = REACTOME_BULK_PATH + ".tmp"
    downloaded = 0
    with open(tmp_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
            downloaded += len(chunk)
    os.replace(tmp_path, REACTOME_BULK_PATH)
    if verbose:
        print(f"[Reactome] Downloaded {downloaded / 1e6:.1f} MB -> {REACTOME_BULK_PATH}")
    return REACTOME_BULK_PATH


def build_reactome_universe(
    taxon_id: str,
    taxon_name: str = "",
    include_iea: bool = True,
    verbose: bool = False,
    force_refresh: bool = False,
) -> dict:
    """Build the Reactome universe from UniProt2Reactome_All_Levels.txt (native Reactome source).

    Columns: UniProt accession | Reactome pathway ID | URL | pathway name | evidence code | organism name
    Filters by organism name (taxon_name, case-insensitive), and by evidence code when include_iea=False.
    Returns {"pathways": {pathway_id: {"name", "genes": set}}, "universe_size", "fetched_at"}.
    """
    cache_db = ("Reactome" if include_iea else "Reactome_noiea") + f"_v{REACTOME_UNIVERSE_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, taxon_id)
        if cached is not None:
            return cached

    tsv_path = _ensure_reactome_file(verbose=verbose)
    target_org = taxon_name.strip().lower()
    pathways: dict = {}
    annotated_genes: set = set()

    with open(tsv_path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            accession    = parts[0].strip()
            pathway_id   = parts[1].strip()
            pathway_name = parts[3].strip()
            evidence     = parts[4].strip()
            organism     = parts[5].strip().lower()

            if not accession or not pathway_id:
                continue
            # Filter by organism name when provided
            if target_org and organism != target_org:
                continue
            # Filter IEA
            if not include_iea and evidence == "IEA":
                continue

            bucket = pathways.setdefault(pathway_id, {"name": pathway_name, "genes": set()})
            bucket["genes"].add(accession)
            annotated_genes.add(accession)

    result = {
        "pathways": pathways,
        "universe_size": len(annotated_genes),
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
    # An empty pathway set is ambiguous - genuinely no Reactome data for this
    # organism, or a fetch/parse problem - so don't cache it and let the next
    # run try again, rather than serving "zero pathways" indefinitely (no TTL).
    if pathways:
        _save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[Reactome universe] not caching: no pathways found for org={taxon_name!r}")
    if verbose:
        print(f"[Reactome universe] org={taxon_name!r}, iea={include_iea}: "
              f"{len(pathways)} pathways, {len(annotated_genes)} annotated genes")
    return result


def load_user_background_list(content: bytes) -> set:
    """Parse an uploaded background gene list (UniProt accessions).

    Accepts a plain one-ID-per-line file, or a CSV/TSV with a recognizable
    UniProt ID column header.
    """
    text = content.decode("utf-8-sig", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return set()

    # Try as delimited table first (comma or tab) with a UniProt-ID-like header.
    first_line_lower = lines[0].lower()
    id_like_headers = ("uniprot", "accession", "uniprot id", "uniprot_id")
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
                str(row.get(id_col, "")).strip()
                for row in reader
                if str(row.get(id_col, "")).strip()
            }

    # Otherwise treat as one accession per line (skip an obvious header line).
    accessions = set()
    for line in lines:
        token = line.split(",")[0].split("\t")[0].strip()
        if token and token.lower() not in id_like_headers:
            accessions.add(token)
    return accessions


# ---------------------------------------------------------------------------
# Custom (user-uploaded) gene-set database
# ---------------------------------------------------------------------------

# Header aliases for the tidy TSV/CSV form.
_CUSTOM_TERM_ALIASES = (
    "term_id", "term id", "termid", "term", "set", "set_id", "setid",
    "pathway", "pathway_id", "category", "gene_set", "geneset", "id",
    "gene_set_id", "gene set id",
)
_CUSTOM_NAME_ALIASES = (
    "term_name", "term name", "termname", "name", "description",
    "label", "pathway_name", "set_name",
    "gene_set_name", "gene set name",
)
_CUSTOM_GENE_ALIASES = (
    "gene", "gene_id", "geneid", "accession", "uniprot", "uniprot id",
    "uniprot_id", "uniprotid", "symbol", "gene symbol", "gene_symbol",
    "uniprot_accession", "uniprot accession",
    "ensembl", "ensembl_gene_id", "ensembl gene id", "ensembl_id", "ensemblid",
)


def load_custom_geneset_db(content: bytes, id_type: str = "accession") -> dict:
    """Parse a user-uploaded custom gene-set database into the ORA universe shape.

    Accepts two formats:
      * GMT: one gene set per line, tab-separated ->
            Term_ID <TAB> Description <TAB> Gene1 <TAB> Gene2 ...
      * Tidy TSV/CSV: a header naming a term column and a gene column (plus an
        optional term-name column); one (term, gene) membership per row.

    id_type: "accession" (UniProt accessions), "symbol" (gene symbols), or
             "ensembl" (Ensembl Gene IDs). Symbol and Ensembl values are stored
             uppercased for case-insensitive matching against the resolved
             input list (mirrors the symbol-native GWAS path).

    Returns {"terms": {term_id: {"name": str, "genes": set}},
             "n_terms": int, "universe_size": int, "id_type": str}.
    Raises ValueError on an empty/unparseable file.
    """
    text = content.decode("utf-8-sig", errors="replace")
    lines = [ln for ln in (l.rstrip("\r\n") for l in text.splitlines()) if ln.strip()]
    if not lines:
        raise ValueError("The uploaded database file is empty.")

    up = id_type in ("symbol", "ensembl")

    def _norm(g: str) -> str:
        g = g.strip()
        return g.upper() if up else g

    terms: dict = {}

    # Detect the tidy table form: a header naming BOTH a term column and a gene column.
    header = lines[0]
    delim = "\t" if "\t" in header else ("," if "," in header else None)
    is_tidy = False
    if delim:
        cols = [c.strip().lower() for c in header.split(delim)]
        has_term = any(c in _CUSTOM_TERM_ALIASES for c in cols)
        has_gene = any(c in _CUSTOM_GENE_ALIASES for c in cols)
        is_tidy = has_term and has_gene

    if is_tidy:
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
            gene = _norm(str(row.get(gene_col, "") or ""))
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
                parts = ln.split(",")   # tolerate comma-separated GMT-ish lines
            if len(parts) < 3:
                continue
            tid = parts[0].strip()
            if not tid:
                continue
            desc = parts[1].strip()
            genes = {_norm(g) for g in parts[2:] if g.strip()}
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

    return {
        "terms": terms,
        "n_terms": len(terms),
        "universe_size": len(universe),
        "id_type": id_type if id_type in ("symbol", "ensembl") else "accession",
    }


# ---------------------------------------------------------------------------
# GWAS Catalog universe builder
# ---------------------------------------------------------------------------

def _ensure_gwas_bulk_file(verbose: bool = False) -> str:
    """Download (if needed) and cache the GWAS Catalog full associations TSV.

    The associations ZIP (~65 MB) is stored at GWAS_BULK_TSV_PATH after
    extraction.  Re-downloaded only when older than GWAS_REFRESH_DAYS.
    Returns the path to the local TSV.
    """
    if os.path.exists(GWAS_BULK_TSV_PATH):
        age_days = (time.time() - os.path.getmtime(GWAS_BULK_TSV_PATH)) / 86400
        if age_days < GWAS_REFRESH_DAYS:
            if verbose:
                print(f"[GWAS] Using cached TSV ({age_days:.1f} days old): {GWAS_BULK_TSV_PATH}")
            return GWAS_BULK_TSV_PATH
        if verbose:
            print(f"[GWAS] Cache expired ({age_days:.1f} days), re-downloading...")

    os.makedirs(CACHE_DIR, exist_ok=True)
    url = GWAS_RELEASE_DIR + GWAS_ASSOC_FILENAME
    if verbose:
        print(f"[GWAS] Downloading {url} ...")
    assert HTTP_SESSION is not None, "requests not installed"
    resp = HTTP_SESSION.get(url, timeout=300, stream=True)
    resp.raise_for_status()

    # Stream to a temp ZIP then extract the TSV member.
    zip_path = GWAS_BULK_TSV_PATH + ".zip"
    downloaded = 0
    with open(zip_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
            downloaded += len(chunk)
    if verbose:
        print(f"[GWAS] Downloaded {downloaded / 1e6:.1f} MB, extracting...")

    with zipfile.ZipFile(zip_path) as zf:
        tsv_members = [m for m in zf.namelist() if m.endswith(".tsv") or m.endswith(".txt")]
        if not tsv_members:
            raise RuntimeError("No TSV/TXT member found in the GWAS Catalog ZIP.")
        # Pick the member with "associations" in the name, else the first TSV.
        member = next((m for m in tsv_members if "association" in m.lower()), tsv_members[0])
        if verbose:
            print(f"[GWAS] Extracting member: {member}")
        # Extract to a temp file, then replace atomically - an interrupted
        # extraction straight to GWAS_BULK_TSV_PATH would leave a truncated
        # file with a fresh mtime, reused as "fresh" for GWAS_REFRESH_DAYS.
        tmp_path = GWAS_BULK_TSV_PATH + ".tmp"
        with zf.open(member) as src, open(tmp_path, "wb") as dst:
            while True:
                chunk = src.read(1024 * 256)
                if not chunk:
                    break
                dst.write(chunk)
        os.replace(tmp_path, GWAS_BULK_TSV_PATH)
    os.remove(zip_path)
    if verbose:
        print(f"[GWAS] Saved TSV to {GWAS_BULK_TSV_PATH}")
    return GWAS_BULK_TSV_PATH


# Accepted column-name variants for GWAS Catalog (columns renamed across releases).
_GWAS_GENE_COLS = ["MAPPED_GENE", "MAPPED GENE", "Mapped Gene"]
_GWAS_TRAIT_NAME_COLS = ["DISEASE/TRAIT", "DISEASE_TRAIT", "Disease/Trait", "MAPPED_TRAIT", "TRAIT"]
_GWAS_TRAIT_URI_COLS = ["MAPPED_TRAIT_URI", "MAPPED TRAIT URI", "Mapped Trait URI"]


def _pick_col(fieldnames: "list[str]", candidates: "list[str]") -> "str | None":
    for c in candidates:
        if c in fieldnames:
            return c
    # Fuzzy fallback: case-insensitive prefix match
    cl = [c.lower() for c in candidates]
    for f in fieldnames:
        if f.lower() in cl:
            return f
    return None


def build_gwas_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """Build the GWAS Catalog annotated universe for a taxon: trait -> gene-set map.

    For non-human taxa (taxon_id != '9606'), the GWAS Catalog has essentially
    no data; a caveat is set and an empty universe is returned rather than
    raising an error.

    Returns {"traits": {trait_key: {"name": str, "genes": set}},
             "universe_size": int, "fetched_at": iso8601 str, "caveat": str|None}.
    """
    if not force_refresh:
        cached = _load_universe_cache("GWAS Catalog", taxon_id)
        if cached is not None:
            return cached

    caveat = None
    if str(taxon_id) != GWAS_HUMAN_TAXON_ID:
        caveat = (
            "GWAS Catalog coverage is overwhelmingly human. "
            f"Taxon {taxon_id} returned few or no associations; "
            "results below are likely empty."
        )
        if verbose:
            print(f"[GWAS] Non-human taxon {taxon_id}: {caveat}")

    tsv_path = _ensure_gwas_bulk_file(verbose=verbose)

    # Symbol-native: store gene symbols directly (uppercase). No UniProt bridge.
    traits: dict = {}
    annotated_genes: set = set()
    rows_parsed = 0

    with open(tsv_path, "r", encoding="utf-8-sig", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames: "list[str]" = list(reader.fieldnames or [])
        gene_col = _pick_col(fieldnames, _GWAS_GENE_COLS)
        trait_name_col = _pick_col(fieldnames, _GWAS_TRAIT_NAME_COLS)
        trait_uri_col = _pick_col(fieldnames, _GWAS_TRAIT_URI_COLS)

        if verbose:
            print(f"[GWAS] Columns found: gene={gene_col!r}, trait_name={trait_name_col!r}, trait_uri={trait_uri_col!r}")
        if gene_col is None or trait_name_col is None:
            raise RuntimeError(
                f"GWAS Catalog TSV missing expected columns. "
                f"Available: {fieldnames[:20]}"
            )

        for row in reader:
            rows_parsed += 1
            raw_genes = str(row.get(gene_col, "") or "").strip()
            if not raw_genes or raw_genes in ("NR", "N/A", "NA", "-"):
                continue

            trait_name = str(row.get(trait_name_col, "") or "").strip()
            if not trait_name:
                continue

            trait_uri_raw = str(row.get(trait_uri_col, "") or "").strip() if trait_uri_col else ""
            if trait_uri_raw:
                trait_key = trait_uri_raw.split(",")[0].strip().rstrip("/").rsplit("/", 1)[-1]
            else:
                trait_key = re.sub(r"[^A-Za-z0-9_-]", "_", trait_name)[:80]

            raw_genes = raw_genes.replace(" - ", ",").replace(";", ",")
            # Normalize to uppercase for case-insensitive matching with the input list
            symbols = [g.strip().upper() for g in raw_genes.split(",") if g.strip()]

            bucket = traits.setdefault(trait_key, {"name": trait_name, "genes": set()})
            for sym in symbols:
                bucket["genes"].add(sym)
                annotated_genes.add(sym)

    if verbose:
        print(f"[GWAS] Parsed {rows_parsed} rows -> {len(traits)} traits, {len(annotated_genes)} unique gene symbols")

    result = {
        "traits": traits,
        "universe_size": len(annotated_genes),
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "caveat": caveat,
    }
    _save_universe_cache("GWAS Catalog", taxon_id, result)
    return result


# ---------------------------------------------------------------------------
# ORA universe cache — one JSON artifact per (database, taxon_id), plus a
# lightweight CSV index so staleness/existence can be checked without
# opening every blob. GO/Reactome cache indefinitely (force_refresh=True to
# bust); GWAS Catalog's bulk file uses its own TTL (see GWAS section).
# ---------------------------------------------------------------------------

def _ensure_universe_cache_dir() -> None:
    os.makedirs(ORA_UNIVERSE_CACHE_DIR, exist_ok=True)


def _universe_cache_file_path(database: str, taxon_id: str) -> str:
    safe_db = re.sub(r"[^A-Za-z0-9_-]", "_", database)
    return os.path.join(ORA_UNIVERSE_CACHE_DIR, f"{safe_db}_{taxon_id}.json")


def _universe_to_json_safe(data: dict) -> dict:
    """Convert internal set-valued universe dicts to JSON-serializable lists."""
    def _convert(obj: object) -> object:
        if isinstance(obj, set):
            return sorted(obj)  # type: ignore[type-var]
        if isinstance(obj, dict):
            return {k: _convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_convert(v) for v in obj]
        return obj
    return _convert(data)  # type: ignore[return-value]


def _universe_from_json_safe(data: dict) -> dict:
    """Reverse of _universe_to_json_safe: turn "genes" lists back into sets."""
    def _convert(obj: object, key: "str | None" = None) -> object:
        if isinstance(obj, dict):
            return {k: _convert(v, key=k) for k, v in obj.items()}
        if isinstance(obj, list) and key == "genes":
            return set(obj)
        return obj
    return _convert(data)  # type: ignore[return-value]


def _load_universe_cache(database: str, taxon_id: str) -> "dict | None":
    path = _universe_cache_file_path(database, taxon_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)
        return _universe_from_json_safe(raw)
    except Exception:
        return None


def _save_universe_cache(database: str, taxon_id: str, data: dict) -> None:
    _ensure_universe_cache_dir()
    path = _universe_cache_file_path(database, taxon_id)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(_universe_to_json_safe(data), fh)
    except Exception:
        return
    _update_universe_index(database, taxon_id, path, data.get("universe_size", 0), data.get("fetched_at", ""))


def _update_universe_index(database: str, taxon_id: str, path: str, universe_size: int, fetched_at: str) -> None:
    _ensure_universe_cache_dir()
    rows = []
    if os.path.exists(ORA_UNIVERSE_INDEX_PATH):
        try:
            existing = pd.read_csv(ORA_UNIVERSE_INDEX_PATH, dtype=str).fillna("")
            rows = [
                row for _, row in existing.iterrows()
                if not (row.get("database") == database and row.get("taxon_id") == str(taxon_id))
            ]
            rows = [dict(row) for row in rows]
        except Exception:
            rows = []
    rows.append({
        "database": database, "taxon_id": str(taxon_id), "file_path": path,
        "universe_size": str(universe_size), "fetched_at": fetched_at,
    })
    pd.DataFrame(rows).to_csv(ORA_UNIVERSE_INDEX_PATH, index=False)


# ---------------------------------------------------------------------------
# ORA orchestrator
# ---------------------------------------------------------------------------

def _run_fisher_for_namespace(
    term_to_genes: dict,
    list_set: set,
    universe_genes: set,
    fdr_alpha: float = FDR_ALPHA_DEFAULT,
    include_bonferroni: bool = True,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
) -> dict:
    """Shared per-test-family Fisher's exact + BH-FDR runner.

    term_to_genes: {term_id: {"name": str, "genes": set}} restricted to terms
                   already known to belong to this namespace/database.
    list_set:      the user's input accessions, NOT YET restricted to the universe.
    universe_genes: the full annotated universe for this namespace/database
                    (already intersected with any user-uploaded background).

    Returns a result_block dict (see run_ora_for_database docstring).
    """
    n_input_in_universe_set = list_set & universe_genes
    n_input_dropped = len(list_set) - len(n_input_in_universe_set)
    n = len(n_input_in_universe_set)
    N = len(universe_genes)

    rows = []
    # Counts terms excluded by the size gate (too small OR, if set, too large) so
    # they never enter the BH-FDR family. Field name kept for API/UI compatibility.
    n_excluded_min_size = 0
    for term_id, info in term_to_genes.items():
        term_genes = info["genes"] & universe_genes
        K = len(term_genes)
        if K < min_term_size or (max_term_size and K > max_term_size):
            n_excluded_min_size += 1
            continue
        overlap = term_genes & n_input_in_universe_set
        k = len(overlap)
        # Short-circuit the k==0 case: a one-sided ("greater") Fisher test with
        # zero overlap is exactly p=1.0 (P(X>=0)=1) and fold enrichment 0, so we
        # skip the scipy call. This is a big speedup when many terms don't
        # overlap the input at all (e.g. per-module enrichment of small modules
        # against the whole GO ontology), and changes no result.
        if n <= 0:
            p_val, fold = float("nan"), float("nan")
        elif k == 0:
            p_val, fold = 1.0, 0.0
        else:
            p_val, fold = fisher_exact_pvalue(k, n, K, N), fold_enrichment(k, n, K, N)
        rows.append({
            "Term ID": term_id,
            "Term Name": info["name"],
            "k": k, "K": K, "n": n, "N": N,
            "Fold Enrichment": fold,
            "p-value": p_val,
            "Overlapping Genes": ", ".join(sorted(overlap)),
        })

    if rows:
        pvalues = [r["p-value"] for r in rows]
        fdrs = benjamini_hochberg(pvalues)
        for row, fdr in zip(rows, fdrs):
            row["FDR"] = fdr
        if include_bonferroni:
            bonf = bonferroni(pvalues)
            for row, b in zip(rows, bonf):
                row["Bonferroni"] = b
        else:
            for row in rows:
                row["Bonferroni"] = float("nan")

    table = pd.DataFrame(rows, columns=ORA_RESULT_COLUMNS)
    if not table.empty:
        table = table.sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)

    return {
        "table": table,
        "n_terms_tested": len(rows),
        "n_terms_excluded_min_size": n_excluded_min_size,
        "universe_size": N,
        "n_input_in_universe": n,
        "n_input_dropped": n_input_dropped,
        "caveat": None,
    }


def run_ora_for_database(
    database: str,
    taxon_id: str,
    list_accessions: list,
    list_symbols: "list | None" = None,
    fallback_symbols: "list | None" = None,
    background_accessions: "set | None" = None,
    fdr_alpha: float = FDR_ALPHA_DEFAULT,
    include_bonferroni: bool = True,
    include_iea: bool = True,
    taxon_name: str = "",
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run ORA for one database against the given gene identifier list.

    Args:
        database:            "GO" | "Reactome" | "GWAS Catalog"
        taxon_id:            NCBI taxon ID for the species choice.
        list_accessions:     flat list of UniProt accessions (used by GO and Reactome).
        list_symbols:        flat list of gene symbols (used by GWAS Catalog — matched
                             in uppercase against MAPPED_GENE). None or empty => GWAS
                             returns an empty result with a caveat.
        fallback_symbols:    gene symbols of input genes that had NO UniProt ID. GO only:
                             recovered to accessions via QuickGO's SYMBOL column so they
                             can still be tested (reported as n_recovered_by_symbol).
        background_accessions: optional background set (accession space; not used for GWAS).
        include_iea:         when False, exclude electronically-inferred (IEA) annotations
                             from GO and Reactome universes.
        taxon_name:          organism name string (used by Reactome to filter its file).

    Returns:
        For database="GO": {"biological_process": <block>, "molecular_function": <block>,
                             "cellular_component": <block>}
        For Reactome/GWAS Catalog: {"<database_key>": <block>}
    """
    list_set = {str(a).strip() for a in list_accessions if str(a).strip()}
    # Uppercase symbol set for GWAS (symbol-native universe)
    symbol_set = {str(s).strip().upper() for s in (list_symbols or []) if str(s).strip()}

    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    if database == "GO":
        _progress(0, 1, "Building GO universe...")
        universe = build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose)

        # Gene-symbol fallback: recover input genes that had no UniProt ID by
        # mapping their symbol -> accession via QuickGO's own SYMBOL column.
        symbol_index = universe.get("symbol_index", {}) or {}
        recovered_accs: set = set()
        n_recovered_by_symbol = 0
        for sym in (fallback_symbols or []):
            acc = symbol_index.get(str(sym).strip().upper())
            if acc and acc not in list_set:
                recovered_accs.add(acc)
                n_recovered_by_symbol += 1
        effective_set = list_set | recovered_accs

        results = {}
        for i, ns in enumerate(GO_NAMESPACES):
            _progress(i, len(GO_NAMESPACES), f"Testing GO:{ns}")
            universe_genes: set = set()
            for info in universe[ns].values():
                universe_genes |= info["genes"]
            if background_accessions:
                universe_genes &= background_accessions
            block = _run_fisher_for_namespace(
                universe[ns], effective_set, universe_genes, fdr_alpha, include_bonferroni,
                min_term_size, max_term_size,
            )
            block["n_recovered_by_symbol"] = n_recovered_by_symbol
            if universe.get("caveat"):
                block["caveat"] = universe["caveat"]
            results[ns] = block
        _progress(len(GO_NAMESPACES), len(GO_NAMESPACES), "Done")
        return results

    if database == "Reactome":
        _progress(0, 1, "Building Reactome universe...")
        universe = build_reactome_universe(
            taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose,
        )
        universe_genes = set()
        for info in universe["pathways"].values():
            universe_genes |= info["genes"]
        if background_accessions:
            universe_genes &= background_accessions
        _progress(0, 1, "Testing Reactome pathways")
        result = _run_fisher_for_namespace(
            universe["pathways"], list_set, universe_genes, fdr_alpha, include_bonferroni,
            min_term_size, max_term_size,
        )
        _progress(1, 1, "Done")
        return {"Reactome": result}

    if database == "GWAS Catalog":
        # GWAS is symbol-native: test in symbol space, not accession space.
        if not symbol_set:
            empty_table = pd.DataFrame(columns=ORA_RESULT_COLUMNS)
            return {"GWAS Catalog": {
                "table": empty_table,
                "n_terms_tested": 0,
                "n_terms_excluded_min_size": 0,
                "universe_size": 0,
                "n_input_in_universe": 0,
                "n_input_dropped": 0,
                "caveat": (
                    "GWAS Catalog requires Gene Symbol input. "
                    "No gene symbols were available for this run."
                ),
            }}
        _progress(0, 1, "Building GWAS Catalog universe...")
        universe = build_gwas_universe(taxon_id, verbose=verbose)
        universe_genes = set()
        for info in universe["traits"].values():
            universe_genes |= info["genes"]
        _progress(0, 1, "Testing GWAS Catalog traits")
        result = _run_fisher_for_namespace(
            universe["traits"], symbol_set, universe_genes, fdr_alpha, include_bonferroni,
            min_term_size, max_term_size,
        )
        result["caveat"] = universe.get("caveat")
        _progress(1, 1, "Done")
        return {"GWAS Catalog": result}

    raise ValueError(f"Unknown ORA database: {database!r}")


def run_ora_for_custom_db(
    term_map: dict,
    id_type: str,
    list_accessions: list,
    list_symbols: "list | None" = None,
    list_ensembl: "list | None" = None,
    db_label: str = "Custom database",
    fdr_alpha: float = FDR_ALPHA_DEFAULT,
    include_bonferroni: bool = True,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run ORA against a user-uploaded gene-set database (see load_custom_geneset_db).

    Reuses the shared Fisher/BH engine. Input identifiers are chosen by id_type:
    accessions for 'accession', uppercased symbols for 'symbol', uppercased
    Ensembl Gene IDs for 'ensembl'. Universe = union of all genes across the
    uploaded sets. The input identifier list must be in the SAME space as
    id_type (e.g. an 'ensembl' custom DB only matches when the resolved input
    table has an Ensembl Gene ID column) — a mismatch surfaces as a caveat.
    """
    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    id_label = {"symbol": "gene symbols", "ensembl": "Ensembl Gene IDs"}.get(id_type, "UniProt accessions")
    if id_type == "symbol":
        input_set = {str(s).strip().upper() for s in (list_symbols or []) if str(s).strip()}
    elif id_type == "ensembl":
        input_set = {str(e).strip().upper() for e in (list_ensembl or []) if str(e).strip()}
    else:
        input_set = {str(a).strip() for a in list_accessions if str(a).strip()}

    universe_genes: set = set()
    for info in term_map.values():
        universe_genes |= info["genes"]

    _progress(0, 1, "Testing custom gene sets")
    block = _run_fisher_for_namespace(
        term_map, input_set, universe_genes, fdr_alpha, include_bonferroni,
        min_term_size, max_term_size,
    )
    if block["n_input_in_universe"] == 0:
        block["caveat"] = (
            f"0 of your {len(input_set)} gene(s) matched this database's identifiers. "
            f"Check that the file uses {id_label} matching your input."
        )
    _progress(1, 1, "Done")
    return {db_label: block}


# ---------------------------------------------------------------------------
# Gene-set Burden Test — competitive Wilcoxon rank-sum, pathway-level analog
# of a per-gene variant burden test. Where ORA tests binary membership with
# Fisher's exact test, this tests whether a term's genes carry systematically
# higher continuous scores (variant severity, GWAS association strength, a
# weighted variant count, ...) than the rest of the scored genes - the same
# competitive framing used by GSEA/CAMERA/MAGMA-style gene-set tests, here
# via a nonparametric rank test rather than a permutation-based statistic.
# Reuses the same universe builders, BH-FDR/Bonferroni helpers, and term-size
# gate as the ORA engine above.
# ---------------------------------------------------------------------------

def _run_wilcoxon_for_namespace(
    term_to_genes: dict,
    score_by_gene: dict,
    universe_genes: set,
    include_bonferroni: bool = True,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
) -> dict:
    """Shared per-test-family competitive Wilcoxon rank-sum + BH-FDR runner.

    term_to_genes: {term_id: {"name": str, "genes": set}} restricted to terms
                   already known to belong to this namespace/database.
    score_by_gene: {gene_id: float} variant/severity score per gene, already
                   oriented so that HIGHER = more burden (the caller negates
                   scores upstream when the user's metric runs the other way,
                   e.g. an untransformed p-value where lower = more signal).
    universe_genes: the full annotated universe for this namespace/database.

    For each term, compares the scores of genes inside the term ("in-set")
    against every other scored gene in the universe ("out-set") with a
    one-sided Mann-Whitney U / Wilcoxon rank-sum test (alternative="greater"
    -> in-set scores are stochastically larger). Terms with fewer than
    min_term_size scored genes, or more than max_term_size, are excluded
    before FDR correction (same size-gate semantics as the ORA engine); a
    term also needs at least one "out-set" gene to compare against.

    Returns a result_block dict shaped like run_ora_for_database's blocks,
    but with the burden-specific columns in BURDEN_RESULT_COLUMNS.
    """
    if _scipy_mannwhitneyu is None:
        raise RuntimeError("scipy is required for the gene-set burden test")

    scored_universe = {g for g in universe_genes if g in score_by_gene}
    n_input_dropped = len(score_by_gene) - len(scored_universe)
    n = len(scored_universe)
    N = len(universe_genes)

    rows = []
    # Counts terms excluded by the size gate (too small OR too large) or for
    # lacking a comparison set, so they never enter the BH-FDR family.
    n_excluded_min_size = 0
    for term_id, info in term_to_genes.items():
        term_genes = info["genes"] & scored_universe
        K = len(term_genes)
        if K < min_term_size or (max_term_size and K > max_term_size):
            n_excluded_min_size += 1
            continue
        out_genes = scored_universe - term_genes
        if not out_genes:
            n_excluded_min_size += 1
            continue
        in_scores = [score_by_gene[g] for g in term_genes]
        out_scores = [score_by_gene[g] for g in out_genes]
        try:
            _stat, p_val = _scipy_mannwhitneyu(in_scores, out_scores, alternative="greater")
            p_val = float(p_val)
        except ValueError:
            # e.g. every score identical (mannwhitneyu can raise on degenerate input)
            p_val = 1.0
        mean_in = sum(in_scores) / len(in_scores)
        mean_out = sum(out_scores) / len(out_scores)
        rows.append({
            "Term ID": term_id,
            "Term Name": info["name"],
            "Genes in set": K,
            "Scored genes (universe)": n,
            "Mean score (set)": mean_in,
            "Mean score (rest)": mean_out,
            "Effect (mean difference)": mean_in - mean_out,
            "p-value": p_val,
            "Genes": ", ".join(sorted(term_genes)),
        })

    if rows:
        pvalues = [r["p-value"] for r in rows]
        fdrs = benjamini_hochberg(pvalues)
        for row, fdr in zip(rows, fdrs):
            row["FDR"] = fdr
        if include_bonferroni:
            bonf = bonferroni(pvalues)
            for row, b in zip(rows, bonf):
                row["Bonferroni"] = b
        else:
            for row in rows:
                row["Bonferroni"] = float("nan")

    table = pd.DataFrame(rows, columns=BURDEN_RESULT_COLUMNS)
    if not table.empty:
        table = table.sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)

    return {
        "table": table,
        "n_terms_tested": len(rows),
        "n_terms_excluded_min_size": n_excluded_min_size,
        "universe_size": N,
        "n_input_in_universe": n,
        "n_input_dropped": n_input_dropped,
        "caveat": None,
    }


def run_burden_for_database(
    database: str,
    taxon_id: str,
    score_by_acc: dict,
    include_bonferroni: bool = True,
    include_iea: bool = True,
    taxon_name: str = "",
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run the competitive Wilcoxon gene-set burden test for one database.

    Args:
        database:      "GO" | "Reactome"
        taxon_id:      NCBI taxon ID for the species choice.
        score_by_acc:  {UniProt accession: float} variant/severity score per
                       gene, already oriented so HIGHER = more burden (the
                       caller negates scores upstream when the user's raw
                       metric runs the other way, e.g. an untransformed
                       p-value).
        include_iea:   when False, exclude electronically-inferred (IEA)
                       annotations from GO and Reactome universes.
        taxon_name:    organism name string (used by Reactome to filter its file).

    Returns:
        For database="GO": {"biological_process": <block>, "molecular_function": <block>,
                             "cellular_component": <block>}
        For Reactome: {"Reactome": <block>}
    """
    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    if database == "GO":
        _progress(0, 1, "Building GO universe...")
        universe = build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose)
        results = {}
        for i, ns in enumerate(GO_NAMESPACES):
            _progress(i, len(GO_NAMESPACES), f"Testing GO:{ns}")
            universe_genes: set = set()
            for info in universe[ns].values():
                universe_genes |= info["genes"]
            block = _run_wilcoxon_for_namespace(
                universe[ns], score_by_acc, universe_genes, include_bonferroni,
                min_term_size, max_term_size,
            )
            if universe.get("caveat"):
                block["caveat"] = universe["caveat"]
            results[ns] = block
        _progress(len(GO_NAMESPACES), len(GO_NAMESPACES), "Done")
        return results

    if database == "Reactome":
        _progress(0, 1, "Building Reactome universe...")
        universe = build_reactome_universe(
            taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose,
        )
        universe_genes = set()
        for info in universe["pathways"].values():
            universe_genes |= info["genes"]
        _progress(0, 1, "Testing Reactome pathways")
        result = _run_wilcoxon_for_namespace(
            universe["pathways"], score_by_acc, universe_genes, include_bonferroni,
            min_term_size, max_term_size,
        )
        _progress(1, 1, "Done")
        return {"Reactome": result}

    raise ValueError(f"Unknown burden-test database: {database!r}")


def run_burden_for_custom_db(
    term_map: dict,
    score_by_id: dict,
    db_label: str = "Custom database",
    include_bonferroni: bool = True,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size: "int | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run the burden test against a user-uploaded gene-set database.

    score_by_id must already be keyed in the SAME identifier space as the
    custom database's id_type (accession/symbol/ensembl) - the caller builds
    that mapping (mirrors run_ora_for_custom_db's id_type handling). Universe
    = union of all genes across the uploaded sets.
    """
    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    universe_genes: set = set()
    for info in term_map.values():
        universe_genes |= info["genes"]

    _progress(0, 1, "Testing custom gene sets")
    block = _run_wilcoxon_for_namespace(
        term_map, score_by_id, universe_genes, include_bonferroni,
        min_term_size, max_term_size,
    )
    if block["n_input_in_universe"] == 0:
        block["caveat"] = (
            f"0 of your {len(score_by_id)} scored gene(s) matched this database's identifiers."
        )
    _progress(1, 1, "Done")
    return {db_label: block}
