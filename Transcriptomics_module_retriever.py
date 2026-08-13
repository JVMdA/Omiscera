# -*- coding: utf-8 -*-
"""
Transcriptomics module retriever for Omiscera.

Turns an uploaded differential-expression (DE) results table (gene identifier +
log2 fold-change + p-value[+ adjusted p-value]) into one clean, analysis-ready
table that every Transcriptomics analysis tile (ORA, GSEA, Enrichment,
TF/Regulon) consumes downstream. This file is standalone: it does not import,
reuse, or depend on any Genomics_* module.

Accepted identifier columns (auto-detected): Gene Symbol, Ensembl Gene ID,
Ensembl Transcript ID, Entrez Gene ID. Species-general: the caller supplies a
taxon ID and identifier mapping is resolved online for any organism, mirroring
the Genomics retriever's organism selector but without any bundled per-species
mapping file.

Identifier resolution uses mygene.info's bulk /v3/query endpoint (any species,
any of the four identifier types in one request, alias/retired-symbol aware),
with the Ensembl REST API (/lookup/id) as a fallback for transcripts mygene
does not resolve directly.

The retriever never fetches from GO/Reactome/MSigDB/TRRUST/DoRothEA/ChEA3 —
those are analysis-engine concerns. It also never thresholds the DE
table or computes a ranking metric: it returns the full, unthresholded table,
and each analysis tile takes its own view of it later.
"""

import io
import os
import re
import time

import numpy as np
import pandas as pd

try:
    import requests as _requests
    HTTP_SESSION = _requests.Session()
except ImportError:
    _requests = None
    HTTP_SESSION = None

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MYGENE_QUERY_URL = "https://mygene.info/v3/query"
MYGENE_CHUNK = 1000  # mygene.info batch query size per request

ENSEMBL_REST_URL = "https://rest.ensembl.org"

OMISCERA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositorio do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "application/json",
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
TRANSCRIPTOMICS_IDMAP_CACHE_PATH = os.path.join(CACHE_DIR, "transcriptomics_idmap_cache.csv")

DEFAULT_TIMEOUT = 20
DEFAULT_RETRIES = 3

RESOLVED_COLUMNS = [
    "input_id", "gene_symbol", "entrez_id", "ensembl_id",
    "log2FC", "pvalue", "padj", "status",
]

# identifier -> mygene.info "scopes" search field(s). "retired" (deprecated/
# withdrawn HGNC symbols) is included for gene_symbol, not just mixed - an
# exact lookup in a curated field, not fuzzy matching; measured to recover a
# meaningful fraction of real-world gene lists that carry older nomenclature.
MYGENE_SCOPES = {
    "gene_symbol": "symbol,alias,retired",
    "ensembl_gene": "ensembl.gene",
    "ensembl_transcript": "ensembl.transcript",
    "entrez": "entrezgene",
    "mixed": "symbol,alias,ensembl.gene,ensembl.transcript,entrezgene,retired",
}

# ---------------------------------------------------------------------------
# Column detection (arbitrary uploaded headers -> canonical fields)
# ---------------------------------------------------------------------------

COLUMN_SYNONYMS = {
    "identifier": {
        "gene", "geneid", "gene_id", "genesymbol", "gene_symbol", "symbol",
        "genename", "gene_name", "ensembl", "ensemblid", "ensembl_id",
        "ensemblgeneid", "ensembl_gene_id", "ensembltranscriptid",
        "ensembl_transcript_id", "transcript", "transcriptid",
        "transcript_id", "targetid", "target_id", "feature", "featureid",
        "entrez", "entrezid", "entrez_id", "entrezgeneid", "entrez_gene_id",
        "ncbigeneid", "ncbi_gene_id",
    },
    "log2FC": {
        "log2foldchange", "log2fc", "logfc", "lfc", "log2_fold_change",
        "log2fold_change", "coef", "beta",
    },
    "pvalue": {
        "pvalue", "pval", "p_value", "p", "pvalraw",
    },
    "padj": {
        "padj", "fdr", "qvalue", "q_value", "adjpvalue", "adj_p_val",
        "adjpval", "p_adj", "bhfdr", "fdrq",
    },
}


def _normalize_col_key(col: str) -> str:
    """Flatten a column name to a lowercase identifier for comparison."""
    return (
        str(col)
        .replace("﻿", "")
        .replace("\xa0", " ")
        .strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
        .replace(".", "")
    )


def detect_columns(df: pd.DataFrame) -> dict:
    """Map arbitrary uploaded headers to canonical identifier/log2FC/pvalue/padj.

    Returns a dict {canonical_field: original_column_name}. Raises ValueError
    if identifier, log2FC, or pvalue cannot be found (padj is optional — the
    Benjamini-Hochberg FDR is computed from pvalue when absent).
    """
    normalized_synonyms = {
        target: {_normalize_col_key(s) for s in syns}
        for target, syns in COLUMN_SYNONYMS.items()
    }
    mapping = {}
    for col in df.columns:
        key = _normalize_col_key(col)
        for target, syns in normalized_synonyms.items():
            if target in mapping:
                continue
            if key in syns:
                mapping[target] = col
                break

    missing = [t for t in ("identifier", "log2FC", "pvalue") if t not in mapping]
    if missing:
        raise ValueError(
            f"DE table is missing required column(s) for: {missing}. "
            f"Detected columns after normalisation: {list(df.columns)}"
        )
    return mapping


# ---------------------------------------------------------------------------
# Identifier type detection
# ---------------------------------------------------------------------------

_ENSG_RE = re.compile(r"^ENS[A-Z]*G\d+(\.\d+)?$", re.IGNORECASE)
_ENST_RE = re.compile(r"^ENS[A-Z]*T\d+(\.\d+)?$", re.IGNORECASE)
_ENTREZ_RE = re.compile(r"^\d+$")


def _strip_version(identifier: str) -> str:
    """Strip an Ensembl version suffix (e.g. 'ENSG00000141510.15' -> '...510')."""
    return re.sub(r"\.\d+$", "", identifier.strip())


def detect_identifier_type(series: pd.Series) -> str:
    """Classify a column of identifiers as one type, or 'mixed' if no clear
    majority (>=90%) agrees."""
    values = series.dropna().astype(str).str.strip()
    values = values[values != ""]
    if values.empty:
        return "gene_symbol"

    counts = {"ensembl_gene": 0, "ensembl_transcript": 0, "entrez": 0, "gene_symbol": 0}
    for val in values:
        if _ENSG_RE.match(val):
            counts["ensembl_gene"] += 1
        elif _ENST_RE.match(val):
            counts["ensembl_transcript"] += 1
        elif _ENTREZ_RE.match(val):
            counts["entrez"] += 1
        else:
            counts["gene_symbol"] += 1

    total = sum(counts.values())
    best_type, best_count = max(counts.items(), key=lambda kv: kv[1])
    if total and (best_count / total) >= 0.9:
        return best_type
    return "mixed"


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def _parse_de_csv_text(text: str) -> pd.DataFrame:
    """Parse CSV text with a multi-delimiter sniff, headers untouched."""
    best_df = None
    for sep in (",", "\t", ";"):
        try:
            candidate = pd.read_csv(io.StringIO(text), sep=sep, dtype=str, keep_default_na=False)
            if len(candidate.columns) < 1:
                continue
            if best_df is None or len(candidate.columns) > len(best_df.columns):
                best_df = candidate
        except Exception:
            continue
    if best_df is None:
        raise ValueError("Could not parse the uploaded CSV file.")
    return best_df


def load_de_table(path: str) -> pd.DataFrame:
    """Load a DE results table from a file path. Headers are left as-is;
    column detection happens in resolve_de_table."""
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin1") as fh:
            text = fh.read()
    return _parse_de_csv_text(text)


def load_de_table_from_bytes(content: bytes) -> pd.DataFrame:
    """Load a DE results table from raw bytes (e.g. an uploaded file)."""
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("utf-8", errors="replace")
    return _parse_de_csv_text(text)


# ---------------------------------------------------------------------------
# Numeric coercion + Benjamini-Hochberg FDR
# ---------------------------------------------------------------------------

def _coerce_numeric(df: pd.DataFrame) -> "tuple[pd.DataFrame, int]":
    """Coerce log2FC/pvalue to numeric, dropping rows where either is unusable.

    Returns (clean_df, n_dropped). padj (if present) is coerced but rows are
    never dropped for a bad padj — it's optional and gets recomputed if absent.
    """
    df = df.copy()
    df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    valid = df["log2FC"].notna() & df["pvalue"].notna()
    n_dropped = int((~valid).sum())
    df = df[valid].reset_index(drop=True)
    if "padj" in df.columns:
        df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    return df, n_dropped


def _benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Standard BH step-up FDR, returned aligned to the input Series' index."""
    n = len(pvalues)
    if n == 0:
        return pvalues.copy()
    sorted_idx = pvalues.sort_values(ascending=True).index
    sorted_p = pvalues.loc[sorted_idx].to_numpy(dtype=float)
    ranks = np.arange(1, n + 1)
    raw_padj = sorted_p * n / ranks
    # enforce monotone non-decreasing from the largest p-value backward
    adjusted = np.minimum.accumulate(raw_padj[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(adjusted, index=sorted_idx).reindex(pvalues.index)


# ---------------------------------------------------------------------------
# mygene.info / Ensembl REST HTTP helpers
# ---------------------------------------------------------------------------

def _mygene_post(ids: list, scopes: str, species: str, verbose: bool = False) -> "list | None":
    """POST a batch identifier query to mygene.info.

    Returns the raw JSON list of hit dicts on success (each carrying 'query',
    and either match fields or 'notfound') - possibly an empty list, which is
    a valid successful response. Returns None if the request could not be
    completed after retries (network error, non-2xx status, unparseable
    response, or the requests library being unavailable): callers MUST treat
    None as "unknown" and never conflate it with "no genes found", since doing
    so would silently mark a whole batch of real genes as unmapped on a
    transient failure.
    """
    if not ids:
        return []
    if HTTP_SESSION is None:
        if verbose:
            print("mygene.info query skipped: requests library not available")
        return None
    data = {
        "q": ",".join(ids),
        "scopes": scopes,
        "species": species,
        "fields": "symbol,entrezgene,ensembl.gene",
    }
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            response = HTTP_SESSION.post(
                MYGENE_QUERY_URL, data=data, headers=OMISCERA_HEADERS, timeout=DEFAULT_TIMEOUT,
            )
            if response.status_code in (200, 201):
                try:
                    result = response.json()
                    return result if isinstance(result, list) else []
                except Exception as exc:
                    if verbose:
                        print(f"mygene.info JSON parse failed | error={exc}")
                    return None
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
        print(f"mygene.info query failed | error={last_error}")
    return None


def _parse_mygene_hit(hit: dict) -> dict:
    symbol = str(hit.get("symbol", "") or "").strip()
    entrez_raw = hit.get("entrezgene", "")
    entrez = str(entrez_raw).strip() if entrez_raw not in (None, "") else ""
    ensembl_field = hit.get("ensembl")
    ensembl_gene = ""
    if isinstance(ensembl_field, dict):
        ensembl_gene = str(ensembl_field.get("gene", "") or "").strip()
    elif isinstance(ensembl_field, list) and ensembl_field:
        first = ensembl_field[0]
        if isinstance(first, dict):
            ensembl_gene = str(first.get("gene", "") or "").strip()
    return {"symbol": symbol, "entrez": entrez, "ensembl": ensembl_gene}


def _ensembl_parent_gene(transcript_id: str, verbose: bool = False) -> str:
    """Fallback ENST -> parent ENSG lookup via the Ensembl REST API."""
    if HTTP_SESSION is None:
        return ""
    url = f"{ENSEMBL_REST_URL}/lookup/id/{transcript_id}"
    try:
        response = HTTP_SESSION.get(
            url, params={"content-type": "application/json"},
            headers={"Accept": "application/json"}, timeout=DEFAULT_TIMEOUT,
        )
        if response.status_code in (200, 201):
            data = response.json()
            return str(data.get("Parent", "") or "").strip()
    except Exception as exc:
        if verbose:
            print(f"Ensembl parent lookup failed | id={transcript_id} | error={exc}")
    return ""


# ---------------------------------------------------------------------------
# On-disk identifier-mapping cache (per taxon + identifier type + id)
# ---------------------------------------------------------------------------

TRANSCRIPTOMICS_IDMAP_CACHE: dict = {}
_IDMAP_CACHE_LOADED = False


def _ensure_cache_loaded() -> None:
    global _IDMAP_CACHE_LOADED
    if _IDMAP_CACHE_LOADED:
        return
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(TRANSCRIPTOMICS_IDMAP_CACHE_PATH):
        try:
            df = pd.read_csv(TRANSCRIPTOMICS_IDMAP_CACHE_PATH, dtype=str).fillna("")
            for _, row in df.iterrows():
                key = row.get("CacheKey", "")
                if key:
                    TRANSCRIPTOMICS_IDMAP_CACHE[key] = {
                        "symbol": row.get("Symbol", ""),
                        "entrez": row.get("Entrez", ""),
                        "ensembl": row.get("Ensembl", ""),
                        "status": row.get("Status", ""),
                    }
        except Exception:
            pass
    _IDMAP_CACHE_LOADED = True


def save_idmap_cache() -> None:
    """Persist the in-memory identifier-mapping cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for key, data in TRANSCRIPTOMICS_IDMAP_CACHE.items():
        rows.append({
            "CacheKey": key,
            "Symbol": data.get("symbol", ""),
            "Entrez": data.get("entrez", ""),
            "Ensembl": data.get("ensembl", ""),
            "Status": data.get("status", ""),
        })
    pd.DataFrame(rows).to_csv(TRANSCRIPTOMICS_IDMAP_CACHE_PATH, index=False)


def _get_cache(key: str):
    _ensure_cache_loaded()
    return TRANSCRIPTOMICS_IDMAP_CACHE.get(key)


def _set_cache(key: str, data: dict) -> None:
    _ensure_cache_loaded()
    # Only cache confident positive matches - caching "unmapped" would poison
    # the cache if a query transiently failed rather than genuinely not existing.
    if data.get("status") in ("mapped", "alias_mapped"):
        TRANSCRIPTOMICS_IDMAP_CACHE[key] = {
            "symbol": data.get("symbol", ""),
            "entrez": data.get("entrez", ""),
            "ensembl": data.get("ensembl", ""),
            "status": data.get("status", ""),
        }


# ---------------------------------------------------------------------------
# Bulk identifier mapping
# ---------------------------------------------------------------------------

def _result_for_query(id_type: str, q: str, parsed: dict) -> dict:
    """Build a {"symbol","entrez","ensembl","status"} result dict for a single
    resolved query, applying the same alias-detection rule everywhere it's
    needed (primary pass, transcript fallback, mixed-scope fallback)."""
    status = "mapped"
    if id_type == "gene_symbol" and parsed["symbol"].lower() != q.lower():
        status = "alias_mapped"
    ensembl_val = parsed["ensembl"] or (q if id_type == "ensembl_gene" else "")
    entrez_val = parsed["entrez"] or (q if id_type == "entrez" else "")
    return {"symbol": parsed["symbol"], "entrez": entrez_val, "ensembl": ensembl_val, "status": status}


def map_identifiers(
    df: pd.DataFrame,
    id_type: str,
    taxon_id: str,
    progress_callback=None,
    verbose: bool = False,
) -> "tuple[pd.DataFrame, dict]":
    """Map every value in df['identifier'] to a gene symbol + Entrez + Ensembl
    gene ID, in bulk (never per-gene), via mygene.info with an Ensembl REST
    fallback for unresolved transcripts and a broadened-scope second pass for
    anything still unresolved (recovers e.g. retired/deprecated symbols).

    Adds columns: gene_symbol, entrez_id, ensembl_id, status
    (status: mapped | alias_mapped | unmapped).

    Returns (df_with_columns, meta) where meta has:
      n_recovered_fallback  - resolved only by the second-pass broadened scope
      n_ambiguous           - queries that matched more than one gene (the
                               first hit mygene returned is kept; deterministic
                               but disclosed, never resolved interactively)
      n_fetch_failures      - mygene.info batch requests that failed outright
                               after retries (network/HTTP error), as opposed
                               to a successful response saying "not found"
      n_ids_affected_by_failure - identifiers still unresolved specifically
                               because their batch failed to fetch, not
                               because mygene confirmed no match

    Raises RuntimeError if EVERY primary-pass batch failed to fetch, rather
    than silently returning a table where everything is marked "unmapped" for
    a reason that has nothing to do with the genes themselves.
    """
    scopes = MYGENE_SCOPES.get(id_type, MYGENE_SCOPES["mixed"])
    strip_versions = id_type in ("ensembl_gene", "ensembl_transcript", "mixed")

    def _query_id(raw: str) -> str:
        raw = str(raw).strip()
        return _strip_version(raw) if strip_versions else raw

    id_to_query = {}
    unique_ids = []
    seen = set()
    for raw in df["identifier"]:
        q = _query_id(raw)
        id_to_query[raw] = q
        qk = q.lower()
        if qk and qk not in seen:
            seen.add(qk)
            unique_ids.append(q)

    results = {}  # query_lower -> {"symbol","entrez","ensembl","status"}; absent/None = unresolved
    to_fetch = []
    for q in unique_ids:
        cache_key = f"{taxon_id}|{id_type}|{q.lower()}"
        cached = _get_cache(cache_key)
        if cached:
            results[q.lower()] = cached
        else:
            to_fetch.append(q)

    n_fetch_failures = 0
    n_primary_chunks = 0
    n_primary_chunk_failures = 0
    failed_query_keys: set = set()
    hit_counts: dict = {}  # query_lower -> number of non-notfound hits seen (ambiguity)

    fetched_count = 0
    for i in range(0, len(to_fetch), MYGENE_CHUNK):
        chunk = to_fetch[i:i + MYGENE_CHUNK]
        n_primary_chunks += 1
        hits = _mygene_post(chunk, scopes, taxon_id, verbose=verbose)

        if hits is None:
            n_fetch_failures += 1
            n_primary_chunk_failures += 1
            failed_query_keys.update(q.lower() for q in chunk)
            fetched_count += len(chunk)
            if progress_callback:
                try:
                    progress_callback(fetched_count, len(to_fetch), "mapping identifiers (retrying failures)")
                except Exception:
                    pass
            continue

        by_query = {}
        for hit in hits:
            q = str(hit.get("query", "")).strip()
            if not q:
                continue
            qk = q.lower()
            if hit.get("notfound"):
                by_query.setdefault(qk, None)
                continue
            hit_counts[qk] = hit_counts.get(qk, 0) + 1
            if by_query.get(qk) is None:
                by_query[qk] = _parse_mygene_hit(hit)

        for q in chunk:
            qk = q.lower()
            parsed = by_query.get(qk)
            if parsed and parsed.get("symbol"):
                result = _result_for_query(id_type, q, parsed)
                results[qk] = result
                _set_cache(f"{taxon_id}|{id_type}|{qk}", result)
            else:
                results[qk] = None

        fetched_count += len(chunk)
        if progress_callback:
            try:
                progress_callback(fetched_count, len(to_fetch), "mapping identifiers")
            except Exception:
                pass

    if to_fetch and n_primary_chunks > 0 and n_primary_chunk_failures == n_primary_chunks:
        raise RuntimeError(
            f"mygene.info identifier mapping failed for all {len(to_fetch)} identifier(s) "
            "after retries - this looks like a network/connectivity problem, not genes that "
            "don't exist. Please try again."
        )

    # Ensembl REST fallback: transcripts mygene didn't resolve directly. This
    # is a per-transcript loop (Ensembl has no bulk lookup for this endpoint),
    # against a service documented to cap at ~15 requests/second - throttled
    # here since a large transcript-level table could otherwise trip that
    # limit (this loop previously had no rate limiting or backoff at all).
    if id_type in ("ensembl_transcript", "mixed"):
        _ensembl_call_count = 0
        for q in to_fetch:
            qk = q.lower()
            if results.get(qk) is not None:
                continue
            if not _ENST_RE.match(q):
                continue
            if _ensembl_call_count > 0:
                time.sleep(0.08)  # ~12/s, safely under Ensembl's ~15/s cap
            _ensembl_call_count += 1
            parent = _ensembl_parent_gene(q, verbose=verbose)
            if not parent:
                continue
            parent_hits = _mygene_post([parent], "ensembl.gene", taxon_id, verbose=verbose)
            if not parent_hits:
                if parent_hits is None:
                    n_fetch_failures += 1
                continue
            parsed = _parse_mygene_hit(parent_hits[0])
            if parsed.get("symbol"):
                result = {
                    "symbol": parsed["symbol"], "entrez": parsed["entrez"],
                    "ensembl": parsed["ensembl"] or parent, "status": "mapped",
                }
                results[qk] = result
                _set_cache(f"{taxon_id}|{id_type}|{qk}", result)

    # Second-pass fallback: broaden the scope for anything STILL unresolved
    # (confirmed not-found under the primary scope, or unresolved because its
    # batch failed to fetch) - catches real cases like a deprecated/retired
    # symbol landing in a non-"gene_symbol" column, or a stray Ensembl/Entrez
    # ID inside a symbol column. Skipped when id_type is already "mixed" (same
    # scope, nothing new to try). Exact lookups only - never fuzzy matching.
    n_recovered_fallback = 0
    if id_type != "mixed":
        still_unresolved = [q for q in to_fetch if results.get(q.lower()) is None]
        for i in range(0, len(still_unresolved), MYGENE_CHUNK):
            chunk = still_unresolved[i:i + MYGENE_CHUNK]
            hits = _mygene_post(chunk, MYGENE_SCOPES["mixed"], taxon_id, verbose=verbose)
            if hits is None:
                n_fetch_failures += 1
                continue

            by_query = {}
            for hit in hits:
                q = str(hit.get("query", "")).strip()
                if not q:
                    continue
                qk = q.lower()
                if hit.get("notfound"):
                    by_query.setdefault(qk, None)
                    continue
                if by_query.get(qk) is None:
                    by_query[qk] = _parse_mygene_hit(hit)

            for q in chunk:
                qk = q.lower()
                parsed = by_query.get(qk)
                if parsed and parsed.get("symbol"):
                    result = _result_for_query(id_type, q, parsed)
                    results[qk] = result
                    _set_cache(f"{taxon_id}|{id_type}|{qk}", result)
                    n_recovered_fallback += 1

    n_ambiguous = sum(1 for c in hit_counts.values() if c > 1)
    n_ids_affected_by_failure = sum(1 for qk in failed_query_keys if results.get(qk) is None)

    save_idmap_cache()

    symbols, entrez_ids, ensembl_ids, statuses = [], [], [], []
    for raw in df["identifier"]:
        r = results.get(id_to_query[raw].lower())
        if r:
            symbols.append(r["symbol"])
            entrez_ids.append(r["entrez"])
            ensembl_ids.append(r["ensembl"])
            statuses.append(r["status"])
        else:
            symbols.append("")
            entrez_ids.append("")
            ensembl_ids.append("")
            statuses.append("unmapped")

    fetch_meta = {
        "n_recovered_fallback": n_recovered_fallback,
        "n_ambiguous": n_ambiguous,
        "n_fetch_failures": n_fetch_failures,
        "n_ids_affected_by_failure": n_ids_affected_by_failure,
    }

    out = df.copy()
    out["gene_symbol"] = symbols
    out["entrez_id"] = entrez_ids
    out["ensembl_id"] = ensembl_ids
    out["status"] = statuses
    return out, fetch_meta


# ---------------------------------------------------------------------------
# Collapsing (transcript -> gene, and duplicate gene symbols)
# ---------------------------------------------------------------------------

def _collapse_by_gene_symbol(df: pd.DataFrame, override_status: "str | None") -> "tuple[pd.DataFrame, int]":
    """Collapse rows sharing the same (non-empty) gene_symbol to one row —
    keeping the most-significant (smallest padj, tie-break largest |log2FC|).
    Rows with an empty gene_symbol (unmapped) pass through untouched."""
    if "gene_symbol" not in df.columns or df.empty:
        return df.copy(), 0

    has_symbol = df["gene_symbol"].astype(str).str.strip() != ""
    mapped = df[has_symbol].copy()
    unmapped = df[~has_symbol].copy()
    if mapped.empty:
        return df.copy(), 0

    def _sort_key(row):
        try:
            padj_val = float(row.get("padj"))
            if pd.isna(padj_val):
                padj_val = float("inf")
        except (TypeError, ValueError):
            padj_val = float("inf")
        try:
            abs_lfc = abs(float(row.get("log2FC", 0.0)))
        except (TypeError, ValueError):
            abs_lfc = 0.0
        return (padj_val, -abs_lfc)

    mapped["_sort_key"] = mapped.apply(_sort_key, axis=1)
    mapped = mapped.sort_values("_sort_key", kind="stable")
    grouped = mapped.groupby(mapped["gene_symbol"].str.lower(), sort=False)

    kept_rows = []
    n_collapsed = 0
    for _, group in grouped:
        if len(group) > 1:
            n_collapsed += len(group) - 1
            best = group.iloc[0].copy()
            if override_status:
                best["status"] = override_status
        else:
            best = group.iloc[0]
        kept_rows.append(best)

    result = pd.DataFrame(kept_rows).drop(columns=["_sort_key"], errors="ignore")
    if not unmapped.empty:
        result = pd.concat([result, unmapped], ignore_index=True)
    else:
        result = result.reset_index(drop=True)
    return result, n_collapsed


def collapse_transcripts_to_genes(df: pd.DataFrame, rule: str = "most_significant") -> "tuple[pd.DataFrame, int]":
    """Collapse multiple transcripts mapping to the same gene into one row."""
    return _collapse_by_gene_symbol(df, override_status="transcript_collapsed")


def deduplicate_genes(df: pd.DataFrame) -> "tuple[pd.DataFrame, int]":
    """Collapse any remaining duplicate gene symbols (would otherwise
    double-count in downstream enrichment)."""
    return _collapse_by_gene_symbol(df, override_status=None)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def resolve_de_table(
    source,
    taxon_id: str,
    identifier_type: "str | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> "tuple[pd.DataFrame, dict]":
    """Resolve an uploaded DE results table into the analysis-ready table every
    Transcriptomics tile consumes.

    Args:
        source:          A DataFrame, or a path to a CSV/TSV file.
        taxon_id:        NCBI/UniProt-style numeric taxon ID (e.g. "9606").
        identifier_type: Force a type instead of auto-detecting
                          ("gene_symbol"|"ensembl_gene"|"ensembl_transcript"|"entrez").
        progress_callback: Optional callable(current, total, label).
        verbose:          Print resolution steps.

    Returns:
        (resolved_df, meta) — resolved_df has RESOLVED_COLUMNS; meta carries
        resolution metrics (see module docstring / plan for the full key list).
    """
    if isinstance(source, pd.DataFrame):
        raw_df = source.copy()
    elif isinstance(source, (str, os.PathLike)):
        raw_df = load_de_table(str(source))
    else:
        raise TypeError("source must be a pandas DataFrame or a file path")

    n_input_rows = len(raw_df)
    column_mapping = detect_columns(raw_df)

    work = pd.DataFrame({
        "identifier": raw_df[column_mapping["identifier"]].astype(str).str.strip(),
        "log2FC": raw_df[column_mapping["log2FC"]],
        "pvalue": raw_df[column_mapping["pvalue"]],
    })
    work["padj"] = raw_df[column_mapping["padj"]] if "padj" in column_mapping else np.nan
    work = work[work["identifier"] != ""].reset_index(drop=True)

    work, n_dropped_nonnumeric = _coerce_numeric(work)

    padj_computed = False
    if "padj" not in column_mapping or work["padj"].isna().all():
        work["padj"] = _benjamini_hochberg(work["pvalue"])
        padj_computed = True

    detected_id_type = identifier_type or detect_identifier_type(work["identifier"])
    input_level = "transcript" if detected_id_type == "ensembl_transcript" else "gene"

    mapped, fetch_meta = map_identifiers(
        work, id_type=detected_id_type, taxon_id=str(taxon_id),
        progress_callback=progress_callback, verbose=verbose,
    )

    n_transcripts_collapsed = 0
    if detected_id_type in ("ensembl_transcript", "mixed"):
        mapped, n_transcripts_collapsed = collapse_transcripts_to_genes(mapped)

    mapped, n_duplicates_collapsed = deduplicate_genes(mapped)

    mapped = mapped.rename(columns={"identifier": "input_id"})
    n_mapped = int((mapped["status"] != "unmapped").sum())
    n_unmapped = int((mapped["status"] == "unmapped").sum())

    resolved_df = mapped[RESOLVED_COLUMNS].reset_index(drop=True)

    meta = {
        "identifier_type": detected_id_type,
        "taxon": str(taxon_id),
        "input_level": input_level,
        "column_mapping": column_mapping,
        "n_input_rows": n_input_rows,
        "n_mapped": n_mapped,
        "n_unmapped": n_unmapped,
        "n_duplicates_collapsed": n_duplicates_collapsed,
        "n_transcripts_collapsed": n_transcripts_collapsed,
        "n_dropped_nonnumeric": n_dropped_nonnumeric,
        "padj_computed": padj_computed,
        # From map_identifiers: recovery via the broadened-scope fallback,
        # ambiguous multi-matches, and mygene.info batch fetch failures (see
        # map_identifiers' docstring for exact semantics of each).
        "n_recovered_fallback": fetch_meta["n_recovered_fallback"],
        "n_ambiguous": fetch_meta["n_ambiguous"],
        "n_fetch_failures": fetch_meta["n_fetch_failures"],
        "n_ids_affected_by_failure": fetch_meta["n_ids_affected_by_failure"],
    }
    return resolved_df, meta
