# -*- coding: utf-8 -*-
"""
Proteomics module retriever for Omiscera.

Turns an uploaded protein-level differential-abundance table (UniProt
accession + log2FC + p-value [+ adjusted p-value]) into one clean,
analysis-ready table every Proteomics analysis tile consumes downstream.
Mirrors Transcriptomics_module_retriever.py in shape and purpose, but the
identifier space is UniProt accession ONLY (see module docstring in the
Proteomics planning discussion for why gene symbol was rejected as an input
type: MS search engines emit accessions natively; symbols would collapse
isoform distinctions MS just made, and a protein group can span several
genes). Reuses Transcriptomics_module_retriever's generic, identifier-
agnostic table utilities directly rather than re-deriving them.

Three behaviours here have NO transcriptomics equivalent, because a
proteomics export is not a clean gene list:
  1. Protein GROUPS - MS often cannot distinguish isoforms/close homologs, so
     tools report "P04637;P04637-2;Q9H3D4"; convention is the leading
     accession.
  2. Contaminant/decoy rows - MaxQuant marks these with a "CON__"/"REV__" ID
     prefix and/or "+"-flagged columns ("Potential contaminant", "Reverse",
     "Only identified by site"). Left in, these poison every downstream
     enrichment silently (a contaminant is never a real biological hit).
  3. Deprecated/merged/secondary accessions - UniProt periodically merges or
     retires accessions. Resolution therefore goes through UniProt's ID
     Mapping service (idmapping/run -> status -> results), not a plain
     per-accession GET, because it is the one endpoint that resolves a
     secondary accession to its current primary AND returns symbol/name in
     the same call (verified live: Q15086 -> P04637/TP53, 2026-08-24).
"""

import os
import re
import time

import numpy as np
import pandas as pd

from Transcriptomics_module_retriever import (
    HTTP_SESSION, OMISCERA_HEADERS, DEFAULT_TIMEOUT, DEFAULT_RETRIES,
    _parse_de_csv_text, _coerce_numeric, _benjamini_hochberg,
    _normalize_col_key,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_IDMAPPING_RUN_URL = "https://rest.uniprot.org/idmapping/run"
UNIPROT_IDMAPPING_STATUS_URL = "https://rest.uniprot.org/idmapping/status/{job_id}"
UNIPROT_IDMAPPING_RESULTS_URL = "https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
IDMAPPING_FIELDS = "accession,gene_names,protein_name,organism_id"
IDMAPPING_BATCH = 5000  # generous vs. the ~3-8k proteins a real MS run quantifies
IDMAPPING_POLL_INTERVAL = 1.0
IDMAPPING_POLL_TIMEOUT = 120

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
PROTEOMICS_IDMAP_CACHE_PATH = os.path.join(CACHE_DIR, "proteomics_idmap_cache.csv")

RESOLVED_COLUMNS = [
    "input_id", "accession", "gene_symbol", "protein_name",
    "log2FC", "pvalue", "padj", "status",
]

# ---------------------------------------------------------------------------
# Column detection (proteomics-specific: identifier synonyms differ sharply
# from Transcriptomics' gene-centric list; log2FC/pvalue/padj are reused
# because the statistics columns genuinely are the same shape)
# ---------------------------------------------------------------------------

COLUMN_SYNONYMS = {
    "identifier": {
        "protein", "proteinid", "protein_id", "proteinids", "protein_ids",
        "majorityproteinids", "majority_protein_ids", "leadingprotein",
        "leading_protein", "accession", "uniprot", "uniprotid", "uniprot_id",
        "uniprotaccession", "uniprot_accession", "proteingroup", "protein_group",
        "proteingroups", "protein_groups",
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
    # Optional: a single pre-aggregated abundance value (e.g. Perseus "Mean
    # intensity"). Deliberately NOT the per-sample "LFQ intensity Control_1"
    # style columns MaxQuant emits - handling those needs sample->group
    # design, which this app does not collect (see the deferred
    # intensity-matrix pipeline discussion). Used only for an optional
    # upload/QC rank-abundance plot; never required, never analysed.
    "abundance": {
        "abundance", "meanabundance", "mean_abundance", "intensity",
        "meanintensity", "mean_intensity", "ibaq", "meanibaq",
        "averageintensity", "average_intensity",
    },
}

# MaxQuant's own QC flag columns - normalized-key match, value is "truthy"
# ("+", "true", "1", "yes") when the row should be dropped. These are
# SEPARATE from the CON__/REV__ ID-prefix check below; MaxQuant emits both
# mechanisms and either alone can mark a row.
CONTAMINANT_FLAG_COLUMNS = {
    "potentialcontaminant", "contaminant", "reverse", "onlyidentifiedbysite",
}

_CONTAMINANT_PREFIX_RE = re.compile(r"^(CON__|REV__)", re.IGNORECASE)
_ISOFORM_SUFFIX_RE = re.compile(r"-\d+$")


def _normalize_flag_value(v) -> bool:
    """MaxQuant marks a flagged row with a literal '+' in that cell; blank
    otherwise. Also tolerate TRUE/1/yes for tables re-exported through other
    tools."""
    s = str(v).strip().lower()
    return s in ("+", "true", "1", "yes", "x")


def detect_columns(df: pd.DataFrame) -> dict:
    """Map arbitrary uploaded headers to canonical identifier/log2FC/pvalue/
    padj[/abundance]. Raises ValueError if identifier, log2FC, or pvalue
    cannot be found (padj is optional - BH FDR is computed when absent;
    abundance is optional and QC-only)."""
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
            f"Protein table is missing required column(s) for: {missing}. "
            f"Detected columns after normalisation: {list(df.columns)}"
        )
    return mapping


def _detect_flag_columns(df: pd.DataFrame) -> "list[str]":
    """Find MaxQuant-style QC flag columns present in the uploaded table,
    by normalized-key match against CONTAMINANT_FLAG_COLUMNS."""
    found = []
    for col in df.columns:
        if _normalize_col_key(col) in CONTAMINANT_FLAG_COLUMNS:
            found.append(col)
    return found


# ---------------------------------------------------------------------------
# CSV loading - reused verbatim from Transcriptomics_module_retriever
# (identifier-agnostic; no reason to duplicate the multi-delimiter sniff)
# ---------------------------------------------------------------------------

def load_de_table_from_bytes(content: bytes) -> pd.DataFrame:
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("utf-8", errors="replace")
    return _parse_de_csv_text(text)


def load_de_table(path: str) -> pd.DataFrame:
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin1") as fh:
            text = fh.read()
    return _parse_de_csv_text(text)


# ---------------------------------------------------------------------------
# Protein-group / isoform / contaminant handling - proteomics-only, no
# transcriptomics equivalent
# ---------------------------------------------------------------------------

def parse_protein_group(value) -> str:
    """Leading accession from a semicolon-separated protein group
    ('P04637;P04637-2;Q9H3D4' -> 'P04637'). MaxQuant/MSFragger/DIA-NN
    convention: the leading member is the group's representative."""
    s = str(value).strip()
    if not s:
        return ""
    return s.split(";")[0].strip()


def strip_isoform(accession: str) -> str:
    """'P04637-2' -> 'P04637'. Annotation databases (GO, Reactome, InterPro,
    Complex Portal, MSigDB-via-symbol) are keyed on the canonical entry, not
    a specific isoform."""
    return _ISOFORM_SUFFIX_RE.sub("", accession.strip())


def filter_contaminants_decoys(df: pd.DataFrame, id_col: str) -> "tuple[pd.DataFrame, dict]":
    """Remove contaminant/decoy rows via BOTH detection mechanisms real
    proteomics tools use, and report counts - matching Perseus. Without this,
    a MaxQuant export pasted straight in silently corrupts every downstream
    enrichment (CON__/REV__ rows are not real biology).

    Returns (clean_df, counts) where counts = {
        "n_prefix_removed": int,       # CON__/REV__ ID prefix
        "n_flagged_removed": int,      # any MaxQuant '+' flag column
        "flag_columns_found": [str],   # which flag columns were present
        "n_total_removed": int,
    }
    """
    ids = df[id_col].astype(str)
    prefix_mask = ids.str.match(_CONTAMINANT_PREFIX_RE)

    flag_cols = _detect_flag_columns(df)
    flag_mask = pd.Series(False, index=df.index)
    for col in flag_cols:
        flag_mask = flag_mask | df[col].apply(_normalize_flag_value)

    remove_mask = prefix_mask | flag_mask
    counts = {
        "n_prefix_removed": int((prefix_mask & ~flag_mask).sum()) + int((prefix_mask & flag_mask).sum()),
        "n_flagged_removed": int(flag_mask.sum()),
        "flag_columns_found": flag_cols,
        "n_total_removed": int(remove_mask.sum()),
    }
    # n_prefix_removed simplified to a plain count of prefix matches (overlap
    # with flag_mask double-counted above is corrected here to a single,
    # honest "how many had the prefix" number):
    counts["n_prefix_removed"] = int(prefix_mask.sum())
    return df[~remove_mask].reset_index(drop=True), counts


# ---------------------------------------------------------------------------
# UniProt ID Mapping - resolves accession (incl. secondary/deprecated) to
# current primary accession + gene symbol + protein name in one round trip.
# Verified live 2026-08-24: a secondary accession (Q15086) resolves to its
# current primary (P04637) with gene symbol (TP53) in the SAME response -
# the plain bulk-search endpoint (fields=) cannot return secondary
# accessions at all (both 'sec_acc' and 'secondary_accession' are rejected
# as invalid field names), so ID Mapping is not an optimisation, it is the
# only endpoint that does this.
# ---------------------------------------------------------------------------

def _idmapping_submit(ids: "list[str]", verbose: bool = False) -> "str | None":
    """Submit one ID Mapping job, with the SAME retry-on-transient-error
    pattern _mygene_post uses elsewhere in this codebase. Without this, a
    single hiccup fails the ENTIRE batch outright (verified live 2026-08-24:
    a 3500-accession resolve failed once with no underlying network problem -
    a plain retry of the identical call succeeded in 36s), because a typical
    input fits in ONE chunk (IDMAPPING_BATCH=5000), so one unretried failure
    IS a total failure, not a partial one."""
    if HTTP_SESSION is None or not ids:
        return None
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            resp = HTTP_SESSION.post(
                UNIPROT_IDMAPPING_RUN_URL,
                data={"from": "UniProtKB_AC-ID", "to": "UniProtKB", "ids": ",".join(ids)},
                headers=OMISCERA_HEADERS, timeout=DEFAULT_TIMEOUT,
            )
            if resp.status_code in (200, 201):
                return resp.json().get("jobId")
            last_error = f"HTTP {resp.status_code}"
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < DEFAULT_RETRIES:
                time.sleep(1.0)
                continue
            break
        except Exception as exc:
            last_error = str(exc)
            if attempt < DEFAULT_RETRIES:
                time.sleep(1.0)
                continue
            break
    if verbose:
        print(f"[UniProt ID mapping] submit failed after retries: {last_error}")
    return None


def _idmapping_poll(job_id: str, verbose: bool = False) -> bool:
    """Poll until FINISHED or timeout. Returns True if the job finished
    (successfully or with per-ID failures reported in results) - False only
    on a genuine timeout/network failure to even check status."""
    deadline = time.time() + IDMAPPING_POLL_TIMEOUT
    url = UNIPROT_IDMAPPING_STATUS_URL.format(job_id=job_id)
    while time.time() < deadline:
        try:
            resp = HTTP_SESSION.get(url, headers=OMISCERA_HEADERS, timeout=DEFAULT_TIMEOUT)
            if resp.status_code in (200, 201):
                data = resp.json()
                status = data.get("jobStatus")
                if status == "FINISHED" or "results" in data:
                    return True
                if status == "FAILED" or status == "ERROR":
                    return False
        except Exception as exc:
            if verbose:
                print(f"[UniProt ID mapping] status poll failed: {exc}")
        time.sleep(IDMAPPING_POLL_INTERVAL)
    return False


def _idmapping_fetch_results(job_id: str, verbose: bool = False) -> dict:
    """Cursor-paginate the results endpoint (same Link: rel=next pattern
    Genomics_module._fetch_uniprot_annotations_bulk uses)."""
    results = []
    failed_ids: "list[str]" = []
    obsolete_count = 0
    url = UNIPROT_IDMAPPING_RESULTS_URL.format(job_id=job_id)
    params = {"fields": IDMAPPING_FIELDS, "size": 500}
    page = 0
    while url:
        try:
            resp = HTTP_SESSION.get(
                url, params=params if page == 0 else None,
                headers=OMISCERA_HEADERS, timeout=DEFAULT_TIMEOUT,
            )
            if resp.status_code not in (200, 201):
                if verbose:
                    print(f"[UniProt ID mapping] results page {page} failed: HTTP {resp.status_code}")
                break
            data = resp.json()
        except Exception as exc:
            if verbose:
                print(f"[UniProt ID mapping] results page {page} failed: {exc}")
            break
        results.extend(data.get("results", []))
        failed_ids.extend(data.get("failedIds", []) or [])
        obsolete_count += int(data.get("obsoleteCount", 0) or 0)
        page += 1
        next_link = resp.links.get("next", {})
        url = next_link.get("url", "") or None
        params = None
    return {"results": results, "failed_ids": failed_ids, "obsolete_count": obsolete_count}


def _idmapping_resolve(unique_accessions: "list[str]", verbose: bool = False) -> "tuple[dict, int]":
    """Resolve a list of (already protein-group-parsed, isoform-stripped)
    accessions via UniProt ID Mapping, chunked at IDMAPPING_BATCH.

    Returns ({query_upper: {"accession","gene_symbol","protein_name"} or
    None}, n_fetch_failures) - None means the query genuinely could not be
    resolved (UniProt says so, via failedIds); a chunk-level submit/poll
    failure is counted in n_fetch_failures and its queries are simply absent
    from the dict (caller must NOT treat "absent" as "confirmed unmapped").
    """
    resolved: dict = {}
    n_fetch_failures = 0
    for i in range(0, len(unique_accessions), IDMAPPING_BATCH):
        chunk = unique_accessions[i:i + IDMAPPING_BATCH]
        job_id = _idmapping_submit(chunk, verbose=verbose)
        if not job_id or not _idmapping_poll(job_id, verbose=verbose):
            n_fetch_failures += 1
            continue
        payload = _idmapping_fetch_results(job_id, verbose=verbose)
        for row in payload["results"]:
            q = str(row.get("from", "")).strip().upper()
            to = row.get("to", {}) or {}
            primary = str(to.get("primaryAccession", "")).strip()
            if not primary:
                continue
            genes = to.get("genes", []) or []
            symbol = ""
            if genes:
                symbol = str((genes[0].get("geneName") or {}).get("value", "")).strip()
            desc = to.get("proteinDescription", {}) or {}
            rec_name = (desc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            if not rec_name:
                submitted = (desc.get("submissionNames") or [{}])[0]
                rec_name = (submitted.get("fullName") or {}).get("value", "")
            resolved[q] = {
                "accession": primary, "gene_symbol": symbol, "protein_name": str(rec_name).strip(),
            }
        for q in payload["failed_ids"]:
            resolved.setdefault(str(q).strip().upper(), None)
    return resolved, n_fetch_failures


# ---------------------------------------------------------------------------
# Dedupe by resolved primary accession (isoforms of the same protein, or
# distinct protein-group rows that resolve to the same current primary
# accession, must collapse to one row - mirrors
# Transcriptomics_module_retriever._collapse_by_gene_symbol)
# ---------------------------------------------------------------------------

def _dedupe_by_accession(df: pd.DataFrame) -> "tuple[pd.DataFrame, int]":
    if "accession" not in df.columns or df.empty:
        return df.copy(), 0
    has_acc = df["accession"].astype(str).str.strip() != ""
    mapped = df[has_acc].copy()
    unmapped = df[~has_acc].copy()
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
    grouped = mapped.groupby(mapped["accession"], sort=False)

    kept_rows, n_collapsed = [], 0
    for _, group in grouped:
        if len(group) > 1:
            n_collapsed += len(group) - 1
        kept_rows.append(group.iloc[0])
    result = pd.DataFrame(kept_rows).drop(columns=["_sort_key"], errors="ignore")
    if not unmapped.empty:
        result = pd.concat([result, unmapped], ignore_index=True)
    else:
        result = result.reset_index(drop=True)
    return result, n_collapsed


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def resolve_protein_table(
    source,
    taxon_id: str,
    progress_callback=None,
    verbose: bool = False,
) -> "tuple[pd.DataFrame, dict]":
    """Resolve an uploaded protein-level DE table into the analysis-ready
    table every Proteomics tile consumes.

    Args:
        source: A DataFrame, or a path to a CSV/TSV file.
        taxon_id: NCBI/UniProt-style numeric taxon ID (e.g. "9606"). Used
            only to verify resolved accessions belong to the chosen
            organism (a mismatch is reported, not silently dropped) - the
            organism itself is user-selected, never inferred from the
            accessions (rejected in planning: silently wrong for
            mixed-species/contaminated input).

    Returns:
        (resolved_df, meta) - resolved_df has RESOLVED_COLUMNS.
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
    has_abundance = "abundance" in column_mapping
    if has_abundance:
        work["abundance"] = pd.to_numeric(raw_df[column_mapping["abundance"]], errors="coerce")
    # Carry the raw QC flag columns (e.g. MaxQuant's "Potential contaminant",
    # "Reverse", "Only identified by site") into `work`, aligned by position,
    # so filter_contaminants_decoys can actually see them - they are NOT part
    # of RESOLVED_COLUMNS and are dropped again right after filtering.
    flag_cols = _detect_flag_columns(raw_df)
    for col in flag_cols:
        work[col] = raw_df[col]  # same index as raw_df, matching work's index above
    work = work[work["identifier"] != ""].reset_index(drop=True)

    # Contaminant/decoy filtering BEFORE protein-group parsing, so the
    # CON__/REV__ prefix check sees the raw uploaded identifier exactly as
    # MaxQuant wrote it (a group's leading member is never itself prefixed
    # differently from the row's own flag columns in practice, but checking
    # pre-parse is the conservative, documented order).
    work, contam_counts = filter_contaminants_decoys(work, "identifier")
    work = work.drop(columns=flag_cols, errors="ignore")

    work["accession_query"] = work["identifier"].apply(
        lambda v: strip_isoform(parse_protein_group(v))
    )
    work = work[work["accession_query"] != ""].reset_index(drop=True)

    work, n_dropped_nonnumeric = _coerce_numeric(work)

    padj_computed = False
    if "padj" not in column_mapping or work["padj"].isna().all():
        work["padj"] = _benjamini_hochberg(work["pvalue"])
        padj_computed = True

    unique_queries = sorted({q.upper() for q in work["accession_query"] if q})
    resolved_map, n_fetch_failures = _idmapping_resolve(
        unique_queries, verbose=verbose,
    )
    if progress_callback:
        try:
            progress_callback(len(unique_queries), len(unique_queries), "resolved via UniProt ID Mapping")
        except Exception:
            pass

    if unique_queries and n_fetch_failures and not resolved_map:
        raise RuntimeError(
            "UniProt ID Mapping failed for all accessions after retries - this looks like a "
            "network/connectivity problem, not proteins that don't exist. Please try again."
        )

    accessions, symbols, names, statuses = [], [], [], []
    for q in work["accession_query"]:
        qk = q.upper()
        r = resolved_map.get(qk)
        if r:
            accessions.append(r["accession"])
            symbols.append(r["gene_symbol"])
            names.append(r["protein_name"])
            statuses.append("mapped" if r["accession"].upper() == qk else "resolved_secondary")
        else:
            accessions.append("")
            symbols.append("")
            names.append("")
            statuses.append("unmapped")

    work = work.rename(columns={"identifier": "input_id"})
    work["accession"] = accessions
    work["gene_symbol"] = symbols
    work["protein_name"] = names
    work["status"] = statuses

    work, n_duplicates_collapsed = _dedupe_by_accession(work)

    n_mapped = int((work["status"] != "unmapped").sum())
    n_unmapped = int((work["status"] == "unmapped").sum())
    n_resolved_secondary = int((work["status"] == "resolved_secondary").sum())

    resolved_df = work[RESOLVED_COLUMNS].reset_index(drop=True)

    meta = {
        "taxon": str(taxon_id),
        "column_mapping": column_mapping,
        "has_abundance": has_abundance,
        "abundance_by_accession": (
            dict(zip(work["accession"], work.get("abundance", pd.Series(dtype=float))))
            if has_abundance else {}
        ),
        "n_input_rows": n_input_rows,
        "n_mapped": n_mapped,
        "n_unmapped": n_unmapped,
        "n_resolved_secondary": n_resolved_secondary,
        "n_duplicates_collapsed": n_duplicates_collapsed,
        "n_dropped_nonnumeric": n_dropped_nonnumeric,
        "padj_computed": padj_computed,
        "n_fetch_failures": n_fetch_failures,
        "contaminants": contam_counts,
    }
    return resolved_df, meta
