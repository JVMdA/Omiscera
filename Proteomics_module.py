# -*- coding: utf-8 -*-
"""
Proteomics statistical analysis module for Omiscera.

Hosts Over-Representation Analysis (ORA) for Proteomics. Identifier space is
UniProt accession throughout (see Proteomics_module_retriever.py), which is
why this file builds on Genomics_module's ORA engine rather than
Transcriptomics_module's: Genomics ORA is already accession-native
(_run_fisher_for_namespace, build_go_universe, build_reactome_universe,
load_custom_geneset_db, load_user_background_list all operate in accession
space already), whereas Transcriptomics ORA is gene-symbol-native. Reused
directly rather than re-derived - this file adds only what Genomics does not
already have: the proteomics-specific databases and the accession<->symbol
bridge MSigDB needs.

Databases
---------
  GO, Reactome          - reuse Genomics_module.build_go_universe /
                           build_reactome_universe verbatim (same accession
                           space, same native sources).
  MSigDB                - gene-symbol-native (reuses
                           Transcriptomics_module.fetch_msigdb_sets); tested
                           via the resolved accession->gene_symbol bridge,
                           the same "symbol-native pathway inside an
                           otherwise accession-native module" pattern
                           Genomics already uses for GWAS Catalog.
  UniProt (disease)      \\
  UniProt (keywords)      |  ALL FOUR sourced from ONE UniProt bulk fetch
  InterPro                |  per taxon (_fetch_uniprot_annotations_bulk),
                           /  verified live 2026-08-24: a single query with
                           fields=accession,xref_interpro,xref_complexportal,
                           cc_disease,keyword returns InterPro IDs, disease
                           annotations (name embedded), and keywords (name
                           embedded) in one response. Only InterPro needs a
                           separate name lookup (entry.list); disease and
                           keyword names arrive embedded in the same call.
  Complex Portal          - NOT sourced via the UniProt cross-reference (it
                           gives membership only). Complex Portal's own
                           ComplexTAB file (verified live 2026-08-24) gives
                           BOTH names and membership in one per-taxon file -
                           a cleaner, more authoritative source.
  Custom database         - reuses Genomics_module.load_custom_geneset_db
                           unchanged (id_type="accession" is its default).

NOT included: PANTHER. Verified live 2026-08-24 that pantherdb.org and
ftp.pantherdb.org are behind Cloudflare bot protection - every request
(varied User-Agent, HEAD, GET) returns 403 or a JS challenge page, including
for what should be static FTP files. This is a live access barrier, not a
licensing issue; revisit if/when a working access path is found.
"""

import csv
import io
import os
import re
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd

try:
    import requests as _requests
    HTTP_SESSION = _requests.Session()
except ImportError:
    _requests = None
    HTTP_SESSION = None

from Genomics_module import (
    OMISCERA_HEADERS, DEFAULT_TIMEOUT, DEFAULT_RETRIES,
    UNIPROT_SEARCH_URL, UNIPROT_PAGE_SIZE,
    MIN_TERM_SIZE, MAX_TERM_SIZE_DEFAULT, FDR_ALPHA_DEFAULT, ORA_RESULT_COLUMNS,
    BURDEN_RESULT_COLUMNS,
    benjamini_hochberg, bonferroni, fisher_exact_pvalue, fold_enrichment,
    _run_fisher_for_namespace, _run_wilcoxon_for_namespace,
    build_go_universe, build_reactome_universe,
    load_custom_geneset_db, load_user_background_list,
    _fetch_uniprot_annotations_bulk,
    _load_universe_cache, _save_universe_cache,
)
from Transcriptomics_module import (
    fetch_msigdb_sets, msigdb_collection_available, MSIGDB_COLLECTIONS,
    _msigdb_selected_keys, GSEA_RANKING_METRICS, GSEA_MIN_P,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PADJ_CUTOFF_DEFAULT = 0.05
LFC_CUTOFF_DEFAULT = 1.0
DIRECTIONS = ("all", "up", "down")

# PANTHER intentionally excluded - see module docstring.
PROTEOMICS_ORA_DATABASES = [
    "GO", "Reactome", "MSigDB", "UniProt (disease)", "UniProt (keywords)",
    "InterPro", "Complex Portal", "DISEASES", "Custom database",
]

# DISEASES (Jensen Lab, https://diseases.jensenlab.org) - added because
# UniProt's own disease annotations, while correct, are structurally
# unusable for enrichment: measured live 2026-08-24, 97.5% of all 6862
# human UniProt disease terms have exactly ONE gene ("Diamond-Blackfan
# anemia 22", "Diamond-Blackfan anemia 19", ...), so only 77 terms (1.1%)
# survive even a minimum-size-3 gate. DISEASES aggregates the same kind of
# knowledge into real disease CATEGORIES instead: 396 of its 1115 curated
# terms (35.5%) clear K>=3, and its largest terms are things like
# "Neurodegenerative disease" (442 genes), "Cancer" (326), "Epilepsy" (324).
# Both are kept - UniProt's narrow monogenic entries are still right for
# someone studying one specific rare disease.
#
# ONLY the manually-curated "knowledge" channel is used. The other two
# channels were deliberately rejected: "experiments" is GWAS-scale (its
# largest term has 2535 genes - far past any sensible max-term-size gate
# and uninformative as a "disease gene set"), and "textmining" is
# lower-confidence still.
#
# Licence: CC BY 4.0 (verified live 2026-08-24) - the same licence family
# as STRING, from the same lab, so it carries no new licensing risk.
# HUMAN ONLY by design (it is a human-disease resource); other organisms
# get an explicit empty universe with a caveat, never a silent zero.
DISEASES_KNOWLEDGE_URL = "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv"
DISEASES_LOCAL_NAME = "diseases_knowledge_filtered.tsv"
DISEASES_REFRESH_DAYS = 30
DISEASES_HUMAN_TAXON_ID = "9606"

# Quantitative Enrichment's default ranking metric. Deliberately NOT
# GSEA_RANKING_METRIC_DEFAULT ("signed_log_pvalue" - the right default for
# permutation-based GSEA, which needs a well-behaved monotonic ranking
# statistic): Perseus' 1D annotation enrichment, the method this module
# mirrors, is defined directly on log2FC, so that is the natural default
# here. GSEA_RANKING_METRICS (the set of valid choices) is still reused
# unchanged - only the DEFAULT differs.
QUANT_METRIC_DEFAULT = "log2fc"

# Proteomics renames the shared engine's "Overlapping Genes" column: in a
# proteomics module these are the measured PROTEINS that fall inside a term,
# not genes. The rename also fixes a real inconsistency - the same column
# silently held two different identifier types depending on the database
# (UniProt accessions for the 6 accession-keyed sources, gene SYMBOLS for
# the 2 symbol-keyed ones, MSigDB and DISEASES). _finalize_member_column()
# below normalises every database to one readable "P05386 (RPLP0)" form.
OVERLAP_COLUMN = "Overlapping Proteins"
PROTEOMICS_ORA_RESULT_COLUMNS = [
    OVERLAP_COLUMN if c == "Overlapping Genes" else c for c in ORA_RESULT_COLUMNS
]
# Quantitative Enrichment reuses Genomics' Gene-set Burden Test engine
# (_run_wilcoxon_for_namespace) - same member-column rename as ORA, applied
# to that engine's own column name ("Genes", not "Overlapping Genes").
PROTEOMICS_BURDEN_RESULT_COLUMNS = [
    OVERLAP_COLUMN if c == "Genes" else c for c in BURDEN_RESULT_COLUMNS
]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")

# InterPro entry.list - global (species-agnostic) ID->name file, verified
# live 2026-08-24: tab-separated ENTRY_AC / ENTRY_TYPE / ENTRY_NAME, ~2.9MB.
INTERPRO_ENTRY_LIST_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
INTERPRO_ENTRY_LIST_LOCAL_NAME = "interpro_entry_list.tsv"
INTERPRO_ENTRY_LIST_PATH = os.path.join(CACHE_DIR, INTERPRO_ENTRY_LIST_LOCAL_NAME)
INTERPRO_REFRESH_DAYS = 30

# Complex Portal ComplexTAB - per-taxon file, names+membership together,
# verified live 2026-08-24 (9606.tsv exists, 5.0MB, tab-separated with a
# '#'-prefixed header; column 5 is pipe-separated ACCESSION(stoichiometry)).
COMPLEXPORTAL_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab"
COMPLEXPORTAL_REFRESH_DAYS = 30

# UniProt accession pattern (used to filter Complex Portal's membership
# column, which can also reference OTHER complexes e.g. "CPX-6" - those are
# not proteins and must not enter the universe as if they were).
_UNIPROT_ACC_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")


# ---------------------------------------------------------------------------
# Significant-protein selection (mirrors Transcriptomics_module.select_degs,
# but keys on 'accession' - Proteomics' primary identifier - not gene_symbol)
# ---------------------------------------------------------------------------

def select_significant_proteins(
    resolved_df: pd.DataFrame, padj_cutoff: float, lfc_cutoff: float, direction: str = "all",
) -> pd.DataFrame:
    """Select significant proteins from the resolved table.

    direction: "all" (either direction), "up" (log2FC >= +lfc_cutoff), or
    "down" (log2FC <= -lfc_cutoff). Rows with no resolved accession (unmapped)
    or a missing padj are never selected (NaN comparisons evaluate False)."""
    df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""].copy()
    mask = (df["padj"] <= padj_cutoff) & (df["log2FC"].abs() >= lfc_cutoff)
    if direction == "up":
        mask &= df["log2FC"] >= lfc_cutoff
    elif direction == "down":
        mask &= df["log2FC"] <= -lfc_cutoff
    return df[mask]


# ---------------------------------------------------------------------------
# InterPro entry.list (global name lookup, not per-taxon)
# ---------------------------------------------------------------------------

def _ensure_interpro_entry_list(verbose: bool = False) -> str:
    """Download entry.list if missing/stale. Global file, not per-taxon."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(INTERPRO_ENTRY_LIST_PATH):
        age_days = (time.time() - os.path.getmtime(INTERPRO_ENTRY_LIST_PATH)) / 86400.0
        if age_days < INTERPRO_REFRESH_DAYS:
            return INTERPRO_ENTRY_LIST_PATH
    if HTTP_SESSION is None:
        if os.path.exists(INTERPRO_ENTRY_LIST_PATH):
            return INTERPRO_ENTRY_LIST_PATH
        raise RuntimeError("requests library not available; cannot fetch InterPro entry.list")
    try:
        resp = HTTP_SESSION.get(INTERPRO_ENTRY_LIST_URL, headers=OMISCERA_HEADERS, timeout=120)
        if resp.status_code == 200:
            tmp_path = INTERPRO_ENTRY_LIST_PATH + ".tmp"
            with open(tmp_path, "wb") as fh:
                fh.write(resp.content)
            os.replace(tmp_path, INTERPRO_ENTRY_LIST_PATH)
            return INTERPRO_ENTRY_LIST_PATH
    except Exception as exc:
        if verbose:
            print(f"[InterPro entry.list] fetch failed: {exc}")
    if os.path.exists(INTERPRO_ENTRY_LIST_PATH):
        return INTERPRO_ENTRY_LIST_PATH  # stale but present beats nothing
    raise RuntimeError("Could not fetch InterPro entry.list and no cached copy exists.")


def _load_interpro_names(verbose: bool = False) -> dict:
    path = _ensure_interpro_entry_list(verbose=verbose)
    names = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header: ENTRY_AC / ENTRY_TYPE / ENTRY_NAME
        for row in reader:
            if len(row) >= 3:
                names[row[0].strip()] = row[2].strip()
    return names


# ---------------------------------------------------------------------------
# UniProt-cross-reference-derived universes: InterPro, UniProt disease,
# UniProt keywords - all three from ONE bulk fetch per taxon.
# ---------------------------------------------------------------------------

_XREF_FIELDS_BY_SOURCE = {
    "InterPro": "accession,xref_interpro",
    "UniProt (disease)": "accession,cc_disease",
    "UniProt (keywords)": "accession,keyword",
}


def _build_xref_universe(source: str, taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """Shared builder for the three UniProt-cross-reference-derived
    databases. Returns {term_id: {"name","genes": set}}, "universe_size",
    "fetched_at"."""
    cache_db = f"Proteomics_{source.replace(' ', '_').replace('(', '').replace(')', '')}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, taxon_id)
        if cached is not None:
            return cached

    fields = _XREF_FIELDS_BY_SOURCE[source]
    entries = _fetch_uniprot_annotations_bulk(taxon_id, fields, verbose=verbose)
    interpro_names = _load_interpro_names(verbose=verbose) if source == "InterPro" else None

    terms: dict = {}
    universe: set = set()
    for entry in entries:
        accession = str(entry.get("primaryAccession", "")).strip()
        if not accession:
            continue
        if source == "InterPro":
            for xref in entry.get("uniProtKBCrossReferences", []) or []:
                if xref.get("database") != "InterPro":
                    continue
                term_id = str(xref.get("id", "")).strip()
                if not term_id:
                    continue
                name = interpro_names.get(term_id, term_id)
                bucket = terms.setdefault(term_id, {"name": name, "genes": set()})
                bucket["genes"].add(accession)
                universe.add(accession)
        elif source == "UniProt (disease)":
            for comment in entry.get("comments", []) or []:
                if comment.get("commentType") != "DISEASE":
                    continue
                disease = comment.get("disease")
                if not disease:
                    continue
                term_id = str(disease.get("diseaseAccession", "")).strip()
                name = str(disease.get("diseaseId", "")).strip()
                if not term_id or not name:
                    continue
                bucket = terms.setdefault(term_id, {"name": name, "genes": set()})
                bucket["genes"].add(accession)
                universe.add(accession)
        elif source == "UniProt (keywords)":
            for kw in entry.get("keywords", []) or []:
                term_id = str(kw.get("id", "")).strip()
                name = str(kw.get("name", "")).strip()
                if not term_id or not name:
                    continue
                bucket = terms.setdefault(term_id, {"name": name, "genes": set()})
                bucket["genes"].add(accession)
                universe.add(accession)

    result = {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
    if entries:
        _save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[{source} universe] not caching: empty fetch (taxon={taxon_id})")
    return result


def build_interpro_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    return _build_xref_universe("InterPro", taxon_id, verbose=verbose, force_refresh=force_refresh)


def build_uniprot_disease_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    return _build_xref_universe("UniProt (disease)", taxon_id, verbose=verbose, force_refresh=force_refresh)


def build_uniprot_keyword_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    return _build_xref_universe("UniProt (keywords)", taxon_id, verbose=verbose, force_refresh=force_refresh)


# ---------------------------------------------------------------------------
# Complex Portal (own ComplexTAB file - names + membership together)
# ---------------------------------------------------------------------------

def _ensure_complexportal_file(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> "str | None":
    os.makedirs(CACHE_DIR, exist_ok=True)
    local_path = os.path.join(CACHE_DIR, f"complexportal_{taxon_id}.tsv")
    if not force_refresh and os.path.exists(local_path):
        age_days = (time.time() - os.path.getmtime(local_path)) / 86400.0
        if age_days < COMPLEXPORTAL_REFRESH_DAYS:
            return local_path
    if HTTP_SESSION is None:
        return local_path if os.path.exists(local_path) else None
    url = f"{COMPLEXPORTAL_BASE_URL}/{taxon_id}.tsv"
    try:
        resp = HTTP_SESSION.get(url, headers=OMISCERA_HEADERS, timeout=120)
        if resp.status_code == 200:
            tmp_path = local_path + ".tmp"
            with open(tmp_path, "wb") as fh:
                fh.write(resp.content)
            os.replace(tmp_path, local_path)
            return local_path
        # A 404 for this taxon means Complex Portal has no data for this
        # organism (verified live: only ~20 species have a file at all) -
        # distinct from a transient error, so not worth retrying blindly.
        if verbose:
            print(f"[Complex Portal] taxon={taxon_id}: HTTP {resp.status_code}")
    except Exception as exc:
        if verbose:
            print(f"[Complex Portal] fetch failed: {exc}")
    return local_path if os.path.exists(local_path) else None


def _ensure_diseases_file(verbose: bool = False, force_refresh: bool = False) -> "str | None":
    """Download DISEASES' curated knowledge channel if missing/stale.
    Single global human file (~600KB), not per-taxon."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    local_path = os.path.join(CACHE_DIR, DISEASES_LOCAL_NAME)
    if not force_refresh and os.path.exists(local_path):
        age_days = (time.time() - os.path.getmtime(local_path)) / 86400.0
        if age_days < DISEASES_REFRESH_DAYS:
            return local_path
    if HTTP_SESSION is None:
        return local_path if os.path.exists(local_path) else None
    try:
        resp = HTTP_SESSION.get(DISEASES_KNOWLEDGE_URL, headers=OMISCERA_HEADERS, timeout=120)
        if resp.status_code == 200 and resp.content:
            tmp_path = local_path + ".tmp"
            with open(tmp_path, "wb") as fh:
                fh.write(resp.content)
            os.replace(tmp_path, local_path)
            return local_path
        if verbose:
            print(f"[DISEASES] fetch failed: HTTP {resp.status_code}")
    except Exception as exc:
        if verbose:
            print(f"[DISEASES] fetch failed: {exc}")
    # A stale cached copy beats nothing, matching every other fetcher here.
    return local_path if os.path.exists(local_path) else None


def build_diseases_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """DISEASES (Jensen Lab) curated gene-disease universe, GENE-SYMBOL keyed.

    File columns (verified live 2026-08-24, tab-separated, no header):
        0 Ensembl protein id | 1 gene symbol | 2 DOID | 3 disease name
        4 source | 5 evidence type | 6 confidence score
    Only columns 1-3 are used; the Ensembl protein id is ignored because
    this module bridges through gene symbols (the same route MSigDB already
    uses here), not through Ensembl.

    Returns {"terms": {DOID: {"name","genes": set[SYMBOL]}}, "universe_size",
    "fetched_at", "caveat": str|None}. Symbols are upper-cased to match the
    accession->symbol map run_ora builds (which upper-cases too).
    """
    cache_db = "Proteomics_DISEASES"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, taxon_id)
        if cached is not None:
            return cached

    # Human-only resource, by design - not a fetch failure, a scope limit.
    if str(taxon_id).strip() != DISEASES_HUMAN_TAXON_ID:
        return {
            "terms": {}, "universe_size": 0,
            "fetched_at": datetime.now(timezone.utc).isoformat(),
            "caveat": (
                f"DISEASES is a human-only resource (taxon {DISEASES_HUMAN_TAXON_ID}); it has no data for "
                f"taxon {taxon_id}, so this database shows no results for this organism."
            ),
        }

    path = _ensure_diseases_file(verbose=verbose, force_refresh=force_refresh)
    terms: dict = {}
    universe: set = set()
    caveat = None
    if path is None:
        caveat = "Could not fetch the DISEASES database for this run (transient service error). Re-run to retry."
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) < 4:
                    continue
                symbol = row[1].strip().upper()
                doid = row[2].strip()
                name = row[3].strip()
                if not symbol or not doid or not name:
                    continue
                bucket = terms.setdefault(doid, {"name": name, "genes": set()})
                bucket["genes"].add(symbol)
                universe.add(symbol)

    result = {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": datetime.now(timezone.utc).isoformat(), "caveat": caveat,
    }
    if terms:
        _save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[DISEASES universe] not caching: empty (taxon={taxon_id})")
    return result


def build_complex_portal_universe(taxon_id: str, verbose: bool = False, force_refresh: bool = False) -> dict:
    """Complex Portal universe from the organism's own ComplexTAB file.
    Returns {"terms": {complex_ac: {"name","genes": set}}, "universe_size",
    "fetched_at", "caveat": str|None}. A missing per-taxon file (organism not
    covered by Complex Portal) yields an empty, explicitly-noted universe -
    never an exception, since "0 species this small database covers this
    organism" is a normal, expected outcome, not a fetch error."""
    cache_db = "Proteomics_ComplexPortal"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, taxon_id)
        if cached is not None:
            return cached

    path = _ensure_complexportal_file(taxon_id, verbose=verbose, force_refresh=force_refresh)
    terms: dict = {}
    universe: set = set()
    caveat = None
    if path is None:
        caveat = (
            f"Complex Portal has no data file for taxon {taxon_id}; it currently covers "
            "only ~20 well-studied species. This database will show no results for this "
            "organism."
        )
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            for row in reader:
                if len(row) < 5:
                    continue
                complex_ac = row[0].strip()
                name = row[1].strip()
                members_field = row[4].strip()
                if not complex_ac or not members_field:
                    continue
                genes = set()
                for token in members_field.split("|"):
                    acc = token.split("(")[0].strip()
                    if _UNIPROT_ACC_RE.match(acc):  # excludes nested "CPX-…" complex refs
                        genes.add(acc)
                if not genes:
                    continue
                terms[complex_ac] = {"name": name or complex_ac, "genes": genes}
                universe |= genes

    result = {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": datetime.now(timezone.utc).isoformat(), "caveat": caveat,
    }
    if terms:
        _save_universe_cache(cache_db, taxon_id, result)
    elif verbose:
        print(f"[Complex Portal universe] not caching: empty (taxon={taxon_id})")
    return result


# ---------------------------------------------------------------------------
# ORA orchestrator
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Shared helper: normalise the overlapping-member column (used by BOTH
# run_ora and run_quantitative_enrichment - see each caller for what column
# name their underlying engine produces).
# ---------------------------------------------------------------------------

def _finalize_member_column(
    tables: dict, accession_to_symbol: dict, source_col: str, target_columns: "list[str]",
) -> None:
    """Normalise the overlapping-member column across every block in
    `tables`, in place.

    Two problems this fixes, both real and user-reported (see Fix commit
    for the ORA case, `52ba651`):
      1. The shared Genomics engines name this column for genes ("Overlapping
         Genes" in the ORA engine, "Genes" in the Burden Test engine); in a
         proteomics module these are PROTEINS. `source_col` names whichever
         one the caller's engine actually produced.
      2. Its contents silently differ in identifier TYPE depending on the
         database - UniProt accessions from the 6 accession-keyed sources,
         bare gene symbols from the 2 symbol-keyed ones (MSigDB, DISEASES) -
         so the same column could mean different things on different rows of
         the same results screen.

    Both are normalised to "ACCESSION (SYMBOL)": unambiguous, and readable
    (P05386 alone tells a reader nothing; RPLP0 alone loses the isoform-level
    identity that is the whole point of using an accession).
    """
    symbol_to_accession: dict = {}
    for acc, sym in accession_to_symbol.items():
        if sym:
            symbol_to_accession.setdefault(sym, acc)

    def _display_member(member: str) -> str:
        m = member.strip()
        if not m:
            return m
        sym = accession_to_symbol.get(m)          # accession-space member
        if sym:
            return f"{m} ({sym})"
        acc = symbol_to_accession.get(m.upper())  # symbol-space member (MSigDB/DISEASES)
        if acc:
            return f"{acc} ({m.upper()})"
        return m   # unresolvable (e.g. a custom-DB gene never quantified) - left as-is

    for block in tables.values():
        tbl = block.get("table")
        if tbl is None:
            continue
        if source_col in tbl.columns:
            if not tbl.empty:
                tbl[source_col] = tbl[source_col].map(
                    lambda cell: ", ".join(_display_member(x) for x in str(cell).split(",") if x.strip())
                )
            tbl = tbl.rename(columns={source_col: OVERLAP_COLUMN})
        # Reindex so an empty block still carries the proteomics column set.
        block["table"] = tbl.reindex(columns=target_columns)


def _mean_score_by_key(keys: "list[str]", values: "np.ndarray") -> dict:
    """Aggregate parallel (key, value) pairs into {key: mean(values)}.

    Needed in two places quantitative enrichment builds a score dict: (a)
    accession-level, where the SAME resolved accession can legitimately
    appear on more than one input row (e.g. two uploaded protein-group rows
    both resolving to the same canonical accession - resolve_protein_table
    does not deduplicate by accession); (b) the accession->symbol bridge for
    MSigDB/DISEASES, where distinct accessions (isoforms, paralogs sharing an
    annotation) can map to the same gene symbol. A naive dict build would
    silently let the LAST row win, an arbitrary, row-order-dependent choice
    - the mean is the same aggregation a protein/gene-level rollup normally
    uses and does not depend on input order.
    """
    out: dict = {}
    counts: dict = {}
    for k, v in zip(keys, values):
        if k in out:
            out[k] += v
            counts[k] += 1
        else:
            out[k] = v
            counts[k] = 1
    return {k: out[k] / counts[k] for k in out}


def _build_protein_scores(resolved_df: pd.DataFrame, metric: str = QUANT_METRIC_DEFAULT) -> "tuple[dict, int]":
    """Build an accession -> score dict for Quantitative Enrichment's
    competitive Wilcoxon test.

    Reuses the SAME four metric formulas as Transcriptomics' GSEA ranking
    (GSEA_RANKING_METRICS / Transcriptomics_module.build_ranked_list) for
    consistency across modules, rather than re-deriving them - but is NOT
    build_ranked_list itself: this is keyed on 'accession' (this module's
    identifier), does not require a resolved gene_symbol, and returns an
    unordered dict (the Wilcoxon engine partitions in-set/out-set itself; no
    ranked ordering is needed).

    Rows with a non-finite score (NaN log2FC/pvalue, or - for
    signed_log_padj only - a missing padj) are dropped. Duplicate accessions
    are mean-aggregated (see _mean_score_by_key) - so an accession with ONE
    valid row and one NaN row still ends up scored, not dropped. n_dropped
    counts DISTINCT accessions left with no score at all, not dropped rows:
    a row-level count would overstate the true number of un-scored proteins
    whenever a duplicate accession has some (or all) rows non-finite, since
    each such row would count separately even though it is the SAME protein.

    Returns (score_by_accession, n_dropped).
    """
    if metric not in GSEA_RANKING_METRICS:
        metric = QUANT_METRIC_DEFAULT
    df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""].copy()
    input_accessions = set(df["accession"].astype(str))

    if metric == "signed_log_padj":
        df = df[df["padj"].notna()]

    log2fc = df["log2FC"].to_numpy(dtype=float)
    # sign(0) is 0 in numpy, which would silently strand a protein with
    # exactly log2FC==0 - define the boundary case as +1 explicitly, same
    # convention as build_ranked_list.
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

    finite = np.isfinite(values)
    accs = df["accession"].astype(str).to_numpy()[finite]
    values = values[finite]

    score_by_acc = _mean_score_by_key(list(accs), values)
    n_dropped = len(input_accessions) - len(score_by_acc)
    return score_by_acc, n_dropped


def run_ora(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    taxon_name: str,
    padj_cutoff: float = PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = LFC_CUTOFF_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    include_iea: bool = True,
    force_refresh: bool = False,
    custom_db: "dict | None" = None,
    background_accessions: "set | None" = None,
    msigdb_collections: "list[str] | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run ORA for the requested databases against the resolved protein
    table. Mirrors Transcriptomics_module.run_ora's shape (per-database,
    per-direction independently-BH-corrected result families under
    "tables"), but built on Genomics_module's accession-native engine - see
    module docstring for why.

    databases: subset of PROTEOMICS_ORA_DATABASES.
    msigdb_collections: required (non-empty) when "MSigDB" is selected.
    custom_db: required when "Custom database" is selected - the dict
        returned by Genomics_module.load_custom_geneset_db (id_type=
        "accession" is that function's own default, matching this module's
        identifier space with no adaptation needed).
    background_accessions: optional UniProt-accession set restricting the
        universe for every database (from load_user_background_list).
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universe_df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""]
    universe_accessions_all = set(universe_df["accession"].astype(str))
    n_universe_before_background = len(universe_accessions_all)
    if background_accessions:
        universe_accessions_all = universe_accessions_all & {str(a).strip() for a in background_accessions}

    accession_to_symbol = dict(zip(
        universe_df["accession"].astype(str), universe_df["gene_symbol"].astype(str).str.upper(),
    ))

    directions = {}
    for direction_key in DIRECTIONS:
        sig_df = select_significant_proteins(resolved_df, padj_cutoff, lfc_cutoff, direction_key)
        raw_set = set(sig_df["accession"].astype(str))
        directions[direction_key] = raw_set & universe_accessions_all

    n_msigdb_selected = len(_msigdb_selected_keys(msigdb_collections))
    # GO increments `step` once per NAMESPACE x DIRECTION (3x3=9), and MSigDB
    # once per SELECTED COLLECTION x DIRECTION, not once per direction like
    # every other database - a plain len(databases)*len(DIRECTIONS) undercounts
    # both, so progress (step/steps_total) overshoots 1.0 and the bar pins at
    # 100% before the run actually finishes. Mirrors the identical, already-
    # documented calculation in Transcriptomics_module.run_ora.
    steps_total = max(
        sum(
            3 if db == "GO" else            # biological_process/molecular_function/cellular_component
            n_msigdb_selected if db == "MSigDB" else
            1
            for db in databases
        ) * len(DIRECTIONS),
        1,
    )
    step = 0
    tables: dict = {}

    def _empty_block(note: str) -> dict:
        return {
            "table": pd.DataFrame(columns=PROTEOMICS_ORA_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
            "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
            "note": note, "caveat": None,
        }

    def _drop_zero_overlap(term_to_genes: dict, sig_set: set, universe_genes: set):
        """Exclude terms with ZERO overlap with the significant-protein set
        from the family entirely, before Fisher/BH - clusterProfiler
        convention. Genomics_module._run_fisher_for_namespace does NOT do
        this itself (verified live 2026-08-24: reusing it unfiltered gave
        1291 "tested" Complex Portal terms for a 4-protein input, which is
        the exact k=0-inflation bug Transcriptomics_module._run_ora_for_sets'
        own docstring documents having fixed - "An earlier version of this
        function did [include k=0 rows]... which measurably worsens every
        other term's FDR without adding information"). That fix was never
        applied to Genomics' function. Filtering HERE, in this module,
        avoids silently reinheriting the bug without touching
        Genomics_module.py (out of scope for this work) - replicates the
        SAME min/max term-size gate _run_fisher_for_namespace applies
        internally so the two counts partition cleanly with no double count.
        Returns (filtered_term_to_genes, n_excluded_min_size, n_no_overlap).
        """
        filtered, n_excluded_min_size, n_no_overlap = {}, 0, 0
        for term_id, info in term_to_genes.items():
            term_genes = info["genes"] & universe_genes
            K = len(term_genes)
            if K < min_term_size or (max_term_size and K > max_term_size):
                n_excluded_min_size += 1
                continue
            if term_genes & sig_set:
                filtered[term_id] = info
            else:
                n_no_overlap += 1
        return filtered, n_excluded_min_size, n_no_overlap

    def _run_family(result_key: str, direction_key: str, term_to_genes: dict, universe_genes: set, caveat=None):
        nonlocal step
        step += 1
        _progress(step, steps_total, result_key)
        sig_set = directions[direction_key]
        if not sig_set:
            tables[result_key] = _empty_block("No significant proteins in this direction at the current thresholds.")
            return
        filtered_terms, n_excl, n_no_overlap = _drop_zero_overlap(term_to_genes, sig_set, universe_genes)
        # min_term_size/max_term_size already applied above - pass through
        # disabled (0 / None) so _run_fisher_for_namespace's own gate is a
        # no-op and n_excl above is the single source of truth for that count.
        block = _run_fisher_for_namespace(
            filtered_terms, sig_set, universe_genes, FDR_ALPHA_DEFAULT, True, 0, None,
        )
        block["n_terms_excluded_min_size"] = n_excl
        block["n_terms_no_overlap"] = n_no_overlap
        if caveat:
            block["caveat"] = caveat
        tables[result_key] = block

    if "GO" in databases:
        go_universe = build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
        for ns in ("biological_process", "molecular_function", "cellular_component"):
            ns_label = ns.replace("_", " ").title()
            ns_genes: set = set()
            for info in go_universe[ns].values():
                ns_genes |= info["genes"]
            # ALWAYS narrow to the proteins actually quantified in this
            # experiment (universe_accessions_all is already narrowed by any
            # background at the top of this function) - matches
            # Transcriptomics_module._run_ora_for_sets' "annotated_universe =
            # universe_symbols & annotated_genes" convention. A conditional
            # narrow (only when a background file was uploaded) was a real
            # bug: without it N was the WHOLE GO-annotated proteome (~150k
            # entries incl. isoforms/unreviewed) instead of the ~3-8k
            # proteins actually measured, inflating every p-value here by
            # many orders of magnitude - verified live 2026-08-24.
            ns_genes &= universe_accessions_all
            for direction_key in DIRECTIONS:
                _run_family(f"GO: {ns_label} ({direction_key})", direction_key, go_universe[ns], ns_genes, go_universe.get("caveat"))

    if "Reactome" in databases:
        if not taxon_name.strip():
            for direction_key in DIRECTIONS:
                step += 1
                tables[f"Reactome ({direction_key})"] = _empty_block(
                    "Organism name unavailable; Reactome's file cannot be filtered by species.")
        else:
            reactome_universe = build_reactome_universe(taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose)
            r_genes: set = set()
            for info in reactome_universe["pathways"].values():
                r_genes |= info["genes"]
            r_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
            for direction_key in DIRECTIONS:
                _run_family(f"Reactome ({direction_key})", direction_key, reactome_universe["pathways"], r_genes)

    if "UniProt (disease)" in databases:
        disease_universe = build_uniprot_disease_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        d_genes: set = set()
        for info in disease_universe["terms"].values():
            d_genes |= info["genes"]
        d_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
        for direction_key in DIRECTIONS:
            _run_family(f"UniProt disease ({direction_key})", direction_key, disease_universe["terms"], d_genes)

    if "UniProt (keywords)" in databases:
        kw_universe = build_uniprot_keyword_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        k_genes: set = set()
        for info in kw_universe["terms"].values():
            k_genes |= info["genes"]
        k_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
        for direction_key in DIRECTIONS:
            _run_family(f"UniProt keywords ({direction_key})", direction_key, kw_universe["terms"], k_genes)

    if "InterPro" in databases:
        ipr_universe = build_interpro_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        i_genes: set = set()
        for info in ipr_universe["terms"].values():
            i_genes |= info["genes"]
        i_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
        for direction_key in DIRECTIONS:
            _run_family(f"InterPro ({direction_key})", direction_key, ipr_universe["terms"], i_genes)

    if "Complex Portal" in databases:
        cpx_universe = build_complex_portal_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        c_genes: set = set()
        for info in cpx_universe["terms"].values():
            c_genes |= info["genes"]
        c_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
        for direction_key in DIRECTIONS:
            _run_family(f"Complex Portal ({direction_key})", direction_key, cpx_universe["terms"], c_genes, cpx_universe.get("caveat"))

    # Symbol-space bridge, shared by the two symbol-keyed databases here
    # (MSigDB and DISEASES): the same accession-space universe, mapped
    # through the resolved accession->gene_symbol pairs - mirrors Genomics'
    # own "GWAS Catalog is symbol-native inside an accession-native module"
    # pattern, not a new architectural choice. Computed once, lazily, only
    # if one of those databases was actually selected.
    symbol_universe_all: set = set()
    directions_symbol: dict = {}
    if "MSigDB" in databases or "DISEASES" in databases:
        symbol_universe_all = {accession_to_symbol[a] for a in universe_accessions_all if accession_to_symbol.get(a)}
        directions_symbol = {
            d: {accession_to_symbol[a] for a in directions[d] if accession_to_symbol.get(a)}
            for d in DIRECTIONS
        }

    if "DISEASES" in databases:
        dis_universe = build_diseases_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        # Narrow to symbols BOTH quantified here AND annotated in DISEASES -
        # the same ALWAYS-narrow rule every other database above follows,
        # applied in symbol space.
        dis_annotated: set = set()
        for info in dis_universe["terms"].values():
            dis_annotated |= info["genes"]
        dis_symbol_universe = symbol_universe_all & dis_annotated
        for direction_key in DIRECTIONS:
            step += 1
            _progress(step, steps_total, f"DISEASES ({direction_key})")
            result_key = f"DISEASES ({direction_key})"
            sig_symbols = directions_symbol[direction_key]
            if not dis_universe["terms"]:
                tables[result_key] = _empty_block(
                    dis_universe.get("caveat") or "No DISEASES data available for this run.")
                continue
            if not sig_symbols:
                tables[result_key] = _empty_block("No significant proteins in this direction at the current thresholds.")
                continue
            filtered_sets, n_excl, n_no_overlap = _drop_zero_overlap(
                dis_universe["terms"], sig_symbols, dis_symbol_universe)
            block = _run_fisher_for_namespace(
                filtered_sets, sig_symbols, dis_symbol_universe, FDR_ALPHA_DEFAULT, True, 0, None,
            )
            block["n_terms_excluded_min_size"] = n_excl
            block["n_terms_no_overlap"] = n_no_overlap
            if dis_universe.get("caveat"):
                block["caveat"] = dis_universe["caveat"]
            tables[result_key] = block

    if "MSigDB" in databases:
        selected = _msigdb_selected_keys(msigdb_collections)
        for coll_key in selected:
            coll_label = MSIGDB_COLLECTIONS.get(coll_key, {}).get("label", coll_key)
            if not msigdb_collection_available(coll_key, taxon_id):
                for direction_key in DIRECTIONS:
                    step += 1
                    tables[f"MSigDB: {coll_label} ({direction_key})"] = _empty_block(
                        f"MSigDB {coll_label} is not available for this organism.")
                continue
            fetched = fetch_msigdb_sets(coll_key, taxon_id, verbose=verbose, force_refresh=force_refresh)
            msigdb_sets = {
                name: {"name": info["name"], "genes": set(info["genes"])}
                for name, info in fetched.get("sets", {}).items()
            }
            # Narrow to symbols BOTH quantified in this experiment AND
            # annotated somewhere in this MSigDB collection - same ALWAYS-
            # narrow fix as every accession-space database above, applied in
            # symbol space. Per-collection (not computed once outside the
            # loop) because each collection annotates a different symbol set.
            msigdb_annotated_symbols: set = set()
            for info in msigdb_sets.values():
                msigdb_annotated_symbols |= info["genes"]
            symbol_universe = symbol_universe_all & msigdb_annotated_symbols
            for direction_key in DIRECTIONS:
                step += 1
                _progress(step, steps_total, f"MSigDB: {coll_label} ({direction_key})")
                sig_symbols = directions_symbol[direction_key]
                result_key = f"MSigDB: {coll_label} ({direction_key})"
                if fetched.get("fetch_failed"):
                    tables[result_key] = _empty_block(f"Could not fetch MSigDB {coll_label} for this run. Re-run to retry.")
                    continue
                if not sig_symbols:
                    tables[result_key] = _empty_block("No significant proteins in this direction at the current thresholds.")
                    continue
                filtered_sets, n_excl, n_no_overlap = _drop_zero_overlap(msigdb_sets, sig_symbols, symbol_universe)
                block = _run_fisher_for_namespace(
                    filtered_sets, sig_symbols, symbol_universe, FDR_ALPHA_DEFAULT, True, 0, None,
                )
                block["n_terms_excluded_min_size"] = n_excl
                block["n_terms_no_overlap"] = n_no_overlap
                tables[result_key] = block

    if "Custom database" in databases and custom_db:
        custom_genes: set = set()
        for info in custom_db["terms"].values():
            custom_genes |= info["genes"]
        custom_genes &= universe_accessions_all  # ALWAYS - see the GO block's comment above
        for direction_key in DIRECTIONS:
            _run_family(f"Custom database ({direction_key})", direction_key, custom_db["terms"], custom_genes)

    # --- Normalise the overlapping-member column across every database -----
    # Applied ONCE here, over every stored block, rather than at each of the
    # several places a block is created (_run_family, the MSigDB loop, the
    # DISEASES loop, _empty_block) - one pass cannot miss a path.
    #
    # Two problems this fixes, both real and user-reported:
    #   1. The shared engine names the column "Overlapping Genes"; in a
    #      proteomics module these are PROTEINS.
    #   2. Its contents silently changed identifier type per database -
    #      accessions from the 6 accession-keyed sources, bare gene symbols
    #      from MSigDB and DISEASES - so the same column meant different
    #      things on different rows of the same results screen.
    # Both are normalised to "ACCESSION (SYMBOL)", which keeps the
    # unambiguous proteomics identifier while staying readable (P05386 alone
    # tells a reader nothing; RPLP0 alone loses the isoform-level identity
    # that is the whole point of an accession).
    _finalize_member_column(tables, accession_to_symbol, "Overlapping Genes", PROTEOMICS_ORA_RESULT_COLUMNS)

    n_deg_before_background = {d: len(select_significant_proteins(resolved_df, padj_cutoff, lfc_cutoff, d)) for d in DIRECTIONS}
    summary = {
        "n_universe": len(universe_accessions_all),
        "n_universe_before_background": n_universe_before_background,
        "background_applied": bool(background_accessions),
        "n_sig_all": len(directions["all"]),
        "n_sig_up": len(directions["up"]),
        "n_sig_down": len(directions["down"]),
        "n_sig_all_before_background": n_deg_before_background["all"],
    }
    return {"tables": tables, "summary": summary}


# ---------------------------------------------------------------------------
# Quantitative Enrichment (1D annotation enrichment / competitive Wilcoxon)
# ---------------------------------------------------------------------------
# Why this exists as a SEPARATE analysis from ORA, and why it is not GSEA:
# ORA only sees proteins that already passed a padj/log2FC cutoff, so a
# coordinated but individually sub-threshold shift across a whole pathway is
# invisible to it. The usual answer is GSEA, but GSEA assumes a long, dense
# ranked list (~20k genes for RNA-seq); shotgun proteomics quantifies far
# fewer proteins (3-8k) AND its missing values are missing-NOT-at-random
# (undetected because low-abundance, not at random), biasing one end of any
# ranking. Perseus - proteomics' de facto standard tool - instead uses "1D
# annotation enrichment": for each term, a Mann-Whitney/Wilcoxon rank-sum
# test asking whether the log2FC values of proteins INSIDE the term are
# shifted relative to every other quantified protein. No thresholding, uses
# every quantified protein, no ranked-list assumption. That is exactly what
# Genomics' Gene-set Burden Test (_run_wilcoxon_for_namespace) already
# computes, extended here to alternative="two-sided" (see that function's
# docstring) because log2FC is signed - a downward shift is exactly as
# interesting as an upward one, and Genomics' own one-sided default would
# silently miss every down-shifted term.
#
# Structural differences from run_ora, both deliberate:
#   - No up/down/all split. The signed score already encodes direction and
#     the test is two-sided, so there is exactly ONE family per database
#     instead of three - a materially smaller multiple-testing burden, and
#     the honest structure (splitting would just be re-adding a threshold
#     through the back door).
#   - No padj/log2FC cutoff. Every quantified protein with a finite score
#     contributes; that is the entire point of a threshold-free test.

def run_quantitative_enrichment(
    resolved_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str,
    taxon_name: str,
    metric: str = QUANT_METRIC_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    include_iea: bool = True,
    force_refresh: bool = False,
    custom_db: "dict | None" = None,
    background_accessions: "set | None" = None,
    msigdb_collections: "list[str] | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Run Quantitative (1D annotation) Enrichment for the requested
    databases against the resolved protein table. See module comment above
    for why this differs from run_ora (Wilcoxon vs Fisher, one family per
    database instead of three, no significance cutoff).

    databases: subset of PROTEOMICS_ORA_DATABASES (same database set as ORA;
        PANTHER is excluded for the same reason - see module docstring).
    metric: one of GSEA_RANKING_METRICS - the per-protein score the Wilcoxon
        test compares inside vs. outside each term. Reuses Transcriptomics'
        GSEA metric definitions (signed_log_pvalue default, log2fc, pi_value,
        signed_log_padj) for cross-module consistency.
    msigdb_collections: required (non-empty) when "MSigDB" is selected.
    custom_db: required when "Custom database" is selected.
    background_accessions: optional UniProt-accession set restricting the
        universe for every database (from load_user_background_list).
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universe_df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""]
    universe_accessions_all = set(universe_df["accession"].astype(str))
    n_universe_before_background = len(universe_accessions_all)
    if background_accessions:
        universe_accessions_all = universe_accessions_all & {str(a).strip() for a in background_accessions}

    accession_to_symbol = dict(zip(
        universe_df["accession"].astype(str), universe_df["gene_symbol"].astype(str).str.upper(),
    ))

    score_by_acc_all, n_score_dropped = _build_protein_scores(universe_df, metric)
    # Scores restricted to the (possibly background-narrowed) universe, same
    # ALWAYS-narrow rule run_ora applies to every database's gene set below.
    score_by_acc = {a: s for a, s in score_by_acc_all.items() if a in universe_accessions_all}

    n_msigdb_selected = len(_msigdb_selected_keys(msigdb_collections))
    # One family per database (no x3 direction multiplier, unlike run_ora) -
    # GO still counts once per NAMESPACE (3), MSigDB once per SELECTED
    # COLLECTION, mirroring run_ora's identical, already-documented reasoning
    # for why a plain len(databases) undercounts both.
    steps_total = max(
        sum(
            3 if db == "GO" else
            n_msigdb_selected if db == "MSigDB" else
            1
            for db in databases
        ),
        1,
    )
    step = 0
    tables: dict = {}

    def _empty_block(note: str) -> dict:
        return {
            "table": pd.DataFrame(columns=PROTEOMICS_BURDEN_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
            "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
            "note": note, "caveat": None,
        }

    def _run_family(result_key: str, term_to_genes: dict, universe_genes: set, scores: dict, caveat=None):
        nonlocal step
        step += 1
        _progress(step, steps_total, result_key)
        if not scores:
            tables[result_key] = _empty_block("No quantified proteins with a valid score at the current metric.")
            return
        block = _run_wilcoxon_for_namespace(
            term_to_genes, scores, universe_genes, True, min_term_size, max_term_size,
            alternative="two-sided",
        )
        if caveat:
            block["caveat"] = caveat
        tables[result_key] = block

    if "GO" in databases:
        go_universe = build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose, force_refresh=force_refresh)
        for ns in ("biological_process", "molecular_function", "cellular_component"):
            ns_label = ns.replace("_", " ").title()
            ns_genes: set = set()
            for info in go_universe[ns].values():
                ns_genes |= info["genes"]
            ns_genes &= universe_accessions_all  # ALWAYS - see run_ora's identical comment
            _run_family(f"GO: {ns_label}", go_universe[ns], ns_genes, score_by_acc, go_universe.get("caveat"))

    if "Reactome" in databases:
        if not taxon_name.strip():
            step += 1
            tables["Reactome"] = _empty_block(
                "Organism name unavailable; Reactome's file cannot be filtered by species.")
        else:
            reactome_universe = build_reactome_universe(taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose)
            r_genes: set = set()
            for info in reactome_universe["pathways"].values():
                r_genes |= info["genes"]
            r_genes &= universe_accessions_all
            _run_family("Reactome", reactome_universe["pathways"], r_genes, score_by_acc)

    if "UniProt (disease)" in databases:
        disease_universe = build_uniprot_disease_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        d_genes: set = set()
        for info in disease_universe["terms"].values():
            d_genes |= info["genes"]
        d_genes &= universe_accessions_all
        _run_family("UniProt disease", disease_universe["terms"], d_genes, score_by_acc)

    if "UniProt (keywords)" in databases:
        kw_universe = build_uniprot_keyword_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        k_genes: set = set()
        for info in kw_universe["terms"].values():
            k_genes |= info["genes"]
        k_genes &= universe_accessions_all
        _run_family("UniProt keywords", kw_universe["terms"], k_genes, score_by_acc)

    if "InterPro" in databases:
        ipr_universe = build_interpro_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        i_genes: set = set()
        for info in ipr_universe["terms"].values():
            i_genes |= info["genes"]
        i_genes &= universe_accessions_all
        _run_family("InterPro", ipr_universe["terms"], i_genes, score_by_acc)

    if "Complex Portal" in databases:
        cpx_universe = build_complex_portal_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        c_genes: set = set()
        for info in cpx_universe["terms"].values():
            c_genes |= info["genes"]
        c_genes &= universe_accessions_all
        _run_family("Complex Portal", cpx_universe["terms"], c_genes, score_by_acc, cpx_universe.get("caveat"))

    # Symbol-space bridge for the 2 symbol-keyed databases (MSigDB,
    # DISEASES) - mirrors run_ora's identical bridge, but aggregates scores
    # (not just accession identity) into per-symbol scores via
    # _mean_score_by_key, since more than one accession can map to the same
    # symbol (see that helper's docstring).
    symbol_universe_all: set = set()
    score_by_symbol: dict = {}
    if "MSigDB" in databases or "DISEASES" in databases:
        symbol_universe_all = {accession_to_symbol[a] for a in universe_accessions_all if accession_to_symbol.get(a)}
        sym_keys, sym_vals = [], []
        for acc, val in score_by_acc.items():
            sym = accession_to_symbol.get(acc)
            if sym:
                sym_keys.append(sym)
                sym_vals.append(val)
        score_by_symbol = _mean_score_by_key(sym_keys, np.array(sym_vals, dtype=float)) if sym_keys else {}

    if "DISEASES" in databases:
        dis_universe = build_diseases_universe(taxon_id, verbose=verbose, force_refresh=force_refresh)
        dis_annotated: set = set()
        for info in dis_universe["terms"].values():
            dis_annotated |= info["genes"]
        dis_symbol_universe = symbol_universe_all & dis_annotated
        step += 1
        _progress(step, steps_total, "DISEASES")
        if not dis_universe["terms"]:
            tables["DISEASES"] = _empty_block(dis_universe.get("caveat") or "No DISEASES data available for this run.")
        elif not score_by_symbol:
            tables["DISEASES"] = _empty_block("No quantified proteins with a valid score at the current metric.")
        else:
            block = _run_wilcoxon_for_namespace(
                dis_universe["terms"], score_by_symbol, dis_symbol_universe, True, min_term_size, max_term_size,
                alternative="two-sided",
            )
            if dis_universe.get("caveat"):
                block["caveat"] = dis_universe["caveat"]
            tables["DISEASES"] = block

    if "MSigDB" in databases:
        selected = _msigdb_selected_keys(msigdb_collections)
        for coll_key in selected:
            coll_label = MSIGDB_COLLECTIONS.get(coll_key, {}).get("label", coll_key)
            result_key = f"MSigDB: {coll_label}"
            if not msigdb_collection_available(coll_key, taxon_id):
                step += 1
                tables[result_key] = _empty_block(f"MSigDB {coll_label} is not available for this organism.")
                continue
            fetched = fetch_msigdb_sets(coll_key, taxon_id, verbose=verbose, force_refresh=force_refresh)
            msigdb_sets = {
                name: {"name": info["name"], "genes": set(info["genes"])}
                for name, info in fetched.get("sets", {}).items()
            }
            msigdb_annotated_symbols: set = set()
            for info in msigdb_sets.values():
                msigdb_annotated_symbols |= info["genes"]
            symbol_universe = symbol_universe_all & msigdb_annotated_symbols
            step += 1
            _progress(step, steps_total, result_key)
            if fetched.get("fetch_failed"):
                tables[result_key] = _empty_block(f"Could not fetch MSigDB {coll_label} for this run. Re-run to retry.")
                continue
            if not score_by_symbol:
                tables[result_key] = _empty_block("No quantified proteins with a valid score at the current metric.")
                continue
            block = _run_wilcoxon_for_namespace(
                msigdb_sets, score_by_symbol, symbol_universe, True, min_term_size, max_term_size,
                alternative="two-sided",
            )
            tables[result_key] = block

    if "Custom database" in databases and custom_db:
        custom_genes: set = set()
        for info in custom_db["terms"].values():
            custom_genes |= info["genes"]
        custom_genes &= universe_accessions_all
        _run_family("Custom database", custom_db["terms"], custom_genes, score_by_acc)

    _finalize_member_column(tables, accession_to_symbol, "Genes", PROTEOMICS_BURDEN_RESULT_COLUMNS)

    summary = {
        "n_universe": len(universe_accessions_all),
        "n_universe_before_background": n_universe_before_background,
        "background_applied": bool(background_accessions),
        "n_scored": len(score_by_acc),
        "n_scored_dropped": n_score_dropped,
        "metric": metric if metric in GSEA_RANKING_METRICS else QUANT_METRIC_DEFAULT,
    }
    return {"tables": tables, "summary": summary}
