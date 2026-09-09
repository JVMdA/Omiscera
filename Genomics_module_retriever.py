# -*- coding: utf-8 -*-
"""
Genomics module gene identifier resolver for Omiscera.

Resolves gene symbols or Ensembl Gene IDs to UniProt accessions using the
UniProt REST API. Resolution happens in three sequential layers:

  Layer 1   : Query reviewed (Swiss-Prot) entries only.
  Layer 1b  : TrEMBL fallback when no reviewed entry exists.
  Layer 2   : Reference-proteome tie-break when multiple candidates remain.
  Layer 3   : Flag as ambiguous when still multiple — never silently picks one.

Resolution statuses
-------------------
  user_provided      UniProt ID was already in the input row; no lookup done.
  resolved           Exactly one reviewed (Swiss-Prot) accession found.
  resolved_unreviewed  No reviewed entry; one unreviewed (TrEMBL) accession found.
  resolved_via_ncbi  UniProt had nothing directly; bridged via an NCBI Entrez
                     Gene ID lookup (_resolve_via_geneid).
  resolved_via_kegg  Was ambiguous under UniProt alone; KEGG's own accession
                     map collapsed the candidates to exactly one gene (see
                     resolve_kegg_ids - runs only when a KEGG organism was
                     chosen, after this module's normal resolution).
  resolved_manual    User manually selected one accession from ambiguous candidates.
  ambiguous          Multiple candidates remain after tie-break; needs user choice.
  not_found          Zero results from UniProt (reviewed and unreviewed).
"""

import io
import json
import os
import re
import time

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

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PROTEOMES_URL = "https://rest.uniprot.org/proteomes/search"
UNIPROT_FIELDS = "accession,id,protein_name,gene_names,xref_proteomes,organism_id"

NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

# NCBI BLAST URL API (Option 3 — Homology path)
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
BLAST_DATABASE = "swissprot"
BLAST_POLL_INTERVAL = 60    # seconds between status polls (NCBI guideline: ≤1/min/RID)
BLAST_SUBMIT_DELAY = 5      # seconds between consecutive new submissions
BLAST_BATCH_DEADLINE = 180  # seconds to wait in-app before returning blast_pending rows
BLAST_RID_TTL = 86400       # seconds before a saved RID is considered expired (24 h)
BLAST_LOW_IDENTITY = 30.0   # percent identity below which Confidence = "low"
BLAST_HIGH_EVALUE = 1e-3    # e-value above which Confidence = "low"

OMISCERA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositorio do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "application/json",
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
GENOMICS_CACHE_PATH = os.path.join(CACHE_DIR, "genomics_uniprot_cache.csv")
GENOMICS_BLAST_CACHE_PATH = os.path.join(CACHE_DIR, "genomics_blast_cache.csv")
GENOMICS_BLAST_PENDING_PATH = os.path.join(CACHE_DIR, "genomics_blast_pending.csv")

DEFAULT_TIMEOUT = 20
DEFAULT_RETRIES = 3
MAX_RESULTS_PER_QUERY = 10  # enough to detect ambiguity; UniProt rarely returns >2-3 exact matches

# "KEGG ID" is filled by resolve_kegg_ids, a separate opt-in pass that runs
# AFTER resolve_gene_table (it needs a chosen KEGG organism code, which is
# independent of - and optional beyond - the taxon used for UniProt
# resolution). resolve_gene_table always writes it as "" so the column shape
# is stable whether or not that pass ever runs.
OUTPUT_COLUMNS_SYMBOL = [
    "Gene Symbol", "UniProt ID", "Entry Name", "Protein Name", "NCBI Gene ID",
    "KEGG ID", "Resolution_status", "Candidates",
]
OUTPUT_COLUMNS_ENSEMBL = [
    "Ensembl Gene ID", "UniProt ID", "Entry Name", "Protein Name", "NCBI Gene ID",
    "KEGG ID", "Resolution_status", "Candidates",
]
OUTPUT_COLUMNS_SEQUENCE = [
    "Sequence", "Sequence Type", "UniProt ID", "Entry Name", "Protein Name",
    "Percent Identity", "E-value", "Confidence", "Resolution_status",
]

# ---------------------------------------------------------------------------
# Cache management
# ---------------------------------------------------------------------------

GENOMICS_CACHE: dict = {}
_CACHE_LOADED = False


def _ensure_cache_loaded() -> None:
    global _CACHE_LOADED
    if _CACHE_LOADED:
        return
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(GENOMICS_CACHE_PATH):
        try:
            df = pd.read_csv(GENOMICS_CACHE_PATH, dtype=str).fillna("")
            for _, row in df.iterrows():
                cache_type = row.get("CacheType", "")
                cache_key = row.get("CacheKey", "")
                if cache_type and cache_key:
                    GENOMICS_CACHE[(cache_type, cache_key)] = {
                        "UniProtID": row.get("UniProtID", ""),
                        "EntryName": row.get("EntryName", ""),
                        "ProteinName": row.get("ProteinName", ""),
                        "Resolution_status": row.get("Resolution_status", ""),
                        "Candidates": row.get("Candidates", ""),
                        "NCBIGeneID": row.get("NCBIGeneID", ""),
                    }
        except Exception:
            pass
    _CACHE_LOADED = True


def save_genomics_cache() -> None:
    """Persist the in-memory gene→UniProt cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for (cache_type, cache_key), data in GENOMICS_CACHE.items():
        rows.append({
            "CacheType": cache_type,
            "CacheKey": cache_key,
            "UniProtID": data.get("UniProtID", ""),
            "EntryName": data.get("EntryName", ""),
            "ProteinName": data.get("ProteinName", ""),
            "Resolution_status": data.get("Resolution_status", ""),
            "Candidates": data.get("Candidates", ""),
            "NCBIGeneID": data.get("NCBIGeneID", ""),
        })
    pd.DataFrame(rows).to_csv(GENOMICS_CACHE_PATH, index=False)


def _get_cache(cache_type: str, cache_key: str) -> dict:
    _ensure_cache_loaded()
    return GENOMICS_CACHE.get((cache_type, cache_key), {})


def _set_cache(cache_type: str, cache_key: str, data: dict) -> None:
    _ensure_cache_loaded()
    # Only cache confident positive resolutions. Caching "ambiguous" would serve
    # stale answers after the user resolves them manually; caching "not_found"
    # would poison the cache when a query transiently fails or returns nothing.
    if data.get("Resolution_status") in ("resolved", "resolved_unreviewed", "resolved_via_ncbi"):
        GENOMICS_CACHE[(cache_type, cache_key)] = {k: data.get(k, "") for k in
            ("UniProtID", "EntryName", "ProteinName", "Resolution_status", "Candidates", "NCBIGeneID")}


# ---------------------------------------------------------------------------
# HTTP helper
# ---------------------------------------------------------------------------

def _safe_uniprot_search(query: str, verbose: bool = False) -> list:
    """Query the UniProt REST API and return the list of raw result dicts."""
    if HTTP_SESSION is None:
        if verbose:
            print("requests library not available; cannot search UniProt")
        return []

    params = {
        "query": query,
        "format": "json",
        "size": MAX_RESULTS_PER_QUERY,
        "fields": UNIPROT_FIELDS,
    }
    last_error = ""
    for attempt in range(DEFAULT_RETRIES + 1):
        try:
            response = HTTP_SESSION.get(
                UNIPROT_SEARCH_URL,
                params=params,
                headers=OMISCERA_HEADERS,
                timeout=DEFAULT_TIMEOUT,
            )
            if response.status_code in (200, 201):
                try:
                    return response.json().get("results", [])
                except Exception as exc:
                    if verbose:
                        print(f"JSON parse failed | query={query} | error={exc}")
                    return []
            if response.status_code == 404:
                return []
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
        print(f"UniProt search failed | query={query} | error={last_error}")
    return []


# ---------------------------------------------------------------------------
# Entry parsing helpers
# ---------------------------------------------------------------------------

def _parse_protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    if rec:
        return str(rec.get("fullName", {}).get("value", "")).strip()
    subs = desc.get("submissionNames", [])
    if subs:
        return str(subs[0].get("fullName", {}).get("value", "")).strip()
    return ""


def _parse_entry(entry: dict) -> dict:
    accession = str(entry.get("primaryAccession", "")).strip()
    entry_name = str(entry.get("uniProtkbId", "")).strip()
    protein_name = _parse_protein_name(entry)

    # Proteome IDs come from the "Proteomes" cross-references (requested via the
    # xref_proteomes field). The entry itself carries NO reference-proteome flag;
    # reference status is resolved separately via _get_reference_proteome_ids().
    proteome_ids = set()
    for xref in entry.get("uniProtKBCrossReferences", []) or []:
        if str(xref.get("database", "")) == "Proteomes":
            upid = str(xref.get("id", "")).strip()
            if upid:
                proteome_ids.add(upid)

    taxon_id = str(entry.get("organism", {}).get("taxonId", "")).strip()

    return {
        "accession": accession,
        "entry_name": entry_name,
        "protein_name": protein_name,
        "proteome_ids": proteome_ids,
        "taxon_id": taxon_id,
    }


# ---------------------------------------------------------------------------
# Reference-proteome lookup (per taxon, cached in-memory per run)
# ---------------------------------------------------------------------------

_REF_PROTEOME_CACHE: dict = {}

# Process-local cache for the "harvest an Entrez ID for an ambiguous gene"
# step below, keyed "symbol|taxon" (or bare ensembl_id). Mirrors
# _REF_PROTEOME_CACHE's pattern: _set_cache deliberately refuses to persist
# "ambiguous" results (a stale answer would survive the user's manual fix),
# so without this an ambiguous gene's NCBI esearch would re-fire on every
# single run against the same table, not just every process.
_NCBI_GENEID_CACHE: dict = {}


def _get_reference_proteome_ids(taxon_id: str, verbose: bool = False) -> set:
    """Return the set of reference-proteome UPIDs for a taxon.

    The UniProtKB entry response does not flag reference proteomes, so this
    queries the Proteomes endpoint and keeps entries whose proteomeType is
    "Reference proteome" (e.g. Homo sapiens taxon 9606 ->{"UP000005640"}).
    """
    taxon_clean = str(taxon_id).strip()
    if not taxon_clean:
        return set()
    if taxon_clean in _REF_PROTEOME_CACHE:
        return _REF_PROTEOME_CACHE[taxon_clean]

    ids: set = set()
    if HTTP_SESSION is not None:
        params = {"query": f"organism_id:{taxon_clean}", "format": "json", "size": 50}
        try:
            response = HTTP_SESSION.get(
                UNIPROT_PROTEOMES_URL,
                params=params,
                headers=OMISCERA_HEADERS,
                timeout=DEFAULT_TIMEOUT,
            )
            if response.status_code in (200, 201):
                for rec in response.json().get("results", []):
                    ptype = str(rec.get("proteomeType", "")).strip().lower()
                    if ptype.startswith("reference"):
                        upid = str(rec.get("id", "")).strip()
                        if upid:
                            ids.add(upid)
        except Exception as exc:
            if verbose:
                print(f"Reference proteome lookup failed | taxon={taxon_clean} | error={exc}")

    _REF_PROTEOME_CACHE[taxon_clean] = ids
    return ids


# ---------------------------------------------------------------------------
# Core selection logic (shared by symbol and Ensembl paths)
# ---------------------------------------------------------------------------

def _select_unambiguous(parsed_entries: list, resolved_status: str, taxon_for_ref: str = "") -> dict:
    """Select one entry unambiguously or flag the list as ambiguous.

    Applies the reference-proteome tie-break when >1 candidate exists.
    resolved_status is either "resolved" (reviewed) or "resolved_unreviewed".
    taxon_for_ref is the taxon used to look up the reference proteome; for the
    gene-symbol path it is the chosen organism, for the Ensembl path it is
    derived from the returned entries.
    """
    if not parsed_entries:
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "Resolution_status": "not_found", "Candidates": "", "NCBIGeneID": "",
        }

    if len(parsed_entries) == 1:
        e = parsed_entries[0]
        return {
            "UniProtID": e["accession"],
            "EntryName": e["entry_name"],
            "ProteinName": e["protein_name"],
            "Resolution_status": resolved_status,
            "Candidates": "", "NCBIGeneID": "",
        }

    # Reference-proteome tie-break: keep entries whose proteome cross-reference
    # is the organism's reference proteome.
    ref_taxon = taxon_for_ref or parsed_entries[0].get("taxon_id", "")
    reference_proteome_ids = _get_reference_proteome_ids(ref_taxon)
    ref_entries = [
        e for e in parsed_entries
        if reference_proteome_ids and (e["proteome_ids"] & reference_proteome_ids)
    ]
    # Use ref-proteome subset if it narrows the list; otherwise keep original.
    candidates = ref_entries if 1 <= len(ref_entries) < len(parsed_entries) else parsed_entries

    if len(candidates) == 1:
        e = candidates[0]
        return {
            "UniProtID": e["accession"],
            "EntryName": e["entry_name"],
            "ProteinName": e["protein_name"],
            "Resolution_status": resolved_status,
            "Candidates": "", "NCBIGeneID": "",
        }

    # Still ambiguous — surface all remaining candidates
    candidates_payload = [
        {"accession": e["accession"], "entry_name": e["entry_name"], "protein_name": e["protein_name"]}
        for e in candidates
    ]
    return {
        "UniProtID": "", "EntryName": "", "ProteinName": "",
        "Resolution_status": "ambiguous",
        "Candidates": json.dumps(candidates_payload), "NCBIGeneID": "",
    }


# ---------------------------------------------------------------------------
# Gene-symbol resolution
# ---------------------------------------------------------------------------

def resolve_gene_symbol(symbol: str, taxon_id: str, verbose: bool = False) -> dict:
    """Resolve a gene symbol + taxon ID to a single UniProt accession.

    Args:
        symbol:   Gene symbol as the user provided it (case is preserved).
        taxon_id: UniProt numeric taxon ID (e.g. "9606" for Homo sapiens).
        verbose:  Print resolution steps.

    Returns:
        Dict with keys UniProtID, EntryName, ProteinName, Resolution_status, Candidates.
    """
    symbol_clean = str(symbol).strip()
    taxon_clean = str(taxon_id).strip()
    if not symbol_clean or not taxon_clean:
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "Resolution_status": "not_found", "Candidates": "", "NCBIGeneID": "",
        }

    cache_key = f"{symbol_clean.lower()}|{taxon_clean}"
    cached = _get_cache("gene_symbol", cache_key)
    if cached and cached.get("Resolution_status"):
        if verbose:
            print(f"[cache hit] gene_symbol | {cache_key} ->{cached['Resolution_status']}")
        return cached

    # Layer 1: reviewed Swiss-Prot entries only
    query_reviewed = (
        f"gene_exact:{symbol_clean} AND organism_id:{taxon_clean} AND reviewed:true"
    )
    entries = _safe_uniprot_search(query_reviewed, verbose=verbose)
    if verbose:
        print(f"[L1 reviewed] {symbol_clean}|{taxon_clean} ->{len(entries)} entries")

    if entries:
        parsed = [_parse_entry(e) for e in entries]
        result = _select_unambiguous(parsed, "resolved", taxon_for_ref=taxon_clean)
    else:
        # Layer 1b: TrEMBL fallback
        query_unreviewed = (
            f"gene_exact:{symbol_clean} AND organism_id:{taxon_clean} AND reviewed:false"
        )
        entries_ur = _safe_uniprot_search(query_unreviewed, verbose=verbose)
        if verbose:
            print(f"[L1b unreviewed] {symbol_clean}|{taxon_clean} ->{len(entries_ur)} entries")
        parsed = [_parse_entry(e) for e in entries_ur]
        result = _select_unambiguous(parsed, "resolved_unreviewed", taxon_for_ref=taxon_clean)

    if result["Resolution_status"] == "not_found":
        result = _ncbi_fallback_symbol(symbol_clean, taxon_clean, verbose=verbose)
    elif result["Resolution_status"] == "ambiguous" and not result.get("NCBIGeneID"):
        # Harvest-only: an Entrez ID here does not resolve the gene (status and
        # Candidates are untouched), but it enables the KEGG bridge in
        # build_recovery_bundle / the retriever's KEGG ID column, and costs
        # nothing extra since _harvest_ncbi_geneid_for_ambiguous is cached.
        harvested = _harvest_ncbi_geneid_for_ambiguous(
            f"{symbol_clean}[Gene Name] AND txid{taxon_clean}[Organism:exp]",
            cache_key, verbose=verbose,
        )
        if harvested:
            result["NCBIGeneID"] = harvested
    _set_cache("gene_symbol", cache_key, result)
    return result


# ---------------------------------------------------------------------------
# Ensembl Gene ID resolution
# ---------------------------------------------------------------------------

def _ensembl_xref_clause(ensembl_id: str) -> str:
    """Build the UniProt xref query clause for an Ensembl Gene ID.

    UniProt cross-reference queries use the form `xref:{database}-{id}`.
    Vertebrate Ensembl gene IDs (ENSG, ENSMUSG, ENSDARG, …) live under the
    "ensembl" database; Ensembl Genomes divisions (plants/fungi/metazoa/
    protists/bacteria) use their own database names, so those are OR-ed.
    """
    if re.match(r"^ENS[A-Z]*G\d+", ensembl_id, re.IGNORECASE):
        return f"xref:ensembl-{ensembl_id}"
    divisions = ["ensemblplants", "ensemblfungi", "ensemblmetazoa", "ensemblprotists", "ensemblbacteria"]
    return "(" + " OR ".join(f"xref:{db}-{ensembl_id}" for db in divisions) + ")"


def resolve_ensembl_id(ensembl_id: str, verbose: bool = False) -> dict:
    """Resolve an Ensembl Gene ID to a single UniProt accession.

    No organism selector is needed because Ensembl IDs are species-specific;
    the reference-proteome tie-break derives the taxon from the matched entries.

    Args:
        ensembl_id: Ensembl Gene ID (e.g. "ENSG00000141510").
        verbose:    Print resolution steps.

    Returns:
        Dict with keys UniProtID, EntryName, ProteinName, Resolution_status, Candidates.
    """
    ensembl_clean = str(ensembl_id).strip()
    if not ensembl_clean:
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "Resolution_status": "not_found", "Candidates": "", "NCBIGeneID": "",
        }

    cache_key = ensembl_clean.lower()
    cached = _get_cache("ensembl", cache_key)
    if cached and cached.get("Resolution_status"):
        if verbose:
            print(f"[cache hit] ensembl | {cache_key} ->{cached['Resolution_status']}")
        return cached

    xref_clause = _ensembl_xref_clause(ensembl_clean)

    # Layer 1: reviewed Swiss-Prot entries only
    query_reviewed = f"{xref_clause} AND reviewed:true"
    entries = _safe_uniprot_search(query_reviewed, verbose=verbose)
    if verbose:
        print(f"[L1 reviewed] ensembl {ensembl_clean} ->{len(entries)} entries")

    if entries:
        parsed = [_parse_entry(e) for e in entries]
        result = _select_unambiguous(parsed, "resolved")
    else:
        # Layer 1b: TrEMBL fallback
        query_unreviewed = f"{xref_clause} AND reviewed:false"
        entries_ur = _safe_uniprot_search(query_unreviewed, verbose=verbose)
        if verbose:
            print(f"[L1b unreviewed] ensembl {ensembl_clean} ->{len(entries_ur)} entries")
        parsed = [_parse_entry(e) for e in entries_ur]
        result = _select_unambiguous(parsed, "resolved_unreviewed")

    if result["Resolution_status"] == "not_found":
        result = _ncbi_fallback_ensembl(ensembl_clean, verbose=verbose)
    elif result["Resolution_status"] == "ambiguous" and not result.get("NCBIGeneID"):
        harvested = _harvest_ncbi_geneid_for_ambiguous(ensembl_clean, cache_key, verbose=verbose)
        if harvested:
            result["NCBIGeneID"] = harvested
    _set_cache("ensembl", cache_key, result)
    return result


# ---------------------------------------------------------------------------
# NCBI Gene fallback (Layer 5)
# ---------------------------------------------------------------------------

def _ncbi_gene_search(term: str, verbose: bool = False) -> list:
    """Query NCBI Gene esearch and return a list of Entrez Gene IDs (strings)."""
    if HTTP_SESSION is None:
        return []
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "retmax": 5,
        "tool": "Omiscera",
        "email": "claudeiajoao@gmail.com",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    # Respect NCBI's rate limit (≤3 req/s without key); Layer 5 fires only for
    # already-failed genes so the per-table volume is small.
    time.sleep(0.34)
    try:
        response = HTTP_SESSION.get(
            NCBI_ESEARCH_URL,
            params=params,
            headers=OMISCERA_HEADERS,
            timeout=DEFAULT_TIMEOUT,
        )
        if response.status_code in (200, 201):
            return response.json().get("esearchresult", {}).get("idlist", [])
    except Exception as exc:
        if verbose:
            print(f"NCBI Gene search failed | term={term} | error={exc}")
    return []


def _harvest_ncbi_geneid_for_ambiguous(term: str, cache_key: str, verbose: bool = False) -> str:
    """Best-effort Entrez ID lookup for an AMBIGUOUS gene - harvest only, never override.

    Returns the Entrez ID string, or "" if there was no single unambiguous
    esearch hit. Deliberately stricter than _ncbi_fallback_symbol/_ensembl
    (which take gene_ids[0] of up to 5): those fire only after UniProt found
    NOTHING, so any hit is better than none, but an ambiguous gene is
    ambiguous precisely because real paralogs likely exist - taking [0] here
    could inject an arbitrary paralog's ID. Requiring exactly one hit makes
    this bridge silent rather than wrong.

    The caller's Resolution_status and Candidates are never touched - this
    only ever adds an NCBIGeneID onto an otherwise-unchanged ambiguous result.
    """
    if cache_key in _NCBI_GENEID_CACHE:
        return _NCBI_GENEID_CACHE[cache_key]
    ids = _ncbi_gene_search(term, verbose=verbose)
    entrez_id = ids[0] if len(ids) == 1 else ""
    _NCBI_GENEID_CACHE[cache_key] = entrez_id
    return entrez_id


def _resolve_via_geneid(entrez_id: str, taxon_for_ref: str = "", verbose: bool = False) -> dict:
    """Bridge an Entrez Gene ID to a UniProt accession via the geneid cross-reference.

    Queries UniProt `xref:geneid-{id}` (reviewed first, then TrEMBL) and runs
    the same reference-proteome tie-break as the direct symbol/Ensembl paths.
    On a single hit the status is overridden to resolved_via_ncbi.

    The Entrez ID itself is now kept UNCONDITIONALLY, even when the UniProt
    bridge fails or stays ambiguous - it is a real, independently-useful piece
    of information the caller already paid an NCBI query for (KEGG gene IDs
    ARE Entrez IDs for RefSeq-sourced organisms, so this is what lets a gene
    with no UniProt match still be matched against KEGG). Previously it was
    thrown away whenever the bridge did not resolve cleanly.
    """
    query_rev = f"xref:geneid-{entrez_id} AND reviewed:true"
    entries = _safe_uniprot_search(query_rev, verbose=verbose)
    if verbose:
        print(f"[L5 bridge reviewed] geneid={entrez_id} ->{len(entries)} entries")

    if entries:
        parsed = [_parse_entry(e) for e in entries]
        result = _select_unambiguous(parsed, "resolved", taxon_for_ref=taxon_for_ref)
    else:
        query_ur = f"xref:geneid-{entrez_id} AND reviewed:false"
        entries_ur = _safe_uniprot_search(query_ur, verbose=verbose)
        if verbose:
            print(f"[L5 bridge unreviewed] geneid={entrez_id} ->{len(entries_ur)} entries")
        parsed = [_parse_entry(e) for e in entries_ur]
        result = _select_unambiguous(parsed, "resolved_unreviewed", taxon_for_ref=taxon_for_ref)

    if result["Resolution_status"] in ("resolved", "resolved_unreviewed"):
        result["Resolution_status"] = "resolved_via_ncbi"
    result["NCBIGeneID"] = entrez_id
    return result


def _ncbi_fallback_symbol(symbol: str, taxon_id: str, verbose: bool = False) -> dict:
    """Layer 5 NCBI Gene fallback for gene symbol mode."""
    term = f"{symbol}[Gene Name] AND txid{taxon_id}[Organism:exp]"
    if verbose:
        print(f"[L5 NCBI symbol] term={term}")
    gene_ids = _ncbi_gene_search(term, verbose=verbose)
    if not gene_ids:
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "Resolution_status": "not_found", "Candidates": "", "NCBIGeneID": "",
        }
    entrez_id = gene_ids[0]
    if verbose:
        print(f"[L5 NCBI symbol] {symbol} ->Entrez Gene ID {entrez_id}")
    return _resolve_via_geneid(entrez_id, taxon_for_ref=taxon_id, verbose=verbose)


def _ncbi_fallback_ensembl(ensembl_id: str, verbose: bool = False) -> dict:
    """Layer 5 NCBI Gene fallback for Ensembl Gene ID mode."""
    if verbose:
        print(f"[L5 NCBI ensembl] term={ensembl_id}")
    gene_ids = _ncbi_gene_search(ensembl_id, verbose=verbose)
    if not gene_ids:
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "Resolution_status": "not_found", "Candidates": "", "NCBIGeneID": "",
        }
    entrez_id = gene_ids[0]
    if verbose:
        print(f"[L5 NCBI ensembl] {ensembl_id} ->Entrez Gene ID {entrez_id}")
    # taxon_for_ref left empty: _resolve_via_geneid derives it from the matched entries
    return _resolve_via_geneid(entrez_id, taxon_for_ref="", verbose=verbose)


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def _norm_gene_key(value) -> str:
    """Normalise a gene identifier for duplicate comparison (case-insensitive)."""
    text = re.sub(r"\s+", " ", str(value).strip()).lower()
    if text in ("nan", "none"):
        return ""
    return text


def deduplicate_gene_table(df: pd.DataFrame, identifier_type: str) -> tuple:
    """Collapse repeated gene identifiers to a single row before resolution.

    Duplicate genes would be counted multiple times and bias downstream
    overrepresentation statistics, so each gene must appear exactly once.
    When the same identifier appears more than once, the most informative
    occurrence is kept (a row carrying a user-provided UniProt ID wins over
    one without; otherwise the first occurrence is kept). Rows with an empty
    identifier are left untouched (not collapsed together).

    Returns:
        (deduplicated_df, stats) where stats has keys
        rows_before, rows_after, duplicates_removed.
    """
    id_col = "Gene Symbol" if identifier_type == "gene_symbol" else "Ensembl Gene ID"
    rows_before = len(df)
    if id_col not in df.columns:
        return df.copy(), {"rows_before": rows_before, "rows_after": rows_before, "duplicates_removed": 0}

    has_uniprot = "UniProt ID" in df.columns
    seen: dict = {}   # normalised key -> position in kept list
    kept: list = []

    for _, row in df.iterrows():
        key = _norm_gene_key(row.get(id_col, ""))
        if not key:
            kept.append(row)
            continue
        if key not in seen:
            seen[key] = len(kept)
            kept.append(row)
        elif has_uniprot:
            # Duplicate: upgrade the kept row only if it lacks a UniProt ID and
            # this occurrence supplies one.
            existing = kept[seen[key]]
            existing_up = str(existing.get("UniProt ID", "")).strip()
            new_up = str(row.get("UniProt ID", "")).strip()
            if not existing_up and new_up and new_up.lower() not in ("nan", "none"):
                kept[seen[key]] = row

    deduped = pd.DataFrame(kept, columns=df.columns).reset_index(drop=True)
    stats = {
        "rows_before": rows_before,
        "rows_after": len(deduped),
        "duplicates_removed": rows_before - len(deduped),
    }
    return deduped, stats


# ---------------------------------------------------------------------------
# Table-level resolution
# ---------------------------------------------------------------------------

def resolve_gene_table(
    df: pd.DataFrame,
    identifier_type: str,
    taxon_id: "str | None" = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Resolve a full gene identifier table to UniProt accessions.

    Args:
        df:              DataFrame with the uploaded gene table.
        identifier_type: "gene_symbol" or "ensembl_gene_id".
        taxon_id:        Required for "gene_symbol" mode (UniProt numeric taxon ID).
        verbose:         Print per-row resolution details.

    Returns:
        DataFrame with columns matching OUTPUT_COLUMNS_SYMBOL or OUTPUT_COLUMNS_ENSEMBL.
    """
    _ensure_cache_loaded()

    if identifier_type == "gene_symbol":
        id_col = "Gene Symbol"
        output_columns = OUTPUT_COLUMNS_SYMBOL
    else:
        id_col = "Ensembl Gene ID"
        output_columns = OUTPUT_COLUMNS_ENSEMBL

    # Defensive dedup: harmless (0 removed) if the caller already deduplicated.
    df, _ = deduplicate_gene_table(df, identifier_type)

    rows = []
    total = len(df)
    for idx, (_, row) in enumerate(df.iterrows(), start=1):
        gene_id = str(row.get(id_col, "")).strip()
        uniprot_provided = str(row.get("UniProt ID", "")).strip()

        if verbose:
            print(f"Processing {idx}/{total}: {gene_id}")

        # Skip resolution when the user already provided a UniProt ID
        if uniprot_provided and uniprot_provided.lower() not in ("nan", "none", ""):
            rows.append({
                id_col: gene_id,
                "UniProt ID": uniprot_provided,
                "Entry Name": "",
                "Protein Name": "",
                "NCBI Gene ID": "",
                "KEGG ID": "",
                "Resolution_status": "user_provided",
                "Candidates": "",
            })
            continue

        if not gene_id or gene_id.lower() in ("nan", "none", ""):
            rows.append({
                id_col: "",
                "UniProt ID": "",
                "Entry Name": "",
                "Protein Name": "",
                "NCBI Gene ID": "",
                "KEGG ID": "",
                "Resolution_status": "not_found",
                "Candidates": "",
            })
            continue

        if identifier_type == "gene_symbol":
            result = resolve_gene_symbol(gene_id, taxon_id or "", verbose=verbose)
        else:
            result = resolve_ensembl_id(gene_id, verbose=verbose)

        rows.append({
            id_col: gene_id,
            "UniProt ID": result.get("UniProtID", ""),
            "Entry Name": result.get("EntryName", ""),
            "Protein Name": result.get("ProteinName", ""),
            "NCBI Gene ID": result.get("NCBIGeneID", ""),
            "KEGG ID": "",
            "Resolution_status": result.get("Resolution_status", "not_found"),
            "Candidates": result.get("Candidates", ""),
        })

    save_genomics_cache()
    return pd.DataFrame(rows, columns=output_columns)


def resolve_kegg_ids(
    resolved_df: pd.DataFrame,
    id_col: str,
    kegg_organism_code: str,
    verbose: bool = False,
) -> pd.DataFrame:
    """Assign a KEGG gene ID to each row of an already-resolved gene table,
    and use KEGG's own accession map to disambiguate rows UniProt alone could
    not (native KEGG source, via Genomics_module.build_kegg_gene_map).

    Must run AFTER resolve_gene_table - it reads "UniProt ID", "NCBI Gene ID",
    id_col and, for ambiguous rows, "Candidates". The KEGG organism is a
    separate, optional choice from the UniProt taxon (KEGG organism codes and
    NCBI taxon IDs cannot be derived from one another - verified live).

    Precedence per gene, first hit wins (measured live on chickpea, a real
    non-model organism, for calibration):
      1. UniProt ID -> KEGG's acc_to_gene (~92% coverage).
      2. NCBI Gene ID -> KEGG's entrez_to_gene (KEGG gene IDs ARE Entrez IDs
         for RefSeq-sourced organisms - an exact, organism-agnostic bridge;
         this is what makes a gene UniProt alone could not place still
         reachable, since _resolve_via_geneid/the ambiguous-harvest step
         above now KEEP the Entrez ID even when the UniProt bridge fails).
      3. Gene symbol -> KEGG's symbol_to_gene (model-organism only - measured
         1.6% coverage for chickpea vs 94.1% for human; kept as a low-value
         last resort, never the primary route).

    AMBIGUOUS rows (empty UniProt ID, JSON candidates in "Candidates") get
    special handling: each candidate accession is checked against
    acc_to_gene. If the candidates collapse to EXACTLY ONE KEGG gene, that is
    not just a KEGG match - it is KEGG's own curated data disambiguating
    which UniProt accession is real for this organism, something the
    reference-proteome tie-break could not do (no reference proteome exists
    for most non-model organisms). The row is promoted: UniProt ID/Entry
    Name/Protein Name are filled from that one candidate, Resolution_status
    becomes "resolved_via_kegg", and Candidates is cleared - so the row flows
    into the normal resolved table through the EXISTING ambiguous/table1_df
    split with no new UI or manual-selection path required. More than one
    matching candidate -> left ambiguous, exactly as before; KEGG never guesses.

    Returns a NEW DataFrame (does not mutate resolved_df) with one added
    column, "KEGG ID" (already present as "" from resolve_gene_table; filled
    in place here).
    """
    import Genomics_module as gm  # local import: keeps this module importable if Genomics_module's heavier deps are unavailable

    df = resolved_df.copy()
    if "KEGG ID" not in df.columns:
        df["KEGG ID"] = ""

    org = str(kegg_organism_code).strip().lower()
    if not org:
        return df

    gene_map = gm.build_kegg_gene_map(org, verbose=verbose)
    acc_to_gene = gene_map.get("acc_to_gene", {})
    entrez_to_gene = gene_map.get("entrez_to_gene", {})
    symbol_to_gene = gene_map.get("symbol_to_gene", {})
    try_symbol_bridge = id_col == "Gene Symbol"

    for idx, row in df.iterrows():
        acc = str(row.get("UniProt ID", "")).strip()
        entrez = str(row.get("NCBI Gene ID", "")).strip()
        status = str(row.get("Resolution_status", "")).strip()
        was_ambiguous = status == "ambiguous"

        # An ambiguous row always has an empty UniProt ID by construction
        # (_select_unambiguous), so the accession bridge can only ever fire
        # here for an ALREADY-resolved row; the entrez and symbol bridges,
        # unlike it, can fire for an ambiguous row too (an ambiguous gene can
        # still carry a harvested NCBI Gene ID, or a gene symbol regardless
        # of status) - which is exactly the case the promotion step below
        # exists to handle.
        kegg_gene = ""
        if acc and acc in acc_to_gene:
            kegg_gene = acc_to_gene[acc]
        elif entrez and entrez in entrez_to_gene:
            kegg_gene = entrez_to_gene[entrez]
        elif try_symbol_bridge:
            symbol = str(row.get(id_col, "")).strip().upper()
            if symbol and symbol in symbol_to_gene:
                kegg_gene = symbol_to_gene[symbol]

        if was_ambiguous and not kegg_gene:
            # No direct bridge fired - fall back to KEGG's own accession map
            # disambiguating the pending candidates. Unlike the direct
            # bridges, this ALSO resolves the UniProt ambiguity itself (it
            # picks which of several real candidate accessions is correct),
            # so it fills UniProt ID/Entry Name/Protein Name too, not only
            # KEGG ID.
            candidates_raw = str(row.get("Candidates", "") or "")
            try:
                candidates = json.loads(candidates_raw) if candidates_raw else []
            except Exception:
                candidates = []

            matched_genes: set = set()
            matched_candidate = None
            for cand in candidates:
                cand_acc = str(cand.get("accession", "")).strip()
                g = acc_to_gene.get(cand_acc)
                if g:
                    matched_genes.add(g)
                    if len(matched_genes) == 1:
                        matched_candidate = cand

            if len(matched_genes) == 1 and matched_candidate is not None:
                df.at[idx, "UniProt ID"] = matched_candidate.get("accession", "")
                df.at[idx, "Entry Name"] = matched_candidate.get("entry_name", "")
                df.at[idx, "Protein Name"] = matched_candidate.get("protein_name", "")
                kegg_gene = next(iter(matched_genes))

        if not kegg_gene:
            continue

        df.at[idx, "KEGG ID"] = kegg_gene
        if was_ambiguous:
            # Whatever bridge fired - the entrez bridge directly, or the
            # accession-candidate disambiguation above - this row must be
            # promoted OUT of "ambiguous" or it would never reach final_df
            # at all (table1_df excludes every ambiguous row; a KEGG ID with
            # no status change would be silently unreachable). When the
            # entrez bridge is what fired, UniProt ID/Entry Name/Protein Name
            # stay empty - KEGG placed the gene but UniProt still could not -
            # which is fine: GO/Reactome/other UniProt-space databases
            # already skip rows with an empty accession, so this promotion's
            # only real effect for them is a no-op, while KEGG (and KEGG
            # alone) can now see the gene.
            df.at[idx, "Resolution_status"] = "resolved_via_kegg"
            df.at[idx, "Candidates"] = ""

    return df


# ---------------------------------------------------------------------------
# Input loading helpers
# ---------------------------------------------------------------------------

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
    )


def _parse_gene_csv_text(text: str) -> tuple:
    """Parse CSV text, normalise column names, detect identifier type.

    Returns:
        (DataFrame, identifier_type) where identifier_type is
        "gene_symbol" or "ensembl_gene_id".
    Raises:
        ValueError if no required column is found after normalisation.
    """
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

    df = best_df.copy()

    rename_map = {}
    for col in df.columns:
        key = _normalize_col_key(col)
        if key in ("genesymbol", "symbol", "gene", "hgncsymbol", "genesymbolid"):
            rename_map[col] = "Gene Symbol"
        elif key in ("ensemblgeneid", "ensemblid", "ensgid", "ensembl", "ensemblgene",
                     "ensembl_gene_id", "ensemblid"):
            rename_map[col] = "Ensembl Gene ID"
        elif key in ("uniprotid", "uniprot", "uniprotaccession", "accession",
                     "uniprotkbid", "entry"):
            rename_map[col] = "UniProt ID"

    df = df.rename(columns=rename_map)

    if "Gene Symbol" in df.columns:
        identifier_type = "gene_symbol"
        if "UniProt ID" not in df.columns:
            df["UniProt ID"] = ""
        keep = ["Gene Symbol", "UniProt ID"]
        extra = [c for c in df.columns if c not in keep]
        return df[keep + extra].copy(), identifier_type

    if "Ensembl Gene ID" in df.columns:
        identifier_type = "ensembl_gene_id"
        if "UniProt ID" not in df.columns:
            df["UniProt ID"] = ""
        keep = ["Ensembl Gene ID", "UniProt ID"]
        extra = [c for c in df.columns if c not in keep]
        return df[keep + extra].copy(), identifier_type

    raise ValueError(
        "Table must contain a 'Gene Symbol' or 'Ensembl Gene ID' column. "
        f"Detected columns after normalisation: {list(df.columns)}"
    )


def load_gene_table(path: str) -> tuple:
    """Load a gene table CSV from a file path.

    Returns:
        (DataFrame, identifier_type)
    """
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin1") as fh:
            text = fh.read()
    return _parse_gene_csv_text(text)


def load_gene_table_from_bytes(content: bytes) -> tuple:
    """Load a gene table CSV from raw bytes (e.g. from a Streamlit file uploader).

    Returns:
        (DataFrame, identifier_type)
    """
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("utf-8", errors="replace")
    return _parse_gene_csv_text(text)


# ---------------------------------------------------------------------------
# BLAST cache (separate from the gene cache — different schema, different path)
# ---------------------------------------------------------------------------

BLAST_CACHE: dict = {}
_BLAST_CACHE_LOADED = False


def _ensure_blast_cache_loaded() -> None:
    global _BLAST_CACHE_LOADED
    if _BLAST_CACHE_LOADED:
        return
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(GENOMICS_BLAST_CACHE_PATH):
        try:
            df = pd.read_csv(GENOMICS_BLAST_CACHE_PATH, dtype=str).fillna("")
            for _, row in df.iterrows():
                key = row.get("CacheKey", "")
                if key:
                    BLAST_CACHE[key] = {
                        "UniProtID": row.get("UniProtID", ""),
                        "EntryName": row.get("EntryName", ""),
                        "ProteinName": row.get("ProteinName", ""),
                        "PercentIdentity": row.get("PercentIdentity", ""),
                        "Evalue": row.get("Evalue", ""),
                        "Confidence": row.get("Confidence", ""),
                        "Resolution_status": row.get("Resolution_status", ""),
                    }
        except Exception:
            pass
    _BLAST_CACHE_LOADED = True


def save_blast_cache() -> None:
    """Persist the in-memory BLAST result cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for key, data in BLAST_CACHE.items():
        rows.append({
            "CacheKey": key,
            "UniProtID": data.get("UniProtID", ""),
            "EntryName": data.get("EntryName", ""),
            "ProteinName": data.get("ProteinName", ""),
            "PercentIdentity": data.get("PercentIdentity", ""),
            "Evalue": data.get("Evalue", ""),
            "Confidence": data.get("Confidence", ""),
            "Resolution_status": data.get("Resolution_status", ""),
        })
    pd.DataFrame(rows).to_csv(GENOMICS_BLAST_CACHE_PATH, index=False)


def _blast_cache_key(seq_type: str, sequence: str) -> str:
    import hashlib
    return f"{seq_type.lower()}|{hashlib.sha1(sequence.upper().encode()).hexdigest()}"


def _get_blast_cache(key: str) -> dict:
    _ensure_blast_cache_loaded()
    return BLAST_CACHE.get(key, {})


def _set_blast_cache(key: str, data: dict) -> None:
    _ensure_blast_cache_loaded()
    if data.get("Resolution_status") == "resolved_homology":
        BLAST_CACHE[key] = {k: data.get(k, "") for k in
            ("UniProtID", "EntryName", "ProteinName", "PercentIdentity", "Evalue", "Confidence", "Resolution_status")}


# ---------------------------------------------------------------------------
# BLAST pending-jobs store (RIDs saved across requests for resumability)
# ---------------------------------------------------------------------------

BLAST_PENDING: dict = {}  # cache_key -> {rid, program, seq_type, submit_epoch}
_BLAST_PENDING_LOADED = False


def _ensure_blast_pending_loaded() -> None:
    global _BLAST_PENDING_LOADED
    if _BLAST_PENDING_LOADED:
        return
    os.makedirs(CACHE_DIR, exist_ok=True)
    if os.path.exists(GENOMICS_BLAST_PENDING_PATH):
        try:
            df = pd.read_csv(GENOMICS_BLAST_PENDING_PATH, dtype=str).fillna("")
            now = time.time()
            for _, row in df.iterrows():
                key = row.get("CacheKey", "")
                epoch = float(row.get("SubmitEpoch", "0") or "0")
                if key and (now - epoch) < BLAST_RID_TTL:
                    BLAST_PENDING[key] = {
                        "rid": row.get("RID", ""),
                        "program": row.get("Program", ""),
                        "seq_type": row.get("SeqType", ""),
                        "submit_epoch": epoch,
                    }
        except Exception:
            pass
    _BLAST_PENDING_LOADED = True


def save_blast_pending() -> None:
    """Persist the in-memory pending-jobs store to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for key, data in BLAST_PENDING.items():
        rows.append({
            "CacheKey": key,
            "RID": data.get("rid", ""),
            "Program": data.get("program", ""),
            "SeqType": data.get("seq_type", ""),
            "SubmitEpoch": str(data.get("submit_epoch", "")),
        })
    pd.DataFrame(rows).to_csv(GENOMICS_BLAST_PENDING_PATH, index=False)


def _get_pending(key: str) -> dict:
    _ensure_blast_pending_loaded()
    entry = BLAST_PENDING.get(key, {})
    if entry:
        age = time.time() - float(entry.get("submit_epoch", 0))
        if age >= BLAST_RID_TTL:
            BLAST_PENDING.pop(key, None)
            return {}
    return entry


def _set_pending(key: str, rid: str, program: str, seq_type: str) -> None:
    _ensure_blast_pending_loaded()
    BLAST_PENDING[key] = {
        "rid": rid,
        "program": program,
        "seq_type": seq_type,
        "submit_epoch": time.time(),
    }


def _clear_pending(key: str) -> None:
    _ensure_blast_pending_loaded()
    BLAST_PENDING.pop(key, None)


# ---------------------------------------------------------------------------
# BLAST HTTP helpers
# ---------------------------------------------------------------------------

def _blast_submit(sequence: str, program: str, verbose: bool = False) -> "tuple[str, int]":
    """Submit a BLAST search. Returns (RID, rtoe) or ("", 0) on failure.

    rtoe is NCBI's estimated time-to-result in seconds (used to seed first poll
    delay instead of a fixed constant).
    """
    if HTTP_SESSION is None:
        return "", 0
    data = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": BLAST_DATABASE,
        "QUERY": sequence,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "1",
        "EXPECT": "10",
        "tool": "Omiscera",
        "email": "claudeiajoao@gmail.com",
    }
    try:
        response = HTTP_SESSION.post(
            NCBI_BLAST_URL,
            data=data,
            headers={"User-Agent": OMISCERA_HEADERS["User-Agent"]},
            timeout=60,
        )
        if response.status_code in (200, 201):
            text = response.text
            block = text
            m = re.search(r"QBlastInfoBegin(.*?)QBlastInfoEnd", text, re.DOTALL)
            if m:
                block = m.group(1)
            rid_match = re.search(r"\bRID\b\s*=\s*(\w+)", block)
            if rid_match:
                rid = rid_match.group(1)
                rtoe = 0
                rtoe_match = re.search(r"\bRTOE\b\s*=\s*(\d+)", block)
                if rtoe_match:
                    rtoe = int(rtoe_match.group(1))
                if verbose:
                    print(f"[BLAST submit] RID={rid} RTOE={rtoe}s")
                return rid, rtoe
    except Exception as exc:
        if verbose:
            print(f"[BLAST submit] error: {exc}")
    return "", 0


def _blast_poll(rid: str, verbose: bool = False) -> str:
    """Poll a BLAST RID for status. Returns 'READY', 'WAITING', or 'FAILED'."""
    if HTTP_SESSION is None:
        return "FAILED"
    params = {
        "CMD": "Get",
        "FORMAT_OBJECT": "SearchInfo",
        "RID": rid,
        "tool": "Omiscera",
        "email": "claudeiajoao@gmail.com",
    }
    try:
        response = HTTP_SESSION.get(
            NCBI_BLAST_URL, params=params,
            headers={"User-Agent": OMISCERA_HEADERS["User-Agent"]},
            timeout=30,
        )
        if response.status_code in (200, 201):
            text = response.text
            block = text
            m = re.search(r"QBlastInfoBegin(.*?)QBlastInfoEnd", text, re.DOTALL)
            if m:
                block = m.group(1)
            if "Status=READY" in block:
                return "READY"
            if "Status=WAITING" in block:
                return "WAITING"
            if "Status=UNKNOWN" in block or "Status=FAILED" in block:
                return "FAILED"
    except Exception as exc:
        if verbose:
            print(f"[BLAST poll] RID={rid} error: {exc}")
    return "WAITING"


def _blast_fetch_xml(rid: str, verbose: bool = False) -> str:
    """Retrieve BLAST results as XML text."""
    if HTTP_SESSION is None:
        return ""
    params = {
        "CMD": "Get",
        "FORMAT_TYPE": "XML",
        "RID": rid,
        "tool": "Omiscera",
        "email": "claudeiajoao@gmail.com",
    }
    try:
        response = HTTP_SESSION.get(
            NCBI_BLAST_URL, params=params,
            headers={"User-Agent": OMISCERA_HEADERS["User-Agent"]},
            timeout=60,
        )
        if response.status_code in (200, 201):
            return response.text
    except Exception as exc:
        if verbose:
            print(f"[BLAST fetch] RID={rid} error: {exc}")
    return ""


def _parse_blast_xml(xml_text: str) -> "dict | None":
    """Parse NCBI BLAST XML and return the best hit, or None if no hits.

    Returns dict with keys: UniProtID, EntryName, ProteinName,
    PercentIdentity (float), Evalue (float).
    """
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None

    # Walk to the first Hit/Hsp
    hits = root.findall(".//Hit")
    if not hits:
        return None
    hit = hits[0]

    hit_id = hit.findtext("Hit_id", "")
    hit_def = hit.findtext("Hit_def", "")
    hit_accession = hit.findtext("Hit_accession", "")

    # swissprot Hit_id format: "sp|P04637|P53_HUMAN" or "sp|ACC|ENTRY Protein Name"
    uniprot_id = ""
    entry_name = ""
    if hit_id.startswith("sp|"):
        parts = hit_id.split("|")
        if len(parts) >= 2:
            uniprot_id = re.sub(r"\.\d+$", "", parts[1].strip())
        if len(parts) >= 3 and parts[2].strip():
            entry_name = parts[2].split()[0].strip()
    if not uniprot_id:
        uniprot_id = re.sub(r"\.\d+$", "", hit_accession.strip())

    protein_name = hit_def.split("[")[0].strip() if hit_def else ""

    hsp = hit.find(".//Hsp")
    if hsp is None:
        return None

    try:
        identity = int(hsp.findtext("Hsp_identity", "0"))
        align_len = int(hsp.findtext("Hsp_align-len", "1"))
        percent_identity = round(identity / max(align_len, 1) * 100, 2)
    except (ValueError, ZeroDivisionError):
        percent_identity = 0.0

    try:
        evalue = float(hsp.findtext("Hsp_evalue", "999"))
    except ValueError:
        evalue = 999.0

    return {
        "UniProtID": uniprot_id,
        "EntryName": entry_name,
        "ProteinName": protein_name,
        "PercentIdentity": percent_identity,
        "Evalue": evalue,
    }


# ---------------------------------------------------------------------------
# Sequence (BLAST homology) resolution
# ---------------------------------------------------------------------------

def _normalise_seq_type(raw: str) -> str:
    """Normalise a Sequence Type cell to 'dna' or 'protein' (or original if unknown)."""
    key = raw.strip().lower().replace("-", "").replace(" ", "")
    if key in ("dna", "nucleotide", "nt", "nucl"):
        return "dna"
    if key in ("protein", "aa", "prot", "peptide", "aminoacid", "aminoacids"):
        return "protein"
    return raw.strip()


def resolve_sequence(sequence: str, seq_type: str, verbose: bool = False) -> dict:
    """Resolve a raw sequence to a UniProt accession via NCBI BLAST (swissprot).

    Args:
        sequence: DNA or protein sequence string.
        seq_type: Normalised type — "dna" (blastx) or "protein" (blastp).
        verbose:  Print resolution steps.

    Returns:
        Dict with keys UniProtID, EntryName, ProteinName, PercentIdentity,
        Evalue, Confidence, Resolution_status.
    """
    failed_result = {
        "UniProtID": "", "EntryName": "", "ProteinName": "",
        "PercentIdentity": "", "Evalue": "", "Confidence": "",
        "Resolution_status": "blast_failed",
    }

    seq_clean = re.sub(r"\s+", "", sequence).upper()
    if not seq_clean:
        return failed_result

    norm_type = _normalise_seq_type(seq_type)
    if norm_type not in ("dna", "protein"):
        result = failed_result.copy()
        result["Resolution_status"] = "blast_failed"
        return result

    cache_key = _blast_cache_key(norm_type, seq_clean)
    cached = _get_blast_cache(cache_key)
    if cached and cached.get("Resolution_status"):
        if verbose:
            print(f"[BLAST cache hit] {norm_type} -> {cached['Resolution_status']}")
        return cached

    program = "blastx" if norm_type == "dna" else "blastp"

    # Reuse a saved pending RID if it's still within TTL
    pending = _get_pending(cache_key)
    if pending and pending.get("rid"):
        rid = pending["rid"]
        if verbose:
            print(f"[BLAST resume] reusing RID={rid}")
    else:
        if verbose:
            print(f"[BLAST submit] program={program} seq_len={len(seq_clean)}")
        rid, rtoe = _blast_submit(seq_clean, program, verbose=verbose)
        if not rid:
            return failed_result
        _set_pending(cache_key, rid, program, norm_type)
        save_blast_pending()
        # Seed the first wait from NCBI's own estimate (clamped to a sane range)
        first_wait = max(BLAST_SUBMIT_DELAY, min(rtoe, BLAST_POLL_INTERVAL))
        time.sleep(first_wait)

    result = _blast_collect_ready(cache_key, rid, deadline=BLAST_BATCH_DEADLINE, verbose=verbose)
    return result


def _blast_collect_ready(
    cache_key: str,
    rid: str,
    deadline: int = BLAST_BATCH_DEADLINE,
    verbose: bool = False,
) -> dict:
    """Poll a single RID until READY (within deadline), then fetch and parse.

    On success: writes to result cache, clears pending entry, returns resolved dict.
    On timeout:  returns blast_pending dict (RID remains in pending store).
    On FAILED:   clears pending, returns blast_failed dict.
    """
    failed_result = {
        "UniProtID": "", "EntryName": "", "ProteinName": "",
        "PercentIdentity": "", "Evalue": "", "Confidence": "",
        "Resolution_status": "blast_failed",
    }
    pending_result = {
        "UniProtID": "", "EntryName": "", "ProteinName": "",
        "PercentIdentity": "", "Evalue": "", "Confidence": "",
        "Resolution_status": "blast_pending",
    }

    elapsed = 0
    status = "WAITING"
    while status == "WAITING" and elapsed < deadline:
        status = _blast_poll(rid, verbose=verbose)
        if verbose:
            print(f"[BLAST poll] RID={rid} status={status} elapsed={elapsed}s")
        if status == "WAITING":
            time.sleep(BLAST_POLL_INTERVAL)
            elapsed += BLAST_POLL_INTERVAL

    if status == "FAILED":
        _clear_pending(cache_key)
        return failed_result

    if status != "READY":
        # Deadline exceeded — leave pending entry intact for Resume
        if verbose:
            print(f"[BLAST] RID={rid} deadline exceeded, returning blast_pending")
        return pending_result

    xml_text = _blast_fetch_xml(rid, verbose=verbose)
    _clear_pending(cache_key)
    if not xml_text:
        return failed_result

    parsed = _parse_blast_xml(xml_text)
    if parsed is None or not parsed.get("UniProtID"):
        return {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "PercentIdentity": "", "Evalue": "", "Confidence": "",
            "Resolution_status": "not_found",
        }

    pct_id = parsed["PercentIdentity"]
    evalue = parsed["Evalue"]
    confidence = "low" if (pct_id < BLAST_LOW_IDENTITY or evalue > BLAST_HIGH_EVALUE) else "high"

    entry_name = parsed["EntryName"]
    protein_name = parsed["ProteinName"]
    try:
        enrich = _safe_uniprot_search(f"accession:{parsed['UniProtID']}", verbose=False)
        if enrich:
            ep = _parse_entry(enrich[0])
            entry_name = ep.get("entry_name", "") or entry_name
            protein_name = ep.get("protein_name", "") or protein_name
    except Exception:
        pass

    result = {
        "UniProtID": parsed["UniProtID"],
        "EntryName": entry_name,
        "ProteinName": protein_name,
        "PercentIdentity": str(pct_id),
        "Evalue": str(evalue),
        "Confidence": confidence,
        "Resolution_status": "resolved_homology",
    }
    _set_blast_cache(cache_key, result)
    return result


def resolve_sequence_table(
    df: pd.DataFrame,
    taxon_id: "str | None" = None,
    progress_callback=None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Resolve a sequence table to UniProt accessions via BLAST.

    Uses a submit-all-then-poll-in-parallel strategy:
      Phase A — submit every unresolved sequence; persist RIDs immediately.
      Phase B — poll all outstanding RIDs round-robin until all READY or deadline.
      Phase C — any still-outstanding rows become blast_pending (resumable).

    Args:
        df:                DataFrame with 'Sequence' and 'Sequence Type' columns.
        taxon_id:          Accepted but NOT used as a query filter.
        progress_callback: Optional callable(current, total, label) for UI updates.
        verbose:           Print per-row resolution details.

    Returns:
        DataFrame with OUTPUT_COLUMNS_SEQUENCE columns.
    """
    _ = taxon_id
    _ensure_blast_cache_loaded()
    _ensure_blast_pending_loaded()
    df, _ = deduplicate_sequence_table(df)

    total = len(df)
    resolved_count = 0

    def _progress(label: str) -> None:
        if progress_callback:
            try:
                progress_callback(resolved_count, total, label)
            except Exception:
                pass

    # Build an ordered list of row descriptors
    row_descs = []
    for _, row in df.iterrows():
        sequence = str(row.get("Sequence", "")).strip()
        seq_type = str(row.get("Sequence Type", "")).strip()
        uniprot_provided = str(row.get("UniProt ID", "")).strip()
        row_descs.append({
            "sequence": sequence,
            "seq_type": seq_type,
            "uniprot_provided": uniprot_provided,
            "result": None,  # filled below
        })

    # -----------------------------------------------------------------------
    # Phase 0: immediate rows (user_provided, blank, result-cache hits)
    # -----------------------------------------------------------------------
    need_blast = []  # indices into row_descs that require BLAST
    for i, rd in enumerate(row_descs):
        seq = rd["sequence"]
        seq_type = rd["seq_type"]
        up = rd["uniprot_provided"]

        if up and up.lower() not in ("nan", "none", ""):
            rd["result"] = {
                "UniProtID": up, "EntryName": "", "ProteinName": "",
                "PercentIdentity": "", "Evalue": "", "Confidence": "",
                "Resolution_status": "user_provided",
            }
            continue

        if not seq or seq.lower() in ("nan", "none"):
            rd["result"] = {
                "UniProtID": "", "EntryName": "", "ProteinName": "",
                "PercentIdentity": "", "Evalue": "", "Confidence": "",
                "Resolution_status": "blast_failed",
            }
            continue

        norm_type = _normalise_seq_type(seq_type)
        if norm_type not in ("dna", "protein"):
            rd["result"] = {
                "UniProtID": "", "EntryName": "", "ProteinName": "",
                "PercentIdentity": "", "Evalue": "", "Confidence": "",
                "Resolution_status": "blast_failed",
            }
            continue

        seq_clean = re.sub(r"\s+", "", seq).upper()
        cache_key = _blast_cache_key(norm_type, seq_clean)
        cached = _get_blast_cache(cache_key)
        if cached and cached.get("Resolution_status"):
            if verbose:
                print(f"[BLAST cache hit] row {i} -> {cached['Resolution_status']}")
            rd["result"] = cached
            continue

        need_blast.append((i, seq_clean, norm_type, cache_key))

    resolved_count = total - len(need_blast)
    _progress(f"Submitted 0/{len(need_blast)}")

    if not need_blast:
        save_blast_cache()
        save_blast_pending()
        return _build_sequence_output(row_descs)

    # -----------------------------------------------------------------------
    # Phase A: submit all (or reuse saved pending RIDs)
    # -----------------------------------------------------------------------
    outstanding = {}  # cache_key -> {idx, seq_clean, norm_type, rid}
    for submit_num, (i, seq_clean, norm_type, cache_key) in enumerate(need_blast):
        program = "blastx" if norm_type == "dna" else "blastp"
        pending = _get_pending(cache_key)
        if pending and pending.get("rid"):
            rid = pending["rid"]
            if verbose:
                print(f"[BLAST resume] row {i} RID={rid}")
        else:
            rid, rtoe = _blast_submit(seq_clean, program, verbose=verbose)
            if not rid:
                row_descs[i]["result"] = {
                    "UniProtID": "", "EntryName": "", "ProteinName": "",
                    "PercentIdentity": "", "Evalue": "", "Confidence": "",
                    "Resolution_status": "blast_failed",
                }
                resolved_count += 1
                _progress(f"Submitted {submit_num + 1}/{len(need_blast)}")
                continue
            _set_pending(cache_key, rid, program, norm_type)
            # Space out new submissions to be polite
            if submit_num < len(need_blast) - 1:
                time.sleep(BLAST_SUBMIT_DELAY)

        outstanding[cache_key] = {"idx": i, "rid": rid, "norm_type": norm_type}
        _progress(f"Submitted {submit_num + 1}/{len(need_blast)}")

    save_blast_pending()

    # -----------------------------------------------------------------------
    # Phase B: poll all outstanding RIDs in round-robin until deadline
    # -----------------------------------------------------------------------
    deadline_start = time.time()
    while outstanding and (time.time() - deadline_start) < BLAST_BATCH_DEADLINE:
        finished_keys = []
        for cache_key, info in list(outstanding.items()):
            rid = info["rid"]
            idx = info["idx"]
            status = _blast_poll(rid, verbose=verbose)
            if verbose:
                elapsed = int(time.time() - deadline_start)
                print(f"[BLAST poll] row {idx} RID={rid} status={status} batch_elapsed={elapsed}s")

            if status == "READY":
                xml_text = _blast_fetch_xml(rid, verbose=verbose)
                _clear_pending(cache_key)
                if not xml_text:
                    row_descs[idx]["result"] = {
                        "UniProtID": "", "EntryName": "", "ProteinName": "",
                        "PercentIdentity": "", "Evalue": "", "Confidence": "",
                        "Resolution_status": "blast_failed",
                    }
                else:
                    parsed = _parse_blast_xml(xml_text)
                    if parsed is None or not parsed.get("UniProtID"):
                        row_descs[idx]["result"] = {
                            "UniProtID": "", "EntryName": "", "ProteinName": "",
                            "PercentIdentity": "", "Evalue": "", "Confidence": "",
                            "Resolution_status": "not_found",
                        }
                    else:
                        pct_id = parsed["PercentIdentity"]
                        evalue = parsed["Evalue"]
                        confidence = (
                            "low" if (pct_id < BLAST_LOW_IDENTITY or evalue > BLAST_HIGH_EVALUE)
                            else "high"
                        )
                        entry_name = parsed["EntryName"]
                        protein_name = parsed["ProteinName"]
                        try:
                            enrich = _safe_uniprot_search(
                                f"accession:{parsed['UniProtID']}", verbose=False)
                            if enrich:
                                ep = _parse_entry(enrich[0])
                                entry_name = ep.get("entry_name", "") or entry_name
                                protein_name = ep.get("protein_name", "") or protein_name
                        except Exception:
                            pass
                        result = {
                            "UniProtID": parsed["UniProtID"],
                            "EntryName": entry_name,
                            "ProteinName": protein_name,
                            "PercentIdentity": str(pct_id),
                            "Evalue": str(evalue),
                            "Confidence": confidence,
                            "Resolution_status": "resolved_homology",
                        }
                        _set_blast_cache(cache_key, result)
                        row_descs[idx]["result"] = result

                finished_keys.append(cache_key)
                resolved_count += 1
                _progress(f"Resolved {resolved_count}/{total}")

            elif status == "FAILED":
                _clear_pending(cache_key)
                row_descs[idx]["result"] = {
                    "UniProtID": "", "EntryName": "", "ProteinName": "",
                    "PercentIdentity": "", "Evalue": "", "Confidence": "",
                    "Resolution_status": "blast_failed",
                }
                finished_keys.append(cache_key)
                resolved_count += 1
                _progress(f"Resolved {resolved_count}/{total}")

        for k in finished_keys:
            outstanding.pop(k, None)

        if outstanding:
            time.sleep(BLAST_POLL_INTERVAL)

    # -----------------------------------------------------------------------
    # Phase C: deadline — remaining outstanding rows become blast_pending
    # -----------------------------------------------------------------------
    for cache_key, info in outstanding.items():
        idx = info["idx"]
        if verbose:
            print(f"[BLAST] row {idx} deadline exceeded -> blast_pending")
        row_descs[idx]["result"] = {
            "UniProtID": "", "EntryName": "", "ProteinName": "",
            "PercentIdentity": "", "Evalue": "", "Confidence": "",
            "Resolution_status": "blast_pending",
        }

    save_blast_cache()
    save_blast_pending()

    if progress_callback:
        try:
            progress_callback(total, total, "Done")
        except Exception:
            pass

    return _build_sequence_output(row_descs)


def _build_sequence_output(row_descs: list) -> pd.DataFrame:
    """Convert resolved row descriptors to the standard OUTPUT_COLUMNS_SEQUENCE DataFrame."""
    rows = []
    for rd in row_descs:
        result = rd.get("result") or {}
        rows.append({
            "Sequence": rd["sequence"],
            "Sequence Type": rd["seq_type"],
            "UniProt ID": result.get("UniProtID", ""),
            "Entry Name": result.get("EntryName", ""),
            "Protein Name": result.get("ProteinName", ""),
            "Percent Identity": result.get("PercentIdentity", ""),
            "E-value": result.get("Evalue", ""),
            "Confidence": result.get("Confidence", ""),
            "Resolution_status": result.get("Resolution_status", "blast_failed"),
        })
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS_SEQUENCE)


# ---------------------------------------------------------------------------
# Sequence table deduplication
# ---------------------------------------------------------------------------

def deduplicate_sequence_table(df: pd.DataFrame) -> tuple:
    """Collapse rows with identical (seq_type, sequence) before BLAST resolution.

    The most informative row (one with a user-supplied UniProt ID) is kept;
    otherwise the first occurrence is kept.

    Returns:
        (deduplicated_df, stats) where stats has keys
        rows_before, rows_after, duplicates_removed.
    """
    rows_before = len(df)
    if "Sequence" not in df.columns or "Sequence Type" not in df.columns:
        return df.copy(), {"rows_before": rows_before, "rows_after": rows_before, "duplicates_removed": 0}

    has_uniprot = "UniProt ID" in df.columns
    seen: dict = {}  # normalised key -> position in kept list
    kept: list = []

    for _, row in df.iterrows():
        seq = re.sub(r"\s+", "", str(row.get("Sequence", ""))).upper()
        stype = _normalise_seq_type(str(row.get("Sequence Type", "")))
        key = f"{stype}|{seq}"
        if not seq:
            kept.append(row)
            continue
        if key not in seen:
            seen[key] = len(kept)
            kept.append(row)
        elif has_uniprot:
            existing = kept[seen[key]]
            existing_up = str(existing.get("UniProt ID", "")).strip()
            new_up = str(row.get("UniProt ID", "")).strip()
            if not existing_up and new_up and new_up.lower() not in ("nan", "none"):
                kept[seen[key]] = row

    deduped = pd.DataFrame(kept, columns=df.columns).reset_index(drop=True)
    stats = {
        "rows_before": rows_before,
        "rows_after": len(deduped),
        "duplicates_removed": rows_before - len(deduped),
    }
    return deduped, stats


# ---------------------------------------------------------------------------
# Sequence input loading helpers
# ---------------------------------------------------------------------------

def _parse_sequence_csv_text(text: str) -> tuple:
    """Parse CSV text, normalise column names for the sequence/BLAST path.

    Returns:
        (DataFrame, "sequence")
    Raises:
        ValueError if required columns are missing.
    """
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

    df = best_df.copy()
    rename_map = {}
    for col in df.columns:
        key = _normalize_col_key(col)
        if key in ("sequence", "seq", "nucleotide", "aminoacid", "dna", "protein"):
            rename_map[col] = "Sequence"
        elif key in ("sequencetype", "seqtype", "type", "moleculetype", "sequencetype"):
            rename_map[col] = "Sequence Type"
        elif key in ("uniprotid", "uniprot", "uniprotaccession", "accession",
                     "uniprotkbid", "entry"):
            rename_map[col] = "UniProt ID"

    df = df.rename(columns=rename_map)

    missing = [c for c in ("Sequence", "Sequence Type") if c not in df.columns]
    if missing:
        raise ValueError(
            f"Table must contain columns: {missing}. "
            f"Detected columns after normalisation: {list(df.columns)}"
        )

    if "UniProt ID" not in df.columns:
        df["UniProt ID"] = ""

    # Normalise Sequence Type synonyms in-place
    if "Sequence Type" in df.columns:
        df["Sequence Type"] = df["Sequence Type"].apply(_normalise_seq_type)

    keep = ["Sequence", "Sequence Type", "UniProt ID"]
    extra = [c for c in df.columns if c not in keep]
    return df[keep + extra].copy(), "sequence"


def load_sequence_table_from_bytes(content: bytes) -> tuple:
    """Load a sequence table CSV from raw bytes (Streamlit file uploader).

    Returns:
        (DataFrame, "sequence")
    """
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("utf-8", errors="replace")
    return _parse_sequence_csv_text(text)


def load_sequence_table(path: str) -> tuple:
    """Load a sequence table CSV from a file path.

    Returns:
        (DataFrame, "sequence")
    """
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin1") as fh:
            text = fh.read()
    return _parse_sequence_csv_text(text)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(
            "Usage: python Genomics_module_retriever.py <input.csv> [taxon_id]\n"
            "  taxon_id is required when the input contains a 'Gene Symbol' column."
        )
        sys.exit(1)

    input_path = sys.argv[1]
    taxon_arg = sys.argv[2] if len(sys.argv) > 2 else None

    input_df, id_type = load_gene_table(input_path)
    print(f"Detected identifier type: {id_type}")
    print(f"Loaded {len(input_df)} rows.")

    result_df = resolve_gene_table(input_df, id_type, taxon_id=taxon_arg, verbose=True)
    print(result_df.to_string(index=False))

    out_path = input_path.replace(".csv", "_resolved.csv")
    result_df.to_csv(out_path, index=False)
    print(f"\nSaved to {out_path}")
