# -*- coding: utf-8 -*-
"""Database catalog and online context utilities for BioRank.

This module provides database metadata and lightweight online context
utilities for the Metabolic Connection workflow, including LOTUS, HMDB,
Rhea + UniProt, and NPAtlas.

It defines:

- supported databases;
- required compound identifiers for each database;
- organism, species, taxon, or study metadata context requirements;
- online context retrieval helpers when safely available;
- standardized lookup result formats.

This file does not perform the full BioRank rank-based analysis.
The rank logic, Table 1 generation, and Table 2 pathway summarization
will be handled by a separate analysis module or by the Streamlit
workflow using the functions provided here.

This file should not download or persist large database files.
No database search should run automatically on import. Any test or
manual execution code must be placed under:

if __name__ == "__main__":
"""

import hashlib
import io
import json
import os
import re
import time
from typing import Any, Dict, Optional, cast

import pandas as pd  # type: ignore[reportMissingModuleSource]
import requests

ENABLE_ONLINE_APIS = os.environ.get("BIORANK_ENABLE_ONLINE_APIS", "1") == "1"
DEFAULT_TIMEOUT = 20
DEFAULT_RETRIES = 3

_MEMORY_CACHE: Dict[str, Any] = {}

BROWSER_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml,application/json,text/plain;q=0.9,*/*;q=0.8",
}

# Used exclusively for Wikidata API and SPARQL calls, per Wikimedia User-Agent Policy.
# Fill in the project URL/repo and contact e-mail before any public deployment.
WIKIDATA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositório do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml,application/json,text/plain;q=0.9,*/*;q=0.8",
}

# Self-identifying User-Agent for Rhea, UniProt, and UniChem calls.
# Fill in the project URL/repo and contact e-mail before any public deployment.
OMISCERA_HEADERS = {
    "User-Agent": (
        "Omiscera/1.0 (<PREENCHER: URL ou repositório do projeto>; "
        "contato: <PREENCHER: seu e-mail>) Python-requests"
    ),
    "Accept": "application/json",
}

RHEA_SPARQL_URL = "https://sparql.rhea-db.org/sparql"
UNIPROT_SPARQL_URL = "https://sparql.uniprot.org/sparql"
UNIPROT_TAXONOMY_URL = "https://rest.uniprot.org/taxonomy"
UNICHEM_INCHIKEY_URL = "https://www.ebi.ac.uk/unichem/rest/inchikey"
NPATLAS_API_URL = "https://www.npatlas.org/api/v1"

SUPPORTED_DATABASES = [
    "LOTUS",
    "HMDB",
    "Rhea + UniProt",
    "NPAtlas",
]

_CANONICAL_NAMES = {
    "lotus": "LOTUS",
    "hmdb": "HMDB",
    "hmdb import": "HMDB",
    "rhea + uniprot": "Rhea + UniProt",
    "rhea+uniprot": "Rhea + UniProt",
    "rhea uniprot": "Rhea + UniProt",
    "npatlas": "NPAtlas",
    "np atlas": "NPAtlas",
    "natural products atlas": "NPAtlas",
}


def get_supported_databases() -> list:
    """Return the supported database names in canonical order."""
    return list(SUPPORTED_DATABASES)


def canonical_database_name(database_name: str) -> str:
    """Normalize a database name to the canonical supported name."""
    if database_name is None:
        raise ValueError("database_name is required")
    key = normalize_text(database_name).lower()
    if key in _CANONICAL_NAMES:
        return _CANONICAL_NAMES[key]
    raise ValueError(f"Unsupported database: {database_name}")


def get_database_metadata() -> Dict[str, Dict[str, Any]]:
    """Return metadata describing each supported database."""
    return {
        "LOTUS": {
            "display_name": "LOTUS",
            "required_identifiers": ["SMILES", "InChIKey"],
            "preferred_identifier_order": ["SMILES", "InChIKey"],
            "context_type": "taxon_or_species",
            "context_label": "Species, taxon, or organism",
            "online_supported": True,
            "local_data_supported": False,
            "notes": "LOTUS compound lookup checks Wikidata occurrence (P703) by InChIKey (P235). Context is a Wikidata taxon Q-id.",
        },
        "HMDB": {
            "display_name": "HMDB",
            "required_identifiers": ["HMDB ID", "Name", "InChIKey"],
            "preferred_identifier_order": ["HMDB ID", "Name", "InChIKey"],
            "context_type": "imported_file",
            "context_label": "Imported HMDB table",
            "online_supported": False,
            "local_data_supported": True,
            "notes": "HMDB lookup checks the user-imported table for the selected rank.",
        },
        "Rhea + UniProt": {
            "display_name": "Rhea + UniProt",
            "required_identifiers": ["InChIKey"],
            "preferred_identifier_order": ["InChIKey"],
            "context_type": "taxon",
            "context_label": "UniProt taxon ID or name",
            "online_supported": True,
            "local_data_supported": False,
            "notes": (
                "Rhea reactions are found via InChIKey→ChEBI (UniChem) then Rhea SPARQL. "
                "Enzyme presence is checked against the full set of Rhea reactions annotated "
                "for the chosen taxon in UniProt (all entries, Swiss-Prot + TrEMBL). "
                "Rows without an InChIKey are marked Not_evaluated."
            ),
        },
        "NPAtlas": {
            "display_name": "NPAtlas",
            "required_identifiers": ["InChIKey"],
            "preferred_identifier_order": ["InChIKey"],
            "context_type": "taxon",
            "context_label": "Producing organism (genus or full species name)",
            "online_supported": True,
            "local_data_supported": False,
            "notes": (
                "NPAtlas records microbial natural products and their producing organisms "
                "(bacteria + fungi). Compound lookup checks whether the InChIKey is recorded "
                "as produced by the selected taxon. Searching a genus name matches all species "
                "under it (descendant-inclusive); typing a full species name matches that "
                "species exactly. No API key required."
            ),
        },
    }

def get_database_retrieval_capabilities() -> Dict[str, Dict[str, Any]]:
    """Return truthful retrieval/lookup implementation status for each database."""
    return {
        "LOTUS": {
            "database": "LOTUS",
            "context_list_status": "not_implemented",
            "context_search_status": "partial",
            "context_retrieval_status": "partial",
            "compound_lookup_status": "implemented",
            "context_source": "Wikidata taxon search used as LOTUS context helper",
            "compound_lookup_source": "Wikidata SPARQL occurrence check (P703 found-in-taxon, matched by InChIKey P235)",
            "context_search_notes": "Wikidata taxon search used as LOTUS context helper; not direct LOTUS occurrence confirmation.",
            "reason": "LOTUS compound occurrence lookup via Wikidata SPARQL is implemented (InChIKey required). Context listing is partial (Wikidata taxon search only).",
        },
        "HMDB": {
            "database": "HMDB",
            "context_list_status": "user_import",
            "context_search_status": "user_import",
            "context_retrieval_status": "user_import",
            "compound_lookup_status": "implemented",
            "context_source": "User-imported HMDB table",
            "compound_lookup_source": "User-imported HMDB table",
            "context_search_notes": "HMDB context is the uploaded per-rank table.",
            "reason": "HMDB lookup checks HMDB ID, Name, then InChIKey in the imported table.",
        },
        "Rhea + UniProt": {
            "database": "Rhea + UniProt",
            "context_list_status": "not_implemented",
            "context_search_status": "partial",
            "context_retrieval_status": "partial",
            "compound_lookup_status": "implemented",
            "context_source": "UniProt taxonomy search API",
            "compound_lookup_source": (
                "Rhea SPARQL (InChIKey→ChEBI via UniChem, then compound→reactions with role). "
                "Enzyme presence: reaction in taxon's UniProt-annotated Rhea set."
            ),
            "context_search_notes": (
                "UniProt taxonomy search by name or numeric taxon ID. "
                "Context listing is not preloadable; search-only."
            ),
            "reason": (
                "Compound reactions and substrate/product roles via Rhea SPARQL. "
                "Enzyme presence checks whether the exact Rhea reaction is annotated "
                "for the selected taxon in UniProt (all entries)."
            ),
        },
        "NPAtlas": {
            "database": "NPAtlas",
            "context_list_status": "not_implemented",
            "context_search_status": "partial",
            "context_retrieval_status": "partial",
            "compound_lookup_status": "implemented",
            "context_source": "NPAtlas /taxon/search API",
            "compound_lookup_source": (
                "NPAtlas /compounds/basicSearch API (InChIKey + taxon name filter). "
                "Found if any record matches in the selected taxon (genus = all species; "
                "full species name = exact match)."
            ),
            "context_search_notes": (
                "NPAtlas taxon search by name/substring. A synthetic row for the raw typed "
                "text allows species-level exact searches not in the taxon node list."
            ),
            "reason": (
                "NPAtlas compound occurrence lookup via REST API (InChIKey required). "
                "Context listing is not preloadable; search-only."
            ),
        },
    }

def get_database_retrieval_capability(database_name: str) -> Dict[str, Any]:
    """Return retrieval/lookup implementation status for one database."""
    canonical = canonical_database_name(database_name)
    return dict(get_database_retrieval_capabilities()[canonical])


def get_required_identifiers(database_name: str) -> list:
    """Return the accepted identifier labels for a database."""
    metadata = get_database_metadata()
    canonical = canonical_database_name(database_name)
    return list(metadata[canonical]["required_identifiers"])


def get_preferred_identifier_order(database_name: str) -> list:
    """Return the preferred identifier order for a database."""
    metadata = get_database_metadata()
    canonical = canonical_database_name(database_name)
    return list(metadata[canonical]["preferred_identifier_order"])


def _as_dict(row: Any) -> Dict[str, Any]:
    if isinstance(row, pd.Series):
        return cast(Dict[str, Any], row.to_dict())
    if isinstance(row, dict):
        return cast(Dict[str, Any], dict(row))
    return {}



def get_best_available_identifier(compound_row: Any, database_name: str) -> Dict[str, Any]:
    """Return the best available identifier for the selected database."""
    canonical = canonical_database_name(database_name)
    row = _as_dict(compound_row)

    if canonical == "HMDB":
        for identifier_type in ("HMDB ID", "Name", "InChIKey"):
            value = normalize_text(row.get(identifier_type, ""))
            if value:
                return {
                    "identifier_type": identifier_type,
                    "identifier_value": value,
                    "status": "available",
                }
        return {
            "identifier_type": None,
            "identifier_value": "",
            "status": "missing",
        }
    if canonical == "LOTUS":
        smiles = normalize_smiles(row.get("SMILES", ""))
        if smiles:
            return {
                "identifier_type": "SMILES",
                "identifier_value": smiles,
                "status": "available",
            }
        inchikey = normalize_inchikey(row.get("InChIKey", ""))
        if inchikey:
            return {
                "identifier_type": "InChIKey",
                "identifier_value": inchikey,
                "status": "available",
            }
        return {
            "identifier_type": None,
            "identifier_value": "",
            "status": "missing",
        }

    if canonical == "Rhea + UniProt":
        inchikey = normalize_inchikey(row.get("InChIKey", ""))
        if inchikey:
            return {
                "identifier_type": "InChIKey",
                "identifier_value": inchikey,
                "status": "available",
            }
        return {
            "identifier_type": None,
            "identifier_value": "",
            "status": "missing",
        }

    if canonical == "NPAtlas":
        inchikey = normalize_inchikey(row.get("InChIKey", ""))
        if inchikey:
            return {
                "identifier_type": "InChIKey",
                "identifier_value": inchikey,
                "status": "available",
            }
        return {
            "identifier_type": None,
            "identifier_value": "",
            "status": "missing",
        }

    return {
        "identifier_type": None,
        "identifier_value": "",
        "status": "missing",
    }


def normalize_text(value: Any) -> str:
    """Normalize text values for safe comparison and display."""
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    text = str(value)
    text = text.replace("\xa0", " ")
    text = text.strip()
    while len(text) >= 2 and text[0] == text[-1] and text[0] in {'"', "'"}:
        text = text[1:-1].strip()
    text = re.sub(r"\s+", " ", text)
    return text



def normalize_inchikey(value: Any) -> str:
    """Normalize InChIKey values by trimming and uppercasing."""
    return normalize_text(value).upper()


def normalize_smiles(value: Any) -> str:
    """Normalize SMILES values without chemical canonicalization."""
    return normalize_text(value)


def make_memory_cache_key(*parts: Any) -> str:
    """Create a stable cache key for in-memory storage."""
    normalized = [normalize_text(part) for part in parts]
    payload = json.dumps(normalized, sort_keys=True, ensure_ascii=True)
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def safe_request_text(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Perform a safe HTTP GET request with retry handling."""
    if not url:
        return {
            "success": False,
            "status_code": None,
            "text": "",
            "error": "Missing URL.",
            "url": url,
        }

    last_error = ""
    for attempt in range(retries + 1):
        try:
            response = requests.get(
                url,
                params=params,
                headers=headers,
                timeout=timeout,
            )
            if response.status_code in {200, 201, 202}:
                return {
                    "success": True,
                    "status_code": response.status_code,
                    "text": response.text,
                    "error": "",
                    "url": response.url,
                }
            if response.status_code in {429, 500, 502, 503, 504}:
                last_error = f"HTTP {response.status_code}"
                if attempt < retries:
                    if verbose:
                        print(f"Retrying request after {last_error}...")
                    time.sleep(1)
                    continue
            else:
                last_error = f"HTTP {response.status_code}"
            return {
                "success": False,
                "status_code": response.status_code,
                "text": response.text,
                "error": last_error,
                "url": response.url,
            }
        except requests.exceptions.Timeout:
            last_error = "Request timeout"
            if attempt < retries:
                time.sleep(1)
                continue
            return {
                "success": False,
                "status_code": None,
                "text": "",
                "error": last_error,
                "url": url,
            }
        except requests.exceptions.RequestException as exc:
            last_error = str(exc)
            if attempt < retries:
                time.sleep(1)
                continue
            return {
                "success": False,
                "status_code": None,
                "text": "",
                "error": last_error,
                "url": url,
            }

    return {
        "success": False,
        "status_code": None,
        "text": "",
        "error": last_error or "Unknown request error",
        "url": url,
    }


def safe_post_json(
    url: str,
    json_body: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    verbose: bool = False,
    extra_params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Perform a safe HTTP POST request with an optional JSON body and/or URL query params."""
    if not url:
        return {"success": False, "status_code": None, "text": "", "error": "Missing URL.", "url": url}

    last_error = ""
    for attempt in range(retries + 1):
        try:
            response = requests.post(url, params=extra_params, json=json_body, headers=headers, timeout=timeout)
            if response.status_code in {200, 201, 202}:
                return {"success": True, "status_code": response.status_code,
                        "text": response.text, "error": "", "url": str(response.url)}
            if response.status_code in {429, 500, 502, 503, 504}:
                last_error = f"HTTP {response.status_code}"
                if attempt < retries:
                    if verbose:
                        print(f"Retrying POST after {last_error}...")
                    time.sleep(1)
                    continue
            else:
                last_error = f"HTTP {response.status_code}"
            return {"success": False, "status_code": response.status_code,
                    "text": response.text, "error": last_error, "url": str(response.url)}
        except requests.exceptions.Timeout:
            last_error = "Request timeout"
            if attempt < retries:
                time.sleep(1)
                continue
            return {"success": False, "status_code": None, "text": "", "error": last_error, "url": url}
        except requests.exceptions.RequestException as exc:
            last_error = str(exc)
            if attempt < retries:
                time.sleep(1)
                continue
            return {"success": False, "status_code": None, "text": "", "error": last_error, "url": url}

    return {"success": False, "status_code": None, "text": "",
            "error": last_error or "Unknown request error", "url": url}


def get_empty_context_dataframe(reason: str = "") -> pd.DataFrame:
    """Return an empty context dataframe with the expected schema."""
    frame = pd.DataFrame(
        columns=[
            "database",
            "context_id",
            "context_name",
            "context_type",
            "source",
            "reason",
        ]
    )
    if reason:
        frame.attrs["reason"] = reason
    return frame




def _build_retrieved_context_dataframe(
    database_name: str,
    rows: list[Dict[str, Any]],
    source: str,
    empty_reason: str,
) -> pd.DataFrame:
    """Build a context dataframe for rows retrieved from an online source."""
    frame = pd.DataFrame(rows)
    if frame.empty:
        return get_empty_context_dataframe(empty_reason)
    expected_columns = [
        "database",
        "context_id",
        "context_name",
        "context_type",
        "source",
        "reason",
    ]
    for column in expected_columns:
        if column not in frame.columns:
            frame[column] = ""
    frame["database"] = database_name
    frame["source"] = source
    frame["reason"] = frame["reason"].fillna("").map(normalize_text)
    frame = frame.drop_duplicates(subset=["context_id", "context_name"]).reset_index(drop=True)
    return cast(pd.DataFrame, frame.loc[:, expected_columns])


def _safe_json_loads(text: str) -> Any:
    try:
        return json.loads(text)
    except Exception:
        return None




def list_lotus_contexts(
    allow_online: bool = True,
    use_memory_cache: bool = True,
    force_refresh: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """Return LOTUS contexts when a lightweight official retrieval path exists."""
    if not allow_online:
        return get_empty_context_dataframe("LOTUS context retrieval requires online access and is not implemented without a lightweight official endpoint.")
    cache_key = make_memory_cache_key("lotus_contexts")
    if use_memory_cache and not force_refresh and cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, pd.DataFrame):
            return cached.copy()
    frame = get_empty_context_dataframe(
        "LOTUS lightweight context retrieval is not implemented because no stable official context endpoint is configured; full dataset download is intentionally avoided."
    )
    if use_memory_cache:
        _MEMORY_CACHE[cache_key] = frame.copy()
    if verbose:
        print(frame.attrs.get("reason", ""))
    return frame


def _require_context_query(database_name: str) -> pd.DataFrame:
    return get_empty_context_dataframe("Search query is required.")


def _cache_context_search(
    database_name: str,
    query: str,
    frame: pd.DataFrame,
    use_memory_cache: bool,
) -> pd.DataFrame:
    if use_memory_cache:
        _MEMORY_CACHE[make_memory_cache_key("context_search", database_name, query)] = frame.copy()
    return frame


def _get_cached_context_search(database_name: str, query: str, use_memory_cache: bool) -> Optional[pd.DataFrame]:
    if not use_memory_cache:
        return None
    cached = _MEMORY_CACHE.get(make_memory_cache_key("context_search", database_name, query))
    if isinstance(cached, pd.DataFrame):
        return cached.copy()
    return None


def search_lotus_contexts_online(
    query: str,
    limit: int = 25,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Search LOTUS contexts only when a lightweight official query path is available."""
    query_text = normalize_text(query)
    if not query_text:
        return _require_context_query("LOTUS")
    if not allow_online:
        return get_empty_context_dataframe("LOTUS context search requires online access.")
    cached = _get_cached_context_search("LOTUS", query_text, use_memory_cache)
    if cached is not None:
        return cached.head(limit).copy()

    # Direct Q-id lookup: if the user typed something like "Q79817", resolve it
    # immediately via wbgetentities and prepend it to results.
    qid_rows: list[Dict[str, Any]] = []
    if re.fullmatch(r"Q\d+", query_text, re.IGNORECASE):
        qid = query_text.upper()
        ge_response = safe_request_text(
            "https://www.wikidata.org/w/api.php",
            params={"action": "wbgetentities", "ids": qid, "languages": "en",
                    "props": "labels|descriptions|claims", "format": "json"},
            headers={**WIKIDATA_HEADERS, "Accept": "application/json"},
            timeout=DEFAULT_TIMEOUT,
            retries=1,
            verbose=verbose,
        )
        if ge_response.get("success"):
            ge_payload = _safe_json_loads(ge_response.get("text", ""))
            entity = (ge_payload or {}).get("entities", {}).get(qid, {})
            if entity.get("id") == qid:
                qid_label = normalize_text(entity.get("labels", {}).get("en", {}).get("value", ""))
                if not qid_label:
                    # Some taxa (e.g. Homo sapiens / Q15978631) carry no English
                    # label; their name lives in P225 (taxon name).
                    p225 = entity.get("claims", {}).get("P225", [])
                    if p225:
                        qid_label = normalize_text(
                            p225[0].get("mainsnak", {}).get("datavalue", {}).get("value", ""))
                qid_desc  = normalize_text(entity.get("descriptions", {}).get("en", {}).get("value", ""))
                if qid_label:
                    qid_rows.append({
                        "database": "LOTUS",
                        "context_id": qid,
                        "context_name": qid_label,
                        "context_type": "taxon_or_species",
                        "source": "Wikidata direct entity lookup (Q-id)",
                        "reason": qid_desc,
                    })

    response = safe_request_text(
        "https://www.wikidata.org/w/api.php",
        params={
            "action": "wbsearchentities",
            "search": query_text,
            "language": "en",
            "format": "json",
            "limit": min(max(limit, 1), 25),
        },
        headers={**WIKIDATA_HEADERS, "Accept": "application/json,text/plain,*/*"},
        timeout=DEFAULT_TIMEOUT,
        retries=1,
        verbose=verbose,
    )
    rows: list[Dict[str, Any]] = []
    if response.get("success"):
        payload = _safe_json_loads(response.get("text", ""))
        if isinstance(payload, dict):
            for item in payload.get("search", []):
                if not isinstance(item, dict):
                    continue
                display_label = item.get("display", {}).get("label", {}).get("value", "")
                label = normalize_text(display_label or item.get("label", ""))
                description = normalize_text(item.get("description", ""))
                entity_id = normalize_text(item.get("id", ""))
                text_for_filter = f"{label} {description}".lower()
                if not label or not entity_id:
                    continue
                if query_text.lower() not in label.lower() and query_text.lower() not in description.lower():
                    continue
                if not any(token in text_for_filter for token in ["taxon", "species", "genus"]):
                    continue
                rows.append(
                    {
                        "database": "LOTUS",
                        "context_id": entity_id,
                        "context_name": label,
                        "context_type": "taxon_or_species",
                        "source": "Wikidata taxon search used as LOTUS context helper",
                        "reason": "",
                    }
                )
    else:
        # Request failed — do NOT cache (caching empty here would poison the
        # session so later identical queries return empty without retrying).
        frame = get_empty_context_dataframe(
            f"LOTUS lightweight taxon/context search failed through Wikidata public route: {response.get('error', '')}"
        )
        return frame.head(limit).copy()

    # Prepend Q-id direct result; deduplicate by context_id.
    qid_ids = {r["context_id"] for r in qid_rows}
    rows = qid_rows + [r for r in rows if r["context_id"] not in qid_ids]

    frame = _build_retrieved_context_dataframe(
        "LOTUS",
        rows,
        "Wikidata taxon search used as LOTUS context helper",
        "LOTUS lightweight taxon/context search returned no real taxon or organism rows.",
    )
    return _cache_context_search("LOTUS", query_text, frame, use_memory_cache).head(limit).copy()


WIKIDATA_SPARQL_URL = "https://query.wikidata.org/sparql"


def get_lotus_occurrences(
    inchikey: str,
    taxon_qid: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[list]:
    """Occurrences of a compound (InChIKey/P235) in taxon_qid or a descendant (P171*),
    each with the matched species label and source DOI.
    None = could not evaluate; [] = none found; non-empty = list of occurrence dicts."""
    key = normalize_inchikey(inchikey)
    qid = normalize_text(taxon_qid)
    if not key or not re.fullmatch(r"Q\d+", qid) or not allow_online or not ENABLE_ONLINE_APIS:
        return None

    cache_key = make_memory_cache_key("lotus_occurrences", key, qid)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return [dict(item) for item in cached]

    query = (
        "PREFIX wd: <http://www.wikidata.org/entity/> "
        "PREFIX wdt: <http://www.wikidata.org/prop/direct/> "
        "PREFIX p: <http://www.wikidata.org/prop/> "
        "PREFIX ps: <http://www.wikidata.org/prop/statement/> "
        "PREFIX wikibase: <http://wikiba.se/ontology#> "
        "PREFIX bd: <http://www.bigdata.com/rdf#> "
        "SELECT DISTINCT ?taxon ?taxonLabel WHERE { "
        f'?c wdt:P235 "{key}" . '
        "?c p:P703 ?stmt . "
        "?stmt ps:P703 ?taxon . "
        f"?taxon wdt:P171* wd:{qid} . "
        'SERVICE wikibase:label { bd:serviceParam wikibase:language "en". } '
        "} LIMIT 500"
    )
    response = safe_request_text(
        WIKIDATA_SPARQL_URL,
        params={"query": query, "format": "json"},
        headers={**WIKIDATA_HEADERS, "Accept": "application/sparql-results+json"},
        timeout=DEFAULT_TIMEOUT,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
    )
    if not response.get("success"):
        return None
    payload = _safe_json_loads(response.get("text", ""))
    if not isinstance(payload, dict):
        return None

    occurrences: list = []
    seen: set = set()
    for binding in payload.get("results", {}).get("bindings", []):
        taxon_uri = normalize_text(binding.get("taxon", {}).get("value", ""))
        species = normalize_text(binding.get("taxonLabel", {}).get("value", ""))
        species_qid = taxon_uri.rsplit("/", 1)[-1] if taxon_uri else ""
        if species_qid in seen:
            continue
        seen.add(species_qid)
        occurrences.append({"species": species, "species_qid": species_qid})

    _MEMORY_CACHE[cache_key] = [dict(item) for item in occurrences]
    return occurrences


def lotus_compound_in_taxon(
    inchikey: str,
    taxon_qid: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[bool]:
    """True if the compound is recorded in taxon_qid or a descendant; None if not evaluable."""
    occurrences = get_lotus_occurrences(inchikey, taxon_qid, allow_online, verbose)
    if occurrences is None:
        return None
    return bool(occurrences)


def get_available_contexts(
    database_name: str,
    force_refresh: bool = False,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Return available contexts for the selected database."""
    canonical = canonical_database_name(database_name)
    if canonical == "LOTUS":
        return list_lotus_contexts(
            allow_online=allow_online,
            use_memory_cache=use_memory_cache,
            force_refresh=force_refresh,
            verbose=verbose,
        )
    if canonical == "Rhea + UniProt":
        # No preloadable full list; context search only.
        return get_empty_context_dataframe(
            "Rhea + UniProt context list is not available; use the search box to find a taxon."
        )
    if canonical == "NPAtlas":
        return get_empty_context_dataframe(
            "NPAtlas context list is not available; use the search box to find a producing organism."
        )
    raise ValueError(f"Unsupported database: {database_name}")


def _build_probe_result(
    database_name: str,
    online_access: str,
    context_listing_possible: bool,
    compound_lookup_possible: bool,
    reason: str,
    tested_url: str = "",
    status_code: Any = None,
) -> Dict[str, Any]:
    """Create a standardized lightweight database access probe result."""
    return {
        "database": canonical_database_name(database_name),
        "online_access": online_access,
        "context_listing_possible": context_listing_possible,
        "compound_lookup_possible": compound_lookup_possible,
        "reason": normalize_text(reason),
        "tested_url": normalize_text(tested_url),
        "status_code": status_code,
    }



def probe_database_access(
    database_name: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Lightly test whether a database has a usable online access path."""
    canonical = canonical_database_name(database_name)

    if canonical == "LOTUS":
        frame = search_lotus_contexts_online("Theobroma", allow_online=allow_online, verbose=verbose)
        return _build_probe_result(
            canonical,
            "partial" if not frame.empty else "not_available",
            not frame.empty,
            False,
            "LOTUS query context search returned retrieved contexts." if not frame.empty else frame.attrs.get("reason", "LOTUS lightweight query context search is not implemented."),
            "",
            None,
        )

    if canonical == "Rhea + UniProt":
        if not allow_online:
            return _build_probe_result(canonical, "not_available", False, False,
                                       "Online probing disabled.", UNIPROT_TAXONOMY_URL, None)
        frame = search_uniprot_taxa_online("Arabidopsis", allow_online=True, verbose=verbose)
        return _build_probe_result(
            canonical,
            "partial" if not frame.empty else "not_available",
            not frame.empty,
            False,
            "UniProt taxonomy search returned results." if not frame.empty
            else frame.attrs.get("reason", "UniProt taxonomy search is not available."),
            UNIPROT_TAXONOMY_URL,
            None,
        )

    if canonical == "HMDB":
        return _build_probe_result(
            canonical, "not_available", False, True,
            "HMDB uses user-imported files; no online probe is needed.", "", None,
        )

    if canonical == "NPAtlas":
        if not allow_online:
            return _build_probe_result(canonical, "not_available", False, False,
                                       "Online probing disabled.", NPATLAS_API_URL, None)
        probe_url = f"{NPATLAS_API_URL}/taxon/search"
        response = safe_post_json(
            probe_url,
            json_body=None,
            headers={**OMISCERA_HEADERS},
            timeout=DEFAULT_TIMEOUT,
            retries=1,
            verbose=verbose,
            extra_params={"taxon": "Streptomyces"},
        )
        ok = response.get("success", False)
        return _build_probe_result(
            canonical,
            "available" if ok else "not_available",
            ok,
            ok,
            "NPAtlas taxon search returned results." if ok
            else response.get("error", "NPAtlas probe failed."),
            probe_url,
            response.get("status_code"),
        )

    raise ValueError(f"Unsupported database: {database_name}")


def probe_all_databases(allow_online: bool = True, verbose: bool = False) -> pd.DataFrame:
    """Return lightweight online feasibility probe results for all databases."""
    rows = [
        probe_database_access(database_name, allow_online=allow_online, verbose=verbose)
        for database_name in SUPPORTED_DATABASES
    ]
    return pd.DataFrame(rows)


def search_contexts(
    database_name: str,
    query: str,
    limit: int = 25,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Search available contexts for a database using case-insensitive matching."""
    canonical = canonical_database_name(database_name)
    if canonical == "LOTUS":
        return search_lotus_contexts_online(
            query,
            limit=limit,
            allow_online=allow_online,
            use_memory_cache=use_memory_cache,
            verbose=verbose,
        )
    if canonical == "Rhea + UniProt":
        return search_uniprot_taxa_online(
            query,
            limit=limit,
            allow_online=allow_online,
            use_memory_cache=use_memory_cache,
            verbose=verbose,
        )
    if canonical == "NPAtlas":
        return search_npatlas_taxa_online(
            query,
            limit=limit,
            allow_online=allow_online,
            use_memory_cache=use_memory_cache,
            verbose=verbose,
        )
    frame = get_available_contexts(
        canonical,
        force_refresh=False,
        allow_online=allow_online,
        use_memory_cache=use_memory_cache,
        verbose=verbose,
    )
    if frame.empty:
        return frame.copy()

    search_term = normalize_text(query).lower()
    if search_term == "":
        subset = frame.head(limit)
    else:
        sortable = frame.copy()
        sortable["_match_rank"] = cast(
            pd.Series,
            sortable.apply(
                lambda row: (
                    0
                    if normalize_text(row.get("context_id", "")).lower() == search_term
                    else 1
                    if normalize_text(row.get("context_name", "")).lower() == search_term
                    else 2
                    if search_term in normalize_text(row.get("context_id", "")).lower()
                    or search_term in normalize_text(row.get("context_name", "")).lower()
                    or search_term in normalize_text(row.get("taxonomy_or_lineage", "")).lower()
                    else 99
                ),
                axis=1,
            ),
        )
        ranked_matches = cast(pd.DataFrame, sortable.loc[sortable["_match_rank"] < 99])
        subset = cast(
            pd.DataFrame,
            pd.concat(
                [
                    ranked_matches.loc[ranked_matches["_match_rank"] == rank]
                    for rank in [0, 1, 2]
                ],
                axis=0,
            )
            .drop(columns=["_match_rank"])
            .head(limit),
        )
    return subset.copy()


def get_context_by_id(
    database_name: str,
    context_id: str,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> Optional[Dict[str, Any]]:
    """Return a matching context dictionary by context ID."""
    frame = get_available_contexts(
        database_name,
        force_refresh=False,
        allow_online=allow_online,
        use_memory_cache=use_memory_cache,
        verbose=verbose,
    )
    if frame.empty:
        return None
    target = normalize_text(context_id).lower()
    for _, row in frame.iterrows():
        if normalize_text(row.get("context_id", "")).lower() == target:
            return cast(Dict[str, Any], row.to_dict())
    return None


def get_context_by_name(
    database_name: str,
    context_name: str,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> Optional[Dict[str, Any]]:
    """Return a matching context dictionary by context name."""
    frame = get_available_contexts(
        database_name,
        force_refresh=False,
        allow_online=allow_online,
        use_memory_cache=use_memory_cache,
        verbose=verbose,
    )
    if frame.empty:
        return None
    target = normalize_text(context_name).lower()
    for _, row in frame.iterrows():
        if normalize_text(row.get("context_name", "")).lower() == target:
            return cast(Dict[str, Any], row.to_dict())
    return None


def build_standard_lookup_result(
    database_name: str,
    status: str,
    reason: str,
    identifier_type: Optional[str] = None,
    identifier_value: Optional[str] = None,
    context_id: Optional[str] = None,
    context_name: Optional[str] = None,
    context_type: Optional[str] = None,
    matched_identifier: Optional[str] = None,
    source: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Create a standardized lookup result dictionary."""
    if extra is None:
        extra = {}
    return {
        "database": canonical_database_name(database_name),
        "status": status,
        "reason": normalize_text(reason),
        "identifier_type": normalize_text(identifier_type),
        "identifier_value": normalize_text(identifier_value),
        "context_id": normalize_text(context_id),
        "context_name": normalize_text(context_name),
        "context_type": normalize_text(context_type),
        "matched_identifier": normalize_text(matched_identifier),
        "source": normalize_text(source),
        "extra": extra,
    }


def _normalize_context_selection(context_selection: Any) -> Optional[Dict[str, Any]]:
    if context_selection is None:
        return None
    if isinstance(context_selection, dict):
        return {
            "context_id": normalize_text(context_selection.get("context_id", "")),
            "context_name": normalize_text(context_selection.get("context_name", "")),
            "context_type": normalize_text(context_selection.get("context_type", "")),
        }
    if isinstance(context_selection, pd.Series):
        return {
            "context_id": normalize_text(context_selection.get("context_id", "")),
            "context_name": normalize_text(context_selection.get("context_name", "")),
            "context_type": normalize_text(context_selection.get("context_type", "")),
        }
    return None


def _context_is_available(database_name: str, context_selection: Any) -> bool:
    context = _normalize_context_selection(context_selection)
    if context is None:
        return False
    if not context["context_id"] and not context["context_name"]:
        return False
    metadata = get_database_metadata()
    required = metadata[canonical_database_name(database_name)]["context_type"]
    if required and context["context_type"]:
        return normalize_text(context["context_type"]).lower() == normalize_text(required).lower()
    return True


def _normalize_lookup_value(value: Any) -> str:
    return normalize_text(value).lower()


def _first_matching_table_value(imported_table: pd.DataFrame, column_name: str, value: str) -> Optional[str]:
    if column_name not in imported_table.columns:
        return None
    target = _normalize_lookup_value(value)
    if not target:
        return None
    matches = imported_table[column_name].apply(_normalize_lookup_value) == target
    if not bool(matches.any()):
        return None
    matched_value = imported_table.loc[matches, column_name].iloc[0]
    return normalize_text(matched_value)


def canonical_hmdb_id_set(value: Any) -> set:
    """Extract and canonicalize all HMDB IDs from a cell value.

    Handles multi-valued cells (semicolons, spaces) and zero-pads numeric parts
    to exactly 7 digits: HMDB97 and HMDB0000097 both become 'HMDB0000097'.
    """
    text = normalize_text(value)
    if not text:
        return set()
    ids = set()
    for match in re.finditer(r"(?i)HMDB(\d+)", text):
        ids.add("HMDB" + match.group(1).zfill(7))
    return ids


def search_hmdb_compound(
    compound_row: Any,
    context_selection: Any,
    imported_table: Optional[pd.DataFrame],
    verbose: bool = False,
    allow_online: bool = True,
) -> Dict[str, Any]:
    """Search a compound in the user-imported HMDB table for the selected rank.

    Match priority: HMDB ID (robust, canonical) → InChIKey → Name (exact) → SMILES.
    """
    del verbose, allow_online
    context = _normalize_context_selection(context_selection)
    row = _as_dict(compound_row)

    ctx_id = context.get("context_id") if context else None
    ctx_name = context.get("context_name") if context else None
    ctx_type = context.get("context_type") if context else None

    if imported_table is None or imported_table.empty:
        return build_standard_lookup_result(
            "HMDB",
            "not_evaluated",
            "No HMDB table imported",
            context_id=ctx_id,
            context_name=ctx_name,
            context_type=ctx_type,
            source="User-imported HMDB table",
        )

    hmdb_id_val = normalize_text(row.get("HMDB ID", ""))
    inchikey_val = normalize_inchikey(row.get("InChIKey", ""))
    name_val = normalize_text(row.get("Name", ""))
    smiles_val = normalize_text(row.get("SMILES", ""))

    has_any = any([hmdb_id_val, inchikey_val, name_val, smiles_val])
    if not has_any:
        return build_standard_lookup_result(
            "HMDB",
            "not_evaluated",
            "No HMDB identifier",
            context_id=ctx_id,
            context_name=ctx_name,
            context_type=ctx_type,
            source="User-imported HMDB table",
        )

    primary_identifier_type = (
        "HMDB ID" if hmdb_id_val
        else "InChIKey" if inchikey_val
        else "Name" if name_val
        else "SMILES"
    )
    primary_identifier_value = (
        hmdb_id_val if hmdb_id_val
        else inchikey_val if inchikey_val
        else name_val if name_val
        else smiles_val
    )

    def _found(field: str, value: str, matched: str) -> Dict[str, Any]:
        return build_standard_lookup_result(
            "HMDB",
            "found",
            f"Matched by {field}",
            identifier_type=field,
            identifier_value=value,
            context_id=ctx_id,
            context_name=ctx_name,
            context_type=ctx_type,
            matched_identifier=matched,
            source="User-imported HMDB table",
            extra={"matched_field": field, "matched_value": matched},
        )

    # Build list of checkable fields: compound must have the value AND table must have the column.
    checkable: list = []
    if hmdb_id_val and "HMDB ID" in imported_table.columns:
        checkable.append(("HMDB ID", hmdb_id_val))
    if inchikey_val and "InChIKey" in imported_table.columns:
        checkable.append(("InChIKey", inchikey_val))
    if name_val and "Name" in imported_table.columns:
        checkable.append(("Name", name_val))
    if smiles_val and "SMILES" in imported_table.columns:
        checkable.append(("SMILES", smiles_val))

    if not checkable:
        return build_standard_lookup_result(
            "HMDB",
            "not_evaluated",
            "Compound has no identifier matching a column in the imported HMDB table",
            identifier_type=primary_identifier_type,
            identifier_value=primary_identifier_value,
            context_id=ctx_id,
            context_name=ctx_name,
            context_type=ctx_type,
            source="User-imported HMDB table",
        )

    # Match in priority order; column existence is guaranteed by checkable construction.
    for field, value in checkable:
        if field == "HMDB ID":
            compound_ids = canonical_hmdb_id_set(value)
            if compound_ids:
                table_ids: set = set()
                for cell in imported_table["HMDB ID"]:
                    table_ids |= canonical_hmdb_id_set(cell)
                overlap = compound_ids & table_ids
                if overlap:
                    return _found("HMDB ID", value, next(iter(overlap)))
        elif field == "InChIKey":
            matches = imported_table["InChIKey"].apply(lambda v: normalize_inchikey(v)) == value
            if bool(matches.any()):
                matched_val = normalize_text(str(imported_table["InChIKey"][matches].iat[0]))
                return _found("InChIKey", value, matched_val)
        elif field == "Name":
            matched = _first_matching_table_value(imported_table, "Name", value)
            if matched is not None:
                return _found("Name", value, matched)
        elif field == "SMILES":
            smiles_stripped = value.strip()
            matches = imported_table["SMILES"].apply(lambda v: normalize_text(v).strip()) == smiles_stripped
            if bool(matches.any()):
                matched_val = normalize_text(str(imported_table["SMILES"][matches].iat[0]))
                return _found("SMILES", value, matched_val)

    return build_standard_lookup_result(
        "HMDB",
        "not_found",
        "No matching HMDB row",
        identifier_type=primary_identifier_type,
        identifier_value=primary_identifier_value,
        context_id=ctx_id,
        context_name=ctx_name,
        context_type=ctx_type,
        source="User-imported HMDB table",
        extra={"checked_fields": [f for f, _ in checkable]},
    )


# ── Rhea + UniProt ─────────────────────────────────────────────────────────────

def _get_chebi_ids_for_inchikey(
    inchikey: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> list:
    """Return ChEBI IDs (e.g. ['CHEBI:4167']) for an InChIKey via UniChem."""
    key = normalize_inchikey(inchikey)
    if not key or not allow_online or not ENABLE_ONLINE_APIS:
        return []

    cache_key = make_memory_cache_key("unichem_chebi", key)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return list(cached)

    response = safe_request_text(
        f"{UNICHEM_INCHIKEY_URL}/{key}",
        headers=OMISCERA_HEADERS,
        timeout=30,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
    )
    if not response.get("success"):
        return []

    try:
        data = json.loads(response.get("text", "[]"))
    except Exception:
        return []

    chebi_ids = [
        normalize_text(item.get("src_compound_id", ""))
        for item in (data if isinstance(data, list) else [])
        if str(item.get("src_id", "")) == "7" and item.get("src_compound_id")
    ]
    _MEMORY_CACHE[cache_key] = chebi_ids
    return chebi_ids


def search_uniprot_taxa_online(
    query: str,
    limit: int = 25,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Search UniProt taxonomy for taxa matching a name or numeric taxon ID."""
    query_text = normalize_text(query)
    if not query_text:
        return _require_context_query("Rhea + UniProt")
    if not allow_online:
        return get_empty_context_dataframe("Rhea + UniProt context search requires online access.")

    cached = _get_cached_context_search("Rhea + UniProt", query_text, use_memory_cache)
    if cached is not None:
        return cached.head(limit).copy()

    rows: list = []

    # Numeric query → direct taxon lookup for exact match
    if re.fullmatch(r"\d+", query_text):
        r = safe_request_text(
            f"{UNIPROT_TAXONOMY_URL}/{query_text}",
            params={"format": "json"},
            headers=OMISCERA_HEADERS,
            timeout=DEFAULT_TIMEOUT,
            retries=1,
            verbose=verbose,
        )
        if r.get("success"):
            try:
                data = json.loads(r.get("text", "{}"))
            except Exception:
                data = {}
            if data.get("taxonId"):
                rows.append({
                    "database": "Rhea + UniProt",
                    "context_id": str(data["taxonId"]),
                    "context_name": normalize_text(data.get("scientificName", "")),
                    "context_type": "taxon",
                    "source": "UniProt taxonomy direct numeric lookup",
                    "reason": normalize_text(data.get("rank", "")),
                })

    # Name-based search
    response = safe_request_text(
        f"{UNIPROT_TAXONOMY_URL}/search",
        params={"query": query_text, "fields": "id,scientific_name,rank", "format": "json",
                "size": min(max(limit, 1), 25)},
        headers=OMISCERA_HEADERS,
        timeout=DEFAULT_TIMEOUT,
        retries=1,
        verbose=verbose,
    )
    if response.get("success"):
        try:
            payload = json.loads(response.get("text", "{}"))
        except Exception:
            payload = {}
        seen_ids = {r["context_id"] for r in rows}
        for item in payload.get("results", []):
            taxon_id = str(item.get("taxonId", "")).strip()
            name = normalize_text(item.get("scientificName", ""))
            if not taxon_id or not name or taxon_id in seen_ids:
                continue
            rows.append({
                "database": "Rhea + UniProt",
                "context_id": taxon_id,
                "context_name": name,
                "context_type": "taxon",
                "source": "UniProt taxonomy search",
                "reason": normalize_text(item.get("rank", "")),
            })
            seen_ids.add(taxon_id)
    elif not rows:
        # Request failed and no prior rows — do NOT cache (avoid poisoning the
        # session cache with an empty result that blocks later retries).
        frame = get_empty_context_dataframe(
            f"UniProt taxonomy search failed: {response.get('error', '')}"
        )
        return frame.head(limit).copy()

    frame = _build_retrieved_context_dataframe(
        "Rhea + UniProt",
        rows,
        "UniProt taxonomy search",
        "UniProt taxonomy search returned no results.",
    )
    return _cache_context_search("Rhea + UniProt", query_text, frame, use_memory_cache).head(limit).copy()


def _rhea_reactions_for_chebi_ids(
    chebi_ids: list,
    allow_online: bool = True,
    verbose: bool = False,
) -> list:
    """Return Rhea reactions for a list of ChEBI IDs.

    Shared core used by both InChIKey and SMILES lookup paths.
    Each entry: {"rhea_id": "RHEA:NNNNN", "ec": "N.N.N.N" or "", "role": "substrate"|"product"|"both"}.
    """
    reactions: Dict[str, Dict[str, Any]] = {}

    for chebi_raw in chebi_ids:
        chebi_num = re.sub(r"(?i)chebi[:\s]*", "", normalize_text(chebi_raw)).strip()
        if not chebi_num:
            continue

        sparql_query = (
            "PREFIX rh: <http://rdf.rhea-db.org/> "
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
            "SELECT DISTINCT ?rheaAcc ?ec ?role WHERE { "
            "  ?reaction rdfs:subClassOf rh:Reaction . "
            "  ?reaction rh:accession ?rheaAcc . "
            "  ?reaction rh:side ?side . "
            "  ?side rh:contains ?part . "
            "  ?part rh:compound ?cpd . "
            f" ?cpd rh:chebi <http://purl.obolibrary.org/obo/CHEBI_{chebi_num}> . "
            "  OPTIONAL { "
            "    ?reaction rh:ec ?ecURI . "
            "    BIND(STRAFTER(STR(?ecURI), \"http://purl.uniprot.org/enzyme/\") AS ?ec) "
            "  } "
            "  BIND(IF(STRENDS(STR(?side), \"_L\"), \"substrate\", "
            "       IF(STRENDS(STR(?side), \"_R\"), \"product\", \"unknown\")) AS ?role) "
            "}"
        )

        response = safe_request_text(
            RHEA_SPARQL_URL,
            params={"query": sparql_query, "format": "json"},
            headers={**OMISCERA_HEADERS, "Accept": "application/sparql-results+json"},
            timeout=DEFAULT_TIMEOUT,
            retries=DEFAULT_RETRIES,
            verbose=verbose,
        )
        if not response.get("success"):
            continue

        try:
            payload = json.loads(response.get("text", "{}"))
        except Exception:
            continue

        for binding in payload.get("results", {}).get("bindings", []):
            rhea_id = normalize_text(binding.get("rheaAcc", {}).get("value", ""))
            ec = normalize_text(binding.get("ec", {}).get("value", ""))
            role = normalize_text(binding.get("role", {}).get("value", "unknown"))
            if not rhea_id:
                continue
            if rhea_id not in reactions:
                reactions[rhea_id] = {"ec": ec, "roles": set()}
            reactions[rhea_id]["roles"].add(role)

    result = []
    for rhea_id, data in reactions.items():
        roles = data["roles"]
        if "substrate" in roles and "product" in roles:
            role_final = "both"
        elif "substrate" in roles:
            role_final = "substrate"
        elif "product" in roles:
            role_final = "product"
        else:
            role_final = "unknown"
        result.append({"rhea_id": rhea_id, "ec": data["ec"], "role": role_final})

    return result


def get_rhea_reactions_for_compound(
    inchikey: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[list]:
    """Return Rhea reactions for a compound identified by InChIKey.

    Each entry: {"rhea_id": "RHEA:NNNNN", "ec": "N.N.N.N" or "", "role": "substrate"|"product"|"both"}.
    None = could not evaluate (offline or lookup error). [] = no reactions found.
    Role uses the Rhea undirected-reaction convention: _L side = substrate, _R side = product.
    """
    key = normalize_inchikey(inchikey)
    if not key or not allow_online or not ENABLE_ONLINE_APIS:
        return None

    cache_key = make_memory_cache_key("rhea_reactions_for_compound", key)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return [dict(item) for item in cached]

    chebi_ids = _get_chebi_ids_for_inchikey(key, allow_online, verbose)
    if not chebi_ids:
        _MEMORY_CACHE[cache_key] = []
        return []

    result = _rhea_reactions_for_chebi_ids(chebi_ids, allow_online, verbose)
    _MEMORY_CACHE[cache_key] = [dict(item) for item in result]
    return result


def _get_chebi_ids_for_smiles(
    smiles: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> list:
    """Return ChEBI IDs matching an exact SMILES string via Rhea's ChEBI snapshot.

    Uses the chebi:smiles predicate in Rhea SPARQL. Matches only when the given
    SMILES is byte-identical to ChEBI's stored SMILES (no local canonicalization).
    Returns [] when no match or offline.
    """
    smi = normalize_smiles(smiles)
    if not smi or not allow_online or not ENABLE_ONLINE_APIS:
        return []

    cache_key = make_memory_cache_key("rhea_chebi_for_smiles", smi)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return list(cached)

    escaped = smi.replace("\\", "\\\\").replace('"', '\\"')
    sparql_query = (
        "PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/> "
        "SELECT DISTINCT ?chebiURI WHERE { "
        f'  ?chebiURI chebi:smiles "{escaped}" . '
        "}"
    )
    response = safe_request_text(
        RHEA_SPARQL_URL,
        params={"query": sparql_query, "format": "json"},
        headers={**OMISCERA_HEADERS, "Accept": "application/sparql-results+json"},
        timeout=DEFAULT_TIMEOUT,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
    )
    if not response.get("success"):
        return []

    try:
        payload = json.loads(response.get("text", "{}"))
    except Exception:
        return []

    chebi_ids = []
    for binding in payload.get("results", {}).get("bindings", []):
        uri = binding.get("chebiURI", {}).get("value", "")
        chebi_num = uri.rsplit("CHEBI_", 1)[-1] if "CHEBI_" in uri else ""
        if chebi_num:
            chebi_ids.append(chebi_num)

    _MEMORY_CACHE[cache_key] = chebi_ids
    return chebi_ids


def get_rhea_reactions_for_smiles(
    smiles: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[list]:
    """Return Rhea reactions for a compound identified by its SMILES string.

    Matches only when the SMILES is byte-identical to ChEBI's stored form.
    None = could not evaluate. [] = no match / no reactions.
    """
    smi = normalize_smiles(smiles)
    if not smi or not allow_online or not ENABLE_ONLINE_APIS:
        return None

    cache_key = make_memory_cache_key("rhea_reactions_for_smiles", smi)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return [dict(item) for item in cached]

    chebi_ids = _get_chebi_ids_for_smiles(smi, allow_online, verbose)
    if not chebi_ids:
        _MEMORY_CACHE[cache_key] = []
        return []

    result = _rhea_reactions_for_chebi_ids(chebi_ids, allow_online, verbose)
    _MEMORY_CACHE[cache_key] = [dict(item) for item in result]
    return result


def get_uniprot_taxon_rhea_reactions(
    taxon_id: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[set]:
    """Return the set of Rhea reaction accessions annotated for a taxon in UniProt.

    Queries all entries (Swiss-Prot + TrEMBL). Paginated in batches of 500.
    None = could not evaluate. set() = empty (no annotated reactions found).
    """
    tid = normalize_text(taxon_id).strip()
    if not tid or not re.fullmatch(r"\d+", tid):
        return None
    if not allow_online or not ENABLE_ONLINE_APIS:
        return None

    cache_key = make_memory_cache_key("uniprot_taxon_rhea_reactions", tid)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, set):
            return set(cached)

    sparql_template = (
        "PREFIX up: <http://purl.uniprot.org/core/> "
        "PREFIX taxon: <http://purl.uniprot.org/taxonomy/> "
        "PREFIX rh: <http://rdf.rhea-db.org/> "
        "SELECT DISTINCT ?rheaAcc WHERE {{ "
        "  ?protein a up:Protein . "
        f" ?protein up:organism taxon:{tid} . "
        "  ?protein up:annotation ?annot . "
        "  ?annot a up:Catalytic_Activity_Annotation . "
        "  ?annot up:catalyticActivity ?ca . "
        "  ?ca up:catalyzedReaction ?reaction . "
        "  ?reaction rh:accession ?rheaAcc . "
        "}} "
        "LIMIT 500 "
        "OFFSET {offset}"
    )

    all_reactions: set = set()
    offset = 0
    _MAX_PAGES = 40  # cap at 20 000 reactions to avoid runaway queries on high-rank taxa
    while offset < _MAX_PAGES * 500:
        response = safe_request_text(
            UNIPROT_SPARQL_URL,
            params={"query": sparql_template.format(offset=offset), "format": "json"},
            headers={**OMISCERA_HEADERS, "Accept": "application/sparql-results+json"},
            timeout=60,
            retries=DEFAULT_RETRIES,
            verbose=verbose,
        )
        if not response.get("success"):
            if not all_reactions:
                return None
            break
        try:
            payload = json.loads(response.get("text", "{}"))
        except Exception:
            break
        batch = [
            normalize_text(b.get("rheaAcc", {}).get("value", ""))
            for b in payload.get("results", {}).get("bindings", [])
            if b.get("rheaAcc", {}).get("value")
        ]
        all_reactions.update(batch)
        if len(batch) < 500:
            break
        offset += 500

    _MEMORY_CACHE[cache_key] = set(all_reactions)
    return all_reactions


def _fetch_ec_names_batch(
    ec_list: list,
    allow_online: bool = True,
    verbose: bool = False,
) -> Dict[str, str]:
    """Return a dict mapping EC number → preferred name via UniProt SPARQL skos:prefLabel.

    Uses _MEMORY_CACHE per EC; only queries uncached entries. Returns {} on failure.
    Never uses external enzyme databases.
    """
    if not ec_list or not allow_online or not ENABLE_ONLINE_APIS:
        return {}

    result: Dict[str, str] = {}
    to_fetch = []
    for ec in ec_list:
        cache_key = make_memory_cache_key("ec_name", ec)
        if cache_key in _MEMORY_CACHE:
            cached = _MEMORY_CACHE[cache_key]
            if isinstance(cached, str):
                result[ec] = cached
        else:
            to_fetch.append(ec)

    if not to_fetch:
        return result

    values_clause = " ".join(
        f"<http://purl.uniprot.org/enzyme/{ec}>" for ec in to_fetch
    )
    sparql_query = (
        "PREFIX skos: <http://www.w3.org/2004/02/skos/core#> "
        "SELECT ?ecURI ?name WHERE { "
        f"  VALUES ?ecURI {{ {values_clause} }} "
        "  ?ecURI skos:prefLabel ?name . "
        "}"
    )
    response = safe_request_text(
        UNIPROT_SPARQL_URL,
        params={"query": sparql_query, "format": "json"},
        headers={**OMISCERA_HEADERS, "Accept": "application/sparql-results+json"},
        timeout=DEFAULT_TIMEOUT,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
    )
    if response.get("success"):
        try:
            payload = json.loads(response.get("text", "{}"))
        except Exception:
            payload = {}
        for binding in payload.get("results", {}).get("bindings", []):
            ec_uri = binding.get("ecURI", {}).get("value", "")
            name = normalize_text(binding.get("name", {}).get("value", ""))
            if ec_uri and name:
                ec_num = ec_uri.rsplit("/", 1)[-1]
                result[ec_num] = name
                _MEMORY_CACHE[make_memory_cache_key("ec_name", ec_num)] = name

    # Cache misses so we don't re-query
    for ec in to_fetch:
        if ec not in result:
            _MEMORY_CACHE[make_memory_cache_key("ec_name", ec)] = ""

    return result


def build_rhea_uniprot_enzyme_rows(
    compound_row: Any,
    taxon_id: str,
    allow_online: bool = True,
    verbose: bool = False,
    taxon_reactions: Optional[set] = None,
) -> list:
    """Return enzyme rows for a compound against a UniProt taxon via Rhea.

    Each row: {"Name", "Enzyme", "Role", "Enzyme presence"}.
    Enzyme is "N.N.N.N - name" (UniProt skos:prefLabel) or bare EC when
    the name is unavailable, or "RHEA:NNNNN" for EC-less reactions.
    Role: Substrate / Product / Both / Unknown.
    Enzyme presence: Found / Not_found / Not_evaluated.

    Pass a precomputed taxon_reactions set to skip the (heavy, paginated) UniProt
    SPARQL fetch — useful when processing multiple compounds for the same taxon.
    taxon_reactions=None means "fetch it now"; pass an empty set explicitly to mean
    "we already know this taxon has no reactions."
    """
    row = _as_dict(compound_row)
    compound_name = normalize_text(row.get("Name", ""))
    inchikey = normalize_inchikey(row.get("InChIKey", ""))
    smiles = normalize_smiles(row.get("SMILES", ""))

    def _not_evaluated(note: str) -> list:
        return [{"Name": compound_name, "Enzyme": note, "Role": "", "Enzyme presence": "Not_evaluated"}]

    if not inchikey and not smiles:
        return _not_evaluated("No InChIKey or SMILES — Rhea lookup skipped")

    reactions: Optional[list] = None
    if inchikey:
        reactions = get_rhea_reactions_for_compound(inchikey, allow_online, verbose)

    # SMILES fallback: try when InChIKey missing or yielded nothing
    if (not reactions) and smiles:
        smiles_reactions = get_rhea_reactions_for_smiles(smiles, allow_online, verbose)
        if smiles_reactions or reactions is None:
            reactions = smiles_reactions

    if reactions is None:
        return _not_evaluated("Rhea lookup failed")
    if not reactions:
        return _not_evaluated("No Rhea reactions found for this compound")

    if taxon_reactions is None:
        taxon_reactions = get_uniprot_taxon_rhea_reactions(taxon_id, allow_online, verbose)

    # Fetch EC names in one batched UniProt SPARQL call
    unique_ecs = list({rxn["ec"] for rxn in reactions if rxn["ec"]})
    ec_names = _fetch_ec_names_batch(unique_ecs, allow_online, verbose)

    # Group reactions by EC (or by rhea_id when no EC), aggregating roles
    ec_groups: Dict[str, Dict[str, Any]] = {}
    for rxn in reactions:
        rhea_id = rxn["rhea_id"]
        ec = rxn["ec"]
        role = rxn["role"]
        key = ec if ec else rhea_id
        if key not in ec_groups:
            ec_groups[key] = {"ec": ec, "rhea_ids": [], "roles": set()}
        ec_groups[key]["rhea_ids"].append(rhea_id)
        ec_groups[key]["roles"].add(role)

    result_rows = []
    for key, data in ec_groups.items():
        roles = data["roles"]
        if "substrate" in roles and "product" in roles:
            role_label = "Both"
        elif "substrate" in roles:
            role_label = "Substrate"
        elif "product" in roles:
            role_label = "Product"
        else:
            role_label = "Unknown"

        ec = data["ec"]
        if ec:
            name = ec_names.get(ec, "")
            enzyme_label = f"{ec} - {name}" if name else ec
        else:
            enzyme_label = key

        if taxon_reactions is None:
            presence = "Not_evaluated"
        else:
            found = any(rhea_id in taxon_reactions for rhea_id in data["rhea_ids"])
            presence = "Found" if found else "Not_found"

        result_rows.append({
            "Name": compound_name,
            "Enzyme": enzyme_label,
            "Role": role_label,
            "Enzyme presence": presence,
        })

    return sorted(result_rows, key=lambda r: r["Enzyme"])


def _npatlas_taxon_has_producers(
    taxon_name: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> bool:
    """Return True if NPAtlas has at least one compound produced by taxon_name.

    Uses basicSearch with taxon as URL query param (server-side organism match).
    Works for genus names (matches any species under it) and full species names.
    Result is cached per name. Returns False when offline or on API failure.
    """
    name = normalize_text(taxon_name).strip()
    if not name or not allow_online or not ENABLE_ONLINE_APIS:
        return False
    # Fast pre-guard: reject obvious non-organism text before hitting the API
    if not re.fullmatch(r"[A-Za-z][A-Za-z\s\-\.]+", name):
        return False

    cache_key = make_memory_cache_key("npatlas_has_producers", name)
    if cache_key in _MEMORY_CACHE:
        return bool(_MEMORY_CACHE[cache_key])

    response = safe_post_json(
        f"{NPATLAS_API_URL}/compounds/basicSearch",
        json_body=None,
        headers={**OMISCERA_HEADERS},
        timeout=DEFAULT_TIMEOUT,
        retries=1,
        verbose=verbose,
        extra_params={"taxon": name, "limit": 1},
    )
    if not response.get("success"):
        return False
    try:
        data = json.loads(response.get("text", "[]"))
    except Exception:
        return False
    result = isinstance(data, list) and len(data) > 0
    _MEMORY_CACHE[cache_key] = result
    return result


def search_npatlas_taxa_online(
    query: str,
    limit: int = 25,
    allow_online: bool = True,
    use_memory_cache: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """Search NPAtlas for taxa matching a name or substring.

    Returns a context DataFrame with validated results only. An empty frame means
    the query matched no real NPAtlas organism (same behavior as LOTUS/Rhea+UniProt).
    """
    query_text = normalize_text(query)
    if not query_text:
        return _require_context_query("NPAtlas")
    if not allow_online:
        return get_empty_context_dataframe("NPAtlas taxon search requires online access.")

    cached = _get_cached_context_search("NPAtlas", query_text, use_memory_cache)
    if cached is not None:
        return cached.head(limit).copy()

    rows: list = []
    # NPAtlas taxon/search requires "taxon" as a URL query parameter (not JSON body)
    response = safe_post_json(
        f"{NPATLAS_API_URL}/taxon/search",
        json_body=None,
        headers={**OMISCERA_HEADERS},
        timeout=DEFAULT_TIMEOUT,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
        extra_params={"taxon": query_text},
    )
    if response.get("success"):
        try:
            data = json.loads(response.get("text", "[]"))
        except Exception:
            data = []
        seen_ids: set = set()
        for item in (data if isinstance(data, list) else []):
            taxon_id = str(item.get("id", "")).strip()
            name = normalize_text(item.get("name", ""))
            rank = normalize_text(item.get("rank", ""))
            if not name:
                continue
            key = taxon_id or name
            if key in seen_ids:
                continue
            seen_ids.add(key)
            rows.append({
                "database": "NPAtlas",
                "context_id": taxon_id,
                "context_name": name,
                "context_type": "taxon",
                "source": "NPAtlas taxon search",
                "reason": rank,
            })

    # Species-level typed row: add only when the query is NOT already a taxon node
    # (NPAtlas taxonomy only goes to genus; species are free-text in original_organism)
    # and the name is validated against real NPAtlas data.
    # The server's basicSearch taxon filter only matches genus-level names, not
    # species strings — so for a multi-word query (e.g. "Curvularia geniculata") we
    # validate the genus part (first word) instead.
    node_names = {normalize_text(r.get("context_name", "")).lower() for r in rows}
    if query_text.lower() not in node_names:
        words = query_text.split()
        validate_name = words[0] if len(words) >= 2 else query_text
        if _npatlas_taxon_has_producers(validate_name, allow_online=allow_online, verbose=verbose):
            rows.append({
                "database": "NPAtlas",
                "context_id": "",
                "context_name": query_text,
                "context_type": "taxon_name_string",
                "source": "NPAtlas organism name (validated)",
                "reason": "Validated producer name",
            })

    frame = _build_retrieved_context_dataframe(
        "NPAtlas",
        rows,
        "NPAtlas taxon search",
        "NPAtlas taxon search returned no results.",
    )
    # Only cache when the taxon/search request actually succeeded. Caching the
    # empty result of a failed/timed-out request would poison the session cache
    # so later identical queries return empty without retrying (the "not finding
    # except the one I chose before" symptom).
    if response.get("success"):
        return _cache_context_search("NPAtlas", query_text, frame, use_memory_cache).head(limit).copy()
    return frame.head(limit).copy()


def get_npatlas_occurrences(
    inchikey: str,
    taxon_name: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[list]:
    """Return NPAtlas records where a compound (InChIKey) is produced by the given taxon.

    The NPAtlas basicSearch API requires inchikey as a URL query parameter.
    Taxon filtering is done client-side: genus name matches any species containing
    it (e.g. 'Streptomyces' matches 'Streptomyces coelicolor'); a full species name
    matches exactly.
    None = could not evaluate; [] = not recorded; non-empty = list of {organism, doi}.
    """
    key = normalize_inchikey(inchikey)
    taxon = normalize_text(taxon_name).strip()
    if not key or not taxon or not allow_online or not ENABLE_ONLINE_APIS:
        return None

    cache_key = make_memory_cache_key("npatlas_occurrences", key, taxon)
    if cache_key in _MEMORY_CACHE:
        cached = _MEMORY_CACHE[cache_key]
        if isinstance(cached, list):
            return [dict(item) for item in cached]

    # Fetch all producers for this InChIKey via URL query param (not JSON body)
    response = safe_post_json(
        f"{NPATLAS_API_URL}/compounds/basicSearch",
        json_body=None,
        headers={**OMISCERA_HEADERS},
        timeout=DEFAULT_TIMEOUT,
        retries=DEFAULT_RETRIES,
        verbose=verbose,
        extra_params={"inchikey": key, "limit": 100},
    )
    if not response.get("success"):
        return None

    try:
        data = json.loads(response.get("text", "[]"))
    except Exception:
        return None

    # Client-side taxon filter: match anywhere in the organism string
    taxon_lower = taxon.lower()
    records = []
    for item in (data if isinstance(data, list) else []):
        organism = normalize_text(item.get("original_organism", ""))
        if taxon_lower not in organism.lower():
            continue
        doi = normalize_text(item.get("original_doi", ""))
        npaid = normalize_text(str(item.get("npaid", "")))
        records.append({"organism": organism, "doi": doi, "npaid": npaid})

    _MEMORY_CACHE[cache_key] = [dict(r) for r in records]
    return records


def npatlas_compound_in_taxon(
    inchikey: str,
    taxon_name: str,
    allow_online: bool = True,
    verbose: bool = False,
) -> Optional[bool]:
    """True if the compound is recorded as produced by the taxon in NPAtlas; None if not evaluable."""
    records = get_npatlas_occurrences(inchikey, taxon_name, allow_online, verbose)
    if records is None:
        return None
    return bool(records)


def search_npatlas_compound(
    compound_row: Any,
    context_selection: Any,
    verbose: bool = False,
    allow_online: bool = True,
) -> Dict[str, Any]:
    """NPAtlas compound-in-taxon lookup via the basicSearch API (InChIKey + taxon name)."""
    row = _as_dict(compound_row)
    context = _normalize_context_selection(context_selection)
    ctx_id = context.get("context_id") if context else None
    ctx_name = context.get("context_name") if context else None
    ctx_type = context.get("context_type") if context else None
    source = "NPAtlas /compounds/basicSearch (InChIKey + taxon name filter)"
    inchikey = normalize_inchikey(row.get("InChIKey", ""))

    def _result(status: str, reason: str) -> Dict[str, Any]:
        return build_standard_lookup_result(
            "NPAtlas", status, reason,
            identifier_type="InChIKey", identifier_value=inchikey,
            context_id=ctx_id, context_name=ctx_name, context_type=ctx_type,
            matched_identifier=inchikey if status == "found" else None,
            source=source,
        )

    if not inchikey:
        return _result("not_evaluated", "NPAtlas lookup needs an InChIKey; none on this compound.")
    if not _context_is_available("NPAtlas", context_selection):
        return _result("not_evaluated", "No valid NPAtlas taxon context selected.")

    # Use context_name as the taxon string (covers both node names and typed species strings)
    taxon_name = normalize_text(ctx_name or "").strip()
    if not taxon_name:
        return _result("not_evaluated", "NPAtlas context has no taxon name.")

    records = get_npatlas_occurrences(inchikey, taxon_name, allow_online=allow_online, verbose=verbose)
    if records is None:
        return _result("not_evaluated", "NPAtlas lookup could not be completed (offline or request failed).")
    if records:
        return _result("found", f"Recorded in NPAtlas for '{taxon_name}' ({len(records)} record(s)).")
    return _result("not_found", f"Not recorded in NPAtlas for '{taxon_name}'.")


def search_compound_in_database(
    database_name: str,
    compound_row: Any,
    context_selection: Any,
    verbose: bool = False,
    allow_online: bool = True,
    imported_table: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """Dispatch compound lookup to the implemented database-specific search."""
    canonical = canonical_database_name(database_name)
    if canonical == "HMDB":
        return search_hmdb_compound(compound_row, context_selection, imported_table, verbose, allow_online)
    if canonical == "LOTUS":
        return search_lotus_compound(compound_row, context_selection, verbose, allow_online)
    if canonical == "Rhea + UniProt":
        return build_standard_lookup_result(
            "Rhea + UniProt",
            "not_evaluated",
            "Rhea + UniProt is evaluated separately via the enzyme presence table",
        )
    if canonical == "NPAtlas":
        return search_npatlas_compound(compound_row, context_selection, verbose, allow_online)
    raise ValueError(f"Unsupported database: {database_name}")


def search_lotus_compound(
    compound_row: Any,
    context_selection: Any,
    verbose: bool = False,
    allow_online: bool = True,
) -> Dict[str, Any]:
    """LOTUS compound-in-taxon lookup via Wikidata (P235 InChIKey + P703 found-in-taxon)."""
    row = _as_dict(compound_row)
    context = _normalize_context_selection(context_selection)
    ctx_id = context.get("context_id") if context else None
    ctx_name = context.get("context_name") if context else None
    ctx_type = context.get("context_type") if context else None
    source = "LOTUS via Wikidata occurrence (P703) matched by InChIKey"
    inchikey = normalize_inchikey(row.get("InChIKey", ""))

    def _result(status: str, reason: str) -> Dict[str, Any]:
        return build_standard_lookup_result(
            "LOTUS", status, reason,
            identifier_type="InChIKey", identifier_value=inchikey,
            context_id=ctx_id, context_name=ctx_name, context_type=ctx_type,
            matched_identifier=inchikey if status == "found" else None,
            source=source,
        )

    if not inchikey:
        return _result("not_evaluated", "LOTUS online lookup needs an InChIKey; none on this compound.")
    if not _context_is_available("LOTUS", context_selection):
        return _result("not_evaluated", "No valid LOTUS taxon context selected.")
    if not re.fullmatch(r"Q\d+", normalize_text(ctx_id or "")):
        return _result("not_evaluated", "LOTUS context is not a Wikidata taxon id.")

    occurrences = get_lotus_occurrences(inchikey, ctx_id or "", allow_online=allow_online, verbose=verbose)
    if occurrences is None:
        return _result("not_evaluated", "LOTUS/Wikidata lookup could not be completed (offline or request failed).")
    if occurrences:
        return _result("found", "Recorded in LOTUS for the selected taxon (or a descendant).")
    return _result("not_found", "Not recorded in LOTUS for the selected taxon.")


def validate_database_catalog_contract(allow_online: bool = True) -> dict:
    """Run lightweight contract checks expected by the Streamlit workflow."""
    errors: list[str] = []
    warnings: list[str] = []
    summary: Dict[str, Any] = {}

    required_functions = [
        "get_supported_databases",
        "get_database_metadata",
        "get_required_identifiers",
        "get_preferred_identifier_order",
        "get_best_available_identifier",
        "get_available_contexts",
        "search_contexts",
        "get_context_by_id",
        "get_context_by_name",
        "search_compound_in_database",
        "probe_all_databases",
        "get_database_retrieval_capabilities",
        "get_database_retrieval_capability",
        "list_lotus_contexts",
        "search_lotus_contexts_online",
        "lotus_compound_in_taxon",
        "get_lotus_occurrences",
        "search_uniprot_taxa_online",
        "get_rhea_reactions_for_compound",
        "get_rhea_reactions_for_smiles",
        "_get_chebi_ids_for_smiles",
        "_rhea_reactions_for_chebi_ids",
        "get_uniprot_taxon_rhea_reactions",
        "build_rhea_uniprot_enzyme_rows",
        "search_npatlas_taxa_online",
        "get_npatlas_occurrences",
        "npatlas_compound_in_taxon",
        "search_npatlas_compound",
        "_npatlas_taxon_has_producers",
    ]
    missing_functions = [
        function_name
        for function_name in required_functions
        if not callable(globals().get(function_name))
    ]
    if missing_functions:
        errors.append(f"Missing or non-callable functions: {', '.join(missing_functions)}")
    summary["required_functions_checked"] = len(required_functions)
    summary["missing_functions"] = missing_functions

    expected_databases = ["LOTUS", "HMDB", "Rhea + UniProt", "NPAtlas"]
    try:
        supported = get_supported_databases()
        summary["supported_databases"] = supported
        if supported != expected_databases:
            errors.append(
                f"Supported databases mismatch: expected {expected_databases}, got {supported}"
            )
    except Exception as exc:
        supported = []
        errors.append(f"get_supported_databases failed: {exc}")

    try:
        metadata = get_database_metadata()
        missing_metadata = [
            database_name for database_name in expected_databases if database_name not in metadata
        ]
        summary["metadata_databases"] = sorted(metadata.keys())
        if missing_metadata:
            errors.append(f"Missing metadata for databases: {', '.join(missing_metadata)}")
    except Exception as exc:
        metadata = {}
        errors.append(f"get_database_metadata failed: {exc}")

    try:
        capabilities = get_database_retrieval_capabilities()
        missing_capabilities = [
            database_name for database_name in expected_databases if database_name not in capabilities
        ]
        summary["retrieval_capability_databases"] = sorted(capabilities.keys())
        if missing_capabilities:
            errors.append(
                f"Missing retrieval capabilities for databases: {', '.join(missing_capabilities)}"
            )
        valid_statuses = {"implemented", "partial", "not_implemented", "failed", "user_import"}
        for database_name in expected_databases:
            capability = capabilities.get(database_name, {})
            list_status = capability.get("context_list_status")
            search_status = capability.get("context_search_status")
            context_status = capability.get("context_retrieval_status")
            compound_status = capability.get("compound_lookup_status")
            if list_status is not None and list_status not in valid_statuses:
                errors.append(f"{database_name} has invalid context list status: {list_status}")
            if search_status is not None and search_status not in valid_statuses:
                errors.append(f"{database_name} has invalid context search status: {search_status}")
            if context_status not in valid_statuses:
                errors.append(f"{database_name} has invalid context retrieval status: {context_status}")
            if compound_status not in valid_statuses:
                errors.append(f"{database_name} has invalid compound lookup status: {compound_status}")
        for database_name in ["LOTUS", "Rhea + UniProt", "NPAtlas"]:
            frame = get_available_contexts(database_name, allow_online=False)
            if not frame.empty:
                errors.append(f"{database_name} returned contexts although real retrieval is not implemented.")
            if not normalize_text(frame.attrs.get("reason", "")):
                errors.append(f"{database_name} empty context dataframe is missing a retrieval reason.")
    except Exception as exc:
        errors.append(f"Retrieval capability contract check failed: {exc}")

    glucose_row = {
        "SMILES": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "InChIKey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
    }
    structure_row = {
        "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O",
        "InChIKey": "KZNIFHPLKGYRTM-UHFFFAOYSA-N",
    }
    identifier_examples = {
        "LOTUS": (structure_row, "SMILES"),
        "Rhea + UniProt": (glucose_row, "InChIKey"),
        "NPAtlas": (glucose_row, "InChIKey"),
    }
    identifier_summary: Dict[str, Any] = {}
    for database_name, example in identifier_examples.items():
        row, expected_identifier = example
        try:
            best = get_best_available_identifier(row, database_name)
            identifier_summary[database_name] = best
            if best.get("identifier_type") != expected_identifier:
                errors.append(
                    f"{database_name} identifier selection mismatch: expected "
                    f"{expected_identifier}, got {best.get('identifier_type')}"
                )
        except Exception as exc:
            errors.append(f"{database_name} identifier selection failed: {exc}")
    summary["identifier_selection"] = identifier_summary

    context_search_examples = {
        "LOTUS": ["arabidopsis", "salvia"],
        "Rhea + UniProt": ["arabidopsis", "salvia"],
        "NPAtlas": ["streptomyces", "penicillium"],
    }
    context_search_summary: Dict[str, int] = {}
    for database_name, queries in context_search_examples.items():
        database_rows = 0
        for query in queries:
            try:
                result = search_contexts(
                    database_name,
                    query,
                    limit=5,
                    allow_online=allow_online,
                )
                database_rows += len(result)
                if allow_online and result.empty:
                    warnings.append(
                        f"{database_name} context search for '{query}' returned no rows: "
                        f"{normalize_text(result.attrs.get('reason', 'No reason provided.'))}"
                    )
                if not allow_online and not result.empty:
                    errors.append(f"{database_name} returned query contexts while online search was disabled.")
            except Exception as exc:
                warnings.append(f"{database_name} context search for '{query}' failed: {exc}")
        context_search_summary[database_name] = database_rows
    summary["context_search_rows_total"] = context_search_summary

    if allow_online:
        try:
            for database_name, queries in context_search_examples.items():
                for query in queries:
                    result = search_contexts(database_name, query, limit=5, allow_online=True)
                    for _, row in result.iterrows():
                        source = normalize_text(row.get("source", "")).lower()
                        if any(marker in source for marker in ["local helper", "manual", "fake"]):
                            errors.append(
                                f"{database_name} context search returned non-retrieved source for '{query}': {source}"
                            )
        except Exception as exc:
            warnings.append(f"Retrieved-source validation failed: {exc}")

    try:
        probe_frame = probe_all_databases(allow_online=allow_online, verbose=True)
        summary["probe_all_databases_rows"] = len(probe_frame)
        if probe_frame.empty:
            warnings.append("probe_all_databases returned an empty DataFrame.")
    except TypeError as exc:
        errors.append(f"probe_all_databases verbose contract failed: {exc}")
    except Exception as exc:
        warnings.append(f"probe_all_databases failed during online probe: {exc}")

    expected_lookup_keys = {
        "database",
        "status",
        "reason",
        "identifier_type",
        "identifier_value",
        "context_id",
        "context_name",
        "context_type",
        "matched_identifier",
        "source",
        "extra",
    }
    allowed_lookup_statuses = {
        "found",
        "not_found",
        "not_evaluated",
        "failed",
        "not_implemented",
    }
    lookup_contract_examples = {
        "LOTUS": (
            structure_row,
            {"context_id": "user_defined_taxon", "context_name": "User-defined taxon", "context_type": "taxon_or_species"},
        ),
        "Rhea + UniProt": (
            glucose_row,
            {"context_id": "3702", "context_name": "Arabidopsis thaliana", "context_type": "taxon"},
        ),
        "NPAtlas": (
            glucose_row,
            {"context_id": "", "context_name": "Streptomyces", "context_type": "taxon"},
        ),
    }
    lookup_status_summary: Dict[str, str] = {}
    try:
        for database_name, example in lookup_contract_examples.items():
            row, context_selection = example
            lookup_result = search_compound_in_database(
                database_name,
                row,
                context_selection,
                allow_online=False,
            )
            missing_lookup_keys = sorted(expected_lookup_keys.difference(lookup_result.keys()))
            if missing_lookup_keys:
                errors.append(
                    f"{database_name} lookup result is missing keys: {', '.join(missing_lookup_keys)}"
                )
            status = normalize_text(lookup_result.get("status", ""))
            lookup_status_summary[database_name] = status
            if status not in allowed_lookup_statuses:
                errors.append(f"{database_name} lookup returned invalid status: {status}")
        summary["standard_lookup_keys"] = sorted(expected_lookup_keys)
        summary["lookup_statuses"] = lookup_status_summary
    except Exception as exc:
        errors.append(f"search_compound_in_database standard result check failed: {exc}")

    return {
        "ok": not errors,
        "errors": errors,
        "warnings": warnings,
        "summary": summary,
    }


def run_demo_tests() -> None:
    """Print compact diagnostics for the catalog layer."""
    print("Supported databases:", get_supported_databases())
    print("Metadata keys:", sorted(get_database_metadata().keys()))

    glucose_row = {
        "Name": "D-Glucose",
        "PubChem ID": "5793",
        "SMILES": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "InChIKey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
        "Formula": "C6H12O6",
        "Mapping_status": "mapped",
    }
    structure_row = {
        "Name": "Apigenin",
        "PubChem ID": "5280443",
        "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(O2)C=C(C=C3O)O)O",
        "InChIKey": "KZNIFHPLKGYRTM-UHFFFAOYSA-N",
        "Formula": "C15H10O5",
        "Mapping_status": "partial",
    }
    example_rows = {
        "LOTUS": structure_row,
        "Rhea + UniProt": glucose_row,
        "NPAtlas": glucose_row,
    }

    print("Best identifier selection:")
    for db_name, row in example_rows.items():
        best = get_best_available_identifier(row, db_name)
        print(
            f"  {db_name}: {best['identifier_type']}={best['identifier_value']} ({best['status']})"
        )

    print("Database online feasibility probes:")
    feasibility = probe_all_databases(allow_online=True)
    if feasibility.empty:
        print("  No probe results.")
    else:
        print(feasibility.to_string(index=False))

    print("Compound lookup checks:")
    lotus_result = search_compound_in_database(
        "LOTUS",
        structure_row,
        {"context_id": "hsa", "context_name": "Homo sapiens", "context_type": "taxon_or_species"},
    )
    print("  LOTUS with context:", lotus_result)

    print("Database catalog contract validation:")
    contract = validate_database_catalog_contract(allow_online=True)
    print(json.dumps(contract, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    run_demo_tests()
