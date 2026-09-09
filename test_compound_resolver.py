# -*- coding: utf-8 -*-
"""Compound resolver for Omiscera.

This script reads a CSV input table, resolves compound identifiers
through PubChem using only essential endpoints,
and writes a result table for testing.

The functions in this module are designed to be reusable and importable
from other modules, such as the future app integration.
"""

import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

try:
    import requests
except ImportError:
    requests = None

# Official API endpoints (constants)
PUBCHEM_NAME_CIDS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/TXT"
PUBCHEM_INCHIKEY_CIDS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{encoded_inchikey}/cids/TXT"
PUBCHEM_SMILES_CIDS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/cids/TXT"
PUBCHEM_TEXT_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={encoded_query}&retmode=json&retmax=1"
PUBCHEM_CID_PROPERTIES_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
PUBCHEM_CID_SYNONYMS_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
PUBCHEM_PUGVIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
# PubChem's "substance" (SID) namespace is DIFFERENT from its "compound" (CID)
# namespace - the two number lines overlap but refer to unrelated entities.
# Needed because KEGG's own DBLINKS "PubChem:" value is an SID, not a CID (a
# legacy KEGG-NIH submission quirk, confirmed live 2026-08-26: KEGG C00031
# (D-glucose) reports "PubChem: 3333", which is SID 3333 - resolving it via
# THIS endpoint correctly gives CID 5793 (real D-glucose), while treating
# 3333 as a CID directly would silently point at an unrelated drug molecule.
PUBCHEM_SID_TO_CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/TXT"
KEGG_FIND_COMPOUND_URL = "https://rest.kegg.jp/find/compound/{encoded_name}"
KEGG_GET_COMPOUND_URL = "https://rest.kegg.jp/get/{kegg_id}"
# KEGG's documented per-request cap on the "get" endpoint - used to batch
# reaction-equation fetches in database_catalog.py rather than one call per
# reaction (a compound's Role column can otherwise need ~10x as many calls).
KEGG_GET_BATCH_SIZE = 10

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache")
PUBCHEM_CACHE_PATH = os.path.join(CACHE_DIR, "pubchem_cache.csv")
KEGG_CACHE_PATH = os.path.join(CACHE_DIR, "kegg_cache.csv")
PUBCHEM_SYNONYMS_CACHE_PATH = os.path.join(CACHE_DIR, "pubchem_synonyms_cache.csv")
HMDB_CACHE_PATH = os.path.join(CACHE_DIR, "hmdb_pubchem_cache.csv")
FINAL_RESULT_CACHE_PATH = os.path.join(CACHE_DIR, "final_result_cache.csv")

PUBCHEM_CACHE = {}
KEGG_CACHE = {}
PUBCHEM_SYNONYMS_CACHE = {}
HMDB_CACHE = {}
FINAL_RESULT_CACHE = {}
CACHE_LOADED = False
CURRENT_ROW_REQUEST_FAILED = False
RUNTIME_DIAGNOSTICS = {
    "final_cache_hits": 0,
    "api_cache_hits": 0,
    "online_requests": 0,
    "pubchem_online_requests": 0,
    "kegg_online_requests": 0,
    "cache_misses": 0,
    "retry_count": 0,
    "total_sleep_time": 0.0,
    "compounds_resolved_online": 0,
    "compounds_resolved_from_final_cache": 0,
}
HTTP_SESSION = requests.Session() if requests is not None else None
OUTPUT_COLUMNS = [
    "Name",
    "KEGG ID",
    "PubChem ID",
    "SMILES",
    "InChIKey",
    "HMDB ID",
    "Formula",
    "Mapping_status",
]


def ensure_cache_loaded():
    """Load lookup caches once per run."""
    global CACHE_LOADED
    if CACHE_LOADED:
        return

    os.makedirs(CACHE_DIR, exist_ok=True)

    if os.path.exists(PUBCHEM_CACHE_PATH):
        df = pd.read_csv(PUBCHEM_CACHE_PATH, dtype=str).fillna("")
        for _, row in df.iterrows():
            cache_key = (row.get("CacheType", ""), row.get("CacheKey", ""))
            PUBCHEM_CACHE[cache_key] = {
                "PubChem ID": row.get("PubChem ID", ""),
                "SMILES": row.get("SMILES", ""),
                "InChIKey": row.get("InChIKey", ""),
                "Formula": row.get("Formula", ""),
            }

    if os.path.exists(KEGG_CACHE_PATH):
        df = pd.read_csv(KEGG_CACHE_PATH, dtype=str).fillna("")
        for _, row in df.iterrows():
            cache_key = (row.get("CacheType", ""), row.get("CacheKey", ""))
            KEGG_CACHE[cache_key] = {
                "KEGG ID": row.get("KEGG ID", ""),
                "Formula": row.get("Formula", ""),
            }

    if os.path.exists(PUBCHEM_SYNONYMS_CACHE_PATH):
        df = pd.read_csv(PUBCHEM_SYNONYMS_CACHE_PATH, dtype=str).fillna("")
        for _, row in df.iterrows():
            cid = clean_identifier(row.get("PubChem ID", ""))
            try:
                PUBCHEM_SYNONYMS_CACHE[cid] = json.loads(row.get("Synonyms", "[]"))
            except Exception:
                PUBCHEM_SYNONYMS_CACHE[cid] = []

    if os.path.exists(HMDB_CACHE_PATH):
        df = pd.read_csv(HMDB_CACHE_PATH, dtype=str).fillna("")
        for _, row in df.iterrows():
            cid = clean_identifier(row.get("PubChem ID", ""))
            if cid:
                HMDB_CACHE[cid] = normalize_hmdb_ids(row.get("HMDB ID", ""))

    if os.path.exists(FINAL_RESULT_CACHE_PATH):
        df = pd.read_csv(FINAL_RESULT_CACHE_PATH, dtype=str).fillna("")
        for _, row in df.iterrows():
            cache_key = str(row.get("CacheKey", "")).strip()
            if cache_key:
                FINAL_RESULT_CACHE[cache_key] = {
                    col: normalize_hmdb_ids(row.get(col, "")) if col == "HMDB ID" else row.get(col, "")
                    for col in OUTPUT_COLUMNS
                }

    CACHE_LOADED = True


def reset_runtime_diagnostics():
    """Reset runtime counters for one resolver run."""
    for key in RUNTIME_DIAGNOSTICS:
        RUNTIME_DIAGNOSTICS[key] = 0


def print_runtime_diagnostics(total_runtime: float = 0.0):
    """Print cache/request counters for verbose runs."""
    print(f"total_runtime: {total_runtime:.2f}")
    for key in [
        "final_cache_hits",
        "api_cache_hits",
        "online_requests",
        "pubchem_online_requests",
        "kegg_online_requests",
        "cache_misses",
        "retry_count",
        "total_sleep_time",
        "compounds_resolved_online",
        "compounds_resolved_from_final_cache",
    ]:
        print(f"{key}: {RUNTIME_DIAGNOSTICS.get(key, 0)}")


def save_pubchem_cache():
    """Persist PubChem cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for (cache_type, cache_key), data in PUBCHEM_CACHE.items():
        rows.append({
            "CacheType": cache_type,
            "CacheKey": cache_key,
            "PubChem ID": data.get("PubChem ID", ""),
            "SMILES": data.get("SMILES", ""),
            "InChIKey": data.get("InChIKey", ""),
            "Formula": data.get("Formula", ""),
        })
    pd.DataFrame(rows).to_csv(PUBCHEM_CACHE_PATH, index=False)


def save_kegg_cache():
    """Persist KEGG cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for (cache_type, cache_key), data in KEGG_CACHE.items():
        rows.append({
            "CacheType": cache_type,
            "CacheKey": cache_key,
            "KEGG ID": data.get("KEGG ID", ""),
            "Formula": data.get("Formula", ""),
        })
    pd.DataFrame(rows).to_csv(KEGG_CACHE_PATH, index=False)


def save_pubchem_synonyms_cache():
    """Persist PubChem synonym cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for cid, synonyms in PUBCHEM_SYNONYMS_CACHE.items():
        rows.append({"PubChem ID": cid, "Synonyms": json.dumps(synonyms)})
    pd.DataFrame(rows).to_csv(PUBCHEM_SYNONYMS_CACHE_PATH, index=False)


def save_hmdb_cache():
    """Persist HMDB IDs recovered from PubChem by CID."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for cid, hmdb_ids in HMDB_CACHE.items():
        rows.append({"PubChem ID": cid, "HMDB ID": normalize_hmdb_ids(hmdb_ids)})
    pd.DataFrame(rows).to_csv(HMDB_CACHE_PATH, index=False)


def save_final_result_cache():
    """Persist final per-input-name result cache to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    rows = []
    for cache_key, data in FINAL_RESULT_CACHE.items():
        row = {"CacheKey": cache_key}
        for col in OUTPUT_COLUMNS:
            row[col] = normalize_hmdb_ids(data.get(col, "")) if col == "HMDB ID" else data.get(col, "")
        rows.append(row)
    pd.DataFrame(rows).to_csv(FINAL_RESULT_CACHE_PATH, index=False)


def make_final_result_cache_key(name, row=None) -> str:
    """Normalize the original input name plus provided identifiers for final-result cache lookup."""
    cleaned = clean_wrapping_quotes(name)
    cleaned = re.sub(r"\s*;\s*", ";", cleaned)
    cleaned = re.sub(r" +", " ", cleaned.strip())
    base_key = cleaned.lower()
    if row is None:
        return base_key

    identifier_parts = []
    for column in ["KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"]:
        if column == "SMILES":
            value = clean_smiles_text(row.get(column, ""))
        elif column == "Formula":
            value = clean_identifier(row.get(column, ""))
        else:
            value = clean_identifier(row.get(column, ""))
        if value:
            identifier_parts.append(f"{column.lower()}={value.lower()}")

    if not identifier_parts:
        return base_key
    return base_key + "|" + "|".join(identifier_parts)


def get_final_result_cache(name, row=None) -> dict:
    ensure_cache_loaded()
    cache_key = make_final_result_cache_key(name, row)
    cached = FINAL_RESULT_CACHE.get(cache_key, {})
    if cached:
        RUNTIME_DIAGNOSTICS["final_cache_hits"] += 1
        return {col: cached.get(col, "") for col in OUTPUT_COLUMNS}
    RUNTIME_DIAGNOSTICS["cache_misses"] += 1
    return {}


def set_final_result_cache(name, data: dict, row=None, user_fields: dict = None):
    """Persist only identifiers the resolution process actually retrieved.

    Any field whose value exactly matches what the user typed in the input
    row is excluded — the cache must never become a copy of the user's own
    PubChem ID/SMILES/InChIKey/HMDB ID/Formula, only of what an online lookup
    genuinely returned. If nothing was retrieved (e.g. the row already had
    every identifier filled in), the row is not cached at all.
    """
    if not data:
        return
    mapping_status = str(data.get("Mapping_status", "")).strip()
    if mapping_status not in ("mapped", "partial", "not_found"):
        return
    if mapping_status == "not_found" and CURRENT_ROW_REQUEST_FAILED:
        return

    ensure_cache_loaded()
    cache_key = make_final_result_cache_key(name, row)
    if not cache_key:
        return

    user_fields = user_fields or {}
    retrieved = {}
    for column in ("KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"):
        if column == "HMDB ID":
            value = normalize_hmdb_ids(data.get(column, ""))
            user_value = normalize_hmdb_ids(user_fields.get(column, ""))
        else:
            value = cell_to_text(data.get(column, ""))
            user_value = cell_to_text(user_fields.get(column, ""))
        if user_value and value == user_value:
            continue  # supplied by the user, not retrieved — do not cache it
        retrieved[column] = value

    if not any(retrieved.values()):
        return  # nothing genuinely retrieved — avoid caching an empty echo row

    FINAL_RESULT_CACHE[cache_key] = {
        "Name": name,
        "KEGG ID": retrieved.get("KEGG ID", ""),
        "PubChem ID": retrieved.get("PubChem ID", ""),
        "SMILES": retrieved.get("SMILES", ""),
        "InChIKey": retrieved.get("InChIKey", ""),
        "HMDB ID": retrieved.get("HMDB ID", ""),
        "Formula": retrieved.get("Formula", ""),
        "Mapping_status": mapping_status,
    }


def get_pubchem_cache(cache_type: str, cache_key: str) -> dict:
    ensure_cache_loaded()
    cached = PUBCHEM_CACHE.get((cache_type, cache_key), {})
    if cached:
        RUNTIME_DIAGNOSTICS["api_cache_hits"] += 1
    else:
        RUNTIME_DIAGNOSTICS["cache_misses"] += 1
    return cached


def set_pubchem_cache(cache_type: str, cache_key: str, data: dict):
    if not data or not data.get("PubChem ID"):
        return
    ensure_cache_loaded()
    PUBCHEM_CACHE[(cache_type, cache_key)] = {
        "PubChem ID": data.get("PubChem ID", ""),
        "SMILES": data.get("SMILES", ""),
        "InChIKey": data.get("InChIKey", ""),
        "Formula": data.get("Formula", ""),
    }



def get_kegg_cache(cache_type: str, cache_key: str) -> dict:
    ensure_cache_loaded()
    cached = KEGG_CACHE.get((cache_type, cache_key), {})
    if cached:
        RUNTIME_DIAGNOSTICS["api_cache_hits"] += 1
    else:
        RUNTIME_DIAGNOSTICS["cache_misses"] += 1
    return cached


def set_kegg_cache(cache_type: str, cache_key: str, data: dict):
    if not data or not data.get("KEGG ID"):
        return
    ensure_cache_loaded()
    KEGG_CACHE[(cache_type, cache_key)] = {
        "KEGG ID": data.get("KEGG ID", ""),
        "Formula": data.get("formula", data.get("Formula", "")),
    }


def safe_http_request(
    url: str,
    timeout: int = 15,
    verbose: bool = False,
    database: str = "",
    query: str = "",
    attempts: int = 3,
) -> tuple:
    """Perform a GET request with retry. Returns (response_text, success_bool)."""
    global CURRENT_ROW_REQUEST_FAILED
    retry_statuses = {429, 500, 502, 503, 504}
    last_error = ""

    for attempt in range(1, attempts + 1):
        RUNTIME_DIAGNOSTICS["online_requests"] += 1
        database_lower = str(database).lower()
        if "pubchem" in database_lower:
            RUNTIME_DIAGNOSTICS["pubchem_online_requests"] += 1
        elif "kegg" in database_lower:
            RUNTIME_DIAGNOSTICS["kegg_online_requests"] += 1

        try:
            if HTTP_SESSION is not None:
                response = HTTP_SESSION.get(url, timeout=timeout, headers={"User-Agent": "Mozilla/5.0"})
                if response.status_code == 404:
                    return "", True
                if response.status_code in retry_statuses:
                    last_error = f"HTTP {response.status_code}: {response.reason}"
                    if attempt < attempts:
                        RUNTIME_DIAGNOSTICS["retry_count"] += 1
                        sleep_seconds = 0.5
                        RUNTIME_DIAGNOSTICS["total_sleep_time"] += sleep_seconds
                        time.sleep(sleep_seconds)
                        continue
                    break
                if response.status_code >= 400:
                    last_error = f"HTTP {response.status_code}: {response.reason}"
                    break
                return response.text, True

            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
                return text, True
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code}: {exc.reason}"
            if exc.code == 404:
                return "", True
            if exc.code in retry_statuses and attempt < attempts:
                RUNTIME_DIAGNOSTICS["retry_count"] += 1
                sleep_seconds = 0.5
                RUNTIME_DIAGNOSTICS["total_sleep_time"] += sleep_seconds
                time.sleep(sleep_seconds)
                continue
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = str(exc)
            break
        except Exception as exc:
            last_error = str(exc)
            break

    if verbose:
        print(f"Request failed | database={database or 'unknown'} | query={query} | error={last_error}")
    CURRENT_ROW_REQUEST_FAILED = True
    return "", False


def load_input_table(path: str, verbose: bool = False) -> pd.DataFrame:
    """Load the input compound table using comma or semicolon separator."""
    def read_header(file_path: str, encoding: str) -> str:
        with open(file_path, "r", encoding=encoding, errors="replace") as handle:
            return handle.readline()

    def detect_separator(header: str) -> str:
        return ";" if header.count(";") > header.count(",") else ","

    def read_single_name_column(file_path: str, encoding: str) -> pd.DataFrame:
        with open(file_path, "r", encoding=encoding, errors="replace") as handle:
            lines = handle.read().splitlines()
        return pd.DataFrame({"Name": lines[1:]})

    try:
        header = read_header(path, "utf-8")
        if header.strip().lower().replace("\ufeff", "") == "name":
            df = read_single_name_column(path, "utf-8")
        else:
            separator = detect_separator(header)
            df = pd.read_csv(path, sep=separator, encoding="utf-8", dtype=str)
    except UnicodeDecodeError:
        header = read_header(path, "latin1")
        if header.strip().lower().replace("\ufeff", "") == "name":
            df = read_single_name_column(path, "latin1")
        else:
            separator = detect_separator(header)
            df = pd.read_csv(path, sep=separator, encoding="latin1", dtype=str)

    # Normalize columns to standard field names
    rename_map = {}
    hmdb_source = ""
    hmdb_priority = ["hmdbid", "hmdb", "hdmbid", "hdmb"]
    normalized_columns = {}
    for column in df.columns:
        normalized = (
            str(column)
            .replace("\ufeff", "")
            .replace("\xa0", " ")
            .replace("\u2007", " ")
            .replace("\u202f", " ")
            .strip()
        )
        lower = normalized.lower().replace(" ", "").replace("_", "").replace("-", "")
        normalized_columns[column] = lower

        if lower == "name":
            rename_map[column] = "Name"
        elif lower in ("keggid", "kegg"):
            rename_map[column] = "KEGG ID"
        elif lower in ("pubchemid", "pubchem"):
            rename_map[column] = "PubChem ID"
        elif lower == "smiles":
            rename_map[column] = "SMILES"
        elif lower in ("inchikey",):
            rename_map[column] = "InChIKey"
        elif lower == "formula":
            rename_map[column] = "Formula"

    for hmdb_name in hmdb_priority:
        for column, lower in normalized_columns.items():
            if lower == hmdb_name:
                hmdb_source = column
                break
        if hmdb_source:
            break

    df = df.rename(columns=rename_map)
    if hmdb_source:
        df["HMDB ID"] = df[hmdb_source]

    if verbose:
        print(f"Detected columns: {list(df.columns)}")

    if "Name" not in df.columns:
        raise ValueError("Input table is missing the required 'Name' column.")

    df["Name"] = df["Name"].apply(clean_wrapping_quotes)

    name_preview = df["Name"].head(5).fillna("").tolist()
    if verbose:
        print(f"First 5 Name values: {name_preview}")

    if df["Name"].fillna("").astype(str).str.strip().eq("").all():
        raise ValueError("Input table 'Name' column is empty for all rows.")

    # Ensure all expected columns are present
    for col in ["Name", "KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"]:
        if col not in df.columns:
            df[col] = ""

    for col in ["KEGG ID", "PubChem ID", "InChIKey", "Formula"]:
        df[col] = df[col].apply(clean_identifier)

    df["HMDB ID"] = df["HMDB ID"].apply(normalize_hmdb_ids)
    df["SMILES"] = df["SMILES"].apply(clean_smiles_text)

    return df


def cell_to_text(value, strip: bool = True) -> str:
    """Convert a dataframe cell to text while keeping missing values empty."""
    if pd.isna(value):
        return ""
    text = str(value)
    return text.strip() if strip else text


def clean_identifier(value) -> str:
    """Clean identifier cells while preserving non-numeric IDs."""
    if pd.isna(value):
        return ""

    text = str(value).replace("\xa0", " ").replace("\u2007", " ").replace("\u202f", " ").strip()
    if not text or text.lower() in ("nan", "none"):
        return ""

    if re.fullmatch(r"\d+\.0+", text):
        return text.split(".", 1)[0]

    return text


def canonicalize_hmdb_id(value) -> str:
    """Convert a single HMDB identifier to HMDB plus seven zero-padded digits."""
    match = re.fullmatch(r"\s*HMDB\s*(\d+)\s*", str(value or ""), flags=re.IGNORECASE)
    if not match:
        return ""
    number_text = match.group(1).lstrip("0") or "0"
    return f"HMDB{int(number_text):07d}"


def normalize_hmdb_ids(value) -> str:
    """Normalize and deduplicate HMDB identifiers while preserving stable order."""
    if pd.isna(value):
        return ""

    normalized_ids = []
    seen = set()
    for match in re.findall(r"\bHMDB\s*\d+\b", str(value), flags=re.IGNORECASE):
        canonical_id = canonicalize_hmdb_id(match)
        if canonical_id and canonical_id not in seen:
            seen.add(canonical_id)
            normalized_ids.append(canonical_id)
    return ";".join(normalized_ids)


def clean_smiles_text(value) -> str:
    """Light cleanup for SMILES strings without altering chemistry characters."""
    if pd.isna(value):
        return ""

    text = str(value).replace("\xa0", " ").replace("\u2007", " ").replace("\u202f", " ").strip()
    if text.lower() in ("nan", "none"):
        return ""
    return text


def clean_wrapping_quotes(value) -> str:
    """Remove accidental matching quotes around a name without touching internal quotes."""
    if value is None or pd.isna(value):
        return ""

    text = str(value).replace("\xa0", " ").replace("\u2007", " ").replace("\u202f", " ").strip()
    if not text or text.lower() in ("nan", "none"):
        return ""

    while len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        text = text[1:-1].strip()

    return text


def get_query_names(name) -> list:
    """Return cleaned semicolon-separated query names in input order."""
    name_text = clean_wrapping_quotes(name)
    if not name_text:
        return []

    parts = [clean_wrapping_quotes(part) for part in name_text.split(";")]
    return [part.strip() for part in parts if part.strip()]


def normalize_for_comparison(value) -> str:
    """Normalize values enough to compare user-provided and online identifiers."""
    return re.sub(r"\s+", "", str(value or "").strip()).lower()


def keep_user_value_or_fill(current_value, online_value, field_name: str, source: str, compound_name: str, verbose: bool):
    """Keep non-empty user values; fill only empty fields from online data."""
    current_text = cell_to_text(current_value)
    online_text = cell_to_text(online_value)
    if not current_text:
        return online_text
    if online_text and normalize_for_comparison(current_text) != normalize_for_comparison(online_text) and verbose:
        print(f"[user value preserved] {compound_name} | {field_name} | source={source} | ignored={online_text}")
    return current_text


def get_user_provided_fields(row: pd.Series, input_name: str) -> dict:
    """Read and normalize user-provided fields from the current input row."""
    return {
        "Name": input_name,
        "KEGG ID": clean_identifier(row.get("KEGG ID", "")),
        "PubChem ID": clean_identifier(row.get("PubChem ID", "")),
        "SMILES": clean_smiles_text(row.get("SMILES", "")),
        "InChIKey": clean_identifier(row.get("InChIKey", "")),
        "HMDB ID": normalize_hmdb_ids(row.get("HMDB ID", "")),
        "Formula": clean_identifier(row.get("Formula", "")),
    }


def calculate_mapping_status(result: dict) -> str:
    """Calculate Mapping_status from the final merged identifier set.

    "mapped" requires KEGG alongside the four identifiers already required
    (PubChem + SMILES + InChIKey + HMDB) - a deliberate change from before
    KEGG was re-added, when HMDB alone filled this slot. This does NOT
    exclude anything from the analysis (API.py passes both "mapped" and
    "partial" rows through to stage 2) - it only changes which label and
    which displayed count a compound gets.
    """
    has_kegg = bool(cell_to_text(result.get("KEGG ID", "")))
    has_pubchem = bool(cell_to_text(result.get("PubChem ID", "")))
    has_smiles = bool(cell_to_text(result.get("SMILES", "")))
    has_inchikey = bool(cell_to_text(result.get("InChIKey", "")))
    has_hmdb = bool(cell_to_text(result.get("HMDB ID", "")))
    useful_identifiers = [
        has_kegg,
        has_pubchem,
        has_smiles,
        has_inchikey,
        has_hmdb,
        bool(cell_to_text(result.get("Formula", ""))),
    ]

    if has_kegg and has_pubchem and has_smiles and has_inchikey and has_hmdb:
        return "mapped"
    if any(useful_identifiers):
        return "partial"
    return "not_found"


def merge_cached_result_with_user_fields(cached_result: dict, user_fields: dict, verbose: bool, source_label: str) -> dict:
    """Merge cached data with current row fields while preserving user-provided values."""
    merged = {col: "" for col in OUTPUT_COLUMNS}
    merged["Name"] = user_fields.get("Name", cached_result.get("Name", ""))

    for field in ["KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"]:
        if field == "HMDB ID":
            user_value = normalize_hmdb_ids(user_fields.get(field, ""))
            cached_value = normalize_hmdb_ids(cached_result.get(field, ""))
        else:
            user_value = cell_to_text(user_fields.get(field, ""))
            cached_value = cell_to_text(cached_result.get(field, ""))
        if user_value:
            if cached_value and normalize_for_comparison(user_value) != normalize_for_comparison(cached_value) and verbose:
                print(
                    f"[user value preserved] {merged['Name']} | {field} | "
                    f"source={source_label} | ignored={cached_value}"
                )
            merged[field] = user_value
        else:
            merged[field] = cached_value

    merged["Mapping_status"] = calculate_mapping_status(merged)
    return merged


def split_compound_names(name_cell) -> list:
    """Split the Name column into alternative query names using semicolon."""
    return get_query_names(name_cell)


def clean_query_name(name) -> str:
    """Clean a compound name for database query use."""
    name_text = clean_wrapping_quotes(name)
    name_text = name_text.replace("\xa0", " ").replace("\u2007", " ").replace("\u202f", " ")
    name_text = re.sub(r" +", " ", name_text.strip())
    return name_text


def make_duplicate_name_key(name) -> str:
    """Build a normalized key for duplicate input Name detection."""
    synonyms = split_compound_names(name)
    if not synonyms:
        return clean_query_name(name).lower()

    cleaned_synonyms = [clean_query_name(synonym).lower() for synonym in synonyms]
    cleaned_synonyms = [synonym for synonym in cleaned_synonyms if synonym]
    return ";".join(sorted(cleaned_synonyms))


def remove_duplicate_input_names(input_df: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicated input Names, keeping the first occurrence."""
    if "Name" not in input_df.columns:
        return input_df

    input_df = input_df.copy()
    input_df["Name"] = input_df["Name"].apply(clean_wrapping_quotes)

    seen_keys = set()
    keep_indices = []
    removed_rows = []

    for row_number, (index, row) in enumerate(input_df.iterrows(), start=1):
        name_value = cell_to_text(row.get("Name", ""), strip=False)
        key = make_duplicate_name_key(name_value)
        if key and key in seen_keys:
            removed_rows.append((row_number, name_value))
        else:
            seen_keys.add(key)
            keep_indices.append(index)

    if removed_rows:
        print("Duplicate input Names detected. Keeping first occurrence and removing duplicated rows:")
        for row_number, name_value in removed_rows:
            print(f"Row {row_number}: {name_value}")

    return input_df.loc[keep_indices].reset_index(drop=True)


def generate_query_variants(name) -> list:
    """Generate internal database query variants without changing output names."""
    cleaned = clean_query_name(name)
    if not cleaned:
        return []

    variants = []

    def add_variant(value):
        value = clean_query_name(value)
        if value and value not in variants:
            variants.append(value)

    add_variant(cleaned)

    spaced_paren_match = re.match(r"^([A-Za-z]+)\s+\((\d+:\d+(?:[/_]\d+:\d+)?)\)$", cleaned)
    compact_paren_match = re.match(r"^([A-Za-z]+)\((\d+:\d+(?:[/_]\d+:\d+)?)\)$", cleaned)
    lipid_match = re.match(r"^([A-Za-z]+)\s+(\d+:\d+(?:[/_]\d+:\d+)?)$", cleaned)

    if spaced_paren_match is not None:
        match = spaced_paren_match
        lipid_class = match.group(1)
        chains = match.group(2)
        slash_chains = chains.replace("_", "/")
        underscore_chains = chains.replace("/", "_")
        space_chains = re.sub(r"[/_]", " ", chains)

        add_variant(f"{lipid_class}({chains})")
        add_variant(f"{lipid_class} ({chains})")
        add_variant(f"{lipid_class} {chains}")
        add_variant(f"{lipid_class} {slash_chains}")
        add_variant(f"{lipid_class}({slash_chains})")
        add_variant(f"{lipid_class} {underscore_chains}")
        add_variant(f"{lipid_class}({underscore_chains})")
        add_variant(f"{lipid_class} {space_chains}")
    elif compact_paren_match is not None:
        match = compact_paren_match
        lipid_class = match.group(1)
        chains = match.group(2)
        slash_chains = chains.replace("_", "/")
        underscore_chains = chains.replace("/", "_")
        space_chains = re.sub(r"[/_]", " ", chains)

        add_variant(f"{lipid_class}({chains})")
        add_variant(f"{lipid_class} ({chains})")
        add_variant(f"{lipid_class} {chains}")
        add_variant(f"{lipid_class} {slash_chains}")
        add_variant(f"{lipid_class}({slash_chains})")
        add_variant(f"{lipid_class} {underscore_chains}")
        add_variant(f"{lipid_class}({underscore_chains})")
        add_variant(f"{lipid_class} {space_chains}")
    elif lipid_match is not None:
        match = lipid_match
        lipid_class = match.group(1)
        chains = match.group(2)
        slash_chains = chains.replace("_", "/")
        underscore_chains = chains.replace("/", "_")
        space_chains = re.sub(r"[/_]", " ", chains)

        add_variant(f"{lipid_class}({chains})")
        add_variant(f"{lipid_class} ({chains})")
        add_variant(f"{lipid_class} {chains}")
        add_variant(f"{lipid_class} {slash_chains}")
        add_variant(f"{lipid_class}({slash_chains})")
        add_variant(f"{lipid_class} {underscore_chains}")
        add_variant(f"{lipid_class}({underscore_chains})")
        add_variant(f"{lipid_class} {space_chains}")
    else:
        if "/" in cleaned:
            add_variant(cleaned.replace("/", "_"))
        if "_" in cleaned:
            add_variant(cleaned.replace("_", "/"))

    return variants[:8]


def parse_first_cid(text: str) -> str:
    """Return the first CID from a plain text CID response."""
    lines = str(text or "").strip().split("\n")
    if lines and lines[0].strip().isdigit():
        return lines[0].strip()
    return ""


def get_pubchem_cid_by_inchikey(inchikey: str) -> str:
    """Search PubChem by InChIKey and return the first CID found."""
    if not inchikey or not str(inchikey).strip():
        return ""

    inchikey_text = clean_identifier(inchikey)
    cached = get_pubchem_cache("inchikey", inchikey_text)
    if cached.get("PubChem ID"):
        return cached.get("PubChem ID", "")

    encoded_inchikey = urllib.parse.quote(inchikey_text, safe="")
    url = PUBCHEM_INCHIKEY_CIDS_URL.format(encoded_inchikey=encoded_inchikey)
    text, success = safe_http_request(url, timeout=20, database="PubChem InChIKey", query=inchikey_text)
    if not success or not text:
        return ""

    cid = parse_first_cid(text)
    if cid:
        set_pubchem_cache("inchikey", inchikey_text, {"PubChem ID": cid})
    return cid


def get_pubchem_cid_by_smiles(smiles: str) -> str:
    """Search PubChem by SMILES and return the first CID found."""
    if not smiles or not str(smiles).strip():
        return ""

    smiles_text = clean_smiles_text(smiles)
    cached = get_pubchem_cache("smiles", smiles_text)
    if cached.get("PubChem ID"):
        return cached.get("PubChem ID", "")

    encoded_smiles = urllib.parse.quote(smiles_text, safe="")
    url = PUBCHEM_SMILES_CIDS_URL.format(encoded_smiles=encoded_smiles)
    text, success = safe_http_request(url, timeout=20, database="PubChem SMILES", query=smiles_text)
    if not success or not text:
        return ""

    cid = parse_first_cid(text)
    if cid:
        set_pubchem_cache("smiles", smiles_text, {"PubChem ID": cid})
    return cid


def get_pubchem_cid_by_text_search(query: str, verbose: bool = False) -> str:
    """Search PubChem broadly through NCBI E-utilities and return the first CID."""
    query_text = clean_query_name(query)
    if not query_text:
        return ""

    cached = get_pubchem_cache("text", query_text.lower())
    if cached.get("PubChem ID"):
        return cached.get("PubChem ID", "")

    encoded_query = urllib.parse.quote(query_text, safe="")
    url = PUBCHEM_TEXT_SEARCH_URL.format(encoded_query=encoded_query)
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="PubChem text", query=query_text)
    if not success or not text:
        return ""

    try:
        data = json.loads(text)
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if id_list and str(id_list[0]).strip().isdigit():
            return str(id_list[0]).strip()
    except Exception as exc:
        global CURRENT_ROW_REQUEST_FAILED
        CURRENT_ROW_REQUEST_FAILED = True
        if verbose:
            print(f"JSON parse failed | database=PubChem text | query={query_text} | error={exc}")
        return ""

    return ""


def get_pubchem_cid_by_name_direct(name: str, verbose: bool = False) -> str:
    """Search PubChem by exact name route and return the first CID found."""
    query_text = clean_query_name(name)
    if not query_text:
        return ""

    cached = get_pubchem_cache("name", query_text.lower())
    if cached.get("PubChem ID"):
        return cached.get("PubChem ID", "")

    encoded_name = urllib.parse.quote(query_text, safe="")
    if verbose and "/" in query_text:
        print(f"Original query: {query_text}")
        print(f"Encoded query: {encoded_name}")
    url = PUBCHEM_NAME_CIDS_URL.format(encoded_name=encoded_name)
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="PubChem name", query=query_text)
    if not success or not text:
        return ""

    cid = parse_first_cid(text)
    if cid:
        set_pubchem_cache("name", query_text.lower(), {"PubChem ID": cid})
    return cid


def normalize_name_for_match(value) -> str:
    """Normalize names for cautious PubChem synonym validation."""
    if pd.isna(value):
        return ""

    text = str(value).replace("\xa0", " ").replace("\u2007", " ").replace("\u202f", " ")
    text = re.sub(r" +", " ", text.strip()).lower()
    text = text.replace("_", "/")
    text = re.sub(r"\s*\(\s*", "(", text)
    text = re.sub(r"\s*\)\s*", ")", text)

    lipid_match = re.match(r"^([a-z]+)\s*\((\d+:\d+(?:/\d+:\d+)?)\)$", text)
    if lipid_match is not None:
        return f"{lipid_match.group(1)} {lipid_match.group(2)}"

    text = re.sub(r"^([a-z]+)\s+(\d+:\d+(?:/\d+:\d+)?)$", r"\1 \2", text)
    return text


def get_pubchem_synonyms_by_cid(cid: str, verbose: bool = False) -> list:
    """Retrieve PubChem synonyms for a CID using PUG-REST."""
    if not cid or not str(cid).strip():
        return []

    cid_text = clean_identifier(cid)
    ensure_cache_loaded()
    if cid_text in PUBCHEM_SYNONYMS_CACHE:
        RUNTIME_DIAGNOSTICS["api_cache_hits"] += 1
        return PUBCHEM_SYNONYMS_CACHE.get(cid_text, [])
    RUNTIME_DIAGNOSTICS["cache_misses"] += 1

    url = PUBCHEM_CID_SYNONYMS_URL.format(cid=urllib.parse.quote(cid_text))
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="PubChem synonyms", query=cid_text)
    if not success or not text:
        return []

    try:
        data = json.loads(text)
        info_list = data.get("InformationList", {}).get("Information", [])
        if not info_list:
            return []
        synonyms = info_list[0].get("Synonym", []) or []
        if synonyms:
            PUBCHEM_SYNONYMS_CACHE[cid_text] = synonyms
        return synonyms
    except Exception as exc:
        global CURRENT_ROW_REQUEST_FAILED
        CURRENT_ROW_REQUEST_FAILED = True
        if verbose:
            print(f"JSON parse failed | database=PubChem synonyms | query={cid_text} | error={exc}")
        return []


def collect_json_strings(value) -> list:
    """Collect all string values recursively from nested JSON-like data."""
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        strings = []
        for nested_value in value.values():
            strings.extend(collect_json_strings(nested_value))
        return strings
    if isinstance(value, list):
        strings = []
        for nested_value in value:
            strings.extend(collect_json_strings(nested_value))
        return strings
    return []


def extract_hmdb_ids_from_text_values(values) -> str:
    """Extract unique HMDB identifiers from PubChem text values."""
    hmdb_ids = []
    seen = set()
    for value in values or []:
        for match in re.findall(r"\bHMDB\s*\d+\b", str(value), flags=re.IGNORECASE):
            hmdb_id = canonicalize_hmdb_id(match)
            if hmdb_id not in seen:
                seen.add(hmdb_id)
                hmdb_ids.append(hmdb_id)
    return ";".join(hmdb_ids)


def fetch_pubchem_json_strings(url: str, database: str, query: str, verbose: bool = False) -> tuple:
    """Fetch a PubChem JSON endpoint and return all nested string values plus success flag."""
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database=database, query=query)
    if not success or not text:
        return [], False

    try:
        data = json.loads(text)
        return collect_json_strings(data), True
    except Exception as exc:
        global CURRENT_ROW_REQUEST_FAILED
        CURRENT_ROW_REQUEST_FAILED = True
        if verbose:
            print(f"JSON parse failed | database={database} | query={query} | error={exc}")
        return [], False


def get_hmdb_ids_by_pubchem_cid(cid: str, verbose: bool = False) -> str:
    """Recover HMDB identifiers from PubChem synonyms and PUG-View data for a CID."""
    cid_text = clean_identifier(cid)
    if not cid_text:
        return ""

    ensure_cache_loaded()
    if cid_text in HMDB_CACHE:
        RUNTIME_DIAGNOSTICS["api_cache_hits"] += 1
        if verbose:
            print(f"[hmdb cache hit] {cid_text}")
        return normalize_hmdb_ids(HMDB_CACHE.get(cid_text, ""))

    RUNTIME_DIAGNOSTICS["cache_misses"] += 1
    all_strings = []

    synonyms_url = PUBCHEM_CID_SYNONYMS_URL.format(cid=urllib.parse.quote(cid_text))
    if verbose:
        print(f"[hmdb pubchem synonyms lookup] {cid_text}")
    synonyms_strings, synonyms_success = fetch_pubchem_json_strings(
        synonyms_url,
        "PubChem synonyms",
        cid_text,
        verbose=verbose,
    )
    all_strings.extend(synonyms_strings)

    pugview_url = PUBCHEM_PUGVIEW_URL.format(cid=urllib.parse.quote(cid_text))
    if verbose:
        print(f"[hmdb pubchem pugview lookup] {cid_text}")
    pugview_strings, pugview_success = fetch_pubchem_json_strings(
        pugview_url,
        "PubChem PUG-View",
        cid_text,
        verbose=verbose,
    )
    all_strings.extend(pugview_strings)

    if verbose:
        print(f"[hmdb strings scanned] {len(all_strings)}")

    hmdb_ids = extract_hmdb_ids_from_text_values(all_strings)
    if hmdb_ids:
        HMDB_CACHE[cid_text] = normalize_hmdb_ids(hmdb_ids)
        if verbose:
            print(f"[hmdb found] {cid_text} -> {HMDB_CACHE[cid_text]}")
        return HMDB_CACHE[cid_text]

    if synonyms_success and pugview_success:
        HMDB_CACHE[cid_text] = ""

    if verbose:
        print(f"[hmdb not found] {cid_text}")
    return ""


def validate_pubchem_name_match(cid: str, query_variants: list, verbose: bool = False) -> tuple:
    """Return (accepted, matched_synonym, reason) for a name-based CID."""
    synonyms = get_pubchem_synonyms_by_cid(cid, verbose=verbose)
    if not synonyms:
        return False, "", "no PubChem synonyms returned"

    normalized_queries = {normalize_name_for_match(value) for value in query_variants}
    normalized_queries = {value for value in normalized_queries if value}

    for synonym in synonyms:
        normalized_synonym = normalize_name_for_match(synonym)
        if normalized_synonym in normalized_queries:
            return True, str(synonym), ""

    return False, "", "no query variant found in PubChem synonyms"


def get_pubchem_cid_by_name(name: str, verbose: bool = False) -> str:
    """Search PubChem by name variants, then broader text search."""
    if not name or not str(name).strip():
        return ""

    query_variants = generate_query_variants(name)
    for query_name in query_variants:
        cid = get_pubchem_cid_by_name_direct(query_name, verbose=verbose)
        if cid:
            if verbose:
                print(f"Candidate CID: {cid}")
                if clean_query_name(query_name) == clean_query_name(name):
                    print("PubChem match source: direct name lookup")
                else:
                    print("PubChem match source: direct variant lookup")
            return cid

        cid = get_pubchem_cid_by_text_search(query_name, verbose=verbose)
        if cid:
            accepted, matched_synonym, reason = validate_pubchem_name_match(cid, query_variants, verbose=verbose)
            if verbose:
                print(f"Candidate CID: {cid}")
                print(f"PubChem match source: broad fallback {'accepted' if accepted else 'rejected'}")
            if accepted:
                if verbose:
                    print("Validation: accepted")
                    print(f"Matched synonym if accepted: {matched_synonym}")
                set_pubchem_cache("text", clean_query_name(query_name).lower(), {"PubChem ID": cid})
                return cid
            if verbose:
                print("Validation: rejected")
                print(f"Reason if rejected: {reason}")

    return ""


def normalize_query_string(value: str) -> str:
    """Normalize a string for name matching and comparison."""
    if value is None:
        return ""
    value = str(value).lower().strip()
    # Replace whitespace and common separators with a single space
    value = re.sub(r"[\s\-_,()\[\]]+", " ", value)
    # Remove any non-alphanumeric characters except spaces
    value = re.sub(r"[^a-z0-9 ]+", "", value)
    # Collapse multiple spaces
    value = re.sub(r" +", " ", value)
    return value.strip()


def get_pubchem_data_by_cid(cid: str, verbose: bool = False) -> dict:
    """Retrieve PubChem compound properties for a given CID."""
    if not cid or not str(cid).strip():
        return {}

    cid_text = clean_identifier(cid)
    result = {
        "PubChem ID": cid_text,
        "SMILES": "",
        "InChIKey": "",
        "Formula": "",
    }

    cached = get_pubchem_cache("cid", cid_text)
    if cached.get("SMILES") or cached.get("InChIKey") or cached.get("Formula"):
        return {
            "PubChem ID": cached.get("PubChem ID", cid_text),
            "SMILES": cached.get("SMILES", ""),
            "InChIKey": cached.get("InChIKey", ""),
            "Formula": cached.get("Formula", ""),
        }

    url = PUBCHEM_CID_PROPERTIES_URL.format(cid=urllib.parse.quote(cid_text))
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="PubChem cid", query=cid_text)
    if not success or not text:
        return result

    try:
        data = json.loads(text)
        props_table = data.get("PropertyTable", {})
        properties_list = props_table.get("Properties", [])
        if not properties_list:
            return result

        props = properties_list[0]

        # Extract fields exactly as specified in the official endpoint
        result["Formula"] = str(props.get("MolecularFormula", "") or "").strip()

        for smiles_key in ("IsomericSMILES", "CanonicalSMILES", "SMILES", "ConnectivitySMILES"):
            selected_smiles = str(props.get(smiles_key, "") or "").strip()
            if selected_smiles:
                result["SMILES"] = selected_smiles
                break
        result["InChIKey"] = str(props.get("InChIKey", "") or "").strip()
        if result.get("SMILES") or result.get("InChIKey") or result.get("Formula"):
            set_pubchem_cache("cid", cid_text, result)

    except Exception as exc:
        global CURRENT_ROW_REQUEST_FAILED
        CURRENT_ROW_REQUEST_FAILED = True
        if verbose:
            print(f"JSON parse failed | database=PubChem cid | query={cid_text} | error={exc}")

    return result


def parse_kegg_get(text: str) -> dict:
    """Parse the raw KEGG get response for basic fields.

    Returns a dict with keys: KEGG ID, name, formula, pubchem_sid. The
    DBLINKS "PubChem:" value KEGG reports is a Substance ID (SID), not a
    Compound ID - see PUBCHEM_SID_TO_CID_URL's comment for why that
    distinction matters.
    """
    result = {"KEGG ID": "", "name": "", "formula": "", "pubchem_sid": ""}
    names = []
    in_dblinks = False

    for line in text.splitlines():
        if line.startswith("ENTRY"):
            parts = line.split()
            if len(parts) >= 2:
                result["KEGG ID"] = parts[1].strip()
            in_dblinks = False
        elif line.startswith("NAME"):
            name_line = line[4:].strip()
            names.append(name_line)
            in_dblinks = False
        elif line.startswith("FORMULA"):
            result["formula"] = line[7:].strip()
            in_dblinks = False
        elif line.startswith("DBLINKS"):
            in_dblinks = True
            rest = line[7:].strip()
            if rest.startswith("PubChem:"):
                result["pubchem_sid"] = rest.split(":", 1)[1].strip()
        elif line.startswith("///"):
            break
        elif line.startswith(" "):
            # continuation line - only DBLINKS continuations matter here,
            # everything else (NAME wraps, pathway/enzyme blocks, ...) is
            # outside what this parser reports.
            cont = line.strip()
            if in_dblinks and cont.startswith("PubChem:"):
                result["pubchem_sid"] = cont.split(":", 1)[1].strip()
            elif names and not in_dblinks:
                names[-1] += " " + cont
        else:
            in_dblinks = False

    if names:
        result["name"] = "; ".join(names)

    return result


def validate_kegg_id(kegg_id: str, verbose: bool = False) -> dict:
    """Validate a KEGG ID and extract basic info (name, formula, PubChem SID)."""
    if not kegg_id or not str(kegg_id).strip():
        return {}

    kegg_id_text = str(kegg_id).strip()
    cached = get_kegg_cache("id", kegg_id_text)
    if cached.get("KEGG ID"):
        return {
            "KEGG ID": cached.get("KEGG ID", kegg_id_text),
            "name": "",
            "formula": cached.get("Formula", ""),
            "pubchem_sid": "",
        }

    url = KEGG_GET_COMPOUND_URL.format(kegg_id=urllib.parse.quote(kegg_id_text))
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="KEGG id", query=kegg_id_text)
    if not success or not text or "ENTRY" not in text:
        return {}

    parsed = parse_kegg_get(text)
    if not parsed.get("KEGG ID"):
        parsed["KEGG ID"] = kegg_id_text
    set_kegg_cache("id", kegg_id_text, parsed)
    return parsed


def get_pubchem_cid_from_kegg_sid(sid: str, verbose: bool = False) -> str:
    """Resolve a PubChem Substance ID (as reported in a KEGG DBLINKS line)
    to its associated Compound ID, via PubChem's substance->CID endpoint.

    Must never be skipped in favour of treating the SID as if it were a
    CID directly - the two ID spaces overlap numerically but refer to
    unrelated entities (verified live 2026-08-26: KEGG C00031/D-glucose
    reports "PubChem: 3333"; SID 3333 correctly resolves here to CID 5793,
    the real D-glucose CID, while CID 3333 taken literally is an unrelated
    dihydropyridine drug).
    """
    sid_text = clean_identifier(sid)
    if not sid_text or not sid_text.isdigit():
        return ""

    cached = get_pubchem_cache("kegg_sid_to_cid", sid_text)
    if cached.get("PubChem ID"):
        return cached.get("PubChem ID", "")

    url = PUBCHEM_SID_TO_CID_URL.format(sid=sid_text)
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="PubChem SID", query=sid_text)
    if not success or not text:
        return ""

    cid = parse_first_cid(text)
    if cid:
        set_pubchem_cache("kegg_sid_to_cid", sid_text, {"PubChem ID": cid})
    return cid


def search_kegg_compound_by_name(name: str, resolved_pubchem_id: str = "", verbose: bool = False) -> dict:
    """Search KEGG compound entries by name and return a verified match.

    Every name-matched candidate is checked against `resolved_pubchem_id`
    (this row's own already-resolved PubChem CID, if any) via KEGG's own
    DBLINKS PubChem cross-reference - see get_pubchem_cid_from_kegg_sid for
    why that check requires an SID->CID resolution step rather than a
    direct comparison. This is what protects against KEGG's name search
    matching a similarly-named but structurally different compound.

    Falls back to the single best unverified name match only when nothing
    could be structurally verified (e.g. no PubChem CID was resolved for
    this compound at all, or the matched KEGG entry has no PubChem DBLINK).
    """
    if not name or not name.strip():
        return {}

    query_name = clean_query_name(name)
    query_key = query_name.lower()
    cached = get_kegg_cache("name", query_key)
    if cached.get("KEGG ID"):
        return {"KEGG ID": cached.get("KEGG ID", ""), "name": ""}

    encoded_name = urllib.parse.quote(query_name, safe="")
    if verbose and "/" in query_name:
        print(f"Original query: {query_name}")
        print(f"Encoded query: {encoded_name}")
    url = KEGG_FIND_COMPOUND_URL.format(encoded_name=encoded_name)
    text, success = safe_http_request(url, timeout=20, verbose=verbose, database="KEGG name", query=query_key)
    if not success or not text:
        return {}

    candidates = []
    for line in text.splitlines():
        line = line.strip()
        if not line or "	" not in line:
            continue
        parts = line.split("	", 1)
        if len(parts) == 2:
            candidates.append({"KEGG ID": parts[0].strip(), "name": parts[1].strip()})

    if not candidates:
        return {}

    # Exact-name candidates first (by normalized string), same ranking the
    # unverified fallback below uses - but here we may check several.
    query_norm = normalize_query_string(query_name)
    exact_candidates = []
    for candidate in candidates:
        candidate_names = [n.strip() for n in candidate["name"].split(";") if n.strip()]
        if any(normalize_query_string(n) == query_norm for n in candidate_names):
            exact_candidates.append(candidate)
    ranked = exact_candidates or candidates

    resolved_pubchem_id_text = clean_identifier(resolved_pubchem_id)
    if resolved_pubchem_id_text:
        for candidate in ranked:
            info = validate_kegg_id(candidate["KEGG ID"], verbose=verbose)
            sid = info.get("pubchem_sid", "")
            if not sid:
                continue
            candidate_cid = get_pubchem_cid_from_kegg_sid(sid, verbose=verbose)
            if candidate_cid and candidate_cid == resolved_pubchem_id_text:
                result = {"KEGG ID": candidate["KEGG ID"], "name": candidate["name"]}
                set_kegg_cache("name", query_key, result)
                if verbose:
                    print(f"KEGG match verified via PubChem CID {candidate_cid}: {candidate['KEGG ID']}")
                return result

    # Nothing could be structurally verified - fall back to the best
    # unverified name match (ranked[0]), matching the original behaviour.
    best = ranked[0]
    result = {"KEGG ID": best["KEGG ID"], "name": best["name"]}
    set_kegg_cache("name", query_key, result)
    return result


def resolve_compound(row: pd.Series, verbose: bool = False) -> dict:
    """Resolve a single compound row and return the mapped result.

    Strategy:
    1. Preserve user-provided identifiers
    2. Use PubChem identifiers to recover missing PubChem fields
    3. Search PubChem by name only for missing IDs
    4. If no KEGG ID was provided, search KEGG by name, VERIFIED against the
       PubChem CID already resolved in step 2/3 (see
       search_kegg_compound_by_name) - this runs after PubChem resolution
       precisely because it needs a resolved CID to verify against.
    5. Determine mapping_status after all recovery attempts
    """
    global CURRENT_ROW_REQUEST_FAILED
    CURRENT_ROW_REQUEST_FAILED = False

    input_name = clean_wrapping_quotes(row.get("Name", ""))
    user_fields = get_user_provided_fields(row, input_name)
    cached_result = get_final_result_cache(input_name, row)
    if cached_result:
        RUNTIME_DIAGNOSTICS["compounds_resolved_from_final_cache"] += 1
        merged_result = merge_cached_result_with_user_fields(cached_result, user_fields, verbose, "final cache")
        if verbose:
            print(f"[final cache hit compatible] {input_name}")
            if user_fields.get("HMDB ID") and merged_result.get("PubChem ID"):
                print(f"[user HMDB preserved] {input_name}: {user_fields.get('HMDB ID')}")
        if merged_result.get("PubChem ID") and not merged_result.get("HMDB ID"):
            merged_result["HMDB ID"] = get_hmdb_ids_by_pubchem_cid(merged_result.get("PubChem ID", ""), verbose=verbose)
        merged_result["HMDB ID"] = normalize_hmdb_ids(merged_result.get("HMDB ID", ""))
        merged_result["Mapping_status"] = calculate_mapping_status(merged_result)
        set_final_result_cache(input_name, merged_result, row, user_fields=user_fields)
        return merged_result

    legacy_cached_result = {}
    if make_final_result_cache_key(input_name, row) != make_final_result_cache_key(input_name):
        legacy_cached_result = get_final_result_cache(input_name)
        if legacy_cached_result and verbose:
            print(f"[legacy cache merged] {input_name}")
        if legacy_cached_result:
            legacy_cached_result = merge_cached_result_with_user_fields(
                legacy_cached_result,
                user_fields,
                verbose,
                "legacy final cache",
            )

    online_requests_before = RUNTIME_DIAGNOSTICS["online_requests"]
    synonyms = get_query_names(input_name)

    # Read existing input identifiers - preserve these
    kegg_id = user_fields["KEGG ID"] or legacy_cached_result.get("KEGG ID", "")
    pubchem_id = user_fields["PubChem ID"] or legacy_cached_result.get("PubChem ID", "")
    smiles = user_fields["SMILES"] or legacy_cached_result.get("SMILES", "")
    inchikey = user_fields["InChIKey"] or legacy_cached_result.get("InChIKey", "")
    hmdb_id = user_fields["HMDB ID"] or legacy_cached_result.get("HMDB ID", "")
    formula = user_fields["Formula"] or legacy_cached_result.get("Formula", "")

    # Step 0 (KEGG): if a KEGG ID was provided, validate it (backfills
    # Formula if still missing) and, if PubChem is still empty, try to
    # recover a PubChem CID from the KEGG entry's own name - the same
    # KEGG-name-seed-to-PubChem enrichment the pre-removal tool did.
    if kegg_id:
        kegg_info = validate_kegg_id(kegg_id, verbose=verbose)
        formula = keep_user_value_or_fill(formula, kegg_info.get("formula", ""), "Formula", "KEGG ID", input_name, verbose)
        if not pubchem_id and kegg_info.get("name"):
            kegg_name_seed = kegg_info["name"].split(";")[0].strip()
            if kegg_name_seed:
                pubchem_id = get_pubchem_cid_by_name(kegg_name_seed, verbose=verbose)
                if verbose and pubchem_id:
                    print(f"PubChem match source: KEGG ID name ({pubchem_id})")

    # Step 1: If PubChem ID is provided, fetch missing properties
    if pubchem_id:
        if verbose:
            print(f"PubChem match source: provided CID ({pubchem_id})")
        pubchem_data = get_pubchem_data_by_cid(pubchem_id, verbose=verbose)
        smiles = keep_user_value_or_fill(smiles, pubchem_data.get("SMILES", ""), "SMILES", "PubChem CID", input_name, verbose)
        inchikey = keep_user_value_or_fill(inchikey, pubchem_data.get("InChIKey", ""), "InChIKey", "PubChem CID", input_name, verbose)
        formula = keep_user_value_or_fill(formula, pubchem_data.get("Formula", ""), "Formula", "PubChem CID", input_name, verbose)

    # Step 2: Search for PubChem ID by provided structure identifiers if not found
    if not pubchem_id and inchikey:
        pubchem_id = get_pubchem_cid_by_inchikey(inchikey)
        if verbose and pubchem_id:
            print(f"PubChem match source: InChIKey ({pubchem_id})")

    if not pubchem_id and smiles:
        pubchem_id = get_pubchem_cid_by_smiles(smiles)
        if verbose and pubchem_id:
            print(f"PubChem match source: SMILES ({pubchem_id})")

    # Step 3: If PubChem ID was recovered, fetch missing properties
    if pubchem_id and (not smiles or not inchikey or not formula):
        pubchem_data = get_pubchem_data_by_cid(pubchem_id, verbose=verbose)
        smiles = keep_user_value_or_fill(smiles, pubchem_data.get("SMILES", ""), "SMILES", "PubChem CID", input_name, verbose)
        inchikey = keep_user_value_or_fill(inchikey, pubchem_data.get("InChIKey", ""), "InChIKey", "PubChem CID", input_name, verbose)
        formula = keep_user_value_or_fill(formula, pubchem_data.get("Formula", ""), "Formula", "PubChem CID", input_name, verbose)

    # Step 4: Search for PubChem ID by synonym if not found
    if not pubchem_id and synonyms:
        for synonym in synonyms:
            cid = get_pubchem_cid_by_name(synonym, verbose=verbose)
            if cid:
                pubchem_id = cid
                break

    # Step 5: If a PubChem ID is available, fetch any missing PubChem properties
    if pubchem_id and (not smiles or not inchikey or not formula):
        pubchem_data = get_pubchem_data_by_cid(pubchem_id, verbose=verbose)
        smiles = keep_user_value_or_fill(smiles, pubchem_data.get("SMILES", ""), "SMILES", "PubChem CID", input_name, verbose)
        inchikey = keep_user_value_or_fill(inchikey, pubchem_data.get("InChIKey", ""), "InChIKey", "PubChem CID", input_name, verbose)
        formula = keep_user_value_or_fill(formula, pubchem_data.get("Formula", ""), "Formula", "PubChem CID", input_name, verbose)

    # Step 6: If PubChem ID is available and HMDB ID is empty, try PubChem synonym enrichment
    if pubchem_id and hmdb_id:
        if user_fields.get("HMDB ID") and verbose:
            print(f"[user HMDB preserved] {input_name}: {hmdb_id}")
    elif pubchem_id and not hmdb_id:
        hmdb_id = get_hmdb_ids_by_pubchem_cid(pubchem_id, verbose=verbose)
    hmdb_id = normalize_hmdb_ids(hmdb_id)

    # Step 7 (KEGG): search by name/synonym only if no KEGG ID exists yet.
    # Runs here, AFTER PubChem resolution, specifically so pubchem_id is
    # available to verify each name-matched KEGG candidate against (see
    # search_kegg_compound_by_name's docstring for why that matters).
    if not kegg_id and synonyms:
        for synonym in synonyms:
            kegg_match = search_kegg_compound_by_name(synonym, resolved_pubchem_id=pubchem_id, verbose=verbose)
            if kegg_match.get("KEGG ID"):
                kegg_id = kegg_match["KEGG ID"]
                if not formula:
                    kegg_info = validate_kegg_id(kegg_id, verbose=verbose)
                    formula = keep_user_value_or_fill(formula, kegg_info.get("formula", ""), "Formula", "KEGG ID", input_name, verbose)
                break

    resolved = {
        "Name": input_name,
        "KEGG ID": kegg_id,
        "PubChem ID": pubchem_id,
        "SMILES": smiles,
        "InChIKey": inchikey,
        "HMDB ID": hmdb_id,
        "Formula": formula,
    }
    resolved["Mapping_status"] = calculate_mapping_status(resolved)
    set_final_result_cache(input_name, resolved, row, user_fields=user_fields)
    if RUNTIME_DIAGNOSTICS["online_requests"] > online_requests_before:
        RUNTIME_DIAGNOSTICS["compounds_resolved_online"] += 1
    return resolved


def resolve_compound_table(df: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
    """Resolve an entire compound input table and return the result dataframe."""
    ensure_cache_loaded()
    reset_runtime_diagnostics()
    start_time = time.time()
    df = df.copy()
    if "Name" in df.columns:
        df["Name"] = df["Name"].apply(clean_wrapping_quotes)
    output_rows = []
    for row_number, (_, row) in enumerate(df.iterrows(), start=1):
        input_name = cell_to_text(row.get("Name", ""), strip=False)
        if verbose:
            print(f"Processing {row_number}/{len(df)}: {input_name[:60]}")

        try:
            resolved = resolve_compound(row, verbose=verbose)
        except Exception as exc:
            fallback_kegg = clean_identifier(row.get("KEGG ID", ""))
            fallback_pubchem = clean_identifier(row.get("PubChem ID", ""))
            fallback_smiles = clean_smiles_text(row.get("SMILES", ""))
            fallback_inchikey = clean_identifier(row.get("InChIKey", ""))
            fallback_hmdb = normalize_hmdb_ids(row.get("HMDB ID", ""))
            fallback_formula = clean_identifier(row.get("Formula", ""))
            # fallback_status is recomputed unconditionally by
            # calculate_mapping_status() below regardless of this value, so
            # this is just for readability of what a raw-passthrough row
            # looks like before that normalization runs.
            if fallback_kegg and fallback_pubchem and fallback_smiles and fallback_inchikey:
                fallback_status = "mapped"
            elif any([fallback_kegg, fallback_pubchem, fallback_smiles, fallback_inchikey, fallback_hmdb, fallback_formula]):
                fallback_status = "partial"
            else:
                fallback_status = "not_found"
            resolved = {
                "Name": input_name,
                "KEGG ID": fallback_kegg,
                "PubChem ID": fallback_pubchem,
                "SMILES": fallback_smiles,
                "InChIKey": fallback_inchikey,
                "HMDB ID": fallback_hmdb,
                "Formula": fallback_formula,
                "Mapping_status": fallback_status,
            }

        # Keep only the required output columns
        resolved["HMDB ID"] = normalize_hmdb_ids(resolved.get("HMDB ID", ""))
        resolved["Mapping_status"] = calculate_mapping_status(resolved)
        filtered_result = {col: resolved.get(col, "") for col in OUTPUT_COLUMNS}
        output_rows.append(filtered_result)

        # (no partial files saved)

    output_df = pd.DataFrame(output_rows, columns=OUTPUT_COLUMNS)
    save_pubchem_cache()
    save_kegg_cache()
    save_pubchem_synonyms_cache()
    save_hmdb_cache()
    save_final_result_cache()
    if verbose:
        print_runtime_diagnostics(time.time() - start_time)
    return output_df


def find_input_file() -> str:
    """Choose the input CSV file path for the resolver."""
    # Look for test_compounds.csv first, fall back to Data/example_input.csv
    if os.path.exists("test_compounds.csv"):
        return "test_compounds.csv"

    example_path = os.path.join("Data", "example_input_2.csv")
    if os.path.exists(example_path):
        return example_path

    raise FileNotFoundError(
        "No input file found. Expected 'test_compounds.csv' in current directory "
        "or 'Data/example_input.csv'."
    )


def test_common_compounds_direct(verbose: bool = True) -> pd.DataFrame:
    """Direct test for common compounds expected to resolve through PubChem."""
    test_names = [
        "Gluconic acid",
        "Uridine",
        "D-glucose 6-phosphate",
        "Apigenin",
        "Linoleic acid",
        "Bilineurine;Choline",
    ]
    input_df = pd.DataFrame({"Name": test_names})
    result_df = resolve_compound_table(input_df, verbose=verbose)
    if verbose:
        print(result_df)
        print(result_df["Mapping_status"].value_counts())
    return result_df


def test_user_provided_identifier_preservation(verbose: bool = True) -> pd.DataFrame:
    """Direct test that user-provided identifiers are preserved during resolution."""
    original_functions = {
        "get_final_result_cache": get_final_result_cache,
        "set_final_result_cache": set_final_result_cache,
        "get_pubchem_data_by_cid": get_pubchem_data_by_cid,
        "get_pubchem_cid_by_name": get_pubchem_cid_by_name,
        "get_pubchem_cid_by_inchikey": get_pubchem_cid_by_inchikey,
        "get_pubchem_cid_by_smiles": get_pubchem_cid_by_smiles,
    }

    def fake_pubchem_data_by_cid(cid, verbose=False):
        return {
            "PubChem ID": str(cid),
            "SMILES": "ONLINE_SMILES",
            "InChIKey": "ONLINE-INCHIKEY",
            "Formula": "C1H1",
        }

    try:
        globals()["get_final_result_cache"] = lambda name, row=None: {}
        globals()["set_final_result_cache"] = lambda name, data, row=None: None
        globals()["get_pubchem_data_by_cid"] = fake_pubchem_data_by_cid
        globals()["get_pubchem_cid_by_name"] = lambda name, verbose=False: ""
        globals()["get_pubchem_cid_by_inchikey"] = lambda inchikey: ""
        globals()["get_pubchem_cid_by_smiles"] = lambda smiles: ""

        input_df = pd.DataFrame(
            [
                {"Name": "No Real ID", "PubChem ID": ""},
                {"Name": "Also Not A Real Name", "PubChem ID": "305"},
                {
                    "Name": "Preserve Provided Structures",
                    "PubChem ID": "305",
                    "SMILES": "USER_SMILES",
                    "InChIKey": "USER-INCHIKEY",
                    "HMDB ID": "HMDB0000001",
                    "Formula": "USERFORMULA",
                },
            ]
        )
        result_df = resolve_compound_table(input_df, verbose=verbose)

        assert result_df.loc[0, "Mapping_status"] in ("partial", "not_found")
        assert result_df.loc[1, "PubChem ID"] == "305"
        assert result_df.loc[2, "PubChem ID"] == "305"
        assert result_df.loc[2, "SMILES"] == "USER_SMILES"
        assert result_df.loc[2, "InChIKey"] == "USER-INCHIKEY"
        assert result_df.loc[2, "HMDB ID"] == "HMDB0000001"
        assert result_df.loc[2, "Formula"] == "USERFORMULA"
        assert list(result_df.columns) == OUTPUT_COLUMNS

        if verbose:
            print("User-provided identifier preservation test passed.")
            print(result_df)
        return result_df
    finally:
        for function_name, function_value in original_functions.items():
            globals()[function_name] = function_value


def test_final_cache_preserves_user_identifiers(verbose: bool = True) -> dict:
    """Direct tests for final-result cache safety with user-provided identifiers."""
    global CACHE_LOADED
    original_functions = {
        "get_pubchem_data_by_cid": get_pubchem_data_by_cid,
        "get_pubchem_cid_by_name": get_pubchem_cid_by_name,
        "get_pubchem_cid_by_inchikey": get_pubchem_cid_by_inchikey,
        "get_pubchem_cid_by_smiles": get_pubchem_cid_by_smiles,
        "save_pubchem_cache": save_pubchem_cache,
        "save_pubchem_synonyms_cache": save_pubchem_synonyms_cache,
        "save_final_result_cache": save_final_result_cache,
    }
    original_final_cache = FINAL_RESULT_CACHE.copy()
    original_cache_loaded = CACHE_LOADED

    def fake_pubchem_data_by_cid(cid, verbose=False):
        return {
            "PubChem ID": str(cid),
            "SMILES": f"SMILES_{cid}",
            "InChIKey": f"INCHIKEY-{cid}",
            "Formula": f"FORMULA{cid}",
        }

    def fake_pubchem_cid_by_name(name, verbose=False):
        return "111" if clean_wrapping_quotes(name) == "Cache Compound" else ""

    try:
        FINAL_RESULT_CACHE.clear()
        CACHE_LOADED = True
        globals()["get_pubchem_data_by_cid"] = fake_pubchem_data_by_cid
        globals()["get_pubchem_cid_by_name"] = fake_pubchem_cid_by_name
        globals()["get_pubchem_cid_by_inchikey"] = lambda inchikey: ""
        globals()["get_pubchem_cid_by_smiles"] = lambda smiles: ""
        globals()["save_pubchem_cache"] = lambda: None
        globals()["save_pubchem_synonyms_cache"] = lambda: None
        globals()["save_final_result_cache"] = lambda: None

        name_only = resolve_compound_table(pd.DataFrame([{"Name": "Cache Compound"}]), verbose=verbose)
        with_pubchem = resolve_compound_table(
            pd.DataFrame([{"Name": "Cache Compound", "PubChem ID": "305"}]),
            verbose=verbose,
        )
        fake_pubchem = resolve_compound_table(
            pd.DataFrame([{"Name": "Fake Unmapped PubChem", "PubChem ID": "777"}]),
            verbose=verbose,
        )

        exact_empty_cache_row = pd.Series({"Name": "Exact Empty Cache Compound", "PubChem ID": ""})
        FINAL_RESULT_CACHE[make_final_result_cache_key("Exact Empty Cache Compound", exact_empty_cache_row)] = {
            "Name": "Exact Empty Cache Compound",
            "PubChem ID": "",
            "SMILES": "",
            "InChIKey": "",
            "HMDB ID": "",
            "Formula": "",
            "Mapping_status": "not_found",
        }
        exact_empty_cache = resolve_compound_table(pd.DataFrame([exact_empty_cache_row.to_dict()]), verbose=verbose)

        assert name_only.loc[0, "PubChem ID"] == "111"
        assert with_pubchem.loc[0, "PubChem ID"] == "305"
        assert with_pubchem.loc[0, "Mapping_status"] in ("partial", "mapped")
        assert fake_pubchem.loc[0, "PubChem ID"] == "777"
        assert fake_pubchem.loc[0, "Mapping_status"] in ("partial", "mapped")
        assert exact_empty_cache.loc[0, "Mapping_status"] == "not_found"

        results = {
            "name_only": name_only,
            "with_pubchem": with_pubchem,
            "fake_pubchem": fake_pubchem,
            "exact_empty_cache": exact_empty_cache,
        }
        if verbose:
            print("Final-result cache user-identifier preservation tests passed.")
        return results
    finally:
        FINAL_RESULT_CACHE.clear()
        FINAL_RESULT_CACHE.update(original_final_cache)
        CACHE_LOADED = original_cache_loaded
        for function_name, function_value in original_functions.items():
            globals()[function_name] = function_value


def test_hmdb_pubchem_enrichment(verbose: bool = True) -> pd.DataFrame:
    """Direct test that HMDB IDs can be enriched from PubChem synonyms."""
    original_functions = {
        "get_final_result_cache": get_final_result_cache,
        "set_final_result_cache": set_final_result_cache,
        "get_pubchem_data_by_cid": get_pubchem_data_by_cid,
        "fetch_pubchem_json_strings": fetch_pubchem_json_strings,
        "get_pubchem_cid_by_name": get_pubchem_cid_by_name,
        "get_pubchem_cid_by_inchikey": get_pubchem_cid_by_inchikey,
        "get_pubchem_cid_by_smiles": get_pubchem_cid_by_smiles,
        "save_hmdb_cache": save_hmdb_cache,
    }
    original_hmdb_cache = HMDB_CACHE.copy()

    def fake_pubchem_data_by_cid(cid, verbose=False):
        return {
            "PubChem ID": str(cid),
            "SMILES": f"SMILES_{cid}",
            "InChIKey": f"INCHIKEY-{cid}",
            "Formula": f"FORMULA{cid}",
        }

    def fake_fetch_pubchem_json_strings(url, database, query, verbose=False):
        if str(query) == "305" and database == "PubChem synonyms":
            return ["Choline", "HMDB0000097", "hmdb0000097"], True
        if str(query) == "5793" and database == "PubChem PUG-View":
            return ["Nested source mentions HMDB0000122"], True
        if str(query) == "145742" and database == "PubChem PUG-View":
            return ["Proline source HMDB0000162; duplicate HMDB00162"], True
        return ["No HMDB here"], True

    try:
        HMDB_CACHE.clear()
        globals()["get_final_result_cache"] = lambda name, row=None: {}
        globals()["set_final_result_cache"] = lambda name, data, row=None: None
        globals()["get_pubchem_data_by_cid"] = fake_pubchem_data_by_cid
        globals()["fetch_pubchem_json_strings"] = fake_fetch_pubchem_json_strings
        globals()["get_pubchem_cid_by_name"] = lambda name, verbose=False: ""
        globals()["get_pubchem_cid_by_inchikey"] = lambda inchikey: ""
        globals()["get_pubchem_cid_by_smiles"] = lambda smiles: ""
        globals()["save_hmdb_cache"] = lambda: None

        input_df = pd.DataFrame(
            [
                {"Name": "HMDB From PubChem", "PubChem ID": "305"},
                {"Name": "HMDB From PUGView", "PubChem ID": "5793"},
                {"Name": "User HMDB Preserved", "PubChem ID": "305", "HMDB ID": "HMDB0097"},
                {"Name": "No HMDB Available", "PubChem ID": "999"},
                {"Name": "Proline HMDB normalization", "PubChem ID": "145742"},
                {"Name": "User HMDB Normalized", "PubChem ID": "305", "HMDB ID": "hmdb162"},
            ]
        )
        result_df = resolve_compound_table(input_df, verbose=verbose)

        assert result_df.loc[0, "HMDB ID"] == "HMDB0000097"
        assert result_df.loc[1, "HMDB ID"] == "HMDB0000122"
        assert result_df.loc[2, "HMDB ID"] == "HMDB0000097"
        assert result_df.loc[3, "HMDB ID"] == ""
        assert result_df.loc[4, "HMDB ID"] == "HMDB0000162"
        assert result_df.loc[5, "HMDB ID"] == "HMDB0000162"
        assert list(result_df.columns) == OUTPUT_COLUMNS

        if verbose:
            print("HMDB PubChem enrichment test passed.")
            print(result_df)
        return result_df
    finally:
        HMDB_CACHE.clear()
        HMDB_CACHE.update(original_hmdb_cache)
        for function_name, function_value in original_functions.items():
            globals()[function_name] = function_value


def test_hmdb_normalization(verbose: bool = True) -> list:
    """Direct test for canonical HMDB ID normalization and deduplication."""
    input_values = ["HMDB0000162", "HMDB00162", "hmdb162"]
    normalized_values = [normalize_hmdb_ids(value) for value in input_values]
    combined = normalize_hmdb_ids(";".join(input_values))

    assert normalized_values == ["HMDB0000162", "HMDB0000162", "HMDB0000162"]
    assert combined == "HMDB0000162"

    if verbose:
        print(f"HMDB normalization test passed: {input_values} -> {combined}")
    return normalized_values


def test_hmdb_pubchem_live_cids(verbose: bool = True) -> pd.DataFrame:
    """Manual live check for HMDB recovery from PubChem CIDs 305 and 5793."""
    input_df = pd.DataFrame(
        [
            {"Name": "Choline live HMDB check", "PubChem ID": "305"},
            {"Name": "Glucose live HMDB check", "PubChem ID": "5793"},
        ]
    )
    result_df = resolve_compound_table(input_df, verbose=verbose)
    if verbose:
        print(result_df[["Name", "PubChem ID", "HMDB ID", "Mapping_status"]])
    return result_df


def run_resolver(verbose: bool = False):
    """Run the resolver script from the command line."""
    try:
        start_time = time.time()
        input_path = find_input_file()
        if verbose:
            print(f"Loading input from: {input_path}")
        df = load_input_table(input_path, verbose=verbose)
        df = remove_duplicate_input_names(df)
        if verbose:
            print(f"Loaded {len(df)} compounds.\n")

        result_df = resolve_compound_table(df, verbose=verbose)

        # Save final results to XLSX only
        xlsx_path = "compound_resolver_test_results.xlsx"
        try:
            import openpyxl
            result_df.to_excel(xlsx_path, index=False, engine="openpyxl")
            print(f"Saved XLSX results to: {xlsx_path}")
        except ImportError:
            print("(openpyxl not available, cannot save XLSX)")
        except PermissionError:
            print(
                "Could not save compound_resolver_test_results.xlsx because the file is locked. "
                "Close it if it is open and run the resolver again."
            )

        total_runtime = time.time() - start_time
        compound_count = len(result_df)
        average_seconds = total_runtime / compound_count if compound_count else 0
        status_counts = result_df["Mapping_status"].value_counts()

        print(f"\nResolution complete. {compound_count} compounds processed.")
        print(f"Total runtime: {total_runtime:.2f} seconds")
        print(f"Average seconds per compound: {average_seconds:.2f}")
        print(f"mapped: {int(status_counts.get('mapped', 0))}")
        print(f"partial: {int(status_counts.get('partial', 0))}")
        print(f"not_found: {int(status_counts.get('not_found', 0))}")
        print(f"Output file path: {xlsx_path}")

        return result_df

    except Exception as exc:
        print(f"Error: {exc}")
        raise


if __name__ == "__main__":
    if os.environ.get("RUN_HMDB_PUBCHEM_TEST") == "1":
        test_hmdb_pubchem_live_cids(verbose=True)
    else:
        run_resolver()
