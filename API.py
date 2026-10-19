# -*- coding: utf-8 -*-
"""
Omiscera FastAPI backend.

Serves the HTML/JS frontend (index.html + static/) and exposes the same
workflows as the Streamlit prototype in main3.py: Metabolic Connection and
Annotation Analytics (Genomics ORA). main3.py itself is not used by this
file - the small set of pure-pandas orchestration helpers it defines have
been ported below so this module has no Streamlit dependency.

Run with:
    uvicorn API:app --reload --port 8000
"""

from __future__ import annotations

import base64
import io
import json
import math
import os
import re
import threading
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Optional

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend; must be set before pyplot import
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False

try:
    from PIL import Image as _PILImage
    _PIL_AVAILABLE = True
except ImportError:
    _PIL_AVAILABLE = False

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, Form, Header, HTTPException, UploadFile
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles

import database_catalog as dc
import test_compound_resolver as resolver_module
import Genomics_module_retriever as gene_resolver
import Genomics_module as genomics_ora
import Genomics_network as genomics_net
import Transcriptomics_network as transcriptomics_network_engine
import Genomics_annotation as genomics_annot
import Transcriptomics_module_retriever as transcriptomics_retriever
import Transcriptomics_module as transcriptomics_engine
import Transcriptomics_regulon as transcriptomics_regulon_engine
import Proteomics_module_retriever as proteomics_retriever
import Proteomics_module as proteomics_engine
import Proteomics_network as proteomics_network_engine
import Proteomics_ptm_retriever as proteomics_ptm_retriever
import Proteomics_ptm as proteomics_ptm_engine
# Metabolomics enrichment. Metabolomics_retriever is a DELIBERATE fork of
# test_compound_resolver (imported above as resolver_module, which backs the
# separate Metabolic Connection feature at /api/mc/*): the fork exists so this
# module can adapt the resolver freely without any risk to that live feature,
# so these two must never be consolidated back into one import.
import Metabolomics_retriever as metabolomics_retriever
import Metabolomics_module as metabolomics_engine
import Multiomics_module as multiomics_engine
import Multiomics_module_retriever as multiomics_retriever
import Multiomics_mapping as multiomics_mapping
import Multiomics_correlation as multiomics_correlation
import Multiomics_network as multiomics_network
import Multiomics_mofa as multiomics_mofa
import Multiomics_diablo as multiomics_diablo

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Constants ported from main3.py
# ---------------------------------------------------------------------------

DISPLAY_COLUMNS = [
    "Name",
    "KEGG ID",
    "PubChem ID",
    "SMILES",
    "InChIKey",
    "HMDB ID",
    "Formula",
    "Mapping_status",
]
EXPECTED_INPUT_COLUMNS = ["Name", "KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"]
USEFUL_DUPLICATE_FIELDS = ["KEGG ID", "PubChem ID", "SMILES", "InChIKey", "HMDB ID", "Formula"]


# ---------------------------------------------------------------------------
# Pure-pandas helpers ported from main3.py (unchanged logic, no Streamlit)
# ---------------------------------------------------------------------------

def normalize_column_name(column_name: str) -> str:
    """Clean a column name by stripping whitespace and removing BOM characters."""
    if column_name is None:
        return ""
    name_text = str(column_name)
    name_text = name_text.strip()
    name_text = name_text.lstrip("﻿").lstrip("￾")
    name_text = name_text.replace("Ã¯Â»Â¿", "")
    name_text = name_text.replace("﻿", "").replace("￾", "")
    return name_text


def standardize_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename dataframe columns to the standard field names when possible."""
    rename_map = {}
    hmdb_source = ""
    hmdb_priority = ["hmdbid", "hmdb", "hdmbid", "hdmb"]
    normalized_columns = {}
    for col in df.columns:
        normalized = (
            normalize_column_name(col)
            .replace("\xa0", " ")
            .replace(" ", " ")
            .replace(" ", " ")
            .strip()
        )
        key = normalized.lower().replace(" ", "").replace("_", "").replace("-", "")
        normalized_columns[col] = key

        if key == "name":
            rename_map[col] = "Name"
        elif key in ("keggid", "kegg"):
            rename_map[col] = "KEGG ID"
        elif key in ("pubchemid", "pubchem"):
            rename_map[col] = "PubChem ID"
        elif key == "smiles":
            rename_map[col] = "SMILES"
        elif key in ("inchikey",):
            rename_map[col] = "InChIKey"
        elif key == "formula":
            rename_map[col] = "Formula"
        else:
            rename_map[col] = normalized

    for hmdb_name in hmdb_priority:
        for col, key in normalized_columns.items():
            if key == hmdb_name:
                hmdb_source = col
                break
        if hmdb_source:
            break

    df = df.rename(columns=rename_map)
    if hmdb_source:
        df["HMDB ID"] = df[hmdb_source]
        final_source_name = rename_map.get(hmdb_source, hmdb_source)
        if final_source_name != "HMDB ID" and final_source_name in df.columns:
            df = df.drop(columns=[final_source_name])
    return df


def clean_wrapping_quotes_for_input(value):
    """Strip spaces and accidental full-value wrapping quotes from uploaded text cells."""
    if pd.isna(value):
        return ""
    text = str(value).replace("\xa0", " ").strip()
    while len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        text = text[1:-1].strip()
    return text


def clean_uploaded_text_cells(df: pd.DataFrame) -> pd.DataFrame:
    """Clean uploaded string cells without changing chemical characters inside values."""
    return df.apply(lambda column: column.map(clean_wrapping_quotes_for_input))


def is_non_empty_input_value(value) -> bool:
    """Return True for uploaded identifier values that contain useful text."""
    if pd.isna(value):
        return False
    text = str(value).strip()
    return bool(text) and text.lower() not in ("nan", "none")


def normalize_name_for_duplicate_detection(value) -> str:
    """Normalize Name values only for duplicate comparison."""
    text = clean_wrapping_quotes_for_input(value)
    text = text.replace("\xa0", " ").replace(" ", " ").replace(" ", " ")
    text = " ".join(text.split())
    return text.lower()


def uploaded_row_value(row: pd.Series, field: str):
    """Return the first scalar value for a possibly duplicated uploaded column."""
    value = row.get(field, "")
    if isinstance(value, pd.Series):
        if value.empty:
            return ""
        return value.iloc[0]
    return value


def count_useful_identifier_fields(row: pd.Series) -> int:
    """Count non-empty useful identifier fields in an uploaded row."""
    count = 0
    for field in USEFUL_DUPLICATE_FIELDS:
        if field in row.index and is_non_empty_input_value(uploaded_row_value(row, field)):
            count += 1
    return count


def first_named_column(df: pd.DataFrame, column_name: str) -> pd.Series:
    """Return the first column with this name, even if duplicate names exist."""
    selected = df.loc[:, column_name]
    if isinstance(selected, pd.DataFrame):
        return pd.Series(selected.to_numpy()[:, 0], index=selected.index, name=column_name)
    if isinstance(selected, pd.Series):
        return selected
    return pd.Series([selected], name=column_name)


def duplicated_column_names(columns) -> list:
    """Return duplicated dataframe column names in display order."""
    seen = set()
    duplicates = []
    for column in columns:
        if column in seen and column not in duplicates:
            duplicates.append(column)
        seen.add(column)
    return duplicates


def score_uploaded_dataframe(df: pd.DataFrame, delimiter_order: int) -> tuple:
    """Score parsed CSV candidates while requiring a Name column."""
    if "Name" not in df.columns:
        return (-1, -delimiter_order)

    expected_count = sum(1 for column in EXPECTED_INPUT_COLUMNS if column in df.columns)
    name_column = first_named_column(df, "Name")
    empty_name_count = int((name_column.astype(str).str.strip() == "").sum())
    duplicate_count = len(duplicated_column_names(df.columns.tolist()))
    unnamed_count = sum(1 for column in df.columns if str(column).lower().startswith("unnamed"))
    non_empty_names = len(df) - empty_name_count
    score = (expected_count * 1000) + (non_empty_names * 10) - (empty_name_count * 10)
    score -= duplicate_count * 100
    score -= unnamed_count * 50
    return (score, -delimiter_order)


def read_uploaded_csv_bytes(raw_bytes: bytes) -> tuple:
    """Read uploaded CSV bytes with normal quote handling and delimiter selection."""
    try:
        text = raw_bytes.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw_bytes.decode("utf-8", errors="replace")

    delimiter_options = [(",", "comma (,)"), ("\t", "tab"), (";", "semicolon (;)")]
    candidates = []
    errors = []

    for delimiter_order, (delimiter, delimiter_label) in enumerate(delimiter_options):
        try:
            parsed_df = pd.read_csv(
                io.StringIO(text),
                sep=delimiter,
                dtype=str,
                quotechar='"',
                doublequote=True,
                keep_default_na=False,
            )
            parsed_df = standardize_dataframe_columns(parsed_df)
            parsed_df = clean_uploaded_text_cells(parsed_df)
            candidates.append((score_uploaded_dataframe(parsed_df, delimiter_order), parsed_df, delimiter_label))
        except Exception as exc:
            errors.append(f"{delimiter_label}: {exc}")

    valid_candidates = [candidate for candidate in candidates if candidate[0][0] >= 0]
    if not valid_candidates:
        error_text = "; ".join(errors) if errors else "No delimiter produced a table with a Name column."
        raise ValueError(f"Unable to parse uploaded CSV with a required Name column. {error_text}")

    _, best_df, detected_delimiter = max(valid_candidates, key=lambda candidate: candidate[0])
    return best_df.copy(), detected_delimiter


def deduplicate_uploaded_compounds(df: pd.DataFrame) -> tuple:
    """Remove duplicated compound Names while preserving the most informative row."""
    rows_before = len(df)
    if "Name" not in df.columns:
        stats = {
            "duplicate_name_count": 0,
            "duplicate_rows_removed": 0,
            "rows_before": rows_before,
            "rows_after": rows_before,
        }
        return df.copy(), pd.DataFrame(), stats

    name_column = first_named_column(df, "Name")
    duplicate_groups = {}
    for position, value in enumerate(name_column.tolist()):
        duplicate_key = normalize_name_for_duplicate_detection(value)
        if duplicate_key:
            duplicate_groups.setdefault(duplicate_key, []).append(position)

    duplicate_groups = {key: positions for key, positions in duplicate_groups.items() if len(positions) > 1}
    if not duplicate_groups:
        stats = {
            "duplicate_name_count": 0,
            "duplicate_rows_removed": 0,
            "rows_before": rows_before,
            "rows_after": rows_before,
        }
        return df.copy(), pd.DataFrame(), stats

    rows_to_remove = set()
    kept_rows = {}
    report_rows = []

    for positions in duplicate_groups.values():
        scored_positions = [
            (count_useful_identifier_fields(df.iloc[position]), -position, position)
            for position in positions
        ]
        _, _, kept_position = max(scored_positions)
        kept_row = df.iloc[kept_position].copy()
        conflict_fields = []

        for position in positions:
            if position == kept_position:
                continue
            rows_to_remove.add(position)
            duplicate_row = df.iloc[position]
            for field in USEFUL_DUPLICATE_FIELDS:
                if field not in df.columns:
                    continue
                kept_value = uploaded_row_value(kept_row, field)
                duplicate_value = uploaded_row_value(duplicate_row, field)
                if not is_non_empty_input_value(kept_value) and is_non_empty_input_value(duplicate_value):
                    kept_row[field] = duplicate_value
                elif (
                    is_non_empty_input_value(kept_value)
                    and is_non_empty_input_value(duplicate_value)
                    and str(kept_value).strip() != str(duplicate_value).strip()
                    and field not in conflict_fields
                ):
                    conflict_fields.append(field)

        kept_rows[kept_position] = kept_row
        display_name = clean_wrapping_quotes_for_input(name_column.iloc[kept_position])
        removed_positions = [position for position in positions if position != kept_position]
        report_rows.append(
            {
                "duplicated Name": display_name,
                "number of occurrences": len(positions),
                "kept row position": kept_position,
                "removed row positions": "; ".join(str(position) for position in removed_positions),
                "conflict found": "Yes" if conflict_fields else "No",
                "conflicting fields": "; ".join(conflict_fields),
            }
        )

    output_rows = []
    for position in range(rows_before):
        if position in rows_to_remove:
            continue
        if position in kept_rows:
            output_rows.append(kept_rows[position])
        else:
            output_rows.append(df.iloc[position])

    deduplicated_df = pd.DataFrame(output_rows, columns=df.columns).reset_index(drop=True)
    report_df = pd.DataFrame(report_rows)
    stats = {
        "duplicate_name_count": len(duplicate_groups),
        "duplicate_rows_removed": len(rows_to_remove),
        "rows_before": rows_before,
        "rows_after": len(deduplicated_df),
    }
    return deduplicated_df, report_df, stats


def order_mapped_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Show expected mapped columns first while preserving any additional resolver columns."""
    ordered_columns = [column for column in DISPLAY_COLUMNS if column in df.columns]
    extra_columns = [column for column in df.columns if column not in ordered_columns]
    return df.reindex(columns=ordered_columns + extra_columns)


def run_connection_analysis(
    mapped_df: pd.DataFrame,
    rank_settings: dict,
    hmdb_import_by_rank: dict,
    allow_online: bool = True,
) -> dict:
    """Run compound presence lookups for mapped and partial compounds."""
    results = {}
    filtered_df = mapped_df[mapped_df["Mapping_status"].isin(["mapped", "partial"])].reset_index(drop=True)

    for rank_name, settings in rank_settings.items():
        rank_letter = str(rank_name).replace("Rank", "").strip()
        hmdb_table = hmdb_import_by_rank.get(rank_letter)
        _enzyme_dbs = {"Rhea + UniProt", "KEGG"}
        table1_databases = [
            db for db in settings if dc.canonical_database_name(db) not in _enzyme_dbs
        ]

        rhea_taxon_id = ""
        kegg_org = ""
        for db, db_settings in settings.items():
            canonical_db = dc.canonical_database_name(db)
            if canonical_db == "Rhea + UniProt":
                rhea_taxon_id = str(db_settings.get("context_id", "")).strip()
            elif canonical_db == "KEGG":
                kegg_org = str(db_settings.get("context_id", "")).strip()
        rhea_taxon_reactions = (
            dc.get_uniprot_taxon_rhea_reactions(rhea_taxon_id, allow_online)
            if rhea_taxon_id else None
        )
        kegg_org_enzymes = (
            dc.get_kegg_organism_enzymes(kegg_org, allow_online)
            if kegg_org else None
        )

        table1_rows = []
        rhea_enzyme_rows = []
        kegg_enzyme_rows = []

        for _, compound in filtered_df.iterrows():
            compound_name = str(compound.get("Name", "")).strip()
            row = {"Name": compound_name}

            for db, db_settings in settings.items():
                if dc.canonical_database_name(db) in _enzyme_dbs:
                    continue
                context_selection = {
                    "database": db,
                    "context_id": db_settings.get("context_id", ""),
                    "context_name": db_settings.get("context_name", ""),
                    "context_type": db_settings.get("context_type", ""),
                }
                imported_table = hmdb_table if dc.canonical_database_name(db) == "HMDB" else None
                result = dc.search_compound_in_database(
                    db,
                    compound,
                    context_selection,
                    allow_online=allow_online,
                    imported_table=imported_table,
                )
                raw_status = result.get("status", "not_evaluated")
                cell = raw_status if raw_status in {"found", "not_found", "unknown_id"} else "not_evaluated"
                row[db] = cell.capitalize()

            table1_rows.append(row)

            if rhea_taxon_id:
                rhea_rows = dc.build_rhea_uniprot_enzyme_rows(
                    compound, rhea_taxon_id, allow_online,
                    taxon_reactions=rhea_taxon_reactions,
                )
                rhea_enzyme_rows.extend(rhea_rows)

            if kegg_org:
                kegg_rows = dc.build_kegg_enzyme_rows(
                    compound, kegg_org, allow_online,
                    org_enzymes=kegg_org_enzymes,
                )
                kegg_enzyme_rows.extend(kegg_rows)

        if table1_databases:
            table1 = pd.DataFrame(table1_rows) if table1_rows else pd.DataFrame(columns=["Name"] + table1_databases)
        else:
            table1 = pd.DataFrame()
        table2 = pd.DataFrame(
            rhea_enzyme_rows,
            columns=["Name", "Enzyme", "Role", "Enzyme presence"],
        )
        if not table2.empty:
            _presence_order = {"Found": 0, "Not_evaluated": 1, "Not_found": 2}
            table2 = (
                table2
                .assign(_rank=table2["Enzyme presence"].map(_presence_order).fillna(99))
                .sort_values(by="_rank", kind="stable")
                .drop(columns="_rank")
                .reset_index(drop=True)
            )

        table3 = pd.DataFrame(
            kegg_enzyme_rows,
            columns=["Name", "Enzyme", "Role", "Enzyme presence", "Pathway"],
        )
        if not table3.empty:
            _presence_order = {"Found": 0, "Not_evaluated": 1, "Not_found": 2}
            table3 = (
                table3
                .assign(_rank=table3["Enzyme presence"].map(_presence_order).fillna(99))
                .sort_values(by="_rank", kind="stable")
                .drop(columns="_rank")
                .reset_index(drop=True)
            )

        results[rank_name] = {"table1": table1, "table2": table2, "table3": table3, "selections": settings}

    return results


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------

def df_payload(df: Optional[pd.DataFrame]) -> dict:
    """Serialize a DataFrame to a JSON-safe {columns, rows} payload."""
    if df is None or df.empty:
        columns = list(df.columns) if df is not None else []
        return {"columns": columns, "rows": []}
    # NB: cast to object BEFORE substituting None. `df.where(cond, None)` on a
    # float64 column silently re-coerces None back to NaN (a float column can't
    # hold None), so the NaN survives into to_dict() and breaks Starlette's
    # strict JSON encoder (allow_nan=False -> "Out of range float ... nan").
    safe_df = df.astype(object).where(pd.notna(df), None)
    return {"columns": list(safe_df.columns), "rows": safe_df.to_dict(orient="records")}


def df_to_csv_response(df: pd.DataFrame, filename: str) -> Response:
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    return Response(
        content=csv_bytes,
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# Session store (mirrors st.session_state, keyed per browser session)
# ---------------------------------------------------------------------------

_SESSIONS_LOCK = threading.Lock()
SESSIONS: dict[str, dict[str, Any]] = {}


def _fresh_ora_state() -> dict:
    """Single source of truth for the 'ora' session sub-dict shape.

    Used by new_session_state() and by both gene/sequence upload endpoints
    (which fully reset this sub-dict on a new upload) so the schema can
    never drift out of sync between the three call sites.
    """
    return {
        "identifier_type": "",
        "gene_df": None,
        "resolved_df": None,
        "table1_df": None,
        "display_cols": [],
        "id_col": "",
        "final_df": None,
        # KEGG organism code (e.g. "hsa") the "KEGG ID" column of resolved_df/
        # table1_df/final_df was actually resolved against, set by
        # _job_resolve_kegg_ids. The single source of truth _job_run_ora/
        # _job_run_burden trust for KEGG's universe organism - never the
        # config screen's taxon_id, which is only synced FROM this value.
        "kegg_organism_code": "",
        "ora_config": None,
        "ora_results": None,
        "custom_db": None,   # parsed user-uploaded gene-set database (Part 10)
        "background_by_db": {},  # {database: [accession, ...]} user-supplied background lists
        "variant_scores": {},    # {gene_id (lowercased): float} Gene-set Burden Test score capture
        "score_direction": "higher",  # "higher" or "lower" = more burden
        "burden_config": None,
        "burden_results": None,
        # Network Propagation / Module Detection (submodule 3, STRING-based).
        # Reuses the shared retrieval (final_df) + the same variant_scores/
        # score_direction capture as Burden (optional seed weighting here).
        "netprop_taxon_id": "",
        "netprop_taxon_name": "",
        "netprop_config": None,
        "netprop_results": None,
        # Gene Function Annotation (submodule 4) — per-gene descriptive GO +
        # UniProt lookup. Reuses the shared retrieval (final_df); single-gene
        # mode also populates final_df (with one row) so the SAME configure/
        # run/results flow below works unchanged for both modes.
        "annotation_config": None,
        "annotation_results": None,
    }


def _fresh_transcriptomics_state() -> dict:
    """Single source of truth for the 'transcriptomics' session sub-dict shape."""
    return {
        "identifier_type": "",
        "de_df": None,
        "column_mapping": None,
        "taxon_id": "",
        "taxon_name": "",
        "resolved_df": None,
        "meta": None,
        # Over-Representation Analysis (first analysis tile). Shares retrieval
        # (resolved_df/taxon_id) above; only this stage is ORA-specific.
        "ora_config": None,
        "ora_results": None,
        # Parsed custom gene-set database (load_custom_geneset_db's return
        # value), when the "Custom database" checkbox is used.
        "custom_db": None,
        # Optional user-uploaded background/reference gene-symbol set
        # (load_user_background_list's return value).
        "background_symbols": None,
        # GSEA (second analysis tile). Shares retrieval AND the custom_db/
        # background_symbols above with ORA - only this stage is GSEA-specific.
        "gsea_config": None,
        "gsea_results": None,
        # Transcription Factor / Regulon (third analysis tile). Shares
        # retrieval above; a separate config/results pair per analysis
        # (target enrichment vs. activity inference - see
        # Transcriptomics_regulon.py), since a run can select either or both.
        "regulon_config": None,
        "regulon_enrichment_results": None,
        "regulon_activity_results": None,
        "network_config": None,
        "network_results": None,
        # KEGG's organism is normally DERIVED from taxon_id/taxon_name above
        # (Transcriptomics_module.resolve_kegg_organism_code), never picked
        # separately. This holds a user's manual confirmation for the rare
        # case that automatic derivation can't taxon-verify a unique KEGG
        # organism (several name-matched KEGG genomes, none carrying this
        # run's NCBI taxon - e.g. a species with only strain-level KEGG
        # genomes). Shape: {"taxon_id", "taxon_name", "organism_code",
        # "organism_name"} when set, else None. The stored taxon_id/
        # taxon_name are what every reader re-checks against the session's
        # CURRENT values before trusting organism_code (see
        # _transcriptomics_kegg_confirmed_code) - this is what stops a
        # confirmation from surviving an organism change and silently
        # applying to a different species.
        "kegg_organism_confirmation": None,
    }


def _fresh_proteomics_state() -> dict:
    """Single source of truth for the 'proteomics' session sub-dict shape.
    Mirrors _fresh_transcriptomics_state() - see that function's comment for
    why a factory function (not an inline literal) is the pattern here."""
    return {
        "de_df": None,
        "column_mapping": None,
        "taxon_id": "",
        "taxon_name": "",
        "resolved_df": None,
        "meta": None,
        "ora_config": None,
        "ora_results": None,
        "quant_config": None,
        "quant_results": None,
        "network_config": None,
        "network_results": None,
        # PTM / Phosphosite Enrichment - a SEPARATE, site-level upload/resolve
        # path (Proteomics_ptm_retriever.py), not the protein-level
        # de_df/resolved_df pair above. Shares taxon_id/taxon_name (an
        # organism selection is still one per session) but nothing else -
        # RESOLVED_SITE_COLUMNS has no protein-level analogue to reuse.
        "ptm_sites_df": None,
        "ptm_column_mapping": None,
        "resolved_sites_df": None,
        "ptm_meta": None,
        "ptm_config": None,
        "ptm_activity_results": None,
        "ptm_ora_results": None,
        "custom_db": None,
        "background_accessions": None,
        # KEGG's organism is normally DERIVED from taxon_id/taxon_name above
        # (Transcriptomics_module.resolve_kegg_organism_code, reused
        # unchanged - see Proteomics_module's module docstring on the KEGG
        # database), never picked separately. This holds a manual
        # confirmation for the rare case automatic resolution can't
        # taxon-verify a unique KEGG organism. Same shape/contract as
        # Transcriptomics' identical field - see _kegg_confirmed_code.
        "kegg_organism_confirmation": None,
    }


def _fresh_metabolomics_state() -> dict:
    """Single source of truth for the 'metabolomics' session sub-dict shape.

    Entirely separate from the existing Metabolic Connection state
    (input_df/mapped_df at the top level, /api/mc/*): that feature answers
    "does this compound occur in organism X", this one runs enrichment. They
    deliberately share no state and no code path - see Metabolomics_module's
    docstring and Metabolomics_retriever's cache-isolation comment.

    TWO independent upload paths, mirroring how _fresh_proteomics_state
    carries the protein-level and PTM site-level paths side by side:
      - the LIST path (a plain compound-name list, no statistics) is SHARED
        by three tiles - Over-Representation Analysis, Chemical Class
        Enrichment, Lipid Class / Structural Enrichment - since all three
        need only ONE resolve+crosswalk-match of the same uploaded list.
      - the CONCENTRATION path (samples x compounds, Omiscera computes the
        statistics) belongs to Quantitative Enrichment alone - it is a
        different upload with a different compound set, so it cannot share
        the list path's resolved table.
    Both paths share one organism selection (taxon_id/taxon_name) and one
    process-cached crosswalk (Metabolomics_module.build_compound_crosswalk
    caches itself; nothing here duplicates that).
    """
    return {
        "taxon_id": "", "taxon_name": "",
        # List path - shared by ORA / Chemical Class / Lipid Class.
        "list_df": None, "list_column_mapping": None,
        "list_resolved_df": None, "list_matched_df": None, "list_meta": None,
        "background_compounds": None, "background_meta": None,
        "ora_config": None, "ora_results": None,
        "chemclass_config": None, "chemclass_results": None,
        "lipidclass_config": None, "lipidclass_results": None,
        # Concentration-table path - Quantitative Enrichment only.
        "conc_values_df": None, "conc_groups": None, "conc_meta": None,
        "conc_resolved_df": None, "conc_matched_df": None,
        "quant_stats_df": None, "quant_stats_meta": None,
        "quant_config": None, "quant_results": None,
        # Pathway Topology Analysis - a 5th tile that reuses EITHER upload
        # path above (its own config says which): "list" reads
        # list_matched_df exactly like ORA; "conc" reads conc_matched_df and
        # computes its OWN per-compound statistics (topology_stats_df/meta)
        # via the same compute_compound_statistics call metabolomics_quant_
        # configure makes, kept separate from quant_stats_df so the two
        # tiles' group/metric choices never need to agree.
        "topology_stats_df": None, "topology_stats_meta": None,
        "topology_config": None, "topology_results": None,
    }


def new_session_state() -> dict:
    return {
        "input_df": None,
        "duplicate_report_df": pd.DataFrame(),
        "mapped_df": None,
        "hmdb_import_by_rank": {},
        "analysis_results": {},
        "ora": _fresh_ora_state(),
        "transcriptomics": _fresh_transcriptomics_state(),
        "proteomics": _fresh_proteomics_state(),
        "metabolomics": _fresh_metabolomics_state(),
        "multiomics": _fresh_multiomics_state(),
    }


def get_session(session_id: str) -> dict:
    with _SESSIONS_LOCK:
        if session_id not in SESSIONS:
            SESSIONS[session_id] = new_session_state()
        return SESSIONS[session_id]


def session_dep(x_session_id: Optional[str] = Header(default=None)) -> dict:
    if not x_session_id:
        raise HTTPException(status_code=400, detail="Missing X-Session-Id header.")
    return get_session(x_session_id)


# ---------------------------------------------------------------------------
# Job store (background work: compound mapping, connection analysis,
# gene resolution, BLAST) - mirrors the Streamlit spinners/progress bars.
# ---------------------------------------------------------------------------

_JOBS_LOCK = threading.Lock()
JOBS: dict[str, dict[str, Any]] = {}
EXECUTOR = ThreadPoolExecutor(max_workers=4)


def create_job() -> str:
    job_id = uuid.uuid4().hex
    with _JOBS_LOCK:
        JOBS[job_id] = {
            "status": "queued",
            "progress": 0.0,
            "message": "",
            "result": None,
            "error": None,
        }
    return job_id


def update_job(job_id: str, **fields) -> None:
    with _JOBS_LOCK:
        if job_id in JOBS:
            JOBS[job_id].update(fields)


def run_job(job_id: str, fn, *args, **kwargs) -> None:
    update_job(job_id, status="running")
    try:
        result = fn(job_id, *args, **kwargs)
        update_job(job_id, status="done", progress=1.0, result=result)
    except Exception as exc:
        update_job(job_id, status="error", error=str(exc))


def submit_job(fn, *args, **kwargs) -> str:
    job_id = create_job()
    EXECUTOR.submit(run_job, job_id, fn, *args, **kwargs)
    return job_id


# ---------------------------------------------------------------------------
# Job worker functions
# ---------------------------------------------------------------------------

def _job_map_compounds(job_id: str, session: dict) -> dict:
    update_job(job_id, message="Running compound mapping...")
    mapped_df = resolver_module.resolve_compound_table(session["input_df"].copy(), verbose=False)
    session["mapped_df"] = mapped_df
    session["analysis_results"] = {}
    counts = mapped_df["Mapping_status"].value_counts().to_dict()
    return {
        "table": df_payload(order_mapped_columns(mapped_df)),
        "counts": {
            "mapped": counts.get("mapped", 0),
            "partial": counts.get("partial", 0),
            "not_found": counts.get("not_found", 0),
        },
    }


def _job_run_analysis(job_id: str, session: dict, rank_settings: dict) -> dict:
    update_job(job_id, message="Running connection analysis...")
    results = run_connection_analysis(
        session["mapped_df"], rank_settings, session["hmdb_import_by_rank"],
    )
    session["analysis_results"] = results
    payload = {}
    for rank_name, result in results.items():
        sel_rows = [
            {
                "Database": db,
                "Context ID": info.get("context_id", ""),
                "Context Name": info.get("context_name", ""),
                "Context Type": info.get("context_type", ""),
            }
            for db, info in result.get("selections", {}).items()
        ]
        payload[rank_name] = {
            "selections": sel_rows,
            "table1": df_payload(result.get("table1")),
            "table2": df_payload(result.get("table2")),
            "table3": df_payload(result.get("table3")),
        }
    return payload


def _job_resolve_genes(job_id: str, session: dict, identifier_type: str, taxon_id: str) -> dict:
    update_job(job_id, message="Resolving gene identifiers via UniProt...")
    gene_df = session["ora"]["gene_df"]
    resolved_df = gene_resolver.resolve_gene_table(
        gene_df,
        identifier_type,
        taxon_id=taxon_id if identifier_type == "gene_symbol" else None,
        verbose=False,
    )
    session["ora"]["resolved_df"] = resolved_df
    session["ora"]["final_df"] = None

    id_col = "Gene Symbol" if identifier_type == "gene_symbol" else "Ensembl Gene ID"
    ambiguous_mask = resolved_df["Resolution_status"] == "ambiguous"
    table1_df = resolved_df[~ambiguous_mask].copy()
    ambiguous_df = resolved_df[ambiguous_mask].copy()
    display_cols = [c for c in resolved_df.columns if c != "Candidates"]

    session["ora"]["table1_df"] = table1_df[display_cols].copy()
    session["ora"]["display_cols"] = display_cols
    session["ora"]["id_col"] = id_col

    status_counts = table1_df["Resolution_status"].value_counts().to_dict()
    metrics = {
        "resolved": status_counts.get("resolved", 0),
        "resolved_unreviewed": status_counts.get("resolved_unreviewed", 0),
        "resolved_via_ncbi": status_counts.get("resolved_via_ncbi", 0),
        "user_provided": status_counts.get("user_provided", 0),
        "not_found": status_counts.get("not_found", 0),
    }

    ambiguous_list = []
    for _, amb_row in ambiguous_df.iterrows():
        gene_id = str(amb_row.get(id_col, ""))
        candidates_raw = str(amb_row.get("Candidates", ""))
        try:
            candidates_list = json.loads(candidates_raw) if candidates_raw else []
        except Exception:
            candidates_list = []
        ambiguous_list.append({"gene_id": gene_id, "candidates": candidates_list})

    auto_final = None
    if ambiguous_df.empty:
        final_df = table1_df[display_cols].copy()
        session["ora"]["final_df"] = final_df
        auto_final = df_payload(final_df)

    return {
        "id_col": id_col,
        "resolved_table": df_payload(table1_df[display_cols]),
        "metrics": metrics,
        "ambiguous": ambiguous_list,
        "auto_final": auto_final,
    }


def _job_resolve_kegg_ids(job_id: str, session: dict, kegg_organism_code: str) -> dict:
    """Optional, opt-in second retrieval pass: assigns each gene a KEGG gene
    ID (Genomics_module_retriever.resolve_kegg_ids), a SEPARATE identifier
    from UniProt ID that only Genomics ORA/Burden's KEGG database matches on
    (see build_kegg_universe's docstring for why KEGG is not accession-keyed).

    Runs on top of whatever /api/ora/genes/resolve already produced - it can
    ALSO resolve genes that were "ambiguous" under UniProt alone (KEGG's own
    accession map disambiguating pending candidates, or an NCBI Entrez ID
    bridge), which is why this re-derives table1_df/ambiguous_df exactly like
    _job_resolve_genes rather than only adding a column.
    """
    update_job(job_id, message=f"Resolving KEGG gene IDs for organism '{kegg_organism_code}'...")
    resolved_df = session["ora"].get("resolved_df")
    id_col = session["ora"].get("id_col")
    if resolved_df is None or not id_col:
        raise ValueError("Resolve gene identifiers first.")

    mapped_df = gene_resolver.resolve_kegg_ids(resolved_df, id_col, kegg_organism_code, verbose=False)
    session["ora"]["resolved_df"] = mapped_df
    session["ora"]["kegg_organism_code"] = kegg_organism_code

    ambiguous_mask = mapped_df["Resolution_status"] == "ambiguous"
    table1_df = mapped_df[~ambiguous_mask].copy()
    ambiguous_df = mapped_df[ambiguous_mask].copy()
    display_cols = [c for c in mapped_df.columns if c != "Candidates"]

    session["ora"]["table1_df"] = table1_df[display_cols].copy()
    session["ora"]["display_cols"] = display_cols

    status_counts = table1_df["Resolution_status"].value_counts().to_dict()
    metrics = {
        "resolved": status_counts.get("resolved", 0),
        "resolved_unreviewed": status_counts.get("resolved_unreviewed", 0),
        "resolved_via_ncbi": status_counts.get("resolved_via_ncbi", 0),
        "resolved_via_kegg": status_counts.get("resolved_via_kegg", 0),
        "user_provided": status_counts.get("user_provided", 0),
        "not_found": status_counts.get("not_found", 0),
    }

    ambiguous_list = []
    for _, amb_row in ambiguous_df.iterrows():
        gene_id = str(amb_row.get(id_col, ""))
        candidates_raw = str(amb_row.get("Candidates", ""))
        try:
            candidates_list = json.loads(candidates_raw) if candidates_raw else []
        except Exception:
            candidates_list = []
        ambiguous_list.append({"gene_id": gene_id, "candidates": candidates_list})

    auto_final = None
    if ambiguous_df.empty:
        final_df = table1_df[display_cols].copy()
        session["ora"]["final_df"] = final_df
        auto_final = df_payload(final_df)
    else:
        # KEGG resolution can change which rows are ambiguous, so any final_df
        # set BEFORE this step is now stale (it could miss newly-recovered
        # rows) - clear it; the user re-applies selections / discard-all on
        # the updated ambiguous set.
        session["ora"]["final_df"] = None

    return {
        "resolved_table": df_payload(table1_df[display_cols]),
        "metrics": metrics,
        "ambiguous": ambiguous_list,
        "auto_final": auto_final,
        "n_kegg_mapped": int((mapped_df["KEGG ID"] != "").sum()),
    }


def _job_resolve_de(job_id: str, session: dict, identifier_type: str, taxon_id: str, taxon_name: str) -> dict:
    """Transcriptomics retrieval worker: resolves the uploaded DE table via
    Transcriptomics_module_retriever.resolve_de_table (bulk mygene.info /
    Ensembl REST identifier mapping - never per-gene, see that module).

    taxon_name is committed HERE, not in the endpoint. It used to be written
    synchronously by /api/transcriptomics/resolve before this job ran, while
    taxon_id was only written below on success - and run_job does no session
    rollback on failure. A failed re-resolve therefore left the session
    holding the NEW organism's name beside the OLD organism's id and resolved
    table, and since taxon_name drives Reactome while taxon_id drives GO,
    MSigDB and the BioMart gene-length correction, the next analysis silently
    mixed two organisms in one result set. Committing both here, together
    with resolved_df/meta, makes a failed resolve a complete no-op.
    """
    update_job(job_id, message="Resolving identifiers via mygene.info...")
    de_df = session["transcriptomics"]["de_df"]

    def _progress(current, total, label):
        frac = current / max(total, 1)
        update_job(job_id, progress=frac, message=f"{label} ({current}/{total})")

    resolved_df, meta = transcriptomics_retriever.resolve_de_table(
        de_df, taxon_id, identifier_type=identifier_type,
        progress_callback=_progress, verbose=False,
    )
    session["transcriptomics"]["taxon_id"] = taxon_id
    session["transcriptomics"]["taxon_name"] = taxon_name
    session["transcriptomics"]["resolved_df"] = resolved_df
    session["transcriptomics"]["meta"] = meta
    # A KEGG organism confirmation (see _fresh_transcriptomics_state) was
    # made against the PREVIOUS organism, if any - committing a new one here
    # invalidates it. _transcriptomics_kegg_confirmed_code would catch this
    # anyway via the stored taxon_id/taxon_name mismatch, but clearing it
    # here too means a stale confirmation never lingers in session state
    # after the organism it was confirmed for is gone.
    session["transcriptomics"]["kegg_organism_confirmation"] = None

    update_job(job_id, progress=1.0, message="Done")
    return {
        "resolved_table": df_payload(resolved_df),
        "metrics": meta,
    }


def _serialize_transcriptomics_ora_tables(tables: dict) -> dict:
    """JSON-safe view of Transcriptomics_module.run_ora's 'tables' dict.

    The un-serialized version (DataFrames) stays in session state; only this
    payload crosses the API boundary. Mirrors _serialize_ora_results' shape
    for the Genomics ORA endpoint, flattened to one level since Transcriptomics
    result keys (e.g. "GO: Biological Process (all)") are already unique.
    """
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block["n_terms_tested"],
            "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
            "n_terms_no_overlap": block.get("n_terms_no_overlap", 0),
            "universe_size": block["universe_size"],
            "n_deg_in_universe": block["n_deg_in_universe"],
            "n_universe_total": block.get("n_universe_total", block["universe_size"]),
            "n_deg_total": block.get("n_deg_total", block["n_deg_in_universe"]),
            "length_correction_used": block.get("length_correction_used", False),
            "length_correction_reason": block.get("length_correction_reason", ""),
            "note": block.get("note"),
        }
    return payload


def _kegg_confirmed_code(session: dict, module: str) -> str:
    """The KEGG organism code a user manually confirmed for THIS session's
    CURRENT organism in the given module ("transcriptomics" or
    "proteomics"), or "" if none applies.

    Belt-and-braces re-check, mirroring the existing kegg_taxon_mismatch
    guard pattern used for Genomics ORA's kegg_organism_code (~API.py:1540):
    never trust stored session state without re-verifying it against the
    run's current organism. _job_resolve_de/_job_resolve_proteins already
    clear kegg_organism_confirmation whenever a new organism is committed,
    but a confirmation could still be stale if that invalidation were ever
    missed - so the stored taxon_id/taxon_name are checked again here, at the
    one point that actually matters (the moment a run is about to use it),
    not only at the one point that currently sets it.

    Generalized from the original _transcriptomics_kegg_confirmed_code (now
    used by Proteomics too - same organism-resolution/confirmation contract,
    see Proteomics_module.build_kegg_accession_sets' docstring for why
    Proteomics needed no picker changes to reuse it, only a different
    session key).
    """
    confirmation = session[module].get("kegg_organism_confirmation")
    if not confirmation:
        return ""
    if (
        confirmation.get("taxon_id") != session[module].get("taxon_id")
        or confirmation.get("taxon_name") != session[module].get("taxon_name")
    ):
        return ""
    return str(confirmation.get("organism_code") or "")


def _job_run_transcriptomics_ora(job_id: str, session: dict) -> dict:
    """Transcriptomics ORA worker: runs Transcriptomics_module.run_ora against
    the retriever's resolved DE table using the stored ora_config."""
    update_job(job_id, message="Starting over-representation analysis...")
    resolved_df = session["transcriptomics"]["resolved_df"]
    if resolved_df is None:
        raise ValueError("No resolved differential-expression table available.")

    config = session["transcriptomics"]["ora_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    result = transcriptomics_engine.run_ora(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["transcriptomics"]["taxon_id"],
        taxon_name=session["transcriptomics"].get("taxon_name", ""),
        padj_cutoff=config["padj_cutoff"],
        lfc_cutoff=config["lfc_cutoff"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        use_length_correction=config["use_length_correction"],
        include_iea=config.get("include_iea", True),
        custom_db=session["transcriptomics"].get("custom_db"),
        background_symbols=(
            session["transcriptomics"].get("background_symbols") if config.get("use_background") else None
        ),
        msigdb_collections=config.get("msigdb_collections"),
        kegg_confirmed_organism_code=_kegg_confirmed_code(session, "transcriptomics"),
        progress_callback=_progress,
        verbose=False,
    )
    session["transcriptomics"]["ora_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_transcriptomics_ora_tables(result["tables"]),
        "summary": result["summary"],
    }


def _job_resolve_proteins(job_id: str, session: dict, taxon_id: str, taxon_name: str) -> dict:
    """Proteomics retrieval worker: resolves the uploaded protein-level DE
    table via Proteomics_module_retriever.resolve_protein_table (UniProt ID
    Mapping - never per-accession, see that module). Mirrors _job_resolve_de's
    commit-both-together discipline: taxon_id/taxon_name/resolved_df/meta are
    all written here, together, on success only, so a failed re-resolve
    cannot leave the session holding a new organism's name beside an old
    organism's resolved table (see _job_resolve_de's docstring for the exact
    failure mode this prevents)."""
    update_job(job_id, message="Resolving accessions via UniProt ID Mapping...")
    de_df = session["proteomics"]["de_df"]

    def _progress(current, total, label):
        frac = current / max(total, 1)
        update_job(job_id, progress=frac, message=f"{label} ({current}/{total})")

    resolved_df, meta = proteomics_retriever.resolve_protein_table(
        de_df, taxon_id, progress_callback=_progress, verbose=False,
    )
    session["proteomics"]["taxon_id"] = taxon_id
    session["proteomics"]["taxon_name"] = taxon_name
    session["proteomics"]["resolved_df"] = resolved_df
    session["proteomics"]["meta"] = meta
    # A KEGG organism confirmation (see _fresh_proteomics_state) was made
    # against the PREVIOUS organism, if any - committing a new one here
    # invalidates it. Mirrors _job_resolve_de's identical reset.
    session["proteomics"]["kegg_organism_confirmation"] = None

    update_job(job_id, progress=1.0, message="Done")
    return {
        "resolved_table": df_payload(resolved_df),
        "metrics": meta,
    }


def _serialize_proteomics_ora_tables(tables: dict) -> dict:
    """JSON-safe view of Proteomics_module.run_ora's 'tables' dict. Mirrors
    _serialize_transcriptomics_ora_tables' shape; no length-correction keys
    since Proteomics ORA has no PWF/gene-length-bias concept in this
    milestone."""
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block["n_terms_tested"],
            "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
            "n_terms_no_overlap": block.get("n_terms_no_overlap", 0),
            "universe_size": block["universe_size"],
            "n_input_in_universe": block.get("n_input_in_universe", 0),
            "n_input_dropped": block.get("n_input_dropped", 0),
            "caveat": block.get("caveat"),
            "note": block.get("note"),
        }
    return payload


def _job_run_proteomics_ora(job_id: str, session: dict) -> dict:
    """Proteomics ORA worker: runs Proteomics_module.run_ora against the
    retriever's resolved protein table using the stored ora_config."""
    update_job(job_id, message="Starting over-representation analysis...")
    resolved_df = session["proteomics"]["resolved_df"]
    if resolved_df is None:
        raise ValueError("No resolved protein table available.")

    config = session["proteomics"]["ora_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    result = proteomics_engine.run_ora(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["proteomics"]["taxon_id"],
        taxon_name=session["proteomics"].get("taxon_name", ""),
        padj_cutoff=config["padj_cutoff"],
        lfc_cutoff=config["lfc_cutoff"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        include_iea=config.get("include_iea", True),
        custom_db=session["proteomics"].get("custom_db"),
        background_accessions=(
            session["proteomics"].get("background_accessions") if config.get("use_background") else None
        ),
        msigdb_collections=config.get("msigdb_collections"),
        kegg_confirmed_organism_code=_kegg_confirmed_code(session, "proteomics"),
        progress_callback=_progress,
        verbose=False,
    )
    session["proteomics"]["ora_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_proteomics_ora_tables(result["tables"]),
        "summary": result["summary"],
    }


def _serialize_proteomics_quant_tables(tables: dict) -> dict:
    """JSON-safe view of Proteomics_module.run_quantitative_enrichment's
    'tables' dict. Mirrors _serialize_proteomics_ora_tables' shape, adjusted
    for the Burden-Test-engine block shape (no n_terms_no_overlap concept -
    the Wilcoxon engine has no analogous zero-overlap exclusion, see that
    engine's docstring)."""
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block["n_terms_tested"],
            "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
            "universe_size": block["universe_size"],
            "n_input_in_universe": block.get("n_input_in_universe", 0),
            "n_input_dropped": block.get("n_input_dropped", 0),
            "caveat": block.get("caveat"),
            "note": block.get("note"),
        }
    return payload


def _job_run_proteomics_quant(job_id: str, session: dict) -> dict:
    """Proteomics Quantitative Enrichment worker: runs
    Proteomics_module.run_quantitative_enrichment against the retriever's
    resolved protein table using the stored quant_config. Mirrors
    _job_run_proteomics_ora's shape."""
    update_job(job_id, message="Starting quantitative enrichment...")
    resolved_df = session["proteomics"]["resolved_df"]
    if resolved_df is None:
        raise ValueError("No resolved protein table available.")

    config = session["proteomics"]["quant_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    result = proteomics_engine.run_quantitative_enrichment(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["proteomics"]["taxon_id"],
        taxon_name=session["proteomics"].get("taxon_name", ""),
        metric=config["metric"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        include_iea=config.get("include_iea", True),
        custom_db=session["proteomics"].get("custom_db"),
        background_accessions=(
            session["proteomics"].get("background_accessions") if config.get("use_background") else None
        ),
        msigdb_collections=config.get("msigdb_collections"),
        kegg_confirmed_organism_code=_kegg_confirmed_code(session, "proteomics"),
        progress_callback=_progress,
        verbose=False,
    )
    session["proteomics"]["quant_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_proteomics_quant_tables(result["tables"]),
        "summary": result["summary"],
    }


# ---------------------------------------------------------------------------
# Metabolomics — Over-Representation Analysis, Chemical Class Enrichment,
# Lipid Class / Structural Enrichment, and Quantitative Enrichment.
#
# TWO upload paths (see _fresh_metabolomics_state's own comment): a plain
# compound-name LIST, shared by the first three tiles, and a CONCENTRATION
# TABLE for Quantitative Enrichment, whose statistics Omiscera computes
# itself (never a user-supplied log2FC/pvalue - that RNA-seq/proteomics
# convention does not fit what a metabolomics researcher actually uploads).
# Both paths funnel through the SAME resolve+crosswalk-match pipeline
# (_resolve_and_match_compound_names) since the underlying problem - turn a
# compound NAME into this module's internal identifier - is identical either
# way; only what supplies the names differs.
# ---------------------------------------------------------------------------

def _compute_user_id_mask(input_df: pd.DataFrame) -> pd.Series:
    """True on rows where the USER'S OWN upload carried a structure
    identifier (InChIKey/PubChem ID/ChEBI ID/LIPID MAPS ID).

    Computed BEFORE resolving, while the table still shows only what the
    user actually supplied: the retriever fills these same columns in from
    name lookups, after which a user-supplied structure is indistinguishable
    from a retriever-guessed one. Only matters for lipid shorthand - see
    Metabolomics_module.match_compounds' user_id_mask parameter.
    """
    user_id_columns = [c for c in ("InChIKey", "PubChem ID", "ChEBI ID", "LIPID MAPS ID")
                       if c in input_df.columns]
    if not user_id_columns:
        return pd.Series([False] * len(input_df))
    return (
        input_df[user_id_columns].astype(str)
        .apply(lambda col: col.str.strip().replace({"nan": "", "None": ""}) != "")
        .any(axis=1)
    )


def _dedupe_metabolomics_names(input_df: pd.DataFrame) -> "tuple[pd.DataFrame, int]":
    """Remove an exact-duplicate compound Name from a plain LIST upload
    (ORA / Chemical Class / Lipid Class), keeping the first occurrence -
    metabolomics_retriever.remove_duplicate_input_names, reused unchanged.

    A repeated Name in a plain list (no numeric column attached) is a
    data-entry duplicate, not a second measurement - unlike the concentration
    table, where a repeated Name genuinely is a second feature/adduct of the
    same compound and must NOT be collapsed (see _job_resolve_metabolites_conc,
    which deliberately never calls this).

    Returns (deduped_df, n_removed).
    """
    deduped = metabolomics_retriever.remove_duplicate_input_names(input_df)
    return deduped, len(input_df) - len(deduped)


def _resolve_and_match_compound_names(
    job_id: str, input_df: pd.DataFrame, progress_start: float = 0.0, progress_span: float = 1.0,
    dedupe_names: bool = True,
) -> "tuple[pd.DataFrame, pd.DataFrame, dict]":
    """Shared by every Metabolomics resolve job: resolve compound Name(+IDs)
    via the retriever, build/reuse the crosswalk, then match into crosswalk
    space. Both the list path and the concentration-table path (whose
    "names" are its own column headers) need this exact three-step
    pipeline - only what gets uploaded differs, never how it gets resolved.

    progress_start/progress_span let a caller reserve part of the job's
    progress bar for work before/after this call (e.g. the concentration
    path's own parse step already ran before this begins).

    dedupe_names: remove an exact-duplicate Name before resolving (the list
    path's default). The concentration path passes False - see
    _dedupe_metabolomics_names' docstring.

    Returns (resolved_df, matched_df, meta) - meta carries match_compounds'
    own tier/skeleton/class-level-only counts plus mapping_status_counts, plus
    n_input_rows (the ORIGINAL upload count, before any dedup) and
    n_duplicate_names_removed.
    """
    def _frac(local):
        return progress_start + progress_span * local

    n_uploaded = len(input_df)
    n_duplicate_names_removed = 0
    if dedupe_names:
        input_df, n_duplicate_names_removed = _dedupe_metabolomics_names(input_df)

    update_job(job_id, progress=_frac(0.0), message="Resolving compound identifiers...")
    user_id_mask = _compute_user_id_mask(input_df)
    resolved_df = metabolomics_retriever.resolve_compound_table(input_df, verbose=False)
    for column in ("ChEBI ID", "LIPID MAPS ID"):
        if column in input_df.columns and column not in resolved_df.columns:
            resolved_df[column] = input_df[column].values

    update_job(job_id, progress=_frac(0.35), message="Building compound crosswalk...")

    def _cw_progress(step, total, label):
        # The crosswalk occupies the middle of this call's progress share.
        update_job(job_id, progress=_frac(0.35 + 0.50 * (step / max(total, 1))), message=label)

    crosswalk = metabolomics_engine.build_compound_crosswalk(
        verbose=False, progress_callback=_cw_progress)

    update_job(job_id, progress=_frac(0.9), message="Matching compounds...")
    matched_df, match_meta = metabolomics_engine.match_compounds(
        resolved_df, crosswalk, user_id_mask=user_id_mask)

    meta = dict(match_meta)
    meta["n_input_rows"] = n_uploaded
    meta["n_duplicate_names_removed"] = n_duplicate_names_removed
    meta["mapping_status_counts"] = (
        resolved_df["Mapping_status"].value_counts().to_dict()
        if "Mapping_status" in resolved_df.columns else {}
    )
    return resolved_df, matched_df, meta


def _job_resolve_metabolites_list(job_id: str, session: dict, taxon_id: str, taxon_name: str) -> dict:
    """List-path retrieval worker, shared by ORA / Chemical Class / Lipid
    Class. Mirrors _job_resolve_proteins' commit-both-together discipline -
    taxon_id/taxon_name/list_resolved_df/list_matched_df/list_meta are
    written together, on success only, so a failed re-resolve can never
    leave the session holding a new organism beside an old organism's
    matched table."""
    input_df = session["metabolomics"]["list_df"]
    if input_df is None:
        raise ValueError("No compound list uploaded.")

    taxon_changed = (
        taxon_id != session["metabolomics"].get("taxon_id", "")
        or taxon_name != session["metabolomics"].get("taxon_name", "")
    )
    resolved_df, matched_df, meta = _resolve_and_match_compound_names(job_id, input_df)
    meta["taxon_id"] = taxon_id
    meta["taxon_name"] = taxon_name

    session["metabolomics"]["taxon_id"] = taxon_id
    session["metabolomics"]["taxon_name"] = taxon_name
    session["metabolomics"]["list_resolved_df"] = resolved_df
    session["metabolomics"]["list_matched_df"] = matched_df
    session["metabolomics"]["list_meta"] = meta
    # A new resolve invalidates any analysis built on the previous one, for
    # all three tiles that share this list - PLUS Pathway Topology, which
    # reuses this same list_matched_df when its own input mode is "list".
    session["metabolomics"]["ora_results"] = None
    session["metabolomics"]["chemclass_results"] = None
    session["metabolomics"]["lipidclass_results"] = None
    session["metabolomics"]["topology_results"] = None
    if taxon_changed:
        # The list and concentration paths share ONE organism. Quant's KEGG /
        # Reactome / SMPDB / PathBank universes are built from that same
        # taxon_id/taxon_name (run_quantitative_enrichment), so committing a
        # NEW one here makes Quant's already-computed results describe an
        # organism that is no longer current - clear them too, not only this
        # path's own. quant_stats_df/meta are untouched: per-compound
        # statistics do not depend on the organism, only the enrichment
        # universes built from them do.
        session["metabolomics"]["quant_results"] = None

    update_job(job_id, progress=1.0, message="Done")
    return {"resolved_table": df_payload(matched_df), "metrics": meta}


def _job_resolve_background_compounds(job_id: str, session: dict, bg_df: pd.DataFrame) -> dict:
    """Resolve an optional background compound list (list tiles only) into
    a set of compound_keys, used to narrow a database's reference universe
    to "compounds that could plausibly have been detected" rather than
    every compound the database has ever annotated."""
    _resolved_df, matched_df, _meta = _resolve_and_match_compound_names(job_id, bg_df)
    bg_keys = set(
        matched_df[matched_df["compound_key"].astype(str) != ""]["compound_key"].astype(str)
    )
    background_meta = {"n_uploaded": len(bg_df), "n_matched": len(bg_keys)}
    session["metabolomics"]["background_compounds"] = bg_keys
    session["metabolomics"]["background_meta"] = background_meta
    update_job(job_id, progress=1.0, message="Done")
    return {"metrics": background_meta}


def _serialize_metabolomics_tables(tables: dict) -> dict:
    """JSON-safe view of run_ora's / run_quantitative_enrichment's /
    run_pathway_topology's 'tables'.

    One serializer for all five tiles: every block differs only in which
    optional keys it carries, and every read here is already a .get() with
    a default. "measure" and "graphs" are Pathway Topology-only (None/{} on
    every other tile's block, harmless).
    """
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block.get("n_terms_tested", 0),
            "n_terms_excluded_min_size": block.get("n_terms_excluded_min_size", 0),
            "n_terms_no_overlap": block.get("n_terms_no_overlap", 0),
            "universe_size": block.get("universe_size", 0),
            "n_input_in_universe": block.get("n_input_in_universe", 0),
            "n_input_dropped": block.get("n_input_dropped", 0),
            "caveat": block.get("caveat"),
            "note": block.get("note"),
            "measure": block.get("measure"),
            "graphs": block.get("graphs") or {},
        }
    return payload


def _job_run_metabolomics_list_tile(job_id: str, session: dict, config_key: str,
                                    results_key: str, databases: "list[str]") -> dict:
    """Shared runner for ORA / Chemical Class / Lipid Class - all three call
    Metabolomics_module.run_ora with a different `databases` list (see that
    function's own docstring for why one function serves all three)."""
    update_job(job_id, message="Starting analysis...")
    matched_df = session["metabolomics"]["list_matched_df"]
    if matched_df is None:
        raise ValueError("No resolved compound list available.")
    config = session["metabolomics"][config_key]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    crosswalk = metabolomics_engine.build_compound_crosswalk(verbose=False)
    result = metabolomics_engine.run_ora(
        matched_df, crosswalk, databases=databases,
        taxon_id=session["metabolomics"].get("taxon_id", ""),
        taxon_name=session["metabolomics"].get("taxon_name", ""),
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        chemical_class_level=config.get("chemical_class_level", "sub_class"),
        lipid_class_level=config.get("lipid_class_level", "lm_main_class"),
        pathbank_primary_only=config.get("pathbank_primary_only", True),
        background_compounds=(
            session["metabolomics"].get("background_compounds") if config.get("use_background") else None
        ),
        kegg_organism_code=config.get("kegg_organism_code", ""),
        kegg_organism_name=config.get("kegg_organism_name", ""),
        progress_callback=_progress,
        verbose=False,
    )
    session["metabolomics"][results_key] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_metabolomics_tables(result["tables"]),
        "summary": result["summary"],
        "universe_caveats": result.get("universe_caveats", {}),
    }


def _job_run_metabolomics_ora(job_id: str, session: dict) -> dict:
    """Over-Representation Analysis worker (KEGG/Reactome/SMPDB/PathBank)."""
    config = session["metabolomics"].get("ora_config") or {}
    return _job_run_metabolomics_list_tile(
        job_id, session, "ora_config", "ora_results",
        config.get("databases", ["KEGG", "Reactome", "SMPDB", "PathBank"]))


def _job_run_metabolomics_chemclass(job_id: str, session: dict) -> dict:
    """Chemical Class Enrichment worker (RefMet taxonomy)."""
    return _job_run_metabolomics_list_tile(
        job_id, session, "chemclass_config", "chemclass_results", ["Chemical Class"])


def _job_run_metabolomics_lipidclass(job_id: str, session: dict) -> dict:
    """Lipid Class / Structural Enrichment worker (LIPID MAPS taxonomy)."""
    return _job_run_metabolomics_list_tile(
        job_id, session, "lipidclass_config", "lipidclass_results", ["Lipid Class"])


def _job_resolve_metabolites_conc(job_id: str, session: dict, taxon_id: str, taxon_name: str) -> dict:
    """Concentration-path retrieval worker (Quantitative Enrichment only).
    The "compound list" here is the concentration table's OWN column
    headers - each column is one compound - resolved through the exact same
    pipeline as the list-upload tiles."""
    values_df = session["metabolomics"]["conc_values_df"]
    if values_df is None:
        raise ValueError("No concentration table uploaded.")
    taxon_changed = (
        taxon_id != session["metabolomics"].get("taxon_id", "")
        or taxon_name != session["metabolomics"].get("taxon_name", "")
    )
    names_df = pd.DataFrame({"Name": [str(c) for c in values_df.columns]})

    # dedupe_names=False: a repeated Name here is two FEATURES of one compound
    # (e.g. two adducts) - each is its own column in values_df, and removing a
    # "duplicate" row would desynchronize matched_df from those columns.
    resolved_df, matched_df, meta = _resolve_and_match_compound_names(
        job_id, names_df, dedupe_names=False)
    meta["taxon_id"] = taxon_id
    meta["taxon_name"] = taxon_name

    session["metabolomics"]["taxon_id"] = taxon_id
    session["metabolomics"]["taxon_name"] = taxon_name
    session["metabolomics"]["conc_resolved_df"] = resolved_df
    session["metabolomics"]["conc_matched_df"] = matched_df
    session["metabolomics"]["conc_meta"] = meta
    session["metabolomics"]["quant_stats_df"] = None
    session["metabolomics"]["quant_stats_meta"] = None
    session["metabolomics"]["quant_results"] = None
    # Same reasoning as quant_stats_df/quant_results above, for Pathway
    # Topology's own concentration-mode statistics/results.
    session["metabolomics"]["topology_stats_df"] = None
    session["metabolomics"]["topology_stats_meta"] = None
    session["metabolomics"]["topology_results"] = None
    if taxon_changed:
        # Chemical Class / Lipid Class never use the organism at all (see
        # Metabolomics_module._build_universes - their universes are built
        # from crosswalk taxonomy only), so they are left alone. ORA's KEGG /
        # Reactome / SMPDB / PathBank universes DO use it, so a taxon commit
        # made HERE (on the concentration path) invalidates ORA's now-stale
        # cached results - the mirror image of the list path's guard above.
        session["metabolomics"]["ora_results"] = None

    update_job(job_id, progress=1.0, message="Done")
    return {"resolved_table": df_payload(matched_df), "metrics": meta}


def _job_run_metabolomics_quant(job_id: str, session: dict) -> dict:
    """Quantitative Enrichment worker (competitive two-sided Wilcoxon),
    scored from the per-compound statistics quant_stats_df already computed
    at configure time - see metabolomics_quant_configure."""
    update_job(job_id, message="Starting quantitative enrichment...")
    matched_df = session["metabolomics"]["conc_matched_df"]
    stats_df = session["metabolomics"]["quant_stats_df"]
    if matched_df is None or stats_df is None:
        raise ValueError("Resolve and configure the analysis first.")
    config = session["metabolomics"]["quant_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    crosswalk = metabolomics_engine.build_compound_crosswalk(verbose=False)
    result = metabolomics_engine.run_quantitative_enrichment(
        matched_df, crosswalk, stats_df,
        databases=config["databases"],
        taxon_id=session["metabolomics"].get("taxon_id", ""),
        taxon_name=session["metabolomics"].get("taxon_name", ""),
        metric=config["metric"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        chemical_class_level=config.get("chemical_class_level", "sub_class"),
        lipid_class_level=config.get("lipid_class_level", "lm_main_class"),
        pathbank_primary_only=config.get("pathbank_primary_only", True),
        background_compounds=(
            session["metabolomics"].get("background_compounds") if config.get("use_background") else None
        ),
        kegg_organism_code=config.get("kegg_organism_code", ""),
        kegg_organism_name=config.get("kegg_organism_name", ""),
        progress_callback=_progress,
        verbose=False,
    )
    session["metabolomics"]["quant_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    # Per-compound log2 fold change keyed by the user's own Name (the same
    # names the results' "Overlapping Compounds" column carries) - feeds the
    # compound-concept network's node colouring, exactly as Proteomics'
    # per-protein log2FC does. Always log2fc regardless of the ranking
    # metric, so the colour keeps its plain biological meaning.
    compound_log2fc = {
        str(name): float(value)
        for name, value in zip(stats_df["Name"], stats_df["log2fc"])
        if pd.notna(value) and math.isfinite(float(value))
    }
    return {
        "tables": _serialize_metabolomics_tables(result["tables"]),
        "summary": result["summary"],
        "universe_caveats": result.get("universe_caveats", {}),
        "compound_log2fc": compound_log2fc,
    }


def _job_run_metabolomics_topology(job_id: str, session: dict) -> dict:
    """Pathway Topology Analysis worker (KEGG only). Reads config["input_mode"]
    to decide which of the two shared upload paths to run against - "list"
    (list_matched_df, exactly like ORA) or "conc" (conc_matched_df plus this
    tile's OWN topology_stats_df, computed at configure time exactly as
    metabolomics_quant_configure computes quant_stats_df)."""
    update_job(job_id, message="Starting pathway topology analysis...")
    mx = session["metabolomics"]
    config = mx.get("topology_config")
    if config is None:
        raise ValueError("Configure the analysis first.")
    input_mode = config["input_mode"]
    if input_mode == "list":
        matched_df = mx["list_matched_df"]
        stats_df, metric = None, metabolomics_engine.QUANT_METRIC_DEFAULT
        if matched_df is None:
            raise ValueError("No resolved compound list available.")
    else:
        matched_df = mx["conc_matched_df"]
        stats_df = mx["topology_stats_df"]
        metric = config.get("metric", metabolomics_engine.QUANT_METRIC_DEFAULT)
        if matched_df is None or stats_df is None:
            raise ValueError("Resolve and configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    crosswalk = metabolomics_engine.build_compound_crosswalk(verbose=False)
    result = metabolomics_engine.run_pathway_topology(
        matched_df, crosswalk, input_mode=input_mode,
        stats_df=stats_df, metric=metric,
        measure=config.get("measure", metabolomics_engine.PATHWAY_TOPOLOGY_MEASURE_DEFAULT),
        taxon_id=mx.get("taxon_id", ""), taxon_name=mx.get("taxon_name", ""),
        min_term_size=config["min_term_size"], max_term_size=config["max_term_size"],
        background_compounds=(mx.get("background_compounds") if config.get("use_background") else None),
        kegg_organism_code=config.get("kegg_organism_code", ""),
        kegg_organism_name=config.get("kegg_organism_name", ""),
        progress_callback=_progress,
        verbose=False,
    )
    mx["topology_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_metabolomics_tables(result["tables"]),
        "summary": result["summary"],
        "universe_caveats": result.get("universe_caveats", {}),
    }


# ---------------------------------------------------------------------------
# Proteomics — PTM / Phosphosite Enrichment
# (a SEPARATE site-level retrieval path from ORA/Quant/Network above - see
# _fresh_proteomics_state's own comment. Shares taxon_id/taxon_name once
# resolved, but resolved_sites_df is independent of resolved_df.)
# ---------------------------------------------------------------------------

def _job_resolve_ptm_sites(job_id: str, session: dict, taxon_id: str, taxon_name: str) -> dict:
    """PTM retrieval worker: resolves the uploaded SITE-level table via
    Proteomics_ptm_retriever.resolve_site_table. Mirrors _job_resolve_proteins'
    commit-both-together discipline. Deliberately does NOT touch resolved_df/
    ora_results/quant_results/network_results - those belong to the
    protein-level pipeline and are untouched by a PTM upload/resolve."""
    update_job(job_id, message="Resolving phosphosites via UniProt ID Mapping...")
    sites_df = session["proteomics"]["ptm_sites_df"]

    def _progress(current, total, label):
        frac = current / max(total, 1)
        update_job(job_id, progress=frac, message=f"{label} ({current}/{total})")

    resolved_sites_df, meta = proteomics_ptm_retriever.resolve_site_table(
        sites_df, taxon_id, progress_callback=_progress, verbose=False,
    )
    session["proteomics"]["taxon_id"] = taxon_id
    session["proteomics"]["taxon_name"] = taxon_name
    session["proteomics"]["resolved_sites_df"] = resolved_sites_df
    session["proteomics"]["ptm_meta"] = meta

    update_job(job_id, progress=1.0, message="Done")
    return {
        "resolved_table": df_payload(resolved_sites_df),
        "metrics": meta,
    }


def _serialize_ptm_activity_block(block: dict) -> dict:
    """JSON-safe view of run_kinase_activity's return - single block (one
    database, iPTMnet), wrapped under a "iPTMnet" key for the SAME {db:
    block} shape _serialize_transcriptomics_regulon_activity_tables already
    uses, so the results screen can reuse that rendering pattern almost
    unchanged rather than inventing a second one."""
    if block.get("table") is None:
        return {"iPTMnet": {
            "table": df_payload(pd.DataFrame(columns=proteomics_ptm_engine.ACTIVITY_RESULT_COLUMNS)),
            "n_kinases_tested": 0, "note": block.get("note"),
            "null_slopes_by_signature": {}, "kinase_regression_detail": {},
            "ranking_values_sample": [], "y_mean": 0.0, "n_sites": 0, "metric": "",
            "lambda_used": 1.0, "null_center": 0.0, "calibration_note": None,
        }}
    return {"iPTMnet": {
        "table": df_payload(block["table"]),
        "n_kinases_tested": block["n_kinases_tested"],
        "note": block.get("note"),
        "null_slopes_by_signature": block.get("null_slopes_by_signature", {}),
        "kinase_regression_detail": block.get("kinase_regression_detail", {}),
        "ranking_values_sample": block.get("ranking_values_sample", []),
        "y_mean": block.get("y_mean", 0.0),
        "n_sites": block.get("n_sites", 0),
        "metric": block.get("metric", ""),
        "lambda_used": block.get("lambda_used", 1.0),
        "null_center": block.get("null_center", 0.0),
        "calibration_note": block.get("calibration_note"),
    }}


def _serialize_ptm_ora_tables(tables: dict) -> dict:
    """JSON-safe view of run_kinase_substrate_ora's 'tables' dict. Mirrors
    _serialize_proteomics_ora_tables' shape (same underlying engine,
    tx._run_ora_for_sets)."""
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block.get("n_terms_tested", 0),
            "n_terms_excluded_min_size": block.get("n_terms_excluded_min_size", 0),
            "n_terms_no_overlap": block.get("n_terms_no_overlap", 0),
            "universe_size": block.get("universe_size", 0),
            "note": block.get("note"),
        }
    return payload


def _job_run_ptm(job_id: str, session: dict) -> dict:
    """Runs kinase activity (A1) then kinase substrate ORA (A2) in one job -
    same reasoning as _job_run_transcriptomics_regulon_both: a single
    progress bar, and both analyses share the same iPTMnet universe fetch so
    running them together costs almost nothing extra once A1 has paid for
    the fetch."""
    update_job(job_id, message="Starting kinase activity inference...")
    resolved_sites_df = session["proteomics"]["resolved_sites_df"]
    if resolved_sites_df is None:
        raise ValueError("No resolved phosphosite table available.")
    config = session["proteomics"]["ptm_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    taxon_id = session["proteomics"]["taxon_id"]
    taxon_name = session["proteomics"].get("taxon_name", "")

    def _progress_activity(step, total, label):
        update_job(job_id, progress=0.5 * step / max(total, 1), message=label)

    activity_block = proteomics_ptm_engine.run_kinase_activity(
        resolved_sites_df, taxon_id, taxon_name=taxon_name,
        metric=config["metric"], min_kinase_sites=config["min_kinase_sites"],
        min_score=config.get("min_score"), n_permutations=config["n_permutations"],
        seed=config["seed"], progress_callback=_progress_activity, verbose=False,
    )
    session["proteomics"]["ptm_activity_results"] = activity_block

    update_job(job_id, progress=0.5, message="Starting kinase substrate enrichment...")
    ora_result = proteomics_ptm_engine.run_kinase_substrate_ora(
        resolved_sites_df, taxon_id, taxon_name=taxon_name,
        padj_cutoff=config["padj_cutoff"], lfc_cutoff=config["lfc_cutoff"],
        min_term_size=config["min_term_size"], max_term_size=config.get("max_term_size"),
        min_score=config.get("min_score"), verbose=False,
    )
    session["proteomics"]["ptm_ora_results"] = ora_result

    update_job(job_id, progress=1.0, message="Done")
    return {
        "activity": {"tables": _serialize_ptm_activity_block(activity_block)},
        "ora": {"tables": _serialize_ptm_ora_tables(ora_result["tables"]), "universe_size": ora_result["universe_size"],
               "caveat": ora_result.get("caveat")},
    }


def _serialize_transcriptomics_gsea_tables(tables: dict) -> dict:
    """JSON-safe view of Transcriptomics_module.run_gsea's 'tables' dict.

    Mirrors _serialize_transcriptomics_ora_tables' shape; GSEA blocks carry
    different summary keys (no length-correction concept - that's an ORA-only
    idea) plus per-term enrichment curves for the frontend's three-panel
    enrichment plot.
    """
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block["n_terms_tested"],
            "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
            "universe_size": block["universe_size"],
            "curves": block.get("curves", {}),
            "null_nes_by_size": block.get("null_nes_by_size", {}),
            "note": block.get("note"),
        }
    return payload


def _job_run_transcriptomics_gsea(job_id: str, session: dict) -> dict:
    """Transcriptomics GSEA worker: runs Transcriptomics_module.run_gsea
    against the retriever's resolved DE table using the stored gsea_config.
    custom_db/background_symbols are shared with ORA (same session keys) -
    both stages read the same upload widgets, which is safe because any new
    file upload (POST /api/transcriptomics/upload) fully resets the whole
    'transcriptomics' session sub-dict before either stage's configure step
    can run again."""
    update_job(job_id, message="Starting GSEA...")
    resolved_df = session["transcriptomics"]["resolved_df"]
    if resolved_df is None:
        raise ValueError("No resolved differential-expression table available.")

    config = session["transcriptomics"]["gsea_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")

    def _progress(step, total, label):
        update_job(job_id, progress=step / max(total, 1), message=label)

    result = transcriptomics_engine.run_gsea(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["transcriptomics"]["taxon_id"],
        taxon_name=session["transcriptomics"].get("taxon_name", ""),
        metric=config["metric"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        n_permutations=config["n_permutations"],
        weight=config["weight"],
        seed=config["seed"],
        include_iea=config.get("include_iea", True),
        custom_db=session["transcriptomics"].get("custom_db"),
        background_symbols=(
            session["transcriptomics"].get("background_symbols") if config.get("use_background") else None
        ),
        msigdb_collections=config.get("msigdb_collections"),
        kegg_confirmed_organism_code=_kegg_confirmed_code(session, "transcriptomics"),
        progress_callback=_progress,
        verbose=False,
    )
    session["transcriptomics"]["gsea_results"] = result
    update_job(job_id, progress=1.0, message="Done")
    return {
        "tables": _serialize_transcriptomics_gsea_tables(result["tables"]),
        "ranking_profile": result["ranking_profile"],
        "ranking_full": result["ranking_full"],
        "summary": result["summary"],
    }


def _serialize_transcriptomics_regulon_enrichment_tables(tables: dict) -> dict:
    """JSON-safe view of run_regulon_enrichment's 'tables' dict - same block
    shape as _serialize_transcriptomics_ora_tables (both wrap
    Transcriptomics_module._run_ora_for_sets), minus the length-correction
    keys this module never sets (see run_regulon_enrichment's docstring)."""
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_terms_tested": block["n_terms_tested"],
            "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
            "n_terms_no_overlap": block.get("n_terms_no_overlap", 0),
            "universe_size": block["universe_size"],
            "n_deg_in_universe": block.get("n_deg_in_universe", 0),
            "note": block.get("note"),
        }
    return payload


def _serialize_transcriptomics_regulon_activity_tables(tables: dict) -> dict:
    """JSON-safe view of run_regulon_activity's 'tables' dict. null_slopes_by_
    signature/tf_regression_detail are already JSON-safe (str keys, plain
    floats/lists) as returned - see run_regulon_activity's docstring for why
    only the top REGULON_MAX_CURVES_STORED TFs per database carry these."""
    payload = {}
    for key, block in tables.items():
        payload[key] = {
            "table": df_payload(block["table"]),
            "n_tfs_tested": block["n_tfs_tested"],
            "note": block.get("note"),
            "null_slopes_by_signature": block.get("null_slopes_by_signature", {}),
            "tf_regression_detail": block.get("tf_regression_detail", {}),
            # The regression figure needs the values the FIT actually used
            # (ranking metric, not log2FC) plus the fit's own centroid - see
            # run_regulon_activity's tf_regression_detail comment.
            "ranking_values_sample": block.get("ranking_values_sample", []),
            "y_mean": block.get("y_mean", 0.0),
            "n_genes": block.get("n_genes", 0),
            "metric": block.get("metric", ""),
            # Empirical-null calibration diagnostics - see
            # run_regulon_activity's docstring and _empirical_null_calibration.
            "lambda_used": block.get("lambda_used", 1.0),
            "null_center": block.get("null_center", 0.0),
            "calibration_note": block.get("calibration_note"),
        }
    return payload


def _regulon_session_inputs(session: dict):
    """Shared guard for every regulon worker below."""
    resolved_df = session["transcriptomics"]["resolved_df"]
    if resolved_df is None:
        raise ValueError("No resolved differential-expression table available.")
    config = session["transcriptomics"]["regulon_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")
    return resolved_df, config


def _regulon_enrichment_payload(job_id: str, session: dict, progress_offset: float = 0.0,
                                 progress_scale: float = 1.0) -> dict:
    """Run TF target enrichment and return the FLAT {tables, summary} payload.

    The custom-db/background widgets are NOT used here (regulons are their own
    gene-set source) - only resolved_df/taxon_id and regulon_config are read.
    Callers wrap this under the "enrichment" key; see
    _job_run_transcriptomics_regulon_enrichment for why the envelope shape
    matters.
    """
    update_job(job_id, message="Starting TF target enrichment...")
    resolved_df, config = _regulon_session_inputs(session)

    def _progress(step, total, label):
        update_job(job_id, progress=progress_offset + progress_scale * step / max(total, 1), message=label)

    result = transcriptomics_regulon_engine.run_regulon_enrichment(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["transcriptomics"]["taxon_id"],
        padj_cutoff=config["padj_cutoff"],
        lfc_cutoff=config["lfc_cutoff"],
        min_term_size=config["min_term_size"],
        max_term_size=config["max_term_size"],
        dorothea_levels=config["dorothea_levels"],
        chea3_libraries=config.get("chea3_libraries"),
        background_symbols=(
            session["transcriptomics"].get("background_symbols") if config.get("use_background") else None
        ),
        progress_callback=_progress,
        verbose=False,
    )
    session["transcriptomics"]["regulon_enrichment_results"] = result
    return {
        "tables": _serialize_transcriptomics_regulon_enrichment_tables(result["tables"]),
        "summary": result["summary"],
    }


def _regulon_activity_payload(job_id: str, session: dict, progress_offset: float = 0.0,
                               progress_scale: float = 1.0) -> dict:
    """Run TF activity inference (ULM) and return the FLAT {tables, summary}
    payload. Callers wrap this under the "activity" key."""
    update_job(job_id, message="Starting TF activity inference...")
    resolved_df, config = _regulon_session_inputs(session)

    def _progress(step, total, label):
        update_job(job_id, progress=progress_offset + progress_scale * step / max(total, 1), message=label)

    result = transcriptomics_regulon_engine.run_regulon_activity(
        resolved_df,
        databases=config["databases"],
        taxon_id=session["transcriptomics"]["taxon_id"],
        metric=config["metric"],
        min_regulon_size=config["min_regulon_size"],
        dorothea_levels=config["dorothea_levels"],
        n_permutations=config["n_permutations"],
        seed=config["seed"],
        progress_callback=_progress,
        verbose=False,
    )
    session["transcriptomics"]["regulon_activity_results"] = result
    return {
        "tables": _serialize_transcriptomics_regulon_activity_tables(result["tables"]),
        "summary": result["summary"],
    }


# The three workers below MUST all return the same envelope shape -
# {"enrichment": {...}} and/or {"activity": {...}} - because the results
# screen reads `results.enrichment` / `results.activity` and renders one pane
# per key it finds (static/app.js, renderTranscriptomicsRegulonResultsScreen).
#
# They did not. The two single-analysis workers used to return a FLAT
# {tables, summary}, so selecting only ONE analysis on the config screen gave
# a completely BLANK results screen - both panes skipped, empty section rail,
# no error anywhere. Every test had configured BOTH analyses, so the path was
# never exercised. Verified against the live API before and after this fix.

def _job_run_transcriptomics_regulon_enrichment(job_id: str, session: dict) -> dict:
    payload = {"enrichment": _regulon_enrichment_payload(job_id, session)}
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _job_run_transcriptomics_regulon_activity(job_id: str, session: dict) -> dict:
    payload = {"activity": _regulon_activity_payload(job_id, session)}
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _job_run_transcriptomics_regulon_both(job_id: str, session: dict) -> dict:
    """Runs enrichment then activity in one job, so the results screen for
    'both' shows a single progress bar instead of two separate job polls -
    each half driving its own 0-50% / 50-100% slice of the bar."""
    payload = {
        "enrichment": _regulon_enrichment_payload(job_id, session, progress_offset=0.0, progress_scale=0.5),
        "activity": _regulon_activity_payload(job_id, session, progress_offset=0.5, progress_scale=0.5),
    }
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _job_run_blast(job_id: str, session: dict) -> dict:
    gene_df = session["ora"]["gene_df"]

    def _progress(current, total, label):
        frac = current / max(total, 1)
        update_job(job_id, progress=frac, message=f"{label} ({current}/{total})")

    update_job(job_id, message="Submitting...")
    resolved_df = gene_resolver.resolve_sequence_table(
        gene_df, taxon_id=None, progress_callback=_progress, verbose=False,
    )
    session["ora"]["resolved_df"] = resolved_df
    session["ora"]["final_df"] = None

    status_counts = resolved_df["Resolution_status"].value_counts().to_dict()
    metrics = {
        "resolved_homology": status_counts.get("resolved_homology", 0),
        "user_provided": status_counts.get("user_provided", 0),
        "not_found": status_counts.get("not_found", 0),
        "blast_failed": status_counts.get("blast_failed", 0),
        "blast_pending": status_counts.get("blast_pending", 0),
    }
    return {"table": df_payload(resolved_df), "metrics": metrics}


def _term_size_bounds(cfg: dict) -> tuple:
    """Optional min/max term-size filter from a per_database payload entry.

    min defaults to the engine's MIN_TERM_SIZE (3). max distinguishes "key
    absent from the payload" (-> genomics_ora.MAX_TERM_SIZE_DEFAULT, 500 -
    clusterProfiler's own default; needed now that GO ancestor propagation
    can otherwise let a near-root term with tens of thousands of genes
    through uncapped) from "explicitly null/empty" (-> no cap, a deliberate
    opt-out), mirroring transcriptomics_ora_configure's identical
    distinction. The running app always sends the key (see
    ORA_DEFAULT_MAX_TERM_SIZE in static/app.js), so this mainly protects a
    direct API caller that omits it.
    Shared by _job_run_ora and _job_run_burden.
    """
    try:
        min_ts = int(cfg.get("min_term_size"))
        if min_ts < 1:
            min_ts = genomics_ora.MIN_TERM_SIZE
    except (TypeError, ValueError):
        min_ts = genomics_ora.MIN_TERM_SIZE
    _MISSING = object()
    max_raw = cfg.get("max_term_size", _MISSING)
    if max_raw is _MISSING:
        max_ts = genomics_ora.MAX_TERM_SIZE_DEFAULT
    elif max_raw in (None, ""):
        max_ts = None
    else:
        try:
            max_ts = int(max_raw)
            if max_ts < min_ts:
                max_ts = None
        except (TypeError, ValueError):
            max_ts = None
    return min_ts, max_ts


def _append_caveat(block: dict, message: str) -> None:
    """Append a caveat message onto a result block, preserving any existing one.

    Shared by _job_run_ora and _job_run_burden's KEGG organism-mismatch and
    dropped-background guards.
    """
    existing = block.get("caveat")
    block["caveat"] = (existing + " " + message) if existing else message


def _job_run_ora(job_id: str, session: dict) -> dict:
    update_job(job_id, message="Starting over-representation analysis...")
    final_df = session["ora"]["final_df"]
    if final_df is None or "UniProt ID" not in final_df.columns:
        raise ValueError("No confirmed UniProt ID list available for ORA.")

    accessions = [a for a in final_df["UniProt ID"].astype(str).tolist() if a.strip()]
    # Gene symbols are available when the input type was Gene Symbol.
    # GWAS Catalog uses these directly (symbol-native, no accession bridge).
    symbols: list = []
    fallback_symbols: list = []
    if "Gene Symbol" in final_df.columns:
        symbols = [s for s in final_df["Gene Symbol"].astype(str).tolist() if s.strip()]
        # Symbols of genes that FAILED UniProt resolution (empty accession) — used by
        # GO as a fallback so hard-to-resolve genes (esp. non-model species) still test.
        uid = final_df["UniProt ID"].astype(str).str.strip()
        sym = final_df["Gene Symbol"].astype(str).str.strip()
        fallback_symbols = [s for s, u in zip(sym.tolist(), uid.tolist()) if s and not u]
    # Ensembl Gene IDs are available when the input type was Ensembl - used only
    # by an "ensembl"-typed custom database (no built-in database is Ensembl-native).
    ensembl_ids: list = []
    if "Ensembl Gene ID" in final_df.columns:
        ensembl_ids = [e for e in final_df["Ensembl Gene ID"].astype(str).tolist() if e.strip()]

    # KEGG gene IDs are available only after Genomics_module_retriever.resolve_kegg_ids
    # has run (an opt-in step after the normal UniProt/Ensembl resolution - see
    # /api/ora/genes/resolve-kegg). KEGG matches on THESE, never on UniProt ID -
    # see build_kegg_universe's docstring for why the universe is KEGG-gene-ID-keyed.
    kegg_ids: list = []
    kegg_id_labels: dict = {}
    if "KEGG ID" in final_df.columns:
        label_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else (
            "Ensembl Gene ID" if "Ensembl Gene ID" in final_df.columns else None
        )
        for _, row in final_df.iterrows():
            kid = str(row.get("KEGG ID", "")).strip()
            if not kid:
                continue
            kegg_ids.append(kid)
            if label_col:
                label = str(row.get(label_col, "")).strip()
                if label:
                    kegg_id_labels[kid] = label

    config = session["ora"]["ora_config"]
    databases = config["databases"]
    per_database = config["per_database"]
    background_by_db = session["ora"].get("background_by_db") or {}

    all_results = {}
    for i, db in enumerate(databases):
        db_cfg = per_database[db]
        min_ts, max_ts = _term_size_bounds(db_cfg)

        # Custom (user-uploaded) gene-set database — reuses the shared engine.
        # Dispatched BEFORE reading taxon_id: its per_database payload has no
        # taxon_id/taxon_name (species-agnostic), only id_type.
        if db == "Custom database":
            custom = session["ora"].get("custom_db")
            if not custom:
                all_results[db] = {db: {
                    "table": pd.DataFrame(columns=genomics_ora.ORA_RESULT_COLUMNS),
                    "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
                    "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
                    "caveat": "No custom gene-set file was uploaded.",
                }}
                continue

            def _cprogress(current, total, label, _db=db, _i=i):
                frac = (_i + current / max(total, 1)) / len(databases)
                update_job(job_id, progress=frac, message=f"{_db}: {label}")

            all_results[db] = genomics_ora.run_ora_for_custom_db(
                term_map=custom["terms"],
                id_type=custom["id_type"],
                list_accessions=accessions,
                list_symbols=symbols,
                list_ensembl=ensembl_ids,
                db_label=db,
                min_term_size=min_ts,
                max_term_size=max_ts,
                progress_callback=_cprogress,
                verbose=False,
            )
            continue

        taxon_id = db_cfg["taxon_id"]
        taxon_name = db_cfg.get("taxon_name", "")
        include_iea = bool(db_cfg.get("include_iea", True))

        # Server-side correctness guard: KEGG's organism must be the one the
        # KEGG gene IDs were actually resolved for (session["ora"]
        # ["kegg_organism_code"], set by /api/ora/genes/resolve-kegg), never
        # whatever the config screen's taxon_id happens to hold. The two used
        # to be independently settable via a second KEGG organism picker on
        # the config screen (removed in static/app.js - see that fix), which
        # could silently diverge from the retrieval-stage organism and test
        # "hsa" gene IDs against an "mmu" universe: an empty result, no error.
        # The client is not the correctness boundary, so this is checked here
        # too, not only prevented in the UI. If kegg_ids is non-empty,
        # kegg_organism_code is guaranteed set - the same job
        # (_job_resolve_kegg_ids) writes both together.
        kegg_taxon_mismatch = False
        if db == "KEGG":
            kegg_org = session["ora"].get("kegg_organism_code", "")
            if kegg_org and kegg_org != taxon_id:
                kegg_taxon_mismatch = True
                taxon_id = kegg_org

        # User-supplied background list (GO/Reactome/KEGG) lives in the session,
        # keyed by database; None -> use the whole annotated-genome universe.
        # Always uploaded/parsed as UniProt accessions (load_user_background_list),
        # so for KEGG it must be mapped into KEGG gene-ID space through the SAME
        # gene map resolve_kegg_ids used - a raw accession background would
        # silently match nothing against a KEGG-gene-ID-keyed universe.
        background = set(background_by_db.get(db) or []) or None
        background_kegg_ids = None
        kegg_background_dropped = False
        if db == "KEGG" and background:
            kegg_org_for_bg = session["ora"].get("kegg_organism_code", "")
            if kegg_org_for_bg:
                bg_gene_map = genomics_ora.build_kegg_gene_map(kegg_org_for_bg)
                bg_acc_to_gene = bg_gene_map.get("acc_to_gene", {})
                background_kegg_ids = {bg_acc_to_gene[a] for a in background if a in bg_acc_to_gene} or None
            else:
                kegg_background_dropped = True

        # Defensive gate: GWAS needs symbols; skip gracefully if unavailable.
        if db == "GWAS Catalog" and not symbols:
            empty_table = pd.DataFrame(columns=genomics_ora.ORA_RESULT_COLUMNS)
            all_results[db] = {"GWAS Catalog": {
                "table": empty_table,
                "n_terms_tested": 0,
                "n_terms_excluded_min_size": 0,
                "universe_size": 0,
                "n_input_in_universe": 0,
                "n_input_dropped": 0,
                "caveat": "GWAS Catalog requires Gene Symbol input. No symbols available for this run.",
            }}
            continue

        # Defensive gate: KEGG needs a "KEGG ID" per gene, a SEPARATE opt-in
        # step from normal resolution (see /api/ora/genes/resolve-kegg). A
        # user who selected KEGG without ever running that step gets a clear,
        # actionable caveat here rather than a silently-empty KEGG table.
        if db == "KEGG" and not kegg_ids:
            empty_table = pd.DataFrame(columns=genomics_ora.ORA_RESULT_COLUMNS)
            all_results[db] = {"KEGG": {
                "table": empty_table,
                "n_terms_tested": 0,
                "n_terms_excluded_min_size": 0,
                "universe_size": 0,
                "n_input_in_universe": 0,
                "n_input_dropped": 0,
                "caveat": (
                    "KEGG requires KEGG gene IDs, resolved separately from UniProt. "
                    "Go back to the retrieval step, choose a KEGG organism and resolve "
                    "KEGG IDs, then re-run."
                ),
            }}
            continue

        def _progress(current, total, label, _db=db, _i=i):
            frac = (_i + current / max(total, 1)) / len(databases)
            update_job(job_id, progress=frac, message=f"{_db}: {label}")

        result = genomics_ora.run_ora_for_database(
            database=db,
            taxon_id=taxon_id,
            list_accessions=accessions,
            list_symbols=symbols,
            fallback_symbols=fallback_symbols,
            background_accessions=background,
            include_iea=include_iea,
            taxon_name=taxon_name,
            min_term_size=min_ts,
            max_term_size=max_ts,
            list_kegg_ids=kegg_ids,
            background_kegg_ids=background_kegg_ids,
            kegg_id_labels=kegg_id_labels,
            progress_callback=_progress,
            verbose=False,
        )
        if db == "KEGG" and "KEGG" in result:
            if kegg_taxon_mismatch:
                _append_caveat(result["KEGG"],
                    f"The configured organism did not match the one KEGG IDs were "
                    f"resolved for; '{taxon_id}' (the resolved organism) was used instead.")
            if kegg_background_dropped:
                _append_caveat(result["KEGG"],
                    "A background gene list was uploaded but could not be mapped to "
                    "KEGG gene space (no KEGG organism on record) - it was ignored "
                    "for this run.")
        all_results[db] = result

    session["ora"]["ora_results"] = all_results
    return _serialize_ora_results(all_results)


def _serialize_ora_results(all_results: dict) -> dict:
    """JSON-safe view of run_ora_for_database's nested output.

    The un-serialized version (DataFrames, etc.) stays in session state;
    only this payload crosses the API boundary.
    """
    payload = {}
    for db, db_result in all_results.items():
        payload[db] = {}
        for key, block in db_result.items():
            payload[db][key] = {
                "table": df_payload(block["table"]),
                "n_terms_tested": block["n_terms_tested"],
                "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
                "universe_size": block["universe_size"],
                "n_input_in_universe": block["n_input_in_universe"],
                "n_input_dropped": block["n_input_dropped"],
                "n_recovered_by_symbol": block.get("n_recovered_by_symbol", 0),
                "caveat": block.get("caveat"),
            }
    return payload


def _serialize_burden_results(all_results: dict, gene_scores_by_db: dict) -> dict:
    """JSON-safe view of the burden run's nested output.

    Same block shape as _serialize_ora_results, plus a per-block "gene_scores"
    map (gene id -> oriented variant score) so the client's gene-concept network
    can color each gene node by its score. The map is the same for every block
    of a database (e.g. GO's three namespaces share score_by_acc); it holds only
    the user's scored genes, so the duplication is negligible.
    """
    payload = {}
    for db, db_result in all_results.items():
        payload[db] = {}
        scores = gene_scores_by_db.get(db, {})
        for key, block in db_result.items():
            payload[db][key] = {
                "table": df_payload(block["table"]),
                "n_terms_tested": block["n_terms_tested"],
                "n_terms_excluded_min_size": block["n_terms_excluded_min_size"],
                "universe_size": block["universe_size"],
                "n_input_in_universe": block["n_input_in_universe"],
                "n_input_dropped": block["n_input_dropped"],
                "caveat": block.get("caveat"),
                "gene_scores": scores,
            }
    return payload


def _job_run_burden(job_id: str, session: dict) -> dict:
    """Gene-set Burden Test worker: mirrors _job_run_ora's shape/dispatch, but
    scores are continuous (variant_scores) rather than binary membership, and
    the statistic is the competitive Wilcoxon rank-sum test (see
    Genomics_module.run_burden_for_database / run_burden_for_custom_db).
    Retrieval (final_df) is fully shared with ORA - only the analysis stage
    differs.
    """
    update_job(job_id, message="Starting gene-set burden test...")
    final_df = session["ora"]["final_df"]
    if final_df is None or "UniProt ID" not in final_df.columns:
        raise ValueError("No confirmed UniProt ID list available for the burden test.")

    variant_scores = session["ora"].get("variant_scores") or {}
    if not variant_scores:
        raise ValueError(
            "No variant score column was captured for this gene table. "
            "Re-upload with a score column selected."
        )
    # Orient scores so higher always means "more burden" for the one-sided
    # Wilcoxon test; a raw metric where lower = more signal (e.g. a p-value)
    # is negated here rather than in the engine, keeping the engine agnostic.
    sign = -1.0 if session["ora"].get("score_direction") == "lower" else 1.0

    id_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else (
        "Ensembl Gene ID" if "Ensembl Gene ID" in final_df.columns else None
    )

    def _aggregate(pairs) -> dict:
        buckets: dict = {}
        for key, raw in pairs:
            buckets.setdefault(key, []).append(raw)
        return {key: sign * (sum(vals) / len(vals)) for key, vals in buckets.items()}

    # Build score maps in every identifier space a database might need:
    # accession (GO/Reactome, and a Custom DB typed "accession"), plus symbol/
    # ensembl for a Custom DB typed to match those input spaces (mirrors how
    # _job_run_ora builds accessions/symbols/ensembl_ids from the same final_df).
    score_by_acc, score_by_symbol, score_by_ensembl = {}, {}, {}
    if id_col:
        uid_series = final_df["UniProt ID"].astype(str).str.strip()
        gid_series = final_df[id_col].astype(str).str.strip()
        acc_pairs, sym_pairs, ens_pairs = [], [], []
        for uid, gid in zip(uid_series.tolist(), gid_series.tolist()):
            gid_key = gid.lower()
            if not gid_key or gid_key not in variant_scores:
                continue
            raw_score = variant_scores[gid_key]
            if uid:
                acc_pairs.append((uid, raw_score))
            if id_col == "Gene Symbol":
                sym_pairs.append((gid.upper(), raw_score))
            elif id_col == "Ensembl Gene ID":
                ens_pairs.append((gid.upper(), raw_score))
        score_by_acc = _aggregate(acc_pairs)
        score_by_symbol = _aggregate(sym_pairs)
        score_by_ensembl = _aggregate(ens_pairs)

    # KEGG gene IDs are available only after resolve_kegg_ids has run (see
    # _job_run_ora's identical block for the full rationale). Built the same
    # way as the accession/symbol/ensembl maps above: one score per KEGG gene,
    # aggregated if the same KEGG gene ID somehow appears on more than one row.
    score_by_kegg_id: dict = {}
    kegg_id_labels: dict = {}
    if "KEGG ID" in final_df.columns and id_col:
        kegg_series = final_df["KEGG ID"].astype(str).str.strip()
        gid_series_for_kegg = final_df[id_col].astype(str).str.strip()
        uid_series_for_kegg = final_df["UniProt ID"].astype(str).str.strip()
        kegg_pairs = []
        for kid, gid, uid in zip(kegg_series.tolist(), gid_series_for_kegg.tolist(), uid_series_for_kegg.tolist()):
            if not kid:
                continue
            gid_key = gid.lower()
            if gid_key and gid_key in variant_scores:
                kegg_pairs.append((kid, variant_scores[gid_key]))
            if gid:
                kegg_id_labels[kid] = gid
        score_by_kegg_id = _aggregate(kegg_pairs)

    config = session["ora"]["burden_config"]
    databases = config["databases"]
    per_database = config["per_database"]

    all_results = {}
    # Per-database gene->oriented-score map, sent to the client so the
    # gene-concept network can color each gene node by its variant score.
    gene_scores_by_db: dict = {}
    for i, db in enumerate(databases):
        db_cfg = per_database[db]
        min_ts, max_ts = _term_size_bounds(db_cfg)

        # Custom (user-uploaded) gene-set database — reuses the shared engine,
        # dispatched before reading taxon_id (species-agnostic, only id_type).
        if db == "Custom database":
            custom = session["ora"].get("custom_db")
            if not custom:
                all_results[db] = {db: {
                    "table": pd.DataFrame(columns=genomics_ora.BURDEN_RESULT_COLUMNS),
                    "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
                    "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
                    "caveat": "No custom gene-set file was uploaded.",
                }}
                continue

            score_by_id = {
                "accession": score_by_acc, "symbol": score_by_symbol, "ensembl": score_by_ensembl,
            }.get(custom["id_type"], {})
            gene_scores_by_db[db] = score_by_id

            def _cprogress(current, total, label, _db=db, _i=i):
                frac = (_i + current / max(total, 1)) / len(databases)
                update_job(job_id, progress=frac, message=f"{_db}: {label}")

            all_results[db] = genomics_ora.run_burden_for_custom_db(
                term_map=custom["terms"],
                score_by_id=score_by_id,
                db_label=db,
                min_term_size=min_ts,
                max_term_size=max_ts,
                progress_callback=_cprogress,
                verbose=False,
            )
            continue

        taxon_id = db_cfg["taxon_id"]
        taxon_name = db_cfg.get("taxon_name", "")
        include_iea = bool(db_cfg.get("include_iea", True))

        # Server-side correctness guard: same reasoning as _job_run_ora's
        # identical block - trust session["ora"]["kegg_organism_code"] (the
        # organism KEGG gene IDs were actually resolved for), never the config
        # screen's taxon_id, which used to be independently settable via a
        # second KEGG organism picker (removed - see static/app.js).
        kegg_taxon_mismatch = False
        if db == "KEGG":
            kegg_org = session["ora"].get("kegg_organism_code", "")
            if kegg_org and kegg_org != taxon_id:
                kegg_taxon_mismatch = True
                taxon_id = kegg_org

        # Defensive gate: KEGG needs a "KEGG ID" per gene (see _job_run_ora's
        # identical gate for the full rationale) - a separate opt-in step from
        # normal resolution.
        if db == "KEGG" and not score_by_kegg_id:
            gene_scores_by_db[db] = {}
            all_results[db] = {"KEGG": {
                "table": pd.DataFrame(columns=genomics_ora.BURDEN_RESULT_COLUMNS),
                "n_terms_tested": 0,
                "n_terms_excluded_min_size": 0,
                "universe_size": 0,
                "n_input_in_universe": 0,
                "n_input_dropped": 0,
                "caveat": (
                    "KEGG requires KEGG gene IDs, resolved separately from UniProt. "
                    "Go back to the retrieval step, choose a KEGG organism and resolve "
                    "KEGG IDs, then re-run."
                ),
            }}
            continue

        # GO/Reactome are accession-space, so their gene-node scores come from
        # score_by_acc; KEGG is KEGG-gene-ID-space, so it uses score_by_kegg_id.
        gene_scores_by_db[db] = score_by_kegg_id if db == "KEGG" else score_by_acc

        def _progress(current, total, label, _db=db, _i=i):
            frac = (_i + current / max(total, 1)) / len(databases)
            update_job(job_id, progress=frac, message=f"{_db}: {label}")

        result = genomics_ora.run_burden_for_database(
            database=db,
            taxon_id=taxon_id,
            score_by_acc=score_by_acc,
            include_iea=include_iea,
            taxon_name=taxon_name,
            min_term_size=min_ts,
            max_term_size=max_ts,
            score_by_kegg_id=score_by_kegg_id,
            kegg_id_labels=kegg_id_labels,
            progress_callback=_progress,
            verbose=False,
        )
        if db == "KEGG" and kegg_taxon_mismatch and "KEGG" in result:
            _append_caveat(result["KEGG"],
                f"The configured organism did not match the one KEGG IDs were "
                f"resolved for; '{taxon_id}' (the resolved organism) was used instead.")
        all_results[db] = result

    session["ora"]["burden_results"] = all_results
    return _serialize_burden_results(all_results, gene_scores_by_db)


def _netprop_param_bounds(cfg: dict) -> dict:
    """Normalize the Network Propagation config into concrete engine params:
    the exposed knobs (STRING confidence, FDR, min module size, # permutations)
    validated against their allowed ranges/choices, plus the fixed-but-reported
    defaults (restart probability, MCL inflation, permutation type, STRING
    version). Returned dict doubles as the 'Parameters used' record."""
    def _int(val, default, choices=None, lo=None):
        try:
            v = int(val)
        except (TypeError, ValueError):
            return default
        if choices is not None and v not in choices:
            return default
        if lo is not None and v < lo:
            return default
        return v

    def _float(val, default, lo, hi):
        try:
            v = float(val)
        except (TypeError, ValueError):
            return default
        return v if (lo < v <= hi) else default

    return {
        "confidence": _int(cfg.get("confidence"), genomics_net.DEFAULT_CONFIDENCE, choices=set(genomics_net.CONFIDENCE_CHOICES)),
        "fdr_alpha": _float(cfg.get("fdr_alpha"), genomics_net.DEFAULT_FDR_ALPHA, 0.0, 1.0),
        "min_module_size": _int(cfg.get("min_module_size"), genomics_net.DEFAULT_MIN_MODULE_SIZE, lo=2),
        "n_permutations": _int(cfg.get("n_permutations"), genomics_net.DEFAULT_N_PERMUTATIONS, choices=set(genomics_net.N_PERMUTATIONS_CHOICES)),
        # Fixed, reported for reproducibility (not user-adjustable).
        "restart_prob": genomics_net.DEFAULT_RESTART_PROB,
        "mcl_inflation": genomics_net.DEFAULT_MCL_INFLATION,
        "permutation_type": "degree-preserving",
        "string_version": genomics_net.STRING_VERSION,
        "enrichment_database": "GO",
        "module_max_term_size": genomics_net.DEFAULT_MODULE_MAX_TERM_SIZE,
    }


def _job_run_netprop(job_id: str, session: dict) -> dict:
    """Network Propagation / Module Detection worker. Retrieval (final_df) is
    shared with ORA/Burden; the optional variant-score column (same capture as
    Burden) sets seed strengths, else seeds are equal-weight. Results are
    ADDITIVE: the propagation block is always produced; the module block only
    when module detection was selected.
    """
    update_job(job_id, message="Starting network propagation...")
    final_df = session["ora"]["final_df"]
    if final_df is None or "UniProt ID" not in final_df.columns:
        raise ValueError("No confirmed UniProt ID list available for network propagation.")

    config = session["ora"]["netprop_config"]
    taxon_id = str(config.get("taxon_id") or session["ora"].get("netprop_taxon_id") or "").strip()
    if not taxon_id:
        raise ValueError("No organism (STRING taxon) selected.")
    do_modules = bool(config.get("module_detection"))
    params = _netprop_param_bounds(config)

    # Build the seed -> UniProt-accession score map from the resolved table.
    # A supplied variant-score column sets seed strength (oriented via the same
    # higher/lower convention as Burden); with no score column, all seeds are
    # equal-weight 1.0 (binary seeds).
    variant_scores = session["ora"].get("variant_scores") or {}
    sign = -1.0 if session["ora"].get("score_direction") == "lower" else 1.0
    id_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else (
        "Ensembl Gene ID" if "Ensembl Gene ID" in final_df.columns else None
    )
    uid_series = final_df["UniProt ID"].astype(str).str.strip()
    seed_scores: dict = {}
    if id_col and variant_scores:
        gid_series = final_df[id_col].astype(str).str.strip()
        buckets: dict = {}
        for uid, gid in zip(uid_series.tolist(), gid_series.tolist()):
            if not uid:
                continue
            gid_key = gid.lower()
            if gid_key in variant_scores:
                buckets.setdefault(uid, []).append(variant_scores[gid_key])
        seed_scores = {uid: sign * (sum(v) / len(v)) for uid, v in buckets.items()}
    if not seed_scores:
        seed_scores = {uid: 1.0 for uid in uid_series.tolist() if uid}
    if not seed_scores:
        raise ValueError("No genes with a UniProt accession to seed the network.")

    def _progress(current, total, label):
        frac = current / max(total, 1)
        update_job(job_id, progress=frac * (0.85 if do_modules else 1.0), message=label)

    update_job(job_id, message=f"Loading STRING network for taxon {taxon_id}...")
    graph = genomics_net.build_string_network(taxon_id, confidence=params["confidence"], verbose=False)

    prop = genomics_net.run_propagation_with_significance(
        graph, seed_scores,
        restart_prob=params["restart_prob"],
        n_permutations=params["n_permutations"],
        fdr_alpha=params["fdr_alpha"],
        progress_callback=_progress,
    )

    payload = {
        "taxon_id": taxon_id,
        "taxon_name": config.get("taxon_name", ""),
        "module_detection": do_modules,
        "parameters": params,
        "propagation": {
            "table": df_payload(prop["table"]),
            "n_seeds_mapped": prop["n_seeds_mapped"],
            "n_seeds_not_in_string": prop["n_seeds_not_in_string"],
            "n_hits": prop["n_hits"],
            "n_recruited": prop["n_recruited"],
            "n_proteins": graph["n_proteins"],
        },
        "modules": None,
    }

    # Edges among the significant hits, for the client-side network figure
    # (only the significant subnetwork is ever drawn, never the ~19k graph).
    hit_ids = prop["table"][prop["table"]["FDR"] < params["fdr_alpha"]]["STRING ID"].tolist()
    payload["propagation"]["edges"] = _netprop_subnetwork_edges(graph, hit_ids)

    if do_modules:
        update_job(job_id, progress=0.88, message="Detecting modules (MCL)...")
        modules = genomics_net.detect_modules(
            hit_ids, graph, min_module_size=params["min_module_size"], inflation=params["mcl_inflation"],
        )
        update_job(job_id, progress=0.92, message="Enriching modules...")
        enr_result = genomics_net.enrich_modules(
            modules, graph, taxon_id, taxon_name=config.get("taxon_name", ""),
            database=params["enrichment_database"], include_iea=True,
            max_term_size=params["module_max_term_size"],
        )
        enr = enr_result["results"]
        module_rows, module_id_map, module_enrichment = [], {}, {}
        for mod in modules:
            genes = [graph["protein_gene"].get(sid, sid) for sid in mod["string_ids"]]
            e = enr.get(mod["module_id"])
            top_term, top_fdr, top_fdr_raw, enr_db = "", None, None, ""
            if e and e["block"] is not None and not e["block"]["table"].empty:
                top = e["block"]["table"].iloc[0]
                top_term = str(top["Term Name"])
                top_fdr_raw = float(top["FDR"])
                # Corrected for having picked the best of N namespaces per
                # module - see Genomics_network.enrich_modules. The raw,
                # within-namespace value is kept alongside for reference.
                top_fdr = float(e["top_fdr_corrected"])
                enr_db = f'{params["enrichment_database"]}:{e["namespace"]}' if e.get("namespace") else params["enrichment_database"]
                # Ship the module's full enrichment table (trimmed to the top 50
                # rows, already FDR-sorted) so the results screen can draw a
                # per-module dot/bar plot - not just the single top term. The
                # untrimmed block has one row per tested GO term (thousands,
                # mostly FDR=1.0); 50 is far more than any dot plot shows and
                # keeps the payload small even across many modules.
                module_enrichment[mod["module_id"]] = {
                    "namespace": e.get("namespace") or "",
                    "database": params["enrichment_database"],
                    "table": df_payload(e["block"]["table"].head(50)),
                }
            else:
                module_enrichment[mod["module_id"]] = None
            module_rows.append({
                "Module ID": mod["module_id"],
                "Size": mod["size"],
                "Genes": ", ".join(sorted(genes)),
                "Top enriched term": top_term,
                "Top enriched FDR": top_fdr,
                "Top enriched FDR (raw, uncorrected)": top_fdr_raw,
                "Enrichment database": enr_db,
            })
            module_id_map[mod["module_id"]] = mod["string_ids"]
        payload["modules"] = {
            "table": df_payload(pd.DataFrame(module_rows, columns=genomics_net.NETPROP_MODULE_COLUMNS)),
            "membership": module_id_map,       # module_id -> [string ids] (for node coloring)
            "enrichment": module_enrichment,   # module_id -> {namespace, database, table} | None
            "caveat": enr_result.get("caveat"),  # applies to every module (shared universe)
        }

    session["ora"]["netprop_results"] = payload
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _job_run_annotation(job_id: str, session: dict) -> dict:
    """Gene Function Annotation worker: per-gene descriptive GO + UniProt
    lookup (no statistics). Retrieval (final_df) is shared with ORA/Burden/
    Netprop; single-gene mode (see /api/annotation/single) also populates
    final_df with one row, so this worker is identical for both modes.
    """
    update_job(job_id, message="Starting gene function annotation...")
    final_df = session["ora"]["final_df"]
    if final_df is None or "UniProt ID" not in final_df.columns:
        raise ValueError("No confirmed UniProt ID list available for annotation.")

    config = session["ora"]["annotation_config"] or {}
    # Initial display preference only — the engine always fetches/returns FULL
    # evidence (see run_gene_function_annotation's docstring); the results
    # screen's "Include IEA" toggle re-filters client-side from the tidy table,
    # with no re-run, starting from this value.
    include_iea_default = bool(config.get("include_iea", True))

    id_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else (
        "Ensembl Gene ID" if "Ensembl Gene ID" in final_df.columns else None
    )
    if id_col is None:
        raise ValueError("No gene identifier column available for annotation.")

    uid_series = final_df["UniProt ID"].astype(str).str.strip()
    gid_series = final_df[id_col].astype(str).str.strip()
    pairs = [(gid, uid) for gid, uid in zip(gid_series.tolist(), uid_series.tolist()) if uid]
    if not pairs:
        raise ValueError("No genes with a UniProt accession to annotate.")

    def _progress(current, total, label):
        update_job(job_id, progress=current / max(total, 1), message=label)

    result = genomics_annot.run_gene_function_annotation(
        pairs, progress_callback=_progress, verbose=False,
    )
    payload = {
        "collapsed": df_payload(result["collapsed"]),
        "tidy": df_payload(result["tidy"]),
        "counts": result["counts"],
        "include_iea_default": include_iea_default,
    }
    session["ora"]["annotation_results"] = payload
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _network_param_bounds(config: dict, module: str = "transcriptomics") -> dict:
    """Validate + clamp the network config. Mirrors _netprop_param_bounds's
    idiom (a bad value falls back to the default rather than erroring), and
    doubles as the reproducibility record shown on the "Parameters used" pane -
    every value here is reported with the results and in the CSV downloads.

    Choice-set parameters are REJECTED to the default rather than clamped: a
    hostile n_permutations of 10**9 must become 500, never "the largest allowed
    value", which would still be a denial-of-service.

    module: "transcriptomics" (default) or "proteomics" - selects the
    identifier-space-specific bits (DEG cutoff defaults, the enrichment-
    database allow-list, whether MSigDB collections apply at all). Everything
    else - confidence/permutations/FDR/weighting/metric/module-size/cap - is
    identical by construction, since both tiles run the SAME
    Genomics_network engine underneath. See _job_run_network, its counterpart.
    """
    txnet = transcriptomics_network_engine
    pnet = proteomics_network_engine
    is_proteomics = module == "proteomics"
    engine = proteomics_engine if is_proteomics else transcriptomics_engine

    def _choice(value, choices, default):
        try:
            v = int(value)
        except (TypeError, ValueError):
            return default
        return v if v in choices else default

    def _float(value, lo, hi, default):
        try:
            v = float(value)
        except (TypeError, ValueError):
            return default
        return v if lo < v <= hi else default

    def _int_min(value, minimum, default):
        try:
            v = int(value)
        except (TypeError, ValueError):
            return default
        return v if v >= minimum else default

    direction = str(config.get("direction") or "all")
    if direction not in engine.DIRECTIONS:
        direction = "all"
    weighting = str(config.get("weighting") or txnet.NETWORK_SEED_WEIGHT_DEFAULT)
    if weighting not in txnet.NETWORK_SEED_WEIGHTINGS:
        weighting = txnet.NETWORK_SEED_WEIGHT_DEFAULT
    metric = str(config.get("metric") or transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT)
    if metric not in transcriptomics_engine.GSEA_RANKING_METRICS:
        metric = transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT

    if is_proteomics:
        db_choices = pnet.PPI_ENRICHMENT_DATABASES
        db_default = pnet.PPI_ENRICHMENT_DATABASE_DEFAULT
    else:
        db_choices = txnet.NETWORK_ENRICHMENT_DATABASES
        db_default = txnet.NETWORK_ENRICHMENT_DATABASE_DEFAULT
    # A run can label modules against several databases at once (each one's
    # namespaces/collections join ONE shared best-of-N-namespaces correction
    # pool - see enrich_modules_by_symbol/enrich_modules_by_accession's own
    # docstrings) - same dedup-and-whitelist idiom msigdb_collections below
    # already uses. Falls back to [db_default] only when the request named
    # NOTHING valid; an empty selection must never silently become "GO
    # anyway" at the CONFIGURE step without the caller knowing - the two
    # /network/configure routes reject an empty list outright instead.
    enrichment_databases = list(dict.fromkeys(
        d for d in (config.get("enrichment_databases") or []) if d in db_choices
    ))

    # MSigDB is Transcriptomics-only - build_string_id_index's identifier
    # bridge is symbol-space, which is what MSigDB's gene sets need; Proteomics
    # never offers it (see Proteomics_network.py's module docstring), so its
    # collections list is always empty regardless of what the client sends.
    msigdb_collections = []
    if not is_proteomics:
        msigdb_collections = list(dict.fromkeys(
            k for k in (config.get("msigdb_collections") or [])
            if k in transcriptomics_engine.MSIGDB_COLLECTIONS
        ))

    return {
        "confidence": _choice(config.get("confidence"), genomics_net.CONFIDENCE_CHOICES,
                              genomics_net.DEFAULT_CONFIDENCE),
        "n_permutations": _choice(config.get("n_permutations"), txnet.NETWORK_N_PERMUTATIONS_CHOICES,
                                  txnet.NETWORK_N_PERMUTATIONS_DEFAULT),
        "fdr_alpha": _float(config.get("fdr_alpha"), 0.0, 1.0, genomics_net.DEFAULT_FDR_ALPHA),
        "padj_cutoff": _float(config.get("padj_cutoff"), 0.0, 1.0, engine.PADJ_CUTOFF_DEFAULT),
        "lfc_cutoff": _float(config.get("lfc_cutoff"), -1e-9, 1e6, engine.LFC_CUTOFF_DEFAULT),
        "min_module_size": _int_min(config.get("min_module_size"), 2,
                                    genomics_net.DEFAULT_MIN_MODULE_SIZE),
        "max_hits_for_modules": _int_min(config.get("max_hits_for_modules"), 50,
                                         txnet.DEFAULT_MAX_HITS_FOR_MODULES),
        "direction": direction,
        "weighting": weighting,
        "metric": metric,
        "enrichment_databases": enrichment_databases,
        "msigdb_collections": msigdb_collections,
        "include_iea": bool(config.get("include_iea", True)),
        "module_detection": bool(config.get("module_detection", True)),
        # Fixed, reported for reproducibility (not user-editable).
        "restart_prob": genomics_net.DEFAULT_RESTART_PROB,
        "mcl_inflation": genomics_net.DEFAULT_MCL_INFLATION,
        "permutation_type": "degree-preserving",
        "string_version": genomics_net.STRING_VERSION,
        "seed_rng": 42,
    }


def _transcriptomics_network_param_bounds(config: dict) -> dict:
    return _network_param_bounds(config, "transcriptomics")


def _job_run_network(job_id: str, session: dict, module: str = "transcriptomics") -> dict:
    """Run Network / Module Detection for either tile. The two differ only in
    which orchestrator hooks apply (run_network_analysis vs.
    Proteomics_network.run_ppi_network_analysis, which itself just calls
    run_network_analysis with accession-space hooks - see that module) and
    which session sub-dict is read/written; the parameter bounding, the
    _graph strip, the edge-list build and the response payload shape are
    identical by construction, so a recruited protein / module / hit table
    looks the same in the client regardless of tile."""
    is_proteomics = module == "proteomics"
    state = session[module]
    resolved_df = state["resolved_df"]
    config = state["network_config"]
    if resolved_df is None:
        raise ValueError("No resolved data table available.")
    if not config:
        raise ValueError("Network analysis is not configured.")

    params = _network_param_bounds(config, module)
    taxon_id = str(state.get("taxon_id") or "")
    taxon_name = str(state.get("taxon_name") or "")

    def _progress(frac, label):
        update_job(job_id, progress=float(frac), message=str(label))

    if is_proteomics:
        result = proteomics_network_engine.run_ppi_network_analysis(
            resolved_df, taxon_id, taxon_name=taxon_name,
            confidence=params["confidence"],
            padj_cutoff=params["padj_cutoff"], lfc_cutoff=params["lfc_cutoff"],
            direction=params["direction"], weighting=params["weighting"], metric=params["metric"],
            restart_prob=params["restart_prob"], n_permutations=params["n_permutations"],
            fdr_alpha=params["fdr_alpha"], module_detection=params["module_detection"],
            min_module_size=params["min_module_size"], mcl_inflation=params["mcl_inflation"],
            max_hits_for_modules=params["max_hits_for_modules"],
            enrichment_databases=params["enrichment_databases"], include_iea=params["include_iea"],
            seed=params["seed_rng"], progress_callback=_progress,
        )
    else:
        result = transcriptomics_network_engine.run_network_analysis(
            resolved_df, taxon_id, taxon_name=taxon_name,
            confidence=params["confidence"],
            padj_cutoff=params["padj_cutoff"], lfc_cutoff=params["lfc_cutoff"],
            direction=params["direction"], weighting=params["weighting"], metric=params["metric"],
            restart_prob=params["restart_prob"], n_permutations=params["n_permutations"],
            fdr_alpha=params["fdr_alpha"], module_detection=params["module_detection"],
            min_module_size=params["min_module_size"], mcl_inflation=params["mcl_inflation"],
            max_hits_for_modules=params["max_hits_for_modules"],
            enrichment_databases=params["enrichment_databases"],
            msigdb_collections=params["msigdb_collections"], include_iea=params["include_iea"],
            seed=params["seed_rng"], progress_callback=_progress,
        )

    # The graph is a ~50-100 MB in-process object; it must never reach the
    # session or the JSON response. Pop it here, use it for the edge list, drop it.
    graph = result.pop("_graph", None)
    hit_ids = result.pop("_hit_ids", [])
    edges = []
    if graph is not None and hit_ids:
        # Rank by the SAME order the client figure draws (propagated score,
        # strongest first) - see _netprop_subnetwork_edges' docstring.
        ranked = result["propagation"]["table"].sort_values(
            by="Propagated score", ascending=False, kind="stable")["STRING ID"].tolist()
        edges = _netprop_subnetwork_edges(graph, hit_ids, rank_by=ranked)

    # Keep the raw DataFrames in the session for CSV download; serialize copies
    # for the client (same split the ORA/GSEA/regulon workers use).
    state["network_results"] = result

    payload = {
        "taxon_id": result["taxon_id"],
        "taxon_name": result["taxon_name"],
        "module_detection": result["module_detection"],
        "parameters": result["parameters"],
        "mapping": {**{k: v for k, v in result["mapping"].items() if k != "table"},
                    "table": df_payload(result["mapping"]["table"].head(500))},
        "seeds": result["seeds"],
        "propagation": {**{k: v for k, v in result["propagation"].items() if k != "table"},
                        "table": df_payload(result["propagation"]["table"].head(2000)),
                        "edges": edges},
        "modules": None,
        "caveats": result["caveats"],
    }
    if result["modules"]:
        mods = result["modules"]
        payload["modules"] = {
            "table": df_payload(mods["table"]),
            "membership": mods["membership"],
            "universe_size": mods["universe_size"],
            "caveat": mods["caveat"],
            "selection_caveat": mods["selection_caveat"],
            "cap_info": mods["cap_info"],
            # "namespace"/"table" describe the WINNER ("best"), same shape the
            # client already reads today; "candidates" additionally carries
            # every OTHER tested namespace/database's own (also capped)
            # table, so the results screen can offer a per-module database
            # selector rather than only ever showing whichever one won. Same
            # head(50) cap discipline as the winner - a handful of small
            # tables per module (at most one per requested database's
            # namespace/collection), not a new round trip per switch.
            "enrichment": {
                mid: ({
                    "namespace": e["best"]["namespace"],
                    "table": df_payload(e["best"]["block"]["table"].head(50)),
                    "candidates": {
                        key: ({"table": df_payload(c["block"]["table"].head(50))} if c else None)
                        for key, c in e["candidates"].items()
                    },
                } if e else None)
                for mid, e in mods["enrichment"].items()
            },
        }
    update_job(job_id, progress=1.0, message="Done")
    return payload


def _job_run_transcriptomics_network(job_id: str, session: dict) -> dict:
    """Run the Transcriptomics Network / Module Detection analysis."""
    return _job_run_network(job_id, session, "transcriptomics")


def _netprop_subnetwork_edges(graph: dict, string_ids: list, max_nodes: int = 400,
                              rank_by: "list | None" = None) -> list:
    """Edges (with confidence weights) among a set of STRING proteins, for the
    client network figure. Caps at max_nodes so an unusually large hit set never
    ships a giant edge list to the browser.

    rank_by: the caller's own preferred ordering of string_ids (strongest
    first). When the cap bites, the kept nodes are taken from THIS order.
    Passing None falls back to top-degree, the historical behaviour.

    Why this parameter exists: the client figure draws the top-N nodes by
    PROPAGATED SCORE (buildNetpropNetwork) and then keeps only nodes that have
    a drawn edge. Selecting the shipped subnetwork by DEGREE instead meant the
    two rankings could disagree - at genomics scale they mostly overlap, but a
    transcriptomics hit set is an order of magnitude larger, where the
    top-scoring nodes can be almost absent from a degree-selected 400, leaving
    the figure empty or, worse, showing a degree-biased subgraph labelled as
    "your top hits". Callers should pass the same ranking the figure uses.
    """
    import numpy as _np
    protein_index = graph["protein_index"]
    raw_adj = graph["raw_adjacency"]
    ids = [s for s in string_ids if s in protein_index]
    if len(ids) > max_nodes:
        if rank_by:
            order = {s: i for i, s in enumerate(rank_by)}
            # Unranked ids sort last, preserving a deterministic result.
            ids = sorted(ids, key=lambda s: order.get(s, len(order)))[:max_nodes]
        else:
            deg = _np.asarray(raw_adj.sum(axis=0)).ravel()
            ids = sorted(ids, key=lambda s: deg[protein_index[s]], reverse=True)[:max_nodes]
    idx_of = {s: protein_index[s] for s in ids}
    id_set = set(ids)
    sub = raw_adj[_np.ix_([idx_of[s] for s in ids], [idx_of[s] for s in ids])].tocoo()
    edges = []
    for r, c, w in zip(sub.row, sub.col, sub.data):
        if r < c:  # undirected: emit each pair once
            edges.append({"source": ids[r], "target": ids[c], "weight": round(float(w), 3)})
    return edges


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------

app = FastAPI(title="Omiscera API")

STATIC_DIR = os.path.join(BASE_DIR, "static")
if os.path.isdir(STATIC_DIR):
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/")
def serve_index():
    """Serve index.html with a cache-busting ?v=<mtime> on the app.js/style.css
    URLs. StaticFiles doesn't set an explicit Cache-Control, so browsers can
    heuristically cache those two assets indefinitely (surviving even a full
    restart); templating in each file's own mtime forces a brand-new URL the
    instant it changes on disk, so a stale copy can never be served."""
    with open(os.path.join(BASE_DIR, "index.html"), "r", encoding="utf-8") as fh:
        html = fh.read()
    v_js = int(os.path.getmtime(os.path.join(STATIC_DIR, "app.js")))
    v_css = int(os.path.getmtime(os.path.join(STATIC_DIR, "style.css")))
    html = html.replace('href="/static/style.css"', f'href="/static/style.css?v={v_css}"')
    html = html.replace('src="/static/app.js"', f'src="/static/app.js?v={v_js}"')
    return Response(content=html, media_type="text/html")


@app.get("/api/assets/{filename}")
def serve_asset(filename: str):
    safe_name = os.path.basename(filename)
    asset_path = os.path.join(BASE_DIR, "Assets", safe_name)
    if not os.path.isfile(asset_path):
        raise HTTPException(status_code=404, detail="Asset not found.")
    return FileResponse(asset_path, headers={"Cache-Control": "no-cache"})


INFO_DIR = os.path.join(BASE_DIR, "info")


@app.get("/api/info")
def list_info_docs():
    """Every module README under info/, for the public "Information" page: {filename, title} pairs, title = the
    file's own first line, in the order the files sort (groups Multi-omics' six together, etc.)."""
    if not os.path.isdir(INFO_DIR):
        return {"docs": []}
    docs = []
    for name in sorted(os.listdir(INFO_DIR)):
        if not name.startswith("README-") or not name.endswith(".txt"):
            continue
        path = os.path.join(INFO_DIR, name)
        with open(path, "r", encoding="utf-8") as fh:
            title = fh.readline().strip() or name
        docs.append({"filename": name, "title": title})
    return {"docs": docs}


@app.get("/api/info/{filename}")
def get_info_doc(filename: str):
    safe_name = os.path.basename(filename)
    if not safe_name.startswith("README-") or not safe_name.endswith(".txt"):
        raise HTTPException(status_code=404, detail="Unknown document.")
    path = os.path.join(INFO_DIR, safe_name)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Unknown document.")
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    return {"filename": safe_name, "text": text}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str):
    with _JOBS_LOCK:
        job = JOBS.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    return job


# ----- Shared / database catalog -----

@app.get("/api/databases")
def list_databases():
    options = [db for db in dc.get_supported_databases() if db != "HMDB"] + ["HMDB Import"]
    return {"databases": options}


@app.get("/api/catalog/validate")
def validate_catalog():
    try:
        contract = dc.validate_database_catalog_contract(allow_online=False)
        return contract
    except Exception as exc:
        return {"ok": False, "errors": [str(exc)], "warnings": [], "summary": {}}


@app.get("/api/contexts/search")
def search_contexts(database: str, query: str, limit: int = 10):
    results_df = dc.search_contexts(database, query, limit=limit, allow_online=True, verbose=False)
    if results_df.empty:
        return {"results": []}
    return {"results": results_df.to_dict(orient="records")}


# ----- Metabolic Connection: Step 1 - Upload -----

@app.post("/api/mc/upload")
async def mc_upload(file: UploadFile = File(...), x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        input_df, detected_delimiter = read_uploaded_csv_bytes(raw_bytes)
        deduplicated_df, duplicate_report_df, stats = deduplicate_uploaded_compounds(input_df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Unable to read CSV file: {exc}")

    session["input_df"] = deduplicated_df.copy()
    session["duplicate_report_df"] = duplicate_report_df.copy()
    session["mapped_df"] = None
    session["analysis_results"] = {}
    session["hmdb_import_by_rank"] = {}

    name_column = first_named_column(deduplicated_df, "Name")
    empty_name_count = int((name_column.astype(str).str.strip() == "").sum())

    return {
        "table": df_payload(deduplicated_df),
        "detected_delimiter": detected_delimiter,
        "empty_name_count": empty_name_count,
        "duplicate_name_count": stats["duplicate_name_count"],
        "duplicate_rows_removed": stats["duplicate_rows_removed"],
        "rows_before_duplicate_removal": stats["rows_before"],
        "rows_after_duplicate_removal": stats["rows_after"],
        "duplicate_report": df_payload(duplicate_report_df),
    }


# ----- Metabolic Connection: Step 2 - Mapping -----

@app.post("/api/mc/map")
def mc_map(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["input_df"] is None:
        raise HTTPException(status_code=400, detail="Please upload a CSV file before running mapping.")
    if "Name" not in session["input_df"].columns:
        raise HTTPException(status_code=400, detail="The uploaded file must contain a 'Name' column.")
    job_id = submit_job(_job_map_compounds, session)
    return {"job_id": job_id}


@app.get("/api/mc/download/mapped.csv")
def mc_download_mapped(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["mapped_df"] is None:
        raise HTTPException(status_code=400, detail="No mapped table available yet.")
    return df_to_csv_response(order_mapped_columns(session["mapped_df"]), "mapped_compounds.csv")


@app.get("/api/mc/download/result")
def mc_download_result(rank: str, table: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session.get("analysis_results") or {}
    result = results.get(rank)
    if result is None:
        raise HTTPException(status_code=404, detail="No results available for this rank.")
    if table not in ("table1", "table2", "table3"):
        raise HTTPException(status_code=400, detail="Invalid table key.")
    df = result.get(table)
    if df is None or df.empty:
        raise HTTPException(status_code=400, detail="This table has no data to download.")
    safe_rank = rank.replace(" ", "_")
    suffix = "compound_hits" if table == "table1" else "enzyme_presence"
    if table == "table3":
        suffix = "kegg_enzyme_presence"
    return df_to_csv_response(df, f"{safe_rank}_{suffix}.csv")


# ----- Metabolic Connection: Step 3 - Rank / database / context selection -----

@app.post("/api/mc/hmdb-import")
async def mc_hmdb_import(
    file: UploadFile = File(...),
    rank: str = Form(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        if file.filename and file.filename.endswith(".xlsx"):
            hmdb_df = pd.read_excel(io.BytesIO(raw_bytes), dtype=str)
        else:
            text = raw_bytes.decode("utf-8-sig", errors="replace")
            hmdb_df = pd.read_csv(io.StringIO(text), dtype=str, sep=None, engine="python")
        hmdb_df = standardize_dataframe_columns(hmdb_df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to read HMDB file: {exc}")

    if not any(c in hmdb_df.columns for c in ("HMDB ID", "Name")):
        raise HTTPException(
            status_code=400,
            detail="File does not contain an HMDB ID or Name column — compounds cannot be matched.",
        )

    key_cols = [c for c in ("HMDB ID", "Name", "InChIKey") if c in hmdb_df.columns]
    session["hmdb_import_by_rank"][rank] = hmdb_df

    return {
        "filename": file.filename,
        "row_count": len(hmdb_df),
        "key_columns": key_cols,
        "preview": df_payload(hmdb_df),
    }


@app.post("/api/mc/analyze")
def mc_analyze(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["mapped_df"] is None:
        raise HTTPException(status_code=400, detail="Run compound mapping before analysis.")

    selections = payload.get("selections", [])
    if not selections:
        raise HTTPException(status_code=400, detail="Select at least one database/context before running the analysis.")

    rank_settings: dict = {}
    for ctx in selections:
        rname = f"Rank {ctx['rank']}"
        rank_settings.setdefault(rname, {})[ctx["database"]] = {
            "context_id": ctx.get("context_id", ""),
            "context_name": ctx.get("context_name", ""),
            "context_type": ctx.get("context_type", ""),
            "status": "selected",
        }

    job_id = submit_job(_job_run_analysis, session, rank_settings)
    return {"job_id": job_id}


# ----- Annotation Analytics: Genomics ORA -----

@app.post("/api/ora/genes/upload")
async def ora_genes_upload(
    file: UploadFile = File(...),
    identifier_type: str = Form(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        gene_df, detected_type = gene_resolver.load_gene_table_from_bytes(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    if detected_type != identifier_type:
        expected_col = "Gene Symbol" if identifier_type == "gene_symbol" else "Ensembl Gene ID"
        raise HTTPException(
            status_code=400,
            detail=(
                f"Uploaded file does not contain a '{expected_col}' column. "
                f"Detected: {detected_type.replace('_', ' ').title()}."
            ),
        )

    gene_df, dup_stats = gene_resolver.deduplicate_gene_table(gene_df, identifier_type)
    session["ora"] = _fresh_ora_state()
    session["ora"]["identifier_type"] = identifier_type
    session["ora"]["gene_df"] = gene_df

    return {
        "table": df_payload(gene_df),
        "duplicates_removed": dup_stats["duplicates_removed"],
    }


@app.post("/api/ora/genes/score-column")
def ora_genes_score_column(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    """Register which already-uploaded column is the Gene-set Burden Test's
    variant score, without re-uploading the file. The gene table is already
    stored in session["ora"]["gene_df"] (extra columns preserved by
    _parse_gene_csv_text) from the shared /api/ora/genes/upload step; the
    resolver later drops extra columns when it rebuilds rows, so the score is
    captured here - keyed by the original (lowercased) gene id - and re-joined
    onto the resolved table at burden-run time (see _job_run_burden).
    """
    session = session_dep(x_session_id)
    gene_df = session["ora"]["gene_df"]
    identifier_type = session["ora"]["identifier_type"]
    if gene_df is None or identifier_type not in ("gene_symbol", "ensembl_gene_id"):
        raise HTTPException(status_code=400, detail="Upload a gene table first.")

    column = str(payload.get("column", "")).strip()
    direction = payload.get("direction", "higher")
    if not column:
        raise HTTPException(status_code=400, detail="No score column specified.")

    col_match = next((c for c in gene_df.columns if c.strip().lower() == column.lower()), None)
    if not col_match:
        raise HTTPException(status_code=400, detail=f"Column '{column}' not found in the uploaded file.")

    id_col = "Gene Symbol" if identifier_type == "gene_symbol" else "Ensembl Gene ID"
    variant_scores: dict = {}
    n_non_numeric = 0
    for _, row in gene_df.iterrows():
        gene_key = str(row.get(id_col, "")).strip().lower()
        if not gene_key:
            continue
        try:
            variant_scores[gene_key] = float(str(row[col_match]).strip())
        except (TypeError, ValueError):
            n_non_numeric += 1
    if not variant_scores:
        raise HTTPException(
            status_code=400,
            detail=f"Column '{column}' has no numeric values to use as the variant score.",
        )

    session["ora"]["variant_scores"] = variant_scores
    session["ora"]["score_direction"] = direction if direction in ("higher", "lower") else "higher"
    return {"n_scored": len(variant_scores), "n_non_numeric": n_non_numeric}


@app.post("/api/ora/seq/upload")
async def ora_seq_upload(file: UploadFile = File(...), x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        gene_df, _detected_type = gene_resolver.load_sequence_table_from_bytes(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    gene_df, dup_stats = gene_resolver.deduplicate_sequence_table(gene_df)
    session["ora"] = _fresh_ora_state()
    session["ora"]["identifier_type"] = "sequence"
    session["ora"]["gene_df"] = gene_df

    n_to_blast = (
        (gene_df["UniProt ID"].str.strip() == "").sum()
        if "UniProt ID" in gene_df.columns else len(gene_df)
    )

    return {
        "table": df_payload(gene_df),
        "duplicates_removed": dup_stats["duplicates_removed"],
        "n_to_blast": int(n_to_blast),
    }


@app.post("/api/ora/genes/resolve")
def ora_genes_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    identifier_type = session["ora"]["identifier_type"]
    if session["ora"]["gene_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a gene table first.")
    if identifier_type == "gene_symbol" and not payload.get("taxon_id"):
        raise HTTPException(status_code=400, detail="Select an organism before resolving gene symbols.")

    job_id = submit_job(_job_resolve_genes, session, identifier_type, payload.get("taxon_id", ""))
    return {"job_id": job_id}


@app.post("/api/ora/genes/resolve-kegg")
def ora_genes_resolve_kegg(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    """Optional second pass on top of /api/ora/genes/resolve: assigns each
    already-resolved gene a KEGG gene ID (a separate identifier from UniProt
    ID; only Genomics ORA/Burden's KEGG database matches on it). The KEGG
    organism is chosen independently from the UniProt taxon - KEGG organism
    codes and NCBI taxon IDs cannot be derived from one another (verified:
    KEGG's own organism hierarchy carries no taxon IDs).
    """
    session = session_dep(x_session_id)
    if session["ora"].get("resolved_df") is None:
        raise HTTPException(status_code=400, detail="Resolve gene identifiers first.")
    kegg_organism_code = str(payload.get("kegg_organism_code", "")).strip()
    if not kegg_organism_code:
        raise HTTPException(status_code=400, detail="Select a KEGG organism before resolving KEGG IDs.")

    job_id = submit_job(_job_resolve_kegg_ids, session, kegg_organism_code)
    return {"job_id": job_id}


@app.post("/api/ora/genes/discard-all")
def ora_genes_discard_all(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    table1_df = session["ora"]["table1_df"]
    if table1_df is None:
        raise HTTPException(status_code=400, detail="No resolved genes available.")
    final_df = table1_df.copy()
    session["ora"]["final_df"] = final_df
    return {"final_table": df_payload(final_df)}


@app.post("/api/ora/genes/apply-selections")
def ora_genes_apply_selections(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    table1_df = session["ora"]["table1_df"]
    display_cols = session["ora"]["display_cols"]
    id_col = session["ora"]["id_col"]
    if table1_df is None:
        raise HTTPException(status_code=400, detail="No resolved genes available.")

    choices = payload.get("choices", [])
    manual_rows = []
    for choice in choices:
        if choice.get("discard"):
            continue
        manual_rows.append({
            id_col: choice.get("gene_id", ""),
            "UniProt ID": choice.get("accession", ""),
            "Entry Name": choice.get("entry_name", ""),
            "Protein Name": choice.get("protein_name", ""),
            "NCBI Gene ID": "",
            "Resolution_status": "resolved_manual",
        })

    manual_df = pd.DataFrame(manual_rows, columns=display_cols) if manual_rows else pd.DataFrame(columns=display_cols)
    final_df = pd.concat([table1_df, manual_df], ignore_index=True)
    session["ora"]["final_df"] = final_df
    return {"final_table": df_payload(final_df)}


@app.get("/api/ora/download/genes.csv")
def ora_download_genes(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    final_df = session["ora"]["final_df"]
    if final_df is None:
        raise HTTPException(status_code=400, detail="No final gene table available.")
    return df_to_csv_response(final_df, "genes_resolved.csv")


@app.post("/api/ora/seq/blast")
def ora_seq_blast(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["gene_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a sequence table first.")
    job_id = submit_job(_job_run_blast, session)
    return {"job_id": job_id}


@app.post("/api/ora/seq/confirm")
def ora_seq_confirm(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    resolved_df = session["ora"]["resolved_df"]
    if resolved_df is None:
        raise HTTPException(status_code=400, detail="Run the BLAST search first.")
    selected_indices = payload.get("selected_indices", [])
    final_df = resolved_df.loc[resolved_df.index.isin(selected_indices)].reset_index(drop=True)
    session["ora"]["final_df"] = final_df
    return {"final_table": df_payload(final_df)}


@app.get("/api/ora/download/blast_all.csv")
def ora_download_blast_all(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    resolved_df = session["ora"]["resolved_df"]
    if resolved_df is None:
        raise HTTPException(status_code=400, detail="No BLAST results available.")
    return df_to_csv_response(resolved_df, "sequences_blast_resolved.csv")


@app.get("/api/ora/download/blast_confirmed.csv")
def ora_download_blast_confirmed(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    final_df = session["ora"]["final_df"]
    if final_df is None:
        raise HTTPException(status_code=400, detail="No confirmed selection available.")
    return df_to_csv_response(final_df, "sequences_blast_confirmed.csv")


# ----- Genomics ORA: statistical analysis (over-representation) -----
# Works identically regardless of which path (Gene Symbol / Ensembl / BLAST)
# produced session["ora"]["final_df"] - only the UniProt ID column is used.

CUSTOM_DB_LABEL = "Custom database"


@app.post("/api/ora/analysis/custom-db-upload")
async def ora_custom_db_upload(
    file: UploadFile = File(...),
    id_type: str = Form("accession"),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    id_type = str(id_type).strip().lower()
    id_type = id_type if id_type in ("symbol", "ensembl") else "accession"
    try:
        parsed = genomics_ora.load_custom_geneset_db(raw_bytes, id_type=id_type)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    parsed["filename"] = file.filename or "custom_database"
    session["ora"]["custom_db"] = parsed

    sample_terms = [
        {"term": tid, "name": info["name"], "n_genes": len(info["genes"])}
        for tid, info in list(parsed["terms"].items())[:5]
    ]
    return {
        "filename": parsed["filename"],
        "id_type": parsed["id_type"],
        "n_terms": parsed["n_terms"],
        "universe_size": parsed["universe_size"],
        "sample_terms": sample_terms,
    }


@app.post("/api/ora/analysis/background-upload")
async def ora_background_upload(
    file: UploadFile = File(...),
    database: str = Form(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """Store a user-supplied background (reference) gene list for one database.

    Background lives in UniProt-accession space, so it is only meaningful for
    the accession-space databases (GO, Reactome, KEGG). GWAS is symbol-native
    and the custom database carries its own universe, so those are rejected
    here. The list restricts the ORA universe to |annotated universe ∩ background|.
    """
    session = session_dep(x_session_id)
    db = str(database).strip()
    if db not in ("GO", "Reactome", "KEGG"):
        raise HTTPException(
            status_code=400,
            detail="A custom background applies only to GO, Reactome and KEGG (UniProt-accession space).",
        )
    raw_bytes = await file.read()
    try:
        accessions = genomics_ora.load_user_background_list(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    if not accessions:
        raise HTTPException(
            status_code=400,
            detail="No UniProt accessions found in the uploaded background file.",
        )
    session["ora"].setdefault("background_by_db", {})[db] = sorted(accessions)
    return {
        "database": db,
        "filename": file.filename or "background",
        "n_background": len(accessions),
    }


# The example ZIP is tailored to the ORA INPUT type (not the custom-DB id-type),
# because the input type is what determines the valid identifier space:
#   - gene_symbol   -> UniProt accession OR gene symbol  (include the symbol GMT)
#   - ensembl_gene_id -> UniProt accession OR Ensembl ID (accession data only)
#   - sequence (BLAST) -> UniProt accession only         (accession data only)
# Accession CSV/GMT are always the format template; the gene-symbol GMT is only
# relevant for gene-symbol input, so it is dropped for ensembl/homology (which
# otherwise made those examples look identical to the gene-symbol one).
_CUSTOM_ACCESSION_DATA = ["example_custom_geneset_db_TIDY.csv", "example_custom_geneset_db_GMT.gmt"]
CUSTOM_DB_EXAMPLE_BY_INPUT_TYPE = {
    "gene_symbol": ("README_custom_geneset_db.txt", _CUSTOM_ACCESSION_DATA + ["example_custom_geneset_db_SYMBOLS.gmt"]),
    "ensembl_gene_id": ("README_custom_geneset_db_ensembl.txt", _CUSTOM_ACCESSION_DATA),
    "sequence": ("README_custom_geneset_db_homology.txt", _CUSTOM_ACCESSION_DATA),
}


@app.get("/api/ora/background-example")
def ora_background_example():
    """Example background (reference gene list) file — plain UniProt accessions,
    the format load_user_background_list expects for the GO/Reactome background upload."""
    path = os.path.join(BASE_DIR, "Data2", "example_ora_background.csv")
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example background file not found.")
    return FileResponse(path, filename="example_ora_background.csv", media_type="text/csv")


BURDEN_EXAMPLE_BY_IDENTIFIER_TYPE = {
    "gene_symbol": "example_burden_input_symbol.csv",
    "ensembl_gene_id": "example_burden_input_ensembl.csv",
}


@app.get("/api/burden/example")
def burden_example(identifier_type: str = "gene_symbol"):
    """Example Gene-set Burden Test retrieval input (identifier column + a
    numeric Variant Score column), tailored to the identifier type selected
    on the retrieval screen - mirrors /api/ora/custom-db-example's per-type
    tailoring, but for the burden test's own upload step (not the custom
    gene-set database)."""
    fname = BURDEN_EXAMPLE_BY_IDENTIFIER_TYPE.get(
        str(identifier_type).strip().lower(), BURDEN_EXAMPLE_BY_IDENTIFIER_TYPE["gene_symbol"]
    )
    path = os.path.join(BASE_DIR, "Data2", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example burden-test input file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


@app.get("/api/ora/custom-db-example")
def ora_custom_db_example(input_type: str = "gene_symbol"):
    data2_dir = os.path.join(BASE_DIR, "Data2")
    readme_name, data_files = CUSTOM_DB_EXAMPLE_BY_INPUT_TYPE.get(
        str(input_type).strip().lower(), CUSTOM_DB_EXAMPLE_BY_INPUT_TYPE["gene_symbol"]
    )
    files_to_zip = [readme_name] + data_files
    missing = [f for f in files_to_zip if not os.path.isfile(os.path.join(data2_dir, f))]
    if missing:
        raise HTTPException(status_code=404, detail=f"Example file(s) not found: {', '.join(missing)}")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname in files_to_zip:
            zf.write(os.path.join(data2_dir, fname), arcname=fname)
    buf.seek(0)
    return Response(
        content=buf.read(),
        media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="custom_geneset_db_example.zip"'},
    )


@app.post("/api/ora/analysis/configure")
def ora_analysis_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["final_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No confirmed gene list available. Complete gene resolution first.",
        )

    databases = payload.get("databases", [])
    per_database = payload.get("per_database", {})
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    for db in databases:
        if db == CUSTOM_DB_LABEL:
            if not session["ora"].get("custom_db"):
                raise HTTPException(status_code=400, detail="Upload a gene-set file for the custom database.")
            continue
        if db not in per_database or not per_database[db].get("taxon_id"):
            raise HTTPException(status_code=400, detail=f"Select a species for {db}.")

    # Prune stale user-supplied backgrounds: keep only the databases whose
    # current payload still flags has_background. A "Start new run" without a
    # fresh upload therefore drops the previous background (falls back to the
    # whole annotated-genome universe).
    bg = session["ora"].get("background_by_db") or {}
    kept = {
        db: bg[db]
        for db in databases
        if db in bg and bool(per_database.get(db, {}).get("has_background"))
    }
    session["ora"]["background_by_db"] = kept

    session["ora"]["ora_config"] = {"databases": databases, "per_database": per_database}
    session["ora"]["ora_results"] = None
    return {"config": session["ora"]["ora_config"]}


@app.post("/api/ora/analysis/run")
def ora_analysis_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["ora_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_ora, session)
    return {"job_id": job_id}


@app.get("/api/ora/analysis/download/{database}.csv")
def ora_analysis_download(
    database: str,
    namespace: Optional[str] = None,
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    ora_results = session["ora"]["ora_results"]
    if not ora_results or database not in ora_results:
        raise HTTPException(status_code=400, detail="No ORA results available for this database.")
    key = namespace if namespace else database
    db_result = ora_results[database]
    if key not in db_result:
        raise HTTPException(status_code=400, detail=f"No result block for '{key}'.")
    table = db_result[key]["table"]
    filename = f"{database}_{namespace}_ora_results.csv" if namespace else f"{database}_ora_results.csv"
    return df_to_csv_response(table, filename)


# ----- Annotation Analytics: Genomics Gene-set Burden Test -----
# Retrieval (gene upload/resolve, organism selection) is fully shared with
# ORA via session["ora"]; only the analysis stage (configure/run/download)
# is separate, mirroring the /api/ora/analysis/* endpoints above.

@app.post("/api/burden/analysis/configure")
def burden_analysis_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["final_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No confirmed gene list available. Complete gene resolution first.",
        )
    if not session["ora"].get("variant_scores"):
        raise HTTPException(
            status_code=400,
            detail="No variant score column was captured. Re-upload your gene table with a score column selected.",
        )

    databases = payload.get("databases", [])
    per_database = payload.get("per_database", {})
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    for db in databases:
        if db == CUSTOM_DB_LABEL:
            if not session["ora"].get("custom_db"):
                raise HTTPException(status_code=400, detail="Upload a gene-set file for the custom database.")
            continue
        if db not in per_database or not per_database[db].get("taxon_id"):
            raise HTTPException(status_code=400, detail=f"Select a species for {db}.")

    session["ora"]["burden_config"] = {"databases": databases, "per_database": per_database}
    session["ora"]["burden_results"] = None
    return {"config": session["ora"]["burden_config"]}


@app.post("/api/burden/analysis/run")
def burden_analysis_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["burden_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_burden, session)
    return {"job_id": job_id}


@app.get("/api/burden/analysis/download/{database}.csv")
def burden_analysis_download(
    database: str,
    namespace: Optional[str] = None,
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    burden_results = session["ora"]["burden_results"]
    if not burden_results or database not in burden_results:
        raise HTTPException(status_code=400, detail="No burden test results available for this database.")
    key = namespace if namespace else database
    db_result = burden_results[database]
    if key not in db_result:
        raise HTTPException(status_code=400, detail=f"No result block for '{key}'.")
    table = db_result[key]["table"]
    filename = f"{database}_{namespace}_burden_results.csv" if namespace else f"{database}_burden_results.csv"
    return df_to_csv_response(table, filename)


# ----- Annotation Analytics: Genomics Network Propagation / Module Detection -----
# Retrieval (gene upload/resolve) is shared with ORA/Burden via session["ora"];
# the organism selector for THIS module is driven by STRING's own species list
# (STRING covers a subset of taxa), and the analysis stage is separate.

NETPROP_EXAMPLE_BY_IDENTIFIER_TYPE = {
    "gene_symbol": "example_netprop_input_symbol.csv",
    "ensembl_gene_id": "example_netprop_input_ensembl.csv",
}


@app.get("/api/netprop/species")
def netprop_species_search(query: str = "", x_session_id: Optional[str] = Header(default=None)):
    """Search STRING's covered-organism list (drives this module's organism
    selector - STRING covers a subset of the taxa ORA/Burden's UniProt search
    does, so picking from STRING guarantees a network exists)."""
    try:
        results = genomics_net.search_string_species(query, limit=15)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Could not load STRING species list: {exc}")
    return {"results": results}


@app.post("/api/netprop/analysis/configure")
def netprop_analysis_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["final_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No confirmed gene list available. Complete gene resolution first.",
        )
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism (STRING network) first.")
    try:
        if not genomics_net.is_taxon_in_string(taxon_id):
            raise HTTPException(status_code=400, detail=f"STRING has no network for taxon {taxon_id}.")
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Could not validate organism against STRING: {exc}")

    # network_propagation always runs; module_detection auto-implies it.
    module_detection = bool(payload.get("module_detection"))
    config = {
        "taxon_id": taxon_id,
        "taxon_name": str(payload.get("taxon_name") or ""),
        "network_propagation": True,
        "module_detection": module_detection,
        "confidence": payload.get("confidence"),
        "fdr_alpha": payload.get("fdr_alpha"),
        "min_module_size": payload.get("min_module_size"),
        "n_permutations": payload.get("n_permutations"),
    }
    session["ora"]["netprop_config"] = config
    session["ora"]["netprop_taxon_id"] = taxon_id
    session["ora"]["netprop_taxon_name"] = config["taxon_name"]
    session["ora"]["netprop_results"] = None
    return {"config": config}


@app.post("/api/netprop/analysis/run")
def netprop_analysis_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["netprop_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_netprop, session)
    return {"job_id": job_id}


@app.get("/api/netprop/analysis/download/{which}.csv")
def netprop_analysis_download(which: str, x_session_id: Optional[str] = Header(default=None)):
    """Download the propagation hit table ('hits') or the module table ('modules')."""
    session = session_dep(x_session_id)
    results = session["ora"]["netprop_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No network propagation results available.")
    if which == "hits":
        payload = results["propagation"]["table"]
        filename = "network_propagation_hits.csv"
    elif which == "modules":
        if not results.get("modules"):
            raise HTTPException(status_code=400, detail="No module results (module detection was not run).")
        payload = results["modules"]["table"]
        filename = "network_modules.csv"
    else:
        raise HTTPException(status_code=400, detail="Unknown table; use 'hits' or 'modules'.")
    df = pd.DataFrame(payload.get("rows", []), columns=payload.get("columns", []))
    return df_to_csv_response(df, filename)


@app.get("/api/netprop/example")
def netprop_example(identifier_type: str = "gene_symbol"):
    """Example network-propagation input (identifier column + an optional
    Gene Score column), tailored to the identifier type on the upload screen."""
    fname = NETPROP_EXAMPLE_BY_IDENTIFIER_TYPE.get(
        str(identifier_type).strip().lower(), NETPROP_EXAMPLE_BY_IDENTIFIER_TYPE["gene_symbol"]
    )
    path = os.path.join(BASE_DIR, "Data2", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example network-propagation input file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


# ----- Annotation Analytics: Genomics Gene Function Annotation -----
# Retrieval (table mode) is shared with ORA/Burden/Netprop via session["ora"];
# single-gene mode (/api/annotation/single) also populates final_df with one
# row, so the SAME configure/run/download flow below serves both modes.

@app.post("/api/annotation/single")
def annotation_single(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    """Single-gene mode: auto-detect whether the typed identifier is a UniProt
    accession (resolution skipped entirely) or a Gene Symbol/Ensembl Gene ID
    (resolved as usual). Populates final_df with one row so the shared
    configure/run flow below is identical for table and single-gene mode."""
    session = session_dep(x_session_id)
    identifier = str(payload.get("identifier") or "").strip()
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not identifier:
        raise HTTPException(status_code=400, detail="Enter a gene identifier.")

    if genomics_annot.looks_like_uniprot_accession(identifier):
        uniprot_id = identifier.upper()
        if genomics_annot.uniprot_accession_exists(uniprot_id):
            session["ora"]["final_df"] = pd.DataFrame({"Gene Symbol": [identifier], "UniProt ID": [uniprot_id]})
            return {"uniprot_id": uniprot_id, "skipped_resolution": True, "status": "user_provided"}
        # Matches the accession pattern but UniProt has no such entry - most
        # likely a typo, or (rarely) a gene symbol that happens to match the
        # pattern. Fall through to symbol resolution instead of silently
        # producing a blank result for a nonexistent accession.
        if not taxon_id:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"'{identifier}' looks like a UniProt accession, but no such entry exists. "
                    "If it's a gene symbol, select an organism and try again."
                ),
            )

    if genomics_annot.looks_like_ensembl_gene_id(identifier):
        result = gene_resolver.resolve_ensembl_id(identifier, verbose=False)
    else:
        if not taxon_id:
            raise HTTPException(
                status_code=400,
                detail="Select an organism first (needed to resolve a Gene Symbol).",
            )
        result = gene_resolver.resolve_gene_symbol(identifier, taxon_id, verbose=False)

    status = result.get("Resolution_status", "not_found")
    uniprot_id = str(result.get("UniProtID") or "").strip()
    if status == "ambiguous" or not uniprot_id:
        raise HTTPException(
            status_code=400,
            detail=(
                f"'{identifier}' matched multiple UniProt entries; try entering its UniProt "
                "accession directly instead."
            ) if status == "ambiguous" else f"'{identifier}' could not be resolved to a UniProt accession.",
        )
    session["ora"]["final_df"] = pd.DataFrame({"Gene Symbol": [identifier], "UniProt ID": [uniprot_id]})
    return {"uniprot_id": uniprot_id, "skipped_resolution": False, "status": status}


@app.post("/api/annotation/analysis/configure")
def annotation_analysis_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["final_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No confirmed gene list available. Complete gene resolution first.",
        )
    config = {"include_iea": bool(payload.get("include_iea", True))}
    session["ora"]["annotation_config"] = config
    session["ora"]["annotation_results"] = None
    return {"config": config}


@app.post("/api/annotation/analysis/run")
def annotation_analysis_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["ora"]["annotation_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_annotation, session)
    return {"job_id": job_id}


@app.get("/api/annotation/analysis/download/{which}.csv")
def annotation_analysis_download(which: str, x_session_id: Optional[str] = Header(default=None)):
    """Download the collapsed one-row-per-gene table ('annotations') or the
    tidy one-row-per-gene-term table ('tidy')."""
    session = session_dep(x_session_id)
    results = session["ora"]["annotation_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No annotation results available.")
    if which == "annotations":
        payload = results["collapsed"]
        filename = "gene_function_annotations.csv"
    elif which == "tidy":
        payload = results["tidy"]
        filename = "gene_function_annotations_tidy.csv"
    else:
        raise HTTPException(status_code=400, detail="Unknown table; use 'annotations' or 'tidy'.")
    df = pd.DataFrame(payload.get("rows", []), columns=payload.get("columns", []))
    return df_to_csv_response(df, filename)


@app.get("/api/annotation/example")
def annotation_example(identifier_type: str = "gene_symbol"):
    """Example gene-list input, reusing the Burden Test's existing example
    files (same required Gene Symbol/Ensembl + UniProt ID columns; the extra
    Variant Score column those files carry is simply irrelevant/ignored here)."""
    fname = BURDEN_EXAMPLE_BY_IDENTIFIER_TYPE.get(
        str(identifier_type).strip().lower(), BURDEN_EXAMPLE_BY_IDENTIFIER_TYPE["gene_symbol"]
    )
    path = os.path.join(BASE_DIR, "Data2", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example gene-list input file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


# ---------------------------------------------------------------------------
# Transcriptomics — retrieval (DE-table upload, identifier resolution)
# ---------------------------------------------------------------------------

TRANSCRIPTOMICS_ID_TYPE_LABELS = {
    "gene_symbol": "Gene Symbol",
    "ensembl_gene": "Ensembl Gene ID",
    "ensembl_transcript": "Ensembl Transcript ID",
    "entrez": "Entrez Gene ID",
    "mixed": "a mix of identifier types",
}

# Each identifier type's "Example file" download. MUST stay in sync with
# static/app.js's TRANSCRIPTOMICS_EXAMPLE_FILENAMES, which names the file the
# browser saves it as - they had drifted (gene_symbol served the 3000-gene
# de_full file but saved it under the name of the unrelated 20-gene
# _symbol.csv), so a consistency test now asserts the two maps match.
#
# gene_symbol/ensembl_gene/entrez all point at ~3000-gene simulated
# transcriptomes (3 genuinely enriched Hallmark sets spiked in), so
# "Example file" -> Resolve -> Run actually produces results worth looking at
# in EVERY downstream module. The small 18-gene ensembl/entrez examples that
# used to be served here could not support TF/Regulon activity inference at
# all (default minimum regulon size 5 -> essentially no TF clears the bar ->
# empty activity table), which is a poor first run to hand a new user.
#
# ensembl_transcript is DELIBERATELY still the small file: it exists to
# demonstrate transcript COLLAPSING (a retriever behaviour), and enlarging it
# would break what it is for. The upload caption says so.
TRANSCRIPTOMICS_EXAMPLE_BY_IDENTIFIER_TYPE = {
    "gene_symbol": "example_transcriptomics_de_full.csv",
    "ensembl_gene": "example_transcriptomics_ensembl_full.csv",
    "ensembl_transcript": "example_transcriptomics_transcript.csv",
    "entrez": "example_transcriptomics_entrez_full.csv",
}


@app.post("/api/transcriptomics/upload")
async def transcriptomics_upload(
    file: UploadFile = File(...),
    identifier_type: str = Form(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        de_df = transcriptomics_retriever.load_de_table_from_bytes(raw_bytes)
        column_mapping = transcriptomics_retriever.detect_columns(de_df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Reject a clear contradiction between what the user selected and what the
    # identifier column actually looks like - e.g. selecting Entrez Gene ID but
    # uploading gene symbols would otherwise query mygene.info with the wrong
    # scope and silently return 100% unmapped, with no error anywhere. A
    # "mixed" detection is allowed through (a real, legitimate case the
    # retriever already handles), so only a confident, contradicting majority
    # blocks the upload. Mirrors /api/ora/genes/upload's identifier_type check.
    detected_type = transcriptomics_retriever.detect_identifier_type(de_df[column_mapping["identifier"]])
    if detected_type != "mixed" and detected_type != identifier_type:
        selected_label = TRANSCRIPTOMICS_ID_TYPE_LABELS.get(identifier_type, identifier_type)
        detected_label = TRANSCRIPTOMICS_ID_TYPE_LABELS.get(detected_type, detected_type)
        raise HTTPException(
            status_code=400,
            detail=(
                f"You selected {selected_label}, but the '{column_mapping['identifier']}' column "
                f"looks like {detected_label}. Choose the matching identifier type and re-upload."
            ),
        )

    session["transcriptomics"] = _fresh_transcriptomics_state()
    session["transcriptomics"]["identifier_type"] = identifier_type
    session["transcriptomics"]["de_df"] = de_df
    session["transcriptomics"]["column_mapping"] = column_mapping

    preview = de_df.head(200)
    return {
        "table": df_payload(preview),
        "column_mapping": column_mapping,
    }


@app.post("/api/transcriptomics/resolve")
def transcriptomics_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["de_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a differential-expression table first.")
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism before resolving identifiers.")
    # taxon_name is NOT written here - it is committed inside the job, beside
    # taxon_id, so a failed resolve leaves the session untouched rather than
    # half-updated. See _job_resolve_de's docstring.
    taxon_name = str(payload.get("taxon_name") or "").strip()

    identifier_type = session["transcriptomics"]["identifier_type"]
    job_id = submit_job(_job_resolve_de, session, identifier_type, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.get("/api/transcriptomics/download/resolved.csv")
def transcriptomics_download_resolved(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    resolved_df = session["transcriptomics"]["resolved_df"]
    if resolved_df is None:
        raise HTTPException(status_code=400, detail="No resolved table available.")
    return df_to_csv_response(resolved_df, "transcriptomics_resolved.csv")


# ---------------------------------------------------------------------------
# KEGG organism resolution. Shared by every module that derives KEGG's
# organism from retrieval's taxon_id/taxon_name rather than offering a
# separate KEGG organism control (Transcriptomics ORA+GSEA, Proteomics
# ORA+Quantitative Enrichment - see resolve_kegg_organism_code's docstring
# for why: a second, independently-settable picker is what caused the
# silent-divergence bug fixed in commit 27ce033). One pair of endpoints per
# module (different session key, same underlying resolver, since
# resolve_kegg_organism_code takes only taxon_id/taxon_name - it has no
# module-specific behaviour at all), backed by these two shared helpers so
# the note-construction/validation logic exists exactly once.
# ---------------------------------------------------------------------------

def _kegg_resolve_organism_for(session: dict, module: str) -> dict:
    """Read-only check: can KEGG's organism be derived automatically for this
    session's current organism (in the given module), and if not, what are
    the candidates?

    Touches no session state - a config screen calls this whenever KEGG is
    checked, purely to decide whether to show a confirmation picker. Returns
    resolve_kegg_organism_code's result as-is; a non-empty "candidates" list
    is the signal to offer a picker."""
    taxon_id = session[module].get("taxon_id", "")
    taxon_name = session[module].get("taxon_name", "")
    return transcriptomics_engine.resolve_kegg_organism_code(taxon_id, taxon_name, verbose=False)


def _kegg_confirm_organism_for(session: dict, module: str, payload: dict) -> dict:
    """Confirm one KEGG organism for the genuinely-ambiguous case (see
    resolve_kegg_organism_code) and persist it for this session's runs in the
    given module.

    Re-derives candidates for the session's CURRENT taxon_id/taxon_name and
    only accepts organism_code if it is really one of them - the same guard
    resolve_kegg_organism_code applies internally, applied again here so a
    stale or forged code from the client is rejected at the API boundary
    too, not just silently ignored three layers down inside a job."""
    taxon_id = session[module].get("taxon_id", "")
    taxon_name = session[module].get("taxon_name", "")
    organism_code = str(payload.get("organism_code") or "").strip().lower()
    if not organism_code:
        raise HTTPException(status_code=400, detail="organism_code is required.")

    result = transcriptomics_engine.resolve_kegg_organism_code(
        taxon_id, taxon_name, verbose=False, confirmed_organism_code=organism_code,
    )
    if result["organism_code"] != organism_code:
        # Not among this run's own candidates - resolve_kegg_organism_code
        # silently fell back to the unconfirmed result rather than raising,
        # so surface that here as a real error instead of appearing to
        # succeed with a different organism than the one requested.
        raise HTTPException(
            status_code=400,
            detail=(
                f"'{organism_code}' is not one of the KEGG organisms matched for "
                f"\"{taxon_name}\". Re-check /api/{module}/kegg/resolve-organism's "
                "candidates and confirm one of those exactly."
            ),
        )
    session[module]["kegg_organism_confirmation"] = {
        "taxon_id": taxon_id, "taxon_name": taxon_name,
        "organism_code": result["organism_code"], "organism_name": result["organism_name"],
    }
    return result


@app.post("/api/transcriptomics/kegg/resolve-organism")
def transcriptomics_kegg_resolve_organism(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_resolve_organism_for(session, "transcriptomics")


@app.post("/api/transcriptomics/kegg/confirm-organism")
def transcriptomics_kegg_confirm_organism(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_confirm_organism_for(session, "transcriptomics", payload)


@app.get("/api/transcriptomics/example")
def transcriptomics_example(identifier_type: str = "gene_symbol"):
    fname = TRANSCRIPTOMICS_EXAMPLE_BY_IDENTIFIER_TYPE.get(
        str(identifier_type).strip().lower(), TRANSCRIPTOMICS_EXAMPLE_BY_IDENTIFIER_TYPE["gene_symbol"]
    )
    path = os.path.join(BASE_DIR, "Data3", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example DE-table input file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


# ---------------------------------------------------------------------------
# Shared guard for every Transcriptomics analysis entry point
# ---------------------------------------------------------------------------

def _require_consistent_taxon(session: dict) -> None:
    """Refuse to configure an analysis when the session's CURRENT organism is
    not the one its resolved table was actually built from.

    Why this exists as a separate invariant rather than just fixing the paths
    that could drift: the organism is consumed by different databases through
    different fields - taxon_id drives GO/QuickGO, MSigDB and the BioMart
    gene-length correction (and keys the retriever's mygene cache), while
    taxon_name drives Reactome - so a drifted session does not fail, it
    silently returns ONE result set mixing two organisms. That is the worst
    possible failure mode for a user who cannot independently check the
    biology, so it is blocked centrally at every entry point instead of
    being chased path by path.

    meta["taxon"] is written by resolve_de_table itself
    (Transcriptomics_module_retriever.py, "taxon": str(taxon_id)) and is the
    only record of what the table was genuinely resolved under.

    No meta yet (nothing resolved) is NOT an inconsistency - the existing
    "resolved table required" checks in each endpoint handle that case, and
    duplicating them here would produce a worse error message.
    """
    tx = session.get("transcriptomics") or {}
    meta = tx.get("meta") or {}
    resolved_taxon = str(meta.get("taxon") or "").strip()
    current_taxon = str(tx.get("taxon_id") or "").strip()
    if not resolved_taxon or not current_taxon:
        return
    if resolved_taxon != current_taxon:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Organism mismatch: the resolved table was built for taxon {resolved_taxon}, "
                f"but the current organism is taxon {current_taxon}. Re-resolve your "
                f"identifiers for the current organism before running an analysis."
            ),
        )


# ---------------------------------------------------------------------------
# Transcriptomics — Over-Representation Analysis (ORA)
# ---------------------------------------------------------------------------

TRANSCRIPTOMICS_ORA_DATABASES = ["GO", "Reactome", "KEGG", "MSigDB", "Custom database"]


@app.post("/api/transcriptomics/ora/custom-db-upload")
async def transcriptomics_custom_db_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """Always gene-symbol space - unlike Genomics' custom database, there is
    no id_type choice here, since Transcriptomics is symbol-native throughout."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        parsed = transcriptomics_engine.load_custom_geneset_db(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    parsed["filename"] = file.filename or "custom_database"
    session["transcriptomics"]["custom_db"] = parsed

    sample_terms = [
        {"term": tid, "name": info["name"], "n_genes": len(info["genes"])}
        for tid, info in list(parsed["terms"].items())[:5]
    ]
    return {
        "filename": parsed["filename"],
        "n_terms": parsed["n_terms"],
        "universe_size": parsed["universe_size"],
        "sample_terms": sample_terms,
    }


@app.get("/api/transcriptomics/ora/custom-db-example")
def transcriptomics_custom_db_example():
    """Example custom gene-set database file, gene-symbol-keyed (GMT format)
    - the only identifier space Transcriptomics' custom database accepts."""
    path = os.path.join(BASE_DIR, "Data2", "example_custom_geneset_db_SYMBOLS.gmt")
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example custom database file not found.")
    return FileResponse(path, filename="example_custom_geneset_db_symbols.gmt", media_type="text/plain")


@app.post("/api/transcriptomics/ora/background-upload")
async def transcriptomics_background_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """Store a user-supplied background (reference) gene-symbol list, applied
    globally across every selected database (unlike Genomics, Transcriptomics
    has one shared universe, not a per-database one)."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        symbols = transcriptomics_engine.load_user_background_list(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    if not symbols:
        raise HTTPException(status_code=400, detail="No gene symbols could be parsed from this file.")
    session["transcriptomics"]["background_symbols"] = symbols
    return {"filename": file.filename or "background", "n_genes": len(symbols)}


@app.post("/api/transcriptomics/ora/background-remove")
def transcriptomics_background_remove(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    session["transcriptomics"]["background_symbols"] = None
    return {"ok": True}


@app.get("/api/transcriptomics/ora/background-example")
def transcriptomics_background_example():
    """Example background (reference gene list) file - plain gene symbols,
    the format load_user_background_list expects (mirrors Genomics'
    /api/ora/background-example, but symbol-space since Transcriptomics is
    symbol-native throughout)."""
    path = os.path.join(BASE_DIR, "Data3", "example_transcriptomics_background.csv")
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example background file not found.")
    return FileResponse(path, filename="example_transcriptomics_background.csv", media_type="text/csv")


@app.post("/api/transcriptomics/ora/configure")
def transcriptomics_ora_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved differential-expression table available. Complete retrieval first.",
        )
    _require_consistent_taxon(session)

    # dict.fromkeys dedupes while preserving order - a repeated database name
    # would otherwise double steps_total in run_ora while the "if db in
    # databases" branch it corresponds to still only runs once, stalling the
    # progress bar at 50%. Not reachable from the checkbox UI, only a direct
    # API caller.
    databases = list(dict.fromkeys(
        d for d in (payload.get("databases") or []) if d in TRANSCRIPTOMICS_ORA_DATABASES
    ))
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    if "Custom database" in databases and not session["transcriptomics"].get("custom_db"):
        raise HTTPException(status_code=400, detail="Upload a custom gene-set database file first.")

    # Same dedup-while-preserving-order idiom as `databases` above; whitelisted
    # against MSIGDB_COLLECTIONS' own keys (defense in depth - the frontend
    # only ever sends real keys, but a direct API caller could send anything).
    msigdb_collections = list(dict.fromkeys(
        c for c in (payload.get("msigdb_collections") or [])
        if c in transcriptomics_engine.MSIGDB_COLLECTIONS
    ))
    if "MSigDB" in databases and not msigdb_collections:
        raise HTTPException(status_code=400, detail="Select at least one MSigDB collection.")

    try:
        padj_cutoff = float(payload.get("padj_cutoff", transcriptomics_engine.PADJ_CUTOFF_DEFAULT))
        lfc_cutoff = float(payload.get("lfc_cutoff", transcriptomics_engine.LFC_CUTOFF_DEFAULT))
        min_term_size = int(payload.get("min_term_size", transcriptomics_engine.MIN_TERM_SIZE_DEFAULT))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")

    # Distinguish "key absent" (-> clusterProfiler default cap) from
    # "explicitly null/empty" (-> caller deliberately wants no cap), so a
    # direct API caller can still opt out of the default the way the config
    # screen's blank "no limit" field does.
    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING:
        max_term_size = transcriptomics_engine.MAX_TERM_SIZE_DEFAULT
    elif max_term_size_raw in (None, ""):
        max_term_size = None
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    config = {
        "databases": databases,
        "msigdb_collections": msigdb_collections,
        "padj_cutoff": padj_cutoff,
        "lfc_cutoff": lfc_cutoff,
        "min_term_size": min_term_size,
        "max_term_size": max_term_size,
        "use_length_correction": bool(payload.get("use_length_correction", True)),
        # Applies to GO and Reactome (both carry an evidence-code column;
        # MSigDB has no IEA concept at all). Default True keeps prior
        # behaviour unchanged for anyone not touching the new toggle.
        "include_iea": bool(payload.get("include_iea", True)),
        # Explicit flag, not "is something stored in the session": mirrors
        # Genomics' has_background pattern, so unchecking the box client-side
        # (without clicking "Remove background") correctly excludes a
        # previously-uploaded background from this run.
        "use_background": bool(payload.get("use_background", False)),
    }
    session["transcriptomics"]["ora_config"] = config
    session["transcriptomics"]["ora_results"] = None
    return {"config": config}


@app.post("/api/transcriptomics/ora/run")
def transcriptomics_ora_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["ora_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_transcriptomics_ora, session)
    return {"job_id": job_id}


@app.get("/api/transcriptomics/ora/download")
def transcriptomics_ora_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    """Download one result table's CSV. `key` is a result key like
    'GO: Biological Process (all)' — passed as a query param (not a path
    segment) since it contains characters (':', '(', ')', spaces) that don't
    survive as a clean URL path."""
    session = session_dep(x_session_id)
    results = session["transcriptomics"]["ora_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No ORA results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"transcriptomics_ora_{safe_name}.csv")


# ---------------------------------------------------------------------------
# Proteomics — upload / resolve / Over-Representation Analysis (ORA)
#
# Identifier space is UniProt accession ONLY (no identifier_type choice,
# unlike Transcriptomics) - see Proteomics_module_retriever.py's module
# docstring for why gene symbol was rejected as an input type.
# ---------------------------------------------------------------------------

PROTEOMICS_EXAMPLE_FILENAME = "example_proteomics_de_full.csv"
PROTEOMICS_EXAMPLE_MESSY_FILENAME = "example_proteomics_maxquant_messy.csv"


@app.post("/api/proteomics/upload")
async def proteomics_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        de_df = proteomics_retriever.load_de_table_from_bytes(raw_bytes)
        column_mapping = proteomics_retriever.detect_columns(de_df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # A brand new upload fully resets every downstream stage - same
    # discipline as transcriptomics_upload, so a stale resolved table /
    # ORA config from a PREVIOUS file can never silently leak into this run.
    session["proteomics"] = _fresh_proteomics_state()
    session["proteomics"]["de_df"] = de_df
    session["proteomics"]["column_mapping"] = column_mapping

    # 200 rows, matching transcriptomics_upload's own preview size exactly
    # (de_df.head(200)) - an earlier version of this endpoint used head(10),
    # an arbitrary, unjustified inconsistency with no reason to differ from
    # the other retriever's preview.
    return {
        "n_rows": len(de_df),
        "columns": list(de_df.columns),
        "column_mapping": column_mapping,
        "preview": df_payload(de_df.head(200)),
    }


@app.post("/api/proteomics/resolve")
def proteomics_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["de_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a protein-level differential-abundance table first.")
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism before resolving identifiers.")
    taxon_name = str(payload.get("taxon_name") or "").strip()

    job_id = submit_job(_job_resolve_proteins, session, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.get("/api/proteomics/download/resolved.csv")
def proteomics_download_resolved(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    resolved_df = session["proteomics"]["resolved_df"]
    if resolved_df is None:
        raise HTTPException(status_code=400, detail="No resolved table available.")
    return df_to_csv_response(resolved_df, "proteomics_resolved.csv")


# KEGG organism resolution for Proteomics (ORA + Quantitative Enrichment) -
# see the shared _kegg_resolve_organism_for/_kegg_confirm_organism_for
# helpers and their surrounding comment (Transcriptomics section above) for
# why this is a thin per-module wrapper, not separate logic.

@app.post("/api/proteomics/kegg/resolve-organism")
def proteomics_kegg_resolve_organism(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_resolve_organism_for(session, "proteomics")


@app.post("/api/proteomics/kegg/confirm-organism")
def proteomics_kegg_confirm_organism(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_confirm_organism_for(session, "proteomics", payload)


@app.get("/api/proteomics/example")
def proteomics_example(messy: bool = False):
    """messy=true serves the deliberately dirty MaxQuant-style example
    (protein groups, CON__/REV__ rows, '+' flag columns) - the only way to
    demonstrate the contaminant-filtering and protein-group-parsing
    behaviours that have no transcriptomics equivalent."""
    fname = PROTEOMICS_EXAMPLE_MESSY_FILENAME if messy else PROTEOMICS_EXAMPLE_FILENAME
    path = os.path.join(BASE_DIR, "Data3", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example protein-table input file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


def _require_consistent_proteomics_taxon(session: dict) -> None:
    """Proteomics analogue of _require_consistent_taxon - see that function's
    docstring for the full rationale (a drifted organism silently mixes two
    species in one result set rather than failing loudly)."""
    px = session.get("proteomics") or {}
    meta = px.get("meta") or {}
    resolved_taxon = str(meta.get("taxon") or "").strip()
    current_taxon = str(px.get("taxon_id") or "").strip()
    if not resolved_taxon or not current_taxon:
        return
    if resolved_taxon != current_taxon:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Organism mismatch: the resolved table was built for taxon {resolved_taxon}, "
                f"but the current organism is taxon {current_taxon}. Re-resolve your "
                f"accessions for the current organism before running an analysis."
            ),
        )


def _require_consistent_proteomics_ptm_taxon(session: dict) -> None:
    """PTM analogue of _require_consistent_proteomics_taxon, reading
    ptm_meta (the SEPARATE site-level resolve's own record) instead of meta
    (the protein-level pipeline's) - see _fresh_proteomics_state's comment
    for why the two resolves are independent."""
    px = session.get("proteomics") or {}
    ptm_meta = px.get("ptm_meta") or {}
    resolved_taxon = str(ptm_meta.get("taxon") or "").strip()
    current_taxon = str(px.get("taxon_id") or "").strip()
    if not resolved_taxon or not current_taxon:
        return
    if resolved_taxon != current_taxon:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Organism mismatch: the resolved phosphosite table was built for taxon "
                f"{resolved_taxon}, but the current organism is taxon {current_taxon}. "
                f"Re-resolve your sites for the current organism before running an analysis."
            ),
        )


@app.post("/api/proteomics/ora/custom-db-upload")
async def proteomics_custom_db_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """UniProt-accession space, matching Genomics' custom database
    (id_type="accession" is load_custom_geneset_db's own default)."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        parsed = genomics_ora.load_custom_geneset_db(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Guard against a silent wrong answer: Proteomics custom databases are
    # accession-only (no symbol/Ensembl option, unlike Genomics), but a
    # symbol-based file still PARSES cleanly here - it just then matches
    # nothing at run time, and the results screen blames the term-size
    # filter instead of the real cause. Reject up front with the real
    # reason if almost nothing in the file looks like a UniProt accession.
    all_genes = {g for info in parsed["terms"].values() for g in info["genes"]}
    if all_genes:
        n_accession_like = sum(1 for g in all_genes if proteomics_engine._UNIPROT_ACC_RE.match(g))
        if n_accession_like / len(all_genes) < 0.5:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"This file parsed, but only {n_accession_like} of {len(all_genes)} gene values look "
                    "like UniProt accessions (e.g. P04637) - Proteomics custom databases must use UniProt "
                    "accessions, not gene symbols or another identifier type. If your gene column has gene "
                    "symbols, convert them to UniProt accessions first (e.g. via UniProt's own ID mapping "
                    "tool at uniprot.org/id-mapping)."
                ),
            )

    parsed["filename"] = file.filename or "custom_database"
    session["proteomics"]["custom_db"] = parsed

    sample_terms = [
        {"term": tid, "name": info["name"], "n_genes": len(info["genes"])}
        for tid, info in list(parsed["terms"].items())[:5]
    ]
    return {
        "filename": parsed["filename"],
        "n_terms": parsed["n_terms"],
        "universe_size": parsed["universe_size"],
        "sample_terms": sample_terms,
    }


@app.get("/api/proteomics/ora/custom-db-example")
def proteomics_custom_db_example():
    """Example custom gene-set database, zipped with a README (mirrors
    /api/ora/custom-db-example's shape) - GMT and tidy CSV formats, built
    from the SAME six real Complex Portal complexes seeded into the bundled
    example_proteomics_de_full.csv DE example, so downloading this and
    running it against that example produces genuine enrichment."""
    data3_dir = os.path.join(BASE_DIR, "Data3")
    files_to_zip = [
        "README_proteomics_custom_geneset_db.txt",
        "example_proteomics_custom_geneset_db.gmt",
        "example_proteomics_custom_geneset_db_TIDY.csv",
    ]
    missing = [f for f in files_to_zip if not os.path.isfile(os.path.join(data3_dir, f))]
    if missing:
        raise HTTPException(status_code=404, detail=f"Example file(s) not found: {', '.join(missing)}")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fname in files_to_zip:
            zf.write(os.path.join(data3_dir, fname), arcname=fname)
    buf.seek(0)
    return Response(
        content=buf.read(),
        media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="proteomics_custom_geneset_db_example.zip"'},
    )


@app.post("/api/proteomics/ora/background-upload")
async def proteomics_background_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """UniProt-accession-space background list, applied globally across
    every selected database (mirrors Transcriptomics' one-shared-universe
    pattern, not Genomics' per-database one)."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        accessions = genomics_ora.load_user_background_list(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    if not accessions:
        raise HTTPException(status_code=400, detail="No UniProt accessions could be parsed from this file.")
    session["proteomics"]["background_accessions"] = accessions
    return {"filename": file.filename or "background", "n_accessions": len(accessions)}


@app.post("/api/proteomics/ora/background-remove")
def proteomics_background_remove(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    session["proteomics"]["background_accessions"] = None
    return {"ok": True}


@app.get("/api/proteomics/ora/background-example")
def proteomics_background_example():
    path = os.path.join(BASE_DIR, "Data3", "example_proteomics_background.csv")
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example background file not found.")
    return FileResponse(path, filename="example_proteomics_background.csv", media_type="text/csv")


@app.post("/api/proteomics/ora/configure")
def proteomics_ora_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved protein table available. Complete retrieval first.",
        )
    _require_consistent_proteomics_taxon(session)

    databases = list(dict.fromkeys(
        d for d in (payload.get("databases") or []) if d in proteomics_engine.PROTEOMICS_ORA_DATABASES
    ))
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    if "Custom database" in databases and not session["proteomics"].get("custom_db"):
        raise HTTPException(status_code=400, detail="Upload a custom gene-set database file first.")

    msigdb_collections = list(dict.fromkeys(
        c for c in (payload.get("msigdb_collections") or [])
        if c in transcriptomics_engine.MSIGDB_COLLECTIONS
    ))
    if "MSigDB" in databases and not msigdb_collections:
        raise HTTPException(status_code=400, detail="Select at least one MSigDB collection.")

    try:
        padj_cutoff = float(payload.get("padj_cutoff", proteomics_engine.PADJ_CUTOFF_DEFAULT))
        lfc_cutoff = float(payload.get("lfc_cutoff", proteomics_engine.LFC_CUTOFF_DEFAULT))
        min_term_size = int(payload.get("min_term_size", genomics_ora.MIN_TERM_SIZE))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")

    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING:
        max_term_size = genomics_ora.MAX_TERM_SIZE_DEFAULT
    elif max_term_size_raw in (None, ""):
        max_term_size = None
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    config = {
        "databases": databases,
        "msigdb_collections": msigdb_collections,
        "padj_cutoff": padj_cutoff,
        "lfc_cutoff": lfc_cutoff,
        "min_term_size": min_term_size,
        "max_term_size": max_term_size,
        "include_iea": bool(payload.get("include_iea", True)),
        "use_background": bool(payload.get("use_background", False)),
    }
    session["proteomics"]["ora_config"] = config
    session["proteomics"]["ora_results"] = None
    return {"config": config}


@app.post("/api/proteomics/ora/run")
def proteomics_ora_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["ora_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_proteomics_ora, session)
    return {"job_id": job_id}


@app.get("/api/proteomics/ora/download")
def proteomics_ora_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["proteomics"]["ora_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No ORA results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"proteomics_ora_{safe_name}.csv")


@app.post("/api/proteomics/quant/configure")
def proteomics_quant_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved protein table available. Complete retrieval first.",
        )
    _require_consistent_proteomics_taxon(session)

    databases = list(dict.fromkeys(
        d for d in (payload.get("databases") or []) if d in proteomics_engine.PROTEOMICS_ORA_DATABASES
    ))
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    if "Custom database" in databases and not session["proteomics"].get("custom_db"):
        raise HTTPException(status_code=400, detail="Upload a custom gene-set database file first.")

    msigdb_collections = list(dict.fromkeys(
        c for c in (payload.get("msigdb_collections") or [])
        if c in transcriptomics_engine.MSIGDB_COLLECTIONS
    ))
    if "MSigDB" in databases and not msigdb_collections:
        raise HTTPException(status_code=400, detail="Select at least one MSigDB collection.")

    metric = payload.get("metric", proteomics_engine.QUANT_METRIC_DEFAULT)
    if metric not in transcriptomics_engine.GSEA_RANKING_METRICS:
        raise HTTPException(status_code=400, detail=f"Unknown ranking metric: {metric!r}")

    try:
        min_term_size = int(payload.get("min_term_size", genomics_ora.MIN_TERM_SIZE))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")

    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING:
        max_term_size = genomics_ora.MAX_TERM_SIZE_DEFAULT
    elif max_term_size_raw in (None, ""):
        max_term_size = None
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    config = {
        "databases": databases,
        "msigdb_collections": msigdb_collections,
        "metric": metric,
        "min_term_size": min_term_size,
        "max_term_size": max_term_size,
        "include_iea": bool(payload.get("include_iea", True)),
        "use_background": bool(payload.get("use_background", False)),
    }
    session["proteomics"]["quant_config"] = config
    session["proteomics"]["quant_results"] = None
    return {"config": config}


@app.post("/api/proteomics/quant/run")
def proteomics_quant_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["quant_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_proteomics_quant, session)
    return {"job_id": job_id}


@app.get("/api/proteomics/quant/download")
def proteomics_quant_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["proteomics"]["quant_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No Quantitative Enrichment results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"proteomics_quant_{safe_name}.csv")


# ---------------------------------------------------------------------------
# Proteomics — PPI Network / Module Detection
# (shares the retriever + session keys resolved_df / taxon_id / taxon_name with
# ORA and Quantitative Enrichment; defines no retrieval endpoints of its own.
# See _network_param_bounds / _job_run_network above, and Proteomics_network.py's
# module docstring for why this tile can use Reactome when Transcriptomics'
# own Network tile cannot.)
# ---------------------------------------------------------------------------

@app.post("/api/proteomics/network/configure")
def proteomics_network_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved protein table available. Complete retrieval first.")
    _require_consistent_proteomics_taxon(session)

    taxon_id = str(session["proteomics"].get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism first.")
    # STRING covers ~12,500 species, far fewer than the UniProt taxonomy the
    # proteomics organism selector searches. Refuse fast and by name here
    # rather than after a long download.
    try:
        in_string = genomics_net.is_taxon_in_string(taxon_id)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Could not reach STRING to check organism coverage: {exc}")
    if not in_string:
        raise HTTPException(
            status_code=400,
            detail=(f"STRING has no protein-interaction network for taxon {taxon_id} "
                    f"({session['proteomics'].get('taxon_name') or 'this organism'}). "
                    f"PPI Network / Module Detection needs one. GO/Reactome/Complex Portal-based "
                    f"analyses (Over-Representation Analysis, Quantitative Enrichment) still work "
                    f"for this organism."))

    config = _network_param_bounds(payload or {}, "proteomics")
    if config["module_detection"] and not config["enrichment_databases"]:
        raise HTTPException(status_code=400, detail="Select at least one database to label modules with.")
    session["proteomics"]["network_config"] = config
    session["proteomics"]["network_results"] = None
    return {"config": config}


@app.post("/api/proteomics/network/run")
def proteomics_network_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["network_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_network, session, "proteomics")
    return {"job_id": job_id}


@app.get("/api/proteomics/network/download")
def proteomics_network_download(which: str = "hits", key: str = "", db: str = "",
                                x_session_id: Optional[str] = Header(default=None)):
    """CSV export. `which` selects the table; `key` names the module for
    which="enrichment". `db` optionally names which tested namespace/
    database's own candidate to export (the exact keys under that module's
    "candidates" in the run's JSON payload, e.g. "biological_process" or
    "Reactome") - omitted downloads the WINNING ("best") one, same as
    before multi-database selection existed. Mirrors
    transcriptomics_network_download exactly - see that function's
    docstring."""
    session = session_dep(x_session_id)
    results = session["proteomics"]["network_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No network results available.")

    if which == "hits":
        return df_to_csv_response(results["propagation"]["table"],
                                  "proteomics_network_hits.csv")
    if which == "mapping":
        return df_to_csv_response(results["mapping"]["table"],
                                  "proteomics_network_protein_mapping.csv")
    if which == "modules":
        if not results.get("modules"):
            raise HTTPException(status_code=400, detail="This run did not include module detection.")
        return df_to_csv_response(results["modules"]["table"],
                                  "proteomics_network_modules.csv")
    if which == "enrichment":
        if not results.get("modules"):
            raise HTTPException(status_code=400, detail="This run did not include module detection.")
        entry = results["modules"]["enrichment"].get(key)
        if entry is None:
            raise HTTPException(status_code=400, detail=f"No enrichment for module {key!r}.")
        if db:
            candidate = entry["candidates"].get(db)
            if candidate is None:
                raise HTTPException(status_code=400, detail=f"No {db!r} result for module {key!r}.")
            block = candidate["block"]
        else:
            block = entry["best"]["block"]
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
        return df_to_csv_response(block["table"],
                                  f"proteomics_network_enrichment_{safe_name}.csv")

    raise HTTPException(status_code=400, detail=f"Unknown table: {which!r}")


# ---------------------------------------------------------------------------
# Proteomics — PTM / Phosphosite Enrichment
# (a SEPARATE site-level upload/resolve path from ORA/Quant/Network above -
# see _fresh_proteomics_state's comment. Shares the organism selector
# (taxon_id/taxon_name) but nothing else.)
# ---------------------------------------------------------------------------

PROTEOMICS_PTM_EXAMPLE_FILENAME = "example_proteomics_phosphosites.csv"


@app.post("/api/proteomics/ptm/upload")
async def proteomics_ptm_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        sites_df = proteomics_ptm_retriever.load_site_table_from_bytes(raw_bytes)
        column_mapping = proteomics_ptm_retriever.detect_columns(sites_df)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Only the PTM-specific fields reset here (NOT session["proteomics"] as a
    # whole, unlike proteomics_upload's protein-level reset) - a new
    # phosphosite file must not silently discard a user's already-completed,
    # independent ORA/Quant/Network results, which use a completely
    # different session sub-state (resolved_df, not resolved_sites_df).
    session["proteomics"]["ptm_sites_df"] = sites_df
    session["proteomics"]["ptm_column_mapping"] = column_mapping
    session["proteomics"]["resolved_sites_df"] = None
    session["proteomics"]["ptm_meta"] = None
    session["proteomics"]["ptm_config"] = None
    session["proteomics"]["ptm_activity_results"] = None
    session["proteomics"]["ptm_ora_results"] = None

    return {
        "n_rows": len(sites_df),
        "columns": list(sites_df.columns),
        "column_mapping": column_mapping,
        "preview": df_payload(sites_df.head(200)),
    }


@app.post("/api/proteomics/ptm/resolve")
def proteomics_ptm_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["ptm_sites_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a phosphosite differential-abundance table first.")
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism before resolving identifiers.")
    taxon_name = str(payload.get("taxon_name") or "").strip()

    job_id = submit_job(_job_resolve_ptm_sites, session, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.get("/api/proteomics/ptm/download/resolved.csv")
def proteomics_ptm_download_resolved(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    resolved_sites_df = session["proteomics"]["resolved_sites_df"]
    if resolved_sites_df is None:
        raise HTTPException(status_code=400, detail="No resolved phosphosite table available.")
    return df_to_csv_response(resolved_sites_df, "proteomics_ptm_resolved.csv")


@app.get("/api/proteomics/ptm/example")
def proteomics_ptm_example():
    path = os.path.join(BASE_DIR, "Data3", PROTEOMICS_PTM_EXAMPLE_FILENAME)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, media_type="text/csv", filename=PROTEOMICS_PTM_EXAMPLE_FILENAME)


@app.post("/api/proteomics/ptm/configure")
def proteomics_ptm_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["resolved_sites_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved phosphosite table available. Complete retrieval first.")
    _require_consistent_proteomics_ptm_taxon(session)

    try:
        padj_cutoff = float(payload.get("padj_cutoff", proteomics_engine.PADJ_CUTOFF_DEFAULT))
        lfc_cutoff = float(payload.get("lfc_cutoff", proteomics_engine.LFC_CUTOFF_DEFAULT))
        min_term_size = int(payload.get("min_term_size", proteomics_ptm_engine.MIN_TERM_SIZE_DEFAULT))
        min_kinase_sites = int(payload.get("min_kinase_sites", proteomics_ptm_engine.MIN_KINASE_SITES_DEFAULT))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")

    metric = str(payload.get("metric") or proteomics_ptm_engine.KINASE_ACTIVITY_METRIC_DEFAULT)
    if metric not in transcriptomics_engine.GSEA_RANKING_METRICS:
        raise HTTPException(status_code=400, detail=f"Unknown ranking metric: {metric!r}")

    n_permutations = payload.get("n_permutations", proteomics_ptm_engine.PTM_N_PERMUTATIONS_DEFAULT)
    try:
        n_permutations = int(n_permutations)
    except (TypeError, ValueError):
        n_permutations = proteomics_ptm_engine.PTM_N_PERMUTATIONS_DEFAULT
    if n_permutations not in proteomics_ptm_engine.PTM_N_PERMUTATIONS_CHOICES:
        n_permutations = proteomics_ptm_engine.PTM_N_PERMUTATIONS_DEFAULT

    min_score_raw = payload.get("min_score")
    min_score = None
    if min_score_raw not in (None, ""):
        try:
            min_score = float(min_score_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Minimum confidence score must be numeric.")

    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING or max_term_size_raw in (None, ""):
        max_term_size = proteomics_ptm_engine.MAX_TERM_SIZE_DEFAULT
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    config = {
        "padj_cutoff": padj_cutoff, "lfc_cutoff": lfc_cutoff,
        "min_term_size": min_term_size, "max_term_size": max_term_size,
        "min_kinase_sites": min_kinase_sites, "metric": metric,
        "n_permutations": n_permutations, "seed": proteomics_ptm_engine.PTM_SEED_DEFAULT,
        "min_score": min_score,
    }
    session["proteomics"]["ptm_config"] = config
    session["proteomics"]["ptm_activity_results"] = None
    session["proteomics"]["ptm_ora_results"] = None
    return {"config": config}


@app.post("/api/proteomics/ptm/run")
def proteomics_ptm_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["proteomics"]["ptm_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_ptm, session)
    return {"job_id": job_id}


@app.get("/api/proteomics/ptm/download")
def proteomics_ptm_download(which: str = "activity", key: str = "",
                            x_session_id: Optional[str] = Header(default=None)):
    """CSV export. which='activity' downloads the kinase activity table;
    which='ora' downloads one direction's kinase-substrate-ORA table (key
    names it, e.g. 'Kinase substrates (up)')."""
    session = session_dep(x_session_id)
    if which == "activity":
        results = session["proteomics"]["ptm_activity_results"]
        if not results or results.get("table") is None:
            raise HTTPException(status_code=400, detail="No kinase activity results available.")
        return df_to_csv_response(results["table"], "proteomics_ptm_kinase_activity.csv")
    if which == "ora":
        results = session["proteomics"]["ptm_ora_results"]
        if not results:
            raise HTTPException(status_code=400, detail="No kinase substrate enrichment results available.")
        block = results["tables"].get(key)
        if block is None:
            raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
        return df_to_csv_response(block["table"], f"proteomics_ptm_ora_{safe_name}.csv")
    raise HTTPException(status_code=400, detail=f"Unknown table: {which!r}")


# ---------------------------------------------------------------------------
# Metabolomics — Over-Representation Analysis, Chemical Class Enrichment,
# Lipid Class / Structural Enrichment, Quantitative Enrichment.
#
# FOUR separate tiles (matching the app's "one tile per analysis" convention
# used by every other module), TWO upload paths:
#   - a plain compound-name LIST, shared by the first three tiles (they
#     differ only in which database(s) they test the list against);
#   - a CONCENTRATION TABLE for Quantitative Enrichment, whose statistics
#     Omiscera computes itself - see Metabolomics_module.compute_compound_
#     statistics. This is a deliberate correction: an earlier version of
#     this module wrongly asked for a user-supplied log2FC/pvalue/padj, an
#     RNA-seq/proteomics differential-expression convention that does not
#     match what a metabolomics researcher actually uploads (MetaboAnalyst
#     takes a compound list or a concentration table, never pre-computed
#     statistics for the latter).
# Identifier handling differs from every other omics module either way: the
# upload is compound NAMES, which a two-step resolve turns into a single
# internal compound key - see Metabolomics_module's docstring on identifier
# bridging, and _resolve_and_match_compound_names.
# ---------------------------------------------------------------------------

METABOLOMICS_LIST_EXAMPLE_FILENAME = "example_metabolomics_list.csv"
METABOLOMICS_QUANT_EXAMPLE_FILENAME = "example_metabolomics_concentrations.csv"
METABOLOMICS_PATHWAY_DATABASES = ["KEGG", "Reactome", "SMPDB", "PathBank"]


def _load_metabolomics_list_upload(raw_bytes: bytes) -> pd.DataFrame:
    """Parse an uploaded compound-NAME list (ORA / Chemical Class / Lipid
    Class). Required: Name. Optional: any identifier column the retriever/
    crosswalk understand. NO statistics columns - a list has no
    significance concept at all (every compound IS the test set; the
    universe is the database's own reference set - see
    Metabolomics_module.run_ora's own docstring). Column matching is
    case/space-insensitive, matching the other modules' mappers.
    """
    try:
        text = raw_bytes.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw_bytes.decode("latin1")
    separator = "\t" if ("\t" in text.split("\n")[0]) else ","
    df = pd.read_csv(io.StringIO(text), sep=separator, dtype=str).fillna("")

    canonical = {
        "name": "Name", "compound": "Name", "compoundname": "Name", "metabolite": "Name",
        "keggid": "KEGG ID", "kegg": "KEGG ID",
        "pubchemid": "PubChem ID", "pubchem": "PubChem ID", "cid": "PubChem ID",
        "inchikey": "InChIKey", "smiles": "SMILES", "formula": "Formula",
        "hmdbid": "HMDB ID", "hmdb": "HMDB ID",
        "chebiid": "ChEBI ID", "chebi": "ChEBI ID",
        "lipidmapsid": "LIPID MAPS ID", "lmid": "LIPID MAPS ID", "lipidmaps": "LIPID MAPS ID",
    }
    rename = {}
    for column in df.columns:
        key = re.sub(r"[\s_\-]+", "", str(column)).lower()
        if key in canonical:
            rename[column] = canonical[key]
    df = df.rename(columns=rename)

    if "Name" not in df.columns:
        raise ValueError(
            "Input table is missing the required 'Name' column. This tile takes a plain "
            "list of compound names - no statistics columns are needed."
        )
    df = df[df["Name"].astype(str).str.strip() != ""].reset_index(drop=True)
    if df.empty:
        raise ValueError("Input table has no rows with a compound Name.")
    return df


@app.post("/api/metabolomics/list/upload")
async def metabolomics_list_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        input_df = _load_metabolomics_list_upload(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # A new list upload resets the list path and every tile built on it, but
    # leaves the concentration path (Quantitative Enrichment) untouched -
    # they are independent uploads, same discipline as PTM's site-level
    # upload not resetting the protein-level path.
    mx = session["metabolomics"]
    mx["list_df"] = input_df
    id_columns = [c for c in ("KEGG ID", "PubChem ID", "InChIKey", "HMDB ID",
                              "ChEBI ID", "LIPID MAPS ID", "SMILES", "Formula")
                  if c in input_df.columns]
    mx["list_column_mapping"] = {"identifier_columns": id_columns}
    mx["list_resolved_df"] = None
    mx["list_matched_df"] = None
    mx["list_meta"] = None
    mx["ora_config"] = None
    mx["ora_results"] = None
    mx["chemclass_config"] = None
    mx["chemclass_results"] = None
    mx["lipidclass_config"] = None
    mx["lipidclass_results"] = None
    return {
        "n_rows": len(input_df),
        "columns": list(input_df.columns),
        "identifier_columns": id_columns,
        "preview": df_payload(input_df.head(200)),
    }


@app.get("/api/metabolomics/list/example")
def metabolomics_list_example():
    path = os.path.join(BASE_DIR, "Data3", METABOLOMICS_LIST_EXAMPLE_FILENAME)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, media_type="text/csv", filename=METABOLOMICS_LIST_EXAMPLE_FILENAME)


@app.post("/api/metabolomics/list/resolve")
def metabolomics_list_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["list_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a compound list first.")
    taxon_id = str(payload.get("taxon_id", "")).strip()
    taxon_name = str(payload.get("taxon_name", "")).strip()
    job_id = submit_job(_job_resolve_metabolites_list, session, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.get("/api/metabolomics/list/download/resolved.csv")
def metabolomics_list_download_resolved(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    matched_df = session["metabolomics"]["list_matched_df"]
    if matched_df is None:
        raise HTTPException(status_code=400, detail="No resolved compound list available.")
    return df_to_csv_response(matched_df, "metabolomics_list_resolved_compounds.csv")


@app.post("/api/metabolomics/list/background-upload")
async def metabolomics_list_background_upload(
    file: UploadFile = File(...),
    x_session_id: Optional[str] = Header(default=None),
):
    """Optional background compound list for the three list-input tiles:
    narrows a database's reference universe to "compounds that could
    plausibly have been detected" (e.g. an LC-MS panel) rather than every
    compound the database has ever annotated."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    try:
        bg_df = _load_metabolomics_list_upload(raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    job_id = submit_job(_job_resolve_background_compounds, session, bg_df)
    return {"job_id": job_id}


# KEGG organism for Metabolomics.
#
# By default it is DERIVED from the organism the user chose (the same shared
# resolver Transcriptomics and Proteomics use - resolve-organism below reports
# what that derivation gives). The user can instead pick a KEGG organism
# explicitly from KEGG's own list, exactly like the organism box in Metabolic
# Connection (/api/contexts/search?database=KEGG); the choice travels in the run
# config (kegg_organism_code) and is validated here. That replaces the
# confirm-a-candidate mechanism the gene modules use: there is no gene-ID
# mapping in Metabolomics that a separately-chosen organism could diverge from
# (the bug of commit 27ce033), because the organism only selects pathway maps.
@app.post("/api/metabolomics/kegg/resolve-organism")
def metabolomics_kegg_resolve_organism(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_resolve_organism_for(session, "metabolomics")


_KEGG_ORGANISM_CODE_RE = re.compile(r"[a-z][a-z0-9]{1,7}")


def _validate_kegg_organism(code) -> dict:
    """{organism_code, organism_name, n_maps} for a KEGG organism code the user
    picked, or HTTP 400 saying why it cannot be used. A code with no pathway maps
    is rejected rather than accepted and silently returning an empty KEGG table."""
    code = str(code or "").strip().lower()
    if not _KEGG_ORGANISM_CODE_RE.fullmatch(code):
        raise HTTPException(status_code=400, detail=f"{code!r} is not a KEGG organism code (e.g. hsa, mmu, sce).")
    index = metabolomics_engine.build_kegg_organism_pathway_index(code)
    if not index.get("names"):
        raise HTTPException(status_code=400, detail=index.get("caveat") or
                            f"KEGG lists no pathway maps for organism '{code}'.")
    name = ""
    try:
        found = dc.search_contexts("KEGG", code, limit=5, allow_online=True, verbose=False)
        if found is not None and not found.empty:
            exact = found[found["context_id"].astype(str).str.lower() == code]
            if not exact.empty:
                name = str(exact.iloc[0].get("context_name") or "")
    except Exception:
        name = ""
    return {"organism_code": code, "organism_name": name or code, "n_maps": len(index["names"])}


@app.post("/api/metabolomics/kegg/check-organism")
def metabolomics_kegg_check_organism(payload: dict):
    return _validate_kegg_organism(payload.get("organism_code"))


def _kegg_config_fields(payload: dict, databases: "list[str]") -> dict:
    """The explicit KEGG organism, validated, for a run whose databases include
    KEGG; otherwise empty (KEGG's organism is then derived, or unused)."""
    code = str(payload.get("kegg_organism_code") or "").strip()
    if not code or "KEGG" not in databases:
        return {"kegg_organism_code": "", "kegg_organism_name": ""}
    org = _validate_kegg_organism(code)
    return {"kegg_organism_code": org["organism_code"], "kegg_organism_name": org["organism_name"]}


@app.post("/api/metabolomics/list/background-remove")
def metabolomics_list_background_remove(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    session["metabolomics"]["background_compounds"] = None
    session["metabolomics"]["background_meta"] = None
    return {"ok": True}


def _metabolomics_max_term_size(payload: dict, default):
    """Blank/absent -> `default`; otherwise a positive int. The default is
    tile-specific: see _list_tile_config."""
    raw = payload.get("max_term_size")
    if raw in (None, "", 0):
        return default
    try:
        value = int(raw)
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Maximum term size must be a whole number.")
    return value if value > 0 else default


def _list_tile_config(payload: dict, databases: "list[str]") -> dict:
    """Shared parameter validation for the three list-input tiles.
    `databases` is FIXED per call site: the ORA route validates a
    user-chosen subset of the 4 pathway databases; Chemical Class / Lipid
    Class each test exactly their own one database, so there is nothing to
    choose there.

    max_term_size defaults to NO LIMIT for the two class tiles and to the
    shared 500 for pathway ORA. The 500 cap is a gene-set convention (a term
    holding a large share of the genome says little); a class taxonomy is
    different - the main-class level users actually want (e.g.
    Glycerophosphocholines, thousands of members in LIPID MAPS) sits far above
    500, so the cap silently emptied exactly the lipid result the tile exists
    for. Fisher's test is exact at any term size; the user can still set one."""
    is_class_tile = list(databases) in (["Chemical Class"], ["Lipid Class"])
    return {
        "databases": databases,
        "min_term_size": max(1, int(payload.get("min_term_size", genomics_ora.MIN_TERM_SIZE))),
        "max_term_size": _metabolomics_max_term_size(
            payload, None if is_class_tile else genomics_ora.MAX_TERM_SIZE_DEFAULT),
        "chemical_class_level": (payload.get("chemical_class_level")
                                 if payload.get("chemical_class_level") in metabolomics_engine.CLASS_LEVELS
                                 else "sub_class"),
        "lipid_class_level": (payload.get("lipid_class_level")
                              if payload.get("lipid_class_level") in metabolomics_engine.LIPID_CLASS_LEVELS
                              else "lm_main_class"),
        "pathbank_primary_only": bool(payload.get("pathbank_primary_only", True)),
        "use_background": bool(payload.get("use_background", False)),
        **_kegg_config_fields(payload, databases),
    }


@app.post("/api/metabolomics/ora/configure")
def metabolomics_ora_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["list_matched_df"] is None:
        raise HTTPException(status_code=400, detail="Resolve the compound list first.")
    requested = payload.get("databases") or METABOLOMICS_PATHWAY_DATABASES
    databases = [d for d in requested if d in METABOLOMICS_PATHWAY_DATABASES]
    if not databases:
        raise HTTPException(
            status_code=400,
            detail="Select at least one database from: " + ", ".join(METABOLOMICS_PATHWAY_DATABASES),
        )
    config = _list_tile_config(payload, databases)
    session["metabolomics"]["ora_config"] = config
    session["metabolomics"]["ora_results"] = None
    return {"config": config}


@app.post("/api/metabolomics/ora/run")
def metabolomics_ora_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["ora_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_metabolomics_ora, session)
    return {"job_id": job_id}


@app.get("/api/metabolomics/ora/download")
def metabolomics_ora_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["metabolomics"]["ora_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No ORA results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"metabolomics_ora_{safe_name}.csv")


@app.post("/api/metabolomics/chemclass/configure")
def metabolomics_chemclass_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["list_matched_df"] is None:
        raise HTTPException(status_code=400, detail="Resolve the compound list first.")
    config = _list_tile_config(payload, ["Chemical Class"])
    session["metabolomics"]["chemclass_config"] = config
    session["metabolomics"]["chemclass_results"] = None
    return {"config": config}


@app.post("/api/metabolomics/chemclass/run")
def metabolomics_chemclass_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["chemclass_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_metabolomics_chemclass, session)
    return {"job_id": job_id}


@app.get("/api/metabolomics/chemclass/download")
def metabolomics_chemclass_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["metabolomics"]["chemclass_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No Chemical Class results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"metabolomics_chemclass_{safe_name}.csv")


@app.post("/api/metabolomics/lipidclass/configure")
def metabolomics_lipidclass_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["list_matched_df"] is None:
        raise HTTPException(status_code=400, detail="Resolve the compound list first.")
    config = _list_tile_config(payload, ["Lipid Class"])
    session["metabolomics"]["lipidclass_config"] = config
    session["metabolomics"]["lipidclass_results"] = None
    return {"config": config}


@app.post("/api/metabolomics/lipidclass/run")
def metabolomics_lipidclass_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["lipidclass_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_metabolomics_lipidclass, session)
    return {"job_id": job_id}


@app.get("/api/metabolomics/lipidclass/download")
def metabolomics_lipidclass_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["metabolomics"]["lipidclass_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No Lipid Class results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"metabolomics_lipidclass_{safe_name}.csv")


# --- Quantitative Enrichment (concentration table) --------------------------

@app.post("/api/metabolomics/quant/upload")
async def metabolomics_quant_upload(
    file: UploadFile = File(...),
    orientation: str = "auto",
    group_column: Optional[str] = None,
    x_session_id: Optional[str] = Header(default=None),
):
    """orientation/group_column are the user's explicit overrides of the
    parser's auto-detection (the UI re-uploads the same file with them
    when the detected layout is wrong) - auto-detection on an ambiguous
    file is exactly where a silent wrong guess would be worst."""
    session = session_dep(x_session_id)
    raw_bytes = await file.read()
    if orientation not in ("auto", "standard", "transposed"):
        raise HTTPException(status_code=400, detail=f"Unknown orientation: {orientation!r}")
    try:
        values_df, groups, meta = metabolomics_engine.parse_concentration_table(
            raw_bytes, orientation=orientation, group_column=group_column or None)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    # Resets only the concentration path - the list path (if any) belongs
    # to different tiles and must not be disturbed by this upload.
    mx = session["metabolomics"]
    mx["conc_values_df"] = values_df
    mx["conc_groups"] = groups
    mx["conc_meta"] = meta
    mx["conc_resolved_df"] = None
    mx["conc_matched_df"] = None
    mx["quant_stats_df"] = None
    mx["quant_stats_meta"] = None
    mx["quant_config"] = None
    mx["quant_results"] = None
    mx["topology_stats_df"] = None
    mx["topology_stats_meta"] = None
    mx["topology_config"] = None
    mx["topology_results"] = None
    return {
        "n_samples": meta["n_samples"], "n_compounds": meta["n_compounds"],
        "orientation": meta["orientation"], "group_column": meta["group_column"],
        "group_levels": meta["group_levels"], "group_counts": meta["group_counts"],
        "candidate_group_columns": meta["candidate_group_columns"],
        "n_repeated_names": meta.get("n_repeated_names", 0),
        "n_empty_dropped": meta.get("n_empty_dropped", 0),
        "n_unnamed_dropped": meta.get("n_unnamed_dropped", 0),
        "preview": df_payload(values_df.reset_index().rename(columns={"index": "Sample"}).head(10)),
    }


@app.get("/api/metabolomics/quant/example")
def metabolomics_quant_example():
    path = os.path.join(BASE_DIR, "Data3", METABOLOMICS_QUANT_EXAMPLE_FILENAME)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, media_type="text/csv", filename=METABOLOMICS_QUANT_EXAMPLE_FILENAME)


@app.post("/api/metabolomics/quant/resolve")
def metabolomics_quant_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["conc_values_df"] is None:
        raise HTTPException(status_code=400, detail="Upload a concentration table first.")
    taxon_id = str(payload.get("taxon_id", "")).strip()
    taxon_name = str(payload.get("taxon_name", "")).strip()
    job_id = submit_job(_job_resolve_metabolites_conc, session, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.get("/api/metabolomics/quant/download/resolved.csv")
def metabolomics_quant_download_resolved(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    matched_df = session["metabolomics"]["conc_matched_df"]
    if matched_df is None:
        raise HTTPException(status_code=400, detail="No resolved compound table available.")
    return df_to_csv_response(matched_df, "metabolomics_quant_resolved_compounds.csv")


@app.post("/api/metabolomics/quant/configure")
def metabolomics_quant_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    """Computes the per-compound statistics (compute_compound_statistics)
    AND stores the enrichment-analysis config in one call: the two-group
    contrast is itself part of "how to score this run", not a separate
    step, and doing both here lets validation errors (too few samples in a
    group, an unknown group name) surface synchronously rather than only
    after a background job starts."""
    session = session_dep(x_session_id)
    mx = session["metabolomics"]
    if mx["conc_matched_df"] is None:
        raise HTTPException(status_code=400, detail="Resolve the concentration table first.")
    values_df, groups = mx["conc_values_df"], mx["conc_groups"]

    group_a = str(payload.get("group_a", "")).strip()
    group_b = str(payload.get("group_b", "")).strip()
    if not group_a or not group_b:
        raise HTTPException(status_code=400, detail="Select two groups to compare.")
    if group_a == group_b:
        raise HTTPException(status_code=400, detail="The two groups must be different.")
    log_transform = bool(payload.get("log_transform", True))
    min_per_group = max(2, int(payload.get("min_per_group", 3)))
    try:
        stats_df, stats_meta = metabolomics_engine.compute_compound_statistics(
            values_df, groups, group_a, group_b,
            log_transform=log_transform, min_per_group=min_per_group,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    if stats_df.empty:
        raise HTTPException(
            status_code=400,
            detail="No compound had enough samples in both groups to test. "
                   "Lower the minimum samples per group or check your group selection.",
        )

    requested = payload.get("databases") or metabolomics_engine.METABOLOMICS_ORA_DATABASES
    databases = [d for d in requested if d in metabolomics_engine.METABOLOMICS_ORA_DATABASES]
    if not databases:
        raise HTTPException(
            status_code=400,
            detail="Select at least one database from: "
                   + ", ".join(metabolomics_engine.METABOLOMICS_ORA_DATABASES),
        )
    metric = payload.get("metric", metabolomics_engine.QUANT_METRIC_DEFAULT)
    if metric not in metabolomics_engine.QUANT_METRICS:
        metric = metabolomics_engine.QUANT_METRIC_DEFAULT

    config = {
        "databases": databases,
        "metric": metric,
        "min_term_size": max(1, int(payload.get("min_term_size", genomics_ora.MIN_TERM_SIZE))),
        "max_term_size": _metabolomics_max_term_size(payload, genomics_ora.MAX_TERM_SIZE_DEFAULT),
        "chemical_class_level": (payload.get("chemical_class_level")
                                 if payload.get("chemical_class_level") in metabolomics_engine.CLASS_LEVELS
                                 else "sub_class"),
        "lipid_class_level": (payload.get("lipid_class_level")
                              if payload.get("lipid_class_level") in metabolomics_engine.LIPID_CLASS_LEVELS
                              else "lm_main_class"),
        "pathbank_primary_only": bool(payload.get("pathbank_primary_only", True)),
        "use_background": bool(payload.get("use_background", False)),
        **_kegg_config_fields(payload, databases),
    }
    mx["quant_stats_df"] = stats_df
    mx["quant_stats_meta"] = stats_meta
    mx["quant_config"] = config
    mx["quant_results"] = None
    return {
        "config": config,
        "stats_meta": stats_meta,
        "stats_preview": df_payload(stats_df.sort_values("pvalue").head(50)),
    }


@app.post("/api/metabolomics/quant/run")
def metabolomics_quant_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["quant_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_metabolomics_quant, session)
    return {"job_id": job_id}


@app.get("/api/metabolomics/quant/download")
def metabolomics_quant_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["metabolomics"]["quant_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No quantitative enrichment results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"metabolomics_quant_{safe_name}.csv")


@app.get("/api/metabolomics/quant/download/statistics.csv")
def metabolomics_quant_download_statistics(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    stats_df = session["metabolomics"]["quant_stats_df"]
    if stats_df is None:
        raise HTTPException(status_code=400, detail="No statistics available - configure the analysis first.")
    return df_to_csv_response(stats_df, "metabolomics_quant_statistics.csv")


# --- Pathway Topology Analysis (KEGG only; list OR concentration input) -----
#
# Unlike the four tiles above, this one has no upload/resolve route of its
# own: its module screen lets the user pick "list" or "concentration table"
# and then drives the EXISTING /list/upload+/list/resolve or /quant/upload+
# /quant/resolve routes unchanged (see _fresh_metabolomics_state's own
# comment). Only configure/run/download are new here.

@app.post("/api/metabolomics/topology/configure")
def metabolomics_topology_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = session["metabolomics"]
    input_mode = payload.get("input_mode")
    if input_mode not in ("list", "conc"):
        raise HTTPException(status_code=400, detail="input_mode must be 'list' or 'conc'.")

    if input_mode == "list":
        if mx["list_matched_df"] is None:
            raise HTTPException(status_code=400, detail="Resolve the compound list first.")
        metric = metabolomics_engine.QUANT_METRIC_DEFAULT
    else:
        if mx["conc_matched_df"] is None:
            raise HTTPException(status_code=400, detail="Resolve the concentration table first.")
        values_df, groups = mx["conc_values_df"], mx["conc_groups"]
        group_a = str(payload.get("group_a", "")).strip()
        group_b = str(payload.get("group_b", "")).strip()
        if not group_a or not group_b:
            raise HTTPException(status_code=400, detail="Select two groups to compare.")
        if group_a == group_b:
            raise HTTPException(status_code=400, detail="The two groups must be different.")
        log_transform = bool(payload.get("log_transform", True))
        min_per_group = max(2, int(payload.get("min_per_group", 3)))
        try:
            stats_df, stats_meta = metabolomics_engine.compute_compound_statistics(
                values_df, groups, group_a, group_b,
                log_transform=log_transform, min_per_group=min_per_group,
            )
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        if stats_df.empty:
            raise HTTPException(
                status_code=400,
                detail="No compound had enough samples in both groups to test. "
                       "Lower the minimum samples per group or check your group selection.",
            )
        metric = payload.get("metric", metabolomics_engine.QUANT_METRIC_DEFAULT)
        if metric not in metabolomics_engine.QUANT_METRICS:
            metric = metabolomics_engine.QUANT_METRIC_DEFAULT
        mx["topology_stats_df"] = stats_df
        mx["topology_stats_meta"] = stats_meta

    measure = payload.get("measure", metabolomics_engine.PATHWAY_TOPOLOGY_MEASURE_DEFAULT)
    if measure not in metabolomics_engine.PATHWAY_TOPOLOGY_MEASURES:
        measure = metabolomics_engine.PATHWAY_TOPOLOGY_MEASURE_DEFAULT

    config = {
        "input_mode": input_mode,
        "metric": metric,
        "measure": measure,
        "min_term_size": max(1, int(payload.get("min_term_size", genomics_ora.MIN_TERM_SIZE))),
        "max_term_size": _metabolomics_max_term_size(payload, genomics_ora.MAX_TERM_SIZE_DEFAULT),
        "use_background": bool(payload.get("use_background", False)),
        # Topology always tests KEGG, so its organism is always validated
        # when a code is given - unlike _kegg_config_fields' conditional
        # check for the multi-database tiles above.
        **_kegg_config_fields(payload, ["KEGG"]),
    }
    mx["topology_config"] = config
    mx["topology_results"] = None
    response = {"config": config}
    if input_mode == "conc":
        response["stats_meta"] = mx["topology_stats_meta"]
        response["stats_preview"] = df_payload(mx["topology_stats_df"].sort_values("pvalue").head(50))
    return response


@app.post("/api/metabolomics/topology/run")
def metabolomics_topology_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["metabolomics"]["topology_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_metabolomics_topology, session)
    return {"job_id": job_id}


@app.get("/api/metabolomics/topology/download")
def metabolomics_topology_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["metabolomics"]["topology_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No pathway topology results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"metabolomics_topology_{safe_name}.csv")


@app.get("/api/metabolomics/topology/download/statistics.csv")
def metabolomics_topology_download_statistics(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    stats_df = session["metabolomics"]["topology_stats_df"]
    if stats_df is None:
        raise HTTPException(status_code=400, detail="No statistics available - configure the analysis first.")
    return df_to_csv_response(stats_df, "metabolomics_topology_statistics.csv")


# ---------------------------------------------------------------------------
# Transcriptomics — Gene Set Enrichment Analysis (GSEA, pre-ranked)
#
# Retrieval (upload/resolve) and the custom-database/background-list widgets
# are FULLY SHARED with ORA above (same session keys: resolved_df, taxon_id/
# taxon_name, custom_db, background_symbols) - only configure/run/download
# are GSEA-specific, so this section reuses the ORA custom-db-upload/example/
# background-upload/remove routes above rather than duplicating them.
# ---------------------------------------------------------------------------

TRANSCRIPTOMICS_GSEA_DATABASES = ["GO", "Reactome", "KEGG", "MSigDB", "Custom database"]


@app.post("/api/transcriptomics/gsea/configure")
def transcriptomics_gsea_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved differential-expression table available. Complete retrieval first.",
        )
    _require_consistent_taxon(session)

    databases = list(dict.fromkeys(
        d for d in (payload.get("databases") or []) if d in TRANSCRIPTOMICS_GSEA_DATABASES
    ))
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")
    if "Custom database" in databases and not session["transcriptomics"].get("custom_db"):
        raise HTTPException(status_code=400, detail="Upload a custom gene-set database file first.")

    msigdb_collections = list(dict.fromkeys(
        c for c in (payload.get("msigdb_collections") or [])
        if c in transcriptomics_engine.MSIGDB_COLLECTIONS
    ))
    if "MSigDB" in databases and not msigdb_collections:
        raise HTTPException(status_code=400, detail="Select at least one MSigDB collection.")

    metric = str(payload.get("metric") or transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT)
    if metric not in transcriptomics_engine.GSEA_RANKING_METRICS:
        metric = transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT

    try:
        min_term_size = int(payload.get("min_term_size", transcriptomics_engine.GSEA_MIN_TERM_SIZE_DEFAULT))
        weight = float(payload.get("weight", transcriptomics_engine.GSEA_WEIGHT_DEFAULT))
        seed = int(payload.get("seed", transcriptomics_engine.GSEA_SEED_DEFAULT))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")
    # The config screen already constrains both (min="2"/min="0"), so this only
    # matters for a direct API call - but the failure mode without it is
    # silent, not loud: np.abs(metric)**weight on a gene whose metric is
    # exactly 0 (reachable - signed_log_pvalue is 0 whenever p=1.0) returns
    # inf for any negative weight, which then poisons ES/NES with NaN instead
    # of raising anywhere.
    min_term_size = max(2, min_term_size)
    weight = max(0.0, weight)

    try:
        n_permutations = int(payload.get("n_permutations", transcriptomics_engine.GSEA_N_PERMUTATIONS_DEFAULT))
    except (TypeError, ValueError):
        n_permutations = transcriptomics_engine.GSEA_N_PERMUTATIONS_DEFAULT
    if n_permutations not in transcriptomics_engine.GSEA_N_PERMUTATIONS_CHOICES:
        n_permutations = transcriptomics_engine.GSEA_N_PERMUTATIONS_DEFAULT

    # Same "key absent" vs "explicit null" distinction as ORA's max_term_size
    # (transcriptomics_ora_configure above) - a direct API caller can opt out
    # of the default the way the config screen's blank field does.
    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING:
        max_term_size = transcriptomics_engine.GSEA_MAX_TERM_SIZE_DEFAULT
    elif max_term_size_raw in (None, ""):
        max_term_size = None
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    config = {
        "databases": databases,
        "msigdb_collections": msigdb_collections,
        "metric": metric,
        "min_term_size": min_term_size,
        "max_term_size": max_term_size,
        "n_permutations": n_permutations,
        "weight": weight,
        "seed": seed,
        "include_iea": bool(payload.get("include_iea", True)),
        "use_background": bool(payload.get("use_background", False)),
    }
    session["transcriptomics"]["gsea_config"] = config
    session["transcriptomics"]["gsea_results"] = None
    return {"config": config}


@app.post("/api/transcriptomics/gsea/run")
def transcriptomics_gsea_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["gsea_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_transcriptomics_gsea, session)
    return {"job_id": job_id}


@app.get("/api/transcriptomics/gsea/download")
def transcriptomics_gsea_download(key: str, x_session_id: Optional[str] = Header(default=None)):
    """Download one result table's CSV. `key` is a result key like
    'GO: Biological Process' - passed as a query param (not a path segment)
    since it contains characters (':', spaces) that don't survive as a clean
    URL path."""
    session = session_dep(x_session_id)
    results = session["transcriptomics"]["gsea_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No GSEA results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"transcriptomics_gsea_{safe_name}.csv")


# ---------------------------------------------------------------------------
# Transcriptomics — Transcription Factor / Regulon
#
# Retrieval (upload/resolve) is FULLY SHARED with ORA/GSEA above (same
# session keys: resolved_df, taxon_id/taxon_name) - this module defines no
# custom-database/background widgets of its own; regulons ARE the gene-set
# source, so those ORA-section routes don't apply here.
# ---------------------------------------------------------------------------

TRANSCRIPTOMICS_REGULON_ANALYSES = ["enrichment", "activity"]


@app.post("/api/transcriptomics/regulon/configure")
def transcriptomics_regulon_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved differential-expression table available. Complete retrieval first.",
        )
    _require_consistent_taxon(session)

    databases = list(dict.fromkeys(
        d for d in (payload.get("databases") or []) if d in transcriptomics_regulon_engine.REGULON_DATABASES
    ))
    if not databases:
        raise HTTPException(status_code=400, detail="Select at least one database.")

    analyses = list(dict.fromkeys(
        a for a in (payload.get("analyses") or TRANSCRIPTOMICS_REGULON_ANALYSES)
        if a in TRANSCRIPTOMICS_REGULON_ANALYSES
    ))
    if not analyses:
        raise HTTPException(
            status_code=400,
            detail="Select at least one analysis (target enrichment and/or activity inference).",
        )

    dorothea_levels = str(payload.get("dorothea_levels") or transcriptomics_regulon_engine.DOROTHEA_CONFIDENCE_DEFAULT)
    if dorothea_levels not in transcriptomics_regulon_engine.DOROTHEA_CONFIDENCE_CHOICES:
        dorothea_levels = transcriptomics_regulon_engine.DOROTHEA_CONFIDENCE_DEFAULT

    chea3_libraries = [
        lib for lib in (payload.get("chea3_libraries") or transcriptomics_regulon_engine.CHEA3_LIBRARIES)
        if lib in transcriptomics_regulon_engine.CHEA3_LIBRARIES
    ] or list(transcriptomics_regulon_engine.CHEA3_LIBRARIES)

    metric = str(payload.get("metric") or transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT)
    if metric not in transcriptomics_engine.GSEA_RANKING_METRICS:
        metric = transcriptomics_engine.GSEA_RANKING_METRIC_DEFAULT

    try:
        padj_cutoff = float(payload.get("padj_cutoff", transcriptomics_engine.PADJ_CUTOFF_DEFAULT))
        lfc_cutoff = float(payload.get("lfc_cutoff", transcriptomics_engine.LFC_CUTOFF_DEFAULT))
        min_term_size = int(payload.get("min_term_size", transcriptomics_regulon_engine.MIN_TERM_SIZE_DEFAULT))
        min_regulon_size = int(payload.get("min_regulon_size", transcriptomics_regulon_engine.MIN_REGULON_SIZE_DEFAULT))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Thresholds must be numeric.")
    # Defensive clamps (same rationale as the ORA/GSEA configure endpoints
    # above): min_term_size <= 0 would make every regulon "large enough" to
    # test regardless of size, and a non-positive min_regulon_size would
    # leave only _ulm_score's own >=2-target guard standing between a bad
    # config and a degenerate activity fit.
    padj_cutoff = max(0.0, min(padj_cutoff, 1.0))
    lfc_cutoff = max(0.0, lfc_cutoff)
    min_term_size = max(1, min_term_size)
    min_regulon_size = max(2, min_regulon_size)

    # Same "key absent" vs "explicit null" distinction as ORA/GSEA's
    # max_term_size above - a direct API caller can opt out of the default.
    _MISSING = object()
    max_term_size_raw = payload.get("max_term_size", _MISSING)
    if max_term_size_raw is _MISSING:
        max_term_size = transcriptomics_regulon_engine.MAX_TERM_SIZE_DEFAULT
    elif max_term_size_raw in (None, ""):
        max_term_size = None
    else:
        try:
            max_term_size = int(max_term_size_raw)
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="Max term size must be numeric.")

    try:
        n_permutations = int(payload.get("n_permutations", transcriptomics_regulon_engine.REGULON_N_PERMUTATIONS_DEFAULT))
    except (TypeError, ValueError):
        n_permutations = transcriptomics_regulon_engine.REGULON_N_PERMUTATIONS_DEFAULT
    if n_permutations not in transcriptomics_regulon_engine.REGULON_N_PERMUTATIONS_CHOICES:
        n_permutations = transcriptomics_regulon_engine.REGULON_N_PERMUTATIONS_DEFAULT
    try:
        seed = int(payload.get("seed", transcriptomics_regulon_engine.REGULON_SEED_DEFAULT))
    except (TypeError, ValueError):
        seed = transcriptomics_regulon_engine.REGULON_SEED_DEFAULT

    config = {
        "databases": databases,
        "analyses": analyses,
        "dorothea_levels": dorothea_levels,
        "chea3_libraries": chea3_libraries,
        "metric": metric,
        "padj_cutoff": padj_cutoff,
        "lfc_cutoff": lfc_cutoff,
        "min_term_size": min_term_size,
        "max_term_size": max_term_size,
        "min_regulon_size": min_regulon_size,
        "n_permutations": n_permutations,
        "seed": seed,
        "use_background": bool(payload.get("use_background", False)),
    }
    session["transcriptomics"]["regulon_config"] = config
    session["transcriptomics"]["regulon_enrichment_results"] = None
    session["transcriptomics"]["regulon_activity_results"] = None
    return {"config": config}


@app.post("/api/transcriptomics/regulon/run")
def transcriptomics_regulon_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    config = session["transcriptomics"]["regulon_config"]
    if config is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    analyses = config["analyses"]
    if analyses == ["enrichment"]:
        job_id = submit_job(_job_run_transcriptomics_regulon_enrichment, session)
    elif analyses == ["activity"]:
        job_id = submit_job(_job_run_transcriptomics_regulon_activity, session)
    else:
        job_id = submit_job(_job_run_transcriptomics_regulon_both, session)
    return {"job_id": job_id}


@app.get("/api/transcriptomics/regulon/download")
def transcriptomics_regulon_download_enrichment(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["transcriptomics"]["regulon_enrichment_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No TF target enrichment results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"transcriptomics_regulon_enrichment_{safe_name}.csv")


# ---------------------------------------------------------------------------
# Transcriptomics - Network / Module Detection
# (shares the retriever + session keys resolved_df / taxon_id / taxon_name with
# ORA, GSEA and TF/Regulon; defines no retrieval endpoints of its own)
# ---------------------------------------------------------------------------

@app.post("/api/transcriptomics/network/configure")
def transcriptomics_network_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["resolved_df"] is None:
        raise HTTPException(
            status_code=400,
            detail="No resolved differential-expression table available. Complete retrieval first.")
    # 4th call site. taxon_id drives STRING/GO/MSigDB while taxon_name drives
    # Reactome, so a drifted session would silently mix two organisms.
    _require_consistent_taxon(session)

    taxon_id = str(session["transcriptomics"].get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism first.")
    # STRING covers ~12,500 species, far fewer than the UniProt taxonomy the
    # transcriptomics organism selector searches. Refuse fast and by name here
    # rather than after a long download.
    try:
        in_string = genomics_net.is_taxon_in_string(taxon_id)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Could not reach STRING to check organism coverage: {exc}")
    if not in_string:
        raise HTTPException(
            status_code=400,
            detail=(f"STRING has no protein-interaction network for taxon {taxon_id} "
                    f"({session['transcriptomics'].get('taxon_name') or 'this organism'}). "
                    f"Network / Module Detection needs one. GO- and Reactome-based analyses "
                    f"(Over-Representation Analysis, Gene Set Enrichment Analysis) still work for "
                    f"this organism."))

    config = _transcriptomics_network_param_bounds(payload or {})
    if config["module_detection"] and not config["enrichment_databases"]:
        raise HTTPException(status_code=400, detail="Select at least one database to label modules with.")
    if "MSigDB" in config["enrichment_databases"] and not config["msigdb_collections"]:
        raise HTTPException(
            status_code=400,
            detail="Select at least one MSigDB collection to label modules with, or choose GO instead.")

    session["transcriptomics"]["network_config"] = config
    session["transcriptomics"]["network_results"] = None
    return {"config": config}


@app.post("/api/transcriptomics/network/run")
def transcriptomics_network_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if session["transcriptomics"]["network_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_transcriptomics_network, session)
    return {"job_id": job_id}


@app.get("/api/transcriptomics/network/download")
def transcriptomics_network_download(which: str = "hits", key: str = "", db: str = "",
                                     x_session_id: Optional[str] = Header(default=None)):
    """CSV export. `which` selects the table; `key` names the module for
    which="enrichment". `db` optionally names which tested namespace/
    database's own candidate to export (the exact keys under that module's
    "candidates" in the run's JSON payload, e.g. "biological_process" or
    "MSigDB: Hallmark") - omitted downloads the WINNING ("best") one, same
    as before multi-database selection existed. Query params (module ids and
    table names are not URL-path-safe), matching the ORA/GSEA/regulon
    download convention. The FULL, uncapped tables are served here - the
    client only ever receives a truncated copy for rendering."""
    session = session_dep(x_session_id)
    results = session["transcriptomics"]["network_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No network results available.")

    if which == "hits":
        return df_to_csv_response(results["propagation"]["table"],
                                  "transcriptomics_network_hits.csv")
    if which == "mapping":
        return df_to_csv_response(results["mapping"]["table"],
                                  "transcriptomics_network_gene_mapping.csv")
    if which == "modules":
        if not results.get("modules"):
            raise HTTPException(status_code=400, detail="This run did not include module detection.")
        return df_to_csv_response(results["modules"]["table"],
                                  "transcriptomics_network_modules.csv")
    if which == "enrichment":
        if not results.get("modules"):
            raise HTTPException(status_code=400, detail="This run did not include module detection.")
        entry = results["modules"]["enrichment"].get(key)
        if entry is None:
            raise HTTPException(status_code=400, detail=f"No enrichment for module {key!r}.")
        if db:
            candidate = entry["candidates"].get(db)
            if candidate is None:
                raise HTTPException(status_code=400, detail=f"No {db!r} result for module {key!r}.")
            block = candidate["block"]
        else:
            block = entry["best"]["block"]
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
        return df_to_csv_response(block["table"],
                                  f"transcriptomics_network_enrichment_{safe_name}.csv")

    raise HTTPException(status_code=400, detail=f"Unknown table: {which!r}")


@app.get("/api/transcriptomics/regulon/download-activity")
def transcriptomics_regulon_download_activity(key: str, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = session["transcriptomics"]["regulon_activity_results"]
    if not results:
        raise HTTPException(status_code=400, detail="No TF activity results available.")
    block = results["tables"].get(key)
    if block is None:
        raise HTTPException(status_code=400, detail=f"Unknown result table: {key!r}")
    safe_name = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
    return df_to_csv_response(block["table"], f"transcriptomics_regulon_activity_{safe_name}.csv")



@app.post("/api/ora/chart-tiff")
def ora_chart_tiff(payload: dict):
    """Convert a client-rendered Plotly chart (PNG data URL) to a 300 DPI TIFF.

    Plotly.js can export SVG/PNG directly in the browser, but not TIFF - so the
    client renders a high-resolution PNG and this endpoint re-encodes it as TIFF
    with proper DPI metadata via Pillow.
    """
    if not _PIL_AVAILABLE:
        raise HTTPException(status_code=501, detail="Pillow is not installed on the server; TIFF export unavailable.")

    image_data_url = payload.get("image", "")
    filename = payload.get("filename") or "chart"
    if not filename.lower().endswith(".tiff"):
        filename = filename.rsplit(".", 1)[0] + ".tiff" if "." in filename else filename + ".tiff"

    if "," not in image_data_url:
        raise HTTPException(status_code=400, detail="Invalid image data URL.")
    try:
        png_bytes = base64.b64decode(image_data_url.split(",", 1)[1])
        img = _PILImage.open(io.BytesIO(png_bytes))
        img.load()
        if img.mode in ("RGBA", "P"):
            img = img.convert("RGB")
        out = io.BytesIO()
        img.save(out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
        out.seek(0)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to convert image to TIFF: {exc}")

    return Response(
        content=out.read(),
        media_type="image/tiff",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


def _df_to_image_response(fig, filename: str, fmt: str) -> Response:
    """Serialize a matplotlib figure to SVG/PNG/PDF bytes and return as Response."""
    buf = io.BytesIO()
    media_map = {"svg": "image/svg+xml", "png": "image/png", "pdf": "application/pdf"}
    media_type = media_map.get(fmt, "application/octet-stream")
    fig.savefig(buf, format=fmt, bbox_inches="tight", dpi=150)
    plt.close(fig)
    buf.seek(0)
    return Response(
        content=buf.read(),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


def _build_dot_plot_fig(df: pd.DataFrame, title: str):
    top = df[df["FDR"].apply(lambda x: pd.notna(x) and str(x) not in ("", "nan"))].head(20).copy()
    if top.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No significant terms", ha="center", va="center")
        ax.set_title(title)
        return fig
    top["gene_ratio"] = top["k"].astype(float) / top["n"].astype(float).clip(lower=1)
    top["neg_log_fdr"] = -top["FDR"].astype(float).apply(lambda x: math.log10(max(x, 1e-300)))
    norm = mcolors.Normalize(vmin=top["neg_log_fdr"].min(), vmax=top["neg_log_fdr"].max())
    cmap = cm.get_cmap("RdYlGn")
    fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.45)))
    sc = ax.scatter(
        top["gene_ratio"],
        range(len(top)),
        s=top["k"].astype(float) * 30 + 30,
        c=top["neg_log_fdr"].values,
        cmap=cmap, norm=norm, alpha=0.85,
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([n[:60] for n in top["Term Name"]], fontsize=8)
    ax.set_xlabel("Gene ratio (k/n)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    cb = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("-log10(FDR)", fontsize=9)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    fig.tight_layout()
    return fig


def _build_bar_plot_fig(df: pd.DataFrame, title: str):
    top = df[df["FDR"].apply(lambda x: pd.notna(x) and str(x) not in ("", "nan"))].head(20).copy()
    if top.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No significant terms", ha="center", va="center")
        ax.set_title(title)
        return fig
    top["neg_log_fdr"] = -top["FDR"].astype(float).apply(lambda x: math.log10(max(x, 1e-300)))
    fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.45)))
    ax.barh(range(len(top)), top["neg_log_fdr"].values, color="#1B639C", alpha=0.85)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([n[:60] for n in top["Term Name"]], fontsize=8)
    ax.set_xlabel("-log10(FDR)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    fig.tight_layout()
    return fig


def _build_volcano_fig(df: pd.DataFrame, title: str):
    rows = df[df["Fold Enrichment"].apply(
        lambda x: pd.notna(x) and str(x) not in ("", "nan", "inf")
    )].copy()
    if rows.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "No terms to plot", ha="center", va="center")
        ax.set_title(title)
        return fig
    rows["log2_fe"] = rows["Fold Enrichment"].astype(float).apply(lambda x: math.log2(max(x, 0.001)))
    rows["neg_log_p"] = rows["p-value"].astype(float).apply(lambda x: -math.log10(max(x, 1e-300)))
    colors = rows["FDR"].astype(float).apply(lambda x: "#c0392b" if x < 0.05 else "#95a5a6")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(rows["log2_fe"], rows["neg_log_p"], c=colors, s=30, alpha=0.7, linewidths=0)
    ax.axhline(-math.log10(0.05), color="#e74c3c", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("log2(Fold Enrichment)", fontsize=10)
    ax.set_ylabel("-log10(p-value)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(linestyle="--", alpha=0.3)
    from matplotlib.lines import Line2D
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=8, label="FDR < 0.05"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#95a5a6", markersize=8, label="FDR ≥ 0.05"),
        ],
        fontsize=8, loc="upper left",
    )
    fig.tight_layout()
    return fig


@app.get("/api/ora/analysis/chart/{kind}")
def ora_chart_export(
    kind: str,
    database: str,
    format: str = "svg",
    namespace: Optional[str] = None,
    x_session_id: Optional[str] = Header(default=None),
):
    if not _MPL_AVAILABLE:
        raise HTTPException(status_code=501, detail="Matplotlib not installed.")
    if format not in ("svg", "png", "pdf"):
        raise HTTPException(status_code=400, detail="format must be svg, png, or pdf.")
    if kind not in ("dot", "bar", "volcano"):
        raise HTTPException(status_code=400, detail="kind must be dot, bar, or volcano.")

    session = session_dep(x_session_id)
    ora_results = session["ora"]["ora_results"]
    if not ora_results or database not in ora_results:
        raise HTTPException(status_code=400, detail="No ORA results for this database.")
    key = namespace if namespace else database
    db_result = ora_results[database]
    if key not in db_result:
        raise HTTPException(status_code=400, detail=f"No result block for '{key}'.")

    table: pd.DataFrame = db_result[key]["table"]
    ns_label = ORA_GO_NAMESPACE_LABELS.get(namespace, namespace) if namespace else None
    title = f"{database}{' — ' + ns_label if ns_label else ''}"

    try:
        if kind == "dot":
            fig = _build_dot_plot_fig(table, title)
        elif kind == "bar":
            fig = _build_bar_plot_fig(table, title)
        else:
            fig = _build_volcano_fig(table, title)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Chart rendering failed: {e}")

    filename = f"{database}_{namespace or database}_{kind}.{format}"
    return _df_to_image_response(fig, filename, format)


# Namespace label map (mirrored from app.js for server-side chart titles)
ORA_GO_NAMESPACE_LABELS = {
    "biological_process": "Biological Process",
    "molecular_function": "Molecular Function",
    "cellular_component": "Cellular Component",
}


# ---------------------------------------------------------------------------
# Multi-omics - Joint Pathway Enrichment
#
# Its own session sub-dict (session["multiomics"]) and its own routes. It never
# reads another tile's session state: every omics is uploaded and resolved
# HERE, through that omics' own existing resolver (Multiomics_module_retriever;
# Metabolomics through the same _resolve_and_match_compound_names helper the
# Metabolomics tiles use). Nothing below edits an existing function - it only
# CALLS the shared helpers (_kegg_resolve_organism_for, _kegg_confirmed_code,
# _load_metabolomics_list_upload, ...).
# ---------------------------------------------------------------------------

MULTIOMICS_EXAMPLE_BY_OMICS = {
    "genomics": "example_multiomics_genes.csv",
    "transcriptomics": "example_transcriptomics_de_full.csv",
    "proteomics": "example_proteomics_de_full.csv",
    "metabolomics": "example_metabolomics_list.csv",
}


def _fresh_multiomics_input() -> dict:
    return {
        "df": None, "kind": None, "identifier_type": None, "column_mapping": None,
        "n_duplicates_removed": 0, "resolved_df": None, "matched_df": None, "meta": None,
    }


def _fresh_multiomics_state() -> dict:
    """Single source of truth for the 'multiomics' session sub-dict shape.

    One upload/resolve slot per omics type, entirely separate from
    session["ora"] / ["transcriptomics"] / ["proteomics"] / ["metabolomics"]:
    the user chose a dedicated upload here, so a re-run of another tile can
    never silently change what Multi-omics sees. ONE organism for the whole run
    (taxon_id/taxon_name) and ONE KEGG organism confirmation - the KEGG
    organism is resolved once per run, never once per omics.
    """
    return {
        "selected_omics": [],
        "taxon_id": "", "taxon_name": "",
        "kegg_organism_confirmation": None,
        "inputs": {o: _fresh_multiomics_input() for o in multiomics_retriever.OMICS_TYPES},
        "joint_config": None, "joint_results": None,
        "mapping": _fresh_multiomics_mapping(),
        "correlation": _fresh_multiomics_correlation(),
        "network": _fresh_multiomics_network(),
        "mofa": _fresh_multiomics_mofa(),
        "diablo": _fresh_multiomics_diablo(),
    }


def _fresh_multiomics_correlation() -> dict:
    """Correlation-based Integration's own slot: the matrices uploaded (samples x features), the gene symbols their
    features resolved to, the sample alignment, the config and the last result. Separate from Joint and Mapping."""
    return {"selected": [], "taxon_id": "", "taxon_name": "", "uploads": {}, "resolved": {}, "aligned": None,
            "config": None, "results": None}


def _fresh_multiomics_mapping() -> dict:
    """Cross-omics Pathway Mapping's own slot: its config (the DE cutoffs), the pickers
    already built per database, and the last mapped pathway (kept for the CSV download)."""
    return {"config": None, "pickers": {}, "results": None}


def _multiomics_clear_results(mx: dict) -> None:
    """Drop every analysis built on the current uploads / resolution (Joint and Mapping)."""
    mx["joint_config"] = None
    mx["joint_results"] = None
    mx["mapping"] = _fresh_multiomics_mapping()
    mx["network"] = _fresh_multiomics_network()


def _multiomics_state(session: dict) -> dict:
    return session["multiomics"]


@app.post("/api/multiomics/select")
def multiomics_select(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    """Choose which omics to combine (2-4). Re-choosing keeps the uploads of any
    omics that stays selected, and drops the joint config/results (they described
    a different combination)."""
    session = session_dep(x_session_id)
    requested = payload.get("omics_types") or []
    chosen = [o for o in multiomics_retriever.OMICS_TYPES if o in requested]
    if len(chosen) < 2:
        raise HTTPException(status_code=400, detail="Select at least two omics to combine.")
    previous = session["multiomics"]
    fresh = _fresh_multiomics_state()
    fresh["selected_omics"] = chosen
    for o in chosen:
        fresh["inputs"][o] = previous["inputs"][o]
    fresh["taxon_id"], fresh["taxon_name"] = previous["taxon_id"], previous["taxon_name"]
    fresh["correlation"] = previous["correlation"]   # a different method's state: choosing omics for Joint must not wipe it
    fresh["mofa"] = previous["mofa"]                 # MOFA reads the Correlation matrices, so it keeps its own result as well
    fresh["diablo"] = previous["diablo"]             # DIABLO reads them too, and keeps its own result the same way
    session["multiomics"] = fresh
    return {"selected": chosen}


@app.get("/api/multiomics/example")
def multiomics_example(omics_type: str = "genomics"):
    fname = MULTIOMICS_EXAMPLE_BY_OMICS.get(str(omics_type).strip().lower())
    if not fname:
        raise HTTPException(status_code=404, detail="Unknown omics type.")
    path = os.path.join(BASE_DIR, "Data3", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


@app.post("/api/multiomics/kegg/resolve-organism")
def multiomics_kegg_resolve_organism(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_resolve_organism_for(session, "multiomics")


@app.post("/api/multiomics/kegg/confirm-organism")
def multiomics_kegg_confirm_organism(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    return _kegg_confirm_organism_for(session, "multiomics", payload)


@app.post("/api/multiomics/resolve")
def multiomics_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    if not mx["selected_omics"]:
        raise HTTPException(status_code=400, detail="Choose the omics to combine first.")
    missing = [multiomics_engine.OMICS_LABELS[o] for o in mx["selected_omics"] if mx["inputs"][o]["df"] is None]
    if missing:
        raise HTTPException(status_code=400, detail="Upload a file for: " + ", ".join(missing) + ".")
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism before resolving identifiers.")
    taxon_name = str(payload.get("taxon_name") or "").strip()
    job_id = submit_job(_job_resolve_multiomics, session, taxon_id, taxon_name)
    return {"job_id": job_id}


@app.post("/api/multiomics/joint/configure")
def multiomics_joint_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    selected = mx["selected_omics"]
    if not selected or any(mx["inputs"][o]["resolved_df"] is None for o in selected):
        raise HTTPException(status_code=400, detail="Resolve every selected omics' identifiers first.")
    requested = payload.get("databases") or list(multiomics_engine.JOINT_DATABASES)
    databases = [d for d in multiomics_engine.JOINT_DATABASES if d in requested]
    if not databases:
        raise HTTPException(status_code=400, detail="Select KEGG and/or Reactome.")
    method = payload.get("method") or multiomics_engine.COMBINE_METHOD_DEFAULT
    if method not in multiomics_engine.COMBINE_METHODS:
        raise HTTPException(status_code=400, detail="Unknown combination method.")
    try:
        min_ts = max(1, int(payload.get("min_term_size", multiomics_engine.MIN_TERM_SIZE_DEFAULT)))
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail="Minimum term size must be a whole number.")
    max_ts = _metabolomics_max_term_size(payload, multiomics_engine.MAX_TERM_SIZE_DEFAULT)
    if max_ts and max_ts < min_ts:
        raise HTTPException(status_code=400, detail="Maximum term size must be at least the minimum.")
    cutoffs = _multiomics_parse_cutoffs(payload, selected)
    config = {"databases": databases, "method": method, "min_term_size": min_ts,
              "max_term_size": max_ts, "de_cutoffs": cutoffs}
    mx["joint_config"] = config
    mx["joint_results"] = None
    return {"config": config}


@app.post("/api/multiomics/joint/run")
def multiomics_joint_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if _multiomics_state(session)["joint_config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_run_multiomics_joint, session)
    return {"job_id": job_id}


@app.get("/api/multiomics/joint/download")
def multiomics_joint_download(database: str, table: str = "combined",
                              x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = _multiomics_state(session)["joint_results"]
    if not results or database not in results["blocks"]:
        raise HTTPException(status_code=400, detail="No joint result available for that database.")
    block = results["blocks"][database]
    df = block["partial_table"] if table == "partial" else block["table"]
    return df_to_csv_response(df, f"multiomics_joint_{database}_{'partial' if table == 'partial' else 'combined'}.csv")


@app.post("/api/multiomics/{omics_type}/upload")
async def multiomics_upload(omics_type: str, file: UploadFile = File(...),
                            x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    if omics_type not in mx["selected_omics"]:
        raise HTTPException(status_code=400, detail=f"'{omics_type}' is not one of the selected omics.")
    raw_bytes = await file.read()
    try:
        if omics_type == "metabolomics":
            df = _load_metabolomics_list_upload(raw_bytes)
            id_columns = [c for c in ("KEGG ID", "PubChem ID", "InChIKey", "HMDB ID", "ChEBI ID",
                                      "LIPID MAPS ID", "SMILES", "Formula") if c in df.columns]
            parsed = {"df": df, "kind": "list", "identifier_type": None,
                      "column_mapping": {"identifier_columns": id_columns}, "n_duplicates_removed": 0}
        else:
            parsed = multiomics_retriever.parse_upload(omics_type, raw_bytes)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    entry = _fresh_multiomics_input()
    entry.update({k: parsed[k] for k in ("df", "kind", "identifier_type", "column_mapping", "n_duplicates_removed")})
    mx["inputs"][omics_type] = entry
    _multiomics_clear_results(mx)
    df = parsed["df"]
    return {
        "omics_type": omics_type, "kind": parsed["kind"], "n_rows": len(df), "columns": list(df.columns),
        "identifier_type": parsed["identifier_type"], "column_mapping": parsed["column_mapping"],
        "n_duplicates_removed": parsed["n_duplicates_removed"], "preview": df_payload(df.head(50)),
    }


def _multiomics_resolve_summary(omics_type: str, entry: dict) -> dict:
    """Small JSON-safe description of one omics' resolved table."""
    resolved, meta = entry["resolved_df"], entry["meta"] or {}
    summary = {"label": multiomics_engine.OMICS_LABELS[omics_type], "kind": entry["kind"],
               "n_input_rows": int(meta.get("n_input_rows", len(entry["df"]))), "meta": meta}
    if omics_type == "metabolomics":
        matched = entry["matched_df"]
        summary["n_usable"] = int((matched["compound_key"].astype(str).str.strip() != "").sum())
    elif omics_type == "genomics":
        summary["n_usable"] = len(multiomics_retriever.usable_genomics_rows(resolved))
    elif omics_type == "transcriptomics":
        summary["n_usable"] = int((resolved["gene_symbol"].astype(str).str.strip() != "").sum())
    else:
        summary["n_usable"] = int((resolved["accession"].astype(str).str.strip() != "").sum())
    return summary


def _job_resolve_multiomics(job_id: str, session: dict, taxon_id: str, taxon_name: str) -> dict:
    """Resolve every selected omics with that omics' OWN resolver. Everything is
    committed together, on success only (same discipline as _job_resolve_de): a
    failure part-way leaves the previous organism and tables untouched."""
    mx = _multiomics_state(session)
    selected = list(mx["selected_omics"])
    staged: dict = {}
    for i, omics_type in enumerate(selected):
        start, span = i / len(selected), 1.0 / len(selected)
        label = multiomics_engine.OMICS_LABELS[omics_type]
        update_job(job_id, progress=start, message=f"Resolving {label} identifiers...")
        entry = mx["inputs"][omics_type]
        if omics_type == "metabolomics":
            resolved_df, matched_df, meta = _resolve_and_match_compound_names(
                job_id, entry["df"], progress_start=start, progress_span=span)
            staged[omics_type] = {"resolved_df": resolved_df, "matched_df": matched_df, "meta": meta}
        else:
            def _progress(current, total, text, _s=start, _p=span, _l=label):
                update_job(job_id, progress=_s + _p * (current / max(total, 1)), message=f"{_l}: {text}")
            out = multiomics_retriever.resolve_upload(
                omics_type, entry["df"], taxon_id, identifier_type=entry["identifier_type"],
                progress_callback=_progress)
            staged[omics_type] = {"resolved_df": out["resolved_df"], "matched_df": None, "meta": out["meta"]}
    for omics_type, got in staged.items():
        mx["inputs"][omics_type].update(got)
    mx["taxon_id"], mx["taxon_name"] = taxon_id, taxon_name
    mx["kegg_organism_confirmation"] = None
    _multiomics_clear_results(mx)
    update_job(job_id, progress=1.0, message="Done")
    return {"omics": {o: _multiomics_resolve_summary(o, mx["inputs"][o]) for o in selected},
            "taxon_id": taxon_id, "taxon_name": taxon_name}


def _serialize_multiomics_results(results: dict) -> dict:
    """JSON-safe view of the joint run; the DataFrames stay in session state."""
    blocks = {}
    for database, block in results["blocks"].items():
        blocks[database] = {
            "table": df_payload(block["table"]), "partial_table": df_payload(block["partial_table"]),
            "per_omics": block["per_omics"], "omics_combined": block["omics_combined"],
            "omics_excluded": block["omics_excluded"], "n_terms_combined": block["n_terms_combined"],
            "n_terms_partial": block["n_terms_partial"], "n_terms_no_overlap": block["n_terms_no_overlap"],
            "note": block["note"], "caveat": block["caveat"], "method": block["method"],
            "columns": block["columns"], "correlation": block.get("correlation"),
        }
    return {"blocks": blocks, "summary": results["summary"]}


def _multiomics_resolve_kegg(session: dict, mx: dict, databases, job_id: "str | None" = None) -> tuple:
    """The KEGG organism is resolved ONCE and the same code goes into every omics'
    KEGG builder (commit 27ce033 was two paths resolving it separately). Returns
    (code, name, note, caveat); a non-empty note means no organism could be used."""
    kegg_code, kegg_name, kegg_note, kegg_caveat = "", "", None, None
    if "KEGG" in databases:
        if job_id:
            update_job(job_id, progress=0.02, message="Resolving the KEGG organism...")
        kegg_org = transcriptomics_engine.resolve_kegg_organism_code(
            mx["taxon_id"], mx["taxon_name"], verbose=False,
            confirmed_organism_code=_kegg_confirmed_code(session, "multiomics"))
        kegg_note, kegg_caveat = kegg_org.get("note"), kegg_org.get("caveat")
        if not kegg_note:
            kegg_code, kegg_name = kegg_org["organism_code"], kegg_org.get("organism_name", "")
    return kegg_code, kegg_name, kegg_note, kegg_caveat


def _multiomics_prepare_omics(mx: dict, selected: list, de_cutoffs: dict, kegg_code: str,
                              job_id: "str | None" = None, assign_kegg_ids: bool = True) -> tuple:
    """Per-omics selection (once, shared by every database): the hit set each omics
    contributes. Used by Joint Pathway Enrichment and by Cross-omics Pathway Mapping,
    so both tiles mean the same thing by "hit". Returns (prep, omics_summary)."""
    prep: dict = {}
    omics_summary: dict = {}
    for o in selected:
        entry = mx["inputs"][o]
        resolved = entry["resolved_df"]
        if o == "genomics":
            table = resolved
            if kegg_code and assign_kegg_ids:
                if job_id:
                    update_job(job_id, progress=0.05, message="Assigning KEGG gene IDs to the Genomics list...")
                table = multiomics_retriever.assign_genomics_kegg_ids(resolved, kegg_code)
            usable = multiomics_retriever.usable_genomics_rows(table)
            prep[o] = {"final_df": usable}
            omics_summary[o] = {"n_uploaded": len(resolved), "n_in_list": len(usable),
                                "n_left_out": len(resolved) - len(usable)}
        elif o in multiomics_retriever.DE_OMICS:
            cut = de_cutoffs[o]
            select = (multiomics_retriever.select_transcriptomics if o == "transcriptomics"
                      else multiomics_retriever.select_proteomics)
            sel = select(resolved, cut["padj"], cut["lfc"])
            prep[o] = sel
            omics_summary[o] = {"n_tested": sel["n_tested"], "n_in_list": sel["n_significant"],
                                "padj_cutoff": cut["padj"], "lfc_cutoff": cut["lfc"]}
        else:
            matched = entry["matched_df"]
            prep[o] = {"matched_df": matched, "crosswalk": metabolomics_engine.build_compound_crosswalk(verbose=False)}
            omics_summary[o] = {"n_uploaded": len(matched),
                                "n_in_list": int((matched["compound_key"].astype(str).str.strip() != "").sum())}
        omics_summary[o]["label"] = multiomics_engine.OMICS_LABELS[o]
    return prep, omics_summary


def _multiomics_database_inputs(mx: dict, selected: list, prep: dict, database: str, kegg_code: str,
                                kegg_name: str, kegg_note: "str | None", kegg_caveat: "str | None") -> dict:
    """One database's per-omics {term_map, list_set, universe_set, ...}, each built by that
    omics' OWN universe builder with the ONE shared KEGG organism code."""
    taxon_id, taxon_name = mx["taxon_id"], mx["taxon_name"]
    inputs: dict = {}
    for o in selected:
        if database == "KEGG" and kegg_note:
            inputs[o] = multiomics_engine._empty_input(kegg_note, kegg_caveat)
            continue
        if o == "genomics":
            inputs[o] = multiomics_engine.build_genomics_input(
                prep[o]["final_df"], database, taxon_id, taxon_name, kegg_code)
        elif o == "transcriptomics":
            inputs[o] = multiomics_engine.build_transcriptomics_input(
                prep[o]["universe_df"], prep[o]["significant"], database, taxon_name, kegg_code)
        elif o == "proteomics":
            inputs[o] = multiomics_engine.build_proteomics_input(
                mx["inputs"][o]["resolved_df"], prep[o]["significant"], database,
                taxon_id, taxon_name, kegg_code)
        else:
            inputs[o] = multiomics_engine.build_metabolomics_input(
                prep[o]["matched_df"], prep[o]["crosswalk"], database, taxon_id, taxon_name,
                kegg_code, kegg_name)
    return inputs


def _job_run_multiomics_joint(job_id: str, session: dict) -> dict:
    update_job(job_id, message="Starting joint pathway enrichment...")
    mx = _multiomics_state(session)
    config = mx["joint_config"]
    if config is None:
        raise ValueError("Configure the analysis first.")
    selected = list(mx["selected_omics"])
    taxon_id, taxon_name = mx["taxon_id"], mx["taxon_name"]
    databases = config["databases"]

    kegg_code, kegg_name, kegg_note, kegg_caveat = _multiomics_resolve_kegg(session, mx, databases, job_id)
    prep, omics_summary = _multiomics_prepare_omics(mx, selected, config["de_cutoffs"], kegg_code, job_id)

    blocks: dict = {}
    for d_index, database in enumerate(databases):
        update_job(job_id, progress=0.1 + 0.85 * d_index / len(databases),
                   message=f"Building {database} pathway sets for each omics...")
        inputs = _multiomics_database_inputs(mx, selected, prep, database, kegg_code, kegg_name,
                                             kegg_note, kegg_caveat)
        update_job(job_id, progress=0.1 + 0.85 * (d_index + 0.6) / len(databases),
                   message=f"Combining {database} p-values across omics...")
        block = multiomics_engine.run_joint_pathway_enrichment(
            inputs, database, config["method"], config["min_term_size"], config["max_term_size"])
        caveats = [c for c in (block.get("caveat"), kegg_caveat if database == "KEGG" else None,
                               multiomics_engine.joint_caveat(config["method"])) if c]
        block["caveat"] = " ".join(dict.fromkeys(caveats))
        blocks[database] = block

    results = {"blocks": blocks, "summary": {
        "omics": selected, "omics_summary": omics_summary, "taxon_id": taxon_id, "taxon_name": taxon_name,
        "kegg_organism_code": kegg_code, "kegg_organism_name": kegg_name,
        "method": config["method"],
        "method_label": multiomics_engine.COMBINE_METHOD_LABELS[config["method"]],
    }}
    mx["joint_results"] = results
    update_job(job_id, progress=1.0, message="Done")
    return _serialize_multiomics_results(results)


# ---------------------------------------------------------------------------
# Multi-omics - Cross-omics Pathway Mapping (descriptive; see Multiomics_mapping).
# Shares the uploads / resolution / organism / hit rule of Joint Pathway
# Enrichment (same Multi-omics session sub-dict), keeps its own config and
# results under session["multiomics"]["mapping"].
# ---------------------------------------------------------------------------

_MAPPING_PATHWAY_ID_RE = {"KEGG": re.compile(r"^[a-z]{2,5}\d{5}$"), "Reactome": re.compile(r"^R-[A-Z]{3}-\d+$")}


def _multiomics_parse_cutoffs(payload: dict, selected: list) -> dict:
    """padj / |log2FC| cutoffs for the differential-expression omics among `selected`."""
    cutoffs = {}
    raw_cutoffs = payload.get("de_cutoffs") or {}
    for o in multiomics_retriever.DE_OMICS:
        if o not in selected:
            continue
        raw = raw_cutoffs.get(o) or {}
        try:
            padj = float(raw.get("padj", multiomics_retriever.PADJ_CUTOFF_DEFAULT))
            lfc = abs(float(raw.get("lfc", multiomics_retriever.LFC_CUTOFF_DEFAULT)))
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="padj and log2FC cutoffs must be numbers.")
        if not (0 < padj <= 1):
            raise HTTPException(status_code=400, detail="The padj cutoff must be above 0 and at most 1.")
        cutoffs[o] = {"padj": padj, "lfc": lfc}
    return cutoffs


@app.post("/api/multiomics/mapping/configure")
def multiomics_mapping_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    selected = mx["selected_omics"]
    if not selected or any(mx["inputs"][o]["resolved_df"] is None for o in selected):
        raise HTTPException(status_code=400, detail="Resolve every selected omics' identifiers first.")
    mx["mapping"] = _fresh_multiomics_mapping()
    mx["mapping"]["config"] = {"de_cutoffs": _multiomics_parse_cutoffs(payload, selected)}
    return {"config": mx["mapping"]["config"]}


def _multiomics_mapping_database(payload: dict) -> str:
    database = str(payload.get("database") or "")
    if database not in multiomics_mapping.MAP_DATABASES:
        raise HTTPException(status_code=400, detail="Choose KEGG or Reactome.")
    return database


@app.post("/api/multiomics/mapping/pathways")
def multiomics_mapping_pathways(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    if mx["mapping"]["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the mapping first.")
    database = _multiomics_mapping_database(payload)
    job_id = submit_job(_job_multiomics_mapping_pathways, session, database)
    return {"job_id": job_id}


def _multiomics_build_picker(session: dict, database: str, job_id: "str | None" = None) -> dict:
    """The pathway picker: every pathway at least one omics hit, ranked by how many omics hit it. Built from the same
    per-omics sets the joint test uses (coverage, no statistics). The rows are kept in the session so an overview
    pathway can show its sub-pathways with the same hit counts without building this twice."""
    mx = _multiomics_state(session)
    cfg = mx["mapping"]["config"]
    if cfg is None:
        raise ValueError("Configure the mapping first.")
    selected = list(mx["selected_omics"])
    kegg_code, kegg_name, kegg_note, kegg_caveat = _multiomics_resolve_kegg(session, mx, [database], job_id)
    prep, omics_summary = _multiomics_prepare_omics(mx, selected, cfg["de_cutoffs"], kegg_code, job_id,
                                                    assign_kegg_ids=(database == "KEGG"))
    if job_id:
        update_job(job_id, progress=0.3, message=f"Building {database} pathway sets for each omics...")
    inputs = _multiomics_database_inputs(mx, selected, prep, database, kegg_code, kegg_name, kegg_note, kegg_caveat)
    joint = (mx["joint_results"] or {}).get("blocks", {}).get(database)
    joint_fdr = {}
    if joint is not None and len(joint["table"]):
        joint_fdr = dict(zip(joint["table"]["Term ID"].astype(str), joint["table"]["Combined FDR"]))
    rows = multiomics_mapping.build_pathway_picker(inputs, selected, joint_fdr)
    notes = [o_in["note"] for o_in in inputs.values() if o_in.get("note")]
    mx["mapping"]["pickers"][database] = {"kegg_code": kegg_code, "rows": rows}
    return {"database": database, "rows": rows, "n_pathways": len(rows), "omics": selected,
            "omics_summary": omics_summary, "kegg_organism_code": kegg_code, "kegg_organism_name": kegg_name,
            "note": " ".join(dict.fromkeys(notes)) or None, "caveat": kegg_caveat if database == "KEGG" else None,
            "has_joint_fdr": bool(joint_fdr)}


def _job_multiomics_mapping_pathways(job_id: str, session: dict, database: str) -> dict:
    update_job(job_id, message=f"Listing {database} pathways hit by your omics...")
    built = _multiomics_build_picker(session, database, job_id)
    columns = multiomics_mapping.picker_columns(built["omics"])
    out = {k: v for k, v in built.items() if k != "rows"}
    out["table"] = df_payload(pd.DataFrame(built["rows"], columns=columns))
    update_job(job_id, progress=1.0, message="Done")
    return out


@app.post("/api/multiomics/mapping/run")
def multiomics_mapping_run(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    if mx["mapping"]["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the mapping first.")
    database = _multiomics_mapping_database(payload)
    pathway_id = str(payload.get("pathway_id") or "").strip()
    if not _MAPPING_PATHWAY_ID_RE[database].match(pathway_id):
        raise HTTPException(status_code=400, detail=f"'{pathway_id}' is not a valid {database} pathway ID.")
    job_id = submit_job(_job_run_multiomics_mapping, session, database, pathway_id,
                        str(payload.get("pathway_name") or ""))
    return {"job_id": job_id}


def _serialize_map_node(n: dict, selected: list) -> dict:
    """JSON view of one mapped box. Only the omics that touched it carry a status (a global map has ~5,000 boxes and
    most are touched by nobody - four 'absent' dicts each would triple the payload); a polyline and a group's
    components travel with it."""
    out = {"id": n["id"], "kind": n["kind"], "label": n["label"], "x": n["x"], "y": n["y"], "w": n["w"], "h": n["h"],
           "status": {o: n["status"][o] for o in selected if n["status"][o]["state"] != "absent"}}
    if "path" in n:
        out["path"] = n["path"]
    if "components" in n:
        out["components"] = n["components"]
    if "cls" in n:
        out["cls"] = n["cls"]
    return out


def _job_run_multiomics_mapping(job_id: str, session: dict, database: str, pathway_id: str,
                                pathway_name: str = "") -> dict:
    update_job(job_id, message="Starting pathway mapping...")
    mx = _multiomics_state(session)
    cfg = mx["mapping"]["config"]
    if cfg is None:
        raise ValueError("Configure the mapping first.")
    selected = list(mx["selected_omics"])
    # The KEGG organism is needed by KEGG itself and, as an identifier bridge, by Reactome.
    kegg_code, kegg_name, kegg_note, kegg_caveat = _multiomics_resolve_kegg(session, mx, ["KEGG"], job_id)
    if database == "KEGG" and kegg_note:
        raise ValueError(kegg_note)

    update_job(job_id, progress=0.1, message=f"Fetching the {database} pathway structure...")
    fetched = multiomics_retriever.fetch_pathway_structure(database, pathway_id)
    result = {"database": database, "pathway_id": pathway_id, "pathway_name": pathway_name,
              "status": fetched["status"], "note": fetched.get("note"),
              "restricted_to": fetched.get("restricted_to"), "omics": selected,
              "nodes": [], "reactions": [], "summary": None, "members": None, "links": None,
              "interactions": None, "relations": [], "decor_lines": [], "line_graphics": False, "classes": [],
              "subpathways": None, "caveat": multiomics_mapping.MAP_CAVEAT}
    if fetched["status"] == "overview":
        # Reactome draws this pathway only as links to its sub-pathways: nothing to put a hit on. List them, with the
        # same hit counts the picker shows, so one click reaches a real diagram.
        result["pathway_name"] = fetched["structure"].get("title") or pathway_name
        picker = mx["mapping"]["pickers"].get(database)
        if not picker or "rows" not in picker:
            update_job(job_id, progress=0.3, message="Counting your hits per sub-pathway...")
            picker = {"rows": _multiomics_build_picker(session, database, job_id)["rows"]}
        rows = multiomics_mapping.overview_rows(fetched["subpathways"], picker["rows"], selected)
        result["subpathways"] = df_payload(pd.DataFrame(rows, columns=multiomics_mapping.picker_columns(selected)))
    if fetched["status"] != "ok":
        mx["mapping"]["results"] = None
        update_job(job_id, progress=1.0, message="Done")
        return result
    structure = fetched["structure"]
    result["pathway_name"] = structure.get("title") or pathway_name

    update_job(job_id, progress=0.35, message="Matching your omics to the pathway's members...")
    prep, _ = _multiomics_prepare_omics(mx, selected, cfg["de_cutoffs"], kegg_code, job_id,
                                        assign_kegg_ids=(database == "KEGG"))
    gene_map = genomics_ora.build_kegg_gene_map(kegg_code) if kegg_code else None
    features: dict = {}
    for o in selected:
        if o == "genomics":
            features[o] = multiomics_mapping.gene_features_genomics(prep[o]["final_df"], gene_map)
        elif o == "transcriptomics":
            features[o] = multiomics_mapping.gene_features_de(o, prep[o]["universe_df"], prep[o]["significant"], gene_map)
        elif o == "proteomics":
            resolved = mx["inputs"][o]["resolved_df"]
            universe_df = resolved[resolved["accession"].astype(str).str.strip() != ""]
            features[o] = multiomics_mapping.gene_features_de(o, universe_df, prep[o]["significant"], gene_map)
        else:
            features[o] = multiomics_mapping.compound_features(prep[o]["matched_df"])
    if database == "KEGG":
        try:
            multiomics_mapping.relabel_kegg_compounds(structure, metabolomics_engine.build_kegg_compound_name_index())
        except Exception:
            pass  # labels stay as KEGG C-numbers; the mapping itself does not depend on them
    skeletons = {}
    if database == "Reactome" and "metabolomics" in selected:
        update_job(job_id, progress=0.6, message="Bridging Reactome's ChEBI IDs to your compounds...")
        skeletons = multiomics_retriever.fetch_chebi_skeletons(multiomics_mapping.structure_chebi_ids(structure))
    mapped = multiomics_mapping.map_omics_to_structure(structure, features, database, skeletons)
    classes = []
    if structure.get("line_graphics") and database == "KEGG":
        # KEGG paints each line of a global map in its pathway CLASS's colour; name the colours from KEGG's own classification
        update_job(job_id, progress=0.8, message="Naming the pathway classes of the map...")
        map_to_class = multiomics_retriever.fetch_pathway_classes()
        universe = genomics_ora.build_kegg_universe(kegg_code)
        gene_to_maps: dict = {}
        for pid, info in (universe.get("pathways") or {}).items():
            for gene in info["genes"]:
                gene_to_maps.setdefault(gene, set()).add(pid[-5:])
        color_class = multiomics_mapping.assign_color_classes(mapped, gene_to_maps, map_to_class)
        for n in mapped:
            if n["kind"] == "gene" and n.get("color"):
                n["cls"] = color_class.get(n["color"], multiomics_mapping.OTHER_CLASS)
            elif n["kind"] == "link" and n.get("pathway") in map_to_class:
                n["cls"] = map_to_class[n["pathway"]]
        classes = multiomics_mapping.class_legend(mapped, color_class)
    links = multiomics_mapping.reaction_links(structure, mapped, selected)
    members = multiomics_mapping.member_rows(mapped, selected)
    interactions = multiomics_mapping.interaction_links(structure, mapped, selected)
    members_df = pd.DataFrame(members, columns=multiomics_mapping.MEMBER_COLUMNS)
    links_df = pd.DataFrame(links, columns=multiomics_mapping.LINK_COLUMNS)
    interactions_df = pd.DataFrame(interactions, columns=multiomics_mapping.INTERACTION_COLUMNS)

    result["nodes"] = [_serialize_map_node(n, selected) for n in mapped]
    result["reactions"] = structure["reactions"]
    result["relations"] = structure.get("relations") or []
    result["classes"] = classes
    result["decor_lines"] = structure.get("decor_lines") or []
    result["line_graphics"] = bool(structure.get("line_graphics"))
    result["summary"] = multiomics_mapping.summarize_mapping(mapped, links, selected, interactions)
    result["members"], result["links"] = df_payload(members_df), df_payload(links_df)
    result["interactions"] = df_payload(interactions_df)
    notes = []
    if database == "Reactome" and "transcriptomics" in selected and not kegg_code:
        notes.append("No KEGG organism was resolved, so Transcriptomics genes could only be matched to Reactome "
                     "entities by gene symbol (less precise than through UniProt).")
    if "metabolomics" in selected and database == "Reactome":
        notes.append("Metabolomics compounds are matched to Reactome chemicals by ChEBI ID, or - where Reactome cites a "
                     "different protonation state of the same compound - by InChIKey skeleton (each match says which).")
    result["note"] = " ".join(notes) or None
    mx["mapping"]["results"] = {"database": database, "pathway_id": pathway_id, "members_df": members_df,
                                "links_df": links_df, "interactions_df": interactions_df}
    update_job(job_id, progress=1.0, message="Done")
    return result


@app.get("/api/multiomics/mapping/download")
def multiomics_mapping_download(table: str = "members", x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    res = _multiomics_state(session)["mapping"]["results"]
    if not res:
        raise HTTPException(status_code=400, detail="No pathway mapping available.")
    table = table if table in ("links", "interactions") else "members"
    df = res[f"{table}_df"]
    return df_to_csv_response(df, f"multiomics_mapping_{res['pathway_id']}_{table}.csv")


# ---------------------------------------------------------------------------
# Multi-omics - Correlation-based Integration (see Multiomics_correlation for the statistics)
#
# Inputs are SAMPLE x FEATURE matrices (the same samples in every file), uploaded here for Transcriptomics, Proteomics
# and Metabolomics; Genomics (a gene list) cannot take part. State: session["multiomics"]["correlation"].
# ---------------------------------------------------------------------------

MULTIOMICS_CORR_EXAMPLES = {
    "transcriptomics": "example_corr_transcriptomics.csv",
    "proteomics": "example_corr_proteomics.csv",
    "metabolomics": "example_corr_metabolomics.csv",
}
_CORR_PAIR_COLUMNS = ["Feature A", "Feature B", "Omics A", "Omics B", "rho", "p-value", "FDR", "n"]
_CORR_GENE_COLUMNS = ["Gene", "rho", "p-value", "FDR", "n", "Transcript id(s)", "Protein id(s)"]


def _corr_state(session: dict) -> dict:
    return session["multiomics"]["correlation"]


def _corr_clear_downstream(cs: dict) -> None:
    """Drop what was built on the current uploads / resolution (the alignment and the result)."""
    cs["aligned"] = None
    cs["results"] = None


@app.post("/api/multiomics/correlation/select")
def multiomics_correlation_select(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    requested = payload.get("omics_types") or []
    if "genomics" in requested:
        raise HTTPException(status_code=400, detail="Genomics is a gene list, not a quantitative table: it cannot be "
                                                    "correlated across samples. Choose Transcriptomics, Proteomics and / or Metabolomics.")
    chosen = [o for o in multiomics_correlation.OMICS_ORDER if o in requested]
    if len(chosen) < 2:
        raise HTTPException(status_code=400, detail="Select at least two omics to correlate.")
    cs["selected"] = chosen
    cs["uploads"] = {o: v for o, v in cs["uploads"].items() if o in chosen}
    cs["resolved"] = {o: v for o, v in cs["resolved"].items() if o in chosen}
    _corr_clear_downstream(cs)
    return {"selected": chosen}


@app.get("/api/multiomics/correlation/example")
def multiomics_correlation_example(omics_type: str = "transcriptomics"):
    fname = MULTIOMICS_CORR_EXAMPLES.get(str(omics_type).strip().lower())
    if not fname:
        raise HTTPException(status_code=404, detail="Unknown omics type.")
    path = os.path.join(BASE_DIR, "Data3", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


@app.post("/api/multiomics/correlation/{omics_type}/upload")
async def multiomics_correlation_upload(omics_type: str, file: UploadFile = File(...), orientation: str = Form("auto"),
                                        x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if omics_type not in cs["selected"]:
        raise HTTPException(status_code=400, detail=f"'{omics_type}' is not one of the selected omics.")
    raw_bytes = await file.read()
    try:
        parsed = multiomics_correlation.parse_matrix(raw_bytes, orientation=orientation)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    cs["uploads"][omics_type] = {"df": parsed["df"], "meta": parsed["meta"], "filename": getattr(file, "filename", "")}
    cs["resolved"].pop(omics_type, None)
    _corr_clear_downstream(cs)
    meta = parsed["meta"]
    return {"omics_type": omics_type, "n_samples": meta["n_samples"], "n_features": meta["n_features"],
            "orientation": meta["orientation"], "orientation_reason": meta["orientation_reason"],
            "set_aside_columns": meta["set_aside_columns"], "first_samples": meta["first_samples"],
            "n_repeated_feature_ids": meta["n_repeated_feature_ids"], "n_empty_dropped": meta["n_empty_dropped"],
            "first_features": [str(c) for c in parsed["df"].columns[:5]]}


def _job_correlation_resolve(job_id: str, session: dict, taxon_id: str, taxon_name: str, case_insensitive: bool) -> dict:
    """Give the transcript / protein features their gene symbols (each omics' own mapping step), then keep only the
    samples every file has. Committed together, on success only."""
    cs = _corr_state(session)
    selected = list(cs["selected"])
    staged: dict = {}
    for i, o in enumerate(selected):
        start, span = i / len(selected), 1.0 / len(selected)
        label = multiomics_correlation.OMICS_LABELS[o]
        update_job(job_id, progress=start, message=f"Resolving {label} identifiers...")
        df = cs["uploads"][o]["df"]
        if o == "metabolomics":
            # nothing here needs a compound's database ids: features are correlated by name as written
            staged[o] = {"gene_of": {}, "n_features": int(df.shape[1]), "n_mapped": None, "n_unmapped": None, "identifier_type": None}
        else:
            def _progress(current, total, text, _s=start, _p=span, _l=label):
                update_job(job_id, progress=_s + _p * (current / max(total, 1)), message=f"{_l}: {text}")
            # a repeated id was renamed "X (2)" on upload; it is resolved by its original id "X"
            original = cs["uploads"][o]["meta"].get("original_ids") or {}
            ids = [str(c) for c in df.columns]
            res = multiomics_retriever.resolve_matrix_features(o, list(dict.fromkeys(original.get(c, c) for c in ids)), taxon_id,
                                                               progress_callback=_progress)
            gene_of = {c: res["gene_of"][original.get(c, c)] for c in ids if original.get(c, c) in res["gene_of"]}
            extra = {k: {c: res[k][original.get(c, c)] for c in ids if original.get(c, c) in res[k]}
                     for k in ("entrez_of", "acc_of") if k in res}
            staged[o] = {**res, **extra, "gene_of": gene_of, "n_features": len(ids), "n_mapped": len(gene_of),
                         "n_unmapped": len(ids) - len(gene_of)}
    update_job(job_id, progress=0.95, message="Matching the samples across the files...")
    mats = {}
    for o in selected:
        m = cs["uploads"][o]["df"].copy()
        m.attrs["alt_samples"] = cs["uploads"][o]["meta"]["alt_samples"]
        mats[o] = m
    aligned = multiomics_correlation.align_samples(mats, case_insensitive=case_insensitive)
    cs["resolved"] = staged
    cs["aligned"] = {"samples": aligned["samples"], "mats": aligned["mats"], "report": aligned["report"],
                     "warnings": aligned["warnings"], "case_insensitive": case_insensitive}
    cs["taxon_id"], cs["taxon_name"] = taxon_id, taxon_name
    cs["results"] = None
    update_job(job_id, progress=1.0, message="Done")
    return _corr_alignment_payload(cs)


def _corr_alignment_payload(cs: dict) -> dict:
    al = cs["aligned"]
    omics = {}
    for o in cs["selected"]:
        up, rs = cs["uploads"][o], cs["resolved"].get(o) or {}
        omics[o] = {"label": multiomics_correlation.OMICS_LABELS[o], "n_samples": al["report"][o]["n_samples"],
                    "n_used": al["report"][o]["n_used"], "dropped": al["report"][o]["dropped"],
                    "n_features": up["meta"]["n_features"], "orientation": up["meta"]["orientation"],
                    "n_mapped": rs.get("n_mapped"), "n_unmapped": rs.get("n_unmapped"),
                    "identifier_type": rs.get("identifier_type")}
    return {"omics": omics, "samples": al["samples"], "n_shared": len(al["samples"]), "warnings": al["warnings"]}


@app.post("/api/multiomics/correlation/resolve")
def multiomics_correlation_resolve(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if len(cs["selected"]) < 2:
        raise HTTPException(status_code=400, detail="Choose the omics to correlate first.")
    missing = [multiomics_correlation.OMICS_LABELS[o] for o in cs["selected"] if o not in cs["uploads"]]
    if missing:
        raise HTTPException(status_code=400, detail="Upload a table for: " + ", ".join(missing) + ".")
    taxon_id = str(payload.get("taxon_id") or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism before resolving identifiers.")
    job_id = submit_job(_job_correlation_resolve, session, taxon_id, str(payload.get("taxon_name") or "").strip(),
                        bool(payload.get("case_insensitive")))
    return {"job_id": job_id}


@app.post("/api/multiomics/correlation/configure")
def multiomics_correlation_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["aligned"] is None:
        raise HTTPException(status_code=400, detail="Resolve the identifiers and match the samples first.")
    cfg = dict(multiomics_correlation.CONFIG_DEFAULTS)
    method = payload.get("method", cfg["method"])
    if method not in multiomics_correlation.METHODS:
        raise HTTPException(status_code=400, detail="Unknown correlation method.")
    cfg["method"] = method
    cfg["log2"] = bool(payload.get("log2", False))

    def _num(key, cast, lo, hi, label):
        try:
            v = cast(payload.get(key, cfg[key]))
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail=f"{label} must be a number.")
        if not (lo <= v <= hi):
            raise HTTPException(status_code=400, detail=f"{label} must be between {lo} and {hi}.")
        return v

    cfg["top_variance"] = _num("top_variance", int, 50, 5000, "The number of top-variance features")
    cfg["fdr"] = _num("fdr", float, 1e-12, 1.0, "The FDR cutoff")
    cfg["rho_cut"] = _num("rho_cut", float, 0.0, 0.999, "The |rho| cutoff")
    cfg["min_samples"] = _num("min_samples", int, multiomics_correlation.MIN_SAMPLES, 1000, "The minimum number of samples per pair")
    cfg["heat_k"] = _num("heat_k", int, 5, 80, "The heatmap size")
    cfg["net_k"] = _num("net_k", int, 5, 150, "The number of network edges")
    if len(cs["aligned"]["samples"]) < cfg["min_samples"]:
        raise HTTPException(status_code=400, detail=f"Only {len(cs['aligned']['samples'])} samples are shared; the minimum per pair "
                                                    f"you set is {cfg['min_samples']}.")
    cs["config"] = cfg
    cs["results"] = None
    return {"config": cfg}


def _corr_rows_df(rows: list, columns: list) -> pd.DataFrame:
    return pd.DataFrame([{c: r[c] for c in columns} for r in rows], columns=columns)


def _job_correlation_run(job_id: str, session: dict) -> dict:
    cs = _corr_state(session)
    cfg = cs["config"]
    if cfg is None or cs["aligned"] is None:
        raise ValueError("Configure the analysis first.")
    update_job(job_id, progress=0.05, message="Cleaning the matrices...")
    clean, feature_notes = {}, {}
    for o in cs["selected"]:
        clean[o], feature_notes[o] = multiomics_correlation.prepare_matrix(
            cs["aligned"]["mats"][o], log2=cfg["log2"], max_missing=multiomics_correlation.MAX_MISSING_FRACTION, top_variance=None)
        if clean[o].shape[1] < 2:
            raise ValueError(f"{multiomics_correlation.OMICS_LABELS[o]}: fewer than two usable features remain.")
    gene_of = {o: (cs["resolved"].get(o) or {}).get("gene_of", {}) for o in cs["selected"]}
    update_job(job_id, progress=0.15, message="Correlating every pair of features (this is the slow step)...")
    res = multiomics_correlation.run_correlation(
        clean, gene_of, cfg,
        progress=lambda f, text: update_job(job_id, progress=0.15 + 0.7 * f, message=f"Correlating every pair of features - {text}"))
    update_job(job_id, progress=0.9, message="Preparing the tables...")
    tables: dict = {}
    out = {"config": res["config"], "samples": res["samples"], "omics": res["omics"], "notes": res["notes"],
           "warnings": cs["aligned"]["warnings"], "alignment": _corr_alignment_payload(cs), "network": res["network"],
           "feature_report": {o: {**feature_notes[o], **res["feature_report"][o]} for o in res["omics"]}, "blocks": {}, "concordance": None}
    for key, blk in res["blocks"].items():
        df = _corr_rows_df(blk["rows"], _CORR_PAIR_COLUMNS)
        tables[key] = pd.concat([df, pd.DataFrame({"Feature A id": [r["a_id"] for r in blk["rows"]], "Feature B id": [r["b_id"] for r in blk["rows"]]})], axis=1)
        out["blocks"][key] = {k: blk[k] for k in ("omics_a", "omics_b", "n_a", "n_b", "n_tests", "n_sig_fdr", "n_reported",
                                                  "truncated", "expected_false", "n_too_few_samples", "heat")}
        out["blocks"][key]["table"] = df_payload(df)
        out["blocks"][key]["ids"] = [[r["a_id"], r["b_id"]] for r in blk["rows"]]
    conc_mats = None
    if res["concordance"] is not None:
        c = res["concordance"]
        conc_mats = c.pop("_matrices")
        df = _corr_rows_df(c["rows"], _CORR_GENE_COLUMNS)
        tables["concordance"] = df
        out["concordance"] = {k: v for k, v in c.items() if k != "rows"}
        out["concordance"]["table"] = df_payload(df)
    cs["results"] = {"payload": out, "mats": clean, "conc_mats": conc_mats, "tables": tables, "labels": res["labels"]}
    update_job(job_id, progress=1.0, message="Done")
    return out


@app.post("/api/multiomics/correlation/run")
def multiomics_correlation_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_correlation_run, session)
    return {"job_id": job_id}


@app.get("/api/multiomics/correlation/scatter")
def multiomics_correlation_scatter(kind: str = "pair", block: str = "", a: str = "", b: str = "", gene: str = "",
                                   x_session_id: Optional[str] = Header(default=None)):
    """One pair's values per sample, for the scatter plot, with the correlation recomputed for the label."""
    session = session_dep(x_session_id)
    res = _corr_state(session)["results"]
    if not res:
        raise HTTPException(status_code=400, detail="No correlation result available.")
    cfg = res["payload"]["config"]
    try:
        if kind == "gene":
            if res["conc_mats"] is None:
                raise KeyError("no concordance")
            Tg, Pg = res["conc_mats"]
            vec = multiomics_correlation.pair_vectors(Tg, gene, Pg, gene)
            labels = {"x_label": "Transcript " + gene, "y_label": "Protein " + gene}
        else:
            oa, _, ob = block.partition("|")
            vec = multiomics_correlation.pair_vectors(res["mats"][oa], a, res["mats"][ob], b)
            def _name(o, f):
                # the gene symbol a protein / transcript resolved to, with its own id when they differ
                lab = (res.get("labels") or {}).get(o, {}).get(f, f)
                return f"{multiomics_correlation.OMICS_LABELS[o]}: {lab}" + (f" ({f})" if lab != f else "")
            labels = {"x_label": _name(oa, a), "y_label": _name(ob, b)}
    except KeyError:
        raise HTTPException(status_code=404, detail="Unknown feature for this result.")
    x = np.array([np.nan if v is None else v for v in vec["x"]])
    y = np.array([np.nan if v is None else v for v in vec["y"]])
    r, n = multiomics_correlation.paired_correlation(x[:, None], y[:, None], cfg["method"], 3)
    p = multiomics_correlation.pvalues_from_r(r, n)
    return {**vec, **labels, "method": cfg["method"],
            "rho": None if not np.isfinite(r[0]) else float(r[0]), "p": None if not np.isfinite(p[0]) else float(p[0]), "n": int(n[0])}


@app.get("/api/multiomics/correlation/download")
def multiomics_correlation_download(table: str = "concordance", x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    res = _corr_state(session)["results"]
    if not res or table not in res["tables"]:
        raise HTTPException(status_code=400, detail="No such result table.")
    return df_to_csv_response(res["tables"][table], f"multiomics_correlation_{table.replace('|', '_')}.csv")


# ---------------------------------------------------------------------------
# Multi-omics - Network-based Integration (see Multiomics_network)
#
# Shares the uploads / resolution / organism of the other Multi-omics tiles (the same session sub-dict, like Cross-omics
# Pathway Mapping); keeps its own config and results under session["multiomics"]["network"]. Two networks from one run:
# the STRING gene / protein network (seeds of every omics, equal weight per omics) and the Reactome pathway-overlap network.
# ---------------------------------------------------------------------------

def _fresh_multiomics_network() -> dict:
    return {"config": None, "results": None}


def _network_int(payload: dict, key: str, default: int, lo: int, hi: int, label: str) -> int:
    raw = payload.get(key, default)
    try:
        v = int(raw)
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail=f"{label} must be a whole number.")
    if not (lo <= v <= hi):
        raise HTTPException(status_code=400, detail=f"{label} must be between {lo} and {hi}.")
    return v


def _network_float(payload: dict, key: str, default: float, lo: float, hi: float, label: str) -> float:
    raw = payload.get(key, default)
    try:
        v = float(raw)
    except (TypeError, ValueError):
        raise HTTPException(status_code=400, detail=f"{label} must be a number.")
    if not (lo < v <= hi):
        raise HTTPException(status_code=400, detail=f"{label} must be above {lo} and at most {hi}.")
    return v


@app.post("/api/multiomics/network/configure")
def multiomics_network_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    mx = _multiomics_state(session)
    selected = mx["selected_omics"]
    if len(selected) < 2 or any(mx["inputs"][o]["resolved_df"] is None for o in selected):
        raise HTTPException(status_code=400, detail="Choose at least two omics and resolve every selected omics' identifiers first.")
    taxon_id = str(mx["taxon_id"] or "").strip()
    if not taxon_id:
        raise HTTPException(status_code=400, detail="Select an organism first.")
    try:
        in_string = genomics_net.is_taxon_in_string(taxon_id)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Could not reach STRING to check organism coverage: {exc}")
    if not in_string:
        raise HTTPException(status_code=400, detail=(
            f"STRING has no protein-interaction network for taxon {taxon_id} ({mx['taxon_name'] or 'this organism'}), and "
            "Network-based Integration needs one. Joint Pathway Enrichment and Cross-omics Pathway Mapping still work for this organism."))
    payload = payload or {}
    confidence = _network_int(payload, "confidence", genomics_net.DEFAULT_CONFIDENCE, 0, 1000, "The STRING confidence cutoff")
    if confidence not in genomics_net.CONFIDENCE_CHOICES:
        raise HTTPException(status_code=400, detail=f"The STRING confidence cutoff must be one of {list(genomics_net.CONFIDENCE_CHOICES)}.")
    n_perm = _network_int(payload, "n_permutations", multiomics_network.CONFIG_DEFAULTS["n_permutations"], 1, 10**7, "The number of permutations")
    if n_perm not in transcriptomics_network_engine.NETWORK_N_PERMUTATIONS_CHOICES:
        raise HTTPException(status_code=400, detail=f"The number of permutations must be one of {list(transcriptomics_network_engine.NETWORK_N_PERMUTATIONS_CHOICES)}.")
    module_db = str(payload.get("module_database") or multiomics_network.CONFIG_DEFAULTS["module_database"])
    if module_db not in proteomics_network_engine.PPI_ENRICHMENT_DATABASES:
        raise HTTPException(status_code=400, detail=f"Label modules with one of {list(proteomics_network_engine.PPI_ENRICHMENT_DATABASES)}.")
    cfg = {
        "confidence": confidence, "restart_prob": genomics_net.DEFAULT_RESTART_PROB, "n_permutations": n_perm,
        "fdr_alpha": _network_float(payload, "fdr_alpha", genomics_net.DEFAULT_FDR_ALPHA, 0.0, 1.0, "The significance threshold"),
        "mcl_inflation": genomics_net.DEFAULT_MCL_INFLATION,
        "min_module_size": _network_int(payload, "min_module_size", genomics_net.DEFAULT_MIN_MODULE_SIZE, 2, 1000, "The minimum module size"),
        "max_hits_for_modules": _network_int(payload, "max_hits_for_modules", transcriptomics_network_engine.DEFAULT_MAX_HITS_FOR_MODULES, 50, 20000, "The maximum number of hits clustered"),
        "module_database": module_db, "include_iea": bool(payload.get("include_iea", True)),
        "max_reactions": _network_int(payload, "max_reactions", multiomics_network.CONFIG_DEFAULTS["max_reactions"], 5, 1000, "The reaction limit for a metabolite"),
        "max_reaction_proteins": _network_int(payload, "max_reaction_proteins", multiomics_network.CONFIG_DEFAULTS["max_reaction_proteins"], 2, 500, "The protein limit for a reaction"),
        "pathway_merge_jaccard": _network_float(payload, "pathway_merge_jaccard", multiomics_network.CONFIG_DEFAULTS["pathway_merge_jaccard"], 0.29, 1.0, "The merge threshold for similar pathways"),
        "pathway_top_n": _network_int(payload, "pathway_top_n", multiomics_network.CONFIG_DEFAULTS["pathway_top_n"], 5, 100, "The number of pathways in the pathway network"),
        "pathway_min_jaccard": _network_float(payload, "pathway_min_jaccard", multiomics_network.CONFIG_DEFAULTS["pathway_min_jaccard"], 0.0, 1.0, "The overlap threshold"),
        "seed": 42, "de_cutoffs": _multiomics_parse_cutoffs(payload, selected),
    }
    mx["network"] = {"config": cfg, "results": None}
    return {"config": cfg}


@app.post("/api/multiomics/network/run")
def multiomics_network_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    if _multiomics_state(session)["network"]["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_multiomics_network, session)
    return {"job_id": job_id}


def _job_multiomics_network(job_id: str, session: dict) -> dict:
    mx = _multiomics_state(session)
    cfg = mx["network"]["config"]
    if cfg is None:
        raise ValueError("Configure the analysis first.")
    selected = list(mx["selected_omics"])
    taxon_id, taxon_name = mx["taxon_id"], mx["taxon_name"]
    inputs: dict = {}
    for o in selected:
        entry = mx["inputs"][o]
        if o == "genomics":
            inputs[o] = {"final_df": multiomics_retriever.usable_genomics_rows(entry["resolved_df"])}
        elif o in multiomics_retriever.DE_OMICS:
            inputs[o] = {"resolved_df": entry["resolved_df"]}
        else:
            inputs[o] = {"features": multiomics_mapping.compound_features(entry["matched_df"])}
    enzyme_maps = None
    if "metabolomics" in selected:
        update_job(job_id, progress=0.01, message="Reading Reactome's reaction files (metabolite -> proteins)...")
        enzyme_maps = multiomics_retriever.fetch_reactome_reaction_maps(taxon_id, taxon_name)

    def _progress(frac, label):
        update_job(job_id, progress=0.02 + 0.80 * float(frac), message=str(label))

    result = multiomics_network.run_string_network(
        inputs, taxon_id, taxon_name, cfg, enzyme_maps=enzyme_maps,
        skeleton_provider=multiomics_retriever.fetch_chebi_skeletons_cached, progress_callback=_progress)
    graph = result.pop("_graph", None)
    hit_ids = result.pop("_hit_ids", [])
    result.pop("_per_omics", None)
    result.pop("_seed_omics", None)
    edges = []
    if graph is not None and hit_ids:
        ranked = result["propagation"]["table"].sort_values(by="Propagated score", ascending=False, kind="stable")["STRING ID"].tolist()
        edges = _netprop_subnetwork_edges(graph, hit_ids, rank_by=ranked)

    update_job(job_id, progress=0.85, message="Building the Reactome pathway-overlap network...")
    pathway = {"nodes": [], "edges": [], "groups": [], "omics": [], "note": None}
    try:
        prep, _ = _multiomics_prepare_omics(mx, selected, cfg["de_cutoffs"], "", job_id, assign_kegg_ids=False)
        pinputs = _multiomics_database_inputs(mx, selected, prep, "Reactome", "", "", None, None)
        block = multiomics_engine.run_joint_pathway_enrichment(
            pinputs, "Reactome", "brown", multiomics_engine.MIN_TERM_SIZE_DEFAULT, multiomics_engine.MAX_TERM_SIZE_DEFAULT)
        pathway = multiomics_network.pathway_overlap_network(block, pinputs, cfg["pathway_top_n"], cfg["pathway_min_jaccard"],
                                                             merge_jaccard=cfg["pathway_merge_jaccard"])
        pathway["method_note"] = block.get("note")
    except Exception as exc:
        pathway["note"] = f"The pathway-overlap network could not be built: {exc}"

    metab_df = pd.DataFrame(result["metabolites"]) if result.get("metabolites") else None
    omics_df = pd.DataFrame(result["omics"])
    nodes_df = pd.DataFrame([{"Term ID": n["id"], "Pathway": n["name"], "Group": n["group"], "Combined p-value": n["combined_p"],
                              "Combined FDR": n["combined_fdr"], "Hit members": n["n_hit_members"],
                              "Similar pathways merged into this one": "; ".join(f"{m['name']} ({m['id']})" for m in n["merged"]),
                              **{f"{multiomics_engine.OMICS_LABELS[o]} hits": n["hits"].get(o) for o in pathway["omics"]}}
                             for n in pathway["nodes"]])
    edges_df = pd.DataFrame([{"Pathway A": e["a"], "Pathway B": e["b"], "Jaccard": e["jaccard"], "Shared hit members": e["shared"]}
                             for e in pathway["edges"]])
    mx["network"]["results"] = {"string": result, "pathway": pathway, "tables": {
        "hits": result["propagation"]["table"], "modules": result["modules"]["table"] if result["modules"] else None,
        "metabolites": metab_df, "omics": omics_df, "pathway_nodes": nodes_df, "pathway_edges": edges_df}}

    mods = result["modules"]
    payload = {
        "taxon_id": result["taxon_id"], "taxon_name": result["taxon_name"], "parameters": result["parameters"],
        "omics": result["omics"], "omics_used": result["omics_used"], "seeds": result["seeds"],
        "metabolites": df_payload(metab_df) if metab_df is not None else None,
        "propagation": {**{k: v for k, v in result["propagation"].items() if k != "table"},
                        "table": df_payload(result["propagation"]["table"].head(2000)), "edges": edges},
        "modules": None, "caveats": result["caveats"], "pathway_network": pathway,
    }
    if mods:
        payload["modules"] = {
            "table": df_payload(mods["table"]), "membership": mods["membership"], "universe_size": mods["universe_size"],
            "caveat": mods["caveat"], "selection_caveat": mods["selection_caveat"], "cap_info": mods["cap_info"], "omics": mods["omics"],
            "enrichment": {
                mid: ({"namespace": e["best"]["namespace"], "table": df_payload(e["best"]["block"]["table"].head(50)),
                       "candidates": {key: ({"table": df_payload(c["block"]["table"].head(50))} if c else None)
                                      for key, c in e["candidates"].items()}} if e else None)
                for mid, e in mods["enrichment"].items()},
        }
    update_job(job_id, progress=1.0, message="Done")
    return payload


@app.get("/api/multiomics/network/download")
def multiomics_network_download(which: str = "hits", key: str = "", db: str = "",
                                x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    results = _multiomics_state(session)["network"]["results"]
    if not results:
        raise HTTPException(status_code=400, detail="No network results available.")
    tables = results["tables"]
    names = {"hits": "hits", "modules": "modules", "metabolites": "metabolite_enzymes", "omics": "omics_seeds",
             "pathway_nodes": "pathway_network_nodes", "pathway_edges": "pathway_network_edges"}
    if which in names:
        df = tables.get(which)
        if df is None:
            raise HTTPException(status_code=400, detail=f"This run has no {which} table.")
        return df_to_csv_response(df, f"multiomics_network_{names[which]}.csv")
    if which == "enrichment":
        mods = results["string"].get("modules")
        if not mods:
            raise HTTPException(status_code=400, detail="This run did not include module detection.")
        entry = mods["enrichment"].get(key)
        if entry is None:
            raise HTTPException(status_code=400, detail=f"No enrichment for module {key!r}.")
        if db:
            candidate = entry["candidates"].get(db)
            if candidate is None:
                raise HTTPException(status_code=400, detail=f"No {db!r} result for module {key!r}.")
            block = candidate["block"]
        else:
            block = entry["best"]["block"]
        safe = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower()
        return df_to_csv_response(block["table"], f"multiomics_network_enrichment_{safe}.csv")
    raise HTTPException(status_code=400, detail=f"Unknown table: {which!r}")


# ---------------------------------------------------------------------------
# Multi-omics - MOFA / Factor-based Integration (see Multiomics_mofa)
#
# Shares the matrices of Correlation-based Integration (uploaded, resolved to gene symbols and matched on their samples there); keeps its own
# config and results under session["multiomics"]["mofa"]. A MOFA result is tied to the exact aligned data it was fitted on (a fingerprint): if
# the Correlation uploads change, the result is dropped, never shown against different data.
# ---------------------------------------------------------------------------

def _fresh_multiomics_mofa() -> dict:
    return {"config": None, "results": None, "sheet": None, "enrich": {}, "metab_match": None, "inputs": {}}


def _mofa_state(session: dict) -> dict:
    return session["multiomics"]["mofa"]


def _mofa_fingerprint(cs: dict) -> "str | None":
    al = cs.get("aligned")
    if not al:
        return None
    parts = [",".join(al["samples"])]
    for o in sorted(al["mats"]):
        m = al["mats"][o]
        parts.append(f"{o}:{m.shape[0]}x{m.shape[1]}:{float(np.nansum(m.to_numpy(dtype=float))):.6f}:{hash(tuple(map(str, m.columns))) & 0xffffffff}")
    return "|".join(parts)


def _mofa_checked(session: dict) -> dict:
    """The MOFA state, with any config / result / enrichment dropped when the Correlation data they were built on has changed."""
    ms = _mofa_state(session)
    fp = _mofa_fingerprint(_corr_state(session))
    if ms.get("fingerprint") != fp:
        ms.update({"config": None, "results": None, "enrich": {}, "metab_match": None, "inputs": {}, "fingerprint": fp})
        if fp is None:
            ms["sheet"] = None
    return ms


def _mofa_group_options(session: dict) -> dict:
    """The sample-group sources on offer: a label column found in a matrix upload (Group, Condition ...) or a column of the sample sheet."""
    cs = _corr_state(session)
    ms = _mofa_state(session)
    samples = list((cs.get("aligned") or {}).get("samples") or [])
    opts = []
    for o in cs["selected"]:
        gc = (cs["uploads"].get(o) or {}).get("meta", {}).get("group_columns") or {}
        for col, mp in gc.items():
            opts.append({"key": f"{o}:{col}", "label": f"{multiomics_correlation.OMICS_LABELS[o]} file, column '{col}'", "mapping": mp})
    for col, mp in ((ms.get("sheet") or {}).get("columns") or {}).items():
        opts.append({"key": f"sheet:{col}", "label": f"Sample sheet, column '{col}'", "mapping": mp})
    out = []
    for op in opts:
        vals = [op["mapping"].get(s) for s in samples]
        n_have = sum(1 for v in vals if v)
        out.append({"key": op["key"], "label": op["label"], "n_samples": n_have, "groups": sorted({v for v in vals if v})})
    return {"options": out, "mappings": {op["key"]: op["mapping"] for op in opts}}


MULTIOMICS_MOFA_EXAMPLES = {"transcriptomics": "example_mofa_transcriptomics.csv", "proteomics": "example_mofa_proteomics.csv",
                            "metabolomics": "example_mofa_metabolomics.csv"}


@app.get("/api/multiomics/mofa/example")
def multiomics_mofa_example(omics_type: str = "transcriptomics"):
    """The synthetic 36-sample example tables (the Correlation examples have only 12 samples, too few for MOFA)."""
    fname = MULTIOMICS_MOFA_EXAMPLES.get(str(omics_type).strip().lower())
    if not fname:
        raise HTTPException(status_code=404, detail="Unknown omics type.")
    path = os.path.join(BASE_DIR, "Data3", fname)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="Example file not found.")
    return FileResponse(path, filename=fname, media_type="text/csv")


@app.get("/api/multiomics/mofa/options")
def multiomics_mofa_options(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    ms = _mofa_checked(session)
    go = _mofa_group_options(session)
    return {"ready": cs["aligned"] is not None, "samples": list((cs["aligned"] or {}).get("samples") or []),
            "omics": list(cs["selected"]), "group_options": go["options"], "sheet": ({"filename": ms["sheet"]["filename"], "n_matched": ms["sheet"]["n_matched"]} if ms.get("sheet") else None),
            "defaults": multiomics_mofa.CONFIG_DEFAULTS, "min_samples": multiomics_mofa.MIN_SAMPLES, "warn_samples": multiomics_mofa.WARN_SAMPLES}


@app.post("/api/multiomics/mofa/groups")
async def multiomics_mofa_groups(file: UploadFile = File(...), x_session_id: Optional[str] = Header(default=None)):
    """A small sample sheet: the first column is the sample id, every other column is a candidate group label (Group, Batch, Sex ...)."""
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["aligned"] is None:
        raise HTTPException(status_code=400, detail="Upload and match the omics tables first (the sample sheet is matched to their samples).")
    raw = await file.read()
    try:
        df = pd.read_csv(io.BytesIO(raw), sep=None, engine="python", dtype=str)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Could not read the sample sheet: {exc}")
    if df.shape[1] < 2 or df.shape[0] < 2:
        raise HTTPException(status_code=400, detail="The sample sheet needs a sample-id column and at least one label column.")
    df.columns = [str(c).strip() for c in df.columns]
    ids = df.iloc[:, 0].astype(str).str.strip()
    samples = set(cs["aligned"]["samples"])
    matched = [s for s in ids if s in samples]
    if len(matched) < 2:
        raise HTTPException(status_code=400, detail=f"Only {len(matched)} of the sheet's sample ids are among the {len(samples)} samples of your tables "
                                                    f"(it must repeat the SAME ids; the first column is read as the sample id).")
    columns = {}
    for col in df.columns[1:]:
        vals = df[col].astype(str).str.strip()
        mp = {s: v for s, v in zip(ids, vals) if s in samples and v not in ("", "nan")}
        distinct = set(mp.values())
        if 2 <= len(distinct) <= max(2, len(mp) // 2):
            columns[col] = mp
    if not columns:
        raise HTTPException(status_code=400, detail="No column of the sample sheet looks like a group label (text with a few repeated values).")
    ms = _mofa_checked(session)
    ms["sheet"] = {"filename": getattr(file, "filename", ""), "columns": columns, "n_matched": len(set(matched))}
    ms["results"] = None
    ms["enrich"] = {}
    return {"filename": ms["sheet"]["filename"], "n_matched": len(set(matched)), "n_samples": len(samples), "columns": list(columns),
            "group_options": _mofa_group_options(session)["options"]}


@app.post("/api/multiomics/mofa/configure")
def multiomics_mofa_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["aligned"] is None:
        raise HTTPException(status_code=400, detail="Upload, resolve and match the omics tables first.")
    n = len(cs["aligned"]["samples"])
    if n < multiomics_mofa.MIN_SAMPLES:
        raise HTTPException(status_code=400, detail=f"MOFA needs at least {multiomics_mofa.MIN_SAMPLES} samples shared by every table; {n} were matched.")
    payload = payload or {}
    d = multiomics_mofa.CONFIG_DEFAULTS
    cfg = dict(d)
    cfg["factors"] = _network_int(payload, "factors", d["factors"], 2, multiomics_mofa.MAX_FACTORS, "The number of factors")
    cfg["convergence"] = str(payload.get("convergence", d["convergence"]))
    if cfg["convergence"] not in multiomics_mofa.CONVERGENCE_MODES:
        raise HTTPException(status_code=400, detail=f"The convergence mode must be one of {list(multiomics_mofa.CONVERGENCE_MODES)}.")
    cfg["top_variance"] = _network_int(payload, "top_variance", d["top_variance"], 50, 5000, "The number of most variable features")
    cfg["n_refits"] = _network_int(payload, "n_refits", d["n_refits"], 0, 5, "The number of stability refits")
    cfg["min_r2"] = _network_float(payload, "min_r2", d["min_r2"], 0.0, 20.0, "The minimum variance explained")
    cfg["stability_cut"] = _network_float(payload, "stability_cut", d["stability_cut"], 0.0, 1.0, "The stability cut-off")
    cfg["log2"] = bool(payload.get("log2", d["log2"]))
    cfg["scale_views"] = bool(payload.get("scale_views", d["scale_views"]))
    group_key = str(payload.get("group_column") or "")
    if group_key:
        keys = {o["key"] for o in _mofa_group_options(session)["options"]}
        if group_key not in keys:
            raise HTTPException(status_code=400, detail="That sample-group column is not available (upload the sample sheet first, or choose one of the listed columns).")
    cfg["group_column"] = group_key
    cfg["seed"] = 42
    ms = _mofa_checked(session)
    ms["config"], ms["results"], ms["enrich"] = cfg, None, {}
    return {"config": cfg, "n_samples": n}


@app.post("/api/multiomics/mofa/run")
def multiomics_mofa_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ms = _mofa_checked(session)
    if ms["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_mofa_run, session)
    return {"job_id": job_id}


def _mofa_labels(cs: dict) -> dict:
    """{omics: {feature id: display label}}: the gene symbol a transcript / protein resolved to, else the id as written."""
    out = {}
    for o in cs["selected"]:
        gene_of = (cs["resolved"].get(o) or {}).get("gene_of") or {}
        out[o] = {str(c): (gene_of.get(c) or str(c)) for c in cs["aligned"]["mats"][o].columns}
    return out


def _job_mofa_run(job_id: str, session: dict) -> dict:
    cs = _corr_state(session)
    ms = _mofa_checked(session)
    cfg = ms["config"]
    if cfg is None or cs["aligned"] is None:
        raise ValueError("Configure the analysis first.")
    gm = _mofa_group_options(session)
    groups = gm["mappings"].get(cfg["group_column"]) if cfg["group_column"] else None
    labels = _mofa_labels(cs)
    res = multiomics_mofa.run_mofa(cs["aligned"]["mats"], cfg, groups=groups, labels=labels,
                                   progress=lambda f, text: update_job(job_id, progress=float(f), message=str(text)))
    K = res["factors"].shape[1]
    names = res["factors"].columns.tolist()
    top = {}
    for o in res["omics"]:
        wdf = res["weights"][o]
        top[o] = {}
        for fn in names:
            col = wdf[fn].to_numpy()
            order = np.argsort(-np.abs(col))[:25]
            top[o][fn] = [{"feature": wdf["Feature"].iloc[i], "label": wdf["Label"].iloc[i], "weight": float(col[i])} for i in order]
    r2 = res["r2"]
    payload = {
        "config": res["config"], "samples": res["samples"], "omics": res["omics"], "warnings": res["warnings"],
        "omics_labels": {o: multiomics_correlation.OMICS_LABELS[o] for o in res["omics"]},
        "reports": res["reports"], "n_factors": K, "n_factors_requested": int(res["config"]["factors"]),
        "converged": bool(res["fit"]["converged"]), "n_iterations": res["fit"]["n_iterations"],
        "r2": {"omics": [multiomics_correlation.OMICS_LABELS[o] for o in res["omics"]], "factors": names, "values": [[float(v) for v in r2.loc[multiomics_correlation.OMICS_LABELS[o]]] for o in res["omics"]]},
        "r2_total": res["r2_total"], "stability": res["stability"],
        "factor_values": {"samples": res["samples"], "names": names, "values": [[float(v) for v in row] for row in res["factors"].to_numpy()]},
        "group_column": cfg["group_column"], "groups": ({s: groups.get(s, "") for s in res["samples"]} if groups else None),
        "group_tests": res["group_tests"], "top_weights": top,
        "n_features": {o: int(res["views"][o].shape[1]) for o in res["omics"]},
    }
    tables = {"factors": res["factors"].reset_index().rename(columns={"index": "Sample"}),
              "r2": r2.reset_index().rename(columns={"index": "Omics"}),
              "stability": pd.DataFrame(res["stability"]["per_factor"]).rename(columns={"factor": "Factor", "stability": "Stability (median |r| over refits)", "min": "Worst refit", "unstable": "Flagged unstable"}),
              "groups": pd.DataFrame(res["group_tests"]["rows"]) if res["group_tests"]["rows"] else pd.DataFrame(columns=["Factor", "Test", "p-value", "FDR", "Effect"])}
    for o in res["omics"]:
        tables[f"weights:{o}"] = res["weights"][o]
    ms["results"] = {"payload": payload, "tables": tables, "views": res["views"], "fit": res["fit"], "labels": labels}
    update_job(job_id, progress=1.0, message="Done")
    return payload


class _MofaEnrichmentUnavailable(Exception):
    pass


def _mofa_enrichment_inputs(session: dict, omics: str, job_id: "str | None" = None) -> tuple:
    """(feature keys in the pathway database's own id space, Reactome term_map) for one omics' features, built with the omics' OWN universe
    builders - the ones Joint Pathway Enrichment calls."""
    cs, ms = _corr_state(session), _mofa_state(session)
    taxon_id, taxon_name = str(cs.get("taxon_id") or ""), str(cs.get("taxon_name") or "")
    feats = [str(c) for c in ms["results"]["views"][omics].columns]
    res = cs["resolved"].get(omics) or {}
    if omics == "transcriptomics":
        gene_of, entrez_of = res.get("gene_of") or {}, res.get("entrez_of") or {}
        rows = [{"gene_symbol": gene_of[f], "entrez_id": entrez_of[f]} for f in feats if f in gene_of and f in entrez_of]
        if not rows:
            raise _MofaEnrichmentUnavailable("None of the transcripts resolved to an Entrez Gene ID, which Reactome's gene sets are keyed by.")
        universe_df = pd.DataFrame(rows).drop_duplicates("gene_symbol")
        built = multiomics_engine.build_transcriptomics_input(universe_df, set(), "Reactome", taxon_name, "")
        keys = [gene_of.get(f) for f in feats]
    elif omics == "proteomics":
        acc_of, gene_of = res.get("acc_of") or {}, res.get("gene_of") or {}
        rows = [{"accession": acc_of[f], "gene_symbol": gene_of.get(f, "")} for f in feats if f in acc_of]
        if not rows:
            raise _MofaEnrichmentUnavailable("None of the proteins resolved to a UniProt accession.")
        built = multiomics_engine.build_proteomics_input(pd.DataFrame(rows).drop_duplicates("accession"), set(), "Reactome", taxon_id, taxon_name, "")
        keys = [acc_of.get(f) for f in feats]
    else:
        if ms.get("metab_match") is None:
            if job_id:
                update_job(job_id, progress=0.1, message="Matching the compound names (once)...")
            names = pd.DataFrame({"Name": feats})
            resolved_df, matched_df, meta = _resolve_and_match_compound_names(job_id or "", names, progress_start=0.1, progress_span=0.5)
            ms["metab_match"] = {"matched": matched_df}
        matched_df = ms["metab_match"]["matched"]
        crosswalk = metabolomics_engine.build_compound_crosswalk(verbose=False)
        built = multiomics_engine.build_metabolomics_input(matched_df, crosswalk, "Reactome", taxon_id, taxon_name, "", "")
        name_to_key = {str(n): str(k) for n, k in zip(matched_df["Name"], matched_df["compound_key"]) if str(k).strip()}
        keys = [name_to_key.get(f) for f in feats]
    if built.get("note") or not built.get("term_map"):
        raise _MofaEnrichmentUnavailable(built.get("note") or "No Reactome pathway could be matched to these features.")
    return keys, built["term_map"]


@app.post("/api/multiomics/mofa/enrich")
def multiomics_mofa_enrich(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ms = _mofa_checked(session)
    if not ms["results"]:
        raise HTTPException(status_code=400, detail="Run the MOFA analysis first.")
    pl = ms["results"]["payload"]
    payload = payload or {}
    omics = str(payload.get("omics") or "")
    if omics not in pl["omics"]:
        raise HTTPException(status_code=400, detail="Choose one of the omics that took part.")
    factor = str(payload.get("factor") or "")
    if factor not in pl["factor_values"]["names"]:
        raise HTTPException(status_code=400, detail="Unknown factor.")
    mode = str(payload.get("mode") or "signed")
    if mode not in ("signed", "absolute"):
        raise HTTPException(status_code=400, detail="The mode must be 'signed' or 'absolute'.")
    min_size = _network_int(payload, "min_size", 10, 3, 200, "The minimum pathway size")
    max_size = _network_int(payload, "max_size", 500, min_size, 5000, "The maximum pathway size")
    job_id = submit_job(_job_mofa_enrich, session, omics, factor, mode, min_size, max_size)
    return {"job_id": job_id}


def _job_mofa_enrich(job_id: str, session: dict, omics: str, factor: str, mode: str, min_size: int, max_size: int) -> dict:
    ms = _mofa_state(session)
    if not ms["results"]:
        raise ValueError("Run the MOFA analysis first.")
    wdf = ms["results"]["tables"][f"weights:{omics}"]
    # the pathway sets depend on the omics only (not on the factor or the test), and building them can take minutes: built once per result
    if omics not in ms["inputs"]:
        update_job(job_id, progress=0.05, message=f"Building the Reactome pathway sets for {multiomics_correlation.OMICS_LABELS[omics]} (the first time only)...")
        try:
            ms["inputs"][omics] = _mofa_enrichment_inputs(session, omics, job_id)
        except _MofaEnrichmentUnavailable as exc:
            return {"omics": omics, "factor": factor, "mode": mode, "table": df_payload(pd.DataFrame()), "note": str(exc), "n_tested": 0}
    keys, term_map = ms["inputs"][omics]
    # the pathway test names its top members by the database key (an InChIKey, an accession): show the feature's own label instead
    key_label = {k: lab for k, lab in zip(keys, wdf["Label"].astype(str)) if k}
    w = wdf[factor].to_numpy(dtype=float)
    w, keys = multiomics_mofa.collapse_by_key(w, keys)
    update_job(job_id, progress=0.6, message="Testing the pathways...")
    table = multiomics_mofa.weight_set_enrichment(w, keys, term_map, min_size, max_size, mode)
    if not table.empty:
        table["Top members"] = [", ".join(key_label.get(k.strip(), k.strip()) for k in str(s).split(", ")) for s in table["Top members"]]
    key_name = f"{omics}:{factor}:{mode}"
    ms["enrich"][key_name] = table
    n_with_key = int(sum(1 for k in keys if k))
    note = None
    if table.empty:
        note = (f"No Reactome pathway had between {min_size} and {max_size} of the {n_with_key} features that could be matched to it, "
                "or too few features could be matched.")
    update_job(job_id, progress=1.0, message="Done")
    return {"omics": omics, "factor": factor, "mode": mode, "table": df_payload(table.head(200)), "note": note, "n_tested": int(len(table)),
            "n_features_matched": n_with_key, "n_features": int(len(wdf))}


@app.get("/api/multiomics/mofa/download")
def multiomics_mofa_download(table: str = "factors", x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ms = _mofa_checked(session)
    if not ms["results"]:
        raise HTTPException(status_code=400, detail="No MOFA result available.")
    if table.startswith("enrichment:"):
        key = table[len("enrichment:"):]
        df = ms["enrich"].get(key)
        if df is None:
            raise HTTPException(status_code=400, detail="That enrichment has not been run.")
        return df_to_csv_response(df, "multiomics_mofa_enrichment_" + re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower() + ".csv")
    tables = ms["results"]["tables"]
    if table not in tables:
        raise HTTPException(status_code=400, detail="No such result table.")
    return df_to_csv_response(tables[table], "multiomics_mofa_" + re.sub(r"[^A-Za-z0-9]+", "_", table).strip("_").lower() + ".csv")


# ===========================================================================
# Multi-omics - DIABLO / Multiblock Integration (see Multiomics_diablo)
#
# Shares the matrices of Correlation-based Integration (uploaded, resolved to gene symbols and matched on their samples there); keeps its own
# config and results under session["multiomics"]["diablo"]. A DIABLO result is tied to the exact aligned data it was fitted on (a fingerprint,
# the same idea as MOFA's): if the Correlation uploads change, the result is dropped, never shown against different data. Unlike MOFA, a sample
# group is REQUIRED (DIABLO is supervised).
# ---------------------------------------------------------------------------

def _fresh_multiomics_diablo() -> dict:
    return {"config": None, "sheet": None, "results": None, "enrich": {}, "metab_match": None, "inputs": {}}


def _diablo_state(session: dict) -> dict:
    return session["multiomics"]["diablo"]


def _diablo_fingerprint(cs: dict) -> "str | None":
    al = cs.get("aligned")
    if not al:
        return None
    parts = [",".join(al["samples"])]
    for o in sorted(al["mats"]):
        m = al["mats"][o]
        parts.append(f"{o}:{m.shape[0]}x{m.shape[1]}:{float(np.nansum(m.to_numpy(dtype=float))):.6f}:{hash(tuple(map(str, m.columns))) & 0xffffffff}")
    return "|".join(parts)


def _diablo_checked(session: dict) -> dict:
    """The DIABLO state, with any config / result / enrichment dropped when the Correlation data they were built on has changed."""
    ds = _diablo_state(session)
    fp = _diablo_fingerprint(_corr_state(session))
    if ds.get("fingerprint") != fp:
        ds.update({"config": None, "results": None, "enrich": {}, "metab_match": None, "inputs": {}, "fingerprint": fp})
        if fp is None:
            ds["sheet"] = None
    return ds


def _diablo_group_options(session: dict) -> dict:
    """The sample-group sources on offer: a label column found in a matrix upload (Group, Condition ...) or a column of the sample sheet.
    Same idea as MOFA's own group-source list, written separately here (DIABLO REQUIRES a group, MOFA's is optional)."""
    cs = _corr_state(session)
    ds = _diablo_state(session)
    samples = list((cs.get("aligned") or {}).get("samples") or [])
    opts = []
    for o in cs["selected"]:
        gc = (cs["uploads"].get(o) or {}).get("meta", {}).get("group_columns") or {}
        for col, mp in gc.items():
            opts.append({"key": f"{o}:{col}", "label": f"{multiomics_correlation.OMICS_LABELS[o]} file, column '{col}'", "mapping": mp})
    for col, mp in ((ds.get("sheet") or {}).get("columns") or {}).items():
        opts.append({"key": f"sheet:{col}", "label": f"Sample sheet, column '{col}'", "mapping": mp})
    out = []
    for op in opts:
        vals = [op["mapping"].get(s) for s in samples]
        counts: dict = {}
        for v in vals:
            if v:
                counts[v] = counts.get(v, 0) + 1
        groups = sorted(g for g, n in counts.items() if n >= multiomics_diablo.MIN_SAMPLES_PER_GROUP)
        out.append({"key": op["key"], "label": op["label"], "n_samples": sum(counts.get(g, 0) for g in groups), "groups": groups,
                    "usable": len(groups) >= 2})
    return {"options": out, "mappings": {op["key"]: op["mapping"] for op in opts}}


@app.get("/api/multiomics/diablo/options")
def multiomics_diablo_options(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    ds = _diablo_checked(session)
    go = _diablo_group_options(session)
    return {"ready": cs["aligned"] is not None, "samples": list((cs["aligned"] or {}).get("samples") or []),
            "omics": list(cs["selected"]), "group_options": go["options"],
            "sheet": ({"filename": ds["sheet"]["filename"], "n_matched": ds["sheet"]["n_matched"]} if ds.get("sheet") else None),
            "defaults": multiomics_diablo.CONFIG_DEFAULTS, "min_samples_per_group": multiomics_diablo.MIN_SAMPLES_PER_GROUP,
            "warn_samples_per_group": multiomics_diablo.WARN_SAMPLES_PER_GROUP, "max_components": multiomics_diablo.MAX_COMPONENTS,
            "distances": list(multiomics_diablo.DISTANCES)}


@app.post("/api/multiomics/diablo/groups")
async def multiomics_diablo_groups(file: UploadFile = File(...), x_session_id: Optional[str] = Header(default=None)):
    """A small sample sheet: the first column is the sample id, every other column is a candidate group label (Group, Batch, Sex ...).
    Same convention as MOFA's sample sheet, kept as DIABLO's own route / state (never shared - a change here must not touch MOFA)."""
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["aligned"] is None:
        raise HTTPException(status_code=400, detail="Upload and match the omics tables first (the sample sheet is matched to their samples).")
    raw = await file.read()
    try:
        df = pd.read_csv(io.BytesIO(raw), sep=None, engine="python", dtype=str)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Could not read the sample sheet: {exc}")
    if df.shape[1] < 2 or df.shape[0] < 2:
        raise HTTPException(status_code=400, detail="The sample sheet needs a sample-id column and at least one label column.")
    df.columns = [str(c).strip() for c in df.columns]
    ids = df.iloc[:, 0].astype(str).str.strip()
    samples = set(cs["aligned"]["samples"])
    matched = [s for s in ids if s in samples]
    if len(matched) < 2:
        raise HTTPException(status_code=400, detail=f"Only {len(matched)} of the sheet's sample ids are among the {len(samples)} samples of your tables "
                                                    f"(it must repeat the SAME ids; the first column is read as the sample id).")
    columns = {}
    for col in df.columns[1:]:
        vals = df[col].astype(str).str.strip()
        mp = {s: v for s, v in zip(ids, vals) if s in samples and v not in ("", "nan")}
        distinct = set(mp.values())
        if 2 <= len(distinct) <= max(2, len(mp) // 2):
            columns[col] = mp
    if not columns:
        raise HTTPException(status_code=400, detail="No column of the sample sheet looks like a group label (text with a few repeated values).")
    ds = _diablo_checked(session)
    ds["sheet"] = {"filename": getattr(file, "filename", ""), "columns": columns, "n_matched": len(set(matched))}
    ds["results"] = None
    ds["enrich"] = {}
    return {"filename": ds["sheet"]["filename"], "n_matched": len(set(matched)), "n_samples": len(samples), "columns": list(columns),
            "group_options": _diablo_group_options(session)["options"]}


@app.post("/api/multiomics/diablo/configure")
def multiomics_diablo_configure(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    cs = _corr_state(session)
    if cs["aligned"] is None:
        raise HTTPException(status_code=400, detail="Upload, resolve and match the omics tables first.")
    payload = payload or {}
    d = multiomics_diablo.CONFIG_DEFAULTS
    cfg = dict(d)
    group_key = str(payload.get("group_column") or "")
    go = _diablo_group_options(session)
    opt = next((o for o in go["options"] if o["key"] == group_key), None)
    if not opt:
        raise HTTPException(status_code=400, detail="Choose a sample-group column (DIABLO needs groups to classify): upload a sample sheet first, "
                                                    "or choose one of the listed columns.")
    if not opt["usable"]:
        raise HTTPException(status_code=400, detail=f"'{opt['label']}' has fewer than two groups with at least {multiomics_diablo.MIN_SAMPLES_PER_GROUP} "
                                                    "samples each; choose another column.")
    cfg["group_column"] = group_key
    cfg["top_variance"] = _network_int(payload, "top_variance", d["top_variance"], 50, 5000, "The number of most variable features")
    cfg["log2"] = bool(payload.get("log2", d["log2"]))
    cfg["design_value"] = _network_float(payload, "design_value", d["design_value"], 0.0, 1.0, "The design value between omics")
    cfg["max_components"] = _network_int(payload, "max_components", d["max_components"], 1, multiomics_diablo.MAX_COMPONENTS, "The maximum number of components")
    grid = payload.get("keepX_grid")
    if grid is not None:
        try:
            grid = sorted({int(v) for v in grid})
        except (TypeError, ValueError):
            raise HTTPException(status_code=400, detail="'keepX_grid' must be a list of whole numbers.")
        if not (1 <= len(grid) <= 8) or any(v < 2 or v > 500 for v in grid):
            raise HTTPException(status_code=400, detail="The keepX grid must have 1 to 8 values, each between 2 and 500.")
        cfg["keepX_grid"] = grid
    cfg["n_folds"] = _network_int(payload, "n_folds", d["n_folds"], 2, 10, "The number of cross-validation folds")
    cfg["n_repeats"] = _network_int(payload, "n_repeats", d["n_repeats"], 1, 50, "The number of cross-validation repeats")
    distance = str(payload.get("distance", d["distance"]))
    if distance not in multiomics_diablo.DISTANCES:
        raise HTTPException(status_code=400, detail=f"'distance' must be one of {list(multiomics_diablo.DISTANCES)}.")
    cfg["distance"] = distance
    cfg["seed"] = 42
    n = len(cs["aligned"]["samples"])
    mapping = go["mappings"][group_key]
    smallest_group = min(sum(1 for v in mapping.values() if v == g) for g in opt["groups"])
    if cfg["n_folds"] > smallest_group:
        raise HTTPException(status_code=400, detail=f"The smallest group has {smallest_group} samples, fewer than the {cfg['n_folds']} "
                                                    "cross-validation folds asked for; choose fewer folds.")
    ds = _diablo_checked(session)
    ds["config"], ds["results"], ds["enrich"] = cfg, None, {}
    return {"config": cfg, "n_samples": n, "groups": opt["groups"]}


@app.post("/api/multiomics/diablo/run")
def multiomics_diablo_run(x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ds = _diablo_checked(session)
    if ds["config"] is None:
        raise HTTPException(status_code=400, detail="Configure the analysis first.")
    job_id = submit_job(_job_diablo_run, session)
    return {"job_id": job_id}


def _diablo_labels(cs: dict) -> dict:
    """{omics: {feature id: display label}}: the gene symbol a transcript / protein resolved to, else the id as written."""
    out = {}
    for o in cs["selected"]:
        gene_of = (cs["resolved"].get(o) or {}).get("gene_of") or {}
        out[o] = {str(c): (gene_of.get(c) or str(c)) for c in cs["aligned"]["mats"][o].columns}
    return out


def _job_diablo_run(job_id: str, session: dict) -> dict:
    cs = _corr_state(session)
    ds = _diablo_checked(session)
    cfg = ds["config"]
    if cfg is None or cs["aligned"] is None:
        raise ValueError("Configure the analysis first.")
    gm = _diablo_group_options(session)
    mapping = gm["mappings"][cfg["group_column"]]
    samples = cs["aligned"]["samples"]
    groups = np.array([mapping.get(s, "") for s in samples], dtype=object)
    labels = multiomics_diablo.check_groups(groups)
    keep_mask = np.isin(groups, labels)
    groups_used = groups[keep_mask]
    mats_used = {o: cs["aligned"]["mats"][o].loc[[s for s, k in zip(samples, keep_mask) if k]] for o in cs["aligned"]["mats"]}
    samples_used = [s for s, k in zip(samples, keep_mask) if k]
    n_dropped = len(samples) - len(samples_used)

    update_job(job_id, progress=0.03, message="Preparing the data...")
    blocks, reports = multiomics_diablo.prepare_blocks(mats_used, cfg)

    def progress_span(lo, hi):
        return lambda f: update_job(job_id, progress=lo + (hi - lo) * max(0.0, min(1.0, f)))

    update_job(job_id, progress=0.05, message="Choosing how many components to keep (cross-validated)...")
    tc = multiomics_diablo.tune_components(blocks, groups_used, labels, cfg, progress=progress_span(0.05, 0.3))
    ncomp = tc["chosen"]
    update_job(job_id, progress=0.3, message=f"Choosing which features to keep per omics (cross-validated, {ncomp} component(s))...")
    tk = multiomics_diablo.tune_keepX(blocks, groups_used, labels, cfg, ncomp, progress=progress_span(0.3, 0.75))
    update_job(job_id, progress=0.78, message="Measuring the cross-validated error with the chosen settings...")
    cv = multiomics_diablo.cv_ber(blocks, groups_used, labels, cfg, ncomp, tk["keepX"], cfg["n_folds"], cfg["n_repeats"], cfg["seed"], cfg["distance"],
                                  progress=progress_span(0.78, 0.86))
    update_job(job_id, progress=0.87, message="Checking the error against random group labels (permutation test)...")
    perm = multiomics_diablo.permutation_test_ber(blocks, groups_used, labels, cfg, ncomp, tk["keepX"], cv["BER"], n_perm=20, seed=cfg["seed"],
                                                  progress=progress_span(0.87, 0.97))
    update_job(job_id, progress=0.98, message="Fitting the final model on every sample...")
    fit = multiomics_diablo.fit_diablo(blocks, groups_used, labels, cfg, ncomp, keepX=tk["keepX"])
    sel = multiomics_diablo.selected_features(fit)
    ave_raw = multiomics_diablo.explained_variance(fit)
    ave = {o: [None if (v != v) else v for v in vals] for o, vals in ave_raw.items()}   # NaN -> null: strict JSON has no NaN
    labels_display = _diablo_labels(cs)

    warnings = []
    if n_dropped:
        warnings.append(f"{n_dropped} sample(s) had no usable group label and were left out of this run.")
    smallest = min(int((groups_used == g).sum()) for g in labels)
    if smallest < multiomics_diablo.WARN_SAMPLES_PER_GROUP:
        warnings.append(f"The smallest group has only {smallest} samples: below {multiomics_diablo.WARN_SAMPLES_PER_GROUP}, the cross-validated "
                        "error and the selected features are both unstable - treat them as exploratory.")

    names = [f"Component{i + 1}" for i in range(ncomp)]
    top = {}
    for o in fit["omics"]:
        w = fit["loadings"][o]
        top[o] = {}
        for h, cname in enumerate(names):
            col = w[:, h]
            order = np.nonzero(col != 0)[0]
            order = order[np.argsort(-np.abs(col[order]))]
            top[o][cname] = [{"feature": fit["feature_names"][o][i], "label": labels_display.get(o, {}).get(str(fit["feature_names"][o][i]), str(fit["feature_names"][o][i])),
                              "weight": float(col[i])} for i in order]
    payload = {
        "config": cfg, "samples": fit["samples"], "omics": fit["omics"], "warnings": warnings,
        "omics_labels": {o: multiomics_correlation.OMICS_LABELS[o] for o in fit["omics"]},
        "n_components_offered": cfg["max_components"], "n_components": ncomp, "component_search": tc["per_ncomp"],
        "component_names": names, "keepX": tk["keepX"], "keepX_search": tk["per_component"],
        "group_column": cfg["group_column"], "group_labels": labels, "groups": {s: str(g) for s, g in zip(fit["samples"], groups_used)},
        "cv": {"BER": cv["BER"], "BER_sd": cv["BER_sd"], "ER": cv["ER"], "ER_sd": cv["ER_sd"], "per_repeat": cv["per_repeat"],
              "n_folds": cfg["n_folds"], "n_repeats": cfg["n_repeats"], "distance": cfg["distance"]},
        "permutation": {"p_value": perm["p_value"], "null_bers": perm["null_bers"], "n_perm": perm["n_perm"]},
        "selected_features": sel, "n_selected": {o: len(sel[o]) for o in sel},
        "n_features": {o: int(blocks[o].shape[1]) for o in blocks}, "n_features_in": {o: int(reports[o]["n_in"]) for o in reports},
        "ave": ave, "top_weights": top,
        "component_values": {o: fit["variates"][o].tolist() for o in fit["omics"]},
        "design": fit["design"].tolist(),
    }
    tables = {
        "component_search": pd.DataFrame(tc["per_ncomp"]).rename(columns={"ncomp": "Components", "BER": "Cross-validated error (BER)", "BER_sd": "+/- s.d."}),
        "cv": pd.DataFrame([{"Metric": "Balanced error rate (BER)", "Value": cv["BER"], "+/- s.d.": cv["BER_sd"]},
                           {"Metric": "Overall error rate", "Value": cv["ER"], "+/- s.d.": cv["ER_sd"]}]),
    }
    for o in fit["omics"]:
        wdf = pd.DataFrame(fit["loadings"][o], columns=names)
        wdf.insert(0, "Label", [labels_display.get(o, {}).get(str(c), str(c)) for c in fit["feature_names"][o]])
        wdf.insert(0, "Feature", [str(c) for c in fit["feature_names"][o]])
        tables[f"weights:{o}"] = wdf
        tdf = pd.DataFrame(fit["variates"][o], index=fit["samples"], columns=names).reset_index().rename(columns={"index": "Sample"})
        tables[f"components:{o}"] = tdf
    ds["results"] = {"payload": payload, "tables": tables, "fit": fit, "blocks": blocks, "labels": labels_display}
    update_job(job_id, progress=1.0, message="Done")
    return payload


class _DiabloEnrichmentUnavailable(Exception):
    pass


def _diablo_enrichment_inputs(session: dict, omics: str, job_id: "str | None" = None) -> tuple:
    """(feature keys in the pathway database's own id space, Reactome term_map, background keys = every feature that
    entered the model for this omics) - written separately from MOFA's own version (same idea, no shared code - a
    change to one tile's enrichment inputs must not silently change the other's)."""
    cs, ds = _corr_state(session), _diablo_state(session)
    taxon_id, taxon_name = str(cs.get("taxon_id") or ""), str(cs.get("taxon_name") or "")
    feats = [str(c) for c in ds["results"]["blocks"][omics].columns]
    res = cs["resolved"].get(omics) or {}
    if omics == "transcriptomics":
        gene_of, entrez_of = res.get("gene_of") or {}, res.get("entrez_of") or {}
        rows = [{"gene_symbol": gene_of[f], "entrez_id": entrez_of[f]} for f in feats if f in gene_of and f in entrez_of]
        if not rows:
            raise _DiabloEnrichmentUnavailable("None of the transcripts resolved to an Entrez Gene ID, which Reactome's gene sets are keyed by.")
        universe_df = pd.DataFrame(rows).drop_duplicates("gene_symbol")
        built = multiomics_engine.build_transcriptomics_input(universe_df, set(), "Reactome", taxon_name, "")
        keys = [gene_of.get(f) for f in feats]
    elif omics == "proteomics":
        acc_of, gene_of = res.get("acc_of") or {}, res.get("gene_of") or {}
        rows = [{"accession": acc_of[f], "gene_symbol": gene_of.get(f, "")} for f in feats if f in acc_of]
        if not rows:
            raise _DiabloEnrichmentUnavailable("None of the proteins resolved to a UniProt accession.")
        built = multiomics_engine.build_proteomics_input(pd.DataFrame(rows).drop_duplicates("accession"), set(), "Reactome", taxon_id, taxon_name, "")
        keys = [acc_of.get(f) for f in feats]
    else:
        if ds.get("metab_match") is None:
            if job_id:
                update_job(job_id, progress=0.1, message="Matching the compound names (once)...")
            names = pd.DataFrame({"Name": feats})
            resolved_df, matched_df, meta = _resolve_and_match_compound_names(job_id or "", names, progress_start=0.1, progress_span=0.5)
            ds["metab_match"] = {"matched": matched_df}
        matched_df = ds["metab_match"]["matched"]
        crosswalk = metabolomics_engine.build_compound_crosswalk(verbose=False)
        built = multiomics_engine.build_metabolomics_input(matched_df, crosswalk, "Reactome", taxon_id, taxon_name, "", "")
        name_to_key = {str(n): str(k) for n, k in zip(matched_df["Name"], matched_df["compound_key"]) if str(k).strip()}
        keys = [name_to_key.get(f) for f in feats]
    if built.get("note") or not built.get("term_map"):
        raise _DiabloEnrichmentUnavailable(built.get("note") or "No Reactome pathway could be matched to these features.")
    return keys, built["term_map"], feats


@app.post("/api/multiomics/diablo/enrich")
def multiomics_diablo_enrich(payload: dict, x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ds = _diablo_checked(session)
    if not ds["results"]:
        raise HTTPException(status_code=400, detail="Run the DIABLO analysis first.")
    pl = ds["results"]["payload"]
    payload = payload or {}
    omics = str(payload.get("omics") or "")
    if omics not in pl["omics"]:
        raise HTTPException(status_code=400, detail="Choose one of the omics that took part.")
    min_size = _network_int(payload, "min_size", 10, 3, 200, "The minimum pathway size")
    max_size = _network_int(payload, "max_size", 500, min_size, 5000, "The maximum pathway size")
    job_id = submit_job(_job_diablo_enrich, session, omics, min_size, max_size)
    return {"job_id": job_id}


def _job_diablo_enrich(job_id: str, session: dict, omics: str, min_size: int, max_size: int) -> dict:
    ds = _diablo_state(session)
    if not ds["results"]:
        raise ValueError("Run the DIABLO analysis first.")
    update_job(job_id, progress=0.05, message=f"Building the Reactome pathway sets for {multiomics_correlation.OMICS_LABELS[omics]}...")
    try:
        if omics not in ds["inputs"]:
            ds["inputs"][omics] = _diablo_enrichment_inputs(session, omics, job_id)
        keys, term_map, feats = ds["inputs"][omics]
    except _DiabloEnrichmentUnavailable as exc:
        return {"omics": omics, "table": df_payload(pd.DataFrame()), "note": str(exc), "n_tested": 0}
    sel_feats = set(ds["results"]["payload"]["selected_features"][omics])
    key_of = dict(zip(feats, keys))
    selected_keys = {key_of[f] for f in sel_feats if key_of.get(f)}
    background_keys = {k for k in keys if k}
    update_job(job_id, progress=0.6, message="Testing the pathways...")
    table = multiomics_diablo.selected_feature_enrichment(selected_keys, background_keys, term_map, min_size, max_size)
    ds["enrich"][omics] = table
    n_matched = len(selected_keys)
    note = None
    if table.empty:
        note = (f"No Reactome pathway had between {min_size} and {max_size} of the {len(background_keys)} background features, "
                f"or too few of the {len(sel_feats)} selected features could be matched (matched: {n_matched}).")
    update_job(job_id, progress=1.0, message="Done")
    return {"omics": omics, "table": df_payload(table.head(200)), "note": note, "n_tested": int(len(table)),
            "n_selected_matched": n_matched, "n_selected": len(sel_feats)}


@app.get("/api/multiomics/diablo/download")
def multiomics_diablo_download(table: str = "cv", x_session_id: Optional[str] = Header(default=None)):
    session = session_dep(x_session_id)
    ds = _diablo_checked(session)
    if not ds["results"]:
        raise HTTPException(status_code=400, detail="No DIABLO result available.")
    if table.startswith("enrichment:"):
        key = table[len("enrichment:"):]
        df = ds["enrich"].get(key)
        if df is None:
            raise HTTPException(status_code=400, detail="That enrichment has not been run.")
        return df_to_csv_response(df, "multiomics_diablo_enrichment_" + re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_").lower() + ".csv")
    tables = ds["results"]["tables"]
    if table not in tables:
        raise HTTPException(status_code=400, detail="No such result table.")
    return df_to_csv_response(tables[table], "multiomics_diablo_" + re.sub(r"[^A-Za-z0-9]+", "_", table).strip("_").lower() + ".csv")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("API:app", host="127.0.0.1", port=8000, reload=True)
