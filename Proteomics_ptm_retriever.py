# -*- coding: utf-8 -*-
"""
Proteomics PTM/phosphosite retriever for Omiscera.

Turns an uploaded SITE-level differential-abundance table (one row per
modified residue, not one row per protein) into a clean, analysis-ready
table for the PTM / Phosphosite Enrichment tile. This is a SEPARATE pipeline
from Proteomics_module_retriever.py, not an extension of it - that pipeline
is irreducibly protein-level (RESOLVED_COLUMNS has no residue/position, and
_dedupe_by_accession keeps exactly one row per accession), so a site-level
input needs its own contract end to end.

THE ONE DECISION THAT MATTERS MOST HERE: isoform suffixes are NEVER
stripped from `accession`. Proteomics_module_retriever.strip_isoform exists
because GO/Reactome/InterPro/Complex Portal are annotated on the canonical
entry - but a phosphosite's POSITION is isoform-specific ("P12345-2" S100 is
not the same residue as "P12345" S100 unless the isoform happens to share
that exact numbering, which is not guaranteed). Stripping here would
silently relabel sites onto the wrong sequence. `canonical_accession`
(isoform-stripped) is computed alongside `accession` and used ONLY for
protein-level joins (symbol/name lookup, kinase-substrate matching against
iPTMnet's own canonical-only accessions) - never for site identity.

Reuses Proteomics_module_retriever's CSV loading, protein-group parsing,
contaminant filtering, and UniProt ID Mapping resolution directly (cross-
module private imports are already this codebase's convention - see that
module's own reuse of Transcriptomics_module_retriever).
"""

import os
import re

import numpy as np
import pandas as pd

from Transcriptomics_module_retriever import (
    _parse_de_csv_text, _coerce_numeric, _benjamini_hochberg, _normalize_col_key,
)
from Proteomics_module_retriever import (
    parse_protein_group, strip_isoform, filter_contaminants_decoys,
    _detect_flag_columns, _idmapping_resolve,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESOLVED_SITE_COLUMNS = [
    "input_id", "accession", "canonical_accession", "isoform", "protein_group",
    "ambiguous_group", "gene_symbol", "protein_name", "residue", "position",
    "site_id", "canonical_site_id", "multiplicity", "localization_prob",
    "log2FC", "pvalue", "padj", "status",
]
# site_id is THIS ROW'S TRUE IDENTITY - built from `accession` (isoform kept),
# so it is the row-level dedup/lookup key everywhere in THIS app.
# canonical_site_id is built from `canonical_accession` and exists ONLY for
# matching against EXTERNAL databases that publish canonical-only
# accessions (iPTMnet, dbPTM) - never used as an identity key here. Conflating
# the two would silently merge two different isoforms' sites onto one
# identity, exactly the mistake this module's docstring says not to make.

LOCALIZATION_PROB_DEFAULT = 0.75  # Class I, the standard phosphoproteomics convention

# ---------------------------------------------------------------------------
# Column detection - site-specific identifiers on top of the log2FC/pvalue/
# padj synonyms Proteomics_module_retriever.COLUMN_SYNONYMS already covers
# (duplicated here rather than imported, since the identifier synonym set
# genuinely differs: "site" and "position" have no protein-level meaning).
# ---------------------------------------------------------------------------

COLUMN_SYNONYMS = {
    "identifier": {
        "protein", "proteinid", "protein_id", "proteinids", "protein_ids",
        "majorityproteinids", "majority_protein_ids", "leadingprotein",
        "leading_protein", "accession", "uniprot", "uniprotid", "uniprot_id",
        "uniprotaccession", "uniprot_accession", "proteingroup", "protein_group",
        "proteingroups", "protein_groups", "proteins",
    },
    # A single combined site token, e.g. "S473", "Ser473", "pS473" - OR, if
    # this column is absent, "residue" + "position" together take its place
    # (detect_columns accepts either shape).
    "site": {
        "site", "siteid", "site_id", "modifiedsite", "modified_site",
        "ptmsite", "ptm_site", "phosphosite", "modifiedresidue", "modified_residue",
    },
    "residue": {
        "residue", "aminoacid", "amino_acid", "aa", "modifiedaminoacid",
    },
    "position": {
        "position", "positionwithinproteins", "position_within_proteins",
        "positionswithinproteins", "positions_within_proteins", "sitepos",
        "sitepositon", "pos", "modifiedresiduenumber", "residueposition",
    },
    "log2FC": {
        "log2foldchange", "log2fc", "logfc", "lfc", "log2_fold_change",
        "log2fold_change", "coef", "beta",
    },
    "pvalue": {"pvalue", "pval", "p_value", "p", "pvalraw"},
    "padj": {
        "padj", "fdr", "qvalue", "q_value", "adjpvalue", "adj_p_val",
        "adjpval", "p_adj", "bhfdr", "fdrq",
    },
    # MaxQuant emits one row per Class (___1/___2/___3) of a multiply-
    # phosphorylated peptide; see _dedupe_by_site for why these must never
    # be averaged together.
    "multiplicity": {"multiplicity"},
    "localization_prob": {
        "localizationprob", "localization_prob", "localizationprobability",
        "localization_probability", "ptmscore", "ptm_score", "score",
    },
    # Optional, QC/display only - not required for any statistic here, but
    # kept for a future dbPTM residue-composition analysis (Stage 3).
    "sequence_window": {
        "sequencewindow", "sequence_window", "sequence", "peptidewindow",
    },
}

_CONTAMINANT_PREFIX_RE = re.compile(r"^(CON__|REV__)", re.IGNORECASE)

# Three-letter -> one-letter amino acid codes, covering every residue any
# PTM type in dbPTM/iPTMnet modifies (not just S/T/Y).
_AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
_VALID_RESIDUES_1LETTER = set(_AA_3TO1.values())
# A combined site token: optional leading 'p' (phospho-notation, "pS473"),
# then a 1- or 3-letter residue code, then digits. Separator (space/_/-)
# between residue and digits is optional and tolerated.
_SITE_TOKEN_RE = re.compile(r"^[pP]?([A-Za-z]{1,3})[\s_-]?(\d+)$")


def normalize_site(raw_site: "str | None" = None, residue: "str | None" = None,
                   position=None) -> "tuple[str, int] | tuple[None, None]":
    """Parse a site into (residue_1letter, position_int), or (None, None) if
    it cannot be parsed. Accepts EITHER a single combined token (raw_site,
    e.g. "S473"/"Ser473"/"pS473") OR separate residue+position (MaxQuant's
    own "Amino acid" + "Positions within proteins" columns - the latter is
    semicolon-aligned with the protein group; only the FIRST position is
    used here, matching the leading-accession convention this codebase
    already applies to protein groups elsewhere)."""
    if raw_site is not None and str(raw_site).strip():
        m = _SITE_TOKEN_RE.match(str(raw_site).strip())
        if not m:
            return None, None
        letters, digits = m.group(1).upper(), m.group(2)
    else:
        if residue is None or position is None:
            return None, None
        letters = str(residue).strip().upper()
        digits = str(position).strip().split(";")[0].strip()
        if not letters or not digits.isdigit():
            return None, None

    if len(letters) == 1:
        res1 = letters
    elif letters in _AA_3TO1:
        res1 = _AA_3TO1[letters]
    else:
        return None, None
    if res1 not in _VALID_RESIDUES_1LETTER or not digits.isdigit():
        return None, None
    return res1, int(digits)


def detect_columns(df: pd.DataFrame) -> dict:
    """Map arbitrary uploaded headers to canonical identifier/site-or-
    residue+position/log2FC/pvalue[/padj/multiplicity/localization_prob].
    Raises ValueError listing exactly what is missing - a site table has
    more required-shape variants than the protein-level table (site XOR
    residue+position), so the message is written to be actionable."""
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
    has_site_shape = "site" in mapping or ("residue" in mapping and "position" in mapping)
    if not has_site_shape:
        missing.append("site (or residue + position)")
    if missing:
        raise ValueError(
            f"Phosphosite table is missing required column(s) for: {missing}. "
            f"Detected columns after normalisation: {list(df.columns)}"
        )
    return mapping


# ---------------------------------------------------------------------------
# CSV loading - reused verbatim (identifier-agnostic multi-delimiter sniff)
# ---------------------------------------------------------------------------

def load_site_table_from_bytes(content: bytes) -> pd.DataFrame:
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = content.decode("utf-8", errors="replace")
    return _parse_de_csv_text(text)


def load_site_table(path: str) -> pd.DataFrame:
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin1") as fh:
            text = fh.read()
    return _parse_de_csv_text(text)


# ---------------------------------------------------------------------------
# Multiplicity collapse - a MaxQuant "Phospho (STY)Sites.txt" export can
# carry ___1/___2/___3 rows for the SAME site (different co-eluting
# multiply-phosphorylated peptide classes). These are collapsed to one row
# per site_id, but NEVER by averaging log2FC - they are different analytes
# (a diphospho peptide's quantification is not interchangeable with a
# monophospho one at the same site) - the smallest-pvalue row is kept
# instead, mirroring Proteomics_module_retriever._dedupe_by_accession's own
# "most significant row wins" convention exactly.
# ---------------------------------------------------------------------------

def _dedupe_by_site(df: pd.DataFrame) -> "tuple[pd.DataFrame, int]":
    if "site_id" not in df.columns or df.empty:
        return df.copy(), 0
    has_site = df["site_id"].astype(str).str.strip() != ""
    mapped = df[has_site].copy()
    unmapped = df[~has_site].copy()
    if mapped.empty:
        return df.copy(), 0

    def _sort_key(row):
        try:
            p = float(row.get("pvalue"))
            if pd.isna(p):
                p = float("inf")
        except (TypeError, ValueError):
            p = float("inf")
        try:
            abs_lfc = abs(float(row.get("log2FC", 0.0)))
        except (TypeError, ValueError):
            abs_lfc = 0.0
        try:
            mult = float(row.get("multiplicity") or 0.0)
        except (TypeError, ValueError):
            mult = 0.0
        return (p, -abs_lfc, mult)

    mapped["_sort_key"] = mapped.apply(_sort_key, axis=1)
    mapped = mapped.sort_values("_sort_key", kind="stable")
    grouped = mapped.groupby(mapped["site_id"], sort=False)

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

def resolve_site_table(
    source,
    taxon_id: str,
    localization_threshold: "float | None" = LOCALIZATION_PROB_DEFAULT,
    drop_ambiguous_groups: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> "tuple[pd.DataFrame, dict]":
    """Resolve an uploaded site-level phosphoproteomics table into
    RESOLVED_SITE_COLUMNS.

    localization_threshold: rows with a localization_prob column present and
    a value BELOW this are dropped (Class I convention, default 0.75);
    None disables the filter. Rows with no localization_prob column at all
    are never dropped for it - `meta["has_localization_prob"]` reports
    whether the filter could even apply, so the caller can show an honest
    caveat instead of a silently-empty filter.

    drop_ambiguous_groups: when True, rows whose protein group has more than
    one member (position validity is only guaranteed for the leading one)
    are dropped instead of kept-with-a-flag; default False, matching
    Proteomics_module_retriever's own "keep and flag" convention.

    Returns (resolved_df, meta). resolved_df has RESOLVED_SITE_COLUMNS.
    """
    if isinstance(source, pd.DataFrame):
        raw_df = source.copy()
    elif isinstance(source, (str, os.PathLike)):
        raw_df = load_site_table(str(source))
    else:
        raise TypeError("source must be a pandas DataFrame or a file path")

    n_input_rows = len(raw_df)
    column_mapping = detect_columns(raw_df)
    has_localization_prob = "localization_prob" in column_mapping

    work = pd.DataFrame({
        "identifier": raw_df[column_mapping["identifier"]].astype(str).str.strip(),
        "log2FC": raw_df[column_mapping["log2FC"]],
        "pvalue": raw_df[column_mapping["pvalue"]],
    })
    work["padj"] = raw_df[column_mapping["padj"]] if "padj" in column_mapping else np.nan
    work["multiplicity"] = raw_df[column_mapping["multiplicity"]] if "multiplicity" in column_mapping else np.nan
    work["localization_prob"] = (
        pd.to_numeric(raw_df[column_mapping["localization_prob"]], errors="coerce")
        if has_localization_prob else np.nan
    )
    if "site" in column_mapping:
        work["_raw_site"] = raw_df[column_mapping["site"]].astype(str)
        work["_raw_residue"] = ""
        work["_raw_position"] = ""
    else:
        work["_raw_site"] = ""
        work["_raw_residue"] = raw_df[column_mapping["residue"]].astype(str)
        work["_raw_position"] = raw_df[column_mapping["position"]].astype(str)

    # Contaminant/decoy filtering, same two mechanisms and same "before
    # group-parsing" ordering as the protein-level pipeline.
    flag_cols = _detect_flag_columns(raw_df)
    for col in flag_cols:
        work[col] = raw_df[col]
    work = work[work["identifier"] != ""].reset_index(drop=True)
    work, contam_counts = filter_contaminants_decoys(work, "identifier")
    work = work.drop(columns=flag_cols, errors="ignore")

    work["protein_group"] = work["identifier"]
    work["accession"] = work["identifier"].apply(parse_protein_group)  # leading member, ISOFORM KEPT
    work["canonical_accession"] = work["accession"].apply(strip_isoform)
    work["isoform"] = work["accession"].apply(lambda a: a.split("-", 1)[1] if "-" in a else "")
    work["_n_group_members"] = work["protein_group"].apply(lambda v: len([m for m in str(v).split(";") if m.strip()]))
    work = work[work["accession"] != ""].reset_index(drop=True)

    parsed = work.apply(
        lambda r: normalize_site(raw_site=r["_raw_site"] or None,
                                 residue=r["_raw_residue"] or None,
                                 position=r["_raw_position"] or None),
        axis=1,
    )
    work["residue"] = [p[0] for p in parsed]
    work["position"] = [p[1] for p in parsed]
    n_unparsed_site = int(work["residue"].isna().sum())
    work = work[work["residue"].notna()].reset_index(drop=True)
    work["position"] = work["position"].astype(int)
    # site_id: this row's TRUE identity, isoform-preserved (see
    # RESOLVED_SITE_COLUMNS' own comment for why this must never be built
    # from canonical_accession). canonical_site_id is the SEPARATE key used
    # only when matching against external canonical-only databases.
    work["site_id"] = work["accession"] + "_" + work["residue"] + work["position"].astype(str)
    work["canonical_site_id"] = work["canonical_accession"] + "_" + work["residue"] + work["position"].astype(str)

    work, n_dropped_nonnumeric = _coerce_numeric(work)

    padj_computed = False
    if "padj" not in column_mapping or work["padj"].isna().all():
        work["padj"] = _benjamini_hochberg(work["pvalue"])
        padj_computed = True

    work["ambiguous_group"] = work["_n_group_members"] > 1
    n_group_ambiguous = int(work["ambiguous_group"].sum())
    if drop_ambiguous_groups:
        work = work[~work["ambiguous_group"]].reset_index(drop=True)

    n_before_loc = len(work)
    if localization_threshold is not None and has_localization_prob:
        keep = work["localization_prob"].isna() | (work["localization_prob"] >= localization_threshold)
        work = work[keep].reset_index(drop=True)
    n_dropped_localization = n_before_loc - len(work)

    # Symbol/name resolution: ONE UniProt ID Mapping call per unique
    # CANONICAL accession (never per site) - reuses
    # Proteomics_module_retriever's own resolver directly rather than
    # re-implementing ID Mapping.
    unique_canonical = sorted({a.upper() for a in work["canonical_accession"] if a})
    resolved_map, n_fetch_failures = _idmapping_resolve(unique_canonical, verbose=verbose)
    if progress_callback:
        try:
            progress_callback(len(unique_canonical), len(unique_canonical), "resolved via UniProt ID Mapping")
        except Exception:
            pass
    if unique_canonical and n_fetch_failures and not resolved_map:
        raise RuntimeError(
            "UniProt ID Mapping failed for all proteins after retries - this looks like a "
            "network/connectivity problem, not proteins that don't exist. Please try again."
        )

    accessions_resolved, symbols, names, statuses = [], [], [], []
    for canon in work["canonical_accession"]:
        r = resolved_map.get(canon.upper())
        if r:
            accessions_resolved.append(r["accession"])
            symbols.append(r["gene_symbol"])
            names.append(r["protein_name"])
            statuses.append("mapped" if r["accession"].upper() == canon.upper() else "resolved_secondary")
        else:
            accessions_resolved.append("")
            symbols.append("")
            names.append("")
            statuses.append("unmapped")
    work["gene_symbol"] = symbols
    work["protein_name"] = names
    work["status"] = statuses

    work = work.rename(columns={"identifier": "input_id"})
    work, n_duplicates_collapsed = _dedupe_by_site(work)

    n_mapped = int((work["status"] != "unmapped").sum())
    n_unmapped = int((work["status"] == "unmapped").sum())
    n_resolved_secondary = int((work["status"] == "resolved_secondary").sum())

    resolved_df = work[RESOLVED_SITE_COLUMNS].reset_index(drop=True)

    meta = {
        "taxon": str(taxon_id),
        "column_mapping": column_mapping,
        "n_input_rows": n_input_rows,
        "n_mapped": n_mapped,
        "n_unmapped": n_unmapped,
        "n_resolved_secondary": n_resolved_secondary,
        "n_duplicates_collapsed": n_duplicates_collapsed,
        "n_dropped_nonnumeric": n_dropped_nonnumeric,
        "n_unparsed_site": n_unparsed_site,
        "n_group_ambiguous": n_group_ambiguous,
        "dropped_ambiguous_groups": bool(drop_ambiguous_groups),
        "has_localization_prob": has_localization_prob,
        "localization_threshold": localization_threshold if has_localization_prob else None,
        "n_dropped_localization": n_dropped_localization,
        "padj_computed": padj_computed,
        "n_fetch_failures": n_fetch_failures,
        "contaminants": contam_counts,
    }
    return resolved_df, meta
