# -*- coding: utf-8 -*-
"""
Metabolomics statistical analysis module for Omiscera.

Hosts Over-Representation Analysis (ORA) and Quantitative Enrichment for
Metabolomics. Built on Genomics_module's engines (_run_fisher_for_namespace,
_run_wilcoxon_for_namespace) exactly as Proteomics_module is - those engines
are genuinely identifier-agnostic (plain sets/dicts of strings), so no fork
of any statistics is needed here. This file adds only what is genuinely new:
the metabolite identifier-bridging layer and the metabolite-specific
databases.

The identifier problem (why build_compound_crosswalk exists)
------------------------------------------------------------
Metabolomics_retriever resolves a compound NAME to PubChem CID, InChIKey,
KEGG ID, HMDB ID and SMILES (via PubChem + KEGG). But the analysis databases
are keyed differently, and two of them cannot be reached from that set at
all:

    KEGG pathways   KEGG C-number      retriever resolves it directly
    SMPDB           KEGG or ChEBI      KEGG direct; ChEBI via crosswalk
    PathBank        KEGG or ChEBI      KEGG direct; ChEBI via crosswalk
    Reactome        ChEBI ONLY         crosswalk required
    RefMet class    name/CID/InChIKey/KEGG/HMDB
    LIPID MAPS      LM_ID ONLY         crosswalk required

So this module builds a local crosswalk from bulk files and joins everything
through ONE internal identifier (see _compound_key): the full InChIKey when
available, else a namespaced fallback. Building a local metabolite ID
crosswalk is standard practice for this problem (cf. Bioconductor's
metaboliteIDmapping); the point of doing it from bulk files is to avoid
thousands of per-compound API calls at run time.

Every crosswalk source is OPTIONAL and independently degradable. If a host
is unreachable the build continues with the remaining sources and names the
missing one (and what capability it costs) in its caveat - never a silent
partial result. This is not defensive padding: ftp.ebi.ac.uk (ChEBI's only
bulk host) was verified completely unreachable on 2026-09-21 from a machine
where www.ebi.ac.uk answered normally, and the same host serves files the
Proteomics module already depends on.

Databases and licences (all verified live 2026-09-21)
-----------------------------------------------------
  RefMet          Metabolomics Workbench, rest/refmet/all_ids - 37,557 rows
                  carrying name, pubchem_cid, chebi_id, kegg_id, hmdb_id,
                  lm_id, inchi_key AND the 3-level chemical taxonomy
                  (super_class/main_class/sub_class). Measured coverage:
                  chebi 66.7%, lm 37.0%, kegg 37.4%, hmdb 26.7%.
                  NOTE: the sibling rest/refmet/all endpoint (43MB) is
                  deliberately NOT fetched - it adds only exact mass and
                  formula, which no enrichment test here uses.
                  LICENCE: Metabolomics Workbench's terms of use say
                  "personal, non-commercial use only". cache/ is gitignored,
                  so this repo ships code, not data: each deployment's own
                  re-fetch binds that user, not this codebase. This applies
                  ONLY to the Chemical Class analysis; every other source
                  below is freely licensed, so skipping that one tile leaves
                  a commercial user unaffected.
  LIPID MAPS      LMSD SDF (structures.sdf inside a ~21MB zip). Supplies the
                  lipid taxonomy (CATEGORY/MAIN_CLASS/SUB_CLASS) AND a full
                  crosswalk (INCHI_KEY, PUBCHEM_CID, HMDB_ID, KEGG_ID,
                  CHEBI_ID). CC BY 4.0 - so Lipid Class enrichment depends on
                  no non-commercial source.
  ChEBI           flat files, crosswalk only (ChEBI ID <-> KEGG/InChI).
                  CC BY 4.0. Optional: see the degradability note above.
  KEGG            link/pathway/compound, via Genomics_module._kegg_safe_get.
  Reactome        ChEBI2Reactome_All_Levels.txt. CC BY 4.0.
  SMPDB           smpdb_metabolites.csv.zip. Human-only. CC BY 4.0.
  PathBank        pathbank_primary_metabolites.csv.zip. Multi-organism. ODbL.

Pathway Topology Analysis (see run_pathway_topology near the end of this
file): adds a topology "impact" score to KEGG pathway enrichment, from each
matched compound's position in that pathway's own reaction graph - relative
betweenness centrality (Freeman 1977, via Brandes' algorithm, Brandes 2001,
J. Math. Sociol. 25:163-177) by default, or out-degree. The impact formula
follows MetPA (Xia & Wishart 2010, Bioinformatics 26:2342); Brandes'
algorithm is this module's own implementation (no networkx dependency -
none is installed here, and KEGG pathway graphs are small). KEGG's KGML
reaction graphs are organism-specific in their EDGES (verified live
2026-09-29: hsa00020/sce00020/eco00020 share the same 23 compounds but have
22/21/21 reactions respectively) even though compound MEMBERSHIP is not
(see build_kegg_organism_pathway_index) - so the graph is fetched and cached
per organism per map, separately from the compound-membership caches above.
"""

import csv
import io
import json
import os
import re
import time
import zipfile
import xml.etree.ElementTree as ET
from collections import deque
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
    OMISCERA_HEADERS,
    MIN_TERM_SIZE, MAX_TERM_SIZE_DEFAULT, FDR_ALPHA_DEFAULT,
    ORA_RESULT_COLUMNS, BURDEN_RESULT_COLUMNS,
    benjamini_hochberg, bonferroni,
    _run_fisher_for_namespace, _run_wilcoxon_for_namespace,
    _load_universe_cache, _save_universe_cache,
    _kegg_safe_get, _parse_kegg_pathway_list, KEGG_LIST_PATHWAY_URL,
    load_user_background_list,
)
# KEGG organism resolution is organism-only (taxon in, KEGG organism code out),
# so it is reused unchanged - exactly as Proteomics_module does - rather than
# re-derived. See that function's docstring for the divergence-bug history
# (commit 27ce033) this design avoids.
from Transcriptomics_module import resolve_kegg_organism_code
# KEGG_OVERVIEW_PATHWAY_NUMBERS: reused, not duplicated - the existing
# Metabolic Connection feature already documents these as "present in nearly
# every organism, so their presence is not evidence of species-specific
# occurrence" for gene/enzyme evidence; the same reasoning applies to compound
# membership (see KEGG_CHEMICAL_STRUCTURE_MAP_NUMBERS below for why compounds
# additionally need a second, lower-numbered range these gene-side maps don't
# cover).
from database_catalog import KEGG_OVERVIEW_PATHWAY_NUMBERS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METABOLOMICS_ORA_DATABASES = [
    "KEGG", "Reactome", "SMPDB", "PathBank", "Chemical Class", "Lipid Class",
]

# Quantitative Enrichment's default score. Perseus' 1D annotation enrichment
# (the method mirrored here) is defined directly on log2FC - but UNLIKE
# Proteomics_module (which reuses Transcriptomics' four log2FC/pvalue-based
# GSEA metric formulas because its input already carries per-protein log2FC/
# pvalue columns), a metabolomics upload here is either a plain compound
# LIST (no statistics at all) or a concentration TABLE whose statistics
# Omiscera computes itself (see compute_compound_statistics). So there is no
# metric menu to reuse - just the two columns compute_compound_statistics
# produces: the log2 fold change between groups, or Welch's t-statistic.
QUANT_METRIC_DEFAULT = "log2fc"
QUANT_METRICS = ("log2fc", "t_stat")

# Metabolomics renames the shared engines' member column: these are measured
# COMPOUNDS, not genes. Mirrors Proteomics_module.OVERLAP_COLUMN.
OVERLAP_COLUMN = "Overlapping Compounds"
METABOLOMICS_ORA_RESULT_COLUMNS = [
    OVERLAP_COLUMN if c == "Overlapping Genes" else c for c in ORA_RESULT_COLUMNS
]
METABOLOMICS_BURDEN_RESULT_COLUMNS = [
    {"Genes": OVERLAP_COLUMN,
     "Genes in set": "Compounds in set",
     "Scored genes (universe)": "Scored compounds (universe)"}.get(c, c)
    for c in BURDEN_RESULT_COLUMNS
]
_BURDEN_COLUMN_RENAMES = {
    "Genes": OVERLAP_COLUMN,
    "Genes in set": "Compounds in set",
    "Scored genes (universe)": "Scored compounds (universe)",
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "cache", "metabolomics")

REFRESH_DAYS = 30

# RefMet - see module docstring for why /all is deliberately not fetched.
REFMET_ALL_IDS_URL = "https://www.metabolomicsworkbench.org/rest/refmet/all_ids"
REFMET_ALL_IDS_LOCAL = "refmet_all_ids.json"

# LIPID MAPS LMSD. The /databases/lmsd/download page is a JS app; this is the
# actual file link it exposes (verified live 2026-09-21, ~21MB zip containing
# a single member "structures.sdf").
LIPIDMAPS_SDF_URL = "https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip"
LIPIDMAPS_SDF_LOCAL = "lipidmaps_lmsd.sdf.zip"
LIPIDMAPS_SDF_MEMBER = "structures.sdf"

# ChEBI flat files. NOTE: ChEBI 2.0 (launched Oct 2025) moved these from
# Flat_file_tab_delimited/ to flat_files/; both spellings are tried because
# the rename is recent and mirrors may lag. OPTIONAL source - the crosswalk
# degrades without it (see module docstring).
CHEBI_FLAT_BASES = [
    "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files",
    "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited",
]
CHEBI_ACCESSION_LOCAL = "chebi_database_accession.tsv"
CHEBI_INCHI_LOCAL = "chebi_id_inchi.tsv"

CHEBI2REACTOME_URL = "https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt"
CHEBI2REACTOME_LOCAL = "ChEBI2Reactome_All_Levels.txt"

SMPDB_METABOLITES_URL = "https://www.smpdb.ca/downloads/smpdb_metabolites.csv.zip"
SMPDB_LOCAL = "smpdb_metabolites.csv.zip"
SMPDB_HUMAN_TAXON_ID = "9606"

PATHBANK_PRIMARY_URL = "https://pathbank.org/downloads/pathbank_primary_metabolites.csv.zip"
PATHBANK_ALL_URL = "https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip"
PATHBANK_PRIMARY_LOCAL = "pathbank_primary_metabolites.csv.zip"
PATHBANK_ALL_LOCAL = "pathbank_all_metabolites.csv.zip"

# Cache-key versions, baked into the universe-cache database string exactly
# like Genomics' KEGG_v2/Reactome_v2 convention, so a parser change
# invalidates old caches instead of silently reading a stale shape.
CROSSWALK_CACHE_VERSION = 2
UNIVERSE_CACHE_VERSION = 2  # bumped: excludes KEGG's 01xxx overview maps now
# The per-organism KEGG map index is a SEPARATE cache artifact from the global
# compound->reference-map links above, so a change to how the organism list is
# parsed costs one small refetch per organism, never a re-download of the big
# link file.
KEGG_ORG_MAPS_CACHE_VERSION = 1
# A THIRD, separate cache artifact: one KGML reaction graph per organism per
# map (Pathway Topology Analysis only). Unlike the compound-membership caches
# above, this cannot be shared across organisms - see the module docstring.
KEGG_KGML_CACHE_VERSION = 2  # v2: nodes carry KEGG's own compound NAME, not just the C-number
# A fourth cache artifact: KEGG's own C-number -> common-name lookup
# (list/compound, ~19,000 compounds, one global 1MB fetch) - used to label
# Pathway Topology's graph-view nodes with real names ("Citrate") instead of
# a bare KEGG ID, since KGML's own graphics name for a compound entry is
# usually just the C-number itself (verified live 2026-09-29).
KEGG_COMPOUND_NAMES_CACHE_VERSION = 1

PATHWAY_TOPOLOGY_MEASURES = ("betweenness", "outdegree")
PATHWAY_TOPOLOGY_MEASURE_DEFAULT = "betweenness"

# KEGG reserves map01000-map01099 ("Overview of biosynthetic pathways",
# "Biosynthesis of plant secondary metabolites", "Biosynthesis of alkaloids
# derived from...", etc.) as GLOBAL reference charts that aggregate compounds
# across huge swaths of chemistry by chemical-structure class, not real
# individual pathways - confirmed live 2026-09-22 against KEGG's own
# list/pathway output. Verified on a real run: "Biosynthesis of terpenoids
# and steroids" (map01062) and "Biosynthesis of alkaloids derived from
# terpenoid and polyketide" (map01066) both outranked the genuinely correct
# "Citrate cycle (TCA cycle)" hit, each with k==K (every one of the map's
# annotated compounds in our universe was also significant) because these
# maps' real content is "most of central metabolism", not a specific route.
# database_catalog.KEGG_OVERVIEW_PATHWAY_NUMBERS (imported above) already
# excludes the SEPARATE 01100-01320 "Global and overview maps" range for gene
# evidence; this covers the lower, compound-specific "Chemical structure
# maps" range that range does not.
KEGG_CHEMICAL_STRUCTURE_MAP_NUMBERS = frozenset(str(n).zfill(5) for n in range(1000, 1100))
KEGG_COMPOUND_EXCLUDED_MAPS = frozenset(
    n.zfill(5) for n in KEGG_OVERVIEW_PATHWAY_NUMBERS
) | KEGG_CHEMICAL_STRUCTURE_MAP_NUMBERS

# The 8 LIPID MAPS categories, used to parse lipid shorthand names (see
# _lipid_shorthand_class) and to sanity-check the SDF's CATEGORY field.
LIPIDMAPS_CATEGORY_CODES = {
    "FA": "Fatty Acyls",
    "GL": "Glycerolipids",
    "GP": "Glycerophospholipids",
    "SP": "Sphingolipids",
    "ST": "Sterol Lipids",
    "PR": "Prenol Lipids",
    "SL": "Saccharolipids",
    "PK": "Polyketides",
}

# Lipid shorthand prefixes -> (LIPID MAPS category, a readable class label).
# Sum-composition names like "PC 34:1" denote a BULK SPECIES covering many
# structures, so they can never resolve to a single InChIKey - but the class
# is unambiguous, which is exactly what Lipid Class enrichment needs. See
# _lipid_shorthand_class for why these rows are deliberately excluded from
# pathway membership.
LIPID_SHORTHAND_PREFIXES = {
    "PC": ("Glycerophospholipids", "Glycerophosphocholines"),
    "LPC": ("Glycerophospholipids", "Glycerophosphocholines"),
    "PE": ("Glycerophospholipids", "Glycerophosphoethanolamines"),
    "LPE": ("Glycerophospholipids", "Glycerophosphoethanolamines"),
    "PS": ("Glycerophospholipids", "Glycerophosphoserines"),
    "PG": ("Glycerophospholipids", "Glycerophosphoglycerols"),
    "PI": ("Glycerophospholipids", "Glycerophosphoinositols"),
    "PA": ("Glycerophospholipids", "Glycerophosphates"),
    "CL": ("Glycerophospholipids", "Glycerophosphoglycerophosphoglycerols"),
    "TG": ("Glycerolipids", "Triradylglycerols"),
    "TAG": ("Glycerolipids", "Triradylglycerols"),
    "DG": ("Glycerolipids", "Diradylglycerols"),
    "DAG": ("Glycerolipids", "Diradylglycerols"),
    "MG": ("Glycerolipids", "Monoradylglycerols"),
    "MAG": ("Glycerolipids", "Monoradylglycerols"),
    "SM": ("Sphingolipids", "Phosphosphingolipids"),
    "CER": ("Sphingolipids", "Ceramides"),
    "HEXCER": ("Sphingolipids", "Neutral glycosphingolipids"),
    "GLCCER": ("Sphingolipids", "Neutral glycosphingolipids"),
    "LACCER": ("Sphingolipids", "Neutral glycosphingolipids"),
    "S1P": ("Sphingolipids", "Sphingoid bases"),
    "SPH": ("Sphingolipids", "Sphingoid bases"),
    "CE": ("Sterol Lipids", "Steryl esters"),
    "FA": ("Fatty Acyls", "Fatty Acids and Conjugates"),
}
# "PC 34:1", "PC(34:1)", "Cer d18:1/16:0", "TG 52:2" ...
_LIPID_SHORTHAND_RE = re.compile(
    r"^([A-Za-z]{2,7})[\s(_-]*[dtoOP]?-?\d{1,3}:\d{1,2}", re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _clean(value) -> str:
    """Normalise one cell to a stripped string ('' for NaN/None)."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in ("nan", "none", "null", "-"):
        return ""
    return text


def _norm_name(name: str) -> str:
    """Normalise a compound name for last-resort name matching."""
    text = _clean(name).lower()
    # Commas become spaces: this string is embedded in "NAME:" compound keys,
    # and the shared Fisher/Wilcoxon engines join and split member lists on
    # ",", so a comma inside a key ("glucose, 6-phosphate") would tear one
    # compound into two on the way to the results table.
    text = text.replace(",", " ")
    text = re.sub(r"[\s_]+", " ", text)
    text = re.sub(r"[^a-z0-9 :/()\[\].+-]", "", text)
    return text.strip()


def inchikey_skeleton(inchikey: str) -> str:
    """First block of an InChIKey - the connectivity/skeleton layer.

    An InChIKey is AAAAAAAAAAAAAA-BBBBBBBBFV-C: block 1 encodes connectivity,
    block 2 stereochemistry and isotopes, block 3 protonation. Matching only
    full keys is strict (exact stereoisomer) and loses a large share of lipid
    matches, because databases routinely disagree on how completely
    stereochemistry is specified - a measured lipid reported without
    sn-position versus LIPID MAPS' fully-specified structure. Matching on
    block 1 alone silently conflates genuine stereoisomers. So callers try
    the full key first and fall back to this, COUNTING the fallback
    separately (see match_compounds' n_skeleton_only) so a user can see how
    much of a result rests on stereochemistry-insensitive matching.
    """
    text = _clean(inchikey).upper()
    if not text:
        return ""
    return text.split("-")[0]


def _compound_key(inchikey: str = "", pubchem_cid: str = "", chebi_id: str = "",
                  kegg_id: str = "", lm_id: str = "", hmdb_id: str = "",
                  name: str = "") -> str:
    """The module's internal identifier, in preference order.

    The full InChIKey is preferred because it is structure-derived and shared
    across every source here. The namespaced fallbacks keep compounds that
    lack one (common for bulk lipid species and for database entries with no
    deposited structure) in ONE identifier space rather than silently
    dropping them.
    """
    if _clean(inchikey):
        return _clean(inchikey).upper()
    if _clean(pubchem_cid):
        return "PUBCHEM:" + _clean(pubchem_cid)
    if _clean(chebi_id):
        return "CHEBI:" + _normalize_chebi(chebi_id)
    if _clean(kegg_id):
        return "KEGG:" + _clean(kegg_id).upper()
    if _clean(lm_id):
        return "LM:" + _clean(lm_id).upper()
    if _clean(hmdb_id):
        return "HMDB:" + _clean(hmdb_id).upper()
    if _norm_name(name):
        return "NAME:" + _norm_name(name)
    return ""


def _normalize_chebi(value) -> str:
    """ChEBI IDs appear bare ("15756") and prefixed ("CHEBI:15756").

    RefMet's chebi_id and Reactome's ChEBI2Reactome column are both bare;
    LIPID MAPS' CHEBI_ID is bare; SMPDB/PathBank use the prefixed form. All
    are normalised to bare digits so the crosswalk has one spelling.
    """
    text = _clean(value).upper().replace("CHEBI:", "").strip()
    return text if text.isdigit() else ""


def _normalize_kegg_compound(value) -> str:
    """KEGG compound IDs: accept "C00249", "cpd:C00249", lowercase."""
    text = _clean(value).upper().replace("CPD:", "").strip()
    return text if re.fullmatch(r"C\d{5}", text) else ""


def _normalize_hmdb(value) -> str:
    """HMDB IDs come 5-digit (HMDB00220) and 7-digit (HMDB0000220)."""
    text = _clean(value).upper()
    match = re.search(r"HMDB(\d+)", text)
    if not match:
        return ""
    return "HMDB" + match.group(1).zfill(7)


def _strip_bracket_code(value: str) -> "tuple[str, str]":
    """LIPID MAPS writes "Fatty Acyls [FA]" / "Straight chain fatty acids
    [FA0101]". Returns (label, code) with the code stripped from the label."""
    text = _clean(value)
    match = re.search(r"\[([A-Z0-9]+)\]\s*$", text)
    if not match:
        return text, ""
    return text[: match.start()].strip(), match.group(1)


def _lipid_shorthand_class(name: str) -> "tuple[str, str] | None":
    """Parse a sum-composition lipid name to (category, class label).

    "PC 34:1" names a BULK SPECIES - one sum composition covering many
    distinct structures - so it cannot resolve to any single InChIKey, and
    structure-level or pathway-level matching is impossible for it by
    definition. Its CLASS, though, is unambiguous and is exactly what Lipid
    Class enrichment tests. Callers must therefore use this tier for the
    Lipid Class analysis ONLY and flag the row class-level-only; letting such
    a row into a pathway universe would assert membership that the name does
    not actually support.
    """
    text = _clean(name)
    if not text:
        return None
    match = _LIPID_SHORTHAND_RE.match(text)
    if not match:
        return None
    prefix = match.group(1).upper()
    return LIPID_SHORTHAND_PREFIXES.get(prefix)


def _ensure_bulk_file(url: str, local_name: str, refresh_days: int = REFRESH_DAYS,
                      verbose: bool = False, force_refresh: bool = False,
                      timeout: int = 300) -> "str | None":
    """Download a bulk file if missing/stale; a stale copy beats nothing.

    Mirrors Proteomics_module._ensure_diseases_file exactly (streaming write
    to .tmp then os.replace, so an interrupted download can never leave a
    truncated file that later parses as a valid-but-short universe).
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    local_path = os.path.join(CACHE_DIR, local_name)
    if not force_refresh and os.path.exists(local_path):
        age_days = (time.time() - os.path.getmtime(local_path)) / 86400.0
        if age_days < refresh_days:
            return local_path
    if HTTP_SESSION is None:
        return local_path if os.path.exists(local_path) else None
    try:
        resp = HTTP_SESSION.get(url, headers=OMISCERA_HEADERS, timeout=timeout, stream=True)
        if resp.status_code == 200:
            tmp_path = local_path + ".tmp"
            with open(tmp_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1024 * 256):
                    if chunk:
                        fh.write(chunk)
            os.replace(tmp_path, local_path)
            if verbose:
                size_mb = os.path.getsize(local_path) / 1e6
                print(f"[metabolomics] fetched {local_name} ({size_mb:.1f} MB)")
            return local_path
        if verbose:
            print(f"[metabolomics] fetch failed for {local_name}: HTTP {resp.status_code}")
    except Exception as exc:
        if verbose:
            print(f"[metabolomics] fetch failed for {local_name}: {exc}")
    return local_path if os.path.exists(local_path) else None


def _iter_sdf_records(path: str, member: str = LIPIDMAPS_SDF_MEMBER):
    """Yield {TAG: value} dicts from an SDF inside a zip, streaming.

    Only the "> <TAG>" data fields are read; the molblock is skipped entirely
    (this module never needs atom coordinates, so no RDKit dependency and no
    need to hold a ~100MB decompressed file in memory).
    """
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        target = member if member in names else (names[0] if names else None)
        if target is None:
            return
        with zf.open(target) as raw:
            stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            record: dict = {}
            tag = None
            buffer: list = []
            for line in stream:
                line = line.rstrip("\n").rstrip("\r")
                if line.startswith("$$$$"):
                    if tag is not None and buffer:
                        record[tag] = "\n".join(buffer).strip()
                    if record:
                        yield record
                    record, tag, buffer = {}, None, []
                    continue
                if line.startswith("> <") and line.endswith(">"):
                    if tag is not None and buffer:
                        record[tag] = "\n".join(buffer).strip()
                    tag = line[3:-1]
                    buffer = []
                    continue
                if tag is not None:
                    if line.strip() == "":
                        # Blank line terminates a tag's value block.
                        if buffer:
                            record[tag] = "\n".join(buffer).strip()
                            tag, buffer = None, []
                    else:
                        buffer.append(line)
            if tag is not None and buffer:
                record[tag] = "\n".join(buffer).strip()
            if record:
                yield record


# ---------------------------------------------------------------------------
# The crosswalk layer - see module docstring for why this exists
# ---------------------------------------------------------------------------

_CROSSWALK_FIELDS = (
    "name", "synonyms", "inchikey", "skeleton", "pubchem_cid", "chebi_id", "kegg_id",
    "hmdb_id", "lm_id", "refmet_id",
    "super_class", "main_class", "sub_class",          # RefMet taxonomy
    "lm_category", "lm_main_class", "lm_sub_class",    # LIPID MAPS taxonomy
)


def _blank_entry() -> dict:
    return {field: "" for field in _CROSSWALK_FIELDS}


def _merge_entry(target: dict, incoming: dict) -> None:
    """Fill empty fields of `target` from `incoming`, never overwriting.

    First writer wins per field. Sources are merged in a deliberate order
    (LIPID MAPS, then RefMet, then ChEBI) so the more specific taxonomy owns
    its own fields, and no source can clobber an identifier another already
    supplied.
    """
    for field in _CROSSWALK_FIELDS:
        if not target.get(field) and incoming.get(field):
            target[field] = incoming[field]


def _parse_lipidmaps(path: str, verbose: bool = False) -> "list[dict]":
    """LIPID MAPS LMSD -> crosswalk entries (lipid taxonomy + IDs).

    Verified live 2026-09-21 against a real record (Palmitic acid,
    LMFA01010001): CATEGORY "Fatty Acyls [FA]", MAIN_CLASS "Fatty Acids and
    Conjugates [FA01]", SUB_CLASS "Straight chain fatty acids [FA0101]",
    plus INCHI_KEY / PUBCHEM_CID / KEGG_ID / HMDB_ID / CHEBI_ID. Not every
    record carries every field (the first record in the file has no
    SUB_CLASS, KEGG_ID or HMDB_ID at all), which is why each is read
    defensively rather than assumed present.
    """
    entries = []
    for rec in _iter_sdf_records(path):
        lm_id = _clean(rec.get("LM_ID"))
        if not lm_id:
            continue
        category, _ = _strip_bracket_code(rec.get("CATEGORY", ""))
        main_class, _ = _strip_bracket_code(rec.get("MAIN_CLASS", ""))
        sub_class, _ = _strip_bracket_code(rec.get("SUB_CLASS", ""))
        inchikey = _clean(rec.get("INCHI_KEY")).upper()
        entry = _blank_entry()
        entry.update({
            "name": _clean(rec.get("NAME")) or _clean(rec.get("SYSTEMATIC_NAME")),
            # LIPID MAPS' SYNONYMS is a semicolon-separated list. Indexing it
            # materially improves name matching: databases disagree on which
            # name is canonical ("Palmitic acid" vs "hexadecanoic acid" vs
            # "C16:0"), and a user's table uses whichever their instrument
            # software emitted.
            "synonyms": _clean(rec.get("SYNONYMS")),
            "inchikey": inchikey,
            "skeleton": inchikey_skeleton(inchikey),
            "pubchem_cid": _clean(rec.get("PUBCHEM_CID")),
            "chebi_id": _normalize_chebi(rec.get("CHEBI_ID")),
            "kegg_id": _normalize_kegg_compound(rec.get("KEGG_ID")),
            "hmdb_id": _normalize_hmdb(rec.get("HMDB_ID")),
            "lm_id": lm_id.upper(),
            "lm_category": category,
            "lm_main_class": main_class,
            "lm_sub_class": sub_class,
        })
        entries.append(entry)
    if verbose:
        print(f"[metabolomics] LIPID MAPS: {len(entries)} records")
    return entries


def _parse_refmet(path: str, verbose: bool = False) -> "list[dict]":
    """RefMet all_ids -> crosswalk entries (chemical taxonomy + IDs).

    Response shape verified live 2026-09-21: a JSON OBJECT keyed "Row1",
    "Row2", ... (not a list), 37,557 rows. Fields are omitted entirely when
    absent rather than being present-but-empty, so every read is a .get().
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        raw = json.load(fh)
    rows = raw.values() if isinstance(raw, dict) else raw
    entries = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        inchikey = _clean(row.get("inchi_key")).upper()
        entry = _blank_entry()
        entry.update({
            "name": _clean(row.get("name")),
            "inchikey": inchikey,
            "skeleton": inchikey_skeleton(inchikey),
            "pubchem_cid": _clean(row.get("pubchem_cid")),
            "chebi_id": _normalize_chebi(row.get("chebi_id")),
            "kegg_id": _normalize_kegg_compound(row.get("kegg_id")),
            "hmdb_id": _normalize_hmdb(row.get("hmdb_id")),
            "lm_id": _clean(row.get("lm_id")).upper(),
            "refmet_id": _clean(row.get("refmet_id")),
            "super_class": _clean(row.get("super_class")),
            "main_class": _clean(row.get("main_class")),
            "sub_class": _clean(row.get("sub_class")),
        })
        entries.append(entry)
    if verbose:
        print(f"[metabolomics] RefMet: {len(entries)} records")
    return entries


def _ensure_chebi_files(verbose: bool = False, force_refresh: bool = False) -> "tuple[str | None, str | None]":
    """ChEBI's two flat files, trying both the new and legacy path spellings.

    OPTIONAL source. ftp.ebi.ac.uk was verified completely unreachable on
    2026-09-21 (no ICMP, connection refused on 80 and 443, direct-IP too)
    from a machine where www.ebi.ac.uk answered normally - and it is the only
    host serving ChEBI bulk files. Returns (accession_path, inchi_path), with
    None for whichever could not be obtained; the crosswalk continues without
    them and says so in its caveat.
    """
    accession_path = inchi_path = None
    for base in CHEBI_FLAT_BASES:
        if accession_path is None:
            accession_path = _ensure_bulk_file(
                f"{base}/database_accession.tsv", CHEBI_ACCESSION_LOCAL,
                verbose=verbose, force_refresh=force_refresh,
            )
        if inchi_path is None:
            inchi_path = _ensure_bulk_file(
                f"{base}/chebiId_inchi.tsv", CHEBI_INCHI_LOCAL,
                verbose=verbose, force_refresh=force_refresh,
            )
        if accession_path and inchi_path:
            break
    return accession_path, inchi_path


def _parse_chebi(accession_path: "str | None", inchi_path: "str | None",
                 verbose: bool = False) -> "list[dict]":
    """ChEBI flat files -> crosswalk entries (identifiers only, no taxonomy).

    database_accession.tsv maps a ChEBI ID to accessions in other resources
    (KEGG COMPOUND, HMDB, LIPID MAPS, PubChem); chebiId_inchi.tsv maps it to
    an InChI, from which the InChIKey's skeleton cannot be derived without a
    cheminformatics toolkit - so only the InChI's presence is used to keep the
    entry, and structure-level joining for ChEBI-only compounds happens
    through the shared IDs instead. Both files are optional.
    """
    by_chebi: dict = {}
    if accession_path and os.path.exists(accession_path):
        try:
            with open(accession_path, "r", encoding="utf-8", errors="replace") as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None) or []
                lowered = [h.strip().lower() for h in header]

                def _col(*names, default=None):
                    for name in names:
                        if name in lowered:
                            return lowered.index(name)
                    return default

                i_chebi = _col("compound_id", "chebi_id", default=1)
                i_type = _col("type", "source", default=3)
                i_acc = _col("accession_number", "accession", default=4)
                for parts in reader:
                    if len(parts) <= max(i_chebi, i_type, i_acc):
                        continue
                    chebi = _normalize_chebi(parts[i_chebi])
                    if not chebi:
                        continue
                    src = _clean(parts[i_type]).upper()
                    acc = _clean(parts[i_acc])
                    entry = by_chebi.setdefault(chebi, _blank_entry())
                    entry["chebi_id"] = chebi
                    if "KEGG" in src and not entry["kegg_id"]:
                        entry["kegg_id"] = _normalize_kegg_compound(acc)
                    elif "HMDB" in src and not entry["hmdb_id"]:
                        entry["hmdb_id"] = _normalize_hmdb(acc)
                    elif "LIPID" in src and not entry["lm_id"]:
                        entry["lm_id"] = acc.upper()
                    elif "PUBCHEM" in src and not entry["pubchem_cid"]:
                        entry["pubchem_cid"] = acc if acc.isdigit() else ""
        except Exception as exc:
            if verbose:
                print(f"[metabolomics] ChEBI accession parse failed: {exc}")
    if verbose:
        print(f"[metabolomics] ChEBI: {len(by_chebi)} records")
    return [e for e in by_chebi.values() if e.get("chebi_id")]


def build_compound_crosswalk(verbose: bool = False, force_refresh: bool = False,
                             progress_callback=None) -> dict:
    """Merge LIPID MAPS + RefMet + ChEBI into one identifier index.

    Returns {"by_key": {compound_key: entry}, plus lookup indexes
    "by_inchikey"/"by_skeleton"/"by_pubchem"/"by_chebi"/"by_kegg"/"by_hmdb"/
    "by_lm_id"/"by_name" (each {value: compound_key}), "sources" (which
    contributed), "universe_size", "fetched_at", "caveat"}.

    Every source is independently optional: a failed fetch subtracts that
    source's contribution and is named in the caveat, never silently reducing
    the result. This is the slowest first-run step in the module (~30MB of
    downloads plus a full SDF parse), hence progress_callback.
    """
    cache_db = f"Metabolomics_Crosswalk_v{CROSSWALK_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, "global")
        if cached is not None and cached.get("by_key"):
            return cached

    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    entries: list = []
    sources: list = []
    missing: list = []

    _progress(0, 4, "Fetching LIPID MAPS structures...")
    lm_path = _ensure_bulk_file(LIPIDMAPS_SDF_URL, LIPIDMAPS_SDF_LOCAL,
                                verbose=verbose, force_refresh=force_refresh)
    if lm_path:
        try:
            entries.extend(_parse_lipidmaps(lm_path, verbose=verbose))
            sources.append("LIPID MAPS")
        except Exception as exc:
            if verbose:
                print(f"[metabolomics] LIPID MAPS parse failed: {exc}")
            missing.append("LIPID MAPS")
    else:
        missing.append("LIPID MAPS")

    _progress(1, 4, "Fetching RefMet classification...")
    refmet_path = _ensure_bulk_file(REFMET_ALL_IDS_URL, REFMET_ALL_IDS_LOCAL,
                                    verbose=verbose, force_refresh=force_refresh)
    if refmet_path:
        try:
            entries.extend(_parse_refmet(refmet_path, verbose=verbose))
            sources.append("RefMet")
        except Exception as exc:
            if verbose:
                print(f"[metabolomics] RefMet parse failed: {exc}")
            missing.append("RefMet")
    else:
        missing.append("RefMet")

    _progress(2, 4, "Fetching ChEBI cross-references...")
    chebi_acc, chebi_inchi = _ensure_chebi_files(verbose=verbose, force_refresh=force_refresh)
    chebi_entries = _parse_chebi(chebi_acc, chebi_inchi, verbose=verbose) if chebi_acc else []
    if chebi_entries:
        entries.extend(chebi_entries)
        sources.append("ChEBI")
    else:
        missing.append("ChEBI")

    _progress(3, 4, "Indexing compounds...")
    by_key: dict = {}
    for entry in entries:
        key = _compound_key(
            inchikey=entry.get("inchikey", ""), pubchem_cid=entry.get("pubchem_cid", ""),
            chebi_id=entry.get("chebi_id", ""), kegg_id=entry.get("kegg_id", ""),
            lm_id=entry.get("lm_id", ""), hmdb_id=entry.get("hmdb_id", ""),
            name=entry.get("name", ""),
        )
        if not key:
            continue
        if key in by_key:
            _merge_entry(by_key[key], entry)
        else:
            by_key[key] = dict(entry)

    index_specs = [
        ("by_inchikey", "inchikey", lambda v: v.upper()),
        ("by_skeleton", "skeleton", lambda v: v.upper()),
        ("by_pubchem", "pubchem_cid", lambda v: v),
        ("by_chebi", "chebi_id", lambda v: v),
        ("by_kegg", "kegg_id", lambda v: v.upper()),
        ("by_hmdb", "hmdb_id", lambda v: v.upper()),
        ("by_lm_id", "lm_id", lambda v: v.upper()),
        ("by_name", "name", _norm_name),
    ]
    indexes = {name: {} for name, _, _ in index_specs}
    for key, entry in by_key.items():
        for index_name, field, transform in index_specs:
            value = _clean(entry.get(field))
            if not value:
                continue
            transformed = transform(value)
            if transformed:
                # First writer wins: a later duplicate must not silently
                # repoint an identifier at a different compound.
                indexes[index_name].setdefault(transformed, key)
        # Synonyms feed the NAME index only, and only as a fallback: a
        # primary name already claimed above must never be displaced by
        # something else's synonym.
        for synonym in _clean(entry.get("synonyms")).split(";"):
            normalized = _norm_name(synonym)
            if normalized:
                indexes["by_name"].setdefault(normalized, key)

    caveat = None
    if missing:
        lost = {
            "LIPID MAPS": "lipid class enrichment and LM_ID matching are unavailable",
            "RefMet": "chemical class enrichment is unavailable",
            "ChEBI": (
                "ChEBI cross-references are reduced to those LIPID MAPS and RefMet "
                "supply, so some compounds may not reach Reactome/SMPDB/PathBank"
            ),
        }
        caveat = "Could not load: " + "; ".join(
            f"{name} ({lost.get(name, 'reduced coverage')})" for name in missing
        ) + ". Results below use the sources that did load."

    result = {
        "by_key": by_key,
        "sources": sources,
        "universe_size": len(by_key),
        "fetched_at": _now_iso(),
        "caveat": caveat,
    }
    result.update(indexes)
    _progress(4, 4, "Done")

    if by_key:
        _save_universe_cache(cache_db, "global", result)
    return result


# ---------------------------------------------------------------------------
# Matching the user's resolved compounds into crosswalk space
# ---------------------------------------------------------------------------

# Join tiers, most specific first. Order matters: a full InChIKey identifies
# one exact stereoisomer, while a name is ambiguous, so a lower tier must
# only ever be consulted when every higher one missed.
_MATCH_TIERS = (
    ("inchikey", "by_inchikey", lambda row: _clean(row.get("InChIKey")).upper()),
    ("pubchem", "by_pubchem", lambda row: _clean(row.get("PubChem ID"))),
    ("chebi", "by_chebi", lambda row: _normalize_chebi(row.get("ChEBI ID"))),
    ("kegg", "by_kegg", lambda row: _normalize_kegg_compound(row.get("KEGG ID"))),
    ("hmdb", "by_hmdb", lambda row: _normalize_hmdb(row.get("HMDB ID"))),
    ("lm_id", "by_lm_id", lambda row: _clean(row.get("LIPID MAPS ID")).upper()),
    ("name", "by_name", lambda row: _norm_name(row.get("Name"))),
)

MATCHED_COLUMNS = [
    "Name", "compound_key", "match_tier", "class_level_only",
    "inchikey", "pubchem_cid", "chebi_id", "kegg_id", "hmdb_id", "lm_id",
    "super_class", "main_class", "sub_class",
    "lm_category", "lm_main_class", "lm_sub_class",
]


def match_compounds(resolved_df: pd.DataFrame, crosswalk: dict,
                    user_id_mask: "pd.Series | None" = None,
                    verbose: bool = False) -> "tuple[pd.DataFrame, dict]":
    """Join the resolved compound table into crosswalk (compound_key) space.

    Tier order: full InChIKey -> PubChem CID -> ChEBI ID -> KEGG ID -> HMDB ID
    -> LM_ID -> normalized name; then an InChIKey-SKELETON fallback, counted
    separately (see inchikey_skeleton for why stereochemistry-insensitive
    matching must be visible rather than silent); then lipid-shorthand class
    parsing, which yields a CLASS but no structure, so those rows are flagged
    class_level_only and must never enter a pathway universe.

    user_id_mask: per-row True where the USER'S OWN upload carried a
    structure identifier (InChIKey / PubChem CID / ChEBI ID / LIPID MAPS ID).
    It matters only for lipid-shorthand names: the retriever also synthesises
    such identifiers FROM a name (PubChem returns one representative
    structure for "PC 34:1"), so the presence of an identifier on the
    resolved row proves nothing about whether a real structure was known.
    Without this mask every shorthand name is treated as class-level-only,
    which is the conservative and correct default; with it, a row where the
    user themselves pinned a structure keeps full pathway capability.

    Returns (matched_df, meta). matched_df has one row per input row, with
    MATCHED_COLUMNS; unmatched rows keep an empty compound_key so they are
    excluded from every universe but still visible to the user.
    """
    by_key = crosswalk.get("by_key", {})
    rows = []
    tier_counts = {name: 0 for name, _, _ in _MATCH_TIERS}
    tier_counts["skeleton"] = 0
    tier_counts["lipid_shorthand"] = 0
    n_unmatched = 0

    for position, (_, src) in enumerate(resolved_df.iterrows()):
        user_pinned = bool(user_id_mask.iloc[position]) if user_id_mask is not None else False
        row = {col: "" for col in MATCHED_COLUMNS}
        row["Name"] = _clean(src.get("Name"))
        row["class_level_only"] = False

        matched_key = ""
        matched_tier = ""
        for tier_name, index_name, extract in _MATCH_TIERS:
            value = extract(src)
            if not value:
                continue
            key = crosswalk.get(index_name, {}).get(value)
            if key:
                matched_key, matched_tier = key, tier_name
                break

        # Skeleton fallback - only after every exact tier missed.
        if not matched_key:
            skeleton = inchikey_skeleton(_clean(src.get("InChIKey")))
            if skeleton:
                key = crosswalk.get("by_skeleton", {}).get(skeleton)
                if key:
                    matched_key, matched_tier = key, "skeleton"

        # Lipid shorthand - a CLASS with no structure. Deliberately last, and
        # deliberately not a pathway-capable match.
        if not matched_key:
            shorthand = _lipid_shorthand_class(row["Name"])
            if shorthand:
                category, class_label = shorthand
                row["lm_category"] = category
                row["lm_main_class"] = class_label
                row["lm_sub_class"] = ""
                row["class_level_only"] = True
                row["compound_key"] = "SHORTHAND:" + _norm_name(row["Name"])
                row["match_tier"] = "lipid_shorthand"
                tier_counts["lipid_shorthand"] += 1
                rows.append(row)
                continue

        if matched_key:
            entry = by_key.get(matched_key, {})
            row["compound_key"] = matched_key
            row["match_tier"] = matched_tier
            tier_counts[matched_tier] = tier_counts.get(matched_tier, 0) + 1
            for field in ("inchikey", "pubchem_cid", "chebi_id", "kegg_id", "hmdb_id",
                          "lm_id", "super_class", "main_class", "sub_class",
                          "lm_category", "lm_main_class", "lm_sub_class"):
                row[field] = _clean(entry.get(field))
            # A shorthand name that matched only BY NAME is still a bulk
            # species, even when the source has an entry for it. RefMet does
            # carry bulk species ("PC 34:1") and assigns them a REPRESENTATIVE
            # ChEBI ID - excellent for classification, but letting that ID
            # through would silently assert that one sum composition covering
            # many structures participates in a specific reaction. So keep the
            # (better) curated class and withhold only the pathway claim.
            # A row whose match came from a structure-specific identifier the
            # user supplied is exempt: they have pinned a single structure.
            # Flag on the NAME, not on the match tier. The retriever
            # synthesises a PubChem CID/InChIKey from a shorthand name too
            # (PubChem returns one representative structure for "PC 34:1"),
            # so an identifier on the resolved row is NOT evidence that a
            # real structure was known - an earlier tier-based exemption let
            # exactly that case slip back into pathway universes. Only a
            # structure the USER supplied earns the exemption.
            shorthand = None if user_pinned else _lipid_shorthand_class(row["Name"])
            if shorthand:
                row["class_level_only"] = True
                # NOT counted as a lipid_shorthand tier hit: the tier counts
                # must partition the matched rows (this row was matched by
                # name), otherwise they sum to more than n_matched and the
                # resolve summary silently double-counts.
                # RefMet classifies bulk species under its own taxonomy but
                # carries no LIPID MAPS fields for them, so fill lm_* from the
                # shorthand prefix - without this a lipidomics upload matched
                # via RefMet produces an entirely empty Lipid Class analysis.
                if not row["lm_category"]:
                    row["lm_category"] = shorthand[0]
                if not row["lm_main_class"]:
                    row["lm_main_class"] = shorthand[1]
        else:
            n_unmatched += 1
        rows.append(row)

    matched_df = pd.DataFrame(rows, columns=MATCHED_COLUMNS)
    meta = {
        "n_input": len(resolved_df),
        "n_matched": int((matched_df["compound_key"].astype(str) != "").sum()),
        "n_unmatched": n_unmatched,
        "n_matched_by_tier": tier_counts,
        "n_skeleton_only": tier_counts.get("skeleton", 0),
        "n_class_level_only": tier_counts.get("lipid_shorthand", 0),
        "crosswalk_sources": crosswalk.get("sources", []),
        "crosswalk_caveat": crosswalk.get("caveat"),
    }
    if verbose:
        print(f"[metabolomics] matched {meta['n_matched']}/{meta['n_input']} "
              f"(tiers: {tier_counts})")
    return matched_df, meta


# ---------------------------------------------------------------------------
# Universe builders
#
# TWO universe rules coexist here, and every builder is told explicitly
# which one is in force via `universe_scope` - a block must never leave a
# caller guessing which rule produced its numbers:
#
#   universe_scope="experiment" - the universe is compounds actually
#     QUANTIFIED in this run that carry an annotation. Correct for the
#     Quantitative Enrichment tile (a concentration table): every compound
#     was measured, so testing against a database's WHOLE annotated set
#     would inflate every p-value by orders of magnitude (the same bug
#     Proteomics_module.run_ora documents having fixed for protein ORA).
#
#   universe_scope="database" - the universe is the DATABASE's own full
#     annotated set (its reference universe, exactly like testing a plain
#     gene list against all of GO/KEGG elsewhere in this app). Correct for
#     the ORA / Chemical Class / Lipid Class tiles: their input is a plain
#     compound LIST with no "measured" concept at all - if the universe
#     were the user's own list, every term would contain 100% of it and
#     Fisher's test is degenerate (k=K, n=N, p=1 for every term).
#
# The four PATHWAY universes translate an external file's own identifiers
# (KEGG compound ID / ChEBI ID) into this module's compound_key space via
# _id_to_key_maps, whose `scope` argument picks which identifier maps to use
# (matched_df's own, vs. the full crosswalk's). The two CLASS universes are
# built directly from either matched_df or the crosswalk's own "by_key"
# table - see _build_class_universe / _build_class_universe_from_crosswalk.
# ---------------------------------------------------------------------------

CLASS_LEVELS = ("super_class", "main_class", "sub_class")
LIPID_CLASS_LEVELS = ("lm_category", "lm_main_class", "lm_sub_class")
CLASS_LEVEL_LABELS = {
    "super_class": "Super class", "main_class": "Main class", "sub_class": "Sub class",
    "lm_category": "Category", "lm_main_class": "Main class", "lm_sub_class": "Sub class",
}


def _build_class_universe(matched_df: pd.DataFrame, level: str, label_prefix: str,
                          exclude_class_level_only: bool = False) -> dict:
    """Group matched compounds into {class label: compound_key set}.
    universe_scope="experiment" builder - see the section comment above."""
    terms: dict = {}
    df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
    if exclude_class_level_only and "class_level_only" in df.columns:
        df = df[~df["class_level_only"].astype(bool)]
    for _, row in df.iterrows():
        class_label = _clean(row.get(level))
        if not class_label:
            continue
        term_id = f"{label_prefix}:{class_label}"
        info = terms.setdefault(term_id, {"name": class_label, "genes": set()})
        info["genes"].add(str(row["compound_key"]))
    universe = set()
    for info in terms.values():
        universe |= info["genes"]
    return {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": _now_iso(), "caveat": None,
    }


def _build_class_universe_from_crosswalk(crosswalk: dict, level: str, label_prefix: str) -> dict:
    """Group the WHOLE crosswalk into {class label: compound_key set}.
    universe_scope="database" builder - the database's own full reference
    set, for a plain compound-list upload with no "measured" concept."""
    terms: dict = {}
    for key, entry in crosswalk.get("by_key", {}).items():
        class_label = _clean(entry.get(level))
        if not class_label:
            continue
        term_id = f"{label_prefix}:{class_label}"
        info = terms.setdefault(term_id, {"name": class_label, "genes": set()})
        info["genes"].add(key)
    universe = set()
    for info in terms.values():
        universe |= info["genes"]
    return {
        "terms": terms, "universe_size": len(universe),
        "fetched_at": _now_iso(), "caveat": None,
    }


def build_chemical_class_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                                  level: str = "sub_class",
                                  universe_scope: str = "experiment") -> dict:
    """RefMet chemical taxonomy universe (super_class/main_class/sub_class)."""
    if level not in CLASS_LEVELS:
        level = "sub_class"
    if universe_scope == "database" and crosswalk:
        result = _build_class_universe_from_crosswalk(crosswalk, level, "CHEMCLASS")
    else:
        result = _build_class_universe(matched_df, level, "CHEMCLASS")
    if not result["terms"]:
        result["caveat"] = (
            "No RefMet chemical classification matched these compounds. RefMet covers "
            "~37,500 commonly measured metabolites; compounds outside it (or matched only "
            "by lipid shorthand) carry no chemical class."
        )
    return result


def build_lipid_class_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                               level: str = "lm_main_class",
                               universe_scope: str = "experiment") -> dict:
    """LIPID MAPS lipid taxonomy universe (CATEGORY/MAIN_CLASS/SUB_CLASS).

    Unlike the pathway universes this one DOES keep class-level-only
    (lipid-shorthand) rows: a sum-composition name like "PC 34:1" cannot
    support a pathway claim, but its lipid class is exactly what this
    analysis tests.

    Under universe_scope="database" the reference universe comes from the
    crosswalk's own LMSD-derived entries - but a shorthand name like
    "PC 34:1" has NO individual LMSD entry at all (it names a bulk species,
    not one structure), so it can never appear in that reference set. Left
    alone, a shorthand-heavy lipidomics list would score its OWN
    best-supported analysis as completely empty. Fix: union the user's own
    shorthand-classified compounds into the term(s) their prefix maps to -
    they become one additional member of their own class, which keeps the
    2x2 table consistent (the hit set stays inside the universe) and barely
    moves N against a ~50k reference set.
    """
    if level not in LIPID_CLASS_LEVELS:
        level = "lm_main_class"
    if universe_scope == "database" and crosswalk:
        result = _build_class_universe_from_crosswalk(crosswalk, level, "LIPIDCLASS")
        # Select on class_level_only, NOT on a "SHORTHAND:" key prefix: a
        # shorthand name that RefMet also lists ("PC 34:1") is keyed on RefMet's
        # representative structure, so it has no such prefix - yet that key
        # carries no LIPID MAPS class in the crosswalk either, so it would
        # silently vanish from both the universe and the hit set.
        shorthand_df = matched_df[
            matched_df["class_level_only"].astype(bool)
            & (matched_df["compound_key"].astype(str).str.strip() != "")
        ] if "class_level_only" in matched_df.columns else matched_df.iloc[0:0]
        if not shorthand_df.empty:
            for _, row in shorthand_df.iterrows():
                class_label = _clean(row.get(level))
                if not class_label:
                    continue
                term_id = f"LIPIDCLASS:{class_label}"
                info = result["terms"].setdefault(term_id, {"name": class_label, "genes": set()})
                info["genes"].add(str(row["compound_key"]))
            universe = set()
            for info in result["terms"].values():
                universe |= info["genes"]
            result["universe_size"] = len(universe)
    else:
        result = _build_class_universe(matched_df, level, "LIPIDCLASS")
    if not result["terms"]:
        result["caveat"] = (
            "No LIPID MAPS classification matched these compounds - expected unless this "
            "is a lipidomics dataset."
        )
    return result


def _pathway_capable_keys(matched_df: pd.DataFrame) -> set:
    """Compound keys eligible for PATHWAY membership.

    Excludes class-level-only (lipid-shorthand) rows: "PC 34:1" names a bulk
    species covering many structures, so asserting it participates in a
    specific reaction would be a claim the name cannot support.
    """
    df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
    if "class_level_only" in df.columns:
        df = df[~df["class_level_only"].astype(bool)]
    return set(df["compound_key"].astype(str))


def _id_to_key_maps(source, scope: str = "experiment") -> "tuple[dict, dict]":
    """{chebi_id: compound_key} and {kegg_id: compound_key}, so a pathway
    file keyed on either can be translated into this module's single
    identifier space.

    scope="experiment": `source` is matched_df, and the maps cover ONLY this
    run's own measured+eligible compounds - correct for the concentration-
    table tile, where the universe must stay narrowed to what was actually
    quantified (see the universe builders' own section comment).

    scope="database": `source` is the crosswalk, and the maps cover EVERY
    compound the crosswalk itself knows (~80k). Correct for the list-input
    tiles: a plain compound list has no "measured" concept, so the pathway
    universe must be the database's own full reference set, translated
    through the WHOLE crosswalk rather than restricted to the user's own
    upload - exactly like testing a gene list against all of KEGG elsewhere
    in this app, not just genes coincidentally already in the list.
    """
    if scope == "database":
        crosswalk = source
        return dict(crosswalk.get("by_chebi", {})), dict(crosswalk.get("by_kegg", {}))
    matched_df = source
    chebi_to_key: dict = {}
    kegg_to_key: dict = {}
    df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
    if "class_level_only" in df.columns:
        df = df[~df["class_level_only"].astype(bool)]
    for _, row in df.iterrows():
        key = str(row["compound_key"])
        chebi = _clean(row.get("chebi_id"))
        kegg = _clean(row.get("kegg_id"))
        if chebi:
            chebi_to_key.setdefault(chebi, key)
        if kegg:
            kegg_to_key.setdefault(kegg, key)
    return chebi_to_key, kegg_to_key


def _crosswalk_id_maps(crosswalk: dict) -> "tuple[dict, dict, dict]":
    """{chebi_id: key}, {kegg_id: key}, {skeleton: key} from the WHOLE
    crosswalk - the database-scope counterpart of the matched_df-derived
    skeleton_to_key map build_reactome_compound_universe uses under
    experiment scope. Used for the ChEBI-skeleton bridge (see that
    function's own docstring) when the universe itself is database-scoped.
    """
    return (
        dict(crosswalk.get("by_chebi", {})),
        dict(crosswalk.get("by_kegg", {})),
        dict(crosswalk.get("by_skeleton", {})),
    )


def build_kegg_organism_pathway_index(organism_code: str, verbose: bool = False,
                                      force_refresh: bool = False) -> dict:
    """Which pathway MAPS one KEGG organism has:
    {"org", "names": {5-digit number: the organism's own pathway name},
     "universe_size", "fetched_at", "caveat"}.

    Verified live against KEGG (2026-09-24): a map's compound content does not
    vary by organism - the KGML compound sets of hsa00020 and sce00020 are
    identical and equal to reference map00020, and KEGG serves no
    organism-specific compound links at all (link/compound/hsa is empty). So
    the organism can only select WHICH maps are tested, never which compounds
    a map holds. That is why this is a small per-organism index beside one big
    shared link file, not a per-organism re-download.
    """
    org = str(organism_code or "").strip().lower()
    if not org:
        return {"org": "", "names": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": "No KEGG organism code was provided."}

    cache_db = f"Metabolomics_KEGG_orgmaps_v{KEGG_ORG_MAPS_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, org)
        if cached is not None and cached.get("names"):
            return cached

    status, text = _kegg_safe_get(KEGG_LIST_PATHWAY_URL.format(org=org), verbose=verbose)
    names: dict = {}
    if status == "ok" and text:
        # The number is the trailing 5 digits: organism codes are NOT always 3
        # characters, so slicing off a fixed prefix would be wrong.
        for pathway_id, name in _parse_kegg_pathway_list(text).items():
            match = re.search(r"(\d{5})$", str(pathway_id).strip())
            if match:
                names[match.group(1)] = name

    caveat = None
    if status == "error":
        caveat = f"Could not reach KEGG to list the pathway maps for organism '{org}'."
    elif not names:
        caveat = f"KEGG lists no pathway maps for organism '{org}'."

    result = {"org": org, "names": names, "universe_size": len(names),
              "fetched_at": _now_iso(), "caveat": caveat}
    # Cache only a complete, non-empty answer - never a transient failure.
    if names and status == "ok":
        _save_universe_cache(cache_db, org, result)
    return result


def _metabolomics_kegg_note(kegg_org: dict, taxon_id: str, taxon_name: str) -> str:
    """Re-word the shared resolver's failure note for Metabolomics.

    resolve_kegg_organism_code is shared with the gene modules, whose notes end
    "GO and MSigDB match on gene symbol and are unaffected" / "try GO or MSigDB
    for this organism" - advice that means nothing here and points at analyses
    that do not exist in this module. The resolver is left untouched (other
    modules and their tests depend on it); its STRUCTURE (candidates) tells us
    which case this is, and the wording is rebuilt to point at the KEGG organism
    box instead.
    """
    note = str(kegg_org.get("note") or "")
    candidates = kegg_org.get("candidates") or []
    if candidates:
        shown = ", ".join(f"{cand['organism_code']} ({cand['organism_name']})" for cand in candidates[:5])
        more = f", and {len(candidates) - 5} more" if len(candidates) > 5 else ""
        return (
            f'"{taxon_name}" matches {len(candidates)} KEGG organisms ({shown}{more}), none of which carries '
            f"NCBI taxon {taxon_id} - typically KEGG has strain-level genomes where this run's organism is a "
            "species. Pick one in the KEGG organism box (search KEGG's organisms by name or code); until then "
            "KEGG is skipped."
        )
    if "has no genome for" in note:
        return (
            f'KEGG has no genome for "{taxon_name}". KEGG covers only organisms with a sequenced, annotated '
            "genome; you can pick a close relative in the KEGG organism box (search KEGG's organisms by name "
            "or code)."
        )
    return note


def _resolve_kegg_for_metabolomics(
    taxon_id: str,
    taxon_name: str,
    kegg_organism_code: str = "",
    kegg_organism_name: str = "",
    verbose: bool = False,
    force_refresh: bool = False,
) -> "tuple[dict, str | None, str | None]":
    """Resolve KEGG's organism and fetch that organism's pathway-map index.

    Two ways to name the organism:
      - EXPLICIT: the user picked a KEGG organism code from KEGG's own list
        (kegg_organism_code). Used as given - nothing to derive or guess. For
        Metabolomics this is safe in a way it is not for genes: there is no
        gene-ID mapping that could be resolved for one organism and tested
        against another (the divergence bug of commit 27ce033), because
        compounds are organism-independent; the organism only selects maps.
      - DERIVED: otherwise the run's own organism (taxon_id/taxon_name) is
        matched to a KEGG genome by the same resolver Transcriptomics and
        Proteomics use (Transcriptomics_module.resolve_kegg_organism_code).
        When that cannot be verified to one genome there is nothing to test,
        and the note says the user can pick a KEGG organism explicitly.

    Returns (org_index, caveat, note), the same contract as
    Proteomics_module._resolve_kegg_for_proteomics:
      - note set (and org_index == {}) means nothing can be tested, and says
        WHICH reason - no organism name, no KEGG genome, ambiguous, KEGG
        unreachable, or no maps listed. Never collapsed into a generic "no
        overlap" a reader would take for a biology conclusion. `caveat` is then
        the organism-resolution caveat ONLY; the caller stores both side by side.
      - note None means success and `caveat` is the organism-resolution caveat
        (e.g. an unverified taxon match), to be combined by the caller.
    There is deliberately NO fallback to KEGG's generic reference maps: a run
    that silently tests pathways the organism does not have is the exact
    behaviour this function replaces.
    """
    explicit = str(kegg_organism_code or "").strip().lower()
    if explicit:
        org_caveat = None
        code = explicit
        organism_name = str(kegg_organism_name or "").strip() or explicit
        label = f'"{organism_name}" ({code})' if organism_name != code else f"({code})"
    else:
        kegg_org = resolve_kegg_organism_code(taxon_id, taxon_name, verbose=verbose, force_refresh=force_refresh)
        org_caveat = kegg_org.get("caveat")
        if kegg_org.get("note"):
            return {}, org_caveat, _metabolomics_kegg_note(kegg_org, taxon_id, taxon_name)
        code = str(kegg_org.get("organism_code") or "").strip().lower()
        if not code:
            return {}, org_caveat, (
                "KEGG's organism could not be resolved for this run, so KEGG's terms are unavailable. "
                "Pick a KEGG organism explicitly."
            )
        organism_name = kegg_org.get("organism_name") or taxon_name
        label = f'"{organism_name}" ({code})'

    index = build_kegg_organism_pathway_index(code, verbose=verbose, force_refresh=force_refresh)
    if not index.get("names"):
        if str(index.get("caveat") or "").startswith("Could not reach"):
            note = (f"Could not reach KEGG to list the pathway maps for {label} in this run "
                    "(transient service error); KEGG's terms are unavailable. Re-run the analysis to retry.")
        elif explicit:
            note = f"KEGG lists no pathway maps for organism {label}, so there is nothing to test."
        else:
            note = f"KEGG has a genome for {label} but lists no pathway maps for it, so there is nothing to test."
        return {}, org_caveat, note
    index["organism_name"] = organism_name
    index["explicit"] = bool(explicit)
    return index, org_caveat, None


def build_kegg_compound_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                                 universe_scope: str = "experiment",
                                 taxon_id: str = "", taxon_name: str = "",
                                 kegg_organism_code: str = "", kegg_organism_name: str = "",
                                 verbose: bool = False,
                                 force_refresh: bool = False) -> dict:
    """KEGG compound -> pathway membership, for the run's ORGANISM.

    KEGG's compound-pathway links are its REFERENCE maps, fetched once globally
    and cached under a single "global" key - that file is organism-independent
    and expensive. The organism then selects which of those maps are tested:
    only maps the organism actually has (KEGG's list/pathway/<org>) are kept,
    and terms are reported under the organism's own ids and names
    (hsa00020, "Citrate cycle (TCA cycle)"). Compound membership WITHIN a map
    is KEGG's reference biochemistry and does not vary by organism (see
    build_kegg_organism_pathway_index), so this is stated in the caveat rather
    than implying species-specific compound content.

    If KEGG has no genome for the organism, or the organism is ambiguous and
    unconfirmed, the result is EMPTY with a `note` naming the reason - never a
    silent fallback to the generic reference maps.
    """
    org_index, org_caveat, note = _resolve_kegg_for_metabolomics(
        taxon_id, taxon_name, kegg_organism_code, kegg_organism_name,
        verbose=verbose, force_refresh=force_refresh,
    )
    if note:
        # The cheap organism check gates the expensive download.
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": org_caveat, "note": note}

    cache_db = f"Metabolomics_KEGG_v{UNIVERSE_CACHE_VERSION}"
    universe = None if force_refresh else _load_universe_cache(cache_db, "global")
    if universe is None:
        status, text = _kegg_safe_get("https://rest.kegg.jp/link/pathway/compound", verbose=verbose)
        if status != "ok" or not text:
            return {
                "terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": org_caveat,
                "note": "Could not reach KEGG to build the compound-pathway map for this run. "
                        "Re-run the analysis to retry.",
            }
        status_names, names_text = _kegg_safe_get("https://rest.kegg.jp/list/pathway", verbose=verbose)
        pathway_names: dict = {}
        if status_names == "ok" and names_text:
            for line in names_text.splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    pathway_names[parts[0].replace("path:", "").strip()] = parts[1].strip()

        links: dict = {}
        for line in text.splitlines():
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            compound = _normalize_kegg_compound(parts[0])
            pathway = parts[1].replace("path:", "").strip()
            if not compound or not pathway:
                continue
            # "map#####" is the reference map; "rn#####"/organism-prefixed
            # duplicates of the same map would double-count the same pathway.
            if not pathway.startswith("map"):
                continue
            # Overview/chemical-structure maps aggregate compounds across huge
            # swaths of chemistry rather than naming a real pathway - see
            # KEGG_COMPOUND_EXCLUDED_MAPS' own comment for the concrete case
            # this fixes (a global map outranking the correct TCA cycle hit).
            # Applied HERE, on the reference id: after organism re-keying
            # ("hsa01100") the "map" prefix this test strips would no longer be
            # there, and every overview map would silently come back.
            if pathway.replace("map", "") in KEGG_COMPOUND_EXCLUDED_MAPS:
                continue
            links.setdefault(pathway, set()).add(compound)

        universe = {
            "links": {p: sorted(c) for p, c in links.items()},
            "names": pathway_names,
            "fetched_at": _now_iso(),
        }
        if links:
            _save_universe_cache(cache_db, "global", universe)

    org = org_index["org"]
    org_names = org_index["names"]
    source = crosswalk if (universe_scope == "database" and crosswalk) else matched_df
    _, kegg_to_key = _id_to_key_maps(source, scope=universe_scope)
    ref_links = universe.get("links") or {}
    terms: dict = {}
    n_ref_kept = 0
    n_with_compounds = 0
    for pathway, compounds in ref_links.items():
        number = pathway.replace("map", "")
        if number not in org_names:
            continue                      # the organism does not have this map
        n_ref_kept += 1
        keys = {kegg_to_key[c] for c in compounds if c in kegg_to_key}
        if keys:
            n_with_compounds += 1
            name = org_names.get(number) or (universe.get("names") or {}).get(pathway) or org + number
            terms[org + number] = {"name": name, "genes": keys}
    covered = set()
    for info in terms.values():
        covered |= info["genes"]

    organism_label = f'{org_index.get("organism_name") or taxon_name} ({org})'
    result = {
        "terms": terms, "universe_size": len(covered),
        "fetched_at": universe.get("fetched_at", _now_iso()),
        "caveat": org_caveat,
    }
    if not terms:
        # Two opposite diagnoses that must not be reported as one:
        any_kegg_compound = any(c in kegg_to_key for cs in ref_links.values() for c in cs)
        result["note"] = (
            "No measured compound carried a KEGG compound ID."
            if not any_kegg_compound else
            f"Your compounds map to KEGG compounds, but none of the KEGG maps that contain them "
            f"exist for {organism_label}."
        )
        return result
    coverage = (
        f"KEGG organism: {organism_label}{' (chosen explicitly)' if org_index.get('explicit') else ''} - "
        f"{n_ref_kept} of {len(ref_links)} reference maps "
        f"exist for this organism and were tested. Which maps are tested depends on the organism; "
        f"compound membership WITHIN a map is KEGG's reference biochemistry and does not vary by organism."
    )
    result["caveat"] = (org_caveat + " " + coverage) if org_caveat else coverage
    return result


CHEBI_BACKEND_COMPOUNDS_URL = "https://www.ebi.ac.uk/chebi/backend/api/public/compounds"
CHEBI_BACKEND_BATCH_SIZE = 150
CHEBI_SKELETON_CACHE_VERSION = 1


def _reactome_referenced_chebi_ids(path: str) -> set:
    """Every distinct bare ChEBI ID Reactome's file references, any species.

    Verified live 2026-09-22: this is a small, CLOSED set (~3,260 ids) -
    Reactome's small-molecule annotation coverage is sparse, not organism-
    specific, so scanning once (not per-run, not per-organism) is enough.
    """
    ids: set = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            chebi = _normalize_chebi(line.split("\t", 1)[0])
            if chebi:
                ids.add(chebi)
    return ids


def _fetch_chebi_skeletons(chebi_ids: set, verbose: bool = False,
                           force_refresh: bool = False, progress_callback=None) -> dict:
    """{chebi_id: InChIKey skeleton} for the given ChEBI IDs, via ChEBI's
    REST API on www.ebi.ac.uk (reachable) rather than its bulk FTP files
    (ftp.ebi.ac.uk, verified unreachable 2026-09-21).

    WHY THIS EXISTS: Reactome's ChEBI2Reactome file and this module's own
    crosswalk (built from RefMet/LIPID MAPS) frequently reference DIFFERENT
    ChEBI IDs for what is chemically the same compound at a different
    protonation state - e.g. Reactome's CHEBI:16947 "citrate(3-)" versus
    RefMet's CHEBI:30769 "citric acid" for the same "Citrate" upload. Exact
    ChEBI ID matching alone therefore misses real matches systematically, not
    occasionally (verified live: only 7 of 45 pathway-capable compounds in a
    realistic run reached Reactome by exact ID). The InChIKey SKELETON (first
    block: connectivity only, ignoring the protonation/stereochemistry blocks
    that differ between these forms) is identical across them and bridges the
    gap - the exact same principle already used for cross-source compound
    matching (see inchikey_skeleton's docstring), applied here to a second
    identifier space.

    Batches requests (the endpoint accepts chebi_ids=id1,id2,...; ~150/request
    keeps individual payloads reasonable) and caches the result like a
    universe, since this set is bounded and organism-agnostic - built once,
    reused by every future Reactome-enabled run regardless of organism.
    """
    cache_db = f"Metabolomics_ChEBISkeletons_v{CHEBI_SKELETON_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, "global")
        if cached is not None and cached.get("skeletons"):
            have = set(cached["skeletons"].keys())
            if chebi_ids <= have:
                return cached["skeletons"]

    if HTTP_SESSION is None:
        return {}
    ids_list = sorted(chebi_ids)
    skeletons: dict = {}
    total_batches = max(1, (len(ids_list) + CHEBI_BACKEND_BATCH_SIZE - 1) // CHEBI_BACKEND_BATCH_SIZE)
    for i in range(0, len(ids_list), CHEBI_BACKEND_BATCH_SIZE):
        batch = ids_list[i:i + CHEBI_BACKEND_BATCH_SIZE]
        if progress_callback:
            try:
                progress_callback(i // CHEBI_BACKEND_BATCH_SIZE, total_batches,
                                  "Bridging Reactome's ChEBI IDs...")
            except Exception:
                pass
        try:
            resp = HTTP_SESSION.get(CHEBI_BACKEND_COMPOUNDS_URL,
                                    params={"chebi_ids": ",".join(batch)},
                                    headers=OMISCERA_HEADERS, timeout=60)
            if resp.status_code != 200:
                continue
            payload = resp.json()
            for chebi_id, entry in payload.items():
                data = (entry or {}).get("data") or {}
                struct = data.get("default_structure") or {}
                ikey = _clean(struct.get("standard_inchi_key")).upper()
                if ikey:
                    skeletons[_normalize_chebi(chebi_id)] = inchikey_skeleton(ikey)
        except Exception as exc:
            if verbose:
                print(f"[metabolomics] ChEBI skeleton batch failed: {exc}")
            continue
    if skeletons:
        result = {"skeletons": skeletons, "fetched_at": _now_iso()}
        # Merge with whatever was already cached rather than overwriting -
        # a partial batch failure must not discard previously-fetched ids.
        prior = _load_universe_cache(cache_db, "global") or {}
        merged = dict(prior.get("skeletons") or {})
        merged.update(skeletons)
        result["skeletons"] = merged
        _save_universe_cache(cache_db, "global", result)
        return merged
    return skeletons


def build_reactome_compound_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                                     universe_scope: str = "experiment",
                                     taxon_name: str = "",
                                     verbose: bool = False, force_refresh: bool = False,
                                     progress_callback=None) -> dict:
    """Reactome pathway membership from ChEBI2Reactome_All_Levels.txt.

    File shape verified live 2026-09-21 - tab-separated, no header:
        0 ChEBI ID (bare) | 1 pathway stable id | 2 url | 3 pathway name
        4 evidence code   | 5 species name
    Column-for-column the same shape Genomics_module.build_reactome_universe
    already parses for UniProt2Reactome, with a ChEBI ID in place of the
    accession.

    Matches a file row to a measured compound two ways: first by EXACT ChEBI
    ID (cheap, and correct whenever Reactome and our crosswalk happen to cite
    the same protonation state), then - for rows that miss - by InChIKey
    SKELETON via _fetch_chebi_skeletons, which is the majority case in
    practice (see that function's docstring). Skeleton-bridged matches are
    counted separately in the returned dict, following the same
    "report the fallback tier honestly" rule as match_compounds' own
    n_skeleton_only, so the bridge's contribution is visible, not silent.
    """
    path = _ensure_bulk_file(CHEBI2REACTOME_URL, CHEBI2REACTOME_LOCAL,
                             verbose=verbose, force_refresh=force_refresh)
    if not path:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": "Could not fetch Reactome's ChEBI-to-pathway mapping for this run."}

    if universe_scope == "database" and crosswalk:
        chebi_to_key, _, skeleton_to_key = _crosswalk_id_maps(crosswalk)
    else:
        chebi_to_key, _ = _id_to_key_maps(matched_df, scope="experiment")
        # Build {skeleton: compound_key} from the measured compounds' own
        # InChIKeys, so a Reactome row can be bridged even when its exact
        # ChEBI ID was never assigned to any measured compound.
        skeleton_to_key = {}
        df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
        if "class_level_only" in df.columns:
            df = df[~df["class_level_only"].astype(bool)]
        for _, row in df.iterrows():
            skel = inchikey_skeleton(_clean(row.get("inchikey")))
            if skel:
                skeleton_to_key.setdefault(skel, str(row["compound_key"]))

    reactome_chebi_ids = _reactome_referenced_chebi_ids(path)
    reactome_skeletons = (
        _fetch_chebi_skeletons(reactome_chebi_ids, verbose=verbose,
                               force_refresh=force_refresh, progress_callback=progress_callback)
        if skeleton_to_key else {}
    )

    wanted_species = _clean(taxon_name).lower()
    terms: dict = {}
    species_seen = False
    bridged_keys: set = set()
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                chebi = _normalize_chebi(parts[0])
                if not chebi:
                    continue
                key = chebi_to_key.get(chebi)
                bridged = False
                if key is None:
                    skel = reactome_skeletons.get(chebi)
                    if skel and skel in skeleton_to_key:
                        key = skeleton_to_key[skel]
                        bridged = True
                if key is None:
                    continue
                species = _clean(parts[5]).lower()
                if wanted_species:
                    if species != wanted_species:
                        continue
                    species_seen = True
                pathway_id, pathway_name = _clean(parts[1]), _clean(parts[3])
                if not pathway_id:
                    continue
                info = terms.setdefault(pathway_id, {"name": pathway_name or pathway_id, "genes": set()})
                info["genes"].add(key)
                if bridged:
                    bridged_keys.add(key)
    except Exception as exc:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": f"Could not parse Reactome's ChEBI mapping: {exc}"}
    n_skeleton_bridged = len(bridged_keys)

    covered = set()
    for info in terms.values():
        covered |= info["genes"]
    caveat = None
    if wanted_species and not species_seen and not terms:
        caveat = (
            f"No Reactome small-molecule annotations matched organism '{taxon_name}'. "
            "Check the organism selection, or that these compounds carry ChEBI IDs."
        )
    elif not terms:
        caveat = "No measured compound reached a Reactome pathway (none carried a ChEBI ID)."
    elif n_skeleton_bridged:
        caveat = (
            f"{n_skeleton_bridged} of {len(covered)} compound(s) reached Reactome only via "
            "InChIKey-skeleton bridging (Reactome cites a different protonation/charge state "
            "of the same compound than this run's own identifiers) - see the README."
        )
    return {"terms": terms, "universe_size": len(covered),
            "fetched_at": _now_iso(), "caveat": caveat, "n_skeleton_bridged": n_skeleton_bridged}


def _parse_pathway_csv_zip(path: str, chebi_to_key: dict, kegg_to_key: dict,
                           taxon_name: str = "", verbose: bool = False) -> "tuple[dict, bool]":
    """Shared parser for SMPDB's and PathBank's metabolite CSVs.

    Both publish the SAME shape (a zipped CSV whose header names a pathway id,
    a pathway name, and per-metabolite KEGG/ChEBI columns), so one parser
    serves both rather than two near-identical copies. Column positions are
    resolved BY HEADER NAME, never by index, because the two files do not
    agree on column order and PathBank's exact spelling was not verifiable
    ahead of time. Returns (terms, saw_species_column).
    """
    terms: dict = {}
    saw_species = False
    with zipfile.ZipFile(path) as zf:
        members = [n for n in zf.namelist() if n.lower().endswith(".csv")]
        if not members:
            return terms, saw_species
        wanted_species = _clean(taxon_name).lower()
        for member in members:
            with zf.open(member) as raw:
                stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                reader = csv.DictReader(stream)
                fieldnames = [f.strip().lower() for f in (reader.fieldnames or [])]

                def _find(*candidates):
                    for cand in candidates:
                        for i, name in enumerate(fieldnames):
                            if name == cand:
                                return reader.fieldnames[i]
                    for cand in candidates:
                        for i, name in enumerate(fieldnames):
                            if cand in name:
                                return reader.fieldnames[i]
                    return None

                col_pw_id = _find("pathbank id", "smpdb id", "pathway id", "pw id")
                col_pw_name = _find("pathway name", "name")
                col_kegg = _find("kegg id", "kegg")
                col_chebi = _find("chebi id", "chebi")
                col_species = _find("species", "organism")
                if col_species:
                    saw_species = True
                if not col_pw_id:
                    continue
                for row in reader:
                    if wanted_species and col_species:
                        if _clean(row.get(col_species)).lower() != wanted_species:
                            continue
                    pathway_id = _clean(row.get(col_pw_id))
                    if not pathway_id:
                        continue
                    key = None
                    chebi = _normalize_chebi(row.get(col_chebi)) if col_chebi else ""
                    if chebi and chebi in chebi_to_key:
                        key = chebi_to_key[chebi]
                    if key is None and col_kegg:
                        kegg = _normalize_kegg_compound(row.get(col_kegg))
                        if kegg and kegg in kegg_to_key:
                            key = kegg_to_key[kegg]
                    if key is None:
                        continue
                    pathway_name = _clean(row.get(col_pw_name)) if col_pw_name else ""
                    info = terms.setdefault(pathway_id, {"name": pathway_name or pathway_id, "genes": set()})
                    info["genes"].add(key)
    return terms, saw_species


def build_smpdb_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                         universe_scope: str = "experiment", taxon_id: str = "",
                         verbose: bool = False, force_refresh: bool = False) -> dict:
    """SMPDB pathway membership. Human-only resource, by design."""
    if _clean(taxon_id) and _clean(taxon_id) != SMPDB_HUMAN_TAXON_ID:
        return {
            "terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
            "caveat": (
                f"SMPDB is a human-only resource (taxon {SMPDB_HUMAN_TAXON_ID}); it has no "
                f"data for taxon {taxon_id}, so this database shows no results for this organism."
            ),
        }
    path = _ensure_bulk_file(SMPDB_METABOLITES_URL, SMPDB_LOCAL,
                             verbose=verbose, force_refresh=force_refresh, timeout=900)
    if not path:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": "Could not fetch SMPDB's metabolite-pathway file for this run."}
    source = crosswalk if (universe_scope == "database" and crosswalk) else matched_df
    chebi_to_key, kegg_to_key = _id_to_key_maps(source, scope=universe_scope)
    try:
        terms, _ = _parse_pathway_csv_zip(path, chebi_to_key, kegg_to_key, verbose=verbose)
    except Exception as exc:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": f"Could not parse SMPDB's metabolite file: {exc}"}
    covered = set()
    for info in terms.values():
        covered |= info["genes"]
    return {"terms": terms, "universe_size": len(covered), "fetched_at": _now_iso(),
            "caveat": None if terms else "No measured compound reached an SMPDB pathway."}


def build_pathbank_universe(matched_df: pd.DataFrame, crosswalk: "dict | None" = None,
                            universe_scope: str = "experiment", taxon_name: str = "",
                            primary_only: bool = True, verbose: bool = False,
                            force_refresh: bool = False) -> dict:
    """PathBank pathway membership. Multi-organism, unlike SMPDB.

    primary_only selects the 4.5MB canonical-pathway file over the 52MB
    all-pathways one; canonical is the better default because PathBank's full
    set fragments the same biology across hundreds of thousands of
    near-duplicate pathways, which would dominate an FDR family.
    """
    url = PATHBANK_PRIMARY_URL if primary_only else PATHBANK_ALL_URL
    local = PATHBANK_PRIMARY_LOCAL if primary_only else PATHBANK_ALL_LOCAL
    path = _ensure_bulk_file(url, local, verbose=verbose, force_refresh=force_refresh, timeout=900)
    if not path:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": "Could not fetch PathBank's metabolite-pathway file for this run."}
    source = crosswalk if (universe_scope == "database" and crosswalk) else matched_df
    chebi_to_key, kegg_to_key = _id_to_key_maps(source, scope=universe_scope)
    try:
        terms, saw_species = _parse_pathway_csv_zip(
            path, chebi_to_key, kegg_to_key, taxon_name=taxon_name, verbose=verbose)
    except Exception as exc:
        return {"terms": {}, "universe_size": 0, "fetched_at": _now_iso(),
                "caveat": f"Could not parse PathBank's metabolite file: {exc}"}
    covered = set()
    for info in terms.values():
        covered |= info["genes"]
    caveat = None
    if taxon_name and not saw_species:
        caveat = ("PathBank's file exposed no species column, so its pathways were used "
                  "without an organism filter.")
    elif not terms:
        caveat = "No measured compound reached a PathBank pathway."
    return {"terms": terms, "universe_size": len(covered), "fetched_at": _now_iso(),
            "caveat": caveat}


# ---------------------------------------------------------------------------
# Shared analysis plumbing
# ---------------------------------------------------------------------------

def _build_universes(matched_df: pd.DataFrame, crosswalk: "dict | None", databases: "list[str]",
                     universe_scope: str = "experiment",
                     taxon_id: str = "", taxon_name: str = "",
                     chemical_class_level: str = "sub_class",
                     lipid_class_level: str = "lm_main_class",
                     pathbank_primary_only: bool = True,
                     kegg_organism_code: str = "", kegg_organism_name: str = "",
                     force_refresh: bool = False, verbose: bool = False,
                     progress_callback=None) -> dict:
    """Build every requested universe once, so ORA and Quant share the work.

    universe_scope: "experiment" (default, concentration-table tile) or
    "database" (list-input tiles) - see the universe builders' own section
    comment for what each means. Passed straight through to every builder.

    Returns {database_label: universe_dict}. A database that could not be
    built is still present, carrying empty terms and its own caveat - never
    omitted, so the results screen can say why rather than silently showing
    one fewer table. Under "database" scope, a successfully-built block's
    caveat gains a scope statement naming its universe explicitly, so a
    block never leaves a reader guessing which rule produced its numbers.
    """
    universes: dict = {}
    total = max(len(databases), 1)
    for i, database in enumerate(databases):
        if progress_callback:
            try:
                progress_callback(i, total, f"Building {database} universe...")
            except Exception:
                pass
        if database == "Chemical Class":
            universes[database] = build_chemical_class_universe(
                matched_df, crosswalk=crosswalk, level=chemical_class_level,
                universe_scope=universe_scope)
        elif database == "Lipid Class":
            universes[database] = build_lipid_class_universe(
                matched_df, crosswalk=crosswalk, level=lipid_class_level,
                universe_scope=universe_scope)
        elif database == "KEGG":
            universes[database] = build_kegg_compound_universe(
                matched_df, crosswalk=crosswalk, universe_scope=universe_scope,
                taxon_id=taxon_id, taxon_name=taxon_name,
                kegg_organism_code=kegg_organism_code, kegg_organism_name=kegg_organism_name,
                verbose=verbose, force_refresh=force_refresh)
        elif database == "Reactome":
            def _reactome_sub_progress(step, total_steps, label, _i=i, _total=total, _db=database):
                if progress_callback:
                    try:
                        progress_callback(_i, _total, f"{_db}: {label}")
                    except Exception:
                        pass
            universes[database] = build_reactome_compound_universe(
                matched_df, crosswalk=crosswalk, universe_scope=universe_scope,
                taxon_name=taxon_name, verbose=verbose, force_refresh=force_refresh,
                progress_callback=_reactome_sub_progress)
        elif database == "SMPDB":
            universes[database] = build_smpdb_universe(
                matched_df, crosswalk=crosswalk, universe_scope=universe_scope,
                taxon_id=taxon_id, verbose=verbose, force_refresh=force_refresh)
        elif database == "PathBank":
            universes[database] = build_pathbank_universe(
                matched_df, crosswalk=crosswalk, universe_scope=universe_scope,
                taxon_name=taxon_name, primary_only=pathbank_primary_only,
                verbose=verbose, force_refresh=force_refresh)
        else:
            continue

        block = universes[database]
        if block.get("terms") and universe_scope == "database":
            scope_note = (
                f"Universe: all {block['universe_size']} compounds annotated in {database} "
                f"(the database's own reference set) - your list contributes the tested subset."
            )
            block["caveat"] = (block["caveat"] + " " + scope_note) if block.get("caveat") else scope_note
    return universes


def _finalize_member_column(tables: dict, key_to_name: dict, source_col: str,
                            target_columns: "list[str]") -> None:
    """Relabel the shared engines' member column to the user's own names.

    The Genomics engines name this column for genes ("Overlapping Genes" in
    the Fisher engine, "Genes" in the Wilcoxon one) and fill it with whatever
    identifier was tested - here an internal compound key (an InChIKey, or a
    namespaced fallback), which is unreadable in a results table. One pass
    over every stored block, rather than at each of the several places a
    block is created, so no path can be missed. Mirrors
    Proteomics_module._finalize_member_column.
    """
    def _display(member: str) -> str:
        m = member.strip()
        return key_to_name.get(m, m) if m else m

    for block in tables.values():
        tbl = block.get("table")
        if tbl is None:
            continue
        if source_col in tbl.columns:
            if not tbl.empty:
                tbl[source_col] = tbl[source_col].map(
                    lambda cell: "; ".join(
                        sorted({_display(x) for x in str(cell).split(",") if x.strip()})
                    )
                )
            tbl = tbl.rename(columns=_BURDEN_COLUMN_RENAMES if source_col == "Genes"
                             else {source_col: OVERLAP_COLUMN})
        else:
            tbl = tbl.rename(columns=_BURDEN_COLUMN_RENAMES)
        block["table"] = tbl.reindex(columns=target_columns)


def _key_to_name_map(matched_df: pd.DataFrame) -> dict:
    """{compound_key: the user's own Name}. First row wins for a duplicated
    key so the displayed name is stable and input-order independent."""
    mapping: dict = {}
    for _, row in matched_df.iterrows():
        key = _clean(row.get("compound_key"))
        name = _clean(row.get("Name"))
        if key and name:
            mapping.setdefault(key, name)
    return mapping


def _eligible_keys(matched_df: pd.DataFrame, database: str) -> set:
    """Which compound keys may participate in THIS database.

    Class databases accept class-level-only rows; pathway databases do not -
    see _pathway_capable_keys.
    """
    if database in ("Chemical Class", "Lipid Class"):
        df = matched_df[matched_df["compound_key"].astype(str).str.strip() != ""]
        return set(df["compound_key"].astype(str))
    return _pathway_capable_keys(matched_df)


def _drop_zero_overlap(term_to_genes: dict, sig_set: set, universe_genes: set,
                       min_term_size: int, max_term_size) -> "tuple[dict, int, int]":
    """Exclude terms with ZERO overlap with the significant set before
    Fisher/BH - clusterProfiler convention, and the same filter
    Proteomics_module applies for the same reason: keeping k=0 rows inflates
    the BH family size m for every real test beside them, which costs real
    discoveries. Replicates the size gate so the two counts partition
    cleanly with no double count."""
    filtered, n_excluded, n_no_overlap = {}, 0, 0
    for term_id, info in term_to_genes.items():
        term_genes = info["genes"] & universe_genes
        K = len(term_genes)
        if K < min_term_size or (max_term_size and K > max_term_size):
            n_excluded += 1
            continue
        if term_genes & sig_set:
            filtered[term_id] = info
        else:
            n_no_overlap += 1
    return filtered, n_excluded, n_no_overlap


# ---------------------------------------------------------------------------
# Over-Representation Analysis (ORA), Chemical Class Enrichment, and
# Lipid Class / Structural Enrichment - all three run through THIS one
# function, called with a different `databases` list per tile
# (METABOLOMICS_ORA_DATABASES' pathway entries for ORA; ["Chemical Class"]
# / ["Lipid Class"] for the other two). They share one function because they
# share every piece of plumbing (_build_universes, _eligible_keys,
# _drop_zero_overlap, _finalize_member_column) and differ only in WHICH
# database(s) get passed in - three separate near-identical copies would be
# pure duplication.
#
# Input is a plain compound-name LIST (no statistics at all). Unlike
# Proteomics/Transcriptomics ORA there is no separate "significant" subset
# to select: every compound the user uploaded IS the test set, the same
# convention a plain gene-list ORA already uses elsewhere in this app. The
# universe is therefore the DATABASE's own full annotated set
# (universe_scope="database" throughout) - see the universe builders'
# section comment for why this differs from the concentration-table tile's
# "measured compounds" rule.
# ---------------------------------------------------------------------------

def run_ora(
    matched_df: pd.DataFrame,
    crosswalk: dict,
    databases: "list[str]",
    taxon_id: str = "",
    taxon_name: str = "",
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    chemical_class_level: str = "sub_class",
    lipid_class_level: str = "lm_main_class",
    pathbank_primary_only: bool = True,
    background_compounds: "set | None" = None,
    kegg_organism_code: str = "",
    kegg_organism_name: str = "",
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """ORA for the requested databases against an UPLOADED COMPOUND LIST.

    One result family per database (no direction split - a plain list has
    no up/down). background_compounds, if given, narrows the database's own
    reference universe to compounds the user also uploaded as a background
    set (e.g. "everything on our LC-MS panel", to test enrichment relative
    to what could have been detected rather than every compound the
    database has ever heard of).
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universes = _build_universes(
        matched_df, crosswalk, databases, universe_scope="database",
        taxon_id=taxon_id, taxon_name=taxon_name,
        chemical_class_level=chemical_class_level, lipid_class_level=lipid_class_level,
        pathbank_primary_only=pathbank_primary_only,
        kegg_organism_code=kegg_organism_code, kegg_organism_name=kegg_organism_name,
        force_refresh=force_refresh,
        verbose=verbose, progress_callback=progress_callback,
    )

    key_to_name = _key_to_name_map(matched_df)
    background = {str(b).strip() for b in background_compounds} if background_compounds else None

    steps_total = max(len(databases), 1)
    step = 0
    tables: dict = {}

    def _empty_block(note: str, caveat: "str | None" = None) -> dict:
        return {
            "table": pd.DataFrame(columns=METABOLOMICS_ORA_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0, "n_terms_no_overlap": 0,
            "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
            "note": note, "caveat": caveat,
        }

    for database in databases:
        step += 1
        _progress(step, steps_total, database)
        universe = universes.get(database) or {}
        terms = universe.get("terms") or {}
        eligible = _eligible_keys(matched_df, database)

        annotated: set = set()
        for info in terms.values():
            annotated |= info["genes"]
        universe_keys = annotated
        if background:
            universe_keys &= background

        if not terms:
            # A universe that explains its own emptiness carries a `note`; the
            # organism-resolution caveat then travels beside it rather than
            # being swallowed.
            tables[database] = _empty_block(
                universe.get("note") or universe.get("caveat")
                or f"No {database} annotations available for this run.",
                caveat=universe.get("caveat") if universe.get("note") else None)
            continue

        list_set = eligible & universe_keys
        if not list_set:
            tables[database] = _empty_block(
                f"None of your uploaded compounds are annotated in {database}.")
            continue

        filtered, n_excl, n_no_overlap = _drop_zero_overlap(
            terms, list_set, universe_keys, min_term_size, max_term_size)
        block = _run_fisher_for_namespace(
            filtered, list_set, universe_keys, FDR_ALPHA_DEFAULT, True, 0, None)
        block["n_terms_excluded_min_size"] = n_excl
        block["n_terms_no_overlap"] = n_no_overlap
        if universe.get("caveat"):
            block["caveat"] = universe["caveat"]
        tables[database] = block

    _finalize_member_column(tables, key_to_name, "Overlapping Genes",
                            METABOLOMICS_ORA_RESULT_COLUMNS)

    all_keys = set(matched_df[matched_df["compound_key"].astype(str) != ""]["compound_key"].astype(str))
    summary = {
        "n_input_compounds": len(matched_df),
        "n_matched": len(all_keys),
        "n_class_level_only": int(matched_df["class_level_only"].astype(bool).sum())
        if "class_level_only" in matched_df.columns else 0,
        "background_applied": bool(background),
        "universe_scope": "database",
    }
    return {"tables": tables, "summary": summary,
            "universe_caveats": {db: (u or {}).get("caveat") for db, u in universes.items()}}


# ---------------------------------------------------------------------------
# Quantitative Enrichment (competitive Wilcoxon / 1D annotation enrichment)
#
# Input is a CONCENTRATION TABLE (samples in rows, a group/class label
# column, remaining columns numeric compound values) - see
# parse_concentration_table. Omiscera computes the per-compound statistics
# itself (compute_compound_statistics); nothing here trusts a user-supplied
# log2FC/pvalue the way an earlier version of this module wrongly did (that
# convention belongs to RNA-seq/proteomics differential expression, not to
# what a metabolomics researcher actually uploads).
# ---------------------------------------------------------------------------

_GROUP_COLUMN_PATTERN = re.compile(
    r"^(group|class|label|condition|phenotype|treatment)s?$", re.IGNORECASE,
)


def parse_concentration_table(raw_bytes: bytes, orientation: str = "auto",
                              group_column: "str | None" = None) -> "tuple[pd.DataFrame, pd.Series, dict]":
    """Parse an uploaded concentration/peak-intensity table.

    Accepted shape: samples in rows, one sample-ID column, one group/class
    label column, remaining columns numeric compound values. Also accepts
    the TRANSPOSED shape (compounds in rows, samples in columns) common to
    lipidomics exports. orientation="auto" detects this by checking whether
    the first column's values are mostly numeric (a strong transposed
    signal - the first column of a standard table is sample IDs, never
    numeric-looking en masse); meta always reports what was detected
    (or used, if given explicitly) so a caller can offer an override rather
    than silently trust a guess on an ambiguous file - auto-detection on an
    ambiguous concentration table is exactly where a silent wrong guess
    would be worst.

    Returns (values_df: samples x compounds, numeric; groups: Series aligned
    to values_df's index; meta: {"orientation", "group_column",
    "group_levels", "group_counts", "n_samples", "n_compounds",
    "candidate_group_columns"}).
    """
    try:
        text = raw_bytes.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw_bytes.decode("latin1")
    sep = "\t" if ("\t" in text.split("\n", 1)[0]) else ","
    raw_df = pd.read_csv(io.StringIO(text), sep=sep, dtype=str)
    # Spreadsheet exports routinely carry fully blank rows and columns (a stray
    # formatted cell far to the right or below the data). They are not data, so
    # they are dropped up front rather than becoming "Unnamed: 7" / "nan"
    # compounds that the name lookup then wastes time on. Whitespace-only cells
    # count as blank.
    raw_df = raw_df.apply(lambda col: col.str.strip()).replace("", np.nan)
    _blank_named_cols = [
        col for col in raw_df.columns[raw_df.isna().all(axis=0).values]
        if not str(col).strip().lower().startswith("unnamed:")
    ]
    raw_df = raw_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if raw_df.shape[1] < 3:
        raise ValueError(
            "Expected at least a sample ID column, a group/class label column, and one "
            "numeric compound column."
        )

    detected_orientation = orientation
    if orientation == "auto":
        # Both layouts have a non-numeric first column (sample IDs in one,
        # compound names in the other), so first-column numeric-ness cannot
        # tell them apart. What DOES differ is where the group label lives:
        # a COLUMN in the standard layout, a ROW in the transposed one - so
        # a first-column cell that itself reads "Group"/"Class"/"Condition"
        # is the reliable transposed signal. Anything else is treated as
        # standard, and meta reports the choice so a caller can offer an
        # explicit override on an ambiguous file.
        first_col = raw_df.iloc[:, 0].astype(str).str.strip()
        has_group_row = first_col.map(lambda v: bool(_GROUP_COLUMN_PATTERN.match(v))).any()
        if has_group_row:
            detected_orientation = "transposed"
        else:
            # The group is not called "Group"/"Class"/... - decide by STRUCTURE
            # instead of the name, so "Sample type", "Cohort", "Arm" work in
            # either layout. A group label is a line of text with 2-10 levels,
            # each repeated: a ROW of that shape means compounds are in rows
            # (transposed); a COLUMN of that shape (beside the sample IDs) means
            # samples are in rows (standard). A numeric line never qualifies, so
            # the rows of a standard table (an ID, a label, then numbers) do not
            # look like a group row.
            def _looks_categorical(values) -> bool:
                values = pd.Series(values).dropna()
                if len(values) < 4:
                    return False
                if pd.to_numeric(values, errors="coerce").notna().mean() > 0.2:
                    return False
                counts = values.astype(str).str.strip().value_counts()
                return 2 <= len(counts) <= 10 and bool((counts >= 2).all())

            row_like = any(_looks_categorical(raw_df.iloc[r, 1:]) for r in range(len(raw_df)))
            col_like = any(_looks_categorical(raw_df.iloc[:, j]) for j in range(1, raw_df.shape[1]))
            detected_orientation = "transposed" if (row_like and not col_like) else "standard"

    if detected_orientation == "transposed":
        raw_df = raw_df.set_index(raw_df.columns[0]).transpose().reset_index()
        raw_df.columns = ["Sample"] + [str(c) for c in raw_df.columns[1:]]
    else:
        raw_df = raw_df.rename(columns={raw_df.columns[0]: "Sample"})

    def _plausible_group_column(col) -> bool:
        series = raw_df[col]
        # A repeated column name (the same compound measured as two features)
        # selects a whole DataFrame, not a Series - and a group label is one
        # column by definition, so a repeated name is never a group column.
        if isinstance(series, pd.DataFrame):
            return False
        if pd.to_numeric(series, errors="coerce").notna().mean() > 0.8:
            return False
        levels = series.dropna().unique()
        return 2 <= len(levels) <= 10 and bool((series.value_counts() >= 2).all())

    # Offered to the user as choices for an override: columns NAMED like a group
    # label first, then any other column that structurally looks like one
    # (a few levels, each with several samples) - so "Cohort" is offered too.
    candidates = [c for c in raw_df.columns[1:] if _GROUP_COLUMN_PATTERN.match(str(c).strip())]
    candidates += [c for c in raw_df.columns[1:] if c not in candidates and _plausible_group_column(c)]

    # An explicit choice is honoured or rejected - never silently replaced by
    # the auto-detected column, which would hand back a result for a grouping
    # the user did not ask for.
    if group_column:
        if group_column == "Sample" or group_column not in raw_df.columns[1:]:
            raise ValueError(
                f"Group column {group_column!r} was not found among the table's columns "
                f"({', '.join(str(c) for c in raw_df.columns[1:6])}...)."
            )
    chosen_group_col = group_column or (candidates[0] if candidates else None)
    if not chosen_group_col:
        raise ValueError(
            "Could not identify a group/class label column. Name it 'Group', 'Class', "
            "'Condition' or similar, or specify it explicitly."
        )

    groups = raw_df[chosen_group_col].astype(str).str.strip()
    # Select by POSITION, not by label: the same compound is often listed twice
    # (two chromatographic features / adducts of one metabolite), and selecting
    # a repeated label returns every match once per occurrence.
    keep = [i for i, c in enumerate(raw_df.columns) if c not in ("Sample", chosen_group_col)]
    values_all = raw_df.iloc[:, keep].apply(pd.to_numeric, errors="coerce")

    def _is_unnamed(name) -> bool:
        text_ = str(name).strip().lower()
        return text_ in ("", "nan", "none") or text_.startswith("unnamed:")

    # A column with no name, or with no numeric value at all, is not a compound:
    # it cannot be looked up or tested. Dropped and COUNTED (never silently), so
    # the meta can tell the user how many were ignored.
    is_unnamed = [_is_unnamed(raw_df.columns[i]) for i in keep]
    is_empty = [bool(values_all.iloc[:, j].isna().all()) for j in range(values_all.shape[1])]
    n_unnamed_dropped = sum(1 for u, e in zip(is_unnamed, is_empty) if u and not e)
    # Plus the NAMED columns that were entirely blank and already removed above
    # (a compound header with no values at all - an empty compound).
    n_empty_dropped = sum(1 for u, e in zip(is_unnamed, is_empty) if e) + len(_blank_named_cols)
    keep = [i for i, u, e in zip(keep, is_unnamed, is_empty) if not (u or e)]
    if not keep:
        raise ValueError("No compound column with a name and numeric values was found.")
    value_cols = [raw_df.columns[i] for i in keep]
    values_df = raw_df.iloc[:, keep].apply(pd.to_numeric, errors="coerce")
    values_df.index = raw_df["Sample"].astype(str)
    n_repeated_names = len(value_cols) - len(set(value_cols))
    groups.index = values_df.index

    levels = sorted(groups.dropna().unique().tolist())
    if len(levels) < 2:
        raise ValueError(
            f"Found only {len(levels)} group(s) ({levels}) in column '{chosen_group_col}' - "
            "at least two are needed to compare."
        )

    meta = {
        "orientation": detected_orientation,
        "group_column": chosen_group_col,
        "group_levels": levels,
        "group_counts": groups.value_counts().to_dict(),
        "n_samples": len(values_df),
        "n_compounds": len(value_cols),
        # Columns whose name repeats an earlier one - kept as separate
        # measurements and averaged per compound after the statistics.
        "n_repeated_names": n_repeated_names,
        # Columns/rows ignored because they were not compounds: no numeric value
        # at all, or values but no name.
        "n_empty_dropped": n_empty_dropped,
        "n_unnamed_dropped": n_unnamed_dropped,
        "candidate_group_columns": candidates,
    }
    return values_df, groups, meta


def compute_compound_statistics(values_df: pd.DataFrame, groups: pd.Series,
                                group_a: str, group_b: str,
                                log_transform: bool = True,
                                min_per_group: int = 3) -> "tuple[pd.DataFrame, dict]":
    """Per-compound Welch's t-test between two groups.

    An explicit A-vs-B contrast is required even when the table has more
    than two groups - never silently collapsed into an ANOVA, because an
    unsigned F-statistic cannot honestly feed the signed, two-sided Wilcoxon
    test downstream (build_scores_from_stats / _run_wilcoxon_for_namespace).

    No imputation of missing values: Welch's test handles unequal n
    honestly, and silently filling gaps would be a much larger, unstated
    assumption for concentration data (a missing value often means "below
    detection limit", not "missing at random"). Compounds with fewer than
    min_per_group finite values in EITHER group are dropped and counted
    (n_dropped_low_n), never silently kept with a degenerate test.

    log_transform (default True) log2-transforms before differencing, with
    an offset of half the smallest positive value in the whole matrix so
    zeros/near-zeros don't blow up - the offset is reported in meta so it can
    be surfaced as a caveat. This is deliberately NOT MetaboAnalyst's full
    normalization suite - one explicit, auditable transform, nothing more.

    Returns (stats_df with columns Name/n_a/n_b/mean_a/mean_b/log2fc/
    t_stat/pvalue/padj, meta).
    """
    from scipy import stats as _scipy_stats

    mask_a = groups == group_a
    mask_b = groups == group_b
    if int(mask_a.sum()) < 2 or int(mask_b.sum()) < 2:
        raise ValueError(
            f"Each group needs at least 2 samples (got {int(mask_a.sum())} in '{group_a}', "
            f"{int(mask_b.sum())} in '{group_b}')."
        )

    work = values_df.copy()
    offset = 0.0
    if log_transform:
        # A log2 of a value <= 0 is undefined. Concentrations and peak areas are
        # never negative, so a negative value means the table was ALREADY
        # log-transformed and/or centred/auto-scaled (e.g. MetaboAnalyst's
        # normalised export). Clipping those to zero would silently corrupt every
        # such value, so refuse and say what to do instead.
        n_negative = int((work < 0).sum().sum())
        if n_negative:
            raise ValueError(
                f"{n_negative} value(s) in the table are negative, so the log2 transform "
                "cannot be applied: concentrations and peak areas are never negative, which "
                "means this table looks already log-transformed and/or scaled (for example a "
                "normalised export). Untick 'Log2-transform' and, if the columns were "
                "auto-scaled, rank by Welch's t rather than the mean difference."
            )
        positive_min = work[work > 0].min().min()
        offset = float(positive_min) / 2.0 if pd.notna(positive_min) and positive_min > 0 else 1e-9
        work = np.log2(work.clip(lower=0) + offset)

    rows = []
    n_dropped_low_n = 0
    for pos, compound in enumerate(work.columns):
        a = work.loc[mask_a].iloc[:, pos].dropna()
        b = work.loc[mask_b].iloc[:, pos].dropna()
        if len(a) < min_per_group or len(b) < min_per_group:
            n_dropped_low_n += 1
            continue
        mean_a, mean_b = float(a.mean()), float(b.mean())
        try:
            t_stat, p_val = _scipy_stats.ttest_ind(a, b, equal_var=False)
            t_stat, p_val = float(t_stat), float(p_val)
        except Exception:
            t_stat, p_val = float("nan"), 1.0
        rows.append({
            "Name": compound, "n_a": len(a), "n_b": len(b),
            "mean_a": mean_a, "mean_b": mean_b,
            "log2fc": mean_a - mean_b, "t_stat": t_stat, "pvalue": p_val,
        })

    stats_df = pd.DataFrame(rows)
    if not stats_df.empty:
        finite = stats_df["pvalue"].notna()
        stats_df["padj"] = float("nan")
        if finite.any():
            stats_df.loc[finite, "padj"] = benjamini_hochberg(stats_df.loc[finite, "pvalue"].tolist())

    meta = {
        "group_a": group_a, "group_b": group_b,
        "n_a": int(mask_a.sum()), "n_b": int(mask_b.sum()),
        "log_transform": log_transform, "log_offset": offset if log_transform else None,
        "n_compounds_tested": len(stats_df), "n_dropped_low_n": n_dropped_low_n,
    }
    return stats_df, meta


def build_scores_from_stats(matched_df: pd.DataFrame, stats_df: pd.DataFrame,
                            metric: str = QUANT_METRIC_DEFAULT) -> "tuple[dict, int]":
    """{compound_key: score} for the competitive Wilcoxon test, from a
    SEPARATELY-COMPUTED per-compound statistics table (compute_compound_
    statistics) - not from log2FC/pvalue/padj columns living on matched_df,
    which this module briefly (and wrongly) required the user to supply
    themselves. metric picks which of stats_df's two signed columns to use:
    "log2fc" (default - Perseus 1D enrichment is defined on exactly this
    mean difference) or "t_stat" (Welch's t). Both are signed, so the
    downstream alternative="two-sided" Wilcoxon call still holds.

    Joins stats_df to matched_df by Name (compute_compound_statistics' rows
    are keyed on the concentration table's own compound column headers,
    which is the same Name space the upload's identifier-resolution ran
    against). Duplicate names mean-aggregate, same convention as every other
    score-building function in this codebase.
    """
    if metric not in QUANT_METRICS:
        metric = QUANT_METRIC_DEFAULT
    name_to_key = dict(zip(matched_df["Name"].astype(str), matched_df["compound_key"].astype(str)))
    scores: dict = {}
    counts: dict = {}
    for _, row in stats_df.iterrows():
        key = name_to_key.get(str(row["Name"]))
        if not key:
            continue
        value = row.get(metric)
        if value is None or not np.isfinite(value):
            continue
        if key in scores:
            scores[key] += value
            counts[key] += 1
        else:
            scores[key] = value
            counts[key] = 1
    scores = {k: scores[k] / counts[k] for k in scores}
    n_dropped = max(len(stats_df) - len(scores), 0)
    return scores, n_dropped


def run_quantitative_enrichment(
    matched_df: pd.DataFrame,
    crosswalk: "dict | None",
    stats_df: pd.DataFrame,
    databases: "list[str]",
    taxon_id: str = "",
    taxon_name: str = "",
    metric: str = QUANT_METRIC_DEFAULT,
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    chemical_class_level: str = "sub_class",
    lipid_class_level: str = "lm_main_class",
    pathbank_primary_only: bool = True,
    background_compounds: "set | None" = None,
    kegg_organism_code: str = "",
    kegg_organism_name: str = "",
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Quantitative Enrichment: competitive two-sided Wilcoxon per term,
    scored from stats_df (see compute_compound_statistics), against the
    experiment-scoped universe (compounds actually measured in this
    concentration table - see the universe builders' section comment for
    why this tile alone keeps that narrower rule).

    This is Perseus-style 1D annotation enrichment, NOT GSEA - see the
    module docstring and Proteomics_module's own rationale: a metabolomics
    panel is far too short and too dropout-biased for a permutation-null
    ranked-list method. One family per database, two-sided because a
    downward shift is as interesting as an upward one.
    """
    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    universes = _build_universes(
        matched_df, crosswalk, databases, universe_scope="experiment",
        taxon_id=taxon_id, taxon_name=taxon_name,
        chemical_class_level=chemical_class_level, lipid_class_level=lipid_class_level,
        pathbank_primary_only=pathbank_primary_only,
        kegg_organism_code=kegg_organism_code, kegg_organism_name=kegg_organism_name,
        force_refresh=force_refresh,
        verbose=verbose, progress_callback=progress_callback,
    )

    key_to_name = _key_to_name_map(matched_df)
    background = {str(b).strip() for b in background_compounds} if background_compounds else None
    score_by_key, n_dropped = build_scores_from_stats(matched_df, stats_df, metric)

    steps_total = max(len(databases), 1)
    step = 0
    tables: dict = {}

    def _empty_block(note: str, caveat: "str | None" = None) -> dict:
        return {
            "table": pd.DataFrame(columns=METABOLOMICS_BURDEN_RESULT_COLUMNS),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0,
            "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
            "note": note, "caveat": caveat,
        }

    for database in databases:
        step += 1
        _progress(step, steps_total, database)
        universe = universes.get(database) or {}
        terms = universe.get("terms") or {}
        eligible = _eligible_keys(matched_df, database)

        annotated: set = set()
        for info in terms.values():
            annotated |= info["genes"]
        universe_keys = annotated & eligible
        if background:
            universe_keys &= background

        if not terms:
            # A universe that explains its own emptiness carries a `note`; the
            # organism-resolution caveat then travels beside it rather than
            # being swallowed.
            tables[database] = _empty_block(
                universe.get("note") or universe.get("caveat")
                or f"No {database} annotations available for this run.",
                caveat=universe.get("caveat") if universe.get("note") else None)
            continue
        if not (universe_keys & set(score_by_key)):
            tables[database] = _empty_block(
                "No scored compound carried an annotation in this database.")
            continue

        block = _run_wilcoxon_for_namespace(
            terms, score_by_key, universe_keys, True, min_term_size, max_term_size,
            alternative="two-sided",
        )
        block["n_input_dropped"] = n_dropped
        if universe.get("caveat"):
            block["caveat"] = universe["caveat"]
        tables[database] = block

    _finalize_member_column(tables, key_to_name, "Genes",
                            METABOLOMICS_BURDEN_RESULT_COLUMNS)

    summary = {
        "metric": metric if metric in QUANT_METRICS else QUANT_METRIC_DEFAULT,
        "n_scored": len(score_by_key),
        "n_dropped": n_dropped,
        "n_matched": int((matched_df["compound_key"].astype(str) != "").sum()),
        "background_applied": bool(background),
        "universe_scope": "experiment",
    }
    return {"tables": tables, "summary": summary,
            "universe_caveats": {db: (u or {}).get("caveat") for db, u in universes.items()}}


# ---------------------------------------------------------------------------
# Pathway Topology Analysis
#
# Adds a topology "impact" score to KEGG pathway enrichment: how central the
# user's matched compounds sit within that pathway's own reaction graph, not
# only how many of them are annotated to it. Method: MetPA (Xia & Wishart
# 2010, Bioinformatics 26:2342) defines impact as
#   sum(importance of matched nodes) / sum(importance of every node)
# for a chosen node-importance measure - here relative betweenness
# centrality (Freeman 1977; Brandes' algorithm, Brandes 2001, J. Math.
# Sociol. 25:163-177 - this module's own implementation, no networkx
# dependency) by default, or out-degree. The ENRICHMENT p-value itself is
# computed by calling the exact same building blocks run_ora / run_
# quantitative_enrichment use for their own KEGG database entry
# (build_kegg_compound_universe, _eligible_keys, _drop_zero_overlap,
# _run_fisher_for_namespace / _run_wilcoxon_for_namespace) - so a Pathway
# Topology run and the corresponding ORA/Quantitative KEGG run report
# IDENTICAL p-values/FDRs for the same input and config; only the extra
# Impact column is new.
# ---------------------------------------------------------------------------

def parse_kegg_kgml(xml_text: str) -> dict:
    """Parse a KGML pathway document into a compound-only reaction graph.

    Nodes are every `type="compound"` KGML <entry> - KEGG C-numbers, plus
    glycan `gl:` and drug `dr:` ids kept as their own graph nodes (only
    `cpd:C#####` ids can later be matched to the user's own compounds, but
    every node still shapes the graph's topology - dropping glycan nodes
    would silently understate a glycan-heavy map like hsa00510's true
    connectivity). An <entry> naming several ids at once (KEGG sometimes
    groups near-duplicate compound tokens on one drawn node) adds each id as
    its own node sharing that entry's coordinates - KEGG draws them all at
    the same spot, so this is the truthful placement, not a guess.

    Edges come from <reaction>: every substrate -> every product; a
    `type="reversible"` reaction adds the reverse direction too - an
    `irreversible` one is one-directional only, matching the biology.

    Returns {"nodes": {id: {"x","y","label"}}, "edges": [[src,dst], ...],
    "n_reactions": int}. A pathway with entries but zero <reaction> elements
    (real - some KEGG maps draw complexes/processes with no reaction graph,
    e.g. hsa00190) returns an empty edge list, never an error; the caller
    reports this as a "no reactions" flag rather than a fetch failure.
    """
    root = ET.fromstring(xml_text)
    nodes: dict = {}
    for entry in root.findall("entry"):
        if entry.get("type") != "compound":
            continue
        names = [n.strip() for n in str(entry.get("name") or "").split() if n.strip()]
        graphics = entry.find("graphics")
        x = float(graphics.get("x")) if graphics is not None and graphics.get("x") else None
        y = float(graphics.get("y")) if graphics is not None and graphics.get("y") else None
        label = graphics.get("name") if graphics is not None else None
        for raw_id in names:
            cid = raw_id.split(":", 1)[1] if ":" in raw_id else raw_id
            if cid not in nodes:
                nodes[cid] = {"x": x, "y": y, "label": label or cid}

    edges: list = []
    n_reactions = 0
    for reaction in root.findall("reaction"):
        n_reactions += 1
        substrates = [s.get("name", "").split(":", 1)[-1] for s in reaction.findall("substrate")]
        products = [p.get("name", "").split(":", 1)[-1] for p in reaction.findall("product")]
        reversible = reaction.get("type") == "reversible"
        for s in substrates:
            for p in products:
                if s == p:
                    continue
                edges.append([s, p])
                if reversible:
                    edges.append([p, s])
    return {"nodes": nodes, "edges": edges, "n_reactions": n_reactions}


def build_kegg_compound_name_index(verbose: bool = False, force_refresh: bool = False) -> dict:
    """{C-number: primary common name} from KEGG's own list/compound - one
    global bulk file (~19,000 compounds, ~1MB), fetched once and cached like
    the compound-pathway link file. KEGG's compound entries carry several
    semicolon-joined synonyms ("H2O; Water"); the FIRST is KEGG's own
    preferred/primary name, which is what this returns.

    Used to label Pathway Topology's graph-view nodes with a real compound
    name instead of a bare KEGG ID - a KGML <entry> for a compound rarely
    carries a human-readable graphics name the way a gene/ortholog entry
    does (verified live 2026-09-29: KEGG's own graphics name for C00022 is
    literally "C00022"). Returns {} on a fetch failure - callers fall back
    to the bare id, a safe degradation, not an error.
    """
    cache_db = f"Metabolomics_KEGG_names_v{KEGG_COMPOUND_NAMES_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, "global")
        if cached is not None and cached.get("names"):
            return cached["names"]

    status, text = _kegg_safe_get("https://rest.kegg.jp/list/compound", verbose=verbose)
    names: dict = {}
    if status == "ok" and text:
        for line in text.splitlines():
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            cid = parts[0].replace("cpd:", "").strip()
            primary = parts[1].split(";")[0].strip()
            if cid and primary:
                names[cid] = primary
    if names:
        _save_universe_cache(cache_db, "global", {"names": names, "fetched_at": _now_iso()})
    return names


def fetch_pathway_graph(org_pathway_id: str, verbose: bool = False,
                        force_refresh: bool = False) -> dict:
    """KEGG's reaction graph for one organism's pathway map (e.g. "hsa00020"),
    cached PER ORGANISM PER MAP - unlike compound membership
    (build_kegg_organism_pathway_index), a map's REACTIONS are organism-
    specific (verified live 2026-09-29: hsa00020 has 22 reactions, sce00020
    21, eco00020 21, though all three share the same 23 compounds), so this
    cannot reuse that function's "one fetch covers every organism" shortcut.

    Each node's "label" is the compound's real common name where KEGG's own
    list/compound knows one (build_kegg_compound_name_index), falling back
    to the bare KGML graphics name (usually just the C-number) otherwise.

    Isolated nodes (no edge at all) are common and REAL, not a parsing gap:
    verified live on a real pathway (hsa00983) that 33 of its 63 compound
    entries participate in ZERO <reaction> elements - KEGG's own diagram
    places compounds it does not connect with a reaction arrow in that
    specific figure (e.g. a cofactor shown beside the drawn chain, or an
    input to a linked pathway box). The frontend explains this beside the
    graph rather than hiding those nodes, which would misrepresent the
    pathway's real annotated membership.

    Returns {"status": "ok"|"nodata"|"error", "nodes", "edges", "n_reactions",
    "fetched_at"}. "nodata" (KEGG's 404 for a map with no KGML diagram) is
    cached like a real answer so it is not re-fetched every run; "error" (a
    transient failure) is never cached.
    """
    cache_db = f"Metabolomics_KEGG_kgml_v{KEGG_KGML_CACHE_VERSION}"
    if not force_refresh:
        cached = _load_universe_cache(cache_db, org_pathway_id)
        if cached is not None and cached.get("status") in ("ok", "nodata"):
            return cached

    url = f"https://rest.kegg.jp/get/{org_pathway_id}/kgml"
    status, text = _kegg_safe_get(url, verbose=verbose)
    if status == "nodata":
        result = {"status": "nodata", "nodes": {}, "edges": [], "n_reactions": 0,
                  "fetched_at": _now_iso()}
        _save_universe_cache(cache_db, org_pathway_id, result)
        return result
    if status != "ok" or not text:
        return {"status": "error", "nodes": {}, "edges": [], "n_reactions": 0,
                "fetched_at": _now_iso()}

    try:
        parsed = parse_kegg_kgml(text)
    except Exception:
        return {"status": "error", "nodes": {}, "edges": [], "n_reactions": 0,
                "fetched_at": _now_iso()}

    compound_names = build_kegg_compound_name_index(verbose=verbose, force_refresh=force_refresh)
    for cid, info in parsed["nodes"].items():
        real_name = compound_names.get(cid)
        if real_name:
            info["label"] = real_name

    result = {"status": "ok", "nodes": parsed["nodes"], "edges": parsed["edges"],
              "n_reactions": parsed["n_reactions"], "fetched_at": _now_iso()}
    _save_universe_cache(cache_db, org_pathway_id, result)
    return result


def _brandes_betweenness(nodes: list, adjacency: dict) -> dict:
    """Relative betweenness centrality via Brandes' algorithm (Brandes 2001)
    for an unweighted DIRECTED graph - this module's own implementation, no
    networkx dependency (none is installed in this environment, and KEGG
    pathway graphs are small: at most a few hundred nodes). The raw shortest-
    path-passing counts are normalised by (n-1)(n-2) so the result is
    comparable across pathways of different sizes, matching MetPA's
    convention - e.g. the middle node of an isolated 3-node chain A->B->C
    lies on the one shortest path (A,C), giving it a raw count of 1, and a
    normalised centrality of 1/((3-1)*(3-2)) = 0.5.
    """
    n = len(nodes)
    centrality = {v: 0.0 for v in nodes}
    if n < 3:
        return centrality
    for s in nodes:
        stack: list = []
        preds = {v: [] for v in nodes}
        sigma = {v: 0.0 for v in nodes}
        sigma[s] = 1.0
        dist = {v: -1 for v in nodes}
        dist[s] = 0
        queue = deque([s])
        while queue:
            v = queue.popleft()
            stack.append(v)
            for w in adjacency.get(v, ()):
                if dist[w] < 0:
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    preds[w].append(v)
        delta = {v: 0.0 for v in nodes}
        while stack:
            w = stack.pop()
            for v in preds[w]:
                if sigma[w] > 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                centrality[w] += delta[w]
    norm = (n - 1) * (n - 2)
    if norm > 0:
        for v in centrality:
            centrality[v] /= norm
    return centrality


def node_importance(graph: dict, measure: str = PATHWAY_TOPOLOGY_MEASURE_DEFAULT) -> dict:
    """{node_id: importance} for a fetch_pathway_graph()-shaped graph.

    "betweenness" (default): relative betweenness centrality, the standard
    topology-impact measure - how often a compound lies on the shortest
    reaction path between two others. "outdegree": raw out-degree - how many
    compounds this one directly feeds, MetPA's other offered measure.
    """
    nodes = list(graph.get("nodes") or {})
    if not nodes:
        return {}
    node_set = set(nodes)
    adjacency: dict = {v: [] for v in nodes}
    for edge in graph.get("edges") or []:
        s, t = edge[0], edge[1]
        if s in node_set and t in node_set:
            adjacency[s].append(t)
    if measure == "outdegree":
        return {v: float(len(adjacency[v])) for v in nodes}
    return _brandes_betweenness(nodes, adjacency)


def pathway_impact(importance: dict, matched_ids: "set | frozenset") -> "tuple[float, bool]":
    """MetPA-style pathway impact: sum(importance of matched nodes) /
    sum(importance of every node) (Xia & Wishart 2010). Returns
    (impact, well_defined) - well_defined is False when the denominator is
    zero (a map with no reactions, so every node's importance is 0 and
    "impact" has no real meaning): impact is then reported as 0.0 but
    FLAGGED by the caller, never silently indistinguishable from a
    genuinely-computed zero impact.
    """
    total = sum(importance.values())
    if total <= 0:
        return 0.0, False
    matched_total = sum(v for k, v in importance.items() if k in matched_ids)
    return matched_total / total, True


def run_pathway_topology(
    matched_df: pd.DataFrame,
    crosswalk: "dict | None",
    input_mode: str,
    stats_df: "pd.DataFrame | None" = None,
    metric: str = QUANT_METRIC_DEFAULT,
    measure: str = PATHWAY_TOPOLOGY_MEASURE_DEFAULT,
    taxon_id: str = "",
    taxon_name: str = "",
    min_term_size: int = MIN_TERM_SIZE,
    max_term_size=MAX_TERM_SIZE_DEFAULT,
    background_compounds: "set | None" = None,
    kegg_organism_code: str = "",
    kegg_organism_name: str = "",
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Pathway Topology Analysis: KEGG enrichment (Fisher for a compound
    LIST, competitive Wilcoxon for a concentration table - see input_mode)
    plus a per-pathway topology Impact score (see the section comment
    above). One result family, database "KEGG" only - Pathway Topology has
    no meaning for Reactome/SMPDB/PathBank (no graph engine built for them)
    or the two class taxonomies (disjoint terms have no internal topology).

    input_mode: "list" runs Fisher's exact test against the DATABASE-scoped
    KEGG universe (matches run_ora's KEGG path exactly). "conc" runs the
    competitive Wilcoxon test against the EXPERIMENT-scoped universe, scored
    from stats_df via build_scores_from_stats (matches run_quantitative_
    enrichment's KEGG path exactly) - stats_df is then required.

    Graphs are fetched ONLY for pathways that survive to the final table
    (the ones actually tested after the size gate and zero-overlap filter),
    never for the organism's whole map list - one independently-cached KGML
    fetch per tested pathway (fetch_pathway_graph).

    Returns {"tables": {"KEGG": block}, "summary", "universe_caveats"},
    the same shape as run_ora/run_quantitative_enrichment, so the API layer
    and frontend can serialize/render it identically. block carries the
    usual ORA/Burden keys plus "measure" and "graphs" (per Term ID: node
    positions/importance/matched flag and edges, for the graph view).
    """
    if measure not in PATHWAY_TOPOLOGY_MEASURES:
        measure = PATHWAY_TOPOLOGY_MEASURE_DEFAULT

    def _progress(step, total, label):
        if progress_callback:
            try:
                progress_callback(step, total, label)
            except Exception:
                pass

    _progress(0, 1, "Building KEGG universe...")
    universe = build_kegg_compound_universe(
        matched_df, crosswalk=crosswalk,
        universe_scope=("database" if input_mode == "list" else "experiment"),
        taxon_id=taxon_id, taxon_name=taxon_name,
        kegg_organism_code=kegg_organism_code, kegg_organism_name=kegg_organism_name,
        verbose=verbose, force_refresh=force_refresh,
    )
    key_to_name = _key_to_name_map(matched_df)
    key_to_kegg: dict = {}
    for _, row in matched_df.iterrows():
        key = _clean(row.get("compound_key"))
        kegg = _normalize_kegg_compound(row.get("kegg_id"))
        if key and kegg:
            key_to_kegg.setdefault(key, kegg)

    topology_extra_columns = ["Impact", "Matched in graph", "Graph nodes", "Topology note"]

    def _empty_block(note: str, caveat: "str | None" = None) -> dict:
        base_cols = METABOLOMICS_ORA_RESULT_COLUMNS if input_mode == "list" else METABOLOMICS_BURDEN_RESULT_COLUMNS
        return {
            "table": pd.DataFrame(columns=base_cols + topology_extra_columns),
            "n_terms_tested": 0, "n_terms_excluded_min_size": 0, "n_terms_no_overlap": 0,
            "universe_size": 0, "n_input_in_universe": 0, "n_input_dropped": 0,
            "note": note, "caveat": caveat, "measure": measure, "graphs": {},
        }

    def _empty_result(block: dict) -> dict:
        return {"tables": {"KEGG": block},
                "summary": {"input_mode": input_mode, "measure": measure, "n_matched": 0},
                "universe_caveats": {"KEGG": universe.get("caveat")}}

    terms = universe.get("terms") or {}
    if not terms:
        return _empty_result(_empty_block(
            universe.get("note") or universe.get("caveat") or "No KEGG annotations available for this run.",
            caveat=universe.get("caveat") if universe.get("note") else None))

    eligible = _eligible_keys(matched_df, "KEGG")
    background = {str(b).strip() for b in background_compounds} if background_compounds else None
    annotated: set = set()
    for info in terms.values():
        annotated |= info["genes"]

    if input_mode == "list":
        universe_keys = annotated & background if background else annotated
        list_set = eligible & universe_keys
        if not list_set:
            return _empty_result(_empty_block("None of your uploaded compounds are annotated in KEGG."))
        filtered, n_excl, n_no_overlap = _drop_zero_overlap(
            terms, list_set, universe_keys, min_term_size, max_term_size)
        block = _run_fisher_for_namespace(filtered, list_set, universe_keys, FDR_ALPHA_DEFAULT, True, 0, None)
        block["n_terms_excluded_min_size"] = n_excl
        block["n_terms_no_overlap"] = n_no_overlap
        member_col = "Overlapping Genes"
        target_columns = METABOLOMICS_ORA_RESULT_COLUMNS
        n_scored = len(list_set)
    else:
        if stats_df is None or stats_df.empty:
            raise ValueError(
                "No per-compound statistics available - configure the group contrast first.")
        universe_keys = annotated & eligible
        if background:
            universe_keys &= background
        score_by_key, n_dropped = build_scores_from_stats(matched_df, stats_df, metric)
        if not (universe_keys & set(score_by_key)):
            return _empty_result(_empty_block("No scored compound carried an annotation in KEGG."))
        block = _run_wilcoxon_for_namespace(
            terms, score_by_key, universe_keys, True, min_term_size, max_term_size,
            alternative="two-sided",
        )
        block["n_input_dropped"] = n_dropped
        member_col = "Genes"
        target_columns = METABOLOMICS_BURDEN_RESULT_COLUMNS
        n_scored = len(score_by_key)

    if universe.get("caveat"):
        block["caveat"] = universe["caveat"]

    table = block["table"]
    n_rows = len(table)
    graphs: dict = {}
    impacts, matched_counts, node_counts, notes = [], [], [], []
    for i in range(n_rows):
        row = table.iloc[i]
        term_id = str(row["Term ID"])
        _progress(i + 1, max(n_rows, 1), f"Fetching pathway graph {i + 1}/{n_rows}: {term_id}")
        graph = fetch_pathway_graph(term_id, verbose=verbose, force_refresh=force_refresh)
        overlap_keys = {k.strip() for k in str(row[member_col]).split(",") if k.strip()}
        matched_cids = {key_to_kegg[k] for k in overlap_keys if k in key_to_kegg} & set(graph["nodes"])
        importance = node_importance(graph, measure)
        impact, well_defined = pathway_impact(importance, matched_cids)
        if graph["status"] == "nodata":
            note = "KEGG has no KGML diagram for this map."
        elif graph["status"] == "error":
            note = "Could not fetch this pathway's graph from KEGG."
        elif not well_defined:
            note = "This map's KGML has no reactions, so impact is undefined (shown as 0)."
        else:
            note = ""
        impacts.append(impact)
        matched_counts.append(len(matched_cids))
        node_counts.append(len(graph["nodes"]))
        notes.append(note)
        graphs[term_id] = {
            "nodes": [
                {"id": cid, "label": info.get("label") or cid, "x": info.get("x"), "y": info.get("y"),
                 "importance": importance.get(cid, 0.0), "matched": cid in matched_cids}
                for cid, info in graph["nodes"].items()
            ],
            "edges": [list(e) for e in graph["edges"]],
            "status": graph["status"], "n_reactions": graph["n_reactions"],
        }

    table = table.assign(Impact=impacts, **{
        "Matched in graph": matched_counts, "Graph nodes": node_counts, "Topology note": notes,
    })
    block["table"] = table
    block["measure"] = measure
    block["graphs"] = graphs

    _finalize_member_column({"KEGG": block}, key_to_name, member_col,
                            target_columns + topology_extra_columns)

    summary = {
        "input_mode": input_mode, "measure": measure, "n_scored": n_scored,
        "n_matched": int((matched_df["compound_key"].astype(str) != "").sum()),
        "background_applied": bool(background),
        "universe_scope": "database" if input_mode == "list" else "experiment",
    }
    return {"tables": {"KEGG": block}, "summary": summary,
            "universe_caveats": {"KEGG": universe.get("caveat")}}
