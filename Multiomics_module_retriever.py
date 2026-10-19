"""Multi-omics retriever - thin glue between an upload and Multiomics_module.

This file holds ONLY multi-omics glue. It contains no identifier-resolution
logic of its own: every resolver it calls is imported, unedited, from the
single-omics retriever for that omics type, so an ID resolves here exactly as it
resolves in that omics' own tile (tests/test_multiomics_retriever.py checks the
routing). What it does add:

  - which existing loader / resolver handles each omics type;
  - the per-omics list rule: Genomics and Metabolomics take a plain list (the
    uploaded features ARE the test set); Transcriptomics and Proteomics take a
    differential-expression table, and Omiscera applies the padj / |log2FC|
    cutoffs with that module's OWN selection function, up- and down-regulated
    pooled (the joint question is "is this pathway perturbed", and a compound
    list has no direction);
  - the reference set for each: the genes / proteins TESTED in the experiment
    (the field's standard for expression data), never the whole database.

Metabolomics resolution (name -> crosswalk -> matched table) is driven from
API.py through the existing _resolve_and_match_compound_names helper, so it is
literally the same code path the Metabolomics tiles use.
"""

from __future__ import annotations

import pandas as pd

import Genomics_module_retriever as gene_resolver
import Proteomics_module_retriever as proteomics_retriever
import Transcriptomics_module_retriever as transcriptomics_retriever

OMICS_TYPES = ("genomics", "transcriptomics", "proteomics", "metabolomics")
DE_OMICS = ("transcriptomics", "proteomics")
PADJ_CUTOFF_DEFAULT = 0.05
LFC_CUTOFF_DEFAULT = 1.0

# Genomics statuses that mean "this gene has no usable identifier". "ambiguous"
# rows are never guessed: Multi-omics has no candidate-picking screen, so they
# are left out and counted in the resolve summary.
_GENOMICS_UNUSABLE = ("ambiguous", "not_found")


def parse_upload(omics_type: str, raw_bytes: bytes) -> dict:
    """Parse an uploaded file for Genomics, Transcriptomics or Proteomics.

    Returns {"df", "kind", "identifier_type", "column_mapping", "n_duplicates_removed"}.
    Metabolomics is parsed by API._load_metabolomics_list_upload (existing,
    shared by the Metabolomics tiles) - not duplicated here.
    """
    if omics_type == "genomics":
        df, identifier_type = gene_resolver.load_gene_table_from_bytes(raw_bytes)
        df, dup = gene_resolver.deduplicate_gene_table(df, identifier_type)
        return {"df": df, "kind": "list", "identifier_type": identifier_type, "column_mapping": None,
                "n_duplicates_removed": int(dup["duplicates_removed"])}
    if omics_type == "transcriptomics":
        df = transcriptomics_retriever.load_de_table_from_bytes(raw_bytes)
        mapping = transcriptomics_retriever.detect_columns(df)
        return {"df": df, "kind": "de", "identifier_type": None, "column_mapping": mapping,
                "n_duplicates_removed": 0}
    if omics_type == "proteomics":
        df = proteomics_retriever.load_de_table_from_bytes(raw_bytes)
        mapping = proteomics_retriever.detect_columns(df)
        return {"df": df, "kind": "de", "identifier_type": None, "column_mapping": mapping,
                "n_duplicates_removed": 0}
    raise ValueError(f"parse_upload does not handle {omics_type!r}")


def resolve_upload(omics_type: str, df: pd.DataFrame, taxon_id: str, identifier_type: "str | None" = None,
                   progress_callback=None) -> dict:
    """Resolve an uploaded table with that omics' own resolver. Returns
    {"resolved_df", "meta"} (Genomics' resolver returns no meta, so one is built)."""
    if omics_type == "genomics":
        resolved = gene_resolver.resolve_gene_table(
            df, identifier_type, taxon_id=taxon_id if identifier_type == "gene_symbol" else None,
            verbose=False)
        status = resolved["Resolution_status"].value_counts().to_dict()
        n_ambiguous = int(status.get("ambiguous", 0))
        n_unusable = n_ambiguous + int(status.get("not_found", 0))
        return {"resolved_df": resolved, "meta": {
            "identifier_type": identifier_type, "n_input_rows": len(df),
            "n_resolved": len(resolved) - n_unusable, "n_ambiguous": n_ambiguous,
            "n_not_found": int(status.get("not_found", 0)), "status_counts": status}}
    if omics_type == "transcriptomics":
        resolved, meta = transcriptomics_retriever.resolve_de_table(
            df, taxon_id, identifier_type=None, progress_callback=progress_callback, verbose=False)
        return {"resolved_df": resolved, "meta": meta}
    if omics_type == "proteomics":
        resolved, meta = proteomics_retriever.resolve_protein_table(
            df, taxon_id, progress_callback=progress_callback, verbose=False)
        return {"resolved_df": resolved, "meta": meta}
    raise ValueError(f"resolve_upload does not handle {omics_type!r}")


def assign_genomics_kegg_ids(resolved_df: pd.DataFrame, kegg_organism_code: str) -> pd.DataFrame:
    """Genomics' own KEGG-gene-ID pass (Genomics_module_retriever.resolve_kegg_ids),
    run with the ONE organism code the whole Multi-omics run shares."""
    id_col = "Gene Symbol" if "Gene Symbol" in resolved_df.columns else "Ensembl Gene ID"
    return gene_resolver.resolve_kegg_ids(resolved_df, id_col, kegg_organism_code, verbose=False)


def usable_genomics_rows(resolved_df: pd.DataFrame) -> pd.DataFrame:
    """Rows with a usable identifier: not ambiguous, not unresolved."""
    return resolved_df[~resolved_df["Resolution_status"].isin(_GENOMICS_UNUSABLE)].copy()


def select_transcriptomics(resolved_df: pd.DataFrame, padj_cutoff: float, lfc_cutoff: float) -> dict:
    """Significant genes (up AND down pooled) via Transcriptomics_module.select_degs
    - that module's own rule - plus the tested-gene reference table."""
    import Transcriptomics_module as tm

    universe_df = resolved_df[resolved_df["gene_symbol"].astype(str).str.strip() != ""]
    degs = tm.select_degs(resolved_df, padj_cutoff, lfc_cutoff, "all")
    symbols = set(degs["gene_symbol"].astype(str).str.upper())
    tested = set(universe_df["gene_symbol"].astype(str).str.upper())
    return {"universe_df": universe_df, "significant": symbols & tested,
            "n_tested": len(tested), "n_significant": len(symbols & tested)}


def select_proteomics(resolved_df: pd.DataFrame, padj_cutoff: float, lfc_cutoff: float) -> dict:
    """Significant proteins (up AND down pooled) via Proteomics_module's own
    select_significant_proteins."""
    import Proteomics_module as pm

    universe_df = resolved_df[resolved_df["accession"].astype(str).str.strip() != ""]
    sig = pm.select_significant_proteins(resolved_df, padj_cutoff, lfc_cutoff, "all")
    tested = set(universe_df["accession"].astype(str))
    significant = set(sig["accession"].astype(str)) & tested
    return {"significant": significant, "n_tested": len(tested), "n_significant": len(significant)}


# ---------------------------------------------------------------------------
# Cross-omics Pathway Mapping: fetch + disk-cache one pathway's structure.
# Parsing lives in Multiomics_mapping (pure); this is the only place that
# touches the network for it. Only complete answers ("ok") and a definite
# "nodata" are cached - a transport "error" never is, so a flaky moment cannot
# poison the cache (the same rule as the other universe builders).
# ---------------------------------------------------------------------------

import json as _json
import time as _time

import Genomics_module as _gm
import Multiomics_mapping as _mapping

KEGG_KGML_URL = "https://rest.kegg.jp/get/{pid}/kgml"
REACTOME_CONTENT = "https://reactome.org/ContentService/data"
REACTOME_DIAGRAM_URL = "https://reactome.org/download/current/diagram/{stid}.json"
REACTOME_GRAPH_URL = "https://reactome.org/download/current/diagram/{stid}.graph.json"
MAP_KEGG_CACHE_DB = "MultiomicsMap_KEGG_kgml_v4"  # v4: line colours + pathway numbers; v3: line graphics, groups, relations; v2: no TITLE box
MAP_REACTOME_CACHE_DB = "MultiomicsMap_Reactome_diagram_v2"  # v2: link boxes keep their sub-pathway, overview status
_HTTP_RETRIES = 2


def _http_get_json(url: str):
    """GET a JSON document. Returns (status, obj): "ok", "nodata" (HTTP 404 - not
    retried) or "error" (obj is the message)."""
    if _gm.HTTP_SESSION is None:
        return "error", "requests not installed"
    headers = dict(_gm.OMISCERA_HEADERS)
    headers["Accept"] = "application/json"
    last = ""
    for attempt in range(_HTTP_RETRIES + 1):
        try:
            resp = _gm.HTTP_SESSION.get(url, headers=headers, timeout=_gm.DEFAULT_TIMEOUT)
            if resp.status_code == 200:
                return "ok", resp.json()
            if resp.status_code in (400, 404):
                return "nodata", None
            last = f"HTTP {resp.status_code}"
            if resp.status_code in (429, 500, 502, 503, 504) and attempt < _HTTP_RETRIES:
                _time.sleep(2.0)
                continue
            break
        except Exception as exc:  # transport or JSON error
            last = str(exc)
            if attempt < _HTTP_RETRIES:
                _time.sleep(1.0)
    return "error", last


def _cached_structure(cache_db: str, key: str, force_refresh: bool):
    if force_refresh:
        return None
    return _gm._load_universe_cache(cache_db, key)


def fetch_kegg_pathway_structure(pathway_id: str, force_refresh: bool = False, verbose: bool = False) -> dict:
    """One KEGG pathway's genes, compounds and reactions (KGML), parsed.
    Returns {"status": "ok"|"nodata"|"error", "structure": dict|None, "note": str|None}."""
    pid = str(pathway_id or "").strip()
    cached = _cached_structure(MAP_KEGG_CACHE_DB, pid, force_refresh)
    if cached is not None:
        return cached
    status, text = _gm._kegg_safe_get(KEGG_KGML_URL.format(pid=pid), verbose=verbose)
    if status == "error":
        return {"status": "error", "structure": None,
                "note": "Could not reach KEGG to fetch this pathway's map. Try again."}
    if status == "nodata" or not text.strip():
        out = {"status": "nodata", "structure": None,
               "note": f"KEGG has no pathway map (KGML) for {pid}."}
    else:
        try:
            out = {"status": "ok", "structure": _mapping.parse_kgml_full(text), "note": None}
        except Exception as exc:
            return {"status": "error", "structure": None, "note": f"KEGG's map for {pid} could not be read: {exc}"}
    _gm._save_universe_cache(MAP_KEGG_CACHE_DB, pid, out)
    return out


def _reactome_reaction_db_ids(contained) -> set:
    """dbIds of the reactions inside a pathway's containedEvents answer (Reactome
    repeats an already-sent object as a bare integer, so both forms count)."""
    ids: set = set()
    for e in contained or []:
        if isinstance(e, int):
            ids.add(e)
        elif isinstance(e, dict) and e.get("dbId") is not None:
            ids.add(int(e["dbId"]))
    return ids


def resolve_overview_subpathways(structure: dict) -> "list | None":
    """The sub-pathways behind an overview pathway's link boxes: [{"stId", "name"}], de-duplicated, in diagram order.
    containedEvents is no help (it repeats an already-sent object as a bare integer, losing its id), so each link box's
    database id is looked up on its own. None when Reactome cannot be reached - the caller must not cache that."""
    subs, seen = [], set()
    for node in structure.get("nodes") or []:
        db_id = node.get("reactome_id") if node.get("kind") == "link" else None
        if db_id is None:
            continue
        st, obj = _http_get_json(f"{REACTOME_CONTENT}/query/{db_id}")
        if st == "error":
            return None
        if st == "ok" and isinstance(obj, dict) and obj.get("stId") and obj["stId"] not in seen:
            seen.add(obj["stId"])
            subs.append({"stId": obj["stId"], "name": str(obj.get("displayName") or node.get("label") or obj["stId"])})
    return subs


def fetch_reactome_pathway_structure(stable_id: str, force_refresh: bool = False, verbose: bool = False) -> dict:
    """One Reactome pathway's diagram + entity graph, parsed. A pathway without a
    diagram of its own (hasDiagram = false) is drawn from its NEAREST ancestor that
    has one, restricted to this pathway's own reactions. Same return shape as the
    KEGG fetcher; "restricted_to" names the ancestor when that happened."""
    stid = str(stable_id or "").strip()
    cached = _cached_structure(MAP_REACTOME_CACHE_DB, stid, force_refresh)
    if cached is not None:
        return cached

    def _bail(status, note):
        return {"status": status, "structure": None, "note": note}

    st, info = _http_get_json(f"{REACTOME_CONTENT}/query/{stid}")
    if st == "error":
        return _bail("error", "Could not reach Reactome to fetch this pathway's diagram. Try again.")
    if st == "nodata" or not isinstance(info, dict):
        out = _bail("nodata", f"Reactome has no pathway {stid}.")
        _gm._save_universe_cache(MAP_REACTOME_CACHE_DB, stid, out)
        return out
    diagram_id, keep, restricted_to = stid, None, None
    if not info.get("hasDiagram"):
        st, anc = _http_get_json(f"{REACTOME_CONTENT}/event/{stid}/ancestors")
        if st == "error":
            return _bail("error", "Could not reach Reactome to find this pathway's diagram. Try again.")
        chain = (anc or [[]])[0] if st == "ok" and anc else []
        parent = next((a for a in chain[1:] if a.get("hasDiagram")), None)
        if parent is None:
            out = _bail("nodata", f"Reactome has no diagram for {stid} or for any pathway above it.")
            _gm._save_universe_cache(MAP_REACTOME_CACHE_DB, stid, out)
            return out
        diagram_id = parent["stId"]
        restricted_to = f"{parent.get('displayName', diagram_id)} ({diagram_id})"
        st, contained = _http_get_json(f"{REACTOME_CONTENT}/pathway/{stid}/containedEvents")
        if st == "error":
            return _bail("error", "Could not reach Reactome to list this pathway's reactions. Try again.")
        keep = _reactome_reaction_db_ids(contained) if st == "ok" else set()
    st_d, diagram = _http_get_json(REACTOME_DIAGRAM_URL.format(stid=diagram_id))
    st_g, graph = _http_get_json(REACTOME_GRAPH_URL.format(stid=diagram_id))
    if "error" in (st_d, st_g):
        return _bail("error", "Could not reach Reactome to download the diagram. Try again.")
    if st_d != "ok" or st_g != "ok":
        out = _bail("nodata", f"Reactome has no downloadable diagram for {diagram_id}.")
        _gm._save_universe_cache(MAP_REACTOME_CACHE_DB, stid, out)
        return out
    try:
        structure = _mapping.parse_reactome_diagram(diagram, graph, keep_reaction_db_ids=keep)
    except Exception as exc:
        return _bail("error", f"Reactome's diagram for {diagram_id} could not be read: {exc}")
    structure["title"] = str(info.get("displayName") or structure["title"])
    if _mapping.is_overview(structure):
        subs = resolve_overview_subpathways(structure)
        if subs is None:
            return _bail("error", "Could not reach Reactome to list this overview pathway's sub-pathways. Try again.")
        out = {"status": "overview", "structure": structure, "subpathways": subs, "restricted_to": restricted_to,
               "note": (f"Reactome draws \"{structure['title']}\" only as links to its {len(subs)} sub-pathways - "
                        "there is no gene or compound on it to put a hit on. Open one of them below.")}
        _gm._save_universe_cache(MAP_REACTOME_CACHE_DB, stid, out)
        return out
    out = {"status": "ok", "structure": structure, "note": None, "restricted_to": restricted_to}
    _gm._save_universe_cache(MAP_REACTOME_CACHE_DB, stid, out)
    return out


def fetch_pathway_structure(database: str, pathway_id: str, force_refresh: bool = False) -> dict:
    if database == "KEGG":
        return fetch_kegg_pathway_structure(pathway_id, force_refresh=force_refresh)
    if database == "Reactome":
        return fetch_reactome_pathway_structure(pathway_id, force_refresh=force_refresh)
    raise ValueError(f"Unsupported database {database!r}")


KEGG_BRITE_CLASSES_URL = "https://rest.kegg.jp/get/br:br08901"
MAP_CLASSES_CACHE_DB = "MultiomicsMap_KEGG_classes_v2"  # v2: metabolism branch only


def fetch_pathway_classes(force_refresh: bool = False, verbose: bool = False) -> dict:
    """{map number: class} from KEGG's pathway classification (br08901), disk-cached. {} when KEGG cannot be reached (the
    map is then drawn with every colour named 'Other'); only a real answer is cached."""
    cached = None if force_refresh else _gm._load_universe_cache(MAP_CLASSES_CACHE_DB, "global")
    if cached is not None and cached.get("classes"):
        return cached["classes"]
    status, text = _gm._kegg_safe_get(KEGG_BRITE_CLASSES_URL, verbose=verbose)
    classes = _mapping.parse_brite_classes(text) if status == "ok" else {}
    if classes:
        _gm._save_universe_cache(MAP_CLASSES_CACHE_DB, "global", {"classes": classes})
    return classes


def fetch_chebi_skeletons(chebi_ids) -> dict:
    """{ChEBI number: InChIKey skeleton} for a pathway's chemicals, from Metabolomics' own cached ChEBI
    lookup (it fetches only what is missing). {} when ChEBI cannot be reached - the mapping then falls back
    to exact ChEBI IDs and says nothing more than it can show."""
    import Metabolomics_module as mm

    ids = {str(c) for c in chebi_ids if str(c)}
    if not ids:
        return {}
    try:
        return mm._fetch_chebi_skeletons(ids, verbose=False)
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Correlation-based Integration: gene symbols for the features of a matrix
#
# The two resolvers' public entry points want a DE table (identifier + log2FC + p-value) and collapse duplicates by that
# statistic, which a matrix has none of. What the correlation needs is only "feature id -> gene symbol" so a transcript
# and a protein of the same gene can be paired, so the SAME lower-level mapping steps those entry points call are used
# directly (mygene.info for transcripts, UniProt ID Mapping for proteins) and every input id keeps its own answer.
# ---------------------------------------------------------------------------

def resolve_matrix_features(omics_type: str, feature_ids: list, taxon_id: str, progress_callback=None) -> dict:
    """{feature id: GENE SYMBOL} for a transcriptomics or proteomics matrix's features (a feature that did not resolve is
    absent). Returns {"gene_of", "n_features", "n_mapped", "n_unmapped", "identifier_type"}."""
    ids = [str(i) for i in feature_ids]
    if omics_type == "transcriptomics":
        id_type = transcriptomics_retriever.detect_identifier_type(pd.Series(ids))
        mapped, _ = transcriptomics_retriever.map_identifiers(
            pd.DataFrame({"identifier": ids}), id_type=id_type, taxon_id=str(taxon_id),
            progress_callback=progress_callback, verbose=False)
        gene_of = {row["identifier"]: str(row["gene_symbol"]).upper() for _, row in mapped.iterrows()
                   if row["status"] != "unmapped" and str(row["gene_symbol"]).strip()}
        # the Entrez Gene ID too: Reactome / KEGG pathways of a transcript are keyed by it (MOFA's pathway enrichment)
        entrez_of = {row["identifier"]: str(row.get("entrez_id", "") or "").strip() for _, row in mapped.iterrows()
                     if row["identifier"] in gene_of and str(row.get("entrez_id", "") or "").strip() not in ("", "nan", "None")}
        key_extra = {"entrez_of": entrez_of}
    elif omics_type == "proteomics":
        id_type = "uniprot_accession"
        query = {i: proteomics_retriever.strip_isoform(proteomics_retriever.parse_protein_group(i)).upper() for i in ids}
        unique = sorted({q for q in query.values() if q})
        resolved, n_fail = proteomics_retriever._idmapping_resolve(unique, verbose=False)
        if progress_callback:
            try:
                progress_callback(len(unique), len(unique), "resolved via UniProt ID Mapping")
            except Exception:
                pass
        if unique and n_fail and not resolved:
            raise RuntimeError("UniProt ID Mapping failed for all accessions after retries - this looks like a "
                               "network problem, not proteins that do not exist. Please try again.")
        gene_of = {i: str(resolved[q]["gene_symbol"]).upper() for i, q in query.items()
                   if q and resolved.get(q) and str(resolved[q].get("gene_symbol") or "").strip()}
        key_extra = {"acc_of": {i: q for i, q in query.items() if q}}          # the cleaned UniProt accession of every feature
    else:
        raise ValueError(f"resolve_matrix_features does not handle {omics_type!r}")
    return {"gene_of": gene_of, "n_features": len(ids), "n_mapped": len(gene_of), "n_unmapped": len(ids) - len(gene_of),
            "identifier_type": id_type, **key_extra}


# ---------------------------------------------------------------------------
# Network-based Integration: which proteins sit in the reactions of a metabolite
#
# Reactome publishes two bulk files for it: ChEBI2ReactomeReactions.txt (ChEBI -> reaction) and
# UniProt2ReactomeReactions.txt (UniProt -> reaction). Both are streamed ONCE, filtered to the organism, and only the
# (small) organism slice is cached. A protein listed for a reaction takes part in it (substrate, complex member,
# regulator or enzyme) - Reactome does not mark the catalyst in these files - which the README says.
# ---------------------------------------------------------------------------

NET_REACTOME_CACHE_DB = "MultiomicsNet_ReactomeReactions_v1"
REACTOME_CHEBI_REACTIONS_URL = "https://reactome.org/download/current/ChEBI2ReactomeReactions.txt"
REACTOME_UNIPROT_REACTIONS_URL = "https://reactome.org/download/current/UniProt2ReactomeReactions.txt"


def _stream_reactome_pairs(url: str, species_lower: str):
    """Yield (identifier, reaction id, reaction name) for the organism's rows of a Reactome 6-column file."""
    resp = _gm.HTTP_SESSION.get(url, headers=_gm.OMISCERA_HEADERS, timeout=180, stream=True)
    resp.raise_for_status()
    for raw in resp.iter_lines():
        if not raw:
            continue
        parts = raw.decode("utf-8", errors="replace").split("\t")
        if len(parts) < 6 or parts[5].strip().lower() != species_lower:
            continue
        ident, rxn = parts[0].strip(), parts[1].strip()
        if ident and rxn:
            yield ident, rxn, parts[3].strip()


def _read_reactome_rows(url: str, species_lower: str, attempts: int = 3) -> list:
    """The organism's rows of one Reactome file, read in full. Reactome's server sometimes answers 5xx or closes the connection
    half-way through a large file (seen: HTTP 522, 'Response ended prematurely'); the whole read is tried again, a few times."""
    last = None
    for attempt in range(attempts):
        try:
            return list(_stream_reactome_pairs(url, species_lower))
        except Exception as exc:
            last = exc
            if attempt + 1 < attempts:
                _time.sleep(2.0 * (attempt + 1))
    raise last


def fetch_reactome_reaction_maps(taxon_id: str, taxon_name: str, force_refresh: bool = False) -> dict:
    """{"status": "ok"|"nodata"|"error", "chebi_reactions": {ChEBI number: [reaction ids]},
    "reaction_proteins": {reaction id: [UniProt accessions]}, "reaction_names": {reaction id: name}, "note": str|None}
    for ONE organism (Reactome's species name = the organism name, as the Reactome pathway universe uses)."""
    name = str(taxon_name or "").strip()
    empty = {"chebi_reactions": {}, "reaction_proteins": {}, "reaction_names": {}}
    if not name:
        return {"status": "nodata", **empty,
                "note": "Reactome's reaction files need an organism name; none was available for this run."}
    cached = _cached_structure(NET_REACTOME_CACHE_DB, str(taxon_id), force_refresh)
    if cached is not None:
        return cached
    if _gm.HTTP_SESSION is None:
        return {"status": "error", **empty, "note": "requests is not installed."}
    species = name.lower()
    chebi_rx: dict = {}
    rx_prot: dict = {}
    names: dict = {}
    try:
        for chebi, rxn, rname in _read_reactome_rows(REACTOME_CHEBI_REACTIONS_URL, species):
            chebi_rx.setdefault(_mapping.norm_chebi(chebi), set()).add(rxn)
            names[rxn] = rname
        for acc, rxn, rname in _read_reactome_rows(REACTOME_UNIPROT_REACTIONS_URL, species):
            rx_prot.setdefault(rxn, set()).add(acc.upper())
            names.setdefault(rxn, rname)
    except Exception as exc:
        return {"status": "error", **empty, "note": f"Could not download Reactome's reaction files ({exc}). Try again."}
    if not chebi_rx or not rx_prot:
        return {"status": "nodata", **empty, "note": f"Reactome's reaction files have no rows for '{name}'."}
    out = {"status": "ok", "note": None,
           "chebi_reactions": {k: sorted(v) for k, v in chebi_rx.items() if k},
           "reaction_proteins": {k: sorted(v) for k, v in rx_prot.items()},
           "reaction_names": names}
    _gm._save_universe_cache(NET_REACTOME_CACHE_DB, str(taxon_id), out)
    return out


NET_SKELETON_CACHE_DB = "MultiomicsNet_ChEBISkeletons_v1"
_SKELETON_NEGATIVE_MIN_COVERAGE = 0.8


def fetch_chebi_skeletons_cached(chebi_ids) -> dict:
    """{ChEBI number: InChIKey skeleton} for the given ids, remembering the ids ChEBI has NO structure for.

    Metabolomics' own lookup (fetch_chebi_skeletons above) treats an id as cached only when it has a structure, so every run
    asks ChEBI again for the hundreds of ids that have none (a ~3,200-id Reactome file took minutes each time). This wrapper
    keeps its own cache that also stores those as '' and only asks for ids it has never seen. A 'no structure' answer is
    stored only when the lookup clearly worked (at least 80% of the asked ids answered); otherwise a failed batch would be
    remembered as 'no structure' for good. Returns only the ids that have a skeleton."""
    ids = {str(c) for c in chebi_ids if str(c)}
    if not ids:
        return {}
    cached = _gm._load_universe_cache(NET_SKELETON_CACHE_DB, "global") or {}
    known: dict = dict(cached.get("skeletons") or {})
    missing = ids - set(known)
    if missing:
        got = fetch_chebi_skeletons(missing)
        for c, sk in got.items():
            known[str(c)] = sk
        if got and len(got) / len(missing) >= _SKELETON_NEGATIVE_MIN_COVERAGE:
            for c in missing - set(got):
                known[c] = ""
        if got:
            _gm._save_universe_cache(NET_SKELETON_CACHE_DB, "global", {"skeletons": known})
    return {c: known[c] for c in ids if known.get(c)}
