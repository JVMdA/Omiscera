# -*- coding: utf-8 -*-
"""
Proteomics statistical analysis module for Omiscera — PPI Network / Module
Detection.

WHAT IT ANSWERS that ORA/Quantitative Enrichment cannot. Those two ask "which
pre-defined gene sets are over-represented among my significant proteins".
This asks two different questions: which proteins sit at the CENTRE of the
response according to the STRING protein-protein interaction network -
including proteins the user never measured ("recruited") - and how do the
hits group into connected functional modules.

THE MATHS IS NOT HERE. This module supplies only the pieces that depend on
Proteomics' own data shape (identifier bridge, seed selection, module
labelling); the whole pipeline (graph build, random walk with restart, the
degree-preserving permutation null, MCL clustering, the honesty caveats) is
Transcriptomics_network.run_network_analysis, called through its hook
parameters (mapper/deg_selector/module_enricher/entity_noun) rather than
duplicated. See that function's docstring for the hook contract. One copy of
the numerics means the Genomics, Transcriptomics and Proteomics tiles can
never drift into giving different answers to the same underlying question.

WHY THE IDENTIFIER BRIDGE IS ITS OWN FUNCTION, NOT graph["uniprot_to_protein"].
Genomics_network.build_string_network already parses STRING's UniProt_AC
aliases into graph["uniprot_to_protein"] - but via a silent
dict(zip(...)) last-wins, which the SAME module's own docstrings flag
elsewhere as sending a seed to an arbitrary protein on a collision. Measured
live on human STRING v12: 230 accessions map to more than one STRING protein,
and 17,224 proteins carry more than one accession (up to 267, mostly
canonical/isoform fan-out). build_accession_string_index below re-parses the
same aliases file independently, through Transcriptomics_network's own
_alias_index_add - the identical explicit collision rule (highest weighted
degree wins, then lexicographically smallest STRING id) that function already
guarantees for gene symbol/Ensembl/Entrez - so accessions get the same
guarantee, not a weaker one.

WHY MODULE LABELLING TESTS IN STRING-ID SPACE, NOT ACCESSION SPACE. GO,
Reactome and Complex Portal are all UniProt-accession-keyed (via
Genomics_module.build_go_universe / build_reactome_universe /
Proteomics_module.build_complex_portal_universe). Testing a module's
enrichment directly in accession space would double-count every protein that
carries more than one accession - the SAME fan-out problem that made Genomics
move its own KEGG universe off accession space (see Genomics_module.
build_kegg_universe's docstring). Instead, each term's accession set is
translated through the SAME accession -> STRING-id index the identifier
bridge uses, so a set of STRING ids counts each protein exactly once and no
"representative accession" ever has to be chosen for a protein. The module's
own members (mod["string_ids"], from MCL) are ALREADY in this space, so
nothing needs translating on that side at all.

WHY REACTOME WORKS HERE, UNLIKE IN TRANSCRIPTOMICS. Transcriptomics_network
deliberately does not offer Reactome for module labelling: its identifier
bridge is built from the UPLOADED table, so a recruited gene (absent from
that table by definition) would be dropped from its own module's label -
exactly the failure this whole feature exists to avoid. Proteomics' bridge is
anchored on STRING's own accessions instead (see above), so a recruited
protein's accession(s) are available regardless of what the user uploaded,
and Reactome is offered alongside GO and Complex Portal.

Also inherited, unchanged, from Transcriptomics_network/Genomics_network: RWR
seed mass must stay strictly positive (direction is expressed by WHICH
proteins seed the walk, never by a signed weight); p-values are a normal
approximation of a right-skewed permutation null, skew-corrected where the
null fits reliably; module enrichment is POST-SELECTION labelling, not a
hypothesis test; and STRING's combined_score folds in text-mining and
cross-organism transfer, so it is not pure experimental evidence and skews
toward well-studied proteins.
"""

import numpy as np
import pandas as pd

import Genomics_module as genomics_ora
import Genomics_network as gnet
import Proteomics_module as pm
import Transcriptomics_module as tx
import Transcriptomics_network as txnet

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reactome works here (see module docstring for why, vs. Transcriptomics)."
# MSigDB is NOT offered: it is gene-symbol-native, and Complex Portal already
# covers the "known complex" question PPI modules are most naturally
# compared against.
PPI_ENRICHMENT_DATABASES = ("GO", "Reactome", "Complex Portal")
PPI_ENRICHMENT_DATABASE_DEFAULT = "GO"

MAPPING_COLUMNS = ["input_id", "accession", "STRING ID", "Mapped via", "Seed", "Seed score"]


# ---------------------------------------------------------------------------
# Identifier bridge: UniProt accession -> STRING protein id
# ---------------------------------------------------------------------------

def build_accession_string_index(graph: dict, aliases_path: "str | None" = None, verbose: bool = False) -> dict:
    """Build a deterministic UniProt-accession -> STRING-id bridge for ONE
    graph, plus a gene-symbol fallback reused unchanged from
    Transcriptomics_network.build_string_id_index (same graph, same
    collision rule - see this module's docstring for why the accession side
    needs its OWN pass rather than graph["uniprot_to_protein"]).

    NEVER cache apart from its graph, for the same reason
    build_string_id_index never does: membership is confidence-dependent.

    Returns {"by_accession": {ACCESSION: string_id}, "by_symbol": {...},
    "collisions": {...} (accession-tier collisions layered onto the symbol
    index's own), "n_graph_proteins": int, "n_accessions": int}.
    """
    protein_index = graph["protein_index"]
    raw_adj = graph["raw_adjacency"]
    degrees = np.asarray(raw_adj.sum(axis=0)).ravel()

    symbol_index = txnet.build_string_id_index(graph, aliases_path=aliases_path, verbose=verbose)
    collisions: dict = dict(symbol_index["collisions"])

    by_accession: dict = {}
    if aliases_path:
        try:
            aliases_df = pd.read_csv(aliases_path, sep="\t", dtype=str)
            aliases_df.columns = [c.lstrip("#").strip() for c in aliases_df.columns]
            uniprot_rows = aliases_df[aliases_df["source"] == "UniProt_AC"]
            for pid, acc in zip(uniprot_rows["string_protein_id"], uniprot_rows["alias"]):
                if not isinstance(acc, str):
                    continue
                txnet._alias_index_add(by_accession, acc, pid, protein_index, degrees, collisions)
        except Exception as exc:  # a missing/odd aliases file must not kill the run
            if verbose:
                print(f"[ProteomicsNetwork] accession alias parsing skipped: {exc}")

    return {
        "by_accession": by_accession,
        "by_symbol": symbol_index["by_symbol"],
        "collisions": collisions,
        "n_graph_proteins": len(protein_index),
        "n_accessions": len(by_accession),
    }


def map_protein_table_to_string(resolved_df: pd.DataFrame, graph: dict, index: "dict | None" = None,
                                aliases_path: "str | None" = None, taxon_id: str = "",
                                verbose: bool = False) -> dict:
    """Map every row of the resolved protein table onto a STRING protein.

    Tiered, first hit wins: accession -> gene_symbol. Unlike Transcriptomics
    (symbol -> Ensembl -> Entrez, three independent identifier spaces
    resolved from scratch), this is a two-tier fallback within data ONE
    resolution step already produced - Proteomics_module_retriever's UniProt
    ID Mapping call returns accession and gene_symbol together for every row.

    Same coverage gates as Transcriptomics_network.map_de_table_to_string
    (imported from that module, not re-declared, so the two thresholds can
    never drift apart), worded for proteins.
    """
    if index is None:
        index = build_accession_string_index(graph, aliases_path=aliases_path, verbose=verbose)

    by_accession = index["by_accession"]
    by_symbol = index["by_symbol"]

    rows = []
    n_by_tier = {"accession": 0, "gene_symbol": 0}
    string_by_row: dict = {}
    for _, r in resolved_df.iterrows():
        acc = str(r.get("accession") or "").strip()
        sym = str(r.get("gene_symbol") or "").strip()
        pid, via = None, ""
        if acc and acc.upper() in by_accession:
            pid, via = by_accession[acc.upper()], "accession"
        elif sym and sym.upper() in by_symbol:
            pid, via = by_symbol[sym.upper()], "gene_symbol"
        if pid is not None:
            n_by_tier[via] += 1
            string_by_row[str(r.get("input_id"))] = pid
        rows.append({
            "input_id": r.get("input_id"),
            "accession": acc,
            "STRING ID": pid or "",
            "Mapped via": via,
            "Seed": False,
            "Seed score": None,
        })

    table = pd.DataFrame(rows, columns=MAPPING_COLUMNS)
    n_input = len(table)
    n_mapped = int((table["STRING ID"] != "").sum())
    fraction = (n_mapped / n_input) if n_input else 0.0

    if fraction < txnet.MIN_MAPPED_FRACTION_ERROR:
        raise ValueError(
            f"Only {n_mapped} of {n_input} proteins ({fraction:.1%}) could be matched to a protein in "
            f"STRING's network for taxon {taxon_id or '?'}. STRING names proteins using each organism's "
            f"own nomenclature, which may not match the identifiers in your table. Network analysis needs "
            f"at least {txnet.MIN_MAPPED_FRACTION_ERROR:.0%} coverage to be meaningful."
        )

    return {
        "table": table,
        "string_by_row": string_by_row,
        "n_input": n_input,
        "n_mapped": n_mapped,
        "n_unmapped": n_input - n_mapped,
        "n_by_tier": n_by_tier,
        "fraction_mapped": fraction,
        "n_collisions_used": len(index["collisions"]),
        "low_mapping_warning": (
            f"Only {fraction:.1%} of your proteins matched a STRING protein. Results are based on that "
            f"subset and may be unrepresentative."
        ) if fraction < txnet.LOW_MAPPED_FRACTION_WARN else None,
    }


# ---------------------------------------------------------------------------
# Per-module functional enrichment, in STRING-protein-ID space
# ---------------------------------------------------------------------------

def enrich_modules_by_accession(modules: list, graph: dict, taxon_id: str,
                                taxon_name: str = "",
                                databases: "list[str] | None" = None,
                                include_iea: bool = True,
                                msigdb_collections: "list[str] | None" = None,
                                min_term_size: int = genomics_ora.MIN_TERM_SIZE,
                                max_term_size=gnet.DEFAULT_MODULE_MAX_TERM_SIZE,
                                verbose: bool = False) -> dict:
    """Label each module with its most enriched term(s), in STRING-protein-ID
    space - see this module's docstring for why that (not accession space)
    is the correct universe: a graph protein can carry more than one UniProt
    accession, and GO/Reactome/Complex Portal's own gene sets are
    accession-keyed.

    databases: any non-empty subset of PPI_ENRICHMENT_DATABASES; defaults to
    [PPI_ENRICHMENT_DATABASE_DEFAULT] when omitted. Every requested
    database's namespaces (GO's 3) or single term set (Reactome, Complex
    Portal) contribute into ONE shared "how many independent tests were
    compared" pool - selecting several databases at once is the SAME
    best-of-N-namespaces correction this function already applied within
    GO alone, now spanning more namespaces, not a new statistical model.

    Each database's term->accessions map is translated through
    build_accession_string_index's accession -> STRING-id index into
    term->STRING-ids (a set, so a protein with several accessions is counted
    once); terms that translate to nothing are dropped, the same "identifier
    gap, not a biology conclusion" discipline used elsewhere in this app.
    The module's own genes (mod["string_ids"]) need no translation at all -
    they are already in this space.

    THE UNIVERSE IS THE WHOLE GRAPH'S proteins (graph["n_proteins"]), not a
    per-database annotated subset and not the DE table's - module members
    are drawn from the graph, so the graph is the sampling frame (same rule
    Transcriptomics_network.enrich_modules_by_symbol documents, for the same
    reason: a recruited protein must not be dropped from its own module's
    label).

    taxon_name is REQUIRED for "Reactome" (build_reactome_universe filters
    by organism NAME, unlike GO which filters by taxon_id) - Reactome is
    offered here specifically because this bridge does not depend on the
    uploaded table (see module docstring for why Transcriptomics excludes it).

    msigdb_collections: accepted and ignored - MSigDB is not in
    PPI_ENRICHMENT_DATABASES, but run_network_analysis calls whichever
    module_enricher it was given through one uniform keyword set.

    Emits Genomics_module.ORA_RESULT_COLUMNS unchanged (the shared Fisher/
    BH-FDR engine every other accession-space Proteomics database uses). The
    "Overlapping Genes" member column is translated to gene SYMBOL for
    readability once each namespace/database's test is done - purely
    cosmetic, the statistics never depended on picking a single
    representative identifier for a protein.

    Returns {"results": {module_id: {"best": {"namespace", "block",
    "top_fdr_corrected"}, "candidates": {namespace_or_database_key: {
    "block"} | None for every tested key}, "n_namespaces_compared": int} |
    None}, "universe_size", "databases", "caveat", "selection_caveat"}.
    "candidates" keeps EVERY tested namespace/database's own result, not
    just the winner ("best") - so a caller can show a module's GO result and
    its Reactome result side by side, not only whichever one happened to win.
    """
    protein_gene = graph["protein_gene"]
    aliases_path = gnet._string_file_paths(taxon_id).get("aliases")
    accession_index = build_accession_string_index(graph, aliases_path=aliases_path, verbose=verbose)
    by_accession = accession_index["by_accession"]

    def _to_string_id_terms(raw_terms: dict) -> dict:
        out: dict = {}
        for term_id, info in raw_terms.items():
            ids = {by_accession[a.upper()] for a in info["genes"] if a.upper() in by_accession}
            if ids:
                out[term_id] = {"name": info["name"], "genes": ids}
        return out

    databases = list(dict.fromkeys(databases)) if databases else [PPI_ENRICHMENT_DATABASE_DEFAULT]
    unsupported = [d for d in databases if d not in PPI_ENRICHMENT_DATABASES]
    if unsupported:
        raise ValueError(f"Unsupported module-enrichment database(s): {unsupported!r}")

    term_maps: dict = {}
    caveats: list = []
    if "GO" in databases:
        universe = genomics_ora.build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose)
        found = False
        for ns in genomics_ora.GO_NAMESPACES:
            translated = _to_string_id_terms(universe.get(ns, {}) or {})
            if translated:
                term_maps[ns] = translated
                found = True
        if not found:
            caveats.append("GO annotations could not be fetched for this organism, so GO-based labels are unavailable.")
    if "Reactome" in databases:
        if not taxon_name.strip():
            caveats.append(
                "Reactome requires a resolved organism name to filter pathways by species; none was "
                "available for this run, so Reactome-based labels are unavailable."
            )
        else:
            universe = genomics_ora.build_reactome_universe(
                taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose)
            translated = _to_string_id_terms(universe.get("pathways", {}) or {})
            if translated:
                term_maps["Reactome"] = translated
            else:
                caveats.append(f"No Reactome pathway data is available for \"{taxon_name}\", so "
                               f"Reactome-based labels are unavailable.")
    if "Complex Portal" in databases:
        universe = pm.build_complex_portal_universe(taxon_id, verbose=verbose)
        translated = _to_string_id_terms(universe.get("terms", {}) or {})
        if translated:
            term_maps["Complex Portal"] = translated
        else:
            caveats.append(universe.get("caveat") or
                           f"Complex Portal has no data for taxon {taxon_id}, so Complex-Portal-based "
                           f"labels are unavailable.")
    caveat = " ".join(caveats) if caveats else None

    # Precompute each namespace's STRING-id universe ONCE - unioning all
    # term members is expensive and identical for every module.
    universe_by_ns = {}
    for key, term_map in term_maps.items():
        u: set = set()
        for info in term_map.values():
            u |= info["genes"]
        universe_by_ns[key] = u

    def _translate_overlap(block: dict) -> dict:
        # Translate a block's member column (STRING ids) to gene symbols for
        # display - never re-tested, purely cosmetic. Applied to every
        # tested candidate, not just the eventual winner, since any of them
        # may be shown on request.
        tbl = block["table"].copy()
        if "Overlapping Genes" in tbl.columns and not tbl.empty:
            tbl["Overlapping Genes"] = tbl["Overlapping Genes"].map(
                lambda cell: ", ".join(
                    str(protein_gene.get(x.strip(), x.strip()))
                    for x in str(cell).split(",") if x.strip()
                )
            )
        return {**block, "table": tbl}

    results: dict = {}
    for mod in modules:
        ids = set(mod["string_ids"])
        if not ids:
            results[mod["module_id"]] = None
            continue
        candidates: dict = {}
        best_block, best_key, best_min_fdr = None, None, None
        n_tested = 0
        for key, term_map in term_maps.items():
            block = genomics_ora._run_fisher_for_namespace(
                term_map, ids, universe_by_ns[key], include_bonferroni=True,
                min_term_size=min_term_size, max_term_size=max_term_size,
            )
            if block["table"].empty:
                candidates[key] = None
                continue
            block = _translate_overlap(block)
            n_tested += 1
            candidates[key] = {"block": block}
            block_min_fdr = float(block["table"]["FDR"].min())
            if best_block is None or block_min_fdr < best_min_fdr:
                best_block, best_key, best_min_fdr = block, key, block_min_fdr
        if best_block is not None:
            # Picking the smallest FDR out of n_tested independent namespaces
            # is itself a selection - an uncorrected "best of N" is
            # optimistically biased low. Bonferroni-style correction for that
            # specific selection (not a claim about the terms WITHIN a
            # namespace, which BH inside _run_fisher_for_namespace already
            # handles).
            results[mod["module_id"]] = {
                "best": {
                    "namespace": best_key,
                    "block": best_block,
                    "top_fdr_corrected": min(best_min_fdr * n_tested, 1.0),
                },
                "candidates": candidates,
                "n_namespaces_compared": n_tested,
            }
        else:
            results[mod["module_id"]] = None

    return {
        "results": results,
        "universe_size": graph["n_proteins"],
        "databases": databases,
        "caveat": caveat,
        "selection_caveat": (
            "Modules were selected for being densely connected and carrying propagated signal, and are "
            "only then tested for enrichment. These labels describe what a module contains; they are not "
            "hypothesis tests. 'Top enriched FDR' additionally corrects for comparing multiple annotation "
            "namespaces/collections per module (min(raw FDR x namespaces compared, 1.0)); the "
            "uncorrected, within-namespace value is also shown for reference."
        ),
    }


# ---------------------------------------------------------------------------
# Orchestrator - a thin wrapper around Transcriptomics_network.run_network_analysis
# ---------------------------------------------------------------------------

def _ppi_mapper(resolved_df: pd.DataFrame, graph: dict, taxon_id: str, verbose: bool = False) -> dict:
    aliases_path = gnet._string_file_paths(taxon_id).get("aliases")
    index = build_accession_string_index(graph, aliases_path=aliases_path, verbose=verbose)
    return map_protein_table_to_string(resolved_df, graph, index=index, taxon_id=taxon_id, verbose=verbose)


def run_ppi_network_analysis(
    resolved_df: pd.DataFrame,
    taxon_id: str,
    taxon_name: str = "",
    confidence: int = gnet.DEFAULT_CONFIDENCE,
    padj_cutoff: float = pm.PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = pm.LFC_CUTOFF_DEFAULT,
    direction: str = "all",
    weighting: str = txnet.NETWORK_SEED_WEIGHT_DEFAULT,
    metric: str = tx.GSEA_RANKING_METRIC_DEFAULT,
    restart_prob: float = gnet.DEFAULT_RESTART_PROB,
    n_permutations: int = txnet.NETWORK_N_PERMUTATIONS_DEFAULT,
    fdr_alpha: float = gnet.DEFAULT_FDR_ALPHA,
    module_detection: bool = True,
    min_module_size: int = gnet.DEFAULT_MIN_MODULE_SIZE,
    mcl_inflation: float = gnet.DEFAULT_MCL_INFLATION,
    max_hits_for_modules: int = txnet.DEFAULT_MAX_HITS_FOR_MODULES,
    enrichment_databases: "list[str] | None" = None,
    include_iea: bool = True,
    seed: int = 42,
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Full PPI Network / Module Detection pipeline for Proteomics: graph ->
    accession identifier bridge -> significant-protein seeds -> propagation
    -> modules -> STRING-space module labels.

    A thin wrapper around Transcriptomics_network.run_network_analysis,
    supplying this module's accession-space hooks (mapper, deg_selector,
    module_enricher, entity_noun="protein") - see that function's docstring
    for the full pipeline, the hook contract, and every parameter below not
    re-documented here. No msigdb_collections parameter: MSigDB is not in
    PPI_ENRICHMENT_DATABASES (see that constant's comment for why).

    Returns the SAME payload shape run_network_analysis does (deliberately -
    the existing Genomics/Transcriptomics client figures work unmodified),
    with "accession" in place of "gene_symbol" wherever the identifier space
    differs (the mapping table's own column - see MAPPING_COLUMNS).
    """
    return txnet.run_network_analysis(
        resolved_df, taxon_id, taxon_name=taxon_name,
        confidence=confidence, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff,
        direction=direction, weighting=weighting, metric=metric,
        restart_prob=restart_prob, n_permutations=n_permutations, fdr_alpha=fdr_alpha,
        module_detection=module_detection, min_module_size=min_module_size,
        mcl_inflation=mcl_inflation, max_hits_for_modules=max_hits_for_modules,
        enrichment_databases=enrichment_databases, msigdb_collections=None,
        include_iea=include_iea, seed=seed, force_refresh=force_refresh,
        mapper=_ppi_mapper, deg_selector=pm.select_significant_proteins,
        module_enricher=enrich_modules_by_accession, entity_noun="protein",
        progress_callback=progress_callback, verbose=verbose,
    )
