# -*- coding: utf-8 -*-
"""
Transcriptomics statistical analysis module for Omiscera — Network / Module
Detection.

WHAT IT ANSWERS that ORA/GSEA/TF-Regulon cannot. Those three all ask "which
pre-defined gene sets are over-represented in my list". This asks two different
questions: which genes sit at the CENTRE of the response according to the
protein-protein interaction network — including genes the user never uploaded
("recruited") — and how do those hits group into connected functional modules.

RELATION TO THE GENOMICS MODULE. The numerics are shared, not copied:
Genomics_network provides STRING parsing, the random walk, the degree-
preserving permutation null and MCL. Only ONE copy of that maths exists, so the
two tiles can never drift into giving different answers to the same question.
Everything that touches the USER'S DATA is different and lives here, because
the two inputs are not comparable:

    Genomics                       Transcriptomics
    ------------------------------ ---------------------------------------
    a curated candidate list       every gene that was TESTED
    ~25 rows                       ~3,000 rows (570 pass default cutoffs)
    UniProt accessions             gene symbol / Ensembl / Entrez
    optional unsigned score        log2FC (SIGNED) + p + padj
    every row is a seed            the table must be THRESHOLDED first
    no "tested but unchanged"      ~2,430 such genes carry real information

THREE THINGS THAT WOULD SILENTLY PRODUCE WRONG SCIENCE, and what stops them:

  1. A SIGNED metric as a seed weight. RWR seed mass is unsigned, and
     Genomics_network._map_seeds_to_indices shifts non-positive weights up by
     |min|+1. Because propagate_rwr then normalizes to sum 1, only RATIOS
     survive — so that shift silently converts "all genes" into an up-weighted
     run (a strongly-down gene can end with a fraction of a strongly-up gene's
     mass). Direction is therefore expressed by SEED SELECTION only, and every
     weight this module produces is strictly positive, so the shift branch is
     unreachable. Asserted, not merely intended.

  2. The module-enrichment UNIVERSE. Module members are drawn from the graph,
     so the graph is the sampling frame. Using the DE table's tested genes
     instead — the obvious default — would silently delete every RECRUITED gene
     from its own module's label, falsifying this module's headline claim while
     every number on screen still looked right.

  3. SCALE. Genomics seeds tens of genes; a DEG list is hundreds to thousands in
     a ~19,700-protein graph. Unconverged MCL on a large dense hit subgraph
     merges everything into one giant "module" with a single GO label —
     plausible-looking and completely wrong. Hence an explicit, reported hit cap
     and a giant-module warning.

Also inherited and surfaced rather than hidden: p-values are a NORMAL
APPROXIMATION of a right-skewed permutation null; module enrichment is
POST-SELECTION labelling (modules are chosen for connectivity and signal, then
tested) and is not a hypothesis test; and STRING's combined_score folds in
text-mining and cross-organism transfer, so it is not pure experimental
evidence and skews toward well-studied genes.
"""

import numpy as np
import pandas as pd

import Genomics_network as gnet
import Transcriptomics_module as tx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Deliberately smaller than Genomics' (500..10000): every permutation is a full
# RWR over the whole organism graph, and API.py's EXECUTOR has 4 workers, so one
# 10,000-permutation job would hold a quarter of the server's capacity for a
# long time. 500 is ample for z-scoring against a per-protein null.
NETWORK_N_PERMUTATIONS_CHOICES = (200, 500, 1000, 2000)
NETWORK_N_PERMUTATIONS_DEFAULT = 500

NETWORK_SEED_WEIGHTINGS = {
    "equal": "Equal weight (recommended)",
    "abs_metric": "Weighted by |ranking metric| (winsorized)",
}
NETWORK_SEED_WEIGHT_DEFAULT = "equal"
SEED_WEIGHT_WINSOR_QUANTILE = 0.99
SEED_WEIGHT_FLOOR = 0.05           # bounds the max/min weight ratio at exactly 20

# Mapping-coverage gates. Below ERROR the run is refused outright; between
# ERROR and WARN it runs with a prominent banner. STRING's preferred_name
# follows each organism's own nomenclature, so a non-model organism can map at
# a few percent and still "succeed" while returning noise - the zero-mapping
# check inherited from Genomics is not enough on its own.
MIN_MAPPED_FRACTION_ERROR = 0.10
LOW_MAPPED_FRACTION_WARN = 0.50

# The permutation null needs >= 2 candidates per degree bin to move a seed at
# all. Genomics_network._degree_bins uses bin_size = max(1, n // 10), so below
# ~20 proteins every bin holds exactly ONE candidate: the null becomes identical
# to the observed, every z-score is 0, and nothing can ever be significant -
# silently. Characterized in tests/test_network_engine_core.py; refused here
# rather than inherited.
MIN_GRAPH_PROTEINS_FOR_NULL = 20 * gnet.DEGREE_BIN_COUNT // 10

DEFAULT_MAX_HITS_FOR_MODULES = 600
GIANT_MODULE_WARN_FRACTION = 0.50

# Symbol-native only. Reactome is deliberately absent: Transcriptomics_module's
# Reactome path builds its Entrez bridge from the DE TABLE, so a recruited gene
# (absent from that table by definition) would be dropped entirely - exactly the
# failure mode this module exists to avoid. Stated in the UI, not silently.
NETWORK_ENRICHMENT_DATABASES = ("GO", "MSigDB")
NETWORK_ENRICHMENT_DATABASE_DEFAULT = "GO"

HIT_COLUMNS = [
    "STRING ID", "Gene Symbol", "Mapped via", "log2FC", "padj",
    "Propagated score", "Z-score", "p-value", "FDR", "Seed", "Seed score",
]
MODULE_COLUMNS = [
    "Module ID", "Size", "Seeds", "Recruited", "Up", "Down", "Mean log2FC",
    "Genes", "Top enriched term", "Top enriched FDR", "Enrichment database",
]
MAPPING_COLUMNS = ["input_id", "gene_symbol", "STRING ID", "Mapped via", "Seed", "Seed score"]

_UNNAMED = "(unnamed)"


# ---------------------------------------------------------------------------
# Identifier bridge: transcriptomics identifiers -> STRING protein ids
# ---------------------------------------------------------------------------

def build_string_id_index(graph: dict, aliases_path: "str | None" = None, verbose: bool = False) -> dict:
    """Build deterministic external-id -> STRING-id bridges for ONE graph.

    Genomics maps via STRING's UniProt_AC aliases; the transcriptomics resolved
    table has no UniProt column at all, so that bridge is unusable here. Gene
    symbol is the primary key instead - graph["protein_gene"] already holds
    STRING id -> preferred_name, so inverting it is nearly free (measured: 99.6%
    of the bundled 3,000-gene example maps this way). Ensembl and Entrez are
    read from the aliases file as fallbacks, which materially rescues organisms
    whose STRING nomenclature differs from mygene's.

    COLLISIONS ARE RESOLVED EXPLICITLY, never by dict(zip(...)) last-wins (the
    pattern used elsewhere in the STRING layer, which silently sends a seed to
    an arbitrary protein). In order:
        1. only proteins actually present in this graph are candidates;
        2. highest weighted degree wins (the better-connected protein is the
           one the network can actually propagate through);
        3. lexicographically smallest STRING id breaks exact ties, so the
           result is reproducible across runs and machines.

    NEVER cache this apart from its graph: membership is confidence-dependent,
    so an index built at one confidence maps ids onto proteins that may not
    exist at another.
    """
    protein_index = graph["protein_index"]
    protein_gene = graph["protein_gene"]
    raw_adj = graph["raw_adjacency"]
    degrees = np.asarray(raw_adj.sum(axis=0)).ravel()

    def _rank(pid):
        # Higher degree first, then lexicographically smallest id.
        return (-float(degrees[protein_index[pid]]), pid)

    def _add(bucket, key, pid, collisions):
        key = str(key).strip().upper()
        if not key or pid not in protein_index:
            return
        prev = bucket.get(key)
        if prev is None:
            bucket[key] = pid
            return
        if prev == pid:
            return
        winner = min((prev, pid), key=_rank)
        loser = pid if winner == prev else prev
        collisions.setdefault(key, []).append({"kept": winner, "dropped": loser})
        bucket[key] = winner

    by_symbol: dict = {}
    collisions: dict = {}
    for pid, name in protein_gene.items():
        _add(by_symbol, name, pid, collisions)

    by_ensembl: dict = {}
    by_entrez: dict = {}
    if aliases_path:
        try:
            aliases_df = pd.read_csv(aliases_path, sep="\t", dtype=str)
            aliases_df.columns = [c.lstrip("#").strip() for c in aliases_df.columns]
            for pid, alias, source in zip(aliases_df["string_protein_id"],
                                          aliases_df["alias"],
                                          aliases_df["source"]):
                if not isinstance(alias, str):
                    continue
                a = alias.strip()
                if a.upper().startswith("ENS") and "G" in a[:6].upper():
                    _add(by_ensembl, a, pid, collisions)
                elif a.isdigit() and isinstance(source, str) and "GeneID" in source:
                    _add(by_entrez, a, pid, collisions)
        except Exception as exc:  # a missing/odd aliases file must not kill the run
            if verbose:
                print(f"[TxNetwork] alias parsing skipped: {exc}")

    return {
        "by_symbol": by_symbol,
        "by_ensembl": by_ensembl,
        "by_entrez": by_entrez,
        "collisions": collisions,
        "n_graph_proteins": len(protein_index),
        "n_symbols": len(by_symbol),
    }


def map_de_table_to_string(resolved_df: pd.DataFrame, graph: dict, id_index: "dict | None" = None,
                           aliases_path: "str | None" = None, taxon_id: str = "",
                           verbose: bool = False) -> dict:
    """Map every row of the resolved DE table onto a STRING protein.

    Tiered, first hit wins: gene_symbol -> ensembl_id -> entrez_id. Later tiers
    are never consulted for an already-mapped row, so the reported per-tier
    counts are a genuine partition of the mapped rows.

    Raises below MIN_MAPPED_FRACTION_ERROR with an organism-specific message,
    and sets "low_mapping_warning" below LOW_MAPPED_FRACTION_WARN, because a
    5%-mapping run otherwise succeeds and returns noise.
    """
    if id_index is None:
        id_index = build_string_id_index(graph, aliases_path=aliases_path, verbose=verbose)

    by_symbol = id_index["by_symbol"]
    by_ensembl = id_index["by_ensembl"]
    by_entrez = id_index["by_entrez"]

    rows = []
    n_by_tier = {"gene_symbol": 0, "ensembl_id": 0, "entrez_id": 0}
    string_by_row: dict = {}
    for _, r in resolved_df.iterrows():
        sym = str(r.get("gene_symbol") or "").strip()
        ens = str(r.get("ensembl_id") or "").strip()
        ent = str(r.get("entrez_id") or "").strip()
        pid, via = None, ""
        if sym and sym.upper() in by_symbol:
            pid, via = by_symbol[sym.upper()], "gene_symbol"
        elif ens and ens.upper() in by_ensembl:
            pid, via = by_ensembl[ens.upper()], "ensembl_id"
        elif ent and ent.upper() in by_entrez:
            pid, via = by_entrez[ent.upper()], "entrez_id"
        if pid is not None:
            n_by_tier[via] += 1
            string_by_row[str(r.get("input_id"))] = pid
        rows.append({
            "input_id": r.get("input_id"),
            "gene_symbol": sym,
            "STRING ID": pid or "",
            "Mapped via": via,
            "Seed": False,
            "Seed score": None,
        })

    table = pd.DataFrame(rows, columns=MAPPING_COLUMNS)
    n_input = len(table)
    n_mapped = int((table["STRING ID"] != "").sum())
    fraction = (n_mapped / n_input) if n_input else 0.0

    if fraction < MIN_MAPPED_FRACTION_ERROR:
        raise ValueError(
            f"Only {n_mapped} of {n_input} genes ({fraction:.1%}) could be matched to a protein in "
            f"STRING's network for taxon {taxon_id or '?'}. STRING names proteins using each organism's "
            f"own nomenclature (often locus tags for less-studied species), which may not match the "
            f"identifiers in your table. Network analysis needs at least "
            f"{MIN_MAPPED_FRACTION_ERROR:.0%} coverage to be meaningful."
        )

    gene_by_string: dict = {}
    for _, r in table[table["STRING ID"] != ""].iterrows():
        gene_by_string.setdefault(r["STRING ID"], r["gene_symbol"])

    return {
        "table": table,
        "string_by_row": string_by_row,
        "gene_by_string": gene_by_string,
        "n_input": n_input,
        "n_mapped": n_mapped,
        "n_unmapped": n_input - n_mapped,
        "n_by_tier": n_by_tier,
        "fraction_mapped": fraction,
        "n_collisions_used": len(id_index["collisions"]),
        "collision_examples": dict(list(id_index["collisions"].items())[:5]),
        "low_mapping_warning": (
            f"Only {fraction:.1%} of your genes matched a STRING protein. Results are based on that "
            f"subset and may be unrepresentative."
        ) if fraction < LOW_MAPPED_FRACTION_WARN else None,
    }


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------

def build_seed_weights(resolved_df: pd.DataFrame, mapping: dict,
                       padj_cutoff: float = tx.PADJ_CUTOFF_DEFAULT,
                       lfc_cutoff: float = tx.LFC_CUTOFF_DEFAULT,
                       direction: str = "all",
                       weighting: str = NETWORK_SEED_WEIGHT_DEFAULT,
                       metric: str = tx.GSEA_RANKING_METRIC_DEFAULT) -> dict:
    """Choose the seeds and their weights.

    Seeds are the DEGs from tx.select_degs (the SAME function ORA uses, so the
    two modules can never disagree about what a DEG is) that also mapped to a
    STRING protein. Direction is expressed here, by selection, and never as a
    sign in the weight - see this module's docstring for why a signed weight
    silently biases the whole run.

    weighting="equal" (default): every seed weighs exactly 1.0. The padj/lfc
    cutoffs already made the evidence call; re-weighting by the same p-values
    double-counts it, and it makes the null test cleanly "which genes are
    seeds", which is the claim a thresholded DEG list actually supports.

    weighting="abs_metric": |ranking metric|, winsorized at the seed set's 99th
    percentile and floored at SEED_WEIGHT_FLOOR, giving a max/min ratio of
    exactly 20. Unbounded weights are not viable: with tx.GSEA_MIN_P = 1e-300,
    -log10(p) spans 1.3 to 300, so one gene would own the propagation and every
    "module" would be its neighbourhood.

    Every returned weight is strictly positive by construction.
    """
    degs = tx.select_degs(resolved_df, padj_cutoff, lfc_cutoff, direction)
    string_by_row = mapping["string_by_row"]

    raw: dict = {}
    n_degs = len(degs)
    n_not_in_string = 0
    for _, r in degs.iterrows():
        pid = string_by_row.get(str(r.get("input_id")))
        if pid is None:
            n_not_in_string += 1
            continue
        if weighting == "abs_metric":
            lfc = float(r.get("log2FC") or 0.0)
            p = float(r.get("pvalue") if pd.notna(r.get("pvalue")) else 1.0)
            padj = float(r.get("padj") if pd.notna(r.get("padj")) else 1.0)
            p = min(max(p, tx.GSEA_MIN_P), 1.0)
            padj = min(max(padj, tx.GSEA_MIN_P), 1.0)
            if metric == "log2fc":
                m = abs(lfc)
            elif metric == "pi_value":
                m = abs(lfc) * abs(np.log10(p))
            elif metric == "signed_log_padj":
                m = abs(np.log10(padj))
            else:  # signed_log_pvalue (default)
                m = abs(np.log10(p))
        else:
            m = 1.0
        raw[pid] = raw.get(pid, 0.0) + float(m)   # two DEGs on one protein sum once

    winsor_q = None
    n_winsorized = 0
    if weighting == "abs_metric" and raw:
        vals = np.asarray(list(raw.values()), dtype=np.float64)
        q = float(np.quantile(vals, SEED_WEIGHT_WINSOR_QUANTILE))
        winsor_q = q
        if q > 0:
            capped = {}
            for k, v in raw.items():
                if v > q:
                    n_winsorized += 1
                capped[k] = max(min(v, q) / q, SEED_WEIGHT_FLOOR)
            raw = capped
        else:
            raw = {k: 1.0 for k in raw}

    if raw:
        lo, hi = min(raw.values()), max(raw.values())
        assert lo > 0, "seed weights must be strictly positive so the engine's shift branch is unreachable"
        ratio = hi / lo
    else:
        ratio = 0.0

    return {
        "weights": raw,
        "n_degs": n_degs,
        "n_seeds": len(raw),
        "n_degs_not_in_string": n_not_in_string,
        "direction": direction,
        "weighting": weighting,
        "metric": metric if weighting == "abs_metric" else None,
        "winsor_q": winsor_q,
        "n_winsorized": n_winsorized,
        "weight_ratio": ratio,
    }


# ---------------------------------------------------------------------------
# Module detection with a scale guard
# ---------------------------------------------------------------------------

def detect_modules_capped(hit_table: pd.DataFrame, graph: dict,
                          min_module_size: int = gnet.DEFAULT_MIN_MODULE_SIZE,
                          inflation: float = gnet.DEFAULT_MCL_INFLATION,
                          max_hits: int = DEFAULT_MAX_HITS_FOR_MODULES) -> tuple:
    """gnet.detect_modules with the transcriptomics scale guard.

    mcl_cluster densifies the hit-induced subgraph (M.toarray()) and reads it
    out in O(n^2); at 10,000 hits that is an ~800 MB allocation. The quieter
    failure matters more: MCL hitting its iteration cap unconverged leaves many
    nonzeros per column, and the union-find read-out then merges everything
    transitively into ONE module of thousands of genes carrying a single GO
    label - entirely plausible-looking and entirely wrong. So the hit set is
    capped by Z-score (strongest signal first) and the cap is reported.
    """
    ranked = hit_table.sort_values(by="Z-score", ascending=False, kind="stable")
    ids = [pid for pid in ranked["STRING ID"].tolist() if pid in graph["protein_index"]]
    n_hits = len(ids)
    capped = n_hits > max_hits
    used = ids[:max_hits] if capped else ids

    modules = gnet.detect_modules(used, graph, min_module_size=min_module_size, inflation=inflation)
    n_clustered = sum(m["size"] for m in modules)
    largest_fraction = (modules[0]["size"] / n_clustered) if modules and n_clustered else 0.0
    giant = largest_fraction > GIANT_MODULE_WARN_FRACTION and len(modules) > 0

    info = {
        "n_hits": n_hits,
        "n_used_for_modules": len(used),
        "n_clustered": n_clustered,
        "capped": capped,
        "max_hits": max_hits,
        "largest_module_fraction": largest_fraction,
        "giant_module_warning": (
            f"One module contains {largest_fraction:.0%} of all clustered genes. That usually means the "
            f"clustering did not separate the hits into distinct groups, so its single functional label "
            f"should not be read as a specific finding."
        ) if giant else None,
    }
    return modules, info


# ---------------------------------------------------------------------------
# Per-module enrichment, against SYMBOL-keyed universes
# ---------------------------------------------------------------------------

def enrich_modules_by_symbol(modules: list, graph: dict, taxon_id: str,
                             database: str = NETWORK_ENRICHMENT_DATABASE_DEFAULT,
                             include_iea: bool = True,
                             msigdb_collections: "list[str] | None" = None,
                             min_term_size: int = tx.MIN_TERM_SIZE_DEFAULT,
                             max_term_size=gnet.DEFAULT_MODULE_MAX_TERM_SIZE,
                             verbose: bool = False) -> dict:
    """Label each module with its most enriched term, via the TRANSCRIPTOMICS
    ORA engine (tx._run_ora_for_sets, uppercase-gene-symbol-keyed). The Genomics
    equivalent is UniProt-keyed and cannot be reused.

    THE UNIVERSE IS THE GRAPH'S GENE SYMBOLS, not the DE table's tested genes.
    Module members are drawn from the graph, so the graph is the sampling frame.
    Using the DE table would silently drop every recruited gene from its own
    module's label - the exact failure this module must not have. This is a
    deliberate divergence from Genomics_network.enrich_modules.

    Emits tx.RESULT_COLUMNS unchanged, so every existing ORA chart builder and
    the "Overlapping Genes" contract work with no frontend change.
    """
    protein_gene = graph["protein_gene"]
    universe_symbols = {
        str(protein_gene.get(pid, "")).strip().upper()
        for pid in graph["protein_index"]
        if str(protein_gene.get(pid, "")).strip()
    }

    term_maps: dict = {}
    caveat = None
    if database == "GO":
        go_sets = tx.fetch_go_sets(taxon_id, include_iea=include_iea, verbose=verbose)
        for ns in tx.GO_NAMESPACES:
            if go_sets.get(ns):
                term_maps[ns] = go_sets[ns]
        if not term_maps:
            caveat = "GO annotations could not be fetched for this organism, so modules are unlabelled."
    elif database == "MSigDB":
        keys = tx._msigdb_selected_keys(msigdb_collections)
        for key in keys:
            if not tx.msigdb_collection_available(key, taxon_id):
                continue
            res = tx.fetch_msigdb_sets(key, taxon_id, verbose=verbose)
            if not res.get("fetch_failed") and res.get("sets"):
                term_maps[f"MSigDB: {tx.MSIGDB_COLLECTIONS[key]['label']}"] = res["sets"]
        if not term_maps:
            caveat = ("No MSigDB collection is available for this organism (MSigDB covers human and "
                      "mouse only), so modules are unlabelled.")
    else:
        raise ValueError(f"Unsupported module-enrichment database: {database!r}")

    results: dict = {}
    for mod in modules:
        symbols = {
            str(protein_gene.get(pid, "")).strip().upper()
            for pid in mod["string_ids"]
            if str(protein_gene.get(pid, "")).strip()
        }
        best = None
        for ns, term_map in term_maps.items():
            block = tx._run_ora_for_sets(term_map, universe_symbols, symbols,
                                         min_term_size, max_term_size, None)
            tbl = block.get("table")
            if tbl is None or tbl.empty:
                continue
            top_fdr = float(tbl.iloc[0]["FDR"])
            if best is None or top_fdr < best["top_fdr"]:
                best = {"namespace": ns, "block": block, "top_fdr": top_fdr}
        results[mod["module_id"]] = (
            {"namespace": best["namespace"], "block": best["block"]} if best else None
        )

    return {
        "results": results,
        "universe_size": len(universe_symbols),
        "database": database,
        "caveat": caveat,
        # Stated once, at the top level, because it applies to every module.
        "selection_caveat": (
            "Modules were selected for being densely connected and carrying propagated signal, and are "
            "only then tested for enrichment. These labels describe what a module contains; they are not "
            "hypothesis tests, and their FDR is not corrected for having compared "
            f"{max(len(term_maps), 1)} annotation source(s) per module."
        ),
    }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_network_analysis(
    resolved_df: pd.DataFrame,
    taxon_id: str,
    taxon_name: str = "",
    confidence: int = gnet.DEFAULT_CONFIDENCE,
    padj_cutoff: float = tx.PADJ_CUTOFF_DEFAULT,
    lfc_cutoff: float = tx.LFC_CUTOFF_DEFAULT,
    direction: str = "all",
    weighting: str = NETWORK_SEED_WEIGHT_DEFAULT,
    metric: str = tx.GSEA_RANKING_METRIC_DEFAULT,
    restart_prob: float = gnet.DEFAULT_RESTART_PROB,
    n_permutations: int = NETWORK_N_PERMUTATIONS_DEFAULT,
    fdr_alpha: float = gnet.DEFAULT_FDR_ALPHA,
    module_detection: bool = True,
    min_module_size: int = gnet.DEFAULT_MIN_MODULE_SIZE,
    mcl_inflation: float = gnet.DEFAULT_MCL_INFLATION,
    max_hits_for_modules: int = DEFAULT_MAX_HITS_FOR_MODULES,
    enrichment_database: str = NETWORK_ENRICHMENT_DATABASE_DEFAULT,
    msigdb_collections: "list[str] | None" = None,
    include_iea: bool = True,
    seed: int = 42,
    force_refresh: bool = False,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Full pipeline: graph -> identifier bridge -> seeds -> propagation ->
    modules -> module labels.

    The returned payload deliberately mirrors the Genomics netprop payload's
    shape for "propagation" and "modules", so the existing client figures
    (buildNetpropNetwork, buildNetpropHistogram, _netpropSignificantHits) work
    unmodified; transcriptomics-only information is added ALONGSIDE rather than
    by changing those keys.
    """
    def _progress(frac, label):
        if progress_callback:
            try:
                progress_callback(frac, label)
            except Exception:
                pass

    caveats = []

    _progress(0.02, f"Loading STRING network for taxon {taxon_id}...")
    graph = gnet.build_string_network(taxon_id, confidence=confidence,
                                      verbose=verbose, force_refresh=force_refresh)
    if graph["n_proteins"] < MIN_GRAPH_PROTEINS_FOR_NULL:
        raise ValueError(
            f"STRING's network for taxon {taxon_id} has only {graph['n_proteins']} proteins at "
            f"confidence {confidence}. The permutation null needs at least "
            f"{MIN_GRAPH_PROTEINS_FOR_NULL} to be able to move a seed at all - below that every "
            f"degree bin holds one protein, the null becomes identical to the observed, and nothing "
            f"could ever be significant. Try a lower confidence cutoff."
        )

    _progress(0.08, "Matching your genes to STRING proteins...")
    paths = gnet._string_file_paths(taxon_id)
    id_index = build_string_id_index(graph, aliases_path=paths.get("aliases"), verbose=verbose)
    mapping = map_de_table_to_string(resolved_df, graph, id_index=id_index,
                                     taxon_id=taxon_id, verbose=verbose)
    if mapping["low_mapping_warning"]:
        caveats.append(mapping["low_mapping_warning"])

    _progress(0.12, "Selecting seed genes...")
    seeds = build_seed_weights(resolved_df, mapping, padj_cutoff=padj_cutoff,
                               lfc_cutoff=lfc_cutoff, direction=direction,
                               weighting=weighting, metric=metric)
    if seeds["n_seeds"] == 0:
        raise ValueError(
            f"No differentially-expressed gene passed padj <= {padj_cutoff} and "
            f"|log2FC| >= {lfc_cutoff} AND matched a STRING protein "
            f"({seeds['n_degs']} DEG(s) found, {seeds['n_degs_not_in_string']} of them not in STRING). "
            f"Loosen the thresholds or check your identifiers."
        )

    # Seeds are already STRING ids, so the identity map is the correct bridge
    # here - _map_seeds_to_indices stays the single shared implementation.
    identity = {pid: pid for pid in seeds["weights"]}
    seed_idx_score, n_unmapped_seeds = gnet._map_seeds_to_indices(
        identity, graph["protein_index"], seeds["weights"])
    if not seed_idx_score:
        raise ValueError("None of the seed genes are present in the STRING network at this confidence.")

    def _core_progress(current, total, label):
        _progress(0.15 + 0.60 * (current / max(total, 1)), label)

    core = gnet._propagation_core(graph, seed_idx_score, restart_prob=restart_prob,
                                  n_permutations=n_permutations, seed_rng=seed,
                                  progress_callback=_core_progress)

    # ----- assemble the hit table, carrying the user's own measurements -----
    _progress(0.78, "Assembling results...")
    protein_gene = graph["protein_gene"]
    index_protein = graph["index_protein"]
    seed_idx_score = core["seed_idx_score"]

    de_by_string: dict = {}
    for _, r in resolved_df.iterrows():
        pid = mapping["string_by_row"].get(str(r.get("input_id")))
        if pid and pid not in de_by_string:
            de_by_string[pid] = (r.get("log2FC"), r.get("padj"))
    via_by_string: dict = {}
    for _, r in mapping["table"].iterrows():
        if r["STRING ID"] and r["STRING ID"] not in via_by_string:
            via_by_string[r["STRING ID"]] = r["Mapped via"]

    n_unnamed = 0
    rows = []
    for i in range(graph["n_proteins"]):
        pid = index_protein[i]
        name = str(protein_gene.get(pid, "")).strip()
        if not name:
            n_unnamed += 1
            name = _UNNAMED
        lfc, padj = de_by_string.get(pid, (None, None))
        is_seed = i in seed_idx_score
        rows.append({
            "STRING ID": pid,
            "Gene Symbol": name,
            "Mapped via": via_by_string.get(pid, ""),
            "log2FC": float(lfc) if lfc is not None and pd.notna(lfc) else None,
            "padj": float(padj) if padj is not None and pd.notna(padj) else None,
            "Propagated score": float(core["observed"][i]),
            "Z-score": float(core["zscores"][i]),
            "p-value": float(core["pvalues"][i]),
            "FDR": float(core["fdrs"][i]),
            "Seed": bool(is_seed),
            "Seed score": float(seed_idx_score[i]) if is_seed else None,
        })
    table = pd.DataFrame(rows, columns=HIT_COLUMNS)
    table = table.sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)

    hits_mask = table["FDR"] < fdr_alpha
    hit_table = table[hits_mask]
    n_hits = int(hits_mask.sum())
    n_recruited = int((hits_mask & ~table["Seed"]).sum())
    n_hits_seed = n_hits - n_recruited

    # A seed clearing the FDR threshold is largely CIRCULAR and must not be
    # presented as a finding. A seed receives restart mass directly, while its
    # null is built from permutations in which it is usually NOT a seed (here:
    # ~3% of the time), so the comparison is structurally rigged in its favour.
    # Measured on the bundled example: 536 of 563 "significant" proteins were
    # the user's own input. The genuinely new information is the RECRUITED set.
    # The Genomics module shares this property and never said so.
    if n_hits_seed:
        caveats.append(
            f"{n_hits_seed} of the {n_hits} significant proteins are genes you supplied as seeds. Seeds "
            f"receive signal directly, and are compared against permutations in which they are usually "
            f"not seeds, so their significance is largely expected rather than a discovery. The "
            f"{n_recruited} RECRUITED protein(s) - significant but never seeded - are the new information "
            f"this analysis adds."
        )
    skew = core["null_skew_at_max"]
    if abs(skew) > 1.0:
        caveats.append(
            f"The permutation null is right-skewed (skewness {skew:.2f} at the top-scoring protein), while "
            f"the p-values assume a normal null. Treat the exact p-values and FDRs as approximate; the "
            f"RANKING is far more trustworthy than the absolute significance."
        )
    if n_unnamed:
        caveats.append(
            f"{n_unnamed} protein(s) in this network have no gene name in STRING and are shown as "
            f"'{_UNNAMED}'."
        )

    # Mark seeds in the mapping table so the Mapping pane is self-contained.
    seed_pids = {index_protein[i]: s for i, s in seed_idx_score.items()}
    map_tbl = mapping["table"].copy()
    map_tbl["Seed"] = map_tbl["STRING ID"].map(lambda p: bool(p) and p in seed_pids)
    map_tbl["Seed score"] = map_tbl["STRING ID"].map(lambda p: seed_pids.get(p))

    payload = {
        "taxon_id": str(taxon_id),
        "taxon_name": taxon_name,
        "module_detection": bool(module_detection),
        "parameters": {
            "string_version": gnet.STRING_VERSION,
            "confidence": int(confidence),
            "restart_prob": float(restart_prob),
            "n_permutations": int(core["n_permutations"]),
            "permutation_type": "degree-preserving",
            "seed_rng": int(seed),
            "fdr_alpha": float(fdr_alpha),
            "padj_cutoff": float(padj_cutoff),
            "lfc_cutoff": float(lfc_cutoff),
            "direction": direction,
            "weighting": weighting,
            "ranking_metric": seeds["metric"],
            "winsor_quantile": SEED_WEIGHT_WINSOR_QUANTILE if weighting == "abs_metric" else None,
            "mcl_inflation": float(mcl_inflation),
            "min_module_size": int(min_module_size),
            "max_hits_for_modules": int(max_hits_for_modules),
            "enrichment_database": enrichment_database,
            "msigdb_collections": list(msigdb_collections or []),
            "include_iea": bool(include_iea),
            "module_max_term_size": gnet.DEFAULT_MODULE_MAX_TERM_SIZE,
        },
        "mapping": {
            "table": map_tbl,
            "n_input": mapping["n_input"],
            "n_mapped": mapping["n_mapped"],
            "n_unmapped": mapping["n_unmapped"],
            "n_by_tier": mapping["n_by_tier"],
            "fraction_mapped": mapping["fraction_mapped"],
            "n_collisions_used": mapping["n_collisions_used"],
            "low_mapping_warning": mapping["low_mapping_warning"],
        },
        "seeds": {k: seeds[k] for k in
                  ("n_degs", "n_seeds", "n_degs_not_in_string", "direction", "weighting",
                   "metric", "winsor_q", "n_winsorized", "weight_ratio")},
        "propagation": {
            "table": table,
            "n_seeds_mapped": len(seed_idx_score),
            "n_seeds_not_in_string": seeds["n_degs_not_in_string"] + n_unmapped_seeds,
            "n_hits": n_hits,
            "n_hits_seed": n_hits_seed,
            "n_recruited": n_recruited,
            "n_proteins": graph["n_proteins"],
            "null_diagnostics": {
                "n_exceed_top": core["n_exceed_top"],
                "null_skew_at_max": core["null_skew_at_max"],
                "n_permutations": core["n_permutations"],
            },
        },
        "modules": None,
        "caveats": caveats,
        "_graph": graph,          # in-process only; API strips before serializing
        "_hit_ids": hit_table["STRING ID"].tolist(),
    }

    if module_detection and n_hits > 0:
        _progress(0.85, "Detecting modules...")
        modules, cap_info = detect_modules_capped(
            hit_table, graph, min_module_size=min_module_size,
            inflation=mcl_inflation, max_hits=max_hits_for_modules)
        if cap_info["giant_module_warning"]:
            caveats.append(cap_info["giant_module_warning"])
        if cap_info["capped"]:
            caveats.append(
                f"Module detection used the {cap_info['n_used_for_modules']} strongest of "
                f"{cap_info['n_hits']} significant genes (ranked by Z-score), to keep clustering "
                f"stable and interpretable."
            )

        _progress(0.90, "Labelling modules...")
        enr = enrich_modules_by_symbol(
            modules, graph, taxon_id, database=enrichment_database, include_iea=include_iea,
            msigdb_collections=msigdb_collections, max_term_size=gnet.DEFAULT_MODULE_MAX_TERM_SIZE,
            verbose=verbose)
        if enr["caveat"]:
            caveats.append(enr["caveat"])

        lfc_by_string = {pid: v[0] for pid, v in de_by_string.items()}
        seed_id_set = set(seed_pids)
        module_rows, membership = [], {}
        for mod in modules:
            pids = mod["string_ids"]
            membership[mod["module_id"]] = pids
            names = [str(protein_gene.get(p, "")).strip() or _UNNAMED for p in pids]
            lfcs = [lfc_by_string.get(p) for p in pids]
            finite = [float(v) for v in lfcs if v is not None and pd.notna(v)]
            n_seed = sum(1 for p in pids if p in seed_id_set)
            e = enr["results"].get(mod["module_id"])
            top_term, top_fdr, db_label = "", None, ""
            if e and e["block"] is not None and not e["block"]["table"].empty:
                top = e["block"]["table"].iloc[0]
                top_term, top_fdr, db_label = top["Term Name"], float(top["FDR"]), e["namespace"]
            module_rows.append({
                "Module ID": mod["module_id"],
                "Size": mod["size"],
                "Seeds": n_seed,
                "Recruited": mod["size"] - n_seed,
                "Up": sum(1 for v in finite if v > 0),
                "Down": sum(1 for v in finite if v < 0),
                "Mean log2FC": (sum(finite) / len(finite)) if finite else None,
                "Genes": ", ".join(names),
                "Top enriched term": top_term,
                "Top enriched FDR": top_fdr,
                "Enrichment database": db_label,
            })
        payload["modules"] = {
            "table": pd.DataFrame(module_rows, columns=MODULE_COLUMNS),
            "membership": membership,
            "enrichment": enr["results"],
            "universe_size": enr["universe_size"],
            "caveat": enr["caveat"],
            "selection_caveat": enr["selection_caveat"],
            "cap_info": cap_info,
        }

    _progress(1.0, "Done")
    return payload
