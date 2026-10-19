"""Multi-omics - Network-based Integration (fourth tile of the Multi-omics submodule).

QUESTIONS ANSWERED
------------------
Two networks over the SAME uploads as Joint Pathway Enrichment and Cross-omics Pathway Mapping:

  1. THE STRING NETWORK (gene / protein level). Which proteins sit at the centre of the COMBINED response of all your
     omics - including proteins nobody measured ("recruited") - and which connected modules do the hits form, and which
     omics feed each module? Seeds from every omics are placed on the STRING protein-protein network together; a random
     walk with restart (Koehler et al. 2008, Am J Hum Genet 82:949) spreads them; a degree-preserving permutation null
     (the numerics of Genomics_network, unchanged) gives each protein a Z-score, p-value and FDR; MCL clustering (van
     Dongen 2000) groups the hits into modules, labelled with GO / Reactome / Complex Portal terms.
  2. THE REACTOME PATHWAY-OVERLAP NETWORK (pathway level). Which pathways hit by your omics overlap in the members that were
     hit, so they form groups of related pathways? Nodes are pathways, an edge joins two whose hit members overlap (Jaccard),
     the idea of the Enrichment Map (Merico et al. 2010, PLoS ONE 5:e13984).

THE THREE DECISIONS THE USER MADE (see info/README-MultiOmics-Network.txt)
---------------------------------------------------------------------------
  - BOTH networks are built.
  - METABOLITES take part THROUGH THEIR ENZYMES: STRING is a protein network, so a hit metabolite seeds the proteins
    Reactome lists in the reactions it takes part in. Ubiquitous compounds (ATP, H2O, NAD+, ...) sit in hundreds of
    reactions and would seed a large part of the proteome, so a compound in more than max_reactions reactions is NOT used
    and is listed with the reason; a reaction that lists more than max_reaction_proteins proteins (a whole receptor or E3 family)
    does not seed either. Reactome's files do not mark the catalyst: a listed protein takes part in the reaction
    (enzyme, complex member, regulator, substrate), so this is "proteins in the reactions of your metabolites", an
    association.
  - EQUAL WEIGHT PER OMICS: each omics' seeds are scaled to sum to 1 and averaged over the omics that have seeds, so an
    omics with 3,000 hits does not drown one with 40, and a protein hit by several omics gets a larger seed. With ONE omics
    every seed weighs exactly 1.0, i.e. the result equals that omics' own network tile (tests/test_multiomics_network.py
    pins this).

WHAT IS REUSED, NOT COPIED
--------------------------
The graph, the identifier bridges (UniProt accession / gene symbol -> STRING id), the DEG rule, the propagation with its
permutation null and FDR, the MCL read-out and the module labelling are the existing functions of Genomics_network,
Transcriptomics_network, Proteomics_network and Proteomics_module, imported and called as they are. What is new here is only
what a single omics cannot do: build seeds from several omics, combine them, report who contributed, and the pathway network.
The hit-table assembly after the propagation core mirrors Transcriptomics_network.run_network_analysis (same columns, same
recruited-family FDR, same circularity caveat) because that function cannot take several seed sources.

Inherited unchanged: a seed is compared with permutations in which it is usually NOT a seed, so a seed's own significance is
largely expected (the recruited family is the discovery); STRING's combined score folds in text-mining and cross-organism
transfer and favours well-studied proteins; module labelling is post-selection description, not a hypothesis test.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

import Genomics_network as gnet
import Proteomics_module as pm
import Proteomics_network as pnet
import Transcriptomics_module as tx
import Transcriptomics_network as txnet

OMICS_ORDER = ("genomics", "transcriptomics", "proteomics", "metabolomics")
OMICS_LABELS = {"genomics": "Genomics", "transcriptomics": "Transcriptomics", "proteomics": "Proteomics",
                "metabolomics": "Metabolomics"}
DE_OMICS = ("transcriptomics", "proteomics")

CONFIG_DEFAULTS = {
    "confidence": gnet.DEFAULT_CONFIDENCE, "restart_prob": gnet.DEFAULT_RESTART_PROB,
    "n_permutations": txnet.NETWORK_N_PERMUTATIONS_DEFAULT, "fdr_alpha": gnet.DEFAULT_FDR_ALPHA,
    "mcl_inflation": gnet.DEFAULT_MCL_INFLATION, "min_module_size": gnet.DEFAULT_MIN_MODULE_SIZE,
    "max_hits_for_modules": txnet.DEFAULT_MAX_HITS_FOR_MODULES, "module_database": pnet.PPI_ENRICHMENT_DATABASE_DEFAULT,
    "include_iea": True, "max_reactions": 60, "max_reaction_proteins": 30, "pathway_top_n": 40, "pathway_min_jaccard": 0.4, "pathway_merge_jaccard": 0.7, "seed": 42,
    "de_cutoffs": {"transcriptomics": {"padj": tx.PADJ_CUTOFF_DEFAULT, "lfc": tx.LFC_CUTOFF_DEFAULT},
                   "proteomics": {"padj": pm.PADJ_CUTOFF_DEFAULT, "lfc": pm.LFC_CUTOFF_DEFAULT}},
}

HIT_COLUMNS = list(txnet.HIT_COLUMNS) + ["Seed omics", "Transcriptomics log2FC", "Transcriptomics padj",
                                           "Proteomics log2FC", "Proteomics padj"]
_UNNAMED = txnet._UNNAMED


# ---------------------------------------------------------------------------
# Combining the omics into one set of seeds
# ---------------------------------------------------------------------------

def combine_seed_weights(seed_sets: dict) -> dict:
    """Equal weight per omics. seed_sets: {omics: set of STRING ids} (empty sets are ignored).

    Each omics' seeds are scaled to sum to 1, the scaled values are averaged over the K omics that have seeds, and the result is
    multiplied by the number of distinct seeds so the mean weight is 1.0:
        weight(p) = n_seeds * (1 / K) * sum over omics o of [p in S_o] / |S_o|
    One omics alone gives 1.0 for every seed (the single-omics tile). All weights are strictly positive, which the
    propagation requires. Returns {"weights", "seed_omics": {p: [omics, in OMICS_ORDER]}, "omics", "n_seeds"}."""
    used = [o for o in OMICS_ORDER if seed_sets.get(o)]
    if not used:
        return {"weights": {}, "seed_omics": {}, "omics": [], "n_seeds": 0}
    union: set = set()
    for o in used:
        union |= set(seed_sets[o])
    k = len(used)
    weights = {p: 0.0 for p in union}
    seed_omics: dict = {p: [] for p in union}
    for o in used:
        share = 1.0 / len(seed_sets[o])
        for p in seed_sets[o]:
            weights[p] += share
            seed_omics[p].append(o)
    n = len(union)
    return {"weights": {p: w * n / k for p, w in weights.items()}, "seed_omics": seed_omics, "omics": used, "n_seeds": n}


# ---------------------------------------------------------------------------
# Seeds per omics
# ---------------------------------------------------------------------------

def genomics_seeds(final_df: pd.DataFrame, acc_index: dict) -> dict:
    """Genomics is a plain list: every usable row is a seed. UniProt ID -> STRING first, gene symbol second (the same
    two-tier rule the protein bridge uses)."""
    by_acc, by_sym = acc_index["by_accession"], acc_index["by_symbol"]
    sym_col = "Gene Symbol" if "Gene Symbol" in final_df.columns else None
    seeds: set = set()
    n_in = n_unmapped = 0
    for _, row in final_df.iterrows():
        n_in += 1
        acc = str(row.get("UniProt ID", "") or "").strip().upper()
        sym = str(row.get(sym_col, "") or "").strip().upper() if sym_col else ""
        pid = by_acc.get(acc) if acc else None
        if pid is None and sym:
            pid = by_sym.get(sym)
        if pid is None:
            n_unmapped += 1
        else:
            seeds.add(pid)
    return {"seeds": seeds, "n_input": n_in, "n_in_list": n_in, "n_unmapped": n_unmapped, "warning": None}


def de_seeds(omics: str, resolved_df: pd.DataFrame, graph: dict, indexes: dict, cutoffs: dict, taxon_id: str) -> dict:
    """Transcriptomics / Proteomics: the DEGs at the cutoffs (the SAME selection function the omics' own ORA and network tiles
    use, up and down pooled) that map to a STRING protein."""
    padj, lfc = float(cutoffs["padj"]), float(cutoffs["lfc"])
    if omics == "transcriptomics":
        mapping = txnet.map_de_table_to_string(resolved_df, graph, id_index=indexes["symbol"], taxon_id=taxon_id)
        selector = tx.select_degs
    else:
        mapping = pnet.map_protein_table_to_string(resolved_df, graph, index=indexes["accession"], taxon_id=taxon_id)
        selector = pm.select_significant_proteins
    seeds = txnet.build_seed_weights(resolved_df, mapping, padj_cutoff=padj, lfc_cutoff=lfc, direction="all",
                                     weighting="equal", deg_selector=selector)
    return {"seeds": set(seeds["weights"]), "n_input": mapping["n_input"], "n_in_list": seeds["n_degs"],
            "n_unmapped": seeds["n_degs_not_in_string"], "warning": mapping["low_mapping_warning"],
            "mapping": mapping, "n_mapped_rows": mapping["n_mapped"]}


def metabolite_enzymes(features: list, chebi_reactions: dict, reaction_proteins: dict, skeletons: dict,
                       max_reactions: int, max_reaction_proteins: int = 30) -> dict:
    """Hit metabolites -> the proteins in the Reactome reactions they take part in.

    features: Multiomics_mapping.compound_features rows ({"label", "chebi": set, "skeleton"}).
    A compound matches on its ChEBI number first; Reactome often cites another protonation state of the same compound
    (citrate(3-) CHEBI:16947 vs citric acid CHEBI:30769), so the file's ChEBI ids that share the compound's InChIKey SKELETON
    are matched too (skeletons: {ChEBI number: skeleton}). A compound whose matched ChEBI ids are in MORE than max_reactions
    reactions altogether is a hub (ATP, H2O ...) and is not used.
    A second guard works on the REACTION: Reactome lists every receptor of a family in a ligand-binding step ("Liganded Gi-
    activating GPCRs bind Gi" lists 212 proteins, the E3 ubiquitin steps list 254), so a metabolite that is a ligand would
    seed hundreds of unrelated proteins from one step. A reaction listing more than max_reaction_proteins proteins is a family
    reaction and does not seed (it is counted per compound). Measured on human data: 959 of 15,950 reactions are above 30.
    Returns {"rows": one dict per compound, "proteins": {accession: [compound labels]}}."""
    skel_to_chebi: dict = {}
    for c, s in (skeletons or {}).items():
        if c in chebi_reactions and s:
            skel_to_chebi.setdefault(s, set()).add(c)
    rows, proteins = [], {}
    for f in features:
        label = f.get("label") or ""
        matched, via = set(), ""
        for c in f.get("chebi") or ():
            if c in chebi_reactions:
                matched.add(c)
                via = "ChEBI ID"
        skel = f.get("skeleton") or ""
        if skel and skel in skel_to_chebi:
            extra = skel_to_chebi[skel] - matched
            if extra:
                matched |= extra
                via = via or "InChIKey skeleton (another protonation state)"
        rxns: set = set()
        for c in matched:
            rxns |= set(chebi_reactions[c])
        accs: set = set()
        n_generic = 0
        for r in rxns:
            plist = reaction_proteins.get(r, ())
            if len(plist) > max_reaction_proteins:
                n_generic += 1
                continue
            accs |= set(plist)
        if not matched:
            status = "left out: not a participant of any Reactome reaction for this organism"
        elif len(rxns) > max_reactions:
            status = f"left out: in {len(rxns)} reactions (more than the limit of {max_reactions}) - a hub such as a cofactor"
        elif not accs:
            status = ("left out: its reactions list no protein" if not n_generic else
                      f"left out: its reactions list no protein, or only large protein families (more than {max_reaction_proteins})")
        else:
            status = "used"
        used = status == "used"
        rows.append({"Compound": label, "Matched ChEBI": ", ".join(sorted(matched)), "Matched via": via,
                     "Reactions": len(rxns), "Family reactions skipped": n_generic, "Proteins": len(accs), "Status": status})
        if used:
            for a in accs:
                proteins.setdefault(a, []).append(label)
    return {"rows": rows, "proteins": proteins}


def metabolomics_seeds(features: list, enzyme_maps: dict, skeletons: dict, acc_index: dict, max_reactions: int,
                       max_reaction_proteins: int = 30) -> dict:
    res = metabolite_enzymes(features, enzyme_maps["chebi_reactions"], enzyme_maps["reaction_proteins"], skeletons,
                             max_reactions, max_reaction_proteins)
    by_acc = acc_index["by_accession"]
    seeds: set = set()
    n_unmapped = 0
    via_compound: dict = {}
    for acc, comps in res["proteins"].items():
        pid = by_acc.get(acc.upper())
        if pid is None:
            n_unmapped += 1
        else:
            seeds.add(pid)
            via_compound.setdefault(pid, set()).update(comps)
    n_used = sum(1 for r in res["rows"] if r["Status"] == "used")
    warn = None
    if features and n_used < len(features):
        warn = (f"{len(features) - n_used} of {len(features)} metabolites could not seed the network (see the metabolite "
                f"table: hubs, no Reactome reaction, no protein).")
    return {"seeds": seeds, "n_input": len(features), "n_in_list": n_used, "n_unmapped": n_unmapped, "warning": warn,
            "table": res["rows"], "via_compound": {p: sorted(c) for p, c in via_compound.items()}}


# ---------------------------------------------------------------------------
# The STRING run
# ---------------------------------------------------------------------------

def run_string_network(inputs: dict, taxon_id: str, taxon_name: str, cfg: dict, enzyme_maps: "dict | None" = None,
                       skeleton_provider=None, progress_callback=None) -> dict:
    """Seeds from every omics -> combined (equal weight per omics) -> propagation -> hit table -> modules.

    inputs: {"genomics": {"final_df"}, "transcriptomics": {"resolved_df"}, "proteomics": {"resolved_df"},
             "metabolomics": {"features"}} (only the omics in the run).
    skeleton_provider: callable(set of ChEBI numbers) -> {ChEBI number: skeleton} (the metabolite bridge).
    Returns the payload shape of Transcriptomics_network.run_network_analysis (so the existing figures draw it) plus
    "omics" (per-omics coverage), "metabolites" (the metabolite -> enzyme table) and "seed_omics" columns."""
    cfg = {**CONFIG_DEFAULTS, **(cfg or {})}

    def _progress(frac, label):
        if progress_callback:
            try:
                progress_callback(frac, label)
            except Exception:
                pass

    caveats: list = []
    _progress(0.02, f"Loading STRING network for taxon {taxon_id}...")
    graph = gnet.build_string_network(taxon_id, confidence=cfg["confidence"], verbose=False)
    if graph["n_proteins"] < txnet.MIN_GRAPH_PROTEINS_FOR_NULL:
        raise ValueError(f"STRING's network for taxon {taxon_id} has only {graph['n_proteins']} proteins at confidence "
                         f"{cfg['confidence']}; the permutation null needs at least {txnet.MIN_GRAPH_PROTEINS_FOR_NULL}. "
                         "Try a lower confidence cutoff.")
    paths = gnet._string_file_paths(taxon_id)
    _progress(0.06, "Building the identifier bridges to STRING...")
    indexes = {"accession": pnet.build_accession_string_index(graph, aliases_path=paths.get("aliases")),
               "symbol": txnet.build_string_id_index(graph, aliases_path=paths.get("aliases"))}

    per_omics: dict = {}
    notes: list = []
    omics_in = [o for o in OMICS_ORDER if o in inputs]
    for i, o in enumerate(omics_in):
        _progress(0.08 + 0.06 * i / max(len(omics_in), 1), f"Selecting {OMICS_LABELS[o]} seeds...")
        try:
            if o == "genomics":
                res = genomics_seeds(inputs[o]["final_df"], indexes["accession"])
            elif o in DE_OMICS:
                res = de_seeds(o, inputs[o]["resolved_df"], graph, indexes, cfg["de_cutoffs"][o], taxon_id)
            else:
                if not enzyme_maps or enzyme_maps.get("status") != "ok":
                    raise ValueError((enzyme_maps or {}).get("note") or "Reactome's reaction files are not available.")
                feats = inputs[o]["features"]
                chebi_needed = set(enzyme_maps["chebi_reactions"])
                _progress(0.09, "Matching metabolites to Reactome's compounds (the first run looks up each compound's structure "
                                "at ChEBI - a few minutes once, remembered afterwards)...")
                skeletons = skeleton_provider(chebi_needed) if skeleton_provider else {}
                res = metabolomics_seeds(feats, enzyme_maps, skeletons, indexes["accession"], int(cfg["max_reactions"]),
                                         int(cfg["max_reaction_proteins"]))
        except ValueError as exc:
            notes.append(f"{OMICS_LABELS[o]} was left out of the network: {exc}")
            continue
        if res["warning"]:
            caveats.append(f"{OMICS_LABELS[o]}: {res['warning']}")
        if not res["seeds"]:
            notes.append(f"{OMICS_LABELS[o]} has no seed that maps to a STRING protein, so it was left out of the network.")
        per_omics[o] = res

    seed_sets = {o: r["seeds"] for o, r in per_omics.items() if r["seeds"]}
    combo = combine_seed_weights(seed_sets)
    if len(combo["omics"]) < 2:
        raise ValueError("Network-based Integration needs seeds from at least two omics. " + " ".join(notes))
    caveats.extend(notes)
    caveats.append("Equal weight per omics: each omics' seeds are scaled to sum to 1 and averaged over the "
                   f"{len(combo['omics'])} omics with seeds ({', '.join(OMICS_LABELS[o] for o in combo['omics'])}); a protein "
                   "hit by several omics receives a larger seed.")
    if "metabolomics" in combo["omics"]:
        caveats.append("Metabolite seeds are proteins Reactome lists in the reactions of your hit metabolites (enzymes, "
                       "complex members, regulators, substrates), not measured proteins: read them as an association.")

    # The permutation null consumes its random numbers seed by seed, so the ORDER of the seeds changes the permutations. A
    # Python set has no stable order (it varies between runs), so the seeds go in graph-index order: the same input and the
    # same random seed always give the same result.
    pindex = graph["protein_index"]
    ordered_weights = dict(sorted(combo["weights"].items(), key=lambda kv: (pindex.get(kv[0], len(pindex)), kv[0])))
    identity = {pid: pid for pid in ordered_weights}
    seed_idx_score, n_unmapped_seeds = gnet._map_seeds_to_indices(identity, pindex, ordered_weights)
    if not seed_idx_score:
        raise ValueError("None of the seeds are present in the STRING network at this confidence.")

    def _core_progress(current, total, label):
        _progress(0.16 + 0.60 * (current / max(total, 1)), label)

    core = gnet._propagation_core(graph, seed_idx_score, restart_prob=cfg["restart_prob"],
                                  n_permutations=cfg["n_permutations"], seed_rng=cfg["seed"],
                                  progress_callback=_core_progress)

    _progress(0.78, "Assembling results...")
    protein_gene, index_protein = graph["protein_gene"], graph["index_protein"]
    seed_idx_score = core["seed_idx_score"]
    de_by_string: dict = {}
    for o in DE_OMICS:
        if o in per_omics and "mapping" in per_omics[o]:
            smap = per_omics[o]["mapping"]["string_by_row"]
            tbl = {}
            for _, r in inputs[o]["resolved_df"].iterrows():
                pid = smap.get(str(r.get("input_id")))
                if pid and pid not in tbl:
                    tbl[pid] = (r.get("log2FC"), r.get("padj"))
            de_by_string[o] = tbl

    def _num(v):
        return float(v) if v is not None and pd.notna(v) else None

    n_unnamed = 0
    rows = []
    for i in range(graph["n_proteins"]):
        pid = index_protein[i]
        name = str(protein_gene.get(pid, "")).strip()
        if not name:
            n_unnamed += 1
            name = _UNNAMED
        is_seed = i in seed_idx_score
        t_lfc, t_padj = de_by_string.get("transcriptomics", {}).get(pid, (None, None))
        p_lfc, p_padj = de_by_string.get("proteomics", {}).get(pid, (None, None))
        rows.append({
            "STRING ID": pid, "Gene Symbol": name, "Mapped via": "", "log2FC": None, "padj": None,
            "Propagated score": float(core["observed"][i]), "Z-score": float(core["zscores"][i]),
            "p-value": float(core["pvalues"][i]), "Empirical p-value": float(core["pvalues_empirical"][i]),
            "FDR": float(core["fdrs"][i]),
            "FDR (recruited family)": None if is_seed else float(core["fdrs_recruited_family"][i]),
            "Seed": bool(is_seed), "Seed score": float(seed_idx_score[i]) if is_seed else None,
            "Seed omics": " + ".join(OMICS_LABELS[o] for o in combo["seed_omics"].get(pid, [])) if is_seed else "",
            "Transcriptomics log2FC": _num(t_lfc), "Transcriptomics padj": _num(t_padj),
            "Proteomics log2FC": _num(p_lfc), "Proteomics padj": _num(p_padj),
        })
    table = pd.DataFrame(rows, columns=HIT_COLUMNS).sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)
    alpha = cfg["fdr_alpha"]
    hits_mask = table["FDR"] < alpha
    hit_table = table[hits_mask]
    n_hits = int(hits_mask.sum())
    n_recruited_all_family = int((hits_mask & ~table["Seed"]).sum())
    n_hits_seed = n_hits - n_recruited_all_family
    n_recruited = int(((table["FDR (recruited family)"] < alpha) & ~table["Seed"]).sum())
    if n_hits_seed:
        caveats.append(f"{n_hits_seed} of the {n_hits} significant proteins are seeds from your omics. Seeds receive signal "
                       "directly, and are compared against permutations in which they are usually not seeds, so their "
                       "significance is largely expected rather than a discovery.")
    if n_recruited_all_family:
        note = (f"{n_recruited} protein(s) that were never seeds are significant against a family containing ONLY non-seed "
                "proteins - a stricter, more honest test. Read this as the discovery.")
        if n_recruited != n_recruited_all_family:
            note += f" ({n_recruited_all_family} would pass the looser all-protein family alone.)"
        caveats.append(note)
    skew = core["null_skew_at_max"]
    skew_caveat = None
    if abs(skew) > 1.0:
        skew_caveat = (f"The permutation null is right-skewed (skewness {skew:.2f} at the top-scoring protein); p-values use a "
                       "skew-aware fit where it is reliable. Treat the RANKING as more reliable than any exact FDR.")
    if n_unnamed:
        caveats.append(f"{n_unnamed} protein(s) in this network have no gene name in STRING and are shown as '{_UNNAMED}'.")

    seed_pids = {index_protein[i]: s for i, s in seed_idx_score.items()}
    omics_rows = []
    for o in omics_in:
        r = per_omics.get(o)
        omics_rows.append({
            "Omics": OMICS_LABELS[o], "Entities uploaded": r["n_input"] if r else None,
            "In the seed list": r["n_in_list"] if r else None,
            "Seeds on STRING": len(r["seeds"]) if r else 0,
            "Not on STRING": r["n_unmapped"] if r else None,
            "Share of the combined weight": (1.0 / len(combo["omics"])) if o in combo["omics"] else 0.0,
            "Used": o in combo["omics"]})

    payload = {
        "taxon_id": str(taxon_id), "taxon_name": taxon_name, "module_detection": True,
        "parameters": {"string_version": gnet.STRING_VERSION, "confidence": int(cfg["confidence"]),
                       "restart_prob": float(cfg["restart_prob"]), "n_permutations": int(core["n_permutations"]),
                       "permutation_type": "degree-preserving", "seed_rng": int(cfg["seed"]), "fdr_alpha": float(alpha),
                       "mcl_inflation": float(cfg["mcl_inflation"]), "min_module_size": int(cfg["min_module_size"]),
                       "max_hits_for_modules": int(cfg["max_hits_for_modules"]),
                       "enrichment_databases": [cfg["module_database"]], "include_iea": bool(cfg["include_iea"]),
                       "module_max_term_size": gnet.DEFAULT_MODULE_MAX_TERM_SIZE, "max_reactions": int(cfg["max_reactions"]),
                       "max_reaction_proteins": int(cfg["max_reaction_proteins"]),
                       "de_cutoffs": cfg["de_cutoffs"], "weighting": "equal weight per omics"},
        "omics": omics_rows,
        "omics_used": combo["omics"],
        "seeds": {"n_seeds": combo["n_seeds"], "n_multi_omics_seeds": sum(1 for v in combo["seed_omics"].values() if len(v) > 1),
                  "weighting": "equal weight per omics", "n_seeds_not_in_string": n_unmapped_seeds},
        "metabolites": per_omics["metabolomics"]["table"] if "metabolomics" in per_omics else None,
        "propagation": {"table": table, "n_seeds_mapped": len(seed_idx_score), "n_seeds_not_in_string": n_unmapped_seeds,
                        "n_hits": n_hits, "n_hits_seed": n_hits_seed, "n_recruited": n_recruited,
                        "n_recruited_all_family": n_recruited_all_family, "n_proteins": graph["n_proteins"],
                        "null_diagnostics": {"n_exceed_top": core["n_exceed_top"], "null_skew_at_max": core["null_skew_at_max"],
                                             "n_gamma_tail": core["n_gamma_tail"], "n_permutations": core["n_permutations"],
                                             "skew_caveat": skew_caveat}},
        "modules": None, "caveats": caveats, "_graph": graph, "_hit_ids": hit_table["STRING ID"].tolist(),
        "_per_omics": per_omics, "_seed_omics": combo["seed_omics"],
    }

    if n_hits > 0:
        _progress(0.85, "Detecting modules...")
        modules, cap_info = txnet.detect_modules_capped(hit_table, graph, min_module_size=cfg["min_module_size"],
                                                        inflation=cfg["mcl_inflation"], max_hits=cfg["max_hits_for_modules"])
        if cap_info["giant_module_warning"]:
            caveats.append(cap_info["giant_module_warning"])
        if cap_info["capped"]:
            caveats.append(f"Module detection used the {cap_info['n_used_for_modules']} strongest of {cap_info['n_hits']} "
                           "significant proteins (ranked by Z-score), to keep clustering stable and interpretable.")
        _progress(0.90, "Labelling modules...")
        enr = pnet.enrich_modules_by_accession(
            modules, graph, taxon_id, taxon_name=taxon_name, databases=[cfg["module_database"]],
            include_iea=cfg["include_iea"], max_term_size=gnet.DEFAULT_MODULE_MAX_TERM_SIZE, verbose=False)
        if enr["caveat"]:
            caveats.append(enr["caveat"])
        module_rows, membership = [], {}
        label_cols = [f"Seeds from {OMICS_LABELS[o]}" for o in combo["omics"]]
        for mod in modules:
            pids = mod["string_ids"]
            membership[mod["module_id"]] = pids
            names = [str(protein_gene.get(p, "")).strip() or _UNNAMED for p in pids]
            seeds_in = [p for p in pids if p in seed_pids]
            e = enr["results"].get(mod["module_id"])
            top_term, top_fdr, top_fdr_raw, db_label = "", None, None, ""
            if e and e["best"]["block"] is not None and not e["best"]["block"]["table"].empty:
                top = e["best"]["block"]["table"].iloc[0]
                top_term, top_fdr_raw, db_label = top["Term Name"], float(top["FDR"]), e["best"]["namespace"]
                top_fdr = float(e["best"]["top_fdr_corrected"])
            row = {"Module ID": mod["module_id"], "Size": mod["size"], "Seeds": len(seeds_in),
                   "Recruited": mod["size"] - len(seeds_in)}
            for o, col in zip(combo["omics"], label_cols):
                row[col] = sum(1 for p in seeds_in if o in combo["seed_omics"].get(p, []))
            row["Seeds from several omics"] = sum(1 for p in seeds_in if len(combo["seed_omics"].get(p, [])) > 1)
            row.update({"Genes": ", ".join(names), "Top enriched term": top_term, "Top enriched FDR": top_fdr,
                        "Top enriched FDR (raw, uncorrected)": top_fdr_raw, "Enrichment database": db_label})
            module_rows.append(row)
        cols = (["Module ID", "Size", "Seeds", "Recruited"] + label_cols +
                ["Seeds from several omics", "Genes", "Top enriched term", "Top enriched FDR",
                 "Top enriched FDR (raw, uncorrected)", "Enrichment database"])
        payload["modules"] = {"table": pd.DataFrame(module_rows, columns=cols), "membership": membership,
                              "enrichment": enr["results"], "universe_size": enr["universe_size"], "caveat": enr["caveat"],
                              "selection_caveat": enr["selection_caveat"], "cap_info": cap_info,
                              "omics": combo["omics"]}
    _progress(1.0, "Done")
    return payload


# ---------------------------------------------------------------------------
# The Reactome pathway-overlap network
# ---------------------------------------------------------------------------

def _hit_labels(spec: dict, term_id: str, omics: str) -> list:
    """A pathway's HIT members of one omics, as display labels: gene symbols upper-cased, compounds as written."""
    info = (spec.get("term_map") or {}).get(term_id)
    if not info:
        return []
    labels = spec.get("labels") or {}
    out = []
    for m in set(info["genes"]) & set(spec.get("list_set") or ()):
        lab = str(labels.get(m, m)).strip()
        if lab:
            out.append(lab if omics == "metabolomics" else lab.upper())
    return sorted(set(out))


def pathway_overlap_network(block: dict, inputs: dict, top_n: int = 40, min_jaccard: float = 0.25, fdr_mark: float = 0.05,
                            merge_jaccard: float = 0.7) -> dict:
    """Nodes = the top_n DISTINCT pathways of a Joint (Reactome) result by combined p-value that have a hit; an edge joins two
    pathways whose pooled HIT members overlap with Jaccard >= min_jaccard (genes pooled across omics by symbol, compounds
    kept apart as 'M:name'). Groups = connected components.
    Reactome pathways sit in a hierarchy, so a parent and its sub-steps often share almost the same hits and would fill the picture with
    near-duplicates (measured on the bundled examples: 40 pathways, 274 links, 16 distinct pathways). Walking the pathways in order of
    combined p-value, one whose hit members overlap an already drawn pathway with Jaccard >= merge_jaccard is MERGED into it (the better-ranked one
    is the representative; the merged ones are listed on it) and the walk continues until top_n distinct pathways are drawn.
    Returns {"nodes", "edges", "groups", "omics", "note", "canvas", "n_merged"}."""
    table = block.get("table")
    out = {"nodes": [], "edges": [], "groups": [], "omics": list(block.get("omics_combined") or []), "note": None}
    if table is None or table.empty:
        out["note"] = ("No pathway was tested in at least two omics, so there is no pathway network. " + (block.get("note") or "")).strip()
        return out
    omics = out["omics"]
    ranked = table.sort_values(by=["Combined p-value", "Term ID"], kind="stable")
    nodes, pooled, n_merged = [], {}, 0
    for _, r in ranked.iterrows():
        if len(nodes) >= int(top_n):
            break
        tid = r["Term ID"]
        per, pool = {}, set()
        for o in omics:
            labs = _hit_labels(inputs[o], tid, o)
            per[o] = labs
            pool |= {("M:" + x if o == "metabolomics" else x) for x in labs}
        if not pool:
            continue
        twin = next((n for n in nodes if len(pool & pooled[n["id"]]) / len(pool | pooled[n["id"]]) >= merge_jaccard), None)
        if twin is not None:
            twin["merged"].append({"id": tid, "name": r["Term Name"]})
            n_merged += 1
            continue
        pooled[tid] = pool
        fdr = r.get("Combined FDR")
        nodes.append({"id": tid, "name": r["Term Name"], "combined_p": float(r["Combined p-value"]),
                      "combined_fdr": None if fdr is None or pd.isna(fdr) else float(fdr),
                      "significant": bool(pd.notna(fdr) and float(fdr) <= fdr_mark),
                      "hits": {o: len(per[o]) for o in omics}, "n_hit_members": len(pool),
                      "members": {o: per[o][:12] for o in omics}, "merged": []})
    edges = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            a, b = pooled[nodes[i]["id"]], pooled[nodes[j]["id"]]
            shared = len(a & b)
            if shared:
                jac = shared / len(a | b)
                if jac >= min_jaccard:
                    edges.append({"a": nodes[i]["id"], "b": nodes[j]["id"], "jaccard": jac, "shared": shared})
    parent = {n["id"]: n["id"] for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for e in edges:
        parent[find(e["a"])] = find(e["b"])
    comp: dict = {}
    for n in nodes:
        comp.setdefault(find(n["id"]), []).append(n["id"])
    groups = sorted(comp.values(), key=lambda g: (-len(g), g[0]))
    gid = {t: k + 1 for k, g in enumerate(groups) for t in g}
    for n in nodes:
        n["group"] = gid[n["id"]]
    out.update({"nodes": nodes, "edges": edges, "n_merged": n_merged, "merge_jaccard": merge_jaccard,
                "groups": [{"id": k + 1, "size": len(g), "terms": g} for k, g in enumerate(groups)]})
    pos, width, height = pathway_canvas_layout([n["id"] for n in nodes], [(e["a"], e["b"], e["jaccard"]) for e in edges])
    out["canvas"] = {"width": width, "height": height, "box_width": LABEL_BOX[0], "box_height": LABEL_BOX[1]}
    for n in nodes:
        n["x"], n["y"] = pos[n["id"]]
    return out


# A pathway's label (up to 3 lines of ~24 characters), its marker and its evidence chips take about this much room on the
# page, in pixels: the layout keeps these boxes from overlapping, so no two labels can cover each other.
LABEL_BOX = (170.0, 92.0)


def separate_boxes(pos: "np.ndarray", box_w: float, box_h: float, iterations: int = 800) -> "np.ndarray":
    """Push apart (in place) any two nodes whose boxes (box_w x box_h, centred on the node) overlap, along the axis where they
    overlap least, until none do. Deterministic: pairs are visited in index order."""
    n = pos.shape[0]
    for _ in range(iterations):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx, dy = pos[i, 0] - pos[j, 0], pos[i, 1] - pos[j, 1]
                ox, oy = box_w - abs(dx), box_h - abs(dy)
                if ox > 0 and oy > 0:
                    if ox <= oy:
                        s = (ox / 2.0 + 0.5) * (1.0 if dx >= 0 else -1.0)
                        pos[i, 0] += s
                        pos[j, 0] -= s
                    else:
                        s = (oy / 2.0 + 0.5) * (1.0 if dy >= 0 else -1.0)
                        pos[i, 1] += s
                        pos[j, 1] -= s
                    moved = True
        if not moved:
            break
    return pos


def pathway_canvas_layout(ids: list, edges: list) -> tuple:
    """Lay the pathway network out directly in PIXELS on a canvas the engine sizes: a force layout (connected pathways close, groups
    apart), then label boxes separated, then the canvas grown to hold every box. Returns ({id: (x, y)}, width, height); y grows downward."""
    n = len(ids)
    if n == 0:
        return {}, 1100.0, 800.0
    width = max(1100.0, math.sqrt(n) * 300.0)
    height = round(width * 0.78)
    unit = force_layout(ids, edges, min_gap=0.0)
    bw, bh = LABEL_BOX
    margin_x, margin_y = bw / 2 + 30, bh / 2 + 40
    pos = np.array([[margin_x + unit[t][0] * (width - 2 * margin_x), margin_y + unit[t][1] * (height - 2 * margin_y)] for t in ids])
    separate_boxes(pos, bw, bh)
    shift_x, shift_y = margin_x - pos[:, 0].min(), margin_y - pos[:, 1].min()
    pos[:, 0] += shift_x
    pos[:, 1] += shift_y
    width = max(width, float(pos[:, 0].max() + margin_x))
    height = max(height, float(pos[:, 1].max() + margin_y))
    return {t: (round(float(pos[i, 0]), 2), round(float(pos[i, 1]), 2)) for i, t in enumerate(ids)}, round(width), round(height)


def force_layout(ids: list, edges: list, iterations: int = 400, seed: int = 7, min_gap: float = 0.09) -> dict:
    """Deterministic Fruchterman-Reingold layout in the unit square, then a pass that pushes apart any two nodes closer
    than min_gap (so no two nodes overlap at a readable chart size). Elementwise numpy only (no BLAS).
    ids: node ids; edges: (a, b, weight). Returns {id: (x, y)}."""
    n = len(ids)
    if n == 0:
        return {}
    if n == 1:
        return {ids[0]: (0.5, 0.5)}
    rng = np.random.RandomState(seed)
    pos = rng.rand(n, 2)
    idx = {t: i for i, t in enumerate(ids)}
    ei = [(idx[a], idx[b], float(w)) for a, b, w in edges if a in idx and b in idx]
    k = math.sqrt(1.0 / n)
    temp = 0.1
    for _ in range(iterations):
        d = pos[:, None, :] - pos[None, :, :]
        dist = np.sqrt((d * d).sum(axis=2)) + 1e-9
        np.fill_diagonal(dist, 1.0)
        rep = (d / dist[:, :, None]) * ((k * k) / dist)[:, :, None]
        for i in range(n):
            rep[i, i, :] = 0.0
        disp = rep.sum(axis=1)
        for a, b, w in ei:
            delta = pos[a] - pos[b]
            dl = math.sqrt(float((delta * delta).sum())) + 1e-9
            f = (dl * dl / k) * (0.5 + w)
            disp[a] -= delta / dl * f
            disp[b] += delta / dl * f
        disp -= 0.3 * (pos - pos.mean(axis=0))           # gravity keeps unconnected nodes near the others
        length = np.sqrt((disp * disp).sum(axis=1)) + 1e-9
        pos += disp / length[:, None] * np.minimum(length, temp)[:, None]
        temp = max(temp * 0.985, 0.002)
    lo, hi = pos.min(axis=0), pos.max(axis=0)
    span = np.where(hi - lo > 1e-9, hi - lo, 1.0)
    pos = (pos - lo) / span
    for _ in range(300 if min_gap > 0 else 0):  # separate AFTER scaling, so the gap is not shrunk by a later rescale
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                delta = pos[i] - pos[j]
                dl = math.sqrt(float((delta * delta).sum()))
                if dl < min_gap:
                    if dl < 1e-9:
                        delta, dl = np.array([1.0, 0.0]), 1.0
                    push = delta / dl * (min_gap * 1.01 - dl) / 2
                    pos[i] += push
                    pos[j] -= push
                    moved = True
        if not moved:
            break
    return {t: (float(pos[idx[t], 0]), float(pos[idx[t], 1])) for t in ids}
